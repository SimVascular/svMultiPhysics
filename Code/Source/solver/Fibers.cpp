// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "Fibers.h"

#include "ComMod.h"
#include "Core/Exception.h"
#include "FE/Common/FEException.h"
#include "Parameters.h"
#include "utils.h"

#include "eigen3/Eigen/Eigenvalues"

#include <mpi.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace {

Eigen::MatrixXd orthonormalize_fiber_basis(
    const Eigen::Ref<const Eigen::MatrixXd>& input,
    const std::string& location)
{
  const int nsd = input.rows();
  const int num_input_fibers = input.cols();

  if (num_input_fibers != nsd - 1 && num_input_fibers != nsd) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
      "An orthonormal fiber basis requires nsd - 1 or nsd directions.");
  }

  Eigen::MatrixXd basis = Eigen::MatrixXd::Zero(nsd, nsd);
  constexpr double linear_dependence_tolerance = 1.0e-12;

  // Modified Gram-Schmidt. Each direction is normalized and every direction
  // after the first is made orthogonal to the previously accepted ones.
  for (int direction = 0; direction < num_input_fibers; ++direction) {
    Eigen::VectorXd candidate = input.col(direction);
    for (int pass = 0; pass < 2; ++pass) {
      for (int previous = 0; previous < direction; ++previous) {
        candidate -= candidate.dot(basis.col(previous)) * basis.col(previous);
      }
    }

    const double norm = candidate.norm();
    if (norm <= linear_dependence_tolerance) {
      svmp::raise<svmp::ParseException>(
        "Fiber direction " + std::to_string(direction + 1) +
        " is linearly dependent " + location + ".");
    }
    basis.col(direction) = candidate / norm;
  }

  // A partial input basis is completed with a right-handed normal. A complete
  // input basis keeps its supplied orientation, including a left-handed one.
  if (num_input_fibers == nsd) {
    return basis;
  }

  if (nsd == 3) {
    const Eigen::Vector3d fiber = basis.col(0);
    const Eigen::Vector3d sheet = basis.col(1);
    basis.col(2) = fiber.cross(sheet);
  } else if (nsd == 2) {
    basis(0, 1) = -basis(1, 0);
    basis(1, 1) = basis(0, 0);
  } else if (nsd != 1) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
      "Fiber basis completion is only supported for spatial dimensions 1, 2, or 3.");
  }

  return basis;
}

// Interpolates one nodal fiber direction to a quadrature point. The dyads are
// interpolated before the principal direction is recovered, so opposite
// representations of the same fiber direction do not cancel. The reference
// direction is the nodal fiber at the node with the largest |N_a|; it is used
// only to choose the sign of the recovered principal eigenvector.
Vector<double> interpolate_dyadic(
    const Array<double>& nodal_fibers,
    const Vector<double>& shape_functions)
{
  const int nsd = nodal_fibers.nrows();
  const int num_nodes = nodal_fibers.ncols();
  if (nsd == 0 || num_nodes == 0) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
      "Nodal fibers must be a non-empty array.");
  }
  if (shape_functions.size() != num_nodes) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
      "The number of shape functions must match the number of nodal fibers.");
  }

  Eigen::MatrixXd interpolated_dyad = Eigen::MatrixXd::Zero(nsd, nsd);
  Eigen::VectorXd reference_direction = Eigen::VectorXd::Zero(nsd);
  double reference_weight = -1.0;

  // Compute the interpolated fiber direction by averaging the nodal fiber
  // directions weighted by the shape functions.
  for (int node = 0; node < num_nodes; ++node) {
    const double shape_value = shape_functions[node];
    if (utils::is_zero(shape_value)) {
      continue;
    }

    Eigen::VectorXd direction(nsd);
    for (int component = 0; component < nsd; ++component) {
      direction[component] = nodal_fibers(component, node);
    }

    interpolated_dyad.noalias() += shape_value * direction * direction.transpose();

    if (std::abs(shape_value) > reference_weight) {
      reference_direction = direction;
      reference_weight = std::abs(shape_value);
    }
  }

  if (reference_weight < 0.0) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
      "At least one shape function must be nonzero.");
  }

  // Compute the principal direction of the interpolated fiber dyad by solving
  // the eigenvalue problem for the interpolated fiber dyad.
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(interpolated_dyad);
  if (eigen_solver.info() != Eigen::Success) {
    svmp::raise<svmp::InternalErrorException>(
      "Failed to compute the principal direction of the interpolated fiber dyad.");
  }

  Eigen::VectorXd interpolated_fiber = eigen_solver.eigenvectors().col(nsd - 1);
  const double interpolated_fiber_norm = interpolated_fiber.norm();

  // Normalize the interpolated fiber direction to have unit length.
  if (!utils::is_zero(interpolated_fiber_norm, 1.0)) {
    interpolated_fiber /= interpolated_fiber_norm;
  }

  // Flip the interpolated fiber direction if it is opposite to the reference direction.
  if (interpolated_fiber.dot(reference_direction) < 0.0) {
    interpolated_fiber = -interpolated_fiber;
  }

  // Return the interpolated fiber direction as a Vector.
  Vector<double> result(nsd);
  for (int component = 0; component < nsd; ++component) {
    result[component] = interpolated_fiber[component];
  }
  return result;
}

// Scatters element-associated data from the master to every rank. On the
// master, global_data stores values_per_element contiguous values for each of
// the mesh.gnEl elements in the original element ordering. Each rank receives
// the values of its mesh.nEl partitioned elements in local_data.
void scatter_element_data(const double* global_data, double* local_data,
    const int values_per_element, const mshType& mesh, const CmMod& cm_mod,
    const cmType& cm)
{
  // Permute the element blocks into the partitioned ordering so that the
  // elements owned by each rank are contiguous.
  std::vector<double> reordered;
  if (cm.mas(cm_mod)) {
    reordered.resize(static_cast<size_t>(values_per_element) * mesh.gnEl);
    for (int e = 0; e < mesh.gnEl; ++e) {
      std::copy_n(global_data + static_cast<size_t>(e) * values_per_element,
        values_per_element,
        reordered.data() + static_cast<size_t>(mesh.otnIEN[e]) * values_per_element);
    }
  }

  // Send each rank its contiguous range of elements.
  const int num_proc = cm.np();
  Vector<int> send_counts(num_proc), displacements(num_proc);
  for (int i = 0; i < num_proc; ++i) {
    displacements[i] = mesh.eDist[i] * values_per_element;
    send_counts[i] = (mesh.eDist[i + 1] - mesh.eDist[i]) * values_per_element;
  }
  MPI_Scatterv(reordered.data(), send_counts.data(), displacements.data(),
    cm_mod::mpreal, local_data, mesh.nEl * values_per_element, cm_mod::mpreal,
    cm_mod.master, cm.com());
}

} // namespace


void Fibers::distribute(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm,
    mshType& mesh)
{
  const int nG = mesh.nG;
  const int nEl = mesh.nEl;

  // Share the fiber metadata read on the master with all ranks.
  int location = static_cast<int>(data_location);
  cm.bcast(cm_mod, &location);
  data_location = static_cast<DataLocation>(location);
  cm.bcast(cm_mod, &num_input_fibers);
  cm.bcast(cm_mod, &basis_input_is_orthogonal);

  if (!has_fibers()) {
    return;
  }

  // Share the number of raw components stored per location.
  int quadrature_components = 0;
  if (cm.mas(cm_mod)) {
    quadrature_components = quadrature_fibers.nrows();
  }
  cm.bcast(cm_mod, &quadrature_components);

  // The input fibers are not distributed: for cell and constant data they are
  // duplicated by quadrature_fibers. The master keeps point input fibers for
  // timestep-zero output.
  if (data_location != DataLocation::point || cm.slv(cm_mod)) {
    input_fibers.clear();
  }

  // Constant fibers are the same everywhere, so every rank gets a copy of the
  // single set of directions.
  if (data_location == DataLocation::constant) {
    if (cm.slv(cm_mod)) {
      quadrature_fibers.resize(quadrature_components, 1, 1);
    }
    MPI_Bcast(quadrature_fibers.data(), quadrature_components, cm_mod::mpreal,
      cm_mod.master, cm.com());
    return;
  }

  if (data_location != DataLocation::cell && data_location != DataLocation::point) {
    svmp::raise<svmp::InternalErrorException>(
      "Fiber input data has no cell, point, or constant direction association.");
  }

  // Cell fibers store one location per element and point fibers store one per
  // quadrature point. Point fibers also keep the nodal values of each element
  // so they can be interpolated with other basis functions later.
  const bool is_point = data_location == DataLocation::point;
  const int locations_per_element = is_point ? nG : 1;

  // Scatter the quadrature fibers.
  const Array3<double> global_quadrature = quadrature_fibers;
  quadrature_fibers.resize(quadrature_components, locations_per_element, nEl);
  scatter_element_data(global_quadrature.data(), quadrature_fibers.data(),
    quadrature_components * locations_per_element, mesh, cm_mod, cm);

  if (!is_point) {
    return;
  }

  // Scatter the nodal fibers of each element. They hold the same components
  // as the raw quadrature fibers.
  const Array3<double> global_nodal = element_nodal_fibers;
  element_nodal_fibers.resize(quadrature_components, mesh.eNoN, nEl);
  scatter_element_data(global_nodal.data(), element_nodal_fibers.data(),
    quadrature_components * mesh.eNoN, mesh, cm_mod, cm);
}


Array<double> Fibers::gather_fibers(
  const ComMod& com_mod, const CmMod& cm_mod, const mshType& mesh,
  const int nsd) const
{
  Array<double> gathered;

  if (data_location == DataLocation::none) {
    return gathered;
  }

  const int num_components = nsd * get_num_directions();
  auto& cm = com_mod.cm;

  if (data_location == DataLocation::point) {
    if (cm.mas(cm_mod)) {
      if (input_fibers.nrows() != nsd * num_input_fibers ||
          input_fibers.ncols() != mesh.gnNo) {
        svmp::raise<svmp::InternalErrorException>(
          "Point fiber dimensions are inconsistent with the mesh.");
      }
      gathered = input_fibers;
    }
    return gathered;
  }

  if (data_location != DataLocation::cell &&
      data_location != DataLocation::constant) {
    return gathered;
  }

  // Cell and constant fibers have one quadrature location per element, so the
  // quadrature arrays hold the output values directly.
  Array<double> local_fibers(num_components, mesh.nEl);
  for (int e = 0; e < mesh.nEl; ++e) {
    const int location = data_location == DataLocation::constant ? 0 : e;
    for (int component = 0; component < num_components; ++component) {
      local_fibers(component, e) = quadrature_fibers(component, 0, location);
    }
  }

  Vector<int> counts(cm.np());
  Vector<int> displacements(cm.np());
  for (int rank = 0; rank < cm.np(); ++rank) {
    displacements[rank] = mesh.eDist[rank] * num_components;
    counts[rank] =
      (mesh.eDist[rank + 1] - mesh.eDist[rank]) * num_components;
  }

  Array<double> partitioned_fibers;
  if (cm.mas(cm_mod)) {
    partitioned_fibers.resize(num_components, mesh.gnEl);
  }
  MPI_Gatherv(local_fibers.data(), mesh.nEl * num_components,
    cm_mod::mpreal, partitioned_fibers.data(), counts.data(),
    displacements.data(), cm_mod::mpreal, cm_mod.master, cm.com());

  if (cm.mas(cm_mod)) {
    gathered.resize(num_components, mesh.gnEl);
    for (int e = 0; e < mesh.gnEl; ++e) {
      gathered.set_col(e, partitioned_fibers.col(mesh.otnIEN[e]));
    }
  }
  return gathered;
}


void Fibers::clear()
{
  input_fibers.clear();
  quadrature_fibers.clear();
  element_nodal_fibers.clear();
  num_input_fibers = 0;
  basis_input_is_orthogonal = true;
  data_location = DataLocation::none;
}


Array<double> Fibers::zero_fibers(const int nsd)
{
  Array<double> fibers(nsd, 1);
  fibers = 0.0;
  return fibers;
}


void Fibers::read(const MeshParameters& mesh_parameters, mshType& mesh, const int nsd)
{
  clear();

  // Count the number of input modes.
  const int num_point_paths =
      mesh_parameters.fiber_direction_file_paths_points.size();
  const int num_cell_paths =
      mesh_parameters.fiber_direction_file_paths_cells.size();
  const int num_dirs = mesh_parameters.fiber_directions.size();

  // If no input modes are specified it means no fibers are needed, return early.
  const int num_input_modes =
      (num_point_paths != 0) + (num_cell_paths != 0) + (num_dirs != 0);
  if (num_input_modes == 0) {
    return;
  }

  // If multiple input modes are specified, throw an error.
  if (num_input_modes > 1) {
    svmp::raise<svmp::ParseException>(
      "Specify fiber directions with only one of Fiber_direction_constant, "
      "Fiber_direction_from_file_points, or Fiber_direction_from_file_cells.");
  }

  // Read the fiber directions from the appropriate source.
  if (num_dirs != 0) {
    read_directions(mesh_parameters, nsd);
  } else if (num_point_paths != 0) {
    read_from_file(mesh_parameters, mesh, nsd, DataLocation::point);
  } else {
    read_from_file(mesh_parameters, mesh, nsd, DataLocation::cell);
  }

  // Check that the input fiber directions are unit vectors.
  check_input_fibers_unit(nsd);
  basis_input_is_orthogonal = input_fibers_are_orthogonal(nsd);

  // Populate the quadrature fiber directions.
  populate_quadrature_fibers(mesh, nsd);
}


void Fibers::read_directions(const MeshParameters& mesh_parameters, const int nsd)
{
  const auto& directions = mesh_parameters.fiber_directions;
  const int num_dirs = directions.size();

  for (int i = 0; i < num_dirs; ++i) {
    const int num_components = directions[i].size();
    if (num_components != nsd) {
      svmp::raise<svmp::ParseException>(
        "Fiber_direction_constant " + std::to_string(i + 1) + " has " +
        std::to_string(num_components) + " components, but exactly " +
        std::to_string(nsd) + " are required for a " +
        std::to_string(nsd) + "D mesh.");
    }
  }

  num_input_fibers = num_dirs;
  data_location = DataLocation::constant;

  input_fibers.resize(nsd * num_input_fibers, 1);
  input_fibers = 0.0;
  for (int i = 0; i < num_input_fibers; ++i) {
    for (int j = 0; j < nsd; ++j) {
      input_fibers(i * nsd + j, 0) = directions[i][j];
    }
  }
}


void Fibers::read_from_file(const MeshParameters& mesh_parameters, mshType& mesh,
  const int nsd, const DataLocation expected_location)
{
  const auto& path_parameter = expected_location == DataLocation::point
      ? mesh_parameters.fiber_direction_file_paths_points
      : mesh_parameters.fiber_direction_file_paths_cells;
  const int num_paths = path_parameter.size();
  num_input_fibers = num_paths;
  data_location = expected_location;

  const auto paths = path_parameter();
  const bool expect_cells = expected_location == DataLocation::cell;
  const std::string association = expect_cells ? "cells" : "points";
  const std::string xml_tag = expect_cells
      ? "Fiber_direction_from_file_cells"
      : "Fiber_direction_from_file_points";

  for (int idx = 0; idx < num_input_fibers; ++idx) {
    const auto& file_name = paths[idx];

    if (FILE* file = fopen(file_name.c_str(), "r")) {
      fclose(file);
    } else {
      svmp::raise<svmp::FileNotFoundException>(file_name);
    }

    auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(file_name.c_str());
    reader->Update();
    auto grid = reader->GetOutput();

    if (mesh.gnEl != grid->GetNumberOfCells()) {
      svmp::raise<svmp::FileFormatException>(file_name,
        "The number of elements (" +
        std::to_string(grid->GetNumberOfCells()) +
        ") is not equal to the number of elements (" +
        std::to_string(mesh.gnEl) + ") for the mesh named '" + mesh.name + "'.");
    }

    vtkDataArray* array = nullptr;
    vtkIdType num_locations = 0;
    if (expect_cells) {
      array = grid->GetCellData()->GetArray("FIB_DIR");
      num_locations = grid->GetNumberOfCells();
    } else {
      array = grid->GetPointData()->GetArray("FIB_DIR");
      num_locations = grid->GetNumberOfPoints();
    }
    if (array == nullptr) {
      svmp::raise<svmp::FileFormatException>(file_name,
        "No 'FIB_DIR' " + association + " data found as required by " +
        xml_tag + ".");
    }
    if (array->GetNumberOfTuples() != num_locations) {
      svmp::raise<svmp::FileFormatException>(file_name,
        "The number of 'FIB_DIR' tuples (" +
        std::to_string(array->GetNumberOfTuples()) +
        ") does not match the number of associated " + association + " (" +
        std::to_string(num_locations) + ").");
    }
    if (array->GetNumberOfComponents() < nsd) {
      svmp::raise<svmp::FileFormatException>(file_name,
        "The 'FIB_DIR' data has fewer than " + std::to_string(nsd) +
        " components.");
    }

    const vtkIdType expected_num_locations = expect_cells ? mesh.gnEl : mesh.gnNo;
    if (num_locations != expected_num_locations) {
      svmp::raise<svmp::FileFormatException>(file_name,
        "The number of " + association + " (" +
        std::to_string(num_locations) + ") does not match the mesh (" +
        std::to_string(expected_num_locations) + ").");
    }

    if (idx == 0) {
      input_fibers.resize(nsd * num_input_fibers, num_locations);
      input_fibers = 0.0;
    }

    const int offset = idx * nsd;
    for (vtkIdType location = 0; location < num_locations; ++location) {
      const auto tuple = array->GetTuple(location);
      for (int component = 0; component < nsd; ++component) {
        input_fibers(offset + component, location) = tuple[component];
      }
    }
  }
}


void Fibers::populate_quadrature_fibers(const mshType& mesh, const int nsd)
{
  const int num_components = nsd * num_input_fibers;

  if (data_location == DataLocation::constant) {
    quadrature_fibers.resize(num_components, 1, 1);
    for (int component = 0; component < num_components; ++component) {
      quadrature_fibers(component, 0, 0) = input_fibers(component, 0);
    }
    return;
  }

  if (data_location == DataLocation::cell) {
    quadrature_fibers.resize(num_components, 1, mesh.gnEl);
    for (int e = 0; e < mesh.gnEl; ++e) {
      for (int component = 0; component < num_components; ++component) {
        quadrature_fibers(component, 0, e) = input_fibers(component, e);
      }
    }
    return;
  }

  if (data_location != DataLocation::point || mesh.nG <= 0) {
    svmp::raise<svmp::InternalErrorException>(
      "Point fiber directions cannot be populated before mesh quadrature points are initialized.");
  }
  quadrature_fibers.resize(num_components, mesh.nG, mesh.gnEl);
  quadrature_fibers = 0.0;
  interpolate_point_fibers_mesh(mesh, nsd);
}


void Fibers::interpolate_point_fibers_mesh(
  const mshType& mesh, const int nsd)
{
  if (mesh.N.nrows() != mesh.eNoN || mesh.N.ncols() != mesh.nG) {
    svmp::raise<svmp::InternalErrorException>(
      "Mesh shape functions are inconsistent with its nodes and quadrature points.");
  }

  element_nodal_fibers.resize(
    nsd * num_input_fibers, mesh.eNoN, mesh.gnEl);
  element_nodal_fibers = 0.0;
  for (int e = 0; e < mesh.gnEl; ++e) {
    for (int a = 0; a < mesh.eNoN; ++a) {
      const int node = mesh.gIEN(a, e);
      for (int direction = 0; direction < num_input_fibers; ++direction) {
        for (int component = 0; component < nsd; ++component) {
          element_nodal_fibers(component + direction * nsd, a, e) =
            input_fibers(component + direction * nsd, node);
        }
      }
    }
  }

  for (int g = 0; g < mesh.nG; ++g) {
    const auto basis_functions = mesh.N.col(g);
    for (int e = 0; e < mesh.gnEl; ++e) {
      const auto fibers =
        get_quadrature_fibers(nsd, basis_functions, e);
      for (int direction = 0; direction < num_input_fibers; ++direction) {
        for (int component = 0; component < nsd; ++component) {
          quadrature_fibers(component + direction * nsd, g, e) =
            fibers(component, direction);
        }
      }
    }
  }
}


void Fibers::check_input_fibers_unit(const int nsd) const
{
  constexpr double tolerance = 1.0e-6;

  // Check that the fiber directions are unit vectors.
  for (int location = 0; location < input_fibers.ncols(); ++location) {
    for (int i = 0; i < num_input_fibers; ++i) {
      double norm_sq = 0.0;
      for (int component = 0; component < nsd; ++component) {
        const double value = input_fibers(i * nsd + component, location);
        norm_sq += value * value;
      }
      const double norm = std::sqrt(norm_sq);
      if (std::abs(norm - 1.0) > tolerance) {
        svmp::raise<svmp::ParseException>(
          "Fiber direction " + std::to_string(i + 1) + " at location " +
          std::to_string(location) + " is not a unit vector (norm = " +
          std::to_string(norm) + ").");
      }
    }
  }
}


bool Fibers::input_fibers_are_orthogonal(const int nsd) const
{
  constexpr double tolerance = 1.0e-6;

  const int nlocations = input_fibers.ncols();
  for (int location = 0; location < nlocations; ++location) {
    for (int i = 0; i < num_input_fibers; ++i) {
      for (int j = 0; j < i; ++j) {
        double dot = 0.0;
        for (int component = 0; component < nsd; ++component) {
          dot += input_fibers(i * nsd + component, location) *
                 input_fibers(j * nsd + component, location);
        }
        if (std::abs(dot) > tolerance) {
          return false;
        }
      }
    }
  }
  return true;
}


Array<double> Fibers::get_quadrature_fibers(
  const int nsd, const int g, const int e,
  const bool orthonormalize) const
{
  if (nsd <= 0) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
      "The spatial dimension must be positive when accessing fiber directions.");
  }

  double* data = nullptr;
  int num_components = 0;

  if (data_location == DataLocation::constant) {
    num_components = quadrature_fibers.nrows();
    data = quadrature_fibers.data();
  } else if (data_location == DataLocation::cell) {
    if (e < 0 || e >= quadrature_fibers.nslices()) {
      svmp::raise<svmp::IndexOutOfRangeException>(
        "The requested fiber element index " + std::to_string(e) +
        " is out of bounds.");
    }
    num_components = quadrature_fibers.nrows();
    data = quadrature_fibers.data() + e * num_components;
  } else if (data_location == DataLocation::point) {
    if (g < 0 || g >= quadrature_fibers.ncols() ||
        e < 0 || e >= quadrature_fibers.nslices()) {
      svmp::raise<svmp::IndexOutOfRangeException>(
        "The requested quadrature fiber indices (" + std::to_string(g) +
        ", " + std::to_string(e) + ") are out of bounds.");
    }
    num_components = quadrature_fibers.nrows();
    const auto offset =
      (e * quadrature_fibers.ncols() + g) * num_components;
    data = quadrature_fibers.data() + offset;
  } else {
    svmp::raise<svmp::InternalErrorException>(
      "Fiber input data has not been initialized.");
  }

  if (num_components == 0 || num_components % nsd != 0) {
    svmp::raise<svmp::InternalErrorException>(
      "The fiber data have " + std::to_string(num_components) +
      " components, which is incompatible with spatial dimension " +
      std::to_string(nsd) + ".");
  }

  if (!orthonormalize) {
    return {nsd, num_components / nsd, data};
  }
  if (num_input_fibers != nsd - 1 && num_input_fibers != nsd) {
    svmp::raise<svmp::ParseException>(
      "An orthonormal fiber basis in " + std::to_string(nsd) +
      "D requires " + std::to_string(nsd - 1) + " or " +
      std::to_string(nsd) + " input directions, but " +
      std::to_string(num_input_fibers) + " were provided.");
  }
  if (!basis_input_is_orthogonal) {
    svmp::raise<svmp::ParseException>(
      "An orthonormal fiber basis requires mutually orthogonal input directions.");
  }

  const Array<double> raw_fibers(nsd, num_components / nsd, data);
  const Eigen::Map<const Eigen::MatrixXd> input(
    raw_fibers.data(), nsd, num_input_fibers);
  const auto basis = orthonormalize_fiber_basis(
    input, "at quadrature point " + std::to_string(g) +
    " of element " + std::to_string(e));

  Array<double> result(nsd, nsd);
  Eigen::Map<Eigen::MatrixXd>(result.data(), nsd, nsd) = basis;
  return result;
}


Array<double> Fibers::get_quadrature_fibers(
  const int nsd, const Vector<double>& basis_functions, const int e,
  const bool orthonormalize) const
{
  if (data_location != DataLocation::point) {
    return get_quadrature_fibers(nsd, 0, e, orthonormalize);
  }

  const auto interpolated_fibers =
    interpolate_point_fibers_element(nsd, basis_functions, e);
  if (!orthonormalize) {
    return interpolated_fibers;
  }
  if (num_input_fibers != nsd - 1 && num_input_fibers != nsd) {
    svmp::raise<svmp::ParseException>(
      "An orthonormal fiber basis in " + std::to_string(nsd) +
      "D requires " + std::to_string(nsd - 1) + " or " +
      std::to_string(nsd) + " input directions, but " +
      std::to_string(num_input_fibers) + " were provided.");
  }
  if (!basis_input_is_orthogonal) {
    svmp::raise<svmp::ParseException>(
      "An orthonormal fiber basis requires mutually orthogonal input directions.");
  }

  const Eigen::Map<const Eigen::MatrixXd> input(
    interpolated_fibers.data(), nsd, num_input_fibers);
  const auto basis = orthonormalize_fiber_basis(
    input, "in element " + std::to_string(e));

  Array<double> result(nsd, nsd);
  Eigen::Map<Eigen::MatrixXd>(result.data(), nsd, nsd) = basis;
  return result;
}


Array<double> Fibers::interpolate_point_fibers_element(
  const int nsd, const Vector<double>& basis_functions, const int e) const
{
  if (nsd <= 0) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
      "The spatial dimension must be positive when accessing fiber directions.");
  }
  if (e < 0 || e >= element_nodal_fibers.nslices()) {
    svmp::raise<svmp::IndexOutOfRangeException>(
      "The requested fiber element index " + std::to_string(e) +
      " is out of bounds.");
  }
  const int num_basis_functions = basis_functions.size();
  if (num_basis_functions <= 0 ||
      num_basis_functions > element_nodal_fibers.ncols()) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
      "The number of basis-function values must be positive and cannot exceed "
      "the number of element nodes.");
  }

  Array<double> interpolated_fibers(nsd, num_input_fibers);
  Array<double> nodal_fibers(nsd, num_basis_functions);

  for (int iFn = 0; iFn < num_input_fibers; ++iFn) {
    for (int a = 0; a < num_basis_functions; ++a) {
      for (int component = 0; component < nsd; ++component) {
        nodal_fibers(component, a) =
          element_nodal_fibers(component + iFn * nsd, a, e);
      }
    }

    const auto interpolated =
      interpolate_dyadic(nodal_fibers, basis_functions);
    for (int component = 0; component < nsd; ++component) {
      interpolated_fibers(component, iFn) = interpolated[component];
    }
  }

  return interpolated_fibers;
}
