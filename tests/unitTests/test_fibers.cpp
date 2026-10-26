// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "ComMod.h"
#include "Core/Exception.h"
#include "Fibers.h"
#include "Parameters.h"

#include <gtest/gtest.h>

#include <vtkDoubleArray.h>
#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <cmath>
#include <filesystem>
#include <numbers>
#include <string>
#include <vector>

namespace {

constexpr auto raw_fibers = Fibers::Representation::raw;
constexpr auto fiber_basis = Fibers::Representation::orthonormal_basis;

mshType make_single_element_mesh()
{
  mshType mesh;
  mesh.gnEl = 1;
  mesh.gnNo = 1;
  mesh.nG = 1;
  mesh.eNoN = 1;
  mesh.name = "unit-test-mesh";
  mesh.N.resize(1, 1);
  mesh.N(0, 0) = 1.0;
  mesh.gIEN.resize(1, 1);
  mesh.gIEN(0, 0) = 0;
  return mesh;
}

mshType make_two_node_line_mesh(const int num_quadrature_points = 1)
{
  mshType mesh;
  mesh.name = "two-node-line";
  mesh.gnEl = 1;
  mesh.gnNo = 2;
  mesh.eNoN = 2;
  mesh.nG = num_quadrature_points;
  mesh.gIEN.resize(mesh.eNoN, mesh.gnEl);
  mesh.gIEN(0, 0) = 0;
  mesh.gIEN(1, 0) = 1;
  mesh.N.resize(mesh.eNoN, mesh.nG);
  for (int g = 0; g < mesh.nG; ++g) {
    mesh.N(0, g) = 0.5;
    mesh.N(1, g) = 0.5;
  }
  return mesh;
}

void write_nodal_fiber_file(
    const std::string& file_name,
    const std::vector<std::vector<double>>& directions)
{
  auto points = vtkSmartPointer<vtkPoints>::New();
  auto fibers = vtkSmartPointer<vtkDoubleArray>::New();
  fibers->SetName("FIB_DIR");
  fibers->SetNumberOfComponents(3);
  for (int node = 0; node < static_cast<int>(directions.size()); ++node) {
    points->InsertNextPoint(static_cast<double>(node), 0.0, 0.0);
    const auto& direction = directions[node];
    fibers->InsertNextTuple3(
      direction[0], direction[1], direction.size() > 2 ? direction[2] : 0.0);
  }

  auto line = vtkSmartPointer<vtkLine>::New();
  line->GetPointIds()->SetId(0, 0);
  line->GetPointIds()->SetId(1, 1);

  auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  grid->SetPoints(points);
  grid->InsertNextCell(line->GetCellType(), line->GetPointIds());
  grid->GetPointData()->AddArray(fibers);

  auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName(file_name.c_str());
  writer->SetInputData(grid);
  ASSERT_EQ(writer->Write(), 1);
}

} // namespace

TEST(Fibers, LoadsRawThenOrthonormalizesConstantDirections)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{1.0, 0.0, 0.0}, false);
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{0.0, 1.0, 0.0}, false);

  auto mesh = make_single_element_mesh();
  Fibers fibers;

  EXPECT_FALSE(fibers.has_fibers());
  EXPECT_EQ(fibers.get_num_input_directions(), 0);
  {
    const auto zeros = Fibers::zero_fibers(3);
    ASSERT_EQ(zeros.nrows(), 3);
    ASSERT_EQ(zeros.ncols(), 1);
    EXPECT_DOUBLE_EQ(zeros(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(zeros(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(zeros(2, 0), 0.0);

    const auto eigen_zeros = Fibers::zero_eigen_fibers(3);
    ASSERT_EQ(eigen_zeros.rows(), 3);
    ASSERT_EQ(eigen_zeros.cols(), 1);
    EXPECT_DOUBLE_EQ(eigen_zeros(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(eigen_zeros(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(eigen_zeros(2, 0), 0.0);
  }
  fibers.read(parameters, mesh, 3);

  ASSERT_TRUE(fibers.has_fibers());
  ASSERT_EQ(fibers.get_num_input_directions(), 2);
  ASSERT_EQ(fibers.get_num_basis_directions(), 0);
  EXPECT_FALSE(fibers.has_orthonormal_basis());
  EXPECT_EQ(fibers.get_default_representation(),
    Fibers::Representation::raw);
  EXPECT_EQ(fibers.get_num_directions(Fibers::Representation::raw), 2);

  auto quadrature_fibers = fibers.get_quadrature_fibers(3, 0, 0, raw_fibers);
  ASSERT_EQ(quadrature_fibers.nrows(), 3);
  ASSERT_EQ(quadrature_fibers.ncols(), 2);
  EXPECT_DOUBLE_EQ(quadrature_fibers(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(2, 0), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(0, 1), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(2, 1), 0.0);

  fibers.initialize_fiber_bases(3);
  EXPECT_TRUE(fibers.has_orthonormal_basis());
  ASSERT_EQ(fibers.get_num_input_directions(), 2);
  ASSERT_EQ(fibers.get_num_basis_directions(), 3);
  EXPECT_EQ(fibers.get_default_representation(),
    Fibers::Representation::orthonormal_basis);
  EXPECT_EQ(fibers.get_num_directions(
    Fibers::Representation::orthonormal_basis), 3);
  const auto raw_after = fibers.get_quadrature_fibers(3, 0, 0, raw_fibers);
  ASSERT_EQ(raw_after.ncols(), 2);
  EXPECT_DOUBLE_EQ(raw_after(0, 0), 1.0);
  const auto basis = fibers.get_quadrature_fibers(3, 0, 0, fiber_basis);
  const auto selected = fibers.get_quadrature_fibers(
    3, 0, 0, fibers.get_default_representation());
  const auto* basis_data = basis.data();
  ASSERT_EQ(basis.ncols(), 3);
  EXPECT_EQ(selected.data(), basis.data());
  EXPECT_DOUBLE_EQ(basis(2, 2), 1.0);
  fibers.initialize_fiber_bases(3);
  EXPECT_EQ(fibers.get_quadrature_fibers(3, 0, 0, fiber_basis).data(), basis_data);

  fibers.clear();
  EXPECT_FALSE(fibers.has_fibers());
  EXPECT_EQ(fibers.get_num_input_directions(), 0);
}

TEST(Fibers, ConstantDirectionsShareStorageAcrossQueries)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{1.0, 0.0, 0.0}, false);
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{0.0, 1.0, 0.0}, false);

  auto mesh = make_single_element_mesh();
  mesh.gnEl = 2;
  mesh.nG = 2;

  Fibers fibers;
  fibers.read(parameters, mesh, 3);
  fibers.initialize_fiber_bases(3);

  const auto reference = fibers.get_quadrature_fibers(3, 0, 0, fiber_basis);
  const auto eigen_fibers = fibers.get_eigen_fibers(3, 1, 1, fiber_basis);
  EXPECT_EQ(eigen_fibers.data(), reference.data());
  EXPECT_EQ(fibers.get_quadrature_fibers(3, 1, 1, fiber_basis).data(), reference.data());
  EXPECT_DOUBLE_EQ(eigen_fibers(0, 2), 0.0);
  EXPECT_DOUBLE_EQ(eigen_fibers(1, 2), 0.0);
  EXPECT_DOUBLE_EQ(eigen_fibers(2, 2), 1.0);
}

TEST(Fibers, AcceptsRawCountRejectedByBasisInitialization)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{1.0, 0.0, 0.0}, false);

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  EXPECT_NO_THROW(fibers.read(parameters, mesh, 3));
  EXPECT_EQ(fibers.get_num_input_directions(), 1);
  EXPECT_THROW(fibers.initialize_fiber_bases(3), svmp::ParseException);
}

TEST(Fibers, AcceptsNoDirections)
{
  MeshParameters parameters;
  auto mesh = make_single_element_mesh();
  Fibers fibers;

  EXPECT_NO_THROW(fibers.read(parameters, mesh, 3));
  EXPECT_FALSE(fibers.has_fibers());
  EXPECT_EQ(fibers.get_num_input_directions(), 0);
  EXPECT_EQ(fibers.get_num_basis_directions(), 0);
}

TEST(Fibers, RejectsConstantDirectionWithWrongComponentCount)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{1.0, 0.0}, false);
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{0.0, 1.0, 0.0}, false);

  auto mesh = make_single_element_mesh();
  Fibers fibers;

  try {
    fibers.read(parameters, mesh, 3);
    FAIL() << "Expected a malformed constant fiber direction to be rejected.";
  } catch (const svmp::ParseException& error) {
    EXPECT_EQ(error.message(),
      "Fiber_direction_constant 1 has 2 components, but exactly 3 are "
      "required for a 3D mesh.");
  }
}

TEST(Fibers, AcceptsArbitraryRawDirectionCount)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{1.0, 0.0, 0.0}, false);
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{0.0, 1.0, 0.0}, false);
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{0.0, 0.0, 1.0}, false);
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{-1.0, 0.0, 0.0}, false);

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  EXPECT_NO_THROW(fibers.read(parameters, mesh, 3));
  EXPECT_EQ(fibers.get_num_input_directions(), 4);
  EXPECT_EQ(fibers.get_quadrature_fibers(3, 0, 0, raw_fibers).ncols(), 4);
  EXPECT_THROW(fibers.initialize_fiber_bases(3), svmp::ParseException);
}

TEST(Fibers, PreservesCompleteLeftHandedBasis)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{1.0, 0.0, 0.0}, false);
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{0.0, 1.0, 0.0}, false);
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{0.0, 0.0, -1.0}, false);

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  fibers.read(parameters, mesh, 3);
  fibers.initialize_fiber_bases(3);

  const auto basis = fibers.get_quadrature_fibers(3, 0, 0, fiber_basis);
  ASSERT_EQ(basis.ncols(), 3);
  EXPECT_DOUBLE_EQ(basis(2, 2), -1.0);
  EXPECT_EQ(fibers.get_quadrature_fibers(3, 0, 0, raw_fibers).ncols(), 3);
}

TEST(Fibers, RejectsNonUnitConstantDirection)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{3.0, 4.0, 0.0}, false);
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{0.0, 1.0, 0.0}, false);

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  EXPECT_THROW(fibers.read(parameters, mesh, 3), svmp::ParseException);
}

TEST(Fibers, RejectsZeroLengthConstantDirection)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{0.0, 0.0, 0.0}, false);
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{0.0, 1.0, 0.0}, false);

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  EXPECT_THROW(fibers.read(parameters, mesh, 3), svmp::ParseException);
}

TEST(Fibers, RejectsNonOrthogonalConstantDirections)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{1.0, 0.0, 0.0}, false);
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{0.6, 0.8, 0.0}, false);

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  EXPECT_NO_THROW(fibers.read(parameters, mesh, 3));
  const auto raw = fibers.get_quadrature_fibers(3, 0, 0, raw_fibers);
  EXPECT_DOUBLE_EQ(raw(0, 1), 0.6);
  EXPECT_DOUBLE_EQ(raw(1, 1), 0.8);
  EXPECT_THROW(fibers.initialize_fiber_bases(3), svmp::ParseException);
}

TEST(Fibers, RejectsAntiparallelConstantDirections)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{1.0, 0.0, 0.0}, false);
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{-1.0, 0.0, 0.0}, false);

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  EXPECT_NO_THROW(fibers.read(parameters, mesh, 3));
  EXPECT_THROW(fibers.initialize_fiber_bases(3), svmp::ParseException);
}

TEST(Fibers, LoadsCellDirectionsFromVtu)
{
  MeshParameters parameters;
  parameters.fiber_direction_file_paths_cells.set(
      "Fiber_direction_from_file_cells", false,
      {std::string(UNIT_TEST_DATA_DIR) + "/fibers_cell_data.vtu",
       std::string(UNIT_TEST_DATA_DIR) + "/fibers_cell_data_sheet.vtu"});

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  fibers.read(parameters, mesh, 3);

  ASSERT_TRUE(fibers.has_fibers());
  ASSERT_EQ(fibers.get_num_input_directions(), 2);
  ASSERT_EQ(fibers.get_num_basis_directions(), 0);

  const auto quadrature_fibers = fibers.get_quadrature_fibers(3, 0, 0, raw_fibers);
  ASSERT_EQ(quadrature_fibers.nrows(), 3);
  ASSERT_EQ(quadrature_fibers.ncols(), 2);
  EXPECT_DOUBLE_EQ(quadrature_fibers(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(2, 0), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(0, 1), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(2, 1), 0.0);
}

TEST(Fibers, CellDirectionsShareStorageAcrossQuadraturePoints)
{
  MeshParameters parameters;
  parameters.fiber_direction_file_paths_cells.set(
      "Fiber_direction_from_file_cells", false,
      {std::string(UNIT_TEST_DATA_DIR) + "/fibers_cell_data.vtu",
       std::string(UNIT_TEST_DATA_DIR) + "/fibers_cell_data_sheet.vtu"});

  auto mesh = make_single_element_mesh();
  mesh.nG = 2;
  Fibers fibers;
  fibers.read(parameters, mesh, 3);

  const auto reference = fibers.get_quadrature_fibers(3, 0, 0, raw_fibers);
  EXPECT_EQ(fibers.get_quadrature_fibers(3, 1, 0, raw_fibers).data(), reference.data());
}

TEST(Fibers, LoadsPointDirectionsFromVtu)
{
  MeshParameters parameters;
  parameters.fiber_direction_file_paths_points.set(
      "Fiber_direction_from_file_points", false,
      {std::string(UNIT_TEST_DATA_DIR) + "/fibers_point_data.vtu",
       std::string(UNIT_TEST_DATA_DIR) + "/fibers_point_data_sheet.vtu"});

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  fibers.read(parameters, mesh, 3);
  fibers.initialize_fiber_bases(3);

  ASSERT_TRUE(fibers.has_fibers());
  ASSERT_EQ(fibers.get_num_input_directions(), 2);
  ASSERT_EQ(fibers.get_num_basis_directions(), 3);

  const auto raw = fibers.get_quadrature_fibers(3, 0, 0, raw_fibers);
  ASSERT_EQ(raw.ncols(), 2);
  const auto quadrature_fibers = fibers.get_quadrature_fibers(3, 0, 0, fiber_basis);
  ASSERT_EQ(quadrature_fibers.nrows(), 3);
  ASSERT_EQ(quadrature_fibers.ncols(), 3);
  EXPECT_DOUBLE_EQ(quadrature_fibers(0, 0), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(1, 0), 1.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(2, 0), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(0, 1), 1.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(1, 1), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(2, 1), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(0, 2), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(1, 2), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(2, 2), -1.0);
}

TEST(Fibers, CompletesTwoDimensionalPointBasisWithPositiveOrientation)
{
  MeshParameters parameters;
  parameters.fiber_direction_file_paths_points.set(
      "Fiber_direction_from_file_points", false,
      {std::string(UNIT_TEST_DATA_DIR) + "/fibers_point_data.vtu"});

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  fibers.read(parameters, mesh, 2);
  fibers.initialize_fiber_bases(2);

  EXPECT_EQ(fibers.get_quadrature_fibers(2, 0, 0, raw_fibers).ncols(), 1);
  const auto quadrature_fibers = fibers.get_quadrature_fibers(2, 0, 0, fiber_basis);
  ASSERT_EQ(quadrature_fibers.nrows(), 2);
  ASSERT_EQ(quadrature_fibers.ncols(), 2);
  EXPECT_DOUBLE_EQ(quadrature_fibers(0, 0), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(1, 0), 1.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(0, 1), -1.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(1, 1), 0.0);
}

TEST(Fibers, InterpolatesNodalFibersAtWeightedQuadraturePoints)
{
  const std::string file_name =
      ::testing::TempDir() + "weighted_nodal_fibers.vtu";
  write_nodal_fiber_file(file_name, {{1.0, 0.0}, {0.0, 1.0}});

  MeshParameters parameters;
  parameters.fiber_direction_file_paths_points.set(
      "Fiber_direction_from_file_points", false, {file_name});

  auto mesh = make_two_node_line_mesh(2);
  mesh.N(0, 0) = 0.6;
  mesh.N(1, 0) = 0.4;
  mesh.N(0, 1) = 0.25;
  mesh.N(1, 1) = 0.75;

  Fibers fibers;
  fibers.read(parameters, mesh, 2);

  ASSERT_EQ(fibers.get_num_input_directions(), 1);
  ASSERT_EQ(fibers.get_num_basis_directions(), 0);
  EXPECT_NEAR(fibers.get_quadrature_fibers(2, 0, 0, raw_fibers)(0, 0), 1.0, 1.0e-12);
  EXPECT_NEAR(fibers.get_quadrature_fibers(2, 0, 0, raw_fibers)(1, 0), 0.0, 1.0e-12);
  EXPECT_NEAR(fibers.get_quadrature_fibers(2, 1, 0, raw_fibers)(0, 0), 0.0, 1.0e-12);
  EXPECT_NEAR(fibers.get_quadrature_fibers(2, 1, 0, raw_fibers)(1, 0), 1.0, 1.0e-12);

  std::filesystem::remove(file_name);
}

TEST(Fibers, RejectsLinearDependenceCreatedByPointInterpolation)
{
  const std::string fiber_file =
      ::testing::TempDir() + "dependent_interpolated_fiber.vtu";
  const std::string sheet_file =
      ::testing::TempDir() + "dependent_interpolated_sheet.vtu";
  write_nodal_fiber_file(
    fiber_file, {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}});
  write_nodal_fiber_file(
    sheet_file, {{0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0}});

  MeshParameters parameters;
  parameters.fiber_direction_file_paths_points.set(
      "Fiber_direction_from_file_points", false,
      {fiber_file, sheet_file});

  auto mesh = make_two_node_line_mesh();
  Fibers fibers;
  ASSERT_NO_THROW(fibers.read(parameters, mesh, 3));

  const auto raw = fibers.get_quadrature_fibers(3, 0, 0, raw_fibers);
  ASSERT_EQ(raw.ncols(), 2);
  EXPECT_NEAR(std::abs(raw.col(0).dot(raw.col(1))), 1.0, 1.0e-12);
  EXPECT_THROW(fibers.initialize_fiber_bases(3),
    svmp::ParseException);

  std::filesystem::remove(fiber_file);
  std::filesystem::remove(sheet_file);
}

TEST(Fibers, DyadicInterpolationDoesNotCancelOppositeDirections)
{
  const std::string file_name =
      ::testing::TempDir() + "opposite_nodal_fibers.vtu";
  write_nodal_fiber_file(file_name, {{1.0, 0.0}, {-1.0, 0.0}});

  MeshParameters parameters;
  parameters.fiber_direction_file_paths_points.set(
      "Fiber_direction_from_file_points", false, {file_name});

  auto mesh = make_two_node_line_mesh();
  Fibers fibers;
  fibers.read(parameters, mesh, 2);

  const auto direction = fibers.get_quadrature_fibers(2, 0, 0, raw_fibers);
  EXPECT_NEAR(std::abs(direction(0, 0)), 1.0, 1.0e-12);
  EXPECT_NEAR(direction(1, 0), 0.0, 1.0e-12);

  std::filesystem::remove(file_name);
}

TEST(Fibers, DyadicInterpolationReturnsPrincipalDirection)
{
  const double sqrt_one_half = std::sqrt(0.5);
  const std::string file_name =
      ::testing::TempDir() + "principal_nodal_fibers.vtu";
  write_nodal_fiber_file(file_name, {
    {1.0, 0.0},
    {sqrt_one_half, sqrt_one_half}
  });

  MeshParameters parameters;
  parameters.fiber_direction_file_paths_points.set(
      "Fiber_direction_from_file_points", false, {file_name});

  auto mesh = make_two_node_line_mesh();
  Fibers fibers;
  fibers.read(parameters, mesh, 2);

  const auto direction = fibers.get_quadrature_fibers(2, 0, 0, raw_fibers);
  EXPECT_NEAR(direction(0, 0), std::cos(std::numbers::pi / 8.0), 1.0e-12);
  EXPECT_NEAR(direction(1, 0), std::sin(std::numbers::pi / 8.0), 1.0e-12);

  std::filesystem::remove(file_name);
}

TEST(Fibers, RejectsCellTagWithPointData)
{
  MeshParameters parameters;
  parameters.fiber_direction_file_paths_cells.set(
      "Fiber_direction_from_file_cells", false,
      {std::string(UNIT_TEST_DATA_DIR) + "/fibers_point_data.vtu",
       std::string(UNIT_TEST_DATA_DIR) + "/fibers_point_data_sheet.vtu"});

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  EXPECT_THROW(fibers.read(parameters, mesh, 3), svmp::FileFormatException);
}

TEST(Fibers, RejectsPointTagWithCellData)
{
  MeshParameters parameters;
  parameters.fiber_direction_file_paths_points.set(
      "Fiber_direction_from_file_points", false,
      {std::string(UNIT_TEST_DATA_DIR) + "/fibers_cell_data.vtu",
       std::string(UNIT_TEST_DATA_DIR) + "/fibers_cell_data_sheet.vtu"});

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  EXPECT_THROW(fibers.read(parameters, mesh, 3), svmp::FileFormatException);
}

TEST(Fibers, RejectsMixedInputModes)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{1.0, 0.0, 0.0}, false);
  parameters.fiber_direction_file_paths_cells.set(
      "Fiber_direction_from_file_cells", false,
      {std::string(UNIT_TEST_DATA_DIR) + "/fibers_cell_data.vtu"});

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  EXPECT_THROW(fibers.read(parameters, mesh, 3), svmp::ParseException);
}
