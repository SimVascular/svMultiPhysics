// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "prestrain.h"

#include "ComMod.h"
#include "Core/Exception.h"
#include "FE/Common/FEException.h"
#include "Parameters.h"
#include "VtkData.h"
#include "consts.h"
#include "nn.h"

#include <algorithm>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "mpi.h"

namespace {

/// @brief Name of the VTU cell array holding Gauss point g.
std::string array_name(const int g)
{
  return "Prestrain_F_g" + std::to_string(g);
}

/**
 * @brief Scatters element-associated data from the master to every rank.
 *
 * On the master, @p global_data stores @p values_per_element contiguous values
 * for each of the @c mesh.gnEl elements in the original element ordering. The
 * blocks are reordered with @c mesh.otnIEN so each rank's elements are
 * contiguous, then @c MPI_Scatterv delivers the local partition into
 * @p local_data (@c mesh.nEl blocks).
 */
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


void Prestrain::read(const MeshParameters& mesh_parameters, mshType& mesh, const int nsd)
{
  if (!mesh_parameters.prestrain_file_path.defined()) {
    return;
  }
  const auto file_name = mesh_parameters.prestrain_file_path.value();

  if (!std::filesystem::is_regular_file(file_name)) {
    svmp::raise<svmp::FileNotFoundException>(file_name);
  }
  std::unique_ptr<VtkData> vtk_data(VtkData::create_reader(file_name));

  svmp::throw_if<svmp::FileFormatException>(vtk_data->num_elems() != mesh.gnEl, file_name,
      "The prestrain file has " + std::to_string(vtk_data->num_elems()) + " elements but the mesh named '" +
      mesh.name + "' has " + std::to_string(mesh.gnEl) + ".");

  const int ncomp = nsd*nsd;
  F0.resize(ncomp, mesh.nG, mesh.gnEl);

  Array<double> Fg(ncomp, mesh.gnEl);
  for (int g = 0; g < mesh.nG; g++) {
    const auto name = array_name(g);
    svmp::throw_if<svmp::FileFormatException>(!vtk_data->has_cell_data(name), file_name,
        "No cell array named '" + name + "': the mesh named '" + mesh.name + "' integrates with " +
        std::to_string(mesh.nG) + " Gauss points per element.");
    vtk_data->copy_cell_data(name, Fg);

    for (int e = 0; e < mesh.gnEl; e++) {
      for (int c = 0; c < ncomp; c++) {
        F0(c, g, e) = Fg(c,e);
      }
    }
  }

  // Elements with linear shape functions are pulled back once, which needs the
  // prestrain to be the same at all their Gauss points.
  if (mesh.lShpF) {
    for (int e = 0; e < mesh.gnEl; e++) {
      for (int g = 1; g < mesh.nG; g++) {
        for (int c = 0; c < ncomp; c++) {
          if (F0(c, g, e) != F0(c, 0, e)) {
            svmp::raise<svmp::FileFormatException>(file_name,
                "The prestrain varies between the Gauss points of element " + std::to_string(e) +
                " of the mesh named '" + mesh.name + "', whose shape functions are linear.");
          }
        }
      }
    }
  }
}


void Prestrain::distribute(const CmMod& cm_mod, const cmType& cm, const mshType& mesh)
{
  // Share what the master read with all ranks.
  bool has = has_prestrain();
  cm.bcast(cm_mod, &has);
  if (!has) {
    return;
  }
  int ncomp = F0.nrows();
  cm.bcast(cm_mod, &ncomp);

  // Scatter the prestrain of each rank's elements.
  const Array3<double> global = F0;
  F0.resize(ncomp, mesh.nG, mesh.nEl);
  scatter_element_data(global.data(), F0.data(), ncomp * mesh.nG, mesh, cm_mod, cm);
}


void Prestrain::set_identity(const int nsd, const int nG, const int nEl)
{
  F0.resize(nsd*nsd, nG, nEl);
  F0 = 0.0;
  for (int e = 0; e < nEl; e++) {
    for (int g = 0; g < nG; g++) {
      for (int k = 0; k < nsd; k++) {
        F0(k*(nsd+1), g, e) = 1.0;
      }
    }
  }
}


void Prestrain::gather(const ComMod& com_mod, const CmMod& cm_mod, const mshType& mesh, Array<double>& gF0) const
{
  const auto& cm = com_mod.cm;

  // The row count follows from the mesh, the same on every rank; the columns
  // are the rank's own elements, none on a rank that holds none of the mesh.
  const int rows = com_mod.nsd * com_mod.nsd * mesh.nG;

  if (cm.seq()) {
    gF0.resize(rows, mesh.nEl);
    std::copy_n(F0.data(), F0.size(), gF0.data());
    return;
  }

  const int np = cm.np();
  Vector<int> sCount(np), disps(np);
  for (int i = 0; i < np; i++) {
    disps(i) = mesh.eDist(i) * rows;
    sCount(i) = mesh.eDist(i+1) * rows - disps(i);
  }

  Array<double> tmp;
  if (cm.mas(cm_mod)) {
    tmp.resize(rows, mesh.gnEl);
  }

  MPI_Gatherv(F0.data(), rows*mesh.nEl, cm_mod::mpreal, tmp.data(), sCount.data(), disps.data(),
      cm_mod::mpreal, cm_mod.master, cm.com());

  // Back to the original element order, the inverse of the partitioning.
  if (cm.mas(cm_mod)) {
    gF0.resize(rows, mesh.gnEl);
    for (int e = 0; e < mesh.gnEl; e++) {
      const int Ec = mesh.otnIEN(e);
      gF0.set_col(e, tmp.col(Ec));
    }
  }
}


namespace prestrain {

void init(ComMod& com_mod)
{
  using namespace consts;
  const int nsd = com_mod.nsd;

  if (com_mod.prestrainEq) {
    svmp::throw_if<svmp::FE::InvalidArgumentException>(com_mod.prestressEq,
        "Prestress and Prestrain can not both be set.");
    svmp::throw_if<svmp::NotImplementedException>(com_mod.stFileFlag,
        "A Prestrain run can not be restarted from a .bin file. "
        "Point Prestrain_file_path at its last VTU file instead.");

    // Every prestrain iteration restarts the solution from rest, and its time is
    // fictitious, so any other physics would be erased or advanced meaninglessly.
    const auto& eq = com_mod.eq[0];
    svmp::throw_if<svmp::FE::InvalidArgumentException>(com_mod.nEq != 1,
        "A Prestrain run solves a single struct equation, but this input has " +
        std::to_string(com_mod.nEq) + " equations. Every prestrain iteration starts the solution again "
        "from rest in fictitious time, which other physics cannot follow.");
    // A ustruct solid carries its whole volumetric stress in its pressure field,
    // which the restart from rest discards, so it could never hold a hydrostatic
    // prestress (see USTRUCT_PRESTRAIN.md for what support would take).
    svmp::throw_if<svmp::NotImplementedException>(eq.phys == EquationType::phys_ustruct,
        "Prestrain is not supported for a ustruct equation: the restart from rest at every iteration "
        "discards the pressure field that carries a ustruct solid's volumetric stress. Compute the "
        "prestrain with a struct equation.");
    svmp::throw_if<svmp::FE::InvalidArgumentException>(eq.phys != EquationType::phys_struct,
        "A Prestrain run solves a single struct equation, but this input's equation is not struct.");

    // A prestrain iteration composes its update into F0; a line search would
    // have to compose every trial step to judge it by the composed residual.
    svmp::throw_if<svmp::FE::InvalidArgumentException>(eq.line_search_enabled,
        "Line_search_experimental is not supported in a Prestrain run. It remains available "
        "in forward runs that read a prestrain through Prestrain_file_path.");
  }

  for (auto& msh : com_mod.msh) {
    // A prestrain run that starts from nothing starts from the identity.
    if (com_mod.prestrainEq && !msh.prestrain.has_prestrain()) {
      msh.prestrain.set_identity(nsd, msh.nG, msh.nEl);
    }

    // The prestrain is held at the mesh's Gauss points, which the momentum and
    // continuity integrations of a Taylor-Hood element do not use.
    svmp::throw_if<svmp::NotImplementedException>(msh.prestrain.has_prestrain() && msh.nFs != 1,
        "A prestrain requires the same quadrature for the momentum and continuity "
        "equations (P1P1) on the mesh named '" + msh.name + "'; Taylor-Hood is not supported.");
  }

  // has_prestrain is the same on every rank once distributed, so the checks
  // below raise on every rank; sstEq and rmsh.isReqd are broadcast.
  bool prestrained = com_mod.prestrainEq;
  for (const auto& msh : com_mod.msh) {
    prestrained = prestrained || msh.prestrain.has_prestrain();
  }

  // A ustruct solid's hydrostatic prestress would have to come from its pressure
  // field, which starts from zero.
  svmp::throw_if<svmp::NotImplementedException>(prestrained && com_mod.sstEq,
      "A prestrain can not be applied to a ustruct (velocity-pressure) solid: its hydrostatic "
      "part would have to be carried by the pressure field, which starts from zero. Use a struct "
      "equation for the prestrained solid.");

  // F0 is held per element of the mesh as read; remeshing would leave it behind.
  svmp::throw_if<svmp::NotImplementedException>(prestrained && com_mod.rmsh.isReqd,
      "A prestrain is not supported with remeshing: the prestrain is stored per element of the "
      "original mesh and is not carried to a new one.");
}

double compose_iteration(ComMod& com_mod, SolutionStates& solutions)
{
  // A prestrain run has a single struct equation (see init), so these are all
  // its rows.
  auto& An = solutions.current.get_acceleration();
  auto& Yn = solutions.current.get_velocity();
  auto& Dn = solutions.current.get_displacement();
  const auto& eq = com_mod.eq[0];
  const int nsd = com_mod.nsd;

  // F0 <- (I + Grad u) F0 at every Gauss point, with u = af * D.
  for (auto& msh : com_mod.msh) {
    if (!msh.prestrain.has_prestrain()) {
      continue;
    }
    const int eNoN = msh.eNoN;
    Array<double> xl(nsd,eNoN), dl(nsd,eNoN), Nx(nsd,eNoN), ksix(nsd,nsd);
    double Jac = 0.0;

    for (int e = 0; e < msh.nEl; e++) {
      for (int a = 0; a < eNoN; a++) {
        const int Ac = msh.IEN(a,e);
        for (int i = 0; i < nsd; i++) {
          xl(i,a) = com_mod.x(i,Ac);
          dl(i,a) = eq.af * Dn(eq.s+i,Ac);
        }
      }

      for (int g = 0; g < msh.nG; g++) {
        // Shape function gradients on the mesh, constant within linear elements
        if (g == 0 || !msh.lShpF) {
          auto Nx_g = msh.Nx.slice(g);
          nn::gnn(eNoN, nsd, nsd, Nx_g, xl, Nx, Jac, ksix);
        }
        if (nsd == 3) {
          const mat_fun::Matrix<3> F_step = mat_fun::Matrix<3>::Identity() +
              mat_fun::eigen_view_rows<3>(dl, 0) * mat_fun::eigen_view<3>(Nx).transpose();
          msh.prestrain.compose<3>(g, e, F_step);
        } else {
          const mat_fun::Matrix<2> F_step = mat_fun::Matrix<2>::Identity() +
              mat_fun::eigen_view_rows<2>(dl, 0) * mat_fun::eigen_view<2>(Nx).transpose();
          msh.prestrain.compose<2>(g, e, F_step);
        }
      }
    }
  }

  // This rank's largest nodal update; a NaN gets through rather than being
  // hidden by max().
  double umax = 0.0;
  for (int a = 0; a < com_mod.tnNo; a++) {
    double u2 = 0.0;
    for (int i = 0; i < nsd; i++) {
      u2 += Dn(eq.s+i,a) * Dn(eq.s+i,a);
    }
    umax = (u2 > umax || std::isnan(u2)) ? u2 : umax;
  }

  An = 0.0;
  Yn = 0.0;
  Dn = 0.0;
  return eq.af * std::sqrt(umax);
}

};
