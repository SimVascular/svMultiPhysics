// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "prestrain.h"

#include "Core/Exception.h"
#include "FE/Common/FEException.h"
#include "VtkData.h"
#include "consts.h"
#include "nn.h"

#include <fstream>
#include <memory>

#include "mpi.h"

namespace prestrain {

/// @brief Name of the VTU cell array holding Gauss point g.
static std::string array_name(const int g)
{
  return "Prestrain_F_g" + std::to_string(g);
}

Array<double> deformation_gradient(const mshType& lM, const int e, const int g, const int nsd)
{
  auto F0 = mat_fun::mat_id(nsd);
  if (lM.F0.size() == 0) {
    return F0;
  }

  // Column-major nsd x nsd block at this Gauss point, as Matrix<nsd> lays it out.
  const double* block = lM.F0.data() + e*lM.F0.nrows() + g*nsd*nsd;
  for (int j = 0; j < nsd; j++) {
    for (int i = 0; i < nsd; i++) {
      F0(i,j) = block[i + j*nsd];
    }
  }
  return F0;
}

template <int nsd>
static double pull_back_gradients(const mshType& lM, const int e, const int g, Array<double>& Nx)
{
  const mat_fun::Matrix<nsd> F0 = deformation_gradient<nsd>(lM, e, g);
  auto Nxm = mat_fun::eigen_view_mutable(Nx);
  Nxm = F0.transpose() * Nxm;
  return F0.determinant();
}

double pull_back(const mshType& lM, const int e, const int g, Array<double>& Nx)
{
  if (Nx.nrows() == 3) {
    return pull_back_gradients<3>(lM, e, g, Nx);
  } else {
    return pull_back_gradients<2>(lM, e, g, Nx);
  }
}

void pull_back_fibers(const mshType& lM, const int e, const int g, Array<double>& fibers)
{
  if (fibers.nrows() == 3) {
    pull_back_fibers<3>(deformation_gradient<3>(lM, e, g), fibers);
  } else {
    pull_back_fibers<2>(deformation_gradient<2>(lM, e, g), fibers);
  }
}

void read(const std::string& file_name, mshType& mesh, const int nsd)
{
  svmp::throw_if<svmp::FileNotFoundException>(!std::ifstream(file_name).good(), file_name);
  std::unique_ptr<VtkData> vtk_data(VtkData::create_reader(file_name));

  svmp::throw_if<svmp::FileFormatException>(vtk_data->num_elems() != mesh.gnEl, file_name,
      "The prestrain file has " + std::to_string(vtk_data->num_elems()) + " elements but the mesh named '" +
      mesh.name + "' has " + std::to_string(mesh.gnEl) + ".");

  const int ncomp = nsd*nsd;
  mesh.F0 = Array<double>(ncomp*mesh.nG, mesh.gnEl);

  Array<double> Fg(ncomp, mesh.gnEl);
  for (int g = 0; g < mesh.nG; g++) {
    const auto name = array_name(g);
    svmp::throw_if<svmp::FileFormatException>(!vtk_data->has_cell_data(name), file_name,
        "No cell array named '" + name + "': the mesh named '" + mesh.name + "' integrates with " +
        std::to_string(mesh.nG) + " Gauss points per element.");
    vtk_data->copy_cell_data(name, Fg);

    for (int e = 0; e < mesh.gnEl; e++) {
      for (int c = 0; c < ncomp; c++) {
        mesh.F0(g*ncomp + c, e) = Fg(c,e);
      }
    }
  }

  // Elements with linear shape functions are pulled back once, which needs the
  // prestrain to be the same at all their Gauss points.
  if (mesh.lShpF) {
    for (int e = 0; e < mesh.gnEl; e++) {
      for (int g = 1; g < mesh.nG; g++) {
        for (int c = 0; c < ncomp; c++) {
          if (mesh.F0(g*ncomp + c, e) != mesh.F0(c, e)) {
            svmp::raise<svmp::FileFormatException>(file_name,
                "The prestrain varies between the Gauss points of element " + std::to_string(e) +
                " of the mesh named '" + mesh.name + "', whose shape functions are linear.");
          }
        }
      }
    }
  }
}

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
    const int rows = nsd*nsd*msh.nG;

    // The prestrain is read with the mesh's own number of Gauss points.
    svmp::throw_if<svmp::InternalErrorException>(msh.F0.size() != 0 && msh.F0.nrows() != rows,
        "The prestrain of the mesh named '" + msh.name + "' holds " + std::to_string(msh.F0.nrows() / (nsd*nsd)) +
        " Gauss points per element but the mesh integrates with " + std::to_string(msh.nG) + ".");

    if (com_mod.prestrainEq && msh.F0.size() == 0) {
      // A prestrain run that starts from nothing starts from the identity.
      msh.F0 = Array<double>(rows, msh.nEl);
      msh.F0 = 0.0;
      for (int e = 0; e < msh.nEl; e++) {
        for (int g = 0; g < msh.nG; g++) {
          for (int k = 0; k < nsd; k++) {
            msh.F0(g*nsd*nsd + k*(nsd+1), e) = 1.0;
          }
        }
      }
    }

    // The prestrain is held at the mesh's Gauss points, which the momentum and
    // continuity integrations of a Taylor-Hood element do not use.
    svmp::throw_if<svmp::NotImplementedException>(msh.F0.size() != 0 && msh.nFs != 1,
        "A prestrain requires the same quadrature for the momentum and continuity "
        "equations (P1P1) on the mesh named '" + msh.name + "'; Taylor-Hood is not supported.");
  }

  // After partitioning every rank's F0 has the mesh's row count, so the checks
  // below raise on every rank; sstEq and rmsh.isReqd are broadcast.
  bool prestrained = com_mod.prestrainEq;
  for (const auto& msh : com_mod.msh) {
    prestrained = prestrained || msh.F0.nrows() != 0;
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

/// @brief F0 <- (I + Grad u) F0 at Gauss point g of element e, in place.
template <int nsd>
static void compose(Array<double>& F0, const int e, const int g, const Array<double>& dl, const Array<double>& Nx)
{
  Eigen::Map<mat_fun::Matrix<nsd>> F(F0.data() + e*F0.nrows() + g*nsd*nsd);
  const mat_fun::Matrix<nsd> F_step = mat_fun::Matrix<nsd>::Identity() +
      mat_fun::eigen_view_rows<nsd>(dl, 0) * mat_fun::eigen_view<nsd>(Nx).transpose();
  F = F_step * F;
}

/// @brief F0 <- (I + Grad u) F0 at every Gauss point, with u = af * D, the
/// rows s to s + nsd - 1 of D: the displacement of the run's one solid equation.
static void accumulate(ComMod& com_mod, const int s, const double af, const Array<double>& D)
{
  const int nsd = com_mod.nsd;

  for (auto& msh : com_mod.msh) {
    if (msh.F0.size() == 0) {
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
          dl(i,a) = af * D(s+i,Ac);
        }
      }

      for (int g = 0; g < msh.nG; g++) {
        // Shape function gradients on the mesh, constant within linear elements
        if (g == 0 || !msh.lShpF) {
          auto Nx_g = msh.Nx.slice(g);
          nn::gnn(eNoN, nsd, nsd, Nx_g, xl, Nx, Jac, ksix);
        }
        if (nsd == 3) {
          compose<3>(msh.F0, e, g, dl, Nx);
        } else {
          compose<2>(msh.F0, e, g, dl, Nx);
        }
      }
    }
  }
}

double compose_iteration(ComMod& com_mod, SolutionStates& solutions)
{
  // A prestrain run has a single struct equation (see init), so these are all
  // its rows.
  auto& An = solutions.current.get_acceleration();
  auto& Yn = solutions.current.get_velocity();
  auto& Dn = solutions.current.get_displacement();
  const auto& eq = com_mod.eq[0];

  accumulate(com_mod, eq.s, eq.af, Dn);

  // This rank's largest nodal update; a NaN gets through rather than being
  // hidden by max().
  double umax = 0.0;
  for (int a = 0; a < com_mod.tnNo; a++) {
    double u2 = 0.0;
    for (int i = 0; i < com_mod.nsd; i++) {
      u2 += Dn(eq.s+i,a) * Dn(eq.s+i,a);
    }
    umax = (u2 > umax || std::isnan(u2)) ? u2 : umax;
  }

  An = 0.0;
  Yn = 0.0;
  Dn = 0.0;
  return eq.af * std::sqrt(umax);
}

void gather(const ComMod& com_mod, const CmMod& cm_mod, const mshType& lM, Array<double>& gF0)
{
  const auto& cm = com_mod.cm;

  if (cm.seq()) {
    gF0 = lM.F0;
    return;
  }

  // The row count follows from the mesh, the same on every rank; the columns
  // are the rank's own elements, none on a rank that holds none of the mesh.
  const int rows = com_mod.nsd * com_mod.nsd * lM.nG;
  const int np = cm.np();

  Vector<int> sCount(np), disps(np);
  for (int i = 0; i < np; i++) {
    disps(i) = lM.eDist(i) * rows;
    sCount(i) = lM.eDist(i+1) * rows - disps(i);
  }

  Array<double> tmp;
  if (cm.mas(cm_mod)) {
    tmp.resize(rows, lM.gnEl);
  }

  MPI_Gatherv(lM.F0.data(), rows*lM.nEl, cm_mod::mpreal, tmp.data(), sCount.data(), disps.data(),
      cm_mod::mpreal, cm_mod.master, cm.com());

  // Back to the original element order, the inverse of the partitioning.
  if (cm.mas(cm_mod)) {
    gF0.resize(rows, lM.gnEl);
    for (int e = 0; e < lM.gnEl; e++) {
      const int Ec = lM.otnIEN(e);
      gF0.set_col(e, tmp.col(Ec));
    }
  }
}

};
