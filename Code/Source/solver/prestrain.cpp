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
#include "solid/mat_models.h"

#include <algorithm>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "mpi.h"

namespace {

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
    const auto name = "Prestrain_F_g" + std::to_string(g);
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


Eigen::MatrixXd Prestrain::get(const int nsd, const int g, const int e) const
{
  return Eigen::Map<const Eigen::MatrixXd>(block(g, e), nsd, nsd);
}


Eigen::MatrixXd Prestrain::get_mean(const int nsd, const int e) const
{
  Eigen::MatrixXd mean = Eigen::MatrixXd::Zero(nsd, nsd);
  for (int g = 0; g < F0.ncols(); g++) {
    mean += get(nsd, g, e);
  }
  return mean / F0.ncols();
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


void Prestrain::compose(const ComMod& com_mod, const mshType& mesh, const Array<double>& u)
{
  const int nsd = com_mod.nsd;
  const int eNoN = mesh.eNoN;
  Array<double> xl(nsd,eNoN), dl(nsd,eNoN), Nx(nsd,eNoN), ksix(nsd,nsd);
  double Jac = 0.0;

  for (int e = 0; e < mesh.nEl; e++) {
    for (int a = 0; a < eNoN; a++) {
      const int Ac = mesh.IEN(a,e);
      for (int i = 0; i < nsd; i++) {
        xl(i,a) = com_mod.x(i,Ac);
        dl(i,a) = u(i,Ac);
      }
    }

    for (int g = 0; g < mesh.nG; g++) {
      // Shape function gradients on the mesh, constant within linear elements
      if (g == 0 || !mesh.lShpF) {
        auto Nx_g = mesh.Nx.slice(g);
        nn::gnn(eNoN, nsd, nsd, Nx_g, xl, Nx, Jac, ksix);
      }
      // F0 <- (I + Grad u) F0
      if (nsd == 3) {
        const mat_fun::Matrix<3> gradu = mat_fun::eigen_view<3>(dl) * mat_fun::eigen_view<3>(Nx).transpose();
        Eigen::Map<mat_fun::Matrix<3>> F0_g(block(g, e));
        F0_g = (mat_fun::Matrix<3>::Identity() + gradu) * F0_g;
      } else {
        const mat_fun::Matrix<2> gradu = mat_fun::eigen_view<2>(dl) * mat_fun::eigen_view<2>(Nx).transpose();
        Eigen::Map<mat_fun::Matrix<2>> F0_g(block(g, e));
        F0_g = (mat_fun::Matrix<2>::Identity() + gradu) * F0_g;
      }
    }
  }
}


namespace prestrain {

template <int nsd>
void correct_composed_tangent(const ComMod& com_mod, const CepMod& cep_mod, const dmnType& dmn,
    const ActiveStress::Evaluator& active_stress_evaluator, const Eigen::Ref<const Eigen::VectorXd>& N,
    const mat_fun::Matrix<nsd>& F0, const mat_fun::Matrix<nsd>& F, const mat_fun::Matrix<nsd>& P,
    const Eigen::Ref<const Eigen::Matrix<double, nsd, Eigen::Dynamic>>& Nx0,
    const mat_models::FiberRef<nsd>& fibers, const double w_afu, const int dof, Array3<double>& lK)
{
  using Matrix = mat_fun::Matrix<nsd>;
  using NodalMatrix = mat_fun::NodalMatrix<nsd>;
  using Direction = Eigen::Matrix<double, nsd, 1>;
  const int eNoN = Nx0.cols();

  // 1. The current configuration does not move with u. The kernel's tangent
  //    carries the change of the current volume and shape function gradients
  //    with u; take it out:
  //        -w [ (tau Grad_x N_a)_i (Grad_x N_b)_j - (tau Grad_x N_b)_i (Grad_x N_a)_j ],
  //    tau = P F^T, Grad_x N = F^-T Grad_0 N.
  const NodalMatrix Nxs = F.inverse().transpose() * Nx0;
  const NodalMatrix tau_Nxs = (P * F.transpose()) * Nxs;
  for (int b = 0; b < eNoN; b++) {
    for (int a = 0; a < eNoN; a++) {
      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          lK(i*dof + j, a, b) -= w_afu * (tau_Nxs(i,a) * Nxs(j,b) - tau_Nxs(i,b) * Nxs(j,a));
        }
      }
    }
  }

  // 2. The fibers turn with u. Composing u pulls the mesh's fibers back anew,
  //    so each fiber f0 the material law sees changes, to first order, by
  //        df0 = -(I - f^ f^T) G f0,   G = F0^-1 Grad u F0,
  //    turning without stretching. The kernel's tangent holds the fibers
  //    fixed; add w F dS Grad_0 N_a, with dS the material law's change along
  //    the turning: its gradient with respect to the fiber, by one-sided
  //    finite differences along the axes, applied to df0.
  const int nfd = fibers.cols();
  if ((fibers.colwise().norm().array() == 0.0).all()) {   // Fibers::zero_fibers
    return;
  }

  // S of the material law at F for given fibers; the active tension, Dm and
  // Ja are recomputed with each call and discarded.
  const auto stress = [&](const mat_models::FiberRef<nsd>& f) {
    Matrix S;
    mat_fun::Matrix<3*(nsd-1)> Dm;
    double Ja;
    const auto Ta = active_stress_evaluator.evaluate<nsd>(N, F, f);
    mat_models::compute_pk2cc<nsd>(com_mod, cep_mod, dmn, F, nfd, f, Ta, S, Dm, Ja);
    return S;
  };
  const Matrix S = stress(fibers);

  constexpr double h = 1.0e-7;   // finite difference step, relative to the fiber's length
  const Matrix F0_inv = F0.inverse();
  Eigen::Matrix<double, nsd, Eigen::Dynamic> f = fibers;   // the fibers, one of them perturbed

  for (int k = 0; k < nfd; k++) {
    const Direction fk = fibers.col(k);
    const double lk = fk.norm();
    if (lk == 0.0) {
      continue;
    }
    const Direction d = fk / lk;

    // Grad_0 N_b . f0 = Grad N_b . F0 f0, the gradient on the mesh along the fiber
    const mat_fun::NodalVector fk_Nx = Nx0.transpose() * fk;

    // The turning, df0 = -(I - d d^T) F0^-1 du_b (Grad_0 N_b . f0): T maps
    // du_b to df0 for node b, up to the factor fk_Nx(b).
    const Matrix T = -(Matrix::Identity() - d * d.transpose()) * F0_inv;

    for (int j = 0; j < nsd; j++) {
      // dS along the axis e_j
      f(j, k) = fk(j) + h * lk;
      const Matrix dS = (stress(f) - S) / (h * lk);
      f(j, k) = fk(j);

      // df0_j = C(:,b) . du_b for every node b, so lK += w_afu (F dS Grad_0 N_a) (x) C(:,b)
      const NodalMatrix dP_Nx = (F * dS) * Nx0;
      const NodalMatrix C = T.row(j).transpose() * fk_Nx.transpose();
      for (int b = 0; b < eNoN; b++) {
        for (int a = 0; a < eNoN; a++) {
          for (int i = 0; i < nsd; i++) {
            for (int l = 0; l < nsd; l++) {
              lK(i*dof + l, a, b) += w_afu * dP_Nx(i,a) * C(l,b);
            }
          }
        }
      }
    }
  }
}

template void correct_composed_tangent<2>(const ComMod&, const CepMod&, const dmnType&,
    const ActiveStress::Evaluator&, const Eigen::Ref<const Eigen::VectorXd>&,
    const mat_fun::Matrix<2>&, const mat_fun::Matrix<2>&, const mat_fun::Matrix<2>&,
    const Eigen::Ref<const Eigen::Matrix<double, 2, Eigen::Dynamic>>&,
    const Eigen::Ref<const Eigen::Matrix<double, 2, Eigen::Dynamic>>&, double, int, Array3<double>&);
template void correct_composed_tangent<3>(const ComMod&, const CepMod&, const dmnType&,
    const ActiveStress::Evaluator&, const Eigen::Ref<const Eigen::VectorXd>&,
    const mat_fun::Matrix<3>&, const mat_fun::Matrix<3>&, const mat_fun::Matrix<3>&,
    const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>&,
    const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>&, double, int, Array3<double>&);

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

};
