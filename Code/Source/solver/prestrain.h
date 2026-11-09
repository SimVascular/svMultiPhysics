// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PRESTRAIN_H
#define PRESTRAIN_H

#include "Array.h"
#include "CmMod.h"
#include "ComMod.h"
#include "SolutionStates.h"
#include "mat_fun.h"

#include <array>
#include <cmath>
#include <string>

/// @brief Prestrain of a body imaged in a loaded configuration.
///
/// The mesh is the imaged, loaded geometry. The prestrain is a deformation
/// gradient F0 at every Gauss point that maps an unknown stress-free
/// configuration to the mesh. It need not be the gradient of any displacement
/// field: the stress-free configuration may be incompatible, as with residual
/// stresses. A displacement u from the mesh then gives the deformation
/// gradient from the stress-free configuration
///
///     F = (I + Grad u) F0,
///
/// with Grad taken on the mesh. The solid kernels integrate on the stress-free
/// configuration: where the assembly loops compute the shape function
/// gradients, they pull them back in place, F0^T Grad N, together with the
/// volume, dV / det F0, and the fibers, which are given on the mesh and pulled
/// back as material lines with F0^{-1}. The kernels then form F = F0 + Grad_0 u
/// and keep their total Lagrangian form. The formulation is objective:
/// rotating the deformed body rigidly rotates F and the internal forces with it.
///
/// Like the shape function gradients, F0 is constant within an element whose
/// shape functions are linear, so those elements are pulled back once. A
/// prestrain run produces such fields, and read() checks fields from file.
///
/// F0 is found incrementally. Every time step starts from rest on the imaged
/// geometry under the full load; the displacement u it reaches is accumulated
/// as F0 <- (I + Grad u) F0 and discarded. At the fixed point the mesh carries
/// the load in Cauchy stress equilibrium without moving. This is the
/// incremental scheme of Gee, Forster and Wall, Int. J. Numer. Meth. Biomed.
/// Engng. 26 (2010) 52-72, section 3, with the increments composed on the
/// imaged mesh, so that equilibrium holds on the imaged geometry rather than
/// on a virtually deformed one. Each step is iterated like any other time
/// step, up to Max_iterations. With Pseudo_transient its time step may be
/// adapted to the residual (see PseudoTransientContinuation.h).
///
/// F0 is stored per mesh as Array<double>(nsd*nsd*nG, nEl): one column per
/// element holding the nG Gauss point tensors back to back, each in the
/// column-major order of Matrix<nsd>. Written to and read from VTU cell arrays
/// named Prestrain_F_g<g>, one per Gauss point.
///
namespace prestrain {

/// @brief Deformation gradient from the stress-free configuration to the mesh
/// at Gauss point g of element e, or the identity when the mesh carries none.
template <int nsd>
mat_fun::Matrix<nsd> deformation_gradient(const mshType& lM, const int e, const int g)
{
  if (lM.F0.size() == 0) {
    return mat_fun::Matrix<nsd>::Identity();
  }
  return Eigen::Map<const mat_fun::Matrix<nsd>>(lM.F0.data() + e*lM.F0.nrows() + g*nsd*nsd);
}

/// @brief Element average of the deformation gradient from the stress-free
/// configuration, or the identity when the mesh carries none. For output whose
/// quadrature differs from the mesh's.
template <int nsd>
mat_fun::Matrix<nsd> mean_deformation_gradient(const mshType& lM, const int e)
{
  mat_fun::Matrix<nsd> F0 = mat_fun::Matrix<nsd>::Zero();
  for (int g = 0; g < lM.nG; g++) {
    F0 += deformation_gradient<nsd>(lM, e, g);
  }
  return F0 / lM.nG;
}

/// @brief The same, for callers whose number of spatial dimensions is not
/// known at compile time.
Array<double> deformation_gradient(const mshType& lM, int e, int g, int nsd);

/// @brief Pull the shape function gradients Nx of Gauss point g of element e
/// back to the stress-free configuration, in place: Nx <- F0^T Nx. Returns
/// det F0, the ratio of mesh to stress-free volume, by which the integration
/// Jacobian is divided.
double pull_back(const mshType& lM, int e, int g, Array<double>& Nx);

/// @brief Whether two directions given on the mesh are perpendicular, both
/// being nonzero.
template <class A, class B>
bool perpendicular_on_mesh(const A& u, const B& v)
{
  constexpr double tol = 1.0e-6;
  const double lu = u.norm();
  const double lv = v.norm();
  return lu > 0.0 && lv > 0.0 && std::abs(u.dot(v)) <= tol * lu * lv;
}

/// @brief Pull fibers given on the mesh back to the stress-free configuration
/// that F0 maps to the mesh, in place, nsd x nFn.
///
/// Each direction is a material line, pulled back with F0^-1 and rescaled to
/// its given length. A sheet given perpendicular to the fiber stays
/// perpendicular to it, within the pulled-back plane of the two, and in 3D a
/// third direction given perpendicular to both is their plane's normal and
/// stays normal to the pulled-back plane. So the orthonormal basis
/// Fibers::get_quadrature_fibers builds for some material models stays
/// orthonormal. fibers must own its storage, as an Array assigned from
/// Fibers::get_quadrature_fibers does.
template <int nsd>
void pull_back_fibers(const mat_fun::Matrix<nsd>& F0, Array<double>& fibers)
{
  using Direction = Eigen::Matrix<double, nsd, 1>;
  auto f0 = mat_fun::eigen_view_mutable(fibers);
  const int nFn = f0.cols();
  if (nFn == 0) {
    return;
  }
  const Eigen::Matrix<double, nsd, Eigen::Dynamic> f = f0;   // as given on the mesh
  const mat_fun::Matrix<nsd> F0_inv = F0.inverse();

  for (int k = 0; k < nFn; k++) {
    const Direction d = F0_inv * f.col(k);
    const double length = d.norm();
    f0.col(k) = (length > 0.0) ? Direction(d * (f.col(k).norm() / length)) : Direction::Zero();
  }

  if (nFn >= 2 && perpendicular_on_mesh(f.col(0), f.col(1))) {
    const Direction a = f0.col(0).normalized();
    const Direction s = f0.col(1) - a.dot(f0.col(1)) * a;
    f0.col(1) = s * (f.col(1).norm() / s.norm());

    if constexpr (nsd == 3) {
      if (nFn >= 3 && perpendicular_on_mesh(f.col(2), f.col(0)) && perpendicular_on_mesh(f.col(2), f.col(1))) {
        const Direction n = a.cross(Direction(f0.col(1).normalized()));
        const double side = Direction(f.col(2)).dot(Direction(f.col(0)).cross(Direction(f.col(1)))) < 0.0 ? -1.0 : 1.0;
        f0.col(2) = (side * f.col(2).norm()) * n;
      }
    }
  }
}

/// @brief The same with the prestrain of Gauss point g of element e, for
/// callers whose number of spatial dimensions is not known at compile time.
void pull_back_fibers(const mshType& lM, int e, int g, Array<double>& fibers);

/// @brief Make a solid kernel's tangent that of a prestrain step.
///
/// A prestrain step composes the displacement u it reaches into F0 and starts
/// again from rest, so the residual the next step sees is that of the Cauchy
/// stress on the mesh, which does not move. The kernel's tangent also carries
/// the change of the current volume and of the current shape function gradients
/// with u; this subtracts that part,
///
///     w (tau Grad_x N_a (x) Grad_x N_b - tau Grad_x N_b (x) Grad_x N_a),
///
/// with tau = P F^T the Kirchhoff stress, Grad_x N = F^-T Grad_0 N, and w the
/// stress-free weight. Without it each step is an inexact Newton step whose
/// error grows with the stress, and the prestrain can stop converging once the
/// pseudo time step is large. The fibers' part is added by
/// correct_step_tangent_fibers.
///
/// Only for steps of one Newton iteration (Max_iterations 1). A step iterated
/// further solves the kernel's own residual, whose tangent the kernel already
/// has, and is composed only once it stops.
///
/// @param Nx0 shape function gradients on the stress-free configuration.
/// @param w_afu the weight times the stiffness scaling of lK.
/// @param lK the element tangent, lK(i*dof + j, a, b).
template <int nsd, class Gradients>
void correct_step_tangent(const mat_fun::Matrix<nsd>& F, const mat_fun::Matrix<nsd>& P, const Gradients& Nx0,
                          const double w_afu, const int dof, Array3<double>& lK)
{
  const mat_fun::NodalMatrix<nsd> Nxs = F.inverse().transpose() * Nx0;
  const mat_fun::NodalMatrix<nsd> tau_Nxs = (P * F.transpose()) * Nxs;
  for (int b = 0; b < Nxs.cols(); b++) {
    for (int a = 0; a < Nxs.cols(); a++) {
      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          lK(i*dof + j, a, b) -= w_afu * (tau_Nxs(i,a) * Nxs(j,b) - tau_Nxs(i,b) * Nxs(j,a));
        }
      }
    }
  }
}

/// @brief An orthonormal basis of the directions perpendicular to the unit
/// vector d.
template <int nsd>
std::array<Eigen::Matrix<double, nsd, 1>, nsd - 1> perpendicular_directions(const Eigen::Matrix<double, nsd, 1>& d)
{
  if constexpr (nsd == 2) {
    return {Eigen::Vector2d(-d(1), d(0))};
  } else {
    // Cross with the axis least aligned with d.
    int k = 0;
    d.cwiseAbs().minCoeff(&k);
    const Eigen::Vector3d t = d.cross(Eigen::Vector3d::Unit(k)).normalized();
    return {t, d.cross(t)};
  }
}

/// @brief Number of leading fibers that form a frame: a fiber, a sheet kept
/// perpendicular to it and, in 3D, their plane's normal (see pull_back_fibers).
/// Inferred from the pulled-back fibers, which are then perpendicular to
/// round-off; 0 when there is no such frame.
template <int nsd, class Fibers>
int fiber_frame_size(const Fibers& fibers)
{
  constexpr double tol = 1.0e-10;
  auto perpendicular = [&](const int i, const int j) {
    const double li = fibers.col(i).norm();
    const double lj = fibers.col(j).norm();
    return li > 0.0 && lj > 0.0 && std::abs(fibers.col(i).dot(fibers.col(j))) <= tol * li * lj;
  };
  const int nFn = fibers.cols();
  if (nFn < 2 || !perpendicular(0, 1)) {
    return 0;
  }
  if (nsd == 3 && nFn >= 3 && perpendicular(2, 0) && perpendicular(2, 1)) {
    return 3;
  }
  return 2;
}

/// @brief Add the fibers' part of a prestrain step's tangent to a solid
/// kernel's tangent.
///
/// Composing a step into F0 pulls the mesh's fibers back anew (see
/// pull_back_fibers), so the fibers the material law sees change with u. To
/// first order, with G = F0^-1 Grad u F0, the gradient on the mesh carried to
/// the stress-free configuration:
///
/// - a fiber pulled back on its own changes by df0 = -(I - f^ f^T) G f0;
/// - a frame of a fiber a^, a sheet s^ kept perpendicular to it and their
///   normal n^ = a^ x s^ turns, without stretching, by the angle w with
///
///       w.a^ = -n^.G s^,    w.s^ = n^.G a^,    w.n^ = -s^.G a^
///
///   in 3D, or w = -p^.G a^ in 2D, p^ being a^ turned a quarter.
///
/// ^ marks unit vectors. The kernel's tangent holds the fibers fixed; this adds
///
///     w F dS Grad_0 N_a,
///
/// taking the change dS of the material law's S by one-sided finite
/// differences: along each axis a frame turns about (three in 3D, one in 2D),
/// and along the nsd - 1 directions perpendicular to each other fiber.
///
/// @param S the material law's stress at F and the given fibers.
/// @param Nx0 shape function gradients on the stress-free configuration.
/// @param fibers the pulled-back fibers, nsd x nFn.
/// @param stress the material law's S at F for other fibers, nsd x nFn.
/// @param w_afu the weight times the stiffness scaling of lK.
/// @param lK the element tangent, lK(i*dof + j, a, b).
template <int nsd, class Gradients, class Fibers, class Stress>
void correct_step_tangent_fibers(const mat_fun::Matrix<nsd>& F0, const mat_fun::Matrix<nsd>& F,
                                 const mat_fun::Matrix<nsd>& S, const Gradients& Nx0, const Fibers& fibers,
                                 const Stress& stress, const double w_afu, const int dof, Array3<double>& lK)
{
  using Direction = Eigen::Matrix<double, nsd, 1>;
  const int nFn = fibers.cols();
  if (nFn == 0) {
    return;
  }

  // Step of the finite differences: an angle, or a length relative to the fiber's
  constexpr double h = 1.0e-7;

  const mat_fun::Matrix<nsd> F0_inv_t = F0.inverse().transpose();
  Eigen::Matrix<double, nsd, Eigen::Dynamic> f = fibers;

  // lK += w_afu (F dS Grad_0 N_a) (x) C(:,b), where the fibers' change dS of S
  // is linear in du through the scalar sum over b of C(:,b) . du_b.
  auto add = [&](const mat_fun::Matrix<nsd>& dS, const mat_fun::NodalMatrix<nsd>& C) {
    const mat_fun::NodalMatrix<nsd> dP_Nx = (F * dS) * Nx0;
    for (int b = 0; b < Nx0.cols(); b++) {
      for (int a = 0; a < Nx0.cols(); a++) {
        for (int i = 0; i < nsd; i++) {
          for (int j = 0; j < nsd; j++) {
            lK(i*dof + j, a, b) += w_afu * dP_Nx(i,a) * C(j,b);
          }
        }
      }
    }
  };

  // A frame turns as a whole, about each of its axes in turn.
  const int frame = fiber_frame_size<nsd>(fibers);
  if (frame > 0) {
    const Direction a = fibers.col(0).normalized();
    const mat_fun::NodalVector a_Nx = Nx0.transpose() * a;
    auto turn = [&](const auto& rotate, const mat_fun::NodalMatrix<nsd>& C) {
      for (int k = 0; k < frame; k++) {
        f.col(k) = fibers.col(k) + h * rotate(Direction(fibers.col(k)));
      }
      add((stress(f) - S) / h, C);
      for (int k = 0; k < frame; k++) {
        f.col(k) = fibers.col(k);
      }
    };
    if constexpr (nsd == 3) {
      const Direction s = fibers.col(1).normalized();
      const Direction n = a.cross(s);
      const mat_fun::NodalVector s_Nx = Nx0.transpose() * s;
      turn([&](const Direction& v) { return Direction(a.cross(v)); }, -(F0_inv_t * n) * s_Nx.transpose());
      turn([&](const Direction& v) { return Direction(s.cross(v)); },  (F0_inv_t * n) * a_Nx.transpose());
      turn([&](const Direction& v) { return Direction(n.cross(v)); }, -(F0_inv_t * s) * a_Nx.transpose());
    } else {
      const Direction p(-a(1), a(0));
      turn([](const Direction& v) { return Direction(-v(1), v(0)); }, -(F0_inv_t * p) * a_Nx.transpose());
    }
  }

  // Every other fiber is pulled back on its own.
  for (int k = frame; k < nFn; k++) {
    const Direction fk = fibers.col(k);
    const double lk = fk.norm();
    if (lk == 0.0) {
      continue;
    }
    const mat_fun::NodalVector fk_Nx = Nx0.transpose() * fk;

    for (const Direction& t : perpendicular_directions<nsd>(Direction(fk / lk))) {
      f.col(k) = fk + (h * lk) * t;
      const mat_fun::Matrix<nsd> dS = (stress(f) - S) / (h * lk);
      f.col(k) = fk;
      add(dS, -(F0_inv_t * t) * fk_Nx.transpose());     // t . df0
    }
  }
}

/// @brief Read the prestrain of a mesh from the Prestrain_F_g<g> cell arrays
/// of a VTU file, in the mesh's original element order. Called on the master
/// before the mesh is partitioned.
void read(const std::string& file_name, mshType& mesh, int nsd);

/// @brief Validate the settings and the prestrains read from file, and seed
/// the identity where a prestrain run has none. Called once the meshes are
/// partitioned.
void init(ComMod& com_mod);

/// @brief Start a prestrain step: accumulate the displacement u the previous
/// step reached into the prestrain, F0 <- (I + Grad u) F0 at every Gauss
/// point, and reset the solution to rest. Called from the Integrator's
/// predictor, once the previous step's solution has become the old one.
///
/// The kernels saw the generalized-alpha level displacement alpha_f*Dn, which
/// is the state the step reached, so that is u; Dn itself would overshoot by
/// 1/alpha_f.
void start_step(ComMod& com_mod, SolutionStates& solutions);

/// @brief Gather a mesh's prestrain onto the master, in the mesh's original
/// element order, for writing. All ranks must call it.
void gather(const ComMod& com_mod, const CmMod& cm_mod, const mshType& lM, Array<double>& gF0);

};

#endif
