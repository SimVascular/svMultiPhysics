// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef STRUCT_H 
#define STRUCT_H 

#include "ActiveStress.h"
#include "ComMod.h"
#include "SolutionStates.h"
#include "solid/mat_models.h"

namespace struct_ns {

void b_struct_2d(const ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, 
    const Array<double>& Nx, const Array<double>& dl, const Vector<double>& hl, const Vector<double>& nV, 
    Array<double>& lR, Array3<double>& lK);

void b_struct_3d(const ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, 
    const Array<double>& Nx, const Array<double>& dl, const Vector<double>& hl, const Vector<double>& nV, 
    Array<double>& lR, Array3<double>& lK);

void construct_dsolid(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const SolutionStates& solutions);

/// @brief Contribution of one Gauss point to the element residual and
/// tangent of the struct equation (displacement-based nonlinear solid).
///
/// Adds the inertia, damping, body force and internal force terms of the
/// momentum balance in total Lagrangian form, at the generalized-alpha time
/// levels.
///
/// @param[in] com_mod Simulation data: the current equation and domain, and
///   the time integration parameters.
/// @param[in] cep_mod Electrophysiology data.
/// @param[in] eNoN Number of nodes of the element.
/// @param[in] w Gauss weight times the Jacobian determinant of the element.
/// @param[in] N Shape functions at the Gauss point, eNoN.
/// @param[in] Nx Shape function gradients at the Gauss point, nsd x eNoN.
/// @param[in] al Nodal accelerations of the element, tDof x eNoN.
/// @param[in] yl Nodal velocities of the element, tDof x eNoN.
/// @param[in] dl Nodal displacements of the element, tDof x eNoN.
/// @param[in] bfl Nodal body forces of the element, nsd x eNoN.
/// @param[in] F0 Prestrain at the Gauss point: the deformation gradient from
///   the unloaded configuration to the mesh, so that F = (I + Grad u) F0
///   (see Prestrain). The identity when the mesh has no prestrain. If present,
///   then w, Nx, and fibers are pulled back to the unloaded configuration.
/// @param[in] fibers Fiber directions at the Gauss point, one per column.
/// @param[in] pS0l Nodal prestress of the element in Voigt order.
/// @param[out] pSl PK2 stress at the Gauss point in Voigt order.
/// @param[in] active_stress_evaluator Active stress model.
/// @param[in,out] lR Element residual, dof x eNoN.
/// @param[in,out] lK Element tangent, indexed lK(i*dof + j, a, b) for the
///   derivative of component i at node a with respect to component j at node b.
/// @param[in] recompute_visc Whether the viscosity model has to recompute the
///   quantities that depend on the shape function gradients.
void struct_2d(ComMod &com_mod, CepMod &cep_mod, const int eNoN,
               const double w, const Vector<double> &N, const Array<double> &Nx,
               const Array<double> &al, const Array<double> &yl,
               const Array<double> &dl, const Array<double> &bfl,
               const Eigen::Matrix<double, 2, 2> &F0,
               const Eigen::Ref<const Eigen::Matrix<double, 2, Eigen::Dynamic>> &fibers,
               const Array<double> &pS0l,
               Vector<double> &pSl, const ActiveStress::Evaluator &active_stress_evaluator,
               Array<double> &lR, Array3<double> &lK, const bool recompute_visc);

/// @brief Contribution of one Gauss point to the element residual and
/// tangent of the struct equation (displacement-based nonlinear solid), in three
/// dimensions.
///
/// Adds the inertia, damping, body force and internal force terms of the
/// momentum balance in total Lagrangian form, at the generalized-alpha time
/// levels.
///
/// @param[in] com_mod Simulation data: the current equation and domain, and
///   the time integration parameters.
/// @param[in] cep_mod Electrophysiology data.
/// @param[in] eNoN Number of nodes of the element.
/// @param[in] w Gauss weight times the Jacobian determinant of the element.
/// @param[in] N Shape functions at the Gauss point, eNoN.
/// @param[in] Nx Shape function gradients at the Gauss point, nsd x eNoN.
/// @param[in] al Nodal accelerations of the element, tDof x eNoN.
/// @param[in] yl Nodal velocities of the element, tDof x eNoN.
/// @param[in] dl Nodal displacements of the element, tDof x eNoN.
/// @param[in] bfl Nodal body forces of the element, nsd x eNoN.
/// @param[in] F0 Prestrain at the Gauss point: the deformation gradient from
///   the unloaded configuration to the mesh, so that F = (I + Grad u) F0
///   (see Prestrain). The identity when the mesh has no prestrain. If present,
///   then w, Nx, and fibers are pulled back to the unloaded configuration.
/// @param[in] fibers Fiber directions at the Gauss point, one per column.
/// @param[in] pS0l Nodal prestress of the element in Voigt order.
/// @param[out] pSl PK2 stress at the Gauss point in Voigt order.
/// @param[in] active_stress_evaluator Active stress model.
/// @param[in,out] lR Element residual, dof x eNoN.
/// @param[in,out] lK Element tangent, indexed lK(i*dof + j, a, b) for the
///   derivative of component i at node a with respect to component j at node b.
/// @param[in] recompute_visc Whether the viscosity model has to recompute the
///   quantities that depend on the shape function gradients.
void struct_3d(ComMod &com_mod, CepMod &cep_mod, const int eNoN,
               const double w, const Vector<double> &N, const Array<double> &Nx,
               const Array<double> &al, const Array<double> &yl,
               const Array<double> &dl, const Array<double> &bfl,
               const Eigen::Matrix<double, 3, 3> &F0,
               const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>> &fibers,
               const Array<double> &pS0l,
               Vector<double> &pSl, const ActiveStress::Evaluator &active_stress_evaluator,
               Array<double> &lR, Array3<double> &lK, const bool recompute_visc);

/// @brief Make a solid kernel's tangent that of a prestrain iteration.
///
/// A prestrain iteration composes the displacement @f$ u @f$ it reaches into
/// @f$ F_0 @f$ and starts again from rest, so the residual the next iteration
/// sees is that of the Cauchy stress on the mesh, which does not move, with
/// the mesh's fibers pulled back anew. The kernel's tangent differs from that
/// residual's in two ways, both corrected here.
///
/// The current configuration: the kernel's tangent carries the change of the
/// current volume and of the current shape function gradients with @f$ u @f$;
/// this subtracts that part,
/// @f[
///   w \left( \tau \nabla_x N_a \otimes \nabla_x N_b
///          - \tau \nabla_x N_b \otimes \nabla_x N_a \right),
/// @f]
/// with @f$ \tau = P F^T @f$ the Kirchhoff stress,
/// @f$ \nabla_x N = F^{-T} \nabla_0 N @f$, and @f$ w @f$ the stress-free
/// weight. Without it each iteration is an inexact Newton step whose error
/// grows with the stress, and the prestrain can stop converging once the
/// pseudo time step is large.
///
/// The fibers: composing the update pulls the mesh's fibers back anew (see
/// Prestrain), so the fibers the material law sees change with
/// @f$ u @f$. To first order, with @f$ G = F_0^{-1} \, \nabla u \, F_0 @f$, the
/// gradient on the mesh carried to the stress-free configuration, each fiber
/// @f$ f_0 @f$ changes by
/// @f[
///   \delta f_0 = -(I - \hat{f} \hat{f}^T) \, G \, f_0,
/// @f]
/// @f$ \hat{f} @f$ being its unit vector: it keeps unit length. The
/// kernel's tangent holds the fibers fixed; this adds
/// @f[
///   w \, F \, \delta S \, \nabla_0 N_a,
/// @f]
/// taking the gradient of the material law's @f$ S @f$ with respect to each
/// fiber by one-sided finite differences along the axes and applying it to
/// @f$ \delta f_0 @f$, so it costs nsd calls of the material law per fiber,
/// plus one for the baseline. A mesh without fibers pays nothing beyond the
/// first part.
///
/// @param[in] com_mod Simulation data.
/// @param[in] cep_mod Electrophysiology data.
/// @param[in] dmn The domain of the element.
/// @param[in] active_stress_evaluator the kernel's evaluator, updated for the element.
/// @param[in] N shape functions at the Gauss point, for the active tension.
/// @param[in] F0 the prestrain at the Gauss point.
/// @param[in] F the deformation gradient from the stress-free configuration,
/// (I + Grad u) F0.
/// @param[in] P the first Piola-Kirchhoff stress the kernel assembled, viscous and
/// prestress parts included.
/// @param[in] Nx0 shape function gradients on the stress-free configuration, nsd x eNoN.
/// @param[in] fibers the pulled-back fibers the material law saw, nsd x nFn.
/// @param[in] w_afu the weight times the stiffness scaling of lK.
/// @param[in] dof the degrees of freedom per node of lK.
/// @param[in,out] lK the element tangent, lK(i*dof + j, a, b); the correction
/// is added to it.
template <int nsd>
void prestrain_tangent_correction(const ComMod& com_mod, const CepMod& cep_mod, const dmnType& dmn,
    const ActiveStress::Evaluator& active_stress_evaluator, const Eigen::Ref<const Eigen::VectorXd>& N,
    const mat_fun::Matrix<nsd>& F0, const mat_fun::Matrix<nsd>& F, const mat_fun::Matrix<nsd>& P,
    const Eigen::Ref<const Eigen::Matrix<double, nsd, Eigen::Dynamic>>& Nx0,
    const mat_models::FiberRef<nsd>& fibers, double w_afu, int dof, Array3<double>& lK);
};

#endif
