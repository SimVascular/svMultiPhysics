// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef STRUCT_H 
#define STRUCT_H 

#include "ActiveStress.h"
#include "ComMod.h"
#include "SolutionStates.h"

namespace struct_ns {

void b_struct_2d(const ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, 
    const Array<double>& Nx, const Array<double>& dl, const Vector<double>& hl, const Vector<double>& nV, 
    Array<double>& lR, Array3<double>& lK);

void b_struct_3d(const ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, 
    const Array<double>& Nx, const Array<double>& dl, const Vector<double>& hl, const Vector<double>& nV, 
    Array<double>& lR, Array3<double>& lK);

void construct_dsolid(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const SolutionStates& solutions);

/// @brief Contribution of one Gauss point to the element residual and
/// tangent of the struct equation (displacement-based nonlinear solid), in two
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
/// @param[in] F0 Prestrain deformation gradient at the Gauss point.
/// @param[in] fibers Fiber directions at the Gauss point, one per column.
/// @param[in] pS0l Nodal prestress of the element in Voigt order.
/// @param[out] pSl PK2 stress at the Gauss point in Voigt order.
/// @param[in] active_stress_evaluator Active stress model.
/// @param[in,out] lR Element residual, dof x eNoN.
/// @param[in,out] lK Element tangent, indexed lK(i*dof + j, a, b) for the
///   derivative of component i at node a with respect to component j at node b.
/// @param[in] recompute_visc Whether the viscosity model has to recompute the
///   quantities that depend on the shape function gradients within an element.
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
/// tangent of the struct equation (displacement-based nonlinear solid), in two
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
/// @param[in] F0 Prestrain deformation gradient at the Gauss point.
/// @param[in] fibers Fiber directions at the Gauss point, one per column.
/// @param[in] pS0l Nodal prestress of the element in Voigt order.
/// @param[out] pSl PK2 stress at the Gauss point in Voigt order.
/// @param[in] active_stress_evaluator Active stress model.
/// @param[in,out] lR Element residual, dof x eNoN.
/// @param[in,out] lK Element tangent, indexed lK(i*dof + j, a, b) for the
///   derivative of component i at node a with respect to component j at node b.
/// @param[in] recompute_visc Whether the viscosity model has to recompute the
///   quantities that depend on the shape function gradients within an element.
void struct_3d(ComMod &com_mod, CepMod &cep_mod, const int eNoN,
               const double w, const Vector<double> &N, const Array<double> &Nx,
               const Array<double> &al, const Array<double> &yl,
               const Array<double> &dl, const Array<double> &bfl,
               const Eigen::Matrix<double, 3, 3> &F0,
               const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>> &fibers,
               const Array<double> &pS0l,
               Vector<double> &pSl, const ActiveStress::Evaluator &active_stress_evaluator,
               Array<double> &lR, Array3<double> &lK, const bool recompute_visc);
};

#endif
