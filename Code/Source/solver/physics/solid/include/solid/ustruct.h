// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef USTRUCT_H 
#define USTRUCT_H 

#include "ActiveStress.h"
#include "ComMod.h"
#include "SolutionStates.h"

namespace ustruct {

void b_ustruct_2d(const ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, 
    const Array<double>& Nx, const Array<double>& dl, const Vector<double>& hl, const Vector<double>& nV, 
    Array<double>& lR, Array3<double>& lK, Array3<double>& lKd);

void b_ustruct_3d(const ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, 
    const Array<double>& Nx, const Array<double>& dl, const Vector<double>& hl, const Vector<double>& nV, 
    Array<double>& lR, Array3<double>& lK, Array3<double>& lKd);

bool constitutive_model_is_valid(consts::ConstitutiveModelType model);

void construct_usolid(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const SolutionStates& solutions);

int get_col_ptr(ComMod& com_mod, const int rowN, const int colN);

/// @brief Contribution of one Gauss point to the element residual and
/// tangent of the continuity (pressure) equation of ustruct, the mixed
/// velocity-pressure nonlinear solid, in two dimensions.
///
/// @param[in] com_mod Simulation data: the current equation and domain, and
///   the time integration parameters.
/// @param[in] cep_mod Electrophysiology data, for the material laws that
///   depend on it.
/// @param[in] vmsFlag Whether the variational multiscale stabilization terms
///   are added.
/// @param[in] eNoNw Number of nodes of the element's velocity function space.
/// @param[in] eNoNq Number of nodes of the element's pressure function space.
/// @param[in] w Gauss weight times the Jacobian determinant of the element.
/// @param[in] Je Jacobian determinant of the element, for the stabilization
///   parameters.
/// @param[in] Nw Velocity shape functions at the Gauss point, eNoNw.
/// @param[in] Nq Pressure shape functions at the Gauss point, eNoNq.
/// @param[in] Nwx Gradients of the velocity shape functions, nsd x eNoNw.
/// @param[in] Nqx Gradients of the pressure shape functions, nsd x eNoNq.
/// @param[in] al Nodal accelerations of the element, tDof x eNoNw.
/// @param[in] yl Nodal velocities and pressures of the element, tDof x eNoNw.
/// @param[in] dl Nodal displacements of the element, tDof x eNoNw.
/// @param[in] bfl Nodal body forces of the element, nsd x eNoNw.
/// @param[in] F0 Prestrain deformation gradient at the Gauss point, from the
///   stress-free configuration to the mesh (see Prestrain), with which the
///   deformation gradient is @f$ F = (I + \nabla u) F_0 @f$. A prestrain is
///   not currently supported for ustruct (see prestrain::init), so this is
///   always the identity.
/// @param[in,out] lR Element residual, dof x eNoNw; the contribution is added.
/// @param[in,out] lK Element tangent with respect to the velocity and
///   pressure, indexed lK(i*dof + j, a, b); the contribution is added.
/// @param[in,out] lKd Element tangent with respect to the displacement; the
///   contribution is added.
void ustruct_2d_c(ComMod& com_mod, CepMod& cep_mod, const bool vmsFlag, const int eNoNw, const int eNoNq,
    const double w, const double Je, const Vector<double>& Nw,  const Vector<double>& Nq,
    const Array<double>& Nwx, const Array<double>& Nqx, const Array<double>& al, const Array<double>& yl,
    const Array<double>& dl, const Array<double>& bfl,
    const Eigen::Matrix<double, 2, 2>& F0, Array<double>& lR, Array3<double>& lK, 
    Array3<double>& lKd);

/// @brief Contribution of one Gauss point to the element residual and
/// tangent of the momentum equation of ustruct, the mixed velocity-pressure
/// nonlinear solid, in two dimensions.
///
/// @param[in] com_mod Simulation data: the current equation and domain, and
///   the time integration parameters.
/// @param[in] cep_mod Electrophysiology data, for the material laws that
///   depend on it.
/// @param[in] vmsFlag Whether the variational multiscale stabilization terms
///   are added.
/// @param[in] eNoNw Number of nodes of the element's velocity function space.
/// @param[in] eNoNq Number of nodes of the element's pressure function space.
/// @param[in] w Gauss weight times the Jacobian determinant of the element.
/// @param[in] Je Jacobian determinant of the element, for the stabilization
///   parameters.
/// @param[in] Nw Velocity shape functions at the Gauss point, eNoNw.
/// @param[in] Nq Pressure shape functions at the Gauss point, eNoNq.
/// @param[in] Nwx Gradients of the velocity shape functions, nsd x eNoNw.
/// @param[in] al Nodal accelerations of the element, tDof x eNoNw.
/// @param[in] yl Nodal velocities and pressures of the element, tDof x eNoNw.
/// @param[in] dl Nodal displacements of the element, tDof x eNoNw.
/// @param[in] bfl Nodal body forces of the element, nsd x eNoNw.
/// @param[in] F0 Prestrain deformation gradient at the Gauss point, from the
///   stress-free configuration to the mesh (see Prestrain), with which the
///   deformation gradient is @f$ F = (I + \nabla u) F_0 @f$. A prestrain is
///   not currently supported for ustruct (see prestrain::init), so this is
///   always the identity.
/// @param[in] fibers Fiber directions at the Gauss point, one per column.
/// @param[in] active_stress_evaluator Active stress model, updated for the
///   element.
/// @param[in,out] lR Element residual, dof x eNoNw; the contribution is added.
/// @param[in,out] lK Element tangent with respect to the velocity and
///   pressure, indexed lK(i*dof + j, a, b); the contribution is added.
/// @param[in,out] lKd Element tangent with respect to the displacement; the
///   contribution is added.
/// @param[in] recompute_visc Whether the viscosity model has to recompute the
///   quantities that depend on the shape function gradients; false when they
///   are unchanged from the element's previous Gauss point.
void ustruct_2d_m(ComMod &com_mod, CepMod &cep_mod, const bool vmsFlag,
                  const int eNoNw, const int eNoNq,
                  const double w, const double Je, const Vector<double> &Nw,
                  const Vector<double> &Nq, const Array<double> &Nwx,
                  const Array<double> &al, const Array<double> &yl,
                  const Array<double> &dl, const Array<double> &bfl,
                  const Eigen::Matrix<double, 2, 2> &F0,
                  const Eigen::Ref<const Eigen::Matrix<double, 2, Eigen::Dynamic>> &fibers,
                  const ActiveStress::Evaluator &active_stress_evaluator,
                  Array<double> &lR, Array3<double> &lK, Array3<double> &lKd,
                  const bool recompute_visc);

/// @brief Three-dimensional counterpart of ustruct_2d_c.
/// @copydetails ustruct_2d_c
void ustruct_3d_c(ComMod& com_mod, CepMod& cep_mod, const bool vmsFlag, const int eNoNw, const int eNoNq,
    const double w, const double Je, const Vector<double>& Nw,  const Vector<double>& Nq,
    const Array<double>& Nwx, const Array<double>& Nqx, const Array<double>& al, const Array<double>& yl, 
    const Array<double>& dl, const Array<double>& bfl, const Eigen::Matrix<double, 3, 3>& F0,
    Array<double>& lR, Array3<double>& lK, 
    Array3<double>& lKd);

/// @brief Three-dimensional counterpart of ustruct_2d_m.
/// @copydetails ustruct_2d_m
void ustruct_3d_m(ComMod &com_mod, CepMod &cep_mod, const bool vmsFlag,
                  const int eNoNw, const int eNoNq,
                  const double w, const double Je, const Vector<double> &Nw,
                  const Vector<double> &Nq, const Array<double> &Nwx,
                  const Array<double> &al, const Array<double> &yl,
                  const Array<double> &dl, const Array<double> &bfl,
                  const Eigen::Matrix<double, 3, 3> &F0,
                  const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>> &fibers,
                  const ActiveStress::Evaluator &active_stress_evaluator,
                  Array<double> &lR, Array3<double> &lK, Array3<double> &lKd,
                  const bool recompute_visc);

void ustruct_do_assem(ComMod& com_mod, const int d, const Vector<int>& eqN, const Array3<double>& lKd, 
    const Array3<double>& lK, const Array<double>& lR);

void ustruct_r(ComMod& com_mod, const SolutionStates& solutions);

};

#endif
