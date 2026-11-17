// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "fils_struct.hpp"

#include "consts.h"
#include <array>

#ifndef FILS_API_H
#define FILS_API_H

namespace fsi_linear_solver {

/// @todo[michelebucelli] Most if not all the functions below would fit well as
///   methods of the FSILS_lhsType class.

/// @brief Create the face of a Dirichlet or Neumann boundary condition in the
/// linear system.
///
/// 'Val' is stored in free_component_mask for a Dirichlet face and in
/// unit_pressure_load for a Neumann face, summed over the processes that share
/// the face. For a Neumann face, flowrate_derivative is set to
/// unit_pressure_load.
///
/// @param[in,out] lhs The left-hand side of the linear system. Only
///   lhs.face[faIn] is modified.
/// @param[in] faIn The index of the face in lhs.face.
/// @param[in] BC_type The type of the boundary condition.
/// @param[in] gNodes The mesh indices of the nodes of the face.
/// @param[in] Val The values at the nodes of the face, nodes ordered as in
///   'gNodes'. Its number of rows is the number of components per node, and
///   its number of columns is the number of nodes of the face on this process.
void fsils_bc_create(FSILS_lhsType &lhs, int faIn, BcType BC_type,
                     const Vector<int> &gNodes, const Array<double> &Val);

/// @brief Create the face of a Dirichlet or Neumann boundary condition in the
/// linear system, with all its nodal values equal to zero.
///
/// For a Dirichlet face, this constrains every component of every node.
///
/// @param[in,out] lhs The left-hand side of the linear system. Only
///   lhs.face[faIn] is modified.
/// @param[in] faIn The index of the face in lhs.face.
/// @param[in] nNo The number of nodes of the face on this process.
/// @param[in] dof The number of components per node.
/// @param[in] BC_type The type of the boundary condition.
/// @param[in] gNodes The mesh indices of the nodes of the face.
void fsils_bc_create(FSILS_lhsType& lhs, int faIn, int nNo, int dof, BcType BC_type, const Vector<int>& gNodes);

void fsils_bc_free(FSILS_lhsType& lhs, int faIn);

/// @brief Update the two vectors of the coupling operator of a Neumann face.
///
/// Both are summed over the processes that share the face. The first replaces
/// lhs.face[faIn].unit_pressure_load. The second replaces
/// lhs.face[faIn].flowrate_derivative, at the nodes of the face that this
/// process owns; entries appended by fsils_bc_append_flowrate_derivative()
/// are dropped.
///
/// @param[in,out] lhs The left-hand side of the linear system. Only
///   lhs.face[faIn] is modified.
/// @param[in] faIn The index of the face in lhs.face.
/// @param[in] unit_pressure_load The integral over the face of N_a * n, for
///   each face node a (N_a: shape function, n: normal). Its size and node
///   order are those of the values that the face was created with.
/// @param[in] flowrate_derivative The derivative of the flow rate through the
///   face with respect to the velocity of each face node, with the same size
///   and node order.
void fsils_bc_update(FSILS_lhsType &lhs, int faIn,
                     const Array<double> &unit_pressure_load,
                     const Array<double> &flowrate_derivative);

/// @brief Add values to the flow rate derivative of a Neumann face.
///
/// The flow rate derivative is a list of (node, value) entries that is only
/// used in scalar products, where two entries of the same node act as one
/// entry holding their sum. So the values are appended as new entries, without
/// looking for an existing entry of the same node.
///
/// Only the nodes that this process owns are appended. The values are not
/// summed over processes: each must be the complete value at its node.
///
/// @param[in] lhs The left-hand side of the linear system that the face
///   belongs to.
/// @param[in,out] face The face, whose flowrate_derivative_node_indices,
///   flowrate_derivative and scaled_flowrate_derivative are extended.
/// @param[in] gNodes The mesh indices of the nodes of the entries.
/// @param[in] values The values at those nodes. Its number of rows is the
///   number of components per node of the face, and its number of columns is
///   the number of nodes.
void fsils_bc_append_flowrate_derivative(const FSILS_lhsType &lhs,
                                         FSILS_faceType &face,
                                         const Vector<int> &gNodes,
                                         const Array<double> &values);

void fsils_commus(const FSILS_lhsType& lhs, Vector<double>& R); 

void fsils_commuv(const FSILS_lhsType& lhs, const int dof, Array<double>& R);

double fsils_cpu_t();

void fsils_ls_create(FSILS_lsType& ls, LinearSolverType LS_type, double relTol = consts::double_inf, 
  double absTol = consts::double_inf, int maxItr = consts::int_inf, int dimKry = consts::int_inf, 
  std::array<double,2> relTolIn = {}, std::array<double,2> absTolIn = {}, std::array<int,2> maxItrIn = {});

void fsils_solve(FSILS_lhsType& lhs, FSILS_lsType& ls, const int dof, Array<double>& Ri, Array<double>& Val,
    const consts::PreconditionerType prec, const Vector<int>& incL, const Vector<double>& res);

};

#endif


