// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "fils_struct.hpp"

namespace add_bc_mul {

using namespace fsi_linear_solver;

/// @brief Add the contribution of the coupled Neumann faces to a
/// matrix-vector product.
///
/// A coupled face contributes the operator
///
///   coef * scaled_unit_pressure_load * scaled_flowrate_derivative^T
///
/// where coef is face.res when adding the contribution to the product, and
/// -face.res / (1 + face.res * face.nS) when computing the contribution to
/// the preconditioner.
///
/// For reference, see
/// Moghadam et al. 2013 eq. 27 (https://doi.org/10.1016/j.jcp.2012.07.035) and
/// Moghadam et al. 2013b (https://doi.org/10.1007/s00466-013-0868-1).
///
/// @param[in] lhs The left-hand side of the linear system, whose faces are
///   read.
/// @param[in] op_Type Whether to compute the contribution to the product or
///   to the preconditioner.
/// @param[in] X The vector that the operator is applied to (components per
///   node x lhs.nNo).
/// @param[in,out] Y The vector that the result is added to, of the same size
///   as 'X'.
void add_bc_mul(const FSILS_lhsType& lhs, const BcopType op_Type, const Array<double>& X, Array<double>& Y);

};
