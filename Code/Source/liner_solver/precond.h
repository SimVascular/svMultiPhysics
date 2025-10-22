// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "fils_struct.hpp"

namespace precond {

void pos_mul(const Array<int>& rowPtr, const Vector<int>& colPtr, const int nNo, const int nnz, const int dof, Array<double>& Val, const Array<double>& W);

void precond_diag(fsi_linear_solver::FSILS_lhsType& lhs, const Array<int>& rowPtr, const Vector<int>& colPtr, const Vector<int>& diagPtr, 
    const int dof, Array<double>& Val, Array<double>& R, Array<double>& W);

void precond_rcs(fsi_linear_solver::FSILS_lhsType& lhs, const Array<int>& rowPtr, const Vector<int>& colPtr,
    const Vector<int>& diagPtr, const int dof, Array<double>& Val, Array<double>& R, Array<double>& W1, Array<double>& W2);

void pre_mul(const Array<int>& rowPtr, const int nNo, const int nnz, const int dof, Array<double>& Val, const Array<double>& W);

};
