// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "fils_struct.hpp"

namespace cgrad {

using namespace fsi_linear_solver;

void cgrad_v(FSILS_lhsType& lhs, FSILS_subLsType& ls, const int dof, const Array<double>& K, Array<double>& R);

void cgrad_s(FSILS_lhsType& lhs, FSILS_subLsType& ls, const Vector<double>& K, Vector<double>& R);

void schur(FSILS_lhsType& lhs, FSILS_subLsType& ls, const int dof, const Array<double>& D,
    const Array<double>& G, const Vector<double>& L, Vector<double>& R);

};
