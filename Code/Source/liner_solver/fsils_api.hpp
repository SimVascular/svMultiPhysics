// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "fils_struct.hpp"

#include "consts.h"
#include <array>

#ifndef FILS_API_H
#define FILS_API_H

namespace fsi_linear_solver {

void fsils_bc_create(FSILS_lhsType& lhs, int faIn, int nNo, int dof, BcType BC_type, const Vector<int>& gNodes,
    const Array<double>& Val);

void fsils_bc_create(FSILS_lhsType& lhs, int faIn, int nNo, int dof, BcType BC_type, const Vector<int>& gNodes);

void fsils_bc_free(FSILS_lhsType& lhs, int faIn);

void fsils_bc_update(FSILS_lhsType& lhs, int faIn, int nNo, int dof, const Array<double>& Val);
    
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


