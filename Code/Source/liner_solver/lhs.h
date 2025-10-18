// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "fsils.hpp"
#include "CmMod.h"

#include "Vector.h"

#ifndef FSI_LINEAR_SOLVER_LHS_H 
#define FSI_LINEAR_SOLVER_LHS_H 

namespace fsi_linear_solver {

void fsils_lhs_create(FSILS_lhsType& lhs, FSILS_commuType& commu, int gnNo, int nNo, int nnz, Vector<int>& gNodes,  
       Vector<int>&rowPtr, Vector<int>&colPtr, int nFaces);

void fsils_lhs_free(FSILS_lhsType& lhs);

};

#endif


