// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "fsils.hpp"
#include "CmMod.h"

#ifndef FSI_LINEAR_SOLVER_COMMU_H 
#define FSI_LINEAR_SOLVER_COMMU_H 

namespace fsi_linear_solver {

void fsils_commu_create(FSILS_commuType& commu, cm_mod::MpiCommWorldType commi);

void fsils_commu_free(FSILS_commuType& commu);

};

#endif


