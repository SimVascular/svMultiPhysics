// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef LS_H 
#define LS_H 

#include "ComMod.h"
#include "Simulation.h"

namespace ls_ns {

void ls_alloc(ComMod& com_mod, eqType& lEq);

void ls_solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res);

//void init_dir_and_coupneu_bc_petsc(ComMod& com_mod, const Vector<int>& incL, const Vector<double>& res);

//void initialize_petsc(ComMod& com_mod);

};

#endif

