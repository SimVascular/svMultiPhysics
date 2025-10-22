// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef BF_H 
#define BF_H 

#include "ComMod.h"
#include "Simulation.h"

namespace bf {

void bf_construct(ComMod& com_mod, const mshType& lM, const int e, const int eNoN, const int idof, Array<double>& xl, 
    const Array<double>& dl, const Array<double>& bfl, const Vector<int>& ptr);

void set_bf(ComMod& com_mode, const Array<double>& Dg);

void set_bf_l(ComMod& com_mod, bfType& lBf, mshType& lM, const Array<double>& Dg);

};

#endif

