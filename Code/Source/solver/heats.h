// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef HEATS_H 
#define HEATS_H 

#include "ComMod.h"

namespace heats {

void b_heats(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const double h, Array<double>& lR);

void construct_heats(ComMod& com_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Dg);

void heats_2d(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx, 
    const Array<double>& al, const Array<double>& yl, Array<double>& lR, Array3<double>& lK);

void heats_3d(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx,
    const Array<double>& al, const Array<double>& yl, Array<double>& lR, Array3<double>& lK);

};

#endif

