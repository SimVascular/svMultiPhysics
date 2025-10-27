// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef HEATF_H 
#define HEATF_H 

#include "ComMod.h"

namespace heatf {

void b_heatf(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Vector<double>& y,
    const double h, const Vector<double>& nV, Array<double>& lR, Array3<double>& lK);

void heatf_2d(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx,
    const Array<double>& al, const Array<double>& yl, const Array<double>& ksix, Array<double>& lR, Array3<double>& lK);

void heatf_3d(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx,
    const Array<double>& al, const Array<double>& yl, const Array<double>& ksix, Array<double>& lR, Array3<double>& lK);

void construct_heatf(ComMod& com_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Yg);

};

#endif

