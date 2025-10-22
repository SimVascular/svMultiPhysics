// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef STOKES_H 
#define STOKES_H 

#include "ComMod.h"

#include "consts.h"

namespace stokes {

void construct_stokes(ComMod& com_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Yg);

void stokes_2d_c(ComMod& com_mod, const int lStab, const int eNoNw, const int eNoNq, const double w,
    const Array<double>& ksix, const Vector<double>& Nw, const Vector<double>& Nq, const Array<double>& Nwx,
    const Array<double>& Nqx, const Array<double>& al, const Array<double>& yl, const Array<double>& bfl, 
    Array<double>& lR, Array3<double>& lK);

void stokes_2d_m(ComMod& com_mod, const int eNoNw, const int eNoNq, const double w,
    const Vector<double>& Nw, const Vector<double>& Nq, const Array<double>& Nwx, const Array<double>& al,
    const Array<double>& yl, const Array<double>& bfl, Array<double>& lR, Array3<double>& lK);

void stokes_3d_c(ComMod& com_mod, const int lStab, const int eNoNw, const int eNoNq, const double w,
    const Array<double>& ksix, const Vector<double>& Nw, const Vector<double>& Nq, const Array<double>& Nwx,
    const Array<double>& Nqx, const Array<double>& al, const Array<double>& yl, const Array<double>& bfl, 
    Array<double>& lR, Array3<double>& lK);

void stokes_3d_m(ComMod& com_mod, const int eNoNw, const int eNoNq, const double w, 
    const Vector<double>& Nw, const Vector<double>& Nq, const Array<double>& Nwx,
    const Array<double>& al, const Array<double>& yl, const Array<double>& bfl, 
    Array<double>& lR, Array3<double>& lK);

};

#endif

