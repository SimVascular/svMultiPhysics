// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SHELLS_H 
#define SHELLS_H 

#include "ComMod.h"
#include "Simulation.h"

namespace shells {

void construct_shell(ComMod& com_mod, const mshType& lM, const Array<double>& Ag,
    const Array<double>& Yg, const Array<double>& Dg);

void shell_3d(ComMod& com_mod, const mshType& lM, const int g, const int eNoN,
    const int nFn, const Array<double>& fN, const Array<double>& al, const Array<double>& yl, 
    const Array<double>& dl, const Array<double>& xl, const Array<double>& bfl, 
    Array<double>& lR, Array3<double>& lK);

void shell_bend_cst(ComMod& com_mod, const mshType& lM, const int e, const Vector<int>& ptr,
    Array<double>& x0, Array<double>& xc, double bb_0[2][2], double bb_x[2][2],
    Array3<double>& Bb, const bool vflag);

void shell_bf(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx,
    const Array<double>& dl, const Array<double>& xl, const Array<double>& tfl, Array<double>& lR, Array3<double>& lK);

void shell_cst(ComMod& com_mod, const mshType& lM, const int e, const int eNoN, const int nFn, const Array<double>& fN,
    const Array<double>& al, const Array<double>& yl, const Array<double>& dl, const Array<double>& xl,
    const Array<double>& bfl, const Vector<int>& ptr);

void shell_fp(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx, 
    const Array<double>& dl, const Array<double>& xl, const Array<double>& tfl, Array<double>& lR, Array3<double>& lK);

void shl_strs_res(const ComMod& com_mod, const dmnType& lDmn, const int nFn, const Array<double>& fNa0,
    const double aa_0[2][2], const double aa_x[2][2], const double bb_0[2][2], const double bb_x[2][2],
    double& lam3, Array<double>& Sm, Array3<double>& Dm);

};

#endif

