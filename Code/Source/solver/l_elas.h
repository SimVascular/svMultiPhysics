// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef L_ELAS_H 
#define L_ELAS_H 

#include "ComMod.h"

namespace l_elas {

void b_l_elas(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const double h, 
    const Vector<double>& nV, Array<double>& lR);

void construct_l_elas(ComMod& com_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Dg);

void l_elas_2d(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, 
    const Array<double>& Nx, const Array<double>& al, const Array<double>& dl, const Array<double>& bfl, 
    const Array<double>& pS0l, Vector<double>& pSl, Array<double>& lR, Array3<double>& lK);

void l_elas_3d(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, 
    const Array<double>& Nx, const Array<double>& al, const Array<double>& dl, const Array<double>& bfl, 
    const Array<double>& pS0l, Vector<double>& pSl, Array<double>& lR, Array3<double>& lK);

};

#endif

