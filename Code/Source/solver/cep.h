// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef CEP_H 
#define CEP_H 

#include "ComMod.h"

namespace cep {

void b_cep(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const double h, Array<double>& lR);

void cep_1d(ComMod& com_mod, CepMod& cep_mod, const int eNoN, const int nFn, const double w,
    const Vector<double>& N, const Array<double>& Nx, const Array<double>& al, const Array<double>& yl,
    Array<double>& lR, Array3<double>& lK);

void cep_2d(ComMod& com_mod, CepMod& cep_mod, const int eNoN, const int nFn, const double w,
    const Vector<double>& N, const Array<double>& Nx, const Array<double>& al, const Array<double>& yl,
    const Array<double>& dl, const Array<double>& fN, Array<double>& lR, Array3<double>& lK);

void cep_3d(ComMod& com_mod, CepMod& cep_mod, const int eNoN, const int nFn, const double w,
    const Vector<double>& N, const Array<double>& Nx, const Array<double>& al, const Array<double>& yl,
    const Array<double>& dl, const Array<double>& fN, Array<double>& lR, Array3<double>& lK);

void construct_cep(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const Array<double>& Ag, 
    const Array<double>& Yg, const Array<double>& Dg);

};

#endif

