// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef USTRUCT_H 
#define USTRUCT_H 

#include "ComMod.h"

namespace ustruct {

void b_ustruct_2d(const ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, 
    const Array<double>& Nx, const Array<double>& dl, const Vector<double>& hl, const Vector<double>& nV, 
    Array<double>& lR, Array3<double>& lK, Array3<double>& lKd);

void b_ustruct_3d(const ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, 
    const Array<double>& Nx, const Array<double>& dl, const Vector<double>& hl, const Vector<double>& nV, 
    Array<double>& lR, Array3<double>& lK, Array3<double>& lKd);

bool constitutive_model_is_valid(consts::ConstitutiveModelType model);

void construct_usolid(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Yg, 
    const Array<double>& Dg);

int get_col_ptr(ComMod& com_mod, const int rowN, const int colN);

void ustruct_2d_c(ComMod& com_mod, CepMod& cep_mod, const bool vmsFlag, const int eNoNw, const int eNoNq,
    const double w, const double Je, const Vector<double>& Nw,  const Vector<double>& Nq,
    const Array<double>& Nwx, const Array<double>& Nqx, const Array<double>& al, const Array<double>& yl,
    const Array<double>& dl, const Array<double>& bfl, Array<double>& lR, Array3<double>& lK, 
    Array3<double>& lKd);

void ustruct_2d_m(ComMod& com_mod, CepMod& cep_mod, const bool vmsFlag, const int eNoNw, const int eNoNq,
    const int nFn, const double w, const double Je, const Vector<double>& Nw,  const Vector<double>& Nq,
    const Array<double>& Nwx, const Array<double>& al, const Array<double>& yl, const Array<double>& dl,
    const Array<double>& bfl, const Array<double>& fN, const Vector<double>& ya_l, Array<double>& lR,
    Array3<double>& lK, Array3<double>& lKd);

void ustruct_3d_c(ComMod& com_mod, CepMod& cep_mod, const bool vmsFlag, const int eNoNw, const int eNoNq,
    const double w, const double Je, const Vector<double>& Nw,  const Vector<double>& Nq,
    const Array<double>& Nwx, const Array<double>& Nqx, const Array<double>& al, const Array<double>& yl, 
    const Array<double>& dl, const Array<double>& bfl, Array<double>& lR, Array3<double>& lK, 
    Array3<double>& lKd);

void ustruct_3d_m(ComMod& com_mod, CepMod& cep_mod, const bool vmsFlag, const int eNoNw, const int eNoNq, 
    const int nFn, const double w, const double Je, const Vector<double>& Nw,  const Vector<double>& Nq, 
    const Array<double>& Nwx, const Array<double>& al, const Array<double>& yl, const Array<double>& dl, 
    const Array<double>& bfl, const Array<double>& fN, const Vector<double>& ya_l, Array<double>& lR, 
    Array3<double>& lK, Array3<double>& lKd);

void ustruct_do_assem(ComMod& com_mod, const int d, const Vector<int>& eqN, const Array3<double>& lKd, 
    const Array3<double>& lK, const Array<double>& lR);

void ustruct_r(ComMod& com_mod, const Array<double>& Yg);

};

#endif

