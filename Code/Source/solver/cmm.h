// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef CMM_H 
#define CMM_H 

#include "ComMod.h"
#include "Simulation.h"

/// @brief These subroutines implement the Coupled Momentum Method (CMM).
namespace cmm {

void cmm_3d(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx,  
    const Array<double>& al, const Array<double>& yl, const Array<double>& bfl, const Array<double>& Kxi,     
    Array<double>& lR, Array3<double>& lK);

void cmm_b(ComMod& com_mod, const faceType& lFa, const int e, const Array<double>& al, const Array<double>& dl, 
    const Array<double>& xl, const Array<double>& bfl, const Vector<double>& pS0l, const Vector<double>& vwp, 
    const Vector<int>& ptr);

void bcmmi(ComMod& com_mod, const int eNoN, const int idof, const double w, const Vector<double>& N, const Array<double>& Nxi,
    const Array<double>& xl, const Array<double>& tfl, Array<double>& lR);

void cmmi(ComMod& com_mod, const mshType& lM, const Array<double>& al, const Array<double>& dl, const Array<double>& xl,
    const Array<double>& bfl, const Array<double>& pS0l, const Vector<double>& vwp, const Vector<int>& ptr);

void cmm_mass(ComMod& com_mod, const double w, const Vector<double>& N, const Array<double>& al, 
    const Array<double>& bfl, const Vector<double>& vwp, Array<double>& lR, Array3<double>& lK);

void cmm_stiffness(ComMod& com_mod, const Array<double>& Nxi, const Array<double>& xl, const Array<double>& dl,                        
    const Vector<double>& pS0l, const Vector<double>& vwp, Vector<double>& pSl, Array<double>& lR, Array3<double>& lK);

void construct_cmm(ComMod& com_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Yg, const Array<double>& Dg);

};

#endif

