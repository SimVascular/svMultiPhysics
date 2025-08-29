/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef FLUID_H 
#define FLUID_H 

#include "ComMod.h"

#include "consts.h"

namespace fluid {

void b_fluid(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Vector<double>& y,   
    const double h, const Vector<double>& nV, Array<double>& lR, Array3<double>& lK);

void bw_fluid_2d(ComMod& com_mod, const int eNoNw, const int eNoNq, const double w, const Vector<double>& Nw, 
    const Vector<double>& Nq, const Array<double>& Nwx, const Array<double>& yl, const Vector<double>& ub, 
    const Vector<double>& nV, const Vector<double>& tauB, Array<double>& lR, Array3<double>& lK);

void bw_fluid_3d(ComMod& com_mod, const int eNoNw, const int eNoNq, const double w, const Vector<double>& Nw, 
    const Vector<double>& Nq, const Array<double>& Nwx, const Array<double>& yl, const Vector<double>& ub, 
    const Vector<double>& nV, const Vector<double>& tauB, Array<double>& lR, Array3<double>& lK);

void construct_fluid(ComMod& com_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Yg);

void fluid_2d_c(ComMod& com_mod, const int vmsFlag, const int eNoNw, const int eNoNq, const double w, const Array<double>& Kxi, 
    const Vector<double>& Nw, const Vector<double>& Nq, const Array<double>& Nwx, const Array<double>& Nqx, 
    const Array<double>& Nwxx, const Array<double>& al, const Array<double>& yl, const Array<double>& bfl, 
    Array<double>& lR, Array3<double>& lK, double K_inverse_darcy_permeabilityx);

void fluid_2d_m(ComMod& com_mod, const int vmsFlag, const int eNoNw, const int eNoNq, const double w, const Array<double>& Kxi, 
    const Vector<double>& Nw, const Vector<double>& Nq, const Array<double>& Nwx, const Array<double>& Nqx, 
    const Array<double>& Nwxx, const Array<double>& al, const Array<double>& yl, const Array<double>& bfl, 
    Array<double>& lR, Array3<double>& lK, double K_inverse_darcy_permeability);

void fluid_3d_c(ComMod& com_mod, const int vmsFlag, const int eNoNw, const int eNoNq, const double w, const Array<double>& Kxi, 
    const Vector<double>& Nw, const Vector<double>& Nq, const Array<double>& Nwx, const Array<double>& Nqx, 
    const Array<double>& Nwxx, const Array<double>& al, const Array<double>& yl, const Array<double>& bfl, 
    Array<double>& lR, Array3<double>& lK, double K_inverse_darcy_permeability, double DDir=0.0);

void fluid_3d_m(ComMod& com_mod, const int vmsFlag, const int eNoNw, const int eNoNq, const double w, const Array<double>& Kxi, 
    const Vector<double>& Nw, const Vector<double>& Nq, const Array<double>& Nwx, const Array<double>& Nqx, 
    const Array<double>& Nwxx, const Array<double>& al, const Array<double>& yl, const Array<double>& bfl, 
    Array<double>& lR, Array3<double>& lK, double K_inverse_darcy_permeability, double DDir=0.0);

void get_viscosity(const ComMod& com_mod, const dmnType& lDmn, double& gamma, double& mu, double& mu_s, double& mu_x);

};

#endif

