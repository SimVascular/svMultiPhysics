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

#ifndef MAT_MODELS_H 
#define MAT_MODELS_H 

#include "Array.h"
#include "CepMod.h"
#include "ComMod.h"
#include "Tensor4.h"

#include "mat_fun.h"

#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"
#include "eigen3/unsupported/Eigen/CXX11/Tensor"

namespace mat_models {

void actv_strain(const ComMod& com_mod, const CepMod& cep_mod, const double gf, 
    const int nfd, const Array<double>& fl, Array<double>& Fa);

void cc_to_voigt(const int nsd, const Tensor4<double>& CC, Array<double>& Dm);

void voigt_to_cc(const int nsd, const Array<double>& Dm, Tensor4<double>& CC);

void compute_fib_stress(const ComMod& com_mod, const CepMod& cep_mod, const fibStrsType& Tfl, double& g);

void compute_pk2cc(const ComMod& com_mod, const CepMod& cep_mod, const dmnType& lDmn, const Array<double>& F, const int nfd,
    const Array<double>& fl, const double ya, Array<double>& S, Array<double>& Dm, double& Ja);

void compute_pk2cc_shlc(const ComMod& com_mod, const dmnType& lDmn, const int nfd, const Array<double>& fNa0,
    const Array<double>& gg_0, const Array<double>& gg_x, double& g33, Vector<double>& Sml, Array<double>& Dml);

void compute_pk2cc_shli(const ComMod& com_mod, const dmnType& lDmn, const int nfd, const Array<double>& fNa0,
    const Array<double>& gg_0, const Array<double>& gg_x, double& g33, Vector<double>& Sml, Array<double>& Dml);

void compute_tau(const ComMod& com_mod, const dmnType& lDmn, const double detF, const double Je, double& tauM, double& tauC);

void compute_svol_p(const ComMod& com_mod, const CepMod& cep_mod, const stModelType& stM, const double J, 
    double& p, double& pl);

void g_vol_pen(const ComMod& com_mod, const dmnType& lDmn, const double p, 
    double& ro, double& bt, double& dro, double& dbt, const double Ja);

void compute_visc_stress_potential(const double mu, const int eNoN, const Array<double>& Nx, const double vx, const double F,
                        Array<double>& Svis, Array3<double>& Kvis_u, Array3<double>& Kvis_v);

void compute_visc_stress_newtonian(const double mu, const int eNoN, const Array<double>& Nx, const Array<double>& vx, const Array<double>& F,
                        Array<double>& Svis, Array3<double>& Kvis_u, Array3<double>& Kvis_v);

void compute_visc_stress_and_tangent(const dmnType& lDmn, const int eNoN, const Array<double>& Nx, const  Array<double>& vx, const  Array<double>& F,
                        Array<double>& Svis, Array3<double>& Kvis_u, Array3<double>& Kvis_v);
};

#endif

