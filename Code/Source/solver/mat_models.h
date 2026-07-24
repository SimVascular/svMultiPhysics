// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

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

/**
 * @brief Compute 2nd Piola-Kirchhoff stress and material stiffness tensors
 * including both dilational and isochoric components.
 *
 * Reproduces the Fortran 'GETPK2CC' subroutine.
 *
 * @param[in] com_mod Object containing global common variables.
 * @param[in] cep_mod Object containing electrophysiology-specific common
 *   variables.
 * @param[in] lDmn Domain object.
 * @param[in] F Deformation gradient tensor.
 * @param[in] nfd Number of fiber directions.
 * @param[in] fl Fiber directions.
 * @param[in] ya_f Active tension along the fiber direction.
 * @param[in] ya_s Active tension along the sheet direction.
 * @param[in] ya_n Active tension along the sheet-normal direction.
 * @param[out] S 2nd Piola-Kirchhoff stress tensor (modified in place).
 * @param[out] Dm Material stiffness tensor (modified in place).
 * @param[out] Ja Jacobian for active strain
 *
 * @return None, but modifies S, Dm, and Ja in place.
 */
void compute_pk2cc(const ComMod &com_mod, const CepMod &cep_mod,
                   const dmnType &lDmn, const Array<double> &F, const int nfd,
                   const Array<double> &fl, const double ya_f,
                   const double ya_s, const double ya_n, Array<double> &S,
                   Array<double> &Dm, double &Ja);

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

