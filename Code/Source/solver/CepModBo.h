// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef CEP_MOD_BO_H 
#define CEP_MOD_BO_H 

#include "Array.h"
#include "Vector.h"
#include <array>

using BoModelParam = Vector<double>;

/// @brief This module defines data structures for Bueno-Orovio cellular
/// activation model for cardiac electrophysiology.
///
/// The classes defined here duplicate the data structures in the Fortran BOMOD module defined 
/// in CEPMOD_BO.f and PARAMS_BO.f files. 
class CepModBo
{
  public:
    CepModBo();
    ~CepModBo();

    // Scaling factors
    /// Voltage scaling
    double Vscale = 85.70;
    /// Time scaling
    double Tscale = 1.0;
    /// Voltage offset parameter
    double Voffset = -84.0;

    /// Model parameters (epi, endo, myo)
    ///
    /// \todo [TODO:DaveP] these guys should be maps map<int,double>.
    ///
    BoModelParam u_o = {0.0, 0.0, 0.0};
    BoModelParam u_u = {1.550, 1.56, 1.61};
    BoModelParam theta_v = {0.30, 0.3, 0.3};
    BoModelParam theta_w = {0.130, 0.13, 0.13};
    BoModelParam thetam_v = {6.E-3, 0.2, 0.1};
    BoModelParam theta_o = {6.E-3, 6.E-3, 5.E-3};
    BoModelParam taum_v1 = {60.0, 75., 80.};
    BoModelParam taum_v2 = {1.15E3, 10., 1.4506};
    BoModelParam taup_v = {1.45060, 1.4506, 1.4506};
    BoModelParam taum_w1 = {60.0, 6., 70.};
    BoModelParam taum_w2 = {15.0, 140., 8.};
    BoModelParam km_w = {65.0, 200., 200.};
    BoModelParam um_w = {3.E-2, 1.6E-2, 1.6E-2};
    BoModelParam taup_w = {200.0, 280., 280.};
    BoModelParam tau_fi = {0.110, 0.1, 0.078};
    BoModelParam tau_o1 = {400.0, 470., 410.};
    BoModelParam tau_o2 = {6.0, 6., 7.};
    BoModelParam tau_so1 = {30.01810, 40., 91.};
    BoModelParam tau_so2 = {0.99570, 1.2, 0.8};
    BoModelParam k_so = {2.04580, 2., 2.1};
    BoModelParam u_so = {0.650, 0.65, 0.6};
    BoModelParam tau_s1 = {2.73420, 2.7342, 2.7342};
    BoModelParam tau_s2 = {16.0, 2., 2.};
    BoModelParam k_s = {2.09940, 2.0994, 2.0994};
    BoModelParam u_s = {0.90870, 0.9087, 0.9087};
    BoModelParam tau_si = {1.88750, 2.9013, 3.3849};
    BoModelParam tau_winf = {7.E-2, 2.73E-2, 1.E-2};
    BoModelParam ws_inf = {0.940, 0.78, 0.5};

    // Electromechanics coupling parameters: active stress model
    /// Resting voltage (mV)
    double Vrest = -84.0;
    /// Critical voltage (mV)
    double Vcrit = -30.0;
    /// Saturation potential
    double eta_T = 5.E-3;
    /// Minimum activation (ms^{-1})
    double eps_0 = 0.10;
    /// Maximum activation (ms^{-1})
    double eps_i = 1.0;
    /// Transition rate (mV^{-1})
    double xi_T = 1.0;

    // Electromechanics coupling parameters: active strain model
    /// Active force of sacromere (-mM^{-2})
    double alFa = -4.E+6;
    /// Resting Ca concentration (mM) := slow inward current variable (s)
    double c0 = 2.155E-4;
    /// Viscous-type constant (ms-mM^{-2})
    double mu_C = 5.E+6;

    // Force-length relationship parameters
    /// Initial length of sacromeres (um)
    double SL0 = 1.950;
    /// Min. length of sacromeres (um)
    double SLmin = 1.70;
    /// Max. length of sacromeres (um)
    double SLmax = 2.60;
    /// Fourier coefficients
    double f0  = -4333.6183355821190;
    double fc1 =  2570.3953553521950;
    double fs1 = -2051.8272789919760;
    double fc2 =  1329.536116891330;
    double fs2 =  302.2167845582220;
    double fc3 =  104.9437703051160;
    double fs3 =  218.3751742294220;

    /// Cm: Cell capacitance per unit surface area
    double Cm  = 1.0;
    /// sV: Surface to volume ratio
    double sV  = 1.0;
    /// rho: Cellular resistivity
    double rho = 1.0;

    void actv_strn(const double c, const double I4f, const double dt, double& gf);
    void actv_strs(const double X, const double dt, double& Tact, double& epsX);

    double delta(const double r);

    void getf(const int i, const int n, const Vector<double>& X, Vector<double>& f, const double fext, 
        Vector<double>& RPAR);
    void getj(const int i, const int n, const Vector<double>& X, Array<double>& JAC);
              
    void init(const int nX, Vector<double> &X);

    void integ_cn2(const int imyo, const int nX, Vector<double>& X, const double Ts, const double Ti,
        const double Istim, const double Ksac, Vector<int>& IPAR, Vector<double>& RPAR);

    void integ_fe(const int imyo, const int nX, Vector<double>& X, const double Ts, const double Ti, 
        const double Istim, const double Ksac, Vector<double>& RPAR);

    void integ_rk(const int imyo, const int nX, Vector<double>& X, const double Ts, const double Ti,
        const double Istim, const double Ksac, Vector<double>& RPAR);

    double step(const double r);

};

#endif

