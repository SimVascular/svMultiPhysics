// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef CEP_MOD_TTP_H 
#define CEP_MOD_TTP_H 

#include "Array.h"
#include "Vector.h"
#include <array>

#include <optional>
#include <functional>

template <class T>
T& make_ref(T&& x) { return x; }

/// @brief This module defines data structures for ten Tusscher-Panfilov
/// epicardial cellular activation model for cardiac electrophysiology
///
/// The classes defined here duplicate the data structures in the Fortran TPPMOD module defined 
/// in CEPMOD_TTP.f and PARAMS_TPP.f files. 
class CepModTtp
{
  public:
    CepModTtp();
    ~CepModTtp();


//--------------------------------------------------------------------
//
//     Constants for TenTusscher-Panfilov Ventricular Myocyte Model.
//
//--------------------------------------------------------------------

//     Default model parameters
      /// Gas constant [J/mol/K]
      double Rc = 8314.472;

      /// Temperature [K]
      double Tc = 310.0;

      /// Faraday constant [C/mmol]
      double Fc = 96485.3415;

      /// Cell capacitance per unit surface area [uF/cm^{2}]
      double Cm = 0.185;

      /// Surface to volume ratio [um^{-1}]
      double sV = 0.2;

      /// Cellular resistivity [\f$\Omega\f$-cm]
      double rho = 162.0;

      /// Cytoplasmic volume [um^{3}]
      double V_c = 16.404E-3;

      /// Sacroplasmic reticulum volume [um^{3}]
      double V_sr = 1.094E-3;

      /// Subspace volume [um^{3}]
      double V_ss = 5.468E-5;

      /// Extracellular K concentration [mM]
      double K_o = 5.4;

      /// Extracellular Na concentration [mM]
      double Na_o = 140.0;

      /// Extracellular Ca concentration [mM]
      double Ca_o = 2.0;

      /// Maximal I_Na conductance [nS/pF]
      double G_Na = 14.838;

      /// Maximal I_K1 conductance [nS/pF]
      double G_K1 = 5.405;

      /// Maximal epicardial I_to conductance [nS/pF]
      Vector<double> G_to = {0.294, 0.073, 0.294};

      /// Maximal I_Kr conductance [nS/pF]
      double G_Kr = 0.153;

//     G_Kr for spiral wave breakup
//      double G_Kr = 0.172;     // units: nS/pF

      /// Maximal epicardial I_Ks conductance [nS/pF]
      Vector<double> G_Ks = {0.392, 0.392, 0.098};

//     G_Ks for spiral wave breakup (epi)
//      double G_Ks(3) = (/0.441, 0.392_RKIND, 0.098_RKIND/)

      /// Relative I_Ks permeability to Na [-]
      double p_KNa = 3.E-2;

      /// Maximal I_CaL conductance [cm^{3}/uF/ms]
      double G_CaL = 3.98E-5;

      /// Maximal I_NaCa [pA/pF]
      double K_NaCa = 1000.;

      /// Voltage dependent parameter of I_NaCa [-]
      double gamma = 0.35;

      /// Ca_i half-saturation constant for I_NaCa [mM]
      double K_mCa = 1.38;

      /// Na_i half-saturation constant for I_NaCa [mM]
      double K_mNai = 87.5;

      /// Saturation factor for I_NaCa [-]
      double K_sat = 0.1;

      /// Factor enhancing outward nature of I_NaCa [-]
      double alpha = 2.5;

      /// Maximal I_NaK [pA/pF]
      double p_NaK = 2.724;

      /// K_o half-saturation constant of I_NaK [mM]
      double K_mK = 1.;

      /// Na_i half-saturation constant of I_NaK [mM]
      double K_mNa = 40.;

      /// Maximal I_pK conductance [nS/pF]
      double G_pK = 1.46E-2;

//     G_pK for spiral wave breakup
//      double G_pK = 2.19E-3;    // units: nS/pF

      /// Maximal I_pCa conductance [pA/pF]
      double G_pCa = 0.1238;

//     G_pCa for spiral wave breakup
//      double G_pCa = 0.8666;    // units: pA/pF

      /// Half-saturation constant of I_pCa [mM]
      double K_pCa = 5.E-4;

      /// Maximal I_bNa conductance [nS/pF]
      double G_bNa = 2.9E-4;

      /// Maximal I_bCa conductance [nS/pF]
      double G_bCa = 5.92E-4;

      /// Maximal I_up conductance [mM/ms]
      double Vmax_up = 6.375E-3;

      /// Half-saturation constant of I_up [mM]
      double K_up = 2.5E-4;

      /// Maximal I_rel conductance [mM/ms]
      double V_rel = 0.102;

      /// R to O and RI to I, I_rel transition rate [mM^{-2}/ms]
      double k1p = 0.15;

      /// O to I and R to RI, I_rel transition rate [mM^{-1}/ms]
      double k2p = 4.5E-2;

      /// O to R and I to RI, I_rel transition rate [ms^{-1}]
      double k3 = 6.E-2;

      /// I to O and Ri to I, I_rel transition rate [ms^{-1}]
      double k4 = 5.E-3;

      /// Ca_sr half-saturation constant of k_casr [mM]
      double EC = 1.5;

      /// Maximum value of k_casr [-]
      double max_sr = 2.5;

      /// Minimum value of k_casr [-]
      double min_sr = 1.;

      /// Maximal I_leak conductance [mM/ms]
      double V_leak = 3.6E-4;

      /// Maximal I_xfer conductance [mM/ms]
      double V_xfer = 3.8E-3;

      /// Total cytoplasmic buffer concentration [mM]
      double Buf_c = 0.2;

      /// Ca_i half-saturation constant for cytplasmic buffer [mM]
      double K_bufc = 1.E-3;

      /// Total sacroplasmic buffer concentration [mM]
      double Buf_sr = 10.;

      /// Ca_sr half-saturation constant for subspace buffer [mM]
      double K_bufsr = 0.3;

      /// Total subspace buffer concentration [mM]
      double Buf_ss = 0.4;

      /// Ca_ss half-saturation constant for subspace buffer [mM]
      double K_bufss = 2.5E-4;

      /// Resting potential [mV]
      double Vrest = -85.23;

//     Electromechanics coupling parameters: active stress model
      /// Resting Ca concentration [mM]
      double Ca_rest = 5.E-5;

      /// Critical Ca concentration [mM]
      double Ca_crit = 8.E-4;

      /// Saturation of concentration [MPa/mM]
      double eta_T = 12.5;

      /// Minimum activation [ms^{-1}]
      double eps_0 = 0.1;

      /// Maximum activation [ms^{-1}]
      double eps_i = 1. ;

      /// Transition rate [mM^{-1}]
      double xi_T = 4.E3;

//     Electromechanics coupling parameters: active strain model
//

      /// Active force of sacromere [-mM^{-2}]
      double alFa = -4.E6;

      /// Resting Ca concentration [mM]
      double c_Ca0 = 2.155E-4;

      /// Viscous-type constant [ms-mM^{-2}]
      double mu_Ca = 5.E6;

//     Force-length relationship parameters
      /// Initial length of sacromeres [um]
      double SL0 = 1.95;

      /// Min. length of sacromeres [um]
      double SLmin = 1.7;

      /// Max. length of sacromeres [um]
      double SLmax = 2.6;

      /// Fourier coefficients
      double f0  = -4333.618335582119;
      double fc1 =  2570.395355352195;
      double fs1 = -2051.827278991976;
      double fc2 =  1329.53611689133;
      double fs2 =  302.216784558222;
      double fc3 =  104.943770305116;
      double fs3 =  218.375174229422;

//-----------------------------------------------------------------------
//     Scaling factors
      /// Voltage scaling
      double Vscale  = 1.;

      /// Time scaling
      double Tscale  = 1.;

      /// Voltage offset parameter
      double Voffset = 0.;

//-----------------------------------------------------------------------
//     Variables
      /// Reverse potentials for Na, K, Ca
      double E_Na;
      double E_K;
      double E_Ca;
      double E_Ks;
//     Cellular transmembrane currents
      /// Fast sodium current
      double I_Na;

      /// inward rectifier outward current
      double I_K1;

      /// transient outward current
      double I_to;

      /// rapid delayed rectifier current
      double I_Kr;

      /// slow delayed rectifier current
      double I_Ks;

      /// L-type Ca current
      double I_CaL;

      /// Na-Ca exchanger current
      double I_NaCa;

      /// Na-K pump current
      double I_NaK;

      /// plateau Ca current
      double I_pCa;

      /// plateau K current
      double I_pK;

      /// background Ca current
      double I_bCa;

      /// background Na current
      double I_bNa;

      /// sacroplasmic reticulum Ca leak current
      double I_leak;

      /// sacroplasmic reticulum Ca pump current
      double I_up;

      /// Ca induced Ca release current
      double I_rel;

      /// diffusive Ca current
      double I_xfer;
//-----------------------------------------------------------------------
//     State variables
      double V;
      double K_i;
      double Na_i;
      double Ca_i;
      double Ca_ss;
      double Ca_sr;
      double R_bar;

//     Gating variables (runtime, steady state)
      double xr1, xr1i;
      double xr2, xr2i;
      double xs, xsi;
      double m, mi;
      double h, hi;
      double j, ji;
      double d, di;
      double f, fi;
      double f2, f2i;
      double fcass, fcassi;
      double s, si;
      double r, ri;

//     Other variables
      double k1;
      double k2;
      double k_casr;
      double O;

//     Jacobian variables
      double E_Na_Nai, E_K_Ki, E_Ca_Cai, E_Ks_Ki, E_Ks_Nai;
      double I_Na_V, I_Na_Nai;
      double I_to_V, I_to_Ki;
      double I_K1_V, I_K1_Ki;
      double I_Kr_V, I_Kr_Ki;
      double I_Ks_V, I_Ks_Ki, I_Ks_Nai;
      double I_CaL_V, I_CaL_Cass;
      double I_NaCa_V, I_NaCa_Nai, I_NaCa_Cai;
      double I_NaK_V, I_NaK_Nai;
      double I_pCa_Cai;
      double I_pK_V, I_pK_Ki;
      double I_bCa_V, I_bCa_Cai;
      double I_bNa_V, I_bNa_Nai;
      double I_leak_Cai, I_leak_Casr;
      double I_up_Cai;
      double I_rel_Cass, I_rel_Casr, I_rel_Rbar;
      double I_xfer_Cai, I_xfer_Cass;
      double k_casr_sr, k1_casr, O_Casr, O_Cass, O_Rbar;

    void actv_strn(const double c_Ca, const double I4f, const double dt, double& gf);
    void actv_strs(const double c_Ca, const double dt, double& Tact, double& epsX);

    void getf(const int i, const int nX, const int nG, const Vector<double>& X, const Vector<double>& Xg, 
        Vector<double>& dX, const double I_stim, const double K_sac, Vector<double>& RPAR);

    void getj(const int i, const int nX, const int nG, const Vector<double>& X, const Vector<double>& Xg, 
        Array<double>& JAC, const double Ksac);

    void init(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg);

    void init(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg,
        Vector<double>& X0, Vector<double>& Xg0);

    void integ_cn2(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg,
        const double Ts, const double dt, const double Istim, const double Ksac, 
        Vector<int>& IPAR, Vector<double>& RPAR);

    void integ_fe(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg, 
        const double Ts, const double dt, const double Istim, const double Ksac, Vector<double>& RPAR);

    void integ_rk(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg, 
        const double Ts, const double dt, const double Istim, const double Ksac, Vector<double>& RPAR);

    void update_g(const int i, const double dt, const int n, const int nG, const Vector<double>& X, 
        Vector<double>& Xg);

};

#endif

