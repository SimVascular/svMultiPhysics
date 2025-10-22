// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "CepModTtp.h"

#include "mat_fun.h"
#include <math.h>

CepModTtp::CepModTtp()
{
}

CepModTtp::~CepModTtp()
{
}

/// @brief Compute macroscopic fiber strain based on sacromere force-length relationship and calcium concentration
void CepModTtp::actv_strn(const double c_Ca, const double I4f, const double dt, double& gf)
{
  // fiber length
  double SL = I4f * SL0;

  //  Sacromere force-length relationship
  if (SL >= SLmin && SL <= SLmax) {
    SL = 0.5*f0 + fc1*cos(SL) + fs1*sin(SL) + fc2*cos(2.0*SL) + fs2*sin(2.0*SL)  + fc3*cos(3.0*SL) + fs3*sin(3.0*SL);
  } else { 
    SL = 0.0;
  } 

  // Active force
  double Fa = alFa * (c_Ca-c_Ca0)*(c_Ca-c_Ca0) * SL;
  double rtmp = 2.0*I4f*(1.0/ pow(1.0+gf,3.0) - 1.0);
  gf = gf + dt*(Fa + rtmp)/(mu_Ca * c_Ca * c_Ca);
}

void CepModTtp::actv_strs(const double c_Ca, const double dt, double& Tact, double& epsX)
{
  epsX = exp(-exp(-xi_T*(c_Ca - Ca_crit)));
  epsX = eps_0 + (eps_i - eps_0)*epsX;
  double nr  = Tact + epsX*dt*eta_T*(c_Ca - Ca_rest);
  Tact = nr / (1.0 + epsX*dt);
}

/// @brief Compute currents and time derivatives of state variables
///
/// Note that is 'i' the myocardium zone id: 1, 2 or 3.
///
/// Reproduces Fortran 'GETF()'.
void CepModTtp::getf(const int i, const int nX, const int nG, const Vector<double>& X, const Vector<double>& Xg, Vector<double>& dX, 
    const double I_stim, const double K_sac, Vector<double>& RPAR)
{
  // Local copies of state variables
  double V     = X(0);
  double K_i   = X(1);
  double Na_i  = X(2);
  double Ca_i  = X(3);
  double Ca_ss = X(4);
  double Ca_sr = X(5);
  double R_bar = X(6);

  // Local copies of gating variables
  double xr1   = Xg(0);
  double xr2   = Xg(1);
  double xs    = Xg(2);
  double m     = Xg(3);
  double h     = Xg(4);
  double j     = Xg(5);
  double d     = Xg(6);
  double f     = Xg(7);
  double f2    = Xg(8);
  double fcass = Xg(9);
  double s     = Xg(10);
  double r     = Xg(11);

  // Stretch-activated currents
  double I_sac = K_sac * (Vrest - V);

  // Diff = 1. / (1.0D1 * rho * Cm * sV)
  double RT   = Rc * Tc / Fc;
  double E_K  = RT * log(K_o/K_i);
  double E_Na = RT * log(Na_o/Na_i);
  double E_Ca = 0.5 * RT * log(Ca_o/Ca_i);
  double E_Ks = RT * log( (K_o + p_KNa*Na_o)/(K_i + p_KNa*Na_i) );

  // I_Na: Fast sodium current
  double I_Na = G_Na * pow(m,3.0) * h * j * (V - E_Na);

  // I_to: transient outward current
  double I_to = G_to[i-1] * r * s * (V - E_K);
  //I_to = G_to(i) * r * s * (V - E_K)

  // I_K1: inward rectifier outward current
  double e1   = exp(0.06*(V - E_K - 200.0));
  double e2   = exp(2.E-4*(V - E_K + 100.0));
  double e3   = exp(0.1*(V - E_K - 10.0));
  double e4   = exp(-0.5*(V - E_K));
  double a    = 0.1/(1.0 + e1);
  double b    = (3.0*e2 + e3) / (1.0 + e4);
  double tau  = a / (a + b);
  double sq5  = sqrt(K_o/5.4);
  double I_K1 = G_K1 * sq5 * tau * (V - E_K);

  // I_Kr: rapid delayed rectifier current
  double I_Kr = G_Kr * sq5 * xr1 * xr2 * (V - E_K);

  // I_Ks: slow delayed rectifier current
  double I_Ks = G_Ks[i-1] * pow(xs,2.0) * (V - E_Ks);

  // I_CaL: L-type Ca current
  a = 2.0*(V-15.)/RT;
  b = 2.0*a*Fc * (0.25*Ca_ss*exp(a) - Ca_o) / (exp(a)-1.0);
  double I_CaL = G_CaL * d * f * f2 * fcass * b;

  // I_NaCa: Na-Ca exchanger current
  e1 = exp(gamma*V/RT);
  e2 = exp((gamma-1.)*V/RT);
  double n1 = e1*pow(Na_i,3.0)*Ca_o - e2*pow(Na_o,3.0)*Ca_i*alpha;
  double d1 = pow(K_mNai,3.0) + pow(Na_o,3.0);
  double d2 = K_mCa + Ca_o;
  double d3 = 1.0 + K_sat*e2;
  I_NaCa = K_NaCa * n1 / (d1*d2*d3);

  // I_NaK: Na-K pump current
  e1 = exp(-0.1*V/RT);
  e2 = exp(-V/RT);
  n1 = p_NaK * K_o * Na_i;
  d1 = K_o + K_mK;
  d2 = Na_i + K_mNa;
  d3 = 1.0 + 0.1245*e1 + 0.0353*e2;
  double I_NaK = n1 / (d1*d2*d3);

  // I_pCa: plateau Ca current
  double I_pCa = G_pCa * Ca_i / (K_pCa + Ca_i);

  // I_pK: plateau K current
  double I_pK  = G_pK * (V-E_K) / (1.0 + exp((25.0-V)/5.98));

  // I_bCa: background Ca current
  double I_bCa = G_bCa * (V - E_Ca);

  // I_bNa: background Na current
  double I_bNa = G_bNa * (V - E_Na);

  // I_leak: Sacroplasmic Reticulum Ca leak current
  double I_leak = V_leak * (Ca_sr - Ca_i);

  // I_up: Sacroplasmic Reticulum Ca pump current
  double I_up  = Vmax_up / (1.0 + pow(K_up/Ca_i,2.0));

  // I_rel: Ca induced Ca current (CICR)
  k_casr = max_sr - ((max_sr-min_sr) / (1.0 + pow(EC/Ca_sr,2.0)));
  k1 = k1p / k_casr;
  O = k1 * R_bar * pow(Ca_ss,2.0) / (k3 + k1*pow(Ca_ss,2.0));
  I_rel  = V_rel * O * (Ca_sr - Ca_ss);

  //  I_xfer: diffusive Ca current between Ca subspae and cytoplasm
  I_xfer = V_xfer * (Ca_ss - Ca_i);

  // Now compute time derivatives
  // dV/dt: rate of change of transmembrane voltage
  //
  dX(0)  = -(I_Na + I_to + I_K1 + I_Kr + I_Ks + I_CaL + I_NaCa + 
      I_NaK + I_pCa + I_pK + I_bCa + I_bNa  + I_stim) + I_sac;

  // dK_i/dt
  dX(1)  = -(Cm/(V_c*Fc)) * (I_K1 + I_to + I_Kr + I_Ks + I_pK - 2.0*I_NaK + I_stim);

  //  dNa_i/dt
  dX(2)  = -(Cm/(V_c*Fc)) * (I_Na + I_bNa + 3.0*(I_NaK + I_NaCa));

  // dCa_i/dt
  n1 = (I_leak - I_up)*V_sr/V_c + I_xfer;
  double n2 = -(Cm/(V_c*Fc)) * (I_bCa + I_pCa - 2.0*I_NaCa) / 2.0;
  d1 = 1.0 + K_bufc*Buf_c/ pow(Ca_i + K_bufc,2.0);
  dX(3)  = (n1 + n2)/d1;

  // dCa_ss: rate of change of Ca_ss
  n1 = (-I_CaL*Cm/(2.0*Fc) + I_rel*V_sr - V_c*I_xfer)/V_ss;
  d1 = 1.0 + K_bufss*Buf_ss/ pow(Ca_ss + K_bufss,2.0);
  dX(4)  = n1 / d1;

  // dCa_sr: rate of change of Ca_sr
  n1 = I_up - I_leak - I_rel;
  d1 = 1. + K_bufsr*Buf_sr/ pow(Ca_sr + K_bufsr,2.0);
  dX(5)  = n1 / d1;

  // Rbar: ryanodine receptor
  k2 = k2p * k_casr;
  dX(6)  = -k2*Ca_ss*R_bar + k4*(1.0 - R_bar);

  // Quantities to be written to file
  RPAR(2)  = I_Na;
  RPAR(3)  = I_K1;
  RPAR(4)  = I_to;
  RPAR(5)  = I_Kr;
  RPAR(6)  = I_Ks;
  RPAR(7)  = I_CaL;
  RPAR(8)  = I_NaCa;
  RPAR(9) = I_NaK;
  RPAR(10) = I_pCa;
  RPAR(11) = I_pK;
  RPAR(12) = I_bCa;
  RPAR(13) = I_bNa;
  RPAR(14) = I_leak;
  RPAR(15) = I_up;
  RPAR(16) = I_rel;
  RPAR(17) = I_xfer;
}

void CepModTtp::getj(const int i, const int nX, const int nG, const Vector<double>& X, const Vector<double>& Xg, 
    Array<double>& JAC, const double Ksac)
{

  double RT, a, b, c, tau, sq5, e1, e2, e3, e4, n1, n2, d1, d2, d3;

  // Local copies of state variables
  double V     = X(0);
  double K_i   = X(1);
  double Na_i  = X(2);
  double Ca_i  = X(3);
  double Ca_ss = X(4);
  double Ca_sr = X(5);
  double R_bar = X(6);

  // Local copies of gating variables
  double xr1   = Xg(0);
  double xr2   = Xg(1);
  double xs    = Xg(2);
  double m     = Xg(3);
  double h     = Xg(4);
  double j     = Xg(5);
  double d     = Xg(6);
  double f     = Xg(7);
  double f2    = Xg(8);
  double fcass = Xg(9);
  double s     = Xg(10);
  double r     = Xg(11);

  RT = Rc * Tc / Fc;
  double E_K  = RT * log(K_o/K_i);
  double E_Na = RT * log(Na_o/Na_i);
  double E_Ca = 0.5 * RT * log(Ca_o/Ca_i);
  double E_Ks = RT * log( (K_o + p_KNa*Na_o)/(K_i + p_KNa*Na_i) );

  E_K_Ki   = -RT / K_i;
  E_Na_Nai = -RT / Na_i;
  E_Ca_Cai = -RT / Ca_i / 2.0;
  E_Ks_Ki  = -RT / (K_i + p_KNa*Na_i);
  E_Ks_Nai = p_KNa * E_Ks_Ki;

  // I_Na: Fast sodium current
  I_Na = G_Na * pow(m,3.0) * h * j * (V - E_Na);
  I_Na_V   = G_Na * pow(m,3.0) * h * j;
  I_Na_Nai = I_Na_V * (-E_Na_Nai);

  // I_to: transient outward current
  I_to = G_to[i-1] * r * s * (V - E_K);
  I_to_V  = G_to[i-1] * r * s;
  I_to_Ki = I_to_V * (-E_K_Ki);

  // I_K1: inward rectifier outward current
  e1   = exp(0.060*(V - E_K - 200.0));
  e2   = exp(2.E-40*(V - E_K + 100.0));
  e3   = exp(0.10*(V - E_K - 10.0));
  e4   = exp(-0.50*(V - E_K));
  a    = 0.10/(1.0 + e1);
  b    = (3.0*e2 + e3) / (1.0 + e4);
  tau  = a / (a + b);
  sq5  = sqrt(K_o/5.40);
  n1   = -6.E-30*e1 / pow(1.0 + e1,2.0);
  n2   = (6.E-40*e2 + 0.10*e3 + 0.50*b*e4) / (1.0 + e4);
  n1   = (a + b)*n1 - a*(n1 + n2);
  d1   = pow(a + b,2.0);
  I_K1 = G_K1 * sq5 * tau * (V - E_K);
  I_K1_V  = G_K1 * sq5 * (tau + (V - E_K)*n1/d1);
  I_K1_Ki = I_K1_V * (-E_K_Ki);

  // I_Kr: rapid delayed rectifier current
  I_Kr = G_Kr * sq5 * xr1 * xr2 * (V - E_K);
  I_Kr_V   = G_Kr * sq5 * xr1 * xr2;
  I_Kr_Ki  = I_Kr_V * (-E_K_Ki);

  // I_Ks: slow delayed rectifier current
  I_Ks = G_Ks[i-1] * pow(xs,2.0) * (V - E_Ks);
  I_Ks_V   = G_Ks[i-1] * pow(xs,2.0);
  I_Ks_Ki  = I_Ks_V * (-E_Ks_Ki);
  I_Ks_Nai = I_Ks_V * (-E_Ks_Nai);

  //  I_CaL: L-type Ca current
  a = 2.0*(V-15.0)/RT;
  b = (0.250*Ca_ss*exp(a) - Ca_o) / (exp(a)-1.0);
  c = G_CaL * d * f * f2 * fcass * (2.0*a*Fc);
  n1 = (exp(a)/RT) / (exp(a)-1.0);
  I_CaL = c * b;
  I_CaL_V = I_CaL/(V-15.0) + n1*(c*0.50*Ca_ss - 2.0*I_CaL);
  I_CaL_Cass  = c * 0.250 * n1 * RT;

  //  I_NaCa: Na-Ca exchanger current
  e1 = exp(gamma*V/RT);
  e2 = exp((gamma-1.0)*V/RT);
  n1 = e1*pow(Na_i,3.0)*Ca_o - e2*pow(Na_o,3.0)*Ca_i*alpha;
  d1 = pow(K_mNai,3.0) + pow(Na_o,3.0);
  d2 = K_mCa + Ca_o;
  d3 = 1.0 + K_sat*e2;
  c  = 1.0/(d1*d2*d3);
  I_NaCa = K_NaCa * n1 * c;

  n1 = K_NaCa * c * ( e1*pow(Na_i,3.0)*Ca_o*(gamma/RT) - e2*pow(Na_o,3.0)*Ca_i*alpha*((gamma-1.0)/RT) );
  n2 = I_NaCa*K_sat*((gamma-1.0)/RT)*e2/d3;
  I_NaCa_V   = n1 - n2;
  I_NaCa_Nai =  K_NaCa * e1 * (3.0*pow(Na_i,2.0)) * Ca_o * c;
  I_NaCa_Cai = -K_NaCa * e2 * pow(Na_o,3.0) * alpha * c;

  // I_NaK: Na-K pump current
  e1 = exp(-0.10*V/RT);
  e2 = exp(-V/RT);
  n1 = p_NaK * K_o * Na_i;
  d1 = K_o + K_mK;
  d2 = Na_i + K_mNa;
  d3 = 1.0 + 0.12450*e1 + 0.03530*e2;
  I_NaK = n1 / (d1*d2*d3);
  n1 = (0.012450*e1 + 0.03530*e2)/RT;
  I_NaK_V = I_NaK * n1 / d3;
  I_NaK_Nai = I_NaK * K_mNa/(Na_i*d2);

  // I_pCa: plateau Ca current
  d1 = (K_pCa + Ca_i);
  I_pCa = G_pCa * Ca_i / d1;
  I_pCa_Cai = G_pCa * K_pCa/(d1*d1);

  //  I_pK: plateau K current
  e1 = exp((25.0-V)/5.980);
  I_pK  = G_pK * (V-E_K) / (1.0 + e1);
  I_pK_V  = (G_pK + I_pK*e1/5.980) / (1.0+e1);
  I_pK_Ki = G_pK * (-E_K_Ki) / (1.0 + e1);

  // I_bCa: background Ca current
  I_bCa = G_bCa * (V - E_Ca);
  I_bCa_V = G_bCa;
  I_bCa_Cai = G_bCa * (-E_Ca_Cai);

  //  I_bNa: background Na current
  I_bNa = G_bNa * (X(0) - E_Na);
  I_bNa_V = G_bNa;
  I_bNa_Nai = G_bNa * (-E_Na_Nai);

  // I_leak: Sacroplasmic Reticulum Ca leak current
  I_leak = V_leak * (Ca_sr - Ca_i);
  I_leak_Cai  = -V_leak;
  I_leak_Casr =  V_leak;

  // I_up: Sacroplasmic Reticulum Ca pump current
  d1 = 1.0 + pow(K_up/Ca_i,2.0);
  I_up  = Vmax_up / d1;
  I_up_Cai = (I_up / d1) * (2.0*pow(K_up,2.0) / pow(Ca_i,3.0));

  // I_rel: Ca induced Ca current (CICR)
  n1 = max_sr - min_sr;
  d1 = 1.0 + pow(EC/Ca_sr,2.0);
  k_casr = max_sr - (n1/d1);
  k1 = k1p / k_casr;
  n2 = Ca_ss*2.0;
  d2 = k3 + k1*n2;
  O = k1 * R_bar * n2 / d2;
  I_rel  = V_rel * O * (Ca_sr - Ca_ss);

  k_casr_sr = (n1 / pow(d1,2.0) ) * (2.0*pow(EC,2.0) / pow(Ca_sr,3.0));
  //k_casr_sr = (n1 / (d1**2.0)) * (2.0*EC**2.0 / Ca_sr**3.0);
  k1_casr   = -k1p * k_casr_sr / pow(k_casr,2.0);
  O_Cass = 2.0 * k3 * O / (Ca_ss * d2);
  O_Casr = k1_casr * n2 * (R_bar - O) / d2;
  O_Rbar = k1 * n2 / d2;

  I_rel_Cass = V_rel * (O_Cass*(Ca_sr - Ca_ss) - O);
  I_rel_Casr = V_rel * (O_Casr*(Ca_sr - Ca_ss) + O);
  I_rel_Rbar = V_rel * O_Rbar *(Ca_sr - Ca_ss);

  //  I_xfer: diffusive Ca current between Ca subspae and cytoplasm
  I_xfer = V_xfer * (Ca_ss - Ca_i);
  I_xfer_Cai  = -V_xfer;
  I_xfer_Cass =  V_xfer;

  // Compute Jacobian matrix
  //
  JAC = 0.0;
  c = Cm/(V_c*Fc);

  //  V
  JAC(0,0)  = -(I_Na_V + I_to_V + I_K1_V + I_Kr_V + I_Ks_V + I_CaL_V + I_NaCa_V + I_NaK_V + I_pK_V + I_bCa_V + I_bNa_V + Ksac);
  JAC(0,1)  = -(I_to_Ki + I_K1_Ki + I_Kr_Ki + I_Ks_Ki + I_pK_Ki);
  JAC(0,2)  = -(I_Na_Nai + I_Ks_Nai + I_NaCa_Nai + I_NaK_Nai + I_bNa_Nai);
  JAC(0,3)  = -(I_NaCa_Cai + I_pCa_Cai + I_bCa_Cai);
  JAC(0,4) = -I_CaL_Cass;

  // K_i
  JAC(1,0)  = -c * (I_K1_V + I_to_V + I_Kr_V + I_Ks_V + I_pK_V - 2.0*I_NaK_V );
  JAC(1,1)  = -c * (I_K1_Ki + I_to_Ki + I_Kr_Ki + I_Ks_Ki + I_pK_Ki);
  JAC(1,2)  = -c * (I_Ks_Nai - 2.0*I_NaK_Nai);

  // Na_i
  JAC(2,0)  = -c * (I_Na_V + I_bNa_V + 3.0*(I_NaK_V + I_NaCa_V));
  JAC(2,2)  = -c * (I_Na_Nai + I_bNa_Nai + 3.0*(I_NaK_Nai + I_NaCa_Nai));
  JAC(2,3)  = -c * (3.0*I_NaCa_Cai);

  //     Ca_i
  n1 = (I_leak - I_up)*V_sr/V_c + I_xfer - 0.50*c*(I_bCa + I_pCa - 2.0*I_NaCa);
  n2 = (I_leak_Cai - I_up_Cai)*V_sr/V_c + I_xfer_Cai - 0.50*c*(I_bCa_Cai + I_pCa_Cai - 2.0*I_NaCa_Cai);
  d1 = 1.0 + K_bufc*Buf_c / pow(Ca_i + K_bufc,2.0);
  d2 = 2.0*K_bufc*Buf_c / pow(Ca_i + K_bufc,3.0);
  JAC(3,0)  = -c * (I_bCa_V - 2.0*I_NaCa_V) / 2.0 / d1;
  JAC(3,2)  = c * I_NaCa_Nai / d1;
  JAC(3,3)  = (n2 + n1*d2/d1) / d1;
  JAC(3,4) = I_xfer_Cass / d1;
  JAC(3,5) = (I_leak_Casr*V_sr/V_c) / d1;

  //     Ca_ss
  a  = Cm/(2.0*Fc*V_ss);
  b  = V_sr / V_ss;
  c  = V_c / V_ss;
  n1 = -a*I_CaL + b*I_rel - c*I_xfer;
  n2 = -a*I_CaL_Cass + b*I_rel_Cass - c*I_xfer_Cass;
  d1 = 1.0 + K_bufss*Buf_ss/ pow(Ca_ss + K_bufss,2.0);
  d2 = 2.0*K_bufss*Buf_ss / pow(Ca_ss + K_bufss,3.0);
  JAC(4,0)  = -a * I_CaL_V / d1;
  JAC(4,3)  = -c * I_xfer_Cai / d1;
  JAC(4,4) = (n2 + n1*d2/d1) / d1;
  JAC(4,5) = b * I_rel_Casr / d1;
  JAC(4,6) = b * I_rel_Rbar / d1;

  // Ca_sr
  n1 = I_up - I_leak - I_rel;
  n2 = -(I_leak_Casr + I_rel_Casr);
  d1 = 1.0 + K_bufsr*Buf_sr / pow(Ca_sr + K_bufsr,2.0);
  d2 = 2.0*K_bufsr*Buf_sr / pow(Ca_sr + K_bufsr,3.0);
  JAC(5,3) = (I_up_Cai - I_leak_Cai) / d1;
  JAC(5,4) = -I_rel_Cass / d1;
  JAC(5,5) = (n2 + n1*d2/d1) / d1;
  JAC(5,6) = -I_rel_Rbar / d1;

  // Rbar: ryanodine receptor
  k2 = k2p * k_casr;
  JAC(6,4) = -k2 * R_bar;
  JAC(6,5) = -(k2p * k_casr_sr) * Ca_ss * R_bar;
  JAC(6,6) = -(k2*Ca_ss + k4);
}

void CepModTtp::init(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg )
{
  switch (imyo) {

    // epi
    case 1:
      // Initialize state variables
      X(0)   = -85.23;      // V      (units: mV)
      X(1)   =  136.89;     // K_i    (units: mM)
      X(2)   =  8.6040;     // Na_i   (units: mM)
      X(3)   =  1.26E-4;    // Ca_i   (units: mM)
      X(4)   =  3.6E-4;     // Ca_ss  (units: mM)
      X(5)   =  3.64;       // Ca_sr  (units: mM)
      X(6)   =  0.9073;     // R'     (dimensionless)

      // Initialize gating variables
      Xg(0)  =  6.21E-3;    // x_r1   (dimensionless)
      Xg(1)  =  0.4712;     // x_r2   (dimensionless)
      Xg(2)  =  9.5E-3;     // x_s    (dimensionless)
      Xg(3)  =  1.72E-3;    // m      (dimensionless)
      Xg(4)  =  0.7444;     // h      (dimensionless)
      Xg(5)  =  0.7045;     // j      (dimensionless)
      Xg(6)  =  3.373E-5;   // d      (dimensionless)
      Xg(7)  =  0.7888;     // f      (dimensionless)
      Xg(8)  =  0.9755;     // f_2    (dimensionless)
      Xg(9)  =  0.9953;     // f_cass (dimensionless)
      Xg(10) =  0.999998;   // s      (dimensionless)
      Xg(11) =  2.42E-8;    // r      (dimensionless)
    break;

    // endo
    case 2:
      // Initialize state variables
      X(0)   = -86.709;     // V      (units: mV)
      X(1)   =  138.4;      // K_i    (units: mM)
      X(2)   =  10.355;     // Na_i   (units: mM)
      X(3)   =  1.3E-4;    // Ca_i   (units: mM)
      X(4)   =  3.6E-4;     // Ca_ss  (units: mM)
      X(5)   =  3.715;      // Ca_sr  (units: mM)
      X(6)   =  0.9068;     // R'     (dimensionless)

      // Initialize gating variables
      Xg(0)  =  4.48E-3;    // x_r1   (dimensionless)
      Xg(1)  =  0.476;      // x_r2   (dimensionless)
      Xg(2)  =  8.7E-3;     // x_s    (dimensionless)
      Xg(3)  =  1.55E-3;    // m      (dimensionless)
      Xg(4)  =  0.7573;     // h      (dimensionless)
      Xg(5)  =  0.7225;     // j      (dimensionless)
      Xg(6)  =  3.164E-5;   // d      (dimensionless)
      Xg(7)  =  0.8009;     // f      (dimensionless)
      Xg(8)  =  0.9778;     // f_2    (dimensionless)
      Xg(9)  =  0.9953;     // f_cass (dimensionless)
      Xg(10) =  0.3212;     // s      (dimensionless)
      Xg(11) =  2.235E-8;   // r      (dimensionless)
    break;

    // mid-myo
    case 3:
      // Initialize state variables
      X(0)   = -85.423;     // V      (units: mV)
      X(1)   =  138.52;     // K_i    (units: mM)
      X(2)   =  10.132;     // Na_i   (units: mM)
      X(3)   =  1.53E-4;    // Ca_i   (units: mM)
      X(4)   =  4.2E-4;     // Ca_ss  (units: mM)
      X(5)   =  4.272;      // Ca_sr  (units: mM)
      X(6)   =  0.8978;     // R'     (dimensionless)

      //     Initialize gating variables
      Xg(0)  =  1.65E-2;    // x_r1   (dimensionless)
      Xg(1)  =  0.4730;     // x_r2   (dimensionless)
      Xg(2)  =  1.74E-2;    // x_s    (dimensionless)
      Xg(3)  =  1.65E-3;    // m      (dimensionless)
      Xg(4)  =  0.7490;     // h      (dimensionless)
      Xg(5)  =  0.6788;     // j      (dimensionless)
      Xg(6)  =  3.288E-5;   // d      (dimensionless)
      Xg(7)  =  0.7026;     // f      (dimensionless)
      Xg(8)  =  0.9526;     // f_2    (dimensionless)
      Xg(9)  =  0.9942;     // f_cass (dimensionless)
      Xg(10) =  0.999998;   // s      (dimensionless)
      Xg(11) =  2.347E-8;   // r      (dimensionless)
    break;
  }

}

void CepModTtp::init(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg, 
    Vector<double>& X0, Vector<double>& Xg0)
{
  init(imyo, nX, nG, X, Xg);

  if (X0.size() != 0) {
    X = X0;
  }

  if (Xg0.size() != 0) {
    Xg = Xg0;
  }
}

/// @brief Time integration performed using Crank-Nicholson method
void CepModTtp::integ_cn2(const int imyo, const int nX, const int nG, Vector<double>& Xn, Vector<double>& Xg, 
    const double Ts, const double dt, const double Istim, const double Ksac, Vector<int>& IPAR, Vector<double>& RPAR)
{
  int itMax = IPAR(0);
  double atol  = RPAR(0);
  double rtol  = RPAR(1);
  auto Im = mat_fun::mat_id(nX);

  Vector<double> fn(nX);
  getf(imyo, nX, nG, Xn, Xg, fn, Istim, Ksac, RPAR);

  int k  = 0;
  auto Xk = Xn;
  bool l1 = false;
  bool l2 = false;
  bool l3 = false;
  double t = Ts + dt;
  double eps = std::numeric_limits<double>::epsilon();

  while (true) {
    k = k + 1;
    Vector<double> fk(nX);
    getf(imyo, nX, nG, Xn, Xg, fn, Istim, Ksac, RPAR);

    auto rK = Xk - Xn - 0.5*dt*(fk + fn);

    double rmsA = 0.0;
    double rmsR = 0.0;

    for (int i = 0; i < nX; i++) {
      rmsA = rmsA + pow(rK(i),2.0);
      rmsR = rmsR + pow(rK(i) / (Xk(i)+eps), 2.0);
    }

    rmsA = sqrt(rmsA / static_cast<double>(nX));
    rmsR = sqrt(rmsR / static_cast<double>(nX));

    l1 = (k > itMax);
    l2 = (rmsA <= atol);
    l3  = (rmsR <= rtol);
    if (l1 || l2 || l3) {
      break;
    }

    Array<double> JAC(nX,nX);
    getj(imyo, nX, nG, Xk, Xg, JAC, Ksac);

    JAC = Im - 0.5*dt*JAC;
    JAC = mat_fun::mat_inv(JAC, nX);
    rK  = mat_fun::mat_mul(JAC, rK);
    Xk  = Xk - rK;
  }

  Xn = Xk;

  update_g(imyo, dt, nX, nG, Xn, Xg);
  getf(imyo, nX, nG, Xn, Xg, fn, Istim, Ksac, RPAR);

  if (!l2 && !l3) {
    IPAR(1) = IPAR(1) + 1;
  }
}

void CepModTtp::integ_fe(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg, 
    const double Ts, const double dt, const double Istim, const double Ksac, Vector<double>& RPAR)
{
  Vector<double> f(nX);

  // Get time derivatives (RHS)
  getf(imyo, nX, nG, X, Xg, f, Istim, Ksac, RPAR);
  //CALL TTP_GETF(imyo, nX, nG, X, Xg, f, Istim, Ksac, RPAR)

  // Update gating variables
  update_g(imyo, dt, nX, nG, X, Xg);
  //CALL TTP_UPDATEG(imyo, dt, nX, nG, X, Xg)

  //  Update state variables
  X = X + dt*f;
}

void CepModTtp::integ_rk(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg, 
    const double Ts, const double dt, const double Istim, const double Ksac, Vector<double>& RPAR)
{
  double dt6 = dt / 6.0;

  Vector<double> frk1(nX), frk2(nX), frk3(nX), frk4(nX);

  // RK4: 1st pass
  auto Xrk = X;
  getf(imyo, nX, nG, Xrk, Xg, frk1, Istim, Ksac, RPAR);

  // Update gating variables by half-dt
  auto Xgr = Xg;
  update_g(imyo, 0.5*dt, nX, nG, X, Xgr);

  // RK4: 2nd pass
  Xrk = X + 0.5*dt*frk1;
  getf(imyo, nX, nG, Xrk, Xgr, frk2, Istim, Ksac, RPAR);

  // RK4: 3rd pass
  Xrk = X + 0.5*dt*frk2;
  getf(imyo, nX, nG, Xrk, Xgr, frk3, Istim, Ksac, RPAR);

  // Update gating variables by full-dt
  Xgr = Xg;
  update_g(imyo, dt, nX, nG, X, Xgr);

  // RK4: 4th pass
  Xrk = X + dt*frk3;
  getf(imyo, nX, nG, Xrk, Xgr, frk4, Istim, Ksac, RPAR);

  X = X + dt6 * (frk1 + 2.0*(frk2 + frk3) + frk4);
  Xg = Xgr;
}

/// @brief Update all the gating variables
void CepModTtp::update_g(const int i, const double dt, const int n, const int nG, const Vector<double>& X, Vector<double>& Xg)
{
  V  = X(0);
  Ca_ss = X(4);

  xr1   = Xg(0);
  xr2   = Xg(1);
  xs    = Xg(2);
  m     = Xg(3);
  h     = Xg(4);
  j     = Xg(5);
  d     = Xg(6);
  f     = Xg(7);
  f2    = Xg(8);
  fcass = Xg(9);
  s     = Xg(10);
  r     = Xg(11);

  double a, b, c, tau;

  // xr1: activation gate for I_Kr
  xr1i = 1.0/(1.0 + exp(-(26.0+V)/7.0));
  a = 450.0/(1.0 + exp(-(45.0+V)/10.0));
  b = 6.0/(1.0 + exp((30.0+V)/11.50));
  tau = a*b;
  Xg(0) = xr1i - (xr1i - xr1)*exp(-dt/tau);

  // xr2: inactivation gate for I_Kr
  xr2i = 1.0 /(1.0 + exp((88.0+V)/24.0));
  a = 3.0 /(1.0 + exp(-(60.0+V)/20.0));
  b = 1.120/(1.0 + exp(-(60.0-V)/20.0));
  tau = a*b;
  Xg(1) = xr2i - (xr2i - xr2)*exp(-dt/tau);

  // xs: activation gate for I_Ks
  xsi = 1.0/(1.0 + exp(-(5.0+V)/14.0));
  a = 1400.0/sqrt(1.0 + exp((5.0-V)/6.0));
  b = 1.0/(1.0 + exp((V-35.0)/15.0));
  tau  = a*b + 80.0;
  Xg(2)  = xsi - (xsi - xs)*exp(-dt/tau);

  // m: activation gate for I_Na
  mi = 1.0 / pow(1.0 + exp(-(56.860+V)/9.030),2.0);
  //mi = 1.0/( (1.0 + exp(-(56.860+V)/9.030))**2.0 )
  a = 1.0/(1.0 + exp(-(60.0+V)/5.0));
  b = 0.10/(1.0 + exp((35.0+V)/5.0)) + 0.10/(1.0 + exp((V-50.0)/200.0));
  tau = a*b;
  Xg(3) = mi - (mi - m)*exp(-dt/tau);

  // h: fast inactivation gate for I_Na
  hi = 1.0 / pow(1.0 + exp((71.550+V)/7.430),2.0);
  //hi = 1.0/( (1.0 + exp((71.550+V)/7.430))**2.0 )

  if (V >= -40.0) {
    a = 0.0;
    b = 0.770/(0.130*(1.0 + exp(-(10.660+V)/11.10)));
  } else {
    a = 5.7E-2*exp(-(80.0+V)/6.80);
    b = 2.70*exp(0.0790*V) + 310000.0*exp(0.34850*V);
  }

  tau  = 1.0 / (a + b);
  Xg(4)  = hi - (hi - h)*exp(-dt/tau);

  // j: slow inactivation gate for I_Na
  ji = 1.0/ pow(1.0 + exp((71.550+V)/7.430),2.0);
  //ji = 1.0/( (1.0 + EXP((71.550+V)/7.430))**2.0 )

  if (V >= -40.0) {
    a = 0.0;
    b = 0.60*exp(5.7E-2*V) / (1.0 + exp(-0.10*(V+32.0)));
  } else {
    a = -(25428.0*exp(0.24440*V) + 6.948E-6*exp(-0.043910*V)) * (V+37.780) / (1.0 + exp(0.3110*(79.230+V)));
    b = 0.024240*exp(-0.010520*V) / (1.0 + exp(-0.13780*(40.140+V)));
  }
  tau = 1.0 / (a + b);
  Xg(5)  = ji - (ji - j)*exp(-dt/tau);

  // d: activation gate for I_CaL
  di = 1.0/(1.0 + exp(-(8.0+V)/7.50));
  a = 1.40/(1.0 + exp(-(35.0+V)/13.0)) + 0.250;
  b = 1.40/(1.0 + exp((5.0+V)/5.0));
  c = 1.0/(1.0 + exp((50.0-V)/20.0));
  tau = a*b + c;
  Xg(6)  = di - (di - d)*exp(-dt/tau);

  // f: slow inactivation gate for I_CaL
  fi = 1.0/(1.0 + exp((20.0+V)/7.0));
  a = 1102.50*exp(-pow(V+27.0,2.0) / 225.0);
  //a = 1102.50*exp(-((V+27.0)**2.0)/225.0);
  b = 200.0/(1.0 + exp((13.0-V)/10.0));
  c = 180.0/(1.0 + exp((30.0+V)/10.0)) + 20.0;
  tau = a + b + c;
  // for spiral wave breakup
  // if (V .GT. 0.0) tau = tau*2.0
  Xg(7)  = fi - (fi - f)*exp(-dt/tau);

  // f2: fast inactivation gate for I_CaL
  f2i = 0.670/(1.0 + exp((35.0+V)/7.0)) + 0.330;
  a = 562.0*exp(-pow(27.0+V,2.0) /240.0);
  //a = 562.0*exp(-((27.0+V)**2.0) /240.0)
  b = 31.0/(1.0 + exp((25.0-V)/10.0));
  c = 80.0/(1.0 + exp((30.0+V)/10.0));
  tau = a + b + c;
  Xg(8)  = f2i - (f2i - f2)*exp(-dt/tau);

  // fCass: inactivation gate for I_CaL into subspace
  // = 1.0/(1.0 + (Ca_ss/0.050)**2.0)
  c = 1.0 / (1.0 + pow(Ca_ss/0.050,2.0));
  fcassi = 0.60*c  + 0.40;
  tau = 80.0*c + 2.0;
  Xg(9) = fcassi - (fcassi - fcass)*exp(-dt/tau);

  // s: inactivation gate for I_to
  if (i == 1 || i == 3) {
     si = 1.0/(1.0 + exp((20.0+V)/5.0));
     tau = 85.0*exp(-pow(V+45.0,2.0) / 320.0) + 5.0/(1.0+exp((V-20.0)/5.0)) + 3.0;
  } else if (i  ==  2) {
     si = 1.0/(1.0 + exp((28.0+V)/5.0));
     tau = 1000.0*exp(-pow(V+67.0,2.0) /1000.0) + 8.0;
     //tau = 1000.0*exp(-((V+67.0)**2.0) /1000.0) + 8.0;
  }
  Xg(10) = si - (si - s)*exp(-dt/tau);

  // r: activation gate for I_to
  ri = 1.0/(1.0 + exp((20.0-V)/6.0));
  tau = 9.50*exp(-pow(V+40.0,2.0) / 1800.0) + 0.80;
  //tau = 9.50*exp(-((V+40.0)**2.0) /1800.0) + 0.80
  Xg(11) = ri - (ri - r)*exp(-dt/tau);
}



