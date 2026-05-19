// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ionic_ttp.h"

void TTP::read_parameters(const IonicModelParameters &params) {
  IonicModel::read_parameters(params);

  Rc = params.get_scalar("Rc");
  Tc = params.get_scalar("Tc");
  Fc = params.get_scalar("Fc");
  Cm = params.get_scalar("Cm");
  sV = params.get_scalar("sV");
  rho = params.get_scalar("rho");
  V_c = params.get_scalar("V_c");
  V_sr = params.get_scalar("V_sr");
  V_ss = params.get_scalar("V_ss");
  K_o = params.get_scalar("K_o");
  Na_o = params.get_scalar("Na_o");
  Ca_o = params.get_scalar("Ca_o");
  G_Na = params.get_scalar("G_Na");
  G_K1 = params.get_scalar("G_K1");
  G_to = params.get_scalar("G_to");
  G_Kr = params.get_scalar("G_Kr");
  G_Ks = params.get_scalar("G_Ks");
  p_KNa = params.get_scalar("p_KNa");
  G_CaL = params.get_scalar("G_CaL");
  K_NaCa = params.get_scalar("K_NaCa");
  gamma = params.get_scalar("gamma");
  K_mCa = params.get_scalar("K_mCa");
  K_mNai = params.get_scalar("K_mNai");
  K_sat = params.get_scalar("K_sat");
  alpha = params.get_scalar("alpha");
  p_NaK = params.get_scalar("p_NaK");
  K_mK = params.get_scalar("K_mK");
  K_mNa = params.get_scalar("K_mNa");
  G_pK = params.get_scalar("G_pK");
  G_pCa = params.get_scalar("G_pCa");
  K_pCa = params.get_scalar("K_pCa");
  G_bNa = params.get_scalar("G_bNa");
  G_bCa = params.get_scalar("G_bCa");
  Vmax_up = params.get_scalar("Vmax_up");
  K_up = params.get_scalar("K_up");
  V_rel = params.get_scalar("V_rel");
  k1p = params.get_scalar("k1p");
  k2p = params.get_scalar("k2p");
  k3 = params.get_scalar("k3");
  k4 = params.get_scalar("k4");
  EC = params.get_scalar("EC");
  max_sr = params.get_scalar("max_sr");
  min_sr = params.get_scalar("min_sr");
  V_leak = params.get_scalar("V_leak");
  V_xfer = params.get_scalar("V_xfer");
  Buf_c = params.get_scalar("Buf_c");
  K_bufc = params.get_scalar("K_bufc");
  Buf_sr = params.get_scalar("Buf_sr");
  K_bufsr = params.get_scalar("K_bufsr");
  Buf_ss = params.get_scalar("Buf_ss");
  K_bufss = params.get_scalar("K_bufss");
}

void TTP::distribute_parameters(const CmMod &cm_mod, const cmType &cm) {
  IonicModel::distribute_parameters(cm_mod, cm);

  cm.bcast(cm_mod, &Rc);
  cm.bcast(cm_mod, &Tc);
  cm.bcast(cm_mod, &Fc);
  cm.bcast(cm_mod, &Cm);
  cm.bcast(cm_mod, &sV);
  cm.bcast(cm_mod, &rho);
  cm.bcast(cm_mod, &V_c);
  cm.bcast(cm_mod, &V_sr);
  cm.bcast(cm_mod, &V_ss);
  cm.bcast(cm_mod, &K_o);
  cm.bcast(cm_mod, &Na_o);
  cm.bcast(cm_mod, &Ca_o);
  cm.bcast(cm_mod, &G_Na);
  cm.bcast(cm_mod, &G_K1);
  cm.bcast(cm_mod, &G_to);
  cm.bcast(cm_mod, &G_Kr);
  cm.bcast(cm_mod, &G_Ks);
  cm.bcast(cm_mod, &p_KNa);
  cm.bcast(cm_mod, &G_CaL);
  cm.bcast(cm_mod, &K_NaCa);
  cm.bcast(cm_mod, &gamma);
  cm.bcast(cm_mod, &K_mCa);
  cm.bcast(cm_mod, &K_mNai);
  cm.bcast(cm_mod, &K_sat);
  cm.bcast(cm_mod, &alpha);
  cm.bcast(cm_mod, &p_NaK);
  cm.bcast(cm_mod, &K_mK);
  cm.bcast(cm_mod, &K_mNa);
  cm.bcast(cm_mod, &G_pK);
  cm.bcast(cm_mod, &G_pCa);
  cm.bcast(cm_mod, &K_pCa);
  cm.bcast(cm_mod, &G_bNa);
  cm.bcast(cm_mod, &G_bCa);
  cm.bcast(cm_mod, &Vmax_up);
  cm.bcast(cm_mod, &K_up);
  cm.bcast(cm_mod, &V_rel);
  cm.bcast(cm_mod, &k1p);
  cm.bcast(cm_mod, &k2p);
  cm.bcast(cm_mod, &k3);
  cm.bcast(cm_mod, &k4);
  cm.bcast(cm_mod, &EC);
  cm.bcast(cm_mod, &max_sr);
  cm.bcast(cm_mod, &min_sr);
  cm.bcast(cm_mod, &V_leak);
  cm.bcast(cm_mod, &V_xfer);
  cm.bcast(cm_mod, &Buf_c);
  cm.bcast(cm_mod, &K_bufc);
  cm.bcast(cm_mod, &Buf_sr);
  cm.bcast(cm_mod, &K_bufsr);
  cm.bcast(cm_mod, &Buf_ss);
  cm.bcast(cm_mod, &K_bufss);
}

void TTP::update_g(const unsigned int zone_id, const double dt,
                   const Vector<double> &X, Vector<double> &Xg) const {
  // Local copies of state and gating variables.
  const double V = X(0);
  const double Ca_ss = X(4);
  const double xr1 = Xg(0);
  const double xr2 = Xg(1);
  const double xs = Xg(2);
  const double m = Xg(3);
  const double h = Xg(4);
  const double j = Xg(5);
  const double d = Xg(6);
  const double f = Xg(7);
  const double f2 = Xg(8);
  const double fcass = Xg(9);
  const double s = Xg(10);
  const double r = Xg(11);

  // xr1: activation gate for I_Kr
  {
    const double xr1i = 1.0 / (1.0 + exp(-(26.0 + V) / 7.0));
    const double a = 450.0 / (1.0 + exp(-(45.0 + V) / 10.0));
    const double b = 6.0 / (1.0 + exp((30.0 + V) / 11.50));
    const double tau = a * b;

    Xg(0) = xr1i - (xr1i - xr1) * exp(-dt / tau);
  }

  // xr2: inactivation gate for I_Kr
  {
    const double xr2i = 1.0 / (1.0 + exp((88.0 + V) / 24.0));
    const double a = 3.0 / (1.0 + exp(-(60.0 + V) / 20.0));
    const double b = 1.120 / (1.0 + exp(-(60.0 - V) / 20.0));
    const double tau = a * b;

    Xg(1) = xr2i - (xr2i - xr2) * exp(-dt / tau);
  }

  // xs: activation gate for I_Ks
  {
    const double xsi = 1.0 / (1.0 + exp(-(5.0 + V) / 14.0));
    const double a = 1400.0 / sqrt(1.0 + exp((5.0 - V) / 6.0));
    const double b = 1.0 / (1.0 + exp((V - 35.0) / 15.0));
    const double tau = a * b + 80.0;

    Xg(2) = xsi - (xsi - xs) * exp(-dt / tau);
  }

  // m: activation gate for I_Na
  {
    const double mi = 1.0 / pow(1.0 + exp(-(56.860 + V) / 9.030), 2.0);
    // mi = 1.0/( (1.0 + exp(-(56.860+V)/9.030))**2.0 )
    const double a = 1.0 / (1.0 + exp(-(60.0 + V) / 5.0));
    const double b = 0.10 / (1.0 + exp((35.0 + V) / 5.0)) +
                     0.10 / (1.0 + exp((V - 50.0) / 200.0));
    const double tau = a * b;

    Xg(3) = mi - (mi - m) * exp(-dt / tau);
  }

  // h: fast inactivation gate for I_Na
  {
    const double hi = 1.0 / pow(1.0 + exp((71.550 + V) / 7.430), 2.0);
    // hi = 1.0/( (1.0 + exp((71.550+V)/7.430))**2.0 )

    double a, b;
    if (V >= -40.0) {
      a = 0.0;
      b = 0.770 / (0.130 * (1.0 + exp(-(10.660 + V) / 11.10)));
    } else {
      a = 5.7E-2 * exp(-(80.0 + V) / 6.80);
      b = 2.70 * exp(0.0790 * V) + 310000.0 * exp(0.34850 * V);
    }

    const double tau = 1.0 / (a + b);

    Xg(4) = hi - (hi - h) * exp(-dt / tau);
  }

  // j: slow inactivation gate for I_Na
  {
    const double ji = 1.0 / pow(1.0 + exp((71.550 + V) / 7.430), 2.0);
    // ji = 1.0/( (1.0 + EXP((71.550+V)/7.430))**2.0 )

    double a, b;
    if (V >= -40.0) {
      a = 0.0;
      b = 0.60 * exp(5.7E-2 * V) / (1.0 + exp(-0.10 * (V + 32.0)));
    } else {
      a = -(25428.0 * exp(0.24440 * V) + 6.948E-6 * exp(-0.043910 * V)) *
          (V + 37.780) / (1.0 + exp(0.3110 * (79.230 + V)));
      b = 0.024240 * exp(-0.010520 * V) / (1.0 + exp(-0.13780 * (40.140 + V)));
    }
    const double tau = 1.0 / (a + b);

    Xg(5) = ji - (ji - j) * exp(-dt / tau);
  }

  // d: activation gate for I_CaL
  {
    const double di = 1.0 / (1.0 + exp(-(8.0 + V) / 7.50));
    const double a = 1.40 / (1.0 + exp(-(35.0 + V) / 13.0)) + 0.250;
    const double b = 1.40 / (1.0 + exp((5.0 + V) / 5.0));
    const double c = 1.0 / (1.0 + exp((50.0 - V) / 20.0));
    const double tau = a * b + c;

    Xg(6) = di - (di - d) * exp(-dt / tau);
  }

  // f: slow inactivation gate for I_CaL
  {
    const double fi = 1.0 / (1.0 + exp((20.0 + V) / 7.0));
    const double a = 1102.50 * exp(-pow(V + 27.0, 2.0) / 225.0);
    // a = 1102.50*exp(-((V+27.0)**2.0)/225.0);
    const double b = 200.0 / (1.0 + exp((13.0 - V) / 10.0));
    const double c = 180.0 / (1.0 + exp((30.0 + V) / 10.0)) + 20.0;
    const double tau = a + b + c;
    // for spiral wave breakup
    // if (V .GT. 0.0) tau = tau*2.0

    Xg(7) = fi - (fi - f) * exp(-dt / tau);
  }

  // f2: fast inactivation gate for I_CaL
  {
    const double f2i = 0.670 / (1.0 + exp((35.0 + V) / 7.0)) + 0.330;
    const double a = 562.0 * exp(-pow(27.0 + V, 2.0) / 240.0);
    // a = 562.0*exp(-((27.0+V)**2.0) /240.0)
    const double b = 31.0 / (1.0 + exp((25.0 - V) / 10.0));
    const double c = 80.0 / (1.0 + exp((30.0 + V) / 10.0));
    const double tau = a + b + c;

    Xg(8) = f2i - (f2i - f2) * exp(-dt / tau);
  }

  // fCass: inactivation gate for I_CaL into subspace
  // = 1.0/(1.0 + (Ca_ss/0.050)**2.0)
  {
    const double c = 1.0 / (1.0 + pow(Ca_ss / 0.050, 2.0));
    const double fcassi = 0.60 * c + 0.40;
    const double tau = 80.0 * c + 2.0;

    Xg(9) = fcassi - (fcassi - fcass) * exp(-dt / tau);
  }

  // s: inactivation gate for I_to
  {
    double si, tau;

    if (zone_id == 1 || zone_id == 3) {
      si = 1.0 / (1.0 + exp((20.0 + V) / 5.0));
      tau = 85.0 * exp(-pow(V + 45.0, 2.0) / 320.0) +
            5.0 / (1.0 + exp((V - 20.0) / 5.0)) + 3.0;
    } else if (zone_id == 2) {
      si = 1.0 / (1.0 + exp((28.0 + V) / 5.0));
      tau = 1000.0 * exp(-pow(V + 67.0, 2.0) / 1000.0) + 8.0;
      // tau = 1000.0*exp(-((V+67.0)**2.0) /1000.0) + 8.0;
    }

    Xg(10) = si - (si - s) * exp(-dt / tau);
  }

  // r: activation gate for I_to
  {
    const double ri = 1.0 / (1.0 + exp((20.0 - V) / 6.0));
    const double tau = 9.50 * exp(-pow(V + 40.0, 2.0) / 1800.0) + 0.80;
    // tau = 9.50*exp(-((V+40.0)**2.0) /1800.0) + 0.80

    Xg(11) = ri - (ri - r) * exp(-dt / tau);
  }
}

Vector<double> TTP::getf(const unsigned int zone_id, const Vector<double> &X,
                         const Vector<double> &Xg, const double I_stim,
                         const double I_sac) const {
  // Local copies of state variables
  const double V = X(0);
  const double K_i = X(1);
  const double Na_i = X(2);
  const double Ca_i = X(3);
  const double Ca_ss = X(4);
  const double Ca_sr = X(5);
  const double R_bar = X(6);
  const double xr1 = Xg(0);
  const double xr2 = Xg(1);
  const double xs = Xg(2);
  const double m = Xg(3);
  const double h = Xg(4);
  const double j = Xg(5);
  const double d = Xg(6);
  const double f = Xg(7);
  const double f2 = Xg(8);
  const double fcass = Xg(9);
  const double s = Xg(10);
  const double r = Xg(11);

  // Diff = 1. / (1.0D1 * rho * Cm * sV)
  const double RT = Rc * Tc / Fc;
  const double E_K = RT * log(K_o / K_i);
  const double E_Na = RT * log(Na_o / Na_i);
  const double E_Ca = 0.5 * RT * log(Ca_o / Ca_i);
  const double E_Ks = RT * log((K_o + p_KNa * Na_o) / (K_i + p_KNa * Na_i));
  const double sq5 = sqrt(K_o / 5.4);
  const double k_casr =
      max_sr - ((max_sr - min_sr) / (1.0 + pow(EC / Ca_sr, 2.0)));

  // I_Na: Fast sodium current
  const double I_Na = G_Na * pow(m, 3.0) * h * j * (V - E_Na);

  // I_to: transient outward current
  const double I_to = G_to * r * s * (V - E_K);

  // I_K1: inward rectifier outward current
  double I_K1;
  {
    const double e1 = exp(0.06 * (V - E_K - 200.0));
    const double e2 = exp(2.E-4 * (V - E_K + 100.0));
    const double e3 = exp(0.1 * (V - E_K - 10.0));
    const double e4 = exp(-0.5 * (V - E_K));
    const double a = 0.1 / (1.0 + e1);
    const double b = (3.0 * e2 + e3) / (1.0 + e4);
    const double tau = a / (a + b);
    I_K1 = G_K1 * sq5 * tau * (V - E_K);
  }

  // I_Kr: rapid delayed rectifier current
  const double I_Kr = G_Kr * sq5 * xr1 * xr2 * (V - E_K);

  // I_Ks: slow delayed rectifier current
  const double I_Ks = G_Ks * pow(xs, 2.0) * (V - E_Ks);

  // I_CaL: L-type Ca current
  double I_CaL;
  {
    const double a = 2.0 * (V - 15.) / RT;
    const double b =
        2.0 * a * Fc * (0.25 * Ca_ss * exp(a) - Ca_o) / (exp(a) - 1.0);
    I_CaL = G_CaL * d * f * f2 * fcass * b;
  }

  // I_NaCa: Na-Ca exchanger current
  double I_NaCa;
  {
    const double e1 = exp(gamma * V / RT);
    const double e2 = exp((gamma - 1.) * V / RT);
    const double n1 =
        e1 * pow(Na_i, 3.0) * Ca_o - e2 * pow(Na_o, 3.0) * Ca_i * alpha;
    const double d1 = pow(K_mNai, 3.0) + pow(Na_o, 3.0);
    const double d2 = K_mCa + Ca_o;
    const double d3 = 1.0 + K_sat * e2;
    I_NaCa = K_NaCa * n1 / (d1 * d2 * d3);
  }

  // I_NaK: Na-K pump current
  double I_NaK;
  {
    const double e1 = exp(-0.1 * V / RT);
    const double e2 = exp(-V / RT);
    const double n1 = p_NaK * K_o * Na_i;
    const double d1 = K_o + K_mK;
    const double d2 = Na_i + K_mNa;
    const double d3 = 1.0 + 0.1245 * e1 + 0.0353 * e2;
    I_NaK = n1 / (d1 * d2 * d3);
  }

  // I_pCa: plateau Ca current
  const double I_pCa = G_pCa * Ca_i / (K_pCa + Ca_i);

  // I_pK: plateau K current
  const double I_pK = G_pK * (V - E_K) / (1.0 + exp((25.0 - V) / 5.98));

  // I_bCa: background Ca current
  const double I_bCa = G_bCa * (V - E_Ca);

  // I_bNa: background Na current
  const double I_bNa = G_bNa * (V - E_Na);

  // I_leak: Sacroplasmic Reticulum Ca leak current
  const double I_leak = V_leak * (Ca_sr - Ca_i);

  // I_up: Sacroplasmic Reticulum Ca pump current
  const double I_up = Vmax_up / (1.0 + pow(K_up / Ca_i, 2.0));

  // I_rel: Ca induced Ca current (CICR)
  double I_rel;
  {
    const double k1 = k1p / k_casr;
    const double O = k1 * R_bar * pow(Ca_ss, 2.0) / (k3 + k1 * pow(Ca_ss, 2.0));
    I_rel = V_rel * O * (Ca_sr - Ca_ss);
  }

  //  I_xfer: diffusive Ca current between Ca subspae and cytoplasm
  const double I_xfer = V_xfer * (Ca_ss - Ca_i);

  // Now compute time derivatives
  Vector<double> dX(X.size());

  // dV/dt: rate of change of transmembrane voltage
  dX(0) = -(I_Na + I_to + I_K1 + I_Kr + I_Ks + I_CaL + I_NaCa + I_NaK + I_pCa +
            I_pK + I_bCa + I_bNa + I_stim) +
          I_sac;

  // dK_i/dt
  dX(1) = -(Cm / (V_c * Fc)) *
          (I_K1 + I_to + I_Kr + I_Ks + I_pK - 2.0 * I_NaK + I_stim);

  //  dNa_i/dt
  dX(2) = -(Cm / (V_c * Fc)) * (I_Na + I_bNa + 3.0 * (I_NaK + I_NaCa));

  // dCa_i/dt
  {
    const double n1 = (I_leak - I_up) * V_sr / V_c + I_xfer;
    const double n2 = -(Cm / (V_c * Fc)) * (I_bCa + I_pCa - 2.0 * I_NaCa) / 2.0;
    const double d1 = 1.0 + K_bufc * Buf_c / pow(Ca_i + K_bufc, 2.0);
    dX(3) = (n1 + n2) / d1;
  }

  // dCa_ss: rate of change of Ca_ss
  {
    const double n1 =
        (-I_CaL * Cm / (2.0 * Fc) + I_rel * V_sr - V_c * I_xfer) / V_ss;
    const double d1 = 1.0 + K_bufss * Buf_ss / pow(Ca_ss + K_bufss, 2.0);
    dX(4) = n1 / d1;
  }

  // dCa_sr: rate of change of Ca_sr
  {
    const double n1 = I_up - I_leak - I_rel;
    const double d1 = 1. + K_bufsr * Buf_sr / pow(Ca_sr + K_bufsr, 2.0);
    dX(5) = n1 / d1;
  }

  // Rbar: ryanodine receptor
  {
    const double k2 = k2p * k_casr;
    dX(6) = -k2 * Ca_ss * R_bar + k4 * (1.0 - R_bar);
  }

  return dX;
}

Array<double> TTP::getj(const unsigned int zone_id, const Vector<double> &X,
                        const Vector<double> &Xg, const double Ksac) const {
  double a, b, c, tau, sq5, e1, e2, e3, e4, n1, n2, d1, d2, d3;

  // Local copies of state and gating variables
  const double V = X(0);
  const double K_i = X(1);
  const double Na_i = X(2);
  const double Ca_i = X(3);
  const double Ca_ss = X(4);
  const double Ca_sr = X(5);
  const double R_bar = X(6);
  const double xr1 = Xg(0);
  const double xr2 = Xg(1);
  const double xs = Xg(2);
  const double m = Xg(3);
  const double h = Xg(4);
  const double j = Xg(5);
  const double d = Xg(6);
  const double f = Xg(7);
  const double f2 = Xg(8);
  const double fcass = Xg(9);
  const double s = Xg(10);
  const double r = Xg(11);

  const double RT = Rc * Tc / Fc;
  const double E_K = RT * log(K_o / K_i);
  const double E_Na = RT * log(Na_o / Na_i);
  const double E_Ca = 0.5 * RT * log(Ca_o / Ca_i);
  const double E_Ks = RT * log((K_o + p_KNa * Na_o) / (K_i + p_KNa * Na_i));

  const double E_K_Ki = -RT / K_i;
  const double E_Na_Nai = -RT / Na_i;
  const double E_Ca_Cai = -RT / Ca_i / 2.0;
  const double E_Ks_Ki = -RT / (K_i + p_KNa * Na_i);
  const double E_Ks_Nai = p_KNa * E_Ks_Ki;

  // I_Na: Fast sodium current
  const double I_Na = G_Na * pow(m, 3.0) * h * j * (V - E_Na);
  const double I_Na_V = G_Na * pow(m, 3.0) * h * j;
  const double I_Na_Nai = I_Na_V * (-E_Na_Nai);

  // I_to: transient outward current
  const double I_to = G_to * r * s * (V - E_K);
  const double I_to_V = G_to * r * s;
  const double I_to_Ki = I_to_V * (-E_K_Ki);

  // I_K1: inward rectifier outward current
  e1 = exp(0.060 * (V - E_K - 200.0));
  e2 = exp(2.E-40 * (V - E_K + 100.0));
  e3 = exp(0.10 * (V - E_K - 10.0));
  e4 = exp(-0.50 * (V - E_K));
  a = 0.10 / (1.0 + e1);
  b = (3.0 * e2 + e3) / (1.0 + e4);
  tau = a / (a + b);
  sq5 = sqrt(K_o / 5.40);
  n1 = -6.E-30 * e1 / pow(1.0 + e1, 2.0);
  n2 = (6.E-40 * e2 + 0.10 * e3 + 0.50 * b * e4) / (1.0 + e4);
  n1 = (a + b) * n1 - a * (n1 + n2);
  d1 = pow(a + b, 2.0);
  const double I_K1 = G_K1 * sq5 * tau * (V - E_K);
  const double I_K1_V = G_K1 * sq5 * (tau + (V - E_K) * n1 / d1);
  const double I_K1_Ki = I_K1_V * (-E_K_Ki);

  // I_Kr: rapid delayed rectifier current
  const double I_Kr = G_Kr * sq5 * xr1 * xr2 * (V - E_K);
  const double I_Kr_V = G_Kr * sq5 * xr1 * xr2;
  const double I_Kr_Ki = I_Kr_V * (-E_K_Ki);

  // I_Ks: slow delayed rectifier current
  const double I_Ks = G_Ks * pow(xs, 2.0) * (V - E_Ks);
  const double I_Ks_V = G_Ks * pow(xs, 2.0);
  const double I_Ks_Ki = I_Ks_V * (-E_Ks_Ki);
  const double I_Ks_Nai = I_Ks_V * (-E_Ks_Nai);

  //  I_CaL: L-type Ca current
  a = 2.0 * (V - 15.0) / RT;
  b = (0.250 * Ca_ss * exp(a) - Ca_o) / (exp(a) - 1.0);
  c = G_CaL * d * f * f2 * fcass * (2.0 * a * Fc);
  n1 = (exp(a) / RT) / (exp(a) - 1.0);
  const double I_CaL = c * b;
  const double I_CaL_V =
      I_CaL / (V - 15.0) + n1 * (c * 0.50 * Ca_ss - 2.0 * I_CaL);
  const double I_CaL_Cass = c * 0.250 * n1 * RT;

  //  I_NaCa: Na-Ca exchanger current
  e1 = exp(gamma * V / RT);
  e2 = exp((gamma - 1.0) * V / RT);
  n1 = e1 * pow(Na_i, 3.0) * Ca_o - e2 * pow(Na_o, 3.0) * Ca_i * alpha;
  d1 = pow(K_mNai, 3.0) + pow(Na_o, 3.0);
  d2 = K_mCa + Ca_o;
  d3 = 1.0 + K_sat * e2;
  c = 1.0 / (d1 * d2 * d3);
  const double I_NaCa = K_NaCa * n1 * c;

  n1 = K_NaCa * c *
       (e1 * pow(Na_i, 3.0) * Ca_o * (gamma / RT) -
        e2 * pow(Na_o, 3.0) * Ca_i * alpha * ((gamma - 1.0) / RT));
  n2 = I_NaCa * K_sat * ((gamma - 1.0) / RT) * e2 / d3;
  const double I_NaCa_V = n1 - n2;
  const double I_NaCa_Nai = K_NaCa * e1 * (3.0 * pow(Na_i, 2.0)) * Ca_o * c;
  const double I_NaCa_Cai = -K_NaCa * e2 * pow(Na_o, 3.0) * alpha * c;

  // I_NaK: Na-K pump current
  e1 = exp(-0.10 * V / RT);
  e2 = exp(-V / RT);
  n1 = p_NaK * K_o * Na_i;
  d1 = K_o + K_mK;
  d2 = Na_i + K_mNa;
  d3 = 1.0 + 0.12450 * e1 + 0.03530 * e2;
  const double I_NaK = n1 / (d1 * d2 * d3);
  n1 = (0.012450 * e1 + 0.03530 * e2) / RT;
  const double I_NaK_V = I_NaK * n1 / d3;
  const double I_NaK_Nai = I_NaK * K_mNa / (Na_i * d2);

  // I_pCa: plateau Ca current
  d1 = (K_pCa + Ca_i);
  const double I_pCa = G_pCa * Ca_i / d1;
  const double I_pCa_Cai = G_pCa * K_pCa / (d1 * d1);

  //  I_pK: plateau K current
  e1 = exp((25.0 - V) / 5.980);
  const double I_pK = G_pK * (V - E_K) / (1.0 + e1);
  const double I_pK_V = (G_pK + I_pK * e1 / 5.980) / (1.0 + e1);
  const double I_pK_Ki = G_pK * (-E_K_Ki) / (1.0 + e1);

  // I_bCa: background Ca current
  const double I_bCa = G_bCa * (V - E_Ca);
  const double I_bCa_V = G_bCa;
  const double I_bCa_Cai = G_bCa * (-E_Ca_Cai);

  //  I_bNa: background Na current
  const double I_bNa = G_bNa * (V - E_Na);
  const double I_bNa_V = G_bNa;
  const double I_bNa_Nai = G_bNa * (-E_Na_Nai);

  // I_leak: Sacroplasmic Reticulum Ca leak current
  const double I_leak = V_leak * (Ca_sr - Ca_i);
  const double I_leak_Cai = -V_leak;
  const double I_leak_Casr = V_leak;

  // I_up: Sacroplasmic Reticulum Ca pump current
  d1 = 1.0 + pow(K_up / Ca_i, 2.0);
  const double I_up = Vmax_up / d1;
  const double I_up_Cai = (I_up / d1) * (2.0 * pow(K_up, 2.0) / pow(Ca_i, 3.0));

  // I_rel: Ca induced Ca current (CICR)
  n1 = max_sr - min_sr;
  d1 = 1.0 + pow(EC / Ca_sr, 2.0);
  const double k_casr = max_sr - (n1 / d1);
  const double k1 = k1p / k_casr;
  n2 = Ca_ss * 2.0;
  d2 = k3 + k1 * n2;
  const double O = k1 * R_bar * n2 / d2;
  const double I_rel = V_rel * O * (Ca_sr - Ca_ss);

  const double k_casr_sr =
      (n1 / pow(d1, 2.0)) * (2.0 * pow(EC, 2.0) / pow(Ca_sr, 3.0));
  // k_casr_sr = (n1 / (d1**2.0)) * (2.0*EC**2.0 / Ca_sr**3.0);
  const double k1_casr = -k1p * k_casr_sr / pow(k_casr, 2.0);
  const double O_Cass = 2.0 * k3 * O / (Ca_ss * d2);
  const double O_Casr = k1_casr * n2 * (R_bar - O) / d2;
  const double O_Rbar = k1 * n2 / d2;

  const double I_rel_Cass = V_rel * (O_Cass * (Ca_sr - Ca_ss) - O);
  const double I_rel_Casr = V_rel * (O_Casr * (Ca_sr - Ca_ss) + O);
  const double I_rel_Rbar = V_rel * O_Rbar * (Ca_sr - Ca_ss);

  //  I_xfer: diffusive Ca current between Ca subspae and cytoplasm
  const double I_xfer = V_xfer * (Ca_ss - Ca_i);
  const double I_xfer_Cai = -V_xfer;
  const double I_xfer_Cass = V_xfer;

  // Compute Jacobian matrix
  Array<double> Jac(X.size(), X.size());

  c = Cm / (V_c * Fc);

  //  V
  Jac(0, 0) = -(I_Na_V + I_to_V + I_K1_V + I_Kr_V + I_Ks_V + I_CaL_V +
                I_NaCa_V + I_NaK_V + I_pK_V + I_bCa_V + I_bNa_V + Ksac);
  Jac(0, 1) = -(I_to_Ki + I_K1_Ki + I_Kr_Ki + I_Ks_Ki + I_pK_Ki);
  Jac(0, 2) = -(I_Na_Nai + I_Ks_Nai + I_NaCa_Nai + I_NaK_Nai + I_bNa_Nai);
  Jac(0, 3) = -(I_NaCa_Cai + I_pCa_Cai + I_bCa_Cai);
  Jac(0, 4) = -I_CaL_Cass;

  // K_i
  Jac(1, 0) = -c * (I_K1_V + I_to_V + I_Kr_V + I_Ks_V + I_pK_V - 2.0 * I_NaK_V);
  Jac(1, 1) = -c * (I_K1_Ki + I_to_Ki + I_Kr_Ki + I_Ks_Ki + I_pK_Ki);
  Jac(1, 2) = -c * (I_Ks_Nai - 2.0 * I_NaK_Nai);

  // Na_i
  Jac(2, 0) = -c * (I_Na_V + I_bNa_V + 3.0 * (I_NaK_V + I_NaCa_V));
  Jac(2, 2) = -c * (I_Na_Nai + I_bNa_Nai + 3.0 * (I_NaK_Nai + I_NaCa_Nai));
  Jac(2, 3) = -c * (3.0 * I_NaCa_Cai);

  //     Ca_i
  n1 = (I_leak - I_up) * V_sr / V_c + I_xfer -
       0.50 * c * (I_bCa + I_pCa - 2.0 * I_NaCa);
  n2 = (I_leak_Cai - I_up_Cai) * V_sr / V_c + I_xfer_Cai -
       0.50 * c * (I_bCa_Cai + I_pCa_Cai - 2.0 * I_NaCa_Cai);
  d1 = 1.0 + K_bufc * Buf_c / pow(Ca_i + K_bufc, 2.0);
  d2 = 2.0 * K_bufc * Buf_c / pow(Ca_i + K_bufc, 3.0);
  Jac(3, 0) = -c * (I_bCa_V - 2.0 * I_NaCa_V) / 2.0 / d1;
  Jac(3, 2) = c * I_NaCa_Nai / d1;
  Jac(3, 3) = (n2 + n1 * d2 / d1) / d1;
  Jac(3, 4) = I_xfer_Cass / d1;
  Jac(3, 5) = (I_leak_Casr * V_sr / V_c) / d1;

  //     Ca_ss
  a = Cm / (2.0 * Fc * V_ss);
  b = V_sr / V_ss;
  c = V_c / V_ss;
  n1 = -a * I_CaL + b * I_rel - c * I_xfer;
  n2 = -a * I_CaL_Cass + b * I_rel_Cass - c * I_xfer_Cass;
  d1 = 1.0 + K_bufss * Buf_ss / pow(Ca_ss + K_bufss, 2.0);
  d2 = 2.0 * K_bufss * Buf_ss / pow(Ca_ss + K_bufss, 3.0);
  Jac(4, 0) = -a * I_CaL_V / d1;
  Jac(4, 3) = -c * I_xfer_Cai / d1;
  Jac(4, 4) = (n2 + n1 * d2 / d1) / d1;
  Jac(4, 5) = b * I_rel_Casr / d1;
  Jac(4, 6) = b * I_rel_Rbar / d1;

  // Ca_sr
  n1 = I_up - I_leak - I_rel;
  n2 = -(I_leak_Casr + I_rel_Casr);
  d1 = 1.0 + K_bufsr * Buf_sr / pow(Ca_sr + K_bufsr, 2.0);
  d2 = 2.0 * K_bufsr * Buf_sr / pow(Ca_sr + K_bufsr, 3.0);
  Jac(5, 3) = (I_up_Cai - I_leak_Cai) / d1;
  Jac(5, 4) = -I_rel_Cass / d1;
  Jac(5, 5) = (n2 + n1 * d2 / d1) / d1;
  Jac(5, 6) = -I_rel_Rbar / d1;

  // Rbar: ryanodine receptor
  const double k2 = k2p * k_casr;
  Jac(6, 4) = -k2 * R_bar;
  Jac(6, 5) = -(k2p * k_casr_sr) * Ca_ss * R_bar;
  Jac(6, 6) = -(k2 * Ca_ss + k4);

  return Jac;
}

REGISTER_IONIC_MODEL("TTP", TTP);