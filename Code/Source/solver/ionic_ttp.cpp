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


REGISTER_IONIC_MODEL("TTP", TTP);