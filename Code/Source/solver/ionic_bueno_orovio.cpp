// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ionic_bueno_orovio.h"

void BuenoOrovio::read_parameters(const IonicModelParameters &params) {
  IonicModel::read_parameters(params);

  u_o = params.get_vector("u_o");
  u_u = params.get_vector("u_u");
  theta_v = params.get_vector("theta_v");
  theta_w = params.get_vector("theta_w");
  thetam_v = params.get_vector("thetam_v");
  theta_o = params.get_vector("theta_o");
  taum_v1 = params.get_vector("taum_v1");
  taum_v2 = params.get_vector("taum_v2");
  taup_v = params.get_vector("taup_v");
  taum_w1 = params.get_vector("taum_w1");
  taum_w2 = params.get_vector("taum_w2");
  km_w = params.get_vector("km_w");
  um_w = params.get_vector("um_w");
  taup_w = params.get_vector("taup_w");
  tau_fi = params.get_vector("tau_fi");
  tau_o1 = params.get_vector("tau_o1");
  tau_o2 = params.get_vector("tau_o2");
  tau_so1 = params.get_vector("tau_so1");
  tau_so2 = params.get_vector("tau_so2");
  k_so = params.get_vector("k_so");
  u_so = params.get_vector("u_so");
  tau_s1 = params.get_vector("tau_s1");
  tau_s2 = params.get_vector("tau_s2");
  k_s = params.get_vector("k_s");
  u_s = params.get_vector("u_s");
  tau_si = params.get_vector("tau_si");
  tau_winf = params.get_vector("tau_winf");
  ws_inf = params.get_vector("ws_inf");
}

void BuenoOrovio::distribute_parameters(const CmMod &cm_mod, const cmType &cm) {
  IonicModel::distribute_parameters(cm_mod, cm);

  cm.bcast(cm_mod, u_o);
  cm.bcast(cm_mod, u_u);
  cm.bcast(cm_mod, theta_v);
  cm.bcast(cm_mod, theta_w);
  cm.bcast(cm_mod, thetam_v);
  cm.bcast(cm_mod, theta_o);
  cm.bcast(cm_mod, taum_v1);
  cm.bcast(cm_mod, taum_v2);
  cm.bcast(cm_mod, taup_v);
  cm.bcast(cm_mod, taum_w1);
  cm.bcast(cm_mod, taum_w2);
  cm.bcast(cm_mod, km_w);
  cm.bcast(cm_mod, um_w);
  cm.bcast(cm_mod, taup_w);
  cm.bcast(cm_mod, tau_fi);
  cm.bcast(cm_mod, tau_o1);
  cm.bcast(cm_mod, tau_o2);
  cm.bcast(cm_mod, tau_so1);
  cm.bcast(cm_mod, tau_so2);
  cm.bcast(cm_mod, k_so);
  cm.bcast(cm_mod, u_so);
  cm.bcast(cm_mod, tau_s1);
  cm.bcast(cm_mod, tau_s2);
  cm.bcast(cm_mod, k_s);
  cm.bcast(cm_mod, u_s);
  cm.bcast(cm_mod, tau_si);
  cm.bcast(cm_mod, tau_winf);
  cm.bcast(cm_mod, ws_inf);
}

void BuenoOrovio::getf(const unsigned int zone_id, const int nX, const int nG,
                       const Vector<double> &X, const Vector<double> &Xg,
                       Vector<double> &f, const double I_stim,
                       const double I_sac) const {
  // Create local copies of the state variables
  const double u = X(0);
  const double v = X(1);
  const double w = X(2);
  const double s = X(3);

  // @todo[michelebucelli] Why minus 1? In getj this seems not to be the case.
  const int i = zone_id - 1;

  // Define step functions
  const double H_uv = step(u - theta_v[i]);
  const double H_uw = step(u - theta_w[i]);
  const double H_umv = step(u - thetam_v[i]);
  const double H_uo = step(u - theta_o[i]);

  // Define additional constants
  const double taum_v = (1.0 - H_umv) * taum_v1[i] + H_umv * taum_v2[i];
  const double taum_w = taum_w1[i] + 0.5 * (taum_w2[i] - taum_w1[i]) *
                                         (1.0 + tanh(km_w[i] * (u - um_w[i])));
  const double tau_so = tau_so1[i] + 0.5 * (tau_so2[i] - tau_so1[i]) *
                                         (1.0 + tanh(k_so[i] * (u - u_so[i])));
  const double tau_s = (1.0 - H_uw) * tau_s1[i] + H_uw * tau_s2[i];
  const double tau_o = (1.0 - H_uo) * tau_o1[i] + H_uo * tau_o2[i];
  const double v_inf = (1.0 - H_umv);
  const double w_inf =
      (1.0 - H_uo) * (1.0 - u / tau_winf[i]) + H_uo * ws_inf[i];

  const double I_fi = -v * H_uv * (u - theta_v[i]) * (u_u[i] - u) / tau_fi[i];
  const double I_so = (u - u_o[i]) * (1.0 - H_uw) / tau_o + H_uw / tau_so;
  const double I_si = -H_uw * w * s / tau_si[i];

  // Compute RHS of state variable equations
  f(0) = -(I_fi + I_so + I_si + I_stim) + I_sac;
  f(1) = (1.0 - H_uv) * (v_inf - v) / taum_v - H_uv * v / taup_v[i];
  f(2) = (1.0 - H_uw) * (w_inf - w) / taum_w - H_uw * w / taup_w[i];
  f(3) = (0.5 * (1.0 + tanh(k_s[i] * (u - u_s[i]))) - s) / tau_s;
}

void BuenoOrovio::getj(const unsigned int zone_id, const int nX, const int nG,
                       const Vector<double> &X, const Vector<double> &Xg,
                       Array<double> &Jac, const double Ksac) const {
  // Create local copies of the state variables
  const double u = X(0);
  const double v = X(1);
  const double w = X(2);
  const double s = X(3);

  //  Define step functions
  const double H_uv = step(u - theta_v[zone_id]);
  const double H_uw = step(u - theta_w[zone_id]);
  const double H_umv = step(u - thetam_v[zone_id]);
  const double H_uo = step(u - theta_o[zone_id]);

  // Define delta functions
  const double D_uw = delta(u - theta_w[zone_id]);
  const double D_uv = delta(u - theta_v[zone_id]);

  // Define additional constants
  const double taum_v =
      (1.0 - H_umv) * taum_v1[zone_id] + H_umv * taum_v2[zone_id];
  const double taum_w =
      taum_w1[zone_id] + 0.50 * (taum_w2[zone_id] - taum_w1[zone_id]) *
                             (1.0 + tanh(km_w[zone_id] * (u - um_w[zone_id])));
  const double tau_so =
      tau_so1[zone_id] + 0.50 * (tau_so2[zone_id] - tau_so1[zone_id]) *
                             (1.0 + tanh(k_so[zone_id] * (u - u_so[zone_id])));
  const double tau_s = (1.0 - H_uw) * tau_s1[zone_id] + H_uw * tau_s2[zone_id];
  const double tau_o = (1.0 - H_uo) * tau_o1[zone_id] + H_uo * tau_o2[zone_id];
  const double v_inf = (1.0 - H_umv);
  const double w_inf =
      (1.0 - H_uo) * (1.0 - u / tau_winf[zone_id]) + H_uo * ws_inf[zone_id];

  // Define Jacobian
  // @todo[michelebucelli] Ksac might need to be added here.
  Jac = 0.0;

  {
    const double n1 = v * H_uv * (u_u[zone_id] + theta_v[zone_id] - 2.0 * u) /
                      tau_fi[zone_id];
    const double n2 = -(1.0 - H_uw) / tau_fi[zone_id];
    const double n3 =
        (-1.0 / tau_so + (theta_w[zone_id] - u_o[zone_id]) / tau_o +
         w * s / tau_si[zone_id]) *
        D_uw;

    Jac(0, 0) = n1 + n2 + n3;
  }

  Jac(0, 1) =
      H_uv * (u - theta_v[zone_id]) * (u_u[zone_id] - u) / tau_fi[zone_id];

  {
    const double n1 = H_uw / tau_si[zone_id];
    Jac(0, 2) = n1 * s;
    Jac(0, 3) = n1 * w;
  }

  {
    const double n1 = -1.0 / taum_v;
    const double n2 = -1.0 / taup_v[zone_id];
    Jac(1, 0) = ((v_inf - v) * n1 + v * n2) * D_uv;
    Jac(1, 1) = (1.0 - H_uv) * n1 + H_uv * n2;
  }

  {
    const double n1 = -1.0 / taum_w;
    const double n2 = -1.0 / taup_w[zone_id];
    Jac(2, 0) = ((w_inf - w) * n1 + w * n2) * D_uw;
    Jac(2, 2) = (1.0 - H_uw) * n1 + H_uw * n2;
  }

  {
    const double n1 = cosh(k_s[zone_id] * (u - u_s[zone_id]));
    const double n2 = 1.0 / (n1 * n1);
    const double n3 = 1.0 / tau_s;
    Jac(3, 0) = 0.50 * k_s[zone_id] * n2 * n3;
    Jac(3, 3) = -n3;
  }
}

REGISTER_IONIC_MODEL("BO", BuenoOrovio);