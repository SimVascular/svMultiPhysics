// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ActiveStressLandNiederer.h"

void ActiveStressLandNiederer::read_model_specific_parameters(
    const ActiveStressModelParameters &params) {
  ActiveStressODE::read_model_specific_parameters(params);

  calcium_scaling_factor = params.get_scalar("calcium_scaling_factor");
  CaRef = params.get_scalar("CaRef");
  eta_Tm = params.get_scalar("eta_Tm");
  k_uw = params.get_scalar("k_uw");
  k_ws = params.get_scalar("k_ws");
  Tref = params.get_scalar("Tref");
  k_TRPN = params.get_scalar("k_TRPN");
  eta_TRPN = params.get_scalar("eta_TRPN");
  k_u = params.get_scalar("k_u");
  TRPN50 = params.get_scalar("TRPN50");
  rw = params.get_scalar("rw");
  rs = params.get_scalar("rs");
  gamma_s = params.get_scalar("gamma_s");
  gamma_w = params.get_scalar("gamma_w");
  phi = params.get_scalar("phi");
  Aeff = params.get_scalar("Aeff");
  beta0 = params.get_scalar("beta0");
  beta1 = params.get_scalar("beta1");

  disable_force_strain_rate_feedback_ =
      params.get_bool("Disable_force_strain_rate_feedback");

  // With the feedback disabled the stretch rate computation can be skipped.
  needs_fiber_stretch_rate_ = !disable_force_strain_rate_feedback_;
}

void ActiveStressLandNiederer::distribute_model_specific_parameters(const CmMod &cm_mod,
                                                        const cmType &cm) {
  ActiveStressODE::distribute_model_specific_parameters(cm_mod, cm);

  cm.bcast(cm_mod, &calcium_scaling_factor);
  cm.bcast(cm_mod, &CaRef);
  cm.bcast(cm_mod, &eta_Tm);
  cm.bcast(cm_mod, &k_uw);
  cm.bcast(cm_mod, &k_ws);
  cm.bcast(cm_mod, &Tref);
  cm.bcast(cm_mod, &k_TRPN);
  cm.bcast(cm_mod, &eta_TRPN);
  cm.bcast(cm_mod, &k_u);
  cm.bcast(cm_mod, &TRPN50);
  cm.bcast(cm_mod, &rw);
  cm.bcast(cm_mod, &rs);
  cm.bcast(cm_mod, &gamma_s);
  cm.bcast(cm_mod, &gamma_w);
  cm.bcast(cm_mod, &phi);
  cm.bcast(cm_mod, &Aeff);
  cm.bcast(cm_mod, &beta0);
  cm.bcast(cm_mod, &beta1);
  cm.bcast(cm_mod, &disable_force_strain_rate_feedback_);
  cm.bcast(cm_mod, &needs_fiber_stretch_rate_);
}

void ActiveStressLandNiederer::init_local(Vector<double> &state) const { 
    state[0] = 1.0; 
    state[1] = 0.0; 
    state[2] = 0.0; 
    state[3] = 0.0; 
    state[4] = 0.0; 
    state[5] = 0.0; 
    const double lambda = 1.0; 
    state[6] = std::max(0.0, 1.0 + beta0 * (lambda + std::min(0.87, lambda) - 1.87)); 
}

void ActiveStressLandNiederer::advance_time_step_local(const double t, const double dt,
                                           const double calcium,
                                           const double fiber_stretch,
                                           const double fiber_stretch_rate,
                                           Vector<double> &state) const {

  Vector<double> ode_state(6); 
  for (int i = 0; i < 6; i++) {
    ode_state[i] = state[i]; 
  }

  ActiveStressODE::advance_time_step_local(t, dt, calcium, fiber_stretch,
                                           fiber_stretch_rate, ode_state);

  for (unsigned int i = 0; i < 6; ++i) {
    state[i] = ode_state[i];
  }

  const double lambda = std::min(1.2, fiber_stretch); 
  state[6] = std::max(0.0, 1.0 + beta0 * (lambda + std::min(0.87, lambda) - 1.87));
}

Vector<double> ActiveStressLandNiederer::getf(const double t, const Vector<double> &state,
                                  const double calcium,
                                  const double fiber_stretch,
                                  const double fiber_stretch_rate) const {
  Vector<double> f(6);

  // State Variables
  const double XB = state[0]; 
  const double XW = std::max(0.0, state[1]); 
  const double CaTRPN = std::max(0.0, state[2]); 
  const double XS = std::max(0.0, state[3]); 
  const double ZETAW = state[4]; 
  const double ZETAS = state[5]; 

  // Bound State Equation
  const double lambda = std::min(1.2, fiber_stretch); 
  const double XU = (1 - XB) - XS - XW; 
  const double k_b = k_u * std::pow(TRPN50, eta_Tm)/ (1 - rs - (1 - rs)*rw);
  f[0] = k_b * std::min(100.0, std::pow(CaTRPN, -(eta_Tm/2.0))) * XU  - k_u * std::pow(CaTRPN, (eta_Tm/2.0)) *  XB;

  // Pre-power Stroke Kinetics 
  const double k_wu = (k_uw * (1.0/rw - 1.0) - k_ws);
  const double gamma_wu = gamma_w * std::abs(ZETAW); 
  f[1] = k_uw*XU - k_wu*XW - k_ws*XW - gamma_wu*XW;

  // Calcium-Troponin Binding Kinetics  
  const double CaT50 = CaRef + beta1*std::min(0.2, lambda - 1.0); 
  f[2] = k_TRPN * (std::pow(((calcium*calcium_scaling_factor)/CaT50), eta_TRPN) * (1.0 - CaTRPN) - CaTRPN); 

  // Post-power Stroke Kinetics                                 
  const double k_su = k_ws * rw * (1.0/rs - 1.0); 
  const double term1 = (ZETAS > 0.0) ? ZETAS : 0.0;
  const double term2 = (ZETAS < -1.0) ? (-ZETAS - 1.0) : 0.0;
  const double gamma_su = gamma_s * std::max(term1, term2);
  f[3] = k_ws*XW - k_su*XS - gamma_su*XS; 

  // Distortion-Decay Kinetics
  const double stretch_rate =
      disable_force_strain_rate_feedback_ ? 0.0 : fiber_stretch_rate;
  const double Aw = Aeff * rs/((1.0 - rs) * rw + rs);
  const double cw = phi * k_uw * (1.0 - rs)*(1.0 - rw) / ((1.0 - rs)*rw);
  f[4] = Aw * stretch_rate - cw*ZETAW;

  const double As = Aw;
  const double cs = phi * k_ws * (1.0 - rs)*rw/rs;
  f[5] =  As * stretch_rate - cs*ZETAS;

  return f;
}

double ActiveStressLandNiederer::compute_active_tension_local(const Vector<double> &state,
                                                   const double fiber_stretch) const {
  const double XW = std::max(0.0, state[1]);
  const double XS = std::max(0.0, state[3]); 
  const double ZETAW = state[4];  
  const double ZETAS = state[5];
  const double LFac = state[6];  
  
  return LFac * (Tref/rs) * ((ZETAS+1.0) * XS + (ZETAW) * XW);
}

REGISTER_ACTIVE_STRESS_MODEL("LandNiederer", ActiveStressLandNiederer);