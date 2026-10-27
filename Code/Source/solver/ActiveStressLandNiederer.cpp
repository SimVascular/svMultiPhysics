// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ActiveStressLandNiederer.h"

void ActiveStressLandNiederer::read_model_specific_parameters(
    const ActiveStressModelParameters &params) {
  ActiveStressODE::read_model_specific_parameters(params);

  // Read a scalar parameter and raise a ParseException if it violates the
  // given condition.
  auto get_checked = [&params](const std::string &label, auto is_valid,
                               const std::string &requirement) {
    const double value = params.get_scalar(label);
    if (!is_valid(value)) {
      std::ostringstream msg;
      msg << params.xml_element_name << ": parameter '" << label
          << "' must be " << requirement << ", got " << value << ".";
      svmp::raise<svmp::ParseException>(msg.str());
    }
    return value;
  };
  auto get_positive = [&](const std::string &label) {
    return get_checked(label, [](double v) { return v > 0.0; }, "> 0");
  };
  auto get_nonnegative = [&](const std::string &label) {
    return get_checked(label, [](double v) { return v >= 0.0; }, ">= 0");
  };
  auto get_nonpositive = [&](const std::string &label) {
    return get_checked(label, [](double v) { return v <= 0.0; }, "<= 0");
  };
  auto get_in_open_unit_interval = [&](const std::string &label) {
    return get_checked(label, [](double v) { return v > 0.0 && v < 1.0; },
                       "in (0, 1)");
  };

  CaRef = get_positive("CaRef");
  eta_Tm = params.get_scalar("eta_Tm");
  k_uw = get_positive("k_uw");
  k_ws = get_positive("k_ws");
  Tref = get_positive("Tref");
  k_TRPN = get_positive("k_TRPN");
  eta_TRPN = params.get_scalar("eta_TRPN");
  k_u = get_positive("k_u");
  TRPN50 = get_positive("TRPN50");
  rw = get_in_open_unit_interval("rw");
  rs = get_in_open_unit_interval("rs");
  gamma_s = params.get_scalar("gamma_s");
  gamma_w = params.get_scalar("gamma_w");
  phi = params.get_scalar("phi");
  Aeff = params.get_scalar("Aeff");
  beta0 = get_nonnegative("beta0");
  beta1 = get_nonpositive("beta1");
}

void ActiveStressLandNiederer::distribute_model_specific_parameters(const CmMod &cm_mod,
                                                        const cmType &cm) {
  ActiveStressODE::distribute_model_specific_parameters(cm_mod, cm);

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
}

void ActiveStressLandNiederer::init_local(Vector<double> &state) const { 
  state[0] = 1.0; // XB: Fraction of Blocked Myosin Binding Sites   
  state[1] = 0.0; // XW: Pre-Power Stroke State 
  state[2] = 0.0; // CaTRPN: Fraction of Calcium Bounded to Troponin
  state[3] = 0.0; // XS: Post-Power Stroke State 
  state[4] = 0.0; // ZETAW: Distortion-Decay of Pre-Power Stroke State W 
  state[5] = 0.0; // ZETAS: Distortion-Decay of Post-Power Stroke State S 
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
  const double CaT50 = CaRef + beta1 * (lambda - 1.0); 
  f[2] = k_TRPN * (std::pow((calcium/CaT50), eta_TRPN) * (1.0 - CaTRPN) - CaTRPN); 

  // Post-power Stroke Kinetics                                 
  const double k_su = k_ws * rw * (1.0/rs - 1.0); 
  const double term1 = (ZETAS > 0.0) ? ZETAS : 0.0;
  const double term2 = (ZETAS < -1.0) ? (-ZETAS - 1.0) : 0.0;
  const double gamma_su = gamma_s * std::max(term1, term2);
  f[3] = k_ws*XW - k_su*XS - gamma_su*XS; 

  // Distortion-Decay Kinetics
  const double A = Aeff * rs/((1.0 - rs) * rw + rs);
  const double cw = phi * k_uw * (1.0 - rs)*(1.0 - rw) / ((1.0 - rs)*rw);
  const double cs = phi * k_ws * (1.0 - rs)*rw/rs;
  f[4] =  A * fiber_stretch_rate - cw*ZETAW;
  f[5] =  A * fiber_stretch_rate - cs*ZETAS;

  return f;
}

double ActiveStressLandNiederer::compute_active_tension_local(const Vector<double> &state,
                                                   const double fiber_stretch) const {
  const double XW = std::max(0.0, state[1]);
  const double XS = std::max(0.0, state[3]); 
  const double ZETAW = state[4];  
  const double ZETAS = state[5];
  const double lambda = std::min(1.2, fiber_stretch); 
  const double LFac = length_dependence_factor(lambda);
  
  return LFac * (Tref/rs) * ((ZETAS+1.0) * XS + (ZETAW) * XW);
}

double ActiveStressLandNiederer::compute_active_tension_derivative_local(
    const Vector<double> &state, const double fiber_stretch) const {
  const double XW = std::max(0.0, state[1]);
  const double XS = std::max(0.0, state[3]);
  const double ZETAW = state[4];
  const double ZETAS = state[5];
  const double lambda = std::min(1.2, fiber_stretch);
  const double d_lambda = (fiber_stretch < 1.2) ? 1.0 : 0.0;
  const double d_LFac = length_dependence_factor_derivative(lambda) * d_lambda;

  return d_LFac * (Tref/rs) * ((ZETAS+1.0) * XS + (ZETAW) * XW);
}

double ActiveStressLandNiederer::length_dependence_factor(double lambda) const {
  return std::max(0.0, 1.0 + beta0 * (lambda + std::min(0.87, lambda) - 1.87));
}

double ActiveStressLandNiederer::length_dependence_factor_derivative(
    double lambda) const {
  if (length_dependence_factor(lambda) <= 0.0) {
    return 0.0;
  }
  return lambda < 0.87 ? 2.0 * beta0 : beta0;
}

REGISTER_ACTIVE_STRESS_MODEL("LandNiederer", ActiveStressLandNiederer);