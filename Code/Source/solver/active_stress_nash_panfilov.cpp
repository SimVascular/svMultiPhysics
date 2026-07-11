// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "active_stress_nash_panfilov.h"

void NashPanfilov::read_model_specific_parameters(
    const ActiveStressModelParameters &params) {
  ActiveStressODE::read_model_specific_parameters(params);

  epsilon_0 = params.get_scalar("epsilon_0");
  epsilon_i = params.get_scalar("epsilon_i");
  xi_T = params.get_scalar("xi_T");
  calcium_rest = params.get_scalar("calcium_rest");
  calcium_crit = params.get_scalar("calcium_crit");
  eta_T = params.get_scalar("eta_T");
}

void NashPanfilov::distribute_model_specific_parameters(const CmMod &cm_mod,
                                                        const cmType &cm) {
  ActiveStressODE::distribute_model_specific_parameters(cm_mod, cm);

  cm.bcast(cm_mod, &epsilon_0);
  cm.bcast(cm_mod, &epsilon_i);
  cm.bcast(cm_mod, &xi_T);
  cm.bcast(cm_mod, &calcium_rest);
  cm.bcast(cm_mod, &calcium_crit);
  cm.bcast(cm_mod, &eta_T);
}

void NashPanfilov::init_local(Vector<double> &state) const { state[0] = 0.0; }

Vector<double> NashPanfilov::getf(const double t, const Vector<double> &state,
                                  const double calcium,
                                  const double fiber_stretch,
                                  const double fiber_stretch_rate) const {
  Vector<double> f(1);

  const double epsilon =
      epsilon_0 + (epsilon_i - epsilon_0) *
                      std::exp(-std::exp(-xi_T * (calcium - calcium_crit)));

  f[0] = epsilon * (eta_T * (calcium - calcium_rest) - state[0]);

  return f;
}

double
NashPanfilov::compute_active_tension_local(const Vector<double> &state) const {
  return state[0];
}

REGISTER_ACTIVE_STRESS_MODEL("NashPanfilov", NashPanfilov);