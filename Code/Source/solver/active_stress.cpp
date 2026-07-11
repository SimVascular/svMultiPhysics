// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "active_stress.h"

bool supports_active_stress(const consts::EquationType eq_type) {
  return eq_type == consts::EquationType::phys_struct ||
         eq_type == consts::EquationType::phys_ustruct ||
         eq_type == consts::EquationType::phys_FSI;
}

void ActiveStress::read_parameters(const ActiveStressParameters &params) {
  eta_f = params.get_eta_f();
  eta_s = params.get_eta_s();
  eta_n = params.get_eta_n();

  read_model_specific_parameters(
      params.get_parameters(params.get_model_name()));
}

void ActiveStress::distribute_parameters(const CmMod &cm_mod,
                                         const cmType &cm) {
  cm.bcast(cm_mod, &eta_f);
  cm.bcast(cm_mod, &eta_s);
  cm.bcast(cm_mod, &eta_n);

  distribute_model_specific_parameters(cm_mod, cm);
}

void ActiveStress::init(const unsigned int tnNo) {
  states.resize(n_states, tnNo);

  if (n_states > 0) {
    Vector<double> state_loc(n_states);
    init_local(state_loc);

    for (unsigned int i = 0; i < tnNo; ++i)
      for (unsigned int j = 0; j < n_states; ++j)
        states(j, i) = state_loc(j);
  }

  active_tension.resize(tnNo);
}

void ActiveStress::advance_time_step(const double t, const double dt,
                                     const Vector<double> &calcium,
                                     const Vector<double> &fiber_stretch,
                                     const Vector<double> &fiber_stretch_rate) {
  time = t;

  for (unsigned int i = 0; i < states.ncols(); ++i) {
    Vector<double> state_loc = states.col(i);
    advance_time_step_local(t, dt, calcium[i], fiber_stretch[i],
                            fiber_stretch_rate[i], state_loc);
    states.set_col(i, state_loc);

    active_tension[i] = compute_active_tension_local(state_loc);
  }
}