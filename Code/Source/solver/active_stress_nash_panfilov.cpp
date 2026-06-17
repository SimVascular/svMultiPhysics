// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "active_stress_nash_panfilov.h"

void NashPanfilov::read_parameters(const ActiveStressModelParameters &params) {
  rate = params.get_scalar("Rate");
}

void NashPanfilov::distribute_parameters(const CmMod &cm_mod,
                                         const cmType &cm) {
  cm.bcast(cm_mod, &rate);
}

void NashPanfilov::init_local(Vector<double> &state) const { state[0] = 0.0; }

Vector<double> NashPanfilov::getf(const double t, const Vector<double> &state,
                                  const double calcium,
                                  const double fiber_stretch,
                                  const double fiber_stretch_rate) const {
  Vector<double> f(1);

  // @todo[michelebucelli] Implement the actual model here.
  f[0] = rate;

  return f;
}

double
NashPanfilov::compute_active_tension_local(const Vector<double> &state) const {
  // @todo[michelebucelli] Implement the actual model here.
  return state[0];
}

REGISTER_ACTIVE_STRESS_MODEL("NashPanfilov", NashPanfilov);