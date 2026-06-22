// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "active_stress_ode.h"

void ActiveStressODE::advance_time_step_local(const double t, const double dt,
                                              const double calcium,
                                              const double fiber_stretch,
                                              const double fiber_stretch_rate,
                                              Vector<double> &state) const {
  // Forward Euler time stepping.
  Vector<double> f = getf(t, state, calcium, fiber_stretch, fiber_stretch_rate);

  state.add(dt, f);
}