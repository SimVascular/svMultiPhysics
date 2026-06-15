// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "active_stress.h"

void ActiveStress::init(const unsigned int tnNo) {
  states.resize(n_states, tnNo);

  if (n_states > 0) {
    Vector<double> state_loc(n_states);
    init_local(state_loc);

    for (unsigned int i = 0; i < tnNo; ++i)
      for (unsigned int j = 0; j < n_states; ++j)
        states(j, i) = state_loc(j);
  }
}

void ActiveStress::advance_time_step(const double t, const double dt) {
  for (unsigned int i = 0; i < states.ncols(); ++i) {
    Vector<double> state_loc = states.col(i);
    advance_time_step_local(t, dt,
                            /* calcium = */ 0.0,
                            /* fiber_stretch = */ 0.0,
                            /* fiber_stretch_rate = */ 0.0, state_loc);
    states.set_col(i, state_loc);
  }
}