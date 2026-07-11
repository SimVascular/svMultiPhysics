// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "active_stress_ode.h"

void ActiveStressODE::read_model_specific_parameters(
    const ActiveStressModelParameters &params) {
  const std::string solver_str = params.get_string("ODE_solver");

  if (solver_str == "FE") {
    ode_solver = ODESolver::ForwardEuler;
  } else {
    svmp::raise<svmp::ParseException>("Unknown ODE solver " + solver_str +
                                      " for active stress models.");
  }
}

void ActiveStressODE::distribute_model_specific_parameters(const CmMod &cm_mod,
                                                           const cmType &cm) {
  cm.bcast_enum(cm_mod, &ode_solver);
}

void ActiveStressODE::advance_time_step_local(const double t, const double dt,
                                              const double calcium,
                                              const double fiber_stretch,
                                              const double fiber_stretch_rate,
                                              Vector<double> &state) const {
  // Forward Euler time stepping.
  Vector<double> f = getf(t, state, calcium, fiber_stretch, fiber_stretch_rate);

  state.add(dt, f);
}