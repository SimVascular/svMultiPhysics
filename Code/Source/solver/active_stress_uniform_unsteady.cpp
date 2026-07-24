// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "active_stress_uniform_unsteady.h"

#include <fstream>
#include <vector>

void UniformUnsteadyActiveStress::init(const unsigned int tnNo) {
  ActiveStress::init(tnNo);

  fourier_interpolation = FourierInterpolation::from_time_series_file(
      temporal_values_file_path, /* n_components = */ 1, ramp);
}

void UniformUnsteadyActiveStress::read_model_specific_parameters(
    const ActiveStressModelParameters &params) {
  ramp = params.get_bool("Ramp");
  temporal_values_file_path = params.get_string("Temporal_values_file_path");
}

void UniformUnsteadyActiveStress::distribute_model_specific_parameters(
    const CmMod &cm_mod, const cmType &cm) {
  cm.bcast(cm_mod, &ramp);
  cm.bcast(cm_mod, temporal_values_file_path);
}

double UniformUnsteadyActiveStress::compute_active_tension_local(
    const Vector<double> &state) const {
  return fourier_interpolation.value(time)[0];
}

REGISTER_ACTIVE_STRESS_MODEL("UniformUnsteady", UniformUnsteadyActiveStress);