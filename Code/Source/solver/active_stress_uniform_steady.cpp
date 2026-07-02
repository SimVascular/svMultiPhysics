// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "active_stress_uniform_steady.h"

void UniformSteadyActiveStress::read_model_specific_parameters(
    const ActiveStressModelParameters &params) {
  value = params.get_scalar("Value");
}

void UniformSteadyActiveStress::distribute_model_specific_parameters(
    const CmMod &cm_mod, const cmType &cm) {
  cm.bcast(cm_mod, &value);
}

REGISTER_ACTIVE_STRESS_MODEL("UniformSteady", UniformSteadyActiveStress);