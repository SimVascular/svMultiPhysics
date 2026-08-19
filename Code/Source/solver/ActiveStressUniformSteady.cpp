// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ActiveStressUniformSteady.h"

void ActiveStressUniformSteady::read_model_specific_parameters(
    const ActiveStressModelParameters &params) {
  value = params.get_scalar("Value");
}

void ActiveStressUniformSteady::distribute_model_specific_parameters(
    const CmMod &cm_mod, const cmType &cm) {
  cm.bcast(cm_mod, &value);
}

REGISTER_ACTIVE_STRESS_MODEL("UniformSteady", ActiveStressUniformSteady);