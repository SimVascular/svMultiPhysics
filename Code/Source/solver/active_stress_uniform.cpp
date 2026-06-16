// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "active_stress_uniform.h"

void UniformActiveStress::read_parameters(
    const ActiveStressModelParameters &params) {
  value = params.get_scalar("Value");
}

void UniformActiveStress::distribute_parameters(const CmMod &cm_mod,
                                                const cmType &cm) {
  cm.bcast(cm_mod, &value);
}

double UniformActiveStress::operator()(const int idx) const { return value; }

REGISTER_ACTIVE_STRESS_MODEL("Uniform", UniformActiveStress);