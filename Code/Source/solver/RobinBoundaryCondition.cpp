// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "RobinBoundaryCondition.h"
#include <iostream>

#define n_debug_robin_bc

RobinBoundaryCondition::RobinBoundaryCondition(double uniform_stiffness, double uniform_damping, bool normal_only, const faceType& face, SimulationLogger& logger)
    : BoundaryCondition({{"Stiffness", uniform_stiffness}, {"Damping", uniform_damping}}, StringBoolMap{{"normal_direction_only", normal_only}}, face, logger)
{
    // Warning if both stiffness and damping are zero
    if (uniform_stiffness == 0.0 && uniform_damping == 0.0) {
        logger_ -> log_message("WARNING [RobinBoundaryCondition] Both stiffness and damping values set to zero. "
                   "This will result in effectively no boundary condition being applied.");
    }
}

