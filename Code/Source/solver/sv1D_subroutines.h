// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SV1D_SUBROUTINES_H
#define SV1D_SUBROUTINES_H

#include "Simulation.h"
#include "consts.h"
#include "sv1D_interface/OneDSolverInterface.h"

namespace sv1D {

/// @brief Initialize the 1D solver and populate the initial cplBC state.
/// Called once from baf_ini() after the BC data structures are set up.
void init_sv1D(ComMod& com_mod, const CmMod& cm_mod);

/// @brief Advance the 1D solver by one time step and update the coupled BC value.
///
/// @param BCFlag  'D' - derivative / perturbation step (state is NOT committed).
///               'L' - last Newton iteration (state IS committed, time advances).
void calc_sv1D(ComMod& com_mod, const CmMod& cm_mod, char BCFlag);

}  // namespace sv1D

#endif  // SV1D_SUBROUTINES_H
