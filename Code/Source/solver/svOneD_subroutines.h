// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SV1D_SUBROUTINES_H
#define SV1D_SUBROUTINES_H

#include "Simulation.h"
#include "consts.h"
#include "svOneD_interface/OneDSolverInterface.h"

namespace svOneD {

/// @brief Initialize the 1D solver and populate the initial cplBC state.
/// Called once from baf_ini() after the BC data structures are set up.
void init_svOneD(ComMod& com_mod, const CmMod& cm_mod);

/// @brief Advance the 1D solver by one time step and update the coupled BC value.
///
/// @param BCFlag  'D' - derivative / perturbation step (state is NOT committed).
///               'L' - last Newton iteration (state IS committed, time advances).
void calc_svOneD(ComMod& com_mod, const CmMod& cm_mod, char BCFlag);

}  // namespace svOneD

#endif  // SV1D_SUBROUTINES_H
