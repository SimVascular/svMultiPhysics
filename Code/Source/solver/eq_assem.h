// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef EQ_ASSEM_H 
#define EQ_ASSEM_H 

#include "ComMod.h"
#include "SolutionStates.h"
#include "Simulation.h"

namespace eq_assem {

void b_assem_neu_bc(ComMod& com_mod, const faceType& lFa, const Vector<double>& hg, const SolutionStates& solutions);

void b_neu_folw_p(ComMod& com_mod, const bcType& lBc, const faceType& lFa, const Vector<double>& hg, const SolutionStates& solutions);

/// @brief Update the vectors that the linear solver uses for the tangent
/// contribution of a resistance or coupled boundary condition, so that they
/// follow the deformed geometry under a follower pressure load.
///
/// Computes sV_a = integral over the face of N_a * n in the configuration at
/// the end of the time step, where N_a is the shape function of node a and n
/// is the normal (see Brown et al. 2024, Eq. 56). The result is stored in the
/// linear solver face com_mod.lhs.face[lBc.lsPtr], as both unit_pressure_load
/// and flowrate_derivative. For a coupled boundary condition with a cap, the
/// contribution of the cap is appended to flowrate_derivative. The linear
/// solver reads the two vectors in add_bc_mul().
///
/// @param[in,out] com_mod The simulation data. Only
///   com_mod.lhs.face[lBc.lsPtr] is modified.
/// @param[in] lBc The boundary condition.
/// @param[in] lFa The face of the boundary condition.
/// @param[in] solutions The solution states, from which the displacement is
///   read.
void fsi_ls_upd(ComMod& com_mod, const bcType& lBc, const faceType& lFa, const SolutionStates& solutions);

void global_eq_assem(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const SolutionStates& solutions);

};

#endif

