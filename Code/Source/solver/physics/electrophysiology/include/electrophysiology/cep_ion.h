// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef CEP_ION_H
#define CEP_ION_H

#include "Array.h"
#include "ComMod.h"
#include "Simulation.h"
#include "SolutionStates.h"

#include "all_fun.h"
#include "consts.h"

#include <string>

namespace cep_ion {

/// @brief Initialize the ionic model state variables and the action potential
/// of the initial solution.
///
/// The state variables of each CEP domain's ionic model are stored in
/// cep_mod.Xion. Its first row, the membrane potential, is also copied into the
/// row of the old velocity reserved for the CEP equation, so that the initial
/// solution agrees with the ionic model.
///
/// @param[in,out] simulation The simulation, whose cep_mod.Xion is filled with
///   the initial ionic model state variables.
/// @param[in,out] solutions The solution states, whose old velocity receives
///   the initial action potential.
void cep_init(Simulation *simulation, SolutionStates &solutions);

void cep_integ(Simulation *simulation, const int iEq, const int iDof,
               SolutionStates &solutions, const Vector<double> &I4f);

void cep_integ_l(CepMod &cep_mod, cepModelType &cep, Vector<double> &X,
                 Vector<double> &Xg, const double t1, const double I4f,
                 const double dt, const Vector<double> &x);
}; // namespace cep_ion

#endif
