// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef FSI_H 
#define FSI_H 

#include "ComMod.h"
#include "SolutionStates.h"
#include "consts.h"

namespace fsi {

void construct_fsi(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const SolutionStates& solutions);

};

#endif

