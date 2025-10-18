// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef MESH_H 
#define MESH_H 

#include "ComMod.h"

#include "consts.h"

namespace mesh {

void construct_mesh(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Dg);

};

#endif

