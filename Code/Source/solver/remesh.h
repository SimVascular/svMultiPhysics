// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef REMESH_H 
#define REMESH_H 

#include "Simulation.h"

namespace remesh {

void remesh_restart(Simulation* simulation);

void set_face_ebc(ComMod& com_mod, CmMod& cm_mod, faceType& lFa, mshType& lM);

};

#endif

