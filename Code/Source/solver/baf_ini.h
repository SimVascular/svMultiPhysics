// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef BAF_INI_H 
#define BAF_INI_H 

#include "ComMod.h"
#include "Simulation.h"

namespace baf_ini_ns {

void baf_ini(Simulation* simulation);

void bc_ini(const ComMod& com_mod, const CmMod& cm_mod, bcType& lBc, faceType& lFa);

void face_ini(Simulation* simulation, mshType& lm, faceType& la);

void fsi_ls_ini(ComMod& com_mod, const CmMod& cm_mod, bcType& lBc, const faceType& lFa, int& lsPtr);

void set_shl_xien(Simulation* simulation, mshType& mesh);

void shl_bc_ini(const ComMod& com_mod, const CmMod& cm_mod, bcType& lBc, faceType& lFa, mshType& lM);

void shl_ini(const ComMod& com_mod, const CmMod& cm_mod, mshType& lM); 

};

#endif

