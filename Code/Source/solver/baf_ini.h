// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef BAF_INI_H 
#define BAF_INI_H 

#include "ComMod.h"
#include "Simulation.h"

namespace baf_ini_ns {

void baf_ini(Simulation* simulation, const Array<double>& Ao, Array<double>& Do, Array<double>& Yo);

void bc_ini(const ComMod& com_mod, const CmMod& cm_mod, bcType& lBc, faceType& lFa, const Array<double>& Do);

void face_ini(Simulation* simulation, mshType& lm, faceType& la, Array<double>& Do);

void fsi_ls_ini(ComMod& com_mod, const CmMod& cm_mod, bcType& lBc, const faceType& lFa, int& lsPtr, const Array<double>& Do);

void set_shl_xien(Simulation* simulation, mshType& mesh);

void shl_bc_ini(const ComMod& com_mod, const CmMod& cm_mod, bcType& lBc, faceType& lFa, mshType& lM, const Array<double>& Do);

void shl_ini(const ComMod& com_mod, const CmMod& cm_mod, mshType& lM); 

};

#endif

