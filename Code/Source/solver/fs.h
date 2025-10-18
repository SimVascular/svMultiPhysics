// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// These routines initialize function spaces.

#ifndef FS_H 
#define FS_H 

#include "Simulation.h"
#include "ComMod.h"

namespace fs {

void alloc_fs(fsType& fs, const int nsd, const int insd);

void get_thood_fs(ComMod& com_mod, std::array<fsType,2>& fs, const mshType& lM, const bool lStab, const int iOpt);

void init_fs(fsType& fs,  const int nsd, const int insd);

void init_fs_face(const ComMod& com_mod, mshType& lM, faceType& lFa);

void init_fs_msh(const ComMod& com_mod, mshType& mesh);

void set_thood_fs(fsType& fs, consts::ElementType eType);

void thood_val_rc(ComMod& com_mod);

};

#endif

