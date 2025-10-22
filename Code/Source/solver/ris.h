// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef RIS_H 
#define RIS_H

#include "ComMod.h"

namespace ris {

void ris_meanq(ComMod& com_mod, CmMod& cm_mod);
void ris_resbc(ComMod& com_mod, const Array<double>& Yg, const Array<double>& Dg);
void setbc_ris(ComMod& com_mod, const bcType& lBc, const mshType& lM, const faceType& lFa, 
    const Array<double>& Yg, const Array<double>& Dg);

void ris_updater(ComMod& com_mod, CmMod& cm_mod);
void ris_status(ComMod& com_mod, CmMod& cm_mod);

void doassem_ris(ComMod& com_mod, const int d, const Vector<int>& eqN, 
    const Array3<double>& lK, const Array<double>& lR); 

void doassem_velris(ComMod& com_mod, const int d, const Array<int>& eqN, 
    const Array3<double>& lK, const Array<double>& lR);

void clean_r_ris(ComMod& com_mod);
void setbcdir_ris(ComMod& com_mod, Array<double>& lA, Array<double>& lY, Array<double>& lD);

// TODO: RIS 0D code
void ris0d_bc(ComMod& com_mod, CmMod& cm_mod, const Array<double>& Yg, const Array<double>& Dg);
void ris0d_status(ComMod& com_mod, CmMod& cm_mod); //, const Array<double>& Yg, const Array<double>& Dg);

};

#endif

