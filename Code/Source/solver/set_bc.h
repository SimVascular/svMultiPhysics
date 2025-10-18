// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SET_BC_H 
#define SET_BC_H 

#include "Simulation.h"
#include "consts.h"
#include "RobinBoundaryCondition.h"

#include <string>

namespace set_bc {

void calc_der_cpl_bc(ComMod& com_mod, const CmMod& cm_mod);

void cplBC_Integ_X(ComMod& com_mod, const CmMod& cm_mod, const bool RCRflag);

void genBC_Integ_X(ComMod& com_mod, const CmMod& cm_mod, const std::string& genFlag);

void rcr_init(ComMod& com_mod, const CmMod& cm_mod);

void RCR_Integ_X(ComMod& com_mod, const CmMod& cm_mod, int istat);

void set_bc_cmm(ComMod& com_mod, const CmMod& cm_mod, const Array<double>& Ag, const Array<double>& Dg);
void set_bc_cmm_l(ComMod& com_mod, const CmMod& cm_mod, const faceType& lFa, const Array<double>& Ag, const Array<double>& Dg );

void set_bc_cpl(ComMod& com_mod, CmMod& cm_mod);

void set_bc_dir(ComMod& com_mod, Array<double>& lA, Array<double>& lY, Array<double>& lD);
void set_bc_dir_l(ComMod& com_mod, const bcType& lBc, const faceType& lFa, Array<double>& lA, Array<double>& lY, int lDof);
void set_bc_dir_w(ComMod& com_mod, const Array<double>& Yg, const Array<double>& Dg);
void set_bc_dir_wl(ComMod& com_mod, const bcType& lBc, const mshType& lM, const faceType& lFa, const Array<double>& Yg, const Array<double>& Dg);

void set_bc_neu(ComMod& com_mod, const CmMod& cm_mod, const Array<double>& Yg, const Array<double>& Dg);
void set_bc_neu_l(ComMod& com_mod, const CmMod& cm_mod, const bcType& lBc, const faceType& lFa, const Array<double>& Yg, const Array<double>& Dg);

void set_bc_rbnl(ComMod& com_mod, const faceType& lFa, const RobinBoundaryCondition& robin_bc,
  const Array<double>& Yg, const Array<double>& Dg);

void set_bc_trac_l(ComMod& com_mod, const CmMod& cm_mod, const bcType& lBc, const faceType& lFa);

void set_bc_undef_neu(ComMod& com_mod);

void set_bc_undef_neu_l(ComMod& com_mod, const bcType& lBc, const faceType& lFa);

};

#endif

