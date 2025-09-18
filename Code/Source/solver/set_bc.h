/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef SET_BC_H 
#define SET_BC_H 

#include "Simulation.h"
#include "consts.h"
#include "VariableRobinBCData.h"

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

void set_bc_rbnl(ComMod& com_mod, const faceType& lFa, const VariableRobinBCData& variable_robin_data, const bool isN, 
  const Array<double>& Yg, const Array<double>& Dg);

void set_bc_trac_l(ComMod& com_mod, const CmMod& cm_mod, const bcType& lBc, const faceType& lFa);

void set_bc_undef_neu(ComMod& com_mod);

void set_bc_undef_neu_l(ComMod& com_mod, const bcType& lBc, const faceType& lFa);

};

#endif

