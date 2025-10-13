// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef EQ_ASSEM_H 
#define EQ_ASSEM_H 

#include "ComMod.h"
#include "Simulation.h"

namespace eq_assem {

void b_assem_neu_bc(ComMod& com_mod, const faceType& lFa, const Vector<double>& hg, const Array<double>& Yg);

void b_neu_folw_p(ComMod& com_mod, const bcType& lBc, const faceType& lFa, const Vector<double>& hg, const Array<double>& Dg);

void fsi_ls_upd(ComMod& com_mod, const bcType& lBc, const faceType& lFa);

void global_eq_assem(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Yg, const Array<double>& Dg);

};

#endif

