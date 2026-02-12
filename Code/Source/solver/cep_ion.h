// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef CEP_ION_H 
#define CEP_ION_H 

#include "Array.h"
#include "ComMod.h"
#include "Simulation.h"

#include "all_fun.h"
#include "consts.h"

#include <string>

namespace cep_ion {

void cep_init(Simulation* simulation);

void cep_init_l(cepModelType& cep, int nX, int nG, Vector<double>& X, Vector<double>& Xg);

void cep_integ(Simulation* simulation, const int iEq, const int iDof, const Array<double>& Dg);

void cep_integ_l(CepMod& cep_mod, cepModelType& cep, int nX, int nG, Vector<double>& X, Vector<double>& Xg,
    const double t1, double& yl, const double I4f, const double dt);

};

#endif

