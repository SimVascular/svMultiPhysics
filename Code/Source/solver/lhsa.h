// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef LHSA_H 
#define LHSA_H 

#include "ComMod.h"
#include "Simulation.h"

namespace lhsa_ns {

  void add_col(const int tnNo, const int rowN, const int colN, int& mnnzeic, Array<int>& uInd);

  void do_assem(ComMod& com_mod, const int d, const Vector<int>& eqN, const Array3<double>& lK, const Array<double>& lR);

  void lhsa(Simulation* simulation, int& nnz);

  void resiz(const int tnNo, int& mnnzeic, Array<int>& uInd);

};

#endif

