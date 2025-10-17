// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "Simulation.h"

#ifndef PIC_H
#define PIC_H

namespace pic {

void picc(Simulation* simulation);

void pic_eth(Simulation* simulation);

void pici(Simulation* simulation, Array<double>& Ag, Array<double>& Yg, Array<double>& Dg);

void picp(Simulation* simulation);

};

#endif

