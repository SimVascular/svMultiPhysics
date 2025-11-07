// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "Simulation.h"

#ifndef INITIALIZE_H
#define INITIALIZE_H

void finalize(Simulation* simulation);

void init_from_bin(Simulation* simulation, const std::string& fName, std::array<double,3>& timeP,
                   Array<double>& Ao, Array<double>& Do, Array<double>& Yo);

void init_from_vtu(Simulation* simulation, const std::string& fName, std::array<double,3>& timeP,
                   Array<double>& Ao, Array<double>& Do, Array<double>& Yo);

void initialize(Simulation* simulation, Vector<double>& timeP);

void zero_init(Simulation* simulation, Array<double>& Ao, Array<double>& Do, Array<double>& Yo);

#endif

