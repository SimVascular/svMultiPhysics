// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "Simulation.h"

#ifndef INITIALIZE_H
#define INITIALIZE_H

void finalize(Simulation* simulation);

void init_from_bin(Simulation* simulation, const std::string& fName, std::array<double,3>& timeP);

void init_from_vtu(Simulation* simulation, const std::string& fName, std::array<double,3>& timeP);

void initialize(Simulation* simulation, Vector<double>& timeP);

void init_ris_data(ComMod& com_mod, std::ifstream& restart_file);
void init_uris_data(ComMod& com_mod, std::ifstream& restart_file);

void zero_init(Simulation* simulation);

#endif

