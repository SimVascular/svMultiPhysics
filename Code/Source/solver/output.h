// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef OUTPUT__H 
#define OUTPUT__H 

#include "Simulation.h"

#include<fstream>
#include<iostream>

namespace output {

void output_result(Simulation* simulation,  std::array<double,3>& timeP, const int co, const int iEq);

void read_restart_header(ComMod& com_mod, std::array<int,7>& tStamp, double& timeP, std::ifstream& restart_file);

void write_restart(Simulation* simulation, std::array<double,3>& timeP);

void write_restart_header(ComMod& com_mod, std::array<double,3>& timeP, std::ofstream& restart_file);

void write_results(ComMod& com_mod, const std::array<double,3>& timeP, const std::string& fName, const bool sstEq);

};

#endif

