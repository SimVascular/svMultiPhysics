// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef POST_H 
#define POST_H 

#include "Simulation.h"
#include "SolutionStates.h"
#include "consts.h"

namespace post {

void all_post(Simulation* simulation, Array<double>& res, const SolutionStates& solutions,
    consts::OutputNameType outGrp, const int iEq);

void bpost(Simulation* simulation, const mshType& lM, Array<double>& res, const SolutionStates& solutions,
    consts::OutputNameType outGrp);

void div_post(Simulation* simulation, const mshType& lM, Array<double>& res, const SolutionStates& solutions, const int iEq);

void fib_algn_post(Simulation* simulation, const mshType& lM, Array<double>& res, const SolutionStates& solutions, const int iEq);

void fib_dir_post(Simulation* simulation, const mshType& lM, const int nFn, Array<double>& res, const SolutionStates& solutions, const int iEq);

void fib_strech(Simulation* simulation, const int iEq, const mshType& lM, const SolutionStates& solutions, Vector<double>& res);

void post(Simulation* simulation, const mshType& lM, Array<double>& res, const SolutionStates& solutions,
    consts::OutputNameType outGrp, const int iEq);

void ppbin2vtk(Simulation* simulation);

void shl_post(Simulation* simulation, const mshType& lM, const int m, Array<double>& res, 
    Vector<double>& resE, const SolutionStates& solutions, const int iEq, consts::OutputNameType outGrp);

void tpost(Simulation* simulation, const mshType& lM, const int m, Array<double>& res, Vector<double>& resE,
    const SolutionStates& solutions, const int iEq, consts::OutputNameType outGrp);

};

#endif

