// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVZEROD_H
#define SVZEROD_H 

#include "Simulation.h"
#include "consts.h"
#include "svZeroD_interface/LPNSolverInterface.h"
#include <vector>

#include <string>

namespace svZeroD {

void get_coupled_QP(ComMod& com_mod, double QCoupled[], double QnCoupled[], double PCoupled[], double PnCoupled[]);

void print_svZeroD(int* nSrfs, const std::vector<int>& surfID, double Q[], double P[]);

/**
 * @brief Set up the svZeroD model and its coupled boundary conditions.
 *
 * Builds the list of svZeroD-coupled boundaries, creates the svZeroD model
 * from the interface data in \c com_mod.cplBC, applies the initial flows and
 * pressures, and writes the header of the svZeroD state output file.
 *
 * @param[in,out] com_mod Simulation data; the coupled boundary conditions and
 *   \c cplBC bookkeeping are initialized here.
 * @param[in] cm_mod MPI communicator data.
 * @param[in] appPath Directory the simulation results are written to. The
 *   svZeroD output files (svZeroD_data, Q_svZeroD, P_svZeroD) are written
 *   there as well.
 */
void init_svZeroD(ComMod& com_mod, const CmMod& cm_mod, const std::string& appPath);

void calc_svZeroD(ComMod& com_mod, const CmMod& cm_mod, char BCFlag);

};

#endif
