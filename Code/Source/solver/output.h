// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef OUTPUT__H 
#define OUTPUT__H 

#include "Simulation.h"
#include "SolutionStates.h"

#include<fstream>
#include<iostream>

namespace output {

/**
 * @brief Write the header of the history table to the standard output and to
 * the history file.
 *
 * The elapsed time reported by the rows of the table is measured from the
 * moment this function is called.
 *
 * @param[in] simulation Simulation whose logger the header is written to.
 * @param[in,out] timeP Timers of the history table: the elapsed time of the
 *   runs preceding this one, replaced by the instant the elapsed time is
 *   measured from, the elapsed time of the previous row and the elapsed time of
 *   the current row.
 */
void output_header(const Simulation *simulation, std::array<double, 3> &timeP);

/**
 * @brief Write one row of the history table to the standard output and to the
 * history file.
 *
 * The row reports the residuals and the linear solver statistics of the
 * nonlinear iteration that has just been completed.
 *
 * @param[in] simulation Simulation whose logger the row is written to.
 * @param[in,out] timeP Timers of the history table: the instant the elapsed
 *   time is measured from, the elapsed time of the previous row, which is
 *   replaced by the one of this row, and the elapsed time of this row.
 * @param[in] save_results True if the results of this time step are written to
 *   a VTU file, which flags the row with an 's'.
 * @param[in] iEq Index of the equation the row reports on.
 */
void output_result(const Simulation *simulation, std::array<double, 3> &timeP,
                   const bool save_results, const int iEq);

void read_restart_header(ComMod& com_mod, std::array<int,7>& tStamp, double& timeP, std::ifstream& restart_file);

void write_restart(Simulation* simulation, std::array<double,3>& timeP, const SolutionStates& solutions);

void write_restart_header(ComMod& com_mod, std::array<double,3>& timeP, std::ofstream& restart_file);

void write_results(ComMod& com_mod, const std::array<double,3>& timeP, const std::string& fName, const bool sstEq, const SolutionStates& solutions);

void write_ris_data(ComMod& com_mod, std::ofstream& restart_file);
void write_uris_data(ComMod& com_mod, std::ofstream& restart_file);

};

#endif

