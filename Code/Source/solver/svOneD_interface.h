// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SV1D_SUBROUTINES_H
#define SV1D_SUBROUTINES_H

#include "Simulation.h"
#include "consts.h"
#include "svOneD_interface/OneDSolverInterface.h"

/**
 * @namespace svOneD
 * @brief 3D-1D coupling subroutines.
 *
 * These routines interface the 3D finite-element solver svMultiPhysics with
 * the 1D blood-flow solver svOneDSolver via a dynamically loaded shared
 * library, such as libsvOneDSolver_interface.so or libsvOneDSolver_interface.dylib.
 *
 * @par Coupling overview
 * NEU coupling, where the 1D inlet is driven by the 3D outflow:
 * - 3D to 1D: flow rate Q, passed through params[3] and params[4].
 * - 1D to 3D: pressure P, returned through cpl_value.
 * - 3D boundary condition: Neumann pressure traction.
 *
 * DIR coupling, where the 1D outlet is driven by the 3D pressure:
 * - 3D to 1D: pressure P, passed through params[3] and params[4].
 * - 1D to 3D: flow rate Q, returned through cpl_value.
 * - 3D boundary condition: Dirichlet velocity profile.
 *
 * @par Parallelism model
 * Unlike the 0D solver, which is solved once on the master rank, each 1D model
 * is independent and has its own input file. Multiple 1D models are therefore
 * read, initialized, and solved in parallel.
 *
 * Initialization in init_svOneD() has two phases:
 * - Parallel initialization:
 *   - Collect all svOneD-coupled faces into a list indexed from 0 to N-1.
 *   - Assign face/model k to MPI rank k % nProcs.
 *   - Each rank reads and initializes only its owned model(s), with no MPI
 *     synchronization, so all ranks work simultaneously.
 * - Batch metadata exchange:
 *   - After all ranks finish initializing their owned model(s), share
 *     system_size and coupled_dof via MPI_Bcast so that all ranks know the
 *     sizes needed for subsequent result broadcasts.
 *
 * Time stepping in calc_svOneD() also has two phases:
 * - Parallel solve:
 *   - Each rank runs the 1D solve for its owned model(s), with no MPI calls,
 *     so different models can run concurrently on different ranks.
 * - Batch result exchange:
 *   - After all ranks finish solving, results are shared via MPI_Bcast so that
 *     all ranks know each result and can update the corresponding coupled
 *     boundary condition value.
 *
 * @par Parameters passed to run_1d_simulation_step_1d_
 * - params[0]: number of time points, currently 2.0.
 * - params[1]: t_old, the time at the start of the step.
 * - params[2]: t_new, the time at the end of the step.
 * - params[3]: BC_val_old, the coupled Q or P value at t_old.
 * - params[4]: BC_val_new, the coupled Q or P value at t_new.
 */
namespace svOneD {

/// @brief Initialize the 1D solver and populate the initial cplBC state.
/// Called once from baf_ini() after the BC data structures are set up.
void init_svOneD(ComMod& com_mod, const CmMod& cm_mod);

/// @brief Advance the 1D solver by one time step and update the coupled BC value.
///
/// @param BCFlag  'D' - derivative / perturbation step (state is NOT committed).
///               'L' - last Newton iteration (state IS committed, time advances).
void calc_svOneD(ComMod& com_mod, const CmMod& cm_mod, char BCFlag);

}  // namespace svOneD

#endif  // SV1D_SUBROUTINES_H
