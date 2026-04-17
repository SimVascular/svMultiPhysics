// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// 3D-1D coupling subroutines.
//
// These routines interface the 3D finite-element solver (svMultiPhysics) with
// the 1D blood-flow solver (svOneDSolver) via a dynamically loaded shared
// library (libsvoned_interface.so/.dylib).
//
// Coupling overview
// -----------------
//   NEU coupling (1D inlet driven by 3D outflow):
//     3D → 1D : flow rate Q  (as params[3..4])
//     1D → 3D : pressure P   (as cpl_value)
//     BC applied on 3D face  : Neumann (pressure traction)
//
//   DIR coupling (1D outlet driven by 3D pressure):
//     3D → 1D : pressure P   (as params[3..4])
//     1D → 3D : flow rate Q  (as cpl_value)
//     BC applied on 3D face  : Dirichlet (velocity profile)
//
// Parallelism model
// -----------------
//   Unlike the 0D solver (which is solved once on the master rank), each
//   1D model is INDEPENDENT and has its own input file.  Multiple 1D models
//   are therefore read, initialized, and solved in parallel:
//
//   Initialization (init_svOneD):
//     Phase 1 – parallel init:
//       - Collect all svOneD-coupled faces into a list indexed 0..N-1.
//       - Assign face (model) k to MPI rank  k % nProcs.
//       - Each rank reads and initializes ONLY its owned model(s) with no
//         MPI synchronization, so all ranks work simultaneously.
//     Phase 2 – batch metadata exchange:
//       - After every rank has finished initializing its own model(s),
//         share system_size and coupled_dof via MPI_Bcast so that all
//         ranks know the sizes needed for subsequent result broadcasts.
//
//   Time-stepping (calc_svOneD):
//     Phase 1 – parallel solve:
//       - Each rank runs run_step() for its owned model(s) with no MPI
//         calls, so model k on rank A and model k+1 on rank B truly run
//         concurrently.
//     Phase 2 – batch result exchange:
//       - After every rank has finished solving, results are shared via
//         MPI_Bcast so that ALL ranks know the result, and each BC's
//         coupled_bc.set_pressure() is updated accordingly.
//
// params array passed to run_1d_simulation_step_1d_:
//   params[0] = 2.0          (number of time points)
//   params[1] = t_old         (time at start of step)
//   params[2] = t_new         (time at end of step)
//   params[3] = BC_val_old    (Q or P at t_old)
//   params[4] = BC_val_new    (Q or P at t_new)

#include "svOneD_subroutines.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "ComMod.h"
#include "consts.h"
#include "utils.h"
#include "svOneD_interface/OneDSolverInterface.h"

#include "mpi.h"

namespace svOneD {

// ---------------------------------------------------------------------------
// Per-model state.  Each entry corresponds to one svOneD-coupled face (one 1D
// model).  Indexed by the sequential order in which coupled faces were found
// in eq[0].bc[].
// ---------------------------------------------------------------------------

struct OneDModelState {
  // Interface object (null on ranks that do not own this model).
  OneDSolverInterface* interface = nullptr;

  // Problem identifier returned by the 1D library.
  int problem_id = 0;

  // Total DOF count (nodes × 2).
  int system_size = 0;

  // "NEU" or "DIR".
  std::string coupling_type;

  // Index in the solution vector that corresponds to the coupled BC DOF.
  int coupled_dof = 0;

  // Authoritative solution vector (only valid on the owning rank).
  std::vector<double> solution;

  // Owning MPI rank for this model.
  int owner_rank = 0;

  // Index into eq[0].bc[] for the BC this model services.
  int iBc = -1;
};

// ---------------------------------------------------------------------------
// Module-level state.
// ---------------------------------------------------------------------------

// One entry per svOneD-coupled face, filled during init_svOneD().
static std::vector<OneDModelState> oned_models;

// Shared library handle (one per process, loaded once).
static OneDSolverInterface* shared_lib_instance = nullptr;

// Simulation time (advanced only on 'L' steps).
static double svOneDTime = 0.0;

// ---------------------------------------------------------------------------
// Helper: resolve the shared-library path (.so / .dylib / as-is).
// ---------------------------------------------------------------------------
static std::string resolve_lib_path(const std::string& lib_base)
{
  if (std::ifstream(lib_base + ".so").good())  return lib_base + ".so";
  if (std::ifstream(lib_base + ".dylib").good()) return lib_base + ".dylib";
  return lib_base;  // already has extension, or will fail at dlopen time
}

// ---------------------------------------------------------------------------
// init_svOneD
// ---------------------------------------------------------------------------
void init_svOneD(ComMod& com_mod, const CmMod& cm_mod)
{
  using namespace consts;

  auto& cplBC    = com_mod.cplBC;
  auto& solver_if = cplBC.sv1d_solver_interface;
  auto& cm       = com_mod.cm;
  const int nProcs = cm.nProcs;
  const int myRank = cm.taskId;

  if (!solver_if.has_data) {
    throw std::runtime_error("[svOneD::init_svOneD] svOneD solver interface data is missing.");
  }

  // Initialize the 1D simulation clock from the 3D solver's current time so
  // that restarts and non-zero start times are handled correctly.
  svOneDTime = com_mod.time;

  // ----- Collect the list of svOneD-coupled faces -----
  // Iterate over eq[0]'s BCs and pick those with iBC_Coupled and a non-empty
  // oned_input_file (stored in coupled_bc).
  {
    const int iEq = 0;
    const auto& eq = com_mod.eq[iEq];
    for (int iBc = 0; iBc < eq.nBc; iBc++) {
      const auto& bc = eq.bc[iBc];
      if (!utils::btest(bc.bType, iBC_Coupled)) continue;
      if (bc.coupled_bc.get_oned_input_file().empty()) continue;

      OneDModelState st;
      st.iBc = iBc;
      st.coupling_type = (bc.coupled_bc.get_bc_type() == BoundaryConditionType::bType_Neu) ? "NEU" : "DIR";
      oned_models.push_back(std::move(st));
    }
  }

  if (oned_models.empty()) {
    throw std::runtime_error("[svOneD::init_svOneD] No svOneD-coupled faces with input files found.");
  }

  // ----- Guard: require at least one MPI rank per 1D model -----
  // Each rank owns exactly one model (owner_rank = k % nProcs).  If nProcs < N
  // a single rank would own multiple models and call shared_lib_instance->initialize()
  // more than once, corrupting the static problem-ID state inside the shared library.
  const int nTotalModels = static_cast<int>(oned_models.size());
  if (nProcs < nTotalModels) {
    throw std::runtime_error(
        "[svOneD::init_svOneD] Number of MPI processes (" + std::to_string(nProcs) +
        ") is less than the number of svOneD-coupled faces (" +
        std::to_string(nTotalModels) +
        ").  Please run with at least " + std::to_string(nTotalModels) +
        " MPI processes.");
  }

  // ----- Load shared library (once per process) -----
  const std::string lib_path = resolve_lib_path(solver_if.solver_library);
  shared_lib_instance = new OneDSolverInterface();
  shared_lib_instance->load_library(lib_path);

  // ----- Assign ranks and initialize owned models (Phase 1: parallel) -----
  // No MPI calls in this loop.  All ranks proceed simultaneously, each
  // reading and initializing only the model(s) it owns.  Rank k owns model k
  // (assigned via k % nProcs), so for N models and N ranks every rank handles
  // exactly one model with no inter-rank synchronization.
  const int iEq = 0;
  auto& eq = com_mod.eq[iEq];
  for (int k = 0; k < nTotalModels; k++) {
    auto& st = oned_models[k];
    st.owner_rank = k % nProcs;

    if (myRank != st.owner_rank) continue;

    // This rank owns model k: read the input file and initialize.
    const std::string& input_file = eq.bc[st.iBc].coupled_bc.get_oned_input_file();
    int problem_id  = 0;
    int system_size = 0;

    shared_lib_instance->initialize(input_file, problem_id, system_size,
                                    st.coupling_type);
    st.problem_id   = problem_id;
    st.system_size  = system_size;
    st.interface    = shared_lib_instance;

    shared_lib_instance->set_external_step_size(problem_id, com_mod.dt);
    shared_lib_instance->extract_coupled_dof(problem_id, st.coupled_dof,
                                             st.coupling_type);

    st.solution.resize(system_size, 0.0);
    shared_lib_instance->return_solution(problem_id, st.solution.data(), system_size);

    // Initial coupled value = 0; first calc_svOneD call sets the real value.
    eq.bc[st.iBc].coupled_bc.set_pressure(0.0);
  }

  // ----- Broadcast metadata for all models (Phase 2: batch exchange) -----
  // All initialization is complete.  Now share system_size and coupled_dof
  // from each owner so that every rank knows the sizes needed for consistent
  // result broadcasts in calc_svOneD.
  for (int k = 0; k < nTotalModels; k++) {
    auto& st = oned_models[k];
    MPI_Bcast(&st.system_size,  1, MPI_INT, st.owner_rank, cm.com());
    MPI_Bcast(&st.coupled_dof,  1, MPI_INT, st.owner_rank, cm.com());
  }

  // Run one 'D' step to populate the initial resistance term bc.r.
  if (cplBC.schm != CplBCType::cplBC_E) {
    calc_svOneD(com_mod, cm_mod, 'D');
  }
}

// ---------------------------------------------------------------------------
// calc_svOneD
// ---------------------------------------------------------------------------
void calc_svOneD(ComMod& com_mod, const CmMod& cm_mod, char BCFlag)
{
  using namespace consts;

  auto& cm     = com_mod.cm;
  const int myRank = cm.taskId;

  const double t_old = svOneDTime;
  const double t_new = svOneDTime + com_mod.dt;
  const int    nTotalModels = static_cast<int>(oned_models.size());

  const int iEq = 0;
  auto& eq = com_mod.eq[iEq];

  // ----- Phase 1: each rank runs its own models without blocking -----
  // All ranks proceed through this loop simultaneously, each executing only
  // the models it owns.  No MPI call here, so model k on rank A and model k+1
  // on rank B truly run at the same time.
  std::vector<double> cpl_values(nTotalModels, 0.0);

  for (int k = 0; k < nTotalModels; k++) {
    auto& st = oned_models[k];
    auto& bc = eq.bc[st.iBc];

    if (myRank != st.owner_rank) continue;

    // Build params = [2, t_old, t_new, BC_val_old, BC_val_new]
    double params[5];
    params[0] = 2.0;
    params[1] = t_old;
    params[2] = t_new;

    if (bc.coupled_bc.get_bc_type() == BoundaryConditionType::bType_Neu) {
      params[3] = bc.coupled_bc.get_Qo();
      params[4] = bc.coupled_bc.get_Qn();
    } else {
      params[3] = bc.coupled_bc.get_Po();
      params[4] = bc.coupled_bc.get_Pn();
    }

    // Working copy of solution so that 'D' steps don't corrupt the
    // committed state.
    std::vector<double> work_sol = st.solution;
    st.interface->update_solution(st.problem_id, work_sol.data(), st.system_size);

    int save_flag  = (BCFlag == 'L') ? 1 : 0;
    int error_code = 0;

    st.interface->run_step(st.problem_id, t_old, save_flag,
                           st.coupling_type, params,
                           work_sol.data(), cpl_values[k], error_code);

    if (error_code != 0) {
      throw std::runtime_error(
          "[svOneD::calc_svOneD] 1D solver step for face '" +
          bc.coupled_bc.get_oned_input_file() + "' failed with error code " +
          std::to_string(error_code));
    }

    // Commit the updated solution only on the final iteration.
    if (BCFlag == 'L') {
      st.solution = work_sol;
    }
  }

  // ----- Phase 2: broadcast all results and update coupled BCs -----
  // After every rank has finished solving its own models, gather the results.
  // Each MPI_Bcast here is a cheap scalar transfer; the expensive 1D solver
  // work has already been done concurrently in Phase 1.
  for (int k = 0; k < nTotalModels; k++) {
    auto& st = oned_models[k];
    MPI_Bcast(&cpl_values[k], 1, MPI_DOUBLE, st.owner_rank, cm.com());
    eq.bc[st.iBc].coupled_bc.set_pressure(cpl_values[k]);
  }

  // Advance the simulation clock after the final iteration.
  if (BCFlag == 'L') {
    svOneDTime += com_mod.dt;
  }
}

}  // namespace svOneD
