// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "svOneD_interface.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <filesystem>

#include "ComMod.h"
#include "consts.h"
#include "utils.h"
#include "svOneD_interface/OneDSolverInterface.h"
#include "Core/Exception.h"

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

  // Pressure ramp for 1D coupling initialization.
  // DIR: over the first ramp_steps committed time steps the pressure sent to
  //      the 1D solver is linearly interpolated from ramp_ref_pressure to the
  //      actual 3D pressure value.
  // NEU: over the first ramp_steps committed time steps the flow rate sent to
  //      the 1D solver is linearly interpolated from 0 to the actual 3D Q.
  // Zero means no ramping.
  int    ramp_steps = 0;
  double ramp_ref_pressure = 0.0;
  int    step_count = 0;  ///< Number of committed (BCFlag=='L') steps taken.

  // Under-relaxation (omega in (0, 1]).  Default 1.0 = no relaxation.
  //
  // DIR coupling:
  //   Input  (P sent to 1D): P_sent  = omega * P_target + (1-omega) * P_prev_sent
  //   Output (Q from 1D)   : Q_relax = omega * Q_raw    + (1-omega) * Q_prev_sent
  //
  // NEU coupling:
  //   Output (P from 1D)   : P_applied = omega * P_target + (1-omega) * P_neu_prev
  //   (P_target already includes the ramp from ramp_ref_pressure to P_raw)
  double relax_factor = 1.0;
  double P_prev_sent_old = 0.0;  ///< Under-relaxed pressure sent at params[3] (t_old) on last 'L' step (DIR).
  double P_prev_sent_new = 0.0;  ///< Under-relaxed pressure sent at params[4] (t_new) on last 'L' step (DIR).
  double Q_prev_sent = 0.0;      ///< Under-relaxed flow rate output on last 'L' step (DIR only).
  double P_neu_prev  = 0.0;      ///< Under-relaxed pressure applied to 3D on last 'L' step (NEU only).
  double Q_prev_sent_old = 0.0;  ///< Under-relaxed Q sent at params[3] (t_old) on last 'L' step (NEU).
  double Q_prev_sent_new = 0.0;  ///< Under-relaxed Q sent at params[4] (t_new) on last 'L' step (NEU).
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
  if (std::filesystem::is_regular_file(lib_base + ".so")) {
    return lib_base + ".so";
  }

  if (std::filesystem::is_regular_file(lib_base + ".dylib")) {
    return lib_base + ".dylib";
  }
  return lib_base;  // already has extension, or will fail at dlopen time
}

// ---------------------------------------------------------------------------
// init_svOneD
// ---------------------------------------------------------------------------
void init_svOneD(ComMod& com_mod, const CmMod& cm_mod)
{
  using namespace consts;

  auto& cplBC    = com_mod.cplBC;
  auto& solver_if = cplBC.svOneD_solver_interface;
  auto& cm       = com_mod.cm;
  const int nProcs = cm.nProcs;
  const int myRank = cm.taskId;

  if (!solver_if.has_data) {
  throw svmp::CoreException(
      "[svOneD::init_svOneD] svOneD solver interface data is missing. Please check the XML input file.",
      svmp::StatusCode::InvalidState,
      __FILE__,
      __LINE__,
      __func__);
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
      st.ramp_steps         = bc.coupled_bc.get_oned_ramp_steps();
      st.ramp_ref_pressure  = bc.coupled_bc.get_oned_ramp_ref_pressure();
      st.relax_factor       = bc.coupled_bc.get_oned_relax_factor();
      oned_models.push_back(std::move(st));
    }
  }

  if (oned_models.empty()) {
  throw svmp::CoreException(
      "[svOneD::init_svOneD] No svOneD-coupled faces with input files found. Please check the XML input file.",
      svmp::StatusCode::InvalidState,
      __FILE__,
      __LINE__,
      __func__);
  }

  // ----- Guard: require at least one MPI rank per 1D model -----
  // Each rank owns exactly one model (owner_rank = k % nProcs).  If nProcs < N
  // a single rank would own multiple models and call shared_lib_instance->initialize()
  // more than once, corrupting the static problem-ID state inside the shared library.
  const int nTotalModels = static_cast<int>(oned_models.size());
  if (nProcs < nTotalModels) {
    throw svmp::CoreException(
        "[svOneD::init_svOneD] Number of MPI processes (" + std::to_string(nProcs) +
        ") is less than the number of svOneD-coupled faces (" +
        std::to_string(nTotalModels) +
        "). Please run with at least " + std::to_string(nTotalModels) +
        " MPI processes.",
        svmp::StatusCode::InvalidState,
        __FILE__,
        __LINE__,
        __func__);
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
      // NEU coupling: apply under-relaxation to the 3D flow rate Q sent to
      // the 1D solver to damp timestep-to-timestep oscillations in the input.
      // No ramping is applied here; ramping is applied to the *output* pressure
      // P that the 1D solver returns (Phase 2 below).
      const double omega = st.relax_factor;
      params[3] = omega * bc.coupled_bc.get_Qo() + (1.0 - omega) * st.Q_prev_sent_new;
      params[4] = omega * bc.coupled_bc.get_Qn() + (1.0 - omega) * st.Q_prev_sent_new;
    } else {
      double raw_P_old = bc.coupled_bc.get_Po();
      double raw_P_new = bc.coupled_bc.get_Pn();

      // Step 1: apply pressure ramp (scales amplitude from ramp_ref_pressure
      //         to actual 3D pressure over the first ramp_steps steps).
      double P_target_old, P_target_new;
      if (st.ramp_steps > 0) {
        double ramp_factor = std::min(1.0, static_cast<double>(st.step_count) / st.ramp_steps);
        double P_ref = st.ramp_ref_pressure;
        P_target_old = P_ref + ramp_factor * (raw_P_old - P_ref);
        P_target_new = P_ref + ramp_factor * (raw_P_new - P_ref);
      } else {
        P_target_old = raw_P_old;
        P_target_new = raw_P_new;
      }

      // Step 2: apply under-relaxation (damps timestep-to-timestep oscillations).
      // P_sent = omega * P_target + (1 - omega) * P_prev_sent
      const double omega = st.relax_factor;
      params[3] = omega * P_target_old + (1.0 - omega) * st.P_prev_sent_new;
      params[4] = omega * P_target_new + (1.0 - omega) * st.P_prev_sent_new;
    }

    // Working copy of solution so that 'D' steps don't corrupt the
    // committed state.
    std::vector<double> work_sol = st.solution;
    st.interface->update_solution(st.problem_id, work_sol.data(), st.system_size);

    int save_incr  = com_mod.saveIncr;
    int error_code = 0;

    try {
      st.interface->run_simulation(st.problem_id, t_old, save_incr,
                                  st.coupling_type, params,
                                  work_sol.data(), cpl_values[k], BCFlag, error_code);
    } catch (const std::exception& e) {
        throw svmp::CoreException(
            std::string("[svOneD::calc_svOneD] 1D solver step failed for svOneD input file '") +
            bc.coupled_bc.get_oned_input_file() + "' with coupling type '" +
            st.coupling_type + "' at time " + std::to_string(t_new) + ". " +
            "Original error: " + e.what(),
            svmp::StatusCode::DependencyError,
            __FILE__,
            __LINE__,
            __func__);
    }

    // Commit the updated solution only on the final iteration.
    if (BCFlag == 'L') {
      st.solution = work_sol;
      // Update the under-relaxation history with the values actually sent.
      // NEU: Q_prev_sent updated here; P history updated in Phase 2.
      // DIR: P_prev_sent tracks the pressure value sent to the 1D solver.
      if (st.coupling_type == "NEU") {
        st.Q_prev_sent_old = params[3];
        st.Q_prev_sent_new = params[4];
      } else {
        st.P_prev_sent_old = params[3];
        st.P_prev_sent_new = params[4];
      }
    }
  }

  // ----- Phase 2: broadcast all results and update coupled BCs -----
  // After every rank has finished solving its own models, gather the results.
  // Each MPI_Bcast here is a cheap scalar transfer; the expensive 1D solver
  // work has already been done concurrently in Phase 1.
  for (int k = 0; k < nTotalModels; k++) {
    auto& st = oned_models[k];
    MPI_Bcast(&cpl_values[k], 1, MPI_DOUBLE, st.owner_rank, cm.com());
    auto& cpl_bc = eq.bc[st.iBc].coupled_bc;
    if (cpl_bc.get_bc_type() == BoundaryConditionType::bType_Dir) {
      // 1D solver returns flow Q for DIR coupling; store it as flowrate so that
      // set_bc can read get_Qn() and build the nodal velocity profile.
      // Negate the sign: the 1D solver returns Q > 0 for inflow, but the 3D
      // code applies velocity as Q * gx * outward_normal, so Q must be negative
      // for an inlet face (outward normal points away from the domain).
      // This matches the svZeroD convention: in_out = -1 for DIR (outlet of 0D
      // = inlet of 3D), giving QCoupled = -1 * lpn_state_y[flow_id].
      double Q_raw = -cpl_values[k];
      // Apply under-relaxation to the Q output to damp timestep-to-timestep oscillations.
      const double omega = st.relax_factor;
      double Qn_relaxed = omega * Q_raw + (1.0 - omega) * st.Q_prev_sent;
      double Qo_prev = cpl_bc.get_Qn();
      cpl_bc.set_flowrates(Qo_prev, Qn_relaxed);
      if (BCFlag == 'L') {
        st.Q_prev_sent = Qn_relaxed;
      }
    } else {
      // 1D solver returns pressure P for NEU coupling.
            //
      // Step 1: apply pressure ramp (scales output P from ramp_ref_pressure
      //         to the actual 1D pressure over the first ramp_steps committed
      //         steps).  This prevents a large sudden pressure jump from being
      //         imposed on the 3D domain at startup, which is the primary cause
      //         of oscillations in Neumann coupling.
      double P_raw = cpl_values[k];
      double P_target;
      if (st.ramp_steps > 0) {
        double ramp_factor = std::min(1.0, static_cast<double>(st.step_count) / st.ramp_steps);
        double P_ref = st.ramp_ref_pressure;
        P_target = P_ref + ramp_factor * (P_raw - P_ref);
      } else {
        P_target = P_raw;
      }
      // Step 2: apply under-relaxation to damp timestep-to-timestep oscillations.
      // P_applied = omega * P_target + (1 - omega) * P_prev_applied
      const double omega = st.relax_factor;
      double P_relaxed = omega * P_target + (1.0 - omega) * st.P_neu_prev;
      cpl_bc.set_pressure(P_relaxed);
      if (BCFlag == 'L') {
        st.P_neu_prev = P_relaxed;
      }
    }
  }

  // Advance the simulation clock after the final iteration.
  if (BCFlag == 'L') {
    svOneDTime += com_mod.dt;
    for (auto& st : oned_models) {
      st.step_count++;
    }
  }
}

}  // namespace svOneD
