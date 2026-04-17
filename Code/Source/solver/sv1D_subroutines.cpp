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
//     1D → 3D : pressure P   (as cplBCvalue)
//     BC applied on 3D face  : Neumann (pressure traction)
//
//   DIR coupling (1D outlet driven by 3D pressure):
//     3D → 1D : pressure P   (as params[3..4])
//     1D → 3D : flow rate Q  (as cplBCvalue)
//     BC applied on 3D face  : Dirichlet (velocity profile)
//
// params array passed to run_1d_simulation_step_1d_:
//   params[0] = 2.0          (number of time points)
//   params[1] = t_old         (time at start of step)
//   params[2] = t_new         (time at end of step)
//   params[3] = BC_val_old    (Q or P at t_old)
//   params[4] = BC_val_new    (Q or P at t_new)

#include "sv1D_subroutines.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "ComMod.h"
#include "consts.h"
#include "utils.h"
#include "sv1D_interface/OneDSolverInterface.h"

namespace sv1D {

// ---------------------------------------------------------------------------
// Module-level state (one 1D solver instance per simulation).
// ---------------------------------------------------------------------------

static OneDSolverInterface* oned_interface = nullptr;

// Problem identifier assigned by the 1D library during initialization.
static int oned_model_id   = 0;

// Total DOF count of the 1D solution vector (nodes * 2: flow + area/pressure).
static int oned_system_size = 0;

// "NEU" or "DIR" — determined from the face BC type at initialization time.
static std::string oned_coupling_type;

// Index within the 1D solution vector that corresponds to the coupled DOF.
static int oned_coupled_dof = 0;

// 1D solution vector maintained by the 3D solver.
// This is the authoritative snapshot used to restore the 1D state before
// each call to run_step().
static std::vector<double> oned_solution;

// Simulation time tracked for the 1D solver (advanced only on 'L' steps).
static double sv1DTime = 0.0;

// ---------------------------------------------------------------------------
// Helper: find the first BC that uses cplBC.fa (iBC_cpl) and sv1D.
// Returns {ptr, bGrp} or throws.
// ---------------------------------------------------------------------------
static std::pair<int, consts::CplBCType> find_1d_face(const ComMod& com_mod)
{
  using namespace consts;
  const int iEq = 0;
  const auto& eq = com_mod.eq[iEq];
  for (int iBc = 0; iBc < eq.nBc; iBc++) {
    const auto& bc = eq.bc[iBc];
    int ptr = bc.cplBCptr;
    if (ptr == -1) continue;
    if (!utils::btest(bc.bType, iBC_cpl)) continue;
    return {ptr, com_mod.cplBC.fa[ptr].bGrp};
  }
  throw std::runtime_error("[sv1D::find_1d_face] No 1D-coupled face found in eq[0].");
}

// ---------------------------------------------------------------------------
// init_sv1D
// ---------------------------------------------------------------------------
void init_sv1D(ComMod& com_mod, const CmMod& cm_mod)
{
  using namespace consts;

  auto& cplBC    = com_mod.cplBC;
  auto& solver_if = cplBC.sv1d_solver_interface;
  auto& cm       = com_mod.cm;

  if (!solver_if.has_data) {
    throw std::runtime_error("[sv1D::init_sv1D] sv1D solver interface data is missing.");
  }

  // Determine the coupling direction (NEU or DIR) from the BC type of the
  // first coupled face.  Only master needs this but we resolve it globally.
  auto [ptr, bGrp] = find_1d_face(com_mod);
  oned_coupling_type = (bGrp == CplBCType::cplBC_Neu) ? "NEU" : "DIR";

  if (cm.mas(cm_mod)) {
    // ----- Load shared library -----
    const std::string lib_base = solver_if.solver_library;
    // Try .so first (Linux), then .dylib (macOS).
    std::string lib_path;
    {
      std::ifstream f_so(lib_base + ".so");
      std::ifstream f_dy(lib_base + ".dylib");
      if (f_so.good()) {
        lib_path = lib_base + ".so";
      } else if (f_dy.good()) {
        lib_path = lib_base + ".dylib";
      } else {
        // Try using the path as-is (already includes extension).
        lib_path = lib_base;
      }
    }

    oned_interface = new OneDSolverInterface();
    oned_interface->load_library(lib_path);

    // ----- Initialize 1D solver -----
    int problem_id  = 0;
    int system_size = 0;
    oned_interface->initialize(solver_if.input_file, problem_id,
                               system_size, oned_coupling_type);
    oned_model_id   = problem_id;
    oned_system_size = system_size;

    // ----- Synchronize time step -----
    oned_interface->set_external_step_size(oned_model_id, com_mod.dt);

    // ----- Retrieve coupled DOF index -----
    oned_interface->extract_coupled_dof(oned_model_id, oned_coupled_dof,
                                        oned_coupling_type);

    // ----- Copy initial 1D solution to local buffer -----
    oned_solution.resize(oned_system_size, 0.0);
    oned_interface->return_solution(oned_model_id, oned_solution.data(),
                                    oned_system_size);

    // ----- Set initial cplBC.fa[ptr].y -----
    // For NEU: initial BC value is the pressure at the coupled node.
    // For DIR: initial BC value is the flow at the coupled node.
    // The coupled DOF index is oned_coupled_dof (0-based) in the solution
    // vector, where even indices = flow, odd indices = area or pressure
    // depending on the 1D formulation.  The shared library returns the
    // correct cplBCvalue via run_step, so we initialise y = 0 here and let
    // the first calc_sv1D call set the real value.
    cplBC.fa[ptr].y = 0.0;
  }

  // Broadcast system_size and coupling metadata to slave processes.
  cm.bcast(cm_mod, &oned_system_size);
  cm.bcast(cm_mod, &oned_coupled_dof);

  // Slave processes also allocate the solution buffer (unused but kept for
  // any future parallel extension).
  if (cm.slv(cm_mod)) {
    oned_solution.resize(oned_system_size, 0.0);
  }

  // Run one 'D' step so that the initial resistance term bc.r is populated
  // before the first Newton iteration (matches the svZeroD pattern).
  if (cplBC.schm != CplBCType::cplBC_E) {
    // Flowrate/pressure are 0 at t=0; calc_sv1D handles this gracefully.
    calc_sv1D(com_mod, cm_mod, 'D');
  }
}

// ---------------------------------------------------------------------------
// calc_sv1D
// ---------------------------------------------------------------------------
void calc_sv1D(ComMod& com_mod, const CmMod& cm_mod, char BCFlag)
{
  using namespace consts;

  auto& cplBC = com_mod.cplBC;
  auto& cm    = com_mod.cm;

  // Only the master process drives the 1D solver.
  if (cm.mas(cm_mod)) {
    const int iEq = 0;
    auto& eq = com_mod.eq[iEq];

    double t_old = sv1DTime;
    double t_new = sv1DTime + com_mod.dt;

    for (int iBc = 0; iBc < eq.nBc; iBc++) {
      auto& bc = eq.bc[iBc];
      int ptr  = bc.cplBCptr;
      if (ptr == -1) continue;
      if (!utils::btest(bc.bType, iBC_cpl)) continue;

      // Build params = [2, t_old, t_new, BC_val_old, BC_val_new]
      // params[1..2] carry the time window; t_old is also passed as the
      // standalone current_time argument required by the 1D library API.
      double params[5];
      params[0] = 2.0;
      params[1] = t_old;
      params[2] = t_new;

      if (cplBC.fa[ptr].bGrp == CplBCType::cplBC_Neu) {
        // NEU: 3D provides flow Q → 1D returns pressure P.
        params[3] = cplBC.fa[ptr].Qo;
        params[4] = cplBC.fa[ptr].Qn;
      } else {
        // DIR: 3D provides pressure P → 1D returns flow Q.
        params[3] = cplBC.fa[ptr].Po;
        params[4] = cplBC.fa[ptr].Pn;
      }

      // Use a working copy of the solution so that perturbation steps
      // ('D') do not corrupt the committed state in oned_solution.
      std::vector<double> work_sol = oned_solution;
      oned_interface->update_solution(oned_model_id, work_sol.data(),
                                      oned_system_size);

      // Save flag: write VTK output only on the final 'L' step.
      int save_flag = (BCFlag == 'L') ? 1 : 0;

      double cpl_value = 0.0;
      int    error_code = 0;

      oned_interface->run_step(oned_model_id, t_old, save_flag,
                               oned_coupling_type, params,
                               work_sol.data(), cpl_value, error_code);

      if (error_code != 0) {
        throw std::runtime_error("[sv1D::calc_sv1D] 1D solver step failed with error code " +
                                 std::to_string(error_code));
      }

      // Store the returned BC value.
      cplBC.fa[ptr].y = cpl_value;

      // On the final iteration only, commit the updated solution so that
      // the next time step starts from the correct state.
      if (BCFlag == 'L') {
        oned_solution = work_sol;
      }
    }

    // Advance the 1D simulation clock on the final iteration.
    if (BCFlag == 'L') {
      sv1DTime += com_mod.dt;
    }
  }

  // Broadcast updated cplBC.fa[i].y to all parallel processes.
  if (!cm.seq()) {
    Vector<double> y(cplBC.nFa);

    if (cm.mas(cm_mod)) {
      for (int i = 0; i < cplBC.nFa; i++) {
        y(i) = cplBC.fa[i].y;
      }
    }

    cm.bcast(cm_mod, y);

    if (cm.slv(cm_mod)) {
      for (int i = 0; i < cplBC.nFa; i++) {
        cplBC.fa[i].y = y(i);
      }
    }
  }
}

}  // namespace sv1D
