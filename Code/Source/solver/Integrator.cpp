// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "Integrator.h"
#include "all_fun.h"
#include "bf.h"
#include "contact.h"
#include "eq_assem.h"
#include "fs.h"
#include "ls.h"
#include "output.h"
#include "pic.h"
#include "ris.h"
#include "set_bc.h"
#include "ustruct.h"

#include <algorithm>
#include <iostream>
#include <set>

using namespace consts;

//------------------------
// Integrator Constructor
//------------------------
Integrator::Integrator(Simulation* simulation)
  : simulation_(simulation), newton_count_(0)
{
  initialize_arrays();
}

//------------------------
// Integrator Destructor
//------------------------
Integrator::~Integrator() {
  // Arrays will be automatically cleaned up
}

//------------------------
// initialize_arrays
//------------------------
void Integrator::initialize_arrays() {
  auto& com_mod = simulation_->com_mod;
  int tDof = com_mod.tDof;
  int tnNo = com_mod.tnNo;
  int nFacesLS = com_mod.nFacesLS;

  Ag_.resize(tDof, tnNo);
  Yg_.resize(tDof, tnNo);
  Dg_.resize(tDof, tnNo);
  res_.resize(nFacesLS);
  incL_.resize(nFacesLS);
}

//------------------------
// step
//------------------------
/// @brief Execute one Newton iteration loop for the current time step
bool Integrator::step() {
  using namespace consts;

  auto& com_mod = simulation_->com_mod;
  auto& cm_mod = simulation_->cm_mod;
  auto& cep_mod = simulation_->get_cep_mod();

  auto& An = com_mod.An;
  auto& Yn = com_mod.Yn;
  auto& Dn = com_mod.Dn;

  int& cTS = com_mod.cTS;
  int& cEq = com_mod.cEq;

  #define n_debug_integrator_step
  #ifdef debug_integrator_step
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  // Newton iteration loop
  newton_count_ = 1;
  int reply;
  int iEqOld;

  // Looping over Newton iterations
  while (true) {
    #ifdef debug_integrator_step
    dmsg << "---------- Newton Iteration " + std::to_string(newton_count_) << " -----------" << std::endl;
    dmsg << "cEq: " << cEq;
    dmsg << "com_mod.eq[cEq].sym: " << com_mod.eq[cEq].sym;
    #endif

    istr_ = "_" + std::to_string(cTS) + "_" + std::to_string(newton_count_);
    iEqOld = cEq;
    auto& eq = com_mod.eq[cEq];

    if (com_mod.cplBC.coupled && cEq == 0) {
      #ifdef debug_integrator_step
      dmsg << "Set coupled BCs " << std::endl;
      #endif
      set_bc::set_bc_cpl(com_mod, cm_mod);
      set_bc::set_bc_dir(com_mod, An, Yn, Dn);
    }

    // Initiator step for Generalized α-Method (quantities at n+am, n+af).
    initiator_step();

    if (com_mod.Rd.size() != 0) {
      com_mod.Rd = 0.0;
      com_mod.Kd = 0.0;
    }

    // Allocate com_mod.R and com_mod.Val arrays
    allocate_linear_system(eq);

    // Compute body forces
    set_body_forces();

    // Assemble equations
    assemble_equations();

    // Treatment of boundary conditions on faces
    apply_boundary_conditions();

    // Synchronize R across processes
    if (!eq.assmTLS) {
      #ifdef debug_integrator_step
      dmsg << "Synchronize R across processes ..." << std::endl;
      #endif
      all_fun::commu(com_mod, com_mod.R);
    }

    // Update residual in displacement equation for USTRUCT phys
    #ifdef debug_integrator_step
    dmsg << "com_mod.sstEq: " << com_mod.sstEq;
    #endif
    if (com_mod.sstEq) {
      ustruct::ustruct_r(com_mod, Yg_);
    }

    // Set the residual of the continuity equation to 0 on edge nodes
    if (std::set<EquationType>{Equation_stokes, Equation_fluid, Equation_ustruct, Equation_FSI}.count(eq.phys) != 0) {
      #ifdef debug_integrator_step
      dmsg << "thood_val_rc ..." << std::endl;
      #endif
      fs::thood_val_rc(com_mod);
    }

    // Treat Neumann boundaries that are not deforming
    #ifdef debug_integrator_step
    dmsg << "set_bc_undef_neu ..." << std::endl;
    #endif
    set_bc::set_bc_undef_neu(com_mod);

    // Update residual and increment arrays
    update_residual_arrays(eq);

    // Solve equation
    solve_linear_system();

    // Solution is obtained, now updating (Corrector) and check for convergence
    bool all_converged = corrector_and_check_convergence();

    // Check if all equations converged
    if (all_converged) {
      #ifdef debug_integrator_step
      dmsg << ">>> All OK" << std::endl;
      dmsg << "iEqOld: " << iEqOld + 1;
      #endif
      return true;
    }

    output::output_result(simulation_, com_mod.timeP, 2, iEqOld);
    newton_count_ += 1;
  } // End of Newton iteration loop

  return false;
}

//------------------------
// initiator_step
//------------------------
void Integrator::initiator_step() {
  #define n_debug_integrator_step
  #ifdef debug_integrator_step
  DebugMsg dmsg(__func__, simulation_->com_mod.cm.idcm());
  dmsg << "Initiator step ..." << std::endl;
  #endif

  pic::pici(simulation_, Ag_, Yg_, Dg_);

  // Debug output
  Ag_.write("Ag_pic" + istr_);
  Yg_.write("Yg_pic" + istr_);
  Dg_.write("Dg_pic" + istr_);
  simulation_->com_mod.Yn.write("Yn_pic" + istr_);
}

//------------------------
// allocate_linear_system
//------------------------
void Integrator::allocate_linear_system(eqType& eq) {
  #define n_debug_integrator_step
  #ifdef debug_integrator_step
  DebugMsg dmsg(__func__, simulation_->com_mod.cm.idcm());
  dmsg << "Allocating the RHS and LHS" << std::endl;
  #endif

  ls_ns::ls_alloc(simulation_->com_mod, eq);

  // Debug output
  simulation_->com_mod.Val.write("Val_alloc" + istr_);
}

//------------------------
// set_body_forces
//------------------------
void Integrator::set_body_forces() {
  #define n_debug_integrator_step
  #ifdef debug_integrator_step
  DebugMsg dmsg(__func__, simulation_->com_mod.cm.idcm());
  dmsg << "Set body forces ..." << std::endl;
  #endif

  bf::set_bf(simulation_->com_mod, Dg_);

  // Debug output
  simulation_->com_mod.Val.write("Val_bf" + istr_);
}

//------------------------
// assemble_equations
//------------------------
void Integrator::assemble_equations() {
  auto& com_mod = simulation_->com_mod;
  auto& cep_mod = simulation_->get_cep_mod();

  #define n_debug_integrator_step
  #ifdef debug_integrator_step
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg << "Assembling equation: " << com_mod.eq[com_mod.cEq].sym;
  #endif

  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    eq_assem::global_eq_assem(com_mod, cep_mod, com_mod.msh[iM], Ag_, Yg_, Dg_);
  }

  // Debug output
  com_mod.R.write("R_as" + istr_);
  com_mod.Val.write("Val_as" + istr_);
}

//------------------------
// apply_boundary_conditions
//------------------------
void Integrator::apply_boundary_conditions() {
  auto& com_mod = simulation_->com_mod;
  auto& cm_mod = simulation_->cm_mod;

  #define n_debug_integrator_step
  #ifdef debug_integrator_step
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg << "Apply boundary conditions ..." << std::endl;
  #endif

  Yg_.write("Yg_vor_neu" + istr_);
  Dg_.write("Dg_vor_neu" + istr_);

  // Apply Neumman or Traction boundary conditions
  set_bc::set_bc_neu(com_mod, cm_mod, Yg_, Dg_);

  // Apply CMM BC conditions
  if (!com_mod.cmmInit) {
    set_bc::set_bc_cmm(com_mod, cm_mod, Ag_, Dg_);
  }

  // Apply weakly applied Dirichlet BCs
  set_bc::set_bc_dir_w(com_mod, Yg_, Dg_);

  if (com_mod.risFlag) {
    ris::ris_resbc(com_mod, Yg_, Dg_);
  }

  if (com_mod.ris0DFlag) {
    ris::ris0d_bc(com_mod, cm_mod, Yg_, Dg_);
  }

  // Apply contact model and add its contribution to residual
  if (com_mod.iCntct) {
    contact::construct_contact_pnlty(com_mod, cm_mod, Dg_);
  }

  // Debug output
  com_mod.Val.write("Val_neu" + istr_);
  com_mod.R.write("R_neu" + istr_);
  Yg_.write("Yg_neu" + istr_);
  Dg_.write("Dg_neu" + istr_);
}

//------------------------
// solve_linear_system
//------------------------
void Integrator::solve_linear_system() {
  auto& com_mod = simulation_->com_mod;
  auto& eq = com_mod.eq[com_mod.cEq];

  #define n_debug_integrator_step
  #ifdef debug_integrator_step
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg << "Solving equation: " << eq.sym;
  #endif

  ls_ns::ls_solve(com_mod, eq, incL_, res_);

  // Debug output
  com_mod.Val.write("Val_solve" + istr_);
  com_mod.R.write("R_solve" + istr_);
}

//------------------------
// corrector_and_check_convergence
//------------------------
bool Integrator::corrector_and_check_convergence() {
  auto& com_mod = simulation_->com_mod;

  #define n_debug_integrator_step
  #ifdef debug_integrator_step
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg << "Update corrector ..." << std::endl;
  #endif

  pic::picc(simulation_);

  // Debug output
  com_mod.Yn.write("Yn_picc" + istr_);

  // Check if all equations converged
  return std::count_if(com_mod.eq.begin(), com_mod.eq.end(),
                       [](eqType& eq) { return eq.ok; }) == com_mod.eq.size();
}

//------------------------
// update_residual_arrays
//------------------------
void Integrator::update_residual_arrays(eqType& eq) {
  auto& com_mod = simulation_->com_mod;
  int nFacesLS = com_mod.nFacesLS;
  double dt = com_mod.dt;

  #define n_debug_integrator_step
  #ifdef debug_integrator_step
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg << "Update res() and incL ..." << std::endl;
  #endif

  incL_ = 0;
  if (eq.phys == Equation_mesh) {
    incL_(nFacesLS - 1) = 1;
  }

  if (com_mod.cmmInit) {
    incL_(nFacesLS - 1) = 1;
  }

  for (int iBc = 0; iBc < eq.nBc; iBc++) {
    int i = eq.bc[iBc].lsPtr;
    if (i != -1) {
      // Resistance term for coupled Neumann BC tangent contribution
      res_(i) = eq.gam * dt * eq.bc[iBc].r;
      incL_(i) = 1;
    }
  }
}
