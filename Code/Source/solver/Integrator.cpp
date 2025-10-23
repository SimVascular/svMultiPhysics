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
  : simulation_(simulation), inner_count_(0)
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

  // Inner loop for Newton iteration
  inner_count_ = 1;
  int reply;
  int iEqOld;

  // Looping over Newton iterations
  while (true) {
    #ifdef debug_integrator_step
    dmsg << "---------- Inner Loop " + std::to_string(inner_count_) << " -----------" << std::endl;
    dmsg << "cEq: " << cEq;
    dmsg << "com_mod.eq[cEq].sym: " << com_mod.eq[cEq].sym;
    #endif

    auto istr = "_" + std::to_string(cTS) + "_" + std::to_string(inner_count_);
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
    #ifdef debug_integrator_step
    dmsg << "Initiator step ..." << std::endl;
    #endif
    initiator_step();
    Ag_.write("Ag_pic" + istr);
    Yg_.write("Yg_pic" + istr);
    Dg_.write("Dg_pic" + istr);
    Yn.write("Yn_pic" + istr);

    if (com_mod.Rd.size() != 0) {
      com_mod.Rd = 0.0;
      com_mod.Kd = 0.0;
    }

    // Allocate com_mod.R and com_mod.Val arrays
    #ifdef debug_integrator_step
    dmsg << "Allocating the RHS and LHS" << std::endl;
    #endif
    allocate_linear_system(eq);
    com_mod.Val.write("Val_alloc" + istr);

    // Compute body forces
    #ifdef debug_integrator_step
    dmsg << "Set body forces ..." << std::endl;
    #endif
    set_body_forces();
    com_mod.Val.write("Val_bf" + istr);

    // Assemble equations
    #ifdef debug_integrator_step
    dmsg << "Assembling equation: " << eq.sym;
    #endif
    assemble_equations();
    com_mod.R.write("R_as" + istr);
    com_mod.Val.write("Val_as" + istr);

    // Treatment of boundary conditions on faces
    #ifdef debug_integrator_step
    dmsg << "Apply boundary conditions ..." << std::endl;
    #endif
    apply_boundary_conditions();
    com_mod.Val.write("Val_neu" + istr);
    com_mod.R.write("R_neu" + istr);
    Yg_.write("Yg_neu" + istr);
    Dg_.write("Dg_neu" + istr);

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
    #ifdef debug_integrator_step
    dmsg << "Update res() and incL ..." << std::endl;
    #endif
    update_residual_arrays(eq);

    // Solve equation
    #ifdef debug_integrator_step
    dmsg << "Solving equation: " << eq.sym;
    #endif
    solve_linear_system();
    com_mod.Val.write("Val_solve" + istr);
    com_mod.R.write("R_solve" + istr);

    // Solution is obtained, now updating (Corrector) and check for convergence
    #ifdef debug_integrator_step
    dmsg << "Update corrector ..." << std::endl;
    #endif
    bool all_converged = corrector_and_check_convergence();
    com_mod.Yn.write("Yn_picc" + istr);

    // Check if all equations converged
    if (all_converged) {
      #ifdef debug_integrator_step
      dmsg << ">>> All OK" << std::endl;
      dmsg << "iEqOld: " << iEqOld + 1;
      #endif
      return true;
    }

    output::output_result(simulation_, com_mod.timeP, 2, iEqOld);
    inner_count_ += 1;
  } // End of inner loop

  return false;
}

//------------------------
// initiator_step
//------------------------
void Integrator::initiator_step() {
  pic::pici(simulation_, Ag_, Yg_, Dg_);
}

//------------------------
// allocate_linear_system
//------------------------
void Integrator::allocate_linear_system(eqType& eq) {
  ls_ns::ls_alloc(simulation_->com_mod, eq);
}

//------------------------
// set_body_forces
//------------------------
void Integrator::set_body_forces() {
  bf::set_bf(simulation_->com_mod, Dg_);
}

//------------------------
// assemble_equations
//------------------------
void Integrator::assemble_equations() {
  auto& com_mod = simulation_->com_mod;
  auto& cep_mod = simulation_->get_cep_mod();

  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    eq_assem::global_eq_assem(com_mod, cep_mod, com_mod.msh[iM], Ag_, Yg_, Dg_);
  }
}

//------------------------
// apply_boundary_conditions
//------------------------
void Integrator::apply_boundary_conditions() {
  auto& com_mod = simulation_->com_mod;
  auto& cm_mod = simulation_->cm_mod;

  Yg_.write("Yg_vor_neu_" + std::to_string(com_mod.cTS) + "_" + std::to_string(inner_count_));
  Dg_.write("Dg_vor_neu_" + std::to_string(com_mod.cTS) + "_" + std::to_string(inner_count_));

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
}

//------------------------
// solve_linear_system
//------------------------
void Integrator::solve_linear_system() {
  auto& com_mod = simulation_->com_mod;
  auto& eq = com_mod.eq[com_mod.cEq];

  ls_ns::ls_solve(com_mod, eq, incL_, res_);
}

//------------------------
// corrector_and_check_convergence
//------------------------
bool Integrator::corrector_and_check_convergence() {
  auto& com_mod = simulation_->com_mod;

  pic::picc(simulation_);

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
