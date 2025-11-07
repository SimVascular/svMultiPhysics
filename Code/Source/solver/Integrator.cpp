// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "Integrator.h"
#include "all_fun.h"
#include "bf.h"
#include "cep_ion.h"
#include "contact.h"
#include "eq_assem.h"
#include "fs.h"
#include "ls.h"
#include "nn.h"
#include "output.h"
#include "ris.h"
#include "set_bc.h"
#include "ustruct.h"
#include "utils.h"

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
  An_.resize(tDof, tnNo);
  Dn_.resize(tDof, tnNo);
  Yn_.resize(tDof, tnNo);
  res_.resize(nFacesLS);
  incL_.resize(nFacesLS);

  // Initialize An_ = Ao (same as what was done in initialize.cpp:760)
  An_ = com_mod.Ao;
  // Initialize Dn_ = Do (same as what was done in initialize.cpp:761)
  Dn_ = com_mod.Do;
  // Initialize Yn_ = Yo (same as what was done in initialize.cpp:760)
  Yn_ = com_mod.Yo;
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

  auto& An = An_;  // Use member variable
  auto& Yn = Yn_;
  auto& Dn = Dn_;

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
      set_bc::set_bc_cpl(com_mod, cm_mod, Yn);
      set_bc::set_bc_dir(com_mod, An, Yn, Dn);
    }

    // Initiator step for Generalized α-Method (quantities at n+am, n+af).
    #ifdef debug_integrator_step
    dmsg << "Initiator step ..." << std::endl;
    #endif
    initiator_step();

    if (com_mod.Rd.size() != 0) {
      com_mod.Rd = 0.0;
      com_mod.Kd = 0.0;
    }

    // Allocate com_mod.R and com_mod.Val arrays
    #ifdef debug_integrator_step
    dmsg << "Allocating the RHS and LHS" << std::endl;
    #endif
    allocate_linear_system(eq);

    // Compute body forces
    #ifdef debug_integrator_step
    dmsg << "Set body forces ..." << std::endl;
    #endif
    set_body_forces();

    // Assemble equations
    #ifdef debug_integrator_step
    dmsg << "Assembling equation: " << eq.sym;
    #endif
    assemble_equations();

    // Treatment of boundary conditions on faces
    #ifdef debug_integrator_step
    dmsg << "Apply boundary conditions ..." << std::endl;
    #endif
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
    #ifdef debug_integrator_step
    dmsg << "Update res() and incL ..." << std::endl;
    #endif
    update_residual_arrays(eq);

    // Solve equation
    #ifdef debug_integrator_step
    dmsg << "Solving equation: " << eq.sym;
    #endif
    solve_linear_system();

    // Solution is obtained, now updating (Corrector) and check for convergence
    #ifdef debug_integrator_step
    dmsg << "Update corrector ..." << std::endl;
    #endif
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
  pici(Ag_, Yg_, Dg_);

  // Debug output
  Ag_.write("Ag_pic" + istr_);
  Yg_.write("Yg_pic" + istr_);
  Dg_.write("Dg_pic" + istr_);
  Yn_.write("Yn_pic" + istr_);
}

//------------------------
// allocate_linear_system
//------------------------
void Integrator::allocate_linear_system(eqType& eq) {
  ls_ns::ls_alloc(simulation_->com_mod, eq);

  // Debug output
  simulation_->com_mod.Val.write("Val_alloc" + istr_);
}

//------------------------
// set_body_forces
//------------------------
void Integrator::set_body_forces() {
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

  auto& Yn = Yn_;

  Yg_.write("Yg_vor_neu" + istr_);
  Dg_.write("Dg_vor_neu" + istr_);

  // Apply Neumman or Traction boundary conditions
  set_bc::set_bc_neu(com_mod, cm_mod, Yg_, Dg_, Yn);

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
    ris::ris0d_bc(com_mod, cm_mod, Yg_, Dg_, Yn);
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

  picc();

  // Debug output
  Yn_.write("Yn_picc" + istr_);

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

//------------------------
// predictor (picp)
//------------------------
/// @brief Predictor step for next time step
///
/// Modifies:
///   pS0
///   Ad
///   Ao
///   Yo
///   Do
///   An
///   Yn
///   Dn
///
void Integrator::predictor()
{
  using namespace consts;

  auto& com_mod = simulation_->com_mod;

  #define n_debug_picp
  #ifdef debug_picp
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "pstEq: " << com_mod.pstEq;
  #endif

  // Variables for prestress calculations
  auto& pS0 = com_mod.pS0;
  auto& pSn = com_mod.pSn;

  // time derivative of displacement
  auto& Ad = com_mod.Ad;

  auto& Ao = com_mod.Ao;
  auto& An = An_;  // Use member variable
  auto& Yo = com_mod.Yo;
  auto& Yn = Yn_;
  auto& Do = com_mod.Do;
  auto& Dn = Dn_;

  // Prestress initialization
  if (com_mod.pstEq) {
     pS0 = pS0 + pSn;
     Ao = 0.0;
     Yo = 0.0;
     Do = 0.0;
  }

  // IB treatment: Set dirichlet BC and update traces. For explicit
  // coupling, compute FSI forcing and freeze it for the time step.
  // For implicit coupling, project IB displacement on background
  // mesh and predict quantities at next time step
  //
  // [NOTE] not implemented.
  /*
  if (ibFlag) {
    // Set IB Dirichlet BCs
    CALL IB_SETBCDIR(ib.Yb, ib.Ubo)

    // Update IB location and tracers
    CALL IB_UPDATE(Do)

    if (ib.cpld == ibCpld_E) {
      // FSI forcing for immersed bodies (explicit coupling)
      CALL IB_CALCFFSI(Ao, Yo, Do, ib.Auo, ib.Ubo)

    } else { if (ib.cpld == ibCpld_I) {
      // Project IB displacement (Ubn) to background mesh
      CALL IB_PRJCTU(Do)

    //  Predictor step for implicit coupling
    CALL IB_PICP()
    }
  }
  */

  const auto& dt = com_mod.dt;
  #ifdef debug_picp
  dmsg << "dt: " << dt;
  dmsg << "dFlag: " << com_mod.dFlag;
  #endif

  for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
    auto& eq = com_mod.eq[iEq];
    int s = eq.s; // start row
    int e = eq.e; // end row

    #ifdef debug_picp
    dmsg << "----- iEq " << iEq << " -----";
    dmsg << "s: " << s;
    dmsg << "e: " << e;
    dmsg << "eq.gam: " << eq.gam;
    dmsg << "coef: " << coef;
    #endif

    // [TODO:DaveP] careful here with s amd e.
    double coef = (eq.gam - 1.0) / eq.gam;
    for (int i = s; i <= e; i++) {
      for (int j = 0; j < Ao.ncols(); j++) {
        // eqn 87 of Bazilevs 2007
        An(i,j) = Ao(i,j) * coef;
      }
    }

    // electrophysiology
    if (eq.phys == Equation_CEP) {
      cep_ion::cep_integ(simulation_, iEq, e, Do);
    }

    // eqn 86 of Bazilevs 2007
    Yn.set_rows(s,e, Yo.rows(s,e));

    if (com_mod.dFlag) {

      // struct, lElas, FSI (struct, mesh)
      if (!com_mod.sstEq) {
        double coef = dt*dt*(0.5*eq.gam - eq.beta) / (eq.gam - 1.0);
        Dn.set_rows(s,e, Do.rows(s,e) + Yn.rows(s,e)*dt + An.rows(s,e)*coef);

      // ustruct, FSI
      //
      } else {

        if (eq.phys == Equation_ustruct || eq.phys == Equation_FSI) {
          double coef = (eq.gam - 1.0) / eq.gam;
          Ad = Ad*coef;
          Dn.set_rows(s,e, Do.rows(s,e));

        } else if (eq.phys == Equation_mesh) {
          double coef = dt*dt*(0.5*eq.gam - eq.beta) / (eq.gam - 1.0);
          Dn.set_rows(s,e, Do.rows(s,e) + Yn.rows(s,e)*dt + An.rows(s,e)*coef);
        }
      }
    } else {
      Dn.set_rows(s,e, Do.rows(s,e));
    }
  }
}

//------------------------
// pici
//------------------------
/// @brief Initiator for Generalized α-Method
///
/// Uses Generalized α− Method for time stepping.
///
/// Modifes Ag from combination of An and Ao defined by coefs from eq.am, eq.af,
///   Ag = (1 - eq.am) * Ao  +  eq.am * An
///   Yg = (1 - eq.af) * Yo  +  eq.af * Yn
///   Dg = (1 - eq.af) * Do  +  eq.af * Dn
///
/// Modifies:
///   Ag - acceleration
///   Yg - velocity
///   Dg - displacement
///
void Integrator::pici(Array<double>& Ag, Array<double>& Yg, Array<double>& Dg)
{
  using namespace consts;

  auto& com_mod = simulation_->com_mod;

  #define n_debug_pici
  #ifdef debug_pici
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  const int cEq = com_mod.cEq;
  const int tnNo = com_mod.tnNo;
  auto& eq = com_mod.eq[cEq];
  auto& dof = com_mod.dof;
  eq.itr = eq.itr + 1;

  // [NOTE] Setting gobal variable 'dof'.
  dof = eq.dof;
  #ifdef debug_pici
  dmsg << "cEq: " << cEq;
  dmsg << "eq.itr: " << eq.itr;
  dmsg << "dof: " << dof;
  dmsg << "tnNo: " << tnNo;
  dmsg << "com_mod.pstEq: " << com_mod.pstEq;
  #endif

  const auto& Ao = com_mod.Ao;
  const auto& An = An_;  // Use member variable
  const auto& Do = com_mod.Do;
  const auto& Dn = Dn_;
  const auto& Yo = com_mod.Yo;
  const auto& Yn = Yn_;

  for (int i = 0; i < com_mod.nEq; i++) {
    auto& eq = com_mod.eq[i];
    int s = eq.s;
    int e = eq.e;
    Vector<double> coef(4);
    coef(0) = 1.0 - eq.am;
    coef(1) = eq.am;
    coef(2) = 1.0 - eq.af;
    coef(3) = eq.af;
    #ifdef debug_pici
    dmsg << "s: " << s;
    dmsg << "e: " << e;
    dmsg << "coef: " << coef[0] << " " << coef[1] << " " << coef[2] << " " << coef[3];
    #endif

    if ((eq.phys == Equation_heatF) && (com_mod.usePrecomp)){
        for (int a = 0; a < tnNo; a++) {
            for (int j = 0; j < com_mod.nsd; j++) {
                //Ag(j, a) = An(j, a);
                //Yg(j, a) = Yn(j, a);
                //Dg(j, a) = Dn(j, a);
                Ag(j, a) = Ao(j, a) * coef(0) + An(j, a) * coef(1);
                Yg(j, a) = Yo(j, a) * coef(2) + Yn(j, a) * coef(3);
                Dg(j, a) = Do(j, a) * coef(2) + Dn(j, a) * coef(3);
            }
        }
        for (int a = 0; a < tnNo; a++) {
            for (int j = s; j <= e; j++) {
                Ag(j, a) = Ao(j, a) * coef(0) + An(j, a) * coef(1);
                Yg(j, a) = Yo(j, a) * coef(2) + Yn(j, a) * coef(3);
                Dg(j, a) = Do(j, a) * coef(2) + Dn(j, a) * coef(3);
            }
        }
    } else {
        for (int a = 0; a < tnNo; a++) {
            for (int j = s; j <= e; j++) {
                // eqn 89 of Bazilevs 2007
                Ag(j, a) = Ao(j, a) * coef(0) + An(j, a) * coef(1);

                // eqn 90 of Bazilevs 2007
                Yg(j, a) = Yo(j, a) * coef(2) + Yn(j, a) * coef(3);

                Dg(j, a) = Do(j, a) * coef(2) + Dn(j, a) * coef(3);
            }
        }
    }
  }

  // prestress
  if (com_mod.pstEq) {
    com_mod.pSn = 0.0;
    com_mod.pSa = 0.0;
  }
}
//------------------------
// picc
//------------------------
/// @brief Corrector with convergence check
///
/// Decision for next eqn is also made here (modifies cEq global).
///
/// Modifies:
/// \code {.cpp}
///   com_mod.Ad
///   An_ (member variable)
///   com_mod.Dn
///   com_mod.Yn
///   cep_mod.Xion
///   com_mod.pS0
///   com_mod.pSa
///   com_mod.pSn
///
///   com_mod.cEq
///   eq.FSILS.RI.iNorm
///   eq.pNorm
/// \endcode
//
void Integrator::picc()
{
  using namespace consts;

  auto& com_mod = simulation_->com_mod;
  auto& cep_mod = simulation_->get_cep_mod();

  #define n_debug_picc
  #ifdef debug_picc
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  const int nsd = com_mod.nsd;
  const int tnNo = com_mod.tnNo;
  const double dt = com_mod.dt;

  const auto& R = com_mod.R;
  const auto& Rd = com_mod.Rd;

  auto& cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];

  auto& An = An_;  // Use member variable
  auto& Ad = com_mod.Ad;
  auto& Dn = Dn_;
  auto& Yn = Yn_;

  auto& pS0 = com_mod.pS0;
  auto& pSa = com_mod.pSa;
  auto& pSn = com_mod.pSn;
  auto& Xion = cep_mod.Xion;

  int s = eq.s;
  int e = eq.e;

  std::array<double,4> coef;
  coef[0] = eq.gam * dt;
  coef[1] = eq.beta*dt*dt;
  coef[2] = 1.0 / eq.am;
  coef[3] = eq.af*coef[0]*coef[2];

  #ifdef debug_picc
  dmsg << "cEq: " << cEq;
  dmsg << "s: " << s;
  dmsg << "e: " << e;
  dmsg << "coef: " << coef[0] << " " << coef[1] << " " << coef[2] << " " << coef[3];
  dmsg << "sstEq: " << com_mod.sstEq;
  dmsg << "An nrows: " << An.nrows_;
  dmsg << "   ncols: " << An.ncols_;
  #endif

  // ustruct, FSI (ustruct)
  //
  if (com_mod.sstEq) {
    if (eq.phys == EquationType::phys_ustruct || eq.phys == EquationType::phys_FSI) {
      Vector<double> dUl(nsd);

      for (int a = 0; a < tnNo; a++) {
        for (int i = 0; i < e-s+1; i++) {
          An(i+s,a) = An(i+s,a) - R(i,a);
          Yn(i+s,a) = Yn(i+s,a) - R(i,a)*coef[0];
        }

        for (int i = 0; i < e-s; i++) {
          dUl(i) = Rd(i,a)*coef[2] + R(i,a)*coef[3];
          Ad(i,a) = Ad(i,a) - dUl(i);
          Dn(i+s,a) = Dn(i+s,a) - dUl(i)*coef[0];
        }
      }

    } else if (eq.phys == EquationType::phys_mesh) {
      for (int a = 0; a < tnNo; a++) {
        for (int i = 0; i < e-s+1; i++) {
          An(i+s,a) = An(i+s,a) - R(i,a);
          Yn(i+s,a) = Yn(i+s,a) - R(i,a)*coef[0];
          Dn(i+s,a) = Dn(i+s,a) - R(i,a)*coef[1];
        }
      }
    }

  } else {
    for (int a = 0; a < tnNo; a++) {
      for (int i = 0; i < e-s+1; i++) {
        // eqn 94 of Bazilevs 2007 // here, -R contains the acceleration update (obtained from Newton solve))?
        An(i+s,a) = An(i+s,a) - R(i,a);

        // eqn 95 of Bazilevs 2007
        Yn(i+s,a) = Yn(i+s,a) - R(i,a)*coef[0];

        Dn(i+s,a) = Dn(i+s,a) - R(i,a)*coef[1];
      }
    }
  }

  if (std::set<EquationType>{Equation_stokes, Equation_fluid, Equation_ustruct, Equation_FSI}.count(eq.phys) != 0) {
    pic_eth();
  }

  if (eq.phys == Equation_FSI) {
    int s = com_mod.eq[1].s;
    int e = com_mod.eq[1].e;
    #ifdef debug_picc
    dmsg << "eq.phys == Equation_FSI ";
    dmsg << "com_mod.eq[1].sym: " << com_mod.eq[1].sym;
    dmsg << "s: " << s;
    dmsg << "e: " << e;
    #endif

    for (int Ac = 0; Ac < tnNo; Ac++) {
      if (all_fun::is_domain(com_mod, eq, Ac, Equation_struct) ||
          all_fun::is_domain(com_mod, eq, Ac, Equation_ustruct) ||
          all_fun::is_domain(com_mod, eq, Ac, Equation_lElas)) {
        for (int i = 0; i < e-s+1; i++) {
          An(i+s,Ac) = An(i,Ac);
          Yn(i+s,Ac) = Yn(i,Ac);
          Dn(i+s,Ac) = Dn(i,Ac);
        }
      }
    }
  }

  // Update Xion for cardiac electrophysiology
  //
  if (eq.phys == Equation_CEP) {
    int s = eq.s;
    for (int a = 0; a < tnNo; a++) {
      Xion(0,a) = Yn(s,a);
    }
  }

  // Update prestress at the nodes and re-initialize
  //
  if (com_mod.pstEq) {
    all_fun::commu(com_mod, pSn);
    all_fun::commu(com_mod, pSa);

    for (int a = 0; a < tnNo; a++) {
      if (!utils::is_zero(pSa(a))) {
        for (int i = 0; i < pSn.nrows(); i++) {
          pSn(i,a) = pSn(i,a) / pSa(a);
        }
      }
    }

    pSa = 0.0;
  }

  // Filter out the non-wall displacements for CMM equation
  //
  if (eq.phys == Equation_CMM && !com_mod.cmmInit) {
    for (int a = 0; a < tnNo; a++) {
      double r1 = static_cast<double>(com_mod.cmmBdry(a));
      for (int i = 0; i < e-s; i++) {
        Dn(i+s,a) = Dn(i+s,a)*r1;
      }
    }
  }

  // IB treatment
  //if (ibFlag) CALL IB_PICC()

  // Computes norms and check for convergence of Newton iterations
  double eps = std::numeric_limits<double>::epsilon();

  if (utils::is_zero(eq.FSILS.RI.iNorm)) {
    eq.FSILS.RI.iNorm = eps;
  }

  if (utils::is_zero(eq.iNorm)) {
    eq.iNorm = eq.FSILS.RI.iNorm;
    #ifdef debug_picc
    dmsg << "eq.iNorm: " << eq.iNorm;
    #endif
  }

  if (eq.itr == 1) {
     eq.pNorm = eq.FSILS.RI.iNorm / eq.iNorm;
    #ifdef debug_picc
    dmsg << "eq.itr: " << eq.itr;
    dmsg << "eq.pNorm: " << eq.pNorm;
    #endif
  }

  double r1 = eq.FSILS.RI.iNorm / eq.iNorm;
  bool l1 = (eq.itr >= eq.maxItr);
  bool l2 = (r1 <= eq.tol);
  bool l3 = (r1 <= eq.tol*eq.pNorm);
  bool l4 = (eq.itr >= eq.minItr);

  #ifdef debug_picc
  dmsg << "eq.itr: " << eq.itr;
  dmsg << "eq.minItr: " << eq.minItr;
  dmsg << "r1: " << r1;
  dmsg << "l1: " << l1;
  dmsg << "l2: " << l2;
  dmsg << "l3: " << l3;
  dmsg << "l4: " << l4;
  #endif

  if (l1 || ((l2 || l3) && l4)) {
    eq.ok = true;
    #ifdef debug_picc
    dmsg << "eq.ok: " << eq.ok;
    dmsg << "com_mod.eq[0].ok: " << com_mod.eq[0].ok;
    dmsg << "com_mod.eq[1].ok: " << com_mod.eq[1].ok;
    #endif
  }

  auto& eqs = com_mod.eq;
  if (std::count_if(eqs.begin(),eqs.end(),[](eqType& eq){return eq.ok;}) == eqs.size()) {
    #ifdef debug_picc
    dmsg << "all ok";
    #endif
    return;
  }
  //if (ALL(eq.ok)) RETURN

  if (eq.coupled) {
    cEq = cEq + 1;
    #ifdef debug_picc
    dmsg << "eq " << " coupled ";
    dmsg << "1st update cEq: " << cEq;
    #endif

    auto& eqs = com_mod.eq;
    if (std::count_if(eqs.begin(),eqs.end(),[](eqType& eq){return !eq.coupled || eq.ok;}) == eqs.size()) {
      while (cEq < com_mod.nEq) {
        if (!eqs[cEq].coupled) {
          break;
        }
        cEq = cEq + 1;
      }

    } else {
      if (cEq >= com_mod.nEq) {
        cEq = 0;
      }

      while (!eqs[cEq].coupled) {
        cEq = cEq + 1;
        if (cEq >= com_mod.nEq) {
          cEq = 0;
        }
      }
    }

  } else {
    if (eq.ok) {
      cEq = cEq + 1;
    }
  }
 #ifdef debug_picc
 dmsg << "eq " << " coupled ";
 dmsg << "2nd update cEq: " << cEq;
 #endif
}

//------------------------
// pic_eth
//------------------------
/// @brief Pressure correction at edge nodes for Taylor-Hood type element
///
/// Here, we interpolate pressure at the edge nodes by interpolating
/// using a reduced basis (such as P1) applied on element vertices
/// (i.e., corner nodes). For e.g., for a P2 element, pressure is
/// interpolated at the edge nodes using P1 vertices.
///
/// Modifies: Yn_
///
void Integrator::pic_eth()
{
  using namespace consts;

  auto& com_mod = simulation_->com_mod;
  auto& cep_mod = simulation_->get_cep_mod();

  const int nsd = com_mod.nsd;
  const int tnNo = com_mod.tnNo;
  const double dt = com_mod.dt;
  const auto& cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];

  auto& cDmn = com_mod.cDmn;
  auto& Yn = Yn_;

  // Check for something ...
  //
  bool THflag = false;

  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    if (com_mod.msh[iM].nFs == 2) {
      THflag = true;
      break;
     }
  }

  if (!THflag) {
    return;
  }

  Vector<double> sA(tnNo), sF(tnNo);
  int s = eq.s;

  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    auto& msh = com_mod.msh[iM];
    if (msh.nFs == 1) {
      continue;
    }

    auto eType = msh.fs[1].eType;
    int eNoN = msh.fs[0].eNoN;
    int eNoNq = msh.fs[1].eNoN;

    Array<double> xl(nsd,eNoN), xql(nsd,eNoNq), Nqx(nsd,eNoNq);
    Vector<double> pl(eNoNq), Nq(eNoNq);

    Vector<double> xp(nsd), xi0(nsd), xi(nsd);
    Array<double> ksix(nsd,nsd);

    for (int g = 0; g < msh.fs[1].nG; g++) {
      for (int i = 0; i < nsd; i++) {
        xi0(i) = xi0(i) + msh.fs[1].xi(i,g);
      }
    }

    xi0 = xi0 / static_cast<double>(msh.fs[1].nG);

    for (int e = 0; e < msh.nEl; e++) {
      cDmn = all_fun::domain(com_mod, msh, cEq, e);       // setting global cDmn
      if ((eq.dmn[cDmn].phys != Equation_stokes) &&
          (eq.dmn[cDmn].phys != Equation_fluid)  &&
          (eq.dmn[cDmn].phys != Equation_ustruct)) {
        continue;
      }

      for (int a = 0; a < eNoN; a++) {
        int Ac = msh.IEN(a,e);
        for (int i = 0; i < nsd; i++) {
          xl(i,a) = com_mod.x(i,Ac);
        }
      }

      for (int a = 0; a < eNoNq; a++) {
        int Ac = msh.IEN(a,e);
        pl(a) = Yn(s+nsd,Ac);
        for (int i = 0; i < nsd; i++) {
          xql(i,a) = xl(i,a);
        }
      }

      double eVol = 0.0;
      double Jac = 0.0;

      for (int g = 0; g < msh.fs[1].nG; g++) {
        if (g == 0 || !msh.fs[1].lShpF) {
          auto Nx = msh.fs[1].Nx.slice(g);
          nn::gnn(eNoNq, nsd, nsd, Nx, xql, Nqx, Jac, ksix);

          if (utils::is_zero(Jac)) {
            throw std::runtime_error("[pic_eth] Jacobian for element " + std::to_string(e) + " is < 0.");
          }
        }

        eVol = eVol + msh.fs[1].w(g)*Jac;
      }

      for (int a = eNoNq; a < eNoN; a++) {
        int Ac = msh.IEN(a,e);
        for (int i = 0; i < nsd; i++) {
          xp(i) = xl(i,a);
        }
        xi = xi0;
        nn::get_nnx(nsd, eType, eNoNq, xql, msh.fs[1].xib, msh.fs[1].Nb, xp, xi, Nq, Nqx);

        double p = 0.0;
        for (int b = 0; b < eNoNq; b++) {
          p = p + pl(b)*Nq(b);
        }

        sF(Ac) = sF(Ac) + p*eVol;
        sA(Ac) = sA(Ac) + eVol;
      }
    } // e-loop
  } // iM-loop

  all_fun::commu(com_mod, sA);
  all_fun::commu(com_mod, sF);

  for (int a = 0; a < tnNo; a++) {
    if (!utils::is_zero(sA(a))) {
      Yn(s+nsd,a) = sF(a) / sA(a);
    }
  }
}
