// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// The functions defined here are used to run a simulation from the command line.
//
// Usage:
//
//   svMultiPhysics XML_FILE_NAME
//
#include "Simulation.h"

#include "all_fun.h"
#include "bf.h"
#include "contact.h"
#include "distribute.h"
#include "eq_assem.h"
#include "fs.h"
#include "initialize.h"
#include "ls.h"
#include "output.h"
#include "pic.h"
#include "read_files.h"
#include "read_msh.h"
#include "remesh.h"
#include "set_bc.h"
#include "txt.h"
#include "ustruct.h"
#include "vtk_xml.h"
#include "ris.h"
#include "uris.h"

#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <fstream>

//------------------------
// add_eq_linear_algebra
//------------------------
// Create a LinearAlgebra object for an equation.
//
void add_eq_linear_algebra(ComMod& com_mod, eqType& lEq)
{
  lEq.linear_algebra = LinearAlgebraFactory::create_interface(lEq.linear_algebra_type);
  lEq.linear_algebra->set_preconditioner(lEq.linear_algebra_preconditioner);
  lEq.linear_algebra->initialize(com_mod, lEq);

  if (lEq.linear_algebra_assembly_type != consts::LinearAlgebraType::none) {
    lEq.linear_algebra->set_assembly(lEq.linear_algebra_assembly_type);
  }
}

void finalize_linear_algebra(eqType& lEq)
{
  lEq.linear_algebra->finalize();
}

/// @brief Read in a solver XML file and all mesh and BC data.  
//
void read_files(Simulation* simulation, const std::string& file_name)
{
  simulation->com_mod.timer.set_time();

  if (simulation->com_mod.cm.slv(simulation->cm_mod)) {
    return;
  }

  read_files_ns::read_files(simulation, file_name);

/*
  try {
    read_files_ns::read_files(simulation, file_name);

  } catch (const std::exception& exception) {
    std::cout << "[svMultiPhysics] ERROR: The svMultiPhysics program has failed." << std::endl;
    std::cout << "[svMultiPhysics] ERROR: " << exception.what() << std::endl;
    exit(1);
  }
*/
  
}



/// @brief Iterate the precomputed state-variables in time using linear interpolation to the current time step size
//
void iterate_precomputed_time(Simulation* simulation) {
  using namespace consts;

  auto& com_mod = simulation->com_mod;
  auto& cm_mod = simulation->cm_mod;
  auto& cm = com_mod.cm;
  auto& cep_mod = simulation->get_cep_mod();

  int nTS = com_mod.nTS;
  int stopTS = nTS;
  int tDof = com_mod.tDof;
  int tnNo = com_mod.tnNo;
  int nFacesLS = com_mod.nFacesLS;
  int nsd = com_mod.nsd;

  auto& Ad = com_mod.Ad;      // Time derivative of displacement
  auto& Rd = com_mod.Rd;      // Residual of the displacement equation
  auto& Kd = com_mod.Kd;      // LHS matrix for displacement equation

  auto& Ao = com_mod.Ao;      // Old time derivative of variables (acceleration)
  auto& Yo = com_mod.Yo;      // Old variables (velocity)
  auto& Do = com_mod.Do;      // Old integrated variables (dissplacement)

  auto& An = com_mod.An;      // New time derivative of variables
  auto& Yn = com_mod.Yn;      // New variables (velocity)
  auto& Dn = com_mod.Dn;      // New integrated variables

  int& cTS = com_mod.cTS;
  int& nITs = com_mod.nITs;
  double& dt = com_mod.dt;

  if (com_mod.usePrecomp) {
#ifdef debug_iterate_solution
        dmsg << "Use precomputed values ..." << std::endl;
#endif
    // This loop is used to interpolate between known time values of the precomputed
    // state-variable solution
    for (int l = 0; l < com_mod.nMsh; l++) {
      auto lM = com_mod.msh[l];
      if (lM.Ys.nslices() > 1) {
        // If there is only one temporal slice, then the solution is assumed constant
        // in time and no interpolation is performed
        // If there are multiple temporal slices, then the solution is linearly interpolated
        // between the known time values and the current time.
        double precompDt = com_mod.precompDt;
        double preTT = precompDt * (lM.Ys.nslices() - 1);
        double cT = cTS * dt;
        double rT = std::fmod(cT, preTT);
        int n1, n2;
        double alpha;
        if (precompDt == dt) {
          alpha =  0.0;
          if (cTS < lM.Ys.nslices()) {
            n1 = cTS - 1;
          } else {
            n1 = cTS % lM.Ys.nslices() - 1;
          }
        } else {
          n1 = static_cast<int>(rT / precompDt) - 1;
          alpha = std::fmod(rT, precompDt);
        }
        n2 = n1 + 1;
        for (int i = 0; i < tnNo; i++) {
          for (int j = 0; j < nsd; j++) {
            if (alpha == 0.0) {
              Yn(j, i) = lM.Ys(j, i, n2);
            } else {
              Yn(j, i) = (1.0 - alpha) * lM.Ys(j, i, n1) + alpha * lM.Ys(j, i, n2);
            }
          }
        }
      } else {
        for (int i = 0; i < tnNo; i++) {
          for (int j = 0; j < nsd; j++) {
            Yn(j, i) = lM.Ys(j, i, 0);
          }
        }
      }
    }
  }
}

/// @brief Iterate the simulation in time.
///
/// Reproduces the outer and inner loops in Fortan MAIN.f.
//
void iterate_solution(Simulation* simulation)
{
  using namespace consts;

  auto& com_mod = simulation->com_mod;
  auto& cm_mod = simulation->cm_mod;
  auto& cm = com_mod.cm;
  auto& cep_mod = simulation->get_cep_mod();
  
  // number of time steps
  int nTS = com_mod.nTS;
  
  int stopTS = nTS;
  
  // total number of degrees of freedom per node
  int tDof = com_mod.tDof;
  
  // total number of nodes across all meshes, but only on current processor
  int tnNo = com_mod.tnNo;
  
  int nFacesLS = com_mod.nFacesLS;
  int nsd = com_mod.nsd;

  std::cout << std::scientific << std::setprecision(16);

  #define n_debug_iterate_solution
  #ifdef debug_iterate_solution
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  #ifdef debug_iterate_solution
  dmsg << "========== iterate_solution ==========" << std::endl;
  dmsg << "tDof: " << tDof;
  dmsg << "tnNo: " << tnNo;
  dmsg << "nFacesLS: " << nFacesLS;
  dmsg << "stopTS: " << stopTS;
  dmsg << "cmmInit: " << com_mod.cmmInit;
  #endif

  Array<double> Ag(tDof,tnNo); 
  Array<double> Yg(tDof,tnNo); 
  Array<double> Dg(tDof,tnNo); 
  Vector<double> res(nFacesLS); 
  Vector<int> incL(nFacesLS);

  // current time step
  int& cTS = com_mod.cTS;
  
  int& nITs = com_mod.nITs;
  
  // time step size
  double& dt = com_mod.dt;
  
  #ifdef debug_iterate_solution
  dmsg;
  dmsg << "cTS: " << cTS;
  dmsg << "nITs: " << nITs;
  dmsg << "dt: " << dt;
  #endif

  if (cTS <= nITs) { 
    dt = dt / 10.0;
  }

  double& time = com_mod.time;
  auto& cEq = com_mod.cEq;

  auto& Ad = com_mod.Ad;      // Time derivative of displacement 
  auto& Rd = com_mod.Rd;      // Residual of the displacement equation
  auto& Kd = com_mod.Kd;      // LHS matrix for displacement equation

  auto& Ao = com_mod.Ao;      // Old time derivative of variables (acceleration)
  auto& Yo = com_mod.Yo;      // Old variables (velocity)
  auto& Do = com_mod.Do;      // Old integrated variables (displacement)

  auto& An = com_mod.An;      // New time derivative of variables (acceleration)
  auto& Yn = com_mod.Yn;      // New variables (velocity)
  auto& Dn = com_mod.Dn;      // New integrated variables (displacement)

  bool l1 = false;
  bool l2 = false;
  bool l3 = false;

  #ifdef debug_iterate_solution
  dmsg << "Start Outer Loop ..." << std::endl;
  #endif

  bool exit_now = false;
  double elapsed_time = 0.0;

  // Uncomment these two lines to enable writting values to a file.
  //Array<double>::write_enabled = true;
  //Array3<double>::write_enabled = true;

  // Outer loop for marching in time. When entering this loop, all old
  // variables are completely set and satisfy BCs.
  // 
  while (true) {
    #ifdef debug_iterate_solution
    dmsg << "========================================= " << std::endl;
    dmsg << "=============== Outer Loop ============== " << std::endl;
    dmsg << "========================================= " << std::endl;
    #endif

    // Adjusting the time step size once initialization stage is over
    //
    if (cTS == nITs) {
      dt = 10.0 * dt;
      #ifdef debug_iterate_solution
      dmsg << "New time step size (dt): " << dt;
      #endif
    }

    // Incrementing time step, hence cTS will be associated with new
    // variables, i.e. An, Yn, and Dn
    //
    cTS = cTS + 1;
    time = time + dt;
    cEq = 0;
    std::string cstr = "_cts_" + std::to_string(cTS);
    #ifdef debug_iterate_solution
    dmsg << "nITs: " << nITs;
    dmsg << "cTS: " << cTS;
    dmsg << "dt: " << dt;
    dmsg << "time: " << time;
    dmsg << "mvMsh: " << com_mod.mvMsh;
    dmsg << "rmsh.isReqd: " << com_mod.rmsh.isReqd;
    #endif

    for (auto& eq : com_mod.eq) {
      eq.itr = 0;
      eq.ok = false;
    }

    // Compute mesh properties to check if remeshing is required
    //
    if (com_mod.mvMsh && com_mod.rmsh.isReqd) {
      read_msh_ns::calc_mesh_props(com_mod, cm_mod, com_mod.nMsh, com_mod.msh);
      if (com_mod.resetSim) {
        #ifdef debug_iterate_solution
        dmsg << "#### resetSim is true " << std::endl;
        dmsg << "#### Breaking out from Outer Loop " << std::endl;
        #endif
        break;
      }
    }

    // Predictor step
    #ifdef debug_iterate_solution
    dmsg << "Predictor step ... " << std::endl;
    #endif
    pic::picp(simulation);

    // Apply Dirichlet BCs strongly
    //
    // Modifes
    //  An - New time derivative of variables
    //  Yn - New variables
    //  Dn -  New integrated variables
    //  com_mod.Ad - Time derivative of displacement
    //
    #ifdef debug_iterate_solution
    dmsg << "Apply Dirichlet BCs strongly ..." << std::endl;
    #endif

    set_bc::set_bc_dir(com_mod, An, Yn, Dn);

    if (com_mod.urisFlag) {uris::uris_calc_sdf(com_mod);}

    iterate_precomputed_time(simulation);

    // Inner loop for Newton iteration
    //
    int inner_count = 1;
    int reply;
    int iEqOld;

    // Looping over Newton iterations
    while (true) { 
      #ifdef debug_iterate_solution
      dmsg << "---------- Inner Loop " + std::to_string(inner_count) << " -----------" << std::endl;
      dmsg << "cEq: " << cEq;
      dmsg << "com_mod.eq[cEq].sym: " << com_mod.eq[cEq].sym;
      //simulation->com_mod.timer.set_time();
      #endif
      //std::cout << "-------------------------------------" << std::endl;
      //std::cout << "inner_count: " << inner_count << std::endl;

      auto istr = "_" + std::to_string(cTS) + "_" + std::to_string(inner_count);
      iEqOld = cEq;
      auto& eq = com_mod.eq[cEq];

      if (com_mod.cplBC.coupled && cEq == 0) {
        #ifdef debug_iterate_solution
        dmsg << "Set coupled BCs " << std::endl;
        #endif
        set_bc::set_bc_cpl(com_mod, cm_mod);

        set_bc::set_bc_dir(com_mod, An, Yn, Dn);
      }

      // Initiator step for Generalized α− Method (quantities at n+am, n+af). 
      //
      // Modifes
      //   Ag((tDof, tnNo) - 
      //   Yg((tDof, tnNo) - 
      //   Dg((tDof, tnNo) - 
      //
      #ifdef debug_iterate_solution
      dmsg << "Initiator step ..." << std::endl;
      #endif
      pic::pici(simulation, Ag, Yg, Dg);
      Ag.write("Ag_pic"+ istr);
      Yg.write("Yg_pic"+ istr);
      Dg.write("Dg_pic"+ istr);
      Yn.write("Yn_pic"+ istr);

      if (Rd.size() != 0) {
        Rd = 0.0;
        Kd = 0.0;
      }

      // Allocate com_mod.R and com_mod.Val arrays.
      //
      // com_mod.R(dof,tnNo)
      // com_mod.Val(dof*dof, com_mod.lhs.nnz)
      //
      // If Trilinos is used then allocate
      //   com_mod.tls.W(dof,tnNo)
      //   com_mod.tls.R(dof,tnNo)
      //
      #ifdef debug_iterate_solution
      dmsg << "Allocating the RHS and LHS"  << std::endl;
      #endif

      ls_ns::ls_alloc(com_mod, eq);
      com_mod.Val.write("Val_alloc"+ istr);

      // Compute body forces. If phys is shells or CMM (init), apply
      // contribution from body forces (pressure) to residual
      //
      // Modifes: com_mod.Bf, Dg
      //
      #ifdef debug_iterate_solution
      dmsg << "Set body forces ..."  << std::endl;
      #endif

      bf::set_bf(com_mod, Dg);
      com_mod.Val.write("Val_bf"+ istr);

      // Assemble equations.
      //
      #ifdef debug_iterate_solution
      dmsg << "Assembling equation:  " << eq.sym;
      #endif

      for (int iM = 0; iM < com_mod.nMsh; iM++) {
        eq_assem::global_eq_assem(com_mod, cep_mod, com_mod.msh[iM], Ag, Yg, Dg);
      }
      com_mod.R.write("R_as"+ istr);
      com_mod.Val.write("Val_as"+ istr);

      // Treatment of boundary conditions on faces
      // Apply Neumman or Traction boundary conditions
      //
      // Modifies: com_mod.R
      //
      #ifdef debug_iterate_solution
      dmsg << "Apply Neumman or Traction BCs ... " << std::endl;
      #endif

      Yg.write("Yg_vor_neu"+ istr);
      Dg.write("Dg_vor_neu"+ istr);

      set_bc::set_bc_neu(com_mod, cm_mod, Yg, Dg);

      com_mod.Val.write("Val_neu"+ istr);
      com_mod.R.write("R_neu"+ istr);
      Yg.write("Yg_neu"+ istr);
      Dg.write("Dg_neu"+ istr);

      // Apply CMM BC conditions
      //
      if (!com_mod.cmmInit) {
        #ifdef debug_iterate_solution
        dmsg << "Apply CMM BC conditions ... " << std::endl;
        #endif
        set_bc::set_bc_cmm(com_mod, cm_mod, Ag, Dg);
      }

      // Apply weakly applied Dirichlet BCs
      //
      #ifdef debug_iterate_solution
      dmsg << "Apply weakly applied Dirichlet BCs ... " << std::endl;
      #endif

      set_bc::set_bc_dir_w(com_mod, Yg, Dg);

      if (com_mod.risFlag) {
        ris::ris_resbc(com_mod, Yg, Dg);
      }

      if (com_mod.ris0DFlag) {
        ris::ris0d_bc(com_mod, cm_mod, Yg, Dg);
      }

      // Apply contact model and add its contribution to residual
      //
      if (com_mod.iCntct) {
        contact::construct_contact_pnlty(com_mod, cm_mod, Dg);

#if 0
        if (cTS <= 2050) {
          Array<double>::write_enabled = true;
          com_mod.R.write("R_"+ std::to_string(cTS));
          //exit(0);
        }
#endif
      }

      // Synchronize R across processes. Note: that it is important
      // to synchronize residual, R before treating immersed bodies as
      // ib.R is already communicated across processes
      //
      if (!eq.assmTLS) {
        #ifdef debug_iterate_solution
        dmsg << "Synchronize R across processes ..." << std::endl;
        #endif
        all_fun::commu(com_mod, com_mod.R);
      }

      // Update residual in displacement equation for USTRUCT phys.
      // Note that this step is done only first iteration. Residual
      // will be 0 for subsequent iterations
      //
      // Modifies com_mod.Rd.
      //
      #ifdef debug_iterate_solution
      dmsg << "com_mod.sstEq: " << com_mod.sstEq;
      #endif
      if (com_mod.sstEq) {
        ustruct::ustruct_r(com_mod, Yg);
      }

      // Set the residual of the continuity equation and its tangent matrix
      // due to variation with pressure to 0 on all the edge nodes.
      //
      if (std::set<EquationType>{Equation_stokes, Equation_fluid, Equation_ustruct, Equation_FSI}.count(eq.phys) != 0) {
        #ifdef debug_iterate_solution
        dmsg << "thood_val_rc ..." << std::endl;
        #endif
        fs::thood_val_rc(com_mod);
      }

      // Treat Neumann boundaries that are not deforming
      //
      #ifdef debug_iterate_solution
      dmsg << "set_bc_undef_neu ..." << std::endl;
      #endif

      set_bc::set_bc_undef_neu(com_mod);

      // IB treatment: for explicit coupling, simply construct residual.
      //
      /* [NOTE] not implemented.
      if (com_mod.ibFlag) {
        if (com_mod.ib.cpld == ibCpld_I) {
          //CALL IB_IMPLICIT(Ag, Yg, Dg)
        }
        // CALL IB_CONSTRUCT()
      }
      */

      #ifdef debug_iterate_solution
      dmsg << "Update res() and incL ..." << std::endl;
      dmsg << "nFacesLS: " << nFacesLS;
      #endif
      incL = 0;
      if (eq.phys == Equation_mesh) {
        incL(nFacesLS-1) = 1;
      }

      if (com_mod.cmmInit) {
        incL(nFacesLS-1) = 1;
      }

      for (int iBc = 0; iBc < eq.nBc; iBc++) {
        int i = eq.bc[iBc].lsPtr;
        if (i != -1) {
          // Resistance term for coupled Neumann BC tangent contribution
          res(i) = eq.gam * dt * eq.bc[iBc].r;
          incL(i) = 1;
        }
      }

      // Solve equation.
      //
      // Modifies: com_mod.R, com_mod.Val 
      //
      #ifdef debug_iterate_solution
      dmsg << "Solving equation: " << eq.sym; 
      #endif

      ls_ns::ls_solve(com_mod, eq, incL, res);

      com_mod.Val.write("Val_solve"+ istr);
      com_mod.R.write("R_solve"+ istr);

      // Solution is obtained, now updating (Corrector) and check for convergence
      //
      // Modifies: com_mod.An com_mod.Dn com_mod.Yn cep_mod.Xion com_mod.pS0 com_mod.pSa
      //            com_mod.pSn com_mod.cEq eq.FSILS.RI.iNorm eq.pNorm 
      //
      #ifdef debug_iterate_solution
      dmsg << "Update corrector ..." << std::endl; 
      #endif

      pic::picc(simulation);
      com_mod.Yn.write("Yn_picc"+ istr);

      // Writing out the time passed, residual, and etc.
      if (std::count_if(com_mod.eq.begin(),com_mod.eq.end(),[](eqType& eq){return eq.ok;}) == com_mod.eq.size()) { 
        #ifdef debug_iterate_solution
        dmsg << ">>> All OK" << std::endl; 
        dmsg << "iEqOld: " << iEqOld+1; 
        #endif
        break;
      } 
      output::output_result(simulation, com_mod.timeP, 2, iEqOld);

      inner_count += 1;
    } // Inner loop

    #ifdef debug_iterate_solution
    dmsg << ">>> End of inner loop " << std::endl; 
    #endif

    // IB treatment: interpolate flow data on IB mesh from background
    // fluid mesh for explicit coupling, update old solution for implicit
    // coupling
    //
    /* [NOTE] Not implemented.
    if (ibFlag) {
      CALL IB_INTERPYU(Yn, Dn)
      if (ib.cpld == ibCpld_I) {
        ib.Auo = ib.Aun
        ib.Ubo = ib.Ubn
      }
    }
    */

    if (com_mod.risFlag) {
      ris::ris_meanq(com_mod, cm_mod);
      ris::ris_status(com_mod, cm_mod);
      if (cm.mas(cm_mod)) {
        std::cout << "Iteration: " << com_mod.cTS << std::endl;
        for (int iProj = 0; iProj < com_mod.ris.nbrRIS; iProj++) {
          std::cout << "Status for RIS projection: " << iProj << std::endl;
          std::cout << "            RIS iteration: " << com_mod.ris.nbrIter(iProj) << std::endl;
          std::cout << "       Is the valve close? " << com_mod.ris.clsFlg[iProj] << std::endl;
          std::cout << "            The status is: " << com_mod.ris.status[iProj] << std::endl;
        }
      }

      if (!std::all_of(com_mod.ris.status.begin(), com_mod.ris.status.end(), [](bool s) { return s; })) {
        if (std::any_of(com_mod.ris.nbrIter.begin(), com_mod.ris.nbrIter.end(), [](int iter) { return iter <= 1; })) {
          if (cm.mas(cm_mod)) {
            std::cout << "Valve status just changed. Do not update" << std::endl;
          }
        } else {
            ris::ris_updater(com_mod, cm_mod);
        }
        // goto label_11;
      }
    }

    // Saving the TXT files containing average and fluxes (or ECGs)
    #ifdef debug_iterate_solution
    dmsg << "Saving the TXT files containing average and fluxes ..." << std::endl;
    dmsg << "Saving the TXT files containing ECGs ..." << std::endl;
    #endif

    txt_ns::txt(simulation, false);

    // If remeshing is required then save current solution.
    //
    if (com_mod.rmsh.isReqd) {
      l1 = ((cTS % com_mod.rmsh.cpVar) == 0);
      if (l1) {
        #ifdef debug_iterate_solution
        dmsg << "Saving last solution for remeshing." << std::endl; 
        #endif
        com_mod.rmsh.rTS = cTS - 1;
        com_mod.rmsh.time = time - dt;
        for (int i = 0; i < com_mod.rmsh.iNorm.size(); i++) {
          com_mod.rmsh.iNorm(i) = com_mod.eq[i].iNorm;
        }

        com_mod.rmsh.A0 = com_mod.Ao;
        com_mod.rmsh.Y0 = com_mod.Yo;
        com_mod.rmsh.D0 = com_mod.Do;
      }
    }

    // Look for a file containg a time step to stop the simulation.
    //
    // stopTrigName = "STOP_SIM"
    //
    auto& stopTrigName = com_mod.stopTrigName;
    bool l1 = false;
    int stopTS = 0;
    int count = -1;

    if (cm.mas(cm_mod)) {
      if (FILE *fp = fopen(stopTrigName.c_str(), "r")) {
        l1 = true;
        count = fscanf(fp, "%d", &stopTS);

        if (count == 0) {
          stopTS = cTS;
        }
        fclose(fp);

      } else {
        stopTS = nTS;
      }
    }

    #ifdef debug_iterate_solution
    dmsg << "cm.bcast(cm_mod, &stopTS)  ..." << std::endl; 
    #endif

    cm.bcast(cm_mod, &stopTS);

    l1 = (cTS >= stopTS);
    l2 = ((cTS % com_mod.stFileIncr) == 0);

    #ifdef debug_iterate_solution
    dmsg; 
    dmsg << "stFileIncr: " << com_mod.stFileIncr; 
    dmsg << "l1: " << l1; 
    dmsg << "l2: " << l2; 
    #endif

    // Saving the result to restart bin file
    if (l1 || l2) {
       output::write_restart(simulation, com_mod.timeP);
    }

    // Writing results into the disk with VTU format
    //
    #ifdef debug_iterate_solution
    dmsg; 
    dmsg << "saveVTK: " << com_mod.saveVTK; 
    #endif

    if (com_mod.saveVTK) {
      l2 = ((cTS % com_mod.saveIncr) == 0);
      l3 = (cTS >= com_mod.saveATS);
      #ifdef debug_iterate_solution
      dmsg << "l2: " << l2; 
      dmsg << "l3: " << l3; 
      #endif

      if (l2 && l3) {
        output::output_result(simulation, com_mod.timeP, 3, iEqOld);
        bool lAvg = false;
        vtk_xml::write_vtus(simulation, An, Yn, Dn, lAvg);
      } else {
        output::output_result(simulation, com_mod.timeP, 2, iEqOld);
      }

    } else {
      output::output_result(simulation, com_mod.timeP, 2, iEqOld);
    }

    // [NOTE] Not implemented.
    //
    if (com_mod.pstEq) {
      //CALL OUTDNORM()
    }

    if (com_mod.ibFlag) {
      //CALL IB_OUTCPUT()
    }

    // [HZ] Part related to RIS0D
    if (cEq == 0 && com_mod.ris0DFlag) {
      ris::ris0d_status(com_mod, cm_mod);
    }

    // [HZ] Part related to unfitted RIS
    // If the valve is active, look at the pressure difference 
    if (com_mod.urisFlag) {
      for (int iUris = 0; iUris < com_mod.nUris; iUris++) {
        com_mod.uris[iUris].cnt++;
        if (com_mod.uris[iUris].clsFlg) {
          uris::uris_meanp(com_mod, cm_mod, iUris);
          // if (com_mod.uris[iUris].cnt == 1) {
          //   // GOTO 11 // The GOTO Statement in the Fortran code
          // }
        } else {
          uris::uris_meanv(com_mod, cm_mod, iUris);
        }
        if (cm.mas(cm_mod)) {
          std::cout << " URIS surface: " << com_mod.uris[iUris].name << ", count: " << com_mod.uris[iUris].cnt << std::endl;
        }
      }

      if (com_mod.mvMsh) {
        uris::uris_update_disp(com_mod, cm_mod);
      }

      if (cm.mas(cm_mod)) {
        if (l2 && l3) {
          uris::uris_write_vtus(com_mod);
        }
      }
    }
    // end RIS/URIS stuff 

    // Exiting outer loop if l1
    if (l1) {
      break;
    }

    // Solution is stored here before replacing it at next time step
    //
    Ao = An;
    Yo = Yn;

    if (com_mod.dFlag) {
      Do = Dn;
    }
    com_mod.cplBC.xo = com_mod.cplBC.xn;

  } // End of outer loop

  #ifdef debug_iterate_solution
  dmsg << "End of outer loop" << std::endl;
  #endif

  //#ifdef debug_iterate_solution
  //dmsg << "=======  Simulation Finished   ========== " << std::endl;
  //#endif
}


void run_simulation(Simulation* simulation)
{
  iterate_solution(simulation);
}


/// @brief Run a simulation from the command line using the name of a solver input 
/// XML file as an argument.
//
int main(int argc, char *argv[])
{
  if (argc < 2) {
    std::cout << "[svMultiPhysics] ERROR: The svMultiPhysics program requires the solver input XML file name as an argument." << std::endl;
    exit(1);
  }

  // Process extra arguments for XML parameter substitution.
  for (int i = 2; i < argc; i++) {
    std::string str(argv[i]);
    int pos = str.find("=");
    auto name = str.substr(0,pos);
    auto value = str.substr(pos+1,str.size());
  }

  std::cout << std::scientific << std::setprecision(16);

  // Initialize MPI.
  //
  int mpi_rank, mpi_size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

#ifdef ENABLE_ARRAY_INDEX_CHECKING
  if (mpi_rank == 0) {
    std::cout << "WARNING: Index checking is enabled" << std::endl;
  }
#endif

  // Create a Simulation object that stores all data structures for a simulation.
  //
  // The MPI prociess rank is set in the cmType::new_cm() method called
  // from the Simulation constructor. 
  //
  auto simulation = new Simulation();
  auto& cm = simulation->com_mod.cm;
  std::string file_name(argv[1]);

  #define n_debug_main
  #ifdef debug_main
  DebugMsg dmsg(__func__, cm.idcm());
  dmsg.banner();
  #endif

  // Iterate for restarting a simulation after remeshing. 
  //
  while (true) {

    // Read in the solver commands .xml file.
    //
    #ifdef debug_main
    dmsg << "Read files " << " ... ";
    #endif
    read_files(simulation, file_name);
    
    // Distribute data to processors.
    #ifdef debug_main
    dmsg << "Distribute data to processors " << " ... ";
    #endif
    distribute(simulation);

    // Initialize simulation data.
    //
    Vector<double> init_time(3);

    #ifdef debug_main
    dmsg << "Initialize " << " ... ";
    #endif
    initialize(simulation, init_time);

    // Create LinearAlgebra objects for each equation.
    //
    for (int iEq = 0; iEq < simulation->com_mod.nEq; iEq++) {
      auto& eq = simulation->com_mod.eq[iEq];
      add_eq_linear_algebra(simulation->com_mod, eq);
    }

    #ifdef debug_main
    for (int iM = 0; iM < simulation->com_mod.nMsh; iM++) {
      dmsg << "---------- iM " << iM;
      dmsg << "msh[iM].nNo: " << simulation->com_mod.msh[iM].nNo;
      dmsg << "msh[iM].gnNo: " << simulation->com_mod.msh[iM].gnNo;
      dmsg << "msh[iM].nEl: " << simulation->com_mod.msh[iM].nEl;
      dmsg << "msh[iM].gnEl: " << simulation->com_mod.msh[iM].gnEl;
    }
    #endif

    // Run the simulation.
    run_simulation(simulation);

    #ifdef debug_main
    dmsg << "resetSim: " << simulation->com_mod.resetSim;
    #endif

    // Remesh and continue the simulation.
    //
    if (simulation->com_mod.resetSim) {
      #ifdef debug_main
      dmsg << "Calling remesh_restart" << " ..."; 
      #endif
      remesh::remesh_restart(simulation);
      #ifdef debug_main
      dmsg << "Continue the simulation " << " ";
      #endif
    } else {
      break;
    }
  }

   for (int iEq = 0; iEq < simulation->com_mod.nEq; iEq++) {
      auto& eq = simulation->com_mod.eq[iEq];
      finalize_linear_algebra(eq);
    }

  MPI_Finalize();
}
