// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// The functions defined here replicate the Fortran functions defined in READFILES.f.

#include "read_files.h"

#include "all_fun.h"
#include "consts.h"
#include "read_msh.h"
#include "fft.h"
#include "vtk_xml.h"

#include "Array.h"
#include "CepMod.h"
#include "VtkData.h"

#include "fsils_api.hpp"
#include "fils_struct.hpp"

#include <fstream>
#include <functional>
#include <math.h>
#include <sstream>
#include <vector>

namespace read_files_ns {

#include "set_equation_props.h"
#include "set_material_props.h"
//#include "set_output_props.h"
#include "set_viscosity_props.h"

/// @brief Match two faces?
///
/// \todo [TODO:DaveP] this has not been tested.
//
void face_match(ComMod& com_mod, faceType& lFa, faceType& gFa, Vector<int>& ptr)
{
  using namespace read_msh_ns;
  int nsd = com_mod.nsd;
  int nBlkd = round(pow(static_cast<double>(gFa.nNo)/1000.0, 0.333));
  if (nBlkd == 0) {
    nBlkd = 1;
  }

  int nBlk = pow(nBlkd, nsd);
  Vector<double> xMin(nsd), xMax(nsd); 
  double eps = std::numeric_limits<double>::epsilon();

  for (int i = 0;  i < nsd; i++) {
    xMin[i] = std::min(lFa.x.min(), gFa.x.min());
    xMax[i] = std::max(lFa.x.max(), gFa.x.max());

    if (xMin[i] <  0.0) {
      xMin[i] = xMin[i] * (1.0 + eps);
    } else { 
      xMin[i] = xMin[i] * (1.0 - eps);
    }

    if (xMax[i] < 0.0) {
      xMax[i] = xMax[i] * (1.0 - eps);
    } else { 
       xMax[i] = xMax[i] * (1.0 + eps);
    }
  }

  auto dx = (xMax - xMin) / static_cast<double>(nBlkd);
  std::vector<bool> nFlt = {true, true, true};

  for (int i = 0; i < nsd; i++) { 
     if (utils::is_zero(dx[i])) {
       nFlt[i] = false;
     }
  }

  std::vector<blkType> blk(nBlk);
  std::vector<int> nodeBlk(gFa.nNo); 

  for (int a = 0; a < gFa.nNo; a++) {
    auto coord = gFa.x.col(a);
    int iBlk = find_blk(nsd, nBlkd, nFlt, xMin, dx, coord);
    nodeBlk[a] = iBlk;
    blk[iBlk].n = blk[iBlk].n + 1;
  }

  for (int iBlk = 0; iBlk < nBlk; iBlk++) {
    blk[iBlk].gN = Vector<int>(blk[iBlk].n);
    blk[iBlk].n = 0;
  }

  for (int a = 0; a < gFa.nNo; a++) {
    int iBlk = nodeBlk[a];
    blk[iBlk].gN[blk[iBlk].n] = a;
    blk[iBlk].n = blk[iBlk].n + 1;
  }

  for (int a = 0; a < lFa.nNo; a++) {
    auto coord = lFa.x.col(a);
    int iBlk = find_blk(nsd, nBlkd, nFlt, xMin, dx, coord);
    auto minS = std::numeric_limits<double>::max();

    for (int i = 0;  i < blk[iBlk].n; i++) {
      int b = blk[iBlk].gN[i];
      auto diff = lFa.x.col(a) - gFa.x.col(b);
      double ds = sqrt(diff*diff);

      if (ds < minS) {
        minS = ds;
        ptr[a] = b;
      }
    }

    if (ptr[a] == -1) { 
      throw std::runtime_error("[face_match] Failed to find matching nodes between faces '" + lFa.name + "' and '" + gFa.name + "'.");
    }
  }
}

//---------
// read_bc
//---------
// Read boundary condition data.
//
void read_bc(Simulation* simulation, EquationParameters* eq_params, eqType& lEq, BoundaryConditionParameters* bc_params, bcType& lBc)
{
  using namespace consts;
  auto bc_type = bc_params->type.value();

  if (std::set<std::string>{"Dirichlet", "Dir"}.count(bc_type)) {
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_Dir)); 

  } else if (std::set<std::string>{"Neumann", "Neu"}.count(bc_type)) {
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_Neu)); 
    if ((lEq.phys == EquationType::phys_fluid) || (lEq.phys == EquationType::phys_FSI)) {
      lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_bfs));
    }

  } else if (std::set<std::string>{"Traction", "Trac"}.count(bc_type)) {
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_trac)); 

    if (bc_params->traction_values_file_path.defined()) { 
      int iM = lBc.iM;
      int iFa = lBc.iFa;
      lBc.gm.nTP = 2;
      auto& com_mod = simulation->com_mod;
      int nsd = com_mod.nsd;
      lBc.gm.dof = nsd;

      lBc.gm.t.resize(2); 
      lBc.gm.d.resize(nsd, com_mod.msh[iM].fa[iFa].nNo, 2);

      lBc.gm.t[0] = 0.0;
      lBc.gm.t[1] = 1.E+10;
      lBc.gm.period = lBc.gm.t[1];
      auto file_name = bc_params->traction_values_file_path.value();

      read_trac_bcff(com_mod, lBc.gm, com_mod.msh[iM].fa[iFa], file_name);

      lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_gen));
      lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_flat));

      if (bc_params->traction_multiplier.defined()) { 
        double rtmp = bc_params->traction_multiplier.value(); 
        lBc.gm.d *= rtmp;
      }

      lBc.eDrn.resize(nsd); 
      lBc.h.resize(nsd);
      lBc.weakDir = false;
      return; 
    }

  } else if (std::set<std::string>{"Robin", "Rbn"}.count(bc_type)) {
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_Robin)); 
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_Neu)); 

  } else if (std::set<std::string>{"Coupled Momentum","CMM"}.count(bc_type)) {
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_CMM)); 

  } else if (std::set<std::string>{"RIS0D", "RIS0D"}.count(bc_type)) {
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_Ris0D));
    auto& com_mod = simulation->com_mod;
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_Neu));
    com_mod.ris0DFlag = true;
    lBc.resistance = bc_params->resistance.value();

  } else {
    throw std::runtime_error("[read_bc] Unknown boundary condition type '" + bc_type + "'.");
  }

  // Allocate traction
  auto& com_mod = simulation->com_mod;
  lBc.h.resize(com_mod.nsd);

  // [NOTE] Direction vectors can only have three values.
  //
  lBc.eDrn.resize(com_mod.nsd);
  auto effective_direction = bc_params->effective_direction();
  for (int i = 0; i <  effective_direction.size(); i++) {
    lBc.eDrn[i] = effective_direction[i];
  }

  auto ctmp = bc_params->time_dependence.value();

  if (ctmp == "Steady") { 
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_std)); 
    if (utils::btest(lBc.bType, enum_int(BoundaryConditionType::bType_trac))) { 
       lBc.h = bc_params->value.value();
    } else { 
      lBc.g = bc_params->value.value();
    }

  } else if (ctmp == "Unsteady") { 
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_ustd)); 

    if (utils::btest(lBc.bType, enum_int(BoundaryConditionType::bType_trac))) { 
      lBc.gt.d = com_mod.nsd;
    } else { 
      lBc.gt.d = 1;
    }

    if (bc_params->temporal_values_file_path.defined()) { 
      lBc.gt.lrmp = bc_params->ramp_function.value();
      auto file_name = bc_params->temporal_values_file_path.value();
      read_temporal_values(file_name, lBc);

    } else { 
      if (!bc_params->fourier_coefficients_file_path.defined()) { 
        throw std::runtime_error("[read_bc] Undefined data for boundary condition type '" + bc_type + "'.");
      }
      lBc.gt.lrmp = false;
      auto file_name = bc_params->fourier_coefficients_file_path.value(); 
      read_fourier_coeff_values_file(file_name, lBc);
    }

  // There are currently two coupling methods: GenBC and svZeroDSolver.
  //
  // A coupling method is defined if com_mod.cplBC.schm is set.
  //
  } else if (ctmp == "Coupled") { 
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_cpl)); 
    com_mod.cplBC.nFa = com_mod.cplBC.nFa + 1;
    lBc.cplBCptr = com_mod.cplBC.nFa - 1;
    auto& face_name = com_mod.msh[lBc.iM].fa[lBc.iFa].name;

    // The svZeroDSolver_interface parameter is defined.
    //
    if (com_mod.cplBC.svzerod_solver_interface.has_data) { 
      if (!bc_params->svzerod_solver_block.defined()) {
        std::string error_msg = std::string("The svZeroDSolver_block parameter must be defined for the 'Coupled' ") + 
            std::string(" boundary condition for the face '") + face_name + "'.";
        throw std::runtime_error(error_msg);
      }
      auto block_name = bc_params->svzerod_solver_block();
      com_mod.cplBC.svzerod_solver_interface.add_block_face(block_name, face_name);

    // Assume coupling with GenBC.
    //
    } else if (com_mod.cplBC.schm != CplBCType::cplBC_NA) {
      if (bc_params->svzerod_solver_block.defined()) {
        std::string error_msg = std::string("The svZeroDSolver_block parameter cannot be defined for the 'Coupled' ") + 
            std::string(" boundary condition for the face '") + face_name + "' when Couple_to_genBC parameters is used.";
        throw std::runtime_error(error_msg);
      }

    } else { 
      std::string error_msg = std::string("A coupling method must be defined for the 'Coupled' ") + 
            std::string(" boundary condition parameter for face '") + face_name + "'.";
      throw std::runtime_error(error_msg);
    }

  } else if (ctmp == "Resistance") { 
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_res)); 
    if (!utils::btest(lBc.bType, enum_int(BoundaryConditionType::bType_Neu))) { 
      throw std::runtime_error("[read_bc] Resistance is only defined for Neu BC.");
    }

    if (utils::btest(lBc.bType, enum_int(BoundaryConditionType::bType_Robin))) { 
      throw std::runtime_error("[read_bc] Resistance is not defined for Robin BC.");
    }

    if (std::set<EquationType>{Equation_fluid,Equation_FSI,Equation_CMM}.count(lEq.phys) == 0) {
      throw std::runtime_error("[read_bc] Resistance is only defined for fluid/CMM/SI equations.");
    }

    lBc.r = bc_params->value.value();

  ////////////////////////
  //       R C R      ////
  ////////////////////////

  } else if (std::set<std::string>{"RCR", "Windkessel"}.count(ctmp)) {
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_RCR)); 
    if (!utils::btest(lBc.bType, enum_int(BoundaryConditionType::bType_Neu))) { 
      throw std::runtime_error("[read_bc] RCR BC is only defined for Neu BC.");
    }

    if (std::set<EquationType>{Equation_fluid,Equation_FSI,Equation_CMM}.count(lEq.phys) == 0) {
      throw std::runtime_error("[read_bc] Resistance is only defined for fluid/CMM/FSI equations.");
    }

    lBc.RCR.Rp = bc_params->rcr.proximal_resistance.value();
    lBc.RCR.C  = bc_params->rcr.capacitance.value();
    lBc.RCR.Rd = bc_params->rcr.distal_resistance.value();
    lBc.RCR.Pd = bc_params->rcr.distal_pressure.value();
    lBc.RCR.Xo = bc_params->rcr.initial_pressure.value();

    if (com_mod.cplBC.schm != CplBCType::cplBC_NA || com_mod.cplBC.xo.size() != 0) {
      throw std::runtime_error("[read_bc] RCR cannot be used in conjunction with cplBC.");
    }
    com_mod.cplBC.nFa = com_mod.cplBC.nFa + 1;
    lBc.cplBCptr = com_mod.cplBC.nFa - 1;

  ////////////////////////
  //   S p a t i a l  ////
  ////////////////////////

  } else if (ctmp == "Spatial") { 
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_gen)); 
    if (!utils::btest(lBc.bType, enum_int(BoundaryConditionType::bType_Neu))) { 
      throw std::runtime_error("[read_bc] Spatial BC is only defined for Neu BC.");
    }

    int iM = lBc.iM;
    int iFa = lBc.iFa;
    lBc.gm.nTP = 2;
    lBc.gm.dof = 1;

    lBc.gm.d.resize(1, com_mod.msh[iM].fa[iFa].nNo, 2);
    lBc.gm.t.resize(2); 
    lBc.gm.t[0] = 0.0;
    lBc.gm.t[1] = 1.E+10;
    lBc.gm.period = lBc.gm.t[1];

    auto file_name = bc_params->spatial_values_file_path.value();
    read_trac_bcff(com_mod, lBc.gm, com_mod.msh[iM].fa[iFa], file_name);

  ////////////////////////
  //   G e n e r a l  ////
  ////////////////////////

  } else if (ctmp == "General") { 
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_gen)); 
    int iM = lBc.iM;
    int iFa = lBc.iFa;

    // [NOTE] This is not implemented.
    if (bc_params->bct_file_path.defined()) {
      auto file_name = bc_params->bct_file_path.value();
      read_bct(com_mod, lBc.gm, com_mod.msh[iM].fa[iFa], file_name);

    } else {
      if (bc_params->temporal_and_spatial_values_file_path.defined()) {
        auto file_name = bc_params->temporal_and_spatial_values_file_path.value(); 
        read_temp_spat_values(com_mod, com_mod.msh[iM], com_mod.msh[iM].fa[iFa], file_name, lBc); 
      } else {
        throw std::runtime_error("[read_bc] No bct.vtp input file provided for General BC.");
      }
    }

  } else {
    throw std::runtime_error("[read_bc] Unknown time dependence type '" + ctmp + "'.");
  }

  // Stiffness and damping parameters for Robin BC
  if (utils::btest(lBc.bType, enum_int(BoundaryConditionType::bType_Robin))) { 
    
    // Read VTP file path for per-node stiffness and damping (optional)
    if (bc_params->spatial_values_file_path.defined()) {
      lBc.robin_bc = RobinBoundaryCondition(bc_params->spatial_values_file_path.value(), 
                                            bc_params->apply_along_normal_direction.value(),
                                            com_mod.msh[lBc.iM].fa[lBc.iFa],
                                            simulation->logger);
    } else {
      lBc.robin_bc = RobinBoundaryCondition(bc_params->stiffness.value(), 
                                            bc_params->damping.value(),
                                            bc_params->apply_along_normal_direction.value(),
                                            com_mod.msh[lBc.iM].fa[lBc.iFa],
                                            simulation->logger);
    }
  }

  // To impose value or flux
  bool ltmp = bc_params->impose_flux.value();
  if (ltmp) {
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_flx)); 
  }

  // To zero-out perimeter or not. Default is .true. for Dir/CMM
  //
  ltmp = false; 
  ltmp = utils::btest(lBc.bType, enum_int(BoundaryConditionType::bType_Dir));
  if (bc_params->zero_out_perimeter.defined()) {
    ltmp = bc_params->zero_out_perimeter.value();
  }
  
  lBc.bType = utils::ibclr(lBc.bType, enum_int(BoundaryConditionType::bType_zp));
  if (ltmp || utils::btest(lBc.bType, enum_int(BoundaryConditionType::bType_CMM))) {
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_zp));
  }

  // Impose BC on the state variable or its integral
  //
  if (std::set<EquationType>{Equation_lElas,Equation_mesh,Equation_struct,Equation_shell}.count(lEq.phys) != 0) {
    ltmp = true;
  } else {
    ltmp = false;
  }

  if (bc_params->impose_on_state_variable_integral.defined()) {
    ltmp = bc_params->impose_on_state_variable_integral.value(); 
  }

  lBc.bType = utils::ibclr(lBc.bType, enum_int(BoundaryConditionType::bType_impD));
  if (ltmp) { 
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_impD));
  }

  // Reading the spatial profile: flat/para/ud.
  //
  ctmp = bc_params->profile.value();

  if (ctmp == "Flat") { 
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_flat));

  } else if (ctmp == "Parabolic") { 
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_para));

  } else if (ctmp == "User_defined") { 
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_ud));
    auto file_name =  bc_params->spatial_profile_file_path.value();

    int iM = lBc.iM;
    int iFa = lBc.iFa;

    read_spatial_values(com_mod, com_mod.msh[iM], com_mod.msh[iM].fa[iFa], file_name, lBc);
  } else {
    throw std::runtime_error("[read_bc] Unknown spatial profile type '" + ctmp + "'.");
  }

  // Weak Dirichlet BC for fluid/FSI equations
  //
  lBc.weakDir = false;
  if (utils::btest(lBc.bType, enum_int(BoundaryConditionType::bType_Dir))) {
    if ((lEq.phys == Equation_fluid) || (lEq.phys == Equation_FSI)) {
      lBc.weakDir = bc_params->weakly_applied.value();
    }
  }

  // Read penalty values for weakly applied Dir BC
  //
  if (lBc.weakDir) {
    lBc.bType = utils::ibclr(lBc.bType, enum_int(BoundaryConditionType::bType_zp));
    if (bc_params->penalty_parameter.defined()) {
      lBc.tauB = bc_params->penalty_parameter.value();
    }
    if (bc_params->penalty_parameter_tangential.defined()) {
      lBc.tauB[0] = bc_params->penalty_parameter_tangential.value();
    }
    if (bc_params->penalty_parameter_normal.defined()) {
      lBc.tauB[1] = bc_params->penalty_parameter_normal.value();
    }
  }

  // Read BCs for shells with triangular elements. Not necessary for
  // NURBS elements
  //
  if (bc_params->cst_shell_bc_type.defined()) { 
    auto ctmp = bc_params->cst_shell_bc_type.value(); 
    if (std::set<std::string>{"Fixed", "fixed", "Clamped", "clamped"}.count(ctmp)) {
      lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_fix));
      if (!utils::btest(lBc.bType, enum_int(BoundaryConditionType::bType_Dir))) {
        throw std::runtime_error("[read_bc] Fixed BC is only defined for Dirichlet boundaries.");
      }
    } else if (std::set<std::string>{"Hinged", "hinged"}.count(ctmp)) {
      lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_hing));
      if (!utils::btest(lBc.bType, enum_int(BoundaryConditionType::bType_Dir))) {
        throw std::runtime_error("[read_bc] Hinged BC is only defined for Dirichlet boundaries.");
      }

    } else if (std::set<std::string>{"Free", "free"}.count(ctmp)) {
      lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_free));
      if (!utils::btest(lBc.bType, enum_int(BoundaryConditionType::bType_Neu))) {
        throw std::runtime_error("[read_bc] Hinged BC is only defined for Neumann boundaries.");
      }
    }
  }

  //  For Neumann BC, is load vector changing with deformation (follower pressure)
  //
  lBc.flwP = false;
  if (utils::btest(lBc.bType, enum_int(BoundaryConditionType::bType_Neu))) {
    if (lEq.phys == Equation_struct || lEq.phys == Equation_ustruct) {
      lBc.flwP = bc_params->follower_pressure_load.value();
    }
  }

  // If a Neumann BC face is undeforming
  //
  lBc.masN = 0;

  if (utils::btest(lBc.bType, enum_int(BoundaryConditionType::bType_Neu))) {
    ltmp = false;
    lBc.bType = utils::ibclr(lBc.bType, enum_int(BoundaryConditionType::bType_undefNeu));
    ltmp = bc_params->undeforming_neu_face.value();

    if (ltmp) {
      if (lEq.phys != Equation_ustruct) {
        throw std::runtime_error("[read_bc] Undeforming Neu face is currently formulated for USTRUCT only.");
      }

      if (utils::btest(lBc.bType, enum_int(BoundaryConditionType::bType_cpl)) ||
          utils::btest(lBc.bType, enum_int(BoundaryConditionType::bType_res))) { 
        throw std::runtime_error("[read_bc] Undeforming Neu face cannot currently bbe used with a resistance or couple BC.");
      }

      // Clear any BC profile
      lBc.bType = utils::ibclr(lBc.bType, enum_int(BoundaryConditionType::bType_flat));
      lBc.bType = utils::ibclr(lBc.bType, enum_int(BoundaryConditionType::bType_para));

      if (utils::btest(lBc.bType, enum_int(BoundaryConditionType::bType_ud)) ||
          utils::btest(lBc.bType, enum_int(BoundaryConditionType::bType_gen))) { 
        throw std::runtime_error("[read_bc] General BC or user defined spatial profile cannot be imposed on an undeforming Neu face.");
      }

      // Clear zero perimeter flag
      lBc.bType = utils::ibclr(lBc.bType, enum_int(BoundaryConditionType::bType_zp));

      // Reset profile to flat and set undeforming Neumann BC flag
      lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_flat));
      lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_undefNeu));

      // Set master-slave node parameters. Set a master node that
      // is not part of any other face (not on perimeter)
      int iM  = lBc.iM;
      int iFa = lBc.iFa;

      Vector<int> ptr(com_mod.gtnNo);
      ptr = 0;

      for (int a = 0; a < com_mod.msh[iM].fa[iFa].nNo; a++) {
        int Ac = com_mod.msh[iM].fa[iFa].gN[a];
        ptr[Ac] = 1;
      }

      for (int j = 0; j < com_mod.msh[iM].nFa; j++) {
        if (j == iFa) {
          continue;
        } 
        for (int a = 0; a < com_mod.msh[iM].fa[j].nNo; a++) {
          int Ac = com_mod.msh[iM].fa[j].gN[a];
          ptr[Ac] = 0;
        }
      }

      for (int a = 0; a < com_mod.msh[iM].fa[iFa].nNo; a++) {
        int Ac = com_mod.msh[iM].fa[iFa].gN[a];
        if (ptr[Ac] == 1) {
          lBc.masN = Ac;
          break;
        }
      }
    }
  }

  // For CMM BC, load wall displacements
  //
  if (utils::btest(lBc.bType, enum_int(BoundaryConditionType::bType_CMM))) { 
    auto cTmp = bc_params->initial_displacements_file_path.value();
    if (!bc_params->initial_displacements_file_path.defined() && (com_mod.Dinit.size() == 0)) {  
      cTmp = bc_params->prestress_file_path.value();
      if (!bc_params->prestress_file_path.defined() && (com_mod.pS0.size() == 0)) {  
        throw std::runtime_error("[read_bc] No wall displacement field or prestress given for CMM."); 
      }

      // Read prestress tensor here
      //
      if (cTmp != "") {
        int iM = lBc.iM;
        int iFa = lBc.iFa;
        auto& face = com_mod.msh[iM].fa[iFa];

         // [NOTE] What's this all about? 
         //
        if (face.x.size() == 0) {
          //face.x.resize(com_mod.nsymd, face.nNo);
          //ALLOCATE(msh(iM).fa(iFa).x(nsymd,msh(iM).fa(iFa).nNo))
          //msh(iM).fa(iFa).x = 0._RKIND
        }

        // Data is stored in face.x so make sure it is the right size.
        face.x.resize(com_mod.nsymd, face.nNo);
        int data_series = 0;
        vtk_xml::read_vtp_pdata(cTmp, "Stress", com_mod.nsd, com_mod.nsymd, data_series, face);

        if (com_mod.pS0.size() == 0) {
          com_mod.pS0.resize(com_mod.nsymd, com_mod.gtnNo);
        }

        for (int a = 0; a < face.nNo; a++) {
          int Ac = face.gN[a];
          Ac = com_mod.msh[iM].gN[Ac];
          for (int i = 0; i < face.x.nrows(); i++) {
            com_mod.pS0(i,Ac) = face.x(i,a);
          }
        }
        face.x.clear();
      }
    }

    //  Read displacement field here
    //
    if (cTmp != "") {
      int iM = lBc.iM;
      int iFa = lBc.iFa;
      auto& face = com_mod.msh[iM].fa[iFa];
      face.x.resize(com_mod.nsd, face.nNo);

      int data_series = 0;
      vtk_xml::read_vtp_pdata(cTmp, "Displacement", com_mod.nsd, com_mod.nsd, data_series, face);

      if (com_mod.Dinit.size() == 0) {
        com_mod.Dinit.resize(com_mod.nsd, com_mod.gtnNo);
      }

      for (int a = 0; a < face.nNo; a++) { 
        int Ac = face.gN[a];
        Ac = com_mod.msh[iM].gN[Ac];
        for (int i = 0; i < face.x.nrows(); i++) {
          com_mod.Dinit(i,Ac) = face.x(i,a);
        }
      }
      face.x.clear();
    }

    // Set cmmBdry vector for wall nodes
    int iM  = lBc.iM;
    int iFa = lBc.iFa;
    auto& face = com_mod.msh[iM].fa[iFa];

    for (int a = 0; a < face.nNo; a++) {
      int Ac = face.gN[a];
      Ac = com_mod.msh[iM].gN[Ac];
      com_mod.cmmBdry[Ac] = 1;
    }
  }
}

//----------
// read_bct
//----------
// Reads general velocity data from bct.vtp.
//
// Reproduces 'SUBROUTINE READBCT(lMB, lFa, fName)' in READFILES.f.
//
void read_bct(ComMod& com_mod, MBType& lMB, faceType& lFa, const std::string& fName)
{
  // The name of the arrays to search for in the bct.vtp file.
  static std::string shdr("velocity_");

  static std::string file_desc("general velocity data file '" + fName + "'"); 

  if (FILE *file = fopen(fName.c_str(), "r")) {
      fclose(file);
  } else {
    throw std::runtime_error("The " + file_desc + " can't be read.");
  }

  // Read the vtp file.
  //
  VtkVtpData vtp_data(fName);
  int nNo = vtp_data.num_points();
  if (nNo == 0) {
    throw std::runtime_error("The " + file_desc + " does not contain any points.");
  }

  int num_elems = vtp_data.num_elems();
  if (num_elems == 0) {
    throw std::runtime_error("The " + file_desc + " does not contain any elements.");
  }

  if (nNo != lFa.nNo) {
    throw std::runtime_error("The number of points (" + std::to_string(nNo) + ") in the " + file_desc + 
        " does not match the number of points (" + std::to_string(nNo) + " for the face '" + lFa.name + "'.");
  }

  // Get all the point data starting with "velocity_"
  //
  auto namesL = vtp_data.get_point_data_names();
  int n = namesL.size();

  int ntime = 0;
  int nj = 1;

  for (int i = 0; i < n; i++) {
    auto stmp = namesL[i];

    for (int j = 0; j < 9; j++) {
      if (shdr[j] != stmp[j]) {
        break; 
      }
      nj += 1;
    }

    if (nj < 9) {
      namesL[i] = "";
      continue; 
    }
    ntime = ntime + 1;
  } 

  // Initialize lMB data structure
  //
  const int nsd = com_mod.nsd;
  lMB.dof = nsd;
  lMB.nTP = ntime;
  int iM = lFa.iM;
  lMB.t.resize(ntime); 
  lMB.d.resize(nsd,nNo,ntime); 
  Vector<int> ptr(com_mod.msh[iM].gnNo);

  ntime = 0;

  for (int i = 0; i < n; i++) {
    auto stmp = namesL[i];
    if (stmp.size() == 0) { 
      continue;
    }

    stmp = stmp.substr(9);
    double t = std::stod(stmp);
    lMB.t(ntime) = t;

    if (ntime == 0) { 
      if (!utils::is_zero(t)) {
        throw std::runtime_error("The first time step in the " + file_desc + " should be zero.");
      }
    } else {
      t = t - lMB.t(ntime-1);
      if (utils::is_zero(t) || t < 0.0) {
        throw std::runtime_error("There is a non-increaing series of times in the " + file_desc + ".");
      }
    }

    ntime = ntime + 1;
  }

  lMB.period = lMB.t(ntime-1);
  auto& msh = com_mod.msh[iM];

  // Prepare pointer array
  //
  for (int a = 0; a < lFa.nNo; a++) {
    int Ac = lFa.gN(a);
    Ac = msh.lN(Ac);

    if (Ac == -1) {
      throw std::runtime_error("Incorrect global node number detected in the " + file_desc + 
          ".  Mesh: '"  + msh.name + "'  Face: '" + lFa.name + "'  Node: " + std::to_string(a+1) + 
          " gN: " + std::to_string(lFa.gN(a)));
    }

    ptr(Ac) = a;
  } 

  // Get GlobalNodeID from vtp file and make sure it is consistent
  // with mesh structure
  //
  Vector<int> gN(nNo);
  vtp_data.copy_point_data("GlobalNodeID", gN);

  for (int a = 0; a < nNo; a++) {
    int Ac = gN(a) - 1;
    if (Ac > msh.gnNo-1 || Ac < 0) {
      throw std::runtime_error("Global node IDs in the '" + file_desc + " are not consistent with the volume mesh '" + msh.name + ".");
    }
    Ac = ptr(Ac);
    if (Ac == -1) { 
      throw std::runtime_error("Global node IDs in the '" + file_desc + " are not consistent with the volume mesh '" + msh.name + ".");
    }
  }

  // Load spatial data for each time point from vtp file
  //
  Array<double> tmpR(nsd,nNo);
  ntime = 0;

  for (int i = 0; i < n; i++) {
    auto stmp = namesL[i];
    if (stmp.size() == 0) { 
      continue;
    }

    tmpR = 0.0;
    vtp_data.copy_point_data(stmp, tmpR);

    for (int a = 0; a < nNo; a++) {
      int Ac = gN(a) - 1;
      Ac = ptr(Ac);
      for (int j = 0; j < nsd; j++) {
        lMB.d(j,Ac,ntime) = tmpR(j,a);
      }
    }

    ntime = ntime + 1;
  }

}

//---------
// read_bf
//---------
// Read body force data.
//
void read_bf(ComMod& com_mod, BodyForceParameters* bf_params, bfType& lBf) 
{
  using namespace consts;

  int iM = lBf.iM;
  lBf.dof = com_mod.nsd;

  // Reading the type: Vol/Neu/Trac
  //
  auto bf_type = bf_params->type.value();
  std::transform(bf_type.begin(), bf_type.end(), bf_type.begin(), ::tolower);

  if (std::set<std::string>{"volumetric","vol","internal","int"}.count(bf_type)) {
    lBf.bType = utils::ibset(lBf.bType, enum_int(BodyForceType::bfType_vol)); 

  } else if (std::set<std::string>{"traction","trac"}.count(bf_type)) {
    lBf.bType = utils::ibset(lBf.bType, enum_int(BodyForceType::bfType_trac)); 

  } else if (std::set<std::string>{"neumann","neu", "pressure"}.count(bf_type)) {
    lBf.bType = utils::ibset(lBf.bType, enum_int(BodyForceType::bfType_Neu)); 
    lBf.dof = 1;

  } else {
    throw std::runtime_error("Unknown body force type '" + bf_params->type.value() + "'.");
  }

  // Time dependence
  //
  auto time_dependence = bf_params->time_dependence.value();
  std::transform(time_dependence.begin(), time_dependence.end(), time_dependence.begin(), ::tolower);

  // Steady //

  if (time_dependence == "steady") {
    lBf.bType = utils::ibset(lBf.bType, enum_int(BodyForceType::bfType_std)); 
    lBf.b.resize(lBf.dof);
    if (lBf.dof == 1) {
      lBf.b[0] = bf_params->value.value();
    } else {
      lBf.b = bf_params->value.value();
    }

  // Unsteady //

  } else if (time_dependence == "unsteady") {
    lBf.bType = utils::ibset(lBf.bType, enum_int(BodyForceType::bfType_ustd)); 
    lBf.bt.d = lBf.dof;
    if (bf_params->temporal_values_file_path.defined()) {
      lBf.bt.lrmp = bf_params->ramp_function.value();
      auto file_name = bf_params->temporal_values_file_path.value();
      read_temporal_values(file_name, lBf);
    } else {
      lBf.bt.lrmp = false;
      auto fTmp = bf_params->fourier_coefficients_file_path.value();
      read_fourier_coeff_values_file(fTmp, lBf);
    }

  // Spatial //

  } else if (time_dependence == "spatial") {
    lBf.bType = utils::ibset(lBf.bType, enum_int(BodyForceType::bfType_spl)); 
    lBf.bx.resize(lBf.dof, com_mod.gtnNo);
    auto cTmp = bf_params->spatial_values_file_path.value();
    com_mod.msh[iM].x.resize(lBf.dof, com_mod.msh[iM].gnNo);
    int data_comp = lBf.dof;
    int data_series = 0;

    bool is_vtu_file;
    if (cTmp.substr(cTmp.find_last_of(".") + 1) == "vtu") {
      is_vtu_file = true;
    } else {
      is_vtu_file = false;
    }

    if (utils::btest(lBf.bType, enum_int(BodyForceType::bfType_vol))) { 
      vtk_xml::read_vtu_pdata(cTmp, "Body_force", com_mod.nsd, data_comp, data_series, com_mod.msh[iM]);

    } else if (utils::btest(lBf.bType, enum_int(BodyForceType::bfType_trac))) { 
      vtk_xml::read_vtu_pdata(cTmp, "Traction", com_mod.nsd, data_comp, data_series, com_mod.msh[iM]);

    } else if (utils::btest(lBf.bType, enum_int(BodyForceType::bfType_Neu))) { 
      vtk_xml::read_vtu_pdata(cTmp, "Pressure", com_mod.nsd, data_comp, data_series, com_mod.msh[iM]);
    }

    for (int a = 0; a <  com_mod.msh[iM].gnNo; a++) {
      int Ac = com_mod.msh[iM].gN[a];
      for (int i = 0; i <  lBf.bx.nrows(); i++) {
        lBf.bx(i,Ac) = com_mod.msh[iM].x(i,a);
      }
    }
    com_mod.msh[iM].x.clear();

  // General //

  } else if (time_dependence == "general") {
    lBf.bType = utils::ibset(lBf.bType, enum_int(BodyForceType::bfType_gen)); 
    auto file_name = bf_params->temporal_and_spatial_values_file_path.value();
    read_temp_spat_values(com_mod, com_mod.msh[iM], file_name, lBf);

  } else { 
    throw std::runtime_error("Unknown body force time dependence type '" + time_dependence + "'.");
  }
}

//----------
// read_cep_domain
//----------
// Set domain-specific parameters for a cardiac electrophysiology model.
//
void read_cep_domain(Simulation* simulation, EquationParameters* eq_params, DomainParameters* domain_params, dmnType& lDmn)
{ 
  auto model_str = domain_params->electrophysiology_model.value();
  std::transform(model_str.begin(), model_str.end(), model_str.begin(), ::tolower);

  // Get the type of electrophysiology model.
  //
  ElectrophysiologyModelType model_type;

  try {
    model_type = cep_model_name_to_type.at(model_str);
  } catch (const std::out_of_range& exception) {
    throw std::runtime_error("[read_cep_domain] Unknown model type '" + model_str + "'.");
  }
  
  // Set model parameters based on model type.
  //
  lDmn.cep.nG = 0;
  lDmn.cep.cepType = model_type;

  switch (model_type) {
    case ElectrophysiologyModelType::AP: 
    case ElectrophysiologyModelType::FN:
      lDmn.cep.nX = 2;
    break;

    case ElectrophysiologyModelType::BO:
      lDmn.cep.nX = 4;
    break;

    case ElectrophysiologyModelType::TTP:
      lDmn.cep.nX = 7;
      lDmn.cep.nG = 12;
    break;

    default: 
    break;
  }
  
  // Set the maximum number of dof for cellular activation model.
  auto& cep_mod = simulation->get_cep_mod();
  if (cep_mod.nXion < lDmn.cep.nX + lDmn.cep.nG) {
    cep_mod.nXion = lDmn.cep.nX + lDmn.cep.nG;
  } 

  // Set conductivity.
  lDmn.cep.Diso = domain_params->isotropic_conductivity();

  // Set the number of fiber directions (nFn) determined by the
  // given number of anisotropic conductivity values.
  lDmn.cep.nFn = domain_params->anisotropic_conductivity.size(); 
  if (lDmn.cep.nFn == 0) { 
    lDmn.cep.nFn = 1;
  }

  // Set the anisotropic conductivity values. 
  lDmn.cep.Dani.resize(lDmn.cep.nFn);
  for (int i = 0; i < domain_params->anisotropic_conductivity.size(); i++) {
    lDmn.cep.Dani[i] = domain_params->anisotropic_conductivity[i];
  }

  // Set myocardium zone id.
  //
  if (domain_params->myocardial_zone.defined()) { 
    auto myocardial_zone = domain_params->myocardial_zone.value();
    std::transform(myocardial_zone.begin(), myocardial_zone.end(), myocardial_zone.begin(), ::tolower);
    if (std::set<std::string>{"epi", "epicardium"}.count(myocardial_zone)) {
      lDmn.cep.imyo = 1;
    } else if (std::set<std::string>{"endo", "endocardium", "pfib", "purkinje"}.count(myocardial_zone)) {
      lDmn.cep.imyo = 2;
    } else if (std::set<std::string>{"myo", "mid-myo", "myocardium"}.count(myocardial_zone)) {
      lDmn.cep.imyo = 3;
    } else {
      throw std::runtime_error("Unknown myocardial zone type '" + myocardial_zone + "'.");
    }
  }

  // Set Ttp parameters.
  //
  if (domain_params->G_Na.defined())  { cep_mod.ttp.G_Na = domain_params->G_Na.value(); }
  if (domain_params->G_Kr.defined())  { cep_mod.ttp.G_Kr = domain_params->G_Kr.value(); }
  if (domain_params->G_Ks.defined())  { cep_mod.ttp.G_Ks[lDmn.cep.imyo - 1] = domain_params->G_Ks.value(); }
  if (domain_params->G_to.defined())  { cep_mod.ttp.G_to[lDmn.cep.imyo - 1] = domain_params->G_to.value(); }
  if (domain_params->G_CaL.defined()) { cep_mod.ttp.G_CaL = domain_params->G_CaL.value(); }

  // Set Bo parameters.
  //
  if (domain_params->tau_si.defined())  { cep_mod.bo.tau_si[lDmn.cep.imyo - 1] = domain_params->tau_si.value(); }
  if (domain_params->tau_fi.defined())  { cep_mod.bo.tau_fi[lDmn.cep.imyo - 1] = domain_params->tau_fi.value(); }

  // Set stimulus parameters. 
  //
  lDmn.cep.Istim.A  = 0.0;
  lDmn.cep.Istim.Ts = 99999.0;
  lDmn.cep.Istim.CL = 99999.0;
  lDmn.cep.Istim.Td = 0.0;

  if (domain_params->stimulus.defined()) { 
    auto& stimulus_params = domain_params->stimulus;
    lDmn.cep.Istim.A = stimulus_params.amplitude.value();
    if (!utils::is_zero(lDmn.cep.Istim.A)) {
      lDmn.cep.Istim.Ts = stimulus_params.start_time.value();
      lDmn.cep.Istim.Td = stimulus_params.duration.value();
      if (stimulus_params.cycle_length.defined()) { 
        lDmn.cep.Istim.CL = stimulus_params.cycle_length.value();
      } else {
        lDmn.cep.Istim.CL = simulation->nTs * simulation->com_mod.dt; 
      }
    }
  }

  // Dual time step for cellular activation model.
  if (domain_params->time_step_for_integration.defined()) {
    lDmn.cep.dt = domain_params->time_step_for_integration.value();
  } else { 
    lDmn.cep.dt = simulation->com_mod.dt;
  }

  // Set time integrator.
  //
  TimeIntegratioType time_integration_type;
  auto ode_solver_str = domain_params->ode_solver.value();
  std::transform(ode_solver_str.begin(), ode_solver_str.end(), ode_solver_str.begin(), ::tolower);
  try {
    time_integration_type = cep_time_int_to_type.at(ode_solver_str);
  } catch (const std::out_of_range& exception) {
    throw std::runtime_error("[read_cep_domain] Unknown ode solver type '" + ode_solver_str + "'.");
  }
  lDmn.cep.odes.tIntType = time_integration_type;

  if ((lDmn.cep.odes.tIntType == TimeIntegratioType::CN2) && (lDmn.cep.cepType == ElectrophysiologyModelType::TTP)) {
    throw std::runtime_error("[read_cep_domain] Implicit time integration for tenTusscher-Panfilov model can give unexpected results. Use FE or RK4 instead");
  }

  if (lDmn.cep.odes.tIntType == TimeIntegratioType::CN2) {
    lDmn.cep.odes.maxItr = 5;
    lDmn.cep.odes.absTol = 1e-8;
    lDmn.cep.odes.relTol = 1e-4;
    lDmn.cep.odes.maxItr = domain_params->maximum_iterations.value();
    lDmn.cep.odes.absTol = domain_params->absolute_tolerance.value();
    lDmn.cep.odes.relTol = domain_params->relative_tolerance.value();
  }

  if (domain_params->feedback_parameter_for_stretch_activated_currents.defined() && cep_mod.cem.cpld) { 
    lDmn.cep.Ksac = domain_params->feedback_parameter_for_stretch_activated_currents.value();
  } else {
    lDmn.cep.Ksac = 0.0;
  }
}

//----------
// read_cep_equation
//----------
// Set parameters that are affecting the whole equation of a cardiac electrophysiology model (for the moment ECG leads only)
//
void read_cep_equation(CepMod* cep_mod, Simulation* simulation, EquationParameters* eq_params)
{
  // Set ECG leads parameters.
  std::vector<double> x_coords, y_coords, z_coords;
  auto& chnl_mod = simulation->get_chnl_mod();
  if (eq_params->ecg_leads.defined()) {

    if (simulation->get_com_mod().nsd != 3) {
      throw std::runtime_error("ECG leads computation is allowed only for 3D geometries");
    }

    std::string line;
    auto& ecg_leads_params = eq_params->ecg_leads;

    std::ifstream x_coords_file;
    std::string x_coords_file_name = ecg_leads_params.x_coords_file_path.value();
    x_coords_file.open(x_coords_file_name);
    if (!x_coords_file.is_open()) {
      throw std::runtime_error("[read_cep_equation] Failed to open the ECG leads x-coordinates file '" + x_coords_file_name + "'.");
    }
    while(std::getline(x_coords_file, line)) { x_coords.push_back(std::stod(line)); }
    x_coords_file.close();

    std::ifstream y_coords_file;
    std::string y_coords_file_name = ecg_leads_params.y_coords_file_path.value();
    y_coords_file.open(y_coords_file_name);
    if (!y_coords_file.is_open()) {
      throw std::runtime_error("[read_cep_equation] Failed to open the ECG leads y-coordinates file '" + y_coords_file_name + "'.");
    }
    while(std::getline(y_coords_file, line)) { y_coords.push_back(std::stod(line)); }
    y_coords_file.close();

    std::ifstream z_coords_file;
    std::string z_coords_file_name = ecg_leads_params.z_coords_file_path.value();
    z_coords_file.open(z_coords_file_name);
    if (!z_coords_file.is_open()) {
      throw std::runtime_error("[read_cep_equation] Failed to open the ECG leads z-coordinates file '" + z_coords_file_name + "'.");
    }
    while(std::getline(z_coords_file, line)) { z_coords.push_back(std::stod(line)); }
    z_coords_file.close();

    if (x_coords.size() != y_coords.size() ||
        y_coords.size() != z_coords.size()) {
      throw std::runtime_error("[read_cep_equation] ECG leads for x,y,z-coordinates must have the same dimension.");
    }
    cep_mod->ecgleads.num_leads = x_coords.size();

    for (int index = 0; index < cep_mod->ecgleads.num_leads; index++) {
      cep_mod->ecgleads.out_files.push_back(simulation->chnl_mod.appPath + "ecglead_" + std::to_string(index + 1) + ".txt");
    }

    cep_mod->ecgleads.x_coords.set_values(x_coords);
    cep_mod->ecgleads.y_coords.set_values(y_coords);
    cep_mod->ecgleads.z_coords.set_values(z_coords);

    cep_mod->ecgleads.pseudo_ECG.resize(cep_mod->ecgleads.num_leads);
    std::fill(cep_mod->ecgleads.pseudo_ECG.begin(), cep_mod->ecgleads.pseudo_ECG.end(), 0.);
  }
}

//--------------------------------
// read_cplbc_initialization_file
//--------------------------------
// Read the float values for a cplBC initialzation.
//
void read_cplbc_initialization_file(const std::string& file_name, cplBCType& cplBC) 
{
  std::ifstream init_file;
  init_file.open(file_name);
  if (!init_file.is_open()) {
    throw std::runtime_error("Failed to open the cplBC initialization file '" + file_name + "'.");
  }

  double value;
  std::string str_value;
  int n = 0;

  while (init_file >> value) {
    if (n == cplBC.nX) { 
      throw std::runtime_error("The number of values in the cplBC initialization file '" + file_name + 
         "' is larger than the number of values given in the Couple_to_cplBC 'Number_of_unknowns' parameter (" + 
         std::to_string(cplBC.nX) + ").");
    }
    cplBC.xo[n] = value;
    n += 1;
  }

  if (n < cplBC.nX) { 
    throw std::runtime_error("The number of values in the cplBC initialization file '" + file_name + 
      "' (" + std::to_string(n) + ") is smaller than the number of values given in the Couple_to_cplBC 'Number_of_unknowns' parameter (" + 
      std::to_string(cplBC.nX) + ").");
  }
}

//-------------
// read_domain
//-------------
// Set all domain properties.
//
void read_domain(Simulation* simulation, EquationParameters* eq_params, eqType& lEq, EquationProps& propList, EquationPhys physI)
{
  using namespace consts;
  auto& com_mod = simulation->get_com_mod();
  int nPhys;
  std::vector<EquationType> phys;

  if (physI.size() != 0) {
    nPhys = physI.size();
    phys = physI;
  } else {
    nPhys = 1;
    phys.push_back(lEq.phys);
  }

  bool flag = false;
  bool domain_defined = true;
  lEq.nDmn = eq_params->domains.size();

  if (lEq.nDmn == 0) {
    lEq.nDmn = 1;
    flag = true;
    domain_defined = false;
  }

  lEq.dmn.resize(lEq.nDmn);

  for (int iDmn = 0; iDmn < lEq.nDmn; iDmn++) {
     auto& domain_params = (eq_params->domains.size() == 0 ? eq_params->default_domain : eq_params->domains[iDmn]);
     if (flag) {
        // If no domain keywork found, we search upper level for data
        lEq.dmn[iDmn].Id = -1;
     } else { 
       // [NOTE] It is not clear exactly how domain IDs are set.

       for (int i = 0; i < iDmn-1; i++) {
         if (lEq.dmn[iDmn].Id == lEq.dmn[i].Id) {
           throw std::runtime_error("Duplicate domain ID " + std::to_string(lEq.dmn[iDmn].Id) +  ".");
         }
       }
       lEq.dmn[iDmn].Id = std::stoi(domain_params->id()); 
     }

     // For a fluid-solid multi-physics system set the fluid and
     // solid equaion types for each domain.
     //
     if (nPhys > 1) {
        if (lEq.phys != EquationType::phys_FSI) {
          throw std::runtime_error("Only FSI problems can have multiple domains.");
        }
        auto eq_type = consts::equation_name_to_type.at(domain_params->equation());
        lEq.dmn[iDmn].phys = eq_type; 

        if (eq_type == EquationType::phys_ustruct) { 
           if (!com_mod.sstEq) com_mod.sstEq = true;
        }

     // Single physics systsm.
     } else { 
        lEq.dmn[iDmn].phys = lEq.phys;
     }

     // Find the index into the current domain's physics. 
     int iPhys;
     for (int i = 0; i < nPhys; i++) {
        if (lEq.dmn[iDmn].phys == phys[i]) { 
          iPhys = i;
          break;
        }
     }

     if (iPhys > nPhys) {
       throw std::runtime_error("Undefined phys is used.");
     }

     // Set domain properties.
     //
     // Get property values from either equation or domain parameters.
     //
     for (int iProp = 0; iProp < maxNProp; iProp++) {
        double rtmp = 0.0;
        auto prop = propList[iProp][iPhys];

        switch (prop) {
          case PhysicalProperyType::backflow_stab:
            rtmp = domain_params->backflow_stabilization_coefficient.value();
          break;

          case PhysicalProperyType::conductivity:
            rtmp = domain_params->conductivity.value();
          break;

          case PhysicalProperyType::ctau_C:
            rtmp = domain_params->continuity_stabilization_coefficient.value(); 
          break;

          case PhysicalProperyType::ctau_M:
            rtmp = domain_params->momentum_stabilization_coefficient.value();
          break;

          case PhysicalProperyType::damping:
            rtmp = domain_params->mass_damping.value();
          break;

          case PhysicalProperyType::elasticity_modulus:
            rtmp = domain_params->elasticity_modulus.value();
          break;

          case PhysicalProperyType::f_x:
            rtmp = domain_params->force_x.value();
          break;

          case PhysicalProperyType::f_y:
            rtmp = domain_params->force_y.value();
          break;

          case PhysicalProperyType::f_z:
            rtmp = domain_params->force_z.value();
          break;

          case PhysicalProperyType::fluid_density:
            if (lEq.phys == EquationType::phys_CMM) {
              rtmp = domain_params->fluid_density.value();
            } else {
              rtmp = domain_params->density.value();
            }
          break;

          case PhysicalProperyType::poisson_ratio:
            rtmp = domain_params->poisson_ratio.value();
          break;

          case PhysicalProperyType::shell_thickness:
            rtmp = domain_params->shell_thickness.value();
          break;

          case PhysicalProperyType::solid_density:
            if (lEq.phys == EquationType::phys_CMM) {
              rtmp = domain_params->solid_density.value();
            } else {
              rtmp = domain_params->density.value();
            }
          break;

          case PhysicalProperyType::source_term:
            rtmp = domain_params->source_term.value();
          break;

          case PhysicalProperyType::inverse_darcy_permeability:
            rtmp = domain_params->inverse_darcy_permeability.value();
          break;
        }

        lEq.dmn[iDmn].prop[prop] = rtmp;
     }

     // Set parameters for a cardiac electrophysiology model.
     if (lEq.dmn[iDmn].phys == EquationType::phys_CEP) {
        read_cep_domain(simulation, eq_params, domain_params, lEq.dmn[iDmn]);
     }

     // Read material/constitutive model parameters for nonlinear
     // elastodynamics simulations (both solids and shells)
     //
     if ( (lEq.dmn[iDmn].phys == EquationType::phys_shell) || 
          (lEq.dmn[iDmn].phys == EquationType::phys_struct) || 
          (lEq.dmn[iDmn].phys == EquationType::phys_ustruct)) { 
        read_mat_model(simulation, eq_params, domain_params, lEq.dmn[iDmn]);
        if (utils::is_zero(lEq.dmn[iDmn].stM.Kpen) && lEq.dmn[iDmn].phys == EquationType::phys_struct) { 
          //err = "Incompressible struct is not allowed. Use "//  "penalty method or ustruct"
          throw std::runtime_error("An incompressible material model is not allowed for 'struct' physics; use penalty method or ustruct."); 
        }
     }

     // Set parameters for a fluid viscosity model.
     if ((lEq.dmn[iDmn].phys == EquationType::phys_fluid) ||  
         (lEq.dmn[iDmn].phys == EquationType::phys_stokes) ||  
         (lEq.dmn[iDmn].phys == EquationType::phys_CMM && !com_mod.cmmInit)) {
       read_fluid_visc_model(simulation, eq_params, domain_params, lEq.dmn[iDmn]);
     }

     // Set parameters for a fluid viscosity model.
     if ((lEq.dmn[iDmn].phys == EquationType::phys_struct) ||  
         (lEq.dmn[iDmn].phys == EquationType::phys_ustruct)) {
       read_solid_visc_model(simulation, eq_params, domain_params, lEq.dmn[iDmn]);
     }
  }
}

//---------
// read_eq
//---------
// Set equation parameters.
//
void read_eq(Simulation* simulation, EquationParameters* eq_params, eqType& lEq)
{
  using namespace consts;

  lEq.useTLS  = false;
  lEq.assmTLS = false;

  lEq.coupled = eq_params->coupled.value();
  lEq.minItr = eq_params->min_iterations.value();
  lEq.maxItr = eq_params->max_iterations.value();
  lEq.tol = eq_params->tolerance.value();

  // Initialize coupled BC.
  //
  auto& chnl_mod = simulation->chnl_mod;
  auto& com_mod = simulation->get_com_mod();
  auto& cplBC = com_mod.cplBC;

  if (cplBC.xo.size() == 0) {
    std::string cplbc_type_str;

    if (eq_params->couple_to_genBC.defined()) {
      cplBC.useGenBC = true;
      cplbc_type_str = eq_params->couple_to_genBC.type.value();

    } else if (eq_params->svzerodsolver_interface_parameters.defined()) {
      cplBC.useSvZeroD = true;
      cplbc_type_str = eq_params->svzerodsolver_interface_parameters.coupling_type.value();
      cplBC.svzerod_solver_interface.set_data(eq_params->svzerodsolver_interface_parameters);

    } else if (eq_params->couple_to_cplBC.defined()) {
      cplbc_type_str = eq_params->couple_to_cplBC.type.value();
    }

    if (eq_params->couple_to_genBC.defined() || eq_params->couple_to_cplBC.defined() || eq_params->svzerodsolver_interface_parameters.defined()) { 
      try {
        cplBC.schm = consts::cplbc_name_to_type.at(cplbc_type_str);
      } catch (const std::out_of_range& exception) {
        throw std::runtime_error("Unknown coupled BC type '" + cplbc_type_str + ".");
      }
    }

    if (cplBC.schm != consts::CplBCType::cplBC_NA) { 
      if (cplBC.useGenBC) {
        if (cplBC.binPath.size() == 0){
          cplBC.binPath = eq_params->couple_to_genBC.zerod_code_file_path.value();
        }
        cplBC.commuName = "GenBC.int";
        cplBC.nX = 0;
        cplBC.xp.resize(cplBC.nX);

      } else if (cplBC.useSvZeroD) {
        cplBC.nX = 0;

      } else {
        auto& cplBC_params = eq_params->couple_to_cplBC;
        cplBC.nX = cplBC_params.number_of_unknowns.value();
        cplBC.xo.resize(cplBC.nX);
        cplBC.binPath = cplBC_params.zerod_code_file_path.value();
        if (cplBC_params.unknowns_initialization_file_path.defined()) { 
          auto file_name = cplBC_params.unknowns_initialization_file_path.value();
          read_cplbc_initialization_file(file_name, cplBC); 
        }

        cplBC.commuName = simulation->chnl_mod.appPath + cplBC_params.file_name_for_0D_3D_communication.value();
        cplBC.saveName = simulation->chnl_mod.appPath + cplBC_params.file_name_for_saving_unknowns.value();

        cplBC.nXp = cplBC_params.number_of_user_defined_outputs.value();
        cplBC.xp.resize(cplBC.nXp);
      }
    }
  }

  // Set equation properties.
  //
  if (eq_params->prestress.defined() && eq_params->prestress.value()) {
    simulation->com_mod.pstEq = true;
  }

  bool THflag = false; 
  if (eq_params->use_taylor_hood_type_basis.defined()) { 
    THflag = eq_params->use_taylor_hood_type_basis.value(); 
  }
  EquationProps propL{consts::PhysicalProperyType::NA};
  EquationOutputs outPuts;
  EquationNdop nDOP;

  set_equation_properties(simulation, eq_params, lEq, propL, outPuts, nDOP);

  // Read VTK files or boundaries. [TODO:DaveP] this is not a correct comment.
  read_outputs(simulation, eq_params, lEq, nDOP, outPuts);

  // Set the number of function spaces
  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    auto& msh = com_mod.msh[iM];
    if (!THflag) {
      msh.nFs = 1;
    } else { 
      if (com_mod.ibFlag) {
        throw std::runtime_error("Taylor-Hood basis is not implemented for immersed boundaries.");
      }
      if (msh.lShl || msh.lFib || (msh.eType == ElementType::NRB)) {
        msh.nFs = 1;
        throw std::runtime_error("Taylor-Hood basis is not allowed for NURBS mesh or shells and fibers.");
      } else { 
        if ( (msh.eType != ElementType::TRI6)  && (msh.eType != ElementType::TET10)) {
          throw std::runtime_error("Taylor-Hood basis is currently applicable for TRI6 (2D) or TET10 (3D) elements only.");
        }
        msh.nFs = 2;
      }
    }
  }

  // Searching for BCs
  //
  int num_bcs = eq_params->boundary_conditions.size();
  lEq.nBc = num_bcs;

  if (num_bcs != 0) { 
    lEq.bc.resize(num_bcs);
  }

  for (int iBc = 0; iBc < num_bcs; iBc++) {
    auto& bc_params = eq_params->boundary_conditions[iBc];
    auto bc_name = bc_params->name.value();
    int iM, iFa;

    all_fun::find_face(com_mod.msh, bc_name, iM, iFa);

    lEq.bc[iBc].iM = iM; 
    lEq.bc[iBc].iFa = iFa;

    read_bc(simulation, eq_params, lEq, bc_params, lEq.bc[iBc]);
  }

  // Initialize cplBC for RCR-type BC
  //
  bool has_rcr_bc = false;
  for (auto& bc : lEq.bc) {
    if (utils::btest(bc.bType, enum_int(BoundaryConditionType::bType_RCR))) { 
      has_rcr_bc = true;
      break;
    }
  }

  if (has_rcr_bc) { 
    if (std::set<EquationType>{Equation_fluid,Equation_FSI,Equation_CMM}.count(lEq.phys) == 0) {
      throw std::runtime_error("RCR-type BC is allowed for fluid/CMM/FSI eq. only.");
    }
    cplBC.schm = CplBCType::cplBC_SI;
    if (lEq.useTLS) {
      cplBC.schm = CplBCType::cplBC_E;
    }
    cplBC.nX = cplBC.nFa;
    cplBC.nXp = cplBC.nFa + 1;

    if (cplBC.xo.size() != 0) {
      throw std::runtime_error("cplBC structure is already initialized. Unexpected behavior.");
    }
    cplBC.xo.resize(cplBC.nX); 
    cplBC.xp.resize(cplBC.nXp);
    cplBC.saveName = chnl_mod.appPath + "/RCR.dat";
    cplBC.initRCR = eq_params->initialize_rcr_from_flow.value();
  }

  // Searching for body forces
  //
  lEq.nBf = eq_params->body_forces.size();
  if (lEq.nBf != 0) { 
    lEq.bf.resize(lEq.nBf);
  }

  for (int iBf = 0; iBf < lEq.nBf; iBf++) {
    auto& bf_params = eq_params->body_forces[iBf];
    auto& bf = lEq.bf[iBf];
    auto bf_mesh_name = bf_params->mesh_name.value();
    bf.mesh_name = bf_mesh_name;

    // Find the mesh ID for the body force.
    all_fun::find_msh(com_mod.msh, bf_mesh_name, bf.iM);

    read_bf(com_mod, bf_params, bf);

    if (utils::btest(bf.bType, enum_int(BoundaryConditionType::bType_Neu)) ||
        utils::btest(bf.bType, enum_int(BoundaryConditionType::bType_trac))) {
      if (!com_mod.shlEq && !com_mod.cmmInit) {
        throw std::runtime_error("Pressure or traction load can be applied only for shells or for initializing CMM.");
      }
    }
  }
}

//---------------------------------
// read_fiber_temporal_values_file
//---------------------------------
// Read fiber reinforcement stress the values from a file.
//
// Note: There is no equivalent Fortran function.
//
void read_fiber_temporal_values_file(FiberReinforcementStressParameters& fiber_params, dmnType& lDmn)
{
  auto file_name = fiber_params.temporal_values_file_path.value();
  std::ifstream temporal_values_file;
  temporal_values_file.open(file_name);
  if (!temporal_values_file.is_open()) {
    throw std::runtime_error("Failed to open the fiber reinforcement stress the values file '" + file_name + "'.");
  }

  int i, j;
  temporal_values_file >> i >> j; 
  if ((i == 0) || (j == 0) || (i < 2)) {
    throw std::runtime_error("Error reading the first line of the fiber reinforcement stress the values file '" + file_name + "' has an incorrect format.");
  }

  lDmn.stM.Tf.gt.d = 1;
  lDmn.stM.Tf.gt.n = j;

  // Read time/value pairs.
  //
  std::vector<std::vector<double>> temporal_values;
  double time, value;
  std::string line;

  while (std::getline(temporal_values_file, line)) { 
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
    if (line == "") {
      continue;
    }
    std::istringstream line_input(line);
    std::vector<double> values;
    while (line_input >> value) {
      values.push_back(value);
    }
    temporal_values.push_back(values);
  }

  if (lDmn.stM.Tf.gt.lrmp) {
    lDmn.stM.Tf.gt.n = 1;
  }

  lDmn.stM.Tf.gt.qi.resize(lDmn.stM.Tf.gt.d);
  lDmn.stM.Tf.gt.qs.resize(lDmn.stM.Tf.gt.d);
  lDmn.stM.Tf.gt.r.resize(lDmn.stM.Tf.gt.d,lDmn.stM.Tf.gt.n);
  lDmn.stM.Tf.gt.i.resize(lDmn.stM.Tf.gt.d,lDmn.stM.Tf.gt.n);

  fft(i, temporal_values, lDmn.stM.Tf.gt);
}

//------------
// read_files
//------------
// Read in solver parameter, mesh, face, etc. files.
//
// This replicates the READFILES subroutine defind in READFILES.f.
//
// The many global COMMOD values set in the Fortan subroutine are set in the
// Simulation::set_module_parameters() method.
//
void read_files(Simulation* simulation, const std::string& file_name)
{
  using namespace consts;

  auto& com_mod = simulation->get_com_mod();

  #define n_debug_read_files
  #ifdef debug_read_files
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  // Read the solver XML file.
  #ifdef debug_read_files
  dmsg << "Read the solver XML file " << " ... ";
  #endif
  if (!com_mod.resetSim) {
    simulation->read_parameters(std::string(file_name));
  }

  auto& chnl_mod = simulation->get_chnl_mod();
  auto& gen_params = simulation->parameters.general_simulation_parameters;

  // Set chnl_mod.appPath.
  //
  // Note that 'com_mod.saveNam' is changed in 'Simulation::set_module_parameters()'.
  //
  if (!com_mod.resetSim) {
    if (gen_params.save_results_in_folder.defined()) {
      com_mod.saveName = gen_params.save_results_in_folder.value();
      chnl_mod.appPath = com_mod.saveName + "/";
    } else {
      com_mod.saveName = "";
      chnl_mod.appPath = std::to_string(com_mod.cm.np()) + "-procs" + "/";
    }

    // [NOTE] not implemented.
    /*
    chnl_mod.std.oTS = gen_params.verbose.value();
    chnl_mod.wrn.oTS = gen_params.warning.value();
    chnl_mod.dbg.oTS = gen_params.debug.value();

    chnl_mod.std.oTF = chnl_mod.std.oTS;
    chnl_mod.wrn.oTF = chnl_mod.wrn.oTS;
    chnl_mod.dbg.oTF = chnl_mod.dbg.oTS;
    */

    com_mod.stFileFlag = gen_params.continue_previous_simulation.value();
  }

  // Setup logging to history file.
  if (!simulation->com_mod.cm.slv(simulation->cm_mod)) {
    std::string hist_file_name;

    if (chnl_mod.appPath != "") { 
      auto mkdir_arg = std::string("mkdir -p ") + chnl_mod.appPath;
      std::system(mkdir_arg.c_str());
      hist_file_name = chnl_mod.appPath + "/" + simulation->history_file_name;
    } else {
      hist_file_name = simulation->history_file_name;
    }

    bool output_to_cout = true;
    simulation->logger.initialize(hist_file_name, output_to_cout);
  }

  // Set simulation and module member data from XML parameters.
  simulation->set_module_parameters();

  // Read mesh and BCs data.
  #ifdef debug_read_files
  dmsg << "Read mesh and BCs data " << " ... ";
  #endif
  read_msh_ns::read_msh(simulation);

  // Reading immersed boundary mesh data.
  //
  // [NOTE] Not implemented.
  //

  /*
  i = list%srch("Add IB")
  if (i .GT. 0) THEN
    ibFlag = .TRUE.
    ALLOCATE(ib)
    CALL IB_READMSH(list)
    CALL IB_READOPTS(list)
  END if
  */

  // Read equations.
  // 
  int nEq = simulation->parameters.equation_parameters.size();
  com_mod.nEq = nEq; 
  com_mod.eq.resize(nEq);
  std::for_each(com_mod.eq.begin(),com_mod.eq.end(),[&](eqType& eq){eq.roInf=simulation->roInf;});
  #ifdef debug_read_files
  dmsg << "Read equations " << " ... ";
  dmsg << "nEq: " << nEq;
  #endif

  for (int iEq = 0; iEq < nEq; iEq++) { 
    auto& eq = com_mod.eq[iEq];
    auto& eq_params = simulation->parameters.equation_parameters[iEq]; 
    read_eq(simulation, eq_params, eq);

    // [TODO:DaveP] IB is not implemented.
    if (com_mod.ibFlag) {
      if ((eq.phys == EquationType::phys_fluid) || (eq.phys == EquationType::phys_FSI)) {
        //CALL IB_READEQ(eq(iEq), lPtr, ctmp)
      }
    }

    // Check equation order.
    //
    if (eq.phys == EquationType::phys_fluid) {
      if (iEq != 0) {
         throw std::runtime_error("fluid equation must come first.");
      }
    }

    if (eq.phys == EquationType::phys_CMM) { 
      if (iEq != 0) {
         throw std::runtime_error("CMM equation must come first.");
      }
      if (com_mod.mvMsh) {
        throw std::runtime_error("Both CMM and ALE cannot be solved together.");
      }
    }

    if (eq.phys == EquationType::phys_FSI) { 
      if (iEq != 0) {
         throw std::runtime_error("FSI equation must come first.");
      }
      for (int i = 0; i < nEq; i++) { 
        if (i != iEq) {
          auto& eq_params = simulation->parameters.equation_parameters[i]; 
          if (eq_params->type.value() != "mesh")  {
            throw std::runtime_error("mesh equation has to be specified after FSI equation.");
          }
        }
      }
    }

    if (com_mod.rmsh.isReqd && com_mod.eq[0].phys != EquationType::phys_FSI) {   
      throw std::runtime_error("Remeshing is applicable only for FSI equation");
    }     

    if (com_mod.iCntct) {   
      if (eq.phys != EquationType::phys_shell) {
        throw std::runtime_error( "Contact model is applicable for shell problems only");
      }
      if (com_mod.nMsh == 0) {
        throw std::runtime_error("More than one mesh is needed to apply contact model");
      }
    }     

    if (eq.phys == EquationType::phys_heatF) {   
      auto& eq1_params = simulation->parameters.equation_parameters[0];
      auto eq1_type = eq1_params->type.value();
      if ((eq1_type != "fluid") && (eq1_type != "FSI") && (!com_mod.usePrecomp)) {
        throw std::runtime_error("heatF equation has to be specified after fluid/FSI equation");
      }
    }     

    if (eq.phys == EquationType::phys_mesh) {   
      if (!com_mod.mvMsh) {
        throw std::runtime_error("mesh equation can only be specified after FSI equation");
      }     
    }     
  }
  #ifdef debug_read_files
  dmsg << "Done Read equations " << " ";
  #endif

  auto& cep_mod = simulation->get_cep_mod();

  // Read equation-specific parameters for a cardiac electrophysiology model (ECG leads only, at the moment)
  if (cep_mod.cepEq) {
    for (int iEq = 0; iEq < nEq; iEq++) {
      auto& eq = com_mod.eq[iEq];

      if (eq.phys == EquationType::phys_CEP) {
        auto& eq_params = simulation->parameters.equation_parameters[iEq];
        read_cep_equation(&cep_mod, simulation, eq_params);
      }
    }
  }

  if (cep_mod.cem.cpld) {
    if (nEq == 1) {
      throw std::runtime_error("Min equations (2) not solved for electro-mechanics coupling");
    }

    int i = 0;
    for (int iEq = 0; iEq < nEq; iEq++) {
      auto& eq = com_mod.eq[iEq];
      if ((eq.phys == EquationType::phys_CEP) || (eq.phys == EquationType::phys_struct) || 
          (eq.phys == EquationType::phys_ustruct)) {
        i = i + 1;
      }
    }

    if (i != 2) {
      throw std::runtime_error("Both electrophysiology and struct have to be solved for electro-mechanics");
    }

    if (cep_mod.cem.aStress &&  cep_mod.cem.aStrain) {
      throw std::runtime_error("Cannot set both active strain and active stress coupling");
    }

    if (cep_mod.cem.aStrain) {
      if (com_mod.nsd != 3) {
        throw std::runtime_error("Active strain coupling is allowed only for 3D bodies");
      }

      for (int iEq = 0; iEq < nEq; iEq++) {
        auto& eq = com_mod.eq[iEq];
        for (int i = 0; i < eq.nDmn; i++) {
          auto& dmn = eq.dmn[i];
          if ((dmn.phys != EquationType::phys_ustruct) && (dmn.phys != EquationType::phys_struct)) { 
            continue; 
          }
          if ((dmn.stM.isoType != ConstitutiveModelType::stIso_HO)) {
            throw std::runtime_error("Active strain is allowed with Holzapfel-Ogden passive constitutive model only");
          }
        }
      }
    }
  }

  // [NOTE] what's going on here?
  if (com_mod.cplBC.xo.size() == 0) {
    com_mod.cplBC.nX = 0;
    com_mod.cplBC.xo.clear();
  }

  #ifdef debug_read_files
  dmsg << "Done" << " ";
  #endif
}

//--------------------------------
// read_fourier_coeff_values_file 
//--------------------------------
// Set boundary condition Fourier coefficients read in from a file. 
//
// Note: There is no equivalent Fortran function.
//
// [TODO:DaveP] this is not fully implemented. 
//
void read_fourier_coeff_values_file(const std::string& file_name, bcType& lBc) 
{
  std::ifstream temporal_values_file;
  temporal_values_file.open(file_name);
  if (!temporal_values_file.is_open()) {
    throw std::runtime_error("Failed to open the Fourier coefficients values file '" + file_name + "'.");
  }

  temporal_values_file >> lBc.gt.ti >> lBc.gt.T; 

  // Read time/value pairs.
  //
  lBc.gt.qi.resize(lBc.gt.d); 
  lBc.gt.qs.resize(lBc.gt.d);

  double time, value;
  std::string line;
  int n = 0;

  while (std::getline(temporal_values_file, line)) { 
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
    if (line == "") {
      continue;
    }
    std::istringstream line_input(line);
    std::vector<double> values;
    while (line_input >> value) {
      values.push_back(value);
    }
    lBc.gt.qi[n] = values[0];
    lBc.gt.qs[n] = values[1];
    n += 1;

    if (n == lBc.gt.d) {
      break;
    }
  }

  temporal_values_file >> lBc.gt.n; 

  lBc.gt.r.resize(lBc.gt.d, lBc.gt.n);
  lBc.gt.i.resize(lBc.gt.d, lBc.gt.n);
  
  int j = 0;

  while (std::getline(temporal_values_file, line)) { 
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
    if (line == "") {
      continue;
    }
    std::istringstream line_input(line);
    std::vector<double> values;
    while (line_input >> value) {
      values.push_back(value);
    }

    for (int i = 0; i < lBc.gt.d; i++) { 
      lBc.gt.r(i,j) = values[i];
      lBc.gt.i(i,j) = values[i + lBc.gt.d];
    }

    j += 1;
    if (j == lBc.gt.n) {
      break;
    }
  }
}

void read_fourier_coeff_values_file(const std::string& file_name, bfType& lBf) 
{
  std::ifstream temporal_values_file;
  temporal_values_file.open(file_name);
  if (!temporal_values_file.is_open()) {
    throw std::runtime_error("Failed to open the Fourier coefficients values file '" + file_name + "'.");
  }

  temporal_values_file >> lBf.bt.ti >> lBf.bt.T; 

  // Read time/value pairs.
  //
  lBf.bt.qi.resize(lBf.bt.d); 
  lBf.bt.qs.resize(lBf.bt.d);

  double time, value;
  std::string line;
  int n = 0;

  while (std::getline(temporal_values_file, line)) { 
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
    if (line == "") {
      continue;
    }
    std::istringstream line_input(line);
    std::vector<double> values;
    while (line_input >> value) {
      values.push_back(value);
    }
    lBf.bt.qi[n] = values[0];
    lBf.bt.qs[n] = values[1];
    n += 1;

    if (n == lBf.bt.d) {
      break;
    }
  }

  temporal_values_file >> lBf.bt.n; 

  lBf.bt.r.resize(lBf.bt.d, lBf.bt.n);
  lBf.bt.i.resize(lBf.bt.d, lBf.bt.n);
  
  int j = 0;

  while (std::getline(temporal_values_file, line)) { 
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
    if (line == "") {
      continue;
    }
    std::istringstream line_input(line);
    std::vector<double> values;
    while (line_input >> value) {
      values.push_back(value);
    }

    for (int i = 0; i < lBf.bt.d; i++) { 
      lBf.bt.r(i,j) = values[i];
      lBf.bt.i(i,j) = values[i + lBf.bt.d];
    }

    j += 1;
    if (j == lBf.bt.n) {
      break;
    }
  }
}

//---------
// read_ls
//---------
// Set solver parameters.
//
void read_ls(Simulation* simulation, EquationParameters* eq_params, consts::SolverType ils_type, eqType& lEq)
{
  using namespace consts;
  using namespace fsi_linear_solver;

  #define n_debug_read_ls
  #ifdef debug_read_ls
  DebugMsg dmsg(__func__, simulation->com_mod.cm.idcm());
  dmsg.banner();
  #endif

  // Map SolverType to LinearSolverType enums.
  static std::map<SolverType,LinearSolverType> solver_to_ls_map = {
    {SolverType::lSolver_NS, LinearSolverType::LS_TYPE_NS},
    {SolverType::lSolver_GMRES, LinearSolverType::LS_TYPE_GMRES},
    {SolverType::lSolver_CG, LinearSolverType::LS_TYPE_CG},
    {SolverType::lSolver_BICGS, LinearSolverType::LS_TYPE_BICGS},
  };

  // Get solver type.
  //
  SolverType solver_type;
  LinearSolverType FSILSType;
  bool solver_type_defined = eq_params->linear_solver.type.defined();
  #ifdef debug_read_ls
  dmsg << "solver_type_defined: " << solver_type_defined;
  #endif

  if (solver_type_defined) {
    auto solver_str = eq_params->linear_solver.type.value();
    std::transform(solver_str.begin(), solver_str.end(), solver_str.begin(), ::tolower);
    try {
     solver_type = solver_name_to_type.at(solver_str);
    } catch (const std::out_of_range& exception) {
      throw std::runtime_error("Unknown solver type '" + solver_str + ".");
    }
  } else {
    solver_type = ils_type;
  }

  FSILSType = solver_to_ls_map[solver_type];
  for (auto entry : solver_name_to_type) {
    if (solver_type == entry.second) {
      #ifdef debug_read_ls
      dmsg << "solver_type: " << entry.first;
      #endif
      break;
    }
  }

  // Set linear solver parameters.
  //
  lEq.ls.LS_type = solver_type;
  fsi_linear_solver::fsils_ls_create(lEq.FSILS, FSILSType);

  // Process linear_algebra parameters.
  //
  auto& linear_algebra = eq_params->linear_solver.linear_algebra;
  #ifdef debug_read_ls
  dmsg << "linear_algebra.defined: " << linear_algebra.defined();
  dmsg << "linear_algebra.type: " << linear_algebra.type();
  dmsg << "linear_algebra.preconditioner: " << linear_algebra.preconditioner();
  #endif

  if (!linear_algebra.defined()) {
    throw std::runtime_error("[svMultiPhysics] ERROR: No <Linear_algebra> section has been defined for equation '" + 
        eq_params->type() + ".");
  }

  lEq.linear_algebra_type = LinearAlgebra::name_to_type.at(linear_algebra.type());
  auto prec_type = consts::preconditioner_name_to_type.at(linear_algebra.preconditioner());
  lEq.linear_algebra_preconditioner = consts::preconditioner_name_to_type.at(linear_algebra.preconditioner());
  lEq.linear_algebra_assembly_type = LinearAlgebra::name_to_type.at(linear_algebra.assembly()); 

  // Check that equation physics is compatible with the LinearAlgebra type. 
  for (auto& domain : lEq.dmn) {
    LinearAlgebra::check_equation_compatibility(domain.phys,  lEq.linear_algebra_type, lEq.linear_algebra_assembly_type);
  }

  if (!solver_type_defined) {
    return;
  } 

  auto& linear_solver = eq_params->linear_solver;

  lEq.ls.mItr = linear_solver.max_iterations.value();
  lEq.FSILS.RI.mItr = lEq.ls.mItr;

  if (linear_solver.tolerance.defined()) {
    lEq.FSILS.RI.relTol = linear_solver.tolerance.value();
  }

  if (linear_solver.absolute_tolerance.defined()) {
    lEq.FSILS.RI.absTol = linear_solver.absolute_tolerance.value(); 
  }

  if (linear_solver.krylov_space_dimension.defined()) {
    lEq.FSILS.RI.sD = linear_solver.krylov_space_dimension.value();
  }

  if (solver_type == SolverType::lSolver_NS) {
    lEq.FSILS.GM.mItr = linear_solver.ns_gm_max_iterations.value();
    lEq.FSILS.CG.mItr = linear_solver.ns_cg_max_iterations.value(); 

    lEq.FSILS.RI.relTol = linear_solver.tolerance.value();
    lEq.FSILS.GM.relTol = linear_solver.ns_gm_tolerance.value(); 
    lEq.FSILS.CG.relTol = linear_solver.ns_cg_tolerance.value(); 

    lEq.FSILS.RI.absTol = linear_solver.absolute_tolerance.value();
    lEq.FSILS.GM.absTol = lEq.FSILS.RI.absTol;
    lEq.FSILS.CG.absTol = lEq.FSILS.RI.absTol;

    lEq.FSILS.GM.sD = lEq.FSILS.RI.sD;
  } 

  #ifdef debug_read_ls
  dmsg << "Done " << " ";
  #endif
}

//----------------
// read_mat_model
//----------------
// Set the material model parameters for the given domain.
//
void read_mat_model(Simulation* simulation, EquationParameters* eq_params, DomainParameters* domain_params, dmnType& lDmn)
{ 
  using namespace consts;

  // Domain properties: elasticity modulus, poisson ratio
  double E = lDmn.prop[PhysicalProperyType::elasticity_modulus];
  double nu = lDmn.prop[PhysicalProperyType::poisson_ratio];

  // Shear modulus
  double mu  = 0.5 * E / (1.0 + nu);

  // Incompressible material
  bool incompFlag = false;
  if (utils::is_zero(nu-0.5)) {
    incompFlag = true;
  }

  // Bulk modulus for compressible case
  double kap  = {0.0};
  double lam  = {0.0};
  if (!incompFlag) {
    kap = E / (1.0 - 2.0*nu) / 3.0;
    lam = E * nu / (1.0+nu) / (1.0 - 2.0*nu);
  } else { 
    kap = 0.0;
  }

  // If no constitutive model was given use a NeoHookean model.
  //
  ConstitutiveModelType cmodel_type;
  std::string cmodel_str;

  if (!domain_params->constitutive_model.defined()) { 
    lDmn.stM.isoType = ConstitutiveModelType::stIso_nHook;
    cmodel_type = ConstitutiveModelType::stIso_nHook;
    cmodel_str = "neoHookean";

  // Get the constitutive model type.
  } else {
    cmodel_str = domain_params->constitutive_model.type.value();
    try {
      cmodel_type = constitutive_model_name_to_type.at(cmodel_str);
    } catch (const std::out_of_range& exception) {
      throw std::runtime_error("Unknown constitutive model type '" + cmodel_str + ".");
    }
  }

  // Set material properties for the domain 'lDmn'.
  try {
    set_material_props[cmodel_type](domain_params, mu, kap, lam, lDmn);
   } catch (const std::bad_function_call& exception) {
    throw std::runtime_error("[read_mat_model] Constitutive model type '" + cmodel_str + "' is not implemented.");
  }

  // Set fiber reinforcement stress.
  if (domain_params->fiber_reinforcement_stress.defined()) { 
    auto& fiber_params = domain_params->fiber_reinforcement_stress;
    auto fiber_stress = fiber_params.type.value();
    std::transform(fiber_stress.begin(), fiber_stress.end(), fiber_stress.begin(), ::tolower);

    if (fiber_stress == "steady") {
      lDmn.stM.Tf.fType = utils::ibset(lDmn.stM.Tf.fType, static_cast<int>(BoundaryConditionType::bType_std));
      lDmn.stM.Tf.g = fiber_params.value.value();

    } else if (fiber_stress == "unsteady") {
      lDmn.stM.Tf.fType = utils::ibset(lDmn.stM.Tf.fType, static_cast<int>(BoundaryConditionType::bType_ustd));
      lDmn.stM.Tf.gt.lrmp = fiber_params.ramp_function.value();
      read_fiber_temporal_values_file(fiber_params, lDmn);
    }

    // Read directional stress distribution parameters
    if (fiber_params.directional_distribution.defined()) {
      // Validate: ensures exactly 3 parameters specified (no empty blocks), sums to 1.0, non-negative
      fiber_params.directional_distribution.validate();
      
      // Read the validated values (validate() ensures all three are defined)
      lDmn.stM.Tf.eta_f = fiber_params.directional_distribution.fiber_direction.value();
      lDmn.stM.Tf.eta_s = fiber_params.directional_distribution.sheet_direction.value();
      lDmn.stM.Tf.eta_n = fiber_params.directional_distribution.sheet_normal_direction.value();
    }
    // Otherwise (no block at all), defaults (eta_f=1.0, eta_s=0.0, eta_n=0.0) from ComMod.h are used
  }

  // Check for shell model
  //
  // ST91 is the default and the only dilational penalty model for
  // compressible shell elements. This is set to avoid any square-
  // root evaulations of the Jacobian during Newton iterations for
  // satisfying plane-stress condition.
  //
  if (lDmn.phys == EquationType::phys_shell) {
    lDmn.stM.Kpen = kap;
    if (!incompFlag) {
      lDmn.stM.volType = ConstitutiveModelType::stVol_ST91;
    }
    return;
  }

  // Look for dilational penalty model. HGO uses quadratic penalty model.
  if (domain_params->dilational_penalty_model.defined()) {
    auto model_str = domain_params->dilational_penalty_model.value();
    try {
      lDmn.stM.volType = constitutive_model_name_to_type.at(model_str);
    } catch (const std::out_of_range& exception) {
      throw std::runtime_error("Unknown constitutive model type '" + model_str + ".");
    }

  } else {
    lDmn.stM.volType = ConstitutiveModelType::stVol_ST91;
  }

  // Default penalty parameter is equal to bulk modulus.
  if (domain_params->penalty_parameter.defined()) {
    lDmn.stM.Kpen = domain_params->penalty_parameter.value();
  } else {
    lDmn.stM.Kpen = kap;
  }
}

//--------------
// read_outputs
//--------------
// Set output parameters. 
//
// nDOP(1) is the total number of outputs 
//
// nDOP(2) is the default number of outputs for VTK file, 
//
// nDOP(3) is for boundaries, and
//
// nDOP(4) is for volume
//
void read_outputs(Simulation* simulation, EquationParameters* eq_params, eqType& lEq, EquationNdop& nDOP,  EquationOutputs& outPuts)
{
  using namespace consts;
  auto& com_mod = simulation->get_com_mod();
  int nsd = com_mod.nsd;
  #include "set_output_props.h"

  lEq.nOutput = nDOP[0];
  lEq.output.resize(nDOP[0]);

  for (int iOut = 0; iOut < nDOP[0]; iOut++) {
    auto& output = lEq.output[iOut]; 
    auto output_type = outPuts.at(iOut);
    try {
      std::tie(output.grp, output.o, output.l, output.name) = output_props_map.at(output_type);
    } catch (const std::out_of_range& exception) { 
      throw std::runtime_error("Unsupporteed output type '" + std::to_string(static_cast<int>(output_type)) + "'.");
    }
  }

  // First reading the outputs for VTK files and then for boundaries and last for the volume
  //
  int nOut = eq_params->outputs.size();
  OutputType output_type;

  for (int iOut = 0; iOut < nOut; iOut++) { 
    auto& output_params = eq_params->outputs[iOut];
    auto output_type_str = output_params->type.value();
    if (output_type_str == "Alias") { 
      continue;
    }

    try {
      output_type = output_type_name_to_type.at(output_type_str);
    } catch (const std::out_of_range& exception) { 
      throw std::runtime_error("Unknown output type '" + output_type_str + "'.");
    }

    auto& output_list = output_params->output_list;

    for (int i = 0; i < lEq.nOutput; i++) {
      lEq.output[i].options.set_option(output_type, output_params->get_output_value(lEq.output[i].name));
      if (lEq.output[i].name == "Vortex") {
        if (nsd != maxNSD) {
          lEq.output[i].options.set_option(output_type, false);
        }
      }
    }
  }

  // Read any alias names for outputs
  //
  for (int iOut = 0; iOut < nOut; iOut++) {
    auto& output_params = eq_params->outputs[iOut];
    auto output_type = output_params->type.value();
    if (output_type != "Alias") {
      continue;
    }
    for (int i = 0; i < lEq.nOutput; i++) {
      auto alias_name = output_params->get_alias_value(lEq.output[i].name);
      if (alias_name.size() != 0) { 
        lEq.output[i].name = alias_name;
      }
    }
  }
}

//---------------------
// read_spatial_values
//---------------------
//
// Read in a file containing spatial values used for a boundary condition.
//
// There is no equivalent Fortran subroutine.
//
// [TODO:DaveP] this is not tested.
//
void read_spatial_values(const ComMod& com_mod, const mshType& msh, const faceType& lFa, 
    const std::string& file_name, bcType& lBc)
{
  std::ifstream file_stream;
  file_stream.open(file_name);
  if (!file_stream.is_open()) {
    throw std::runtime_error("Failed to open the spatial values file '" + file_name + "'.");
  }

  lBc.gx.resize(lFa.nNo); 

  // Preparing the pointer array
  //
  Vector<int> ptr(msh.gnNo);
  ptr = -1;
  for (int a = 0; a < lFa.nNo; a++) {
    lBc.gx[a] = 0.0;
    int Ac = lFa.gN[a];
    Ac = msh.lN[Ac];
    if (Ac == -1) {
      throw std::runtime_error("Incorrect global node number detected for BC for mesh '" + msh.name + 
          "' and face '" + lFa.name + "' for node " + std::to_string(a) + ".");
    }
    ptr[Ac] = a;
  }

  for (int b = 0; b < lFa.nNo; b++) {
    int Ac;
    double rtmp;
    file_stream >> Ac >> rtmp;
    if ((Ac >= msh.gnNo) || (Ac < 0)) {
      throw std::runtime_error("The node number " + std::to_string(Ac) + 
            " in the spatial values file '" + file_name + " is larger than the number of nodes in the mesh.");
    }
    int a = ptr[Ac];
    if (a == -1) {
      throw std::runtime_error("The node number " + std::to_string(Ac) + 
            " from the spatial values file '" + file_name + " does not belong to the face '" + lFa.name + "'."); 
    }     

    lBc.gx[a] = rtmp;
  }
}

//-----------------------
// read_temp_spat_values 
//-----------------------
// Read in a file containing temporal and spatial values
// used for a boundary condition.
//
// There is no equivalent Fortran subroutine.
//
void read_temp_spat_values(const ComMod& com_mod, const mshType& msh, const faceType& lFa, 
    const std::string& file_name, bcType& lBc)
{
  std::ifstream file_stream;
  file_stream.open(file_name);
  if (!file_stream.is_open()) {
    throw std::runtime_error("Failed to open the temporal and spatial values file '" + file_name + "'.");
  }

  int ndof, num_ts, num_nodes;
  file_stream >> ndof >> num_ts >> num_nodes;
  if ((ndof == 0) || (num_ts == 0) || (num_nodes == 0)) {
    throw std::runtime_error("Error reading the first line of the temporal and spatial values file '" + file_name + "'.");
  }

  if (num_nodes != lFa.nNo) {
    throw std::runtime_error("The number of nodes (" + std::to_string(num_nodes) + ") in the temporal and spatial values file '" + 
        file_name + "' are not equal to the number of nodes (" + std::to_string(lFa.nNo) + ") in face '" + lFa.name + "'");
  }

  if ((ndof < 1) || (ndof > com_mod.nsd)) { 
    throw std::runtime_error("The number of degrees of freedom (" + std::to_string(ndof) + ") in the temporal and spatial values file '" +
        file_name + "' don't agree with the number of degrees of freedom (" + std::to_string(com_mod.nsd) + ") of the simulation.");
  }

  lBc.gm.t.resize(num_ts); 
  lBc.gm.d.resize(ndof,num_nodes,num_ts); 
  lBc.gm.dof = ndof;
  lBc.gm.nTP = num_ts;

  Vector<int> ptr(msh.gnNo);
  ptr = -1;

  // Preparing the pointer array
  for (int a = 0; a < lFa.nNo; a++) {
    int Ac = lFa.gN[a];
    Ac = msh.lN[Ac];
    if (Ac == -1) {
      throw std::runtime_error("Incorrect global node number detected for BC for mesh '" + msh.name + 
          "' and face '" + lFa.name + "' for node " + std::to_string(a) + ".");
    }
    ptr[Ac] = a;
  }

  // Read time sequence.
  //
  for (int i = 0; i < num_ts; i++) {
    double rtmp;
    file_stream >> rtmp;
    lBc.gm.t[i] = rtmp;

    if (i == 0) {
      if (!utils::is_zero(rtmp)) { 
        throw std::runtime_error("The first time step (" + std::to_string(rtmp) + 
            ") in the temporal and spatial values file '" + file_name + " is not zero.");
      }

    } else { 
      rtmp = rtmp - lBc.gm.t[i-1];
      if (utils::is_zero(rtmp) || (rtmp < 0.0)) { 
        throw std::runtime_error("A non-increasing time step was found in the temporal and spatial values file '" + 
            file_name + ".");
      }
    }
  }

  lBc.gm.period = lBc.gm.t[num_ts-1];

  // Read in data.
  //
  // Note: This file contains node IDs so be careful
  // to subbtract 1 from them.
  //
  for (int b = 0; b < lFa.nNo; b++) {
    int Ac;
    file_stream >> Ac;
    Ac -= 1;

    if ((Ac >= msh.gnNo) || (Ac < 0)) {
      throw std::runtime_error("The node number " + std::to_string(Ac) + 
            " in the temporal and spatial values file '" + file_name + " is larger than the number of nodes in the mesh.");
    }     

    int a = ptr[Ac];
    if (a == -1) {
      throw std::runtime_error("The node number " + std::to_string(Ac) + 
            " from the temporal and spatial values file '" + file_name + " does not belong to the face '" + lFa.name + "'."); 
    }     

    for (int i = 0; i < num_ts; i++) { 
      double value;
      for (int k = 0; k < ndof; k++) { 
        file_stream >> value;
        lBc.gm.d(k,a,i) = value;
      }
    } 
  } 
}

//-----------------------
// read_temp_spat_values 
//-----------------------
// Read body force temporal and spatial data.
//
void read_temp_spat_values(const ComMod& com_mod, const mshType& msh, const std::string& file_name, bfType& lBf)
{
  #define n_debug_read_ts_values_bf
  #ifdef debug_read_ts_values_bf 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "msh.gnNo: " << msh.gnNo;
  dmsg << "lBf.dof: " << lBf.dof;
  dmsg << "file_name: " << file_name;
  #endif
  lBf.file_name = file_name;

  std::ifstream file_stream;
  file_stream.open(file_name);
  if (!file_stream.is_open()) {
    throw std::runtime_error("Failed to open the temporal and spatial values file '" + file_name + "'.");
  }

  // Read number dof (dimension), number of time steps and the number of nodes.
  //
  int ndof, num_ts, num_nodes;
  file_stream >> ndof >> num_ts >> num_nodes;
  if ((ndof == 0) || (num_ts == 0) || (num_nodes == 0)) {
    throw std::runtime_error("Error reading the first line of the temporal and spatial values file '" + file_name + "'.");
  }
  #ifdef debug_read_ts_values_bf 
  dmsg << "ndof: " << ndof;
  dmsg << "num_ts: " << num_ts;
  dmsg << "num_nodes: " << num_nodes;
  #endif

  if (num_nodes > msh.gnNo) {
    throw std::runtime_error("The number of nodes (" + std::to_string(num_nodes) + ") in the temporal and spatial values file '" + 
        file_name + "' are not equal to the number of nodes (" + std::to_string(msh.gnNo) + ") in face '" + msh.name + "'");
  }

  if (ndof != lBf.dof) { 
    throw std::runtime_error("The number of degrees of freedom (" + std::to_string(ndof) + ") in the temporal and spatial values file '" +
        file_name + "' don't agree with the number of degrees of freedom (" + std::to_string(lBf.dof) + ") of the simulation.");
  }

  lBf.bm.dof = lBf.dof;
  lBf.bm.nTP = num_ts;

  lBf.bm.t.resize(num_ts); 
  lBf.bm.d.resize(lBf.dof, com_mod.gtnNo, num_ts);

  // Read time sequence.
  //
  #ifdef debug_read_ts_values_bf 
  dmsg << "Read time sequence ...";
  #endif
  for (int i = 0; i < lBf.bm.nTP; i++) {
    double rtmp;
    file_stream >> rtmp;
    #ifdef debug_read_ts_values_bf 
    dmsg << "----- i " << i << " -----";
    dmsg << "rtmp: " << rtmp;
    #endif
    lBf.bm.t[i] = rtmp;
    if (i == 0) {
      if (!utils::is_zero(rtmp)) { 
        throw std::runtime_error("The first time step (" + std::to_string(rtmp) + 
            ") in the temporal and spatial values file '" + file_name + " is not zero.");
      }
    } else { 
      rtmp = rtmp - lBf.bm.t[i-1];
      if (utils::is_zero(rtmp) || (rtmp < 0.0)) { 
        throw std::runtime_error("A non-increasing time step was found in the temporal and spatial values file '" + 
            file_name + ".");
      }
    }
  }

  lBf.bm.period = lBf.bm.t[lBf.bm.nTP-1];

  // Read in data.
  //
  // Note: This file contains node IDs so be careful
  // to subbtract 1 from them.
  //
  #ifdef debug_read_ts_values_bf 
  dmsg << "Read data ...";
  #endif
  for (int i = 0; i < msh.gnNo; i++) {
    int Ac;
    file_stream >> Ac;
    //dmsg << "Node : " << Ac;
    Ac -= 1;
    if ((Ac >= msh.gnNo) || (Ac < 0)) {
      throw std::runtime_error("The node number " + std::to_string(Ac) + 
            " in the temporal and spatial values file '" + file_name + " is larger than the number of nodes in the mesh.");
    }     

    Ac = msh.gN(Ac);

    for (int j = 0; j < lBf.bm.nTP; j++) { 
      double value;
      for (int k = 0; k < lBf.bm.dof; k++) { 
        file_stream >> value;
        lBf.bm.d(k,Ac,j) = value;
      }
    } 
  } 
}

//----------------------
// read_temporal_values
//----------------------
// Set boundary condition temporal values read in from a file. 
//
// Data modified:
//
//  Set lBc.gt Fourier coefficients data (fcType) 
//     lBc.gt.n 
//     lBc.gt.qi
//     lBc.gt.qs
//     lBc.gt.r
//     lBc.gt.i
//
// Note: There is no equivalent Fortran function.
//
void read_temporal_values(const std::string& file_name, bcType& lBc) 
{
  std::ifstream temporal_values_file;
  temporal_values_file.open(file_name);
  if (!temporal_values_file.is_open()) {
    throw std::runtime_error("Failed to open the temporal values file '" + file_name + "'.");
  }

  int i, j;
  temporal_values_file >> i >> j; 
  if ((i == 0) || (j == 0) || (i < 2)) {
    throw std::runtime_error("Error reading the first line of the temporal values file '" + file_name + "'.");
  }
  lBc.gt.n = j;

  // Read time/value pairs.
  //
  std::vector<std::vector<double>> temporal_values;
  double time, value;
  std::string line;
  int line_number = 1;
  int num_values_per_line = lBc.gt.d + 1;

  while (std::getline(temporal_values_file, line)) { 
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
    if (line == "") {
      continue;
    }
    // Remove leading and trailing spaces.
    auto cleaned_line = std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1");
    std::istringstream line_input(cleaned_line);
    std::vector<double> values;

    while (!line_input.eof()) {
      line_input >> value;

      if (line_input.fail()) { 
        throw std::runtime_error("Error reading values for the temporal values file '" + file_name + "' for line " +
            std::to_string(line_number) + ": '" + line + "'; value number " + std::to_string(values.size()+1) + " is not a double.");
      }
      values.push_back(value);
    }

    if (values.size() != num_values_per_line) { 
      throw std::runtime_error("Error reading values for the temporal values file '" + file_name + "' for line " +
          std::to_string(line_number) + ": '" + line + "'; expected " + std::to_string(num_values_per_line) + " values per line.");
    }

    temporal_values.push_back(values);
    line_number += 1;
  }

  if (lBc.gt.lrmp) {
    lBc.gt.n = 1;
  }

  lBc.gt.qi.resize(lBc.gt.d); 
  lBc.gt.qs.resize(lBc.gt.d);
  lBc.gt.r.resize(lBc.gt.d,lBc.gt.n);
  lBc.gt.i.resize(lBc.gt.d,lBc.gt.n);

  fft(i, temporal_values, lBc.gt);
}

// [NOTE] This is a hack, should really just have a single function for this..
//
void read_temporal_values(const std::string& file_name, bfType& lBf) 
{
  std::ifstream temporal_values_file;
  temporal_values_file.open(file_name);
  if (!temporal_values_file.is_open()) {
    throw std::runtime_error("Failed to open the temporal values file '" + file_name + "'.");
  }

  int i, j;
  temporal_values_file >> i >> j; 
  if (i < 2) {
    throw std::runtime_error("The temporal values file '" + file_name + "' has an incorrect format.");
  }
  lBf.bt.n = j;

  // Read time/value pairs.
  //
  std::vector<std::vector<double>> temporal_values;
  double time, value;
  std::string line;

  while (std::getline(temporal_values_file, line)) { 
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
    if (line == "") {
      continue;
    }
    std::istringstream line_input(line);
    std::vector<double> values;
    while (line_input >> value) {
      values.push_back(value);
    }
    temporal_values.push_back(values);
  }

  if (lBf.bt.lrmp) {
    lBf.bt.n = 1;
  }

  lBf.bt.qi.resize(lBf.bt.d); 
  lBf.bt.qs.resize(lBf.bt.d);
  lBf.bt.r.resize(lBf.bt.d, lBf.bt.n);
  lBf.bt.i.resize(lBf.bt.d, lBf.bt.n);

  fft(i, temporal_values, lBf.bt);
}

//----------------
// read_trac_bcff
//----------------
// Reads pressure/traction data from a vtk vtp or vtu file and stores in moving BC data structure.
//
// Reproduces 'SUBROUTINE READTRACBCFF(lMB, lFa, fName)' defined in READFILES.f.
//
void read_trac_bcff(ComMod& com_mod, MBType& lMB, faceType& lFa, const std::string& fName)
{
  if (FILE *file = fopen(fName.c_str(), "r")) {
      fclose(file);
  } else {
    throw std::runtime_error("The VTK VTP traction data file '" + fName + "' can't be read.");
  }

  // Read the vtk file.
  //
  auto vtk_data = VtkData::create_reader(fName);
  int num_points = vtk_data->num_points();
  if (num_points == 0) {
    throw std::runtime_error("The VTK VTP traction data file '" + fName + "' does not contain any points.");
  }

  int num_elems = vtk_data->num_elems();
  if (num_elems == 0) {
    throw std::runtime_error("The VTK VTP traction data file '" + fName + "' does not contain any elements.");
  }

  // Set the name of the data to get from the file.
  //
  std::string data_name;
  if (lMB.dof == 1) {
    data_name = "Pressure";
  } else { 
    data_name = "Traction";
  }

  // Check that the vtk file has traction data.
  if (!vtk_data->has_point_data(data_name)) {
    throw std::runtime_error("No PointData DataArray named '" + data_name + "' found in the VTK VTP traction data file '" + fName +
        "' for the '" + lFa.name + "' face.");
  }

  Vector<double> tmpX1(num_points); 
  Array<double> tmpX2(consts::maxNSD, num_points);

  // Set coordinates.
  //
  faceType gFa;
  gFa.nNo = num_points;
  gFa.x.resize(com_mod.nsd, gFa.nNo);
  vtk_data->copy_points(tmpX2);

  for (int i = 0; i < num_points; i++) {
    for (int j = 0; j < com_mod.nsd; j++) {
      gFa.x(j, i) = tmpX2(j, i); 
    }
  }

  if (data_name == "Pressure") {
    vtk_data->copy_point_data(data_name, tmpX1);
  } else { 
    vtk_data->copy_point_data(data_name, tmpX2);
  }

  // Project traction from gFa to lFa. First prepare lFa%x, lFa%IEN
  //
  lFa.x.resize(com_mod.nsd, lFa.nNo);

  for (int a = 0; a < lFa.nNo; a++) {
    int Ac = lFa.gN[a];
    lFa.x.set_col(a, com_mod.x.col(Ac));
  }

  Vector<int> ptr(lFa.nNo);
  ptr = -1;
  gFa.name = "face from traction vtp";
  face_match(com_mod, lFa, gFa, ptr);

  // Copy pressure/traction data to MB data structure
  //
  if (lMB.dof == 1) {
    for (int a = 0; a < lFa.nNo; a++) {
      int Ac = ptr[a];
      lMB.d(0,a,0) = -tmpX1[Ac];
      lMB.d(0,a,1) = -tmpX1[Ac];
    }

  } else { 
    for (int a = 0; a < lFa.nNo; a++) {
      int Ac = ptr[a];
      for (int i = 0; i < lMB.dof; i++) {
        lMB.d(i,a,0) = tmpX2(i,Ac);
        lMB.d(i,a,1) = tmpX2(i,Ac);
      }
    }
  }
}

//-----------
// read_rmsh
//-----------
// Set remesher parameters.
//
// Replicates Fortran 'SUBROUTINE READRMSH(list)'.
//
void read_rmsh(Simulation* simulation, EquationParameters* eq_param)
{ 
  using namespace consts;

  #define n_debug_read_rmsh 
  #ifdef debug_read_rmsh 
  DebugMsg dmsg(__func__, simulation->com_mod.cm.idcm());
  dmsg.banner();
  #endif

  auto& remesher = eq_param->remesher;
  if (!remesher.defined()) { 
    return;
  }

  auto& com_mod = simulation->com_mod;
  auto& rmsh = com_mod.rmsh;
  int nMsh = com_mod.nMsh;

  auto mesh_gen_str = remesher.type.value();
  #ifdef debug_read_rmsh 
  dmsg << "Remesh type: " << mesh_gen_str;
  #endif

  try {
    rmsh.method = mesh_generator_name_to_type.at(mesh_gen_str);
  } catch (const std::out_of_range& exception) {
    throw std::runtime_error("Unknown mesh generator '" + mesh_gen_str + ".");
  }

  rmsh.maxEdgeSize.resize(nMsh);
  rmsh.maxEdgeSize = 0.5;

  for (int i = 0; i < nMsh; i++) {
    if (remesher.has_edge_size(com_mod.msh[i].name)) {
      rmsh.maxEdgeSize[i] = remesher.get_edge_size(com_mod.msh[i].name);
    }
    #ifdef debug_read_rmsh 
    dmsg << "mesh: " << com_mod.msh[i].name + "  edge size: " + std::to_string(rmsh.maxEdgeSize[i]);
    #endif
  }

  rmsh.minDihedAng = remesher.min_dihedral_angle.value();
  rmsh.maxRadRatio = remesher.max_radius_ratio.value();
  rmsh.freq = remesher.remesh_frequency.value();
  rmsh.cpVar = remesher.frequency_for_copying_data.value();

  #ifdef debug_read_rmsh 
  dmsg << "rmsh.minDihedAng: " << rmsh.minDihedAng; 
  dmsg << "rmsh.maxRadRatio: " << rmsh.maxRadRatio; 
  dmsg << "rmsh.freq: " << rmsh.freq;
  dmsg << "rmsh.cpVar: " << rmsh.cpVar;
  #endif
}

//-----------------
// read_fluid_visc_model
//-----------------
// Set the fluid viscosity material model parameters for the given domain.
//
void read_fluid_visc_model(Simulation* simulation, EquationParameters* eq_params, DomainParameters* domain_params, dmnType& lDmn)
{ 
  using namespace consts;

  // Get viscosity model.
  //
  FluidViscosityModelType vmodel_type;
  std::string vmodel_str; 

  if (domain_params->fluid_viscosity.model.defined()) {
    vmodel_str = domain_params->fluid_viscosity.model.value();
    std::transform(vmodel_str.begin(), vmodel_str.end(), vmodel_str.begin(), ::tolower);

    try {
      vmodel_type = fluid_viscosity_model_name_to_type.at(vmodel_str);
    } catch (const std::out_of_range& exception) {
      throw std::runtime_error("Unknown fluid viscosity model '" + vmodel_str + "'.");
    }
  } else {
    vmodel_type = FluidViscosityModelType::viscType_Const;
  }

  // Set the parameters for the given viscosity model.
  //
  auto& viscosity_params = domain_params->fluid_viscosity;

  try {
    set_fluid_viscosity_props[vmodel_type](simulation, viscosity_params, lDmn);
   } catch (const std::bad_function_call& exception) {
    throw std::runtime_error("[read_fluid_visc_model] Viscosity model '" + vmodel_str + "' is not supported.");
  }

  if ((lDmn.phys == EquationType::phys_stokes) && (lDmn.fluid_visc.viscType != FluidViscosityModelType::viscType_Const)) {
    throw std::runtime_error("Only constant viscosity is allowed for Stokes flow.");
  }
}

//-----------------
// read_solid_visc_model
//-----------------
// Set the solid viscosity material model parameters for the given domain.
//
void read_solid_visc_model(Simulation* simulation, EquationParameters* eq_params, DomainParameters* domain_params, dmnType& lDmn)
{ 
  using namespace consts;

  // Get viscosity model.
  //
  SolidViscosityModelType vmodel_type;
  std::string vmodel_str; 

  if (domain_params->solid_viscosity.model.defined()) {
    vmodel_str = domain_params->solid_viscosity.model.value();
    std::transform(vmodel_str.begin(), vmodel_str.end(), vmodel_str.begin(), ::tolower);

    try {
      vmodel_type = solid_viscosity_model_name_to_type.at(vmodel_str);
    } catch (const std::out_of_range& exception) {
      throw std::runtime_error("Unknown solid viscosity model '" + vmodel_str + "'.");
    }
  } else {
    vmodel_type = SolidViscosityModelType::viscType_Newtonian;
  }

  // Set the parameters for the given viscosity model.
  //
  auto& viscosity_params = domain_params->solid_viscosity;

  try {
    set_solid_viscosity_props[vmodel_type](simulation, viscosity_params, lDmn);
   } catch (const std::bad_function_call& exception) {
    throw std::runtime_error("[read_solid_visc_model] Viscosity model '" + vmodel_str + "' is not supported.");
  }
}

//--------------------
// read_wall_props_ff
//--------------------
// Read CMM variable wall properties from a file.
//
// Modifies:
//   com_mod.msh[iM].x - seems to use this as a scratch array.
//
void read_wall_props_ff(ComMod& com_mod, const std::string& file_name, const int iM, const int iFa)
{
  #define n_debug_read_wall_props_ff
  #ifdef debug_read_wall_props_ff
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "com_mod.cmmInit: " << com_mod.cmmInit;
  dmsg << "com_mod.varWallProps.nrows(): " << com_mod.varWallProps.nrows();
  dmsg << "com_mod.varWallProps.ncols(): " << com_mod.varWallProps.ncols();
  #endif

  if (com_mod.cmmInit) {
    auto& mesh = com_mod.msh[iM];
    mesh.x.resize(1, mesh.gnNo);

    // Read thickness
    int data_comp = 1;
    int data_series = 1;
    vtk_xml::read_vtu_pdata(file_name, "Thickness", com_mod.nsd, data_comp, data_series, mesh);

    for (int a = 0; a < mesh.gnNo; a++) {
      int Ac = mesh.gN[a];
      com_mod.varWallProps(0,Ac) = mesh.x(0,a);
    }

    // Read elasticity modulus
    mesh.x = 0.0; 
    data_comp = 1;
    data_series = 1;
    #ifdef debug_read_wall_props_ff
    dmsg << "read_vtp_pdata Elasticity_modulus ..." << com_mod.cmmInit;
    #endif
    vtk_xml::read_vtu_pdata(file_name, "Elasticity_modulus", com_mod.nsd, data_comp, data_series, mesh);

    #ifdef debug_read_wall_props_ff
    dmsg << "Set com_mod.varWallProps ... " << com_mod.cmmInit;
    #endif
    for (int a = 0; a < mesh.gnNo; a++) {
      int Ac = mesh.gN[a];
      com_mod.varWallProps(1,Ac) = mesh.x(0,a);
    }

  } else { 
    auto& mesh = com_mod.msh[iM];
    auto& face = mesh.fa[iFa];
    face.x.resize(1, face.nNo);

    // Read thickness
    int data_comp = 1;
    int data_series = 0;
    #ifdef debug_read_wall_props_ff
    dmsg << "read_vtp_pdata Thickness " << " ...";
    #endif
    vtk_xml::read_vtp_pdata(file_name, "Thickness", com_mod.nsd, data_comp, data_series, face);

    #ifdef debug_read_wall_props_ff
    dmsg << "Set com_mod.varWallProps(0,Ac) " << " ...";
    #endif

    for (int a = 0; a < face.nNo; a++) {
      int Ac = face.gN[a];
      Ac = mesh.gN[Ac];
      com_mod.varWallProps(0,Ac) = face.x(0,a);
    }

    #ifdef debug_read_wall_props_ff
    dmsg << "read_vtp_pdata Elasticity_modulus " << "...";
    #endif
    vtk_xml::read_vtp_pdata(file_name, "Elasticity_modulus", com_mod.nsd, data_comp, data_series, face);

    #ifdef debug_read_wall_props_ff
    dmsg << "Set com_mod.varWallProps(1,Ac) " << " ...";
    #endif

    for (int a = 0; a < face.nNo; a++) {
      int Ac = face.gN[a];
      Ac = mesh.gN[Ac];
      com_mod.varWallProps(1,Ac) = face.x(0,a);
    }
  }
}

//--------------
// set_cmm_bdry
//--------------
// Find boundary edge nodes for CMM initialization.
//
// Modifies:
//  lM.eAdj.nnz 
//  lM.eAdj.prow
//  lM.eAdj.pcol
//
void set_cmm_bdry(mshType& lM, Vector<int>& bNds)
{
  // Get mesh adjacency.
  //
  Vector<int> incL(lM.gnNo);

  for (int e = 0; e < lM.gnEl; e++) {
    for (int a = 0; a < lM.eNoN; a++) {
      int Ac = lM.gIEN(a,e);
      incL[Ac] = incL[Ac] + 1;
    }
  }

  int nAdj = incL.max();
  Array<int> tmpI(nAdj, lM.gnNo);
  incL = 0;
  tmpI = -1;

  for (int e = 0; e < lM.gnEl; e++) {
    for (int a = 0; a < lM.eNoN; a++) {
      int Ac = lM.gIEN(a,e);
      tmpI(incL[Ac], Ac) = e;
      incL[Ac] = incL[Ac] + 1;
    }
  }

  int b = 2*nAdj;
  bool finshed = false;
  Array<int> adjL;

  while (!finshed) {
    b = b + nAdj;
    incL.resize(lM.gnEl);
    adjL.resize(b, lM.gnEl);
    adjL = -1;
    finshed = true;

    for (int e = 0; e < lM.gnEl; e++) {
      for (int a = 0; a < lM.eNoN; a++) {
        int Ac = lM.gIEN(a,e);
        for (int i = 0; i < nAdj; i++) {
          if (tmpI(i,Ac) == -1) {
            break;
          }

          bool flag = true;
          for (int j = 0; j < incL(e); j++) {
            if (adjL(j,e) == tmpI(i,Ac)) {
              flag = false;
              break;
            }
          }

          if (flag) {
            if (incL(e) > b) {
              finshed = false;
              continue;
            }
            adjL(incL[e], e) = tmpI(i,Ac);
            incL[e] = incL[e] + 1;
          }
        }
      } 
    }
  }

  nAdj = incL.max();

  lM.eAdj.nnz = 0;

  for (int e = 0; e < lM.gnEl; e++) {
    for (int i = 0; i < nAdj; i++) {
      if (adjL(i,e) != -1) {
        lM.eAdj.nnz = lM.eAdj.nnz + 1;
      } else { 
        break;
      }
    }
  }

  lM.eAdj.prow.resize(lM.gnEl+1); 
  lM.eAdj.pcol.resize(lM.eAdj.nnz);

  int j = 0;
  lM.eAdj.prow[0] = j;

  for (int e = 0; e < lM.gnEl; e++) {
    for (int i = 0; i < nAdj; i++) {
      if (adjL(i,e) != -1) {
        lM.eAdj.pcol(j) = adjL(i,e);
        j = j + 1;
      } else { 
        break;
      }
    }
    lM.eAdj.prow[e+1] = j;
  }

  // Loop over all elements to find edge nodes
  //
  for (int e = 0; e < lM.gnEl; e++) {
    // Create a list of neighboring elements.
    int nAdj = lM.eAdj.prow[e+1] - lM.eAdj.prow[e];
    incL.resize(nAdj);
    int i = 0;

    for (int a = lM.eAdj.prow(e); a < lM.eAdj.prow[e+1]; a++) {
      incL[i] = lM.eAdj.pcol[a];
      i = i + 1;
    }

    // Select an edge pair
    for (int a = 0; a < lM.eNoN; a++) {
      int b = a + 1;
      if (a == lM.eNoN-1) {
        b = 0;
      }
      int Ac = lM.gIEN(a,e);
      int Bc = lM.gIEN(b,e);
      if (Ac > Bc) {
        utils::swap(Ac, Bc);
      }

      // Loop over all neighbors and check if this edge is found
      bool flag = true;
      for (int i = 0; i < nAdj; i++) {
        int e1 = incL[i];
        if (e1 == e) {
          continue; 
        }
        for (int a1 = 0; a1 < lM.eNoN; a1++) {
          int b1 = a1 + 1;
          if (a1 == lM.eNoN-1) {
            b1 = 0;
          }
          int Ac1 = lM.gIEN(a1,e1);
          int Bc1 = lM.gIEN(b1,e1);
          if (Ac1 > Bc1) {
            utils::swap(Ac1, Bc1);
          }
          if (Ac == Ac1 && Bc == Bc1) {
            flag = false;
            break;
          }
        }
        if (!flag) {
          break;
        }
      }

      if (flag) {
        Ac = lM.gN[Ac];
        Bc = lM.gN[Bc];
        bNds[Ac] = 1;
        bNds[Bc] = 1;
      }
    }
  }
}

//-------------------------
// set_equation_properties
//-------------------------
// Set the parameter values for the equation.
//
// This replaces the 'SELECT CASE (eqName)' statement in the Fortran 'READEQ()' subroutine.
//
void set_equation_properties(Simulation* simulation, EquationParameters* eq_params, eqType& lEq, EquationProps& propL, EquationOutputs& outPuts, 
  EquationNdop& nDOP)
{
  using namespace consts;

  EquationType eq_type;
  auto eq_type_str = eq_params->type();

  try {
    eq_type = equation_name_to_type.at(eq_type_str);
  } catch (const std::out_of_range& exception) {
    throw std::runtime_error("Unknown equation type '" + eq_type_str + ".");
  }

  try {
    set_equation_props[eq_type](simulation, eq_params, lEq, propL, outPuts, nDOP);
   } catch (const std::bad_function_call& exception) {
    throw std::runtime_error("[set_equation_properties] Equation type " + eq_type_str + " is not supported..");
  }
}

};

