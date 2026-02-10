// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "svZeroD_subroutines.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include "cmm.h"
#include "consts.h"
#include "svZeroD_interface/LPNSolverInterface.h"
#include "ComMod.h"
#include "utils.h"

#include <map>
#include <algorithm>
#include <iterator>

static std::map<int,LPNSolverInterface*> interfaces;

namespace svZeroD {

static int numCoupledSrfs;
static bool writeSvZeroD = true;
static double svZeroDTime = 0.0;

int num_output_steps;
int system_size;
int model_id;

std::vector<int> nsrflistCoupled(numCoupledSrfs);
std::vector<std::string> svzd_blk_names(numCoupledSrfs);
std::vector<int> svzd_blk_name_len(numCoupledSrfs);
std::vector<double> in_out_sign(numCoupledSrfs);
std::vector<double> lpn_times(num_output_steps);
std::vector<double> lpn_solutions((num_output_steps * system_size));
std::vector<double> lpn_state_y(system_size);
std::vector<double> last_state_y(system_size);
std::vector<double> last_state_ydot(system_size);
std::vector<int> sol_IDs(2 * numCoupledSrfs);

void create_svZeroD_model(std::string lpn_library_name, std::string lpn_json_file)
{
    // Load library
  auto interface = new LPNSolverInterface();
  
  // Get correct library name based on operating system
  std::string interface_lib_path = lpn_library_name.substr(0, lpn_library_name.find("libsvzero_interface"));
  std::string interface_lib_so = interface_lib_path + "libsvzero_interface.so";
  std::string interface_lib_dylib = interface_lib_path + "libsvzero_interface.dylib";
  std::ifstream lib_so_exists(interface_lib_so);
  std::ifstream lib_dylib_exists(interface_lib_dylib);
  if (lib_so_exists) {
    interface->load_library(interface_lib_so);
  } else if (lib_dylib_exists) {
    interface->load_library(interface_lib_dylib);
  } else {
    throw std::runtime_error("Could not find shared libraries " + interface_lib_so + " or " + interface_lib_dylib);
  }

  // Initialize model
  interface->initialize(std::string(lpn_json_file));
  model_id = interface->problem_id_;
  interfaces[model_id] = interface;
  
  // Save model parameters
  num_output_steps = interface->num_output_steps_;
  system_size = interface->system_size_;
}

void get_svZeroD_variable_ids(std::string block_name, int* blk_ids, double* inlet_or_outlet)
{
  auto interface = interfaces[model_id];
  std::vector<int> IDs;
  interface->get_block_node_IDs(block_name, IDs);
  // IDs in the above function stores info in the following format:
  // {num inlet nodes, inlet flow[0], inlet pressure[0],..., num outlet nodes, outlet flow[0], outlet pressure[0],...}
  int num_inlet_nodes = IDs[0];
  int num_outlet_nodes = IDs[1+num_inlet_nodes*2];
  if ((num_inlet_nodes == 0) && (num_outlet_nodes == 1)) {
    blk_ids[0] = IDs[1+num_inlet_nodes*2+1]; // Outlet flow
    blk_ids[1] = IDs[1+num_inlet_nodes*2+2]; // Outlet pressure
    *inlet_or_outlet = -1.0; // Signifies inlet to LPN
  } else if ((num_inlet_nodes == 1) && (num_outlet_nodes == 0)) {
    blk_ids[0] = IDs[1]; // Inlet flow
    blk_ids[1] = IDs[2]; // Inlet pressure
    *inlet_or_outlet = 1.0; // Signifies outlet to LPN
  } else {
    std::runtime_error("ERROR: [lpn_interface_get_variable_ids] Not a flow/pressure block.");
  }
}


void update_svZeroD_block_params(std::string block_name, double* time, double* params)
{
  auto interface = interfaces[model_id];
  int param_len = 2; // Usually 2 for this use case
  std::vector<double> new_params(1+2*param_len);
  // Format of new_params for flow/pressure blocks: 
  // [N, time_1, time_2, ..., time_N, value1, value2, ..., value_N]
  // where N is number of time points and value* is flow/pressure
  new_params[0] = (double) param_len;
  for (int i = 0; i < param_len; i++) {
    new_params[1+i] = time[i];
    new_params[1+param_len+i] = params[i];
  }
  interface->update_block_params(block_name, new_params);
}


void write_svZeroD_solution(const double* lpn_time, std::vector<double>& lpn_solution, int* flag)
{
  auto interface = interfaces[model_id];
  if (*flag == 0) { // Initialize output file: Write header with variable names
    std::vector<std::string> variable_names;
    variable_names = interface->variable_names_;
    std::ofstream out_file;
    out_file.open("svZeroD_data", std::ios::out | std::ios::app);
    out_file<<system_size<<" ";
    for (int i = 0; i < system_size; i++) {
      out_file<<static_cast<std::string>(variable_names[i])<<" ";
    }
    out_file<<'\n';
  } else {
    std::ofstream out_file;
    out_file.open("svZeroD_data", std::ios::out | std::ios::app);
    out_file<<*lpn_time<<" ";
    for (int i = 0; i < system_size; i++) {
      out_file<<lpn_solution[i]<<" ";
    }
    out_file<<'\n';
    out_file.close();
  }
}

void get_coupled_QP(ComMod& com_mod, const CmMod& cm_mod, double QCoupled[], double QnCoupled[], double PCoupled[], double PnCoupled[]){
  using namespace consts;
  
  auto& cplBC = com_mod.cplBC;
  int ind = 0;
  
  // Get Q/P from cplBC.fa for Dir BCs
  for (int iFa = 0; iFa < cplBC.nFa; iFa++) {
    auto& fa = cplBC.fa[iFa];
    if (fa.bGrp == consts::CplBCType::cplBC_Dir) {
      QCoupled[ind] = fa.Qo;
      QnCoupled[ind] = fa.Qn;
      PCoupled[ind] = fa.Po;
      PnCoupled[ind] = fa.Pn;
      ind = ind + 1;
    }
  }
  
  // Get Q/P from cplBC.fa for standard Neu BCs
  for (int iFa = 0; iFa < cplBC.nFa; iFa++) {
    auto& fa = cplBC.fa[iFa];
    if (fa.bGrp == consts::CplBCType::cplBC_Neu) {
      QCoupled[ind] = fa.Qo;
      QnCoupled[ind] = fa.Qn;
      PCoupled[ind] = fa.Po;
      PnCoupled[ind] = fa.Pn;
      ind = ind + 1;
    }
  }
  
  // Get Q/P from ZeroDBoundaryCondition for ZeroD BCs by iterating through eq.bc[]
  int iBC_ZeroD = static_cast<int>(BoundaryConditionType::bType_ZeroD);
  for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
    auto& eq = com_mod.eq[iEq];
    for (int iBc = 0; iBc < eq.nBc; iBc++) {
      auto& bc = eq.bc[iBc];
      if (utils::btest(bc.bType, iBC_ZeroD)) {
        // Only process Neumann ZeroD BCs here
        if (bc.zerod_bc.get_bc_type() == consts::BoundaryConditionType::bType_Neu) {
          QCoupled[ind] = bc.zerod_bc.get_Qo();
          QnCoupled[ind] = bc.zerod_bc.get_Qn();
          PCoupled[ind] = 0.0;  // Pressure not used for Neumann BC input
          PnCoupled[ind] = 0.0;
          ind = ind + 1;
        } else if (bc.zerod_bc.get_bc_type() == consts::BoundaryConditionType::bType_Dir) {
          QCoupled[ind] = bc.zerod_bc.get_Qo();
          QnCoupled[ind] = bc.zerod_bc.get_Qn();
          PCoupled[ind] = bc.zerod_bc.get_Po();
          PnCoupled[ind] = bc.zerod_bc.get_Pn();
          ind = ind + 1;
        } else {
          throw std::runtime_error("ERROR: [get_coupled_QP] Invalid ZeroD BC type.");
        }
      }
    }
  }
}

void print_svZeroD(int* nSrfs, std::vector<int>& surfID, double Q[], double P[]) {
  int nParam = 2;
  const char* fileNames[2] = {"Q_svZeroD", "P_svZeroD"};
  std::vector<std::vector<double>> R(nParam, std::vector<double>(*nSrfs));

  if (*nSrfs == 0) return;

  for (int i = 0; i < *nSrfs; ++i) {
    R[0][i] = Q[i];
    R[1][i] = P[i];
  }

  // Set formats
  std::string myFMT1 = "(" + std::to_string(*nSrfs) + "(E13.5))";
  std::string myFMT2 = "(" + std::to_string(*nSrfs) + "(I13))";

  for (int i = 0; i < nParam; ++i) {
    std::ifstream file(fileNames[i]);
    if (!file) {
      std::ofstream newFile(fileNames[i], std::ios::app);
      for (int j = 0; j < *nSrfs; ++j) {
        newFile << std::scientific << std::setprecision(5) << R[i][j] << std::endl;
      }
    } else {
      std::ofstream newFile(fileNames[i]);
      for (int j = 0; j < *nSrfs; ++j) {
        newFile << std::setw(13) << surfID[j] << std::endl;
      }
      for (int j = 0; j < *nSrfs; ++j) {
        newFile << std::scientific << std::setprecision(5) << R[i][j] << std::endl;
      }
    }
  }
}

//--------------
// init_svZeroD
//--------------
//
void init_svZeroD(ComMod& com_mod, const CmMod& cm_mod) 
{
  using namespace consts;
  
  #define n_debug_init_svZeroD
  #ifdef debug_init_svZeroD
  DebugMsg dmsg(__func__, com_mod.cm.idcm()); 
  dmsg.banner();
  #endif

  auto& cplBC = com_mod.cplBC;
  auto& solver_interface = cplBC.svzerod_solver_interface;
  auto& cm = com_mod.cm;
  double dt = com_mod.dt;

  // Count ZeroD BCs by iterating through eq.bc[]
  int iBC_ZeroD = static_cast<int>(BoundaryConditionType::bType_ZeroD);
  int nZeroD_count = 0;
  for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
    auto& eq = com_mod.eq[iEq];
    for (int iBc = 0; iBc < eq.nBc; iBc++) {
      auto& bc = eq.bc[iBc];
      if (utils::btest(bc.bType, iBC_ZeroD)) {
        nZeroD_count++;
      }
    }
  }
  
  // Total coupled surfaces = cplBC.nFa (Dir + Neu) + ZeroD BCs
  numCoupledSrfs = cplBC.nFa + nZeroD_count;
  
  int nDir = 0;
  int nNeu = 0;
  int nZeroD = 0;
  
  // If this process is the master process on the communicator
  if (cm.mas(cm_mod)) {
    // Count Dir BCs from cplBC
    for (int iFa = 0; iFa < cplBC.nFa; iFa++) {
      auto& fa = cplBC.fa[iFa];

      if (fa.bGrp == consts::CplBCType::cplBC_Dir) {
        nsrflistCoupled.push_back(iFa);
        nDir = nDir + 1;
      }
    }
    
    // Count Neu BCs from cplBC
    for (int iFa = 0; iFa < cplBC.nFa; iFa++) {
      auto& fa = cplBC.fa[iFa];

      if (fa.bGrp == consts::CplBCType::cplBC_Neu) {
        nsrflistCoupled.push_back(iFa);
        nNeu = nNeu + 1;
      }
    }
    
    // Count ZeroD BCs by iterating through eq.bc[]
    // Use a special offset for ZeroD surface IDs to distinguish from cplBC indices
    for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
      auto& eq = com_mod.eq[iEq];
      for (int iBc = 0; iBc < eq.nBc; iBc++) {
        auto& bc = eq.bc[iBc];
        if (utils::btest(bc.bType, iBC_ZeroD)) {
          nsrflistCoupled.push_back(cplBC.nFa + nZeroD);
          nZeroD = nZeroD + 1;
        }
      }
    }
  }

  #ifdef debug_init_svZeroD
  dmsg << "nDir: " << nDir;
  dmsg << "nNeu: " << nNeu;
  dmsg << "numCoupledSrfs: " << numCoupledSrfs;
  for (int surf_id : nsrflistCoupled) {
    dmsg << "surface id: " << surf_id;
  }
  #endif

  std::string svzerod_library;
  std::string svzerod_file;
  std::string buffer;
  int ids[2];
  
  // For Dir and Neu BCs (from cplBC)
  int nCplBCFaces = cplBC.nFa;
  std::vector<std::string> svzd_blk_names_unsrtd(nCplBCFaces);
  std::vector<int> svzd_blk_ids(nCplBCFaces);
  int init_flow_flag, init_press_flag;
  double init_flow, init_press, in_out;
  
  if (cm.mas(cm_mod)) {

    if (solver_interface.has_data) { 
      #ifdef debug_init_svZeroD
      dmsg << "#### Use XML data #### " << " ";
      #endif

      svzerod_library = solver_interface.solver_library;
      svzerod_file = solver_interface.configuration_file;

      // Process block_surface_map for Dir and Neu BCs (not ZeroD)
      int i = 0;
      for (const auto& pair : solver_interface.block_surface_map) {
        #ifdef debug_init_svZeroD
        dmsg << "block_surface_map: '" + pair.first << "'";
        #endif
        
        // Skip if this block is handled by a ZeroD BC (check by iterating through eq.bc[])
        bool is_zerod = false;
        for (int iEq = 0; iEq < com_mod.nEq && !is_zerod; iEq++) {
          auto& eq = com_mod.eq[iEq];
          for (int iBc = 0; iBc < eq.nBc && !is_zerod; iBc++) {
            auto& bc = eq.bc[iBc];
            if (utils::btest(bc.bType, iBC_ZeroD) && bc.zerod_bc.get_block_name() == pair.first) {
              is_zerod = true;
            }
          }
        }
        
        if (is_zerod) {
          continue;  // Skip - this block is handled by ZeroDBoundaryCondition
        }
        
        if (i >= nCplBCFaces) {
          break;  // No more space in cplBC arrays
        }
        
        svzd_blk_ids[i] = -1;
        svzd_blk_names_unsrtd[i] = pair.first;
        for (int j = 0; j < cplBC.nFa; j++) {
          auto& fa = cplBC.fa[j];
          if (fa.name == pair.second) { 
            svzd_blk_ids[i] = j;
          }
        }

        if (svzd_blk_ids[i] == -1) { 
          throw std::runtime_error("ERROR: Did not find a coupled boundary condition for block '" + 
              pair.first + "' and surface '" + pair.second + "'; check the Block_to_surface_map solver XML parameter.");
        }
 
        i += 1;
      }

      init_flow_flag = solver_interface.have_initial_flows;
      init_flow = solver_interface.initial_flows;

      init_press_flag = solver_interface.have_initial_pressures;
      init_press = solver_interface.initial_pressures;

    } 

    // Arrange svzd_blk_names in the same order as surface IDs in nsrflistCoupled
    #ifdef debug_init_svZeroD
    dmsg << "Arrange svzd_blk_names ... " << " ";
    #endif

    for (int s = 0; s < numCoupledSrfs; ++s) {
      int found = 0;
      int surf_id = nsrflistCoupled[s];
      
      #ifdef debug_init_svZeroD
      dmsg << ">>> s " << s;
      dmsg << "  nsrflistCoupled[s]: " << surf_id;
      #endif

      // Check if this is a ZeroD BC (surface ID >= cplBC.nFa)
      if (surf_id >= cplBC.nFa) {
        // This is a ZeroD BC - find it by iterating through eq.bc[]
        int zerod_idx = surf_id - cplBC.nFa;
        int current_idx = 0;
        for (int iEq = 0; iEq < com_mod.nEq && !found; iEq++) {
          auto& eq = com_mod.eq[iEq];
          for (int iBc = 0; iBc < eq.nBc && !found; iBc++) {
            auto& bc = eq.bc[iBc];
            if (utils::btest(bc.bType, iBC_ZeroD)) {
              if (current_idx == zerod_idx) {
                std::string blk_name = bc.zerod_bc.get_block_name();
                svzd_blk_names.push_back(blk_name);
                svzd_blk_name_len.push_back(blk_name.length());
                found = 1;
                #ifdef debug_init_svZeroD
                dmsg << "    Found ZeroD block: '" << blk_name << "'";
                #endif
              }
              current_idx++;
            }
          }
        }
      } else {
        // This is a Dir or Neu BC - search in svzd_blk_ids
        for (int t = 0; t < nCplBCFaces; ++t) {
          #ifdef debug_init_svZeroD
          dmsg << "  >>> t " << t;
          dmsg << "    svzd_blk_ids[t]: " << svzd_blk_ids[t];
          #endif
          if (svzd_blk_ids[t] == surf_id) {
            #ifdef debug_init_svZeroD
            dmsg << "    Found  " << " " ;
            dmsg << "    svzd_blk_names_unsrtd[t]: '" << svzd_blk_names_unsrtd[t];
            #endif
            found = 1;
            svzd_blk_names.push_back(svzd_blk_names_unsrtd[t]);
            svzd_blk_name_len.push_back(svzd_blk_names_unsrtd[t].length());
            break;
          }
        }
      }

      if (found == 0) {
        throw std::runtime_error("ERROR: Did not find block name for surface ID: " + std::to_string(surf_id));
      }
    }

    // Create the svZeroD model
    create_svZeroD_model(svzerod_library, svzerod_file);
    auto interface = interfaces[model_id];
    interface->set_external_step_size(dt);

    // Save IDs of relevant variables in the solution vector
    sol_IDs.assign(2 * numCoupledSrfs, 0);
    for (int s = 0; s < numCoupledSrfs; ++s) {
      int len = svzd_blk_name_len[s];
      get_svZeroD_variable_ids(svzd_blk_names[s], ids, &in_out);
      sol_IDs[2 * s] = ids[0];
      sol_IDs[2 * s + 1] = ids[1];
      in_out_sign.push_back(in_out);
      
      // For ZeroD BCs, store the solution IDs in ZeroDBoundaryCondition
      int surf_id = nsrflistCoupled[s];
      if (surf_id >= cplBC.nFa) {
        int zerod_idx = surf_id - cplBC.nFa;
        int current_idx = 0;
        for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
          auto& eq = com_mod.eq[iEq];
          for (int iBc = 0; iBc < eq.nBc; iBc++) {
            auto& bc = eq.bc[iBc];
            if (utils::btest(bc.bType, iBC_ZeroD)) {
              if (current_idx == zerod_idx) {
                bc.zerod_bc.set_solution_ids(ids[0], ids[1], in_out);
              }
              current_idx++;
            }
          }
        }
      }
    }

    // Initialize lpn_state variables corresponding to external coupling blocks
    lpn_times.assign(num_output_steps, 0.0);
    lpn_solutions.assign(num_output_steps*system_size, 0.0);
    lpn_state_y.assign(system_size, 0.0);
    last_state_y.assign(system_size, 0.0);
    last_state_ydot.assign(system_size, 0.0);

    interface->return_y(lpn_state_y);
    interface->return_ydot(last_state_ydot);
    
    for (int s = 0; s < numCoupledSrfs; ++s) {
      int surf_id = nsrflistCoupled[s];
      
      if (surf_id >= cplBC.nFa) {
        // ZeroD BC - initialize by iterating through eq.bc[]
        int zerod_idx = surf_id - cplBC.nFa;
        int current_idx = 0;
        for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
          auto& eq = com_mod.eq[iEq];
          for (int iBc = 0; iBc < eq.nBc; iBc++) {
            auto& bc = eq.bc[iBc];
            if (utils::btest(bc.bType, iBC_ZeroD)) {
              if (current_idx == zerod_idx) {
                if (bc.zerod_bc.get_bc_type() == consts::BoundaryConditionType::bType_Neu) {
                  // Neumann ZeroD BC: initialize flow, pressure comes from 0D solver
                  if (init_flow_flag == 1) {
                    lpn_state_y[sol_IDs[2 * s]] = init_flow;
                  }
                  if (init_press_flag == 1) {
                    lpn_state_y[sol_IDs[2 * s + 1]] = init_press;
                    bc.zerod_bc.set_pressure(lpn_state_y[sol_IDs[2 * s + 1]]);
                  }
                } else if (bc.zerod_bc.get_bc_type() == consts::BoundaryConditionType::bType_Dir) {
                  // Dirichlet ZeroD BC: initialize pressure, flow comes from 0D solver
                  if (init_press_flag == 1) {
                    lpn_state_y[sol_IDs[2 * s + 1]] = init_press;
                    bc.zerod_bc.set_pressure(lpn_state_y[sol_IDs[2 * s + 1]]);
                  }
                  if (init_flow_flag == 1) {
                    lpn_state_y[sol_IDs[2 * s]] = init_flow;
                  }
                } else {
                  throw std::runtime_error("ERROR: [init_svZeroD] Invalid ZeroD BC type.");
                }
              }
              current_idx++;
            }
          }
        }
      } else {
        // Dir or Neu BC - use cplBC.fa
        if (init_flow_flag == 1) {
          lpn_state_y[sol_IDs[2 * s]] = init_flow;
          cplBC.fa[surf_id].y = lpn_state_y[sol_IDs[2 * s]];
        }
        if (init_press_flag == 1) {
          lpn_state_y[sol_IDs[2 * s + 1]] = init_press;
          cplBC.fa[surf_id].y = lpn_state_y[sol_IDs[2 * s + 1]];
        }
      }
    }
    std::copy(lpn_state_y.begin(), lpn_state_y.end(), last_state_y.begin());

    if (writeSvZeroD == 1) {
      // Initialize output file
      int flag = 0;
      write_svZeroD_solution(&svZeroDTime, lpn_state_y, &flag);
    }
  }

  // Broadcast initial values to follower processes
  if (!cm.seq()) {
    // For cplBC.fa (Dir and Neu BCs) - only if there are any
    if (cplBC.nFa > 0) {
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
    
    // For ZeroD BCs - broadcast pressure values (only for Neumann, Dirichlet has pressure as input)
    if (nZeroD_count > 0) {
      // Count Neumann ZeroD BCs for broadcasting
      int nZeroD_Neu = 0;
      for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
        auto& eq = com_mod.eq[iEq];
        for (int iBc = 0; iBc < eq.nBc; iBc++) {
          auto& bc = eq.bc[iBc];
          if (utils::btest(bc.bType, iBC_ZeroD) && 
              bc.zerod_bc.get_bc_type() == consts::BoundaryConditionType::bType_Neu) {
            nZeroD_Neu++;
          }
        }
      }
      
      if (nZeroD_Neu > 0) {
        Vector<double> p(nZeroD_Neu);
        
        if (cm.mas(cm_mod)) {
          int idx = 0;
          for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
            auto& eq = com_mod.eq[iEq];
            for (int iBc = 0; iBc < eq.nBc; iBc++) {
              auto& bc = eq.bc[iBc];
              if (utils::btest(bc.bType, iBC_ZeroD) && 
                  bc.zerod_bc.get_bc_type() == consts::BoundaryConditionType::bType_Neu) {
                p(idx) = bc.zerod_bc.get_pressure();
                idx++;
              }
            }
          }
        }
        
        cm.bcast(cm_mod, p);
        
        if (cm.slv(cm_mod)) {
          int idx = 0;
          for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
            auto& eq = com_mod.eq[iEq];
            for (int iBc = 0; iBc < eq.nBc; iBc++) {
              auto& bc = eq.bc[iBc];
              if (utils::btest(bc.bType, iBC_ZeroD) && 
                  bc.zerod_bc.get_bc_type() == consts::BoundaryConditionType::bType_Neu) {
                bc.zerod_bc.set_pressure(p(idx));
                idx++;
              }
            }
          }
        }
      }
    }
  }
}


void calc_svZeroD(ComMod& com_mod, const CmMod& cm_mod, char BCFlag) {
  using namespace consts;
  
  int nDir = 0;
  int nNeu = 0;
  double dt = com_mod.dt;
  auto& cplBC = com_mod.cplBC;
  auto& cm = com_mod.cm;
  
  // Count ZeroD BCs by iterating through eq.bc[]
  int iBC_ZeroD = static_cast<int>(BoundaryConditionType::bType_ZeroD);
  int nZeroD = 0;
  for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
    auto& eq = com_mod.eq[iEq];
    for (int iBc = 0; iBc < eq.nBc; iBc++) {
      auto& bc = eq.bc[iBc];
      if (utils::btest(bc.bType, iBC_ZeroD)) {
        nZeroD++;
      }
    }
  }

  // If this process is the master process on the communicator
  if (cm.mas(cm_mod)) {
    // Count Dir and Neu BCs from cplBC
    for (int iFa = 0; iFa < cplBC.nFa; iFa++) {
      auto& fa = cplBC.fa[iFa];

      if (fa.bGrp == consts::CplBCType::cplBC_Dir) {
        nDir = nDir + 1;
      } else if (fa.bGrp == consts::CplBCType::cplBC_Neu) {
        nNeu = nNeu + 1;
      }
    }

    double QCoupled[numCoupledSrfs], QnCoupled[numCoupledSrfs], PCoupled[numCoupledSrfs], PnCoupled[numCoupledSrfs];
    double total_flow;
    double params[2];
    double times[2];
    int error_code;
    
    get_coupled_QP(com_mod, cm_mod, QCoupled, QnCoupled, PCoupled, PnCoupled);

    if (writeSvZeroD == 1) {
      if (BCFlag == 'L') {
        int i = numCoupledSrfs;
        print_svZeroD(&i, nsrflistCoupled, QCoupled, PCoupled);
      }
    }

    auto interface = interfaces[model_id];
    
    if (BCFlag != 'I') {
      // Set initial condition from the previous state
      interface->update_state(last_state_y, last_state_ydot);

      times[0] = svZeroDTime;
      times[1] = svZeroDTime + com_mod.dt;

      total_flow = 0.0;

      // Update pressure and flow in the zeroD model
      for (int i = 0; i < numCoupledSrfs; ++i) {
        int surf_id = nsrflistCoupled[i];
        bool is_dirichlet = false;
        
        if (i < nDir) {
          // Standard Dirichlet BC from cplBC
          is_dirichlet = true;
        } else if (surf_id >= cplBC.nFa) {
          // ZeroD BC - check actual type
          int zerod_idx = surf_id - cplBC.nFa;
          int current_idx = 0;
          for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
            auto& eq = com_mod.eq[iEq];
            for (int iBc = 0; iBc < eq.nBc; iBc++) {
              auto& bc = eq.bc[iBc];
              if (utils::btest(bc.bType, iBC_ZeroD)) {
                if (current_idx == zerod_idx) {
                  is_dirichlet = (bc.zerod_bc.get_bc_type() == consts::BoundaryConditionType::bType_Dir);
                  break;
                }
                current_idx++;
              }
            }
            if (is_dirichlet || current_idx > zerod_idx) break;
          }
        }
        // Otherwise it's a standard Neumann BC from cplBC
        
        if (is_dirichlet) {
          // Dirichlet BC: use pressure as input
          params[0] = PCoupled[i];
          params[1] = PnCoupled[i];
        } else {
          // Neumann BC: use flow as input
          // For ZeroD BCs, use the sign from ZeroDBoundaryCondition
          double sign = in_out_sign[i];
          if (surf_id >= cplBC.nFa) {
            // ZeroD BC - get sign from ZeroDBoundaryCondition
            int zerod_idx = surf_id - cplBC.nFa;
            int current_idx = 0;
            for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
              auto& eq = com_mod.eq[iEq];
              for (int iBc = 0; iBc < eq.nBc; iBc++) {
                auto& bc = eq.bc[iBc];
                if (utils::btest(bc.bType, iBC_ZeroD)) {
                  if (current_idx == zerod_idx) {
                    sign = bc.zerod_bc.get_in_out_sign();
                    break;
                  }
                  current_idx++;
                }
              }
            }
          }
          params[0] = sign * QCoupled[i];
          params[1] = sign * QnCoupled[i];
          total_flow += QCoupled[i];
        }
        update_svZeroD_block_params(svzd_blk_names[i], times, params);
      }

      // Run zeroD simulation
      interface->run_simulation(svZeroDTime, lpn_times, lpn_solutions, error_code);

      // Extract pressure and flow from zeroD solution
      std::copy(lpn_solutions.begin() + (num_output_steps-1)*system_size, lpn_solutions.end(), lpn_state_y.begin());
      
      for (int i = 0; i < numCoupledSrfs; ++i) {
        int surf_id = nsrflistCoupled[i];
        
        if (i < nDir) {
          // Dirichlet BC - get flow from 0D
          QCoupled[i] = in_out_sign[i] * lpn_state_y[sol_IDs[2 * i]];
          cplBC.fa[surf_id].y = QCoupled[i];
        } else if (surf_id >= cplBC.nFa) {
          // ZeroD BC - handle based on type
          int zerod_idx = surf_id - cplBC.nFa;
          int current_idx = 0;
          for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
            auto& eq = com_mod.eq[iEq];
            for (int iBc = 0; iBc < eq.nBc; iBc++) {
              auto& bc = eq.bc[iBc];
              if (utils::btest(bc.bType, iBC_ZeroD)) {
                if (current_idx == zerod_idx) {
                  // Use solution IDs stored in ZeroDBoundaryCondition
                  int flow_id = bc.zerod_bc.get_flow_sol_id();
                  int pressure_id = bc.zerod_bc.get_pressure_sol_id();
                  double in_out = bc.zerod_bc.get_in_out_sign();
                  
                  // Validate solution IDs
                  if (flow_id < 0 || flow_id >= system_size || pressure_id < 0 || pressure_id >= system_size) {
                    throw std::runtime_error("ERROR: [calc_svZeroD] Invalid solution IDs for ZeroD BC: flow_id=" + 
                                            std::to_string(flow_id) + ", pressure_id=" + std::to_string(pressure_id) + 
                                            ", system_size=" + std::to_string(system_size));
                  }
                  
                  if (bc.zerod_bc.get_bc_type() == consts::BoundaryConditionType::bType_Neu) {
                    // Neumann ZeroD BC: get pressure from 0D solver
                    PCoupled[i] = lpn_state_y[pressure_id];
                    bc.zerod_bc.set_pressure(PCoupled[i]);
                  } else if (bc.zerod_bc.get_bc_type() == consts::BoundaryConditionType::bType_Dir) {
                    // Dirichlet ZeroD BC: get flow from 0D solver
                    QCoupled[i] = in_out * lpn_state_y[flow_id];
                    // Store flow in ZeroDBoundaryCondition (use previous Qn as Qo, current as Qn)
                    double Qo_prev = bc.zerod_bc.get_Qn();
                    bc.zerod_bc.set_flowrates(Qo_prev, QCoupled[i]);
                  } else {
                    throw std::runtime_error("ERROR: [calc_svZeroD] Invalid ZeroD BC type.");
                  }
                  break;
                }
                current_idx++;
              }
            }
          }
        } else {
          // Standard Neu BC - use cplBC.fa
          PCoupled[i] = lpn_state_y[sol_IDs[2 * i + 1]];
          cplBC.fa[surf_id].y = PCoupled[i];
        }
      }

      if (BCFlag == 'L') {
        // Save state and update time only after the last inner iteration
        interface->return_ydot(last_state_ydot);
        std::copy(lpn_state_y.begin(), lpn_state_y.end(), last_state_y.begin());

        if (writeSvZeroD == 1) {
          // Write the state vector to a file
          int arg = 1;
          write_svZeroD_solution(&svZeroDTime, lpn_state_y, &arg);
        }

        // Keep track of the current time
        svZeroDTime = svZeroDTime + com_mod.dt;
      }
    }
  }

  // If there are multiple procs (not sequential), broadcast outputs to follower procs
  if (!cm.seq()) {
    // Broadcast cplBC.fa values (Dir and standard Neu BCs) - only if there are any
    if (cplBC.nFa > 0) {
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
    
    // Broadcast ZeroD BC pressure values (only for Neumann, Dirichlet has pressure as input)
    if (nZeroD > 0) {
      // Count Neumann ZeroD BCs for broadcasting
      int nZeroD_Neu = 0;
      for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
        auto& eq = com_mod.eq[iEq];
        for (int iBc = 0; iBc < eq.nBc; iBc++) {
          auto& bc = eq.bc[iBc];
          if (utils::btest(bc.bType, iBC_ZeroD) && 
              bc.zerod_bc.get_bc_type() == consts::BoundaryConditionType::bType_Neu) {
            nZeroD_Neu++;
          }
        }
      }
      
      if (nZeroD_Neu > 0) {
        Vector<double> p(nZeroD_Neu);
        
        if (cm.mas(cm_mod)) {
          int idx = 0;
          for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
            auto& eq = com_mod.eq[iEq];
            for (int iBc = 0; iBc < eq.nBc; iBc++) {
              auto& bc = eq.bc[iBc];
              if (utils::btest(bc.bType, iBC_ZeroD) && 
                  bc.zerod_bc.get_bc_type() == consts::BoundaryConditionType::bType_Neu) {
                p(idx) = bc.zerod_bc.get_pressure();
                idx++;
              }
            }
          }
        }
        
        cm.bcast(cm_mod, p);
        
        if (cm.slv(cm_mod)) {
          int idx = 0;
          for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
            auto& eq = com_mod.eq[iEq];
            for (int iBc = 0; iBc < eq.nBc; iBc++) {
              auto& bc = eq.bc[iBc];
              if (utils::btest(bc.bType, iBC_ZeroD) && 
                  bc.zerod_bc.get_bc_type() == consts::BoundaryConditionType::bType_Neu) {
                bc.zerod_bc.set_pressure(p(idx));
                idx++;
              }
            }
          }
        }
      }
    }
  }
}
};
