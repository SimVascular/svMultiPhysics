// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "svZeroD_interface.h"

#include <fstream>
#include <string>
#include <iomanip>
#include "consts.h"
#include "svZeroD_interface/LPNSolverInterface.h"
#include "ComMod.h"
#include "utils.h"

#include <map>
#include <algorithm>

static std::map<int,LPNSolverInterface*> interfaces;

namespace svZeroD {

namespace {

/// Fills \c com_mod.cplBC.nSvZeroD_coupled_bc and \c svZeroD_coupled_bc_idxs (same order as nested eq/bc loops).
static void build_svzero_coupled_bc_idxs(ComMod& com_mod)
{
  auto& cpl = com_mod.cplBC;
  cpl.svZeroD_coupled_bc_idxs.clear();
  for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
    for (int iBc = 0; iBc < com_mod.eq[iEq].nBc; iBc++) {
      if (utils::btest(com_mod.eq[iEq].bc[iBc].bType, consts::iBC_Coupled)) {
        cpl.svZeroD_coupled_bc_idxs.emplace_back(iEq, iBc);
      }
    }
  }
  cpl.nSvZeroD_coupled_bc = static_cast<int>(cpl.svZeroD_coupled_bc_idxs.size());
}

static bool nth_coupled_bc(ComMod& com_mod, int n, bcType** out_bc)
{
  auto& cpl = com_mod.cplBC;
  if (n < 0 || n >= cpl.nSvZeroD_coupled_bc) {
    return false;
  }
  const auto pr = cpl.svZeroD_coupled_bc_idxs[static_cast<size_t>(n)];
  *out_bc = &com_mod.eq[pr.first].bc[pr.second];
  return true;
}

}  // namespace

static int numCoupledSrfs;
static bool writeSvZeroD = true;
static double svZeroDTime = 0.0;

int num_output_steps;
int system_size;
int model_id;

std::vector<int> nsrflistCoupled(numCoupledSrfs);
std::vector<std::string> svzd_blk_names(numCoupledSrfs);
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
    *inlet_or_outlet = 1.0; // Signifies inlet to LPN
  } else if ((num_inlet_nodes == 1) && (num_outlet_nodes == 0)) {
    blk_ids[0] = IDs[1]; // Inlet flow
    blk_ids[1] = IDs[2]; // Inlet pressure
    *inlet_or_outlet = -1.0; // Signifies outlet to LPN
  } else {
    throw std::runtime_error("ERROR: [lpn_interface_get_variable_ids] Not a flow/pressure block.");
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

void get_coupled_QP(ComMod& com_mod, double QCoupled[], double QnCoupled[], double PCoupled[], double PnCoupled[])
{
  using namespace consts;

  const auto& cpl = com_mod.cplBC;
  for (int ind = 0; ind < cpl.nSvZeroD_coupled_bc; ind++) {
    const size_t i_path = static_cast<size_t>(ind);
    const std::pair<int, int>& pr = cpl.svZeroD_coupled_bc_idxs[i_path];
    auto& bc = com_mod.eq[pr.first].bc[pr.second];
    if (bc.coupled_bc.get_bc_type() == consts::BoundaryConditionType::bType_Neu) {
      QCoupled[ind] = bc.coupled_bc.get_Qo();
      QnCoupled[ind] = bc.coupled_bc.get_Qn();
      PCoupled[ind] = 0.0;
      PnCoupled[ind] = 0.0;
    } else if (bc.coupled_bc.get_bc_type() == consts::BoundaryConditionType::bType_Dir) {
      QCoupled[ind] = bc.coupled_bc.get_Qo();
      QnCoupled[ind] = bc.coupled_bc.get_Qn();
      PCoupled[ind] = bc.coupled_bc.get_Po();
      PnCoupled[ind] = bc.coupled_bc.get_Pn();
    } else {
      throw std::runtime_error("ERROR: [get_coupled_QP] Invalid Coupled BC type.");
    }
  }
}

void print_svZeroD(int* nSrfs, const std::vector<int>& surfID, double Q[], double P[]) {
  int nParam = 2;
  const char* fileNames[2] = {"Q_svZeroD", "P_svZeroD"};
  std::vector<std::vector<double>> R(nParam, std::vector<double>(*nSrfs));

  if (*nSrfs == 0) return;

  for (int i = 0; i < *nSrfs; ++i) {
    R[0][i] = Q[i];
    R[1][i] = P[i];
  }

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

  build_svzero_coupled_bc_idxs(com_mod);
  if (cplBC.nSvZeroD_coupled_bc == 0) {
    throw std::runtime_error(
        "ERROR: [init_svZeroD] svZeroDSolver is enabled but no Time_dependence Coupled boundaries with "
        "<Coupling_interface> were found.");
  }
  numCoupledSrfs = cplBC.nSvZeroD_coupled_bc;

#ifdef debug_init_svZeroD
  dmsg << "numCoupledSrfs: " << numCoupledSrfs;
#endif

  std::string svzerod_library;
  std::string svzerod_file;
  int ids[2];
  double in_out;

  if (cm.mas(cm_mod)) {
    nsrflistCoupled.clear();
    svzd_blk_names.clear();
    in_out_sign.clear();

    for (int k = 0; k < numCoupledSrfs; k++) {
      nsrflistCoupled.push_back(k);
    }

    if (!solver_interface.has_data) {
      throw std::runtime_error("ERROR: [init_svZeroD] svZeroDSolver interface data is missing.");
    }

#ifdef debug_init_svZeroD
    dmsg << "#### Use XML data #### ";
#endif

    svzerod_library = solver_interface.solver_library;
    svzerod_file = solver_interface.configuration_file;

    const bool have_init_flow = solver_interface.have_initial_flows;
    const double init_flow = solver_interface.initial_flows;
    const bool have_init_press = solver_interface.have_initial_pressures;
    const double init_press = solver_interface.initial_pressures;

    for (int s = 0; s < numCoupledSrfs; ++s) {
      bcType* bc = nullptr;
      if (!nth_coupled_bc(com_mod, s, &bc)) {
        throw std::runtime_error("ERROR: [init_svZeroD] Internal error resolving Coupled BC index.");
      }
      const std::string blk_name = bc->coupled_bc.get_block_name();
      svzd_blk_names.push_back(blk_name);
#ifdef debug_init_svZeroD
      dmsg << "  coupled surface s=" << s << " block='" << blk_name << "'";
#endif
    }

    // Create the svZeroD model
    create_svZeroD_model(svzerod_library, svzerod_file);
    auto interface = interfaces[model_id];
    interface->set_external_step_size(dt);

    // Save IDs of relevant variables in the solution vector
    sol_IDs.assign(2 * numCoupledSrfs, 0);
    for (int s = 0; s < numCoupledSrfs; ++s) {
      get_svZeroD_variable_ids(svzd_blk_names[s], ids, &in_out);
      sol_IDs[2 * s] = ids[0];
      sol_IDs[2 * s + 1] = ids[1];
      in_out_sign.push_back(in_out);

      bcType* bc = nullptr;
      if (!nth_coupled_bc(com_mod, s, &bc)) {
        throw std::runtime_error("ERROR: [init_svZeroD] Internal error resolving Coupled BC for solution IDs.");
      }
      bc->coupled_bc.set_solution_ids(ids[0], ids[1], in_out);
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
      bcType* bc = nullptr;
      if (!nth_coupled_bc(com_mod, s, &bc)) {
        throw std::runtime_error("ERROR: [init_svZeroD] Internal error resolving Coupled BC for initialization.");
      }
      if (bc->coupled_bc.get_bc_type() == consts::BoundaryConditionType::bType_Neu) {
        if (have_init_flow) {
          lpn_state_y[sol_IDs[2 * s]] = init_flow;
        }
        if (have_init_press) {
          lpn_state_y[sol_IDs[2 * s + 1]] = init_press;
          bc->coupled_bc.set_pressure(lpn_state_y[sol_IDs[2 * s + 1]]);
        }
      } else if (bc->coupled_bc.get_bc_type() == consts::BoundaryConditionType::bType_Dir) {
        if (have_init_press) {
          lpn_state_y[sol_IDs[2 * s + 1]] = init_press;
          bc->coupled_bc.set_pressure(lpn_state_y[sol_IDs[2 * s + 1]]);
        }
        if (have_init_flow) {
          lpn_state_y[sol_IDs[2 * s]] = init_flow;
        }
      } else {
        throw std::runtime_error("ERROR: [init_svZeroD] Invalid Coupled BC type.");
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
    // Coupled BCs - broadcast Neumann pressures (one scalar bcast per BC).
    if (cplBC.nSvZeroD_coupled_bc > 0) {
      for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
        auto& eq = com_mod.eq[iEq];
        for (int iBc = 0; iBc < eq.nBc; iBc++) {
          auto& bc = eq.bc[iBc];
          if (utils::btest(bc.bType, iBC_Coupled) &&
              bc.coupled_bc.get_bc_type() == BoundaryConditionType::bType_Neu) {
            bc.coupled_bc.bcast_coupled_neumann_pressure(cm_mod, cm);
          }
        }
      }
    }
  }
}


void calc_svZeroD(ComMod& com_mod, const CmMod& cm_mod, char BCFlag)
{
  using namespace consts;

  auto& cplBC = com_mod.cplBC;
  auto& cm = com_mod.cm;

  const int nCoupledBc = cplBC.nSvZeroD_coupled_bc;

  if (cm.mas(cm_mod)) {
    double QCoupled[numCoupledSrfs], QnCoupled[numCoupledSrfs], PCoupled[numCoupledSrfs], PnCoupled[numCoupledSrfs];
    double params[2];
    double times[2];
    int error_code;
    
    get_coupled_QP(com_mod, QCoupled, QnCoupled, PCoupled, PnCoupled);

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

      // Update pressure and flow in the zeroD model
      for (int i = 0; i < numCoupledSrfs; ++i) {
        bcType* bc = nullptr;
        if (!nth_coupled_bc(com_mod, i, &bc)) {
          throw std::runtime_error("ERROR: [calc_svZeroD] Internal error resolving Coupled BC.");
        }
        const bool is_dirichlet =
            (bc->coupled_bc.get_bc_type() == consts::BoundaryConditionType::bType_Dir);
        double sign = bc->coupled_bc.get_in_out_sign();

        if (is_dirichlet) {
          params[0] = PCoupled[i];
          params[1] = PnCoupled[i];
        } else {
          params[0] = sign * QCoupled[i];
          params[1] = sign * QnCoupled[i];
        }
        update_svZeroD_block_params(svzd_blk_names[i], times, params);
      }

      // Run zeroD simulation
      interface->run_simulation(svZeroDTime, lpn_times, lpn_solutions, error_code);

      // Extract pressure and flow from zeroD solution
      std::copy(lpn_solutions.begin() + (num_output_steps-1)*system_size, lpn_solutions.end(), lpn_state_y.begin());
      
      for (int i = 0; i < numCoupledSrfs; ++i) {
        bcType* bc = nullptr;
        if (!nth_coupled_bc(com_mod, i, &bc)) {
          throw std::runtime_error("ERROR: [calc_svZeroD] Internal error resolving Coupled BC for output.");
        }
        int flow_id = bc->coupled_bc.get_flow_sol_id();
        int pressure_id = bc->coupled_bc.get_pressure_sol_id();
        double in_out = bc->coupled_bc.get_in_out_sign();

        if (flow_id < 0 || flow_id >= system_size || pressure_id < 0 || pressure_id >= system_size) {
          throw std::runtime_error("ERROR: [calc_svZeroD] Invalid solution IDs for Coupled BC: flow_id=" +
                                   std::to_string(flow_id) + ", pressure_id=" + std::to_string(pressure_id) +
                                   ", system_size=" + std::to_string(system_size));
        }

        if (bc->coupled_bc.get_bc_type() == consts::BoundaryConditionType::bType_Neu) {
          PCoupled[i] = lpn_state_y[pressure_id];
          bc->coupled_bc.set_pressure(PCoupled[i]);
        } else if (bc->coupled_bc.get_bc_type() == consts::BoundaryConditionType::bType_Dir) {
          QCoupled[i] = in_out * lpn_state_y[flow_id];
          double Qo_prev = bc->coupled_bc.get_Qn();
          bc->coupled_bc.set_flowrates(Qo_prev, QCoupled[i]);
        } else {
          throw std::runtime_error("ERROR: [calc_svZeroD] Invalid Coupled BC type.");
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

  // If there are multiple procs (not sequential), broadcast coupled Neumann CBC pressures to followers.
  if (!cm.seq()) {
    // Broadcast Coupled BC Neumann pressures (one scalar bcast per BC).
    if (nCoupledBc > 0) {
      for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
        auto& eq = com_mod.eq[iEq];
        for (int iBc = 0; iBc < eq.nBc; iBc++) {
          auto& bc = eq.bc[iBc];
          if (utils::btest(bc.bType, iBC_Coupled) &&
              bc.coupled_bc.get_bc_type() == BoundaryConditionType::bType_Neu) {
            bc.coupled_bc.bcast_coupled_neumann_pressure(cm_mod, cm);
          } else if (utils::btest(bc.bType, iBC_Coupled) &&
                     bc.coupled_bc.get_bc_type() == BoundaryConditionType::bType_Dir) {
            bc.coupled_bc.bcast_coupled_dir_flowrate(cm_mod, cm);
          }
        }
      }
    }
  }
}

}
