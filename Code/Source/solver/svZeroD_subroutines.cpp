/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "svZeroD_subroutines.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include "cmm.h"
#include "consts.h"
#include "svZeroD_interface/LPNSolverInterface.h"

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
  if ((num_inlet_nodes == 0) && (num_outlet_nodes = 1)) {
    blk_ids[0] = IDs[1+num_inlet_nodes*2+1]; // Outlet flow
    blk_ids[1] = IDs[1+num_inlet_nodes*2+2]; // Outlet pressure
    *inlet_or_outlet = 1.0; // Signifies inlet to LPN
  } else if ((num_inlet_nodes == 1) && (num_outlet_nodes == 0)) {
    blk_ids[0] = IDs[1]; // Inlet flow
    blk_ids[1] = IDs[2]; // Inlet pressure
    *inlet_or_outlet = -1.0; // Signifies outlet to LPN
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
  auto& cplBC = com_mod.cplBC;
  int ind = 0;
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
  #define n_debug_init_svZeroD
  #ifdef debug_init_svZeroD
  DebugMsg dmsg(__func__, com_mod.cm.idcm()); 
  dmsg.banner();
  #endif

  auto& cplBC = com_mod.cplBC;
  auto& solver_interface = cplBC.svzerod_solver_interface;
  auto& cm = com_mod.cm;
  double dt = com_mod.dt;

  numCoupledSrfs = cplBC.nFa;
  int nDir = 0;
  int nNeu = 0;
  
  // If this process is the master process on the communicator
  if (cm.mas(cm_mod)) {
    for (int iFa = 0; iFa < cplBC.nFa; iFa++) {
      auto& fa = cplBC.fa[iFa];

      if (fa.bGrp == consts::CplBCType::cplBC_Dir) {
        nsrflistCoupled.push_back(iFa);
        nDir = nDir + 1;
      }
    }
    for (int iFa = 0; iFa < cplBC.nFa; iFa++) {
      auto& fa = cplBC.fa[iFa];

      if (fa.bGrp == consts::CplBCType::cplBC_Neu) {
        nsrflistCoupled.push_back(iFa);
        nNeu = nNeu + 1;
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
  std::vector<std::string> svzd_blk_names_unsrtd(numCoupledSrfs);
  std::vector<int> svzd_blk_ids(numCoupledSrfs);
  int init_flow_flag, init_press_flag;
  double init_flow, init_press, in_out;
  
  if (cm.mas(cm_mod)) {

    if (solver_interface.has_data) { 
      #ifdef debug_init_svZeroD
      dmsg << "#### Use XML data #### " << " ";
      #endif

      svzerod_library = solver_interface.solver_library;
      svzerod_file = solver_interface.configuration_file;

      int i = 0;
      for (const auto& pair : solver_interface.block_surface_map) {
        #ifdef debug_init_svZeroD
        dmsg << "block_surface_map: '" + pair.first << "'";
        #endif
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
              pair.first + "' and surface '" + pair.second + "'; check th Block_to_surface_map solver XML parameter.");
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
      #ifdef debug_init_svZeroD
      dmsg << ">>> s " << s;
      dmsg << "  nsrflistCoupled[s]: " << nsrflistCoupled[s];
      #endif

      for (int t = 0; t < numCoupledSrfs; ++t) {
        #ifdef debug_init_svZeroD
        dmsg << "  >>> t " << t;
        dmsg << "    svzd_blk_ids[t]: " << svzd_blk_ids[t];
        #endif
        if (svzd_blk_ids[t] == nsrflistCoupled[s]) {
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

      if (found == 0) {
        throw std::runtime_error("ERROR: Did not find block name for surface ID: " + nsrflistCoupled[s]);
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
      if (init_flow_flag == 1) {
        lpn_state_y[sol_IDs[2 * s]] = init_flow;
        cplBC.fa[s].y = lpn_state_y[sol_IDs[2 * s]];
      }
      if (init_press_flag == 1) {
        lpn_state_y[sol_IDs[2 * s + 1]] = init_press;
        cplBC.fa[s].y = lpn_state_y[sol_IDs[2 * s + 1]];
      }
    }
    std::copy(lpn_state_y.begin(), lpn_state_y.end(), last_state_y.begin());

    if (writeSvZeroD == 1) {
      // Initialize output file
      int flag = 0;
      write_svZeroD_solution(&svZeroDTime, lpn_state_y, &flag);
    }
  }

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


void calc_svZeroD(ComMod& com_mod, const CmMod& cm_mod, char BCFlag) {
  int nDir = 0;
  int nNeu = 0;
  double dt = com_mod.dt;
  auto& cplBC = com_mod.cplBC;
  auto& cm = com_mod.cm;

  // If this process is the master process on the communicator
  if (cm.mas(cm_mod)) {
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
        if (i < nDir) {
          params[0] = PCoupled[i];
          params[1] = PnCoupled[i];
        } else {
          params[0] = in_out_sign[i] * QCoupled[i];
          params[1] = in_out_sign[i] * QnCoupled[i];
          total_flow += QCoupled[i];
        }
        update_svZeroD_block_params(svzd_blk_names[i], times, params);
      }

      // Run zeroD simulation
      interface->run_simulation(svZeroDTime, lpn_times, lpn_solutions, error_code);

      // Extract pressure and flow from zeroD solution
      std::copy(lpn_solutions.begin() + (num_output_steps-1)*system_size, lpn_solutions.end(), lpn_state_y.begin());
      
      for (int i = 0; i < numCoupledSrfs; ++i) {
        if (i < nDir) {
          QCoupled[i] = in_out_sign[i] * lpn_state_y[sol_IDs[2 * i]];
          cplBC.fa[i].y = QCoupled[i];
        } else {
          PCoupled[i] = lpn_state_y[sol_IDs[2 * i + 1]];
          cplBC.fa[i].y = PCoupled[i];
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
};
