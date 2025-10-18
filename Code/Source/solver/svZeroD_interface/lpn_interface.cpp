// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "LPNSolverInterface.h"

#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <string>
#include <iterator>

static std::map<int,LPNSolverInterface*> interfaces;


//-------------------------
// Load the svZeroD library and initialize the 0D model
//-------------------------
//
void lpn_interface_add_model_(const char* lpn_library_name, int* lib_filename_len, const char* lpn_json_file, int* json_filename_len, int* interface_id, int* num_output_steps, int* system_size)
{
  // Load library
  std::string lpn_library(lpn_library_name,0,*lib_filename_len);
  std::cout<<"lpn_library: "<<lpn_library<<std::endl;
  auto interface = new LPNSolverInterface();
  interface->load_library(lpn_library);
  std::cout << "[lpn_interface_init_] Loaded library." << std::endl;

  // Initialize model
  std::string lpn_file(lpn_json_file,0,*json_filename_len);
  interface->initialize(std::string(lpn_file));
  *interface_id = interface->problem_id_;
  std::cout << "[lpn_interface_add_model_] interface->problem_id_:" <<interface->problem_id_<< std::endl;
  interfaces[*interface_id] = interface;
  
  // Save model parameters
  *num_output_steps = interface->num_output_steps_;
  *system_size = interface->system_size_;
}

//----------------------------
// Set the time step size for the external solver in the 0D solver
//----------------------------
//
extern "C" void lpn_interface_set_external_step_size_(const int* interface_id, double* step_size)
{
  auto interface = interfaces[*interface_id];
  interface->set_external_step_size(*step_size);
}

//----------------------------
// Get the IDs of the inlet/outlet variables of a given block in the state vector
//----------------------------
//
extern "C" void lpn_interface_get_variable_ids_(const int* interface_id, const char* block_name, int* block_name_len, int* blk_ids, double* inlet_or_outlet)
{
  auto interface = interfaces[*interface_id];
  //std::vector<std::string> variable_names;
  //variable_names = interface->variable_names_;
  std::vector<int> IDs;
  std::string block_name_cpp(block_name,0,*block_name_len);
  interface->get_block_node_IDs(std::string(block_name_cpp), IDs);
  // IDs in the above function stores info in the following format:
  // {num inlet nodes, inlet flow[0], inlet pressure[0],..., num outlet nodes, outlet flow[0], outlet pressure[0],...}
  int num_inlet_nodes = IDs[0];
  int num_outlet_nodes = IDs[1+num_inlet_nodes*2];
  if ((num_inlet_nodes == 0) && (num_outlet_nodes = 1)) {
    //std::cout<<"Only outlet nodes"<<std::endl;
    blk_ids[0] = IDs[1+num_inlet_nodes*2+1]; // Outlet flow
    blk_ids[1] = IDs[1+num_inlet_nodes*2+2]; // Outlet pressure
    *inlet_or_outlet = 1.0; // Signifies inlet to LPN
  } else if ((num_inlet_nodes == 1) && (num_outlet_nodes == 0)) {
    //std::cout<<"Only inlet nodes"<<std::endl;
    blk_ids[0] = IDs[1]; // Inlet flow
    blk_ids[1] = IDs[2]; // Inlet pressure
    *inlet_or_outlet = -1.0; // Signifies outlet to LPN
  } else {
    std::runtime_error("ERROR: [lpn_interface_get_variable_ids] Not a flow/pressure block.");
  }
}

//----------------------------
// Increment the 0D solver by one time step
//----------------------------
//
extern "C" void lpn_interface_get_solution_(const int* interface_id, const double* time, double* solution)
{
  auto interface = interfaces[*interface_id];
  std::vector<double> lpn_solution(interface->system_size_);
  interface->increment_time(*time, lpn_solution);

  for (int i = 0; i < lpn_solution.size(); i++) {
    solution[i] = lpn_solution[i];
  }
}

//----------------------------
// Overwrite the y and ydot state vectors in the 0D solver
//----------------------------
//
extern "C" void lpn_interface_update_state_(const int* interface_id, std::vector<double>& y, std::vector<double>& ydot)
{
  auto interface = interfaces[*interface_id];
  /*std::vector<double> state_y(interface->system_size_);
  std::vector<double> state_ydot(interface->system_size_);
  // Convert arrays to std::vector
  for (int i = 0; i < interface->system_size_; i++) {
    state_y[i] = y[i];
    state_ydot[i] = ydot[i];
  }*/
  //interface->update_state(state_y, state_ydot);
  interface->update_state(y, ydot);
}

//----------------------------
// Return the y state vector
//----------------------------
//
extern "C" void lpn_interface_return_y_(const int* interface_id, std::vector<double>& y)
{
  auto interface = interfaces[*interface_id];
  //std::vector<double> state_y(interface->system_size_);
  interface->return_y(y);
  // Convert std::vector to array for fortran
  /*for (int i = 0; i < interface->system_size_; i++) {
    y[i] = state_y[i];
  }*/
}

//----------------------------
// Return the ydot state vector
//----------------------------
//
extern "C" void lpn_interface_return_ydot_(const int* interface_id, std::vector<double>& ydot)
{
  auto interface = interfaces[*interface_id];
  //std::vector<double> state_ydot(interface->system_size_);
  interface->return_ydot(ydot);
  // Convert std::vector to array for fortran
  /*for (int i = 0; i < interface->system_size_; i++) {
    ydot[i] = state_ydot[i];
  }*/
}

//----------------------------
// Update the parameters of a particular 0D block
//----------------------------
//
extern "C" void lpn_interface_update_block_params_(const int* interface_id, const char* block_name, int* block_name_len, double* time, double* params, int* num_time_pts)
{
  auto interface = interfaces[*interface_id];
  int param_len = *num_time_pts; // Usually 2 for this use case
  std::vector<double> new_params(1+2*param_len);
  // Format of new_params for flow/pressure blocks: 
  // [N, time_1, time_2, ..., time_N, value1, value2, ..., value_N]
  // where N is number of time points and value* is flow/pressure
  new_params[0] = (double) param_len;
  for (int i = 0; i < param_len; i++) {
    new_params[1+i] = time[i];
    new_params[1+param_len+i] = params[i];
  }
  std::string block_name_cpp(block_name,0,*block_name_len);
  //std::cout<<"[lpn_interface_update_block_params_] new_params = "<<new_params[0]<<" "<<new_params[1]<<" "<<new_params[2]<<" "<<new_params[3]<<" "<<new_params[4]<<std::endl;
  interface->update_block_params(std::string(block_name_cpp), new_params);
}

//----------------------------
// Run a full 0D simulation
//----------------------------
//
extern "C" void lpn_interface_run_simulation_(const int* interface_id, const double* time, std::vector<double>& lpn_times, std::vector<double>& lpn_solutions, int* error_code)
{
  auto interface = interfaces[*interface_id];
  int solutions_vec_size = interface->system_size_*interface->num_output_steps_;
  //std::vector<double> lpn_solutions_vec(solutions_vec_size);
  //std::vector<double> lpn_times_vec(interface->num_output_steps_);
  int error_code_ret = 0;
  //interface->run_simulation(*time, lpn_times_vec, lpn_solutions_vec, error_code_ret);
  interface->run_simulation(*time, lpn_times, lpn_solutions, error_code_ret);
  *error_code = error_code_ret;

  /*for (int i = 0; i < interface->num_output_steps_; i++) {
    lpn_times[i] = lpn_times_vec[i];
  }
  for (int i = 0; i < solutions_vec_size; i++) {
    lpn_solutions[i] = lpn_solutions_vec[i];
  }*/
}

//----------------------------
// Write out the solution vector to a runing file
//----------------------------
//
extern "C" void lpn_interface_write_solution_(const int* interface_id, const double* lpn_time, std::vector<double>& lpn_solution, int* flag)
{
  auto interface = interfaces[*interface_id];
  int sys_size = interface->system_size_;
  if (*flag == 0) { // Initialize output file: Write header with variable names
    std::vector<std::string> variable_names;
    variable_names = interface->variable_names_;
    std::ofstream out_file;
    out_file.open("svZeroD_data", std::ios::out | std::ios::app);
    out_file<<sys_size<<" ";
    for (int i = 0; i < sys_size; i++) {
      out_file<<static_cast<std::string>(variable_names[i])<<" ";
    }
    out_file<<'\n';
  } else {
    std::ofstream out_file;
    out_file.open("svZeroD_data", std::ios::out | std::ios::app);
    out_file<<*lpn_time<<" ";
    for (int i = 0; i < sys_size; i++) {
      out_file<<lpn_solution[i]<<" ";
    }
    out_file<<'\n';
    out_file.close();
  }
}
