// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "OneDSolverInterface.h"

#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

OneDSolverInterface::~OneDSolverInterface()
{
  if (library_handle_) {
    dlclose(library_handle_);
    library_handle_ = nullptr;
  }
}

void OneDSolverInterface::load_library(const std::string& interface_lib)
{
  library_handle_ = dlopen(interface_lib.c_str(), RTLD_LAZY);
  if (!library_handle_) {
    throw std::runtime_error(std::string("[OneDSolverInterface] Could not load shared library '") +
                             interface_lib + "': " + dlerror());
  }

  // Clear any existing error.
  dlerror();

  auto load_sym = [&](const char* name) -> void* {
    void* sym = dlsym(library_handle_, name);
    const char* err = dlerror();
    if (err) {
      throw std::runtime_error(std::string("[OneDSolverInterface] Could not load symbol '") +
                               name + "': " + err);
    }
    return sym;
  };

  *(void**)(&initialize_1d_)              = load_sym("initialize_1d");
  *(void**)(&set_external_step_size_1d_)  = load_sym("set_external_step_size_1d");
  *(void**)(&return_1d_solution_)         = load_sym("return_1d_solution");
  *(void**)(&update_1d_solution_)         = load_sym("update_1d_solution");
  *(void**)(&run_1d_simulation_step_1d_)  = load_sym("run_1d_simulation_step_1d");
  *(void**)(&extract_coupled_dof_)        = load_sym("extract_coupled_dof");
}

void OneDSolverInterface::initialize(const std::string& input_file,
                                     int& problem_id,
                                     int& system_size,
                                     const std::string& coupling_type)
{
  if (!initialize_1d_) {
    throw std::runtime_error("[OneDSolverInterface] initialize_1d not loaded");
  }
  initialize_1d_(input_file.c_str(), problem_id, system_size,
                 coupling_type.c_str());
  problem_id_  = problem_id;
  system_size_ = system_size;
}

void OneDSolverInterface::set_external_step_size(int problem_id, double dt)
{
  if (!set_external_step_size_1d_) {
    throw std::runtime_error("[OneDSolverInterface] set_external_step_size_1d not loaded");
  }
  set_external_step_size_1d_(problem_id, dt);
}

void OneDSolverInterface::return_solution(int problem_id, double* solution, int size)
{
  if (!return_1d_solution_) {
    throw std::runtime_error("[OneDSolverInterface] return_1d_solution not loaded");
  }
  return_1d_solution_(problem_id, solution, size);
}

void OneDSolverInterface::update_solution(int problem_id, double* solution, int size)
{
  if (!update_1d_solution_) {
    throw std::runtime_error("[OneDSolverInterface] update_1d_solution not loaded");
  }
  update_1d_solution_(problem_id, solution, size);
}

void OneDSolverInterface::run_step(int problem_id, double current_time,
                                   int save_incr,
                                   const std::string& coupling_type,
                                   double* params, double* solution,
                                   double& cpl_value, char last_flag,
                                   int& error_code)
{
  if (!run_1d_simulation_step_1d_) {
    throw std::runtime_error("[OneDSolverInterface] run_1d_simulation_step_1d not loaded");
  }
  // Copy coupling_type into a mutable buffer (shared-library uses char*).
  std::vector<char> ctype_buf(coupling_type.begin(), coupling_type.end());
  ctype_buf.push_back('\0');
  // Copy last_flag into a mutable single-character buffer.
  char flag_buf[2] = { last_flag, '\0' };
  run_1d_simulation_step_1d_(problem_id, current_time, save_incr,
                              ctype_buf.data(), params, solution,
                              cpl_value, flag_buf, error_code);
}

void OneDSolverInterface::extract_coupled_dof(int problem_id, int& coupled_dof,
                                              const std::string& coupling_type)
{
  if (!extract_coupled_dof_) {
    throw std::runtime_error("[OneDSolverInterface] extract_coupled_dof not loaded");
  }
  // Copy into a mutable buffer; the shared-library function signature uses
  // char* (not const char*) so we must pass a writable copy.
  std::vector<char> buf(coupling_type.begin(), coupling_type.end());
  buf.push_back('\0');
  extract_coupled_dof_(problem_id, coupled_dof, buf.data());
}
