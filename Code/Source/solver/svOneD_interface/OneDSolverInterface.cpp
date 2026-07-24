// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "OneDSolverInterface.h"

#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>
#include "Core/Exception.h"

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
    throw svmp::CoreException(
        std::string("[OneDSolverInterface] Could not load shared library '") +
        interface_lib + "': " + dlerror(),
        svmp::StatusCode::DependencyError,
        __FILE__,
        __LINE__,
        __func__);
  }

  // Clear any existing error.
  dlerror();

  auto load_sym = [&](const char* name) -> void* {
    void* sym = dlsym(library_handle_, name);
    const char* err = dlerror();
    if (err) {
      throw svmp::CoreException(
          std::string("[OneDSolverInterface] Could not find the 1D interface function '") +
          name + "' in the shared library '" + interface_lib + "': " + err,
          svmp::StatusCode::DependencyError,
          __FILE__,
          __LINE__,
          __func__);
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
    throw svmp::CoreException(
        "[OneDSolverInterface] Cannot initialize the 1D solver because the "
        "1D interface library has not been loaded. Call load_library() with "
        "the svOneDSolver interface shared library before initialize(), and "
        "check that the library path is correct.",
        svmp::StatusCode::DependencyError,
        __FILE__,
        __LINE__,
        __func__);
  }
  initialize_1d_(input_file.c_str(), problem_id, system_size,
                 coupling_type.c_str());
  problem_id_  = problem_id;
  system_size_ = system_size;
}

void OneDSolverInterface::set_external_step_size(int problem_id, double dt)
{
  if (!set_external_step_size_1d_) {
    throw svmp::CoreException(
        "[OneDSolverInterface] set_external_step_size_1d not loaded",
        svmp::StatusCode::DependencyError,
        __FILE__,
        __LINE__,
        __func__);
  }
  set_external_step_size_1d_(problem_id, dt);
}

void OneDSolverInterface::return_solution(int problem_id, double* solution, int size)
{
  if (!return_1d_solution_) {
    throw svmp::CoreException(
        "[OneDSolverInterface] return_1d_solution not loaded",
        svmp::StatusCode::DependencyError,
        __FILE__,
        __LINE__,
        __func__);
  }
  return_1d_solution_(problem_id, solution, size);
}

void OneDSolverInterface::update_solution(int problem_id, double* solution, int size)
{
  if (!update_1d_solution_) {
    throw svmp::CoreException(
        "[OneDSolverInterface] update_1d_solution not loaded",
        svmp::StatusCode::DependencyError,
        __FILE__,
        __LINE__,
        __func__);
  }
  update_1d_solution_(problem_id, solution, size);
}

void OneDSolverInterface::run_simulation(int problem_id, double current_time,
                                   int save_incr,
                                   const std::string& coupling_type,
                                   double* params, double* solution,
                                   double& cpl_value, char last_flag,
                                   int& error_code)
{
  if (!run_1d_simulation_step_1d_) {
    throw svmp::CoreException(
        "[OneDSolverInterface] Cannot run the 1D solver because the "
        "run_1d_simulation_step_1d interface function is not available. "
        "Check that the svOneDSolver interface shared library was loaded correctly.",
        svmp::StatusCode::DependencyError,
        __FILE__,
        __LINE__,
        __func__);
  }

  std::vector<char> ctype_buf(coupling_type.begin(), coupling_type.end());
  ctype_buf.push_back('\0');

  char flag_buf[2] = { last_flag, '\0' };

  error_code = 0;

  run_1d_simulation_step_1d_(problem_id, current_time, save_incr,
                             ctype_buf.data(), params, solution,
                             cpl_value, flag_buf, error_code);

  if (error_code != 0) {
    throw svmp::CoreException(
        "[OneDSolverInterface] svOneDSolver failed while advancing the 1D model. "
        "The 1D solver returned error code " + std::to_string(error_code) + ".",
        svmp::StatusCode::DependencyError,
        __FILE__,
        __LINE__,
        __func__);
  }
}

void OneDSolverInterface::extract_coupled_dof(int problem_id, int& coupled_dof,
                                              const std::string& coupling_type)
{
  if (!extract_coupled_dof_) {
    throw svmp::CoreException(
        "[OneDSolverInterface] extract_coupled_dof not loaded",
        svmp::StatusCode::DependencyError,
        __FILE__,
        __LINE__,
        __func__);
  }
  // Copy into a mutable buffer; the shared-library function signature uses
  // char* (not const char*) so we must pass a writable copy.
  std::vector<char> buf(coupling_type.begin(), coupling_type.end());
  buf.push_back('\0');
  extract_coupled_dof_(problem_id, coupled_dof, buf.data());
}
