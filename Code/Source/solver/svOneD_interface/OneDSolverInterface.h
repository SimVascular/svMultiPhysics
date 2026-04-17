// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef ONEDSOLVER_INTERFACE_H
#define ONEDSOLVER_INTERFACE_H

#include <dlfcn.h>
#include <string>
#include <stdexcept>

/// @brief Wrapper class for dynamically loading and calling the 1D solver shared library.
///
/// This class uses dlopen/dlsym to load the libsvoned_interface shared library at runtime
/// and provides C++ method wrappers for each of its exported C functions.
///
/// Shared library functions:
///   - initialize_1d(input_file, problem_id, system_size, coupling_type)
///   - set_external_step_size_1d(problem_id, dt)
///   - return_1d_solution(problem_id, solution, size)
///   - update_1d_solution(problem_id, solution, size)
///   - run_1d_simulation_step_1d(problem_id, time, save_flag, coupling_type,
///                               params, solution, cpl_value, error_code)
///   - extract_coupled_dof(problem_id, coupled_dof, coupling_type)
//
class OneDSolverInterface {
 public:
  OneDSolverInterface() = default;
  ~OneDSolverInterface();

  /// @brief Load the 1D solver shared library from the given path.
  void load_library(const std::string& interface_lib);

  /// @brief Initialize the 1D solver from an input file.
  /// @param input_file Path to the 1D solver .in file.
  /// @param problem_id Output: problem identifier assigned by the solver.
  /// @param system_size Output: total number of DOFs (nodes * 2: flow + area).
  /// @param coupling_type "NEU" or "DIR" coupling direction.
  void initialize(const std::string& input_file, int& problem_id,
                  int& system_size, const std::string& coupling_type);

  /// @brief Synchronize the 1D solver's internal time step with the 3D solver.
  void set_external_step_size(int problem_id, double dt);

  /// @brief Copy the current 1D solution into the caller-provided buffer.
  void return_solution(int problem_id, double* solution, int size);

  /// @brief Push a solution vector into the 1D solver as the current state.
  void update_solution(int problem_id, double* solution, int size);

  /// @brief Advance the 1D solver by one time step.
  /// @param problem_id   Problem identifier.
  /// @param current_time Current simulation time (start of the step).
  /// @param save_flag    Non-zero to write VTK output for this step.
  /// @param coupling_type "NEU" or "DIR".
  /// @param params       Array [N, t1, t2, ..., val1, val2, ...] where N=2.
  /// @param solution     In/out: solution vector updated after the step.
  /// @param cpl_value    Output: the BC value returned by the 1D solver
  ///                     (pressure for NEU, flow for DIR).
  /// @param error_code   Output: non-zero on failure.
  void run_step(int problem_id, double current_time, int save_flag,
                const std::string& coupling_type, double* params,
                double* solution, double& cpl_value, int& error_code);

  /// @brief Retrieve the index within the solution vector that corresponds to
  ///        the coupled boundary DOF.
  void extract_coupled_dof(int problem_id, int& coupled_dof,
                           const std::string& coupling_type);

  // Public data members set after initialize().
  int problem_id_ = 0;
  int system_size_ = 0;

 private:
  void* library_handle_ = nullptr;

  // Function pointers to shared-library symbols.
  void (*initialize_1d_)(const char*, int&, int&, const char*) = nullptr;
  void (*set_external_step_size_1d_)(int, double) = nullptr;
  void (*return_1d_solution_)(int, double*, int) = nullptr;
  void (*update_1d_solution_)(int, double*, int) = nullptr;
  void (*run_1d_simulation_step_1d_)(int, double, int, const char*, double*,
                                     double*, double&, int&) = nullptr;
  void (*extract_coupled_dof_)(int, int&, char*) = nullptr;
};

#endif  // ONEDSOLVER_INTERFACE_H
