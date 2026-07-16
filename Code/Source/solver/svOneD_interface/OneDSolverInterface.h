// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef ONEDSOLVER_INTERFACE_H
#define ONEDSOLVER_INTERFACE_H

#include <dlfcn.h>
#include <string>
#include <stdexcept>

/// @brief Wrapper class for dynamically loading and calling the 1D solver shared library.
///
/// This class loads the svOneDSolver shared library and
/// provides a C++ wrapper around its exported C interface.
///
/// The interface is used to couple 3D Navier-Stokes simulations in
/// svMultiPhysics with reduced-order 1D blood flow models in
/// svOneDSolver.
///
/// Supported coupling modes:
///   - NEU coupling: 3D flow rate -> 1D model, 1D pressure -> 3D solver
///   - DIR coupling: 3D pressure -> 1D model, 1D flow rate -> 3D solver
///
/// Each initialized 1D model is identified by a problem_id assigned
/// by svOneDSolver and used in subsequent library calls.
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
  ///
  /// @param[in] input_file Path to the 1D solver .in file.
  /// @param[out] problem_id Problem identifier assigned by the solver.
  /// @param[out] system_size Total number of DOFs (nodes * 2: flow + area).
  /// @param[in] coupling_type "NEU" or "DIR" coupling direction.
  void initialize(const std::string& input_file, int& problem_id,
                  int& system_size, const std::string& coupling_type);

  /// @brief Synchronize the 1D solver's internal time step with the 3D solver.
  void set_external_step_size(int problem_id, double dt);

  /// @brief Copy the current 1D solution into the caller-provided buffer.
  void return_solution(int problem_id, double* solution, int size);

  /// @brief Push a solution vector into the 1D solver as the current state.
  void update_solution(int problem_id, double* solution, int size);

  /// @brief Advance the 1D solver by one time step.
  ///
  /// @param[in] problem_id   Problem identifier.
  /// @param[in] save_incr    VTK output interval (Increment_in_saving_VTK_files from solver.xml);
  ///                     the 1D library decides internally whether to write output.
  /// @param[in] coupling_type "NEU" or "DIR".
  /// @param[in] params       Array [N, t1, t2, ..., val1, val2, ...] where N=2.
  /// @param[in,out] solution  Solution vector updated after the step.
  /// @param[out] cpl_value    The BC value returned by the 1D solver
  ///                     (pressure for NEU, flow for DIR).
  /// @param[in] last_flag    'L' for the final (committed) iteration, 'D' for
  ///                     derivative / predictor steps.
  /// @param[out] error_code   Non-zero on failure.
  ///                          0 indicates success. Any non-zero value is
  ///                          interpreted by this wrapper as a 1D solver failure.
  void run_simulation(int problem_id, double current_time, int save_incr,
                const std::string& coupling_type, double* params,
                double* solution, double& cpl_value, char last_flag,
                int& error_code);

  /// @brief Retrieve the index within the solution vector that corresponds to
  ///        the coupled boundary DOF.
  void extract_coupled_dof(int problem_id, int& coupled_dof,
                           const std::string& coupling_type);


 private:
  /// @brief Problem identifier assigned by the 1D solver during initialization.
  int problem_id_ = 0;
  /// @brief Total number of DOFs (nodes * 2: flow + area) assigned by the 1D solver during initialization.
  int system_size_ = 0;

  void* library_handle_ = nullptr;

  // Function pointers to shared-library symbols.

  /// @brief Initialize the 1D solver
  void (*initialize_1d_)(const char*, int&, int&, const char*) = nullptr;
  
  /// @brief Set the 1D solver's internal time step size
  void (*set_external_step_size_1d_)(int, double) = nullptr;

  /// @brief Return the current 1D solution vector
  void (*return_1d_solution_)(int, double*, int) = nullptr;

  /// @brief Update the 1D solver's current solution vector
  void (*update_1d_solution_)(int, double*, int) = nullptr;

  /// @brief Advance the 1D solver by one 3D time step
  void (*run_1d_simulation_step_1d_)(int, double, int, const char*, double*,
                                     double*, double&, char*, int&) = nullptr;
  
  /// @brief Retrieve the index within the solution vector that corresponds to
  ///        the coupled boundary DOF.
  void (*extract_coupled_dof_)(int, int&, char*) = nullptr;
};

#endif  // ONEDSOLVER_INTERFACE_H
