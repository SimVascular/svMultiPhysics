// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "Array.h"
#include "Vector.h"
#include "Simulation.h"

/**
 * @brief Integrator class encapsulates the Newton iteration loop for time integration
 *
 * This class handles the nonlinear Newton iteration scheme for solving coupled
 * multi-physics equations in svMultiPhysics. It manages:
 * - Solution variables (Ag, Yg, Dg) at generalized-alpha time levels
 * - Newton iteration loop with convergence checking
 * - Linear system assembly and solve
 * - Boundary condition application
 *
 * Related to GitHub issue #442: Encapsulate the Newton iteration in main.cpp
 */
class Integrator {

public:
  /**
   * @brief Construct a new Integrator object
   *
   * @param simulation Pointer to the Simulation object containing problem data
   */
  Integrator(Simulation* simulation);

  /**
   * @brief Destroy the Integrator object
   */
  ~Integrator();

  /**
   * @brief Execute one time step with Newton iteration loop
   *
   * Performs the complete Newton iteration sequence including initialization,
   * assembly, boundary condition application, linear solve, and convergence check.
   *
   * @return True if all equations converged, false otherwise
   */
  bool step();

  /**
   * @brief Get reference to solution variable Ag (time derivative of variables)
   *
   * @return Reference to Ag array (acceleration in structural mechanics)
   */
  Array<double>& get_Ag() { return Ag_; }

  /**
   * @brief Get reference to solution variable Yg (variables)
   *
   * @return Reference to Yg array (velocity in structural mechanics)
   */
  Array<double>& get_Yg() { return Yg_; }

  /**
   * @brief Get reference to solution variable Dg (integrated variables)
   *
   * @return Reference to Dg array (displacement in structural mechanics)
   */
  Array<double>& get_Dg() { return Dg_; }

private:
  /** @brief Pointer to the simulation object */
  Simulation* simulation_;

  /** @brief Time derivative of variables (acceleration in structural mechanics) */
  Array<double> Ag_;

  /** @brief Variables (velocity in structural mechanics) */
  Array<double> Yg_;

  /** @brief Integrated variables (displacement in structural mechanics) */
  Array<double> Dg_;

  /** @brief Residual vector for face-based quantities */
  Vector<double> res_;

  /** @brief Increment flag for faces in linear solver */
  Vector<int> incL_;

  /** @brief Newton iteration counter for current time step */
  int inner_count_;

  /**
   * @brief Initialize solution arrays for Ag, Yg, Dg based on problem size
   */
  void initialize_arrays();

  /**
   * @brief Perform initiator step for Generalized-alpha Method
   *
   * Computes quantities at intermediate time levels (n+alpha_m, n+alpha_f)
   */
  void initiator_step();

  /**
   * @brief Allocate right-hand side (RHS) and left-hand side (LHS) arrays
   *
   * @param eq Reference to the equation being solved
   */
  void allocate_linear_system(eqType& eq);

  /**
   * @brief Set body forces for the current time step
   */
  void set_body_forces();

  /**
   * @brief Assemble global equations for all meshes
   */
  void assemble_equations();

  /**
   * @brief Apply all boundary conditions (Neumann, Dirichlet, CMM, contact, etc.)
   */
  void apply_boundary_conditions();

  /**
   * @brief Solve the assembled linear system
   */
  void solve_linear_system();

  /**
   * @brief Perform corrector step and check convergence of all equations
   *
   * @return True if all equations converged, false otherwise
   */
  bool corrector_and_check_convergence();

  /**
   * @brief Update residual and increment arrays for linear solver
   *
   * @param eq Reference to the equation being solved
   */
  void update_residual_arrays(eqType& eq);
};

#endif // INTEGRATOR_H
