// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "Array.h"
#include "Vector.h"
#include "Simulation.h"

/**
 * @brief Solution state variables container
 *
 * Contains solution arrays at two time levels (n and n+1) for time integration
 */
struct soluStateVars {
  Array<double> An, Dn, Yn;  // New (n+1): acceleration, displacement, velocity at next time step
  Array<double> Ao, Do, Yo;  // Old (n): acceleration, displacement, velocity at current time step
};

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
   * @param Ao Old acceleration array (takes ownership via move)
   * @param Do Old displacement array (takes ownership via move)
   * @param Yo Old velocity array (takes ownership via move)
   */
  Integrator(Simulation* simulation, Array<double>&& Ao, Array<double>&& Do, Array<double>&& Yo);

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
   * @brief Perform predictor step for next time step
   *
   * Performs predictor step using generalized-alpha method to estimate
   * solution at n+1 time level based on current solution at n time level.
   * This should be called once per time step before the Newton iteration loop.
   */
  void predictor();

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

  /**
   * @brief Get reference to An (new time derivative of variables at n+1)
   *
   * @return Reference to An array (acceleration at next time step)
   */
  Array<double>& get_An() { return solu_state_vars_.An; }

  /**
   * @brief Get reference to Dn (new integrated variables at n+1)
   *
   * @return Reference to Dn array (displacement at next time step)
   */
  Array<double>& get_Dn() { return solu_state_vars_.Dn; }

  /**
   * @brief Get reference to Yn (new variables at n+1)
   *
   * @return Reference to Yn array (velocity at next time step)
   */
  Array<double>& get_Yn() { return solu_state_vars_.Yn; }

  /**
   * @brief Get reference to Ao (old time derivative of variables at n)
   *
   * @return Reference to Ao array (acceleration at current time step)
   */
  Array<double>& get_Ao() { return solu_state_vars_.Ao; }

  /**
   * @brief Get reference to Do (old integrated variables at n)
   *
   * @return Reference to Do array (displacement at current time step)
   */
  Array<double>& get_Do() { return solu_state_vars_.Do; }

  /**
   * @brief Get reference to Yo (old variables at n)
   *
   * @return Reference to Yo array (velocity at current time step)
   */
  Array<double>& get_Yo() { return solu_state_vars_.Yo; }

  /**
   * @brief Get reference to solution state variables struct
   *
   * @return Reference to soluStateVars struct containing all solution arrays
   */
  soluStateVars& get_solu_state_vars() { return solu_state_vars_; }

private:
  /** @brief Pointer to the simulation object */
  Simulation* simulation_;

  /** @brief Time derivative of variables (acceleration in structural mechanics) */
  Array<double> Ag_;

  /** @brief Variables (velocity in structural mechanics) */
  Array<double> Yg_;

  /** @brief Integrated variables (displacement in structural mechanics) */
  Array<double> Dg_;

  /** @brief Solution state variables (An, Dn, Yn, Ao, Do, Yo) */
  soluStateVars solu_state_vars_;

  /** @brief Residual vector for face-based quantities */
  Vector<double> res_;

  /** @brief Increment flag for faces in linear solver */
  Vector<int> incL_;

  /** @brief Newton iteration counter for current time step */
  int newton_count_;

  /** @brief Debug output suffix string combining time step and iteration number */
  std::string istr_;

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

  /**
   * @brief Initiator function for generalized-alpha method (initiator)
   *
   * Computes solution variables at intermediate time levels using
   * generalized-alpha parameters (am, af) for time integration.
   * Updates Ag, Yg, Dg based on An, Ao, Yn, Yo, Dn, Do.
   *
   * @param Ag Time derivative array at generalized-alpha level
   * @param Yg Solution variable array at generalized-alpha level
   * @param Dg Integrated variable array at generalized-alpha level
   */
  void initiator(Array<double>& Ag, Array<double>& Yg, Array<double>& Dg);

  /**
   * @brief Corrector function with convergence check (corrector)
   *
   * Updates solution at n+1 time level and checks convergence of Newton
   * iterations. Also handles equation switching for coupled problems.
   */
  void corrector();

  /**
   * @brief Pressure correction for Taylor-Hood elements (corrector_taylor_hood)
   *
   * Interpolates pressure at edge nodes using reduced basis applied
   * on element vertices for Taylor-Hood type elements.
   */
  void corrector_taylor_hood();
};

#endif // INTEGRATOR_H
