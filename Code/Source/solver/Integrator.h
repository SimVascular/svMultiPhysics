// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "Array.h"
#include "SolutionStates.h"
#include "Vector.h"
#include "Simulation.h"

#include <functional>

/// @brief Newmark time integration utilities.
///
/// Compute consistent state variables from a prescribed displacement or
/// velocity using the Newmark-beta / generalized-alpha relationships:
///   Dn = Do + dt*Yo + dt^2*((0.5-beta)*Ao + beta*An)
///   Yn = Yo + dt*((1-gamma)*Ao + gamma*An)
namespace newmark {

/// Compute acceleration and velocity from prescribed displacement.
inline void state_from_displacement(
    double d_new, double d_old, double v_old, double a_old,
    double dt, double beta, double gam,
    double& a_new, double& v_new)
{
  a_new = (d_new - d_old - dt * v_old) / (beta * dt * dt)
        - (0.5 - beta) / beta * a_old;
  v_new = v_old + dt * ((1.0 - gam) * a_old + gam * a_new);
}

/// Compute acceleration from prescribed velocity.
inline void state_from_velocity(
    double v_new, double v_old, double a_old,
    double dt, double gam,
    double& a_new)
{
  a_new = (v_new - v_old) / (gam * dt) - (1.0 - gam) / gam * a_old;
}

} // namespace newmark

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
   * @param solutions Solution states containing old time level arrays (takes ownership via move)
   */
  Integrator(Simulation* simulation, SolutionStates&& solutions);

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
   * @brief Solve a single equation to convergence (for partitioned coupling)
   *
   * Runs Newton iterations for only the specified equation, without cycling
   * to other equations. Used by partitioned FSI to solve fluid, solid, and
   * mesh equations independently.
   *
   * @param iEq Index of the equation to solve
   * @param post_assembly Optional callback invoked after boundary condition
   *        application but before the linear solve. Used by partitioned FSI
   *        to inject interface traction into the residual.
   * @return True if the equation converged, false if max iterations reached
   */
  bool step_equation(int iEq, std::function<void()> post_assembly = nullptr);

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
  Array<double>& get_Ag() { return solutions_.intermediate.get_acceleration(); }
  const Array<double>& get_Ag() const { return solutions_.intermediate.get_acceleration(); }

  /**
   * @brief Get reference to solution variable Yg (variables)
   *
   * @return Reference to Yg array (velocity in structural mechanics)
   */
  Array<double>& get_Yg() { return solutions_.intermediate.get_velocity(); }
  const Array<double>& get_Yg() const { return solutions_.intermediate.get_velocity(); }

  /**
   * @brief Get reference to solution variable Dg (integrated variables)
   *
   * @return Reference to Dg array (displacement in structural mechanics)
   */
  Array<double>& get_Dg() { return solutions_.intermediate.get_displacement(); }
  const Array<double>& get_Dg() const { return solutions_.intermediate.get_displacement(); }

  /**
   * @brief Get reference to solution states struct
   *
   * Provides access to all solution arrays at old (n) and current (n+1) time levels.
   * Use this to access An, Dn, Yn (current) and Ao, Do, Yo (old) via:
   *   auto& solutions = integrator.get_solutions();
   *   auto& An = solutions.current.get_acceleration();
   *   auto& Do = solutions.old.get_displacement();
   *
   * @return Reference to SolutionStates struct containing all solution arrays
   */
  SolutionStates& get_solutions() { return solutions_; }
  const SolutionStates& get_solutions() const { return solutions_; }

private:
  /** @brief Pointer to the simulation object */
  Simulation* simulation_;

  /** @brief Solution states at old, current, and intermediate time levels */
  SolutionStates solutions_;

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
   * Updates solutions.intermediate (Ag, Yg, Dg) based on solutions.current
   * (An, Yn, Dn) and solutions.old (Ao, Yo, Do).
   *
   * @param solutions Solution states containing old, current, and intermediate levels
   */
  void initiator(SolutionStates& solutions);

  /**
   * @brief Update solution and check convergence of current equation
   *
   * Performs the corrector step: updates An, Yn, Dn from the linear solve
   * result and checks convergence norms. Sets eq.ok if converged.
   * Does NOT handle equation cycling -- that is done separately in corrector().
   */
  void update_solution();

  /**
   * @brief Corrector function with convergence check (corrector)
   *
   * Calls update_solution(), then handles equation switching for coupled
   * problems (modifies cEq global).
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
