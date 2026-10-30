// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "Array.h"
#include "SolutionStates.h"
#include "Vector.h"
#include "Simulation.h"

/**
 * @brief Integrator class that encapsulates the Newton iteration loop for time
 * integration.
 *
 * This class handles the nonlinear Newton iteration scheme for solving coupled
 * multi-physics equations in svMultiPhysics. It manages:
 * - Solution variables (Ag, Yg, Dg) at generalized-alpha time levels
 * - Newton iteration loop with convergence checking
 * - Linear system assembly and solve
 * - Boundary condition application
 *
 * ## Line search
 *
 * The nonlinear iterations of an equation can be globalized by a backtracking
 * line search, enabled by the `Line_search` XML element of the equation (see
 * @ref LineSearchParameters). Without it, every iteration applies the full
 * increment of the linear solve.
 *
 * At every nonlinear iteration of an equation with line search enabled, the
 * increment @f$\delta u@f$ of the linear solve is scaled by a step length
 * @f$\alpha@f$ before being applied, so that the solution is updated to
 * @f$u + \alpha\,\delta u@f$. The step lengths @f$\alpha = 1, \rho, \rho^2,
 * \dots@f$ are tried in turn, with @f$\rho@f$ the step reduction factor, and
 * the first one satisfying
 * @f[
 *   \| R(u + \alpha\,\delta u) \|_W \leq (1 - c\,\alpha)\,\| R(u) \|_W
 * @f]
 * is applied, with @f$c@f$ the minimum decrease. Once the step length reaches
 * the minimum step length, it is applied regardless of the residual decrease.
 *
 * The norm @f$\|\cdot\|_W@f$ weights every residual entry by the inverse
 * square root of the magnitude of the corresponding diagonal entry of the
 * tangent at @f$u@f$, as the Jacobi preconditioner of the linear solver does,
 * so that it is dimensionally homogeneous across the degrees of freedom of the
 * equation. Degrees of freedom constrained by a Dirichlet condition have zero
 * weight, since their residual is the reaction of the constraint. The weights
 * are computed once per iteration and kept fixed while step lengths are tried,
 * so that the norms of different step lengths are comparable (see
 * compute_residual_weights).
 *
 * The linear system is solved once per iteration, whatever the number of step
 * lengths tried. Measuring a step length costs one assembly of the linear
 * system at the trial solution, which re-advances implicitly coupled active
 * stress models as any other assembly does. The convergence test that follows
 * (see finalize_iteration) is the same with or without line search.
 *
 * Line search is not implemented for prestress simulations: the nodal
 * prestress is accumulated over the assemblies of a time step, which would
 * include those of the trial step lengths unless they were excluded.
 *
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
   * One line is written to the standard output and to the history file for every
   * Newton iteration, including the converged one.
   *
   * @param[in] save_results True if the results of this time step are written to a
   *   VTU file, which flags the converged iteration with an 's' in the standard
   *   output.
   *
   * @return True if all equations converged, false otherwise
   */
  bool step(bool save_results = false);

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
   * @brief Search for a step length that decreases the residual.
   *
   * Measures the increment of the linear solve at decreasing step lengths,
   * starting from one, until the residual assembled at the resulting solution
   * falls below \p reference_norm by the margin the settings of the equation
   * ask for, or until the shortest step length is reached, which is returned
   * without being measured. Measuring a step length costs one assembly of the
   * linear system, while the increment itself is reused, so the linear system
   * is solved once per nonlinear iteration however many step lengths are
   * tried.
   *
   * The solution and the increment are left as they were found, so that the
   * caller can apply the returned step length to them.
   *
   * @param[in,out] eq Equation whose increment is measured, and whose settings
   *   configure the search.
   * @param[in] weights Weight of every residual entry in the norm, as computed
   *   by compute_residual_weights at the solution the search starts from.
   * @param[in] reference_norm Norm of the residual at the solution the search
   *   starts from, measured with \p weights.
   *
   * @return Step length to apply to the increment.
   */
  double line_search_step_length(eqType &eq, const Array<double> &weights,
                                 const double reference_norm);

  /**
   * @brief Update residual and increment arrays for linear solver
   *
   * @param eq Reference to the equation being solved
   */
  void update_residual_arrays(eqType& eq);

  /**
   * @brief Assemble the linear system of an equation
   *
   * Evaluates the boundary conditions of the equation, forms the solution at
   * the intermediate generalized-alpha time levels, and assembles the residual
   * into com_mod.R and the tangent into com_mod.Val, including the
   * contributions of the coupled surfaces.
   *
   * @param[in,out] eq Equation whose linear system is assembled, and whose
   *   active stress models are re-advanced when their state is coupled
   *   implicitly.
   */
  void assemble_linear_system(eqType &eq);

  /**
   * @brief Compute the weight of every residual entry in the residual norm.
   *
   * The weight of a degree of freedom is the inverse square root of the
   * magnitude of the corresponding diagonal entry of the tangent, assembled
   * across processes, which is what the Jacobi preconditioner of the linear
   * solver scales the rows of the linear system by. The weight of a degree of
   * freedom constrained by a Dirichlet condition is zero, so that the residual
   * norm ignores the rows the linear solve eliminates, whose residual is a
   * reaction rather than an equation left to be satisfied.
   *
   * The resulting norm is the one the convergence of the nonlinear iterations
   * is tested with, and it is dimensionally homogeneous across the degrees of
   * freedom of an equation, which a norm of the residual entries themselves is
   * not.
   *
   * The weights are read from the tangent assembled last, which the assembly
   * of a further solution overwrites, so a line search computes them once at
   * the solution it starts from and reuses them for every step length it
   * tries. Keeping them fixed is what makes the norm a function of the
   * solution alone, and so makes the norms of two step lengths comparable.
   *
   * @param[out] weights Weight of every degree of freedom of every node, in
   *   the node numbering of the linear solver.
   */
  void compute_residual_weights(Array<double> &weights) const;

  /**
   * @brief Compute the weighted norm of the assembled residual
   *
   * Sums the squares of the weighted entries of com_mod.R belonging to the
   * nodes owned by each process and reduces them across processes, so that the
   * result does not depend on the partitioning. It is available as soon as the
   * residual has been assembled, before the linear system is solved.
   *
   * com_mod.R holds the residual of the equation assembled last, and its rows
   * are the degrees of freedom of that equation, numbered from zero.
   *
   * @param[in] weights Weight of every degree of freedom of every node, as
   *   computed by compute_residual_weights.
   *
   * @return Weighted Euclidean norm of the assembled residual.
   */
  double residual_norm(const Array<double> &weights) const;

  /**
   * @brief Compute the fiber stretch and fiber stretch rate at every node.
   *
   * The stretch is computed from the current displacement, so that calling this
   * function within the nonlinear iterations yields the stretch of the current
   * displacement iterate. Both vectors are left empty if no equation needs
   * them, and are filled with the neutral values (1 for the stretch, 0 for the
   * stretch rate) if they are needed but no equation solves for the
   * displacement.
   *
   * @param[out] fiber_stretch Fiber stretch at every node.
   * @param[out] fiber_stretch_rate Fiber stretch rate at every node.
   */
  void compute_fiber_stretch(Vector<double>& fiber_stretch,
                             Vector<double>& fiber_stretch_rate);

  /**
   * @brief Whether any domain uses an active stress model with implicit state
   * coupling, i.e. one whose state is updated within the nonlinear iterations.
   *
   * @return True if at least one such domain exists, false otherwise.
   */
  bool has_implicit_active_stress_state_coupling() const;

  /**
   * @brief Advance active stress models of an equation to the next time step.
   */
  void time_advance_active_stress(eqType &eq);

  /**
   * @brief Update the active stress models of an equation.
   *
   * @param[in,out] eq Equation whose domains carry the active stress models.
   * @param[in] fiber_stretch Fiber stretch at every node.
   * @param[in] fiber_stretch_rate Fiber stretch rate at every node.
   * @param[in] within_nonlinear_iterations True when called within the
   *   nonlinear iterations, in which case only the models with implicit
   *   coupling are advanced again, from the state stored at the beginning of
   *   the time step. False when called once per time step from the predictor,
   *   in which case all models are advanced from that state.
   */
  void update_active_stress(eqType &eq, const Vector<double> &fiber_stretch,
                            const Vector<double> &fiber_stretch_rate,
                            const bool within_nonlinear_iterations);

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
   * @brief Apply a multiple of the increment of the linear solve to the
   * solution
   *
   * Advances the solution of an equation at the n+1 time level by @c
   * step_length times the increment held in com_mod.R (and com_mod.Rd for the
   * velocity-based formulation), and updates the quantities that are derived
   * from the resulting solution: the Taylor-Hood pressure at edge nodes, the
   * copy of the solution onto the solid subdomain of an FSI equation, the
   * ionic state of an electrophysiology equation, and the wall filter of a CMM
   * equation.
   *
   * The solution change is proportional to step_length, so calling this
   * function from the same starting solution with different step lengths traces
   * the segment between that solution and the one the full increment produces.
   * Every quantity derived from the solution is recomputed, so any step length
   * leaves a consistent state.
   *
   * Modifies:
   * \code {.cpp}
   *   com_mod.Ad
   *   solutions_.current.A
   *   solutions_.current.D
   *   solutions_.current.Y
   *   cep_mod.Xion
   * \endcode
   *
   * @param[in] eq Equation whose solution is advanced, which selects the
   *   degrees of freedom the increment is applied to and the coefficients of
   *   the time integration scheme.
   * @param[in] step_length Step length multiplying the increment. One applies
   * the full increment of the linear solve.
   */
  void apply_increment(const eqType &eq, const double step_length);

  /**
   * @brief Close a nonlinear iteration of the current equation
   *
   * Normalizes the nodal prestress accumulated during the assembly, tests the
   * residual norm of the current equation for convergence, and selects the
   * equation to be solved by the next iteration.
   *
   * The nodal quantities this function reduces are accumulated once per
   * assembly and its convergence bookkeeping counts one iteration, so it is
   * called once per iteration, after the solution has been advanced by
   * apply_increment.
   *
   * Modifies:
   * \code {.cpp}
   *   com_mod.pSa
   *   com_mod.pSn
   *   solutions_.current.A
   *   solutions_.current.D
   *   solutions_.current.Y
   *
   *   com_mod.cEq
   *   eq.FSILS.RI.iNorm
   *   eq.iNorm
   *   eq.ok
   *   eq.pNorm
   * \endcode
   */
  void finalize_iteration();

  /**
   * @brief Pressure correction for Taylor-Hood elements (corrector_taylor_hood)
   *
   * Interpolates pressure at edge nodes using reduced basis applied
   * on element vertices for Taylor-Hood type elements.
   */
  void corrector_taylor_hood();
};

#endif // INTEGRATOR_H
