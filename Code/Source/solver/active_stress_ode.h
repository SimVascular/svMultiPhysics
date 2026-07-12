// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef ACTIVE_STRESS_ODE_H
#define ACTIVE_STRESS_ODE_H

#include "active_stress.h"

/**
 * @brief Abstract ODE-based active stress model.
 *
 * This class provides an interface for defining active stress models for which
 * the evolution of the state vector @f$\astressstate@f$ is governed by a system
 * of ordinary differential equations (ODEs):
 * @f[ \begin{aligned}
 *   \dv{\astressstate}{t} &=
 *     \mathbf{F}_\text{AS}(t, \astressstate, \calcium, \fiberstretch,
 *                          \fiberstretchrate)\;, \\
 *   \Tact &= \Tact(\astressstate)\;.
 * \end{aligned} @f]
 *
 * ### Numerical methods
 *
 * The ODE system is advanced with one of the methods described in
 * @ref ODESolver. After that, the active tension is computed for every node
 * @f$i@f$ as:
 * @f[
 *   {\Tact}_{i}^{n+1} = \Tact(\astressstate_i^{n+1})\;.
 * @f]
 *
 * ### Implementing derived models
 *
 * To implement a new ODE-based active stress model, the following steps need to
 * be taken:
 *
 * 1. Create a new class derived from @ref ActiveStressODE.
 * 2. Override the method init_local (from the base class @ref ActiveStress) to
 *    define the initial condition for the state vector at a single node.
 * 3. Override the method @ref getf to define the right-hand side function
 *    @f$\mathbf{F}_\text{AS}@f$ of the ODE system.
 * 4. Create a new class derived from @ref ActiveStressODE::Parameters to store
 *    the parameters specific to the new active stress model.
 * 5. Override the methods @ref get_parameters, @ref read_parameters and
 *    @ref distribute_parameters to manage the parameters of the new active
 *    stress model. Both read_parameters and distribute_parameters need to call
 *    the corresponding method of ActiveStressODE.
 * 6. Register the new class into the active stress model factory by using the
 *    macro @ref REGISTER_ACTIVE_STRESS_MODEL. The macro should be called in a
 *    `.cpp` file, not in a header file.
 */
class ActiveStressODE : public ActiveStress {
public:
  /**
   * @brief Enumeration of supported ODE solvers.
   *
   * In the documentation of the individual methods below, the subscript @f$i@f$
   * is used to denote the value of variables at the @f$i@f$-th node, and the
   * superscript @f$n@f$ is used to denote the time step.
   */
  enum class ODESolver {
    /**
     * @brief Forward Euler.
     *
     * The state vector is updated as follows:
     * @f[
     *  \astressstate_i^{n+1} = \astressstate_i^n
     *    + \Delta t \mathbf{F}_\text{AS}(t^n, \astressstate_i^n, \calcium_i^n,
     *                                    \fiberstretch_i^n,
     *                                    \fiberstretchrate_i^n)\;.
     * @f]
     */
    ForwardEuler
  };

  /**
   * @brief Model parameters class.
   *
   * Declares parameters that are common to all ODE-based active stress models.
   * Classes derived from ActiveStressODE one should declare their own
   * Parameters class, deriving from ActiveStressODE::parameters.
   */
  class Parameters : public ActiveStressModelParameters {
  public:
    Parameters(const std::string &label) : ActiveStressModelParameters(label) {
      constexpr bool required = true;

      add_parameter("ODE_solver", std::string("FE"), required);
    }
  };

  /**
   * @brief Constructor.
   *
   * @param n_states Number of state variables for this model.
   */
  ActiveStressODE(const unsigned int n_states) : ActiveStress(n_states) {}

protected:
  /**
   * @brief Read model parameters from a parameter object.
   */
  virtual void read_model_specific_parameters(
      const ActiveStressModelParameters &params) override;

  /**
   * @brief Distribute model parameters to all parallel processes.
   */
  virtual void distribute_model_specific_parameters(const CmMod &cm_mod,
                                                    const cmType &cm) override;

  /**
   * @brief Advance in time for a single node.
   *
   * Solves one forward Euler time step for the ODE system.
   *
   * @param[in] t Current time (i.e. the time instant being advanced to).
   * @param[in] dt Time step size.
   * @param[in] calcium Calcium concentration at the current node.
   * @param[in] fiber_stretch Fiber stretch at the current node.
   * @param[in] fiber_stretch_rate Fiber stretch rate at the current node.
   * @param[in,out] state State vector for a single node, to be updated by
   *   this function.
   *
   * @todo[michelebucelli] It might be necessary or useful to implement other
   *   timestepping schemes, e.g. Runge-Kutta. In that case, we might want to
   *   expand the interface to support implicit time stepping too, e.g. by
   *   adding a method to evaluate the Jacobian matrix of the system, as in
   *   @ref IonicModel.
   */
  virtual void advance_time_step_local(const double t, const double dt,
                                       const double calcium,
                                       const double fiber_stretch,
                                       const double fiber_stretch_rate,
                                       Vector<double> &state) const override;

  /**
   * @brief Compute the rate of change in the state variables.
   *
   * @param[in] t Current time (i.e. the time instant being advanced to).
   * @param[in] dt Time step size.
   * @param[in] calcium Calcium concentration at the current node.
   * @param[in] fiber_stretch Fiber stretch at the current node.
   * @param[in] fiber_stretch_rate Fiber stretch rate at the current node.
   *
   * @return A vector containing the rate of change for each state variable.
   */
  virtual Vector<double> getf(const double t, const Vector<double> &state,
                              const double calcium, const double fiber_stretch,
                              const double fiber_stretch_rate) const = 0;

  /// ODE solver.
  ODESolver ode_solver;
};

#endif