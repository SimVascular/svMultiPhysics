// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef ACTIVE_STRESS_H
#define ACTIVE_STRESS_H

#include "Array.h"
#include "Vector.h"

/**
 * @brief Abstract active stress class.
 *
 * @todo Detailed documentation.
 */
class ActiveStress {
public:
  /**
   * @brief Constructor.
   *
   * @param n_states_ Number of state variables for this model.
   */
  ActiveStress(const unsigned int n_states_) : n_states(n_states_) {}

  /**
   * @brief Virtual destructor.
   */
  virtual ~ActiveStress() = default;

  /**
   * @brief Evaluate the active stress at a given point.
   *
   * @param[in] idx Index of the degree of freedom for which the active stress
   *   must be returned.
   *
   * @return Active stress value at the given point.
   * @todo[michelebucelli] Document the unit of measure.
   */
  virtual double operator()(const int idx) const = 0;

  /**
   * @brief Initialize the model.
   *
   * Allocates the internal state vector and initializes it with the model's
   * initial conditions.
   *
   * @param[in] tnNo Total number of mesh nodes for the current rank.
   *
   * @todo[michelebucelli] Double check that tnNo is the right number to pass
   *   here (as opposed to e.g. the number of nodes in the domain where this
   *   model is defined).
   */
  void init(const unsigned int tnNo);

  /**
   * @brief Advance in time.
   *
   * @param[in] t Current time (i.e. the time instant being advanced to).
   * @param[in] dt Time step size.
   */
  virtual void advance_time_step(const double t, const double dt);

  /// Number of state variables for this model.
  const unsigned int n_states;

protected:
  /**
   * @brief Initialize the state vector for a single node.
   *
   * @param[out] state State vector for a single node, to be initialized by
   *   this function.ß
   */
  virtual void init_local(Vector<double> &state) const = 0;

  /**
   * @brief Advance in time for a single node.
   *
   * @param[in] t Current time (i.e. the time instant being advanced to).
   * @param[in] dt Time step size.
   * @param[in] calcium Calcium concentration at the current node.
   * @param[in] fiber_stretch Fiber stretch at the current node.
   * @param[in] fiber_stretch_rate Fiber stretch rate at the current node.
   * @param[in,out] state State vector for a single node, to be updated by
   *   this function.
   */
  virtual void advance_time_step_local(const double t, const double dt,
                                       const double calcium,
                                       const double fiber_stretch,
                                       const double fiber_stretch_rate,
                                       Vector<double> &state) const = 0;

  /// State variables for the model.
  Array<double> states;
};

/**
 * @brief Uniform active stress model.
 *
 * Returns an active tension value constant in space and in time, provided at
 * construction.
 */
class UniformActiveStress : public ActiveStress {
public:
  /**
   * @brief Constructor.
   *
   * @param g_ Active tension value.
   */
  UniformActiveStress(const double g_)
      : ActiveStress(/* n_states = */ 0), g(g_) {}

  /**
   * @brief Evaluate the active stress at a given point.
   */
  virtual double operator()(const int idx) const override { return g; }

protected:
  /**
   * @brief Initialize the state vector for a single node.
   *
   * This model has no states, so this function does nothing.
   */
  virtual void init_local(Vector<double> &state) const override {}

  /**
   * @brief Advance in time for a single node.
   *
   * This model has no states, so this function does nothing.
   */
  virtual void advance_time_step_local(const double t, const double dt,
                                       const double calcium,
                                       const double fiber_stretch,
                                       const double fiber_stretch_rate,
                                       Vector<double> &state) const override {}

  /// Active tension value.
  /// @todo[michelebucelli] Document unit of measure.
  double g;
};

#endif