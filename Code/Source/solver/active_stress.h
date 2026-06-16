// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef ACTIVE_STRESS_H
#define ACTIVE_STRESS_H

#include "Array.h"
#include "Vector.h"
#include "factory.h"

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
 * @brief Alias for the active stress model factory.
 *
 * See the documentation for @ref Factory for more details on how this works.
 */
using ActiveStressFactory = Factory<ActiveStress>;

/**
 * @brief Macro to register an active stress model in the factory.
 */
#define REGISTER_ACTIVE_STRESS_MODEL(name, type)                               \
  REGISTER_IN_FACTORY(ActiveStress, type, name)

#endif