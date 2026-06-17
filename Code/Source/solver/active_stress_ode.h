// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef ACTIVE_STRESS_ODE_H
#define ACTIVE_STRESS_ODE_H

#include "active_stress.h"

/**
 * @brief Abstract ODE-based active stress model.
 *
 * @todo[michelebucelli] Detailed documentation.
 */
class ActiveStressODE : public ActiveStress {
public:
  /**
   * @brief Constructor.
   *
   * @param n_states Number of state variables for this model.
   */
  ActiveStressODE(const unsigned int n_states) : ActiveStress(n_states) {}

protected:
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
                                       Vector<double> &state) const override;

  /**
   * @brief Compute the rate of change in the state variables.
   */
  virtual Vector<double> getf(const double t, const Vector<double> &state,
                              const double calcium, const double fiber_stretch,
                              const double fiber_stretch_rate) const = 0;
};

#endif