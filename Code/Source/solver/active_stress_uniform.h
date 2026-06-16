// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef ACTIVE_STRESS_UNIFORM_H
#define ACTIVE_STRESS_UNIFORM_H

#include "active_stress.h"

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
   */
  UniformActiveStress() : ActiveStress(/* n_states = */ 0), g(4.9875) {}

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