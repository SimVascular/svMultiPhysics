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
  /// Model label.
  static inline const std::string label = "Uniform";

  /// Model parameters class.
  class Parameters : public ActiveStressModelParameters {
  public:
    Parameters() : ActiveStressModelParameters(label) {
      constexpr bool required = true;

      add_parameter("Value", 0.0, required);
    }
  };

  /**
   * @brief Constructor.
   */
  UniformActiveStress() : ActiveStress(/* n_states = */ 0) {}

  /**
   * @brief Construct an instance of model parameters.
   */
  virtual std::unique_ptr<ActiveStressModelParameters>
  get_parameters() const override {
    return std::make_unique<Parameters>();
  }

  /**
   * @brief Read model parameters from a parameter object.
   */
  virtual void
  read_parameters(const ActiveStressModelParameters &params) override;

  /**
   * @brief Distribute model parameters to all parallel processes.
   */
  virtual void distribute_parameters(const CmMod &cm_mod,
                                     const cmType &cm) override;

  /**
   * @brief Evaluate the active stress at a given point.
   */
  virtual double operator()(const int idx) const override;

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
  double value;
};

#endif