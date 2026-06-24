// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef ACTIVE_STRESS_UNIFORM_STEADY_H
#define ACTIVE_STRESS_UNIFORM_STEADY_H

#include "active_stress.h"

/**
 * @brief Uniform and steady active stress model.
 *
 * Defines an active tension that is constant in space and time, i.e.
 * @f[
 *   \Tact(t, \calcium, \fiberstretch, \fiberstretchrate, \astressstate) = g\;,
 * @f]
 * where @f$g@f$ is a user-defined constant value.
 */
class UniformSteadyActiveStress : public ActiveStress {
public:
  /// Model label.
  static inline const std::string label = "UniformSteady";

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
  UniformSteadyActiveStress() : ActiveStress(/* n_states = */ 0) {}

  /**
   * @brief Construct an instance of model parameters.
   */
  virtual std::unique_ptr<ActiveStressModelParameters>
  get_parameters() const override {
    return std::make_unique<Parameters>();
  }

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

  /**
   * @brief Compute the active tension for a single node.
   */
  virtual double
  compute_active_tension_local(const Vector<double> &state) const override {
    return value;
  }

  /// Active tension value.
  double value;
};

#endif