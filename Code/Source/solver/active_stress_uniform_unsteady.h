// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef ACTIVE_STRESS_UNIFORM_UNSTEADY_H
#define ACTIVE_STRESS_UNIFORM_UNSTEADY_H

#include "FourierInterpolation.h"

#include "active_stress.h"

/**
 * @brief Uniform and time dependent active stress model.
 *
 * Defines an active tension that is constant in space but time dependent, i.e.
 * @f[
 *   \Tact(t, \calcium, \fiberstretch, \fiberstretchrate, \astressstate) =
 *     g(t)\;,
 * @f]
 * where @f$g(t)@f$ is a user-defined function of time.
 */
class UniformUnsteadyActiveStress : public ActiveStress {
public:
  /// Model label.
  static inline const std::string label = "UniformUnsteady";

  /// Model parameters class.
  class Parameters : public ActiveStressModelParameters {
  public:
    Parameters() : ActiveStressModelParameters(label) {
      constexpr bool required = true;

      add_parameter("Ramp", false, required);
      add_parameter("Temporal_values_file_path", std::string(""), required);
    }
  };

  /**
   * @brief Constructor.
   */
  UniformUnsteadyActiveStress() : ActiveStress(/* n_states = */ 0) {}

  /**
   * @brief Construct an instance of model parameters.
   */
  virtual std::unique_ptr<ActiveStressModelParameters>
  get_parameters() const override {
    return std::make_unique<Parameters>();
  }

  /**
   * @brief Initialization.
   *
   * Calls the parent class initialization method, and reads the Fourier
   * coefficient from file.
   */
  virtual void init(const unsigned int tnNo) override;

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
  compute_active_tension_local(const Vector<double> &state) const override;

  /// Toggle between ramp or Fourier transform.
  bool ramp;

  /// Name of the file containing the temporal values.
  std::string temporal_values_file_path;

  /// Fourier interpolation of the time dependent data.
  FourierInterpolation fourier_interpolation;
};

#endif