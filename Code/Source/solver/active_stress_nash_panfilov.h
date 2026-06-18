// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef ACTIVE_STRESS_NASH_PANFILOV_H
#define ACTIVE_STRESS_NASH_PANFILOV_H

#include "active_stress_ode.h"

/**
 * @brief Nash-Panfilov active stress model.
 *
 * @todo[michelebucelli] Detailed documentation, with reference.
 */
class NashPanfilov : public ActiveStressODE {
public:
  /// Model label.
  static inline const std::string label = "NashPanfilov";

  /// Model parameters class.
  class Parameters : public ActiveStressModelParameters {
  public:
    Parameters() : ActiveStressModelParameters(label) {
      constexpr bool required = true;

      add_parameter("epsilon_0", 1.0, required);
      add_parameter("epsilon_i", 1.0, required);
      add_parameter("xi_T", 1.0, required);
      add_parameter("calcium_rest", 1.0, required);
      add_parameter("calcium_crit", 1.0, required);
      add_parameter("eta_T", 1.0, required);
    }
  };

  /**
   * @brief Constructor.
   */
  NashPanfilov() : ActiveStressODE(1) {}

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

protected:
  /**
   * @brief Initialize the state vector for a single node.
   *
   * @param[out] state State vector for a single node, to be initialized by
   *   this function.
   */
  virtual void init_local(Vector<double> &state) const override;

  /**
   * @brief Compute the rate of change in the state variables.
   */
  virtual Vector<double> getf(const double t, const Vector<double> &state,
                              const double calcium, const double fiber_stretch,
                              const double fiber_stretch_rate) const override;

  /**
   * @brief Compute the active tension for a single node.
   */
  virtual double
  compute_active_tension_local(const Vector<double> &state) const override;

  /// @name Model parameters.
  /// @todo[michelebucelli] Document meaning and units of measure for all
  /// parameters.
  /// @{

  /// @f$\varepsilon_0@f$.
  double epsilon_0;

  /// @f$\varepsilon_i@f$.
  double epsilon_i;

  /// @f$\xi_T@f$.
  double xi_T;

  /// Resting calcium value.
  double calcium_rest;

  /// Critical calcium value.
  double calcium_crit;

  /// @f$\eta_T@f$.
  double eta_T;

  /// @}
};

#endif