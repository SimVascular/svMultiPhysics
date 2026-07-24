// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef ACTIVE_STRESS_NASH_PANFILOV_H
#define ACTIVE_STRESS_NASH_PANFILOV_H

#include "active_stress_ode.h"

/**
 * @brief Nash-Panfilov active stress model.
 *
 * This class implements the Nash-Panfilov active stress model [1], with the
 * modifications introduced by Goktepe and Kuhl [2].
 *
 * The model equations are the following:
 * @f[ \begin{aligned}
 * \dv{\Tact}{t} &= \varepsilon(\calcium)(
 *   \eta_\text{T} (\calcium - \calcium_\text{rest}) - \Tact)\;, \\
 * \varepsilon(\calcium) &=
 *   \varepsilon_0 + (\varepsilon_i - \varepsilon_0)
 *   \exp(-\exp(-\xi_T (\calcium - \calcium_\text{crit})))\;,
 * \end{aligned} @f]
 * where @f$\eta_\text{T}@f$, @f$\calcium_\text{rest}@f$,
 * @f$\calcium_\text{crit}@f$, @f$\varepsilon_0@f$, @f$\varepsilon_i@f$ and
 * @f$\xi_T@f$ are user-defined model parameters. The function
 * @f$\varepsilon(\calcium)@f$ is a sigmoidal-shaped calcium-dependent time
 * constant (see Figure 3 in [2] for more details).
 *
 * @note The sensitivity of the model to calcium is controlled by the paramter
 * @f$\eta_\text{T}@f$, which has the same units of active tension over calcium.
 * Therefore, if the ionic model providing the calcium is phenomenological (see
 * @ref IonicModel) and calcium is non-dimensional, this parameter may need to
 * be rescaled as well. Similar considerations apply to @f$\xi_T@f$.
 *
 * **References**:
 * 1. [Nash, Panfilov (2004)](https://doi.org/10.1016/j.pbiomolbio.2004.01.016)
 * 2. [Goktepe, Kuhl (2009)](https://doi.org/10.1007/s00466-009-0434-z)
 */
class NashPanfilov : public ActiveStressODE {
public:
  /// Model label.
  static inline const std::string label = "NashPanfilov";

  /// Model parameters class.
  class Parameters : public ActiveStressODE::Parameters {
  public:
    Parameters() : ActiveStressODE::Parameters(label) {
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
  /// @{

  /// Minimum time constant @f$\varepsilon_0@f$. The unit of measure for this
  /// parameter must be the inverse of the unit of measure for time.
  double epsilon_0;

  /// Maximum time constant @f$\varepsilon_i@f$. The unit of measure for this
  /// parameter must be the inverse of the unit of measure for time.
  double epsilon_i;

  /// Sigmoidal function steepness @f$\xi_T@f$. The unit of measure for this
  /// parameter must be the inverse of the unit of measure for calcium
  /// concentration.
  double xi_T;

  /// Resting calcium value. Active tension will increase if calcium is above
  /// this value.
  double calcium_rest;

  /// Critical calcium value, i.e. the threshold value for switching between
  /// minimum and maximum time constant.
  double calcium_crit;

  /// @f$\eta_T@f$. The unit of measure for this parameter must be the ratio of
  /// the unit for tension and the unit for calcium concentration.
  double eta_T;

  /// @}
};

#endif