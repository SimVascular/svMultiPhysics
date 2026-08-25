// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef ACTIVE_STRESS_REGAZZONI_H
#define ACTIVE_STRESS_REGAZZONI_H

#include "ActiveStress.h"

#include <array>

/**
 * @brief Mean-field active stress model (implements the RDQ20-MF formulation).
 *
 * This class implements the mean-field RDQ20-MF sarcomere model of
 * cardiomyocyte force generation of Regazzoni, Dede', and Quarteroni (2020),
 * described in [1] and validated against the authors' reference implementation
 * [2]. The node-local state has 20 variables: 16 regulatory-unit (RU)
 * probabilities (entries 0-15) describing the tropomyosin/troponin
 * configuration of a triplet of neighbouring units, and 4 crossbridge (XB)
 * moments (entries 16-19). The RU probabilities are advanced with an explicit
 * forward-Euler substepping scheme — for every macro time step, a number of
 * smaller sub-steps are taken to update the RU states — and the XB moments with
 * one implicit-Euler step per time step; the active tension is then
 * reconstructed from the XB first moments.
 *
 * The returned scalar active tension is
 * @f[
 *   \Tact = a_\text{XB} \, (\mu_P^1 + \mu_N^1) \, \phi(SL)\;,
 * @f]
 * where @f$\mu_P^1@f$ and @f$\mu_N^1@f$ are the permissive and non-permissive
 * first XB moments (state entries 17 and 19), @f$\phi(SL)@f$ is the
 * single-overlap fraction of the sarcomere at sarcomere length @f$SL = SL_0 \,
 * \fiberstretch@f$ (with @f$\fiberstretch@f$ the fiber stretch), and
 * @f$a_\text{XB}@f$ is the tension upscaling factor. Because @f$\mu_P^1 +
 * \mu_N^1@f$ and @f$\phi(SL)@f$ are dimensionless, @f$a_\text{XB}@f$ sets the
 * units of the returned active tension.
 *
 * **References**:
 * 1. [Regazzoni, Dede', Quarteroni
 * (2020)](https://doi.org/10.1371/journal.pcbi.1008294)
 * 2. [F. Regazzoni, cardiac-activation reference
 * implementation](https://github.com/FrancescoRegazzoni/cardiac-activation)
 *
 * @note Although this model is governed by a system of ODEs, it inherits from
 * @c ActiveStress rather than @c ActiveStressODE because it requires a
 * customized time-stepping scheme to handle the stiffness of the model.
 *
 * @todo Force-strain-rate feedback requires a stabilization strategy for robust
 * use in coupled electromechanics. This will be addressed in a follow-up PR.
 */
class ActiveStressRegazzoni : public ActiveStress {
public:
  /// Model label, used for factory registration and XML selection.
  static inline const std::string label = "Regazzoni";

  /// @name State vector layout
  /// @{

  /// Number of regulatory-unit (RU) probability states (entries 0-15).
  static constexpr unsigned int n_ru_states = 16;

  /// Number of crossbridge (XB) moment states (entries 16-19).
  static constexpr unsigned int n_xb_states = 4;

  /// Total number of state variables.
  static constexpr unsigned int n_state_variables = n_ru_states + n_xb_states;

  /**
   * @brief Flat index of the RU probability state P(TL, TC, TR, CC).
   *
   * Each argument is 0 or 1 and denotes the state of, respectively, the left
   * tropomyosin unit, the central tropomyosin unit, the right tropomyosin unit
   * and the central troponin (calcium unbound/bound). The ordering matches the
   * reference implementation's serialization (TL outermost, CC innermost) and
   * spans [0, 15].
   */
  static constexpr unsigned int ru_index(unsigned int TL, unsigned int TC,
                                         unsigned int TR, unsigned int CC) {
    return 8 * TL + 4 * TC + 2 * TR + CC;
  }

  /// Flat index of the XB moment state @p i (in [0, 3]), spanning [16, 19].
  static constexpr unsigned int xb_index(unsigned int i) {
    return n_ru_states + i;
  }

  /// @}

  /**
   * @brief Model parameters class.
   *
   * Declares the parameters required by the model. All parameters are
   * marked as required, and omitting a parameter will cause a parse error.
   */
  class Parameters : public ActiveStressModelParameters {
  public:
    Parameters() : ActiveStressModelParameters(label) {
      constexpr bool required = true;

      // Reference values: Regazzoni 2020 human body-temperature calibration,
      // expressed consistently with the unit system used by this parameter set.
      add_parameter("Kbasic", 0.013, required);
      add_parameter("Koff", 0.1, required);
      add_parameter("Q", 2.0, required);
      add_parameter("mu", 10.0, required);
      add_parameter("gamma", 12.0, required);
      add_parameter("Kd0", 3.81e-4, required);
      add_parameter("alphaKd", -5.71e-4, required);
      add_parameter("SL0", 2.2, required);
      add_parameter("ru_substep", 2.5e-2, required);
      add_parameter("kd_reference_sarcomere_length", 2.15, required);

      add_parameter("r0", 0.13431, required);
      add_parameter("alpha", 25.184, required);
      add_parameter("mu0_fP", 0.032653, required);
      add_parameter("mu1_fP", 7.78e-4, required);

      add_parameter("LA", 1.25, required);
      add_parameter("LM", 1.65, required);
      add_parameter("LB", 0.18, required);
      add_parameter("a_XB", 22.894, required);

      add_parameter("Disable_force_strain_rate_feedback", false, !required);
    }
  };

  /**
   * @brief Constructor.
   */
  ActiveStressRegazzoni()
      : ActiveStress(/* n_state_variables = */ n_state_variables,
                     /* needs_fiber_stretch = */ true,
                     /* needs_fiber_stretch_rate = */ true) {}

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
   * Sets the state to (1, 0, ..., 0), i.e. all probability mass in the RU state
   * P(0, 0, 0, 0) and all crossbridge moments equal to zero.
   *
   * @param[out] state State vector for a single node, to be initialized by
   *   this function.
   */
  virtual void init_local(Vector<double> &state) const override;

  /**
   * @brief Advance in time for a single node.
   *
   * Advances the RU probabilities (entries 0-15) with the forward-Euler
   * substepping scheme and then the XB moments (entries 16-19) with one
   * implicit-Euler step, using the calcium, fiber stretch and fiber-stretch
   * rate at the node.
   */
  virtual void advance_time_step_local(const double t, const double dt,
                                       const double calcium,
                                       const double fiber_stretch,
                                       const double fiber_stretch_rate,
                                       Vector<double> &state) const override;

  /**
   * @brief Compute the scalar active tension for a single node.
   *
   * Evaluates @f$\Tact@f$ as defined in the class description, using
   * @p fiber_stretch to compute the sarcomere length
   * @f$SL = SL_0 \, \fiberstretch@f$. The returned value has the stress
   * units of @f$a_\text{XB}@f$.
   */
  virtual double
  compute_active_tension_local(const Vector<double> &state,
                               const double fiber_stretch) const override;

private:
  /// Array indexed over the four binary RU configuration variables (TL, TC, TR,
  /// CC).
  using RUArray =
      std::array<std::array<std::array<std::array<double, 2>, 2>, 2>, 2>;

  /// Array indexed over a pair of binary state variables.
  using BinaryPairArray = std::array<std::array<double, 2>, 2>;

  /// Array of the four crossbridge moment state variables.
  using XBArray = std::array<double, 4>;

  /// @name Regulatory-unit (RU) dynamics helpers
  /// @{

  /**
   * @brief Compute the central-tropomyosin transition rate for each local RU
   * configuration.
   *
   * Returns an @c RUArray where entry @c [TL][TC][TR][CC] is the rate at which
   * the central tropomyosin changes state for that configuration. Because the
   * rate depends on the neighbour states TL and TR, nearest-neighbour
   * cooperativity is retained through the tracked TL-TC-TR configuration.
   * These rates depend only on the model parameters, not on calcium or stretch.
   *
   * @return Central-tropomyosin transition rates, indexed @c [TL][TC][TR][CC].
   */
  RUArray ru_transition_rates_tropomyosin() const;

  /**
   * @brief Advance the 16 RU-state probabilities by one forward-Euler substep.
   *
   * Computes the probability fluxes caused by central-state transitions and the
   * effective boundary-neighbour transitions from the mean-field closure, then
   * updates @p state_RU in place.
   *
   * @param[in] dt Substep size [time].
   * @param[in] rates_T Central-tropomyosin transition rates,
   *   indexed @c rates_T[TL][TC][TR][CC].
   * @param[in] rates_C Troponin transition rates, indexed @c rates_C[CC][TC].
   * @param[in,out] state_RU The 16 RU-state probabilities,
   *   indexed @c state_RU[TL][TC][TR][CC].
   */
  void ru_forward_euler_substep(double dt, const RUArray &rates_T,
                                const BinaryPairArray &rates_C,
                                RUArray &state_RU) const;

  /**
   * @brief Advance the four crossbridge moments by one implicit-Euler step.
   *
   * Computes the permissivity and the effective permissive/non-permissive
   * transition rates from the updated RU probabilities, forms the 4x4 linear
   * system for the implicit update, and returns the updated moments.
   *
   * @param[in] dt Outer time step [time].
   * @param[in] velocity Shortening velocity @f$-\dot{SL}/SL_0@f$ [1/time].
   * @param[in] rates_T Central-tropomyosin transition rates,
   *   indexed @c rates_T[TL][TC][TR][CC].
   * @param[in] state_RU The updated 16 RU-state probabilities,
   *   indexed @c state_RU[TL][TC][TR][CC].
   * @param[in] state_XB The four crossbridge moments (input), ordered
   *   @f$[\mu_P^0, \mu_P^1, \mu_N^0, \mu_N^1]@f$.
   * @return Updated crossbridge moments, ordered
   *   @f$[\mu_P^0, \mu_P^1, \mu_N^0, \mu_N^1]@f$.
   */
  XBArray xb_implicit_update(double dt, double velocity, const RUArray &rates_T,
                             const RUArray &state_RU,
                             const XBArray &state_XB) const;

  /**
   * @brief Single-overlap fraction of the sarcomere at a given length.
   *
   * Returns the fraction @f$\phi(SL) \in [0, 1]@f$ of the sarcomere over which
   * thin and thick filaments overlap exactly once, a piecewise-linear function
   * of the sarcomere length built from the filament geometry (LA, LM, LB).
   *
   * @param[in] sarcomere_length Sarcomere length @f$SL@f$ [length].
   */
  double fraction_single_overlap(double sarcomere_length) const;

  /// @}

  /// @name RU model parameters
  /// @{

  double Kbasic; ///< Basic tropomyosin transition rate [1/time].
  double Koff;   ///< Troponin unbinding rate [1/time].
  double Q;      ///< Tropomyosin transition-rate asymmetry factor [-].
  double mu;     ///< Calcium-binding cooperativity factor [-].
  double gamma;  ///< Nearest-neighbour cooperativity factor [-].
  double Kd0; ///< Calcium dissociation constant at reference length [calcium].
  double alphaKd; ///< Length dependence of the dissociation constant
                  ///< [calcium/length].
  double SL0; ///< Reference sarcomere length [length]; maps stretch to length.
  double ru_substep; ///< RU forward-Euler substep size [time].

  /// Reference sarcomere length [length] used in the length-dependent
  /// dissociation constant (distinct from the parameter SL0).
  double kd_reference_sarcomere_length;

  double r0; ///< Combined attachment-detachment rate at zero velocity [1/time].
  double alpha;  ///< Coefficient of |v| in r(v) = r0 + alpha * |v| [-].
  double mu0_fP; ///< Permissive influx into the zeroth-moment crossbridge state
                 ///< [1/time].
  double mu1_fP; ///< Permissive influx into the first-moment crossbridge state
                 ///< [1/time].

  double LA; ///< Thin-filament (actin) length [length].
  double LM; ///< Thick-filament (myosin) length [length].
  double LB; ///< Length of the myosin bare zone [length].

  /// Tension upscaling factor [stress].
  ///
  /// Because the crossbridge moments and the overlap fraction are
  /// dimensionless, a_XB sets the stress unit of the returned active tension.
  /// It must be expressed in the same stress unit as the mechanical
  /// configuration.
  double a_XB;

  /// Controls force–strain-rate feedback in the XB update.
  /// - @c false (default): use @f$v = -\dot{\lambda}@f$ as shortening velocity.
  /// - @c true: set shortening velocity to zero, disabling the feedback.
  bool disable_force_strain_rate_feedback_ = false;

  /// @}
};

#endif
