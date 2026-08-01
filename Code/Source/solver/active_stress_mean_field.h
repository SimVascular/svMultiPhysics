// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef ACTIVE_STRESS_MEAN_FIELD_H
#define ACTIVE_STRESS_MEAN_FIELD_H

#include "active_stress.h"

/**
 * @brief Mean-field active stress model (implements the RDQ20-MF formulation).
 *
 * This class implements the mean-field RDQ20-MF sarcomere model of cardiomyocyte
 * force generation of Regazzoni, Dedè, and Quarteroni (2020), described in [1]
 * and validated against the authors' reference implementation [2]. The node-local state has 20 variables: 16
 * regulatory-unit (RU) probabilities (entries 0-15) describing the
 * tropomyosin/troponin configuration of a triplet of neighbouring units, and 4
 * crossbridge (XB) moments (entries 16-19). The RU probabilities are advanced
 * with an explicit forward-Euler substepping scheme and the XB moments with one
 * implicit-Euler step per time step; the active tension is then reconstructed
 * from the XB first moments.
 *
 * The active tension is
 * @f[
 *   T_\text{act} = a_\text{XB} \, (\mu_P^1 + \mu_N^1) \, \phi(SL)\;,
 * @f]
 * where @f$\mu_P^1@f$ and @f$\mu_N^1@f$ are the permissive and non-permissive
 * first XB moments (state entries 17 and 19), @f$\phi(SL)@f$ is the single-overlap
 * fraction of the sarcomere at sarcomere length @f$SL = SL_0 \, \lambda@f$ (with
 * @f$\lambda@f$ the fiber stretch), and @f$a_\text{XB}@f$ is the tension
 * upscaling factor. Because @f$\mu_P^1 + \mu_N^1@f$ and @f$\phi(SL)@f$ are
 * dimensionless, @f$a_\text{XB}@f$ sets the units of the returned active tension;
 * it is stored in the stress units of the simulation, so no separate output
 * conversion is applied.
 *
 * The model is calibrated in a fixed unit system (calcium in [uM], time in [s],
 * length in [um]); the svMultiPhysics inputs are converted to these units at the
 * interface (see the conversion members below).
 *
 * **References**:
 * 1. [Regazzoni, Dede', Quarteroni (2020)](https://doi.org/10.1371/journal.pcbi.1008294)
 * 2. [F. Regazzoni, cardiac-activation reference implementation](https://github.com/FrancescoRegazzoni/cardiac-activation)
 */
class MeanFieldActiveStress : public ActiveStress {
public:
  /// Model label, used for factory registration and XML selection.
  static inline const std::string label = "MeanFieldActiveStress";

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
  static constexpr int ru_index(int TL, int TC, int TR, int CC) {
    return 8 * TL + 4 * TC + 2 * TR + CC;
  }

  /// Flat index of the XB moment state @p i (in [0, 3]), spanning [16, 19].
  static constexpr int xb_index(int i) {
    return static_cast<int>(n_ru_states) + i;
  }

  /// @}

  /**
   * @brief Model parameters class.
   *
   * Declares the regulatory-unit (RU), crossbridge (XB), geometry, and tension
   * parameters. The values registered below correspond to the published human
   * body-temperature calibration of the reference implementation, expressed in
   * the units documented for each member. The registered value of @c a_XB,
   * 22.894, expresses the reference calibration in MPa; the value supplied in
   * solver.xml must instead use the stress unit of the simulation's mechanical
   * configuration.
   *
   * All parameters are required. A complete @c MeanFieldActiveStress parameter block
   * containing every parameter must be provided in solver.xml. Any value may be changed to
   * use a different calibration, but omitting a parameter causes a parse error;
   * the registered reference value is not used as an automatic default.
   */
  class Parameters : public ActiveStressModelParameters {
  public:
    Parameters() : ActiveStressModelParameters(label) {
      constexpr bool required = true;

      add_parameter("Kbasic", 13.0, required);
      add_parameter("Koff", 100.0, required);
      add_parameter("Q", 2.0, required);
      add_parameter("mu", 10.0, required);
      add_parameter("gamma", 12.0, required);
      add_parameter("Kd0", 0.381, required);
      add_parameter("alphaKd", -0.571, required);
      add_parameter("SL0", 2.2, required);

      add_parameter("r0", 134.31, required);
      add_parameter("alpha", 25.184, required);
      add_parameter("mu0_fP", 32.653, required);
      add_parameter("mu1_fP", 0.778, required);

      add_parameter("LA", 1.25, required);
      add_parameter("LM", 1.65, required);
      add_parameter("LB", 0.18, required);
      add_parameter("a_XB", 22.894, required);
    }
  };

  /**
   * @brief Constructor.
   */
  MeanFieldActiveStress() : ActiveStress(n_state_variables) {}

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
   * implicit-Euler step, using the calcium, fiber stretch and fiber-stretch rate
   * at the node.
   */
  virtual void advance_time_step_local(const double t, const double dt,
                                       const double calcium,
                                       const double fiber_stretch,
                                       const double fiber_stretch_rate,
                                       Vector<double> &state) const override;

  /**
   * @brief Compute the PK2 active-tension coefficient for a single node.
   *
   * Computes the paper tension
   * @f$\widetilde{T}_\text{act} = a_\text{XB} (\mu_P^1 + \mu_N^1) \phi(SL)@f$
   * from the XB first moments (state entries 17 and 19) and the single-overlap
   * fraction at @f$SL = SL_0 \lambda_f@f$, then returns
   * @f$\widetilde{T}_\text{act} / \lambda_f@f$.
   *
   * The division by @f$\lambda_f@f$ is necessary because svMultiPhysics
   * assembles the active stress as
   * @f$\mathbf{S}_\text{act} = T_\text{code}\,\mathbf{f}_0\otimes\mathbf{f}_0@f$
   * (second-PK), so the resulting first-PK active stress satisfies
   * @f$\|\mathbf{P}_\text{act}\| = \widetilde{T}_\text{act}@f$, matching the
   * paper's definition of active tension.
   *
   * The returned value has the stress units of @f$a_\text{XB}@f$.
   *
   * @note The slab regression test uses a (0.7, 0.2, 0.1) fiber/sheet/normal
   *   directional distribution, which is an svMultiPhysics extension of the
   *   paper's fiber-only active stress formulation.
   *
   * @throws svmp::FE::InvalidArgumentException if @p fiber_stretch is zero
   *   or near zero (degenerate/collapsed element).
   */
  virtual double
  compute_active_tension_local(const Vector<double> &state,
                               const double fiber_stretch) const override;

private:
  /// @name Regulatory-unit (RU) dynamics helpers
  /// @{

  /**
   * @brief Compute the central-tropomyosin transition rate for each local RU
   * configuration.
   *
   * Fills @p rates_T, indexed as @c rates_T[TL][TC][TR][CC], where TL, TC and TR
   * are the binary left, central and right tropomyosin states and CC is the
   * central troponin calcium-binding state. Each entry is the rate at which the
   * central tropomyosin changes state for that configuration. Because the rate
   * depends on the neighbour states TL and TR, nearest-neighbour cooperativity
   * is retained through the tracked TL-TC-TR configuration. These rates depend
   * only on the model parameters, not on calcium or stretch.
   */
  void ru_transition_rates_tropomyosin(double (&rates_T)[2][2][2][2]) const;

  /**
   * @brief Advance the 16 RU-state probabilities by one forward-Euler substep.
   *
   * Computes the probability fluxes caused by central-state transitions and the
   * effective boundary-neighbour transitions from the mean-field closure, then
   * updates @p state_RU in place.
   *
   * @param[in] dt Substep size [s].
   * @param[in] rates_T Central-tropomyosin transition rates,
   *   indexed @c rates_T[TL][TC][TR][CC].
   * @param[in] rates_C Troponin transition rates, indexed @c rates_C[CC][TC].
   * @param[in,out] state_RU The 16 RU-state probabilities,
   *   indexed @c state_RU[TL][TC][TR][CC].
   */
  void ru_forward_euler_substep(double dt,
                                const double (&rates_T)[2][2][2][2],
                                const double (&rates_C)[2][2],
                                double (&state_RU)[2][2][2][2]) const;

  /**
   * @brief Advance the four crossbridge moments by one implicit-Euler step.
   *
   * Computes the permissivity and the effective permissive/non-permissive
   * transition rates from the updated RU probabilities, forms the 4x4 linear
   * system for the implicit update and solves it in place for @p state_XB.
   *
   * @param[in] dt Outer time step [s].
   * @param[in] velocity Shortening velocity @f$-\dot{SL}/SL_0@f$ [s^-1].
   * @param[in] rates_T Central-tropomyosin transition rates,
   *   indexed @c rates_T[TL][TC][TR][CC].
   * @param[in] state_RU The updated 16 RU-state probabilities,
   *   indexed @c state_RU[TL][TC][TR][CC].
   * @param[in,out] state_XB The four crossbridge moments, ordered
   *   @f$[\mu_P^0, \mu_P^1, \mu_N^0, \mu_N^1]@f$.
   */
  void xb_implicit_update(double dt, double velocity,
                          const double (&rates_T)[2][2][2][2],
                          const double (&state_RU)[2][2][2][2],
                          double (&state_XB)[4]) const;

  /**
   * @brief Single-overlap fraction of the sarcomere at a given length.
   *
   * Returns the fraction @f$\phi(SL) \in [0, 1]@f$ of the sarcomere over which
   * thin and thick filaments overlap exactly once, a piecewise-linear function
   * of the sarcomere length built from the filament geometry (LA, LM, LB).
   *
   * @param[in] sarcomere_length Sarcomere length @f$SL@f$ [um].
   */
  double fraction_single_overlap(double sarcomere_length) const;

  /// @}

  /// @name Interface unit conversions (svMultiPhysics EM units to reference units)
  /// @{

  /// Calcium conversion, millimolar [mM] to micromolar [uM].
  static constexpr double calcium_mM_to_microM = 1.0e3;

  /// Time conversion, milliseconds [ms] to seconds [s].
  static constexpr double time_ms_to_s = 1.0e-3;

  /// @}

  /// RU forward-Euler substep size [s].
  static constexpr double ru_substep = 2.5e-5;

  /// Fixed reference sarcomere length [um] in the length-dependent dissociation
  /// constant (distinct from the parameter SL0).
  static constexpr double kd_reference_sarcomere_length = 2.15;

  /// @name RU model parameters (reference units)
  /// @{

  double Kbasic;  ///< Basic tropomyosin transition rate [s^-1].
  double Koff;    ///< Troponin unbinding rate [s^-1].
  double Q;       ///< Tropomyosin transition-rate asymmetry factor [-].
  double mu;      ///< Calcium-binding cooperativity factor [-].
  double gamma;   ///< Nearest-neighbour cooperativity factor [-].
  double Kd0;     ///< Calcium dissociation constant at reference length [uM].
  double alphaKd; ///< Length dependence of the dissociation constant [uM/um].
  double SL0;     ///< Reference sarcomere length [um]; maps stretch to length.

  double r0;     ///< Combined attachment-detachment rate at zero velocity [s^-1].
  double alpha;  ///< Coefficient of |v| in r(v) = r0 + alpha * |v| [-].
  double mu0_fP; ///< Permissive influx into the zeroth-moment crossbridge state [s^-1].
  double mu1_fP; ///< Permissive influx into the first-moment crossbridge state [s^-1].

  double LA; ///< Thin-filament (actin) length [um].
  double LM; ///< Thick-filament (myosin) length [um].
  double LB; ///< Length of the myosin bare zone [um].

  /// Tension upscaling factor.
  ///
  /// Because the crossbridge moments and the overlap fraction are dimensionless,
  /// a_XB is the only quantity carrying stress units, so the returned active
  /// tension has the same stress unit as a_XB and no stress-unit conversion is
  /// performed. a_XB must therefore be expressed in the same stress unit as the
  /// mechanical configuration. The reference calibration value 22.894 is
  /// expressed in MPa, consistent with the coupled electromechanics slab case;
  /// the equivalent values are 22.894e3 in kPa and 22.894e6 in Pa. Provide the
  /// value matching the case's stress unit.
  double a_XB;

  /// @}
};

#endif
