// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef IONIC_TTP_H
#define IONIC_TTP_H

#include "ionic_model.h"

#include "Parameters.h"

#include "Vector.h"
#include "utils.h"

/**
 * @brief Ten Tusscher-Panfilov ionic model.
 *
 * This class implements the 2006 version of the model, described in the second
 * reference below.
 *
 * **References**:
 * 1. Ten Tusscher, Noble, Noble, Panfilov. A model for human ventricular
 * tissue. American Journal of Physiology - Heart and Circulatory Physiology
 * (2004).
 * 2. Ten Tusscher, Panfilov. Alternans and spiral breakup in a human
 * ventricular tissue model. American Journal of Physiology - Heart and
 * Circulatory Physiology (2006).
 *
 * Model parameters are from reference 2 above. Default parameters are for
 * epicardium state (source: https://models.cellml.org/e/80d)
 */
class TTP : public IonicModel {
public:
  /// Model label.
  static inline const std::string label = "TTP";

  /// State variables.
  static inline const InitialStates initial_X = {
      {"V", -85.23}, {"K_i", 136.89}, {"Na_i", 8.6040}, {"Ca_i", 1.26e-4},
      {"Ca_ss", 3.6e-4}, {"Ca_sr", 3.64}, {"R_bar", 0.9073}};

  /// Gating variables.
  static inline const InitialStates initial_Xg = {
      {"x_r1_rectifier", 6.21e-3}, {"x_r2_rectifier", 0.4712},
      {"x_s_rectifier", 9.5e-3}, {"m_fast_Na", 1.72e-3}, {"h_fast_Na", 0.7444},
      {"j_fast_Na", 0.7045}, {"d_slow_in", 3.373e-5}, {"f_slow_in", 0.7888},
      {"f2_slow_in", 0.9755}, {"fcass_slow_in", 0.9953}, {"s_out", 0.999998},
      {"r_out", 2.42e-8}};

  /// Index of the intracellular calcium concentration (Ca_i) in the state
  /// vector.
  static constexpr unsigned int calcium_index = 3;

  /// Model parameters class.
  class Parameters : public IonicModelParameters {
  public:
    Parameters() : IonicModelParameters(label, initial_X, initial_Xg) {
      constexpr bool required = true;

      add_parameter("Rc", 8314.472, required);
      add_parameter("Tc", 310.0, required);
      add_parameter("Fc", 96485.3415, required);
      add_parameter("Cm", 0.185, required);
      add_parameter("sV", 0.2, required);
      add_parameter("rho", 162.0, required);
      add_parameter("V_c", 16.404E-3, required);
      add_parameter("V_sr", 1.094E-3, required);
      add_parameter("V_ss", 5.468E-5, required);
      add_parameter("K_o", 5.4, required);
      add_parameter("Na_o", 140.0, required);
      add_parameter("Ca_o", 2.0, required);
      add_parameter("G_Na", 14.838, required);
      add_parameter("G_K1", 5.405, required);
      add_parameter("G_to", 0.294, required);
      add_parameter("G_Kr", 0.153, required);
      add_parameter("G_Ks", 0.392, required);
      add_parameter("p_KNa", 3.E-2, required);
      add_parameter("G_CaL", 3.98E-5, required);
      add_parameter("K_NaCa", 1000., required);
      add_parameter("gamma", 0.35, required);
      add_parameter("K_mCa", 1.38, required);
      add_parameter("K_mNai", 87.5, required);
      add_parameter("K_sat", 0.1, required);
      add_parameter("alpha", 2.5, required);
      add_parameter("p_NaK", 2.724, required);
      add_parameter("K_mK", 1., required);
      add_parameter("K_mNa", 40., required);
      add_parameter("G_pK", 1.46E-2, required);
      add_parameter("G_pCa", 0.1238, required);
      add_parameter("K_pCa", 5.E-4, required);
      add_parameter("G_bNa", 2.9E-4, required);
      add_parameter("G_bCa", 5.92E-4, required);
      add_parameter("Vmax_up", 6.375E-3, required);
      add_parameter("K_up", 2.5E-4, required);
      add_parameter("V_rel", 0.102, required);
      add_parameter("k1p", 0.15, required);
      add_parameter("k2p", 4.5E-2, required);
      add_parameter("k3", 6.E-2, required);
      add_parameter("k4", 5.E-3, required);
      add_parameter("EC", 1.5, required);
      add_parameter("max_sr", 2.5, required);
      add_parameter("min_sr", 1., required);
      add_parameter("V_leak", 3.6E-4, required);
      add_parameter("V_xfer", 3.8E-3, required);
      add_parameter("Buf_c", 0.2, required);
      add_parameter("K_bufc", 1.E-3, required);
      add_parameter("Buf_sr", 10., required);
      add_parameter("K_bufsr", 0.3, required);
      add_parameter("Buf_ss", 0.4, required);
      add_parameter("K_bufss", 2.5E-4, required);
    }
  };

  /// Constructor.
  TTP()
      : IonicModel(initial_X, initial_Xg,
                   /* Vrest_ = */ -85.23, /* Vscale_ = */ 1.0,
                   /* Tscale_ = */ 1.0, /* Voffset_ = */ 0.0) {}

  /// Construct an instance of model parameters.
  virtual std::unique_ptr<IonicModelParameters>
  get_parameters() const override {
    return std::make_unique<Parameters>();
  }

  /// Read model parameters from a parameter object.
  virtual void read_parameters(const IonicModelParameters &params) override;

  /// Distribute model parameters to all parallel processes.
  virtual void distribute_parameters(const CmMod &cm_mod,
                                     const cmType &cm) override;

  /// Get the index of Ca_i in the state vector.
  virtual unsigned int get_calcium_index() const override {
    return calcium_index;
  }

protected:
  /// @name Model parameters
  /// @{

  /// Gas constant [J/mol/K]
  double Rc = 8314.472;

  /// Temperature [K]
  double Tc = 310.0;

  /// Faraday constant [C/mmol]
  double Fc = 96485.3415;

  /// Cell capacitance per unit surface area [uF/cm^{2}]
  double Cm = 0.185;

  /// Surface to volume ratio [um^{-1}]
  double sV = 0.2;

  /// Cellular resistivity [\f$\Omega\f$-cm]
  double rho = 162.0;

  /// Cytoplasmic volume [um^{3}]
  double V_c = 16.404E-3;

  /// Sacroplasmic reticulum volume [um^{3}]
  double V_sr = 1.094E-3;

  /// Subspace volume [um^{3}]
  double V_ss = 5.468E-5;

  /// Extracellular K concentration [mM]
  double K_o = 5.4;

  /// Extracellular Na concentration [mM]
  double Na_o = 140.0;

  /// Extracellular Ca concentration [mM]
  double Ca_o = 2.0;

  /// Maximal I_Na conductance [nS/pF]
  double G_Na = 14.838;

  /// Maximal I_K1 conductance [nS/pF]
  double G_K1 = 5.405;

  /// Maximal I_to conductance [nS/pF]
  double G_to = 0.294;

  /// Maximal I_Kr conductance [nS/pF]
  double G_Kr = 0.153;

  /// Maximal I_Ks conductance [nS/pF]
  double G_Ks = 0.392;

  /// Relative I_Ks permeability to Na [-]
  double p_KNa = 3.E-2;

  /// Maximal I_CaL conductance [cm^{3}/uF/ms]
  double G_CaL = 3.98E-5;

  /// Maximal I_NaCa [pA/pF]
  double K_NaCa = 1000.;

  /// Voltage dependent parameter of I_NaCa [-]
  double gamma = 0.35;

  /// Ca_i half-saturation constant for I_NaCa [mM]
  double K_mCa = 1.38;

  /// Na_i half-saturation constant for I_NaCa [mM]
  double K_mNai = 87.5;

  /// Saturation factor for I_NaCa [-]
  double K_sat = 0.1;

  /// Factor enhancing outward nature of I_NaCa [-]
  double alpha = 2.5;

  /// Maximal I_NaK [pA/pF]
  double p_NaK = 2.724;

  /// K_o half-saturation constant of I_NaK [mM]
  double K_mK = 1.;

  /// Na_i half-saturation constant of I_NaK [mM]
  double K_mNa = 40.;

  /// Maximal I_pK conductance [nS/pF]
  double G_pK = 1.46E-2;

  /// Maximal I_pCa conductance [pA/pF]
  double G_pCa = 0.1238;

  /// Half-saturation constant of I_pCa [mM]
  double K_pCa = 5.E-4;

  /// Maximal I_bNa conductance [nS/pF]
  double G_bNa = 2.9E-4;

  /// Maximal I_bCa conductance [nS/pF]
  double G_bCa = 5.92E-4;

  /// Maximal I_up conductance [mM/ms]
  double Vmax_up = 6.375E-3;

  /// Half-saturation constant of I_up [mM]
  double K_up = 2.5E-4;

  /// Maximal I_rel conductance [mM/ms]
  double V_rel = 0.102;

  /// R to O and RI to I, I_rel transition rate [mM^{-2}/ms]
  double k1p = 0.15;

  /// O to I and R to RI, I_rel transition rate [mM^{-1}/ms]
  double k2p = 4.5E-2;

  /// O to R and I to RI, I_rel transition rate [ms^{-1}]
  double k3 = 6.E-2;

  /// I to O and Ri to I, I_rel transition rate [ms^{-1}]
  double k4 = 5.E-3;

  /// Ca_sr half-saturation constant of k_casr [mM]
  double EC = 1.5;

  /// Maximum value of k_casr [-]
  double max_sr = 2.5;

  /// Minimum value of k_casr [-]
  double min_sr = 1.;

  /// Maximal I_leak conductance [mM/ms]
  double V_leak = 3.6E-4;

  /// Maximal I_xfer conductance [mM/ms]
  double V_xfer = 3.8E-3;

  /// Total cytoplasmic buffer concentration [mM]
  double Buf_c = 0.2;

  /// Ca_i half-saturation constant for cytplasmic buffer [mM]
  double K_bufc = 1.E-3;

  /// Total sacroplasmic buffer concentration [mM]
  double Buf_sr = 10.;

  /// Ca_sr half-saturation constant for subspace buffer [mM]
  double K_bufsr = 0.3;

  /// Total subspace buffer concentration [mM]
  double Buf_ss = 0.4;

  /// Ca_ss half-saturation constant for subspace buffer [mM]
  double K_bufss = 2.5E-4;

  /// @}

  /**
   * @brief Update gating variables.
   *
   * The evolution rate of all the gating variables for this model has an
   * expression in the form:
   * @f[
   *   \frac{\text{d}y}{\text{d}t} = \frac{y_\infty(v) - y}{\tau(v)}\;,
   * @f]
   * with @f$v@f$ denotes the transmembrane potential. Assuming @f$v@f$ to be
   * constant over the integration step, the above equation admits the following
   * exact solution:
   * @f[
   *   y^{n+1} = y_\infty(v) - (y_\infty(v) - y^n) e^{-\Delta t / \tau(v)}\;.
   * @f]
   * Accordingly, this function goes through all gating variables, computes
   * the values of @f$y_\infty(v)\f$ and @f$\tau(v)\f$, and updates the
   * gating variables using the above formula.
   *
   * @param[in] zone_id Identifier for the transmural zone (epicardium,
   *   endocardium, myocardium).
   * @param[in] dt Time step.
   * @param[in] X Vector of state variables.
   * @param[out] Xg Vector of gating variables.
   */
  virtual void update_g(const unsigned int zone_id, const double dt,
                        const Vector<double> &X,
                        Vector<double> &Xg) const override;

  /// Model right-hand side.
  virtual Vector<double> getf(const unsigned int zone_id,
                              const Vector<double> &X, const Vector<double> &Xg,
                              const double I_stim,
                              const double I_sac) const override;

  /// Model jacobian.
  virtual Array<double> getj(const unsigned int zone_id,
                             const Vector<double> &X, const Vector<double> &Xg,
                             const double Ksac) const override;
};

#endif
