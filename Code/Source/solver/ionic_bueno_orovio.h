// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef IONIC_BUENO_OROVIO_H
#define IONIC_BUENO_OROVIO_H

#include "ionic_model.h"

#include "Vector.h"
#include "utils.h"

/**
 * @brief Bueno-Orovio ionic model.
 *
 * **Reference**: Bueno-Orovio, Cherry, Fenton. Minimal model for human
 * ventricular action potentials in tissue. Journal of Theoretical Biology
 * (2008)
 */
class BuenoOrovio : public IonicModel {
public:
  /// Model label.
  static inline const std::string label = "BO";

  /// State variables.
  static inline const InitialStates initial_X = {
      {"u", -84.0}, {"v", 1.0}, {"w", 1.0}, {"s", 0.0}};

  /// Gating variables.
  static inline const InitialStates initial_Xg = {};

  /// Index of the slow inward current gate (s), used as calcium proxy for
  /// electromechanical coupling.
  static constexpr unsigned int calcium_index = 3;

  /// Model parameters class.
  class Parameters : public IonicModelParameters {
  public:
    Parameters() : IonicModelParameters(label, initial_X, initial_Xg) {
      constexpr bool required = true;

      add_parameter("u_o", {0.0, 0.0, 0.0}, required);
      add_parameter("u_u", {1.550, 1.56, 1.61}, required);
      add_parameter("theta_v", {0.30, 0.3, 0.3}, required);
      add_parameter("theta_w", {0.130, 0.13, 0.13}, required);
      add_parameter("thetam_v", {6.E-3, 0.2, 0.1}, required);
      add_parameter("theta_o", {6.E-3, 6.E-3, 5.E-3}, required);
      add_parameter("taum_v1", {60.0, 75., 80.}, required);
      add_parameter("taum_v2", {1.15E3, 10., 1.4506}, required);
      add_parameter("taup_v", {1.45060, 1.4506, 1.4506}, required);
      add_parameter("taum_w1", {60.0, 6., 70.}, required);
      add_parameter("taum_w2", {15.0, 140., 8.}, required);
      add_parameter("km_w", {65.0, 200., 200.}, required);
      add_parameter("um_w", {3.E-2, 1.6E-2, 1.6E-2}, required);
      add_parameter("taup_w", {200.0, 280., 280.}, required);
      add_parameter("tau_fi", {0.110, 0.1, 0.078}, required);
      add_parameter("tau_o1", {400.0, 470., 410.}, required);
      add_parameter("tau_o2", {6.0, 6., 7.}, required);
      add_parameter("tau_so1", {30.01810, 40., 91.}, required);
      add_parameter("tau_so2", {0.99570, 1.2, 0.8}, required);
      add_parameter("k_so", {2.04580, 2., 2.1}, required);
      add_parameter("u_so", {0.650, 0.65, 0.6}, required);
      add_parameter("tau_s1", {2.73420, 2.7342, 2.7342}, required);
      add_parameter("tau_s2", {16.0, 2., 2.}, required);
      add_parameter("k_s", {2.09940, 2.0994, 2.0994}, required);
      add_parameter("u_s", {0.90870, 0.9087, 0.9087}, required);
      add_parameter("tau_si", {1.88750, 2.9013, 3.3849}, required);
      add_parameter("tau_winf", {7.E-2, 2.73E-2, 1.E-2}, required);
      add_parameter("ws_inf", {0.940, 0.78, 0.5}, required);
    }
  };

  /// Constructor.
  BuenoOrovio()
      : IonicModel(initial_X, initial_Xg,
                   /* Vrest_ = */ -84.0, /* Vscale_ = */ 85.70,
                   /* Tscale_ = */ 1.0, /* Voffset_ = */ -84.0) {}

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

  /// Get the calcium proxy index.
  virtual unsigned int get_calcium_index() const override {
    return calcium_index;
  }

protected:
  /// @name Model parameters
  /// @{

  /// Alias for model parameters container. The three entries in each of these
  /// correspond to epicardium, endocardium and myocardium, respectively (see
  /// also table 1 in the reference paper).
  /// \todo [TODO:DaveP] these guys should be maps map<int,double>.
  /// @todo[michelebucelli] Would it make sense to treat these in the same way
  /// as the epi/endo/myo zone are treated in the TTP model?
  using ModelParam = Vector<double>;

  ModelParam u_o = {0.0, 0.0, 0.0};              ///< [1].
  ModelParam u_u = {1.550, 1.56, 1.61};          ///< [1].
  ModelParam theta_v = {0.30, 0.3, 0.3};         ///< [1].
  ModelParam theta_w = {0.130, 0.13, 0.13};      ///< [1].
  ModelParam thetam_v = {6.E-3, 0.2, 0.1};       ///< [1].
  ModelParam theta_o = {6.E-3, 6.E-3, 5.E-3};    ///< [1].
  ModelParam taum_v1 = {60.0, 75., 80.};         ///< [1/s].
  ModelParam taum_v2 = {1.15E3, 10., 1.4506};    ///< [1/s].
  ModelParam taup_v = {1.45060, 1.4506, 1.4506}; ///< [1/s].
  ModelParam taum_w1 = {60.0, 6., 70.};          ///< [1/s].
  ModelParam taum_w2 = {15.0, 140., 8.};         ///< [1/s].
  ModelParam km_w = {65.0, 200., 200.};          ///< [1].
  ModelParam um_w = {3.E-2, 1.6E-2, 1.6E-2};     ///< [1].
  ModelParam taup_w = {200.0, 280., 280.};       ///< [1/s].
  ModelParam tau_fi = {0.110, 0.1, 0.078};       ///< [1/s].
  ModelParam tau_o1 = {400.0, 470., 410.};       ///< [1/s].
  ModelParam tau_o2 = {6.0, 6., 7.};             ///< [1/s].
  ModelParam tau_so1 = {30.01810, 40., 91.};     ///< [1/s].
  ModelParam tau_so2 = {0.99570, 1.2, 0.8};      ///< [1/s].
  ModelParam k_so = {2.04580, 2., 2.1};          ///< [1].
  ModelParam u_so = {0.650, 0.65, 0.6};          ///< [1].
  ModelParam tau_s1 = {2.73420, 2.7342, 2.7342}; ///< [1/s].
  ModelParam tau_s2 = {16.0, 2., 2.};            ///< [1/s].
  ModelParam k_s = {2.09940, 2.0994, 2.0994};    ///< [1].
  ModelParam u_s = {0.90870, 0.9087, 0.9087};    ///< [1].
  ModelParam tau_si = {1.88750, 2.9013, 3.3849}; ///< [1/s].
  ModelParam tau_winf = {7.E-2, 2.73E-2, 1.E-2}; ///< [1/s].
  ModelParam ws_inf = {0.940, 0.78, 0.5};        ///< [1].

  /// @}

  /// Update variable with analytical solution. This model has none, so this
  /// method does nothing.
  virtual void update_g(const unsigned int zone_id, const double dt,
                        const Vector<double> &X,
                        Vector<double> &Xg) const override {}

  /// Model right-hand side.
  virtual Vector<double> getf(const unsigned int zone_id,
                              const Vector<double> &X, const Vector<double> &Xg,
                              const double I_stim,
                              const double I_sac) const override;

  /// Model jacobian.
  virtual Array<double> getj(const unsigned int zone_id,
                             const Vector<double> &X, const Vector<double> &Xg,
                             const double Ksac) const override;

  /// Step function.
  inline double step(const double r) const { return r < 0.0 ? 0.0 : 1.0; }

  /// Delta function.
  inline double delta(const double r) const {
    return utils::is_zero(r) ? 1.0 : 0.0;
  }
};

#endif
