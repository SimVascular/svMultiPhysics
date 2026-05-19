// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef IONIC_ALIEV_PANFILOV_H
#define IONIC_ALIEV_PANFILOV_H

#include "ionic_model.h"

#include "Vector.h"

/**
 * @brief Aliev-Panfilov ionic model.
 *
 * **Reference**: Aliev, Panfilov. A simple two-variable model of cardiac
 * excitation. Chaos, Solitons and Fractals (1996).
 */
class AlievPanfilov : public IonicModel {
public:
  /// Model label.
  static inline const std::string label = "AP";

  /// State variables.
  static inline const InitialStates initial_X = {{"V", -80.0}, {"w", 1.0e-3}};

  /// Gating variables.
  static inline const InitialStates initial_Xg = {};

  /// Index of the recovery variable (w), used as calcium proxy for
  /// electromechanical coupling.
  static constexpr unsigned int calcium_index = 1;

  /// Model parameters class.
  class Parameters : public IonicModelParameters {
  public:
    Parameters() : IonicModelParameters(label, initial_X, initial_Xg) {
      constexpr bool required = true;

      add_parameter("alpha", 1.0e-2, required);
      add_parameter("a", 2.0e-3, required);
      add_parameter("b", 0.15, required);
      add_parameter("c", 8.0, required);
      add_parameter("mu1", 0.20, required);
      add_parameter("mu2", 0.30, required);
    }
  };

  /// Constructor.
  AlievPanfilov()
      : IonicModel(initial_X, initial_Xg,
                   /* Vrest_ = */ -80.0, /* Vscale_ = */ 100.0,
                   /* Tscale_ = */ 12.90, /* Voffset_ = */ -80.0) {}

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

  /// Corresponding to parameter a in Aliev-Panfilov paper [1].
  double alpha = 1.0e-2;

  /// Corresponding to parameter epsilon0 in Aliev-Panfilov paper [1].
  double a = 2.0e-3;

  /// Corresponding to parameter a in Aliev-Panfilov paper [1].
  double b = 0.15;

  /// Corresponding to parameter k in Aliev-Panfilov paper [1].
  double c = 8.0;

  double mu1 = 0.20; ///< [1].
  double mu2 = 0.30; ///< [1].

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
};

#endif
