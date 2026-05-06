// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef IONIC_FITZHUGH_NAGUMO_H
#define IONIC_FITZHUGH_NAGUMO_H

#include "ionic_model.h"

#include "Vector.h"
#include "utils.h"

/**
 * @brief FitzHugh-Nagumo ionic model.
 *
 * **References**:
 * - FitzHugh, Impulses and physiological states in theoretical models of nerve
 *   membrane. Biophysical Journal (1961)
 * - Nagumo, Arimoto, Yoshizawa. An active pulse transmission line simulating
 *   nerve axon. Proceedings of the IRE (1962).
 */
class FitzHughNagumo : public IonicModel {
public:
  /// Model label.
  static inline const std::string label = "FN";

  /// State variables.
  static inline const InitialStates initial_X = {{"V", 1.0e-3}, {"w", 1.0e-3}};

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

      add_parameter("alpha", -0.50, required);
      add_parameter("a", 0.0, required);
      add_parameter("b", -0.60, required);
      add_parameter("c", 50.0, required);
    }
  };

  /// Constructor.
  FitzHughNagumo()
      : IonicModel(initial_X, initial_Xg,
                   /* Vrest_ = */ 0.0, /* Vscale_ = */ 1.0,
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

protected:
  /// @name Model parameters
  /// @todo Document units of measure.
  /// @{

  double alpha = -0.50;
  double a = 0.0;
  double b = -0.60;
  double c = 50.0;

  /// @}

  /// Update variable with analytical solution. This model has none, so this
  /// method does nothing.
  virtual void update_g(const unsigned int zone_id, const double dt,
                        const int nX, const int nG, const Vector<double> &X,
                        Vector<double> &Xg) const override {}

  /// Model right-hand side.
  virtual void getf(const unsigned int zone_id, const int nX, const int nG,
                    const Vector<double> &X, const Vector<double> &Xg,
                    Vector<double> &f, const double I_stim,
                    const double I_sac) const override;

  /// Model jacobian.
  virtual void getj(const unsigned int zone_id, const int nX, const int nG,
                    const Vector<double> &X, const Vector<double> &Xg,
                    Array<double> &Jac, const double Ksac) const override;
};

#endif
