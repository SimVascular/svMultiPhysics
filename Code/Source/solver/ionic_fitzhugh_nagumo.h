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
  /// Constructor.
  FitzHughNagumo()
      : IonicModel(/* initial_X_ = */ {{"V", 1.0e-3}, {"w", 1.0e-3}},
                   /* initial_Xg_ = */ {},
                   /* Vrest_ = */ 0.0, /* Vscale_ = */ 1.0,
                   /* Tscale_ = */ 1.0, /* Voffset_ = */ 0.0) {}

protected:
  /// @name Model parameters
  /// @todo Document units of measure.
  /// @{

  const double alpha = -0.50;
  const double a = 0.0;
  const double b = -0.60;
  const double c = 50.0;

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
