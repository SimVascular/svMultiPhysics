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
  /// Constructor.
  AlievPanfilov()
      : IonicModel(/* states_X_ = */ {{"V", -80.0}, {"w", 1.0e-3}},
                   /* states_Xg_ = */ {},
                   /* Vrest_ = */ -80.0, /* Vscale_ = */ 100.0,
                   /* Tscale_ = */ 12.90, /* Voffset_ = */ -80.0) {}

protected:
  /// @name Model parameters
  /// @todo Document units of measure.
  /// @{

  /// Corresponding to parameter a in Aliev-Panfilov paper.
  const double alpha = 1.0e-2;

  /// Corresponding to parameter epsilon0 in Aliev-Panfilov paper.
  const double a = 2.0e-3;

  /// Corresponding to parameter a in Aliev-Panfilov paper.
  const double b = 0.15;

  /// Corresponding to parameter k in Aliev-Panfilov paper.
  const double c = 8.0;

  const double mu1 = 0.20;
  const double mu2 = 0.30;

  /// Cell capacitance per unit surface area.
  const double Cm = 1.0;

  /// Membrane surface to volume ratio.
  const double sV = 1.0;

  /// Cellular resistivity.
  const double rho = 1.0;

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
