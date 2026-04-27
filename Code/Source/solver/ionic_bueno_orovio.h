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
  /// Constructor.
  BuenoOrovio()
      : IonicModel(/* Vrest_ = */ -84.0, /* Vscale_ = */ 85.70,
                   /* Tscale_ = */ 1.0, /* Voffset_ = */ -84.0) {}

  /// Setup of initial conditions.
  virtual void init(const int nX, Vector<double> &X) const override;

protected:
  /// @name Model parameters
  /// @todo Document units of measure.
  /// @{

  /// Alias for model parameters container. The three entries in each of these
  /// correspond to epicardium, endocardium and myocardium, respectively (see
  /// also table 1 in the reference paper).
  /// \todo [TODO:DaveP] these guys should be maps map<int,double>.
  using ModelParam = Vector<double>;

  const ModelParam u_o = {0.0, 0.0, 0.0};
  const ModelParam u_u = {1.550, 1.56, 1.61};
  const ModelParam theta_v = {0.30, 0.3, 0.3};
  const ModelParam theta_w = {0.130, 0.13, 0.13};
  const ModelParam thetam_v = {6.E-3, 0.2, 0.1};
  const ModelParam theta_o = {6.E-3, 6.E-3, 5.E-3};
  const ModelParam taum_v1 = {60.0, 75., 80.};
  const ModelParam taum_v2 = {1.15E3, 10., 1.4506};
  const ModelParam taup_v = {1.45060, 1.4506, 1.4506};
  const ModelParam taum_w1 = {60.0, 6., 70.};
  const ModelParam taum_w2 = {15.0, 140., 8.};
  const ModelParam km_w = {65.0, 200., 200.};
  const ModelParam um_w = {3.E-2, 1.6E-2, 1.6E-2};
  const ModelParam taup_w = {200.0, 280., 280.};
  const ModelParam tau_fi = {0.110, 0.1, 0.078};
  const ModelParam tau_o1 = {400.0, 470., 410.};
  const ModelParam tau_o2 = {6.0, 6., 7.};
  const ModelParam tau_so1 = {30.01810, 40., 91.};
  const ModelParam tau_so2 = {0.99570, 1.2, 0.8};
  const ModelParam k_so = {2.04580, 2., 2.1};
  const ModelParam u_so = {0.650, 0.65, 0.6};
  const ModelParam tau_s1 = {2.73420, 2.7342, 2.7342};
  const ModelParam tau_s2 = {16.0, 2., 2.};
  const ModelParam k_s = {2.09940, 2.0994, 2.0994};
  const ModelParam u_s = {0.90870, 0.9087, 0.9087};
  const ModelParam tau_si = {1.88750, 2.9013, 3.3849};
  const ModelParam tau_winf = {7.E-2, 2.73E-2, 1.E-2};
  const ModelParam ws_inf = {0.940, 0.78, 0.5};

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

  /// Step function.
  inline double step(const double r) const { return r < 0.0 ? 0.0 : 1.0; }

  /// Delta function.
  inline double delta(const double r) const {
    return utils::is_zero(r) ? 1.0 : 0.0;
  }
};

#endif
