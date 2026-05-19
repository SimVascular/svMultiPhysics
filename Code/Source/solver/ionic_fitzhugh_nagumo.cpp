// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ionic_fitzhugh_nagumo.h"

void FitzHughNagumo::read_parameters(const IonicModelParameters &params) {
  IonicModel::read_parameters(params);

  alpha = params.get_scalar("alpha");
  a = params.get_scalar("a");
  b = params.get_scalar("b");
  c = params.get_scalar("c");
}

void FitzHughNagumo::distribute_parameters(const CmMod &cm_mod,
                                           const cmType &cm) {
  IonicModel::distribute_parameters(cm_mod, cm);

  cm.bcast(cm_mod, &alpha);
  cm.bcast(cm_mod, &a);
  cm.bcast(cm_mod, &b);
  cm.bcast(cm_mod, &c);
}

Vector<double> FitzHughNagumo::getf(const unsigned int zone_id,
                                    const Vector<double> &X,
                                    const Vector<double> &Xg,
                                    const double I_stim,
                                    const double I_sac) const {
  Vector<double> f(X.size());

  f(0) = c * (X(0) * (X(0) - alpha) * (1.0 - X(0)) - X(1)) + I_stim + I_sac;
  f(1) = X(0) - b * X(1) + a;

  return f;
}

Array<double> FitzHughNagumo::getj(const unsigned int zone_id,
                                   const Vector<double> &X,
                                   const Vector<double> &Xg,
                                   const double Ksac) const {
  Array<double> Jac(X.size(), X.size());

  double n1 = -3.0 * pow(X(0), 2.0);
  double n2 = 2.0 * (1.0 + alpha) * X(0);

  Jac(0, 0) = c * (n1 + n2 - alpha) - Ksac;
  Jac(0, 1) = -c;
  Jac(1, 0) = 1.0;
  Jac(1, 1) = -b;
}

REGISTER_IONIC_MODEL("FN", FitzHughNagumo);