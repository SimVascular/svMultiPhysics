// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ionic_aliev_panfilov.h"

void AlievPanfilov::read_parameters(const IonicModelParameters &params) {
  IonicModel::read_parameters(params);

  alpha = params.get_scalar("alpha");
  a = params.get_scalar("a");
  b = params.get_scalar("b");
  c = params.get_scalar("c");
  mu1 = params.get_scalar("mu1");
  mu2 = params.get_scalar("mu2");
}

void AlievPanfilov::distribute_parameters(const CmMod &cm_mod,
                                          const cmType &cm) {
  IonicModel::distribute_parameters(cm_mod, cm);

  cm.bcast(cm_mod, &alpha);
  cm.bcast(cm_mod, &a);
  cm.bcast(cm_mod, &b);
  cm.bcast(cm_mod, &c);
  cm.bcast(cm_mod, &mu1);
  cm.bcast(cm_mod, &mu2);
}

Vector<double> AlievPanfilov::getf(const unsigned int zone_id,
                                   const Vector<double> &X,
                                   const Vector<double> &Xg,
                                   const double I_stim,
                                   const double I_sac) const {
  Vector<double> f(X.size());

  f(0) = X(0) * (c * (X(0) - alpha) * (1.0 - X(0)) - X(1)) - I_stim + I_sac;
  f(1) =
      (a + mu1 * X(1) / (mu2 + X(0))) * (-X(1) - c * X(0) * (X(0) - b - 1.0));

  return f;
}

Array<double> AlievPanfilov::getj(const unsigned int zone_id,
                                  const Vector<double> &X,
                                  const Vector<double> &Xg,
                                  const double Ksac) const {
  Array<double> Jac(X.size(), X.size());

  double n1 = X(0) - alpha;
  double n2 = 1.0 - X(0);

  Jac(0, 0) = c * (n1 * n2 + X(0) * (n2 - n1)) - X(1) - Ksac;
  Jac(0, 1) = -X(0);

  n1 = mu1 * X(1) / (mu2 + X(0));
  n2 = n1 / (mu2 + X(0));
  double n3 = X(1) + c * X(0) * (X(0) - b - 1.0);

  Jac(1, 0) = n2 * n3 - c * (a + n1) * (2.0 * X(0) - b - 1.0);

  n1 = mu1 / (mu2 + X(0));
  n2 = a + n1 * X(1);
  n3 = -n3;
  Jac(1, 1) = n1 * n3 - n2;

  return Jac;
}

REGISTER_IONIC_MODEL("AP", AlievPanfilov);