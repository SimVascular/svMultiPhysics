// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ionic_aliev_panfilov.h"

void AlievPanfilov::init(const int nX, Vector<double> &X) const {
  X = 1.0e-3;
  X(0) = Voffset;
}

void AlievPanfilov::getf(const unsigned int zone_id, const int nX, const int nG,
                         const Vector<double> &X, const Vector<double> &Xg,
                         Vector<double> &f, const double I_stim,
                         const double I_sac) const {
  f(0) = X(0) * (c * (X(0) - alpha) * (1.0 - X(0)) - X(1)) - I_stim + I_sac;
  f(1) =
      (a + mu1 * X(1) / (mu2 + X(0))) * (-X(1) - c * X(0) * (X(0) - b - 1.0));
}

void AlievPanfilov::getj(const unsigned int zone_id, const int nX, const int nG,
                         const Vector<double> &X, const Vector<double> &Xg,
                         Array<double> &Jac, const double Ksac) const {
  Jac = 0.0;

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
}