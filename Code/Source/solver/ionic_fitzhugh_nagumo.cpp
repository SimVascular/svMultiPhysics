// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ionic_fitzhugh_nagumo.h"

void FitzHughNagumo::init(const int nX, Vector<double> &X) const { X = 1.e-3; }

void FitzHughNagumo::getf(const unsigned int zone_id, const int nX,
                          const int nG, const Vector<double> &X,
                          const Vector<double> &Xg, Vector<double> &f,
                          const double fext) const {
  f(0) = c * (X(0) * (X(0) - alpha) * (1.0 - X(0)) - X(1)) + fext;
  f(1) = X(0) - b * X(1) + a;
}

void FitzHughNagumo::getj(const unsigned int zone_id, const int nX,
                          const int nG, const Vector<double> &X,
                          const Vector<double> &Xg, Array<double> &Jac,
                          const double Ksac) const {
  Jac = 0.0;

  double n1 = -3.0 * pow(X(0), 2.0);
  double n2 = 2.0 * (1.0 + alpha) * X(0);

  // @todo This should probably also account for Ksac.
  Jac(0, 0) = c * (n1 + n2 - alpha);
  Jac(0, 1) = -c;
  Jac(1, 0) = 1.0;
  Jac(1, 1) = -b;
}