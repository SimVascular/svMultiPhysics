// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef FC_TYPE_H
#define FC_TYPE_H

#include "Array.h"
#include "Vector.h"

/// @brief Fourier coefficients that are used to specify time-dependent
/// functions.
///
/// These are used for boundary conditions and for time-dependent prescribed
/// active stress.
class fcType {
public:
  bool defined() { return n != 0; };

  /// Toggle whether this is a ramp function or not.
  bool lrmp = false;

  /// Number of Fourier coefficients.
  int n = 0;

  /// Numbero of dimensions (scalar or vector).
  int d = 0;

  /// Initial value.
  Vector<double> qi;

  /// Time derivative of linear part.
  Vector<double> qs;

  /// Period.
  double T = 0.0;

  /// Initial time.
  double ti = 0.0;

  /// Imaginary part of coefficients.
  Array<double> i;

  /// Real part of coefficients.
  Array<double> r;
};

#endif