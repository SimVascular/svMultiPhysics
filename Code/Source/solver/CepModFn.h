// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef CEP_MOD_FN_H 
#define CEP_MOD_FN_H 

#include "Array.h"
#include "Vector.h"

/// @brief This module defines data structures for Fitzhugh-Nagumo cellular
/// activation model for cardiac electrophysiology.
///
/// The classes defined here duplicate the data structures in the Fortran FNMOD module defined in CEPMOD_FN.f 
/// and PARAMS_FN.f files. 
class CepModFn
{
  public:
    CepModFn();
    ~CepModFn();

    // Scaling factors
    /// Voltage scaling
    double Vscale  = 1.0;
    /// Time scaling
    double Tscale  = 1.0;
    /// Voltage offset parameter
    double Voffset = 0.0;

    /// Model parameters
    double alpha = -0.50;
    double a = 0.0;
    double b = -0.60;
    double c = 50.0;

    void getf(const int n, const Vector<double>& X, Vector<double>& f, const double fext);
    void getj(const int n, const Vector<double>& X, Array<double>& JAC);

    void init(const int nX, Vector<double> &X);
    void init(const int nX, Vector<double> &X, double X0);
    void init(const int nX, Vector<double> &X, Vector<double>& X0);

    void integ_cn2(const int nX, Vector<double>& Xn, const double Ts, const double Ti, const double Istim, 
        Vector<int>& IPAR, Vector<double>& RPAR);
    void integ_fe(const int nX, Vector<double>& X, const double Ts, const double Ti, const double Istim);
    void integ_rk(const int nX, Vector<double>& X, const double Ts, const double Ti, const double Istim);

};

#endif

