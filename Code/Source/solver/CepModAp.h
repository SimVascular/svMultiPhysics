// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef CEP_MOD_AP_H 
#define CEP_MOD_AP_H 

#include "Array.h"
#include "Vector.h"

/// @brief This module defines data structures for Aliev-Panfilov cellular
/// activation model for cardiac electrophysiology.
///
/// The classes defined here duplicate the data structures in the Fortran APMOD module defined in CEPMOD_AP.f 
/// and PARAMS_AP.f files. 
class CepModAp
{
  public:
    CepModAp();
    ~CepModAp();

    // Scaling factors
    /// Voltage scaling
    double Vscale  = 100.0;

    /// Time scaling
    double Tscale  = 12.90;

    /// Voltage offset parameter
    double Voffset = -80.0;

    /// Model parameters
    double alpha = 1.E-2;
    double a     = 2.E-3;
    double b     = 0.150;
    double c     = 8.0;
    double mu1   = 0.20;
    double mu2   = 0.30;

    // Electromechanics coupling parameters: active stress model
    /// Resting voltage (mV)
    double Vrest = -80.0;

    /// Critical voltage (mV)
    double Vcrit = -30.0;

    /// Saturation potential
    double eta_T = 5.E-3;

    /// Minimum activation (ms^{-1})
    double eps_0 = 0.10;

    /// Maximum activation (ms^{-1})
    double eps_i = 1.0;

    /// Transition rate (mV^{-1})
    double xi_T  = 1.0;

    /// Cm: Cell capacitance per unit surface area
    double Cm  = 1.0;

    /// sV: Surface to volume ratio
    double sV  = 1.0;

    /// rho: Cellular resistivity
    double rho = 1.0;

    void actv_strs(const double X, const double dt, double& Tact, double& epsX);

    void getf(const int n, const Vector<double>& X, Vector<double>& f, const double fext);
    void getj(const int n, const Vector<double>& X, Array<double>& Jac, const double Ksac);

    void init(const int nX, Vector<double> &X);
    void init(const int nX, Vector<double> &X, double X0);
    void init(const int nX, Vector<double> &X, Vector<double>& X0);

    void integ_cn2(const int nX, Vector<double>& Xn, const double Ts, const double Ti, 
        const double Istim, const double Ksac, Vector<int>& IPAR, Vector<double>& RPAR);
    void integ_fe(const int nX, Vector<double>& X, const double Ts, const double Ti, 
        const double Istim, const double Ksac);
    void integ_rk(const int nX, Vector<double>& X, const double Ts, const double Ti, 
        const double Istim, const double Ksac);
};

#endif

