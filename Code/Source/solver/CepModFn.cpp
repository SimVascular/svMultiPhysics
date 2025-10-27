// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "CepModFn.h"

#include "mat_fun.h"
#include <math.h>

CepModFn::CepModFn()
{
}

CepModFn::~CepModFn()
{
}

void CepModFn::getf(const int n, const Vector<double>& X, Vector<double>& f, const double fext)
{
  f(0) = c * ( X(0)*(X(0)-alpha)*(1.0-X(0)) - X(1) ) + fext;
  f(1) = X(0) - b*X(1) + a;
}

void CepModFn::getj(const int n, const Vector<double>& X, Array<double>& JAC)
{
  JAC = 0.0;

  double n1 = -3.0*pow(X(0),2.0);
  double n2 = 2.0*(1.0+alpha)*X(0);

  JAC(0,0) = c * (n1 + n2 - 1.0);
  JAC(0,1) = -c;
  JAC(1,1) = 1.0;
  JAC(1,1) = -b;
}

/// @brief SUBROUTINE FN_INIT0(nX, X)
void CepModFn::init(const int nX, Vector<double> &X)
{
  X = 1.e-3;
}

/// @brief SUBROUTINE FN_INITS(nX, X, X0)
void CepModFn::init(const int nX, Vector<double> &X, double X0)
{
  X = X0;
}

/// @brief SUBROUTINE FN_INITV(nX, X, X0)
void CepModFn::init(const int nX, Vector<double> &X, Vector<double>& X0)
{
  X = X0;
}

/// @brief Time integration performed using Crank-Nicholson method
void CepModFn::integ_cn2(const int nX, Vector<double>& Xn, const double Ts, const double Ti, const double Istim, 
    Vector<int>& IPAR, Vector<double>& RPAR)
{
  int itMax = IPAR(0);
  double atol  = RPAR(0);
  double rtol  = RPAR(1);

  double t = Ts / Tscale;
  double dt = Ti / Tscale;
  double fext  = Istim * Tscale / Vscale;
  Xn(0) = (Xn(0) - Voffset) / Vscale;
  auto Im = mat_fun::mat_id(nX);

  Vector<double> fn(nX);
  getf(nX, Xn, fn, fext);

  int k  = 0;
  auto Xk = Xn;
  bool l1 = false;
  bool l2 = false;
  bool l3 = false;
  t = Ts + dt;
  double eps = std::numeric_limits<double>::epsilon();

  while (true) {
    k = k + 1;
    Vector<double> fk(nX);
    getf(nX, Xk, fk, fext);
    auto rK = Xk - Xn - 0.5*dt*(fk + fn);

    double rmsA = 0.0;
    double rmsR = 0.0;

    for (int i = 0; i < nX; i++) {
      rmsA = rmsA + pow(rK(i),2.0);
      rmsR = rmsR + pow(rK(i) / (Xk(i)+eps), 2.0);
    }

    rmsA = sqrt(rmsA / static_cast<double>(nX));
    rmsR = sqrt(rmsR / static_cast<double>(nX));

    l1 = (k > itMax);
    l2 = (rmsA <= atol);
    l3  = (rmsR <= rtol);
    if (l1 || l2 || l3) {
      break;
    }

    Array<double> JAC(nX,nX);
    getj(nX, Xk, JAC);

    JAC = Im - 0.50 * dt * JAC;
    JAC = mat_fun::mat_inv(JAC, nX);
    rK = mat_fun::mat_mul(JAC, rK);
    Xk = Xk - rK;
  }
}

/// @brief Time integration performed using Forward Euler method
void CepModFn::integ_fe(const int nX, Vector<double>& X, const double Ts, const double Ti, const double Istim)
{
  double t = Ts / Tscale;
  double dt= Ti / Tscale;
  double fext = Istim * Tscale / Vscale;

  X(0) = (X(0) - Voffset)/Vscale;

  Vector<double> f(nX);
  //CALL FN_GETF(nX, X, f, fext)
  X = X + dt*f;
  X(0) = X(0)*Vscale + Voffset;
}

/// @brief Time integration performed using 4th order Runge-Kutta method
///
/// Replicates 'SUBROUTINE AP_INTEGRK(nX, X, Ts, Ti, Istim, Ksac)' defined in 'CEPMOD_AP.f'.
void CepModFn::integ_rk(const int nX, Vector<double>& X, const double Ts, const double Ti, const double Istim)
{
  double t = Ts / Tscale;
  double dt = Ti / Tscale;
  double dt6  = dt / 6.0;

  double fext = Istim * Tscale / Vscale; 
  X(0) = (X(0) - Voffset) / Vscale;

  Array<double> frk(nX,4);
  Vector<double> Xrk(nX);
  Vector<double> frk1(nX), frk2(nX), frk3(nX), frk4(nX);

  // RK4: 1st pass
  Xrk = X;
  getf(nX, Xrk, frk1, fext);

  // RK4: 2nd pass
  Xrk  = X + 0.5 * dt * frk1;
  getf(nX, Xrk, frk2, fext);

  // RK4: 3rd pass
  Xrk  = X + 0.5 * dt * frk2;
  getf(nX, Xrk, frk3, fext);

  // RK4: 4th pass
  Xrk  = X + dt*frk3;
  getf(nX, Xrk, frk4, fext);

  X = X + dt6 * (frk1 + 2.0*(frk2 + frk3) + frk4);
  X(0) = X(0)*Vscale + Voffset;
}


