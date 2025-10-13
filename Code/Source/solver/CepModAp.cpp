// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "CepModAp.h"

#include "mat_fun.h"
#include <math.h>

CepModAp::CepModAp()
{
}

CepModAp::~CepModAp()
{
}

/// @brief Compute activation force for electromechanics based on active stress model
///
/// Replicates 'SUBROUTINE AP_ACTVSTRS(X, dt, Tact, epsX)' defined in 'CEPMOD_AP.f'.
void CepModAp::actv_strs(const double X, const double dt, double& Tact, double& epsX)
{
  epsX = exp(-exp(-xi_T*(X - Vcrit)));
  epsX = eps_0 + (eps_i - eps_0)*epsX;
  double nr = Tact + epsX*dt*eta_T*(X - Vrest);
  Tact = nr / (1.0 + epsX*dt);
}

void CepModAp::getf(const int n, const Vector<double>& X, Vector<double>& f, const double fext)
{
  f(0) = X(0)*(c*(X(0)-alpha)*(1.0-X(0)) - X(1)) + fext;
  f(1) = (a + mu1*X(1)/(mu2 + X(0))) * (-X(1) - c*X(0)*(X(0) - b - 1.0));
}

void CepModAp::getj(const int n, const Vector<double>& X, Array<double>& JAC, const double Ksac)
{
  JAC = 0.0;

  double n1 = X(0) - alpha;
  double n2 = 1.0 - X(0);

  JAC(0,0) = c * (n1*n2 + X(0)*(n2 - n1)) - X(1) - Ksac;
  JAC(0,1) = -X(0);

  n1 = mu1*X(1) / (mu2 + X(0));
  n2 = n1 / (mu2 + X(0));
  double n3 = X(1) + c*X(0)*(X(0) - b - 1.0);

  JAC(1,0) = n2*n3 - c*(a+n1)*(2.0*X(0) - b - 1.0);

  n1 = mu1 / (mu2 + X(0));
  n2 = a + n1*X(1);
  n3 = -n3;
  JAC(1,1) = n1*n3 - n2;
}

/// @brief SUBROUTINE AP_INIT0(nX, X)
void CepModAp::init(const int nX, Vector<double> &X)
{
  X = 1.e-3;
  X(0) = Voffset;
}

/// @brief SUBROUTINE AP_INITS(nX, X, X0)
void CepModAp::init(const int nX, Vector<double> &X, double X0)
{
  X = X0;
}

/// @brief SUBROUTINE AP_INITV(nX, X, X0)
void CepModAp::init(const int nX, Vector<double> &X, Vector<double>& X0)
{
  X = X0;
}

/// @brief Time integration performed using Crank-Nicholson method
///
/// Replicates 'SUBROUTINE AP_INTEGCN2(nX, Xn, Ts, Ti, Istim, Ksac, IPAR, RPAR)' defined in 'CEPMOD_AP.f'.
void CepModAp::integ_cn2(const int nX, Vector<double>& Xn, const double Ts, const double Ti, const double Istim, const double Ksac,
    Vector<int>& IPAR, Vector<double>& RPAR)
{
  int itMax = IPAR(0);
  double atol  = RPAR(0);
  double rtol  = RPAR(1);

  double t = Ts / Tscale;
  double dt = Ti / Tscale;
  double Isac  = Ksac * (Vrest - Xn(0));
  double fext  = (Istim + Isac) * Tscale / Vscale;
  Xn(0) = (Xn(0) - Voffset) / Vscale;
  auto Im = mat_fun::mat_id(nX);

  Vector<double> fn(nX);

  getf(nX, Xn, fn, fext);
  //CALL AP_GETF(nX, Xn, fn, fext)

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
    //CALL AP_GETF(nX, Xk, fk, fext)
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
    getj(nX, Xk, JAC, Ksac*Tscale);
    //CALL AP_GETJ(nX, Xk, JAC, Ksac*Tscale)

    JAC = Im - 0.50 * dt * JAC;
    JAC = mat_fun::mat_inv(JAC, nX);
    rK = mat_fun::mat_mul(JAC, rK);
    Xk = Xk - rK;
  }

  Xn = Xk;

  getf(nX, Xn, fn, fext);
  //CALL AP_GETF(nX, Xn, fn, fext)

  Xn(0) = Xn(0)*Vscale + Voffset;

  if (!l2 && !l3) {
    IPAR(1) = IPAR(1) + 1;
  }
}

/// @brief Time integration performed using Forward Euler method
///
/// Replicates 'SUBROUTINE AP_INTEGFE(nX, X, Ts, Ti, Istim, Ksac)' defined in 'CEPMOD_AP.f'.
void CepModAp::integ_fe(const int nX, Vector<double>& X, const double Ts, const double Ti, const double Istim, const double Ksac)
{
  double t = Ts / Tscale;
  double dt = Ti / Tscale;
  double Isac = Ksac * (Vrest - X(0));
  double fext = (Istim + Isac) * Tscale / Vscale;

  X(0) = (X(0) - Voffset) / Vscale;
  Vector<double> f(nX);
  getf(nX, X, f, fext);
  //CALL AP_GETF(nX, X, f, fext)
  X = X + dt*f;
  X(0) = X(0)*Vscale + Voffset;
}

/// @brief Time integration performed using 4th order Runge-Kutta method
///
/// Replicates 'SUBROUTINE AP_INTEGRK(nX, X, Ts, Ti, Istim, Ksac)' defined in 'CEPMOD_AP.f'.
void CepModAp::integ_rk(const int nX, Vector<double>& X, const double Ts, const double Ti, const double Istim, const double Ksac)
{
  double t = Ts / Tscale;
  double dt = Ti / Tscale;
  double dt6  = dt / 6.0;

  double Isac = Ksac * (Vrest - X(0));
  double fext = (Istim + Isac) * Tscale / Vscale;
  X(0) = (X(0) - Voffset) / Vscale;

  Array<double> frk(nX,4);
  Vector<double> Xrk(nX);
  Vector<double> frk1(nX), frk2(nX), frk3(nX), frk4(nX);

  // RK4: 1st pass
  Xrk = X;
  getf(nX, Xrk, frk1, fext);
  //CALL AP_GETF(nX, Xrk, frk(:,1), fext)

  // RK4: 2nd pass
  Xrk  = X + 0.5 * dt * frk1;
  getf(nX, Xrk, frk2, fext);
  //CALL AP_GETF(nX, Xrk, frk(:,2), fext)

  // RK4: 3rd pass
  Xrk  = X + 0.5 * dt * frk2;
  getf(nX, Xrk, frk3, fext);
  //CALL AP_GETF(nX, Xrk, frk(:,3), fext)

  // RK4: 4th pass
  Xrk  = X + dt*frk3;
  getf(nX, Xrk, frk4, fext);
  //CALL AP_GETF(nX, Xrk, frk(:,4), fext)

  X = X + dt6 * (frk1 + 2.0*(frk2 + frk3) + frk4);
  //X = X + dt6*(frk(:,1) + 2._RKIND*(frk(:,2) + frk(:,3)) + frk(:,4))

  X(0) = X(0)*Vscale + Voffset;

}


