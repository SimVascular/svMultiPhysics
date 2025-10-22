// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "CepModBo.h"

#include "mat_fun.h"
#include "utils.h"
#include <math.h>

CepModBo::CepModBo()
{
}

CepModBo::~CepModBo()
{
}

/// @brief Compute macroscopic fiber strain based on sacromere force-length
/// relationship and slow inward current variable (s)
void CepModBo::actv_strn(const double c, const double I4f, const double dt, double& gf)
{
  //  fiber length
  double SL = I4f * SL0;

  // Sacromere force-length relationship
  if (SL >= SLmin && SL <= SLmax) {
     SL = 0.5*f0 + fc1*cos(SL) + fs1*sin(SL) + fc2*cos(2.0*SL) + fs2*sin(2.0*SL)  + fc3*cos(3.0*SL) + fs3*sin(3.0*SL);
  } else { 
     SL = 0.0;
  }

  // Active force
  double Fa = alFa * (c-c0)*(c-c0) * SL;

  double rtmp = 2.0*I4f*(1.0/ pow(1.0+gf,3.0) - 1.0);
  gf = gf + dt*(Fa + rtmp)/(mu_C * c * c);
}

/// @brief Compute activation force for electromechanics based on active stress model
void CepModBo::actv_strs(const double X, const double dt, double& Tact, double& epsX)
{
  epsX = exp(-exp(-xi_T*(X - Vcrit)));
  epsX = eps_0 + (eps_i - eps_0)*epsX;
  double nr   = Tact + epsX*dt*eta_T*(X - Vrest);
  Tact = nr / (1.0 + epsX*dt);
}

double CepModBo::delta(const double r)
{
  double result{0.0};

  if (utils::is_zero(r)) {
    result = 1.0;
  }

  return result;
}

/// @brief The 'zone_id' parameter is the myocardium zone id: 1, 2 or 3.
void CepModBo::getf(const int zone_id, const int n, const Vector<double>& X, Vector<double>& f, const double fext, Vector<double>& RPAR)
{
  // Create local copies of the 4 state variables
  double u = X(0);
  double v = X(1);
  double w = X(2);
  double s = X(3);
  int i = zone_id - 1;

  // Define step functions
  double H_uv = step(u - theta_v[i]);
  double H_uw = step(u - theta_w[i]);
  double H_umv = step(u - thetam_v[i]);
  double H_uo = step(u - theta_o[i]);

  // Define additional constants
  double taum_v = (1.0-H_umv)*taum_v1[i] + H_umv*taum_v2[i];
  double taum_w = taum_w1[i] + 0.5*(taum_w2[i]-taum_w1[i])* (1.0 + tanh(km_w[i]*(u-um_w[i])));
  double tau_so = tau_so1[i] + 0.5*(tau_so2[i]-tau_so1[i])* (1.0 + tanh(k_so[i]*(u-u_so[i])));
  double tau_s  = (1.0-H_uw)*tau_s1[i] + H_uw*tau_s2[i];
  double tau_o  = (1.0-H_uo)*tau_o1[i] + H_uo*tau_o2[i];
  double v_inf  = (1.0-H_umv);
  double w_inf  = (1.0-H_uo)*(1.0 - u/tau_winf[i]) + H_uo*ws_inf[i];

  //  Compute RHS of state variable equations
  double I_fi = -v*H_uv*(u-theta_v[i])*(u_u[i] - u)/tau_fi[i];
  double I_so =  (u-u_o[i])*(1.0-H_uw)/tau_o + H_uw/tau_so;
  double I_si = -H_uw*w*s/tau_si[i];

  f(0) = -(I_fi + I_so + I_si + fext);

  f(1) = (1.0-H_uv)*(v_inf-v)/taum_v - H_uv*v/taup_v[i];

  f(2) = (1.0-H_uw)*(w_inf-w)/taum_w - H_uw*w/taup_w[i];

  f(3) = (0.5*(1.0 + tanh(k_s[i]*(u-u_s[i])))-s)/tau_s;

  RPAR(2) = I_fi;
  RPAR(3) = I_so;
  RPAR(3) = I_si;
}

void CepModBo::getj(const int i, const int n, const Vector<double>& X, Array<double>& JAC)
{

  // Create local variables
  double u = X(0);
  double v = X(1);
  double w = X(2);
  double s = X(3);

  //  Define step functions
  double H_uv  = step(u - theta_v[i]);
  double H_uw  = step(u - theta_w[i]);
  double H_umv = step(u - thetam_v[i]);
  double H_uo  = step(u - theta_o[i]);

  // Define delta functions
  double D_uw = delta(u - theta_w[i]);
  double D_uv = delta(u - theta_v[i]);

  // Define additional constants
  double taum_v = (1.0-H_umv)*taum_v1[i] + H_umv*taum_v2[i];
  double taum_w = taum_w1[i] + 0.50*(taum_w2[i]-taum_w1[i])* (1.0 + tanh(km_w[i]*(u-um_w[i])));
  double tau_so = tau_so1[i] + 0.50*(tau_so2[i]-tau_so1[i])* (1.0+tanh(k_so[i]*(u-u_so[i])));
  double tau_s  = (1.0-H_uw)*tau_s1[i] + H_uw*tau_s2[i];
  double tau_o  = (1.0-H_uo)*tau_o1[i] + H_uo*tau_o2[i];
  double v_inf  = (1.0-H_umv);
  double w_inf  = (1.0-H_uo)*(1.0 - u/tau_winf[i]) + H_uo*ws_inf[i];

  //  Define Jacobian
  JAC = 0.0;

  double n1 = v*H_uv*(u_u[i] + theta_v[i] - 2.0*u)/tau_fi[i];
  double n2 = -(1.0-H_uw)/tau_fi[i];
  double n3 = (-1.0/tau_so + (theta_w[i]-u_o[i])/tau_o + w*s/tau_si[i])*D_uw;

  JAC(0,0) = n1 + n2 + n3;
  JAC(0,1) = H_uv*(u-theta_v[i])*(u_u[i]-u)/tau_fi[i];

  n1 = H_uw/tau_si[i];
  JAC(0,2) = n1*s;
  JAC(0,3) = n1*w;

  n1 = -1.0/taum_v;
  n2 = -1.0/taup_v[i];
  JAC(1,0) = ((v_inf-v)*n1 + v*n2)*D_uv;
  JAC(1,1) = (1.0-H_uv)*n1 + H_uv*n2;

  n1 = -1.0/taum_w;
  n2 = -1.0/taup_w[i];
  JAC(2,0) = ((w_inf-w)*n1 + w*n2)*D_uw;
  JAC(2,2) = (1.0-H_uw)*n1 + H_uw*n2;

  n1 = cosh(k_s[i]*(u-u_s[i]));
  n2 = 1.0/(n1*n1);
  n3 = 1.0/tau_s;
  JAC(3,0) = 0.50*k_s[i]*n2*n3;
  JAC(3,3) = -n3;
}

void CepModBo::init(const int nX, Vector<double> &X)
{
  X(0) = Voffset;
  X(1) = 1.0;
  X(2) = 1.0;
  X(3) = 0.0;
}

void CepModBo::integ_cn2(const int imyo, const int nX, Vector<double>& Xn, const double Ts, const double Ti, 
    const double Istim, const double Ksac, Vector<int>& IPAR, Vector<double>& RPAR)
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
  getf(imyo, nX, Xn, fn, fext, RPAR);

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
    getf(imyo, nX, Xk, fk, fext, RPAR);
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
    getj(imyo, nX, Xk, JAC);

    JAC = Im - 0.50 * dt * JAC;
    JAC = mat_fun::mat_inv(JAC, nX);
    rK = mat_fun::mat_mul(JAC, rK);
    Xk = Xk - rK;
  }

  Xn = Xk;

  getf(imyo, nX, Xn, fn, fext, RPAR);

  Xn(0) = Xn(0)*Vscale + Voffset;

  if (!l2 && !l3) {
    IPAR(1) = IPAR(1) + 1;
  }
}

void CepModBo::integ_fe(const int imyo, const int nX, Vector<double>& X, const double Ts, const double Ti, 
    const double Istim, const double Ksac, Vector<double>& RPAR)
{
  double t = Ts / Tscale;
  double dt = Ti / Tscale;

  double Isac = Ksac * (Vrest - X(0));
  double fext = (Istim + Isac) * Tscale / Vscale;

  X(0) = (X(0) - Voffset) / Vscale;

  Vector<double> f(nX);
  getf(imyo, nX, X, f, fext, RPAR);
  //CALL BO_GETF(imyo, nX, X, f, fext, RPAR);

  X = X + dt*f;
  X(0) = X(0)*Vscale + Voffset;
}

void CepModBo::integ_rk(const int imyo, const int nX, Vector<double>& X, const double Ts, const double Ti, 
    const double Istim, const double Ksac, Vector<double>& RPAR)
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
  getf(imyo, nX, Xrk, frk1, fext, RPAR);

  // RK4: 2nd pass
  Xrk  = X + 0.5 * dt * frk1;
  getf(imyo, nX, Xrk, frk2, fext, RPAR);

  // RK4: 3rd pass
  Xrk  = X + 0.5 * dt * frk2;
  getf(imyo, nX, Xrk, frk3, fext, RPAR);

  // RK4: 4th pass
  Xrk  = X + dt*frk3;
  getf(imyo, nX, Xrk, frk4, fext, RPAR);

  X = X + dt6 * (frk1 + 2.0*(frk2 + frk3) + frk4);

  X(0) = X(0)*Vscale + Voffset;

}

double CepModBo::step(const double r)
{
  double result;

  if (r < 0.0) {
    result = 0.0;
  } else { 
    result = 1.0;
  } 

  return result;
}



