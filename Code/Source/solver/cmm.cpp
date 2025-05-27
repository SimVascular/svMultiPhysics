/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// These subroutines implement the Coupled Momentum Method (CMM).

#include "cmm.h"

#include "all_fun.h"
#include "fluid.h"
#include "lhsa.h"
#include "mat_fun.h"
#include "nn.h"
#include "utils.h"
#include <math.h>

namespace cmm {

void cmm_3d(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx, 
    const Array<double>& al, const Array<double>& yl, const Array<double>& bfl, const Array<double>& Kxi, 
    Array<double>& lR, Array3<double>& lK)
{
  using namespace consts;

  #define n_debug_cmm_3d 
  #ifdef debug_cmm_3d 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "w: " << w ;
  #endif

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  const auto cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  const double ctM = 1.0;
  const double ctC = 36.0;

  double rho = dmn.prop.at(PhysicalProperyType::fluid_density);
  Vector<double> f({dmn.prop.at(PhysicalProperyType::f_x), 
                    dmn.prop.at(PhysicalProperyType::f_y), 
                    dmn.prop.at(PhysicalProperyType::f_z)});

  double T1 = eq.af * eq.gam * dt;
  double amd = eq.am/T1;
  double wl = w*T1;
  double wr = w*rho;
  #ifdef debug_cmm_3d 
  dmsg << "rho: " << rho ;
  dmsg << "T1: " << T1 ;
  dmsg << "amd: " << amd ;
  #endif

  // Indices are not selected based on the equation only
  // because fluid equation always come first
  //
  double p = 0.0;
  Vector<double> u(3), px(3);
  Array<double> ux(3,3);
  auto ud = -f;

  for (int a = 0; a < eNoN; a++) {
    p = p + N(a)*yl(3,a);

    ud(0) = ud(0) + N(a)*(al(0,a)-bfl(0,a));
    ud(1) = ud(1) + N(a)*(al(1,a)-bfl(1,a));
    ud(2) = ud(2) + N(a)*(al(2,a)-bfl(2,a));

    px(0) = px(0) + Nx(0,a)*yl(3,a);
    px(1) = px(1) + Nx(1,a)*yl(3,a);
    px(2) = px(2) + Nx(2,a)*yl(3,a);

    u(0) = u(0) + N(a)*yl(0,a);
    u(1) = u(1) + N(a)*yl(1,a);
    u(2) = u(2) + N(a)*yl(2,a);

    ux(0,0) = ux(0,0) + Nx(0,a)*yl(0,a);
    ux(1,0) = ux(1,0) + Nx(1,a)*yl(0,a);
    ux(2,0) = ux(2,0) + Nx(2,a)*yl(0,a);

    ux(0,1) = ux(0,1) + Nx(0,a)*yl(1,a);
    ux(1,1) = ux(1,1) + Nx(1,a)*yl(1,a);
    ux(2,1) = ux(2,1) + Nx(2,a)*yl(1,a);

    ux(0,2) = ux(0,2) + Nx(0,a)*yl(2,a);
    ux(1,2) = ux(1,2) + Nx(1,a)*yl(2,a);
    ux(2,2) = ux(2,2) + Nx(2,a)*yl(2,a);
  }

  double divU = ux(0,0) + ux(1,1) + ux(2,2);

  if (com_mod.mvMsh) {
    for (int a = 0; a < eNoN; a++) {
      u(0) = u(0) - N(a)*yl(4,a);
      u(1) = u(1) - N(a)*yl(5,a);
      u(2) = u(2) - N(a)*yl(6,a);
    }
  }

  // Strain rate tensor 2*e_ij := (u_ij + u_ji)
  //
  Array<double> es(3,3);
  es(0,0) = ux(0,0) + ux(0,0);
  es(1,0) = ux(1,0) + ux(0,1);
  es(2,0) = ux(2,0) + ux(0,2);

  es(0,1) = es(1,0);
  es(1,1) = ux(1,1) + ux(1,1);
  es(2,1) = ux(2,1) + ux(1,2);

  es(0,2) = es(2,0);
  es(1,2) = es(2,1);
  es(2,2) = ux(2,2) + ux(2,2);

  Array<double> es_x(3,eNoN);
  for (int a = 0; a < eNoN; a++) {
    es_x(0,a) = es(0,0)*Nx(0,a) + es(1,0)*Nx(1,a) + es(2,0)*Nx(2,a);
    es_x(1,a) = es(0,1)*Nx(0,a) + es(1,1)*Nx(1,a) + es(2,1)*Nx(2,a);
    es_x(2,a) = es(0,2)*Nx(0,a) + es(1,2)*Nx(1,a) + es(2,2)*Nx(2,a);
  }

  // Shear-rate := (1*e_ij*e_ij)^.5
  //
  double gam = es(0,0)*es(0,0) + es(1,0)*es(1,0) + es(2,0)*es(2,0) + 
               es(0,1)*es(0,1) + es(1,1)*es(1,1) + es(2,1)*es(2,1) + 
               es(0,2)*es(0,2) + es(1,2)*es(1,2) + es(2,2)*es(2,2);
  gam = sqrt(0.50*gam);

  // Compute viscosity based on shear-rate and chosen viscosity model
  //
  double mu, mu_s, mu_x;
  fluid::get_viscosity(com_mod, dmn, gam, mu, mu_s, mu_x);

  if (utils::is_zero(gam)) {
     mu_x = 0.0;
  } else { 
     mu_x = mu_x / gam;
  }

  double kT = 4.0*pow(ctM/dt,2.0);

  double kU = u(0)*u(0)*Kxi(0,0) + u(1)*u(0)*Kxi(1,0) + u(2)*u(0)*Kxi(2,0) +
              u(0)*u(1)*Kxi(0,1) + u(1)*u(1)*Kxi(1,1) + u(2)*u(1)*Kxi(2,1) + 
              u(0)*u(2)*Kxi(0,2) + u(1)*u(2)*Kxi(1,2) + u(2)*u(2)*Kxi(2,2);

  double kS = Kxi(0,0)*Kxi(0,0) + Kxi(1,0)*Kxi(1,0) + Kxi(2,0)*Kxi(2,0) + 
              Kxi(0,1)*Kxi(0,1) + Kxi(1,1)*Kxi(1,1) + Kxi(2,1)*Kxi(2,1) + 
              Kxi(0,2)*Kxi(0,2) + Kxi(1,2)*Kxi(1,2) + Kxi(2,2)*Kxi(2,2);

  kS = ctC * kS * pow(mu/rho,2.0);
  double tauM = 1.0 / (rho * sqrt( kT + kU + kS ));
  double tauC = 1.0 / (tauM * (Kxi(0,0) + Kxi(1,1) + Kxi(2,2)));

  Vector<double> rV(3);
  rV(0) = ud(0) + u(0)*ux(0,0) + u(1)*ux(1,0) + u(2)*ux(2,0);
  rV(1) = ud(1) + u(0)*ux(0,1) + u(1)*ux(1,1) + u(2)*ux(2,1);
  rV(2) = ud(2) + u(0)*ux(0,2) + u(1)*ux(1,2) + u(2)*ux(2,2);

  Vector<double> up(3);
  up(0) = -tauM*(rho*rV(0) + px(0));
  up(1) = -tauM*(rho*rV(1) + px(1));
  up(2) = -tauM*(rho*rV(2) + px(2));

  double tauB = up(0)*up(0)*Kxi(0,0) + up(1)*up(0)*Kxi(1,0) + 
                up(2)*up(0)*Kxi(2,0) + up(0)*up(1)*Kxi(0,1) + 
                up(1)*up(1)*Kxi(1,1) + up(2)*up(1)*Kxi(2,1) + 
                up(0)*up(2)*Kxi(0,2) + up(1)*up(2)*Kxi(1,2) + 
                up(2)*up(2)*Kxi(2,2);

  if (utils::is_zero(tauB)) {
    tauB = std::numeric_limits<double>::epsilon();
  }
  tauB = rho / sqrt(tauB);

  Vector<double> ua(3);
  ua(0) = u(0) + up(0);
  ua(1) = u(1) + up(1);
  ua(2) = u(2) + up(2);
  double pa = p - tauC*divU;

  rV(0) = tauB*(up(0)*ux(0,0) + up(1)*ux(1,0) + up(2)*ux(2,0));
  rV(1) = tauB*(up(0)*ux(0,1) + up(1)*ux(1,1) + up(2)*ux(2,1));
  rV(2) = tauB*(up(0)*ux(0,2) + up(1)*ux(1,2) + up(2)*ux(2,2));

  Array<double> rM(3,3);
  rM(0,0) = mu*es(0,0) - rho*up(0)*ua(0) + rV(0)*up(0) - pa;
  rM(1,0) = mu*es(1,0) - rho*up(0)*ua(1) + rV(0)*up(1);
  rM(2,0) = mu*es(2,0) - rho*up(0)*ua(2) + rV(0)*up(2);

  rM(0,1) = mu*es(0,1) - rho*up(1)*ua(0) + rV(1)*up(0);
  rM(1,1) = mu*es(1,1) - rho*up(1)*ua(1) + rV(1)*up(1) - pa;
  rM(2,1) = mu*es(2,1) - rho*up(1)*ua(2) + rV(1)*up(2);

  rM(0,2) = mu*es(0,2) - rho*up(2)*ua(0) + rV(2)*up(0);
  rM(1,2) = mu*es(1,2) - rho*up(2)*ua(1) + rV(2)*up(1);
  rM(2,2) = mu*es(2,2) - rho*up(2)*ua(2) + rV(2)*up(2) - pa;

  rV(0) = ud(0) + ua(0)*ux(0,0) + ua(1)*ux(1,0) + ua(2)*ux(2,0);
  rV(1) = ud(1) + ua(0)*ux(0,1) + ua(1)*ux(1,1) + ua(2)*ux(2,1);
  rV(2) = ud(2) + ua(0)*ux(0,2) + ua(1)*ux(1,2) + ua(2)*ux(2,2);

  Vector<double> uNx(eNoN), upNx(eNoN), uaNx(eNoN);

  for (int a = 0; a < eNoN; a++) {
    uNx(a) = u(0)*Nx(0,a)  + u(1)*Nx(1,a)  + u(2)*Nx(2,a);
    upNx(a) = up(0)*Nx(0,a) + up(1)*Nx(1,a) + up(2)*Nx(2,a);
    uaNx(a) = uNx(a) + upNx(a);

    lR(0,a) = lR(0,a) + wr*N(a)*rV(0) + w*(Nx(0,a)*rM(0,0) + Nx(1,a)*rM(1,0) + Nx(2,a)*rM(2,0));
    lR(1,a) = lR(1,a) + wr*N(a)*rV(1) + w*(Nx(0,a)*rM(0,1) + Nx(1,a)*rM(1,1) + Nx(2,a)*rM(2,1));
    lR(2,a) = lR(2,a) + wr*N(a)*rV(2) + w*(Nx(0,a)*rM(0,2) + Nx(1,a)*rM(1,2) + Nx(2,a)*rM(2,2));
    lR(3,a) = lR(3,a) + w*(N(a)*divU - upNx(a));
  }

  for (int a = 0; a < eNoN; a++) {
    for (int b = 0; b < eNoN; b++) {
      rM(0,0) = Nx(0,a)*Nx(0,b);
      rM(1,0) = Nx(1,a)*Nx(0,b);
      rM(2,0) = Nx(2,a)*Nx(0,b);
      rM(0,1) = Nx(0,a)*Nx(1,b);
      rM(1,1) = Nx(1,a)*Nx(1,b);
      rM(2,1) = Nx(2,a)*Nx(1,b);
      rM(0,2) = Nx(0,a)*Nx(2,b);
      rM(1,2) = Nx(1,a)*Nx(2,b);
      rM(2,2) = Nx(2,a)*Nx(2,b);

      double NxNx = Nx(0,a)*Nx(0,b) + Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b);
      double T1 = mu*NxNx + tauB*upNx(a)*upNx(b) + rho*( N(a)*(amd*N(b) + uaNx(b)) + 
                  rho*tauM*uaNx(a)*(uNx(b) + amd*N(b)) );
      double T2 = rho*tauM*uNx(a);
      double T3 = rho*tauM*(amd*N(b) + uNx(b));

      // dM/dU
      //
      lK(0,a,b)  = lK(0,a,b)  + wl*((mu + tauC)*rM(0,0) + T1 + mu_x*es_x(0,a)*es_x(0,b));
      lK(1,a,b)  = lK(1,a,b)  + wl*(mu*rM(1,0) + tauC*rM(0,1) + mu_x*es_x(0,a)*es_x(1,b));
      lK(2,a,b)  = lK(2,a,b)  + wl*(mu*rM(2,0) + tauC*rM(0,2) + mu_x*es_x(0,a)*es_x(2,b));
      lK(4,a,b)  = lK(4,a,b)  + wl*(mu*rM(0,1) + tauC*rM(1,0) + mu_x*es_x(1,a)*es_x(0,b));
      lK(5,a,b)  = lK(5,a,b)  + wl*((mu + tauC)*rM(1,1) + T1 + mu_x*es_x(1,a)*es_x(1,b));
      lK(6,a,b)  = lK(6,a,b)  + wl*(mu*rM(2,1) + tauC*rM(1,2) + mu_x*es_x(1,a)*es_x(2,b));
      lK(8,a,b)  = lK(8,a,b)  + wl*(mu*rM(0,2) + tauC*rM(2,0) + mu_x*es_x(2,a)*es_x(0,b));
      lK(9,a,b) = lK(9,a,b) + wl*(mu*rM(1,2) + tauC*rM(2,1) + mu_x*es_x(2,a)*es_x(1,b));
      lK(10,a,b) = lK(10,a,b) + wl*((mu + tauC)*rM(2,2) + T1 + mu_x*es_x(2,a)*es_x(2,b));

      // dM/dP
      lK(3,a,b)  = lK(3,a,b)  - wl*(Nx(0,a)*N(b) - Nx(0,b)*T2);
      lK(7,a,b)  = lK(7,a,b)  - wl*(Nx(1,a)*N(b) - Nx(1,b)*T2);
      lK(11,a,b) = lK(11,a,b) - wl*(Nx(2,a)*N(b) - Nx(2,b)*T2);

      // dC/dU
      lK(12,a,b) = lK(12,a,b) + wl*(N(a)*Nx(0,b) + Nx(0,a)*T3);
      lK(13,a,b) = lK(13,a,b) + wl*(N(a)*Nx(1,b) + Nx(1,a)*T3);
      lK(14,a,b) = lK(14,a,b) + wl*(N(a)*Nx(2,b) + Nx(2,a)*T3);

      // dC/dP
      lK(15,a,b) = lK(15,a,b) + wl*(tauM*NxNx);
    }
  }
}


void cmm_b(ComMod& com_mod, const faceType& lFa, const int e, const Array<double>& al, const Array<double>& dl, 
    const Array<double>& xl, const Array<double>& bfl, const Vector<double>& pS0l, const Vector<double>& vwp, 
    const Vector<int>& ptr) 
{
  const int nsd  = com_mod.nsd;
  const int dof = com_mod.dof;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];

  Array<double> lR(dof,3);
  Array3<double> lK(dof*dof,3,3);

  // Internal stresses (stiffness) contribution
  Vector<double> pSl(6);
  auto Nx = lFa.Nx.slice(0);
  cmm_stiffness(com_mod, Nx, xl, dl, pS0l, vwp, pSl, lR, lK);

  // Inertia and body forces (mass) contribution
  //
  for (int g = 0; g < lFa.nG; g++) {
    Vector<double> nV(nsd);
    auto Nx = lFa.Nx.slice(g);
    nn::gnnb(com_mod, lFa, e, g, nsd, nsd-1, 3, Nx, nV);
    double Jac = sqrt(utils::norm(nV));
    nV = nV / Jac;
    double w = lFa.w(g)*Jac;
    auto N = lFa.N.col(g);
    cmm_mass(com_mod, w, N, al, bfl, vwp, lR, lK);
  }

  eq.linear_algebra->assemble(com_mod, 3, ptr, lK, lR);
}

void bcmmi(ComMod& com_mod, const int eNoN, const int idof, const double w, const Vector<double>& N, const Array<double>& Nxi, 
    const Array<double>& xl, const Array<double>& tfl, Array<double>& lR)
{
  #define n_debug_bcmmi 
  #ifdef debug_bcmmi 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "idof: " << idof;
  #endif

  // Get traction vector
  Vector<double> tfn(idof);

  for (int a = 0; a < eNoN; a++) {
    tfn = tfn + N(a)*tfl.col(a);
  }

  // Get surface normal vector (reference configuration)
  //
  Array<double> xXi(3,2); 

  for (int a = 0; a < eNoN; a++) {
    xXi.set_col(0, xXi.col(0) + xl.col(a)*Nxi(0,a));
    xXi.set_col(1, xXi.col(1) + xl.col(a)*Nxi(1,a));
  }
  auto nV = utils::cross(xXi);
  #ifdef debug_bcmmi 
  dmsg << "nV: " << nV ;
  #endif

  // Local residual
  //
  if (idof == 1) {
    double wl = w * tfn(0);
    for (int a = 0; a < eNoN; a++) {
      lR(0,a) = lR(0,a) - wl*N(a)*nV(0);
      lR(1,a) = lR(1,a) - wl*N(a)*nV(1);
      lR(2,a) = lR(2,a) - wl*N(a)*nV(2);
    }

  } else {
    double wl = w * sqrt(utils::norm(nV));
    for (int a = 0; a < eNoN; a++) {
      lR(0,a) = lR(0,a) - wl*N(a)*tfn(0);
      lR(1,a) = lR(1,a) - wl*N(a)*tfn(1);
      lR(2,a) = lR(2,a) - wl*N(a)*tfn(2);
    }
  }
}


/// @brief CMM initialization (interior).
///
/// Reproduces Fortran 'CMMI'.
//
void cmmi(ComMod& com_mod, const mshType& lM, const Array<double>& al, const Array<double>& dl, const Array<double>& xl,
    const Array<double>& bfl, const Vector<double>& pS0l, const Vector<double>& vwp, const Vector<int>& ptr)
{
  #define n_debug_cmmi 
  #ifdef debug_cmmi 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  const int dof = com_mod.dof;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  const bool pstEq = com_mod.pstEq; 
  auto& pSn = com_mod.pSn;
  auto& pSa = com_mod.pSa;

  auto Nxi = lM.Nx.slice(0);
  Array<double> xXi(3,2); 

  for (int a = 0; a < 3; a++) {
    for (int i = 0; i < 3; i++) {
      xXi(i,0) = xXi(i,0) + xl(i,a)*Nxi(0,a);
      xXi(i,1) = xXi(i,1) + xl(i,a)*Nxi(1,a);
    }
  }

  auto nV = utils::cross(xXi);
  auto Jac = sqrt(utils::norm(nV));
  nV = nV / Jac;

  Array<double> lR(dof,3); 
  Array3<double> lK(dof*dof,3,3);
  Vector<double> pSl(6);

  // Internal stresses (stiffness) contribution
  cmm_stiffness(com_mod, Nxi, xl, dl, pS0l, vwp, pSl, lR, lK);

  // Inertia and body forces (mass) contribution
  //
  for (int g = 0; g < lM.nG; g++) {
    auto N = lM.N.col(g);
    double w = lM.w(g)*Jac;
    cmm_mass(com_mod, w, N, al, bfl, vwp, lR, lK);

    //  Prestress
    if (pstEq) {
      for (int a = 0; a < 3; a++) {
        int Ac = ptr(a);
        for (int i = 0; i < pSn.nrows(); i++) {
          pSn(i,Ac) = pSn(i,Ac) + w*N(a)*pSl(i);
        }
        pSa(Ac) = pSa(Ac) + w*N(a);
      }
    }
  }

  eq.linear_algebra->assemble(com_mod, 3, ptr, lK, lR);
}

void cmm_mass(ComMod& com_mod, const double w, const Vector<double>& N, const Array<double>& al, 
    const Array<double>& bfl, const Vector<double>& vwp, Array<double>& lR, Array3<double>& lK)
{
  using namespace consts;

  #define n_debug_cmm_mass 
  #ifdef debug_cmm_mass 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  const int nsd  = com_mod.nsd;
  const int dof = com_mod.dof;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  const double dt = com_mod.dt;
  const auto cDmn = com_mod.cDmn;
  #ifdef debug_cmm_mass
  dmsg << "nsd: " << nsd;
  dmsg << "dof: " << dof;
  #endif

  Vector<double> f(3);
  double rho = eq.dmn[cDmn].prop.at(PhysicalProperyType::solid_density);
  f(0) = eq.dmn[cDmn].prop.at(PhysicalProperyType::f_x);
  f(1) = eq.dmn[cDmn].prop.at(PhysicalProperyType::f_y);
  f(2) = eq.dmn[cDmn].prop.at(PhysicalProperyType::f_z);
  #ifdef debug_cmm_mass
  dmsg << "rho: " << rho ;
  dmsg << "f: " << f ;
  dmsg << "cmmVarWall: " << com_mod.cmmVarWall ;
  #endif

  double ht{0.0};

  if (com_mod.cmmVarWall) { 
    ht = vwp(0);
  } else { 
    ht = eq.dmn[cDmn].prop.at(PhysicalProperyType::shell_thickness);
  }

  double wl = w * ht * rho;
  double am = eq.am;
  int i = eq.s;
  int j = i + 1;
  int k = j + 1;
  auto ud = -f;

  for (int a = 0; a < 3; a++) {
    ud(0) = ud(0) + N(a)*(al(i,a)-bfl(0,a));
    ud(1) = ud(1) + N(a)*(al(j,a)-bfl(1,a));
    ud(2) = ud(2) + N(a)*(al(k,a)-bfl(2,a));
  }

  for (int a = 0; a < 3; a++) {
    lR(0,a) = lR(0,a) + wl*N(a)*ud(0);
    lR(1,a) = lR(1,a) + wl*N(a)*ud(1);
    lR(2,a) = lR(2,a) + wl*N(a)*ud(2);
  }

  for (int b = 0; b < 3; b++) {
    for (int a = 0; a < 3; a++) {
      double T1 = wl*am*N(a)*N(b);
      lK(0,a,b) = lK(0,a,b) + T1;
      lK(dof+1,a,b) = lK(dof+1,a,b) + T1;
      lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + T1;
    }
  }
}


void cmm_stiffness(ComMod& com_mod, const Array<double>& Nxi, const Array<double>& xl, const Array<double>& dl,
    const Vector<double>& pS0l, const Vector<double>& vwp, Vector<double>& pSl, Array<double>& lR, Array3<double>& lK)
{
  static const double kT = 5.0 /6.0;

  using namespace consts;

  #define n_debug_cmm_stiffness 
  #ifdef debug_cmm_stiffness 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  const int nsd  = com_mod.nsd;
  const int dof = com_mod.dof;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  const double dt = com_mod.dt;
  const auto cDmn = com_mod.cDmn;

  double nu = eq.dmn[cDmn].prop.at(PhysicalProperyType::poisson_ratio);
  double ht, elM; 

  if (com_mod.cmmVarWall) {
    // thickness
    ht = vwp[0];
    
    // elasticity modulus
    elM = vwp[1];
  } else { 
    ht = eq.dmn[cDmn].prop.at(PhysicalProperyType::shell_thickness);
    elM = eq.dmn[cDmn].prop.at(PhysicalProperyType::elasticity_modulus);
  }

  double lam = elM /(1.0 - nu*nu);
  double mu = 0.50 * elM / (1.0 + nu);

  double afl = eq.af * eq.beta * dt * dt;
  int i = eq.s;
  int j = i + 1;
  int k = j + 1;

  #ifdef debug_cmm_stiffness 
  dmsg << "cmmVarWall: " << com_mod.cmmVarWall;
  dmsg << "nu: " << nu;
  dmsg << "ht: " << ht;
  dmsg << "elM: " << elM;
  dmsg << "lam: " << lam;
  dmsg << "mu: " << mu;
  dmsg << "i: " << i;
  #endif

  Array<double> xXi(3,2);

  for (int a = 0; a < 3; a++) { 
    for (int i = 0; i < 3; i++) { 
      xXi(i,0) = xXi(i,0) + xl(i,a)*Nxi(0,a);
      xXi(i,1) = xXi(i,1) + xl(i,a)*Nxi(1,a);
    }
  }

  auto nV = utils::cross(xXi);
  double Jac = sqrt(utils::norm(nV));
  nV = nV / Jac;

  //  Rotation matrix
  //
  Array<double> thet(3,3);
  thet.set_row(0, xXi.col(0) / sqrt(utils::norm(xXi.col(0))));
  thet.set_row(2, nV);

  thet(1,0) = thet(2,1)*thet(0,2) - thet(2,2)*thet(0,1);
  thet(1,1) = thet(2,2)*thet(0,0) - thet(2,0)*thet(0,2);
  thet(1,2) = thet(2,0)*thet(0,1) - thet(2,1)*thet(0,0);
  auto thetT = mat_fun::transpose(thet);

  // Define phi matrix
  //
  Array<double> phi(9,9);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      phi(i,j) = thet(i,j);
      phi(i+3,j+3) = thet(i,j);
      phi(i+6,j+6) = thet(i,j);
    }
  }

  // Transform the global coordinates into local frame (planar). Copy
  // displacements into a vector.
  //
  Array<double> xloc(2,3);
  Vector<double> ul(9);

  for (int a = 0; a < 3; a++) {
    xloc(0,a) = thet(0,0)*xl(0,a) + thet(0,1)*xl(1,a) + thet(0,2)*xl(2,a);
    xloc(1,a) = thet(1,0)*xl(0,a) + thet(1,1)*xl(1,a) + thet(1,2)*xl(2,a);

    int b = 3 * a;
    ul(b+0) = dl(i,a);
    ul(b+1) = dl(j,a);
    ul(b+2) = dl(k,a);
  }

  // Transformation Jacobian from local element to parent element
  //
  Array<double> xlXi(2,2);

  for (int a = 0; a < 3; a++) {
    for (int i = 0; i < 2; i++) {
      xlXi(i,0) = xlXi(i,0) + xloc(i,a)*Nxi(0,a);
      xlXi(i,1) = xlXi(i,1) + xloc(i,a)*Nxi(1,a);
    }
  }

  // Shape function derivatives in local coordinates
  //
  Array<double> Nxl(2,3);

  for (int a = 0; a < 3; a++) {
    Nxl(0,a) = (Nxi(0,a)*xlXi(1,1) - Nxi(1,a)*xlXi(1,0)) / Jac;
    Nxl(1,a) = (Nxi(1,a)*xlXi(0,0) - Nxi(0,a)*xlXi(0,1)) / Jac;
  }

  // B matrix
  //
  Array<double> Bm(5,9);

  for (int a = 0; a < 3; a++) {
    int b = a*3;
    Bm(0,b+0) = Nxl(0,a);
    Bm(1,b+1) = Nxl(1,a);

    Bm(2,b+0) = Bm(1,b+1);
    Bm(2,b+1) = Bm(0,b+0);
    Bm(3,b+2) = Bm(0,b+0);
    Bm(4,b+2) = Bm(1,b+1);
  }

  // Transform B using phi and compute its transpose
  Bm  = mat_fun::mat_mul(Bm, phi);
  auto BmT = mat_fun::transpose(Bm);

  // Material tensor, D
  Array<double> Dm(5,5);
  Dm(0,0) = lam;
  Dm(0,1) = lam*nu;
  Dm(2,2) = mu;
  Dm(3,3) = mu*kT;

  Dm(1,0) = Dm(0,1);
  Dm(1,1) = Dm(0,0);
  Dm(4,4) = Dm(3,3);

  //  D*Bm
  auto DBm = mat_fun::mat_mul(Dm, Bm);

  // Stress tensor in local frame
  Vector<double> Sl(5);

  for (int a = 0; a < 9; a++) {
    Sl(0) = Sl(0) + DBm(0,a)*ul(a);
    Sl(1) = Sl(1) + DBm(1,a)*ul(a);
    Sl(2) = Sl(2) + DBm(2,a)*ul(a);
    Sl(3) = Sl(3) + DBm(3,a)*ul(a);
    Sl(4) = Sl(4) + DBm(4,a)*ul(a);
  }

  // If prestress is present, convert into full matrix, and transform
  // into local coordinates, convert back into voigt notation
  //
  Array<double> pSm(3,3);
  pSm(0,0) = pS0l(0);
  pSm(1,1) = pS0l(1);
  pSm(2,2) = pS0l(2);
  pSm(0,1) = pS0l(3);
  pSm(0,2) = pS0l(4);
  pSm(1,2) = pS0l(5);

  pSm(1,0) = pSm(0,1);
  pSm(2,0) = pSm(0,2);
  pSm(2,1) = pSm(1,2);

  pSm = mat_fun::mat_mul(pSm, thetT);
  pSm = mat_fun::mat_mul(thet, pSm);

  Vector<double> S0l(5);
  S0l(0) = pSm(0,0);
  S0l(1) = pSm(1,1);
  S0l(2) = pSm(0,1);
  S0l(3) = pSm(0,2);
  S0l(4) = pSm(1,2);

  // Prestress is updated with new stress and pulled to global frame
  //
  pSl = 0.0;

  if (com_mod.pstEq) {
    pSm = 0.0;
    pSm(0,0) = Sl(0);
    pSm(1,1) = Sl(1);
    pSm(0,1) = Sl(2);
    pSm(0,2) = Sl(3);
    pSm(1,2) = Sl(4);

    pSm(1,0) = pSm(0,1);
    pSm(2,0) = pSm(0,2);
    pSm(2,1) = pSm(1,2);

    pSm = mat_fun::mat_mul(pSm,  thet);
    pSm = mat_fun::mat_mul(thetT, pSm);

    pSl(0) = pSm(0,0);
    pSl(1) = pSm(1,1);
    pSl(2) = pSm(2,2);
    pSl(3) = pSm(0,1);
    pSl(4) = pSm(0,2);
    pSl(5) = pSm(1,2);
  }

  // Internal stress contribution to residual together with prestress
  //
  Sl = Sl + S0l;
  Vector<double> BtS(9);

  for (int a = 0; a < 5; a++) {
    BtS(0) = BtS(0) + BmT(0,a)*Sl(a);
    BtS(1) = BtS(1) + BmT(1,a)*Sl(a);
    BtS(2) = BtS(2) + BmT(2,a)*Sl(a);
    BtS(3) = BtS(3) + BmT(3,a)*Sl(a);
    BtS(4) = BtS(4) + BmT(4,a)*Sl(a);
    BtS(5) = BtS(5) + BmT(5,a)*Sl(a);
    BtS(6) = BtS(6) + BmT(6,a)*Sl(a);
    BtS(7) = BtS(7) + BmT(7,a)*Sl(a);
    BtS(8) = BtS(8) + BmT(8,a)*Sl(a);
  }

  // Now all the required tensors are defined, contributions to
  // residual and stiffness can be computed
  //
  double T1 = ht * Jac * 0.50;

  for (int a = 0; a < 3; a++) { 
    int b = a*3;
    lR(0,a) = lR(0,a) + T1*BtS(b+0);
    lR(1,a) = lR(1,a) + T1*BtS(b+1);
    lR(2,a) = lR(2,a) + T1*BtS(b+2);
  }

  // Compute element level global stiffness matrix and remapping
  auto Ke = mat_fun::mat_mul(BmT, DBm);
  T1 = T1*afl;

  for (int a = 0; a < 3; a++) { 
    int i = a*3;

    for (int b = 0; b < 3; b++) { 
      int j = b*3;

      lK(0,a,b) = lK(0,a,b) + T1*Ke(i+0,j+0);
      lK(1,a,b) = lK(1,a,b) + T1*Ke(i+0,j+1);
      lK(2,a,b) = lK(2,a,b) + T1*Ke(i+0,j+2);

      lK(dof+0,a,b) = lK(dof+0,a,b) + T1*Ke(i+1,j+0);
      lK(dof+1,a,b) = lK(dof+1,a,b) + T1*Ke(i+1,j+1);
      lK(dof+2,a,b) = lK(dof+2,a,b) + T1*Ke(i+1,j+2);

      lK(2*dof+0,a,b) = lK(2*dof+0,a,b) + T1*Ke(i+2,j+0);
      lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + T1*Ke(i+2,j+1);
      lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + T1*Ke(i+2,j+2);
    }
  }
}


/// @brief Reproduces Fortran 'CONSTRUCT_CMM'.
//
void construct_cmm(ComMod& com_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Yg, const Array<double>& Dg)
{
  using namespace consts;

  #define n_debug_construct_cmm 
  #ifdef debug_construct_cmm 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  const int nsd  = com_mod.nsd;
  const int tDof = com_mod.tDof;
  const int dof = com_mod.dof;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  auto& cDmn = com_mod.cDmn;
  const int nsymd = com_mod.nsymd;
  auto& pS0 = com_mod.pS0;
  auto& pSn = com_mod.pSn;
  auto& pSa = com_mod.pSa;
  bool pstEq = com_mod.pstEq;

  const bool cmmInit = com_mod.cmmInit;
  const bool cmmVarWall = com_mod.cmmVarWall;
  auto& varWallProps = com_mod.varWallProps;

  int eNoN = lM.eNoN;
  #ifdef debug_construct_cmm 
  dmsg << "cmmVarWall: " << cmmVarWall ;
  dmsg << "cmmInit: " << cmmInit ;
  dmsg << "lM.nEl: " << lM.nEl ;
  dmsg << "eNoN: " << eNoN ;
  #endif

  // CMM: dof = nsd+1
  // CMM init: dof = nsd
  //
  Vector<int> ptr(eNoN);
  Vector<double> pSl(nsymd), N(eNoN);
  Array<double> xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN), vwpl(2,eNoN),
                bfl(nsd,eNoN), pS0l(nsymd,eNoN), Nx(nsd,eNoN), lR(dof,eNoN);
  Array3<double> lK(dof*dof,eNoN,eNoN);


  for (int e = 0; e < lM.nEl; e++) {
    // Change the current domain which will be used in later function calls.
    cDmn = all_fun::domain(com_mod, lM, cEq, e);
    auto cPhys = eq.dmn[cDmn].phys;
    if (cPhys != EquationType::phys_CMM) {
      continue;
    }

    if (cmmInit) {
      if (lM.eType != ElementType::TRI3) {
        throw std::runtime_error("[construct_cmm] CMM initialization is allowed for triangular meshes only");
      }
    } else { 
      if (lM.eType != ElementType::TET4) {
        throw std::runtime_error("[construct_cmm] CMM equation is allowed for tetrahedral meshes only");
      }
    }

    // Create local copies
    pS0l = 0.0;
    vwpl = 0.0;

    for (int a = 0; a < eNoN; a++) {
      int Ac = lM.IEN(a,e);
      ptr(a) = Ac;

      for (int i = 0; i < nsd; i++) {
        xl(i,a) = com_mod.x(i,Ac);
        bfl(i,a) = com_mod.Bf(i,Ac);
      }

      for (int i = 0; i < tDof; i++) {
        al(i,a) = Ag(i,Ac);
        dl(i,a) = Dg(i,Ac);
        yl(i,a) = Yg(i,Ac);
      }

      if (pS0.size() != 0) {
        pS0l.set_col(a, pS0.col(Ac));
      }

      if (cmmVarWall) {
        vwpl.set_col(a, varWallProps.col(Ac));
      }
    }

    if (cmmInit) {
      pSl = 0.0;
      
      // vwp = [thickness, elastic modulus], but averaged over the nodes on the element
      Vector<double> vwp(2);
      
      for (int a = 0; a < eNoN; a++) {
        pSl = pSl + pS0l.col(a);
        vwp = vwp + vwpl.col(a);
      }
      pSl = pSl / static_cast<double>(eNoN);
      vwp = vwp / static_cast<double>(eNoN);
      cmmi(com_mod, lM, al, dl, xl, bfl, pSl, vwp, ptr);

    // Gauss integration
    //
    } else {

      lR = 0.0;
      lK = 0.0;

      double Jac{0.0};
      Array<double> ksix(nsd,nsd);

      for (int g = 0; g < lM.nG; g++) {
        if (g == 0 || !lM.lShpF) {
          auto Nx_g = lM.Nx.slice(g);
          nn::gnn(eNoN, nsd, nsd, Nx_g, xl, Nx, Jac, ksix);
          if (utils::is_zero(Jac)) {
            throw std::runtime_error("[construct_dsolid] Jacobian for element " + std::to_string(e) + " is < 0.");
          }
        }
        double w = lM.w(g) * Jac;
        N = lM.N.col(g);

        cmm_3d(com_mod, eNoN, w, N, Nx, al, yl, bfl, ksix, lR, lK);
      }

      eq.linear_algebra->assemble(com_mod, eNoN, ptr, lK, lR);
    }
  }
}

};
