// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "fluid.h"

#include "all_fun.h"
#include "consts.h"
#include "fs.h"
#include "lhsa.h"
#include "nn.h"
#include "utils.h"
#include "ris.h"

#include <array>
#include <iomanip>
#include <math.h>

namespace fluid {

void b_fluid(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Vector<double>& y, 
    const double h, const Vector<double>& nV, Array<double>& lR, Array3<double>& lK)
{
  using namespace consts;

  #define n_debug_b_fluid
  #ifdef debug_b_fluid
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  const int nsd  = com_mod.nsd;
  const int tDof = com_mod.tDof;
  const int dof = com_mod.dof;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  const int cDmn = com_mod.cDmn;
  const auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  double wl  = w * eq.af * eq.gam * dt;
  double udn = 0.0;
  #ifdef debug_b_fluid
  dmsg << "nsd: " << nsd;
  dmsg << "w: " << w;
  dmsg << "h: " << h;
  dmsg << "nV: " << nV(0) << " " << nV(1) << " " << nV(2);
  dmsg << "wl: " << wl;
  dmsg << "com_mod.mvMsh: " << com_mod.mvMsh;
  #endif 

  Vector<double> u(nsd);

  if (com_mod.mvMsh) {
    for (int i = 0; i < nsd; i++) {
      int j = i + nsd + 1;
      u(i) = y(i) - y(j);
      udn  = udn + u(i)*nV(i);
    }

  } else {
    for (int i = 0; i < nsd; i++) {
      u(i) = y(i);
      udn  = udn + u(i)*nV(i);
    }
  }
  // Compute u dot n for backflow stabilization
  udn = 0.50 * dmn.prop.at(PhysicalProperyType::backflow_stab) * dmn.prop.at(PhysicalProperyType::fluid_density) * (udn - fabs(udn));
  auto hc  = h*nV + udn*u;
  #ifdef debug_b_fluid
  dmsg << "udn: " << udn;
  dmsg << "u: " << u(0) << " " << u(1) << " " << u(2);
  dmsg << "hc: " << hc(0) << " " << hc(1) << " " << hc(2);
  #endif

  // Here the loop is started for constructing left and right hand side
  // Note, if the boundary is a coupled or resistance boundary, the boundary
  // pressure is included in the residual here, but the corresponding tangent
  // contribution is not explicit included in the tangent here. Instead, the 
  // tangent contribution is accounted for by the ADDBCMUL() function within the
  // linear solver
  if (nsd == 2) {
    for (int a = 0; a < eNoN; a++) {
      lR(0,a) = lR(0,a) - w*N(a)*hc(0);
      lR(1,a) = lR(1,a) - w*N(a)*hc(1);

      for (int b = 0; b < eNoN; b++) {
        double T1 = wl*N(a)*N(b)*udn;
        lK(0,a,b) = lK(0,a,b) - T1;
        lK(4,a,b) = lK(4,a,b) - T1;
      }
    }

  } else {
    for (int a = 0; a < eNoN; a++) {
      lR(0,a) = lR(0,a) - w*N(a)*hc(0);
      lR(1,a) = lR(1,a) - w*N(a)*hc(1);
      lR(2,a) = lR(2,a) - w*N(a)*hc(2);

      for (int b = 0; b < eNoN; b++) {
        double T1 = wl*N(a)*N(b)*udn;
        lK(0,a,b)  = lK(0,a,b)  - T1;
        lK(5,a,b)  = lK(5,a,b)  - T1;
        lK(10,a,b) = lK(10,a,b) - T1;
      }
    }
  }
} 


void bw_fluid_2d(ComMod& com_mod, const int eNoNw, const int eNoNq, const double w, const Vector<double>& Nw, 
    const Vector<double>& Nq, const Array<double>& Nwx, const Array<double>& yl, const Vector<double>& ub, 
    const Vector<double>& nV, const Vector<double>& tauB, Array<double>& lR, Array3<double>& lK)
{
  using namespace consts;

  #define n_debug_bw_fluid_2d
  #ifdef debug_bw_fluid_2d 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  double rho = dmn.prop.at(PhysicalProperyType::fluid_density);
  double tauT = tauB(0);
  double tauN = tauB(1);

  double T1 = eq.af * eq.gam * dt;
  double wl = w * T1;

  #ifdef debug_bw_fluid_2d 
  dmsg << "tauT: " << tauT;
  dmsg << "tauN: " << tauN;
  dmsg << "T1: " << T1;
  dmsg << "yl: " << yl;
  dmsg << "ub: " << ub;
  #endif

  Vector<double> u(2), Nxn(eNoNw); 
  Array<double> ux(2,2);

  for (int a = 0; a < eNoNw; a++) {
    u(0) = u(0) + Nw(a)*yl(0,a);
    u(1) = u(1) + Nw(a)*yl(1,a);

    ux(0,0) = ux(0,0) + Nwx(0,a)*yl(0,a);
    ux(1,0) = ux(1,0) + Nwx(1,a)*yl(0,a);
    ux(0,1) = ux(0,1) + Nwx(0,a)*yl(1,a);
    ux(1,1) = ux(1,1) + Nwx(1,a)*yl(1,a);

    Nxn(a) = Nwx(0,a)*nV(0) + Nwx(1,a)*nV(1);
  }

  double p = 0.0;
  for (int a = 0; a < eNoNq; a++) {
    p = p + Nq(a)*yl(2,a);
  }

  Vector<double> uh(2); 

  if (com_mod.mvMsh) {
    for (int a = 0; a < eNoNw; a++) {
      uh(0) = uh(0) + Nw(a)*yl(3,a);
      uh(1) = uh(1) + Nw(a)*yl(4,a);
    }
  }

  double un = (u(0)-uh(0))*nV(0) + (u(1)-uh(1))*nV(1);
  un = (fabs(un) - un) * 0.50;

  u = u - ub;
  double ubn  = u(0)*nV(0) + u(1)*nV(1);

  //  Strain rate tensor 2*e_ij := (u_i,j + u_j,i)
  Array<double> es(2,2);
  es(0,0) = ux(0,0) + ux(0,0);
  es(1,1) = ux(1,1) + ux(1,1);
  es(1,0) = ux(1,0) + ux(0,1);
  es(0,1) = es(1,0);

  // Shear-rate := (2*e_ij*e_ij)^.5
  double gam = es(0,0)*es(0,0) + es(1,0)*es(1,0) + es(0,1)*es(0,1) + es(1,1)*es(1,1);
  gam = sqrt(0.50*gam);

  // Compute viscosity based on shear-rate and chosen viscosity model
  // The returned mu_g := (d\mu / d\gamma)
  double mu, mu_s, mu_g;
  get_viscosity(com_mod, dmn, gam, mu, mu_s, mu_g);

  // sigma.n (deviatoric)
  Vector<double> sgmn(2);
  sgmn(0) = mu*(es(0,0)*nV(0) + es(1,0)*nV(1));
  sgmn(1) = mu*(es(0,1)*nV(0) + es(1,1)*nV(1));

  Vector<double> rV(2);
  rV(0) = p*nV(0) - sgmn(0) + (tauT + rho*un)*u(0) + (tauN-tauT)*ubn*nV(0);
  rV(1) = p*nV(1) - sgmn(1) + (tauT + rho*un)*u(1) + (tauN-tauT)*ubn*nV(1);

  Array<double> rM(2,2);
  rM(0,0) = -mu*( u(0)*nV(0) + u(0)*nV(0) );
  rM(1,0) = -mu*( u(0)*nV(1) + u(1)*nV(0) );
  rM(0,1) = -mu*( u(1)*nV(0) + u(0)*nV(1) );
  rM(1,1) = -mu*( u(1)*nV(1) + u(1)*nV(1) );

  // Local residual (momentum)
  //
  for (int a = 0; a < eNoNw; a++) {
    lR(0,a) = lR(0,a) + w*(Nw(a)*rV(0) + Nwx(0,a)*rM(0,0) + Nwx(1,a)*rM(1,0));
    lR(1,a) = lR(1,a) + w*(Nw(a)*rV(1) + Nwx(0,a)*rM(0,1) + Nwx(1,a)*rM(1,1));
  }

  // Local residual (continuity)
  for (int a = 0; a < eNoNq; a++) {
    lR(2,a) = lR(2,a) - w*Nq(a)*ubn;
  }

  // Tangent (stiffness) matrices
  //
  for (int b = 0; b < eNoNw; b++) {
    for (int a = 0; a < eNoNw; a++) {
      double T3 = Nw(a)*Nw(b);
      double T1 = (tauT + rho*un)*T3 - mu*(Nw(a)*Nxn(b) + Nxn(a)*Nw(b));
      double T2{0};

     // dRm_a1/du_b1
     T2 = (tauN - tauT)*T3*nV(0)*nV(0) - mu*(Nw(a)*Nwx(0,b)*nV(0) + Nw(b)*Nwx(0,a)*nV(0));
     lK(0,a,b) = lK(0,a,b) + wl*(T2 + T1);

     // dRm_a1/du_b2
     T2 = (tauN - tauT)*T3*nV(0)*nV(1) - mu*(Nw(a)*Nwx(0,b)*nV(1) + Nw(b)*Nwx(1,a)*nV(0));
     lK(1,a,b) = lK(1,a,b) + wl*T2;

     // dRm_a2/du_b1
     T2 = (tauN - tauT)*T3*nV(1)*nV(0) - mu*(Nw(a)*Nwx(1,b)*nV(0) + Nw(b)*Nwx(0,a)*nV(1));
     lK(3,a,b) = lK(3,a,b) + wl*T2;

     // dRm_a2/du_b2
     T2 = (tauN - tauT)*T3*nV(1)*nV(1) - mu*(Nw(a)*Nwx(1,b)*nV(1) + Nw(b)*Nwx(1,a)*nV(1));
     lK(4,a,b) = lK(4,a,b)  + wl*(T2 + T1);
    }
  }

  for (int b = 0; b < eNoNq; b++) {
    for (int a = 0; a < eNoNw; a++) {
      double T1 = wl*Nw(a)*Nq(b);

      // dRm_a1/dp_b
      lK(2,a,b) = lK(2,a,b) + T1*nV(0);

      // dRm_a2/dp_b
      lK(5,a,b) = lK(5,a,b) + T1*nV(1);
    }
  }

  for (int b = 0; b < eNoNw; b++) {
    for (int a = 0; a < eNoNq; a++) {
      double T1 = wl*Nq(a)*Nw(b);

      // dRc_a/du_b1
      lK(6,a,b) = lK(6,a,b) - T1*nV(0);

      // dRc_a/du_b2
      lK(7,a,b) = lK(7,a,b) - T1*nV(1);
    }
  }
}


void bw_fluid_3d(ComMod& com_mod, const int eNoNw, const int eNoNq, const double w, const Vector<double>& Nw, 
    const Vector<double>& Nq, const Array<double>& Nwx, const Array<double>& yl, const Vector<double>& ub, 
    const Vector<double>& nV, const Vector<double>& tauB, Array<double>& lR, Array3<double>& lK)
{
  using namespace consts;

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  double rho = dmn.prop.at(PhysicalProperyType::fluid_density);
  double tauT = tauB(0);
  double tauN = tauB(1);

  double T1 = eq.af * eq.gam * dt;
  double wl = w * T1;

  Vector<double> u(3), Nxn(eNoNw); 
  Array<double> ux(3,3);

  for (int a = 0; a < eNoNw; a++) {
    u(0) = u(0) + Nw(a)*yl(0,a);
    u(1) = u(1) + Nw(a)*yl(1,a);
    u(2) = u(2) + Nw(a)*yl(2,a);

    ux(0,0) = ux(0,0) + Nwx(0,a)*yl(0,a);
    // ux(1,0) = ux(1,1) + Nwx(1,a)*yl(0,a);
    ux(1,0) = ux(1,0) + Nwx(1,a)*yl(0,a);
    ux(2,0) = ux(2,0) + Nwx(2,a)*yl(0,a);
    ux(0,1) = ux(0,1) + Nwx(0,a)*yl(1,a);
    ux(1,1) = ux(1,1) + Nwx(1,a)*yl(1,a);
    ux(2,1) = ux(2,1) + Nwx(2,a)*yl(1,a);
    ux(0,2) = ux(0,2) + Nwx(0,a)*yl(2,a);
    ux(1,2) = ux(1,2) + Nwx(1,a)*yl(2,a);
    ux(2,2) = ux(2,2) + Nwx(2,a)*yl(2,a);

    Nxn(a)  = Nwx(0,a)*nV(0) + Nwx(1,a)*nV(1) + Nwx(2,a)*nV(2);
  }

  double p = 0.0;
  for (int a = 0; a < eNoNq; a++) {
    p = p + Nq(a)*yl(3,a);
  }

  Vector<double> uh(3); 

  if (com_mod.mvMsh) {
    for (int a = 0; a < eNoNw; a++) {
      uh(0) = uh(0) + Nw(a)*yl(4,a);
      uh(1) = uh(1) + Nw(a)*yl(5,a);
      uh(2) = uh(2) + Nw(a)*yl(6,a);
    }
  }

  double un = (u(0)-uh(0))*nV(0) + (u(1)-uh(1))*nV(1) + (u(2)-uh(2))*nV(2);
  un = (fabs(un) - un) * 0.50;

  u = u - ub;
  double ubn  = u(0)*nV(0) + u(1)*nV(1) + u(2)*nV(2);

  //  Strain rate tensor 2*e_ij := (u_i,j + u_j,i)
  Array<double> es(3,3);
  es(0,0) = ux(0,0) + ux(0,0);
  es(1,1) = ux(1,1) + ux(1,1);
  es(2,2) = ux(2,2) + ux(2,2);
  es(1,0) = ux(1,0) + ux(0,1);
  es(2,1) = ux(2,1) + ux(1,2);
  es(0,2) = ux(0,2) + ux(2,0);
  es(0,1) = es(1,0);
  es(1,2) = es(2,1);
  es(2,0) = es(0,2);

  // Shear-rate := (2*e_ij*e_ij)^.5
  double gam = es(0,0)*es(0,0) + es(1,0)*es(1,0) + es(2,0)*es(2,0) + 
               es(0,1)*es(0,1) + es(1,1)*es(1,1) + es(2,1)*es(2,1) + 
               es(0,2)*es(0,2) + es(1,2)*es(1,2) + es(2,2)*es(2,2);
  gam = sqrt(0.50*gam);

  // Compute viscosity based on shear-rate and chosen viscosity model
  // The returned mu_g := (d\mu / d\gamma)
  double mu, mu_s, mu_g;
  get_viscosity(com_mod, dmn, gam, mu, mu_s, mu_g);

  // sigma.n (deviatoric)
  Vector<double> sgmn(3);
  sgmn(0) = mu*(es(0,0)*nV(0) + es(1,0)*nV(1) + es(2,0)*nV(2));
  sgmn(1) = mu*(es(0,1)*nV(0) + es(1,1)*nV(1) + es(2,1)*nV(2));
  sgmn(2) = mu*(es(0,2)*nV(0) + es(1,2)*nV(1) + es(2,2)*nV(2));

  Vector<double> rV(3);
  rV(0) = p*nV(0) - sgmn(0) + (tauT + rho*un)*u(0) + (tauN-tauT)*ubn*nV(0);
  rV(1) = p*nV(1) - sgmn(1) + (tauT + rho*un)*u(1) + (tauN-tauT)*ubn*nV(1);
  rV(2) = p*nV(2) - sgmn(2) + (tauT + rho*un)*u(2) + (tauN-tauT)*ubn*nV(2);

  Array<double> rM(3,3);
  rM(0,0) = -mu*( u(0)*nV(0) + u(0)*nV(0) );
  rM(1,0) = -mu*( u(0)*nV(1) + u(1)*nV(0) );
  rM(2,0) = -mu*( u(0)*nV(2) + u(2)*nV(0) );

  rM(0,1) = -mu*( u(1)*nV(0) + u(0)*nV(1) );
  rM(1,1) = -mu*( u(1)*nV(1) + u(1)*nV(1) );
  rM(2,1) = -mu*( u(1)*nV(2) + u(2)*nV(1) );

  rM(0,2) = -mu*( u(2)*nV(0) + u(0)*nV(2) );
  rM(1,2) = -mu*( u(2)*nV(1) + u(1)*nV(2) );
  rM(2,2) = -mu*( u(2)*nV(2) + u(2)*nV(2) );

  // Local residual (momentum)
  //
  for (int a = 0; a < eNoNw; a++) {
    lR(0,a) = lR(0,a) + w*(Nw(a)*rV(0) + Nwx(0,a)*rM(0,0) + Nwx(1,a)*rM(1,0) + Nwx(2,a)*rM(2,0));
    lR(1,a) = lR(1,a) + w*(Nw(a)*rV(1) + Nwx(0,a)*rM(0,1) + Nwx(1,a)*rM(1,1) + Nwx(2,a)*rM(2,1));
    lR(2,a) = lR(2,a) + w*(Nw(a)*rV(2) + Nwx(0,a)*rM(0,2) + Nwx(1,a)*rM(1,2) + Nwx(2,a)*rM(2,2));
  }

  // Local residual (continuity)
  // for (int a = 0; a < eNoNq; a++) {
  //   lR(2,a) = lR(2,a) - w*Nq(a)*ubn;
  // }
  for (int a = 0; a < eNoNq; a++) {
    lR(3,a) = lR(3,a) - w*Nq(a)*ubn;
  }

  // Tangent (stiffness) matrices
  //
  for (int b = 0; b < eNoNw; b++) {
    for (int a = 0; a < eNoNw; a++) {
      double T3 = Nw(a)*Nw(b);
      double T1 = (tauT + rho*un)*T3 - mu*(Nw(a)*Nxn(b) + Nxn(a)*Nw(b));
      double T2{0.0};

      // dRm_a1/du_b1
      T2 = (tauN - tauT)*T3*nV(0)*nV(0) - mu*(Nw(a)*Nwx(0,b)*nV(0) + Nw(b)*Nwx(0,a)*nV(0));
      lK(0,a,b) = lK(0,a,b) + wl*(T2 + T1);

      // dRm_a1/du_b2
      T2 = (tauN - tauT)*T3*nV(0)*nV(1) - mu*(Nw(a)*Nwx(0,b)*nV(1) + Nw(b)*Nwx(1,a)*nV(0));
      lK(1,a,b) = lK(1,a,b) + wl*T2;

      // dRm_a1/du_b3
      T2 = (tauN - tauT)*T3*nV(0)*nV(2) - mu*(Nw(a)*Nwx(0,b)*nV(2) + Nw(b)*Nwx(2,a)*nV(0));
      lK(2,a,b) = lK(2,a,b) + wl*T2;

      // dRm_a2/du_b1
      T2 = (tauN - tauT)*T3*nV(1)*nV(0) - mu*(Nw(a)*Nwx(1,b)*nV(0) + Nw(b)*Nwx(0,a)*nV(1));
      lK(4,a,b) = lK(4,a,b) + wl*T2;

      // dRm_a2/du_b2
      T2 = (tauN - tauT)*T3*nV(1)*nV(1) - mu*(Nw(a)*Nwx(1,b)*nV(1) + Nw(b)*Nwx(1,a)*nV(1));
      lK(5,a,b) = lK(5,a,b)  + wl*(T2 + T1);

      // dRm_a2/du_b3
      T2 = (tauN - tauT)*T3*nV(1)*nV(2) - mu*(Nw(a)*Nwx(1,b)*nV(2) + Nw(b)*Nwx(2,a)*nV(1));
      lK(6,a,b) = lK(6,a,b)  + wl*T2;

      // dRm_a3/du_b1
      T2 = (tauN - tauT)*T3*nV(2)*nV(0) - mu*(Nw(a)*Nwx(2,b)*nV(0) + Nw(b)*Nwx(0,a)*nV(2));
      lK(8,a,b) = lK(8,a,b) + wl*T2;

      // dRm_a3/du_b2
      T2 = (tauN - tauT)*T3*nV(2)*nV(1) - mu*(Nw(a)*Nwx(2,b)*nV(1) + Nw(b)*Nwx(1,a)*nV(2));
      lK(9,a,b) = lK(9,a,b)  + wl*T2;

      // dRm_a3/du_b3
      T2 = (tauN - tauT)*T3*nV(2)*nV(2) - mu*(Nw(a)*Nwx(2,b)*nV(2) + Nw(b)*Nwx(2,a)*nV(2));
      lK(10,a,b) = lK(10,a,b)  + wl*(T2 + T1);
    }
  }

  for (int b = 0; b < eNoNq; b++) {
    for (int a = 0; a < eNoNw; a++) {
      double T1 = wl*Nw(a)*Nq(b);

      // dRm_a1/dp_b
      lK(3,a,b)  = lK(3,a,b)  + T1*nV(0);

      // dRm_a2/dp_b
      lK(7,a,b) = lK(7,a,b)  + T1*nV(1);

      // dRm_a3/dp_b
      lK(11,a,b) = lK(11,a,b)  + T1*nV(2);
    }
  }

  for (int b = 0; b < eNoNw; b++) {
    for (int a = 0; a < eNoNq; a++) {
      double T1 = wl*Nq(a)*Nw(b);

      // dRc_a/du_b1
      lK(12,a,b) = lK(13,a,b) - T1*nV(0);

      // dRc_a/du_b2
      lK(13,a,b) = lK(14,a,b) - T1*nV(1);

      // dRc_a/du_b3
      lK(14,a,b) = lK(15,a,b) - T1*nV(2);
    }
  }
}

/// @brief This is for solving fluid transport equation solving Navier-Stokes
/// equations. Dirichlet boundary conditions are either treated
/// strongly or weakly.
//
void construct_fluid(ComMod& com_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Yg)
{
  #define n_debug_construct_fluid
  #ifdef debug_construct_fluid
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  double start_time = utils::cput();
  #endif

  using namespace consts;

  const int eNoN = lM.eNoN;
  bool vmsStab = false;

  if (lM.nFs == 1) {
     vmsStab = true;
  } else {
     vmsStab = false;
  }

  // l = 3, if nsd==2 ; else 6;
  const int l = com_mod.nsymd;
  const int nsd  = com_mod.nsd;
  const int tDof = com_mod.tDof;
  const int dof = com_mod.dof;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  auto& cDmn = com_mod.cDmn;

  #ifdef debug_construct_fluid
  dmsg << "cEq: " << cEq;
  dmsg << "eq.sym: " << eq.sym;
  dmsg << "eq.dof: " << eq.dof;
  dmsg << "eNoN: " << eNoN;
  dmsg << "vmsStab: " << vmsStab;
  dmsg << "l: " << l;
  dmsg << "tDof: " << tDof;
  dmsg << "dof: " << dof;
  dmsg << "lM.nEl: " <<  lM.nEl;
  dmsg << "nsd: " <<  nsd;
  #endif

  // FLUID: dof = nsd+1
  Vector<int> ptr(eNoN); 
  Array<double> xl(nsd,eNoN); 
  
  // local acceleration vector (for a single element)
  Array<double> al(tDof,eNoN);
  
  // local velocity vector (for a single element)
  Array<double> yl(tDof,eNoN);
  Array<double> bfl(nsd,eNoN);
  
  // local (weak form) residual vector (for a single element) 
  Array<double> lR(dof,eNoN);
  
  // local tangent matrix (for a single element)
  Array3<double> lK(dof*dof,eNoN,eNoN);

  double DDir = 0.0;

  // Loop over all elements of mesh
  //
  int num_c = lM.nEl / 10;

  for (int e = 0; e < lM.nEl; e++) {
    #ifdef debug_construct_fluid
    dmsg << "---------- e: " << e+1;
    #endif
    cDmn = all_fun::domain(com_mod, lM, cEq, e);
    auto cPhys = eq.dmn[cDmn].phys;

    if (cPhys != EquationType::phys_fluid) {
      continue;
    }
    
    double K_inverse_darcy_permeability = eq.dmn[cDmn].prop.at(PhysicalProperyType::inverse_darcy_permeability);

    //  Update shape functions for NURBS
    if (lM.eType == ElementType::NRB) {
      //CALL NRBNNX(lM, e)
    }

    // Create local copies
    for (int a = 0; a < eNoN; a++) {
      int Ac = lM.IEN(a,e);
      ptr(a) = Ac;

      for (int i = 0; i < xl.nrows(); i++) {
        xl(i,a) = com_mod.x(i,Ac);
        bfl(i,a) = com_mod.Bf(i,Ac);
     }
      for (int i = 0; i < al.nrows(); i++) {
        al(i,a) = Ag(i,Ac);
        yl(i,a) = Yg(i,Ac);
      }
    }

    // Initialize residual and tangents
    lR = 0.0;
    lK = 0.0;
    std::array<fsType,2> fs;

    // Set function spaces for velocity and pressure.
    fs::get_thood_fs(com_mod, fs, lM, vmsStab, 1);

    // Define element coordinates appropriate for function spaces
    Array<double> xwl(nsd,fs[0].eNoN); 
    Array<double> Nwx(nsd,fs[0].eNoN); 
    Array<double> Nwxx(l,fs[0].eNoN);

    Array<double> xql(nsd,fs[1].eNoN); 
    Array<double> Nqx(nsd,fs[1].eNoN);

    #ifdef debug_construct_fluid
    dmsg;
    dmsg << "l: " << l;
    dmsg << "fs[0].eNoN: " << fs[0].eNoN;
    dmsg << "fs[1].eNoN: " << fs[1].eNoN;
    #endif

    xwl = xl;

    for (int i = 0; i < xql.nrows(); i++) { 
      for (int j = 0; j < fs[1].eNoN; j++) { 
        xql(i,j) = xl(i,j);
      }
    }

    // Gauss integration 1
    //
    #ifdef debug_construct_fluid
    dmsg;
    dmsg << "Gauss integration 1 ... " << "";
    dmsg << "fs[1].nG: " << fs[0].nG;
    dmsg << "fs[1].lShpF: " << fs[0].lShpF;
    dmsg << "fs[2].nG: " << fs[1].nG;
    dmsg << "fs[2].lShpF: " << fs[1].lShpF;
    #endif

    double Jac{0.0};
    Array<double> ksix(nsd,nsd);

    for (int g = 0; g < fs[0].nG; g++) {
      #ifdef debug_construct_fluid
      dmsg << "===== g: " << g+1;
      #endif
      if (g == 0 || !fs[1].lShpF) {
        auto Nx = fs[1].Nx.rslice(g);
        nn::gnn(fs[1].eNoN, nsd, nsd, Nx, xql, Nqx, Jac, ksix);
        if (utils::is_zero(Jac)) {
           throw std::runtime_error("[construct_fluid] Jacobian for element " + std::to_string(e) + " is < 0.");
        }
      }

      if (g == 0 || !fs[0].lShpF) {
        auto Nx = fs[0].Nx.rslice(g);
        nn::gnn(fs[0].eNoN, nsd, nsd, Nx, xwl, Nwx, Jac, ksix);
        if (utils::is_zero(Jac)) {
           throw std::runtime_error("[construct_fluid] Jacobian for element " + std::to_string(e) + " is < 0.");
        }

        auto Nxx = fs[0].Nxx.rslice(g);
        nn::gn_nxx(l, fs[0].eNoN, nsd, nsd, Nx, Nxx, xwl, Nwx, Nwxx); 
      }

      double w = fs[0].w(g) * Jac;
      #ifdef debug_construct_fluid
      dmsg << "Jac: " << Jac;
      dmsg << "w: " << w;
      #endif

      // Plot the coordinates of the quad point in the current configuration
      if (com_mod.urisFlag) {
        Vector<double> distSrf(com_mod.nUris);
        distSrf = 0.0;
        for (int a = 0; a < eNoN; a++) {
          int Ac = lM.IEN(a,e);
          for (int iUris = 0; iUris < com_mod.nUris; iUris++) {
            distSrf(iUris) += fs[0].N(a,g) * std::fabs(com_mod.uris[iUris].sdf(Ac));
          }
        }

        DDir = 0.0;
        double DDirTmp = 0.0;
        double sdf_deps_temp = 0.0;
        for (int iUris = 0; iUris < com_mod.nUris; iUris++) {
          if (com_mod.uris[iUris].clsFlg) {
            sdf_deps_temp = com_mod.uris[iUris].sdf_deps_close;
          } else {
            sdf_deps_temp = com_mod.uris[iUris].sdf_deps;
          }
          if (distSrf(iUris) <= sdf_deps_temp) {
            DDirTmp = (1 + cos(pi*distSrf(iUris)/sdf_deps_temp))/
                      (2*sdf_deps_temp*sdf_deps_temp);
            if (DDirTmp > DDir) {DDir = DDirTmp;}
          }
        }

        if (!com_mod.urisActFlag) {DDir = 0.0;}
      }

      // Compute momentum residual and tangent matrix.
      //
      if (nsd == 3) {
        auto N0 = fs[0].N.rcol(g); 
        auto N1 = fs[1].N.rcol(g); 
        fluid_3d_m(com_mod, vmsStab, fs[0].eNoN, fs[1].eNoN, w, ksix, N0, N1, 
            Nwx, Nqx, Nwxx, al, yl, bfl, lR, lK, K_inverse_darcy_permeability, DDir);

      } else if (nsd == 2) {
        auto N0 = fs[0].N.rcol(g); 
        auto N1 = fs[1].N.rcol(g); 
        fluid_2d_m(com_mod, vmsStab, fs[0].eNoN, fs[1].eNoN, w, ksix, N0, N1, 
            Nwx, Nqx, Nwxx, al, yl, bfl, lR, lK, K_inverse_darcy_permeability);
      }
    } // g: loop

    // Set function spaces for velocity and pressure.
    //
    fs::get_thood_fs(com_mod, fs, lM, vmsStab, 2);

    // Gauss integration 2
    //
    #ifdef debug_construct_fluid
    dmsg;
    dmsg << "Gauss integration 2 ... " << "";
    dmsg << "fs[1].nG: " << fs[0].nG;
    dmsg << "fs[1].lShpF: " << fs[0].lShpF;
    dmsg << "fs[2].nG: " << fs[1].nG;
    dmsg << "fs[2].lShpF: " << fs[1].lShpF;
    #endif

    for (int g = 0; g < fs[1].nG; g++) {
      if (g == 0 || !fs[0].lShpF) {
        auto Nx = fs[0].Nx.rslice(g);
        nn::gnn(fs[0].eNoN, nsd, nsd, Nx, xwl, Nwx, Jac, ksix);

        if (utils::is_zero(Jac)) {
           throw std::runtime_error("[construct_fluid] Jacobian for element " + std::to_string(e) + " is < 0.");
        }
      }

      if (g == 0 || !fs[1].lShpF) {
        auto Nx = fs[1].Nx.rslice(g);
        nn::gnn(fs[1].eNoN, nsd, nsd, Nx, xql, Nqx, Jac, ksix);

        if (utils::is_zero(Jac)) {
           throw std::runtime_error("[construct_fluid] Jacobian for element " + std::to_string(e) + " is < 0.");
        }
      }
      double w = fs[1].w(g) * Jac;

      // Compute continuity residual and tangent matrix.
      //
      if (nsd == 3) {
        auto N0 = fs[0].N.rcol(g); 
        auto N1 = fs[1].N.rcol(g); 
        fluid_3d_c(com_mod, vmsStab, fs[0].eNoN, fs[1].eNoN, w, ksix, N0, N1, Nwx, Nqx, Nwxx, al, yl, bfl, lR, lK, K_inverse_darcy_permeability, DDir);

      } else if (nsd == 2) {
        auto N0 = fs[0].N.rcol(g); 
        auto N1 = fs[1].N.rcol(g); 
        fluid_2d_c(com_mod, vmsStab, fs[0].eNoN, fs[1].eNoN, w, ksix, N0, N1, Nwx, Nqx, Nwxx, al, yl, bfl, lR, lK, K_inverse_darcy_permeability);
      }

    } // g: loop

    eq.linear_algebra->assemble(com_mod, eNoN, ptr, lK, lR);
    if (com_mod.risFlag) {
      if (!std::all_of(com_mod.ris.clsFlg.begin(), com_mod.ris.clsFlg.end(), [](bool v) { return v; })) {
        ris::doassem_ris(com_mod, eNoN, ptr, lK, lR);
      }
    }

  } // e: loop

  #ifdef debug_construct_fluid
  double end_time = utils::cput();
  double etime = end_time - start_time;
  #endif
}

/// @brief Reproduces Fortran 'FLUID2D_C()'.
//
void fluid_2d_c(ComMod& com_mod, const int vmsFlag, const int eNoNw, const int eNoNq, const double w, 
    const Array<double>& Kxi, const Vector<double>& Nw, const Vector<double>& Nq, const Array<double>& Nwx, 
    const Array<double>& Nqx, const Array<double>& Nwxx, const Array<double>& al, const Array<double>& yl, 
    const Array<double>& bfl, Array<double>& lR, Array3<double>& lK, double K_inverse_darcy_permeability)
{
  using namespace consts;

  #define n_debug_fluid_2d_c
  #ifdef debug_fluid_2d_c
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "vmsFlag: " << vmsFlag;
  dmsg << "eNoNw: " << eNoNw;
  dmsg << "eNoNq: " << eNoNq;
  dmsg << "Nw: " << Nw;
  dmsg << "Nq: " << Nq;
  dmsg << "Nwx: " << Nwx;
  dmsg << "Nqx: " << Nqx;
  dmsg << "w: " << w;
  dmsg << "bfl: " << bfl;
  #endif

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  const double ctM = 1.0;
  const double ctC = 36.0;

  double rho = dmn.prop[PhysicalProperyType::fluid_density];
  Vector<double> f(2);
  // f_x is internal force in x-direction; what is internal force?
  f[0] = dmn.prop[PhysicalProperyType::f_x];
  
  f[1] = dmn.prop[PhysicalProperyType::f_y];

  double T1 = eq.af * eq.gam * dt;
  double amd = eq.am / T1;
  double wl = w*T1;
  #ifdef debug_fluid_2d_c 
  dmsg << "T1: " << T1;
  dmsg << "amd: " << amd;
  dmsg << "wl: " << wl;
  #endif

  // Note that indices are not selected based on the equation because
  // fluid equation always come first
  // Velocity and its gradients, inertia (acceleration & body force)
  //
  // ud[j] = jth component of (acceleration - body force)
  Vector<double> ud{-f[0], -f[1]};
  
  // u[j] = jth component of velocity, u
  Vector<double> u(2);
  
  // ux[i, j] = derivative of jth component of velocity with respect to ith component of x = du_j/dx_i
  Array<double> ux(2,2);
  
  // uxx[i, j, k] = 2nd derivative of jth component of velocity with respect to ith component of x and kth component of x = d2(u_j)/(dx_i*dx_k)
  Array3<double> uxx(2,2,2);

  // eNoNw is number of basis functions for this element, where fluid_2d_c is called for each element individually
  for (int a = 0; a < eNoNw; a++) {
    // a_x - f_x // bfl is body force. why is body force being multiplied by the shape function and summed over all shape functions?
    ud(0) = ud(0) + Nw(a)*(al(0,a)-bfl(0,a));
    
    // a_y - f_y
    ud(1) = ud(1) + Nw(a)*(al(1,a)-bfl(1,a));

    // u_x
    u(0) = u(0) + Nw(a)*yl(0,a);
    
    // u_y
    u(1) = u(1) + Nw(a)*yl(1,a);
    
    // du_x/dx
    ux(0,0) = ux(0,0) + Nwx(0,a)*yl(0,a);
    
    // du_x/dy
    ux(1,0) = ux(1,0) + Nwx(1,a)*yl(0,a);
    
    // du_y/dx
    ux(0,1) = ux(0,1) + Nwx(0,a)*yl(1,a);
    
    // du_y/dy
    ux(1,1) = ux(1,1) + Nwx(1,a)*yl(1,a);

    // d2(u_x)/dx2
    uxx(0,0,0) = uxx(0,0,0) + Nwxx(0,a)*yl(0,a);
    
    // d2(u_x)/dy2
    uxx(1,0,1) = uxx(1,0,1) + Nwxx(1,a)*yl(0,a);
    
    // d2(u_x)/dydx
    uxx(1,0,0) = uxx(1,0,0) + Nwxx(2,a)*yl(0,a);

    // d2(u_y)/dx2
    uxx(0,1,0) = uxx(0,1,0) + Nwxx(0,a)*yl(1,a);
    
    // d2(u_y)/dy2
    uxx(1,1,1) = uxx(1,1,1) + Nwxx(1,a)*yl(1,a);
    
    // d2(u_y)/dydx
    uxx(1,1,0) = uxx(1,1,0) + Nwxx(2,a)*yl(1,a);
  }

  // divergence of velocity
  double divU = ux(0,0) + ux(1,1);
  
  // d2(u_x)/dxdy
  uxx(0,0,1) = uxx(1,0,0);
  
  // d2(u_y)/dxdy
  uxx(0,1,1) = uxx(1,1,0);

  // d2u2[j] = laplacian of jth component of velocity
  Vector<double> d2u2(2);
  d2u2(0) = uxx(0,0,0) + uxx(1,0,1);
  d2u2(1) = uxx(0,1,0) + uxx(1,1,1);


  // Pressure and its gradient
  // px[i] = derivative of pressure with respect to ith component of x = dp/dx_i
  Vector<double> px(2);
  for (int a = 0; a < eNoNq; a++) {
    px(0) = px(0) + Nqx(0,a)*yl(2,a);
    px(1) = px(1) + Nqx(1,a)*yl(2,a);
  }

  //  Update convection velocity relative to mesh velocity
  if (com_mod.mvMsh) {
    for (int a = 0; a < eNoNw; a++) {
      u(0) = u(0) - Nw(a)*yl(3,a);
      u(1) = u(1) - Nw(a)*yl(4,a);
     }
  }

  // 2 * strain rate tensor = 2*e_ij = du_i/dx_j + du_j/dx_i)
  Array<double> es(2,2);
  es(0,0) = ux(0,0) + ux(0,0);
  es(1,0) = ux(1,0) + ux(0,1);
  es(1,1) = ux(1,1) + ux(1,1);
  es(0,1) = es(1,0);

  #ifdef debug_fluid_2d_c 
  dmsg << "ud: " << ud;
  dmsg << "u: " << u;
  dmsg << "ux: " << ux;
  dmsg << "uxx: " << uxx;
  dmsg << "px: " << px;
  dmsg << "es: " << es;
  #endif

  // esNx[i, a] = jth column of (ith row of 2 * strain rate tensor) * derivative of shape function (for element node, a) with respect to x_j, with sum over j = 2 * e_ij * dN_a/dx_j
  Array<double> esNx(2,eNoNw);
  for (int a = 0; a < eNoNw; a++) {
    esNx(0,a) = es(0,0)*Nwx(0,a) + es(1,0)*Nwx(1,a);
    esNx(1,a) = es(0,1)*Nwx(0,a) + es(1,1)*Nwx(1,a);
  }
  
  // es_x[i, j, k] = derivative of 2 * strain rate tensor with respect to kth component of x = d(2 * e_ij)/d(x_k)
  Array3<double> es_x(2,2,2);
  for (int k = 0; k < 2; k++) { 
    // d2(u_0)/(dx_0*dx_k) * 2
    es_x(0,0,k) = uxx(0,0,k) + uxx(0,0,k);
    
    // d2(u_1)/(dx_1*dx_k) * 2
    es_x(1,1,k) = uxx(1,1,k) + uxx(1,1,k);
    
    // d2(u_0)/(dx_1*dx_k) + d2(u_1)/(dx_0*dx_k)
    es_x(1,0,k) = uxx(1,0,k) + uxx(0,1,k);
    
    // d2(u_1)/(dx_0*dx_k) + d2(u_0)/(dx_1*dx_k)
    es_x(0,1,k) = es_x(1,0,k);
  }
  
  // mu_x[j] = gamma * derivative of gamma with respect to jth component of x
  Vector<double> mu_x(2);
  mu_x(0) = (es_x(0,0,0)*es(0,0) + es_x(1,1,0)*es(1,1))*0.50 +  es_x(1,0,0)*es(1,0);
  mu_x(1) = (es_x(0,0,1)*es(0,0) + es_x(1,1,1)*es(1,1))*0.50 +  es_x(1,0,1)*es(1,0);

  // shear rate = gamma = (2*e_ij*e_ij)^0.5
  double gam = es(0,0)*es(0,0) + es(1,0)*es(1,0) + es(0,1)*es(0,1) + es(1,1)*es(1,1);
  gam = sqrt(0.50*gam);

  // Compute viscosity based on shear-rate and chosen viscosity model
  // The returned mu_g := (d\mu / d\gamma)
  // mu_g = derivative of effective dynamic viscosity with respect to gamma
  double mu, mu_s, mu_g;
  get_viscosity(com_mod, dmn, gam, mu, mu_s, mu_g);

  if (utils::is_zero(gam)) {
     mu_g = 0.0;
  } else {
     // mu_g = derivative of effective dynamic viscosity with respect to gamma, then divided by gamma
     mu_g = mu_g / gam;
  }
  
  // mu_x[j] = derivative of effective dynamic viscosity with respect to jth component of x
  std::transform(mu_x.begin(), mu_x.end(), mu_x.begin(), [mu_g](double &v){return mu_g*v;});
  //mu_x(:) = mu_g * mu_x(:)

  #ifdef debug_fluid_2d_c 
  dmsg;
  dmsg << "esNx: " << esNx;
  dmsg << "es_x: " << es_x;
  dmsg << "mu_x: " << mu_x;
  dmsg << "gam: " << gam;
  #endif

  // Stabilization parameters
  //
  // updu[j, i,:] = negative of derivative of ith component of momentum PDE residual (not weak form residual) with respect to jth component of velocity
  Array3<double> updu(2,2,eNoNw);
  
  // up[i] = ith component of u_prime (where u_prime = fine-scale velocity in VMS) = -tau_M / rho * ith component of momentum PDE residual (not weak form residual)
  Vector<double> up(2);
  
  // tauM = tau_M / rho = tau_SUPS / rho
  double tauM{0.0};

  if (vmsFlag) {
    // Stabilization parameters
    double kT = 4.0 * pow(ctM/dt,2.0);
    
    // If we consider the NSB model, we need to add an extra term inside the computation for the stab parameter 
    kT = kT + pow(K_inverse_darcy_permeability*mu/rho, 2.0);  
    
    double kU = u(0)*u(0)*Kxi(0,0) + u(1)*u(0)*Kxi(1,0) + u(0)*u(1)*Kxi(0,1) + u(1)*u(1)*Kxi(1,1);
    double kS = Kxi(0,0)*Kxi(0,0) + Kxi(1,0)*Kxi(1,0) + Kxi(0,1)*Kxi(0,1) + Kxi(1,1)*Kxi(1,1);

    kS = ctC * kS * pow(mu/rho,2.0);
    tauM = 1.0 / (rho * sqrt( kT + kU + kS ));

    // rV[i] = ith component of (acceleration + convective term - body force)
    Vector<double> rV(2);
    rV(0) = ud(0) + u(0)*ux(0,0) + u(1)*ux(1,0);
    rV(1) = ud(1) + u(0)*ux(0,1) + u(1)*ux(1,1);

    // rS[i] = ith component of divergence of (2 * effective dynamic viscosity * strain rate tensor) = d(2 * mu * e_ij)/dx_j
    Vector<double> rS(2);
    rS(0) = mu_x(0)*es(0,0) + mu_x(1)*es(1,0) + mu*d2u2(0);
    rS(1) = mu_x(0)*es(0,1) + mu_x(1)*es(1,1) + mu*d2u2(1);

    up(0) = -tauM*(rho*rV(0) + px(0) - rS(0) + mu*K_inverse_darcy_permeability*u(0));
    up(1) = -tauM*(rho*rV(1) + px(1) - rS(1) + mu*K_inverse_darcy_permeability*u(1));

    for (int a = 0; a < eNoNw; a++) {
      double uNx = u(0)*Nwx(0,a) + u(1)*Nwx(1,a);
      T1 = -rho*uNx + mu*(Nwxx(0,a) + Nwxx(1,a)) + mu_x(0)*Nwx(0,a) + mu_x(1)*Nwx(1,a) - mu*K_inverse_darcy_permeability*Nw(a);

      updu(0,0,a) = mu_x(0)*Nwx(0,a) + d2u2(0)*mu_g*esNx(0,a) + T1;
      updu(1,0,a) = mu_x(1)*Nwx(0,a) + d2u2(1)*mu_g*esNx(0,a);

      updu(0,1,a) = mu_x(0)*Nwx(1,a) + d2u2(0)*mu_g*esNx(1,a);
      updu(1,1,a) = mu_x(1)*Nwx(1,a) + d2u2(1)*mu_g*esNx(1,a) + T1;
    } 

    #ifdef debug_fluid_2d_c 
    dmsg;
    dmsg << "kT: " << kT;
    dmsg << "kS: " << kS;
    dmsg << "tauM: " << tauM;
    dmsg << "rV: " << rV;
    dmsg << "rS: " << rS;
    dmsg << "up: " << up;
    dmsg << "updu: " << updu;
    #endif

  } else {
    tauM = 0.0;
    up = 0.0;
    updu = 0.0;
  }

  // Local residual
  //
  for (int a = 0; a < eNoNq; a++) {
    double upNx = up(0)*Nqx(0,a) + up(1)*Nqx(1,a);
    
    // continuity (weak form) residual
    lR(2,a) = lR(2,a) + w*(Nq(a)*divU - upNx);
  }

  // Tangent (stiffness) matrices
  //
  for (int b = 0; b < eNoNw; b++) {
    T1 = rho*amd*Nw(b);

    for (int a = 0; a < eNoNq; a++) {
      double T2 = Nqx(0,a)*(updu(0,0,b) - T1) + Nqx(1,a)*updu(0,1,b);
      
      // dRc_a/dU_b1
      // derivative of continuity (weak form) residual with respect to the x-component of (the acceleration at the next time step)
      lK(6,a,b) = lK(6,a,b) + wl*(Nq(a)*Nwx(0,b) - tauM*T2);

      T2 = Nqx(0,a)*updu(1,0,b) + Nqx(1,a)*(updu(1,1,b) - T1);
      
      // dRc_a/dU_b2
      // derivative of continuity (weak form) residual with respect to the y-component of (the acceleration at the next time step)
      lK(7,a,b) = lK(7,a,b) + wl*(Nq(a)*Nwx(1,b) - tauM*T2);
    }
  }

  if (vmsFlag) {
    for (int b = 0; b < eNoNq; b++) {
      for (int a = 0; a < eNoNq; a++) {
        double NxNx = Nqx(0,a)*Nqx(0,b) + Nqx(1,a)*Nqx(1,b);
        
        // dRc_a/dP_b
        // derivative of continuity (weak form) residual with respect to (time derivative of pressure at the next time step)
        lK(8,a,b) = lK(8,a,b) + wl*tauM*NxNx;
      }
    }
  }
}


void fluid_2d_m(ComMod& com_mod, const int vmsFlag, const int eNoNw, const int eNoNq, const double w, 
    const Array<double>& Kxi, const Vector<double>& Nw, const Vector<double>& Nq, const Array<double>& Nwx, 
    const Array<double>& Nqx, const Array<double>& Nwxx, const Array<double>& al, const Array<double>& yl, 
    const Array<double>& bfl, Array<double>& lR, Array3<double>& lK, double K_inverse_darcy_permeability)
{
  using namespace consts;

  #define n_debug_fluid_2d_m
  #ifdef debug_fluid_2d_m
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "vmsFlag: " << vmsFlag;
  dmsg << "eNoNw: " << eNoNw;
  dmsg << "eNoNq: " << eNoNq;
  #endif

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  double ctM = 1.0;
  double ctC = 36.0;

  double rho = dmn.prop[PhysicalProperyType::fluid_density];
  Vector<double> f(2);
  
  // f_x is internal force in x-direction; what is internal force?
  f[0] = dmn.prop[PhysicalProperyType::f_x];
  
  f[1] = dmn.prop[PhysicalProperyType::f_y];

  double T1 = eq.af * eq.gam * dt;
  double amd = eq.am / T1;
  double wl = w*T1;
  double wr = w*rho;

  // Note that indices are not selected based on the equation because
  // fluid equation always come first
  // Velocity and its gradients, inertia (acceleration & body force)
  //
  // ud[j] = jth component of (acceleration - body force)
  Vector<double> ud{-f[0], -f[1]};
  
  // u[j] = jth component of velocity, u
  Vector<double> u(2);
  
  // ux[i, j] = derivative of jth component of velocity with respect to ith component of x = du_j/dx_i
  Array<double> ux(2,2);
  
  // uxx[i, j, k] = 2nd derivative of jth component of velocity with respect to ith component of x and kth component of x = d2(u_j)/(dx_i*dx_k)
  Array3<double> uxx(2,2,2);

  // eNoNw is number of basis functions for this element, where fluid_2d_m is called for each element individually
  for (int a = 0; a < eNoNw; a++) {
    
    // a_x - f_x // bfl is body force. why is body force being multiplied by the shape function and summed over all shape functions?
    ud(0) = ud(0) + Nw(a)*(al(0,a)-bfl(0,a));
    
    // a_y - f_y
    ud(1) = ud(1) + Nw(a)*(al(1,a)-bfl(1,a));

    // u_x
    u(0) = u(0) + Nw(a)*yl(0,a);
    
    // u_y
    u(1) = u(1) + Nw(a)*yl(1,a);
    
    // du_x/dx
    ux(0,0) = ux(0,0) + Nwx(0,a)*yl(0,a);
    
    // du_x/dy
    ux(1,0) = ux(1,0) + Nwx(1,a)*yl(0,a);
    
    // du_y/dx
    ux(0,1) = ux(0,1) + Nwx(0,a)*yl(1,a);
    
    // du_y/dy
    ux(1,1) = ux(1,1) + Nwx(1,a)*yl(1,a);
    
    // d2(u_x)/dx2
    uxx(0,0,0) = uxx(0,0,0) + Nwxx(0,a)*yl(0,a);
    
    // d2(u_x)/dy2
    uxx(1,0,1) = uxx(1,0,1) + Nwxx(1,a)*yl(0,a);
    
    // d2(u_x)/dydx
    uxx(1,0,0) = uxx(1,0,0) + Nwxx(2,a)*yl(0,a);
    
    // d2(u_y)/dx2
    uxx(0,1,0) = uxx(0,1,0) + Nwxx(0,a)*yl(1,a);
    
    // d2(u_y)/dy2
    uxx(1,1,1) = uxx(1,1,1) + Nwxx(1,a)*yl(1,a);
    
    // d2(u_y)/dydx
    uxx(1,1,0) = uxx(1,1,0) + Nwxx(2,a)*yl(1,a);
  }

  // divergence of velocity
  double divU = ux(0,0) + ux(1,1);
  
  // d2(u_x)/dxdy
  uxx(0,0,1) = uxx(1,0,0);
  
  // d2(u_y)/dxdy
  uxx(0,1,1) = uxx(1,1,0);
  
  // d2u2[j] = laplacian of jth component of velocity
  Vector<double> d2u2(2);
  d2u2(0) = uxx(0,0,0) + uxx(1,0,1);
  d2u2(1) = uxx(0,1,0) + uxx(1,1,1);

  // Pressure and its gradient
  
  // pressure
  double p = 0.0;
  
  // px[i] = derivative of pressure with respect to ith component of x = dp/dx_i
  Vector<double> px(2);
  
  for (int a = 0; a < eNoNq; a++) {
    p = p + Nq(a)*yl(2,a);
    px(0) = px(0) + Nqx(0,a)*yl(2,a);
    px(1) = px(1) + Nqx(1,a)*yl(2,a);
  }

  //  Update convection velocity relative to mesh velocity
  if (com_mod.mvMsh) {
    for (int a = 0; a < eNoNw; a++) {
      u(0) = u(0) - Nw(a)*yl(3,a);
      u(1) = u(1) - Nw(a)*yl(4,a);
     }
  }

  // 2 * strain rate tensor = 2*e_ij = du_i/dx_j + du_j/dx_i)
  Array<double> es(2,2);
  es(0,0) = ux(0,0) + ux(0,0);
  es(1,0) = ux(1,0) + ux(0,1);
  es(1,1) = ux(1,1) + ux(1,1);
  es(0,1) = es(1,0);

  // esNx[i, a] = jth column of (ith row of 2 * strain rate tensor) * derivative of shape function (for element node, a) with respect to x_j, with sum over j = 2 * e_ij * dN_a/dx_j
  Array<double> esNx(2,eNoNw);
  for (int a = 0; a < eNoNw; a++) {
    esNx(0,a) = es(0,0)*Nwx(0,a) + es(1,0)*Nwx(1,a);
    esNx(1,a) = es(0,1)*Nwx(0,a) + es(1,1)*Nwx(1,a);
  }
  
  // es_x[i, j, k] = derivative of 2 * strain rate tensor with respect to kth component of x = d(2 * e_ij)/d(x_k)
  Array3<double> es_x(2,2,2);
  for (int k = 0; k < 2; k++) { 
    // d2(u_0)/(dx_0*dx_k) * 2
    es_x(0,0,k) = uxx(0,0,k) + uxx(0,0,k);
    
    // d2(u_1)/(dx_1*dx_k) * 2
    es_x(1,1,k) = uxx(1,1,k) + uxx(1,1,k);
    
    // d2(u_0)/(dx_1*dx_k) + d2(u_1)/(dx_0*dx_k)
    es_x(1,0,k) = uxx(1,0,k) + uxx(0,1,k);
    
    // d2(u_1)/(dx_0*dx_k) + d2(u_0)/(dx_1*dx_k)
    es_x(0,1,k) = es_x(1,0,k);
  }
  
  // mu_x[j] = gamma * derivative of gamma with respect to jth component of x
  Vector<double> mu_x(2);
  mu_x(0) = (es_x(0,0,0)*es(0,0) + es_x(1,1,0)*es(1,1))*0.50 +  es_x(1,0,0)*es(1,0);
  mu_x(1) = (es_x(0,0,1)*es(0,0) + es_x(1,1,1)*es(1,1))*0.50 +  es_x(1,0,1)*es(1,0);

  // shear rate = gamma = (2*e_ij*e_ij)^0.5
  double gam = es(0,0)*es(0,0) + es(1,0)*es(1,0) + es(0,1)*es(0,1) + es(1,1)*es(1,1);
  gam = sqrt(0.50*gam);

  // Compute viscosity based on shear-rate and chosen viscosity model
  // The returned mu_g := (d\mu / d\gamma)
  // mu_g = derivative of effective dynamic viscosity with respect to gamma
  double mu, mu_s, mu_g;
  get_viscosity(com_mod, dmn, gam, mu, mu_s, mu_g);

  if (utils::is_zero(gam)) {
     mu_g = 0.0;
  } else {
     // mu_g = derivative of effective dynamic viscosity with respect to gamma, then divided by gamma
     mu_g = mu_g / gam;
  }
  
  // mu_x[j] = derivative of effective dynamic viscosity with respect to jth component of x
  std::transform(mu_x.begin(), mu_x.end(), mu_x.begin(), [mu_g](double &v){return mu_g*v;});
  //mu_x(:) = mu_g * mu_x(:)

  //  Stabilization parameters
  double kT = 4.0 * pow(ctM/dt,2.0);
  
  // If we consider the NSB model, we need to add an extra term inside the computation for the stab parameter 
  kT = kT + pow(K_inverse_darcy_permeability*mu/rho, 2.0);

  double kU = u(0)*u(0)*Kxi(0,0) + u(1)*u(0)*Kxi(1,0) + u(0)*u(1)*Kxi(0,1) + u(1)*u(1)*Kxi(1,1);
  double kS = Kxi(0,0)*Kxi(0,0) + Kxi(1,0)*Kxi(1,0) + Kxi(0,1)*Kxi(0,1) + Kxi(1,1)*Kxi(1,1);
  
  kS = ctC * kS * pow(mu/rho,2.0);
  
  // tauM = tau_M / rho = tau_SUPS / rho
  double tauM = 1.0 / (rho * sqrt( kT + kU + kS ));

  // rV[i] = ith component of (acceleration + convective term - body force)
  Vector<double> rV(2);
  rV(0) = ud(0) + u(0)*ux(0,0) + u(1)*ux(1,0);
  rV(1) = ud(1) + u(0)*ux(0,1) + u(1)*ux(1,1);

  // rS[i] = ith component of divergence of (2 * effective dynamic viscosity * strain rate tensor) = d(2 * mu * e_ij)/dx_j
  Vector<double> rS(2);
  rS(0) = mu_x(0)*es(0,0) + mu_x(1)*es(1,0) + mu*d2u2(0);
  rS(1) = mu_x(0)*es(0,1) + mu_x(1)*es(1,1) + mu*d2u2(1);

  // up[i] = ith component of u_prime (where u_prime = fine-scale velocity in VMS) = -tau_M / rho * ith component of momentum PDE residual (not weak form residual)
  Vector<double> up(2);
  up(0) = -tauM*(rho*rV(0) + px(0) - rS(0) + mu*K_inverse_darcy_permeability*u(0));
  up(1) = -tauM*(rho*rV(1) + px(1) - rS(1) + mu*K_inverse_darcy_permeability*u(1));
  
  // tauC = rho * tau_C; tauB = rho * tau_bar; pa = pressure - rho * tau_C * divergence of velocity
  double tauC, tauB, pa;
  double eps = std::numeric_limits<double>::epsilon();
  
  // ua[i] = ith component of (velocity - tau_M / rho * momentum PDE residual)
  Vector<double> ua(2); 

  if (vmsFlag) {
     tauC = 1.0 / (tauM * (Kxi(0,0) + Kxi(1,1)));
     tauB = up(0)*up(0)*Kxi(0,0) + up(1)*up(0)*Kxi(1,0) + up(0)*up(1)*Kxi(0,1) + up(1)*up(1)*Kxi(1,1);

    if (utils::is_zero(tauB)) {
      tauB = eps;
    }
     tauB = rho / sqrt(tauB);

     ua(0) = u(0) + up(0);
     ua(1) = u(1) + up(1);
     pa = p - tauC*divU;
  } else {
     tauC = 0.0;
     tauB = 0.0;
     ua = u;
     pa = p;
  }
  
  // rV[i] = -tau_bar * tau_M * jth component of momentum PDE residual * du_i/dx_j
  rV(0) = tauB*(up(0)*ux(0,0) + up(1)*ux(1,0));
  rV(1) = tauB*(up(0)*ux(0,1) + up(1)*ux(1,1));
  
  // rM[k, i] = r_Mk * r_Mj * du_i/dx_j + (u_k * r_Mi + r_Mk * r_Mi) + (mu*es(0,0) + pa)
  Array<double> rM(2,2);
  rM(0,0) = mu*es(0,0) - rho*up(0)*ua(0) + rV(0)*up(0) - pa;
  rM(1,0) = mu*es(1,0) - rho*up(0)*ua(1) + rV(0)*up(1);

  rM(0,1) = mu*es(0,1) - rho*up(1)*ua(0) + rV(1)*up(0);
  rM(1,1) = mu*es(1,1) - rho*up(1)*ua(1) + rV(1)*up(1) - pa;

  rV(0) = ud(0) + ua(0)*ux(0,0) + ua(1)*ux(1,0);
  rV(1) = ud(1) + ua(0)*ux(0,1) + ua(1)*ux(1,1);

  // Local residual
  
  // updu[j, i,:] = negative of derivative of ith component of momentum PDE residual (not weak form residual) with respect to jth component of velocity
  Array3<double> updu(2,2,eNoNw);
  
  // uNx = u_i * dN_a/dx_i; upNx = -tau_M / rho * ith component of momentum PDE residual * dN_a/dx_i; uaNx = -tau_M / rho * (u_i + ith component of momentum PDE residual) * dN_a/dx_i
  Vector<double> uNx(eNoNw), upNx(eNoNw), uaNx(eNoNw);

  for (int a = 0; a < eNoNw; a++) {
    lR(0,a) = lR(0,a) + wr*Nw(a)*rV(0) + w*(Nwx(0,a)*rM(0,0) + Nwx(1,a)*rM(1,0));
    lR(1,a) = lR(1,a) + wr*Nw(a)*rV(1) + w*(Nwx(0,a)*rM(0,1) + Nwx(1,a)*rM(1,1));

    // Quantities used for stiffness matrix
    uNx(a) = u(0)*Nwx(0,a) + u(1)*Nwx(1,a);
    upNx(a) = up(0)*Nwx(0,a) + up(1)*Nwx(1,a);

    if (vmsFlag) {
      uaNx(a) = uNx(a) + upNx(a);
    } else {
      uaNx(a) = uNx(a);
    }

    T1 = -rho*uNx(a) + mu*(Nwxx(0,a) + Nwxx(1,a)) + mu_x(0)*Nwx(0,a) + mu_x(1)*Nwx(1,a) - mu*K_inverse_darcy_permeability*Nw(a);

    updu(0,0,a) = mu_x(0)*Nwx(0,a) + d2u2(0)*mu_g*esNx(0,a) + T1;
    updu(1,0,a) = mu_x(1)*Nwx(0,a) + d2u2(1)*mu_g*esNx(0,a);

    updu(0,1,a) = mu_x(0)*Nwx(1,a) + d2u2(0)*mu_g*esNx(1,a);
    updu(1,1,a) = mu_x(1)*Nwx(1,a) + d2u2(1)*mu_g*esNx(1,a) + T1;
  }

  // Tangent (stiffness) matrices
  for (int b = 0; b < eNoNw; b++) {
    for (int a = 0; a < eNoNw; a++) {
      rM(0,0) = Nwx(0,a)*Nwx(0,b);
      rM(1,0) = Nwx(1,a)*Nwx(0,b);
      rM(0,1) = Nwx(0,a)*Nwx(1,b);
      rM(1,1) = Nwx(1,a)*Nwx(1,b);

      double NxNx = Nwx(0,a)*Nwx(0,b) + Nwx(1,a)*Nwx(1,b);
      T1 = mu*NxNx + rho*amd*Nw(b)*(Nw(a) + rho*tauM*uaNx(a)) + rho*Nw(a)*(uNx(b)+upNx(b)) + tauB*upNx(a)*upNx(b);

      double T2 = (mu + tauC)*rM(0,0) + esNx(0,a)*mu_g*esNx(0,b) - rho*tauM*uaNx(a)*updu(0,0,b);
      lK(0,a,b) = lK(0,a,b)  + wl*(T2 + T1);
      
      // dRm_a1/du_b1
      // derivative of x-component of momentum (weak form) residual with respect to the x-component of (the acceleration at the next time step)
      lK(0,a,b) = lK(0,a,b)  + mu*K_inverse_darcy_permeability*wl*Nw(b)*Nw(a);

      T2 = mu*rM(1,0) + tauC*rM(0,1) + esNx(0,a)*mu_g*esNx(1,b) - rho*tauM*uaNx(a)*updu(1,0,b);
      
      // dRm_a1/du_b2
      // derivative of x-component of momentum (weak form) residual with respect to the y-component of (the acceleration at the next time step)
      lK(1,a,b) = lK(1,a,b) + wl*(T2);

      T2 = mu*rM(0,1) + tauC*rM(1,0) + esNx(1,a)*mu_g*esNx(0,b) - rho*tauM*uaNx(a)*updu(0,1,b);
      
      // dRm_a2/du_b1
      // derivative of y-component of momentum (weak form) residual with respect to the x-component of (the acceleration at the next time step)
      lK(3,a,b) = lK(3,a,b) + wl*(T2);

      T2 = (mu + tauC)*rM(1,1) + esNx(1,a)*mu_g*esNx(1,b) - rho*tauM*uaNx(a)*updu(1,1,b);
      
      // dRm_a2/du_b2
      // derivative of y-component of momentum (weak form) residual with respect to the y-component of (the acceleration at the next time step)
      lK(4,a,b) = lK(4,a,b) + wl*(T2 + T1);
      lK(4,a,b) = lK(4,a,b)  + mu*K_inverse_darcy_permeability*wl*Nw(b)*Nw(a);
    }
  }

  for (int b = 0; b < eNoNq; b++) {
    for (int a = 0; a < eNoNw; a++) {
      T1 = rho*tauM*uaNx(a);

      // dRm_a1/dp_b
      // derivative of x-component of momentum (weak form) residual with respect to (time derivative of pressure at the next time step)
      lK(2,a,b)  = lK(2,a,b)  - wl*(Nwx(0,a)*Nq(b) - Nqx(0,b)*T1);

      // dRm_a2/dp_b
      // derivative of y-component of momentum (weak form) residual with respect to (time derivative of pressure at the next time step)
      lK(5,a,b) = lK(5,a,b)  - wl*(Nwx(1,a)*Nq(b) - Nqx(1,b)*T1);
    }
  }
  
  // Residual contribution Birkman term 
  // Local residue
  for (int a = 0; a < eNoNw; a++) {
      lR(0,a) = lR(0,a) + mu*K_inverse_darcy_permeability*w*Nw(a)*(u(0)+up(0));
      lR(1,a) = lR(1,a) + mu*K_inverse_darcy_permeability*w*Nw(a)*(u(1)+up(1));
  }
}


/// @brief Element continuity residual.
//
void fluid_3d_c(ComMod& com_mod, const int vmsFlag, const int eNoNw, const int eNoNq, const double w, 
    const Array<double>& Kxi, const Vector<double>& Nw, const Vector<double>& Nq, const Array<double>& Nwx, 
    const Array<double>& Nqx, const Array<double>& Nwxx, const Array<double>& al, const Array<double>& yl, 
    const Array<double>& bfl, Array<double>& lR, Array3<double>& lK, double K_inverse_darcy_permeability, double DDir)
{
  #define n_debug_fluid3d_c
  #ifdef debug_fluid3d_c
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "vmsFlag: " << vmsFlag;
  dmsg << "eNoNw: " << eNoNw;
  dmsg << "eNoNq: " << eNoNq;
  double start_time = utils::cput();
  #endif
  
  // Maximum size of arrays sized by (3,eNoNw) -> (3,MAX_SIZE).
  const int MAX_SIZE = 27;

  using namespace consts;

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  const double ctM  = 1.0;
  const double ctC  = 36.0;

  double Res;
  if (!com_mod.urisFlag) {
    Res = 0.0;
  } else {
    Res = com_mod.urisRes;
    if (com_mod.uris[0].clsFlg) {Res = com_mod.urisResClose;}
  }

  double rho = dmn.prop[PhysicalProperyType::fluid_density];
  double f[3];
  f[0] = dmn.prop[PhysicalProperyType::f_x];
  f[1] = dmn.prop[PhysicalProperyType::f_y];
  f[2] = dmn.prop[PhysicalProperyType::f_z];

  double T1 = eq.af * eq.gam * dt;
  double amd = eq.am / T1;
  double wl = w*T1;
  double wr = w*rho;

  // Note that indices are not selected based on the equation because
  // fluid equation always come first
  // Velocity and its gradients, inertia (acceleration & body force)
  //
  double ud[3] = {-f[0], -f[1], -f[2]};
  double u[3] = {};
  double ux[3][3] = {};
  double uxx[3][3][3] = {};

  for (int a = 0; a < eNoNw; a++) {
    ud[0] = ud[0] + Nw(a)*(al(0,a)-bfl(0,a));
    ud[1] = ud[1] + Nw(a)*(al(1,a)-bfl(1,a));
    ud[2] = ud[2] + Nw(a)*(al(2,a)-bfl(2,a));

    u[0] = u[0] + Nw(a)*yl(0,a);
    u[1] = u[1] + Nw(a)*yl(1,a);
    u[2] = u[2] + Nw(a)*yl(2,a);

    ux[0][0] = ux[0][0] + Nwx(0,a)*yl(0,a);
    ux[1][0] = ux[1][0] + Nwx(1,a)*yl(0,a);
    ux[2][0] = ux[2][0] + Nwx(2,a)*yl(0,a);
    ux[0][1] = ux[0][1] + Nwx(0,a)*yl(1,a);
    ux[1][1] = ux[1][1] + Nwx(1,a)*yl(1,a);
    ux[2][1] = ux[2][1] + Nwx(2,a)*yl(1,a);
    ux[0][2] = ux[0][2] + Nwx(0,a)*yl(2,a);
    ux[1][2] = ux[1][2] + Nwx(1,a)*yl(2,a);
    ux[2][2] = ux[2][2] + Nwx(2,a)*yl(2,a);

    uxx[0][0][0] += Nwxx(0,a)*yl(0,a);
    uxx[1][0][1] += Nwxx(1,a)*yl(0,a);
    uxx[2][0][2] += Nwxx(2,a)*yl(0,a);
    uxx[1][0][0] += Nwxx(3,a)*yl(0,a);
    uxx[2][0][1] += Nwxx(4,a)*yl(0,a);
    uxx[0][0][2] += Nwxx(5,a)*yl(0,a);

    uxx[0][1][0] += Nwxx(0,a)*yl(1,a);
    uxx[1][1][1] += Nwxx(1,a)*yl(1,a);
    uxx[2][1][2] += Nwxx(2,a)*yl(1,a);
    uxx[1][1][0] += Nwxx(3,a)*yl(1,a);
    uxx[2][1][1] += Nwxx(4,a)*yl(1,a);
    uxx[0][1][2] += Nwxx(5,a)*yl(1,a);

    uxx[0][2][0] += Nwxx(0,a)*yl(2,a);
    uxx[1][2][1] += Nwxx(1,a)*yl(2,a);
    uxx[2][2][2] += Nwxx(2,a)*yl(2,a);
    uxx[1][2][0] += Nwxx(3,a)*yl(2,a);
    uxx[2][2][1] += Nwxx(4,a)*yl(2,a);
    uxx[0][2][2] += Nwxx(5,a)*yl(2,a);
  }

  double divU = ux[0][0] + ux[1][1] + ux[2][2];

  uxx[0][0][1] = uxx[1][0][0];
  uxx[1][0][2] = uxx[2][0][1];
  uxx[2][0][0] = uxx[0][0][2];

  uxx[0][1][1] = uxx[1][1][0];
  uxx[1][1][2] = uxx[2][1][1];
  uxx[2][1][0] = uxx[0][1][2];

  uxx[0][2][1] = uxx[1][2][0];
  uxx[1][2][2] = uxx[2][2][1];
  uxx[2][2][0] = uxx[0][2][2];

  double d2u2[3] = {};
  d2u2[0] = uxx[0][0][0] + uxx[1][0][1] + uxx[2][0][2];
  d2u2[1] = uxx[0][1][0] + uxx[1][1][1] + uxx[2][1][2];
  d2u2[2] = uxx[0][2][0] + uxx[1][2][1] + uxx[2][2][2];

  // Pressure and its gradient
  //
  double px[3] = {};

  for (int a = 0; a < eNoNq; a++) {
    px[0] = px[0] + Nqx(0,a)*yl(3,a);
    px[1] = px[1] + Nqx(1,a)*yl(3,a);
    px[2] = px[2] + Nqx(2,a)*yl(3,a);
  }

  // Update convection velocity relative to mesh velocity
  //
  if (com_mod.mvMsh) {
    for (int a = 0; a < eNoNw; a++) {
      u[0] = u[0] - Nw(a)*yl(4,a);
      u[1] = u[1] - Nw(a)*yl(5,a);
      u[2] = u[2] - Nw(a)*yl(6,a);
    }
  }

  // Strain rate tensor 2*e_ij := (u_i,j + u_j,i)
  //
  double es[3][3] = {};
  es[0][0] = ux[0][0] + ux[0][0];
  es[1][1] = ux[1][1] + ux[1][1];
  es[2][2] = ux[2][2] + ux[2][2];
  es[1][0] = ux[1][0] + ux[0][1];
  es[2][1] = ux[2][1] + ux[1][2];
  es[0][2] = ux[0][2] + ux[2][0];
  es[0][1] = es[1][0];
  es[1][2] = es[2][1];
  es[2][0] = es[0][2];

  double esNx[3][MAX_SIZE];

  for (int a = 0; a < eNoNw; a++) {
    esNx[0][a] = es[0][0]*Nwx(0,a) + es[1][0]*Nwx(1,a) + es[2][0]*Nwx(2,a);
    esNx[1][a] = es[0][1]*Nwx(0,a) + es[1][1]*Nwx(1,a) + es[2][1]*Nwx(2,a);
    esNx[2][a] = es[0][2]*Nwx(0,a) + es[1][2]*Nwx(1,a) + es[2][2]*Nwx(2,a);
  }

  double es_x[3][3][3] = {};

  for (int k = 0; k < 3; k++) { 
     es_x[0][0][k] = uxx[0][0][k] + uxx[0][0][k];
     es_x[1][1][k] = uxx[1][1][k] + uxx[1][1][k];
     es_x[2][2][k] = uxx[2][2][k] + uxx[2][2][k];
     es_x[1][0][k] = uxx[1][0][k] + uxx[0][1][k];
     es_x[2][1][k] = uxx[2][1][k] + uxx[1][2][k];
     es_x[0][2][k] = uxx[0][2][k] + uxx[2][0][k];

     es_x[0][1][k] = es_x[1][0][k];
     es_x[1][2][k] = es_x[2][1][k];
     es_x[2][0][k] = es_x[0][2][k];
  }

  double mu_x[3];

  mu_x[0] = (es_x[0][0][0]*es[0][0] + es_x[1][1][0]*es[1][1]
          +  es_x[2][2][0]*es[2][2])*0.5
          +  es_x[1][0][0]*es[1][0] + es_x[2][1][0]*es[2][1]
          +  es_x[0][2][0]*es[0][2];

  mu_x[1] = (es_x[0][0][1]*es[0][0] + es_x[1][1][1]*es[1][1]
          +  es_x[2][2][1]*es[2][2])*0.5
          +  es_x[1][0][1]*es[1][0] + es_x[2][1][1]*es[2][1]
          +  es_x[0][2][1]*es[0][2];

  mu_x[2] = (es_x[0][0][2]*es[0][0] + es_x[1][1][2]*es[1][1]
          +  es_x[2][2][2]*es[2][2])*0.5
          +  es_x[1][0][2]*es[1][0] + es_x[2][1][2]*es[2][1]
          +  es_x[0][2][2]*es[0][2];

  // Shear-rate := (2*e_ij*e_ij)^.5
  double gam = es[0][0]*es[0][0] + es[1][0]*es[1][0] + es[2][0]*es[2][0]
             + es[0][1]*es[0][1] + es[1][1]*es[1][1] + es[2][1]*es[2][1]
             + es[0][2]*es[0][2] + es[1][2]*es[1][2] + es[2][2]*es[2][2];
  gam = sqrt(0.5*gam);

  // Compute viscosity based on shear-rate and chosen viscosity model
  // The returned mu_g := (d\mu / d\gamma)
  double mu, mu_s, mu_g;
  get_viscosity(com_mod, dmn, gam, mu, mu_s, mu_g);

  if (utils::is_zero(gam)) {
     mu_g = 0.0;
  } else {
     mu_g = mu_g / gam;
  }

  for (int i = 0; i < 3; i++) {
    mu_x[i] = mu_g * mu_x[i];
  }

  // Stabilization parameters
  //
  double up[3] = {};
  double updu[3][3][MAX_SIZE] = {};
  double tauM = 0.0;

  if (vmsFlag) {
    // Stabilization parameters
    double kT = 4.0 * pow(ctM/dt,2.0);
    
    // If we consider the NSB model, we need to add an extra term inside the computation for the stab parameter 
    kT = kT + pow(K_inverse_darcy_permeability*mu/rho, 2.0); 
    
    // In case of unfitted RIS, compute the delta function at the quad point,
    // add the additional value to the stabilization param 
    kT = kT + pow(Res*DDir, 2.0);

    double kU = u[0]*u[0]*Kxi(0,0) + u[1]*u[0]*Kxi(1,0) + u[2]*u[0]*Kxi(2,0)
              + u[0]*u[1]*Kxi(0,1) + u[1]*u[1]*Kxi(1,1) + u[2]*u[1]*Kxi(2,1)
              + u[0]*u[2]*Kxi(0,2) + u[1]*u[2]*Kxi(1,2) + u[2]*u[2]*Kxi(2,2);

    double kS = Kxi(0,0)*Kxi(0,0) + Kxi(1,0)*Kxi(1,0) + Kxi(2,0)*Kxi(2,0)
              + Kxi(0,1)*Kxi(0,1) + Kxi(1,1)*Kxi(1,1) + Kxi(2,1)*Kxi(2,1)
              + Kxi(0,2)*Kxi(0,2) + Kxi(1,2)*Kxi(1,2) + Kxi(2,2)*Kxi(2,2);

    kS = ctC * kS * pow(mu/rho,2.0);
    tauM = 1.0 / (rho * sqrt( kT + kU + kS ));

    double rV[3];
    rV[0] = ud[0] + u[0]*ux[0][0] + u[1]*ux[1][0] + u[2]*ux[2][0];
    rV[1] = ud[1] + u[0]*ux[0][1] + u[1]*ux[1][1] + u[2]*ux[2][1];
    rV[2] = ud[2] + u[0]*ux[0][2] + u[1]*ux[1][2] + u[2]*ux[2][2];

    double rS[3];
    rS[0] = mu_x[0]*es[0][0] + mu_x[1]*es[1][0] + mu_x[2]*es[2][0] + mu*d2u2[0];
    rS[1] = mu_x[0]*es[0][1] + mu_x[1]*es[1][1] + mu_x[2]*es[2][1] + mu*d2u2[1];
    rS[2] = mu_x[0]*es[0][2] + mu_x[1]*es[1][2] + mu_x[2]*es[2][2] + mu*d2u2[2];

    // up[0] = -tauM*(rho*rV[0] + px[0] - rS[0] + mu*K_inverse_darcy_permeability*u[0]);
    // up[1] = -tauM*(rho*rV[1] + px[1] - rS[1] + mu*K_inverse_darcy_permeability*u[1]);
    // up[2] = -tauM*(rho*rV[2] + px[2] - rS[2] + mu*K_inverse_darcy_permeability*u[2]);

    up[0] = -tauM*(rho*rV[0] + px[0] - rS[0] + mu*K_inverse_darcy_permeability*u[0]
                   + (Res*DDir)*u[0]);
    up[1] = -tauM*(rho*rV[1] + px[1] - rS[1] + mu*K_inverse_darcy_permeability*u[1]
                   + (Res*DDir)*u[1]);
    up[2] = -tauM*(rho*rV[2] + px[2] - rS[2] + mu*K_inverse_darcy_permeability*u[2]
                   + (Res*DDir)*u[2]);

    for (int a = 0; a < eNoNw; a++) {
      double uNx = u[0]*Nwx(0,a) + u[1]*Nwx(1,a) + u[2]*Nwx(2,a);
      // T1 = -rho*uNx + mu*(Nwxx(0,a) + Nwxx(1,a) + Nwxx(2,a)) + mu_x[0]*Nwx(0,a) + mu_x[1]*Nwx(1,a) + mu_x[2]*Nwx(2,a) - mu*K_inverse_darcy_permeability*Nw(a);

      T1 = -rho*uNx + mu*(Nwxx(0,a) + Nwxx(1,a) + Nwxx(2,a)) 
           + mu_x[0]*Nwx(0,a) + mu_x[1]*Nwx(1,a) + mu_x[2]*Nwx(2,a) 
           - mu*K_inverse_darcy_permeability*Nw(a)
           - (Res*DDir)*Nw(a);

      updu[0][0][a] = mu_x[0]*Nwx(0,a) + d2u2[0]*mu_g*esNx[0][a] + T1;
      updu[1][0][a] = mu_x[1]*Nwx(0,a) + d2u2[1]*mu_g*esNx[0][a];
      updu[2][0][a] = mu_x[2]*Nwx(0,a) + d2u2[2]*mu_g*esNx[0][a];
  
      updu[0][1][a] = mu_x[0]*Nwx(1,a) + d2u2[0]*mu_g*esNx[1][a];
      updu[1][1][a] = mu_x[1]*Nwx(1,a) + d2u2[1]*mu_g*esNx[1][a] + T1;
      updu[2][1][a] = mu_x[2]*Nwx(1,a) + d2u2[2]*mu_g*esNx[1][a];
  
      updu[0][2][a] = mu_x[0]*Nwx(2,a) + d2u2[0]*mu_g*esNx[2][a];
      updu[1][2][a] = mu_x[1]*Nwx(2,a) + d2u2[1]*mu_g*esNx[2][a];
      updu[2][2][a] = mu_x[2]*Nwx(2,a) + d2u2[2]*mu_g*esNx[2][a] + T1;
    }

  } else {
    tauM = 0.0;
    std::memset(up, 0, sizeof up);
    std::memset(updu, 0, sizeof updu);
  }

  //  Local residual
  //
  for (int a = 0; a < eNoNq; a++) {
    double upNx = up[0]*Nqx(0,a) + up[1]*Nqx(1,a) + up[2]*Nqx(2,a);
    lR(3,a) = lR(3,a) + w*(Nq(a)*divU - upNx);
  }

  // Tangent (stiffness) matrices
  //
  for (int b = 0; b < eNoNw; b++) {
    T1 = rho*amd*Nw(b);

    for (int a = 0; a < eNoNq; a++) {
      // dRc_a/dU_b1
      double T2 = Nqx(0,a)*(updu[0][0][b] - T1) + Nqx(1,a)*updu[0][1][b] + Nqx(2,a)*updu[0][2][b];
      lK(12,a,b) = lK(12,a,b) + wl*(Nq(a)*Nwx(0,b) - tauM*T2);

      // dRc_a/dU_b2
      T2 = Nqx(0,a)*updu[1][0][b] + Nqx(1,a)*(updu[1][1][b] - T1) + Nqx(2,a)*updu[1][2][b];
      lK(13,a,b) = lK(13,a,b) + wl*(Nq(a)*Nwx(1,b) - tauM*T2);

      // dRc_a/dU_b3
      T2 = Nqx(0,a)*updu[2][0][b] + Nqx(1,a)*updu[2][1][b] + Nqx(2,a)*(updu[2][2][b] - T1);
      lK(14,a,b) = lK(14,a,b) + wl*(Nq(a)*Nwx(2,b) - tauM*T2);
    }
  }

  if (vmsFlag) {
    for (int b = 0; b < eNoNq; b++) {
      for (int a = 0; a < eNoNq; a++) {
        // dC/dP
        double NxNx = Nqx(0,a)*Nqx(0,b) + Nqx(1,a)*Nqx(1,b) + Nqx(2,a)*Nqx(2,b);
        lK(15,a,b) = lK(15,a,b) + wl*tauM*NxNx;
      }
    }
  }
}

/// @brief Element momentum residual.
///
///  Modifies:
///    lR(dof,eNoN)  - Residual
///    lK(dof*dof,eNoN,eNoN) - Stiffness matrix
//
void fluid_3d_m(ComMod& com_mod, const int vmsFlag, const int eNoNw, const int eNoNq, const double w,
    const Array<double>& Kxi, const Vector<double>& Nw, const Vector<double>& Nq, const Array<double>& Nwx,
    const Array<double>& Nqx, const Array<double>& Nwxx, const Array<double>& al, const Array<double>& yl,
    const Array<double>& bfl, Array<double>& lR, Array3<double>& lK, double K_inverse_darcy_permeability, double DDir)
{
  #define n_debug_fluid_3d_m
  #ifdef debug_fluid_3d_m
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "vmsFlag: " << vmsFlag;
  dmsg << "eNoNw: " << eNoNw;
  dmsg << "eNoNq: " << eNoNq;
  dmsg << "w: " << w;
  double start_time = utils::cput();
  #endif

  // Maximum size of arrays sized by (3,eNoNw) -> (3,MAX_SIZE).
  const int MAX_SIZE = 27;

  using namespace consts;

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  double Res;
  if (!com_mod.urisFlag) {
    Res = 0.0;
  } else {
    Res = com_mod.urisRes;
    if (com_mod.uris[0].clsFlg) {Res = com_mod.urisResClose;}
  }

  double ctM  = 1.0;
  double ctC  = 36.0;

  double rho = dmn.prop[PhysicalProperyType::fluid_density];
  std::array<double,3> f;
  f[0] = dmn.prop[PhysicalProperyType::f_x];
  f[1] = dmn.prop[PhysicalProperyType::f_y];
  f[2] = dmn.prop[PhysicalProperyType::f_z];

  double T1 = eq.af * eq.gam * dt;
  double amd = eq.am / T1;
  double wl = w*T1;
  double wr = w*rho;

  #ifdef debug_fluid_3d_m
  dmsg << "rho: " << rho;
  dmsg << "T1: " << T1;
  dmsg << "wl: " << wl;
  dmsg << "wr: " << wr;
  #endif

  // Note that indices are not selected based on the equation because
  // fluid equation always come first
  // Velocity and its gradients, inertia (acceleration & body force)
  //
  std::array<double,3> ud{-f[0], -f[1], -f[2]};
  //ud  = -f
  double u[3] = {};
  double ux[3][3] = {};
  double uxx[3][3][3] = {};

  for (int a = 0; a < eNoNw; a++) {
    ud[0] = ud[0] + Nw(a)*(al(0,a)-bfl(0,a));
    ud[1] = ud[1] + Nw(a)*(al(1,a)-bfl(1,a));
    ud[2] = ud[2] + Nw(a)*(al(2,a)-bfl(2,a));

    u[0] = u[0] + Nw(a)*yl(0,a);
    u[1] = u[1] + Nw(a)*yl(1,a);
    u[2] = u[2] + Nw(a)*yl(2,a);
    
    ux[0][0] += Nwx(0,a)*yl(0,a);
    ux[1][0] += Nwx(1,a)*yl(0,a);
    ux[2][0] += Nwx(2,a)*yl(0,a);
    ux[0][1] += Nwx(0,a)*yl(1,a);
    ux[1][1] += Nwx(1,a)*yl(1,a);
    ux[2][1] += Nwx(2,a)*yl(1,a);
    ux[0][2] += Nwx(0,a)*yl(2,a);
    ux[1][2] += Nwx(1,a)*yl(2,a);
    ux[2][2] += Nwx(2,a)*yl(2,a);

    uxx[0][0][0] += Nwxx(0,a)*yl(0,a);
    uxx[1][0][1] += Nwxx(1,a)*yl(0,a);
    uxx[2][0][2] += Nwxx(2,a)*yl(0,a);
    uxx[1][0][0] += Nwxx(3,a)*yl(0,a);
    uxx[2][0][1] += Nwxx(4,a)*yl(0,a);
    uxx[0][0][2] += Nwxx(5,a)*yl(0,a);

    uxx[0][1][0] += Nwxx(0,a)*yl(1,a);
    uxx[1][1][1] += Nwxx(1,a)*yl(1,a);
    uxx[2][1][2] += Nwxx(2,a)*yl(1,a);
    uxx[1][1][0] += Nwxx(3,a)*yl(1,a);
    uxx[2][1][1] += Nwxx(4,a)*yl(1,a);
    uxx[0][1][2] += Nwxx(5,a)*yl(1,a);

    uxx[0][2][0] += Nwxx(0,a)*yl(2,a);
    uxx[1][2][1] += Nwxx(1,a)*yl(2,a);
    uxx[2][2][2] += Nwxx(2,a)*yl(2,a);
    uxx[1][2][0] += Nwxx(3,a)*yl(2,a);
    uxx[2][2][1] += Nwxx(4,a)*yl(2,a);
    uxx[0][2][2] += Nwxx(5,a)*yl(2,a);
  }

  double divU = ux[0][0] + ux[1][1] + ux[2][2];
  #ifdef debug_fluid_3d_m
  dmsg << "divU: " << divU;
  #endif

  uxx[0][0][1] = uxx[1][0][0];
  uxx[1][0][2] = uxx[2][0][1];
  uxx[2][0][0] = uxx[0][0][2];

  uxx[0][1][1] = uxx[1][1][0];
  uxx[1][1][2] = uxx[2][1][1];
  uxx[2][1][0] = uxx[0][1][2];

  uxx[0][2][1] = uxx[1][2][0];
  uxx[1][2][2] = uxx[2][2][1];
  uxx[2][2][0] = uxx[0][2][2];

  std::array<double,3> d2u2{0.0};
  d2u2[0] = uxx[0][0][0] + uxx[1][0][1] + uxx[2][0][2];
  d2u2[1] = uxx[0][1][0] + uxx[1][1][1] + uxx[2][1][2];
  d2u2[2] = uxx[0][2][0] + uxx[1][2][1] + uxx[2][2][2];

  // Pressure and its gradient
  //
  double p = 0.0;
  double px[3] = {};

  for (int a = 0; a < eNoNq; a++) {
    p  = p + Nq(a)*yl(3,a);
    px[0] = px[0] + Nqx(0,a)*yl(3,a);
    px[1] = px[1] + Nqx(1,a)*yl(3,a);
    px[2] = px[2] + Nqx(2,a)*yl(3,a);
  }
  #ifdef debug_fluid_3d_m
  dmsg << "u: " << u[0] << " " << u[1] << " " << u[2];
  dmsg << "p: " << p;
  dmsg << "px: " << px[0] << " " << px[1] << " " << px[2];
  #endif

  // Update convection velocity relative to mesh velocity
  //
  if (com_mod.mvMsh) {
    for (int a = 0; a < eNoNw; a++) {
      u[0] = u[0] - Nw(a)*yl(4,a);
      u[1] = u[1] - Nw(a)*yl(5,a);
      u[2] = u[2] - Nw(a)*yl(6,a);
    }
  }

  // Strain rate tensor 2*e_ij := (u_i,j + u_j,i)
  //
  double es[3][3] = {};
  es[0][0] = ux[0][0] + ux[0][0];
  es[1][1] = ux[1][1] + ux[1][1];
  es[2][2] = ux[2][2] + ux[2][2];
  es[1][0] = ux[1][0] + ux[0][1];
  es[2][1] = ux[2][1] + ux[1][2];
  es[0][2] = ux[0][2] + ux[2][0];
  es[0][1] = es[1][0];
  es[1][2] = es[2][1];
  es[2][0] = es[0][2];

  double esNx[3][MAX_SIZE];

  for (int a = 0; a < eNoNw; a++) {
    esNx[0][a] = es[0][0]*Nwx(0,a) + es[1][0]*Nwx(1,a) + es[2][0]*Nwx(2,a);
    esNx[1][a] = es[0][1]*Nwx(0,a) + es[1][1]*Nwx(1,a) + es[2][1]*Nwx(2,a);
    esNx[2][a] = es[0][2]*Nwx(0,a) + es[1][2]*Nwx(1,a) + es[2][2]*Nwx(2,a);
  }

  double es_x[3][3][3];

  for (int k = 0; k < 3; k++) { 
     es_x[0][0][k] = uxx[0][0][k] + uxx[0][0][k];
     es_x[1][1][k] = uxx[1][1][k] + uxx[1][1][k];
     es_x[2][2][k] = uxx[2][2][k] + uxx[2][2][k];
     es_x[1][0][k] = uxx[1][0][k] + uxx[0][1][k];
     es_x[2][1][k] = uxx[2][1][k] + uxx[1][2][k];
     es_x[0][2][k] = uxx[0][2][k] + uxx[2][0][k];

     es_x[0][1][k] = es_x[1][0][k];
     es_x[1][2][k] = es_x[2][1][k];
     es_x[2][0][k] = es_x[0][2][k];
  }

  std::array<double,3> mu_x{0.0};

  mu_x[0] = (es_x[0][0][0]*es[0][0] + es_x[1][1][0]*es[1][1]
          +  es_x[2][2][0]*es[2][2])*0.5
          +  es_x[1][0][0]*es[1][0] + es_x[2][1][0]*es[2][1]
          +  es_x[0][2][0]*es[0][2];

  mu_x[1] = (es_x[0][0][1]*es[0][0] + es_x[1][1][1]*es[1][1]
          +  es_x[2][2][1]*es[2][2])*0.5
          +  es_x[1][0][1]*es[1][0] + es_x[2][1][1]*es[2][1]
          +  es_x[0][2][1]*es[0][2];

  mu_x[2] = (es_x[0][0][2]*es[0][0] + es_x[1][1][2]*es[1][1]
          +  es_x[2][2][2]*es[2][2])*0.5
          +  es_x[1][0][2]*es[1][0] + es_x[2][1][2]*es[2][1]
          +  es_x[0][2][2]*es[0][2];

  #ifdef debug_fluid_3d_m
  dmsg << "mu_x: " << mu_x[0] << " " << mu_x[1] << " " << mu_x[2];
  #endif

  // Shear-rate := (2*e_ij*e_ij)^.5
  //
  //dmsg << "Compute shear rate ... ";
  double gam = es[0][0]*es[0][0] + es[1][0]*es[1][0] + es[2][0]*es[2][0]
             + es[0][1]*es[0][1] + es[1][1]*es[1][1] + es[2][1]*es[2][1]
             + es[0][2]*es[0][2] + es[1][2]*es[1][2] + es[2][2]*es[2][2];
  gam = sqrt(0.5*gam);
  #ifdef debug_fluid_3d_m
  dmsg << "gam: " << gam;
  #endif

  // Compute viscosity based on shear-rate and chosen viscosity model
  // The returned mu_g := (d\mu / d\gamma)
  //
  double mu, mu_s, mu_g;
  fluid::get_viscosity(com_mod, dmn, gam, mu, mu_s, mu_g);

  if (utils::is_zero(gam)) {
     mu_g = 0.0;
  } else {
     mu_g = mu_g / gam;
  }
  std::transform(mu_x.begin(), mu_x.end(), mu_x.begin(), [mu_g](double &v){return mu_g*v;});
  //mu_x(:) = mu_g * mu_x(:)

  // Stabilization parameters
  //
  double kT = 4.0 * pow(ctM/dt,2.0);
  
  // If we consider the NSB model, we need to add an extra term inside the computation for the stab parameter 
  kT = kT + pow(K_inverse_darcy_permeability*mu/rho, 2.0);   

  // In case of unfitted RIS, compute the delta function at the quad point,
  // add the additional value to the stabilization param 
  kT = kT + pow(Res*DDir, 2.0);

  double kU = u[0]*u[0]*Kxi(0,0) + u[1]*u[0]*Kxi(1,0) + u[2]*u[0]*Kxi(2,0)
            + u[0]*u[1]*Kxi(0,1) + u[1]*u[1]*Kxi(1,1) + u[2]*u[1]*Kxi(2,1)
            + u[0]*u[2]*Kxi(0,2) + u[1]*u[2]*Kxi(1,2) + u[2]*u[2]*Kxi(2,2);

  double kS = Kxi(0,0)*Kxi(0,0) + Kxi(1,0)*Kxi(1,0) + Kxi(2,0)*Kxi(2,0)
            + Kxi(0,1)*Kxi(0,1) + Kxi(1,1)*Kxi(1,1) + Kxi(2,1)*Kxi(2,1)
            + Kxi(0,2)*Kxi(0,2) + Kxi(1,2)*Kxi(1,2) + Kxi(2,2)*Kxi(2,2);
  kS = ctC * kS * pow(mu/rho,2.0);
  double tauM = 1.0 / (rho * sqrt( kT + kU + kS ));

  #ifdef debug_fluid_3d_m
  dmsg << "kT: " << kT;
  dmsg << "kU: " << kU;
  dmsg << "kS: " << kS;
  dmsg << "tauM: " << tauM;
  #endif

  double rV[3] = {};
  rV[0] = ud[0] + u[0]*ux[0][0] + u[1]*ux[1][0] + u[2]*ux[2][0];
  rV[1] = ud[1] + u[0]*ux[0][1] + u[1]*ux[1][1] + u[2]*ux[2][1];
  rV[2] = ud[2] + u[0]*ux[0][2] + u[1]*ux[1][2] + u[2]*ux[2][2];

  double rS[3] = {};
  rS[0] = mu_x[0]*es[0][0] + mu_x[1]*es[1][0] + mu_x[2]*es[2][0] + mu*d2u2[0];
  rS[1] = mu_x[0]*es[0][1] + mu_x[1]*es[1][1] + mu_x[2]*es[2][1] + mu*d2u2[1];
  rS[2] = mu_x[0]*es[0][2] + mu_x[1]*es[1][2] + mu_x[2]*es[2][2] + mu*d2u2[2];

  double up[3] = {};
  // up[0] = -tauM*(rho*rV[0] + px[0] - rS[0] + mu*K_inverse_darcy_permeability * u[0]);
  // up[1] = -tauM*(rho*rV[1] + px[1] - rS[1] + mu*K_inverse_darcy_permeability * u[1]);
  // up[2] = -tauM*(rho*rV[2] + px[2] - rS[2] + mu*K_inverse_darcy_permeability * u[2]);

  up[0] = -tauM*(rho*rV[0] + px[0] - rS[0] + mu*K_inverse_darcy_permeability * u[0]
                 + (Res*DDir)*u[0]);
  up[1] = -tauM*(rho*rV[1] + px[1] - rS[1] + mu*K_inverse_darcy_permeability * u[1]
                 + (Res*DDir)*u[1]);
  up[2] = -tauM*(rho*rV[2] + px[2] - rS[2] + mu*K_inverse_darcy_permeability * u[2]
                 + (Res*DDir)*u[2]);

  double tauC, tauB, pa;
  double eps = std::numeric_limits<double>::epsilon();
  double ua[3] = {};

  if (vmsFlag) {
    tauC = 1.0 / (tauM * (Kxi(0,0) + Kxi(1,1) + Kxi(2,2)));
    tauB = up[0]*up[0]*Kxi(0,0) + up[1]*up[0]*Kxi(1,0)
         + up[2]*up[0]*Kxi(2,0) + up[0]*up[1]*Kxi(0,1)
         + up[1]*up[1]*Kxi(1,1) + up[2]*up[1]*Kxi(2,1)
         + up[0]*up[2]*Kxi(0,2) + up[1]*up[2]*Kxi(1,2)
         + up[2]*up[2]*Kxi(2,2);

    if (utils::is_zero(tauB)) {
      tauB = eps;
    }
    tauB = rho / sqrt(tauB);

    ua[0] = u[0] + up[0];
    ua[1] = u[1] + up[1];
    ua[2] = u[2] + up[2];
    pa = p - tauC*divU;

  } else {
    tauC = 0.0;
    tauB = 0.0;
    for (int i = 0; i < 3; i++) {
      ua[i] = u[i];
    }
    pa = p;
  }

  rV[0] = tauB*(up[0]*ux[0][0] + up[1]*ux[1][0] + up[2]*ux[2][0]);
  rV[1] = tauB*(up[0]*ux[0][1] + up[1]*ux[1][1] + up[2]*ux[2][1]);
  rV[2] = tauB*(up[0]*ux[0][2] + up[1]*ux[1][2] + up[2]*ux[2][2]);

  double rM[3][3];
  rM[0][0] = mu*es[0][0] - rho*up[0]*ua[0] + rV[0]*up[0] - pa;
  rM[1][0] = mu*es[1][0] - rho*up[0]*ua[1] + rV[0]*up[1];
  rM[2][0] = mu*es[2][0] - rho*up[0]*ua[2] + rV[0]*up[2];

  rM[0][1] = mu*es[0][1] - rho*up[1]*ua[0] + rV[1]*up[0];
  rM[1][1] = mu*es[1][1] - rho*up[1]*ua[1] + rV[1]*up[1] - pa;
  rM[2][1] = mu*es[2][1] - rho*up[1]*ua[2] + rV[1]*up[2];

  rM[0][2] = mu*es[0][2] - rho*up[2]*ua[0] + rV[2]*up[0];
  rM[1][2] = mu*es[1][2] - rho*up[2]*ua[1] + rV[2]*up[1];
  rM[2][2] = mu*es[2][2] - rho*up[2]*ua[2] + rV[2]*up[2] - pa;

  rV[0] = ud[0] + ua[0]*ux[0][0] + ua[1]*ux[1][0] + ua[2]*ux[2][0];
  rV[1] = ud[1] + ua[0]*ux[0][1] + ua[1]*ux[1][1] + ua[2]*ux[2][1];
  rV[2] = ud[2] + ua[0]*ux[0][2] + ua[1]*ux[1][2] + ua[2]*ux[2][2];

  //  Local residual
  //
  double updu[3][3][MAX_SIZE] = {};
  double uNx[MAX_SIZE] = {};
  double upNx[MAX_SIZE] = {}; 
  double uaNx[MAX_SIZE] = {}; 

  for (int a = 0; a < eNoNw; a++) {
    lR(0,a) = lR(0,a) + wr*Nw(a)*rV[0] + w*(Nwx(0,a)*rM[0][0] + Nwx(1,a)*rM[1][0] + Nwx(2,a)*rM[2][0]);
    lR(1,a) = lR(1,a) + wr*Nw(a)*rV[1] + w*(Nwx(0,a)*rM[0][1] + Nwx(1,a)*rM[1][1] + Nwx(2,a)*rM[2][1]);
    lR(2,a) = lR(2,a) + wr*Nw(a)*rV[2] + w*(Nwx(0,a)*rM[0][2] + Nwx(1,a)*rM[1][2] + Nwx(2,a)*rM[2][2]);

    // Quantities used for stiffness matrix

    uNx[a]  = u[0]*Nwx(0,a)  + u[1]*Nwx(1,a)  + u[2]*Nwx(2,a);
    upNx[a] = up[0]*Nwx(0,a) + up[1]*Nwx(1,a) + up[2]*Nwx(2,a);

    if (vmsFlag) {
       uaNx[a] = uNx[a] + upNx[a];
    } else {
       uaNx[a] = uNx[a];
    }

    // T1 = -rho*uNx[a] + mu*(Nwxx(0,a) + Nwxx(1,a) + Nwxx(2,a)) + mu_x[0]*Nwx(0,a) + mu_x[1]*Nwx(1,a) + mu_x[2]*Nwx(2,a) - mu*K_inverse_darcy_permeability*Nw(a);

    T1 = -rho*uNx[a] + mu*(Nwxx(0,a) + Nwxx(1,a) + Nwxx(2,a)) 
         + mu_x[0]*Nwx(0,a) + mu_x[1]*Nwx(1,a) + mu_x[2]*Nwx(2,a) 
         - mu*K_inverse_darcy_permeability*Nw(a)
         - (Res*DDir)*Nw(a);

    updu[0][0][a] = mu_x[0]*Nwx(0,a) + d2u2[0]*mu_g*esNx[0][a] + T1;
    updu[1][0][a] = mu_x[1]*Nwx(0,a) + d2u2[1]*mu_g*esNx[0][a];
    updu[2][0][a] = mu_x[2]*Nwx(0,a) + d2u2[2]*mu_g*esNx[0][a];

    updu[0][1][a] = mu_x[0]*Nwx(1,a) + d2u2[0]*mu_g*esNx[1][a];
    updu[1][1][a] = mu_x[1]*Nwx(1,a) + d2u2[1]*mu_g*esNx[1][a] + T1;
    updu[2][1][a] = mu_x[2]*Nwx(1,a) + d2u2[2]*mu_g*esNx[1][a];

    updu[0][2][a] = mu_x[0]*Nwx(2,a) + d2u2[0]*mu_g*esNx[2][a];
    updu[1][2][a] = mu_x[1]*Nwx(2,a) + d2u2[1]*mu_g*esNx[2][a];
    updu[2][2][a] = mu_x[2]*Nwx(2,a) + d2u2[2]*mu_g*esNx[2][a] + T1;
  }

  // Tangent (stiffness) matrices
  //
  for (int b = 0; b < eNoNw; b++) {
    for (int a = 0; a < eNoNw; a++) {
      rM[0][0] = Nwx(0,a)*Nwx(0,b);
      rM[1][0] = Nwx(1,a)*Nwx(0,b);
      rM[2][0] = Nwx(2,a)*Nwx(0,b);
      rM[0][1] = Nwx(0,a)*Nwx(1,b);
      rM[1][1] = Nwx(1,a)*Nwx(1,b);
      rM[2][1] = Nwx(2,a)*Nwx(1,b);
      rM[0][2] = Nwx(0,a)*Nwx(2,b);
      rM[1][2] = Nwx(1,a)*Nwx(2,b);
      rM[2][2] = Nwx(2,a)*Nwx(2,b);

      double NxNx = Nwx(0,a)*Nwx(0,b) + Nwx(1,a)*Nwx(1,b) + Nwx(2,a)*Nwx(2,b);
      T1 = mu*NxNx + rho*amd*Nw(b)*(Nw(a) + rho*tauM*uaNx[a]) + rho*Nw(a)*(uNx[b]+upNx[b]) + tauB*upNx[a]*upNx[b];

      // dRm_a1/du_b1
      double T2 = (mu + tauC)*rM[0][0] + esNx[0][a]*mu_g*esNx[0][b] - rho*tauM*uaNx[a]*updu[0][0][b];
      lK(0,a,b)  = lK(0,a,b)  + wl*(T2 + T1);
      // lK(0,a,b)  = lK(0,a,b)  + mu*K_inverse_darcy_permeability*wl*Nw(b)*Nw(a);
      lK(0,a,b)  = lK(0,a,b)  + mu*K_inverse_darcy_permeability*wl*Nw(b)*Nw(a)
                              + (Res*DDir)*wl*Nw(b)*Nw(a);

      // dRm_a1/du_b2
      T2 = mu*rM[1][0] + tauC*rM[0][1] + esNx[0][a]*mu_g*esNx[1][b] - rho*tauM*uaNx[a]*updu[1][0][b];
      lK(1,a,b)  = lK(1,a,b)  + wl*(T2);

      // dRm_a1/du_b3
      T2 = mu*rM[2][0] + tauC*rM[0][2] + esNx[0][a]*mu_g*esNx[2][b] - rho*tauM*uaNx[a]*updu[2][0][b];
      lK(2,a,b)  = lK(2,a,b)  + wl*(T2);

      // dRm_a2/du_b1
      T2 = mu*rM[0][1] + tauC*rM[1][0] + esNx[1][a]*mu_g*esNx[0][b] - rho*tauM*uaNx[a]*updu[0][1][b];
      lK(4,a,b) = lK(4,a,b)  + wl*(T2);

      // dRm_a2/du_b2
      T2 = (mu + tauC)*rM[1][1] + esNx[1][a]*mu_g*esNx[1][b] - rho*tauM*uaNx[a]*updu[1][1][b];
      lK(5,a,b)  = lK(5,a,b)  + wl*(T2 + T1);
      // lK(5,a,b)  = lK(5,a,b)  + mu*K_inverse_darcy_permeability*wl*Nw(b)*Nw(a);
      lK(5,a,b)  = lK(5,a,b)  + mu*K_inverse_darcy_permeability*wl*Nw(b)*Nw(a)
                              + (Res*DDir)*wl*Nw(b)*Nw(a);

      // dRm_a2/du_b3
      T2 = mu*rM[2][1] + tauC*rM[1][2] + esNx[1][a]*mu_g*esNx[2][b] - rho*tauM*uaNx[a]*updu[2][1][b];
      lK(6,a,b)  = lK(6,a,b)  + wl*(T2);

      // dRm_a3/du_b1
      T2 = mu*rM[0][2] + tauC*rM[2][0] + esNx[2][a]*mu_g*esNx[0][b] - rho*tauM*uaNx[a]*updu[0][2][b];
      lK(8,a,b)  = lK(8,a,b)  + wl*(T2);

      // dRm_a3/du_b2
      T2 = mu*rM[1][2] + tauC*rM[2][1] + esNx[2][a]*mu_g*esNx[1][b] - rho*tauM*uaNx[a]*updu[1][2][b];
      lK(9,a,b) = lK(9,a,b) + wl*(T2);

      // dRm_a3/du_b3;
      T2 = (mu + tauC)*rM[2][2] + esNx[2][a]*mu_g*esNx[2][b] - rho*tauM*uaNx[a]*updu[2][2][b];
      lK(10,a,b) = lK(10,a,b) + wl*(T2 + T1);
      // lK(10,a,b) = lK(10,a,b) + mu*K_inverse_darcy_permeability*wl*Nw(b)*Nw(a);
      lK(10,a,b) = lK(10,a,b) + mu*K_inverse_darcy_permeability*wl*Nw(b)*Nw(a)
                              + (Res*DDir)*wl*Nw(b)*Nw(a);
      //dmsg << "lK(10,a,b): " << lK(10,a,b);
    }
  }

  //exit(0);

  for (int b = 0; b < eNoNq; b++) {
    for (int a = 0; a < eNoNw; a++) {
      T1 = rho*tauM*uaNx[a];

      // dRm_a1/dp_b
      lK(3,a,b)  = lK(3,a,b)  - wl*(Nwx(0,a)*Nq(b) - Nqx(0,b)*T1);

      // dRm_a2/dp_b
      lK(7,a,b)  = lK(7,a,b)  - wl*(Nwx(1,a)*Nq(b) - Nqx(1,b)*T1);

      // dRm_a3/dp_b
      lK(11,a,b) = lK(11,a,b) - wl*(Nwx(2,a)*Nq(b) - Nqx(2,b)*T1);
    }
  }
  
  // Residual contribution Birkman term 
  // Local residue
  for (int a = 0; a < eNoNw; a++) {
      lR(0,a) = lR(0,a) + mu*K_inverse_darcy_permeability*w*Nw(a)*(u[0]+up[0])
                        + Res*DDir*w*Nw(a)*(u[0]+up[0]);
      lR(1,a) = lR(1,a) + mu*K_inverse_darcy_permeability*w*Nw(a)*(u[1]+up[1])
                        + Res*DDir*w*Nw(a)*(u[1]+up[1]);
      lR(2,a) = lR(2,a) + mu*K_inverse_darcy_permeability*w*Nw(a)*(u[2]+up[2])
                        + Res*DDir*w*Nw(a)*(u[2]+up[2]);
  }

}


void get_viscosity(const ComMod& com_mod, const dmnType& lDmn, double& gamma, double& mu, double& mu_s, double& mu_x)
{
  using namespace consts;
  
  // effective dynamic viscosity
  mu = 0.0;
  
  mu_s = 0.0;
  
  // derivative of effective dynamic viscosity with respect to gamma
  mu_x = 0.0;

  double mu_i, mu_o, lam, a, n, T1, T2;

  switch (lDmn.fluid_visc.viscType) {

    case FluidViscosityModelType::viscType_Const:
      mu = lDmn.fluid_visc.mu_i;
      mu_s = mu;
      mu_x = 0.0;
    break;
    
    // Carreau-Yasuda
    case FluidViscosityModelType::viscType_CY:
      mu_i = lDmn.fluid_visc.mu_i;
      mu_o = lDmn.fluid_visc.mu_o;
      lam = lDmn.fluid_visc.lam;
      a = lDmn.fluid_visc.a;
      n = lDmn.fluid_visc.n;

      T1 = 1.0 + pow(lam*gamma, a);
      T2 = pow(T1,((n-1.0)/a));
      mu = mu_i + (mu_o - mu_i)*T2;
      mu_s = mu_i;

      T1 = T2 / T1;
      T2 = pow(lam,a) * pow(gamma,(a-1.0)) * T1;
      mu_x = (mu_o - mu_i) * (n - 1.0) * T2;
    break;

    // Casson
    case FluidViscosityModelType::viscType_Cass:
      mu_i = lDmn.fluid_visc.mu_i;
      mu_o = lDmn.fluid_visc.mu_o;
      lam  = lDmn.fluid_visc.lam;

      if (gamma < lam) { 
         mu_o = mu_o / sqrt(lam);
         gamma = lam;
      } else { 
         mu_o = mu_o / sqrt(gamma);
      }

      mu  = (mu_i + mu_o) * (mu_i + mu_o);
      mu_s = mu_i * mu_i;
      mu_x = 2.0 * mu_o * (mu_o + mu_i) / gamma;
    break;
  } 
}

};
