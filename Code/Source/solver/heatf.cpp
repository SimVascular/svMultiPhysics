// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "heatf.h"

#include "all_fun.h"
#include "lhsa.h"
#include "mat_fun.h"
#include "nn.h"
#include "utils.h"

#include <math.h>

namespace heatf {

void b_heatf(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Vector<double>& y,
    const double h, const Vector<double>& nV, Array<double>& lR, Array3<double>& lK)
{
  const int nsd = com_mod.nsd;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  const double dt = com_mod.dt;

  double wl = w * eq.af * eq.gam * dt;
  double T = y(eq.s);
  double udn = 0.0;

  for (int i = 0; i < nsd; i++) {
    udn = udn + y(i)*nV(i);
  }

  udn = 0.5 * (udn - fabs(udn));
  double T1  = h - udn*T;

  for (int a = 0; a < eNoN; a++) { 
    lR(0,a) = lR(0,a) + w*N(a)*T1;

    for (int b = 0; b < eNoN; b++) { 
      lK(0,a,b) = lK(0,a,b) - wl*N(a)*N(b)*udn;
    }
  }
}


void construct_heatf(ComMod& com_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Yg)
{
  #define n_debug_construct_heatf 
  #ifdef debug_construct_heatf
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  using namespace consts;

  const int nsd  = com_mod.nsd;
  const int tDof = com_mod.tDof;
  const int dof = com_mod.dof;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  auto& cDmn = com_mod.cDmn;

  int eNoN = lM.eNoN;
  #ifdef debug_construct_heatf
  dmsg << "cEq: " << cEq;
  dmsg << "cDmn: " << cDmn;
  #endif

  Vector<int> ptr(eNoN);
  Vector<double> N(eNoN);
  Array<double> xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN), Nx(nsd,eNoN), lR(dof,eNoN);
  Array3<double> lK(dof*dof,eNoN,eNoN);
  Array<double> ksix(nsd,nsd);

  for (int e = 0; e < lM.nEl; e++) {
    // Update domain and proceed if domain phys and eqn phys match
    cDmn = all_fun::domain(com_mod, lM, cEq, e);
    auto cPhys = eq.dmn[cDmn].phys;
    if (cPhys != EquationType::phys_heatF) {
      continue;
    }

    // Update shape functions for NURBS
    if (lM.eType == ElementType::NRB) {
      //CALL NRBNNX(lM, e)
    }

    // Create local copies
    //
    for (int a = 0; a < eNoN; a++) {
      int Ac = lM.IEN(a,e);
      ptr(a) = Ac;

      for (int i = 0; i < nsd; i++) {
        xl(i,a) = com_mod.x(i,Ac);
      }

      for (int i = 0; i < tDof; i++) {
        al(i,a) = Ag(i,Ac);
        yl(i,a) = Yg(i,Ac);
      }
    }

    // Gauss integration
    //
    lR = 0.0;
    lK = 0.0;
    double Jac{0.0};

    for (int g = 0; g < lM.nG; g++) {
      if (g == 0 || !lM.lShpF) {
        auto Nx_g = lM.Nx.slice(g);
        nn::gnn(eNoN, nsd, nsd, Nx_g, xl, Nx, Jac, ksix);
        if (utils::is_zero(Jac)) {
          throw std::runtime_error("[construct_heatf] Jacobian for element " + std::to_string(e) + " is < 0.");
        }
      }

      double w = lM.w(g) * Jac;
      N = lM.N.col(g);

      if (nsd == 3) {
        heatf_3d(com_mod, eNoN, w, N, Nx, al, yl, ksix, lR, lK);

      } else if (nsd == 2) {
        heatf_2d(com_mod, eNoN, w, N, Nx, al, yl, ksix, lR, lK);
      }
    } // for g = 0

    eq.linear_algebra->assemble(com_mod, eNoN, ptr, lK, lR);

  } // for e = 0
}

void heatf_2d(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx,
    const Array<double>& al, const Array<double>& yl, const Array<double>& ksix, Array<double>& lR, Array3<double>& lK)
{
  #define n_debug_heatf_2d 
  #ifdef debug_heatf_2d 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "w: " << w;
  #endif

  static Vector<double> ct({4.0, 1.0, 3.0, 1.0});

  using namespace consts;
  using namespace mat_fun;

  const int nsd = com_mod.nsd;
  const int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  const int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;
  const int i = eq.s;

  double T1 = eq.af * eq.gam * dt;
  double amd = eq.am / T1;
  double nu = dmn.prop.at(PhysicalProperyType::conductivity);
  double s = dmn.prop.at(PhysicalProperyType::source_term);
  double wl = w * T1;

  #ifdef debug_heats_2d 
  dmsg << "nu: " << nu;
  dmsg << "s: " << s;
  dmsg << "T1: " << T1;
  dmsg << "i: " << i;
  dmsg << "wl: " << wl;
  #endif

  double Td = -s;
  Vector<double> Tx(nsd), u(nsd), udNx(eNoN);

  for (int a = 0; a < eNoN; a++) {
    u(0) = u(0) + N(a)*yl(0,a);
    u(1) = u(1) + N(a)*yl(1,a);

    Td = Td + N(a)*al(i,a);

    Tx(0) = Tx(0) + Nx(0,a)*yl(i,a);
    Tx(1) = Tx(1) + Nx(1,a)*yl(i,a);
  }

  if (com_mod.mvMsh) {
    for (int a = 0; a < eNoN; a++) {
      u(0) = u(0) - N(a)*yl(4,a);
      u(1) = u(1) - N(a)*yl(5,a);
      u(2) = u(2) - N(a)*yl(6,a);
    }
  } 

  double kU = u(0)*u(0)*ksix(0,0) + u(1)*u(0)*ksix(1,0) + u(0)*u(1)*ksix(0,1) + u(1)*u(1)*ksix(1,1);
  double kS = ksix(0,0)*ksix(0,0) + ksix(1,0)*ksix(1,0) + ksix(0,1)*ksix(0,1) + ksix(1,1)*ksix(1,1);
  double nTx = ksix(0,0)*Tx(0)*Tx(0) + ksix(1,1)*Tx(1)*Tx(1) + (ksix(0,1) + ksix(1,0))*Tx(0)*Tx(1);

  if (utils::is_zero(nTx)) {
    nTx = std::numeric_limits<double>::epsilon();
  }

  double udTx = u(0)*Tx(0) + u(1)*Tx(1);
  double Tp = fabs(Td + udTx);
  nu  = nu + 0.5 * Tp / sqrt(nTx);
  double tauM = ct(3) / sqrt((ct(0)/(dt*dt)) + ct(1)*kU + ct(2)*nu*nu*kS);
  Tp = -tauM*(Td + udTx);

  for (int a = 0; a < eNoN; a++) {
    udNx(a) = u(0)*Nx(0,a) + u(1)*Nx(1,a);
  }

  for (int a = 0; a < eNoN; a++) {
    lR(0,a) = lR(0,a) + w*(N(a)*(Td + udTx) + (Nx(0,a)*Tx(0) + Nx(1,a)*Tx(1))*nu - udNx(a)*Tp);

    for (int b = 0; b < eNoN; b++) {
      lK(0,a,b) = lK(0,a,b) + wl*(nu*(Nx(0,a)*Nx(0,b) + Nx(1,a)*Nx(1,b)) + (N(a) + tauM*udNx(a))*(N(b)*amd + udNx(b)));
    }
  }
}


void heatf_3d(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx,
    const Array<double>& al, const Array<double>& yl, const Array<double>& ksix, Array<double>& lR, Array3<double>& lK)
{
  #define n_debug_heatf_3d 
  #ifdef debug_heatf_3d 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "w: " << w;
  dmsg << "N: " << N;
  dmsg << "Nx: " << Nx;
  dmsg << "yl: " << yl;
  dmsg << "al: " << al;
  #endif

  static Vector<double> ct({4.0, 1.0, 3.0, 1.0});

  using namespace consts;
  using namespace mat_fun;

  const int nsd = com_mod.nsd;
  const int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  const int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;
  const int i = eq.s;

  double T1 = eq.af * eq.gam * dt;
  double amd = eq.am / T1;
  double nu = dmn.prop.at(PhysicalProperyType::conductivity);
  double s = dmn.prop.at(PhysicalProperyType::source_term);
  double wl = w * T1;

  #ifdef debug_heatf_3d 
  dmsg;
  dmsg << "eNoN: " << eNoN;
  dmsg << "nu: " << nu;
  dmsg << "s: " << s;
  dmsg << "T1: " << T1;
  dmsg << "i: " << i;
  dmsg << "wl: " << wl;
  dmsg << "ct: " << ct;
  #endif

  double Td = -s;
  Vector<double> Tx(nsd), u(nsd), udNx(eNoN);

  for (int a = 0; a < eNoN; a++) {
    u(0) = u(0) + N(a)*yl(0,a);
    u(1) = u(1) + N(a)*yl(1,a);
    u(2) = u(2) + N(a)*yl(2,a);

    Td = Td + N(a)*al(i,a);

    Tx(0) = Tx(0) + Nx(0,a)*yl(i,a);
    Tx(1) = Tx(1) + Nx(1,a)*yl(i,a);
    Tx(2) = Tx(2) + Nx(2,a)*yl(i,a);
  }

  if (com_mod.mvMsh) {
    for (int a = 0; a < eNoN; a++) {
      u(0) = u(0) - N(a)*yl(4,a);
      u(1) = u(1) - N(a)*yl(5,a);
      u(2) = u(2) - N(a)*yl(6,a);
    }
  }

  double kU = u(0)*u(0)*ksix(0,0) + u(1)*u(0)*ksix(1,0) + u(2)*u(0)*ksix(2,0) + u(0)*u(1)*ksix(0,1) + 
              u(1)*u(1)*ksix(1,1) + u(2)*u(1)*ksix(2,1) + u(0)*u(2)*ksix(0,2) + u(1)*u(2)*ksix(1,2) + 
              u(2)*u(2)*ksix(2,2);

  double kS = ksix(0,0)*ksix(0,0) + ksix(1,0)*ksix(1,0) + ksix(2,0)*ksix(2,0) + ksix(0,1)*ksix(0,1) + 
              ksix(1,1)*ksix(1,1) + ksix(2,1)*ksix(2,1) + ksix(0,2)*ksix(0,2) + ksix(1,2)*ksix(1,2) + 
              ksix(2,2)*ksix(2,2);

  double nTx = ksix(0,0)*Tx(0)*Tx(0) + ksix(1,1)*Tx(1)*Tx(1) + ksix(2,2)*Tx(2)*Tx(2) + 
               (ksix(0,1) + ksix(1,0))*Tx(0)*Tx(1) + (ksix(0,2) + ksix(2,0))*Tx(0)*Tx(2) + 
               (ksix(1,2) + ksix(2,1))*Tx(1)*Tx(2);

  #ifdef debug_heatf_3d 
  dmsg << "u: " << u;
  dmsg << "Tx: " << Tx;
  dmsg << "kU: " << kU;
  dmsg << "kS: " << kS;
  dmsg << "nTx: " << nTx<< std::endl;
  #endif

  if (utils::is_zero(nTx)) {
    nTx = std::numeric_limits<double>::epsilon();
  }

  double udTx = u(0)*Tx(0) + u(1)*Tx(1) + u(2)*Tx(2);
  double Tp = fabs(Td + udTx);
  nu = nu + 0.5 * Tp / sqrt(nTx);
  double tauM = ct(3) / sqrt((ct(0)/(dt*dt)) + ct(1)*kU + ct(2)*nu*nu*kS);
  Tp = -tauM*(Td + udTx);

  for (int a = 0; a < eNoN; a++) {
    udNx(a) = u(0)*Nx(0,a) + u(1)*Nx(1,a) + u(2)*Nx(2,a);
  }

  for (int a = 0; a < eNoN; a++) {
    lR(0,a) = lR(0,a) + w*(N(a)*(Td + udTx) + (Nx(0,a)*Tx(0) + Nx(1,a)*Tx(1) + Nx(2,a)*Tx(2))*nu - udNx(a)*Tp);

    for (int b = 0; b < eNoN; b++) {
      lK(0,a,b) = lK(0,a,b) + wl*(nu*(Nx(0,a)*Nx(0,b) + Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b)) + 
                  (N(a) + tauM*udNx(a))*(N(b)*amd + udNx(b)));
    }
  }
}

};
