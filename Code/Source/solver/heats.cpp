// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "heats.h"

#include "all_fun.h"
#include "lhsa.h"
#include "mat_fun.h"
#include "nn.h"
#include "utils.h"

namespace heats {

void b_heats(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const double h, Array<double>& lR)
{
  for (int a = 0; a < eNoN; a++) {
    lR(0,a) = lR(0,a) + w * N(a) * h;
  } 
}


void construct_heats(ComMod& com_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Yg)
{
  #define n_debug_construct_heats 
  #ifdef debug_construct_heats
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
  #ifdef debug_construct_heats
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
    if (cPhys != EquationType::phys_heatS) {
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
          throw std::runtime_error("[construct_heats] Jacobian for element " + std::to_string(e) + " is < 0.");
        }
      }

      double w = lM.w(g) * Jac;
      N = lM.N.col(g);

      if (nsd == 3) {
        heats_3d(com_mod, eNoN, w, N, Nx, al, yl, lR, lK);
      } else if (nsd == 2) {
        heats_2d(com_mod, eNoN, w, N, Nx, al, yl, lR, lK);
      }
    }

    eq.linear_algebra->assemble(com_mod, eNoN, ptr, lK, lR);
  }
}

void heats_2d(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx, 
    const Array<double>& al, const Array<double>& yl, Array<double>& lR, Array3<double>& lK)
{
  #define n_debug_heats_2d 
  #ifdef debug_heats_2d 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  using namespace consts;
  using namespace mat_fun;

  const int nsd = com_mod.nsd;
  const int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  const int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;
  const int i = eq.s;

  double nu = dmn.prop.at(PhysicalProperyType::conductivity);
  double s = dmn.prop.at(PhysicalProperyType::source_term);
  double rho = dmn.prop.at(PhysicalProperyType::solid_density);

  double T1 = eq.af * eq.gam * dt;
  double amd = eq.am * rho / T1;
  double wl = w * T1;

  #ifdef debug_heats_2d 
  dmsg;
  dmsg << "nu: " << nu;
  dmsg << "s: " << s;
  dmsg << "rho: " << rho ;
  dmsg << "T1: " << T1;
  dmsg << "i: " << i;
  dmsg << "wl: " << wl;
  #endif

  double Td = -s;
  Vector<double> Tx(nsd);

  for (int a = 0; a < eNoN; a++) {
    Td = Td + N(a)*al(i,a);
    Tx(0) = Tx(0) + Nx(0,a)*yl(i,a);
    Tx(1) = Tx(1) + Nx(1,a)*yl(i,a);
  }

  Td = Td * rho;
  #ifdef debug_heats_2d 
  dmsg;
  dmsg << "Td: " << Td;
  dmsg << "Tx: " << Tx;
  #endif

  for (int a = 0; a < eNoN; a++) {
    lR(0,a) = lR(0,a) + w*(N(a)*Td + (Nx(0,a)*Tx(0) + Nx(1,a)*Tx(1))*nu);

    for (int b = 0; b < eNoN; b++) {
      lK(0,a,b) = lK(0,a,b) + wl*(N(a)*N(b)*amd + nu*(Nx(0,a)*Nx(0,b) + Nx(1,a)*Nx(1,b)));
    }
  }
}


void heats_3d(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx,
    const Array<double>& al, const Array<double>& yl, Array<double>& lR, Array3<double>& lK)
{
  #define n_debug_heats_3d 
  #ifdef debug_heats_3d 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "w: " << w;
  #endif

  using namespace consts;
  using namespace mat_fun;

  const int nsd = com_mod.nsd;
  const int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  const int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;
  const int i = eq.s;

  double nu = dmn.prop.at(PhysicalProperyType::conductivity);
  double s = dmn.prop.at(PhysicalProperyType::source_term);
  double rho = dmn.prop.at(PhysicalProperyType::solid_density);

  double T1 = eq.af * eq.gam * dt;
  double amd = eq.am * rho / T1;
  double wl = w * T1;

  #ifdef debug_heats_3d 
  dmsg;
  dmsg << "nu: " << nu;
  dmsg << "s: " << s;
  dmsg << "rho: " << rho ;
  dmsg << "T1: " << T1;
  dmsg << "i: " << i;
  dmsg << "wl: " << wl;
  #endif

  double Td = -s;
  Vector<double> Tx(nsd);

  for (int a = 0; a < eNoN; a++) {
    Td = Td + N(a)*al(i,a);
    Tx(0) = Tx(0) + Nx(0,a)*yl(i,a);
    Tx(1) = Tx(1) + Nx(1,a)*yl(i,a);
    Tx(2) = Tx(2) + Nx(2,a)*yl(i,a);
  }

  Td = Td * rho;

  for (int a = 0; a < eNoN; a++) {
    lR(0,a) = lR(0,a) + w*(N(a)*Td + (Nx(0,a)*Tx(0) + Nx(1,a)*Tx(1) + Nx(2,a)*Tx(2))*nu);

    for (int b = 0; b < eNoN; b++) {
      lK(0,a,b) = lK(0,a,b) + wl*(N(a)*N(b)*amd + nu*(Nx(0,a)*Nx(0,b) +Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b)));
    }
  }
}

};
