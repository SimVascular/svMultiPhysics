// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "stokes.h"

#include "all_fun.h"
#include "consts.h"
#include "fs.h"
#include "lhsa.h"
#include "nn.h"
#include "utils.h"
//
// These routines are for solving the Stokes equations.
//
#include <array>
#include <iomanip>
#include <math.h>

namespace stokes {

void construct_stokes(ComMod& com_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Yg)
{
  #define n_debug_construct_stokes
  #ifdef debug_construct_stokes
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "lM.nFs: " << lM.nFs;
  #endif

  using namespace consts;

  const int eNoN = lM.eNoN;
  bool lStab = false;

  if (lM.nFs == 1) {
     lStab = true;
  } else {
     lStab = false;
  }

  // l = 3, if nsd==2 ; else 6;
  const int l = com_mod.nsymd;
  const int nsd  = com_mod.nsd;
  const int tDof = com_mod.tDof;
  const int dof = com_mod.dof;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  auto& cDmn = com_mod.cDmn;

  #ifdef debug_construct_stokes
  dmsg << "cEq: " << cEq;
  dmsg << "eq.sym: " << eq.sym;
  dmsg << "eq.dof: " << eq.dof;
  dmsg << "eq.assmTLS: " << eq.assmTLS;

  dmsg << "eNoN: " << eNoN;
  dmsg << "lStab: " << lStab;
  dmsg << "tDof: " << tDof;
  dmsg << "dof: " << dof;
  dmsg << "lM.nEl: " <<  lM.nEl;
  dmsg << "nsd: " <<  nsd;
  #endif

  // FLUID: dof = nsd+1
  Vector<int> ptr(eNoN); 
  Array<double> xl(nsd,eNoN); 
  Array<double> al(tDof,eNoN); 
  Array<double> yl(tDof,eNoN); 
  Array<double> bfl(nsd,eNoN); 
  Array<double> lR(dof,eNoN); 
  Array3<double> lK(dof*dof,eNoN,eNoN);

  // Loop over all elements of mesh
  //
  int num_c = lM.nEl / 10;

  for (int e = 0; e < lM.nEl; e++) {
    //Update domain and proceed if domain phys and eqn phys match
    cDmn = all_fun::domain(com_mod, lM, cEq, e);
    auto cPhys = eq.dmn[cDmn].phys;
    if (cPhys != EquationType::phys_stokes) {
      continue;
    }
    bool pr_active = false;

    //  Update shape functions for NURBS
    if (lM.eType == ElementType::NRB) {
      //CALL NRBNNX(lM, e)
    }

    // Create local copies
    for (int a = 0; a < eNoN; a++) {
      int Ac = lM.IEN(a,e);
      ptr(a) = Ac;

      for (int i = 0; i < nsd; i++) {
        xl(i,a) = com_mod.x(i,Ac);
        bfl(i,a) = com_mod.Bf(i,Ac);
     }
      for (int i = 0; i < tDof; i++) {
        al(i,a) = Ag(i,Ac);
        yl(i,a) = Yg(i,Ac);
      }
    }

    // Initialize residual and tangents
    lR = 0.0;
    lK = 0.0;
    std::array<fsType,2> fs;

    // Set function spaces for velocity and pressure.
    fs::get_thood_fs(com_mod, fs, lM, lStab, 1);

    // Define element coordinates appropriate for function spaces
    Array<double> xwl(nsd,fs[0].eNoN); 
    Array<double> Nwx(nsd,fs[0].eNoN); 

    Array<double> xql(nsd,fs[1].eNoN); 
    Array<double> Nqx(nsd,fs[1].eNoN);

    xwl = xl;

    for (int i = 0; i < nsd; i++) { 
      for (int j = 0; j < fs[1].eNoN; j++) { 
        xql(i,j) = xl(i,j);
      }
    }

    // Gauss integration 1
    //
    double Jac{0.0};
    Array<double> ksix(nsd,nsd);

    for (int g = 0; g < fs[0].nG; g++) {
      if (g == 0 || !fs[0].lShpF) {
        auto Nx = fs[0].Nx.slice(g);
        nn::gnn(fs[0].eNoN, nsd, nsd, Nx, xwl, Nwx, Jac, ksix);
        if (utils::is_zero(Jac)) {
           throw std::runtime_error("[construct_stokes] Jacobian for element " + std::to_string(e) + " is < 0.");
        }
      }

      double w = fs[0].w(g) * Jac;

      // Compute momentum residual and tangent matrix.
      //
      if (nsd == 3) {
        auto N0 = fs[0].N.col(g); 
        auto N1 = fs[1].N.col(g); 
        stokes_3d_m(com_mod, fs[0].eNoN, fs[1].eNoN, w, N0, N1, Nwx, al, yl, bfl, lR, lK);

      } else if (nsd == 2) {
        auto N0 = fs[0].N.col(g); 
        auto N1 = fs[1].N.col(g); 
        stokes_2d_m(com_mod, fs[0].eNoN, fs[1].eNoN, w, N0, N1, Nwx, al, yl, bfl, lR, lK);
      }

    } // g: loop

    // Set function spaces for velocity and pressure.
    //
    fs::get_thood_fs(com_mod, fs, lM, lStab, 2);

    // Gauss integration 2
    //
    for (int g = 0; g < fs[1].nG; g++) {
      if (g == 0 || !fs[0].lShpF) {
        auto Nx = fs[0].Nx.slice(g);
        nn::gnn(fs[0].eNoN, nsd, nsd, Nx, xwl, Nwx, Jac, ksix);

        if (utils::is_zero(Jac)) {
           throw std::runtime_error("[construct_stokes] Jacobian for element " + std::to_string(e) + " is < 0.");
        }
      }

      if (g == 0 || !fs[1].lShpF) {
        auto Nx = fs[1].Nx.slice(g);
        nn::gnn(fs[1].eNoN, nsd, nsd, Nx, xql, Nqx, Jac, ksix);

        if (utils::is_zero(Jac)) {
           throw std::runtime_error("[construct_stokes] Jacobian for element " + std::to_string(e) + " is < 0.");
        }
      }

      double w = fs[1].w(g) * Jac;

      // Compute continuity residual and tangent matrix.
      //
      if (nsd == 3) {
        auto N0 = fs[0].N.col(g); 
        auto N1 = fs[1].N.col(g); 
        stokes_3d_c(com_mod, lStab, fs[0].eNoN, fs[1].eNoN, w, ksix, N0, N1, Nwx, Nqx, al, yl, bfl, lR, lK);

      } else if (nsd == 2) {
        auto N0 = fs[0].N.col(g); 
        auto N1 = fs[1].N.col(g); 
        stokes_2d_c(com_mod, lStab, fs[0].eNoN, fs[1].eNoN, w, ksix, N0, N1, Nwx, Nqx,  al, yl, bfl, lR, lK);
      }

    } // g: loop

    eq.linear_algebra->assemble(com_mod, eNoN, ptr, lK, lR);

  } // e: loop

}

/// @brief Reproduces Fortran 'STOKES2D_C()'.
//
void stokes_2d_c(ComMod& com_mod, const int lStab, const int eNoNw, const int eNoNq, const double w, 
    const Array<double>& ksix, const Vector<double>& Nw, const Vector<double>& Nq, const Array<double>& Nwx, 
    const Array<double>& Nqx, const Array<double>& al, const Array<double>& yl, const Array<double>& bfl, 
    Array<double>& lR, Array3<double>& lK)
{
  using namespace consts;

  #define n_debug_stokes_2d_c
  #ifdef debug_stokes_2d_c
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "w: " << w;
  dmsg << "eNoNw: " << eNoNw;
  dmsg << "eNoNq: " << eNoNq;
  #endif

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;
  double mu = dmn.fluid_visc.mu_i;
  double ctM = dmn.prop[PhysicalProperyType::ctau_M];

  Vector<double> fb(2);
  fb[0] = dmn.prop[PhysicalProperyType::f_x];
  fb[1] = dmn.prop[PhysicalProperyType::f_y];
  double wm = w * eq.am;
  double wf = w * eq.af * eq.gam * dt;

  // {i,j} := velocity dofs; {k} := pressure dof
  int i = eq.s;
  int j = i + 1;
  int k = j + 1;

  #ifdef debug_stokes_2d_c
  dmsg << "dt: " << dt;
  dmsg << "mu: " << mu;
  dmsg << "wf: " << wf;
  dmsg << "wm: " << wm;
  dmsg << "ctM: " << ctM;
  dmsg << "fb: " << fb;
  dmsg;
  #endif

  // Note that indices are not selected based on the equation because
  // fluid equation always come first
  // Velocity and its gradients, inertia (acceleration & body force)
  //
  Vector<double> vd{-fb[0], -fb[1]};
  Vector<double> v(2);
  Array<double> vx(2,2);

  double div = 0.0;

  for (int a = 0; a < eNoNw; a++) {
    vd(0) = vd(0) + Nw(a)*(al(i,a)-bfl(0,a));
    vd(1) = vd(1) + Nw(a)*(al(j,a)-bfl(1,a));
    div = div + Nwx(0,a)*yl(i,a) + Nwx(1,a)*yl(j,a);
  }

  // Pressure and its gradient
  Vector<double> px(2);
  for (int a = 0; a < eNoNq; a++) {
    px(0) = px(0) + Nqx(0,a)*yl(k,a);
    px(1) = px(1) + Nqx(1,a)*yl(k,a);
  }

  double kS = ksix(0,0)*ksix(0,0) + ksix(1,0)*ksix(1,0) + ksix(0,1)*ksix(0,1) + ksix(1,1)*ksix(1,1);
  kS = sqrt(kS);

  double tauM{0.0};

  if (lStab) {
    tauM = ctM / (2.0*mu*kS);
  } else {
    tauM = 0.0;
  }

  // Local residual
  //
  auto rMv = vd + px;

  #ifdef debug_stokes_2d_c
  dmsg << "div: " << div;
  dmsg << "px: " << px;
  dmsg << "kS: " << kS;
  dmsg << "tauM: " << tauM;
  dmsg << "rMv: " << rMv;
  #endif

  for (int a = 0; a < eNoNq; a++) {
    double rM = rMv(0)*Nqx(0,a) + rMv(1)*Nqx(1,a);
    lR(2,a) = lR(2,a) + w*(Nq(a)*div + tauM*rM);
  }

  // Tangent (stiffness) matrices
  //
  wm = wm * tauM;

  for (int b = 0; b < eNoNw; b++) {
    for (int a = 0; a < eNoNq; a++) {
      lK(6,a,b) = lK(6,a,b) + wm*Nqx(0,a)*Nw(b) + wf*Nq(a)*Nwx(0,b);
      lK(7,a,b) = lK(7,a,b) + wm*Nqx(1,a)*Nw(b) + wf*Nq(a)*Nwx(1,b);
    }
  }

  if (lStab) {
    wf = wf * tauM;

    for (int b = 0; b < eNoNq; b++) {
      for (int a = 0; a < eNoNq; a++) {
        double NxNx = Nqx(0,a)*Nqx(0,b) + Nqx(1,a)*Nqx(1,b);
        lK(8,a,b) = lK(8,a,b) + wf*NxNx;
      }
    }
  }
}

//-------------
// stokes_2d_m
//-------------
//
void stokes_2d_m(ComMod& com_mod, const int eNoNw, const int eNoNq, const double w, 
    const Vector<double>& Nw, const Vector<double>& Nq, const Array<double>& Nwx, const Array<double>& al, 
    const Array<double>& yl, const Array<double>& bfl, Array<double>& lR, Array3<double>& lK)
{
  using namespace consts;

  #define n_debug_stokes_2d_m
  #ifdef debug_stokes_2d_m
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "w: " << w;
  dmsg << "eNoNw: " << eNoNw;
  dmsg << "eNoNq: " << eNoNq;
  #endif

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  double mu = dmn.fluid_visc.mu_i;
  Vector<double> fb(2);
  fb[0] = dmn.prop[PhysicalProperyType::f_x];
  fb[1] = dmn.prop[PhysicalProperyType::f_y];

  double af = eq.af * eq.gam * dt;
  double wf = w * af;
  double wm = w * eq.am;

  // {i,j} := velocity dofs; {k} := pressure dof
  int i = eq.s;
  int j = i + 1;
  int k = j + 1;

  #ifdef debug_stokes_2d_m
  dmsg << "dt: " << dt;
  dmsg << "mu: " << mu;
  dmsg << "af: " << af;
  dmsg << "fb: " << fb;
  dmsg;
  #endif

  // Note that indices are not selected based on the equation because
  // fluid equation always come first
  // Velocity and its gradients, inertia (acceleration & body force)
  //
  Vector<double> vd{-fb[0], -fb[1]};
  Vector<double> v(2);
  Array<double> vx(2,2);

  for (int a = 0; a < eNoNw; a++) {
    v(0) = v(0) + Nw(a)*yl(i,a);
    v(1) = v(1) + Nw(a)*yl(j,a);

    vx(0,0) = vx(0,0) + Nwx(0,a)*yl(i,a);
    vx(0,1) = vx(0,1) + Nwx(1,a)*yl(i,a);

    vx(1,0) = vx(1,0) + Nwx(0,a)*yl(j,a);
    vx(1,1) = vx(1,1) + Nwx(1,a)*yl(j,a);

    vd(0) = vd(0) + Nw(a)*(al(i,a)-bfl(0,a));
    vd(1) = vd(1) + Nw(a)*(al(j,a)-bfl(1,a));
  }

  // Pressure 
  double p = 0.0;
  for (int a = 0; a < eNoNq; a++) {
    p = p + Nq(a)*yl(k,a);
  }

  // Viscous stress tensor 
  Array<double> es(2,2);
  es(0,0) = mu * (vx(0,0) + vx(0,0));
  es(1,1) = mu * (vx(1,1) + vx(1,1));
  es(0,1) = mu * (vx(0,1) + vx(1,0));
  es(1,0) = mu * es(0,1);

  for (int a = 0; a < eNoNw; a++) {
    double rM = Nwx(0,a)*(-p + es(0,0)) + Nwx(1,a)*es(0,1);
    lR(0,a) = lR(0,a) + w*(Nw(a)*vd(0) + rM);
    rM = Nwx(0,a)*es(1,0) + Nwx(1,a)*(-p + es(1,1));
    lR(1,a) = lR(1,a) + w*(Nw(a)*vd(1) + rM);
  }

  // Tangent (stiffness) matrices
  for (int b = 0; b < eNoNw; b++) {
    for (int a = 0; a < eNoNw; a++) {
      double NxNx = Nwx(0,a)*Nwx(0,b) + Nwx(1,a)*Nwx(1,b);

      lK(0,a,b) = lK(0,a,b) + wf*mu*(Nwx(0,a)*Nwx(0,b) + NxNx) + wm*Nw(a)*Nw(b);
      lK(1,a,b) = lK(1,a,b) + wf*mu*(Nwx(1,a)*Nwx(0,b));
      lK(3,a,b) = lK(3,a,b) + wf*mu*(Nwx(0,a)*Nwx(1,b));
      lK(4,a,b) = lK(4,a,b) + wf*mu*(Nwx(1,a)*Nwx(1,b) + NxNx) + wm*Nw(a)*Nw(b);
    }
  }

  for (int b = 0; b < eNoNq; b++) {
    for (int a = 0; a < eNoNw; a++) {
      lK(2,a,b) = lK(2,a,b) - wf*Nwx(0,a)*Nq(b);
      lK(5,a,b) = lK(5,a,b) - wf*Nwx(1,a)*Nq(b);
    }
  }
}

//-------------
// stokes_3d_c
//-------------
// Element continuity residual.
//
void stokes_3d_c(ComMod& com_mod, const int lStab, const int eNoNw, const int eNoNq, const double w, 
    const Array<double>& ksix, const Vector<double>& Nw, const Vector<double>& Nq, const Array<double>& Nwx, 
    const Array<double>& Nqx, const Array<double>& al, const Array<double>& yl, 
    const Array<double>& bfl, Array<double>& lR, Array3<double>& lK)
{
  #define n_debug_stokes3d_c
  #ifdef debug_stokes3d_c
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "vmsFlag: " << vmsFlag;
  dmsg << "eNoNw: " << eNoNw;
  dmsg << "eNoNq: " << eNoNq;
  #endif

  using namespace consts;

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;
  double mu = dmn.fluid_visc.mu_i;
  double ctM = dmn.prop[PhysicalProperyType::ctau_M];

  Vector<double> fb(3);
  fb[0] = dmn.prop[PhysicalProperyType::f_x];
  fb[1] = dmn.prop[PhysicalProperyType::f_y];
  fb[2] = dmn.prop[PhysicalProperyType::f_z];

  double wm = w * eq.am;
  double wf = w * eq.af * eq.gam * dt;

  // {i,j} := velocity dofs; {k} := pressure dof
  int i = eq.s;
  int j = i + 1;
  int k = j + 1;
  int l = k + 1;

  // Note that indices are not selected based on the equation because
  // fluid equation always come first
  // Velocity and its gradients, inertia (acceleration & body force)
  //
  Vector<double> vd{-fb[0], -fb[1], -fb[2]};
  Vector<double> v(3);
  Array<double> vx(3,3);

  double div = 0.0;

  for (int a = 0; a < eNoNw; a++) {
    vd(0) = vd(0) + Nw(a)*(al(i,a)-bfl(0,a));
    vd(1) = vd(1) + Nw(a)*(al(j,a)-bfl(1,a));
    vd(2) = vd(2) + Nw(a)*(al(k,a)-bfl(2,a));
    div = div + Nwx(0,a)*yl(i,a) + Nwx(1,a)*yl(j,a) + Nwx(2,a)*yl(k,a);
  }

  // Pressure and its gradient
  Vector<double> px(3);
  for (int a = 0; a < eNoNq; a++) {
    px(0) = px(0) + Nqx(0,a)*yl(l,a);
    px(1) = px(1) + Nqx(1,a)*yl(l,a);
    px(2) = px(2) + Nqx(2,a)*yl(l,a);
  }

  double kS = ksix(0,0)*ksix(0,0) + ksix(1,0)*ksix(1,0) + ksix(2,0)*ksix(2,0) + ksix(0,1)*ksix(0,1) + 
      ksix(1,1)*ksix(1,1) + ksix(2,1)*ksix(2,1) + ksix(0,2)*ksix(0,2) + ksix(1,2)*ksix(1,2) + 
      ksix(2,2)*ksix(2,2);
  kS = sqrt(kS);

  double tauM{0.0};

  if (lStab) {
    tauM = ctM / (2.0*mu*kS);
  } else {
    tauM = 0.0;
  }

  // Local residual
  //
  auto rMv = vd + px;

  for (int a = 0; a < eNoNq; a++) {
    double rM = rMv(0)*Nqx(0,a) + rMv(1)*Nqx(1,a) + rMv(2)*Nqx(2,a);
    lR(3,a) = lR(3,a) + w*(Nq(a)*div + tauM*rM);
  }

  // Tangent (stiffness) matrices
  //
  wm = wm * tauM;

  for (int b = 0; b < eNoNw; b++) {
    for (int a = 0; a < eNoNq; a++) {
      lK(13,a,b) = lK(13,a,b) + wm*Nqx(0,a)*Nw(b) + wf*Nq(a)*Nwx(0,b);
      lK(14,a,b) = lK(14,a,b) + wm*Nqx(1,a)*Nw(b) + wf*Nq(a)*Nwx(1,b);
      lK(15,a,b) = lK(15,a,b) + wm*Nqx(2,a)*Nw(b) + wf*Nq(a)*Nwx(2,b);
    }
  }

  if (lStab) {
    wf = wf * tauM;

    for (int b = 0; b < eNoNq; b++) {
      for (int a = 0; a < eNoNq; a++) {
        double NxNx = Nqx(0,a)*Nqx(0,b) + Nqx(1,a)*Nqx(1,b) + Nqx(2,a)*Nqx(2,b);
        lK(16,a,b) = lK(16,a,b) + wf*NxNx;
      }
    }
  }
}

/// @brief Element momentum residual
///
///  Modifies:
///    lR(dof,eNoN)  - Residual
///    lK(dof*dof,eNoN,eNoN) - Stiffness matrix
//
void stokes_3d_m(ComMod& com_mod, const int eNoNw, const int eNoNq, const double w, 
    const Vector<double>& Nw, const Vector<double>& Nq, const Array<double>& Nwx, 
    const Array<double>& al, const Array<double>& yl, const Array<double>& bfl, 
    Array<double>& lR, Array3<double>& lK)
{
  using namespace consts;

  #define n_debug_stokes_3d_m
  #ifdef debug_stokes_3d_m
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "w: " << w;
  dmsg << "eNoNw: " << eNoNw;
  dmsg << "eNoNq: " << eNoNq;
  #endif

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  double mu = dmn.fluid_visc.mu_i;
  Vector<double> fb(3);
  fb[0] = dmn.prop[PhysicalProperyType::f_x];
  fb[1] = dmn.prop[PhysicalProperyType::f_y];
  fb[2] = dmn.prop[PhysicalProperyType::f_z];

  double af = eq.af * eq.gam * dt;
  double wf = w * af;
  double wm = w * eq.am;

  // {i,j} := velocity dofs; {k} := pressure dof
  int i = eq.s;
  int j = i + 1;
  int k = j + 1;
  int l = k + 1;

  #ifdef debug_stokes_3d_m
  dmsg << "dt: " << dt;
  dmsg << "mu: " << mu;
  dmsg << "af: " << af;
  dmsg << "fb: " << fb;
  #endif

  // Note that indices are not selected based on the equation because
  // fluid equation always come first
  // Velocity and its gradients, inertia (acceleration & body force)
  //
  Vector<double> vd{-fb[0], -fb[1], -fb[2]};
  Vector<double> v(3);
  Array<double> vx(3,3);

  for (int a = 0; a < eNoNw; a++) {
    v(0) = v(0) + Nw(a)*yl(i,a);
    v(1) = v(1) + Nw(a)*yl(j,a);
    v(2) = v(2) + Nw(a)*yl(k,a);

    vx(0,0) = vx(0,0) + Nwx(0,a)*yl(i,a);
    vx(0,1) = vx(0,1) + Nwx(1,a)*yl(i,a);
    vx(0,2) = vx(0,2) + Nwx(2,a)*yl(i,a);

    vx(1,0) = vx(1,0) + Nwx(0,a)*yl(j,a);
    vx(1,1) = vx(1,1) + Nwx(1,a)*yl(j,a);
    vx(1,2) = vx(1,2) + Nwx(2,a)*yl(j,a);

    vx(2,0) = vx(2,0) + Nwx(0,a)*yl(k,a);
    vx(2,1) = vx(2,1) + Nwx(1,a)*yl(k,a);
    vx(2,2) = vx(2,2) + Nwx(2,a)*yl(k,a);

    vd(0) = vd(0) + Nw(a)*(yl(i,a)-bfl(0,a));
    vd(1) = vd(1) + Nw(a)*(yl(j,a)-bfl(1,a));
    vd(2) = vd(2) + Nw(a)*(yl(k,a)-bfl(2,a));
  }

  // Pressure 
  double p = 0.0;
  for (int a = 0; a < eNoNq; a++) {
    p = p + Nq(a)*yl(l,a);
  }

  // Viscous stress tensor 
  Array<double> es(3,3);
  es(0,0) = mu*(vx(0,0) + vx(0,0));
  es(1,1) = mu*(vx(1,1) + vx(1,1));
  es(2,2) = mu*(vx(2,2) + vx(2,2));

  es(0,1) = mu*(vx(0,1) + vx(1,0));
  es(0,2) = mu*(vx(0,2) + vx(2,0));
  es(1,2) = mu*(vx(1,2) + vx(2,1));

  es(1,0) = es(0,1);
  es(2,0) = es(0,2);
  es(2,1) = es(1,2);

  // Local residual
  for (int a = 0; a < eNoNw; a++) {
    double rM = Nwx(0,a)*(-p + es(0,0)) + Nwx(1,a)*es(0,1) + Nwx(2,a)*es(0,2);
    lR(0,a) = lR(0,a) + w*(Nw(a)*vd(0) + rM);

    rM = Nwx(0,a)*es(1,0) + Nwx(1,a)*(-p + es(1,1)) + Nwx(2,a)*es(1,2);
    lR(1,a) = lR(1,a) + w*(Nw(a)*vd(1) + rM);

    rM = Nwx(0,a)*es(2,0) + Nwx(1,a)*es(2,1) + Nwx(2,a)*(-p + es(2,2));
    lR(2,a) = lR(2,a) + w*(Nw(a)*vd(2) + rM);
  }

  // Tangent (stiffness) matrices
  //
  for (int b = 0; b < eNoNw; b++) {
    for (int a = 0; a < eNoNw; a++) {
      double NxNx = Nwx(0,a)*Nwx(0,b) + Nwx(1,a)*Nwx(1,b) +  Nwx(2,a)*Nwx(2,b);

      lK(0,a,b) = lK(0,a,b) + wf*mu*(Nwx(0,a)*Nwx(0,b) + NxNx) + wm*Nw(a)*Nw(b);
      lK(1,a,b) = lK(1,a,b) + wf*mu*(Nwx(1,a)*Nwx(0,b));
      lK(2,a,b) = lK(2,a,b) + wf*mu*(Nwx(2,a)*Nwx(0,b));

      lK(4,a,b) = lK(4,a,b) + wf*mu*(Nwx(0,a)*Nwx(1,b));
      lK(5,a,b) = lK(5,a,b) + wf*mu*(Nwx(1,a)*Nwx(1,b) + NxNx) + wm*Nw(a)*Nw(b);
      lK(6,a,b) = lK(6,a,b) + wf*mu*(Nwx(2,a)*Nwx(1,b));

      lK(8,a,b) = lK(8,a,b)   + wf*mu*(Nwx(0,a)*Nwx(2,b));
      lK(9,a,b) = lK(9,a,b)   + wf*mu*(Nwx(1,a)*Nwx(2,b));
      lK(10,a,b) = lK(10,a,b) + wf*mu*(Nwx(2,a)*Nwx(2,b) + NxNx) + wm*Nw(a)*Nw(b);
    }
  }

  for (int b = 0; b < eNoNq; b++) {
    for (int a = 0; a < eNoNw; a++) {
      lK(3,a,b) = lK(3,a,b) - wf*Nwx(0,a)*Nq(b);
      lK(7,a,b) = lK(7,a,b) - wf*Nwx(1,a)*Nq(b);
      lK(11,a,b) = lK(11,a,b) - wf*Nwx(2,a)*Nq(b);
    }
  }
}

};
