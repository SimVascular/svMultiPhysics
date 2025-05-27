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

// This is for solving linear elasticity.
//
// Reproduces Fortran code in LELAS.f.

#include "l_elas.h"

#include "all_fun.h"
#include "lhsa.h"
#include "nn.h"
#include "utils.h"

namespace l_elas {

/// @brief This is for boundary elasticity equation
///
/// Reproduces 'SUBROUTINE BLELAS (eNoN, w, N, h, nV, lR)'
//
void b_l_elas(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const double h, 
    const Vector<double>& nV, Array<double>& lR)
{
  const int nsd = com_mod.nsd;
  auto hc = h*nV;

  for (int a = 0; a < eNoN; a++) {
    for (int i = 0; i < nsd; i++) {
      lR(i,a) = lR(i,a) - w*N(a)*hc(i);
    }
  }
}

/// @brief Reproduces Fortran 'CONSTRUCT_LELAS'.
//
void construct_l_elas(ComMod& com_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Dg)
{
  using namespace consts;

  #define debug_construct_l_elas 
  #ifdef debug_construct_l_elas 
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

  int eNoN = lM.eNoN;

  Vector<int> ptr(eNoN);
  Vector<double> pSl(nsymd), N(eNoN);
  Array<double> xl(nsd,eNoN), al(tDof,eNoN), dl(tDof,eNoN), bfl(nsd,eNoN), 
      pS0l(nsymd,eNoN), Nx(nsd,eNoN), lR(dof,eNoN);
  Array3<double> lK(dof*dof,eNoN,eNoN);

  for (int e = 0; e < lM.nEl; e++) {
    // Change the current domain which will be used in later function calls. 
    cDmn = all_fun::domain(com_mod, lM, cEq, e);
    auto cPhys = eq.dmn[cDmn].phys;
    if (cPhys != EquationType::phys_lElas) {
      continue;
    }

    // Update shape functions for NURBS
    if (lM.eType == ElementType::NRB) {
      //CALL NRBNNX(lM, e)
    }

    // Create local copies
    pS0l = 0.0;

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
      }

      if (pS0.size() != 0) {
        pS0l.set_col(a, pS0.col(Ac));
      }
    }

    // Gauss integration
    //
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
      pSl = 0.0;

      if (nsd == 3) {
        l_elas_3d(com_mod, eNoN, w, N, Nx, al, dl, bfl, pS0l, pSl, lR, lK);
      } else if (nsd == 2) {
        l_elas_2d(com_mod, eNoN, w, N, Nx, al, dl, bfl, pS0l, pSl, lR, lK);
      }

      // Prestress
      if (pstEq) {
        for (int a = 0; a < eNoN; a++) {
          int Ac = ptr(a);
          pSa(Ac) = pSa(Ac) + w*N(a);
          for (int i = 0; i < pSn.nrows(); i++) {
            pSn(i,Ac) = pSn(i,Ac) + w*N(a)*pSl(i);
          }
        }
      }
    }

    eq.linear_algebra->assemble(com_mod, eNoN, ptr, lK, lR);
  }
}

/// @brief Reproducs Fortran 'LELAS2D()'.
//
void l_elas_2d(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N,
    const Array<double>& Nx, const Array<double>& al, const Array<double>& dl, const Array<double>& bfl,
    const Array<double>& pS0l, Vector<double>& pSl, Array<double>& lR, Array3<double>& lK)
{
  using namespace consts;

  const int dof = com_mod.dof;
  const int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  const auto cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];

  #define n_debug_l_elas_2d 
  #ifdef debug_l_elas_2d 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "w: " << w;
  dmsg << "cEq: " << cEq;
  dmsg << "cDmn: " << cDmn;
  dmsg << "dmn.phys: " << dmn.phys;
  #endif

  const double dt = com_mod.dt;
  double rho = dmn.prop.at(PhysicalProperyType::solid_density);
  double elM = dmn.prop.at(PhysicalProperyType::elasticity_modulus);
  double nu = dmn.prop.at(PhysicalProperyType::poisson_ratio);

  Vector<double> f({dmn.prop.at(PhysicalProperyType::f_x),
                    dmn.prop.at(PhysicalProperyType::f_y),});

  int i = eq.s;
  int j = i + 1;

  double lambda = elM*nu / (1.0+nu) / (1.0-2.0*nu);
  double mu = elM * 0.5 / (1.0+nu);
  double lDm = lambda / mu;
  double T1 = eq.af*eq.beta*dt*dt;
  double amd = eq.am / T1*rho;
  double wl = w * T1 * mu;

  #ifdef debug_l_elas_2d 
  dmsg << "rho: " << rho;
  dmsg << "elM: " << elM;
  dmsg << "nu: " << nu;
  dmsg << "i: " << i;
  dmsg << "j: " << j;
  dmsg << "f: " << f;
  dmsg << "lambda: " << lambda;
  #endif

  Vector<double> ed(3), ud(2), S0(3), S(3);

  ud = -f;

  for (int a = 0; a < eNoN; a++) {
    ud(0) = ud(0) + N(a)*(al(i,a)-bfl(0,a));
    ud(1) = ud(1) + N(a)*(al(j,a)-bfl(1,a));

    ed(0) = ed(0) + Nx(0,a)*dl(i,a);
    ed(1) = ed(1) + Nx(1,a)*dl(j,a);
    ed(2) = ed(2) + Nx(1,a)*dl(i,a) + Nx(0,a)*dl(j,a);

    S0(0) = S0(0) + N(a)*pS0l(0,a);
    S0(1) = S0(1) + N(a)*pS0l(1,a);
    S0(2) = S0(2) + N(a)*pS0l(2,a);
  }

  double divD = lambda * (ed(0) + ed(1));

  // Stress in Voigt notation
  S(0) = divD + 2.0*mu*ed(0);
  S(1) = divD + 2.0*mu*ed(1);
  S(2) = mu*ed(2);  // 2*eps_12
  pSl = S;

  // Add prestress contribution
  S = pSl + S0;

  for (int a = 0; a < eNoN; a++) {
    lR(0,a) = lR(0,a) + w*(rho*N(a)*ud(0) + Nx(0,a)*S(0) + Nx(1,a)*S(2)); 
    lR(1,a) = lR(1,a) + w*(rho*N(a)*ud(1) + Nx(0,a)*S(2) + Nx(1,a)*S(1));

    for (int b = 0; b < eNoN; b++) {
      double NxdNx = Nx(0,a)*Nx(0,b) + Nx(1,a)*Nx(1,b);
      double T1 = amd*N(a)*N(b) / mu + NxdNx;

      lK(0,a,b) = lK(0,a,b) + wl*(T1 + (1.0 + lDm)*Nx(0,a)*Nx(0,b));
      lK(1,a,b) = lK(1,a,b) + wl*(lDm*Nx(0,a)*Nx(1,b) + Nx(1,a)*Nx(0,b));

      lK(dof+0,a,b) = lK(dof+0,a,b) + wl*(lDm*Nx(1,a)*Nx(0,b) + Nx(0,a)*Nx(1,b));
      lK(dof+1,a,b) = lK(dof+1,a,b) + wl*(T1 + (1.0 + lDm)*Nx(1,a)*Nx(1,b));
    }
  }
}

//-----------
// l_elas_3d
//-----------
//
// Reproducs Fortran 'LELAS3D()'.
//
void l_elas_3d(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, 
    const Array<double>& Nx, const Array<double>& al, const Array<double>& dl, const Array<double>& bfl, 
    const Array<double>& pS0l, Vector<double>& pSl, Array<double>& lR, Array3<double>& lK)
{
  using namespace consts;

  const int dof = com_mod.dof;
  const int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  const auto cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];

  #define n_debug_l_elas_3d 
  #ifdef debug_l_elas_3d 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "w: " << w;
  dmsg << "cEq: " << cEq;
  dmsg << "cDmn: " << cDmn;
  dmsg << "dmn.phys: " << dmn.phys;
  #endif

  const double dt = com_mod.dt;
  double rho = dmn.prop.at(PhysicalProperyType::solid_density);
  double elM = dmn.prop.at(PhysicalProperyType::elasticity_modulus);
  double nu = dmn.prop.at(PhysicalProperyType::poisson_ratio);

  Vector<double> f({dmn.prop.at(PhysicalProperyType::f_x),
                    dmn.prop.at(PhysicalProperyType::f_y),
                    dmn.prop.at(PhysicalProperyType::f_z)});

  int i = eq.s;
  int j = i + 1;
  int k = j + 1;

  double lambda = elM*nu / (1.0+nu) / (1.0-2.0*nu);
  double mu = elM * 0.5 / (1.0+nu);
  double lDm = lambda / mu;
  double T1 = eq.af*eq.beta*dt*dt;
  double amd = eq.am / T1*rho;
  double wl = w * T1 * mu;

  #ifdef debug_l_elas_3d 
  dmsg << "rho: " << rho;
  dmsg << "elM: " << elM;
  dmsg << "nu: " << nu;
  dmsg << "i: " << i;
  dmsg << "j: " << j;
  dmsg << "k: " << k;
  dmsg << "f: " << f;
  dmsg << "lambda: " << lambda;
  #endif

  Vector<double> ed(6), ud(3), S0(6), S(6);

  ud = -f;

  for (int a = 0; a < eNoN; a++) {
    ud(0) = ud(0) + N(a)*(al(i,a)-bfl(0,a));
    ud(1) = ud(1) + N(a)*(al(j,a)-bfl(1,a));
    ud(2) = ud(2) + N(a)*(al(k,a)-bfl(2,a));

    ed(0) = ed(0) + Nx(0,a)*dl(i,a);
    ed(1) = ed(1) + Nx(1,a)*dl(j,a);
    ed(2) = ed(2) + Nx(2,a)*dl(k,a);
    ed(3) = ed(3) + Nx(1,a)*dl(i,a) + Nx(0,a)*dl(j,a);
    ed(4) = ed(4) + Nx(2,a)*dl(j,a) + Nx(1,a)*dl(k,a);
    ed(5) = ed(5) + Nx(0,a)*dl(k,a) + Nx(2,a)*dl(i,a);

    S0(0) = S0(0) + N(a)*pS0l(0,a);
    S0(1) = S0(1) + N(a)*pS0l(1,a);
    S0(2) = S0(2) + N(a)*pS0l(2,a);
    S0(3) = S0(3) + N(a)*pS0l(3,a);
    S0(4) = S0(4) + N(a)*pS0l(4,a);
    S0(5) = S0(5) + N(a)*pS0l(5,a);
  }

  double divD = lambda * (ed(0) + ed(1) + ed(2));

  // Stress in Voigt notation
  S(0) = divD + 2.0*mu*ed(0);
  S(1) = divD + 2.0*mu*ed(1);
  S(2) = divD + 2.0*mu*ed(2);
  S(3) = mu*ed(3);  // 2*eps_12
  S(4) = mu*ed(4);  // 2*eps_23
  S(5) = mu*ed(5);  // 2*eps_13
  pSl = S;

  // Add prestress contribution
  S = pSl + S0;

  for (int a = 0; a < eNoN; a++) {
    lR(0,a) = lR(0,a) + w*(rho*N(a)*ud(0) + Nx(0,a)*S(0) + Nx(1,a)*S(3) + Nx(2,a)*S(5));
    lR(1,a) = lR(1,a) + w*(rho*N(a)*ud(1) + Nx(0,a)*S(3) + Nx(1,a)*S(1) + Nx(2,a)*S(4));
    lR(2,a) = lR(2,a) + w*(rho*N(a)*ud(2) + Nx(0,a)*S(5) + Nx(1,a)*S(4) + Nx(2,a)*S(2));

    for (int b = 0; b < eNoN; b++) {
      double NxdNx = Nx(0,a)*Nx(0,b) + Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b);
      double T1 = amd*N(a)*N(b) / mu + NxdNx;

      lK(0,a,b) = lK(0,a,b) + wl*(T1 + (1.0 + lDm)*Nx(0,a)*Nx(0,b));
      lK(1,a,b) = lK(1,a,b) + wl*(lDm*Nx(0,a)*Nx(1,b) + Nx(1,a)*Nx(0,b));
      lK(2,a,b) = lK(2,a,b) + wl*(lDm*Nx(0,a)*Nx(2,b) + Nx(2,a)*Nx(0,b));

      lK(dof+0,a,b) = lK(dof+0,a,b) + wl*(lDm*Nx(1,a)*Nx(0,b) + Nx(0,a)*Nx(1,b));
      lK(dof+1,a,b) = lK(dof+1,a,b) + wl*(T1 + (1.0 + lDm)*Nx(1,a)*Nx(1,b));
      lK(dof+2,a,b) = lK(dof+2,a,b) + wl*(lDm*Nx(1,a)*Nx(2,b) + Nx(2,a)*Nx(1,b));

      lK(2*dof+0,a,b) = lK(2*dof+0,a,b) + wl*(lDm*Nx(2,a)*Nx(0,b) + Nx(0,a)*Nx(2,b));
      lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + wl*(lDm*Nx(2,a)*Nx(1,b) + Nx(1,a)*Nx(2,b));
      lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + wl*(T1 + (1.0 + lDm)*Nx(2,a)*Nx(2,b));
    } 
  } 
}

};

