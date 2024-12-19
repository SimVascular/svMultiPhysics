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

// Functions for solving nonlinear structural mechanics
// problems (pure displacement-based formulation).
//
// Replicates the Fortran functions in 'STRUCT.f'. 

#include "sv_struct.h"

#include "all_fun.h"
#include "consts.h"
#include "lhsa.h"
#include "mat_fun.h"
#include "mat_models.h"
#include "nn.h"
#include "utils.h"

namespace struct_ns {

void b_struct_2d(const ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, 
    const Array<double>& Nx, const Array<double>& dl, const Vector<double>& hl, const Vector<double>& nV, 
    Array<double>& lR, Array3<double>& lK)
{
  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  double dt = com_mod.dt;
  int dof = com_mod.dof;

  double af = eq.af * eq.gam*dt;
  double afm = af / eq.am;
  int i = eq.s;
  int j = i + 1;

  Vector<double> nFi(2);  
  Array<double> NxFi(2,eNoN);

  Array<double> F(2,2); 
  F(0,0) = 1.0;
  F(1,1) = 1.0;

  double h = 0.0;

  for (int a = 0; a  < eNoN; a++) {
    h  = h + N(a)*hl(a);
    F(0,0) = F(0,0) + Nx(0,a)*dl(i,a);
    F(0,1) = F(0,1) + Nx(1,a)*dl(i,a);
    F(1,0) = F(1,0) + Nx(0,a)*dl(j,a);
    F(1,1) = F(1,1) + Nx(1,a)*dl(j,a);
  }

  double Jac = F(0,0)*F(1,1) - F(0,1)*F(1,0);
  auto Fi = mat_fun::mat_inv(F, 2);

  for (int a = 0; a  < eNoN; a++) {
    NxFi(0,a) = Nx(0,a)*Fi(0,0) + Nx(1,a)*Fi(1,0);
    NxFi(1,a) = Nx(0,a)*Fi(0,1) + Nx(1,a)*Fi(1,1);
  }

  nFi(0) = nV(0)*Fi(0,0) + nV(1)*Fi(1,0);
  nFi(1) = nV(0)*Fi(0,1) + nV(1)*Fi(1,1);
  double wl = w * Jac * h;

  for (int a = 0; a  < eNoN; a++) {
    lR(0,a) = lR(0,a) - wl*N(a)*nFi(0);
    lR(1,a) = lR(1,a) - wl*N(a)*nFi(1);

    for (int b = 0; b < eNoN; b++) {
      double Ku = wl*af*N(a)*(nFi(1)*NxFi(0,b) - nFi(0)*NxFi(1,b));
      lK(1,a,b) = lK(1,a,b) + Ku;
      lK(dof,a,b) = lK(dof,a,b) - Ku;
    }
  }
}

/// @brief Add follower pressure load contributions to the local residual and stiffness matrix.
/// @param com_mod 
/// @param eNoN 
/// @param w  Gauss point weight times reference configuration area
/// @param N  Shape function values at the Gauss point
/// @param Nx Shape function derivatives at the Gauss point
/// @param dl Displacement vector
/// @param hl Magnitude of pressure
/// @param nV Normal vector (in reference configuration)
/// @param lR Local residual
/// @param lK Local stiffness matrix
void b_struct_3d(const ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, 
    const Array<double>& Nx, const Array<double>& dl, const Vector<double>& hl, const Vector<double>& nV, 
    Array<double>& lR, Array3<double>& lK)
{
  #define n_debug_b_struct_3d 
  #ifdef debug_b_struct_3d 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  double dt = com_mod.dt;
  int dof = com_mod.dof;

  double af = eq.af * eq.beta*dt*dt;
  int i = eq.s;
  int j = i + 1;
  int k = j + 1;

  #ifdef debug_b_struct_3d 
  debug << "af: " << af;
  debug << "i: " << i;
  debug << "j: " << j;
  debug << "k: " << k;
  #endif

  Vector<double> nFi(3);  
  Array<double> NxFi(3,eNoN);

  Array<double> F(3,3); 
  F(0,0) = 1.0;
  F(1,1) = 1.0;
  F(2,2) = 1.0;

  double h = 0.0;

  // Compute deformation gradient tensor F
  for (int a = 0; a  < eNoN; a++) {
    h  = h + N(a)*hl(a);
    F(0,0) = F(0,0) + Nx(0,a)*dl(i,a);
    F(0,1) = F(0,1) + Nx(1,a)*dl(i,a);
    F(0,2) = F(0,2) + Nx(2,a)*dl(i,a);
    F(1,0) = F(1,0) + Nx(0,a)*dl(j,a);
    F(1,1) = F(1,1) + Nx(1,a)*dl(j,a);
    F(1,2) = F(1,2) + Nx(2,a)*dl(j,a);
    F(2,0) = F(2,0) + Nx(0,a)*dl(k,a);
    F(2,1) = F(2,1) + Nx(1,a)*dl(k,a);
    F(2,2) = F(2,2) + Nx(2,a)*dl(k,a);
  }

  double Jac = mat_fun::mat_det(F, 3);
  auto Fi = mat_fun::mat_inv(F, 3);

  for (int a = 0; a  < eNoN; a++) {
    NxFi(0,a) = Nx(0,a)*Fi(0,0) + Nx(1,a)*Fi(1,0) + Nx(2,a)*Fi(2,0);
    NxFi(1,a) = Nx(0,a)*Fi(0,1) + Nx(1,a)*Fi(1,1) + Nx(2,a)*Fi(2,1);
    NxFi(2,a) = Nx(0,a)*Fi(0,2) + Nx(1,a)*Fi(1,2) + Nx(2,a)*Fi(2,2);
  }
  // Compute N.F^-1, used for Nanson' formula da.n = J*dA*N.F^-1
  nFi(0) = nV(0)*Fi(0,0) + nV(1)*Fi(1,0) + nV(2)*Fi(2,0);
  nFi(1) = nV(0)*Fi(0,1) + nV(1)*Fi(1,1) + nV(2)*Fi(2,1);
  nFi(2) = nV(0)*Fi(0,2) + nV(1)*Fi(1,2) + nV(2)*Fi(2,2);

  double wl = w * Jac * h;

  #ifdef debug_b_struct_3d 
  debug;
  debug << "Jac: " << Jac;
  debug << "h: " << h;
  debug << "wl: " << wl;
  #endif

  // Compute follower pressure load contributions to the local residual and stiffness matrix
  for (int a = 0; a  < eNoN; a++) {
    lR(0,a) = lR(0,a) - wl*N(a)*nFi(0);
    lR(1,a) = lR(1,a) - wl*N(a)*nFi(1);
    lR(2,a) = lR(2,a) - wl*N(a)*nFi(2);

    for (int b = 0; b < eNoN; b++) {
      double Ku = wl * af * N(a) * (nFi(1)*NxFi(0,b) - nFi(0)*NxFi(1,b));
      lK(1,a,b) = lK(1,a,b) + Ku;
      lK(dof,a,b) = lK(dof,a,b) - Ku;

      Ku = wl*af*N(a)*(nFi(2)*NxFi(0,b) - nFi(0)*NxFi(2,b));
      lK(2,a,b) = lK(2,a,b) + Ku;
      lK(2*dof,a,b) = lK(2*dof,a,b) - Ku;

      Ku = wl*af*N(a)*(nFi(2)*NxFi(1,b) - nFi(1)*NxFi(2,b));
      lK(dof+2,a,b) = lK(dof+2,a,b) + Ku;
      lK(2*dof+1,a,b) = lK(2*dof+1,a,b) - Ku;
    }
  }
}

/// @brief Replicates the Fortan 'CONSTRUCT_dSOLID' subroutine.
//
void construct_dsolid(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Yg, const Array<double>& Dg)
{
  using namespace consts;

  #define n_debug_construct_dsolid
  #ifdef debug_construct_dsolid
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  auto& cem = cep_mod.cem;
  const int nsd  = com_mod.nsd;
  const int tDof = com_mod.tDof;
  const int dof = com_mod.dof;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  auto cDmn = com_mod.cDmn;
  const int nsymd = com_mod.nsymd;
  auto& pS0 = com_mod.pS0;
  auto& pSn = com_mod.pSn;
  auto& pSa = com_mod.pSa;
  bool pstEq = com_mod.pstEq;

  int eNoN = lM.eNoN;
  int nFn = lM.nFn;
  if (nFn == 0) {
    nFn = 1;
  }

  #ifdef debug_construct_dsolid
  dmsg << "lM.nEl: " << lM.nEl;
  dmsg << "eNoN: " << eNoN;
  dmsg << "nsymd: " << nsymd;
  dmsg << "nFn: " << nFn;
  dmsg << "lM.nG: " << lM.nG;
  #endif

  // STRUCT: dof = nsd

  Vector<int> ptr(eNoN);
  Vector<double> pSl(nsymd), ya_l(eNoN), N(eNoN);
  Array<double> xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN), 
                bfl(nsd,eNoN), fN(nsd,nFn), pS0l(nsymd,eNoN), Nx(nsd,eNoN), lR(dof,eNoN);
  Array3<double> lK(dof*dof,eNoN,eNoN);

  // Loop over all elements of mesh

  for (int e = 0; e < lM.nEl; e++) {
    // Update domain and proceed if domain phys and eqn phys match
    cDmn = all_fun::domain(com_mod, lM, cEq, e);
    auto cPhys = eq.dmn[cDmn].phys;
    if (cPhys != EquationType::phys_struct) {
      continue; 
    }

    // Update shape functions for NURBS
    if (lM.eType == ElementType::NRB) {
      //CALL NRBNNX(lM, e)
    }

    // Create local copies
    fN  = 0.0;
    pS0l = 0.0;
    ya_l = 0.0;

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

      if (lM.fN.size() != 0) {
        for (int iFn = 0; iFn < nFn; iFn++) {
          for (int i = 0; i < nsd; i++) {
            fN(i,iFn) = lM.fN(i+nsd*iFn,e);
          }
        }
      }

      if (pS0.size() != 0) { 
        pS0l.set_col(a, pS0.col(Ac));
      }

      if (cem.cpld) {
        ya_l(a) = cem.Ya(Ac);
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
        struct_3d(com_mod, cep_mod, eNoN, nFn, w, N, Nx, al, yl, dl, bfl, fN, pS0l, pSl, ya_l, lR, lK);

#if 0
        if (e == 0 && g == 0) {
          Array3<double>::write_enabled = true;
          Array<double>::write_enabled = true;
          lR.write("lR");
          lK.write("lK");
          exit(0);
        }
#endif

      } else if (nsd == 2) {
        struct_2d(com_mod, cep_mod, eNoN, nFn, w, N, Nx, al, yl, dl, bfl, fN, pS0l, pSl, ya_l, lR, lK);
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

/// @brief Reproduces Fortran 'STRUCT2D' subroutine.
//
void struct_2d(ComMod& com_mod, CepMod& cep_mod, const int eNoN, const int nFn, const double w, 
    const Vector<double>& N, const Array<double>& Nx, const Array<double>& al, const Array<double>& yl, 
    const Array<double>& dl, const Array<double>& bfl, const Array<double>& fN, const Array<double>& pS0l, 
    Vector<double>& pSl, const Vector<double>& ya_l, Array<double>& lR, Array3<double>& lK) 
{
  using namespace consts;
  using namespace mat_fun;

  #define n_debug_struct_2d 
  #ifdef debug_struct_2d 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  const int dof = com_mod.dof;
  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  // Set parameters
  //
  double rho = dmn.prop.at(PhysicalProperyType::solid_density);
  double dmp = dmn.prop.at(PhysicalProperyType::damping);
  Vector<double> fb({dmn.prop.at(PhysicalProperyType::f_x), dmn.prop.at(PhysicalProperyType::f_y)});
  double afu = eq.af * eq.beta*dt*dt;
  double afv = eq.af * eq.gam*dt;
  double amd = eq.am * rho  +  eq.af * eq.gam * dt * dmp;
  double afl = eq.af * eq.beta * dt * dt;

  int i = eq.s;
  int j = i + 1;
  #ifdef debug_struct_2d 
  debug << "i: " << i;
  debug << "j: " << j;
  debug << "amd: " << amd;
  debug << "afl: " << afl;
  debug << "w: " << w;
  #endif

  // Inertia, body force and deformation tensor (F)
  //
  Array<double> F(2,2), S0(2,2), vx(2,2);
  Vector<double> ud(2);

  ud = -rho*fb;
  F = 0.0;
  F(0,0) = 1.0;
  F(1,1) = 1.0;
  S0 = 0.0;
  double ya_g = 0.0;

  for (int a = 0; a < eNoN; a++) {
    ud(0) = ud(0) + N(a)*(rho*(al(i,a)-bfl(0,a)) + dmp*yl(i,a));
    ud(1) = ud(1) + N(a)*(rho*(al(j,a)-bfl(1,a)) + dmp*yl(j,a));

    vx(0,0) = vx(0,0) + Nx(0,a)*yl(i,a);
    vx(0,1) = vx(0,1) + Nx(1,a)*yl(i,a);
    vx(1,0) = vx(1,0) + Nx(0,a)*yl(j,a);
    vx(1,1) = vx(1,1) + Nx(1,a)*yl(j,a);

    F(0,0) = F(0,0) + Nx(0,a)*dl(i,a);
    F(0,1) = F(0,1) + Nx(1,a)*dl(i,a);
    F(1,0) = F(1,0) + Nx(0,a)*dl(j,a);
    F(1,1) = F(1,1) + Nx(1,a)*dl(j,a);

    S0(0,0) = S0(0,0) + N(a)*pS0l(0,a);
    S0(1,1) = S0(1,1) + N(a)*pS0l(1,a);
    S0(0,1) = S0(0,1) + N(a)*pS0l(2,a);

    ya_g = ya_g + N(a)*ya_l(a);
  }
  #ifdef debug_struct_2d 
  debug << "ud: " << ud(0) << " " << ud(1);
  debug << "F: " << F(0,0);
  debug << "ya_g: " << ya_g;
  #endif

  S0(1,0) = S0(0,1);

  // 2nd Piola-Kirchhoff stress (S) and material stiffness tensor in Voight notation (Dm)
  Array<double> S(2,2), Dm(3,3);
  double Ja;
  mat_models::get_pk2cc(com_mod, cep_mod, dmn, F, nFn, fN, ya_g, S, Dm, Ja);

  // Viscous 2nd Piola-Kirchhoff stress and tangent contributions
  Array<double> Svis(2,2);
  Array3<double> Kvis_u(4, eNoN, eNoN);
  Array3<double> Kvis_v(4, eNoN, eNoN);

  mat_models::get_visc_stress_and_tangent(dmn, eNoN, Nx, vx, F, Svis, Kvis_u, Kvis_v);

  // Elastic + Viscous stresses
  S = S + Svis;

  // Prestress
  pSl(0) = S(0,0);
  pSl(1) = S(1,1);
  pSl(2) = S(0,1);

  // Total 2nd Piola-Kirchhoff stress
  S = S + S0;

  // 1st Piola-Kirchhoff tensor (P)
  //
  Array<double> P(2,2), DBm(3,2);
  Array3<double> Bm(3,2,eNoN);
  P = mat_fun::mat_mul(F, S);
  #ifdef debug_struct_2d 
  debug << "P: " << P(0,0) << " " << P(0,1);
  debug << "   " << P(1,0) << " " << P(1,1);
  #endif

  // Local residual
  for (int a = 0; a < eNoN; a++) {
    lR(0,a) = lR(0,a) + w*(N(a)*ud(0) + Nx(0,a)*P(0,0) + Nx(2,a)*P(0,1));
    lR(1,a) = lR(1,a) + w*(N(a)*ud(1) + Nx(0,a)*P(1,0) + Nx(1,a)*P(1,1));
  }

  // Auxilary quantities for computing stiffness tensor
  //
  for (int a = 0; a < eNoN; a++) {
    Bm(0,0,a) = Nx(0,a)*F(0,0);
    Bm(0,1,a) = Nx(0,a)*F(1,0);

    Bm(1,0,a) = Nx(1,a)*F(0,1);
    Bm(1,1,a) = Nx(1,a)*F(1,1);

    Bm(2,0,a) = (Nx(0,a)*F(0,1) + F(0,0)*Nx(1,a));
    Bm(2,1,a) = (Nx(0,a)*F(1,1) + F(1,0)*Nx(1,a));
  }

  Array<double> NxFi(2,eNoN), DdNx(2,eNoN), VxNx(2,eNoN);


  //     Local stiffness tensor
  double T1, NxNx, NxSNx, BmDBm;

  for (int b = 0; b < eNoN; b++) { 
    for (int a = 0; a < eNoN; a++) { 

      // Geometric stiffness
      NxSNx = Nx(0,a)*S(0,0)*Nx(0,b) + Nx(1,a)*S(1,0)*Nx(0,b) +
              Nx(0,a)*S(0,1)*Nx(1,b) + Nx(1,a)*S(1,1)*Nx(1,b);
      T1 = amd*N(a)*N(b) + afu*NxSNx;

      // Material stiffness (Bt*D*B)
      DBm(0,0) = Dm(0,0)*Bm(0,0,b) + Dm(0,1)*Bm(1,0,b) + Dm(0,2)*Bm(2,0,b);
      DBm(0,1) = Dm(0,0)*Bm(0,1,b) + Dm(0,1)*Bm(1,1,b) + Dm(0,2)*Bm(2,1,b);

      DBm(1,0) = Dm(1,0)*Bm(0,0,b) + Dm(1,1)*Bm(1,0,b) + Dm(1,2)*Bm(2,0,b);
      DBm(1,1) = Dm(1,0)*Bm(0,1,b) + Dm(1,1)*Bm(1,1,b) + Dm(1,2)*Bm(2,1,b);

      DBm(2,0) = Dm(2,0)*Bm(0,0,b) + Dm(2,1)*Bm(1,0,b) + Dm(2,2)*Bm(2,0,b);
      DBm(2,1) = Dm(2,0)*Bm(0,1,b) + Dm(2,1)*Bm(1,1,b) + Dm(2,2)*Bm(2,1,b);


      // dM1/du1
      // Material stiffness: Bt*D*B
      BmDBm = Bm(0,0,a)*DBm(0,0) + Bm(1,0,a)*DBm(1,0) + Bm(2,0,a)*DBm(2,0);

      lK(0,a,b) = lK(0,a,b) + w*( T1 + afu*(BmDBm + Kvis_u(0,a,b)) + afv*Kvis_v(0,a,b) );

      // dM1/du2
      // Material stiffness: Bt*D*B
      BmDBm = Bm(0,0,a)*DBm(0,1) + Bm(1,0,a)*DBm(1,1) + Bm(2,0,a)*DBm(2,1);

      lK(1,a,b) = lK(1,a,b) + w*( afu*(BmDBm + Kvis_u(1,a,b)) + afv*Kvis_v(1,a,b) );

      // dM2/du1
      // Material stiffness: Bt*D*B
      BmDBm = Bm(0,1,a)*DBm(0,0) + Bm(1,1,a)*DBm(1,0) + Bm(2,1,a)*DBm(2,0);

      lK(dof+0,a,b) = lK(dof+0,a,b) + w*( afu*(BmDBm + Kvis_u(2,a,b)) + afv*Kvis_v(2,a,b) );

      // dM2/du2
      // Material stiffness: Bt*D*B
      BmDBm = Bm(0,1,a)*DBm(0,1) + Bm(1,1,a)*DBm(1,1) + Bm(2,1,a)*DBm(2,1);

      lK(dof+1,a,b) = lK(dof+1,a,b) + w*( T1 + afu*(BmDBm + Kvis_u(3,a,b)) + afv*Kvis_v(3,a,b) );
    }
  }
}

/// @brief Reproduces Fortran 'STRUCT3D' subroutine.
void struct_3d(ComMod& com_mod, CepMod& cep_mod, const int eNoN, const int nFn, const double w, 
    const Vector<double>& N, const Array<double>& Nx, const Array<double>& al, const Array<double>& yl, 
    const Array<double>& dl, const Array<double>& bfl, const Array<double>& fN, const Array<double>& pS0l, 
    Vector<double>& pSl, const Vector<double>& ya_l, Array<double>& lR, Array3<double>& lK) 
{
  using namespace consts;
  using namespace mat_fun;
  // std::cout << "\n==================== struct_3d ===============" << std::endl;

  #define n_debug_struct_3d 
  #ifdef debug_struct_3d 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "eNoN: " << eNoN;
  dmsg << "nFn: " << nFn;
  #endif

  const int dof = com_mod.dof;
  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  // Set parameters
  //
  double rho = dmn.prop.at(PhysicalProperyType::solid_density);
  double dmp = dmn.prop.at(PhysicalProperyType::damping);
  Vector<double> fb({dmn.prop.at(PhysicalProperyType::f_x), 
                     dmn.prop.at(PhysicalProperyType::f_y), 
                     dmn.prop.at(PhysicalProperyType::f_z)});

  double afu = eq.af * eq.beta*dt*dt;
  double afv = eq.af * eq.gam*dt;
  double amd = eq.am * rho  +  eq.af * eq.gam * dt * dmp;

  #ifdef debug_struct_3d 
  dmsg << "rho: " << rho;
  dmsg << "dmp: " << dmp;
  dmsg << "afu: " << afu;
  dmsg << "afv: " << afv;
  dmsg << "amd: " << amd;
  #endif

  int i = eq.s;
  int j = i + 1;
  int k = j + 1;

  // Inertia, body force and deformation tensor (F)
  //
  Array<double> F(3,3), S0(3,3), vx(3,3);
  Vector<double> ud(3);

  double F_f[3][3]={}; 
  F_f[0][0] = 1.0;
  F_f[1][1] = 1.0;
  F_f[2][2] = 1.0;

  ud = -rho*fb;
  F = 0.0;
  F(0,0) = 1.0;
  F(1,1) = 1.0;
  F(2,2) = 1.0;
  S0 = 0.0;
  double ya_g = 0.0;

  for (int a = 0; a < eNoN; a++) {
    ud(0) = ud(0) + N(a)*(rho*(al(i,a)-bfl(0,a)) + dmp*yl(i,a));
    ud(1) = ud(1) + N(a)*(rho*(al(j,a)-bfl(1,a)) + dmp*yl(j,a));
    ud(2) = ud(2) + N(a)*(rho*(al(k,a)-bfl(2,a)) + dmp*yl(k,a));

    vx(0,0) = vx(0,0) + Nx(0,a)*yl(i,a);
    vx(0,1) = vx(0,1) + Nx(1,a)*yl(i,a);
    vx(0,2) = vx(0,2) + Nx(2,a)*yl(i,a);
    vx(1,0) = vx(1,0) + Nx(0,a)*yl(j,a);
    vx(1,1) = vx(1,1) + Nx(1,a)*yl(j,a);
    vx(1,2) = vx(1,2) + Nx(2,a)*yl(j,a);
    vx(2,0) = vx(2,0) + Nx(0,a)*yl(k,a);
    vx(2,1) = vx(2,1) + Nx(1,a)*yl(k,a);
    vx(2,2) = vx(2,2) + Nx(2,a)*yl(k,a);

    F(0,0) = F(0,0) + Nx(0,a)*dl(i,a);
    F(0,1) = F(0,1) + Nx(1,a)*dl(i,a);
    F(0,2) = F(0,2) + Nx(2,a)*dl(i,a);
    F(1,0) = F(1,0) + Nx(0,a)*dl(j,a);
    F(1,1) = F(1,1) + Nx(1,a)*dl(j,a);
    F(1,2) = F(1,2) + Nx(2,a)*dl(j,a);
    F(2,0) = F(2,0) + Nx(0,a)*dl(k,a);
    F(2,1) = F(2,1) + Nx(1,a)*dl(k,a);
    F(2,2) = F(2,2) + Nx(2,a)*dl(k,a);

    S0(0,0) = S0(0,0) + N(a)*pS0l(0,a);
    S0(1,1) = S0(1,1) + N(a)*pS0l(1,a);
    S0(2,2) = S0(2,2) + N(a)*pS0l(2,a);
    S0(0,1) = S0(0,1) + N(a)*pS0l(3,a);
    S0(1,2) = S0(1,2) + N(a)*pS0l(4,a);
    S0(2,0) = S0(2,0) + N(a)*pS0l(5,a);

    ya_g = ya_g + N(a)*ya_l(a);
  }

  S0(1,0) = S0(0,1);
  S0(2,1) = S0(1,2);
  S0(0,2) = S0(2,0);

  // 2nd Piola-Kirchhoff tensor (S) and material stiffness tensor in
  // Voigt notationa (Dm)
  //
  Array<double> S(3,3), Dm(6,6); 
  double Ja;
  mat_models::get_pk2cc(com_mod, cep_mod, dmn, F, nFn, fN, ya_g, S, Dm, Ja);

  // Viscous 2nd Piola-Kirchhoff stress and tangent contributions
  Array<double> Svis(3,3);
  Array3<double> Kvis_u(9, eNoN, eNoN);
  Array3<double> Kvis_v(9, eNoN, eNoN);
  
  mat_models::get_visc_stress_and_tangent(dmn, eNoN, Nx, vx, F, Svis, Kvis_u, Kvis_v);

  // Elastic + Viscous stresses
  S = S + Svis;

  #ifdef debug_struct_3d 
  dmsg << "Jac: " << Jac;
  dmsg << "Fi: " << Fi;
  dmsg << "VxFi: " << VxFi;
  dmsg << "ddev: " << ddev;
  dmsg << "S: " << S;
  #endif

  // Prestress
  pSl(0) = S(0,0);
  pSl(1) = S(1,1);
  pSl(2) = S(2,2);
  pSl(3) = S(0,1);
  pSl(4) = S(1,2);
  pSl(5) = S(2,0);

  // Total 2nd Piola-Kirchhoff stress
  S += S0;

  // 1st Piola-Kirchhoff tensor (P)
  //
  Array<double> P(3,3);
  Array3<double> Bm(6,3,eNoN); 
  mat_fun::mat_mul(F, S, P);

  // Local residual
  for (int a = 0; a < eNoN; a++) {
    lR(0,a) = lR(0,a) + w*(N(a)*ud(0) + Nx(0,a)*P(0,0) + Nx(1,a)*P(0,1) + Nx(2,a)*P(0,2));
    lR(1,a) = lR(1,a) + w*(N(a)*ud(1) + Nx(0,a)*P(1,0) + Nx(1,a)*P(1,1) + Nx(2,a)*P(1,2));
    lR(2,a) = lR(2,a) + w*(N(a)*ud(2) + Nx(0,a)*P(2,0) + Nx(1,a)*P(2,1) + Nx(2,a)*P(2,2));
  }

  // Auxilary quantities for computing stiffness tensor
  //
  for (int a = 0; a < eNoN; a++) {
    Bm(0,0,a) = Nx(0,a)*F(0,0);
    Bm(0,1,a) = Nx(0,a)*F(1,0);
    Bm(0,2,a) = Nx(0,a)*F(2,0);

    Bm(1,0,a) = Nx(1,a)*F(0,1);
    Bm(1,1,a) = Nx(1,a)*F(1,1);
    Bm(1,2,a) = Nx(1,a)*F(2,1);

    Bm(2,0,a) = Nx(2,a)*F(0,2);
    Bm(2,1,a) = Nx(2,a)*F(1,2);
    Bm(2,2,a) = Nx(2,a)*F(2,2);

    Bm(3,0,a) = (Nx(0,a)*F(0,1) + F(0,0)*Nx(1,a));
    Bm(3,1,a) = (Nx(0,a)*F(1,1) + F(1,0)*Nx(1,a));
    Bm(3,2,a) = (Nx(0,a)*F(2,1) + F(2,0)*Nx(1,a));

    Bm(4,0,a) = (Nx(1,a)*F(0,2) + F(0,1)*Nx(2,a));
    Bm(4,1,a) = (Nx(1,a)*F(1,2) + F(1,1)*Nx(2,a));
    Bm(4,2,a) = (Nx(1,a)*F(2,2) + F(2,1)*Nx(2,a));

    Bm(5,0,a) = (Nx(2,a)*F(0,0) + F(0,2)*Nx(0,a));
    Bm(5,1,a) = (Nx(2,a)*F(1,0) + F(1,2)*Nx(0,a));
    Bm(5,2,a) = (Nx(2,a)*F(2,0) + F(2,2)*Nx(0,a));
  }

  // Local stiffness tensor
  double NxSNx, T1, NxNx, BmDBm, Tv;

  Array<double> DBm(6,3);

  for (int b = 0; b < eNoN; b++) {

    for (int a = 0; a < eNoN; a++) {

      // Geometric stiffness
      NxSNx = Nx(0,a)*S(0,0)*Nx(0,b) + Nx(1,a)*S(1,0)*Nx(0,b) +
              Nx(2,a)*S(2,0)*Nx(0,b) + Nx(0,a)*S(0,1)*Nx(1,b) +
              Nx(1,a)*S(1,1)*Nx(1,b) + Nx(2,a)*S(2,1)*Nx(1,b) +
              Nx(0,a)*S(0,2)*Nx(2,b) + Nx(1,a)*S(1,2)*Nx(2,b) +
              Nx(2,a)*S(2,2)*Nx(2,b);

      T1 = amd*N(a)*N(b) + afu*NxSNx;

      // Material Stiffness (Bt*D*B)
      mat_mul(Dm, Bm.rslice(b), DBm);

      // dM1/du1
      // Material stiffness: Bt*D*B
      BmDBm = Bm(0,0,a)*DBm(0,0) + Bm(1,0,a)*DBm(1,0) +
              Bm(2,0,a)*DBm(2,0) + Bm(3,0,a)*DBm(3,0) +
              Bm(4,0,a)*DBm(4,0) + Bm(5,0,a)*DBm(5,0);

      lK(0,a,b) = lK(0,a,b) + w*( T1 + afu*(BmDBm + Kvis_u(0,a,b)) + afv*Kvis_v(0,a,b) );

      // dM1/du2
      // Material stiffness: Bt*D*B
      BmDBm = Bm(0,0,a)*DBm(0,1) + Bm(1,0,a)*DBm(1,1) +
              Bm(2,0,a)*DBm(2,1) + Bm(3,0,a)*DBm(3,1) +
              Bm(4,0,a)*DBm(4,1) + Bm(5,0,a)*DBm(5,1);


      lK(1,a,b) = lK(1,a,b) + w*( afu*(BmDBm + Kvis_u(1,a,b)) + afv*(Kvis_v(1,a,b)) );

      // dM1/du3
      // Material stiffness: Bt*D*B
      BmDBm = Bm(0,0,a)*DBm(0,2) + Bm(1,0,a)*DBm(1,2) +
              Bm(2,0,a)*DBm(2,2) + Bm(3,0,a)*DBm(3,2) +
              Bm(4,0,a)*DBm(4,2) + Bm(5,0,a)*DBm(5,2);

      lK(2,a,b) = lK(2,a,b) + w*( afu*(BmDBm + Kvis_u(2,a,b)) + afv*Kvis_v(2,a,b) );

      // dM2/du1
      // Material stiffness: Bt*D*B
      BmDBm = Bm(0,1,a)*DBm(0,0) + Bm(1,1,a)*DBm(1,0) +
              Bm(2,1,a)*DBm(2,0) + Bm(3,1,a)*DBm(3,0) +
              Bm(4,1,a)*DBm(4,0) + Bm(5,1,a)*DBm(5,0);

      lK(dof+0,a,b) = lK(dof+0,a,b) + w*( afu*(BmDBm + Kvis_u(3,a,b)) + afv*Kvis_v(3,a,b) );

      // dM2/du2
      // Material stiffness: Bt*D*B
      BmDBm = Bm(0,1,a)*DBm(0,1) + Bm(1,1,a)*DBm(1,1) +
              Bm(2,1,a)*DBm(2,1) + Bm(3,1,a)*DBm(3,1) +
              Bm(4,1,a)*DBm(4,1) + Bm(5,1,a)*DBm(5,1);

      lK(dof+1,a,b) = lK(dof+1,a,b) + w*(T1 + afu*(BmDBm + Kvis_u(4,a,b)) + afv*Kvis_v(4,a,b) );

      // dM2/du3
      // Material stiffness: Bt*D*B
      BmDBm = Bm(0,1,a)*DBm(0,2) + Bm(1,1,a)*DBm(1,2) +
              Bm(2,1,a)*DBm(2,2) + Bm(3,1,a)*DBm(3,2) +
              Bm(4,1,a)*DBm(4,2) + Bm(5,1,a)*DBm(5,2);

      lK(dof+2,a,b) = lK(dof+2,a,b) + w*( afu*(BmDBm + Kvis_u(5,a,b)) + afv*Kvis_v(5,a,b) );

      // dM3/du1
      // Material stiffness: Bt*D*B
      BmDBm = Bm(0,2,a)*DBm(0,0) + Bm(1,2,a)*DBm(1,0) +
              Bm(2,2,a)*DBm(2,0) + Bm(3,2,a)*DBm(3,0) +
              Bm(4,2,a)*DBm(4,0) + Bm(5,2,a)*DBm(5,0);

      lK(2*dof+0,a,b) = lK(2*dof+0,a,b) + w*( afu*(BmDBm + Kvis_u(6,a,b)) + afv*Kvis_v(6,a,b) );

      // dM3/du2
      // Material stiffness: Bt*D*B
      BmDBm = Bm(0,2,a)*DBm(0,1) + Bm(1,2,a)*DBm(1,1) +
              Bm(2,2,a)*DBm(2,1) + Bm(3,2,a)*DBm(3,1) +
              Bm(4,2,a)*DBm(4,1) + Bm(5,2,a)*DBm(5,1);

     lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + w*( afu*(BmDBm + Kvis_u(7,a,b)) + afv*Kvis_v(7,a,b) );

      // dM3/du3
      // Material stiffness: Bt*D*B
      BmDBm = Bm(0,2,a)*DBm(0,2) + Bm(1,2,a)*DBm(1,2) +
              Bm(2,2,a)*DBm(2,2) + Bm(3,2,a)*DBm(3,2) +
              Bm(4,2,a)*DBm(4,2) + Bm(5,2,a)*DBm(5,2);

      lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + w*( T1 + afu*(BmDBm + Kvis_u(8,a,b)) + afv*Kvis_v(8,a,b) );
    }
  }
}

};

