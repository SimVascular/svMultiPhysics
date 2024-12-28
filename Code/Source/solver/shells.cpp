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

// This routines is for solving nonlinear shell mechanics problem
// using linear triangle finite elements and IGA.

#include "shells.h"

#include "all_fun.h"
#include "consts.h"
#include "lhsa.h"
#include "mat_fun.h"
#include "mat_models.h"
#include "nn.h"
#include "utils.h"

namespace shells {

/// @brief This routines is for solving nonlinear shell mechanics problem
/// using finite elements and IGA.
/// 
/// Reproduces Fortran CONSTRUCT_SHELL
//
void construct_shell(ComMod& com_mod, const mshType& lM, const Array<double>& Ag,
    const Array<double>& Yg, const Array<double>& Dg)
{
  using namespace consts;

  #define n_debug_construct_shell
  #ifdef debug_construct_shell
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "lM.nFn: " << lM.nFn;
  dmsg << "lM.nFs: " << lM.nFs;
  dmsg << "lM.eNoN: " << lM.eNoN;
  #endif

  const int nsd  = com_mod.nsd;
  const int tDof = com_mod.tDof;
  const int dof = com_mod.dof;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  auto cDmn = com_mod.cDmn;

  int eNoN = lM.eNoN;

  if (lM.eType == ElementType::TRI3) { 
    eNoN = 2 * eNoN;
  }

  int nFn = lM.nFn;
  if (nFn == 0) {
    nFn = 1;
  }

  // Initialize tensor operations
  mat_fun::ten_init(2);

  // SHELLS: dof = nsd
  //
  Vector<int> ptr(eNoN);
  Array<double> xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN),
                bfl(nsd,eNoN), fN(3,nFn), lR(dof,eNoN);
  Array3<double> lK(dof*dof,eNoN,eNoN); 

  // Loop over all elements of mesh
  //
  for (int e = 0; e < lM.nEl; e++) {
    // Update domain and proceed if domain phys and eqn phys match
    cDmn = all_fun::domain(com_mod, lM, cEq, e);
    auto cPhys = eq.dmn[cDmn].phys;
    if (cPhys != EquationType::phys_shell) {
      continue;
    }

#if 0
    dmsg << " " << " ";
    dmsg << " " << " ";
    dmsg << "-------------------------------------" << " ";
    dmsg << "--------------- e: " << e+1;
    dmsg << "-------------------------------------" << " ";
#endif

    //  Create local copies
    xl  = 0.0;
    al  = 0.0;
    yl  = 0.0;
    dl  = 0.0;
    bfl = 0.0;

    for (int a = 0; a < eNoN; a++) {
      //dmsg << "----- a: " << a;
      int Ac = -1;

      if (a < lM.eNoN) {
        Ac = lM.IEN(a,e);
        ptr(a) = Ac;
      } else {
        int b = a - lM.eNoN;
        //dmsg << "b: " << b;
        Ac = lM.eIEN(b,e);
        ptr(a) = Ac;
        //dmsg << "Ac: " << Ac;
        if (Ac == -1) {
          continue;
        }
      } 

      //dmsg << "Ac: " << Ac;

      for (int i = 0; i < nsd; i++) {
        xl(i,a) = com_mod.x(i,Ac);
        bfl(i,a) = com_mod.Bf(i,Ac);
      }

      for (int i = 0; i < tDof; i++) {
        al(i,a) = Ag(i,Ac);
        dl(i,a) = Dg(i,Ac);
        yl(i,a) = Yg(i,Ac);
      }
    }

    //dmsg << "lM.fN.size(): " << lM.fN.size();
    if (lM.fN.size() != 0) {
      for (int iFn = 0; iFn < nFn; iFn++) {
        for (int i = 0; i < 3; i++) {
          fN(i,iFn) = lM.fN(i+nsd*iFn,e);
        }
      }
    } else {
      fN = 0.0;
    }

    #ifdef debug_construct_shell
    dmsg << "-------------------" << "-------------------";
    dmsg << "ptr: " << ptr;
    dmsg << "-------------------" << "-------------------";
    #endif

    //  Constant strain triangles, no numerical integration
    //
    if (lM.eType == ElementType::TRI3) {
      shell_cst(com_mod, lM, e, eNoN, nFn, fN, al, yl, dl, xl, bfl, ptr);

    } else {
      lR = 0.0;
      lK = 0.0;

      // Update shape functions for NURBS elements
      //if (lM.eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

      // Gauss integration
      for (int g = 0; g < lM.nG; g++) {
        shell_3d(com_mod, lM, g, eNoN, nFn, fN, al, yl, dl, xl, bfl, lR, lK);
      }
    }

    eq.linear_algebra->assemble(com_mod, eNoN, ptr, lK, lR);

  } // e: loop

}

//----------
// shell_3d
//----------
// Construct shell mechanics for higher order elements/NURBS
//
void shell_3d(ComMod& com_mod, const mshType& lM, const int g, const int eNoN, 
    const int nFn, const Array<double>& fN,
    const Array<double>& al, const Array<double>& yl, const Array<double>& dl, const Array<double>& xl,
    const Array<double>& bfl, Array<double>& lR, Array3<double>& lK)
{
  std::cout << "========== shell_3d ==========" << std::endl;
  std::cout << "[shell_3d] g: " << g << std::endl;

  using namespace consts;
  using namespace mat_fun;

  const int nsd = com_mod.nsd;
  const int dof = com_mod.dof;
  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  auto cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  // Define parameters
  double rho = eq.dmn[cDmn].prop.at(PhysicalProperyType::solid_density);
  double dmp = dmn.prop.at(PhysicalProperyType::damping);
  double ht = eq.dmn[cDmn].prop.at(PhysicalProperyType::shell_thickness);
  Vector<double> fb({dmn.prop.at(PhysicalProperyType::f_x), dmn.prop.at(PhysicalProperyType::f_y), 
      dmn.prop.at(PhysicalProperyType::f_z)});
  double amd = eq.am * rho  +  eq.af * eq.gam * dt * dmp;
  double afl = eq.af * eq.beta * dt * dt;

  int i = eq.s;
  int j = i + 1;
  int k = j + 1;

  // Get the reference configuration
  auto x0 = xl;

  // Get the current configuration
  //
  Array<double> xc(3,eNoN);

  for (int a = 0; a < eNoN; a++) {
    xc(0,a) = x0(0,a) + dl(i,a);
    xc(1,a) = x0(1,a) + dl(j,a);
    xc(2,a) = x0(2,a) + dl(k,a);
  }

  // Define shape functions and their derivatives at Gauss point
  //
  Vector<double> N; 
  Array<double> Nx, Nxx;

  /* [TODO] Nurbs are not supported.
  if (lM.eType == ElementType::eType_NRB) {
    N = lM.N.rcol(g);
    Nx = lM.Nx.rslice(g);
    Nxx = lM.Nxx.rslice(g);
  } else {
    N = lM.fs(0).N.rcol(g);
    Nx = lM.fs(0).Nx.rslice(g);
    Nxx = lM.fs(0).Nxx.rslice(g);
  }
  */
  N = lM.fs[0].N.rcol(g);
  Nx = lM.fs[0].Nx.rslice(g);
  Nxx = lM.fs[0].Nxx.rslice(g);

//=====================================================================
//    TODO: Might have to call GNNxx for Jacobian transformation. Check
//    formulation again.
//
//     CALL GNNxx(2, eNoN, 2, lM.fs(0).Nx(:,:,g), lM.fs(0).Nxx(:,:,g),
//    xl, Nx, Nxx)
//
//=====================================================================

  // Compute preliminaries on the reference configuration
  // Covariant and contravariant bases (reference config)
  //
  Array<double> aCov0(3,2), aCnv0(3,2);
  Vector<double> nV0(3);
  nn::gnns(nsd, lM.eNoN, Nx, x0, nV0, aCov0, aCnv0);
  double Jac0 = sqrt(utils::norm(nV0));
  nV0 = nV0 / Jac0;

  // Second derivatives for computing curvature coeffs. (ref. config)
  //
  Array3<double> r0_xx(2,2,3);

  for (int a = 0; a < eNoN; a++) {
    for (int i = 0; i < 3; i++) {
      r0_xx(0,0,i) = r0_xx(0,0,i) + Nxx(0,a)*x0(i,a);
      r0_xx(1,1,i) = r0_xx(1,1,i) + Nxx(1,a)*x0(i,a);
      r0_xx(0,1,i) = r0_xx(0,1,i) + Nxx(2,a)*x0(i,a);
    }
  }

  for (int i = 0; i < 3; i++) {
    r0_xx(1,0,i) = r0_xx(0,1,i);
  }

  //  Compute metric tensor and curvature coefficients (ref. config)
  //
  double aa_0[2][2]{}, bb_0[2][2]{};

  for (int l = 0; l < nsd; l++) {
    aa_0[0][0] = aa_0[0][0] + aCov0(l,0)*aCov0(l,0);
    aa_0[0][1] = aa_0[0][1] + aCov0(l,0)*aCov0(l,1);
    aa_0[1][0] = aa_0[1][0] + aCov0(l,1)*aCov0(l,0);
    aa_0[1][1] = aa_0[1][1] + aCov0(l,1)*aCov0(l,1);

    bb_0[0][0] = bb_0[0][0] + r0_xx(0,0,l)*nV0(l);
    bb_0[0][1] = bb_0[0][1] + r0_xx(0,1,l)*nV0(l);
    bb_0[1][0] = bb_0[1][0] + r0_xx(1,0,l)*nV0(l);
    bb_0[1][1] = bb_0[1][1] + r0_xx(1,1,l)*nV0(l);
  }

  // Compute fiber orientation in curvature coordinates
  //
  Array<double> fNa0(2,nFn);

  for (int iFn = 0; iFn < nFn; iFn++) {
    for (int l = 0; l < 3; l++) { 
      fNa0(0,iFn) = fNa0(0,iFn) + fN(l,iFn)*aCnv0(l,0);
      fNa0(1,iFn) = fNa0(1,iFn) + fN(l,iFn)*aCnv0(l,1);
    }
  }

  // Now compute preliminaries on the current configuration
  // Covariant and contravariant bases (current/spatial config)
  //
  Array<double> aCov(3,2), aCnv(3,2);
  Vector<double> nV(3);
  nn::gnns(nsd, eNoN, Nx, xc, nV, aCov, aCnv);
  double Jac = sqrt(utils::norm(nV));
  nV = nV / Jac;

  // Second derivatives for computing curvature coeffs. (cur. config)
  Array3<double> r_xx(2,2,3);

  for (int a = 0; a < eNoN; a++) {
     r_xx(0,0,i) = r_xx(0,0,i) + Nxx(0,a)*xc(i,a);
     r_xx(1,1,i) = r_xx(1,1,i) + Nxx(1,a)*xc(i,a);
     r_xx(0,1,i) = r_xx(0,1,i) + Nxx(2,a)*xc(i,a);
  }

  for (int i = 0; i < 3; i++) {
    r_xx(1,0,i) = r_xx(0,1,i);
  }

  // Compute metric tensor and curvature coefficients (cur. config)
  //
  // [TODO] are these dimensions correct? is nsd = 3?
  //
  double aa_x[2][2]{}, bb_x[2][2]{};

  for (int l = 0; l < nsd; l++) {
     aa_x[0][0] = aa_x[0][0] + aCov(l,0)*aCov(l,0);
     aa_x[0][1] = aa_x[0][1] + aCov(l,0)*aCov(l,1);
     aa_x[1][0] = aa_x[1][0] + aCov(l,1)*aCov(l,0);
     aa_x[1][1] = aa_x[1][1] + aCov(l,1)*aCov(l,1);

     bb_x[0][0] = bb_x[0][0] + r_xx(0,0,l)*nV(l);
     bb_x[0][1] = bb_x[0][1] + r_xx(0,1,l)*nV(l);
     bb_x[1][0] = bb_x[1][0] + r_xx(1,0,l)*nV(l);
     bb_x[1][1] = bb_x[1][1] + r_xx(1,1,l)*nV(l);
  }

  // Compute stress resultants by integrating 2nd Piola Kirchhoff
  // stress and elasticity tensors through the shell thickness. These
  // resultants are computed in Voigt notation.
  //
  Array3<double> Dm(3,3,3);
  Array<double> Sm(3,2);
  double lam3;

  shl_strs_res(com_mod, dmn, nFn, fNa0, aa_0, aa_x, bb_0, bb_x, lam3, Sm, Dm);

  // Variation in the membrane strain
  //
  Array3<double> Bm(3,3,eNoN);

  for (int a = 0; a < eNoN; a++) {
    Bm(0,0,a) = Nx(0,a)*aCov(0,0);
    Bm(0,1,a) = Nx(0,a)*aCov(1,0);
    Bm(0,2,a) = Nx(0,a)*aCov(2,0);

    Bm(1,0,a) = Nx(1,a)*aCov(0,1);
    Bm(1,1,a) = Nx(1,a)*aCov(1,1);
    Bm(1,2,a) = Nx(1,a)*aCov(2,1);

    Bm(2,0,a) = Nx(1,a)*aCov(0,0) + Nx(0,a)*aCov(0,1);
    Bm(2,1,a) = Nx(1,a)*aCov(1,0) + Nx(0,a)*aCov(1,1);
    Bm(2,2,a) = Nx(1,a)*aCov(2,0) + Nx(0,a)*aCov(2,1);
  }

  // Variation in the bending strain
  // dB = -(B1 + B2) du; B1 = N_xx * n;
  // B2 = (r_xx Nm M1 Nx - r_xx N M2 Nx)
  //
  //     Second derivatives of the position vector (current)
  //
  Array<double> Kc(3,3);

  for (int i = 0; i < 3; i++) {
    Kc(0,i) = r_xx(0,0,i);
    Kc(1,i) = r_xx(1,1,i);
    Kc(2,i) = r_xx(0,1,i)  + r_xx(1,0,i);
  }

  // N matrix
  auto Nm = mat_id(3) - mat_dyad_prod(nV, nV, 3);
  Nm = Nm / Jac;

  // M1, M2 matrices
  //
  Array<double> Mm(3,3);
  Array3<double> KNmMm(3,3,2);

  for (int l = 0; l < 2; l++) { 
    Mm = 0.0;
    Mm(0,1) = -aCov(2,l);
    Mm(0,2) =  aCov(1,l);
    Mm(1,2) = -aCov(0,l);

    // Skew-symmetric
    Mm(1,0) = -Mm(0,1);
    Mm(2,0) = -Mm(0,2);
    Mm(2,1) = -Mm(1,2);

    KNmMm.set_slice(l, mat_mul(Kc, mat_mul(Nm, Mm)));
  }

  // Define variation in bending strain tensor (Bb), Voigt notation
  //
  Array3<double> Bb(3,3,eNoN);

  for (int a = 0; a < eNoN; a++) {
    for (int i = 0; i < 3; i++) {
      Bb(0,i,a) = -Nxx(0,a)*nV(i);
      Bb(1,i,a) = -Nxx(1,a)*nV(i);
      Bb(2,i,a) = -Nxx(2,a)*nV(i)*2.0;
    }
  }

  for (int a = 0; a < eNoN; a++) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        Bb(i,j,a) = Bb(i,j,a) + Nx(0,a)*KNmMm(i,j,1) - Nx(1,a)*KNmMm(i,j,0);
      }
    }
  }

  //  Contribution to tangent matrices: Dm * Bm, Dm*Bb
  //
  Array3<double> D0Bm(3,3,eNoN), D1Bm(3,3,eNoN), D1Bb(3,3,eNoN), D2Bb(3,3,eNoN);

  for (int a = 0; a < eNoN; a++) {
    D0Bm.set_slice(a, mat_mul(Dm.rslice(0), Bm.rslice(a)));
    D1Bm.set_slice(a, mat_mul(Dm.rslice(1), Bm.rslice(a)));
    D1Bb.set_slice(a, mat_mul(Dm.rslice(1), Bb.rslice(a)));
    D2Bb.set_slice(a, mat_mul(Dm.rslice(2), Bb.rslice(a)));
  }

  // Acceleration and mass damping at the integration point
  //
  auto ud = -fb;

  for (int a = 0; a < eNoN; a++) {
    ud(0) = ud(0) + N(a)*(rho*(al(i,a)-bfl(0,a)) + dmp*yl(i,a));
    ud(1) = ud(1) + N(a)*(rho*(al(j,a)-bfl(1,a)) + dmp*yl(j,a));
    ud(2) = ud(2) + N(a)*(rho*(al(k,a)-bfl(2,a)) + dmp*yl(k,a));
  }

  // Local residual
  //
  auto w = lM.w(g)*Jac0;
  auto wh = w*ht;

  for (int a = 0; a < eNoN; a++) {
     double BmS = Bm(0,0,a)*Sm(0,0) + Bm(1,0,a)*Sm(1,0) + Bm(2,0,a)*Sm(2,0);
     double BbS = Bb(0,0,a)*Sm(0,1) + Bb(1,0,a)*Sm(1,1) + Bb(2,0,a)*Sm(2,1);
     lR(0,a) = lR(0,a) + wh*N(a)*ud(0) + w*(BmS + BbS);

     BmS = Bm(0,1,a)*Sm(0,0) + Bm(1,1,a)*Sm(1,0) + Bm(2,1,a)*Sm(2,0);
     BbS = Bb(0,1,a)*Sm(0,1) + Bb(1,1,a)*Sm(1,1) + Bb(2,1,a)*Sm(2,1);
     lR(1,a) = lR(1,a) + wh*N(a)*ud(1) + w*(BmS + BbS);

     BmS = Bm(0,2,a)*Sm(0,0) + Bm(1,2,a)*Sm(1,0) + Bm(2,2,a)*Sm(2,0);
     BbS = Bb(0,2,a)*Sm(0,1) + Bb(1,2,a)*Sm(1,1) + Bb(2,2,a)*Sm(2,1);
     lR(2,a) = lR(2,a) + wh*N(a)*ud(2) + w*(BmS + BbS);
  }

  // Local stiffness
  //
  amd = wh*amd;
  afl = w*afl;

  for (int b = 0; b < eNoN; b++) {
    for (int a = 0; a < eNoN; a++) {
      // Contribution from inertia and geometric stiffness
      double NxSNx = Nx(0,a)*Nx(0,b)*Sm(0,0) + Nx(1,a)*Nx(1,b)*Sm(1,0) + Nx(0,a)*Nx(1,b)*Sm(2,0) + 
          Nx(1,a)*Nx(0,b)*Sm(2,0);
      double T1 = amd*N(a)*N(b) + afl*NxSNx;

      lK(0,a,b) = lK(0,a,b) + T1;
      lK(dof+1,a,b) = lK(dof+1,a,b) + T1;
      lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + T1;

      // Contribution from material stiffness
      double BmDBm = Bm(0,0,a)*D0Bm(0,0,b) + Bm(1,0,a)*D0Bm(1,0,b) + Bm(2,0,a)*D0Bm(2,0,b);
      double BmDBb = Bm(0,0,a)*D1Bb(0,0,b) + Bm(1,0,a)*D1Bb(1,0,b) + Bm(2,0,a)*D1Bb(2,0,b);
      double BbDBm = Bb(0,0,a)*D1Bm(0,0,b) + Bb(1,0,a)*D1Bm(1,0,b) + Bb(2,0,a)*D1Bm(2,0,b);
      double BbDBb = Bb(0,0,a)*D2Bb(0,0,b) + Bb(1,0,a)*D2Bb(1,0,b) + Bb(2,0,a)*D2Bb(2,0,b);
      lK(0,a,b) = lK(0,a,b) + afl*(BmDBm + BmDBb + BbDBm + BbDBb);

      BmDBm = Bm(0,0,a)*D0Bm(0,1,b) + Bm(1,0,a)*D0Bm(1,1,b) + Bm(2,0,a)*D0Bm(2,1,b);
      BmDBb = Bm(0,0,a)*D1Bb(0,1,b) + Bm(1,0,a)*D1Bb(1,1,b) + Bm(2,0,a)*D1Bb(2,1,b);
      BbDBm = Bb(0,0,a)*D1Bm(0,1,b) + Bb(1,0,a)*D1Bm(1,1,b) + Bb(2,0,a)*D1Bm(2,1,b);
      BbDBb = Bb(0,0,a)*D2Bb(0,1,b) + Bb(1,0,a)*D2Bb(1,1,b) + Bb(2,0,a)*D2Bb(2,1,b);
      lK(1,a,b) = lK(1,a,b) + afl*(BmDBm + BmDBb + BbDBm + BbDBb);

      BmDBm = Bm(0,0,a)*D0Bm(0,2,b) + Bm(1,0,a)*D0Bm(1,2,b) + Bm(2,0,a)*D0Bm(2,2,b);
      BmDBb = Bm(0,0,a)*D1Bb(0,2,b) + Bm(1,0,a)*D1Bb(1,2,b) + Bm(2,0,a)*D1Bb(2,2,b);
      BbDBm = Bb(0,0,a)*D1Bm(0,2,b) + Bb(1,0,a)*D1Bm(1,2,b) + Bb(2,0,a)*D1Bm(2,2,b);
      BbDBb = Bb(0,0,a)*D2Bb(0,2,b) + Bb(1,0,a)*D2Bb(1,2,b) + Bb(2,0,a)*D2Bb(2,2,b);
      lK(2,a,b) = lK(2,a,b) + afl*(BmDBm + BmDBb + BbDBm + BbDBb);

      BmDBm = Bm(0,1,a)*D0Bm(0,0,b) + Bm(1,1,a)*D0Bm(1,0,b) + Bm(2,1,a)*D0Bm(2,0,b);
      BmDBb = Bm(0,1,a)*D1Bb(0,0,b) + Bm(1,1,a)*D1Bb(1,0,b) + Bm(2,1,a)*D1Bb(2,0,b);
      BbDBm = Bb(0,1,a)*D1Bm(0,0,b) + Bb(1,1,a)*D1Bm(1,0,b) + Bb(2,1,a)*D1Bm(2,0,b);
      BbDBb = Bb(0,1,a)*D2Bb(0,0,b) + Bb(1,1,a)*D2Bb(1,0,b) + Bb(2,1,a)*D2Bb(2,0,b);
      lK(dof+0,a,b) = lK(dof+0,a,b) + afl*(BmDBm + BmDBb + BbDBm + BbDBb);

      BmDBm = Bm(0,1,a)*D0Bm(0,1,b) + Bm(1,1,a)*D0Bm(1,1,b) + Bm(2,1,a)*D0Bm(2,1,b);
      BbDBm = Bb(0,1,a)*D1Bm(0,1,b) + Bb(1,1,a)*D1Bm(1,1,b) + Bb(2,1,a)*D1Bm(2,1,b);
      BbDBb = Bb(0,1,a)*D2Bb(0,1,b) + Bb(1,1,a)*D2Bb(1,1,b) + Bb(2,1,a)*D2Bb(2,1,b);
      lK(dof+1,a,b) = lK(dof+1,a,b) + afl*(BmDBm + BmDBb + BbDBm + BbDBb);

      BmDBm = Bm(0,1,a)*D0Bm(0,2,b) + Bm(1,1,a)*D0Bm(1,2,b) + Bm(2,1,a)*D0Bm(2,2,b);
      BmDBb = Bm(0,1,a)*D1Bb(0,2,b) + Bm(1,1,a)*D1Bb(1,2,b) + Bm(2,1,a)*D1Bb(2,2,b);
      BbDBm = Bb(0,1,a)*D1Bm(0,2,b) + Bb(1,1,a)*D1Bm(1,2,b) + Bb(2,1,a)*D1Bm(2,2,b);
      BbDBb = Bb(0,1,a)*D2Bb(0,2,b) + Bb(1,1,a)*D2Bb(1,2,b) + Bb(2,1,a)*D2Bb(2,2,b);
      lK(dof+2,a,b) = lK(dof+2,a,b) + afl*(BmDBm + BmDBb + BbDBm + BbDBb);

      BmDBm = Bm(0,2,a)*D0Bm(0,0,b) + Bm(1,2,a)*D0Bm(1,0,b) + Bm(2,2,a)*D0Bm(2,0,b);
      BmDBb = Bm(0,2,a)*D1Bb(0,0,b) + Bm(1,2,a)*D1Bb(1,0,b) + Bm(2,2,a)*D1Bb(2,0,b);
      BbDBm = Bb(0,2,a)*D1Bm(0,0,b) + Bb(1,2,a)*D1Bm(1,0,b) + Bb(2,2,a)*D1Bm(2,0,b);
      BbDBb = Bb(0,2,a)*D2Bb(0,0,b) + Bb(1,2,a)*D2Bb(1,0,b) + Bb(2,2,a)*D2Bb(2,0,b);
      lK(2*dof+0,a,b) = lK(2*dof+0,a,b) + afl*(BmDBm + BmDBb + BbDBm + BbDBb);

      BmDBm = Bm(0,2,a)*D0Bm(0,1,b) + Bm(1,2,a)*D0Bm(1,1,b) + Bm(2,2,a)*D0Bm(2,1,b);
      BmDBb = Bm(0,2,a)*D1Bb(0,1,b) + Bm(1,2,a)*D1Bb(1,1,b) + Bm(2,2,a)*D1Bb(2,1,b);
      BbDBm = Bb(0,2,a)*D1Bm(0,1,b) + Bb(1,2,a)*D1Bm(1,1,b) + Bb(2,2,a)*D1Bm(2,1,b);
      BbDBb = Bb(0,2,a)*D2Bb(0,1,b) + Bb(1,2,a)*D2Bb(1,1,b) + Bb(2,2,a)*D2Bb(2,1,b);
      lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + afl*(BmDBm + BmDBb + BbDBm + BbDBb);

      BmDBm = Bm(0,2,a)*D0Bm(0,2,b) + Bm(1,2,a)*D0Bm(1,2,b) + Bm(2,2,a)*D0Bm(2,2,b);
      BmDBb = Bm(0,2,a)*D1Bb(0,2,b) + Bm(1,2,a)*D1Bb(1,2,b) + Bm(2,2,a)*D1Bb(2,2,b);
      BbDBm = Bb(0,2,a)*D1Bm(0,2,b) + Bb(1,2,a)*D1Bm(1,2,b) + Bb(2,2,a)*D1Bm(2,2,b);
      BbDBb = Bb(0,2,a)*D2Bb(0,2,b) + Bb(1,2,a)*D2Bb(1,2,b) + Bb(2,2,a)*D2Bb(2,2,b);
      lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + afl*(BmDBm + BmDBb + BbDBm + BbDBb);
    }
  }
}

//----------------
// shell_bend_cst
//----------------
// This subroutine computes bending strain, Eb, and its variational
// derivative, Bb, for CST elements
//
// Reproduces Fortran SHELLBENDCST.
//
void shell_bend_cst(ComMod& com_mod, const mshType& lM, const int e, const Vector<int>& ptr, 
    Array<double>& x0, Array<double>& xc, double bb_0[2][2], double bb_x[2][2], 
    Array3<double>& Bb, const bool vflag)
{
  using namespace consts;
  using namespace mat_fun;
  using namespace utils;

  #define n_debug_shell_bend_cst
  #ifdef debug_shell_bend_cst 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "e: " << e;
  #endif

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  auto cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];

  // Define parameters
  double rho = eq.dmn[cDmn].prop.at(PhysicalProperyType::solid_density);
  double dmp = dmn.prop.at(PhysicalProperyType::damping);

  int nsd = com_mod.nsd;
  int eNoN = 2 * lM.eNoN;

  // Boundary element check
  bool bFlag = false;

  for (int j = lM.eNoN; j < eNoN; j++) {
    if (ptr(j) == -1) {
      bFlag = true;
      break;
    }
  }

  //  Edge vectors of the main element (reference config)
  //
  Array<double> a0(3,6); 

  for (int i = 0; i < 3; i++) {
    a0(i,0) = x0(i,2) - x0(i,1);
    a0(i,1) = x0(i,0) - x0(i,2);
    a0(i,2) = x0(i,1) - x0(i,0);
  }

  // Edge vectors of the main element (current config)
  //
  Array<double> a(3,6); 

  for (int i = 0; i < 3; i++) {
    a(i,0) = xc(i,2) - xc(i,1);
    a(i,1) = xc(i,0) - xc(i,2);
    a(i,2) = xc(i,1) - xc(i,0);
  }

  // Covariant and contravariant bases in reference config
  Array<double> tmpA(3,3);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      tmpA(i,j) = x0(i,j);
    }
  }

  Array<double> aCov0(3,2), aCnv0(3,2);
  Vector<double> nV0(3);
  nn::gnns(nsd, lM.eNoN, lM.Nx.rslice(0), tmpA, nV0, aCov0, aCnv0);
  double Jac0 = sqrt(utils::norm(nV0));
  nV0 = nV0 / Jac0;

  // Covariant and contravariant bases in current config
  //
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      tmpA(i,j) = xc(i,j);
    }
  }

  Array<double> aCov(3,2), aCnv(3,2);
  Vector<double> nV(3);
  nn::gnns(nsd, lM.eNoN, lM.Nx.rslice(0), tmpA, nV, aCov, aCnv);
  double Jac = sqrt(norm(nV));
  nV = nV / Jac;

  // Update the position vector of the `artificial' or `ghost' nodes
  // depending on the boundary condition.
  //
  if (bFlag) {
    for (int j = lM.eNoN; j < eNoN; j++) {
      if (ptr(j) != -1) {
        continue;
      }
      int i = j - lM.eNoN;
      int p = i - 1;

      if (i == 0) {
        p = 2;
      }

      // Reference config
      // eI = eI0 = aI0/|aI0| (reference config)
      //
      auto aIi = 1.0 / sqrt(norm(a0.col(i)));
      auto eI = a0.col(i) * aIi;

      // nI = nI0 = eI0 x n0 (reference config)
      //
      Vector<double> nI(3);
      nI(0) = eI(1)*nV0(2) - eI(2)*nV0(1);
      nI(1) = eI(2)*nV0(0) - eI(0)*nV0(2);
      nI(2) = eI(0)*nV0(1) - eI(1)*nV0(0);

      // xJ = xI + 2(nI \ctimes nI)aP
      //
      auto nInI = mat_dyad_prod(nI, nI, 3);
      x0(0,j) = 2.0 * (nInI(0,0)*a0(0,p) + nInI(0,1)*a0(1,p) + nInI(0,2)*a0(2,p)) + x0(0,i);
      x0(1,j) = 2.0 * (nInI(1,0)*a0(0,p) + nInI(1,1)*a0(1,p) + nInI(1,2)*a0(2,p)) + x0(1,i);
      x0(2,j) = 2.0 * (nInI(2,0)*a0(0,p) + nInI(2,1)*a0(1,p) + nInI(2,2)*a0(2,p)) + x0(2,i);

      // Current config
      // eI = aI/|aI| (current config)
      //
      aIi = 1.0 / sqrt(utils::norm(a.col(i)));
      eI = a.col(i)*aIi;

      // nI = eI x n (currnt config)
      //
      nI(0) = eI(1)*nV(2) - eI(2)*nV(1);
      nI(1) = eI(2)*nV(0) - eI(0)*nV(2);
      nI(2) = eI(0)*nV(1) - eI(1)*nV(0);

      // xJ = xI + 2(nI \ctimes nI)aP
      //
      nInI = mat_dyad_prod(nI, nI, 3);
      xc(0,j) = 2.0*(nInI(0,0)*a(0,p) + nInI(0,1)*a(1,p) + nInI(0,2)*a(2,p)) + xc(0,i);
      xc(1,j) = 2.0*(nInI(1,0)*a(0,p) + nInI(1,1)*a(1,p) + nInI(1,2)*a(2,p)) + xc(1,i);
      xc(2,j) = 2.0*(nInI(2,0)*a(0,p) + nInI(2,1)*a(1,p) + nInI(2,2)*a(2,p)) + xc(2,i);

      if (utils::btest(lM.sbc(i,e),enum_int(BoundaryConditionType::bType_fix))) {
        for (int i = 0; i < 3; i++) {
          xc(i,j) = x0(i,j);
        }
      }
    }
  }

  // Edge vector of surrounding nodes (reference config)
  //
  for (int i = 0; i < 3; i++) {
    a0(i,3) = x0(i,3) - x0(i,1);
    a0(i,4) = x0(i,4) - x0(i,2);
    a0(i,5) = x0(i,5) - x0(i,0);
  }

  // Edge vector of surrounding nodes (current config)
  //
  for (int i = 0; i < 3; i++) {
    a(i,3) = xc(i,3) - xc(i,1);
    a(i,4) = xc(i,4) - xc(i,2);
    a(i,5) = xc(i,5) - xc(i,0);
  }

  // a.gCnv (reference config)
  //
  Array<double> adg0(3,3);
  adg0(0,0) = norm(a0.col(3), aCnv0.col(0));    // xi_4
  adg0(0,1) = norm(a0.col(4), aCnv0.col(0));    // xi_5
  adg0(0,2) = norm(a0.col(5), aCnv0.col(0));    // xi_6

  adg0(1,0) = norm(a0.col(3), aCnv0.col(1));    // eta_4
  adg0(1,1) = norm(a0.col(4), aCnv0.col(1));    // eta_5
  adg0(1,2) = norm(a0.col(5), aCnv0.col(1));    // eta_6

  adg0(2,0) = norm(a0.col(3), nV0);        // z_4
  adg0(2,1) = norm(a0.col(4), nV0);        // z_5
  adg0(2,2) = norm(a0.col(5), nV0);        // z_6

  // a.gCnv (current config)
  //
  Array<double> adg(3,3);
  adg(0,0) = norm(a.col(3), aCnv.col(0));   // xi_4
  adg(0,1) = norm(a.col(4), aCnv.col(0));   // xi_5
  adg(0,2) = norm(a.col(5), aCnv.col(0));   // xi_6

  adg(1,0) = norm(a.col(3), aCnv.col(1));   // eta_4
  adg(1,1) = norm(a.col(4), aCnv.col(1));   // eta_5
  adg(1,2) = norm(a.col(5), aCnv.col(1));   // eta_6

  adg(2,0) = norm(a.col(3), nV);        // z_4
  adg(2,1) = norm(a.col(4), nV);        // z_5
  adg(2,2) = norm(a.col(5), nV);        // z_6

  // Xi matrix (reference config)
  //
  auto xi0 = adg0;
  xi0(0,2) = adg0(0,2) + 1.0;    // xi_6
  xi0(1,0) = adg0(1,0) + 1.0;    // eta_4

  //  Xi matrix (current config)
  //
  auto xi = adg;
  xi(0,2) = adg(0,2) + 1.0;     // xi_6
  xi(1,0) = adg(1,0) + 1.0;     // eta_4

  // Tmat and inverse (reference config)
  //
  Array<double> Tm0(3,3);

  for (int i = 0; i < 3; i++) {
    Tm0(i,0) = xi0(0,i)*(xi0(0,i) - 1.0);   // xi**2 - xi
    Tm0(i,1) = xi0(1,i)*(xi0(1,i) - 1.0);   // eta**2 - eta
    Tm0(i,2) = xi0(0,i)*xi0(1,i);           // xi * eta
  }

  Tm0 = mat_inv(Tm0, 3);

  // Tmat and inverse (current config)
  //
  Array<double> Tm(3,3);

  for (int i = 0; i < 3; i++) {
    Tm(i,0) = xi(0,i)*(xi(0,i) - 1.0);   // xi**2 - xi
    Tm(i,1) = xi(1,i)*(xi(1,i) - 1.0);   // eta**2 - eta
    Tm(i,2) = xi(0,i)*xi(1,i);           // xi * eta
  }

  Tm = mat_inv(Tm, 3);

  // v = Inv(T) * z (reference config)
  //
  Vector<double> v0(3);

  v0(0) = Tm0(0,0)*xi0(2,0) + Tm0(0,1)*xi0(2,1) + Tm0(0,2)*xi0(2,2);
  v0(1) = Tm0(1,0)*xi0(2,0) + Tm0(1,1)*xi0(2,1) + Tm0(1,2)*xi0(2,2);
  v0(2) = Tm0(2,0)*xi0(2,0) + Tm0(2,1)*xi0(2,1) + Tm0(2,2)*xi0(2,2);

  // Curvature coefficients (ref. config)
  bb_0[0][0] = 2.0 * v0(0);
  bb_0[1][1] = 2.0 * v0(1);
  bb_0[0][1] = v0(2);
  bb_0[1][0] = bb_0[0][1];

  // v = Inv(T) * z (current config)
  //
  double v[3];
  v[0] = Tm(0,0)*xi(2,0) + Tm(0,1)*xi(2,1) + Tm(0,2)*xi(2,2);
  v[1] = Tm(1,0)*xi(2,0) + Tm(1,1)*xi(2,1) + Tm(1,2)*xi(2,2);
  v[2] = Tm(2,0)*xi(2,0) + Tm(2,1)*xi(2,1) + Tm(2,2)*xi(2,2);

  // Curvature coefficients (current config)
  //
  bb_x[0][0] = 2.0*v[0];
  bb_x[1][1] = 2.0*v[1];
  bb_x[0][1] = v[2];
  bb_x[1][0] = bb_x[0][1];

  if (!vflag) {
    Bb = 0.0;
    return;
  }

  // Now compute variation in bending strain
  //
  // B1 bar
  //
  Array<double> B1b(3,6);

  for (int i = 0; i < 3; i++) {
    B1b(i,0) = -Tm(i,0) * ((2.0*xi(0,0)-1.0)*v[0] + xi(1,0)*v[2]);
    B1b(i,1) = -Tm(i,0) * ((2.0*xi(1,0)-1.0)*v[1] + xi(0,0)*v[2]);

    B1b(i,2) = -Tm(i,1) * ((2.0*xi(0,1)-1.0)*v[0] + xi(1,1)*v[2]);
    B1b(i,3) = -Tm(i,1) * ((2.0*xi(1,1)-1.0)*v[1] + xi(0,1)*v[2]);

    B1b(i,4) = -Tm(i,2) * ((2.0*xi(0,2)-1.0)*v[0] + xi(1,2)*v[2]);
    B1b(i,5) = -Tm(i,2) * ((2.0*xi(1,2)-1.0)*v[1] + xi(0,2)*v[2]);
  }

  //  H1
  //
  Array<double> H1(6,18);

  for (int i = 0; i < 3; i++) {
    H1(0, i) = aCnv(i,0)*adg(1,0);
    H1(1, i) = aCnv(i,1)*adg(1,0);
    H1(2, i) = aCnv(i,0)*adg(1,1);
    H1(3, i) = aCnv(i,1)*adg(1,1);
    H1(4, i) = aCnv(i,0)*adg(1,2);
    H1(5, i) = aCnv(i,1)*adg(1,2);

    H1(0, i+3 ) = -aCnv(i,0)*adg(0,0);
    H1(1, i+3 ) = -aCnv(i,1)*adg(0,0);
    H1(2, i+3 ) = -aCnv(i,0)*adg(0,1);
    H1(3, i+3 ) = -aCnv(i,1)*adg(0,1);
    H1(4, i+3 ) = -aCnv(i,0)*adg(0,2);
    H1(5, i+3 ) = -aCnv(i,1)*adg(0,2);

    H1(0,i+9) =  aCnv(i,0);
    H1(1,i+9) =  aCnv(i,1);

    H1(2,i+12) =  aCnv(i,0);
    H1(3,i+12) =  aCnv(i,1);

    H1(4,i+15) =  aCnv(i,0);
    H1(5,i+15) =  aCnv(i,1);
  }

  // H2
  Array<double> H2(18,18);
  H2( 0, 3) = -1.0;
  H2( 1, 4) = -1.0;
  H2( 2, 5) = -1.0;
  H2( 3, 6) = -1.0;
  H2( 4, 7) = -1.0;
  H2( 5, 8) = -1.0;
  H2( 6, 0) = -1.0;
  H2( 7, 1) = -1.0;
  H2( 8, 2) = -1.0;

  H2( 0, 6) =  1.0;
  H2( 1, 7) =  1.0;
  H2( 2, 8) =  1.0;
  H2( 3, 0) =  1.0;
  H2( 4, 1) =  1.0;
  H2( 5, 2) =  1.0;
  H2( 6, 3) =  1.0;
  H2( 7, 4) =  1.0;
  H2( 8, 5) =  1.0;

  H2(9, 3) = -1.0;
  H2(10, 4) = -1.0;
  H2(11, 5) = -1.0;
  H2(12, 6) = -1.0;
  H2(13, 7) = -1.0;
  H2(14, 8) = -1.0;
  H2(15, 0) = -1.0;
  H2(16, 1) = -1.0;
  H2(17, 2) = -1.0;

  for (int i = 9; i < 18; i++) {
    H2(i,i) = 1.0;
  }

  // N matrix
  auto Nm = mat_id(3) - mat_dyad_prod(nV, nV, 3);

  // M1, M2 matrices
  //
  Array3<double> Mm(3,3,2);

  for (int i = 0; i < 2; i++) {
    Mm(0,1,i) = -aCov(2,i);
    Mm(0,2,i) =  aCov(1,i);
    Mm(1,2,i) = -aCov(0,i);

    // Skew-symmetric
    Mm(1,0,i) = -Mm(0,1,i);
    Mm(2,0,i) = -Mm(0,2,i);
    Mm(2,1,i) = -Mm(1,2,i);
  }

  // H3 matrix
  //
  Array<double> H3(3,18);
  tmpA = mat_mul(Nm, Mm.rslice(0));
  tmpA = -tmpA / Jac;

  H3(0,0) = a(0,3)*tmpA(0,0) + a(1,3)*tmpA(1,0) + a(2,3)*tmpA(2,0);
  H3(0,1) = a(0,3)*tmpA(0,1) + a(1,3)*tmpA(1,1) + a(2,3)*tmpA(2,1);
  H3(0,2) = a(0,3)*tmpA(0,2) + a(1,3)*tmpA(1,2) + a(2,3)*tmpA(2,2);

  H3(1,0) = a(0,4)*tmpA(0,0) + a(1,4)*tmpA(1,0) + a(2,4)*tmpA(2,0);
  H3(1,1) = a(0,4)*tmpA(0,1) + a(1,4)*tmpA(1,1) + a(2,4)*tmpA(2,1);
  H3(1,2) = a(0,4)*tmpA(0,2) + a(1,4)*tmpA(1,2) + a(2,4)*tmpA(2,2);

  H3(2,0) = a(0,5)*tmpA(0,0) + a(1,5)*tmpA(1,0) + a(2,5)*tmpA(2,0);
  H3(2,1) = a(0,5)*tmpA(0,1) + a(1,5)*tmpA(1,1) + a(2,5)*tmpA(2,1);
  H3(2,2) = a(0,5)*tmpA(0,2) + a(1,5)*tmpA(1,2) + a(2,5)*tmpA(2,2);

  tmpA = mat_mul(Nm, Mm.rslice(1));
  tmpA = -tmpA / Jac;

  H3(0,3) = a(0,3)*tmpA(0,0) + a(1,3)*tmpA(1,0) + a(2,3)*tmpA(2,0);
  H3(0,4) = a(0,3)*tmpA(0,1) + a(1,3)*tmpA(1,1) + a(2,3)*tmpA(2,1);
  H3(0,5) = a(0,3)*tmpA(0,2) + a(1,3)*tmpA(1,2) + a(2,3)*tmpA(2,2);

  H3(1,3) = a(0,4)*tmpA(0,0) + a(1,4)*tmpA(1,0) + a(2,4)*tmpA(2,0);
  H3(1,4) = a(0,4)*tmpA(0,1) + a(1,4)*tmpA(1,1) + a(2,4)*tmpA(2,1);
  H3(1,5) = a(0,4)*tmpA(0,2) + a(1,4)*tmpA(1,2) + a(2,4)*tmpA(2,2);

  H3(2,3) = a(0,5)*tmpA(0,0) + a(1,5)*tmpA(1,0) + a(2,5)*tmpA(2,0);
  H3(2,4) = a(0,5)*tmpA(0,1) + a(1,5)*tmpA(1,1) + a(2,5)*tmpA(2,1);
  H3(2,5) = a(0,5)*tmpA(0,2) + a(1,5)*tmpA(1,2) + a(2,5)*tmpA(2,2);

  for (int i = 0; i < 3; i++) {
    H3(0,i+9) = nV(i);
    H3(1,i+12) = nV(i);
    H3(2,i+15) = nV(i);
  }

  // Variation in bending strain (Bb = -2*(B1b*H1*H2 + Tinv*H3*H2))
  auto H1H2 = mat_mul(H1, H2);
  auto Bb1 = mat_mul(B1b, H1H2);
  auto H3H2 = mat_mul(H3, H2);
  Bb1 = -2.0 * (Bb1 + mat_mul(Tm, H3H2));
  Bb.set_values(Bb1);
  //Bb   = RESHAPE(Bb1, SHAPE(Bb))

  //  Update Bb for boundary elements
  if (bFlag) {
    std::vector<bool> lFix = {false, false, false};
    auto Im = mat_id(3);

    for (int j = lM.eNoN; j < eNoN; j++) {
      if (ptr(j) != -1) {
        continue; 
      }

      int i = j - lM.eNoN;
      int p = i - 1;
      int f = i + 1;

      if (i == 0) {
        p = 2;
      }

      if (i == 2) {
        f = 0;
      }

      double aIi, cI;
      Vector<double> eI(3), nI(3);
      Array<double> nInI(3,3), eIeI(3,3), eIaP(3,3); 

      if (utils::btest(lM.sbc(i,e),enum_int(BoundaryConditionType::bType_fix))) {
        // eI = eI0 = aI0/|aI0| (reference config)
        aIi = 1.0 / sqrt(norm(a0.col(i)));
        eI = a0.col(i) * aIi;

        // nI = nI0 = eI0 x n0 (reference config)
        nI(0) = eI(1)*nV0(2) - eI(2)*nV0(1);
        nI(1) = eI(2)*nV0(0) - eI(0)*nV0(2);
        nI(2) = eI(0)*nV0(1) - eI(1)*nV0(0);
        nInI = mat_dyad_prod(nI, nI, 3);

      } else { 
        // eI = aI/|aI| (current config)
        aIi = 1.0 / sqrt(norm(a.col(i)));
        eI = a.col(i)*aIi;

        // nI = eI x n (currnt config)
        //
        Vector<double> nI(3);

        nI(0) = eI(1)*nV(2) - eI(2)*nV(1);
        nI(1) = eI(2)*nV(0) - eI(0)*nV(2);
        nI(2) = eI(0)*nV(1) - eI(1)*nV(0);

        cI = norm(a.col(i),a.col(p))*aIi*aIi;
        nInI = mat_dyad_prod(nI, nI, 3);
        eIeI = mat_dyad_prod(eI, eI, 3);
        eIaP = mat_dyad_prod(eI, a.col(p), 3);
      }

      // Update Bb now
      // Free boundary conditions: assumed that the `artificial'
      // triangle is always located in the plane of the main element
      // of the patch
      //
      if (utils::btest(lM.sbc(i,e),enum_int(BoundaryConditionType::bType_free))) {

        // E_I
        //
        if (!lFix[i]) {
          tmpA = -Im + 2.0 * eIeI;
          Bb.rslice(i) = Bb.rslice(i) + mat_mul(Bb.rslice(j), tmpA);
        }

        // E_P
        //
        if (!lFix[p]) {
          tmpA = -2.0 *(cI*Im - 2.0 *cI*eIeI + aIi*eIaP);
          Bb.rslice(p) = Bb.rslice(p) + mat_mul(Bb.rslice(j), tmpA);
        }

        // E_F
        //
        if (!lFix[f]) {
          tmpA = 2.0 * ((1.0-cI)*Im - (1.0-2.0*cI)*eIeI -aIi*eIaP);
          Bb.rslice(f) = Bb.rslice(f) + mat_mul(Bb.rslice(j), tmpA);
        }

        Bb.rslice(j) = 0.0;

        // Hinged boundary conditions: a special case of simple support
        // in which no translation displacements are allowed.

      } else if (utils::btest(lM.sbc(i,e),enum_int(BoundaryConditionType::bType_hing))) {

        // E_I
        if (!lFix[i]) {
          tmpA = -Im + 2.0*eIeI;
          Bb.rslice(i) = Bb.rslice(i) + mat_mul(Bb.rslice(j), tmpA);
        }

	lFix[p] = true; 
        lFix[f] = true; 
        Bb.rslice(p) = 0.0;
        Bb.rslice(f) = 0.0;
        Bb.rslice(j) = 0.0;

      // Fixed boundary condition: no displacements and no rotations
      // are allowed.
      //
      } else if (utils::btest(lM.sbc(i,e),enum_int(BoundaryConditionType::bType_fix))) {

        if (!lFix[i]) {
          tmpA = Im - 2.0*nInI;
          Bb.rslice(i) = Bb.rslice(i) + mat_mul(Bb.rslice(j), tmpA);
        }

        lFix[p] = true;
        lFix[f] = true;
        Bb.rslice(f) = 0.0;
        Bb.rslice(p) = 0.0;

      // Symmetric BCs (need to be verified)
      //
      } else if (utils::btest(lM.sbc(i,e),enum_int(BoundaryConditionType::bType_symm))) {
        if (!lFix[i]) {
          tmpA = Im - 2.0*nInI;
          Bb.rslice(i) = Bb.rslice(i) + mat_mul(Bb.rslice(j), tmpA);
        }

        tmpA.rcol(0) = eI;
        tmpA.rcol(1) = nV;
        tmpA.rcol(2) = nI;

        Bb.rslice(j) = 0.0;
        Bb.rslice(f) = mat_mul(Bb.rslice(f), tmpA);
        Bb.rslice(p) = mat_mul(Bb.rslice(p), tmpA);

        for (int i = 0; i < 3; i++) {
          Bb(i,2,f) = 0.0;
          Bb(i,2,p) = 0.0;
        }
      }
    }
  }
}
 
//----------
// shell_bf
//----------
// Set follower pressure load/net traction on shells. The traction
// on shells is treated as body force and the subroutine is called
// from BF.f
//
// Reproduces Fortran SHELLBF.
//
void shell_bf(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx, 
    const Array<double>& dl, const Array<double>& xl, const Array<double>& tfl, Array<double>& lR, Array3<double>& lK)
{
  using namespace consts;

  const int nsd = com_mod.nsd;
  const int dof = com_mod.dof;
  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  auto cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  double afl = eq.af * eq.beta * dt * dt;

  int i = eq.s;
  int j = i + 1;
  int k = j + 1;

  // Get the current configuration and traction vector
  //
  double tfn = 0.0;
  Array<double> xc(3,eNoN);
  // [NOTE] This is a hack for enabling 'tfl' to be used as a vector in the Fortran.
  auto tfl_data = tfl.data();

  for (int a = 0; a < eNoN; a++) {
    xc(0,a) = xl(0,a) + dl(i,a);
    xc(1,a) = xl(1,a) + dl(j,a);
    xc(2,a) = xl(2,a) + dl(k,a);

    tfn = tfn + N(a)*tfl_data[a];
  }

  double wl = w * tfn;

  // Covariant and contravariant bases in current config
  Array<double> gCov(3,2), gCnv(3,2);
  Vector<double> nV(3);
  nn::gnns(nsd, eNoN, Nx, xc, nV, gCov, gCnv);

  //  Local residual
  for (int a = 0; a < eNoN; a++) {
    lR(0,a) = lR(0,a) - wl*N(a)*nV(0);
    lR(1,a) = lR(1,a) - wl*N(a)*nV(1);
    lR(2,a) = lR(2,a) - wl*N(a)*nV(2);
  }

  // Local stiffness: mass matrix and stiffness contribution due to
  // follower traction load
  //
  double T1 = afl*wl*0.50;

  for (int b = 0; b < eNoN; b++) {
    for (int a = 0; a < eNoN; a++) {
      auto lKp = gCov.rcol(0)*(N(b)*Nx(1,a) - N(a)*Nx(1,b)) - 
          gCov.rcol(1)*(N(b)*Nx(0,a) - N(a)*Nx(0,b));

      lK(1,a,b) = lK(1,a,b) - T1*lKp(2);
      lK(2,a,b) = lK(2,a,b) + T1*lKp(1);

      lK(dof+0,a,b) = lK(dof+0,a,b) + T1*lKp(2);
      lK(dof+2,a,b) = lK(dof+2,a,b) - T1*lKp(0);

      lK(2*dof+0,a,b) = lK(2*dof+0,a,b) - T1*lKp(1);
      lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + T1*lKp(0);
    }
  }
}

//-----------
// shell_cst
//-----------
// Construct shell mechanics for constant strain triangle elements
//
// Note that for triangular elements, eNoN=6 and lM.eNoN=3
//
// Reproduces Fortran SHELLCST.
//
void shell_cst(ComMod& com_mod, const mshType& lM, const int e, const int eNoN, const int nFn, const Array<double>& fN,  
    const Array<double>& al, const Array<double>& yl, const Array<double>& dl, const Array<double>& xl, 
    const Array<double>& bfl, const Vector<int>& ptr)
{
  using namespace consts;

  const int nsd = com_mod.nsd;
  const int dof = com_mod.dof;
  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  auto cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  #define n_debug_shell_cst
  #ifdef debug_shell_cst 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "lM.nFn: " << lM.nFn;
  dmsg << "lM.nFs: " << lM.nFs;
  dmsg << "dof: " << dof;
  dmsg << "eNoN: " << eNoN;
  #endif

  // Define parameters
  double rho = eq.dmn[cDmn].prop.at(PhysicalProperyType::solid_density);
  double dmp = dmn.prop.at(PhysicalProperyType::damping);
  double ht = eq.dmn[cDmn].prop.at(PhysicalProperyType::shell_thickness);
  Vector<double> fb({dmn.prop.at(PhysicalProperyType::f_x), 
                     dmn.prop.at(PhysicalProperyType::f_y), 
                     dmn.prop.at(PhysicalProperyType::f_z)});
  double amd = eq.am * rho  +  eq.af * eq.gam * dt * dmp;
  double afl = eq.af * eq.beta * dt * dt;

  int i = eq.s;
  int j = i + 1;
  int k = j + 1;

  #ifdef debug_shell_cst 
  dmsg << "rho: " << rho;
  dmsg << "dmp: " << dmp;
  dmsg << "ht: " << ht;
  dmsg << "i: " << i;
  dmsg << "j: " << j;
  dmsg << "k: " << k;
  #endif

  //  Get the reference configuration
  auto x0 = xl;

  // Get the current configuration
  //
  Array<double> xc(3,eNoN);

  for (int a = 0; a < eNoN; a++) {
    xc(0,a) = x0(0,a) + dl(i,a);
    xc(1,a) = x0(1,a) + dl(j,a);
    xc(2,a) = x0(2,a) + dl(k,a);
  }

  auto Nx = lM.Nx.slice(0);

  // Covariant and contravariant bases in reference config
  //
  Array<double> tmpX(nsd,lM.eNoN);

  for (int i = 0; i < nsd; i++) {
    for (int j = 0; j < lM.eNoN; j++) {
      tmpX(i,j) = x0(i,j);
    }
  }
 
  Array<double> aCov0(3,2), aCnv0(3,2);
  Vector<double> nV0(3);

  nn::gnns(nsd, lM.eNoN, Nx, tmpX, nV0, aCov0, aCnv0);
  //CALL GNNS(lM%eNoN, Nx, tmpX, nV0, aCov0, aCnv0)

  double Jac0 = sqrt(utils::norm(nV0));
  nV0 = nV0 / Jac0;

  // Covariant and contravariant bases in current config

  for (int i = 0; i < nsd; i++) {
    for (int j = 0; j < lM.eNoN; j++) {
      tmpX(i,j) = xc(i,j);
    }
  }

  Array<double> aCov(3,2), aCnv(3,2);
  Vector<double> nV(3);

  nn::gnns(nsd, lM.eNoN, Nx, tmpX, nV, aCov, aCnv);

  double Jac = sqrt(utils::norm(nV));
  nV = nV / Jac;

  // Compute metric tensor in reference and current config
  //
  double aa_0[2][2] = {};
  double aa_x[2][2] = {};

  for (int g = 0; g < nsd; g++) {
    aa_0[0][0] = aa_0[0][0] + aCov0(g,0)*aCov0(g,0);
    aa_0[0][1] = aa_0[0][1] + aCov0(g,0)*aCov0(g,1);
    aa_0[1][0] = aa_0[1][0] + aCov0(g,1)*aCov0(g,0);
    aa_0[1][1] = aa_0[1][1] + aCov0(g,1)*aCov0(g,1);

    aa_x[0][0] = aa_x[0][0] + aCov(g,0)*aCov(g,0);
    aa_x[0][1] = aa_x[0][1] + aCov(g,0)*aCov(g,1);
    aa_x[1][0] = aa_x[1][0] + aCov(g,1)*aCov(g,0);
    aa_x[1][1] = aa_x[1][1] + aCov(g,1)*aCov(g,1);
  }

  // Compute fiber orientation in curvature coordinates
  //
  Array<double> fNa0(2,nFn);
  //dmsg << "nFn: " << nFn;
  //dmsg << "fN: " << fN;
  //dmsg << "aCnv0: " << aCnv0;

  for (int iFn = 0; iFn < nFn; iFn++) {
    for (int l = 0; l < 3; l++) {
      fNa0(0,iFn) = fNa0(0,iFn) + fN(l,iFn)*aCnv0(l,0);
      fNa0(1,iFn) = fNa0(1,iFn) + fN(l,iFn)*aCnv0(l,1);
    }
  }
  //dmsg << "fNa0: " << fNa0;

  // Define variation in membrane strain only for the main element
  //
  Array3<double> Bm(3,3,lM.eNoN);

  for (int a = 0; a < lM.eNoN; a++) {
     Bm(0,0,a) = Nx(0,a)*aCov(0,0);
     Bm(0,1,a) = Nx(0,a)*aCov(1,0);
     Bm(0,2,a) = Nx(0,a)*aCov(2,0);

     Bm(1,0,a) = Nx(1,a)*aCov(0,1);
     Bm(1,1,a) = Nx(1,a)*aCov(1,1);
     Bm(1,2,a) = Nx(1,a)*aCov(2,1);

     Bm(2,0,a) = Nx(1,a)*aCov(0,0) + Nx(0,a)*aCov(0,1);
     Bm(2,1,a) = Nx(1,a)*aCov(1,0) + Nx(0,a)*aCov(1,1);
     Bm(2,2,a) = Nx(1,a)*aCov(2,0) + Nx(0,a)*aCov(2,1);
  }
  //dmsg << " " << " ";
  //dmsg << "Bm: " << Bm;
  //exit(0);

  // For the boundary elements, zero-out Bm for fixed/clamped BC.
  //
  std::vector<bool> setIt = {false, false, false};
  int a = lM.eNoN;

  while (a < eNoN) {
    //dmsg << "----- a: " << a+1;
    if (ptr(a) == -1) {
      int b = a - lM.eNoN;
      //dmsg << "b: " << b+1;
      if (utils::btest(lM.sbc(b,e),enum_int(BoundaryConditionType::bType_fix))) {
        //dmsg << "set to true " << b+1;
        setIt[b] = true;
      }
    }

    a = a + 1;
  }

  for (int a = 0; a < lM.eNoN; a++) {
    if (setIt[a]) {
      for (int b = 0; b < lM.eNoN; b++) {
        if (a == b) {
          continue; 
        }
        Bm.rslice(b) = 0.0;
        //dmsg << "Zero slice b: " << b;
      }
    }
  }
  //dmsg << " " << " ";
  //dmsg << "Bm: " << Bm;
  //exit(0);

  // Compute curvature coefficients for bending strain and its
  // variation for CST elements
  //
  double bb_0[2][2], bb_x[2][2];
  Array3<double> Bb(3,3,eNoN);
  shell_bend_cst(com_mod, lM, e, ptr, x0, xc, bb_0, bb_x, Bb, true);
  //CALL SHELLBENDCST(lM, e, ptr, x0, xc, bb_0, bb_x, Bb, .TRUE.)
  //dmsg << "Bb: " << Bb;

  // Compute stress resultants by integrating 2nd Piola Kirchhoff
  // stress and elasticity tensors through the shell thickness. These
  // resultants are computed in Voigt notation.
  //
  /*
  dmsg << "aa_0: " << aa_0[0][0];
  dmsg << "    : " << aa_0[0][1];
  dmsg << "    : " << aa_0[1][0];
  dmsg << "    : " << aa_0[1][1];
  dmsg << "aa_x: " << aa_x[0][0];
  dmsg << "    : " << aa_x[0][1];
  dmsg << "    : " << aa_x[1][0];
  dmsg << "    : " << aa_x[1][1];
  */
  Array3<double> Dm(3,3,3); 
  Array<double> Sm(3,2);
  double lam3;
  shl_strs_res(com_mod, dmn, nFn, fNa0, aa_0, aa_x, bb_0, bb_x, lam3, Sm, Dm);

  /*
  dmsg << " " << " ";
  dmsg << "Bm: " << Bm;
  dmsg << "Bb: " << Bb;

  dmsg << " " << " ";
  dmsg << "Sm: " << Sm;
  dmsg << "Dm: " << Dm;
  exit(0);
  */

  // Contribution to tangent matrices: Dm * Bm, Dm*Bb
  //
  Array3<double> D0Bm(3,3,lM.eNoN), D1Bm(3,3,lM.eNoN);

  for (int a = 0; a < lM.eNoN; a++) {
    D0Bm.rslice(a) = mat_fun::mat_mul(Dm.rslice(0), Bm.rslice(a));
    D1Bm.rslice(a) = mat_fun::mat_mul(Dm.rslice(1), Bm.rslice(a));
  }

  Array3<double> D1Bb(3,3,eNoN), D2Bb(3,3,eNoN);

  for (int a = 0; a < eNoN; a++) {
    D1Bb.rslice(a) = mat_fun::mat_mul(Dm.rslice(1), Bb.rslice(a));
    D2Bb.rslice(a) = mat_fun::mat_mul(Dm.rslice(2), Bb.rslice(a));
  }

  // Contribution to residual and stiffness matrices due to inertia and
  // body forces
  //
  Array<double> lR(dof,eNoN);
  Array3<double> lK(dof*dof,eNoN,eNoN);

  for (int g = 0; g < lM.nG; g++) {
    auto N = lM.N.rcol(g);
    double w = lM.w(g) * Jac0 * ht;

    // Acceleration and mass damping at the integration point
    //
    auto ud = -fb;

    for (int a = 0; a < lM.eNoN; a++) {
      ud(0) = ud(0) + N(a)*(rho*(al(i,a)-bfl(0,a)) + dmp*yl(i,a));
      ud(1) = ud(1) + N(a)*(rho*(al(j,a)-bfl(1,a)) + dmp*yl(j,a));
      ud(2) = ud(2) + N(a)*(rho*(al(k,a)-bfl(2,a)) + dmp*yl(k,a));
     }

     // Local residual
     for (int a = 0; a < lM.eNoN; a++) {
       lR(0,a) = lR(0,a) + N(a)*w*ud(0);
       lR(1,a) = lR(1,a) + N(a)*w*ud(1);
       lR(2,a) = lR(2,a) + N(a)*w*ud(2);
     }

     //  Local stiffness contribution from mass matrix
     for (int b = 0; b < lM.eNoN; b++) {
       for (int a = 0; a < lM.eNoN; a++) {
         double BtS = w*amd*N(a)*N(b);
         lK(0,a,b) = lK(0,a,b) + BtS;
         lK(dof+1,a,b) = lK(dof+1,a,b) + BtS;
         lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + BtS;
       }
     }
  }
  //std::cout << "[shell_cst]  " << " " << std::endl;
  //std::cout << "[shell_cst] lK: " << lK << std::endl;
  //exit(0);

  // Contribution to residual from membrane strain
  //
  double w = Jac0 * 0.5;

  for (int a = 0; a < lM.eNoN; a++) {
    double BtS = Bm(0,0,a)*Sm(0,0) + Bm(1,0,a)*Sm(1,0) + Bm(2,0,a)*Sm(2,0);
    lR(0,a) = lR(0,a) + w*BtS;

    BtS = Bm(0,1,a)*Sm(0,0) + Bm(1,1,a)*Sm(1,0) + Bm(2,1,a)*Sm(2,0);
    lR(1,a) = lR(1,a) + w*BtS;

    BtS = Bm(0,2,a)*Sm(0,0) + Bm(1,2,a)*Sm(1,0) + Bm(2,2,a)*Sm(2,0);
    lR(2,a) = lR(2,a) + w*BtS;
  }

  // Contribution to residual from bending strain
  //
  for (int a = 0; a < eNoN; a++) {
     double BtS = Bb(0,0,a)*Sm(0,1) + Bb(1,0,a)*Sm(1,1) + Bb(2,0,a)*Sm(2,1);
     lR(0,a) = lR(0,a) + w*BtS;

     BtS = Bb(0,1,a)*Sm(0,1) + Bb(1,1,a)*Sm(1,1) + Bb(2,1,a)*Sm(2,1);
     lR(1,a) = lR(1,a) + w*BtS;

     BtS = Bb(0,2,a)*Sm(0,1) + Bb(1,2,a)*Sm(1,1) + Bb(2,2,a)*Sm(2,1);
     lR(2,a) = lR(2,a) + w*BtS;
  }

  // Contribution to stiffness from membrane-membrane interactions
  w = afl * Jac0 * 0.5;
  //dmsg << "w: " << w;

  for (int b = 0; b < lM.eNoN; b++) {
    for (int a = 0; a < lM.eNoN; a++) {
      double NxSNx = Nx(0,a)*Nx(0,b)*Sm(0,0) + Nx(1,a)*Nx(1,b)*Sm(1,0) + 
                     Nx(0,a)*Nx(1,b)*Sm(2,0) + Nx(1,a)*Nx(0,b)*Sm(2,0);

      double BtDB = Bm(0,0,a)*D0Bm(0,0,b) + Bm(1,0,a)*D0Bm(1,0,b) + Bm(2,0,a)*D0Bm(2,0,b);
      lK(0,a,b) = lK(0,a,b) + w*(BtDB + NxSNx);

      BtDB = Bm(0,0,a)*D0Bm(0,1,b) + Bm(1,0,a)*D0Bm(1,1,b) + Bm(2,0,a)*D0Bm(2,1,b);
      lK(1,a,b) = lK(1,a,b) + w*BtDB;

      BtDB = Bm(0,0,a)*D0Bm(0,2,b) + Bm(1,0,a)*D0Bm(1,2,b) + Bm(2,0,a)*D0Bm(2,2,b);
      lK(2,a,b) = lK(2,a,b) + w*BtDB;

      BtDB = Bm(0,1,a)*D0Bm(0,0,b) + Bm(1,1,a)*D0Bm(1,0,b) + Bm(2,1,a)*D0Bm(2,0,b);
      lK(dof+0,a,b) = lK(dof+0,a,b) + w*BtDB;

      BtDB = Bm(0,1,a)*D0Bm(0,1,b) + Bm(1,1,a)*D0Bm(1,1,b) + Bm(2,1,a)*D0Bm(2,1,b);
      lK(dof+1,a,b) = lK(dof+1,a,b) + w*(BtDB + NxSNx);

      BtDB = Bm(0,1,a)*D0Bm(0,2,b) + Bm(1,1,a)*D0Bm(1,2,b) + Bm(2,1,a)*D0Bm(2,2,b);
      lK(dof+2,a,b) = lK(dof+2,a,b) + w*BtDB;

      BtDB = Bm(0,2,a)*D0Bm(0,0,b) + Bm(1,2,a)*D0Bm(1,0,b) + Bm(2,2,a)*D0Bm(2,0,b);
      lK(2*dof+0,a,b) = lK(2*dof+0,a,b) + w*BtDB;

      BtDB = Bm(0,2,a)*D0Bm(0,1,b) + Bm(1,2,a)*D0Bm(1,1,b) + Bm(2,2,a)*D0Bm(2,1,b);
      lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + w*BtDB;

      BtDB = Bm(0,2,a)*D0Bm(0,2,b) + Bm(1,2,a)*D0Bm(1,2,b) + Bm(2,2,a)*D0Bm(2,2,b);
      lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + w*(BtDB + NxSNx);
    }
  }
  //std::cout << "[shell_cst]  " << " " << std::endl;
  //std::cout << "[shell_cst] lK: " << lK << std::endl;
  //exit(0);

  //  Contribution to stiffness from membrane-bending interactions
  //
  for (int b = 0; b < eNoN; b++) {
    for (int a = 0; a < lM.eNoN; a++) {
      double BtDB = Bm(0,0,a)*D1Bb(0,0,b) + Bm(1,0,a)*D1Bb(1,0,b) + Bm(2,0,a)*D1Bb(2,0,b);
      lK(0,a,b) = lK(0,a,b) + w*BtDB;

      BtDB = Bm(0,0,a)*D1Bb(0,1,b) + Bm(1,0,a)*D1Bb(1,1,b) + Bm(2,0,a)*D1Bb(2,1,b);
      lK(1,a,b) = lK(1,a,b) + w*BtDB;

      BtDB = Bm(0,0,a)*D1Bb(0,2,b) + Bm(1,0,a)*D1Bb(1,2,b) + Bm(2,0,a)*D1Bb(2,2,b);
      lK(2,a,b) = lK(2,a,b) + w*BtDB;

      BtDB = Bm(0,1,a)*D1Bb(0,0,b) + Bm(1,1,a)*D1Bb(1,0,b) + Bm(2,1,a)*D1Bb(2,0,b);
      lK(dof+0,a,b) = lK(dof+0,a,b) + w*BtDB;

      BtDB = Bm(0,1,a)*D1Bb(0,1,b) + Bm(1,1,a)*D1Bb(1,1,b) + Bm(2,1,a)*D1Bb(2,1,b);
      lK(dof+1,a,b) = lK(dof+1,a,b) + w*BtDB;

      BtDB = Bm(0,1,a)*D1Bb(0,2,b) + Bm(1,1,a)*D1Bb(1,2,b) + Bm(2,1,a)*D1Bb(2,2,b);
      lK(dof+2,a,b) = lK(dof+2,a,b) + w*BtDB;

      BtDB = Bm(0,2,a)*D1Bb(0,0,b) + Bm(1,2,a)*D1Bb(1,0,b) + Bm(2,2,a)*D1Bb(2,0,b);
      lK(2*dof+0,a,b) = lK(2*dof+0,a,b) + w*BtDB;

      BtDB = Bm(0,2,a)*D1Bb(0,1,b) + Bm(1,2,a)*D1Bb(1,1,b) + Bm(2,2,a)*D1Bb(2,1,b);
      lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + w*BtDB;

      BtDB = Bm(0,2,a)*D1Bb(0,2,b) + Bm(1,2,a)*D1Bb(1,2,b) + Bm(2,2,a)*D1Bb(2,2,b);
      lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + w*BtDB;
    }
  }

  // Contribution to stiffness from bending-membrane interactions
  //
  for (int b = 0; b < lM.eNoN; b++) {
    for (int a = 0; a < eNoN; a++) {
      double BtDB = Bb(0,0,a)*D1Bm(0,0,b) + Bb(1,0,a)*D1Bm(1,0,b) + Bb(2,0,a)*D1Bm(2,0,b);
      lK(0,a,b) = lK(0,a,b) + w*BtDB;

      BtDB = Bb(0,0,a)*D1Bm(0,1,b) + Bb(1,0,a)*D1Bm(1,1,b) + Bb(2,0,a)*D1Bm(2,1,b);
      lK(1,a,b) = lK(1,a,b) + w*BtDB;

      BtDB = Bb(0,0,a)*D1Bm(0,2,b) + Bb(1,0,a)*D1Bm(1,2,b) + Bb(2,0,a)*D1Bm(2,2,b);
      lK(2,a,b) = lK(2,a,b) + w*BtDB;

      BtDB = Bb(0,1,a)*D1Bm(0,0,b) + Bb(1,1,a)*D1Bm(1,0,b) + Bb(2,1,a)*D1Bm(2,0,b);
      lK(dof+0,a,b) = lK(dof+0,a,b) + w*BtDB;

      BtDB = Bb(0,1,a)*D1Bm(0,1,b) + Bb(1,1,a)*D1Bm(1,1,b) + Bb(2,1,a)*D1Bm(2,1,b);
      lK(dof+1,a,b) = lK(dof+1,a,b) + w*BtDB;

      BtDB = Bb(0,1,a)*D1Bm(0,2,b) + Bb(1,1,a)*D1Bm(1,2,b) + Bb(2,1,a)*D1Bm(2,2,b);
      lK(dof+2,a,b) = lK(dof+2,a,b) + w*BtDB;

      BtDB = Bb(0,2,a)*D1Bm(0,0,b) + Bb(1,2,a)*D1Bm(1,0,b) + Bb(2,2,a)*D1Bm(2,0,b);
      lK(2*dof+0,a,b) = lK(2*dof+0,a,b) + w*BtDB;

      BtDB = Bb(0,2,a)*D1Bm(0,1,b) + Bb(1,2,a)*D1Bm(1,1,b) + Bb(2,2,a)*D1Bm(2,1,b);
      lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + w*BtDB;

      BtDB = Bb(0,2,a)*D1Bm(0,2,b) + Bb(1,2,a)*D1Bm(1,2,b) + Bb(2,2,a)*D1Bm(2,2,b);
      lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + w*BtDB;
     }
  }

  // Contribution to stiffness from bending-bending interactions
  //
  for (int b = 0; b < eNoN; b++) {
    for (int a = 0 ; a < eNoN; a++) {
      double BtDB = Bb(0,0,a)*D2Bb(0,0,b) + Bb(1,0,a)*D2Bb(1,0,b) + Bb(2,0,a)*D2Bb(2,0,b);
      lK(0,a,b) = lK(0,a,b) + w*BtDB;

      BtDB = Bb(0,0,a)*D2Bb(0,1,b) + Bb(1,0,a)*D2Bb(1,1,b) + Bb(2,0,a)*D2Bb(2,1,b);
      lK(1,a,b) = lK(1,a,b) + w*BtDB;

      BtDB = Bb(0,0,a)*D2Bb(0,2,b) + Bb(1,0,a)*D2Bb(1,2,b) + Bb(2,0,a)*D2Bb(2,2,b);
      lK(2,a,b) = lK(2,a,b) + w*BtDB;

      BtDB = Bb(0,1,a)*D2Bb(0,0,b) + Bb(1,1,a)*D2Bb(1,0,b) + Bb(2,1,a)*D2Bb(2,0,b);
      lK(dof+0,a,b) = lK(dof+0,a,b) + w*BtDB;

      BtDB = Bb(0,1,a)*D2Bb(0,1,b) + Bb(1,1,a)*D2Bb(1,1,b) + Bb(2,1,a)*D2Bb(2,1,b);
      lK(dof+1,a,b) = lK(dof+1,a,b) + w*BtDB;

      BtDB = Bb(0,1,a)*D2Bb(0,2,b) + Bb(1,1,a)*D2Bb(1,2,b) + Bb(2,1,a)*D2Bb(2,2,b);
      lK(dof+2,a,b) = lK(dof+2,a,b) + w*BtDB;

      BtDB = Bb(0,2,a)*D2Bb(0,0,b) + Bb(1,2,a)*D2Bb(1,0,b) + Bb(2,2,a)*D2Bb(2,0,b);
      lK(2*dof+0,a,b) = lK(2*dof+0,a,b) + w*BtDB;

      BtDB = Bb(0,2,a)*D2Bb(0,1,b) + Bb(1,2,a)*D2Bb(1,1,b) + Bb(2,2,a)*D2Bb(2,1,b);
      lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + w*BtDB;

      BtDB = Bb(0,2,a)*D2Bb(0,2,b) + Bb(1,2,a)*D2Bb(1,2,b) + Bb(2,2,a)*D2Bb(2,2,b);
      lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + w*BtDB;
    }
  }

  eq.linear_algebra->assemble(com_mod, eNoN, ptr, lK, lR);
}

//----------
// shell_fp
//----------
//
void shell_fp(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx, 
    const Array<double>& dl, const Array<double>& xl, const Array<double>& tfl, Array<double>& lR, Array3<double>& lK)
{
  int nsd = com_mod.nsd;
  int dof = com_mod.dof;
  double dt = com_mod.dt;
  auto& cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];

  double afl = eq.af * eq.beta*dt*dt;
  int i = eq.s;
  int j = i + 1;
  int k = j + 1;

  // Get the current configuration and traction vector
  //
  Array<double> xc(3,eNoN);
  double tfn = 0.0;
  // [NOTE] This is a hack enabling 'tfl' to be used
  // like a vector in the Fortan.
  auto tfl_data = tfl.data();

  for (int a = 0; a < eNoN; a++) {
    xc(0,a) = xl(0,a) + dl(i,a);
    xc(1,a) = xl(1,a) + dl(j,a);
    xc(2,a) = xl(2,a) + dl(k,a);
    tfn = tfn + N(a)*tfl_data[a];
  }

  double wl = w * tfn;

  // Covariant and contravariant bases in current config
  Vector<double> nV(3); 
  Array<double> gCov(3,2);
  Array<double> gCnv(3,2);
  nn::gnns(nsd, eNoN, Nx, xc, nV, gCov, gCnv);

  // Local residual
  for (int a = 0; a < eNoN; a++) {
    lR(0,a) = lR(0,a) - wl*N(a)*nV(0);
    lR(1,a) = lR(1,a) - wl*N(a)*nV(1);
    lR(2,a) = lR(2,a) - wl*N(a)*nV(2);
  }

  // Local stiffness: mass matrix and stiffness contribution due to
  // follower traction load
  //
  double T1 = afl*wl*0.5;

  for (int b = 0; b < eNoN; b++) {
    for (int a = 0; a < eNoN; a++) {
      auto lKp = gCov.col(0) * (N(b)*Nx(1,a) - N(a)*Nx(1,b)) - gCov.col(1)*(N(b)*Nx(0,a) - N(a)*Nx(0,b));

      lK(1,a,b) = lK(1,a,b) - T1*lKp(2);
      lK(2,a,b) = lK(2,a,b) + T1*lKp(1);

      lK(dof+0,a,b) = lK(dof+0,a,b) + T1*lKp(2);
      lK(dof+2,a,b) = lK(dof+2,a,b) - T1*lKp(0);

      lK(2*dof,a,b) = lK(2*dof,a,b) - T1*lKp(1);
      lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + T1*lKp(0);
    }
  }
}

//--------------
// shl_strs_res
//--------------
// Compute stress resultants for shell elements
//
// Reproduces Fortan SHL_STRS_RES.
//
void shl_strs_res(const ComMod& com_mod, const dmnType& lDmn, const int nFn, const Array<double>& fNa0, 
    const double aa_0[2][2], const double aa_x[2][2], const double bb_0[2][2], const double bb_x[2][2], 
    double& lam3, Array<double>& Sm, Array3<double>& Dm)
{
  using namespace consts;

  #define n_debug_shl_strs_res
  #ifdef debug_shl_strs_res
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  // Set shell thickness
  double ht = lDmn.prop.at(PhysicalProperyType::shell_thickness); 
  double nu = lDmn.prop.at(PhysicalProperyType::poisson_ratio); 

  // Check for incompressibility
  bool flag = false;
  if (utils::is_zero(nu-0.50)) {
    flag = true;
  }
  //dmsg << "flag: " << flag;

  // Set integration parameters (Gauss coordinates and weights)
  double wh[3] = { 5.0/9.0, 5.0/9.0, 8.0/9.0 };
  double xi[3] = { -sqrt(0.60), sqrt(0.60), 0.0 };

  //  Initialize stress-resultants
  Sm = 0.0;
  Dm = 0.0;

  Vector<double> Sml(3);
  Array<double> Dml(3,3);

  // Averaged SQRT(g33) over the thickness
  lam3 = 0.0;

  // Gauss integration through shell thickness
  //
  Array<double> gg_0(2,2), gg_x(2,2);

  for (int g = 0; g < 3; g++) { 
    //dmsg << "---------- g: " << g+1;
    // Local shell thickness coordinate
    double xis = 0.50 * ht * xi[g];

    // Metric coefficients in shell continuum (ref/cur)
    for (int i = 0; i < 2; i++) { 
      for (int j = 0; j < 2; j++) { 
        gg_0(i,j) = aa_0[i][j] - 2.0*xis*bb_0[i][j];
        gg_x(i,j) = aa_x[i][j] - 2.0*xis*bb_x[i][j];
      }
    }
    //dmsg << "gg_0: " << gg_0;
    //dmsg << "gg_x: " << gg_x;

    // Get 2nd Piola-Kirchhoff and elasticity tensors
    //
    double g33;

    // For incompressible materials
    if (flag) {
      mat_models::compute_pk2cc_shli(com_mod, lDmn, nFn, fNa0, gg_0, gg_x, g33, Sml, Dml);

    // For compressible materials
    } else { 
      mat_models::compute_pk2cc_shlc(com_mod, lDmn, nFn, fNa0, gg_0, gg_x, g33, Sml, Dml);
    }

    //dmsg << "      " << " ";
    //dmsg << "g33: " << g33;
    //dmsg << "Sml: " << Sml;
    //dmsg << "Dml: " << Dml;

    double wl[3];
    wl[0] = 0.50 * wh[g] * ht;
    wl[1] = wl[0] * xis;
    wl[2] = wl[1] * xis;

    lam3  = lam3 + wl[0] * sqrt(g33);

    for (int i = 0; i < 3; i++) { 
      Sm(i,0) = Sm(i,0) + wl[0]*Sml(i);
      Sm(i,1) = Sm(i,1) + wl[1]*Sml(i);

      for (int j = 0; j < 3; j++) { 
        Dm(i,j,0) = Dm(i,j,0) + wl[0]*Dml(i,j);
        Dm(i,j,1) = Dm(i,j,1) + wl[1]*Dml(i,j);
        Dm(i,j,2) = Dm(i,j,2) + wl[2]*Dml(i,j);
      }
    }
  }

  lam3 = lam3 / ht;

  //dmsg << "Sm: " << Sm;
  //dmsg << "Dm: " << Dm;
  //dmsg << "lam3: " << lam3;
}

};

