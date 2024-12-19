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

/**
 * @file ustruct.cpp
 * @brief Structural mechanics implementation based on the following reference:
 *
 * Ju Liu, Alison L. Marsden,
 * A unified continuum and variational multiscale formulation for fluids, solids, and fluid–structure interaction,
 * Computer Methods in Applied Mechanics and Engineering,
 * Volume 337,
 * 2018,
 * Pages 549-597,
 * ISSN 0045-7825,
 * https://doi.org/10.1016/j.cma.2018.03.045.
 *
 * This paper describes a unified framework for fluid, solids, and FSI. The code
 * in this file is based on the solid mechanics portion of the paper, which can be
 * found in Section 4.
 */

#include "ustruct.h"

#include "all_fun.h"
#include "fs.h"
#include "mat_fun.h"
#include "mat_models.h"
#include "nn.h"
#include "utils.h"

#include <math.h>

namespace ustruct {

void b_ustruct_2d(const ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, 
    const Array<double>& Nx, const Array<double>& dl, const Vector<double>& hl, const Vector<double>& nV, 
    Array<double>& lR, Array3<double>& lK, Array3<double>& lKd)
{
  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  double dt = com_mod.dt;

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
      lKd(1,a,b) = lKd(1,a,b) + Ku;
      lK(1,a,b) = lK(1,a,b) + afm*Ku;

      lKd(2,a,b) = lKd(2,a,b) - Ku;
      lK(3,a,b)  = lK(3,a,b)  - afm*Ku;
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
/// @param lKd Local stiffness matrix (displacement)
void b_ustruct_3d(const ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, 
    const Array<double>& Nx, const Array<double>& dl, const Vector<double>& hl, const Vector<double>& nV, 
    Array<double>& lR, Array3<double>& lK, Array3<double>& lKd)
{
  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  double dt = com_mod.dt;

  double af = eq.af * eq.gam*dt;
  double afm = af / eq.am;
  int i = eq.s;
  int j = i + 1;
  int k = j + 1;

  Vector<double> nFi(3);  
  Array<double> NxFi(3,eNoN);

  Array<double> F(3,3); 
  F(0,0) = 1.0;
  F(1,1) = 1.0;
  F(2,2) = 1.0;

  double h = 0.0;

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

  nFi(0) = nV(0)*Fi(0,0) + nV(1)*Fi(1,0) + nV(2)*Fi(2,0);
  nFi(1) = nV(0)*Fi(0,1) + nV(1)*Fi(1,1) + nV(2)*Fi(2,1);
  nFi(2) = nV(0)*Fi(0,2) + nV(1)*Fi(1,2) + nV(2)*Fi(2,2);

  double wl = w * Jac * h;

  for (int a = 0; a  < eNoN; a++) {
    lR(0,a) = lR(0,a) - wl*N(a)*nFi(0);
    lR(1,a) = lR(1,a) - wl*N(a)*nFi(1);
    lR(2,a) = lR(2,a) - wl*N(a)*nFi(2);

    for (int b = 0; b < eNoN; b++) {
      double Ku = wl*af*N(a)*(nFi(1)*NxFi(0,b) - nFi(0)*NxFi(1,b));
      lKd(1,a,b) = lKd(1,a,b) + Ku;
      lK(1,a,b) = lK(1,a,b) + afm*Ku;

      lKd(3,a,b) = lKd(3,a,b) - Ku;
      lK(4,a,b)  = lK(4,a,b)  - afm*Ku;

      Ku = wl*af*N(a)*(nFi(2)*NxFi(0,b) - nFi(0)*NxFi(2,b));
      lKd(2,a,b) = lKd(2,a,b) + Ku;
      lK(2,a,b) = lK(2,a,b)  + afm*Ku;

      lKd(6,a,b) = lKd(6,a,b) - Ku;
      lK(8,a,b)  = lK(8,a,b)  - afm*Ku;

      Ku = wl*af*N(a)*(nFi(2)*NxFi(1,b) - nFi(1)*NxFi(2,b));
      lKd(5,a,b) = lKd(5,a,b) + Ku;
      lK(6,a,b)  = lK(6,a,b)  + afm*Ku;

      lKd(7,a,b) = lKd(7,a,b) - Ku;
      lK(9,a,b) = lK(9,a,b) - afm*Ku;
    }
  }
}

/// @brief Reproduces Fortran CONSTRUCT_uSOLID.
//
void construct_usolid(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const Array<double>& Ag, 
    const Array<double>& Yg, const Array<double>& Dg)
{
  using namespace consts;

  #define n_debug_construct_usolid
  #ifdef debug_construct_usolid
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "lM.nFn: " << lM.nFn;
  dmsg << "lM.nFs: " << lM.nFs;
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

  bool vmsStab = false;

  if (lM.nFs == 1) {
    vmsStab = true;
  } else {
    vmsStab = false;
  }

  #ifdef debug_construct_usolid
  dmsg << "nsd: " << nsd ;
  dmsg << "vmsStab: " << vmsStab;
  #endif

  // USTRUCT: dof = nsd+1
  Vector<int> ptr(eNoN);
  Vector<double> pSl(nsymd), ya_l(eNoN), N(eNoN);
  Array<double> xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN),
                bfl(nsd,eNoN), fN(nsd,nFn), pS0l(nsymd,eNoN), Nx(nsd,eNoN), lR(dof,eNoN);
  Array3<double> lK(dof*dof,eNoN,eNoN), lKd(dof*nsd,eNoN,eNoN);

  for (int e = 0; e < lM.nEl; e++) {
    // Update domain and proceed if domain phys and eqn phys match
    cDmn = all_fun::domain(com_mod, lM, cEq, e);
    auto cPhys = eq.dmn[cDmn].phys;
    if (cPhys != EquationType::phys_ustruct) {
      continue;
    }

    // Create local copies
    fN  = 0.0;
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

      if (cem.cpld) {
        ya_l(a) = cem.Ya(Ac);
      }
    }

    // Initialize residual and tangents
    lR = 0.0;
    lK = 0.0;
    lKd = 0.0;
    std::array<fsType,2> fs;

    // Set function spaces for velocity and pressure.
    fs::get_thood_fs(com_mod, fs, lM, vmsStab, 1);

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
           throw std::runtime_error("[construct_usolid] Jacobian for element " + std::to_string(e) + " is < 0.");
        }
      }

      double w = fs[0].w(g) * Jac;

      if (nsd == 3) {
        auto N0 = fs[0].N.col(g);
        auto N1 = fs[1].N.col(g);
        ustruct_3d_m(com_mod, cep_mod, vmsStab, fs[0].eNoN, fs[1].eNoN, nFn, w, Jac, N0, N1, Nwx, al, yl, dl, bfl, fN, ya_l, lR, lK, lKd);

      } else if (nsd == 2) {
        auto N0 = fs[0].N.col(g);
        auto N1 = fs[1].N.col(g);
        ustruct_2d_m(com_mod, cep_mod, vmsStab, fs[0].eNoN, fs[1].eNoN, nFn, w, Jac, N0, N1, Nwx, al, yl, dl, bfl, fN, ya_l, lR, lK, lKd);
      }

    } // for g = 0 to fs[0].nG

    // Set function spaces for velocity/displacement and pressure.
    fs::get_thood_fs(com_mod, fs, lM, vmsStab, 2);

    // Gauss integration 2
    //
    for (int g = 0; g < fs[1].nG; g++) {
      if (g == 0 || !fs[0].lShpF) {
        auto Nx = fs[0].Nx.slice(g);
        nn::gnn(fs[0].eNoN, nsd, nsd, Nx, xwl, Nwx, Jac, ksix);
        if (utils::is_zero(Jac)) {
           throw std::runtime_error("[construct_usolid] Jacobian for element " + std::to_string(e) + " is < 0.");
        }
      }

      if (g == 0 || !fs[1].lShpF) {
        auto Nx = fs[1].Nx.slice(g);
        nn::gnn(fs[1].eNoN, nsd, nsd, Nx, xql, Nqx, Jac, ksix);
        if (utils::is_zero(Jac)) {
           throw std::runtime_error("[construct_usolid] Jacobian for element " + std::to_string(e) + " is < 0.");
        }
      }

      double w = fs[1].w(g) * Jac;

      if (nsd == 3) {
        auto N0 = fs[0].N.col(g);
        auto N1 = fs[1].N.col(g);
        ustruct_3d_c(com_mod, cep_mod, vmsStab, fs[0].eNoN, fs[1].eNoN, w, Jac, N0, N1, Nwx, 
            Nqx, al, yl, dl, bfl, lR, lK, lKd);

      } else if (nsd == 2) {
        auto N0 = fs[0].N.col(g);
        auto N1 = fs[1].N.col(g);
        ustruct_2d_c(com_mod, cep_mod, vmsStab, fs[0].eNoN, fs[1].eNoN, w, Jac, N0, N1, Nwx, 
            Nqx, al, yl, dl, bfl, lR, lK, lKd);
      }

    } // for g = 0 to fs[1].nG

    ustruct_do_assem(com_mod, eNoN, ptr, lKd, lK, lR);

  } // for e = 0 to lM.nEl

}

int get_col_ptr(ComMod& com_mod, const int rowN, const int colN)
{
  auto& rowPtr = com_mod.rowPtr;
  auto& colPtr = com_mod.colPtr;

  int left = rowPtr(rowN);
  int right = rowPtr(rowN+1);
  int ptr = (right + left) / 2;
 
  while (colN != colPtr(ptr)) {
    if (colN > colPtr(ptr)) {
      left  = ptr;
    } else { 
      right = ptr;
    }
    ptr = (right + left) / 2;
  }

  return ptr;
} 

/// @brief Reproduces Fortran USTRUCT2D_C.
//
void ustruct_2d_c(ComMod& com_mod, CepMod& cep_mod, const bool vmsFlag, const int eNoNw, const int eNoNq,
    const double w, const double Je, const Vector<double>& Nw,  const Vector<double>& Nq,
    const Array<double>& Nwx, const Array<double>& Nqx, const Array<double>& al, const Array<double>& yl, 
    const Array<double>& dl, const Array<double>& bfl, Array<double>& lR, Array3<double>& lK, 
    Array3<double>& lKd)
{
  using namespace consts;
  using namespace mat_fun;

  #define n_debug_ustruct_2d_c
  #ifdef debug_ustruct_2d_c
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "w: " << w;
  dmsg << "eNoNw: " << eNoNw;
  dmsg << "eNoNq: " << eNoNq;
  #endif

  // Define parameters
  //
  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  Vector<double> fb(2);
  fb[0] = dmn.prop[PhysicalProperyType::f_x];
  fb[1] = dmn.prop[PhysicalProperyType::f_y];
  fb[2] = dmn.prop[PhysicalProperyType::f_z];

  double am = eq.am;
  double af = eq.af * eq.gam * dt;
  double afm = af / am;

  // {i,j} := velocity dofs; {k} := pressure dof
  int i = eq.s;
  int j = i + 1;
  int k = j + 1;

  #ifdef debug_ustruct_2d_c
  dmsg << "am: " << am;
  dmsg << "af: " << af;
  dmsg << "afm: " << afm;
  dmsg << "i: " << i;
  #endif

  // Inertia (velocity and acceleration), body force, fiber directions,
  // and deformation tensor (F) at integration point
  //
  Vector<double> vd{-fb[0], -fb[1]};
  Vector<double> v(2);
  Array<double> vx(2,2), F(2,2);
  F(0,0) = 1.0;
  F(1,1) = 1.0;

  for (int a = 0; a < eNoNw; a++) {
    v(0) = v(0) + Nw(a)*yl(i,a);
    v(1) = v(1) + Nw(a)*yl(j,a);

    vd(0) = vd(0) + Nw(a)*(al(i,a)-bfl(0,a));
    vd(1) = vd(1) + Nw(a)*(al(j,a)-bfl(1,a));

    vx(0,0) = vx(0,0) + Nwx(0,a)*yl(i,a);
    vx(0,1) = vx(0,1) + Nwx(1,a)*yl(i,a);
    vx(1,0) = vx(1,0) + Nwx(0,a)*yl(j,a);
    vx(1,1) = vx(1,1) + Nwx(1,a)*yl(j,a);

    F(0,0) = F(0,0) + Nwx(0,a)*dl(i,a);
    F(0,1) = F(0,1) + Nwx(1,a)*dl(i,a);
    F(1,0) = F(1,0) + Nwx(0,a)*dl(j,a);
    F(1,1) = F(1,1) + Nwx(1,a)*dl(j,a);
  }

  double Jac = mat_fun::mat_det(F, 2);
  auto Fi = mat_fun::mat_inv(F, 2);

  // Pressure and its gradients 
  //
  double p = 0.0;
  double pd = 0.0;
  Vector<double> px(2);

  for (int a = 0; a < eNoNq; a++) {
    p = p + Nq(a)*yl(k,a);
    pd = pd + Nq(a)*al(k,a);
    px(0) = px(0) + Nqx(0,a)*yl(k,a);
    px(1) = px(1) + Nqx(1,a)*yl(k,a);
  }

  // Compute rho and beta depending on the volumetric penalty model
  double rho= 0;
  double beta = 0;
  double drho= 0;
  double dbeta = 0;
  mat_models::g_vol_pen(com_mod, eq.dmn[cDmn], p, rho, beta, drho, dbeta, 1.0);

  // Compute stabilization parameters
  //
  double tauM = 0.0;
  double tauC = 0.0;

  if (vmsFlag) {
    mat_models::get_tau(com_mod, eq.dmn[cDmn], Jac, Je, tauM, tauC);

  } else {
    tauM = 0.0;
    tauC = 0.0;
  }

  Array<double> NwxFi(2,eNoNw);

  for (int a = 0; a < eNoNw; a++) {
    NwxFi(0,a) = Nwx(0,a)*Fi(0,0) + Nwx(1,a)*Fi(1,0);
    NwxFi(1,a) = Nwx(0,a)*Fi(0,1) + Nwx(1,a)*Fi(1,1);
  }

  Array<double> NqxFi(2,eNoNw);

  for (int a = 0; a < eNoNq; a++) {
    NqxFi(0,a) = Nqx(0,a)*Fi(0,0) + Nqx(1,a)*Fi(1,0);
    NqxFi(1,a) = Nqx(0,a)*Fi(0,1) + Nqx(1,a)*Fi(1,1);
  }

  Array<double> VxFi(2,2);

  VxFi(0,0) = vx(0,0)*Fi(0,0) + vx(0,1)*Fi(1,0);
  VxFi(0,1) = vx(0,0)*Fi(0,1) + vx(0,1)*Fi(1,1);
  VxFi(1,0) = vx(1,0)*Fi(0,0) + vx(1,1)*Fi(1,0);
  VxFi(1,1) = vx(1,0)*Fi(0,1) + vx(1,1)*Fi(1,1);

  Vector<double> PxFi(2);
  PxFi(0) = px(0)*Fi(0,0) + px(1)*Fi(1,0);
  PxFi(1) = px(0)*Fi(0,1) + px(1)*Fi(1,1);

  double rC  = beta*pd + VxFi(0,0) + VxFi(1,1);

  Vector<double> rM(2);
  rM(0) = rho*vd(0) + PxFi(0);
  rM(1) = rho*vd(1) + PxFi(1);

  // Local residual
  //
  Vector<double> rMNqx(eNoNq);

  for (int a = 0; a < eNoNq; a++) {
    rMNqx(a) = rM(0)*NqxFi(0,a) + rM(1)*NqxFi(1,a);
    lR(2,a) = lR(2,a) + w*Jac*(Nq(a)*rC + tauM*rMNqx(a));
  }

  Vector<double> rMNwx(eNoNw);
  Array<double> VxNwx(3,eNoNw);

  for (int a = 0; a < eNoNw; a++) {
    rMNwx(a) = rM(0)*NwxFi(0,a) + rM(1)*NwxFi(1,a);
    VxNwx(0,a) = VxFi(0,0)*NwxFi(0,a) + VxFi(1,0)*NwxFi(1,a);
    VxNwx(1,a) = VxFi(0,1)*NwxFi(0,a) + VxFi(1,1)*NwxFi(1,a);
  }

  // Tangent (stiffness) matrices
  //
  double NxNx, T0, T1, T2, Ku;

  for (int b = 0; b < eNoNw; b++) {
    for (int a = 0; a < eNoNq; a++) {
      NxNx = NqxFi(0,a)*NwxFi(0,b) + NqxFi(1,a)*NwxFi(1,b);

      // dC/dV_1 + af/am *dC/dU_1 
      //
      T0 = Nq(a)*(rC*NwxFi(0,b) - VxNwx(0,b));
      T1 = tauM*(rMNqx(a)*NwxFi(0,b) - rMNwx(b)*NqxFi(0,a));
      T2 = -tauM*NxNx*PxFi(0);
      Ku = w*af*Jac*(T0 + T1 + T2);
      lKd(4,a,b) = lKd(4,a,b) + Ku;

      T1  = (am*tauM*rho)*NqxFi(0,a)*Nw(b) + af*Nq(a)*NwxFi(0,b);
      lK(6,a,b) = lK(6,a,b) + w*Jac*T1 + afm*Ku;

      // dC/dV_2 + af/am *dC/dU_2 
      T0 = Nq(a)*(rC*NwxFi(1,b) - VxNwx(1,b));
      T1 = tauM*(rMNqx(a)*NwxFi(1,b) - rMNwx(b)*NqxFi(1,a));
      T2 = -tauM*NxNx*PxFi(1);
      Ku = w*af*Jac*(T0 + T1 + T2);
      lKd(5,a,b) = lKd(5,a,b) + Ku;

      T1 = (am*tauM*rho)*NqxFi(1,a)*Nw(b) + af*Nq(a)*NwxFi(1,b);
      lK(8,a,b) = lK(8,a,b) + w*Jac*T1 + afm*Ku;
    }
  }

  for (int b = 0; b < eNoNq; b++) {
    for (int a = 0; a < eNoNq; a++) {
      // dC/dP
      NxNx = NqxFi(0,a)*NqxFi(0,b) + NqxFi(1,a)*NqxFi(1,b);
      T0 = (am*beta + af*dbeta*pd)*Nq(a)*Nq(b);
      T1 = NqxFi(0,a)*vd(0) + NqxFi(1,a)*vd(1);
      T2 = T0 + af*tauM*(NxNx + drho*T1*Nq(b));
      lK(9,a,b) = lK(9,a,b) + w*Jac*T2;
    }
  }
}

/// @brief Reproduces Fortran USTRUCT3D_C.
//
void ustruct_3d_c(ComMod& com_mod, CepMod& cep_mod, const bool vmsFlag, const int eNoNw, const int eNoNq,
    const double w, const double Je, const Vector<double>& Nw,  const Vector<double>& Nq,
    const Array<double>& Nwx, const Array<double>& Nqx, const Array<double>& al, const Array<double>& yl, 
    const Array<double>& dl, const Array<double>& bfl, Array<double>& lR, Array3<double>& lK, 
    Array3<double>& lKd)
{
  using namespace consts;
  using namespace mat_fun;

  #define n_debug_ustruct_3d_c
  #ifdef debug_ustruct_3d_c
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "w: " << w;
  dmsg << "eNoNw: " << eNoNw;
  dmsg << "eNoNq: " << eNoNq;
  #endif

  // Define parameters
  //
  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  Vector<double> fb(3);
  fb[0] = dmn.prop[PhysicalProperyType::f_x];
  fb[1] = dmn.prop[PhysicalProperyType::f_y];
  fb[2] = dmn.prop[PhysicalProperyType::f_z];

  double am = eq.am;
  double af = eq.af * eq.gam * dt;
  double afm = af / am;

  // {i,j} := velocity dofs; {k} := pressure dof
  int i = eq.s;
  int j = i + 1;
  int k = j + 1;
  int l = k + 1;

  #ifdef debug_ustruct_3d_c
  dmsg << "am: " << am;
  dmsg << "af: " << af;
  dmsg << "afm: " << afm;
  dmsg << "i: " << i;
  #endif

  // Inertia (velocity and acceleration), body force, fiber directions,
  // and deformation tensor (F) at integration point
  //
  Vector<double> vd{-fb[0], -fb[1], -fb[2]};
  Vector<double> v(3);
  Array<double> vx(3,3), F(3,3);
  F(0,0) = 1.0;
  F(1,1) = 1.0;
  F(2,2) = 1.0;

  for (int a = 0; a < eNoNw; a++) {
    v(0) = v(0) + Nw(a)*yl(i,a);
    v(1) = v(1) + Nw(a)*yl(j,a);
    v(2) = v(2) + Nw(a)*yl(k,a);

    vd(0) = vd(0) + Nw(a)*(al(i,a)-bfl(0,a));
    vd(1) = vd(1) + Nw(a)*(al(j,a)-bfl(1,a));
    vd(2) = vd(2) + Nw(a)*(al(k,a)-bfl(2,a));

    vx(0,0) = vx(0,0) + Nwx(0,a)*yl(i,a);
    vx(0,1) = vx(0,1) + Nwx(1,a)*yl(i,a);
    vx(0,2) = vx(0,2) + Nwx(2,a)*yl(i,a);

    vx(1,0) = vx(1,0) + Nwx(0,a)*yl(j,a);
    vx(1,1) = vx(1,1) + Nwx(1,a)*yl(j,a);
    vx(1,2) = vx(1,2) + Nwx(2,a)*yl(j,a);

    vx(2,0) = vx(2,0) + Nwx(0,a)*yl(k,a);
    vx(2,1) = vx(2,1) + Nwx(1,a)*yl(k,a);
    vx(2,2) = vx(2,2) + Nwx(2,a)*yl(k,a);

    F(0,0) = F(0,0) + Nwx(0,a)*dl(i,a);
    F(0,1) = F(0,1) + Nwx(1,a)*dl(i,a);
    F(0,2) = F(0,2) + Nwx(2,a)*dl(i,a);

    F(1,0) = F(1,0) + Nwx(0,a)*dl(j,a);
    F(1,1) = F(1,1) + Nwx(1,a)*dl(j,a);
    F(1,2) = F(1,2) + Nwx(2,a)*dl(j,a);

    F(2,0) = F(2,0) + Nwx(0,a)*dl(k,a);
    F(2,1) = F(2,1) + Nwx(1,a)*dl(k,a);
    F(2,2) = F(2,2) + Nwx(2,a)*dl(k,a);
  }

  double Jac = mat_fun::mat_det(F, 3);
  auto Fi = mat_fun::mat_inv(F, 3);

  // Pressure and its gradients 
  //
  double p = 0.0;
  double pd = 0.0;
  Vector<double> px(3);

  for (int a = 0; a < eNoNq; a++) {
    p = p + Nq(a)*yl(l,a);
    pd = pd + Nq(a)*al(l,a);
    px(0) = px(0) + Nqx(0,a)*yl(l,a);
    px(1) = px(1) + Nqx(1,a)*yl(l,a);
    px(2) = px(2) + Nqx(2,a)*yl(l,a);
  }

  // Compute rho and beta depending on the volumetric penalty model
  double rho= 0;
  double beta = 0;
  double drho= 0;
  double dbeta = 0;
  mat_models::g_vol_pen(com_mod, eq.dmn[cDmn], p, rho, beta, drho, dbeta, 1.0);

  // Compute stabilization parameters
  //
  double tauM = 0.0;
  double tauC = 0.0;

  if (vmsFlag) {
    mat_models::get_tau(com_mod, eq.dmn[cDmn], Jac, Je, tauM, tauC);

  } else {
    tauM = 0.0;
    tauC = 0.0;
  }

  Array<double> NwxFi(3,eNoNw);

  for (int a = 0; a < eNoNw; a++) {
    NwxFi(0,a) = Nwx(0,a)*Fi(0,0) + Nwx(1,a)*Fi(1,0) + Nwx(2,a)*Fi(2,0);
    NwxFi(1,a) = Nwx(0,a)*Fi(0,1) + Nwx(1,a)*Fi(1,1) + Nwx(2,a)*Fi(2,1);
    NwxFi(2,a) = Nwx(0,a)*Fi(0,2) + Nwx(1,a)*Fi(1,2) + Nwx(2,a)*Fi(2,2);
  }

  Array<double> NqxFi(3,eNoNw);

  for (int a = 0; a < eNoNq; a++) {
    NqxFi(0,a) = Nqx(0,a)*Fi(0,0) + Nqx(1,a)*Fi(1,0) + Nqx(2,a)*Fi(2,0);
    NqxFi(1,a) = Nqx(0,a)*Fi(0,1) + Nqx(1,a)*Fi(1,1) + Nqx(2,a)*Fi(2,1);
    NqxFi(2,a) = Nqx(0,a)*Fi(0,2) + Nqx(1,a)*Fi(1,2) + Nqx(2,a)*Fi(2,2);
  }

  Array<double> VxFi(3,3);

  VxFi(0,0) = vx(0,0)*Fi(0,0) + vx(0,1)*Fi(1,0) + vx(0,2)*Fi(2,0);
  VxFi(0,1) = vx(0,0)*Fi(0,1) + vx(0,1)*Fi(1,1) + vx(0,2)*Fi(2,1);
  VxFi(0,2) = vx(0,0)*Fi(0,2) + vx(0,1)*Fi(1,2) + vx(0,2)*Fi(2,2);

  VxFi(1,0) = vx(1,0)*Fi(0,0) + vx(1,1)*Fi(1,0) + vx(1,2)*Fi(2,0);
  VxFi(1,1) = vx(1,0)*Fi(0,1) + vx(1,1)*Fi(1,1) + vx(1,2)*Fi(2,1);
  VxFi(1,2) = vx(1,0)*Fi(0,2) + vx(1,1)*Fi(1,2) + vx(1,2)*Fi(2,2);

  VxFi(2,0) = vx(2,0)*Fi(0,0) + vx(2,1)*Fi(1,0) + vx(2,2)*Fi(2,0);
  VxFi(2,1) = vx(2,0)*Fi(0,1) + vx(2,1)*Fi(1,1) + vx(2,2)*Fi(2,1);
  VxFi(2,2) = vx(2,0)*Fi(0,2) + vx(2,1)*Fi(1,2) + vx(2,2)*Fi(2,2);

  Vector<double> PxFi(3);
  PxFi(0) = px(0)*Fi(0,0) + px(1)*Fi(1,0) + px(2)*Fi(2,0);
  PxFi(1) = px(0)*Fi(0,1) + px(1)*Fi(1,1) + px(2)*Fi(2,1);
  PxFi(2) = px(0)*Fi(0,2) + px(1)*Fi(1,2) + px(2)*Fi(2,2);

  double rC  = beta*pd + VxFi(0,0) + VxFi(1,1) + VxFi(2,2);

  Vector<double> rM(3);
  rM(0) = rho*vd(0) + PxFi(0);
  rM(1) = rho*vd(1) + PxFi(1);
  rM(2) = rho*vd(2) + PxFi(2);

  // Local residual
  //
  Vector<double> rMNqx(eNoNq);

  for (int a = 0; a < eNoNq; a++) {
    rMNqx(a) = rM(0)*NqxFi(0,a) + rM(1)*NqxFi(1,a) + rM(2)*NqxFi(2,a);
    lR(3,a) = lR(3,a) + w*Jac*(Nq(a)*rC + tauM*rMNqx(a));
  }

  Vector<double> rMNwx(eNoNw);
  Array<double> VxNwx(3,eNoNw);

  for (int a = 0; a < eNoNw; a++) {
    rMNwx(a) = rM(0)*NwxFi(0,a) + rM(1)*NwxFi(1,a) + rM(2)*NwxFi(2,a);
    VxNwx(0,a) = VxFi(0,0)*NwxFi(0,a) + VxFi(1,0)*NwxFi(1,a) + VxFi(2,0)*NwxFi(2,a);
    VxNwx(1,a) = VxFi(0,1)*NwxFi(0,a) + VxFi(1,1)*NwxFi(1,a) + VxFi(2,1)*NwxFi(2,a);
    VxNwx(2,a) = VxFi(0,2)*NwxFi(0,a) + VxFi(1,2)*NwxFi(1,a) + VxFi(2,2)*NwxFi(2,a);
  }

  // Tangent (stiffness) matrices
  //
  double NxNx, T0, T1, T2, Ku;

  for (int b = 0; b < eNoNw; b++) {
    for (int a = 0; a < eNoNq; a++) {
      NxNx = NqxFi(0,a)*NwxFi(0,b) + NqxFi(1,a)*NwxFi(1,b) + NqxFi(2,a)*NwxFi(2,b);

      // dC/dV_1 + af/am *dC/dU_1 
      //
      T0 = Nq(a)*(rC*NwxFi(0,b) - VxNwx(0,b));
      T1 = tauM*(rMNqx(a)*NwxFi(0,b) - rMNwx(b)*NqxFi(0,a));
      T2 = -tauM*NxNx*PxFi(0);
      Ku = w*af*Jac*(T0 + T1 + T2);
      lKd(9,a,b) = lKd(9,a,b) + Ku;

      T1  = (am*tauM*rho)*NqxFi(0,a)*Nw(b) + af*Nq(a)*NwxFi(0,b);
      lK(12,a,b) = lK(12,a,b) + w*Jac*T1 + afm*Ku;

      // dC/dV_2 + af/am *dC/dU_2 
      T0 = Nq(a)*(rC*NwxFi(1,b) - VxNwx(1,b));
      T1 = tauM*(rMNqx(a)*NwxFi(1,b) - rMNwx(b)*NqxFi(1,a));
      T2 = -tauM*NxNx*PxFi(1);
      Ku = w*af*Jac*(T0 + T1 + T2);
      lKd(10,a,b) = lKd(10,a,b) + Ku;

      T1 = (am*tauM*rho)*NqxFi(1,a)*Nw(b) + af*Nq(a)*NwxFi(1,b);
      lK(13,a,b) = lK(13,a,b) + w*Jac*T1 + afm*Ku;

      // dC/dV_3 + af/am *dC/dU_3 
      //
      T0 = Nq(a)*(rC*NwxFi(2,b) - VxNwx(2,b));
      T1 = tauM*(rMNqx(a)*NwxFi(2,b) - rMNwx(b)*NqxFi(2,a));
      T2 = -tauM*NxNx*PxFi(2);
      Ku = w*af*Jac*(T0 + T1 + T2);
      lKd(11,a,b) = lKd(11,a,b) + Ku;

      T1 = (am*tauM*rho)*NqxFi(2,a)*Nw(b) + af*Nq(a)*NwxFi(2,b);
      lK(14,a,b) = lK(14,a,b) + w*Jac*T1 + afm*Ku;
    }
  }

  for (int b = 0; b < eNoNq; b++) {
    for (int a = 0; a < eNoNq; a++) {
      // dC/dP
      NxNx = NqxFi(0,a)*NqxFi(0,b) + NqxFi(1,a)*NqxFi(1,b) + NqxFi(2,a)*NqxFi(2,b);
      T0 = (am*beta + af*dbeta*pd)*Nq(a)*Nq(b);
      T1 = NqxFi(0,a)*vd(0) + NqxFi(1,a)*vd(1) + NqxFi(2,a)*vd(2);
      T2 = T0 + af*tauM*(NxNx + drho*T1*Nq(b));
      lK(15,a,b) = lK(15,a,b) + w*Jac*T2;
    }
  }
}

/// @brief Replicates Fortran USTRUCT2D_M.
//
void ustruct_2d_m(ComMod& com_mod, CepMod& cep_mod, const bool vmsFlag, const int eNoNw, const int eNoNq, 
    const int nFn, const double w, const double Je, const Vector<double>& Nw,  const Vector<double>& Nq, 
    const Array<double>& Nwx, const Array<double>& al, const Array<double>& yl, const Array<double>& dl, 
    const Array<double>& bfl, const Array<double>& fN, const Vector<double>& ya_l, Array<double>& lR, 
    Array3<double>& lK, Array3<double>& lKd)
{
  using namespace consts;
  using namespace mat_fun;

  #define n_debug_ustruct_2d_m
  #ifdef debug_ustruct_2d_m
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

  // Define parameters
  //
  Vector<double> fb(2);
  fb[0] = dmn.prop[PhysicalProperyType::f_x];
  fb[1] = dmn.prop[PhysicalProperyType::f_y];

  double am = eq.am;
  double af = eq.af * eq.gam * dt;
  double afm = af / am;

  // {i,j} := velocity dofs; {k} := pressure dof
  int i = eq.s;
  int j = i + 1;
  int k = j + 1;

  #ifdef debug_ustruct_2d_m
  dmsg << "am: " << am;
  dmsg << "af: " << af;
  dmsg << "afm: " << afm;
  dmsg << "i: " << i;
  #endif

  // Inertia (velocity and acceleration), body force, fiber directions,
  // and deformation tensor (F) at integration point
  //
  Vector<double> vd{-fb[0], -fb[1]};
  Vector<double> v(2);
  Array<double> vx(2,2), F(2,2);
  double ya_g = 0.0;
  F(0,0) = 1.0;
  F(1,1) = 1.0;

  for (int a = 0; a < eNoNw; a++) {
    v(0) = v(0) + Nw(a)*yl(i,a);
    v(1) = v(1) + Nw(a)*yl(j,a);

    vd(0) = vd(0) + Nw(a)*(al(i,a)-bfl(0,a));
    vd(1) = vd(1) + Nw(a)*(al(j,a)-bfl(1,a));

    vx(0,0) = vx(0,0) + Nwx(0,a)*yl(i,a);
    vx(0,1) = vx(0,1) + Nwx(1,a)*yl(i,a);
    vx(1,0) = vx(1,0) + Nwx(0,a)*yl(j,a);
    vx(1,1) = vx(1,1) + Nwx(1,a)*yl(j,a);

    F(0,0) = F(0,0) + Nwx(0,a)*dl(i,a);
    F(0,1) = F(0,1) + Nwx(1,a)*dl(i,a);
    F(1,0) = F(1,0) + Nwx(0,a)*dl(j,a);
    F(1,1) = F(1,1) + Nwx(1,a)*dl(j,a);

    ya_g = ya_g + Nw(a)*ya_l(a);
  }

  double Jac = mat_fun::mat_det(F, 2);
  auto Fi = mat_fun::mat_inv(F, 2);

  // Pressure and its time derivative
  //
  double p = 0.0;
  double pd = 0.0;

  for (int a = 0; a < eNoNq; a++) {
    p = p + Nq(a)*yl(k,a);
    pd = pd + Nq(a)*al(k,a);
  }

  // Compute deviatoric 2nd Piola-Kirchhoff stress tensor (Siso) and
  // isochoric elasticity tensor in Voigt notation (Dm)
  Array<double> Siso(2,2), Dm(3,3);
  double Ja = 0;
  mat_models::get_pk2cc(com_mod, cep_mod, eq.dmn[cDmn], F, nFn, fN, ya_g, Siso, Dm, Ja);

   // Viscous 2nd Piola-Kirchhoff stress and tangent contributions
  Array<double> Svis(2,2);
  Array3<double> Kvis_u(4, eNoNw, eNoNw);
  Array3<double> Kvis_v(4, eNoNw, eNoNw);
  
  mat_models::get_visc_stress_and_tangent(dmn, eNoNw, Nwx, vx, F, Svis, Kvis_u, Kvis_v);

  // Compute rho and beta depending on the volumetric penalty model
  //
  double rho= 0;
  double beta = 0;
  double drho= 0;
  double dbeta = 0;
  mat_models::g_vol_pen(com_mod, eq.dmn[cDmn], p, rho, beta, drho, dbeta, Ja);

  // Compute stabilization parameters
  double tauM = 0.0;
  double tauC = 0.0;

  if (vmsFlag) {
    mat_models::get_tau(com_mod, eq.dmn[cDmn], Jac, Je, tauM, tauC);
  } else {
    tauM = 0.0;
    tauC = 0.0;
  }

  // Total isochoric 2nd Piola-Kirchhoff stress (Elastic + Viscous)
  Siso = Siso + Svis;

  // Deviatoric 1st Piola-Kirchhoff tensor (P)
  //
  auto Pdev = mat_fun::mat_mul(F, Siso);


  // Shape function gradients in the current configuration
  //
  Array<double> NxFi(2,eNoNw);

  for (int a = 0; a < eNoNw; a++) {
    NxFi(0,a) = Nwx(0,a)*Fi(0,0) + Nwx(1,a)*Fi(1,0);
    NxFi(1,a) = Nwx(0,a)*Fi(0,1) + Nwx(1,a)*Fi(1,1);
  }

   // Velocity gradient in current configuration
  auto VxFi = mat_mul(vx, Fi);
  double rC  = beta*pd + VxFi(1,1) + VxFi(2,2);
  double rCl = -p + tauC*rC;

  // Local residual
  //
  for (int a = 0; a < eNoNw; a++) {
    double T1 = Jac*rho*vd(0)*Nw(a);
    double T2 = Pdev(0,0)*Nwx(0,a) + Pdev(0,1)*Nwx(1,a);
    double T3 = Jac*rCl*NxFi(0,a);
    lR(0,a) = lR(0,a) + w*(T1 + T2 + T3);

    T1 = Jac*rho*vd(1)*Nw(a);
    T2 = Pdev(1,0)*Nwx(0,a) + Pdev(1,1)*Nwx(1,a);
    T3 = Jac*rCl*NxFi(1,a);
    lR(1,a) = lR(1,a) + w*(T1 + T2 + T3);
  }

  // Auxilary quantities for computing stiffness tensors
  //
  Array3<double> Bm(3,2,eNoNw);

  for (int a = 0; a < eNoNw; a++) {
    Bm(0,0,a) = Nwx(0,a)*F(0,0);
    Bm(0,1,a) = Nwx(0,a)*F(1,0);

    Bm(1,0,a) = Nwx(1,a)*F(0,1);
    Bm(1,1,a) = Nwx(1,a)*F(1,1);

    Bm(2,0,a) = Nwx(2,a)*F(0,2) + F(0,0)*Nwx(1,a);
    Bm(2,1,a) = Nwx(2,a)*F(1,2) + F(1,0)*Nwx(1,a);
  }

  Array<double> VxNx(2,eNoNw);

  for (int a = 0; a < eNoNw; a++) {
    VxNx(0,a) = VxFi(0,0)*NxFi(0,a) + VxFi(1,0)*NxFi(1,a);
    VxNx(1,a) = VxFi(0,1)*NxFi(0,a) + VxFi(1,1)*NxFi(1,a);
  }

  // Tangent (stiffness) matrices
  //
  double NxSNx{0.0}, BtDB{0.0}, 
      T1{0.0}, T2{0.0}, T3{0.0},
      Tv{0.0}, Ku{0.0};

  Array<double> DBm(3,2);

  for (int b = 0; b < eNoNw; b++) {
    for (int a = 0; a < eNoNw; a++) {
      NxSNx = Nwx(0,a)*Siso(0,0)*Nwx(0,b)
            + Nwx(0,a)*Siso(0,1)*Nwx(1,b)
            + Nwx(1,a)*Siso(1,0)*Nwx(0,b)
            + Nwx(1,a)*Siso(1,1)*Nwx(1,b);

      DBm(0,0) = Dm(0,0)*Bm(0,0,b) + Dm(0,1)*Bm(1,0,b) + Dm(0,2)*Bm(2,0,b);
      DBm(0,1) = Dm(0,0)*Bm(0,1,b) + Dm(0,1)*Bm(1,1,b) + Dm(0,2)*Bm(2,1,b);

      DBm(1,0) = Dm(1,0)*Bm(0,0,b) + Dm(1,1)*Bm(1,0,b) + Dm(1,2)*Bm(2,0,b);
      DBm(1,1) = Dm(1,0)*Bm(0,1,b) + Dm(1,1)*Bm(1,1,b) + Dm(1,2)*Bm(2,1,b);

      DBm(2,0) = Dm(2,0)*Bm(0,0,b) + Dm(2,1)*Bm(1,0,b) + Dm(2,2)*Bm(2,0,b);
      DBm(2,1) = Dm(2,0)*Bm(0,1,b) + Dm(2,1)*Bm(1,1,b) + Dm(2,2)*Bm(2,1,b);


      // dM1_dV1 + af/am *dM_1/dU_1
      //
      BtDB = Bm(0,0,a)*DBm(0,0) + Bm(1,0,a)*DBm(1,0) + Bm(2,0,a)*DBm(2,0);
      T1   = Jac*rho*vd(0)*Nw(a)*NxFi(0,b);
      T2   = -tauC*Jac*NxFi(0,a)*VxNx(0,b);

      Ku   = w*af*(T1 + T2 + BtDB + NxSNx + Kvis_u(0,a,b));
      lKd(0,a,b) = lKd(0,a,b) + Ku;

      T1   = am*Jac*rho*Nw(a)*Nw(b);
      T2   = T1 + af*Jac*tauC*rho*NxFi(0,a)*NxFi(0,b);
      Tv   = af*Kvis_v(0,a,b);
      lK(0,a,b)  = lK(0,a,b) + w*(T2 + Tv) + afm*Ku;

      // dM_1/dV_2 + af/am *dM_1/dU_2
      //
      BtDB = Bm(0,0,a)*DBm(0,1) + Bm(1,0,a)*DBm(1,1) + Bm(2,0,a)*DBm(2,1);
      T1   = Jac*rho*vd(0)*Nw(a)*NxFi(1,b);
      T2   = -tauC*Jac*NxFi(0,a)*VxNx(1,b);
      T3   = Jac*rCl*(NxFi(0,a)*NxFi(1,b) - NxFi(1,a)*NxFi(0,b));

      Ku   = w*af*(T1 + T2 + T3 + BtDB + Kvis_u(1,a,b));
      lKd(1,a,b) = lKd(1,a,b) + Ku;

      T2   = af*Jac*tauC*rho*NxFi(0,a)*NxFi(1,b);
      Tv   = af*Kvis_v(1,a,b);
      lK(1,a,b) = lK(1,a,b) + w*(T2 + Tv) + afm*Ku;

      // dM_2/dV_1 + af/am *dM_2/dU_1
      //
      BtDB = Bm(0,1,a)*DBm(0,0) + Bm(1,1,a)*DBm(1,0) + Bm(2,1,a)*DBm(2,0);
      T1   = Jac*rho*vd(1)*Nw(a)*NxFi(0,b);
      T2   = -tauC*Jac*NxFi(1,a)*VxNx(0,b);
      T3   = Jac*rCl*(NxFi(1,a)*NxFi(0,b) - NxFi(0,a)*NxFi(1,b));

      Ku   = w*af*(T1 + T2 + T3 + BtDB + Kvis_u(2,a,b));
      lKd(2,a,b) = lKd(2,a,b) + Ku;

      T2   = af*Jac*tauC*rho*NxFi(1,a)*NxFi(0,b);
      Tv   = af*Kvis_v(2,a,b);
      lK(3,a,b) = lK(3,a,b) + w*(T2 + Tv) + afm*Ku;

      // dM_2/dV_2 + af/am *dM_2/dU_2
      //
      BtDB = Bm(0,1,a)*DBm(0,1) + Bm(1,1,a)*DBm(1,1) + Bm(2,1,a)*DBm(2,1);
      T1   = Jac*rho*vd(1)*Nw(a)*NxFi(1,b);
      T2   = -tauC*Jac*NxFi(1,a)*VxNx(1,b);

      Ku   = w*af*(T1 + T2 + BtDB + NxSNx + Kvis_u(3,a,b));
      lKd(3,a,b) = lKd(3,a,b) + Ku;

      T1   = am*Jac*rho*Nw(a)*Nw(b);
      T2   = T1 + af*Jac*tauC*rho*NxFi(1,a)*NxFi(1,b);
      Tv   = af*Kvis_v(3,a,b);
      lK(4,a,b) = lK(4,a,b) + w*(T2 + Tv) + afm*Ku;
    }
  }

  for (int b = 0; b < eNoNq; b++) {
    for (int a = 0; a < eNoNw; a++) {
      double T0, T1;

      // dM_0/dP
      T0 = am*tauC*beta + af*(tauC*dbeta*pd - 1.0);
      T1 = T0*NxFi(0,a)*Nq(b) + af*drho*vd(0)*Nw(a)*Nq(b);
      lK(2,a,b) = lK(2,a,b) + w*Jac*T1;

      // dM_1/dP
      T1 = T0*NxFi(1,a)*Nq(b) + af*drho*vd(1)*Nw(a)*Nq(b);
      lK(6,a,b) = lK(6,a,b) + w*Jac*T1;
    }
  }
}

/// @brief Reproduces Fortran USTRUCT3D_M.
//
void ustruct_3d_m(ComMod& com_mod, CepMod& cep_mod, const bool vmsFlag, const int eNoNw, const int eNoNq, 
    const int nFn, const double w, const double Je, const Vector<double>& Nw,  const Vector<double>& Nq, 
    const Array<double>& Nwx, const Array<double>& al, const Array<double>& yl, const Array<double>& dl, 
    const Array<double>& bfl, const Array<double>& fN, const Vector<double>& ya_l, Array<double>& lR, 
    Array3<double>& lK, Array3<double>& lKd)
{
  using namespace consts;
  using namespace mat_fun;

  #define n_debug_ustruct_3d_m
  #ifdef debug_ustruct_3d_m
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

  // Define parameters

  Vector<double> fb(3);
  fb[0] = dmn.prop[PhysicalProperyType::f_x];
  fb[1] = dmn.prop[PhysicalProperyType::f_y];
  fb[2] = dmn.prop[PhysicalProperyType::f_z];

  double am = eq.am;
  double af = eq.af * eq.gam * dt;
  double afm = af / am;

  // {i,j} := velocity dofs; {k} := pressure dof
  int i = eq.s;
  int j = i + 1;
  int k = j + 1;
  int l = k + 1;

  #ifdef debug_ustruct_3d_m
  dmsg << "fb: " << fb;
  dmsg << "am: " << am;
  dmsg << "af: " << af;
  dmsg << "afm: " << afm;
  dmsg << "i: " << i;
  #endif

  // Inertia (velocity and acceleration), body force, fiber directions,
  // and deformation tensor (F) at integration point
  //
  Vector<double> vd{-fb[0], -fb[1], -fb[2]};
  Vector<double> v(3);
  Array<double> vx(3,3), F(3,3);
  double ya_g = 0.0;
  F(0,0) = 1.0;
  F(1,1) = 1.0;
  F(2,2) = 1.0;

  for (int a = 0; a < eNoNw; a++) {
    v(0) = v(0) + Nw(a)*yl(i,a);
    v(1) = v(1) + Nw(a)*yl(j,a);
    v(2) = v(2) + Nw(a)*yl(k,a);

    vd(0) = vd(0) + Nw(a)*(al(i,a)-bfl(0,a));
    vd(1) = vd(1) + Nw(a)*(al(j,a)-bfl(1,a));
    vd(2) = vd(2) + Nw(a)*(al(k,a)-bfl(2,a));

    vx(0,0) = vx(0,0) + Nwx(0,a)*yl(i,a);
    vx(0,1) = vx(0,1) + Nwx(1,a)*yl(i,a);
    vx(0,2) = vx(0,2) + Nwx(2,a)*yl(i,a);

    vx(1,0) = vx(1,0) + Nwx(0,a)*yl(j,a);
    vx(1,1) = vx(1,1) + Nwx(1,a)*yl(j,a);
    vx(1,2) = vx(1,2) + Nwx(2,a)*yl(j,a);

    vx(2,0) = vx(2,0) + Nwx(0,a)*yl(k,a);
    vx(2,1) = vx(2,1) + Nwx(1,a)*yl(k,a);
    vx(2,2) = vx(2,2) + Nwx(2,a)*yl(k,a);

    F(0,0) = F(0,0) + Nwx(0,a)*dl(i,a);
    F(0,1) = F(0,1) + Nwx(1,a)*dl(i,a);
    F(0,2) = F(0,2) + Nwx(2,a)*dl(i,a);

    F(1,0) = F(1,0) + Nwx(0,a)*dl(j,a);
    F(1,1) = F(1,1) + Nwx(1,a)*dl(j,a);
    F(1,2) = F(1,2) + Nwx(2,a)*dl(j,a);

    F(2,0) = F(2,0) + Nwx(0,a)*dl(k,a);
    F(2,1) = F(2,1) + Nwx(1,a)*dl(k,a);
    F(2,2) = F(2,2) + Nwx(2,a)*dl(k,a);

    ya_g = ya_g + Nw(a)*ya_l(a);
  }

  double Jac = mat_fun::mat_det(F, 3);
  auto Fi = mat_fun::mat_inv(F, 3);

  // Pressure and its time derivative
  //
  double p = 0.0;
  double pd = 0.0;

  for (int a = 0; a < eNoNq; a++) {
    p = p + Nq(a)*yl(l,a);
    pd = pd + Nq(a)*al(l,a);
  }

  // Compute deviatoric 2nd Piola-Kirchhoff stress tensor (Siso) and
  // isochoric elasticity tensor in Voigt notation (Dm)
  //
  Array<double> Siso(3,3), Dm(6,6);
  double Ja = 0;
  mat_models::get_pk2cc(com_mod, cep_mod, eq.dmn[cDmn], F, nFn, fN, ya_g, Siso, Dm, Ja);

  // Viscous 2nd Piola-Kirchhoff stress and tangent contributions
  Array<double> Svis(3,3);
  Array3<double> Kvis_u(9, eNoNw, eNoNw);
  Array3<double> Kvis_v(9, eNoNw, eNoNw);
  
  mat_models::get_visc_stress_and_tangent(dmn, eNoNw, Nwx, vx, F, Svis, Kvis_u, Kvis_v);


  // Compute rho and beta depending on the volumetric penalty model
  //
  double rho= 0;
  double beta = 0;
  double drho= 0;
  double dbeta = 0;
  mat_models::g_vol_pen(com_mod, eq.dmn[cDmn], p, rho, beta, drho, dbeta, Ja);

  // Compute stabilization parameters
  double tauM = 0.0;
  double tauC = 0.0;

  if (vmsFlag) {
    mat_models::get_tau(com_mod, eq.dmn[cDmn], Jac, Je, tauM, tauC);
  } else { 
    tauM = 0.0;
    tauC = 0.0;
  }

  // Total isochoric 2nd Piola-Kirchhoff stress (Elastic + Viscous)
  Siso = Siso + Svis;

  // Deviatoric 1st Piola-Kirchhoff tensor (P)
  //
  auto Pdev = mat_fun::mat_mul(F, Siso);

  // Shape function gradients in the current configuration
  //
  Array<double> NxFi(3,eNoNw);

  for (int a = 0; a < eNoNw; a++) {
    NxFi(0,a) = Nwx(0,a)*Fi(0,0) + Nwx(1,a)*Fi(1,0) + Nwx(2,a)*Fi(2,0);
    NxFi(1,a) = Nwx(0,a)*Fi(0,1) + Nwx(1,a)*Fi(1,1) + Nwx(2,a)*Fi(2,1);
    NxFi(2,a) = Nwx(0,a)*Fi(0,2) + Nwx(1,a)*Fi(1,2) + Nwx(2,a)*Fi(2,2);
  } 

  // Velocity gradient in current configuration
  auto VxFi = mat_mul(vx, Fi);
  double rC  = beta*pd + VxFi(0,0) + VxFi(1,1) + VxFi(2,2);
  double rCl = -p + tauC*rC;

  // Local residual
  //
  double T1, T2, T3;

  for (int a = 0; a < eNoNw; a++) {
    T1 = Jac*rho*vd(0)*Nw(a);
    T2 = Pdev(0,0)*Nwx(0,a) + Pdev(0,1)*Nwx(1,a) + Pdev(0,2)*Nwx(2,a);
    T3 = Jac*rCl*NxFi(0,a);
    lR(0,a) = lR(0,a) + w*(T1 + T2 + T3);

    T1 = Jac*rho*vd(1)*Nw(a);
    T2 = Pdev(1,0)*Nwx(0,a) + Pdev(1,1)*Nwx(1,a) + Pdev(1,2)*Nwx(2,a);
    T3 = Jac*rCl*NxFi(1,a);
    lR(1,a) = lR(1,a) + w*(T1 + T2 + T3);

    T1 = Jac*rho*vd(2)*Nw(a);
    T2 = Pdev(2,0)*Nwx(0,a) + Pdev(2,1)*Nwx(1,a) + Pdev(2,2)*Nwx(2,a);
    T3 = Jac*rCl*NxFi(2,a);
    lR(2,a) = lR(2,a) + w*(T1 + T2 + T3);
  }

  // Auxilary quantities for computing stiffness tensors
  //
  Array3<double> Bm(6,3,eNoNw);

  for (int a = 0; a < eNoNw; a++) {
    Bm(0,0,a) = Nwx(0,a)*F(0,0);
    Bm(0,1,a) = Nwx(0,a)*F(1,0);
    Bm(0,2,a) = Nwx(0,a)*F(2,0);

    Bm(1,0,a) = Nwx(1,a)*F(0,1);
    Bm(1,1,a) = Nwx(1,a)*F(1,1);
    Bm(1,2,a) = Nwx(1,a)*F(2,1);

    Bm(2,0,a) = Nwx(2,a)*F(0,2);
    Bm(2,1,a) = Nwx(2,a)*F(1,2);
    Bm(2,2,a) = Nwx(2,a)*F(2,2);

    Bm(3,0,a) = (Nwx(0,a)*F(0,1) + F(0,0)*Nwx(1,a));
    Bm(3,1,a) = (Nwx(0,a)*F(1,1) + F(1,0)*Nwx(1,a));
    Bm(3,2,a) = (Nwx(0,a)*F(2,1) + F(2,0)*Nwx(1,a));

    Bm(4,0,a) = (Nwx(1,a)*F(0,2) + F(0,1)*Nwx(2,a));
    Bm(4,1,a) = (Nwx(1,a)*F(1,2) + F(1,1)*Nwx(2,a));
    Bm(4,2,a) = (Nwx(1,a)*F(2,2) + F(2,1)*Nwx(2,a));

    Bm(5,0,a) = (Nwx(2,a)*F(0,0) + F(0,2)*Nwx(0,a));
    Bm(5,1,a) = (Nwx(2,a)*F(1,0) + F(1,2)*Nwx(0,a));
    Bm(5,2,a) = (Nwx(2,a)*F(2,0) + F(2,2)*Nwx(0,a));
  }

  Array<double> VxNx(3,eNoNw);

  for (int a = 0; a < eNoNw; a++) {
    VxNx(0,a) = VxFi(0,0)*NxFi(0,a) + VxFi(1,0)*NxFi(1,a) + VxFi(2,0)*NxFi(2,a);
    VxNx(1,a) = VxFi(0,1)*NxFi(0,a) + VxFi(1,1)*NxFi(1,a) + VxFi(2,1)*NxFi(2,a);
    VxNx(2,a) = VxFi(0,2)*NxFi(0,a) + VxFi(1,2)*NxFi(1,a) + VxFi(2,2)*NxFi(2,a);
  }

  // Tangent (stiffness) matrices
  //
  double r13 = 1.0 / 3.0;
  double r23 = 2.0 / 3.0;
  double NxSNx{0.0}, BtDB{0.0};
  double Tv{0.0}, Ku{0.0};

  for (int b = 0; b < eNoNw; b++) {
    for (int a = 0; a < eNoNw; a++) {
      NxSNx = Nwx(0,a)*Siso(0,0)*Nwx(0,b)
       + Nwx(0,a)*Siso(0,1)*Nwx(1,b) + Nwx(0,a)*Siso(0,2)*Nwx(2,b)
       + Nwx(1,a)*Siso(1,0)*Nwx(0,b) + Nwx(1,a)*Siso(1,1)*Nwx(1,b)
       + Nwx(1,a)*Siso(1,2)*Nwx(2,b) + Nwx(2,a)*Siso(2,0)*Nwx(0,b)
       + Nwx(2,a)*Siso(2,1)*Nwx(1,b) + Nwx(2,a)*Siso(2,2)*Nwx(2,b);

      auto DBm = mat_mul(Dm, Bm.rslice(b));

      // dM1_dV1 + af/am *dM_1/dU_1
      BtDB = Bm(0,0,a)*DBm(0,0) + Bm(1,0,a)*DBm(1,0) +
             Bm(2,0,a)*DBm(2,0) + Bm(3,0,a)*DBm(3,0) +
             Bm(4,0,a)*DBm(4,0) + Bm(5,0,a)*DBm(5,0);
      T1   = Jac*rho*vd(0)*Nw(a)*NxFi(0,b);
      T2   = -tauC*Jac*NxFi(0,a)*VxNx(0,b);
 
      Ku   = w*af*(T1 + T2 + BtDB + NxSNx + Kvis_u(0,a,b));
      lKd(0,a,b) = lKd(0,a,b) + Ku;
 
      T1   = am*Jac*rho*Nw(a)*Nw(b);
      T2   = T1 + af*Jac*tauC*rho*NxFi(0,a)*NxFi(0,b);
      Tv   = af*Kvis_v(0,a,b);
      lK(0,a,b)  = lK(0,a,b) + w*(T2 + Tv) + afm*Ku;

      // dM_1/dV_2 + af/am *dM_1/dU_2
      BtDB = Bm(0,0,a)*DBm(0,1) + Bm(1,0,a)*DBm(1,1) +
             Bm(2,0,a)*DBm(2,1) + Bm(3,0,a)*DBm(3,1) +
             Bm(4,0,a)*DBm(4,1) + Bm(5,0,a)*DBm(5,1);
      T1   = Jac*rho*vd(0)*Nw(a)*NxFi(1,b);
      T2   = -tauC*Jac*NxFi(0,a)*VxNx(1,b);
      T3   = Jac*rCl*(NxFi(0,a)*NxFi(1,b) - NxFi(1,a)*NxFi(0,b));
 
      Ku   = w*af*(T1 + T2 + T3 + BtDB + Kvis_u(1,a,b));
      lKd(1,a,b) = lKd(1,a,b) + Ku;
 
      T2   = af*Jac*tauC*rho*NxFi(0,a)*NxFi(1,b);
      Tv   = af*Kvis_v(1,a,b);
      lK(1,a,b) = lK(1,a,b) + w*(T2 + Tv) + afm*Ku;

      // dM_1/dV_3 + af/am *dM_1/dU_3
      //
      BtDB = Bm(0,0,a)*DBm(0,2) + Bm(1,0,a)*DBm(1,2) +
             Bm(2,0,a)*DBm(2,2) + Bm(3,0,a)*DBm(3,2) +
             Bm(4,0,a)*DBm(4,2) + Bm(5,0,a)*DBm(5,2);
      T1   = Jac*rho*vd(0)*Nw(a)*NxFi(2,b);
      T2   = -tauC*Jac*NxFi(0,a)*VxNx(2,b);
      T3   = Jac*rCl*(NxFi(0,a)*NxFi(2,b) - NxFi(2,a)*NxFi(0,b));
 
      Ku   = w*af*(T1 + T2 + T3 + BtDB + Kvis_u(2,a,b));
      lKd(2,a,b) = lKd(2,a,b) + Ku;
 
      T2   = af*Jac*tauC*rho*NxFi(0,a)*NxFi(2,b);
      Tv   = af*Kvis_v(2,a,b);
      lK(2,a,b) = lK(2,a,b) + w*(T2 + Tv) + afm*Ku;

      // dM_2/dV_1 + af/am *dM_2/dU_1
      //
      BtDB = Bm(0,1,a)*DBm(0,0) + Bm(1,1,a)*DBm(1,0) +
             Bm(2,1,a)*DBm(2,0) + Bm(3,1,a)*DBm(3,0) +
             Bm(4,1,a)*DBm(4,0) + Bm(5,1,a)*DBm(5,0);

      T1   = Jac*rho*vd(1)*Nw(a)*NxFi(0,b);
      T2   = -tauC*Jac*NxFi(1,a)*VxNx(0,b);
      T3   = Jac*rCl*(NxFi(1,a)*NxFi(0,b) - NxFi(0,a)*NxFi(1,b));
 
      Ku   = w*af*(T1 + T2 + T3 + BtDB + Kvis_u(3,a,b));
      lKd(3,a,b) = lKd(3,a,b) + Ku;
 
      T2   = af*Jac*tauC*rho*NxFi(1,a)*NxFi(0,b);
      Tv   = af*Kvis_v(3,a,b);

      lK(4,a,b) = lK(4,a,b) + w*(T2 + Tv) + afm*Ku;

      // dM_2/dV_2 + af/am *dM_2/dU_2
      //
      BtDB = Bm(0,1,a)*DBm(0,1) + Bm(1,1,a)*DBm(1,1) +
             Bm(2,1,a)*DBm(2,1) + Bm(3,1,a)*DBm(3,1) +
             Bm(4,1,a)*DBm(4,1) + Bm(5,1,a)*DBm(5,1);

      T1   = Jac*rho*vd(1)*Nw(a)*NxFi(1,b);

      T2   = -tauC*Jac*NxFi(1,a)*VxNx(1,b);

 
      Ku   = w*af*(T1 + T2 + BtDB + NxSNx + Kvis_u(4,a,b));
      lKd(4,a,b) = lKd(4,a,b) + Ku;
 
      T1   = am*Jac*rho*Nw(a)*Nw(b);
      T2   = T1 + af*Jac*tauC*rho*NxFi(1,a)*NxFi(1,b);
      Tv   = af*Kvis_v(4,a,b);
      lK(5,a,b) = lK(5,a,b) + w*(T2 + Tv) + afm*Ku;


      // dM_2/dV_3 + af/am *dM_2/dU_3
      //
      BtDB = Bm(0,1,a)*DBm(0,2) + Bm(1,1,a)*DBm(1,2) +
             Bm(2,1,a)*DBm(2,2) + Bm(3,1,a)*DBm(3,2) +
             Bm(4,1,a)*DBm(4,2) + Bm(5,1,a)*DBm(5,2);

      T1   = Jac*rho*vd(1)*Nw(a)*NxFi(2,b);
      T2   = -tauC*Jac*NxFi(1,a)*VxNx(2,b);
      T3   = Jac*rCl*(NxFi(1,a)*NxFi(2,b) - NxFi(2,a)*NxFi(1,b));

 
      Ku   = w*af*(T1 + T2 + T3 + BtDB + Kvis_u(5,a,b));
      lKd(5,a,b) = lKd(5,a,b) + Ku;
 
      T2   = af*Jac*tauC*rho*NxFi(1,a)*NxFi(2,b);
      Tv   = af*Kvis_v(5,a,b);
      lK(6,a,b) = lK(6,a,b) + w*(T2 + Tv) + afm*Ku;

      // dM_3/dV_1 + af/am *dM_3/dU_1
      //
      BtDB = Bm(0,2,a)*DBm(0,0) + Bm(1,2,a)*DBm(1,0) +
             Bm(2,2,a)*DBm(2,0) + Bm(3,2,a)*DBm(3,0) +
             Bm(4,2,a)*DBm(4,0) + Bm(5,2,a)*DBm(5,0);

      T1   = Jac*rho*vd(2)*Nw(a)*NxFi(0,b);
      T2   = -tauC*Jac*NxFi(2,a)*VxNx(0,b);
      T3   = Jac*rCl*(NxFi(2,a)*NxFi(0,b) - NxFi(0,a)*NxFi(2,b));
 
      Ku   = w*af*(T1 + T2 + T3 + BtDB + Kvis_u(6,a,b));
      lKd(6,a,b) = lKd(6,a,b) + Ku;
 
      T2   = af*Jac*tauC*rho*NxFi(2,a)*NxFi(0,b);
      Tv   = af*Kvis_v(6,a,b);
      lK(8,a,b) = lK(8,a,b) + w*(T2 + Tv) + afm*Ku;

      // dM_3/dV_2 + af/am *dM_3/dU_2
      //
      BtDB = Bm(0,2,a)*DBm(0,1) + Bm(1,2,a)*DBm(1,1) +
             Bm(2,2,a)*DBm(2,1) + Bm(3,2,a)*DBm(3,1) +
             Bm(4,2,a)*DBm(4,1) + Bm(5,2,a)*DBm(5,1);

      T1   = Jac*rho*vd(2)*Nw(a)*NxFi(1,b);
      T2   = -tauC*Jac*NxFi(2,a)*VxNx(1,b);
      T3   = Jac*rCl*(NxFi(2,a)*NxFi(1,b) - NxFi(1,a)*NxFi(2,b));
 
      Ku   = w*af*(T1 + T2 + T3 + BtDB + Kvis_u(7,a,b));
      lKd(7,a,b) = lKd(7,a,b) + Ku;
 
      T2   = af*Jac*tauC*rho*NxFi(2,a)*NxFi(1,b);
      Tv   = af*Kvis_v(7,a,b);

      lK(9,a,b) = lK(9,a,b) + w*(T2 + Tv) + afm*Ku;

      // dM_3/dV_3 + af/am *dM_3/dU_3
      //
      BtDB = Bm(0,2,a)*DBm(0,2) + Bm(1,2,a)*DBm(1,2) +
             Bm(2,2,a)*DBm(2,2) + Bm(3,2,a)*DBm(3,2) +
             Bm(4,2,a)*DBm(4,2) + Bm(5,2,a)*DBm(5,2);

      T1   = Jac*rho*vd(2)*Nw(a)*NxFi(2,b);
      T2   = -tauC*Jac*NxFi(2,a)*VxNx(2,b);
 
      Ku   = w*af*(T1 + T2 + BtDB + NxSNx + Kvis_u(8,a,b));
      lKd(8,a,b) = lKd(8,a,b) + Ku;
 
      T1   = am*Jac*rho*Nw(a)*Nw(b);
      T2   = T1 + af*Jac*tauC*rho*NxFi(2,a)*NxFi(2,b);
      Tv   = af*Kvis_v(8,a,b);

      lK(10,a,b) = lK(10,a,b) + w*(T2 + Tv) + afm*Ku;
    }
  }

  for (int b = 0; b < eNoNq; b++) {
    for (int a = 0; a < eNoNw; a++) {
      double T0, T1;

      // dM_0/dP
      T0 = am*tauC*beta + af*(tauC*dbeta*pd - 1.0);
      T1 = T0*NxFi(0,a)*Nq(b) + af*drho*vd(0)*Nw(a)*Nq(b);
      lK(3,a,b) = lK(3,a,b) + w*Jac*T1;

      // dM_1/dP
      T1 = T0*NxFi(1,a)*Nq(b) + af*drho*vd(1)*Nw(a)*Nq(b);
      lK(7,a,b) = lK(7,a,b) + w*Jac*T1;

      // dM_2/dP
      T1 = T0*NxFi(2,a)*Nq(b) + af*drho*vd(2)*Nw(a)*Nq(b);
      lK(11,a,b) = lK(11,a,b) + w*Jac*T1;
    }
  }

}

/// @brief Replicates 'SUBROUTINE USTRUCT_DOASSEM(d, eqN, lKd, lK, lR)'
//
void ustruct_do_assem(ComMod& com_mod, const int d, const Vector<int>& eqN, const Array3<double>& lKd, 
    const Array3<double>& lK, const Array<double>& lR)
{
  const int nsd = com_mod.nsd;
  const auto& idMap = com_mod.idMap;
  auto& R = com_mod.R;
  auto& Kd = com_mod.Kd;
  auto& Val = com_mod.Val;

  for (int a = 0; a < d; a++) {
    // Momentum equation residual is assembled at mapped rows
    int rowN = idMap(eqN(a));
    for (int i = 0; i < nsd; i++) {
      R(i,rowN) = R(i,rowN) + lR(i,a);
    }

    // Continuity equation residual is assembled at unmapped rows
    rowN = eqN(a);
    R(nsd,rowN) = R(nsd,rowN) + lR(nsd,a);
  }

  if (nsd == 3) {
    // Stiffness matrix (A) is assembled using mapped rows and columns
    // Gradient matrix (B) is assembled using mapped row but unmapped
    // column
    //
    for (int a = 0; a < d; a++) {
      int rowN = idMap(eqN(a));

      // A - matrix
      for (int b = 0; b < d; b++) {
        int colN = idMap(eqN(b));
        int ptr = get_col_ptr(com_mod, rowN, colN);

        for (int i = 0; i < 9; i++) {
          Kd (i,ptr) = Kd(i,ptr) + lKd(i,a,b);
        }

        for (int i = 0; i < 3; i++) {
          Val(i,ptr) = Val(i,ptr) + lK (i,a,b);
          Val(i+4,ptr) = Val(i+4,ptr) + lK(i+4,a,b);
          Val(i+8,ptr) = Val(i+8,ptr) + lK(i+8,a,b);
        }
      }

      // B - matrix
      for (int b = 0; b < d; b++) {
        int colN = eqN(b);
        int ptr = get_col_ptr(com_mod, rowN, colN);
        Val(3 ,ptr) = Val(3 ,ptr) + lK(3 ,a,b);
        Val(7 ,ptr) = Val(7 ,ptr) + lK(7 ,a,b);
        Val(11,ptr) = Val(11,ptr) + lK(11,a,b);
      }
    }

    // Divergence matrix (C) is assembled using unmapped rows but
    // mapped columns. Mass matrix (D) is assembled using unmapped
    // rows and columns
    //
    for (int a = 0; a < d; a++) {
      int rowN = eqN(a);

      // C - matrix
      for (int b = 0; b < d; b++) {
        int colN = idMap(eqN(b));
        int ptr = get_col_ptr(com_mod, rowN, colN);

        for (int i = 0; i < 3; i++) {
          Kd(i+9,ptr) = Kd(i+9,ptr) + lKd(i+9,a,b);
          Val(i+12,ptr) = Val(i+12,ptr) + lK(i+12,a,b);
        }
      }

      // D - matrix
      for (int b = 0; b < d; b++) {
        int colN = eqN(b);
        int ptr = get_col_ptr(com_mod, rowN, colN);

        Val(15,ptr) = Val(15,ptr) + lK(15,a,b);
      }
    }

  } else { 

    // Stiffness matrix (A) is assembled using mapped rows and columns
    // Gradient matrix (B) is assembled using mapped row but unmapped
    // column
    //
    for (int a = 0; a < d; a++) {
      int rowN = idMap(eqN(a));

      // A - matrix
      for (int b = 0; b < d; b++) {
        int colN = idMap(eqN(b));
        int ptr = get_col_ptr(com_mod, rowN, colN);

        for (int i = 0; i < 4; i++) {
          Kd(i,ptr) = Kd(i,ptr) + lKd(i,a,b);
        }

        for (int i = 0; i < 2; i++) {
          Val(i,ptr) = Val(i,ptr) + lK(i,a,b);
          Val(i+3,ptr) = Val(i+3,ptr) + lK(i+3,a,b);
        }
      }

      // B - matrix
      for (int b = 0; b < d; b++) {
        int colN = eqN(b);
        int ptr = get_col_ptr(com_mod, rowN, colN);

        Val(2,ptr) = Val(2,ptr) + lK(2,a,b);
        Val(5,ptr) = Val(5,ptr) + lK(5,a,b);
      }
    }

    // Divergence matrix (C) is assembled using unmapped rows but
    // mapped columns. Mass matrix (D) is assembled using unmapped
    // rows and columns
    //
    for (int a = 0; a < d; a++) {
      int rowN = eqN(a);

      // C - matrix
      for (int b = 0; b < d; b++) {
        int colN = idMap(eqN(b));
        int ptr = get_col_ptr(com_mod, rowN, colN);

        for (int i = 0; i < 2; i++) {
          Kd (i+4,ptr) = Kd(i+4,ptr) + lKd(i+4,a,b);
          Val(i+6,ptr) = Val(i+6,ptr) + lK(i+6,a,b);
        }
      }

      // D - matrix
      for (int b = 0; b < d; b++) {
        int colN = eqN(b);
        int ptr = get_col_ptr(com_mod, rowN, colN);
        Val(8,ptr) = Val(8,ptr) + lK(8,a,b);
      }
    }
  }
}

/// Modifies:
///   com_mod.Rd
//
void ustruct_r(ComMod& com_mod, const Array<double>& Yg)
{
  using namespace consts;

  #define n_debug_ustruct_r 
  #ifdef debug_ustruct_r 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];

  if (eq.phys != EquationType::phys_ustruct && eq.phys != EquationType::phys_FSI) {
    return; 
  }

  const int tnNo = com_mod.tnNo;
  const int nsd = com_mod.nsd;

  const auto& Ad = com_mod.Ad;
  const auto& Kd = com_mod.Kd;
  const auto& rowPtr = com_mod.rowPtr;
  const auto& colPtr = com_mod.colPtr;

  auto& R = com_mod.R;
  auto& Rd = com_mod.Rd;

  int s = eq.s;
  double amg = (eq.gam - eq.am) / (eq.gam - 1.0);
  double ami = 1.0 / eq.am;
  #ifdef debug_ustruct_r 
  dmsg << "nsd: " << nsd;
  dmsg << "s: " << s;
  dmsg << "eq.itr: " << eq.itr;
  dmsg << "amg: " << amg;
  dmsg << "ami: " << ami;
  #endif

  if (eq.itr > 1) {
     Rd = 0.0;
  } else {
    for (int a = 0; a < tnNo; a++) { 
      if (!all_fun::is_domain(com_mod, eq, a, EquationType::phys_ustruct)) {
        continue;
      }
      for (int i = 0; i < nsd; i++) {
        Rd(i,a) = amg*Ad(i,a) - Yg(s+i,a);
      }
    }

    if (nsd == 3) {
      Array<double> KU(4,tnNo);

      for (int a = 0; a < tnNo; a++) { 
        if (!all_fun::is_domain(com_mod, eq, a, EquationType::phys_ustruct)) {
          continue;
        }

        for (int i = rowPtr(a); i <= rowPtr(a+1)-1; i++) {
          int c = colPtr(i);

          KU(0,a) = KU(0,a) + Kd(0 ,i)*Rd(0,c) + Kd(1 ,i)*Rd(1,c) + Kd(2 ,i)*Rd(2,c);
          KU(1,a) = KU(1,a) + Kd(3 ,i)*Rd(0,c) + Kd(4 ,i)*Rd(1,c) + Kd(5 ,i)*Rd(2,c);
          KU(2,a) = KU(2,a) + Kd(6 ,i)*Rd(0,c) + Kd(7 ,i)*Rd(1,c) + Kd(8 ,i)*Rd(2,c);
          KU(3,a) = KU(3,a) + Kd(9,i)*Rd(0,c) + Kd(10,i)*Rd(1,c) + Kd(11,i)*Rd(2,c);
        }
      }

      all_fun::commu(com_mod, KU);

      for (int a = 0; a < tnNo; a++) { 
        R(0,a) = R(0,a) - ami*KU(0,a);
        R(1,a) = R(1,a) - ami*KU(1,a);
        R(2,a) = R(2,a) - ami*KU(2,a);
        R(3,a) = R(3,a) - ami*KU(3,a);
      }
    } else {
      Array<double> KU(3,tnNo);

      for (int a = 0; a < tnNo; a++) { 
        if (!all_fun::is_domain(com_mod, eq, a, EquationType::phys_ustruct)) {
          continue;
        }

        for (int i = rowPtr(a); i <= rowPtr(a+1); i++) {
          int c = colPtr(i);
          KU(0,a) = KU(0,a) + Kd(0,i)*Rd(0,c) + Kd(1,i)*Rd(1,c);
          KU(1,a) = KU(1,a) + Kd(2,i)*Rd(0,c) + Kd(3,i)*Rd(1,c);
          KU(2,a) = KU(2,a) + Kd(4,i)*Rd(0,c) + Kd(5,i)*Rd(1,c);
        }
      }

      all_fun::commu(com_mod, KU);

      for (int a = 0; a < tnNo; a++) { 
        R(0,a) = R(0,a) - ami*KU(0,a);
        R(1,a) = R(1,a) - ami*KU(1,a);
        R(2,a) = R(2,a) - ami*KU(2,a);
      }
    }
  } 
} 

};

