// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

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
#include "DebugMsg.h"
#include <array>

namespace struct_ns {

void b_struct_2d(const ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, 
    const Array<double>& Nx, const Array<double>& dl, const Vector<double>& hl, const Vector<double>& nV, 
    Array<double>& lR, Array3<double>& lK)
{
  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  double dt = com_mod.dt;
  int dof = com_mod.dof;

  double af = eq.af * eq.beta * dt *dt;
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

/// @brief Assemble the residual and tangent contributions of one solid mesh.
///
/// @param[in,out] com_mod Global common variables.
/// @param[in] cep_mod Electrophysiology variables, supplying the active stress.
/// @param[in] lM Mesh whose elements are assembled.
/// @param[in] solutions Acceleration, velocity and displacement.
void construct_dsolid(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const SolutionStates& solutions)
{
  const auto& Ag = solutions.intermediate.get_acceleration();
  const auto& Yg = solutions.intermediate.get_velocity();
  const auto& Dg = solutions.intermediate.get_displacement();
  using namespace consts;

  #define n_debug_construct_dsolid
  #ifdef debug_construct_dsolid
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
  Vector<double> pSl(nsymd), N(eNoN);
  ActiveStress::Evaluator active_stress_evaluator;
  Array<double> xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN), 
                bfl(nsd,eNoN), fN(nsd,nFn), pS0l(nsymd,eNoN), Nx(nsd,eNoN), lR(dof,eNoN);
  Array3<double> lK(dof*dof,eNoN,eNoN);

  // Loop over all elements of mesh

  for (int e = 0; e < lM.nEl; e++) {
    // Change the current domain which will be used in later function calls.
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

    if (lM.fN.size() != 0) {
      for (int iFn = 0; iFn < nFn; iFn++) {
        for (int i = 0; i < nsd; i++) {
          fN(i,iFn) = lM.fN(i+nsd*iFn,e);
        }
      }
    }

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
    }

    if (eq.dmn[cDmn].active_stress != nullptr) {
      active_stress_evaluator.update(*eq.dmn[cDmn].active_stress, ptr);
    } else {
      active_stress_evaluator.clear();
    }

    // Gauss integration
    //
    lR = 0.0;
    lK = 0.0;

    double Jac{0.0};
    Array<double> ksix(nsd,nsd);

    for (int g = 0; g < lM.nG; g++) {
      // Shape function gradients and the viscous response are constant
      // within linear triangles and tetrahedra.
      const bool recompute_visc = (g == 0 || !lM.lShpF);

      if (recompute_visc) {
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
        struct_3d(com_mod, cep_mod, eNoN, nFn, w, N, Nx, al, yl, dl, bfl, fN,
                  pS0l, pSl, active_stress_evaluator, lR, lK, recompute_visc);

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
        struct_2d(com_mod, cep_mod, eNoN, nFn, w, N, Nx, al, yl, dl, bfl, fN,
                  pS0l, pSl, active_stress_evaluator, lR, lK, recompute_visc);
      }

      // Prestress
      if (pstEq) {
        for (int a = 0; a < eNoN; a++) {
          int Ac = ptr(a);
          pSa(Ac) += w*N(a);
          for (int i = 0; i < pSn.nrows(); i++) {
            pSn(i,Ac) += w*N(a)*pSl(i);
          }
        }
      }
    } 

    eq.linear_algebra->assemble(com_mod, eNoN, ptr, lK, lR);
  } 
}

/// @brief Reproduces Fortran 'STRUCT2D' subroutine.
//
void struct_2d(ComMod &com_mod, CepMod &cep_mod, const int eNoN, const int nFn,
               const double w, const Vector<double> &N, const Array<double> &Nx,
               const Array<double> &al, const Array<double> &yl,
               const Array<double> &dl, const Array<double> &bfl,
               const Array<double> &fN, const Array<double> &pS0l,
               Vector<double> &pSl, const ActiveStress::Evaluator &active_stress_evaluator,
               Array<double> &lR, Array3<double> &lK, const bool recompute_visc) {
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
  const int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  // Set parameters
  //
  double rho = dmn.prop.at(PhysicalPropertyType::solid_density);
  double dmp = dmn.prop.at(PhysicalPropertyType::damping);
  const Eigen::Vector2d fb{dmn.prop.at(PhysicalPropertyType::f_x),
                           dmn.prop.at(PhysicalPropertyType::f_y)};
  double afu = eq.af * eq.beta*dt*dt;
  double afv = eq.af * eq.gam*dt;
  double amd = eq.am * rho  +  eq.af * eq.gam * dt * dmp;
  double afl = eq.af * eq.beta * dt * dt;

  int i = eq.s;
  #ifdef debug_struct_2d 
  dmsg << "i: " << i;
  dmsg << "amd: " << amd;
  dmsg << "afl: " << afl;
  dmsg << "w: " << w;
  #endif

  // This element's nodal fields, as Eigen views over the caller's storage
  const auto Nxm  = eigen_view<2>(Nx);                    // grad(N_a) per column
  const auto Nm   = eigen_view(N);                        // shape functions
  const auto disp = eigen_view_rows<2>(dl, i);            // nodal displacements
  const auto vel  = eigen_view_rows<2>(yl, i);            // nodal velocities
  const auto acc  = eigen_view_rows<2>(al, i);            // nodal accelerations
  const auto bfm  = eigen_view<2>(bfl);                   // nodal body force
  auto       lRv  = eigen_view_mutable(lR).topRows<2>();  // rows this kernel adds to

  // Inertia, damping and body force: the term the residual weights with N
  const Eigen::Vector2d ud = (rho*(acc - bfm) + dmp*vel) * Nm - rho * fb;

  // Prestress at this Gauss point, in Voigt order [11, 22, 12]
  const Eigen::Vector<double,3> pS0g = eigen_view<3>(pS0l) * Nm;

  Matrix<2> S0;
  S0 << pS0g(0), pS0g(2),
        pS0g(2), pS0g(1);

  // Velocity and deformation gradients: Grad(v) and F = I + Grad(u)
  const Matrix<2> vx = vel * Nxm.transpose();
  const Matrix<2> F  = Matrix<2>::Identity() + disp * Nxm.transpose();

  // Active tension, evaluated here from the fiber stretch of F.
  const auto Ta = active_stress_evaluator.evaluate<2>(Nm, F, eigen_view<2>(fN));
  
  #ifdef debug_struct_2d 
  dmsg << "ud: " << ud(0) << " " << ud(1);
  dmsg << "F: " << F(0,0);
  dmsg << "Ta.fibers: " << Ta.fibers;
  dmsg << "Ta.sheets: " << Ta.sheets;
  dmsg << "Ta.sheet_normals: " << Ta.sheet_normals;
  #endif

  // 2nd Piola-Kirchhoff stress (S) and material stiffness tensor in Voight notation (Dm)
  Matrix<2> S;
  Matrix<3> Dm;
  double Ja;
  mat_models::compute_pk2cc<2>(com_mod, cep_mod, dmn, F, nFn, eigen_view<2>(fN), Ta, S, Dm, Ja);

  // Viscous 2nd Piola-Kirchhoff stress and tangent contributions.
  // Reuse from the previous Gauss point when shape function gradients
  // are constant within an element (e.g. linear triangles, tetrahedra).
  static mat_models::ViscousResponse<2> visc;
  visc.update(dmn, eNoN, Nx, vx, F, recompute_visc);

  // Elastic + Viscous stresses
  S = S + visc.S();

  // Prestress
  pSl(0) = S(0,0);
  pSl(1) = S(1,1);
  pSl(2) = S(0,1);

  // Total 2nd Piola-Kirchhoff stress
  S = S + S0;

  // 1st Piola-Kirchhoff tensor (P)
  //
  const Matrix<2> P = F * S;
  #ifdef debug_struct_2d 
  dmsg << "P: " << P(0,0) << " " << P(0,1);
  dmsg << "   " << P(1,0) << " " << P(1,1);
  #endif

  // Local residual: inertia and body force, plus div P
  lRv += w * (ud * Nm.transpose() + P * Nxm);

  // Strain-displacement matrix; Bm[a] maps node a to Voigt strain
  //
  std::array<Eigen::Matrix<double, 3, 2>, consts::maxNoN> Bm;
  const Matrix<2> Ft = F.transpose();

  for (int a = 0; a < eNoN; a++) {
    const auto g = Nxm.col(a);   // grad(N_a)

    Bm[a].row(0) = g(0) * Ft.row(0);                     // dE_11
    Bm[a].row(1) = g(1) * Ft.row(1);                     // dE_22
    Bm[a].row(2) = g(0) * Ft.row(1) + g(1) * Ft.row(0);  // 2 dE_12
  }

  // Local stiffness tensor
  double T1, NxSNx, BmDBm;

  for (int b = 0; b < eNoN; b++) {

    // Material stiffness for node b
    const Eigen::Matrix<double, 3, 2> DBm = Dm * Bm[b];

    // Geometric stiffness: S*grad(N_b)
    const Eigen::Vector2d SNx = S * Nxm.col(b);

    for (int a = 0; a < eNoN; a++) { 

      // Geometric stiffness
      NxSNx = Nxm.col(a).dot(SNx);
      T1 = amd*N(a)*N(b) + afu*NxSNx;

      // dM1/du1
      BmDBm = Bm[a].col(0).dot(DBm.col(0));
      lK(0,a,b) += w*( T1 + afu*(BmDBm + visc.du(0,a,b)) + afv*visc.dv(0,a,b) );

      // dM1/du2
      BmDBm = Bm[a].col(0).dot(DBm.col(1));
      lK(1,a,b) += w*( afu*(BmDBm + visc.du(1,a,b)) + afv*visc.dv(1,a,b) );

      // dM2/du1
      BmDBm = Bm[a].col(1).dot(DBm.col(0));
      lK(dof+0,a,b) += w*( afu*(BmDBm + visc.du(2,a,b)) + afv*visc.dv(2,a,b) );

      // dM2/du2
      BmDBm = Bm[a].col(1).dot(DBm.col(1));
      lK(dof+1,a,b) += w*( T1 + afu*(BmDBm + visc.du(3,a,b)) + afv*visc.dv(3,a,b) );
    }
  }
}

/// @brief Reproduces Fortran 'STRUCT3D' subroutine.
void struct_3d(ComMod &com_mod, CepMod &cep_mod, const int eNoN, const int nFn,
               const double w, const Vector<double> &N, const Array<double> &Nx,
               const Array<double> &al, const Array<double> &yl,
               const Array<double> &dl, const Array<double> &bfl,
               const Array<double> &fN, const Array<double> &pS0l,
               Vector<double> &pSl, const ActiveStress::Evaluator &active_stress_evaluator,
               Array<double> &lR, Array3<double> &lK, const bool recompute_visc) {          
  using namespace consts;
  using namespace mat_fun;

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
  const int cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  // Set parameters
  //
  double rho = dmn.prop.at(PhysicalPropertyType::solid_density);
  double dmp = dmn.prop.at(PhysicalPropertyType::damping);
  const Eigen::Vector3d fb{dmn.prop.at(PhysicalPropertyType::f_x),
                           dmn.prop.at(PhysicalPropertyType::f_y),
                           dmn.prop.at(PhysicalPropertyType::f_z)};

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

  // This element's nodal fields, as Eigen views over the caller's storage
  const auto Nxm  = eigen_view<3>(Nx);                    // grad(N_a) per column
  const auto Nm   = eigen_view(N);                        // shape functions
  const auto disp = eigen_view_rows<3>(dl, i);            // nodal displacements
  const auto vel  = eigen_view_rows<3>(yl, i);            // nodal velocities
  const auto acc  = eigen_view_rows<3>(al, i);            // nodal accelerations
  const auto bfm  = eigen_view<3>(bfl);                   // nodal body force
  auto       lRv  = eigen_view_mutable(lR).topRows<3>();  // rows this kernel adds to

  // Inertia, damping and body force.
  const Eigen::Vector3d ud = (rho*(acc - bfm) + dmp*vel) * Nm - rho * fb;

  // Prestress at this Gauss point, in Voigt order [11, 22, 33, 12, 23, 31]
  const Eigen::Vector<double,6> pS0g = eigen_view<6>(pS0l) * Nm;

  Matrix<3> S0;
  S0 << pS0g(0), pS0g(3), pS0g(5),
        pS0g(3), pS0g(1), pS0g(4),
        pS0g(5), pS0g(4), pS0g(2);

  // Velocity and deformation gradients: Grad(v) and F = I + Grad(u)
  const Matrix<3> vx = vel * Nxm.transpose();
  const Matrix<3> F  = Matrix<3>::Identity() + disp * Nxm.transpose();

  // Active tension, evaluated here from the fiber stretch of F.
  const auto Ta = active_stress_evaluator.evaluate<3>(Nm, F, eigen_view<3>(fN));

  // 2nd Piola-Kirchhoff tensor (S) and material stiffness tensor in
  // Voigt notation (Dm)
  //
  Matrix<3> S;
  Matrix<6> Dm;
  double Ja;
  mat_models::compute_pk2cc<3>(com_mod, cep_mod, dmn, F, nFn, eigen_view<3>(fN), Ta, S, Dm, Ja);

  // Viscous 2nd Piola-Kirchhoff stress and tangent contributions.
  // Reuse from the previous Gauss point when shape function gradients
  // are constant within an element (e.g. linear triangles, tetrahedra).
  static mat_models::ViscousResponse<3> visc;
  visc.update(dmn, eNoN, Nx, vx, F, recompute_visc);

  // Elastic + Viscous stresses
  S = S + visc.S();

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
  const Matrix<3> P = F * S;

  // Local residual: inertia and body force, plus div P
  lRv += w * (ud * Nm.transpose() + P * Nxm);

  // Strain-displacement matrix; Bm[a] maps node a to Voigt strain
  //
  std::array<Eigen::Matrix<double, 6, 3>, consts::maxNoN> Bm;
  const Matrix<3> Ft = F.transpose();

  for (int a = 0; a < eNoN; a++) {
    const auto g = Nxm.col(a);   // grad(N_a)

    Bm[a].row(0) = g(0) * Ft.row(0);                     // dE_11
    Bm[a].row(1) = g(1) * Ft.row(1);                     // dE_22
    Bm[a].row(2) = g(2) * Ft.row(2);                     // dE_33
    Bm[a].row(3) = g(0) * Ft.row(1) + g(1) * Ft.row(0);  // 2 dE_12
    Bm[a].row(4) = g(1) * Ft.row(2) + g(2) * Ft.row(1);  // 2 dE_23
    Bm[a].row(5) = g(2) * Ft.row(0) + g(0) * Ft.row(2);  // 2 dE_31
  }

  // Local stiffness tensor
  double NxSNx, T1, BmDBm;

  for (int b = 0; b < eNoN; b++) {

    // Material stiffness for node b
    const Eigen::Matrix<double, 6, 3> DBm = Dm * Bm[b];

    // Geometric stiffness: S*grad(N_b)
    const Eigen::Vector3d SNx = S * Nxm.col(b);

    for (int a = 0; a < eNoN; a++) {

      NxSNx = Nxm.col(a).dot(SNx);
      T1 = amd*N(a)*N(b) + afu*NxSNx;

      // dM1/du1
      BmDBm = Bm[a].col(0).dot(DBm.col(0));
      lK(0,a,b) += w*( T1 + afu*(BmDBm + visc.du(0,a,b)) + afv*visc.dv(0,a,b) );

      // dM1/du2
      BmDBm = Bm[a].col(0).dot(DBm.col(1));
      lK(1,a,b) += w*( afu*(BmDBm + visc.du(1,a,b)) + afv*visc.dv(1,a,b) );

      // dM1/du3
      BmDBm = Bm[a].col(0).dot(DBm.col(2));
      lK(2,a,b) += w*( afu*(BmDBm + visc.du(2,a,b)) + afv*visc.dv(2,a,b) );

      // dM2/du1
      BmDBm = Bm[a].col(1).dot(DBm.col(0));
      lK(dof+0,a,b) += w*( afu*(BmDBm + visc.du(3,a,b)) + afv*visc.dv(3,a,b) );

      // dM2/du2
      BmDBm = Bm[a].col(1).dot(DBm.col(1));
      lK(dof+1,a,b) += w*(T1 + afu*(BmDBm + visc.du(4,a,b)) + afv*visc.dv(4,a,b) );

      // dM2/du3
      BmDBm = Bm[a].col(1).dot(DBm.col(2));
      lK(dof+2,a,b) += w*( afu*(BmDBm + visc.du(5,a,b)) + afv*visc.dv(5,a,b) );

      // dM3/du1
      BmDBm = Bm[a].col(2).dot(DBm.col(0));
      lK(2*dof+0,a,b) += w*( afu*(BmDBm + visc.du(6,a,b)) + afv*visc.dv(6,a,b) );

      // dM3/du2
      BmDBm = Bm[a].col(2).dot(DBm.col(1));
      lK(2*dof+1,a,b) += w*( afu*(BmDBm + visc.du(7,a,b)) + afv*visc.dv(7,a,b) );

      // dM3/du3
      BmDBm = Bm[a].col(2).dot(DBm.col(2));
      lK(2*dof+2,a,b) += w*( T1 + afu*(BmDBm + visc.du(8,a,b)) + afv*visc.dv(8,a,b) );
    }
  }
}
};

