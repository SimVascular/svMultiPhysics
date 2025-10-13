// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "cep.h"

#include "all_fun.h"
#include "lhsa.h"
#include "mat_fun.h"
#include "nn.h"
#include "utils.h"

#include <math.h>

namespace cep {

void b_cep(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const double h, Array<double>& lR)
{
  double f = w*h;

  // Here the loop is started for constructing left and right hand side
  for (int a = 0; a < eNoN; a++) {
    lR(0,a) = lR(0,a) + N(a)*f;
  }

}

/// @brief This is for solving 1D electrophysiology diffusion equation for Purkinje fibers
///
/// Reproduces Fortran 'CEP1D' subroutine.
// 
void cep_1d(ComMod& com_mod, CepMod& cep_mod, const int eNoN, const int nFn, const double w,
    const Vector<double>& N, const Array<double>& Nx, const Array<double>& al, const Array<double>& yl,
    Array<double>& lR, Array3<double>& lK)
{
  #define n_debug_cep_1d 
  #ifdef debug_cep_1d 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "w: " << w;
  #endif

  using namespace consts;
  using namespace mat_fun;

  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  auto& cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;

  double T1 = eq.af * eq.gam * dt;
  double amd = eq.am / T1;
  double Diso = dmn.cep.Diso;
  int i = eq.s;
  double wl = w*T1;
  #ifdef debug_cep_1d 
  dmsg << "T1: " << T1;
  dmsg << "amd: " << amd;
  dmsg << "Diso: " << Diso;
  dmsg << "i: " << i;
  #endif

  double Td = 0.0;
  double Tx = 0.0;
  Vector<double> DNx(eNoN);

  for (int a = 0; a < eNoN; a++) {
    Td = Td + N(a)*al(i,a);
    Tx = Tx + Nx(0,a)*yl(i,a);
    DNx(a) = Diso*Nx(0,a);
  }

  for (int a = 0; a < eNoN; a++) {
    lR(0,a) = lR(0,a) + w*(N(a)*Td + Nx(0,a)*Diso*Tx);

    for (int b = 0; b < eNoN; b++) {
      lK(0,a,b) = lK(0,a,b) + wl*(N(a)*N(b)*amd + Nx(0,a)*DNx(b));
    }
  }
}

//--------
// cep_2d
//--------
// Reproduces Fortran 'CEP2D' subroutine.
// 
void cep_2d(ComMod& com_mod, CepMod& cep_mod, const int eNoN, const int nFn, const double w,
    const Vector<double>& N, const Array<double>& Nx, const Array<double>& al, const Array<double>& yl,
    const Array<double>& dl, const Array<double>& fN, Array<double>& lR, Array3<double>& lK)
{
  #define n_debug_cep_2d 
  #ifdef debug_cep_2d 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  using namespace consts;
  using namespace mat_fun;

  const int nsd = com_mod.nsd;
  const int tDof = com_mod.tDof;
  const int dof = com_mod.dof;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  auto& cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;
  const auto& cem = cep_mod.cem;

  Vector<double> Dani(nFn), Vx(2), Ls(nFn), DVx(2);
  Array<double> F(2,2), C(2,2), fl(2,nFn), D(2,2), DNx(2,eNoN);

  if (nFn < dmn.cep.nFn) { 
    throw std::runtime_error("[cep_2d] No. of anisotropic conductivies exceed mesh fibers.");
  }

  double T1 = eq.af * eq.gam * dt;
  double amd = eq.am / T1;
  double wl = w * T1;
  double Diso = dmn.cep.Diso;
  #ifdef debug_cep_2d 
  dmsg << "Diso: " << Diso;
  #endif

  for (int i = 0; i < nFn; i++) {
    if (i+1 <= dmn.cep.nFn) {
      Dani(i) = dmn.cep.Dani(i);
    } else { 
      Dani(i) = Dani(i-1);
    }
  }

  // Compute the isotropic part of diffusion tensor based on spatial
  // isotropy for electromechanics. This models stretch induced changes
  // in conduction velocities
  //
  Ls = 1.0;
  int i;

  if (cem.cpld) { 
    for (int a = 0; a < com_mod.nEq; a++) { 
      if (com_mod.eq[a].phys == EquationType::phys_struct || 
          com_mod.eq[a].phys == EquationType::phys_ustruct) {
        i = com_mod.eq[a].s;
        break;
      }
    }

    // Compute deformation gradient tensor
    //
    F(0,0) = 1.0;
    F(1,1) = 1.0;

    for (int a = 0; a < eNoN; a++) {
      F(0,0) = F(0,0) + Nx(0,a)*dl(i,a);
      F(0,1) = F(0,1) + Nx(1,a)*dl(i,a);
      F(1,0) = F(1,0) + Nx(0,a)*dl(i+1,a);
      F(1,1) = F(1,1) + Nx(1,a)*dl(i+1,a);
    }

    // Jacobian
    double Jac = mat_fun::mat_det(F, 2);

    // Compute Cauchy-Green tensor and its inverse
    C = mat_mul(transpose(F), F);
    C = mat_inv(C, 2);

    // Compute fiber stretch
    for (int i = 0; i < nFn; i++) {
      Ls(i) = sqrt(utils::norm(fN.col(i), mat_mul(C, fN.col(i))));
      for (int j = 0; j < 2; j++) {
        fl(j,i) = fN(j,i) / Ls(i);
      }
    }

    if (Ls(0) <= 1.0) {
      Ls(0) = 1.0;
    }

    // Diffusion tensor - spatial isotropy
    //
    Diso = Diso * Jac;
    Dani = Dani * Jac;
    D = Diso * C;

  } else { 
    D  = 0.0;
    D(0,0) = Diso;
    D(1,1) = Diso;
    fl = fN;
  }

  for (int i = 0 ; i < nFn; i++) {
     D(0,0) = D(0,0) + Dani(i)*fl(0,i)*fl(0,i);
     D(0,1) = D(0,1) + Dani(i)*fl(0,i)*fl(1,i);

     D(1,0) = D(1,0) + Dani(i)*fl(1,i)*fl(0,i);
     D(1,1) = D(1,1) + Dani(i)*fl(1,i)*fl(1,i);
  }

  i = eq.s;
  double Vd = 0.0;
  Vx = 0.0;

  for (int a = 0; a < eNoN; a++) {
    Vd = Vd + N(a)*al(i,a);

    Vx(0) = Vx(0) + Nx(0,a)*yl(i,a);
    Vx(1) = Vx(1) + Nx(1,a)*yl(i,a);

    DNx(0,a) = D(0,0)*Nx(0,a) + D(0,1)*Nx(1,a);
    DNx(1,a) = D(1,0)*Nx(0,a) + D(1,1)*Nx(1,a);
  }

  DVx(0) = D(0,0)*Vx(0) + D(0,1)*Vx(1);
  DVx(1) = D(1,0)*Vx(0) + D(1,1)*Vx(1);

  for (int a = 0; a < eNoN; a++) {
    lR(0,a) = lR(0,a) + w*(N(a)*Vd + Nx(0,a)*DVx(0) + Nx(1,a)*DVx(1));

    for (int b = 0; b < eNoN; b++) {
      lK(0,a,b) = lK(0,a,b) + wl*(N(a)*N(b)*amd + Nx(0,a)*DNx(0,b) + Nx(1,a)*DNx(1,b));
    }
  }
}

//--------
// cep_3d
//--------
// Reproduces Fortran 'CEP3D' subroutine.
// 
void cep_3d(ComMod& com_mod, CepMod& cep_mod, const int eNoN, const int nFn, const double w,
    const Vector<double>& N, const Array<double>& Nx, const Array<double>& al, const Array<double>& yl,
    const Array<double>& dl, const Array<double>& fN, Array<double>& lR, Array3<double>& lK)
{
  #define n_debug_cep_3d 
  #ifdef debug_cep_3d 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "nFn: " << nFn;
  dmsg << "w: " << w;
  #endif

  using namespace consts;
  using namespace mat_fun;

  const int nsd = com_mod.nsd;
  const int tDof = com_mod.tDof;
  const int dof = com_mod.dof;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  auto& cDmn = com_mod.cDmn;
  auto& dmn = eq.dmn[cDmn];
  const double dt = com_mod.dt;
  const auto& cem = cep_mod.cem;

  Vector<double> Dani(nFn), Vx(3), Ls(nFn), DVx(3);
  Array<double> F(3,3), C(3,3), fl(3,nFn), D(3,3), DNx(3,eNoN);

  double T1 = eq.af * eq.gam * dt;
  double amd = eq.am / T1;
  double wl = w * T1;
  double Diso = dmn.cep.Diso;
  #ifdef debug_cep_3d 
  dmsg << "T1: " << T1;
  dmsg << "amd: " << amd;
  dmsg << "wl: " << wl;
  dmsg << "Diso: " << Diso;
  #endif

  for (int i = 0; i < nFn; i++) {
    if (i+1 <= dmn.cep.nFn) {
      Dani(i) = dmn.cep.Dani(i);
    } else {
      Dani(i) = Dani(i-1);
    }
  }

  // Compute the isotropic part of diffusion tensor based on spatial
  // isotropy for electromechanics. This models stretch induced changes
  // in conduction velocities
  //
  Ls = 1.0;
  int i = 0;

  if (cem.cpld) {
    // Get the displacement degrees of freedom
    for (int a = 0; a < com_mod.nEq; a++) {
      if (com_mod.eq[a].phys == EquationType::phys_struct || 
          com_mod.eq[a].phys == EquationType::phys_ustruct) {
        i = com_mod.eq[a].s;
        break; 
      }
    }

    // Compute deformation gradient tensor
    //
    F(0,0) = 1.0;
    F(1,1) = 1.0;
    F(2,2) = 1.0;

    for (int a = 0; a < eNoN; a++) {
      F(0,0) = F(0,0) + Nx(0,a)*dl(i,a);
      F(0,1) = F(0,1) + Nx(1,a)*dl(i,a);
      F(0,2) = F(0,2) + Nx(2,a)*dl(i,a);
      F(1,0) = F(1,0) + Nx(0,a)*dl(i+1,a);
      F(1,1) = F(1,1) + Nx(1,a)*dl(i+1,a);
      F(1,2) = F(1,2) + Nx(2,a)*dl(i+1,a);
      F(2,0) = F(2,0) + Nx(0,a)*dl(i+2,a);
      F(2,1) = F(2,1) + Nx(1,a)*dl(i+2,a);
      F(2,2) = F(2,2) + Nx(2,a)*dl(i+2,a);
    }

    // Jacobian
    double Jac = mat_fun::mat_det(F, 3);

    // Compute Cauchy-Green tensor and its inverse
    C = mat_mul(transpose(F), F);
    C = mat_inv(C, 3);

    // Compute fiber stretch
    for (int i = 0; i < nFn; i++) {
      Ls(i) = sqrt(utils::norm(fN.col(i), mat_mul(C, fN.col(i))));
      for (int j = 0; j < 3; j++) {
        fl(j,i) = fN(j,i) / Ls(i);
      }
    }

    if (Ls(0) <= 1.0) {
      Ls(0) = 1.0;
    }

    // Diffusion tensor - spatial isotropy
    //
    Diso = Diso * Jac;
    Dani = Dani * Jac;
    D = Diso * C;

  } else {
    D(0,0)  = Diso;
    D(1,1)  = Diso;
    D(2,2)  = Diso;
    fl= fN;
  }

  // Compute anisotropic components of diffusion tensor
  //
  for (int i = 0 ; i < nFn; i++) {
    D(0,0) = D(0,0) + Dani(i)*fl(0,i)*fl(0,i);
    D(0,1) = D(0,1) + Dani(i)*fl(0,i)*fl(1,i);
    D(0,2) = D(0,2) + Dani(i)*fl(0,i)*fl(2,i);

    D(1,0) = D(1,0) + Dani(i)*fl(1,i)*fl(0,i);
    D(1,1) = D(1,1) + Dani(i)*fl(1,i)*fl(1,i);
    D(1,2) = D(1,2) + Dani(i)*fl(1,i)*fl(2,i);

    D(2,0) = D(2,0) + Dani(i)*fl(2,i)*fl(0,i);
    D(2,1) = D(2,1) + Dani(i)*fl(2,i)*fl(1,i);
    D(2,2) = D(2,2) + Dani(i)*fl(2,i)*fl(2,i);
  }

  i = eq.s;
  double Vd = 0.0;
  Vx = 0.0;

  for (int a = 0; a < eNoN; a++) {
     Vd = Vd + N(a)*al(i,a);

     Vx(0) = Vx(0) + Nx(0,a)*yl(i,a);
     Vx(1) = Vx(1) + Nx(1,a)*yl(i,a);
     Vx(2) = Vx(2) + Nx(2,a)*yl(i,a);

     DNx(0,a) = D(0,0)*Nx(0,a) + D(0,1)*Nx(1,a) + D(0,2)*Nx(2,a);
     DNx(1,a) = D(1,0)*Nx(0,a) + D(1,1)*Nx(1,a) + D(1,2)*Nx(2,a);
     DNx(2,a) = D(2,0)*Nx(0,a) + D(2,1)*Nx(1,a) + D(2,2)*Nx(2,a);
  }

  DVx(0) = D(0,0)*Vx(0) + D(0,1)*Vx(1) + D(0,2)*Vx(2);
  DVx(1) = D(1,0)*Vx(0) + D(1,1)*Vx(1) + D(1,2)*Vx(2);
  DVx(2) = D(2,0)*Vx(0) + D(2,1)*Vx(1) + D(2,2)*Vx(2);

  for (int a = 0; a < eNoN; a++) {
    lR(0,a) = lR(0,a) + w*(N(a)*Vd + Nx(0,a)*DVx(0) + Nx(1,a)*DVx(1) + Nx(2,a)*DVx(2));

    for (int b = 0; b < eNoN; b++) {
      lK(0,a,b) = lK(0,a,b) + wl*(N(a)*N(b)*amd + Nx(0,a)*DNx(0,b) + Nx(1,a)*DNx(1,b) + Nx(2,a)*DNx(2,b));
    }
  }
}

//---------------
// construct_cep
//---------------
//
void construct_cep(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const Array<double>& Ag, 
    const Array<double>& Yg, const Array<double>& Dg)
{
  #define n_debug_construct_cep 
  #ifdef debug_construct_cep 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  using namespace consts;

  const int nsd = com_mod.nsd;
  const int tDof = com_mod.tDof;
  const int dof = com_mod.dof;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  auto& cDmn = com_mod.cDmn;

  int insd = nsd;
  const int eNoN = lM.eNoN;
  int nFn  = lM.nFn;

  if (lM.lFib) insd = 1;
  if (nFn == 0) nFn = 1;
  #ifdef debug_construct_cep 
  dmsg << "lM.nEl: " << lM.nEl;
  dmsg << "insd: " << insd;
  dmsg << "tDof: " << tDof;
  dmsg << "eNoN: " << eNoN;
  dmsg << "Dg.nrows: " << Dg.nrows();
  dmsg << "Dg.ncols: " << Dg.ncols();
  #endif

  // CEP: dof = 1
  Vector<int> ptr(eNoN); 
  Array<double> xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN), 
      fN(nsd,nFn), Nx(insd,eNoN), lR(dof,eNoN);
  Array3<double> lK(dof*dof,eNoN,eNoN);
  Vector<double>  N(eNoN); 
  
  // ECG computation
  Vector<double> pseudo_ECG_proc(cep_mod.ecgleads.num_leads);
  Vector<double> Vx(3);
  double x_coords;
  double y_coords;
  double z_coords;

  pseudo_ECG_proc = 0.0;

  // Loop over all elements of mesh
  for (int e = 0; e < lM.nEl; e++) {
    cDmn = all_fun::domain(com_mod, lM, cEq, e);
    auto cPhys = eq.dmn[cDmn].phys;
    if (cPhys != EquationType::phys_CEP) {
      continue;
    }

    // Update shape functions for NURBS
    //if (lM.eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

    // Create local copies
    fN = 0.0;

    for (int a = 0; a < eNoN; a++) {
      int Ac = lM.IEN(a,e);
      ptr(a) = Ac;

      for (int i = 0; i < nsd; i++) {
        xl(i,a) = com_mod.x(i,Ac);
      }

      for (int i = 0; i < tDof; i++) {
        al(i,a) = Ag(i,Ac);
        yl(i,a) = Yg(i,Ac);
        dl(i,a) = Dg(i,Ac);
      }

      if (lM.fN.size() != 0) {
        for (int iFn = 0; iFn < nFn; iFn++) { 
          for (int i = 0; i < nsd; i++) {
            fN(i,iFn) = lM.fN(i+nsd*iFn,e);
          }
        }
      }
    }

    // Gauss integration
    lR = 0.0;
    lK = 0.0;
    double Jac{0.0};
    Array<double> ksix(nsd,nsd);

    for (int g = 0; g < lM.nG; g++) {
      if (g == 0 || !lM.lShpF) {
        auto Nx_g = lM.Nx.slice(g);
        nn::gnn(eNoN, nsd, insd, Nx_g, xl, Nx, Jac, ksix);
        if (utils::is_zero(Jac)) {
          throw std::runtime_error("[construct_cep] Jacobian for element " + std::to_string(e) + " is < 0.");
        }
      }

      double w = lM.w(g) * Jac;
      N = lM.N.col(g);

      #ifdef debug_construct_cep 
      dmsg << "   " << " ";
      dmsg << "g: " << g+1;
      dmsg << "lM.Nx.slice(g): " << lM.Nx.slice(g);
      dmsg << "xl: " << xl;
      dmsg << "Nx: " << Nx;
      dmsg << "Jac: " << Jac;
      dmsg << "lM.w(g): " << lM.w(g);
      #endif

      if (insd == 3) {
        cep_3d(com_mod, cep_mod, eNoN, nFn, w, N, Nx, al, yl, dl, fN, lR, lK);

      } else if (insd == 2) {
        cep_2d(com_mod, cep_mod, eNoN, nFn, w, N, Nx, al, yl, dl, fN, lR, lK);

      } else if (insd == 1) {
        cep_1d(com_mod, cep_mod, eNoN, nFn, w, N, Nx, al, yl, lR, lK);
      }

      // ECG computation
      if (cep_mod.ecgleads.num_leads) {
        // Compute transmembrane gauss points location and potential space derivative
        x_coords = 0.0;
        y_coords = 0.0;
        z_coords = 0.0;
        Vx       = 0.0;
        for (int a = 0; a < eNoN; a++) {
          x_coords += N(a) * xl(0,a);
          y_coords += N(a) * xl(1,a);
          z_coords += N(a) * xl(2,a);
          Vx(0) += Nx(0,a) * yl(0,a);
          Vx(1) += Nx(1,a) * yl(0,a);
          Vx(2) += Nx(2,a) * yl(0,a);
        }

        // Compute integral from Equation (8) in Costabal, Yao, Kuhl 2018
        for (int index = 0; index < cep_mod.ecgleads.num_leads; index++) {
          double r_sq = (x_coords * x_coords + y_coords * y_coords + z_coords * z_coords
                        - 2 * (x_coords * cep_mod.ecgleads.x_coords[index] +
                               y_coords * cep_mod.ecgleads.y_coords[index] +
                               z_coords * cep_mod.ecgleads.z_coords[index])
                        + cep_mod.ecgleads.x_coords[index] * cep_mod.ecgleads.x_coords[index]
                        + cep_mod.ecgleads.y_coords[index] * cep_mod.ecgleads.y_coords[index]
                        + cep_mod.ecgleads.z_coords[index] * cep_mod.ecgleads.z_coords[index]);

          double drinv_x = pow(r_sq, -3./2.) * (cep_mod.ecgleads.x_coords[index] - x_coords);
          double drinv_y = pow(r_sq, -3./2.) * (cep_mod.ecgleads.y_coords[index] - y_coords);
          double drinv_z = pow(r_sq, -3./2.) * (cep_mod.ecgleads.z_coords[index] - z_coords);

          pseudo_ECG_proc(index) += w * (-Vx(0) * drinv_x
                                         -Vx(1) * drinv_y
                                         -Vx(2) * drinv_z);
        }
      }

    } 

    // Assembly
    eq.linear_algebra->assemble(com_mod, eNoN, ptr, lK, lR);
  }

  // Communications among processors for ECG leads computation
  if (cep_mod.ecgleads.num_leads) {
    MPI_Reduce(pseudo_ECG_proc.data(),
               cep_mod.ecgleads.pseudo_ECG.data(),
               cep_mod.ecgleads.num_leads,
               cm_mod::mpreal,
               MPI_SUM,
               0,
               com_mod.cm.com());
  }
}

};
