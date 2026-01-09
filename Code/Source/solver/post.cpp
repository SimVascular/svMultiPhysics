// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "post.h"

#include "all_fun.h"
#include "fluid.h"
#include "fs.h"
#include "initialize.h"
#include "mat_fun.h"
#include "mat_models.h"
#include "nn.h"
#include "shells.h"
#include "utils.h"
#include "vtk_xml.h"
#include <math.h>

namespace post {

void all_post(Simulation* simulation, Array<double>& res, const Array<double>& lY, const Array<double>& lD, 
    consts::OutputNameType outGrp, const int iEq) 
{
  using namespace consts;

  auto& com_mod = simulation->com_mod;
  auto& cm = com_mod.cm;
  auto& cm_mod = simulation->cm_mod;

  #define n_dbug_all_post
  #ifdef dbug_all_post
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "outGrp: " << outGrp;
  #endif

  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    auto& msh = com_mod.msh[iM];
    Array<double> tmpV(maxNSD,msh.nNo);

    if (outGrp == OutputNameType::outGrp_WSS ||  outGrp == OutputNameType::outGrp_trac) {
      bpost(simulation, msh,  tmpV, lY, lD, outGrp);
      for (int a = 0; a < com_mod.msh[iM].nNo; a++) {
        int Ac = msh.gN(a);
        res.set_col(Ac, tmpV.col(a));
      }

    } else if (outGrp == OutputNameType::outGrp_J) {
      Array<double> tmpV(1,msh.nNo); 
      Vector<double> tmpVe(msh.nEl);
      tpost(simulation, msh, 1, tmpV, tmpVe, lD, lY, iEq, outGrp);
      res = 0.0;
      for (int a = 0; a < com_mod.msh[iM].nNo; a++) {
        int Ac = msh.gN(a);
         res(0,Ac) = tmpV(0,a);
      }

     } else if (outGrp == OutputNameType::outGrp_mises) {
       Array<double> tmpV(1,msh.nNo); 
       Vector<double> tmpVe(msh.nEl);
       tpost(simulation, msh, 1, tmpV, tmpVe, lD, lY, iEq, outGrp);
       res = 0.0;
       for (int a = 0; a < com_mod.msh[iM].nNo; a++) {
         int Ac = msh.gN(a);
         res(0,Ac) = tmpV(0,a);
       }

     } else if (outGrp ==  OutputNameType::outGrp_divV) {
       Array<double> tmpV(1,msh.nNo); 
       div_post(simulation, msh, tmpV, lY, lD, iEq);
       res = 0.0;
       for (int a = 0; a < com_mod.msh[iM].nNo; a++) {
         int Ac = msh.gN(a);
         res(0,Ac) = tmpV(0,a);
       }

     } else {
       post(simulation, msh, tmpV, lY, lD, outGrp, iEq);
       for (int a = 0; a < com_mod.msh[iM].nNo; a++) {
         int Ac = msh.gN(a);
         res.set_col(Ac, tmpV.col(a));
        }
     }
  }
}

/// @brief General purpose routine for post processing outputs at the
/// faces. Currently this calculates WSS, which is t.n - (n.t.n)n
/// Here t is stress tensor: t = \mu (grad(u) + grad(u)^T)
//
void bpost(Simulation* simulation, const mshType& lM, Array<double>& res, const Array<double>& lY, const Array<double>& lD, 
    consts::OutputNameType outGrp)
{
  using namespace consts;

  auto& com_mod = simulation->com_mod;
  auto& cm = com_mod.cm;
  auto& cm_mod = simulation->cm_mod;

  #define n_debug_bpost
  #ifdef debug_bpost
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "outGrp: " << outGrp;
#endif

  if ((outGrp != OutputNameType::outGrp_WSS) && (outGrp != OutputNameType::outGrp_trac)) {
    throw std::runtime_error("Invalid output group. Correction is required in BPOST");
  }

  int iEq = 0;
  int eNoN = lM.eNoN;

  auto& eq = com_mod.eq[iEq];
  bool FSIeq = false; 

  if (eq.phys == EquationType::phys_FSI) {
    FSIeq = true; 
  }

  const int tnNo = com_mod.tnNo;
  const int nsd = com_mod.nsd;

  Vector<double> sA(tnNo); 
  Array<double> sF(maxNSD,tnNo); 
  Array<double> xl(nsd,eNoN); 
  Array<double> ul(nsd,eNoN); 
  Array<double> gnV(nsd,tnNo); 
  Array<double> lnV(nsd,eNoN); 
  Vector<double> N(eNoN); 
  Array<double> Nx(nsd,eNoN);

  // First creating the norm field
  //
  for (int iFa = 0; iFa < lM.nFa; iFa++) {
    auto& fa = lM.fa[iFa];

    for (int a = 0; a < fa.nNo; a++) {
      int Ac = fa.gN(a);
      for (int i = 0; i < nsd; i++) {
        gnV(i,Ac) = fa.nV(i,a);
      }
    }
  }

  // Update pressure function spaces
  //
  fsType fsP;

  if (lM.nFs == 1) {
    fsP.nG = lM.fs[0].nG;
    fsP.eType = lM.fs[0].eType;
    fsP.eNoN = lM.fs[0].eNoN;

    fs::alloc_fs(fsP, nsd, nsd);
 
    fsP.w = lM.fs[0].w;
    fsP.xi = lM.fs[0].xi;
    fsP.N = lM.fs[0].N;
    fsP.Nx = lM.fs[0].Nx;

  } else {
    fsP.nG = lM.fs[0].nG;
    fsP.eType = lM.fs[1].eType;
    fsP.eNoN = lM.fs[1].eNoN;

    fs::alloc_fs(fsP, nsd, nsd);

    fsP.xi = lM.fs[0].xi;

    for (int g = 0; g < fsP.nG; g++) {
      nn::get_gnn(nsd, fsP.eType, fsP.eNoN, g, fsP.xi, fsP.N, fsP.Nx);
    }
  }

  Vector<double> pl(fsP.eNoN);

  for (int iFa = 0; iFa < lM.nFa; iFa++) {
    auto& fa = lM.fa[iFa];

    for (int e = 0; e < fa.nEl; e++) {
      int Ec = fa.gE(e);
      int cDmn = all_fun::domain(com_mod, lM, iEq, Ec);
      if (cDmn == -1) {
        continue;
      }
      if (lM.eType == ElementType::NRB) {
        // [TODO:DaveP] not implemented.
        //CALL NRBNNX(lM, Ec)
      }

      // Finding the norm for all the nodes of this element, including
      // those that don't belong to this face, which will be inerpolated
      // from the nodes of the face
      //
      Vector<double> nV(nsd);

      for (int a = 0; a < eNoN; a++) {
        int Ac = lM.IEN(a,Ec);

        for (int i = 0; i < nsd; i++) {
          lnV(i,a) = gnV(i,Ac);
          nV(i) = nV(i) + lnV(i,a);
          xl(i,a) = com_mod.x(i,Ac);

          if (FSIeq) {
            xl(i,a) = xl(i,a) + lD(i+nsd+1,Ac);
            ul(i,a) = lY(i,Ac) - lY(i+nsd+1,Ac);
          } else {
            ul(i,a) = lY(i,Ac);
          }
        }
      }

      for (int a = 0; a < fsP.eNoN; a++) {
        int Ac = lM.IEN(a,Ec);
        pl(a) = lY(nsd,Ac);
      }

      nV = nV / fa.eNoN;
      double Jac;

      for (int a = 0; a < eNoN; a++) {
        if (utils::is_zero(utils::norm(lnV.col(a)))) {
          lnV.set_col(a, nV);
        }
      }

      Array<double> ks(nsd,nsd);

      for (int g = 0; g < lM.nG; g++) {
        if (g == 0 || !lM.lShpF) {
          auto lM_Nx = lM.Nx.slice(g);
          nn::gnn(eNoN, nsd, nsd, lM_Nx, xl, Nx, Jac, ks);
        }

        double w = lM.w(g)*Jac;
        auto N = lM.N.col(g);

        // Calculating ux = grad(u) and nV at a Gauss point
        //
        Vector<double> Tdn(nsd); 
        Vector<double> taue(nsd);         
        Array<double> ux(nsd,nsd);
        nV = 0.0;

        for (int a = 0; a < eNoN; a++) {
          nV = nV + N(a)*lnV.col(a);
          for (int i = 0; i < nsd; i++) {
            for (int j = 0; j < nsd; j++) {
              ux(i,j) = ux(i,j) + Nx(i,a)*ul(j,a);
            }
          }
        }

        double p = 0.0;
        for (int a = 0; a < fsP.eNoN; a++) {
          p = p + fsP.N(a,g)*pl(a);
        }

        // Shear rate, gam := (2*e_ij*e_ij)^0.5
        //
        double gam = 0.0;
        for (int i = 0; i < nsd; i++) {
          for (int j = 0; j < nsd; j++) {
            gam = gam + (ux(i,j)+ux(j,i))*(ux(i,j)+ux(j,i));
          }
        }
        gam = sqrt(0.5*gam);

        // Compute viscosity
        double mu, mu_s;
        fluid::get_viscosity(com_mod, eq.dmn[cDmn], gam, mu, mu_s, mu_s);

        // Now finding grad(u).n and n.grad(u).n
        //
        double ndTdn = 0.0;

        for (int i = 0; i < nsd; i++) {
          for (int j = 0; j < nsd; j++) {
            Tdn(i) = Tdn(i) + mu*(ux(i,j) + ux(j,i))*nV(j);
          }

          ndTdn = ndTdn + Tdn(i)*nV(i);
        }

        taue = Tdn - ndTdn*nV;
        Vector<double> lRes(maxNSD);

        if (outGrp == OutputNameType::outGrp_WSS) {
          for (int i = 0; i < nsd; i++) {
            lRes(i) = -taue(i);
          }

        } else if (outGrp == OutputNameType::outGrp_trac) {
          for (int i = 0; i < nsd; i++) {
            lRes(i) = p*nV(i) - Tdn(i);
          } 
        }

        // Mapping Tau into the nodes by assembling it into a local vector
        //
        for (int a = 0; a < eNoN; a++) {
          int Ac = lM.IEN(a,Ec);
          sA(Ac) = sA(Ac) + w*N(a);

          for (int i = 0; i < maxNSD; i++) {
            sF(i,Ac) = sF(i,Ac) + w*N(a)*lRes(i);
          }
        }
      }
    }
  }

  all_fun::commu(com_mod, sF);
  all_fun::commu(com_mod, sA);

  for (int iFa = 0; iFa < lM.nFa; iFa++) {
    for (int a = 0; a < lM.fa[iFa].nNo; a++) {
      int Ac = lM.fa[iFa].gN(a);
      if (!utils::is_zero(sA(Ac))) {
        for (int i = 0; i < maxNSD; i++) {
          sF(i,Ac) = sF(i,Ac) / sA(Ac);
        }
        sA(Ac) = 1.0;
      }
    }
  }

  res = 0.0;

  for (int a = 0; a < lM.nNo; a++) {
    int Ac = lM.gN(a);
    for (int i = 0; i < maxNSD; i++) {
      res(i,a) = sF(i,Ac);
    }
  }
}

void div_post(Simulation* simulation, const mshType& lM, Array<double>& res, const Array<double>& lY, const Array<double>& lD,
    const int iEq)
{
  using namespace consts;

  #define n_debug_div_post
  #ifdef debug_div_post
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  auto& com_mod = simulation->com_mod;
  auto& cm = com_mod.cm;
  auto& cm_mod = simulation->cm_mod;
  auto& eq = com_mod.eq[iEq];

  // [NOTE] Setting gobal variable 'dof'.
  com_mod.dof = eq.dof;

  int eNoN = lM.eNoN;
  int i = eq.s;
  int j = i + 1;
  int k = j + 1;

  const int nsd = com_mod.nsd; 
  const int tnNo = com_mod.tnNo; 
  const int tDof = com_mod.tDof; 

  Vector<double> sA(tnNo); 
  Vector<double> sF(tnNo); 
  Array<double> xl(nsd,eNoN); 
  Array<double> yl(tDof,eNoN); 
  Array<double> dl(tDof,eNoN); 
  Vector<double> N(eNoN); 
  Array<double> Nx(nsd,eNoN);

  for (int e = 0; e < lM.nEl; e++) {
    if (lM.eType == ElementType::NRB) {
      //CALL NRBNNX(lM, e)
    }
    int cDmn = all_fun::domain(com_mod, lM, iEq, e);
    auto cPhys = eq.dmn[cDmn].phys;

    for (int a = 0; a < eNoN; a++) {
      int Ac = lM.IEN(a,e);
      for (int i = 0; i < nsd; i++) {
        xl(i,a) = com_mod.x(i,Ac);
      }
      for (int i = 0; i < tDof; i++) {
        yl(i,a) = lY(i,Ac);
        dl(i,a) = lD(i,Ac);
      }
    }

    Array<double> F;
    double divV = 0.0;
    double Jac = 0.0;
    Array<double> ksix(nsd,nsd);

    for (int g = 0; g < lM.nG; g++) {
      if (g == 0 || !lM.lShpF) {
        auto Nx_g = lM.Nx.slice(g);
        nn::gnn(eNoN, nsd, nsd, Nx_g, xl, Nx, Jac, ksix);
      }

      double w = lM.w(g)*Jac;
      N = lM.N.col(g);

      if ((cPhys == EquationType::phys_fluid) || (cPhys == EquationType::phys_CMM)) {
        Array<double> vx(nsd,nsd);

        if (nsd == 3) {
          for (int a = 0; a < eNoN; a++) {
            vx(0,0) = vx(0,0) + Nx(0,a)*yl(0,a);
            vx(1,1) = vx(1,1) + Nx(1,a)*yl(1,a);
            vx(2,2) = vx(2,2) + Nx(2,a)*yl(2,a);
          }
          divV = vx(0,0) + vx(1,1) + vx(2,2);

        } else { 
          for (int a = 0; a < eNoN; a++) {
            vx(0,0) = vx(0,0) + Nx(0,a)*yl(0,a);
            vx(1,1) = vx(1,1) + Nx(1,a)*yl(1,a);
          }
          divV = vx(0,0) + vx(1,1);
        }

      } else if (cPhys == EquationType::phys_ustruct) {
        Array<double> vx(nsd,nsd);
        auto F = mat_fun::mat_id(nsd);
        Vector<double> VxFi(nsd);

        if (nsd == 3) {
          for (int a = 0; a < eNoN; a++) {
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
          }

          auto Fi = mat_fun::mat_inv(F,3);

          VxFi(0) = vx(0,0)*Fi(0,0) + vx(0,1)*Fi(1,0) + vx(0,2)*Fi(2,0);
          VxFi(1) = vx(1,0)*Fi(0,1) + vx(1,1)*Fi(1,1) + vx(1,2)*Fi(2,1);
          VxFi(2) = vx(2,0)*Fi(0,2) + vx(2,1)*Fi(1,2) + vx(2,2)*Fi(2,2);
          divV = VxFi(0) + VxFi(1) + VxFi(2);

        } else { 
          for (int a = 0; a < eNoN; a++) {
            vx(0,0) = vx(0,0) + Nx(0,a)*yl(i,a);
            vx(0,1) = vx(0,1) + Nx(1,a)*yl(i,a);
            vx(1,0) = vx(1,0) + Nx(0,a)*yl(j,a);
            vx(1,1) = vx(1,1) + Nx(1,a)*yl(j,a);

            F(0,0) = F(0,0) + Nx(0,a)*dl(i,a);
            F(0,1) = F(0,1) + Nx(1,a)*dl(i,a);
            F(1,0) = F(1,0) + Nx(0,a)*dl(j,a);
            F(1,1) = F(1,1) + Nx(1,a)*dl(j,a);
          }

          auto Fi =  mat_fun::mat_inv(F,2);
          VxFi(0) = vx(0,0)*Fi(0,0) + vx(0,1)*Fi(1,0);
          VxFi(1) = vx(1,0)*Fi(0,1) + vx(1,1)*Fi(1,1);
          divV = VxFi(0) + VxFi(1);
        }
      }

      for (int a = 0; a < eNoN; a++) {
        int Ac = lM.IEN(a,e);
        sA(Ac) = sA(Ac) + w*N(a);
        sF(Ac) = sF(Ac) + w*N(a)*divV;
      }
    }
  }

  all_fun::commu(com_mod, sF);
  all_fun::commu(com_mod, sA);

  for (int a = 0; a < lM.nNo; a++) {
    int Ac = lM.gN(a);
    if (!utils::is_zero(sA(Ac))) {
      res(0,a) = res(0,a) + sF(Ac) / sA(Ac);
    } 
  }
}

/// @brief Routine for post processing fiber alignment
//
void fib_algn_post(Simulation* simulation, const mshType& lM, Array<double>& res, const Array<double>& lD, const int iEq)
{
  using namespace consts;

  auto& com_mod = simulation->com_mod;
  auto& cm = com_mod.cm;
  auto& cm_mod = simulation->cm_mod;
  auto& eq = com_mod.eq[iEq];

  int eNoN = lM.eNoN;
  int dof  = eq.dof;

  int i = eq.s;
  int j = i + 1;
  int k = j + 1;

  const int nsd = com_mod.nsd; 
  const int tnNo = com_mod.tnNo; 
  const int tDof = com_mod.tDof; 

  Vector<double> sA(tnNo); 
  Vector<double> sF(tnNo); 
  Array<double> xl(nsd,eNoN); 
  Array<double> dl(tDof,eNoN); 
  Array<double> fN(nsd,2); 
  Array<double> fl(nsd,2); 
  Array<double> Nx(nsd,eNoN); 
  Vector<double> N(eNoN);

  for (int e = 0; e < lM.nEl; e++) {
    int cDmn = all_fun::domain(com_mod, lM, iEq, e);
    auto cPhys = eq.dmn[cDmn].phys;
    if (cPhys != EquationType::phys_struct && cPhys != EquationType::phys_ustruct) {
      continue; 
    } 
    if (lM.eType == ElementType::NRB) {
      // CALL NRBNNX(lM, e)
    }

    for (int a = 0; a < eNoN; a++) {
      int Ac = lM.IEN(a,e);
      for (int i = 0; i < nsd; i++) {
        xl(i,a) = com_mod.x(i,Ac);
        dl(i,a) = lD(i,Ac);
      }
    }

    for (int i = 0; i < nsd; i++) {
      fN(i,0) = lM.fN(i,e);
      fN(i,0) = lM.fN(i+nsd,e);
    }

    Array<double> F;

    for (int g = 0; g < lM.nG; g++) {
      double Jac = 0.0;

      if (g == 0 || !lM.lShpF) {
        auto Nx = lM.Nx.slice(g);
        nn::gnn(eNoN, nsd, nsd, Nx, xl, Nx, Jac, F);
      }

      double w = lM.w(g)*Jac;
      auto F = mat_fun::mat_id(nsd); 

      for (int a = 0; a < eNoN; a++) {
        if (nsd == 3) {
          F(0,0) = F(0,0) + Nx(0,a)*dl(i,a);
          F(0,1) = F(0,1) + Nx(1,a)*dl(i,a);
          F(0,2) = F(0,2) + Nx(2,a)*dl(i,a);
          F(1,0) = F(1,0) + Nx(0,a)*dl(j,a);
          F(1,1) = F(1,1) + Nx(1,a)*dl(j,a);
          F(1,2) = F(1,2) + Nx(2,a)*dl(j,a);
          F(2,0) = F(2,0) + Nx(0,a)*dl(k,a);
          F(2,1) = F(2,1) + Nx(1,a)*dl(k,a);
          F(2,2) = F(2,2) + Nx(2,a)*dl(k,a);
        } else {
          F(0,0) = F(0,0) + Nx(0,a)*dl(i,a);
          F(0,1) = F(0,1) + Nx(1,a)*dl(i,a);
          F(1,0) = F(1,0) + Nx(0,a)*dl(j,a);
          F(1,1) = F(1,1) + Nx(1,a)*dl(j,a);
        }
      }
      for (int iFn = 0; iFn < 2; iFn++) {
        for (int i = 0; i < nsd; i++) {
          auto fN_col = fN.col(iFn);
          auto F_fN = mat_fun::mat_mul(F, fN_col);
          fl.set_col(iFn, F_fN / sqrt(utils::norm(F_fN)));
        }
      }

      double sHat = utils::norm(fl.col(0), fl.col(1));

      for (int a = 0; a < eNoN; a++) {
        int Ac = lM.IEN(a,e);
        sA(Ac) = sA(Ac) + w*N(a);
        sF(Ac) = sF(Ac) + w*N(a)*sHat;
      }
    }
  }

  all_fun::commu(com_mod, sF);
  all_fun::commu(com_mod, sA);

  for (int a = 0; a < lM.nNo; a++) {
    int Ac = lM.gN(a);
    if (!utils::is_zero(sA(Ac))) {
      res(0,a) = res(0,a) + sF(Ac) / sA(Ac);
    }
  }
}

/// @brief Routine for post processing fiber directions.
//
void fib_dir_post(Simulation* simulation, const mshType& lM, const int nFn, Array<double>& res, const Array<double>& lD, const int iEq)
{
  using namespace consts;

  auto& com_mod = simulation->com_mod;
  auto& cm = com_mod.cm;
  auto& cm_mod = simulation->cm_mod;
  auto& eq = com_mod.eq[iEq];

  // [NOTE] Setting gobal variable 'dof'.
  com_mod.dof = eq.dof;

  int eNoN = lM.eNoN;
  int i = eq.s;
  int j = i + 1;
  int k = j + 1;

  const int nsd = com_mod.nsd; 
  const int tnNo = com_mod.tnNo; 
  const int tDof = com_mod.tDof; 

  Vector<double> sA(tnNo); 
  Array<double> sF(nFn*nsd,tnNo); 
  Array<double> xl(nsd,eNoN); 
  Array<double> dl(tDof,eNoN); 
  Array<double> fN(nsd,lM.nFn); 
  Array<double> fl(nsd,lM.nFn); 
  Array<double> Nx(nsd,eNoN); 
  Vector<double> N(eNoN);

  for (int e = 0; e < lM.nEl; e++) {
    int cDmn = all_fun::domain(com_mod, lM, iEq, e);
    auto cPhys = eq.dmn[cDmn].phys;
    if (lM.eType == ElementType::NRB) {
      // CALL NRBNNX(lM, e)
    }

    for (int a = 0; a < eNoN; a++) {
      int Ac = lM.IEN(a,e);
      for (int i = 0; i < nsd; i++) {
        xl(i,a) = com_mod.x(i,Ac);
      }
      for (int i = 0; i < tDof; i++) {
        dl(i,a) = lD(i,Ac);
      }
    }

    for (int iFn = 0; iFn < lM.nFn; iFn++) {
      for (int i = 0; i < nsd; i++) {
        fN(i,iFn) = lM.fN(i+iFn*nsd,e);
      }
    }

    Array<double> F(nsd,nsd);
    Array<double> Nx(nsd,eNoN);
    double Jac = 0.0;

    for (int g = 0; g < lM.nG; g++) {
      if (g == 0 || !lM.lShpF) {
        auto Nxi = lM.Nx.slice(g);
        nn::gnn(eNoN, nsd, nsd, Nxi, xl, Nx, Jac, F);
      }

      double w = lM.w(g) * Jac;
      N = lM.N.col(g);
      F = mat_fun::mat_id(nsd); 

      for (int a = 0; a < eNoN; a++) {
        if (nsd == 3) {
          F(0,0) = F(0,0) + Nx(0,a)*dl(i,a);
          F(0,1) = F(0,1) + Nx(1,a)*dl(i,a);
          F(0,2) = F(0,2) + Nx(2,a)*dl(i,a);
          F(1,0) = F(1,0) + Nx(0,a)*dl(j,a);
          F(1,1) = F(1,1) + Nx(1,a)*dl(j,a);
          F(1,2) = F(1,2) + Nx(2,a)*dl(j,a);
          F(2,0) = F(2,0) + Nx(0,a)*dl(k,a);
          F(2,1) = F(2,1) + Nx(1,a)*dl(k,a);
          F(2,2) = F(2,2) + Nx(2,a)*dl(k,a);
        } else {
          F(0,0) = F(0,0) + Nx(0,a)*dl(i,a);
          F(0,1) = F(0,1) + Nx(1,a)*dl(i,a);
          F(1,0) = F(1,0) + Nx(0,a)*dl(j,a);
          F(1,1) = F(1,1) + Nx(1,a)*dl(j,a);
        }
      }

      for (int iFn = 0; iFn < lM.nFn; iFn++) {
        for (int i = 0; i < nsd; i++) {
          auto fN_col = fN.col(iFn);
          auto F_fN = mat_fun::mat_mul(F, fN_col);
          fl.set_col(iFn, F_fN / sqrt(utils::norm(F_fN)));
        }
      }

      for (int a = 0; a < eNoN; a++) {
        int Ac = lM.IEN(a,e);
        sA(Ac) = sA(Ac) + w*N(a);
        for (int iFn = 0; iFn < lM.nFn; iFn++) {
          int b = iFn*nsd;
          for (int l = 0; l < nsd; l++) {
            sF(b+l,Ac) = sF(b+l,Ac) + w*N(a)*fl(l,iFn);
          }
        }
      }
    }
  }

  all_fun::commu(com_mod, sF);
  all_fun::commu(com_mod, sA);

  for (int a = 0; a < lM.nNo; a++) {
    int Ac = lM.gN(a);
    if (!utils::is_zero(sA(Ac))) {
      for (int i = 0; i < nFn*nsd; i++) {
        res(i,a) = res(i,a) + sF(i,Ac) / sA(Ac);
      }
    } 
  }

}

/// @brief Compute fiber stretch based on 4th invariant: I_{4,f}
//
void fib_strech(Simulation* simulation, const int iEq, const mshType& lM, const Array<double>& lD, Vector<double>& res)
{
  using namespace consts;

  auto& com_mod = simulation->com_mod;
  auto& cm = com_mod.cm;
  auto& cm_mod = simulation->cm_mod;
  auto& eq = com_mod.eq[iEq];

  int nsd = com_mod.nsd;
  int tnNo = com_mod.tnNo;
  int tDof = com_mod.tDof;

  // [NOTE] Setting gobal variable 'dof'.
  com_mod.dof = eq.dof;

  int eNoN = lM.eNoN;
  int i = eq.s;
  int j = i + 1;
  int k = j + 1;

  Vector<double> sA(tnNo); 
  Vector<double> sF(tnNo); 
  Array<double> xl(nsd,eNoN); 
  Array<double> dl(tDof,eNoN);
  Array<double> Nx(nsd,eNoN); 
  Vector<double> N(eNoN);

  for (int e = 0; e < lM.nEl; e++) {
    int cDmn  = all_fun::domain(com_mod, lM, iEq, e);
    auto cPhys = eq.dmn[cDmn].phys;
    if (lM.eType == ElementType::NRB) {
      //CALL NRBNNX(lM, e)
    }

    for (int a = 0; a < eNoN; a++) { 
      int Ac = lM.IEN(a,e);
      xl.set_col(a, com_mod.x.col(Ac));
      dl.set_col(a, lD.col(Ac));
    }

    for (int g = 0; g < lM.nG; g++) {
      double Jac = 0.0;
      Array<double> F(nsd,nsd);
      if (g == 0 || !lM.lShpF) {
        auto Nx = lM.Nx.slice(g);
        nn::gnn(eNoN, nsd, nsd, Nx, xl, Nx, Jac, F);
      }

      double w = lM.w(g)*Jac;
      auto N = lM.N.col(g);
      F = mat_fun::mat_id(nsd);

      for (int a = 0; a < eNoN; a++) { 
        if (nsd == 3) {
          F(0,0) = F(0,0) + Nx(0,a)*dl(i,a);
          F(0,1) = F(0,1) + Nx(1,a)*dl(i,a);
          F(0,2) = F(0,2) + Nx(2,a)*dl(i,a);
          F(1,0) = F(1,0) + Nx(0,a)*dl(j,a);
          F(1,1) = F(1,1) + Nx(1,a)*dl(j,a);
          F(1,2) = F(1,2) + Nx(2,a)*dl(j,a);
          F(2,0) = F(2,0) + Nx(0,a)*dl(k,a);
          F(2,1) = F(2,1) + Nx(1,a)*dl(k,a);
          F(2,2) = F(2,2) + Nx(2,a)*dl(k,a);
        } else {
          F(0,0) = F(0,0) + Nx(0,a)*dl(i,a);
          F(0,1) = F(0,1) + Nx(1,a)*dl(i,a);
          F(1,0) = F(1,0) + Nx(0,a)*dl(j,a);
          F(1,1) = F(1,1) + Nx(1,a)*dl(j,a);
        }
      }

      auto fl = mat_fun::mat_mul(F, lM.fN.rows(0,nsd-1,e));
      double I4f = utils::norm(fl);

      for (int a = 0; a < eNoN; a++) { 
        int Ac = lM.IEN(a,e);
        sA(Ac) = sA(Ac) + w*N(a);
        sF(Ac) = sF(Ac) + w*N(a)*I4f;
      }
    }
  }

  all_fun::commu(com_mod, sF);
  all_fun::commu(com_mod, sA);

  res = 0.0;

  for (int a = 0; a < lM.nNo; a++) {
    int Ac = lM.gN(a);
    if (!utils::is_zero(sA(Ac))) {
      res(a) = res(a) + sF(Ac) / sA(Ac);
    }
  }

}

void post(Simulation* simulation, const mshType& lM, Array<double>& res, const Array<double>& lY, const Array<double>& lD, 
    consts::OutputNameType outGrp, const int iEq)
{
  using namespace consts;
  auto& com_mod = simulation->com_mod;
  auto& cm = com_mod.cm;
  auto& cm_mod = simulation->cm_mod;

  #define n_debug_post
  #ifdef debug_post
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "outGrp: " << outGrp;
  dmsg << "iEq: " << iEq;
  #endif

  bool FSIeq = false;
  auto& eq = com_mod.eq[iEq];

  if (eq.phys == EquationType::phys_FSI) {
    FSIeq = true;
  }

  // Since energy flux can be directly calculated from nodal values.
  // Note that if there are more than one domain, we need to rely on
  // element based calculations
  //
  int nsd = com_mod.nsd;

  if ((outGrp == OutputNameType::outGrp_eFlx) && (com_mod.dmnId.size() == 0)) {
    double rho = eq.dmn[0].prop[PhysicalProperyType::fluid_density];
    for (int a = 0; a < lM.nNo; a++) {
      int Ac = lM.gN(a);
      double p  = lY(nsd,Ac);

      auto u  = lY.col(Ac, {0,nsd-1});
      double unorm = utils::norm(u);
      for (int i = 0; i < nsd; i++) {
        res(i,Ac) = (p + 0.5 * rho * unorm) * u(i);
      }
     }
     return; 
  }

  // Other outputs require more calculations
  //
  int eNoN  = lM.eNoN;
  int tnNo = com_mod.tnNo;
  int tDof = com_mod.tDof;
  Vector<double> sA(tnNo); 
  Array<double> sF(maxNSD,tnNo); 
  Array<double> xl(nsd,eNoN); 
  Array<double> yl(tDof,eNoN); 
  Array<double> Nx(nsd,eNoN); 
  Vector<double> N(eNoN);

  int insd = nsd;
  if (lM.lFib) {
    insd = 1;
  }

  for (int e = 0; e < lM.nEl; e++) {
    int cDmn = all_fun::domain(com_mod, lM, iEq, e);
    if (cDmn == -1) {
      continue;
    } 
    if (lM.eType == ElementType::NRB) {
      // CALL NRBNNX(lM, e)
    }

    // Finding the norm for all the nodes of this element, including
    // those that don't belong to this face, which will be interpolated
    // from the nodes of the face
    //
    for (int a = 0; a < eNoN; a++) {
      int Ac = lM.IEN(a,e);
      for (int i = 0; i < nsd; i++) {
        xl(i,a) = com_mod.x(i,Ac);
      }
      for (int i = 0; i < tDof; i++) {
        yl(i,a) = lY(i,Ac);
      }

      if (FSIeq) {
        for (int i = 0; i < nsd; i++) {
          xl(i,a) = xl(i,a) + lD(i+nsd+1,Ac);
          yl(i,a) = yl(i,a) - lY(i+nsd+1,Ac);
        }
      }
    }

    Array<double> ksix(nsd,nsd);
    double Jac = 0.0;

    for (int g = 0; g < lM.nG; g++) {
      if (g == 0 || !lM.lShpF) {
        auto Nx_g = lM.Nx.slice(g);
        nn::gnn(eNoN, nsd, insd, Nx_g, xl, Nx, Jac, ksix);
      }

      double w = lM.w(g) * Jac;
      auto N = lM.N.col(g);
      Vector<double> lRes(maxNSD);

      // Vorticity calculation 
      //
      if (outGrp == OutputNameType::outGrp_vort) {
        for (int a = 0; a < eNoN; a++) {
          if (nsd == 2) {
            lRes(2) = lRes(2) + Nx(0,a)*yl(1,a)- Nx(1,a)*yl(0,a);
          } else {
            lRes(0) = lRes(0)+ Nx(1,a)*yl(2,a)- Nx(2,a)*yl(1,a);
            lRes(1) = lRes(1)+ Nx(2,a)*yl(0,a)- Nx(0,a)*yl(2,a);
            lRes(2) = lRes(2)+ Nx(0,a)*yl(1,a)- Nx(1,a)*yl(0,a);
          }
        }

      // Vortex Identification Criterion (lamda_ci)
      } else if (outGrp == OutputNameType::outGrp_vortex) {
        Array<double> ux(nsd,nsd);
        for (int a = 0; a < eNoN; a++) {
          for (int i = 0; i < nsd; i++) {
            for (int j = 0; j < nsd; j++) {
              ux(i,j) = ux(i,j) + Nx(i,a)*yl(j,a);
            }
          }
        }

        // [NOTE] Not implemented.
        // eig = MAT_EIG(ux, nsd)
        // lRes(1:nsd) = AIMAG(eig)
        // lRes(1) = MAXVAL(lRes(1:nsd))
        // lRes(2:maxnsd) = 0.

      //  Energy flux calculation   
      //
      } else if (outGrp == OutputNameType::outGrp_eFlx) {
        double rho = eq.dmn[cDmn].prop[PhysicalProperyType::fluid_density];
        double p = 0.0;
        Vector<double> u(nsd);
        Vector<double> lRes(maxNSD);

        for (int a = 0; a < eNoN; a++) {
          p = p + N(a)*yl(nsd,a);
          for (int i = 0; i < nsd; i++) {
            u(i) = u(i) + N(a)*yl(i,a);     
          }
        }
        double unorm = utils::norm(u);

        for (int i = 0; i < nsd; i++) {
          lRes(i) = (p + 0.5 * rho * unorm) * u(i);
        }

      // Heat flux calculation   
      //
      } else if (outGrp == OutputNameType::outGrp_hFlx) {
        double kappa = eq.dmn[cDmn].prop[PhysicalProperyType::conductivity];
        int i = eq.s;

        if (eq.phys == EquationType::phys_heatF) {
          Vector<double> u(nsd);
          Vector<double> T(nsd);
          Vector<double> q(nsd);

          for (int a = 0; a < eNoN; a++) {
            for (int j = 0; j < nsd; j++) {
              q(j) = q(j) + Nx(j,a)*yl(i,a);
              u(j) = u(j) + N(a)*yl(j,a);
              T(j) = T(j) + N(a)*yl(i,a);
            }
          }
          for (int j = 0; j < nsd; j++) {
            lRes(j) = u(j)*T(j) - kappa*q(j);
          }

        } else {
          Vector<double> q(nsd);
          for (int a = 0; a < eNoN; a++) {
            for (int j = 0; j < nsd; j++) {
              q(j) = q(j) + Nx(j,a)*yl(i,a);
            }
          }
          for (int j = 0; j < nsd; j++) {
            lRes(j) = -kappa * q(j);
          }
        }

      // Strain tensor invariants calculation   
      //
      } else if (outGrp == OutputNameType::outGrp_stInv) {
        Array<double> ksix(nsd,nsd);
        Vector<double> lRes(maxNSD);

        for (int a = 0; a < eNoN; a++) {
          ksix(0,0) = ksix(0,0) + Nx(0,a)*yl(0,a);
          ksix(1,0) = ksix(1,0) +(Nx(0,a)*yl(1,a) + Nx(1,a)*yl(0,a))*0.5;
          ksix(1,1) = ksix(1,1) + Nx(1,a)*yl(1,a);

          if (nsd == 3) {
            ksix(2,0) = ksix(2,0) +(Nx(0,a)*yl(2,a) + Nx(2,a)*yl(0,a))*0.5;
            ksix(2,1) = ksix(2,1) +(Nx(1,a)*yl(2,a) + Nx(2,a)*yl(1,a))*0.5;
            ksix(2,2) = ksix(2,2) + Nx(1,a)*yl(1,a);
          }
        }

        if (nsd == 2) {
          lRes(0) = ksix(0,0) + ksix(1,1);
          lRes(1) = ksix(0,0)*ksix(1,1) - ksix(1,0)*ksix(1,0);
        } else {
          lRes(0) = ksix(0,0) + ksix(1,1) + ksix(2,2);
          lRes(1) = ksix(0,0)*ksix(1,1) + ksix(1,1)*ksix(2,2) + ksix(2,2)*ksix(0,0) - ksix(1,0)*ksix(1,0) - 
              ksix(2,0)*ksix(2,0) - ksix(2,1)*ksix(2,1);
          lRes(2) = ksix(0,0)*ksix(1,1)*ksix(2,2) + ksix(1,0)*ksix(2,1)*ksix(2,0)*2.0 - ksix(0,0)*ksix(2,1)*ksix(2,1) - 
              ksix(2,0)*ksix(1,1)*ksix(2,0) - ksix(1,0)*ksix(1,0)*ksix(2,2);
        }
        for (int j = 0; j < maxNSD; j++) {
          lRes(j) = fabs(lRes(j));
        }

      // Viscosity
      //
      } else if (outGrp == OutputNameType::outGrp_Visc) {
        Array<double> ux(nsd,nsd);
        Vector<double> lRes(maxNSD);

        for (int a = 0; a < eNoN; a++) {
          for (int i = 0; i < nsd; i++) {
            for (int j = 0; j < nsd; j++) {
              ux(i,j) = ux(i,j) + Nx(i,a)*yl(j,a);
            }
          }
        }

        // Shear rate, gam := (2*e_ij*e_ij)^0.5
        double gam = 0.0;
        for (int i = 0; i < nsd; i++) {
          for (int j = 0; j < nsd; j++) {
            gam = gam + (ux(i,j)+ux(j,i))*(ux(i,j)+ux(j,i));
          }
        }

        double mu = 0.0;
        double mu_s = 0.0;
        gam = sqrt(0.5*gam);
        // Compute viscosity
        fluid::get_viscosity(com_mod, eq.dmn[cDmn], gam, mu, mu_s, mu_s);
        lRes(0) = mu;
      } else {
        throw std::runtime_error("Error in the post() function.");
      }

      // Mapping Tau into the nodes by assembling it into a local vector
      for (int a = 0; a < eNoN; a++) {
        int Ac = lM.IEN(a,e);
        sA(Ac) = sA(Ac) + w*N(a);
        for (int i = 0; i < maxNSD; i++) {
          sF(i,Ac) = sF(i,Ac) + w*N(a)*lRes(i);
        }
      }
    }
  }

  all_fun::commu(com_mod, sF);
  all_fun::commu(com_mod, sA);

  for (int a = 0; a < lM.nNo; a++) {
    int Ac = lM.gN(a);
    for (int i = 0; i < maxNSD; i++) {
      res(i,a) = sF(i,Ac) / sA(Ac);
    }
  }
}

/// @brief Postprocessing function - used to convert restart bin files into VTK format.
///
/// Reproducces 'SUBROUTINE PPBIN2VTK()' defined in POST.f.
///
/// \todo [NOTE] This is not fully implemeted and is not tested, there are
/// no tests in 'svFSI-Tests' for this.
//
void ppbin2vtk(Simulation* simulation)
{
  auto& com_mod = simulation->com_mod;
  auto& cm = com_mod.cm;
  auto& cm_mod = simulation->cm_mod;

  int stFileIncr = com_mod.stFileIncr;
  auto saveName = com_mod.saveName;
  auto stFileName = com_mod.stFileName;

  for (int iTS = 0; iTS < com_mod.nTS; iTS++) {
    if (iTS % stFileIncr == 0) {
      auto stmp = std::to_string(iTS);
      std::string fName;
      bool flag = false;
      //WRITE(stmp,'(I3.3)') iTS
      //if (iTS .GE. 1000) stmp = STR(iTS)

      // Ignore if vtu file already exists
      if (cm.mas(cm_mod)) {
        fName = saveName + "_" + stmp + ".vtu";
        if (FILE *file = fopen(fName.c_str(), "r")) {
          fclose(file);
          flag = true;
        } else {
          flag = false;
        }
        //fName = TRIM(saveName)//"_"//TRIM(ADJUSTL(stmp))//".vtu"
        //INQUIRE(FILE=TRIM(fName), EXIST=flag)
      }

      cm.bcast(cm_mod, &flag);

      if (flag) {
        continue;
      }

      // Ignore if bin file does not exist
      fName = stFileName + "_" + stmp + ".bin";
      if (FILE *file = fopen(fName.c_str(), "r")) {
        fclose(file);
        flag = true;
      } else {
        flag = false;
      }

      if (!flag) {
        continue;
      }

      std::array<double,3> rtmp;
      init_from_bin(simulation, fName, rtmp);

      bool lAve = false;

      vtk_xml::write_vtus(simulation, com_mod.Ao, com_mod.Yo, com_mod.Do, lAve);
    }
  }

  finalize(simulation);
  
  MPI_Finalize();
}

//----------
// shl_post
//----------
// Routine for post processing shell-based quantities
//
// Reproduces Fortran SHLPOST.
//
void shl_post(Simulation* simulation, const mshType& lM, const int m, Array<double>& res, 
    Vector<double>& resE, const Array<double>& lD, const int iEq, consts::OutputNameType outGrp)
{
  using namespace consts;
  using namespace mat_fun;
  using namespace utils;

  #define n_debug_shl_post
  #ifdef debug_shl_post
  DebugMsg dmsg(__func__, 0);
  dmsg.banner();
  #endif

  auto& com_mod = simulation->com_mod;
  auto& cm = com_mod.cm;
  auto& cm_mod = simulation->cm_mod;
  auto& eq = com_mod.eq[iEq];

  const int nsd = com_mod.nsd;
  const int tnNo = com_mod.tnNo;
  const int tDof = com_mod.tDof;

  // [NOTE] Setting gobal variable 'dof'.
  com_mod.dof = eq.dof;

  int i = eq.s;
  int j = i + 1;
  int k = j + 1;

  int nFn = lM.nFn;
  if (nFn == 0) {
    nFn = 1;
  }

  // Set shell dimension := 2
  int insd = nsd-1;

  // Initialize tensor operations
  ten_init(insd);

  // Set eNoN (number of nodes per element)
  //
  int eNoN = lM.eNoN;

  if (lM.eType == ElementType::TRI3) {
    eNoN = 2*eNoN;
  }

  Vector<double> sA(tnNo), sE(lM.nEl), resl(m), N(lM.eNoN);
  Vector<int> ptr(eNoN);
  Array<double> sF(m,tnNo), dl(tDof,eNoN), x0(3,eNoN), xc(3,eNoN), 
      fN(3,nFn), fNa0(2,eNoN), Nx(2,lM.eNoN);
  Array3<double> Bb(3,3,6);

  // Initialize arrays
  sA = 0.0;
  sF = 0.0;
  sE = 0.0;
  Bb = 0.0;

  // Compute quantities at the element level and project them to nodes
  //
  for (int e = 0; e < lM.nEl; e++) {
    int cDmn = all_fun::domain(com_mod, lM, iEq, e);
    auto cPhys = eq.dmn[cDmn].phys;
    //dmsg << "========== e: " << e+1;

    if (cPhys != EquationType::phys_shell) {
      continue;
    }
    //if (lM.eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

    // Get shell properties
    double nu = eq.dmn[cDmn].prop.at(PhysicalProperyType::poisson_ratio);
    double ht = eq.dmn[cDmn].prop.at(PhysicalProperyType::shell_thickness);

    // Check for incompressibility
    //
    bool incompFlag = false;
    if (is_zero(nu-0.50)) {
      incompFlag = true;
    }
    //dmsg << "incompFlag: " << incompFlag;

    // Get the reference configuration and displacement field
    x0 = 0.0;
    dl = 0.0;

    for (int a = 0; a < eNoN; a++) {
      int Ac;
      if (a < lM.eNoN) {
        Ac = lM.IEN(a,e);
        ptr(a) = Ac;
      } else {
        int b  = a - lM.eNoN;
        Ac = lM.eIEN(b,e);
        ptr(a) = Ac;
        if (Ac == -1) {
          continue; 
        }
      }

      for (int i = 0; i < 3; i++) {
        x0(i,a) = com_mod.x(i,Ac);
      }

      for (int i = 0; i < tDof; i++) {
        dl(i,a) = lD(i,Ac);
      }
    }

    // Get the current configuration
    xc = 0.0;

    for (int a = 0; a < eNoN; a++) {
      xc(0,a) = x0(0,a) + dl(i,a);
      xc(1,a) = x0(1,a) + dl(j,a);
      xc(2,a) = x0(2,a) + dl(k,a);
    }

    // Get fiber directions
    //
    fN = 0.0;

    if (lM.fN.size() != 0) {
      for (int iFn = 0; iFn < nFn; iFn++) {
        for (int i = 0; i < 3; i++) {
          fN(i,iFn) = lM.fN(i+nsd*iFn,e);
        }
      }
    }
    
    // Set number of integration points.
    // Note: Gauss integration performed for NURBS elements
    // Not required for constant-strain triangle elements
    //
    int nwg;

    if (lM.eType == ElementType::TRI3) {
      nwg = 1;
    } else {
      nwg = lM.nG;
    }

    // Update shapefunctions for NURBS elements
    //
    //if (lM.eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

    double Je = 0.0;
    resl = 0.0;
    Vector<double> N;
    Array<double> Nxx, Nx;

    for (int g = 0; g < nwg; g++) {
      double w = 0.0;
      double lam3 = 0.0;
      double aa_0[2][2]{}, bb_0[2][2]{};
      double aa_x[2][2]{}, bb_x[2][2]{};

      Array<double> aCov0(3,2), aCnv0(3,2), aCov(3,2), aCnv(3,2);
      Vector<double> nV0(3), nV(3);
      //dmsg << "---------- g: " << g;

      // [TODO] This is not fully implemented
      //
      if (lM.eType != ElementType::TRI3) {
        // Set element shape functions and their derivatives
        if (lM.eType ==  ElementType::NRB) {
          //N   = lM.N(:,g)
          //Nx  = lM.Nx(:,:,g)
          //Nxx = lM.Nxx(:,:,g)
        } else {
          N = lM.fs[0].N.rcol(g);
          Nx = lM.fs[0].Nx.rslice(g);
          Nxx = lM.fs[0].Nxx.rslice(g);
        }

        // Covariant and contravariant bases (ref. config.)
        //
        nn::gnns(nsd, eNoN, Nx, x0, nV0, aCov0, aCnv0);
        auto Jac0 = sqrt(norm(nV0));
        nV0 = nV0 / Jac0;

        // Covariant and contravariant bases (spatial config.)
        nn::gnns(nsd, eNoN, Nx, xc, nV, aCov, aCnv);
        auto Jac = sqrt(norm(nV));
        nV = nV/Jac;

        // Second derivatives for curvature coeffs. (ref. config)
#if  0
        r0_xx(:,:,:) = 0.0
        r_xx(:,:,:)  = 0.0

        for (int a = 0; a < eNoN; a++) {
          r0_xx(0,0,:) = r0_xx(0,0,:) + Nxx(0,a)*x0(:,a)
          r0_xx(1,1,:) = r0_xx(1,1,:) + Nxx(1,a)*x0(:,a)
          r0_xx(0,1,:) = r0_xx(0,1,:) + Nxx(2,a)*x0(:,a)

          r_xx(0,0,:) = r_xx(0,0,:) + Nxx(0,a)*xc(:,a)
          r_xx(1,1,:) = r_xx(1,1,:) + Nxx(1,a)*xc(:,a)
          r_xx(0,1,:) = r_xx(0,1,:) + Nxx(2,a)*xc(:,a)
        }

        r0_xx(1,0,:) = r0_xx(0,1,:)
        r_xx(1,0,:)  = r_xx(0,1,:)

        // Compute metric tensor (aa) and curvature coefficients(bb)
        //
        double aa_0[2][2]{}, bb_0[2][2]{};
        double aa_x[2][2]{}, bb_x[2][2]{};

        for (int l = 0; l < nsd; l++) {
          aa_0(0,0) = aa_0(0,0) + aCov0(l,0)*aCov0(l,0)
          aa_0(0,1) = aa_0(0,1) + aCov0(l,0)*aCov0(l,1)
          aa_0(1,0) = aa_0(1,0) + aCov0(l,1)*aCov0(l,0)
          aa_0(1,1) = aa_0(1,1) + aCov0(l,1)*aCov0(l,1)

          aa_x(0,0) = aa_x(0,0) + aCov(l,0)*aCov(l,0)
          aa_x(0,1) = aa_x(0,1) + aCov(l,0)*aCov(l,1)
          aa_x(1,0) = aa_x(1,0) + aCov(l,1)*aCov(l,0)
          aa_x(1,1) = aa_x(1,1) + aCov(l,1)*aCov(l,1)

          bb_0(0,0) = bb_0(0,0) + r0_xx(0,0,l)*nV0(l)
          bb_0(0,1) = bb_0(0,1) + r0_xx(0,1,l)*nV0(l)
          bb_0(1,0) = bb_0(1,0) + r0_xx(1,0,l)*nV0(l)
          bb_0(1,1) = bb_0(1,1) + r0_xx(1,1,l)*nV0(l)

          bb_x(0,0) = bb_x(0,0) + r_xx(0,0,l)*nV(l)
          bb_x(0,1) = bb_x(0,1) + r_xx(0,1,l)*nV(l)
          bb_x(1,0) = bb_x(1,0) + r_xx(1,0,l)*nV(l)
          bb_x(1,1) = bb_x(1,1) + r_xx(1,1,l)*nV(l)
        }

        //  Set weight of the Gauss point
        w  = lM.w(g)*Jac0

#endif

      // for constant strain triangles
      } else {      

        // Set element shape functions and their derivatives
        N = lM.N.rcol(g);
        Nx = lM.Nx.rslice(0);

        // Covariant and contravariant bases (ref. config.)
        //
        Array<double> tmpX(nsd,lM.eNoN);

        for (int i = 0; i < nsd; i++) {
          for (int j = 0; j < lM.eNoN; j++) {
            tmpX(i,j) = x0(i,j);
          }
        }

        nn::gnns(nsd, lM.eNoN, Nx, tmpX, nV0, aCov0, aCnv0);
        auto Jac0 = sqrt(norm(nV0));
        nV0  = nV0 / Jac0;

        // Covariant and contravariant bases (spatial config.)
        //
        for (int i = 0; i < nsd; i++) {
          for (int j = 0; j < lM.eNoN; j++) {
            tmpX(i,j) = xc(i,j);
          }
        }

        nn::gnns(nsd, lM.eNoN, Nx, tmpX, nV, aCov, aCnv);
        auto Jac = sqrt(norm(nV));
        nV = nV / Jac;

        // Compute metric tensor (aa)
        //
        for (int l = 0; l < nsd; l++) {
          aa_0[0][0] = aa_0[0][0] + aCov0(l,0)*aCov0(l,0);
          aa_0[0][1] = aa_0[0][1] + aCov0(l,0)*aCov0(l,1);
          aa_0[1][0] = aa_0[1][0] + aCov0(l,1)*aCov0(l,0);
          aa_0[1][1] = aa_0[1][1] + aCov0(l,1)*aCov0(l,1);

          aa_x[0][0] = aa_x[0][0] + aCov(l,0)*aCov(l,0);
          aa_x[0][1] = aa_x[0][1] + aCov(l,0)*aCov(l,1);
          aa_x[1][0] = aa_x[1][0] + aCov(l,1)*aCov(l,0);
          aa_x[1][1] = aa_x[1][1] + aCov(l,1)*aCov(l,1);
        }

        shells::shell_bend_cst(com_mod, lM, e, ptr, x0, xc, bb_0, bb_x, Bb, false);

        // Set weight of the Gauss point
        w = Jac0*0.50;
      }

      // Compute fiber direction in curvature coordinates
      //
      fNa0 = 0.0;

      for (int iFn = 0; iFn < nFn; iFn++) {
        for (int l = 0; l < nsd; l++) {
          fNa0(0,iFn) = fNa0(0,iFn) + fN(l,iFn)*aCnv0(l,0);
          fNa0(1,iFn) = fNa0(1,iFn) + fN(l,iFn)*aCnv0(l,1);
        }
      }

      // Compute stress resultants and lambda3 (integrated through
      //       the shell thickness)
      //
      Array<double> Sm(3,2);       
      Array3<double> Dm(3,3,3);
      shells::shl_strs_res(com_mod, eq.dmn[cDmn], nFn, fNa0, aa_0, aa_x, bb_0, bb_x, lam3, Sm, Dm);

      // Shell in-plane deformation gradient tensor
      //
      auto F = mat_dyad_prod(aCov.rcol(0), aCnv0.rcol(0), 3) + 
               mat_dyad_prod(aCov.rcol(1), aCnv0.rcol(1), 3);

      // D deformation gradient tensor in shell continuum
      auto F3d = F + lam3*mat_dyad_prod(nV, nV0, 3);
      auto detF = mat_det(F3d, nsd);
      Je = Je + w;
      auto Im = mat_id(nsd);

      switch (outGrp) {
        //dmsg << "outGrp: " << outGrp;
        case OutputNameType::outGrp_J: {
          // Jacobian := determinant of deformation gradient tensor
          resl(0) = detF;
          sE(e) = sE(e) + w*detF;
        } break;

        case OutputNameType::outGrp_F:
          // 3D deformation gradient tensor (F)
          resl(0) = F3d(0,0);
          resl(1) = F3d(0,1);
          resl(2) = F3d(0,2);
          resl(3) = F3d(1,0);
          resl(4) = F3d(1,1);
          resl(5) = F3d(1,2);
          resl(6) = F3d(2,0);
          resl(7) = F3d(2,1);
          resl(8) = F3d(2,2);
        break; 

        case OutputNameType::outGrp_strain:
        case OutputNameType::outGrp_C:
        case OutputNameType::outGrp_I1: {
          // In-plane Cauchy-Green deformation tensor
          auto C = mat_mul(transpose(F), F);

          // In-plane Green-Lagrange strain tensor
          auto Eg = 0.50 * (C - Im);

          if (outGrp == OutputNameType::outGrp_strain) {
            // resl is used to remap Eg
            resl(0) = Eg(0,0);
            resl(1) = Eg(1,1);
            resl(2) = Eg(2,2);
            resl(3) = Eg(0,1);
            resl(4) = Eg(1,2);
            resl(5) = Eg(2,0);

          } else if (outGrp == OutputNameType::outGrp_C) {
            // resl is used to remap C
            resl(0) = C(0,0);
            resl(1) = C(1,1);
            resl(2) = C(2,2);
            resl(3) = C(0,1);
            resl(4) = C(1,2);
            resl(5) = C(2,0);

          } else if (outGrp == OutputNameType::outGrp_I1) {
            resl(0) = mat_trace(C, 3);
            sE(e) = sE(e) + w*resl(0);
          }
        } break;

        case OutputNameType::outGrp_stress: {
          //dmsg << "outGrp: " << " outGrp_stress";
          Array<double> S(3,3);

          // 2nd Piola-Kirchhoff stress
          S(0,0) = Sm(0,0);
          S(1,1) = Sm(1,0);
          S(0,1) = Sm(2,0);
          S(1,0) = S(0,1);

          //  Normalizing stress by thickness
          S = S / ht;

          // 2nd Piola-Kirchhoff stress tensor
          resl(0) = S(0,0);
          resl(1) = S(1,1);
          resl(2) = S(2,2);
          resl(3) = S(0,1);
          resl(4) = S(1,2);
          resl(5) = S(2,0);
          //dmsg << "resl: " << resl;
        } break;
      }

      for (int a = 0; a < lM.eNoN; a++) {
        int Ac = lM.IEN(a,e);
        sA(Ac)= sA(Ac) + w*N(a);
        for (int i = 0; i < m; i++) {
          sF(i,Ac) = sF(i,Ac) + w*N(a)*resl(i);
        }
      }
    }

    if (!is_zero(Je)) {
      sE(e) = sE(e) / Je;
    }
  }

  resE = sE;

  // Exchange data at the shared nodes across processes
  all_fun::commu(com_mod, sF);
  all_fun::commu(com_mod, sA);

  for (int a = 0; a < lM.nNo; a++) {
    int Ac = lM.gN(a);
    if (!is_zero(sA(Ac))) {
      for (int i = 0; i < m; i++) {
        res(i,a) = res(i,a) + sF(i,Ac) / sA(Ac);
      }
    }
  }
}

//-------
// tpost
//-------
// Routine for post processing stress tensor
//
void tpost(Simulation* simulation, const mshType& lM, const int m, Array<double>& res, Vector<double>& resE, const Array<double>& lD, 
    const Array<double>& lY, const int iEq, consts::OutputNameType outGrp)
{
  using namespace consts;
  using namespace mat_fun;

  auto& com_mod = simulation->com_mod;
  auto& cep_mod = simulation->cep_mod; 
  auto& cm = com_mod.cm;
  auto& cm_mod = simulation->cm_mod;
  auto& eq = com_mod.eq[iEq];

  #define n_debug_tpost
  #ifdef debug_tpost
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "outGrp: " << outGrp;
  dmsg << "m: " << m;
  #endif

  // [NOTE] Setting gobal variable 'dof'.
  com_mod.dof = eq.dof;

  int i = eq.s;     //  Pointer to start of unknown Yo(:,s:e) so 'i' is an index.
  int j = i + 1;
  int k = j + 1;

  int nFn = lM.nFn;
  if (nFn == 0) {
    nFn = 1;
  }

  #ifdef debug_tpost
  dmsg << "i: " << i;
  dmsg << "j: " << j;
  dmsg << "k: " << k;
  dmsg << "nFn: " << nFn;
  #endif

  // For higher order elements, we lower the order of shape functions
  // to compute the quantities at corner nodes of elements. We then
  // use these lower order shape functions to interpolate the values
  // at edge nodes and elements centers (if applicable)
  //
  bool flag = false; 
  fsType fs;

  if (lM.eType == ElementType::TRI6  || lM.eType == ElementType::QUD8   || 
      lM.eType == ElementType::QUD9  || lM.eType == ElementType::TET10  || 
      lM.eType == ElementType::HEX20 || lM.eType == ElementType::HEX27) {
    flag = true; 
    fs::set_thood_fs(fs, lM.eType);
  } else { 
    fs.eType = lM.eType;
    fs.lShpF = lM.lShpF;
    fs.eNoN  = lM.eNoN;
    fs.nG    = lM.nG;
  }

  #ifdef debug_tpost
  dmsg << "fs.eType: " << fs.eType;
  dmsg << "fs.eNoN: " << fs.eNoN;
  dmsg << "fs.nG: " << fs.nG;
  #endif

  int tnNo = com_mod.tnNo;
  int nsd = com_mod.nsd;
  int tDof = com_mod.tDof;
  int nsymd = com_mod.nsymd;

  #ifdef debug_tpost
  dmsg;
  dmsg << "tnNo: " << tnNo;
  dmsg << "tDof: " << tDof;
  dmsg << "nsymd: " << nsymd;
  #endif

  fs::init_fs(fs, nsd, nsd);

  Vector<double> sA(tnNo); 
  Array<double> sF(m,tnNo); 
  Vector<double> sE(lM.nEl); 
  Array<double> xl(nsd,fs.eNoN); 
  Array<double> dl(tDof,fs.eNoN); 
  Array<double> yl(tDof,fs.eNoN); 
  Array<double> fN(nsd,nFn); 
  Vector<double> resl(m); 
  Array<double> Nx(nsd,fs.eNoN); 
  Vector<double> N(fs.eNoN);

  double ya = 0.0;

  int insd = nsd;
  if (lM.lFib) {
    insd = 1;
  }

  Array<double> Im(nsd, nsd);
  double Je = 0.0; 

  for (int e = 0; e < lM.nEl; e++) {
    int cDmn = all_fun::domain(com_mod, lM, iEq, e);
    auto cPhys = eq.dmn[cDmn].phys;
    if (cPhys != EquationType::phys_struct && cPhys != EquationType::phys_ustruct && cPhys != EquationType::phys_lElas) {
      continue; 
    } 

    double elM = 0.0; 
    double nu = 0.0; 
    double lambda = 0.0; 
    double mu = 0.0; 
    double w = 0.0; 

    if (cPhys == EquationType::phys_lElas) {
      elM = eq.dmn[cDmn].prop[PhysicalProperyType::elasticity_modulus];
      nu = eq.dmn[cDmn].prop[PhysicalProperyType::poisson_ratio];
      lambda = elM*nu / (1.0 + nu) / (1.0 - 2.0*nu);
      mu = 0.5*elM / (1.0 + nu);
    }

    if (lM.eType == ElementType::NRB) {
      // CALL NRBNNX(lM, e)
    }

    fN = 0.0;

    if (lM.fN.size() != 0) {
      for (int l = 0; l < nFn; l++) {
        for (int i = 0; i < nsd; i++) {
          fN(i,l) = lM.fN(i+l*nsd,e);
        }
      }
    }

    dl = 0.0;
    yl = 0.0;

    for (int a = 0; a < fs.eNoN; a++) {
      int Ac = lM.IEN(a,e);
      for (int i = 0; i < nsd; i++) {
        xl(i,a) = com_mod.x(i,Ac);
      }
      for (int i = 0; i < tDof; i++) {
        dl(i,a) = lD(i,Ac);
        yl(i,a) = lY(i,Ac);
      }
    }

    Je = 0.0;
    double Jac = 0.0;

    for (int g = 0; g < fs.nG; g++) {
      if (g == 0  ||  !fs.lShpF) {
        auto Nx_g = fs.Nx.slice(g);
        nn::gnn(fs.eNoN, nsd, insd, Nx_g, xl, Nx, Jac, Im);
      }

      w  = fs.w(g) * Jac;
      N  = fs.N.col(g);
      Je = Je + w;

      auto Im = mat_fun::mat_id(nsd); 
      auto F = Im;

      for (int a = 0; a < fs.eNoN; a++) {
        if (nsd == 3) {
          F(0,0) = F(0,0) + Nx(0,a)*dl(i,a);
          F(0,1) = F(0,1) + Nx(1,a)*dl(i,a);
          F(0,2) = F(0,2) + Nx(2,a)*dl(i,a);
          F(1,0) = F(1,0) + Nx(0,a)*dl(j,a);
          F(1,1) = F(1,1) + Nx(1,a)*dl(j,a);
          F(1,2) = F(1,2) + Nx(2,a)*dl(j,a);
          F(2,0) = F(2,0) + Nx(0,a)*dl(k,a);
          F(2,1) = F(2,1) + Nx(1,a)*dl(k,a);
          F(2,2) = F(2,2) + Nx(2,a)*dl(k,a);
        } else {
          F(0,0) = F(0,0) + Nx(0,a)*dl(i,a);
          F(0,1) = F(0,1) + Nx(1,a)*dl(i,a);
          F(1,0) = F(1,0) + Nx(0,a)*dl(j,a);
          F(1,1) = F(1,1) + Nx(1,a)*dl(j,a);
        }
      }

      double detF = mat_fun::mat_det(F, nsd);

      Vector<double> ed(com_mod.nsymd);

      if (cPhys == EquationType::phys_lElas) {
        for (int a = 0; a < fs.eNoN; a++) {
          if (nsd ==  3) {
            ed(0) = ed(0) + Nx(0,a)*dl(i,a);
            ed(1) = ed(1) + Nx(1,a)*dl(j,a);
            ed(2) = ed(2) + Nx(2,a)*dl(k,a);
            ed(3) = ed(3) + Nx(1,a)*dl(i,a) + Nx(0,a)*dl(j,a);
            ed(4) = ed(4) + Nx(2,a)*dl(j,a) + Nx(1,a)*dl(k,a);
            ed(5) = ed(5) + Nx(0,a)*dl(k,a) + Nx(2,a)*dl(i,a);
          } else { 
            ed(0) = ed(0) + Nx(0,a)*dl(i,a);
            ed(1) = ed(1) + Nx(1,a)*dl(j,a);
            ed(2) = ed(2) + Nx(1,a)*dl(i,a) + Nx(1,a)*dl(j,a);
          }
        }
      }

      switch (outGrp) {

        // Jacobian := determinant of deformation gradient tensor
        case OutputNameType::outGrp_J:
          resl(0) = detF;
          sE(e) = sE(e) + w*detF;
        break;

        //  Deformation gradient tensor (F)
        case OutputNameType::outGrp_F:
          if (nsd == 3) {
            resl(0) = F(0,0);
            resl(1) = F(0,1);
            resl(2) = F(0,2);
            resl(3) = F(1,0);
            resl(4) = F(1,1);
            resl(5) = F(1,2);
            resl(6) = F(2,0);
            resl(7) = F(2,1);
            resl(8) = F(2,2);
          } else { 
            resl(0) = F(0,0);
            resl(1) = F(0,1);
            resl(2) = F(1,0);
            resl(3) = F(1,1);
          }   
        break;

        // Green-Lagrange strain tensor
        case OutputNameType::outGrp_strain:
          if (cPhys == EquationType::phys_lElas) {
            resl = ed;
          } else { 
            auto C = mat_fun::mat_mul(mat_fun::transpose(F), F);
            auto Eg = 0.5 * (C - Im);

            // resl is used to remap Eg
            if (nsd  ==  3) {
              resl(0) = Eg(0,0);
              resl(1) = Eg(1,1);
              resl(2) = Eg(2,2);
              resl(3) = Eg(0,1);
              resl(4) = Eg(1,2);
              resl(5) = Eg(2,0);
            } else {  
              resl(0) = Eg(0,0);
              resl(1) = Eg(1,1);
              resl(2) = Eg(0,1);
            }
          }
        break;

        case OutputNameType::outGrp_stress:
        case OutputNameType::outGrp_cauchy: 
        case OutputNameType::outGrp_mises:
          Array<double> sigma(nsd,nsd);
          Array<double> S(nsd,nsd);

          if (cPhys == EquationType::phys_lElas) {
            if (nsd == 3) {
              double detF = lambda*(ed(0) + ed(1) + ed(2));
              sigma(0,0) = detF + 2.0*mu*ed(0);
              sigma(1,1) = detF + 2.0*mu*ed(1);
              sigma(2,2) = detF + 2.0*mu*ed(2);

              sigma(0,1) = mu*ed(3);
              sigma(1,2) = mu*ed(4);
              sigma(2,0) = mu*ed(5);

              sigma(1,0) = sigma(0,1);
              sigma(2,1) = sigma(1,2);
              sigma(0,2) = sigma(2,0);
            } else {  
              double detF = lambda*(ed(0) + ed(1));
              sigma(0,0) = detF + 2.0*mu*ed(0);
              sigma(1,1) = detF + 2.0*mu*ed(1);
              sigma(0,1) = mu*ed(2);
              sigma(1,0) = sigma(0,1);
            }

          } else if (cPhys == EquationType::phys_ustruct) {
            double p = 0.0;
            for (int a = 0; a < fs.eNoN; a++) {
              p = p + N(a)*yl(k+1,a);
            }
            p = (-p) * detF;

            Array<double> Dm(nsymd,nsymd);
            double Ja;
            
            mat_models::compute_pk2cc(com_mod, cep_mod, eq.dmn[cDmn], F, nFn, fN, ya, S, Dm, Ja);

            // TODO: Add viscous stress

            // Add pressure
            auto C = mat_mul(transpose(F), F);
            S = S + p*mat_inv(C, nsd);

            auto P1 = mat_mul(F, S);
            sigma = mat_mul(P1, transpose(F));

            if (!utils::is_zero(detF)) {
              sigma = sigma / detF;
            }

          } else if (cPhys == EquationType::phys_struct) {
            Array<double> Dm(nsymd,nsymd);
            double Ja;
            mat_models::compute_pk2cc(com_mod, cep_mod, eq.dmn[cDmn], F, nFn, fN, ya, S, Dm, Ja);

            // TODO: Add viscous stress

            auto P1 = mat_mul(F, S);
            sigma = mat_mul(P1, transpose(F));

            if (!utils::is_zero(detF)) {
              sigma = sigma / detF;
            }
          }

          // 2nd Piola-Kirchhoff stress tensor
          if (outGrp == OutputNameType::outGrp_stress) {
            if (nsd == 3) {
              resl(0) = S(0,0);
              resl(1) = S(1,1);
              resl(2) = S(2,2);
              resl(3) = S(0,1);
              resl(4) = S(1,2);
              resl(5) = S(2,0);
            } else { 
              resl(0) = S(0,0);
              resl(1) = S(1,1);
              resl(2) = S(0,1);
            }

          //Cauchy stress tensor
          } else if (outGrp == OutputNameType::outGrp_cauchy) {
            if (nsd == 3) {
              resl(0) = sigma(0,0);
              resl(1) = sigma(1,1);
              resl(2) = sigma(2,2);
              resl(3) = sigma(0,1);
              resl(4) = sigma(1,2);
              resl(5) = sigma(2,0);
            } else { 
              resl(0) = sigma(0,0);
              resl(1) = sigma(1,1);
              resl(2) = sigma(0,1);
            }

          // Von Mises stress
          } else if (outGrp == OutputNameType::outGrp_mises) {
            double trS = mat_trace(sigma, nsd) / static_cast<double>(nsd);
            for (int l = 0; l < nsd; l++) {
              sigma(l,l) = sigma(l,l) - trS;
            }
            double vmises = sqrt(1.5 * mat_ddot(sigma, sigma, nsd));
            resl(0) = vmises;
            sE(e) = sE(e) + w*vmises;
          }
        break;

      } // switch

      for (int a = 0; a < fs.eNoN; a++) { 
        int Ac = lM.IEN(a,e);
        sA(Ac) = sA(Ac) + w*N(a);
        for (int i = 0; i < sF.nrows(); i++) { 
          sF(i,Ac) = sF(i,Ac) + w*N(a)*resl(i);
        }
      }
    }

    if (!utils::is_zero(Je)) {
      sE(e) = sE(e) / Je;
    }

  }

  resE = sE;

  all_fun::commu(com_mod, sF);
  all_fun::commu(com_mod, sA);

  for (int a = 0; a < lM.nNo; a++) {
    int Ac = lM.gN(a);
    if (!utils::is_zero(sA(Ac))) {
      for (int i = 0; i < res.nrows(); i++) { 
        res(i,a) = res(i,a) + sF(i,Ac) / sA(Ac);
      }
    } 
  }

  // For higher order elements, values are interpolated at the edge
  // nodes and element centers using values computed at corners and
  // low-order shape functions
  //
  if (flag) {
    sF = 0.0;
    sA = 0.0;

    Vector<double> xi0(nsd);
    for (int g = 0; g < fs.nG; g++) {
      for (int i = 0; i < fs.xi.nrows(); i++) {
        xi0 = xi0 + fs.xi(i,g);
      }
    }
    xi0 = xi0 / static_cast<double>(fs.nG);

    Array<double> yl(m,fs.eNoN); 
    Vector<double> eNds(tnNo);

    for (int e = 0; e < lM.nEl; e++) {
      int cDmn = all_fun::domain(com_mod, lM, iEq, e);
      auto cPhys = eq.dmn[cDmn].phys;
      if ((cPhys != EquationType::phys_struct) && (cPhys != EquationType::phys_ustruct) && 
          (cPhys != EquationType::phys_lElas)) {
        continue;
      }

      yl = 0.0;
      for (int a = 0; a < fs.eNoN; a++) {
        int Ac = lM.IEN(a,e);
        for (int i = 0; i < nsd; i++) {
          xl(i,a) = com_mod.x(i,Ac);
        }
        for (int i = 0; i < m; i++) {
          yl(i,a) = res(i,lM.lN(Ac));
        }
      }

      double Jac;
      Je = 0.0;
      for (int g = 0; g < fs.nG; g++) {
        if (g == 0 || !fs.lShpF) {
          auto fsNx_g = fs.Nx.slice(g);
          nn::gnn(fs.eNoN, nsd, insd, fsNx_g, xl, Nx, Jac, Im);
        }
      Je = Je + fs.w(g)*Jac;
      }

      for (int a = fs.eNoN; a < lM.eNoN; a++) {
        int Ac = lM.IEN(a,e);
        auto xp = com_mod.x.col(Ac);
        auto xi = xi0;
        nn::get_nnx(nsd, fs.eType, fs.eNoN, xl, fs.xib, fs.Nb, xp, xi, N, Nx);

        resl = 0.0;
        for (int i = 0; i < fs.eNoN; i++) {
          resl = resl + N(i)*yl.col(i);
        }
        i = fs.eNoN - 1;

        if (eNds(Ac) == 0) {
          eNds(Ac) = 1;
        }

        for (int j = 0; j < sF.nrows(); j++) {
          sF(j,Ac) = sF(j,Ac) + resl(j)*Je;
        }
        sA(Ac) = sA(Ac) + Je;
      }

    } // for e = 0; e < lM.nEl

    all_fun::commu(com_mod, sF);
    all_fun::commu(com_mod, sA);

    for (int a = 0; a < lM.nNo; a++) {
      int Ac = lM.gN(a);
      if (eNds(Ac) == 1 && !utils::is_zero(sA(Ac))) {
        for (int i = 0; i < res.nrows(); i++) {
          res(i,a) = res(i,a) + sF(i,Ac) / sA(Ac);
        }
      }
    }
  }
}

};
