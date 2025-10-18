// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "contact.h"

#include "all_fun.h"
#include "consts.h"
#include "lhsa.h"
#include "nn.h"
#include "utils.h"
#include <math.h>

namespace contact {

/// @brief This routine applies penalty-based contact model for possible
/// contacting shell surfaces.
///
/// Reproduces Fortran CONSTRUCT_CONTACTPNLTY.
//
void construct_contact_pnlty(ComMod& com_mod, CmMod& cm_mod, const Array<double>& Dg)
{
  using namespace consts;

  const int nsd  = com_mod.nsd;
  const int tDof = com_mod.tDof;
  const int tnNo = com_mod.tnNo;
  const int dof = com_mod.dof;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  const auto& cntctM = com_mod.cntctM;

  #define n_debug_construct_contact_pnlty
  #ifdef debug_construct_contact_pnlty
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  if (eq.phys != EquationType::phys_shell) {
    return;
  }

  int i = eq.s;
  int j = i + 1;
  int k = j + 1;
  double kl = cntctM.k;
  double hl = cntctM.h;
  #ifdef debug_construct_contact_pnlty
  //dmsg << "kl: " << kl;
  //dmsg << "hl: " << hl;
  //dmsg << "cntctM.c: " << cntctM.c;
  #endif

  // Compute normal vectors at each node in the current configuration
  //
  Array<double> sF(nsd,tnNo); 
  Vector<double> sA(tnNo);

  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    auto& msh = com_mod.msh[iM];
    if (!msh.lShl) {
      continue;
    }
    int eNoN = msh.eNoN;
    int insd = nsd - 1;

    Array<double> Nx(insd,eNoN); 
    Vector<double> N(eNoN); 
    Array<double> xl(nsd,eNoN), gCov(nsd,insd), gCnv(nsd,insd);
    Nx = msh.Nx.slice(0);

    for (int e = 0; e < msh.nEl; e++) {
      for (int a = 0; a < eNoN; a++) {
        int Ac = msh.IEN(a,e);
        xl(0,a) = com_mod.x(0,Ac) + Dg(i,Ac);
        xl(1,a) = com_mod.x(1,Ac) + Dg(j,Ac);
        xl(2,a) = com_mod.x(2,Ac) + Dg(k,Ac);
      }

      Vector<double> nV1(nsd);
      nn::gnns(nsd, eNoN, Nx, xl, nV1, gCov, gCnv);
      double Jac = sqrt(utils::norm(nV1));
      nV1 = nV1 / Jac;

      for (int g = 0; g < msh.nG; g++) {
        double w = msh.w(g) * Jac;
        N = msh.N.col(g);

        for (int a = 0; a < eNoN; a++) {
          int Ac = msh.IEN(a,e);
          sA(Ac) = sA(Ac) + w*N(a);
          for (int i = 0; i < nsd; i++) {
            sF(i,Ac) = sF(i,Ac) + w*N(a)*nV1(i);
          }
        }
      }
    }
  }

  all_fun::commu(com_mod, sF);
  all_fun::commu(com_mod, sA);

  for (int Ac = 0; Ac < tnNo; Ac++) {
    if (!utils::is_zero(sA(Ac))) {
      for (int i = 0; i < nsd; i++) {
        sF(i,Ac) = sF(i,Ac) / sA(Ac);
      }
    }

    double Jac = sqrt(sF.col(Ac) * sF.col(Ac));

    if (!utils::is_zero(Jac)) {
      for (int i = 0; i < nsd; i++) {
        sF(i,Ac) = sF(i,Ac) / Jac;
      }
    }
  }

  // Create a bounding box around possible region of contact and bin
  // the box with neighboring nodes
  //
  int maxNnb = 15;
  label_101: maxNnb = maxNnb + 5;

  Array<int> bBox(maxNnb,tnNo);
  bBox = -1;

  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    auto& msh = com_mod.msh[iM];

    if (!msh.lShl) {
      continue;
    }

    Vector<double> x1(nsd), xmin(nsd), xmax(nsd);

    for (int a = 0; a < msh.nNo; a++) {
      int Ac = msh.gN(a);
      x1(0) = com_mod.x(0,Ac) + Dg(i,Ac);
      x1(1) = com_mod.x(1,Ac) + Dg(j,Ac);
      x1(2) = com_mod.x(2,Ac) + Dg(k,Ac);

      // Box limits for each node
      xmin = x1 - cntctM.c;
      xmax = x1 + cntctM.c;

      // Load the box with neighboring nodes lying within it
      //
      for (int jM = 0; jM < com_mod.nMsh; jM++) {
        auto& msh_jM = com_mod.msh[jM];
        if (iM == jM || !msh_jM.lShl) {
          continue;
        }

        Vector<double> x2(nsd); 

        for (int b = 0; b < msh_jM.nNo; b++) {
          int Bc = msh_jM.gN(b);
          x2(0) = com_mod.x(0,Bc) + Dg(i,Bc);
          x2(1) = com_mod.x(1,Bc) + Dg(j,Bc);
          x2(2) = com_mod.x(2,Bc) + Dg(k,Bc);

          if ((x2(0) >= xmin(0)) && (x2(0) <= xmax(0)) && (x2(1) >= xmin(1)) && (x2(1) <= xmax(1)) && 
              (x2(2) >= xmin(2)) && (x2(2) <= xmax(2))) {

            int l = 0;

            for (int i = 0; i < maxNnb; i++, l++) {
              if (bBox(l,Ac) == -1) {
                bBox(l,Ac) = Bc;
                break; 
              }
              if (Bc > bBox(l,Ac)) {
                continue;
              }

              if (Bc == bBox(l,Ac)) {
                break;
              } 

              if (bBox(maxNnb-1,Ac) != -1) goto label_101;

              for (int m = maxNnb-1; m >= l; m--) {
                bBox(m,Ac) = bBox(m-1,Ac);
              }
              bBox(l,Ac) = Bc;
              break; 
            }
            if (l > maxNnb-1) goto label_101;
          }
        } // b
      } // jM
    } // a
  } // iM

  // Check if any node is strictly involved in contact and compute
  // corresponding penalty forces assembled to the residual
  //
  Array<double> lR(dof,tnNo); 
  Vector<int> incNd(tnNo);
  Vector<double> x1(nsd), x2(nsd);

  for (int Ac = 0; Ac < tnNo; Ac++) {
    if (bBox(0,Ac) == -1) {
      continue; 
    }
    x1(0) = com_mod.x(0,Ac) + Dg(i,Ac);
    x1(1) = com_mod.x(1,Ac) + Dg(j,Ac);
    x1(2) = com_mod.x(2,Ac) + Dg(k,Ac);
    auto nV1 = sF.rcol(Ac);
    int nNb = 0;

    for (int a = 0; a < maxNnb; a++) {
      int Bc = bBox(a,Ac);
      if (Bc == -1) {
        continue; 
      }
      x2(0) = com_mod.x(0,Bc) + Dg(i,Bc);
      x2(1) = com_mod.x(1,Bc) + Dg(j,Bc);
      x2(2) = com_mod.x(2,Bc) + Dg(k,Bc);
      auto nV2 = sF.rcol(Bc);

      auto x12 = x1 - x2;
      double c = sqrt(utils::norm(x12));
      double al = sqrt(fabs(utils::norm(nV1, nV2)));

      if (c <= cntctM.c && al >= cntctM.al) {
        double d = utils::norm(x12, nV2);
        bool flag = false;
        double pk{0.0};
   
        if (d >= -cntctM.h && d < 0.0) {
          pk = 0.5 * (kl / hl) * pow(d+hl, 2.0);
          flag = true;
        } else if (d >= 0.0) {
          pk = 0.5 * kl * (hl + d);
          flag = true;
        } else {
          pk = 0.0;
        }

        if (flag) {
          incNd(Ac) = 1;
          nNb = nNb + 1;
          for (int i = 0; i < nsd; i++) {
            lR(i,Ac) = lR(i,Ac) - pk*nV1(i);
          }
          #ifdef debug_construct_contact_pnlty
          dmsg << "          " << " ";
          dmsg << "Ac: " << Ac+1;
          dmsg << "Bc: " << Bc+1;
          dmsg << "nV1: " << nV1;
          dmsg << "nV2: " << nV2;
          dmsg << "pk: " << pk;
          dmsg << "x12: " << x12;
          dmsg << "d: " << d;
          #endif
        }
      }
    }

    if (nNb != 0) {
      for (int i = 0; i < dof; i++) {
        lR(i,Ac) = lR(i,Ac) / static_cast<double>(nNb);
      }
    }

  }

  // Return if no penalty forces are to be added
  if (incNd.sum() == 0) {
    return;
  }

  // Global assembly
  //
  for (int Ac = 0; Ac < tnNo; Ac++) {
    if (incNd(Ac) == 0) {
      continue;
    }
    for (int i = 0; i < dof; i++) {
      com_mod.R(i,Ac) = com_mod.R(i,Ac) + lR(i,Ac);
    }
  }

}

};
