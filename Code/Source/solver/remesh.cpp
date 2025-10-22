// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "remesh.h"

#include "all_fun.h"
#include "mat_fun.h"
#include "nn.h"
#include "output.h"
#include "read_msh.h"
#include "remeshTet.h"
#include "vtk_xml.h"

#include <array>
#include<iostream>
#include <filesystem>
#include<fstream>

namespace remesh {

/// @brief Distribute the new mesh elements to all processors
//
void distre(ComMod& com_mod, CmMod& cm_mod, mshType& lM, int& nEl, Vector<int>& gE)
{
  auto& cm = com_mod.cm;
  int gnEl = lM.gnEl;
  int eNoN = lM.eNoN;

  #define n_debug_distre
  #ifdef debug_distre
  DebugMsg dmsg(__func__, cm.idcm());
  dmsg.banner();
  dmsg << "eNoN: " << eNoN;
  #endif

  Vector<int> part(gnEl);
  nEl = 0;

  for (int e = 0; e < gnEl; e++) {
    int i = 0;
    for (int a = 0; a < eNoN; a++) {
      int Ac = lM.gIEN(a,e);
      if (lM.gpN(Ac) == cm.tF(cm_mod)) {
        i = i + 1;
      }
    }

    if (i == eNoN) { 
      nEl = nEl + 1;
      part(e) = 1;
    }
  }

  gE.resize(nEl);
  nEl = 0;

  for (int e = 0; e < gnEl; e++) {
    if (part(e) > 0) {
      gE(nEl) = e;
      nEl = nEl + 1;
    }
  }
}

/// @brief Reproduces Fortran 'SUBROUTINE DISTMSHSRF(lFa, lM, iOpt)'
//
void dist_msh_srf(ComMod& com_mod, ChnlMod& chnl_mod, faceType& lFa, mshType& lM, const int iOpt)
{
  const int nsd = com_mod.nsd;
  auto& rmsh = com_mod.rmsh;

  #define n_debug_dist_msh_srf 
  #ifdef debug_dist_msh_srf 
  auto& cm = com_mod.cm;
  DebugMsg dmsg(__func__, cm.idcm());
  dmsg.banner();
  #endif

  std::string sTmp = chnl_mod.appPath + "/" + ".remesh_tmp.dir";
  auto mkdir_arg = std::string("mkdir -p ") + sTmp;
  std::system(mkdir_arg.c_str());

  for (int e = 0; e < lFa.nEl; e++) {
    for (int a = 0; a < lFa.eNoN; a++) {
      int Ac = lFa.IEN(a,e);
      Ac = lFa.gN(Ac);
      lFa.IEN(a,e) = Ac;
    }
  }

  for (int iFa = 0; iFa < lM.nFa; iFa++) {
    lM.fa[iFa].nEl = lM.fa[iFa].gnEl;
    #ifdef debug_dist_msh_srf 
    dmsg << "---------- iFa: " << iFa;
    dmsg << "lM.fa[iFa].name: " << lM.fa[iFa].name;
    dmsg << "lM.fa[iFa].nEl: " << lM.fa[iFa].nEl;
    #endif

    lM.fa[iFa].IEN.resize(lM.fa[iFa].eNoN,lM.fa[iFa].nEl);

    lM.fa[iFa].gE.resize(lM.fa[iFa].nEl);

    int eoff = 0;
    if (iFa > 0) {
      for (int i = 0; i < iFa; i++) {
        eoff += lM.fa[i].gnEl;
      }
    }
    #ifdef debug_dist_msh_srf 
    dmsg << "eoff: " << eoff;
    #endif

    for (int e = 0; e < lM.fa[iFa].gnEl; e++) {
      lM.fa[iFa].gE(e) = lFa.gE(eoff+e);
      for (int i = 0; i < lM.fa[iFa].eNoN; i++) {
        lM.fa[iFa].IEN(i,e) = lFa.IEN(i,eoff+e);
      }
    }

    read_msh_ns::calc_nbc(lM, lM.fa[iFa]);

    lM.fa[iFa].x.resize(nsd, lM.fa[iFa].nNo);

    for (int a = 0; a < lM.fa[iFa].nNo; a++) {
      int Ac = lM.fa[iFa].gN(a);
      for (int i = 0; i < nsd; i++) {
        lM.fa[iFa].x(i,a) = lM.x(i,Ac);
      }
    }

    for (int i = 0; i < lM.fa[iFa].gebc.ncols(); i++) {
      lM.fa[iFa].gebc(0,i) = lM.fa[iFa].gE(i);
    }

    for (int i = 0; i < lM.fa[iFa].eNoN; i++) {
      for (int j = 0; j < lM.fa[iFa].nEl; j++) {
        lM.fa[iFa].gebc(i+1,j) = lM.fa[iFa].IEN(i,j);
      }
    }

    if (iOpt == 1) {
      std::string fTmp = sTmp + "/" + lM.fa[iFa].name + "_" + std::to_string(rmsh.rTS) + ".vtp";
      Vector<int> incNd(lM.gnNo);

      for (int a = 0; a < lM.fa[iFa].nNo; a++) {
        int Ac = lM.fa[iFa].gN(a);
        incNd(Ac) = a;
      }

      for (int e = 0; e < lM.fa[iFa].nEl; e++) {
        for (int a = 0; a < lM.fa[iFa].eNoN; a++) {
          int Ac = lM.fa[iFa].IEN(a,e);
          lM.fa[iFa].IEN(a,e) = incNd(Ac);
        }
      }

      vtk_xml::write_vtp(com_mod, lM.fa[iFa], fTmp);

      for (int e = 0; e < lM.fa[iFa].nEl; e++) {
        for (int a = 0; a < lM.fa[iFa].eNoN; a++) {
          int Ac = lM.fa[iFa].IEN(a,e);
          Ac = lM.fa[iFa].gN(Ac);
          lM.fa[iFa].IEN(a,e) = Ac;
        }
      }
    }
  }
}

/// @brief Modifies
//   gN(nNo) - list of node indices 0  
//   lM.gpN - processor ID (1, 2, 3, ...) for each node
//
void distrn(ComMod& com_mod, CmMod& cm_mod, const int iM, mshType& lM, Array<double>& Dg, int& nNo, Vector<int>& gN)
{
  const int nsd = com_mod.nsd;
  auto& rmsh = com_mod.rmsh;
  auto& cm = com_mod.cm;

  #define n_debug_distrn
  #ifdef debug_distrn
  DebugMsg dmsg(__func__, cm.idcm());
  dmsg.banner();
  #endif

  int gnNo = lM.gnNo;
  int i = 0;
  double f = 2.5e-2;
  nNo = 0;

  Vector<int> part(gnNo), tmpI(gnNo);

  #ifdef debug_distrn
  dmsg << "iM: " << iM;
  dmsg << "lM.nNo: " << lM.nNo;
  dmsg << "lM.gnNo: " << lM.gnNo;
  dmsg << "lM.nEl: " << lM.nEl;
  dmsg << "lM.IEN.nrows(): " << lM.IEN.nrows();
  dmsg << "lM.IEN.ncols(): " << lM.IEN.ncols();
  dmsg << "lM.gnEl: " << lM.gnEl;
  dmsg << "lM.gIEN.nrows(): " << lM.gIEN.nrows();
  dmsg << "lM.gIEN.ncols(): " << lM.gIEN.ncols();
  dmsg << "gnNo: " << gnNo;
  dmsg << "msh(iM).nNo: " << com_mod.msh[iM].nNo;
  dmsg << "cm.tF(cm_mod): " << cm.tF(cm_mod);
  #endif

  // Set nodes for each partition.
  //  
  // lM.x are remeshed nodes (size 3 x gnNo).
  //
  // x are original nodes (size 3 x msh(iM).nNo).
  // 
  while (true) {
    part = 0;
    tmpI = 0;
    nNo = 0;
    f = 2.0 * f;
    double tol = (1.0 + f) * rmsh.maxEdgeSize(iM);
    i = i+1;
    //dmsg << "---------- pass " << std::to_string(i) + " ----------";
    //dmsg << "tol: " << tol;

    for (int a = 0; a < gnNo; a++) {
      //dmsg << "---------- a " << std::to_string(a) + " ----------";
      if (part(a) != 0) {
        continue; 
      }
      double minS = std::numeric_limits<double>::max();

      for (int b = 0; b < com_mod.msh[iM].nNo; b++) {
        int Ac = com_mod.msh[iM].gN(b);
        double dS = 0.0;

        for (int i = 0; i < nsd; i++) {
          double diff = com_mod.x(i,Ac) + Dg(i,Ac) - lM.x(i,a);
          dS += diff * diff;
        }
        dS = sqrt(dS);
        //dS = SQRT(SUM((x(:,Ac)+Dg(:,Ac)-lM.x(:,a))**2._RKIND))

        if (minS > dS) minS = dS;
      }

      if (minS < tol) {
        nNo = nNo + 1;
        part(a) = cm.tF(cm_mod);
      }
    }

    MPI_Allreduce(part.data(), tmpI.data(), gnNo, cm_mod::mpint, MPI_MAX, cm.com());

    int b = 0;

    for (int a = 0; a < gnNo; a++) {
      if (tmpI(a) > 0) b = b + 1;
    }
    //dmsg << "b: " << b;

    if (b == gnNo) {
      break;
    } else { 
      std::cout << "[distrn] Found only " + std::to_string(b) + " nodes in pass " + std::to_string(i) + 
        " out of " + std::to_string(gnNo) + " nodes." << std::endl;

      if (i > 6) {
        throw std::runtime_error("[distrn] Could not distribute all nodes in " + std::to_string(i) + " passes.");
      }
    }
  }

  gN.resize(nNo);
  lM.gpN.resize(gnNo);
  nNo = 0;

  for (int a = 0; a < gnNo; a++) {
    if (part(a) == cm.tF(cm_mod)) {
      gN(nNo) = a;
      nNo = nNo + 1;
    }
  }

  lM.gpN = part;
}

//--------
// find_n
//--------
// Find element in the old mesh for each node in the new mesh
//
// Reproduces Fortran 'SUBROUTINE FINDN(xp, iM, Dg, eList, Ec, Nsf)'
//
void find_n(ComMod& com_mod, const Vector<double>& Xp, const int iM, const Array<double>& Dg, 
    const Vector<int>& eList, int& Ec, Vector<double>& Nsf, int pcount)
{
  const int nsd = com_mod.nsd;
  auto& msh = com_mod.msh[iM];
  int ne = eList.size();
  Array<double> Amat(nsd+1,msh.eNoN);

  Nsf = 0.0;

  for (int e = 0; e < ne; e++) {
    //std::cout << "[find_n] --------------- e " << e+1 << std::endl;
    Ec = eList(e);

    if (Ec == -1) {
      break;
    }

    Amat = 1.0;

    for (int a = 0; a <  msh.eNoN; a++) {
      int Ac = msh.IEN(a,Ec);
      for (int i = 0; i < nsd; i++) {
        Amat(i,a) = com_mod.x(i,Ac) + Dg(i,Ac);
      }
    }

    Amat = mat_fun::mat_inv(Amat, msh.eNoN);

    double tol = 1e-14;
    int a = 0;

    for (int i = 0; i < nsd+1; i++) {
      Nsf(i) = 0.0;

      for (int j = 0; j < msh.eNoN; j++) {
        Nsf(i) = Nsf(i) + Amat(i,j)*Xp(j);
      }

      if ( (Nsf(i) > -tol) && (Nsf(i) < (1.0+tol))) {
        a = a + 1;
      }
    }

    if (a == nsd+1) {
      return;
    }
  }

  Ec = -1;
  Nsf = 0.0;
}

/// @brief Create list of connected/adjacent elements for old/source mesh
///
/// Reproduces Fortran 'SUBROUTINE GETADJESRC(lM, kneList)'
//
void get_adj_esrc(ComMod& com_mod, mshType& lM, Array<int>& kneList)
{
  #define n_debug_get_adj_esrc
  #ifdef debug_get_adj_esrc
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  auto& rmsh = com_mod.rmsh;

  // First get elements around all nodes
  //
  Vector<int> nL(lM.nNo);

  for (int e = 0; e < lM.nEl; e++) {
    for (int a = 0; a < lM.eNoN; a++) {
      int Ac = lM.IEN(a,e);
      Ac = lM.lN(Ac);
      nL(Ac) = nL(Ac) + 1;
    }
  }

  int maxKNE = 0;
  if (nL.size() != 0) {
    maxKNE = std::max(nL.max(), 0);
  }

  Array<int> tmpList(maxKNE, lM.nNo);
  tmpList = -1;  // Stores element IDs
  nL = 0;        // Stores element counts

  for (int e = 0; e < lM.nEl; e++) {
    for (int a = 0; a < lM.eNoN; a++) {
      int Ac = lM.IEN(a,e);
      Ac = lM.lN(Ac);
      tmpList(nL(Ac), Ac) = e;
      nL(Ac) = nL(Ac) + 1;
    }
  }

  int b = 2 * maxKNE;

  // Now get elements around each element and avoid redundancy
  //
  get_elems:
    b = b + maxKNE;
    nL.resize(lM.nEl);
    kneList.resize(b, lM.nEl);
    kneList = -1;
    bool flag = false;

    for (int e = 0; e < lM.nEl; e++) {
      for (int a = 0; a < lM.eNoN; a++) {
        int Ac = lM.IEN(a,e);
        Ac = lM.lN(Ac);

        for (int i = 0; i < maxKNE; i++) {
          if (tmpList(i,Ac) != -1) {
            flag = true;

            for (int j = 0; j < nL(e); j++) {
              if (kneList(j,e) == tmpList(i,Ac)) {
                flag = false;
                break;
              }
            }

            if (flag) {
              nL(e) = nL(e) + 1;
              if (nL(e) >= b) goto get_elems;
              kneList(nL(e)-1, e) = tmpList(i,Ac);
            }
          } else {
            break;
          }

        } // for i

      } // for a

    } // for e

  maxKNE = 0;
  if (nL.size() != 0) {
    maxKNE = std::max(nL.max(), 0);
  }

  tmpList.resize(b,lM.nEl);
  tmpList = kneList;
  kneList.resize(maxKNE, lM.nEl);
  kneList = -1;

  for (int e = 0; e < lM.nEl; e++) {
    for (int i = 0; i < maxKNE; i++) {
      kneList(i, e) = tmpList(i, e);
    }
  }
}

/// @brief Create list of connected/adjacent nodes for new/target mesh
///
/// Reproduces Fortran 'SUBROUTINE GETADJNTGT(lM, nNo, nEl, gN, gE, knnList)'
//
void get_adj_ntgt(ComMod& com_mod, mshType& lM, const int nNo, const int nEl, const Vector<int>& gN, 
    const Vector<int>& gE, Array<int>& knnList)
{
  Vector<int> lN(lM.gnNo);

  for (int a = 0; a < nNo; a++) {
    int Ac = gN(a);
    lN(Ac) = a;
  }

  Vector<int> nL(nNo);

  for (int e = 0; e < nEl; e++) {
    int Ec = gE(e);

    for (int a = 0; a < lM.eNoN; a++) {
      int Ac = lM.gIEN(a,Ec);
      int i = lN(Ac);

      for (int b = 0; b < lM.eNoN; b++) {
        if (b == a) {
          continue;
        }
        int Bc = lM.gIEN(b,Ec);
        nL(i) = nL(i) + 1;
      }
    }
  }

  int maxKNN = 0;
  if (nL.size() != 0) {
    maxKNN = std::max(nL.max(), 0);
  }

  Array<int> tmpList(maxKNN, nNo);

  tmpList = -1;
  nL = 0;

  for (int e = 0; e < nEl; e++) {
    int Ec = gE(e);

    for (int a = 0; a < lM.eNoN; a++) {
      int Ac = lM.gIEN(a,Ec);
      int i = lN(Ac);

      for (int b = 0; b < lM.eNoN; b++) {
        if (a == b) {
          continue;
        }

        int Bc = lM.gIEN(b,Ec);
        int j = lN(Bc);
        bool flag = true;

        for (int n = 0; n < nL(i); n++) {
          if (tmpList(n,i) ==  j) {
            flag = false;
            break;
          }
        }

        if (flag) {
          tmpList(nL(i), i) = j;
          nL(i) = nL(i) + 1;
        }
      }
    }
  }

  maxKNN = 0;
  if (nL.size() != 0) {
    maxKNN = std::max(nL.max(), 0);
  }
  knnList.resize(maxKNN, nNo);
  knnList = -1;

  for (int a = 0; a < nNo; a++) {
    //std::cout << "[get_adj_ntgt] " << a+1 << std::endl;
    for (int i = 0; i < maxKNN; i++) {
      knnList(i,a) = tmpList(i,a);
      //std::cout << knnList(i,a)+1 << " ";
    }
    //std::cout << std::endl;
  }

}

/// @brief Interpolation of data variables from source mesh to target mesh
//
void interp(ComMod& com_mod, CmMod& cm_mod, const int lDof, const int iM, mshType& tMsh, Array<double>& sD, Array<double>& tgD)
{
  const int nsd = com_mod.nsd;
  const int tnNo = com_mod.tnNo;
  auto& msh = com_mod.msh;
  auto& rmsh = com_mod.rmsh;
  auto& cm = com_mod.cm;
  #define n_debug_interp
  #define n_debug_interp_1
  #ifdef debug_interp
  DebugMsg dmsg(__func__, cm.idcm());
  dmsg.banner();
  #endif

  if (nsd+1 != msh[iM].eNoN) {
    throw std::runtime_error("[interp] Inconsistent element type for interpolation. Can support 2D Tri or 3D Tet elements only.");
  }

  int gnNo = tMsh.gnNo;
  int eNoN = tMsh.eNoN;
  int gnEl = tMsh.gnEl;
  int nNo = 0;
  #ifdef debug_interp
  dmsg << "tMsh.gnNo: " << tMsh.gnNo;
  dmsg << "tMsh.eNoN: " << tMsh.eNoN;
  dmsg << "tMsh.gnEl: " << tMsh.gnEl;
  #endif

  // Distribute the new mesh nodes among all the processors
  // Need to transfer mesh displacement
  //
  Vector<int> gN;
  Array<double> Dg(nsd,tnNo);
  int i = nsd+1;

  for (int j = 0; j < nsd; j++) { 
    for (int k = 0; k < tnNo; k++) { 
      Dg(j,k) = rmsh.D0(j+nsd+1,k);
    }
  }

  #ifdef debug_interp
  dmsg << "Distribute nodes ..." << "";
  #endif

  distrn(com_mod, cm_mod, iM, tMsh, Dg, nNo, gN);
  //CALL DISTRN(iM, tMsh, Dg, nNo, gN)
  #ifdef debug_interp
  dmsg << "nNo: " << nNo;
  #endif

  // Distribute elements of the new mesh to all processors
  #ifdef debug_interp
  dmsg << "Distribute elements ..." << "";
  #endif
  int nEl = 0;
  Vector<int> gE;
  distre(com_mod, cm_mod, tMsh, nEl, gE);
  //CALL DISTRE(tMsh, nEl, gE)
  #ifdef debug_interp
  dmsg << "nEl: " << nEl;
  #endif

  // Setup data structures for octree search
  // Get adjacent cells for source (old) mesh
  //
  #ifdef debug_interp
  dmsg << "Setup data structures for octree search ... " << "";
  #endif
  Array<int> srcAdjEl;
  get_adj_esrc(com_mod, msh[iM], srcAdjEl);
  //CALL GETADJESRC(msh(iM), srcAdjEl)
  int maxKNE = srcAdjEl.nrows();
  #ifdef debug_interp
  dmsg << "maxKNE: " << maxKNE;
  //dmsg << "srcAdjEl: " << srcAdjEl;
  #endif

  // Get adjacent nodes for each node on the new mesh
  //
  // tgtAdjNd stores node IDs.
  //
  #ifdef debug_interp
  dmsg << "Get adjacent nodes " << " ...";
  #endif
  Array<int> tgtAdjNd;
  get_adj_ntgt(com_mod, tMsh, nNo, nEl, gN, gE, tgtAdjNd);
  //CALL GETADJNTGT(tMsh, nNo, nEl, gN, gE, tgtAdjNd)

  int maxKNN = tgtAdjNd.nrows();
  //DEALLOCATE(gE)
  #ifdef debug_interp
  dmsg << "maxKNN: " << maxKNN;
  dmsg << "nNo: " << nNo;
  //dmsg << "tgtAdjNd: " << tgtAdjNd;
  #endif

  Vector<double> Xp(nsd+1), Nsf(eNoN); 
  Array<double> gNsf(eNoN,nNo); 
  Vector<int> tagNd(gnNo), rootEl(nNo), masEList(msh[iM].nEl);
  //ALLOCATE(Xp(nsd+1), Nsf(eNoN), gNsf(eNoN,nNo), tagNd(gnNo), gE(nNo), rootEl(nNo), chckNp(nNo), masEList(msh(iM).nEl))
  rootEl = -1;

  std::vector<bool> chckNp(nNo);
  std::fill(chckNp.begin(), chckNp.end(), false);

  for (int e = 0; e < msh[iM].nEl; e++) {
    masEList(e) = e;
  }

  // Determine boundary nodes on the new mesh, where interpolation is
  // not needed, or boundary search is performed
  //
  Vector<int> tmpL(gnNo);
  //tmpL = -1;
  #ifdef debug_interp
  dmsg << "gnNo: " << gnNo;
  #endif

  for (int e = 0; e < tMsh.fa[0].nEl; e++) {
    for (int a = 0; a < tMsh.fa[0].eNoN; a++) {
      int Ac = tMsh.fa[0].IEN(a,e);
      if (tmpL(Ac) == 0) {
        tmpL(Ac) = 1;
      }
    }
  }

  // srfNds is a really a bool array.
  //
  Vector<int>srfNds(nNo);
  srfNds = 0;
  #ifdef debug_interp
  dmsg << "nNo: " << nNo;
  #endif

  for (int a = 0; a < nNo; a++) {
    int Ac = gN(a);
    if (tmpL(Ac) > 0) {
      srfNds(a) = 1;
    }
  }

  //DEALLOCATE(tmpL)

  // tagNd stores procesors IDs ?
  int bTag = 2*cm.np();

  for (int a = 0; a < nNo; a++) {
    int Ac = gN(a);

    if (srfNds(a) > 0) {
      chckNp[a] = true;
      tagNd(Ac) = bTag;
      gE(a) = -1;

      for (int i = 0; i < eNoN; i++) {
        gNsf(i,a) = 0.0;
      }
    }
  }
  
  // Find a starting element through brute search and the starting node
  // is used to initialize the queue (local numbering)
  //
  #ifdef debug_interp
  dmsg << "Find a starting element ... " << "";
  #endif
  utils::queueType rootNdQ;

  for (int a = 0; a < nNo; a++) {
    #ifdef debug_interp_1
    dmsg << "---------- a " << std::to_string(a+1) + " ----------";
    #endif
    if (chckNp[a]) {
      continue; 
    }
    int Ac = gN(a);
    Xp = 1.0;
    int Ec = -1;
    for (int i = 0; i < nsd; i++) {
      Xp(i) = tMsh.x(i,Ac);
    }

    find_n(com_mod, Xp, iM, Dg, masEList, Ec, Nsf, 0);
    #ifdef debug_interp_1
    dmsg << "Ec: " << Ec+1; 
    #endif
    //CALL FINDN(Xp, iM, Dg, masEList, Ec, Nsf)
    chckNp[a] = true;

    // Once an element in the source mesh is found, set 'this' element as
    // 'root' element for all the neighbor nodes of the target mesh. Then
    // add these neighbor nodes to the queue to form the 'front'
    //
    if (Ec > -1) {
      gE(a) = Ec;
      tagNd(Ac) = cm.tF(cm_mod);
      #ifdef debug_interp_s
      if (Ac == 8864-1) {
        dmsg << "#### 1 Ac: " << Ac+1; 
      }
      #endif

      for (int i = 0; i < eNoN; i++) {
        gNsf(i,a) = Nsf(i);
      }
      //gNsf(:,a) = Nsf(:)

      for (int nn = 0; nn < maxKNN; nn++) {
        int b = tgtAdjNd(nn,a);

        if (b > -1) {
          #ifdef debug_interp_1
          dmsg << "nn: " << std::to_string(nn+1) +  "  b: " + std::to_string(b+1); 
          #endif
          utils::enqueue(rootNdQ, b);
          rootEl(b) = Ec;
        }
      }

      if (rootNdQ.n > 1) {
        break;
      }
    }
  }

  // Node-Cell search begins here. Uses fastest grid-to-grid algorithm
  // based on advancing front in the nearest vicinity
  //
  #ifdef debug_interp
  dmsg << " " << "";
  dmsg << "Node-Cell search begins ... " << "";
  #endif
  Vector<int> eList(maxKNE); 
  int probe;
  int pcount = 0;
 
  while (dequeue(rootNdQ, probe)) {
    #ifdef debug_interp_1
    dmsg << " " << " ";
    dmsg << "---------- probe " << std::to_string(probe+1) + " ----------";
    dmsg << "pcount: " << pcount+1;
    #endif
    pcount += 1;

    if (std::all_of(chckNp.begin(), chckNp.end(), [](bool v) { return v; })) {
      break;
    }

    if (chckNp[probe]) {
      continue;
    }

    int Ac = gN(probe);
    Xp = 1.0;
    for (int i = 0; i < nsd; i++) {
      Xp(i) = tMsh.x(i,Ac);
    }
    eList = -1;
    bool flag = false;
    #ifdef debug_interp_1
    dmsg << "Ac: " << Ac+1;
    dmsg << "rootEl(probe): " << rootEl(probe)+1;
    #endif

    for (int e = 0; e < maxKNE; e++) {
      int Ec = srcAdjEl(e,rootEl(probe));
      if (Ec > -1) {
        eList(e) = Ec;
        #ifdef debug_interp_1
        //dmsg << "srcAdjEl:  e Ec: " << std::to_string(e+1) + " " + std::to_string(Ec+1);
        #endif
      }
    }

    for (int itry = 1; itry <= 2; itry++) {
      #ifdef debug_interp_1
      dmsg << ">>>>> try: " << itry; 
      dmsg << "pcount: " << pcount; 
      dmsg << "Xp: " << Xp; 
      dmsg << "eList: " << eList; 
      #endif
      int Ec;
      find_n(com_mod, Xp, iM, Dg, eList, Ec, Nsf, pcount);
      #ifdef debug_interp_1
      dmsg << "Ec: " << Ec+1; 
      #endif

      #ifdef debug_interp_s
      if (pcount == 19195+1) {
        dmsg << "eList: " << eList; 
        exit(0);
      }
      #endif

      // If an element is found, the neighbors of the node are initialized
      // with this element as root element and are added to the queue
      // (front propagation)

      if (Ec > -1) {
        gE(probe) = Ec;
        tagNd(Ac) = cm.tF(cm_mod);

        for (int i = 0; i < Nsf.size(); i++) {
          gNsf(i,probe) = Nsf(i);
        }

        for (int nn = 0; nn < maxKNN; nn++) {
          int a = tgtAdjNd(nn,probe);

          if (a > -1) {
            #ifdef debug_interp_1
            //dmsg << "tgtAdjNd:   nn a Ec: " << std::to_string(nn+1) + " " + std::to_string(a+1) + " " + std::to_string(Ec+1); 
            #endif
            utils::enqueue(rootNdQ, a);
            rootEl(a) = Ec;
          }
        }

        eList.resize(maxKNE); 
        eList = -1;
        break;
      }

      // If failed in first attempt, try again with increased search area
      // but add neighbors uniquely without repetition
      //
      if (itry == 1) {
        tmpL.resize(nEl);
        tmpL = 0;
        int ne = 0;
        flag = false; 

        for (int e = 0; e < maxKNE; e++) {
          int Ec = eList(e);

          if (Ec > -1) {
            for (int a = 0; a < maxKNE; a++) {
              int Bc = srcAdjEl(a,Ec);

              if (Bc != -1 && Bc != Ec) {
                flag = true; 

                for (int i = 0; i < ne; i++) {
                  if (Bc == tmpL(i)) {
                    flag = false; 
                    break;
                  }
                }

                if (flag) {
                  tmpL(ne) = Bc;
                  ne = ne + 1;
                }
              }
            }
          }
        }

        eList.resize(ne);
        eList = -1;
        for (int i = 0; i < ne; i++) {
          eList(i) = tmpL(i);
        }

      } else if (itry == 2) {
        eList.resize(maxKNE); 
        eList = -1;
      }
    } // itry

    chckNp[probe] = true;
  }

  tmpL.resize(gnNo);
  tmpL = 0;

  #ifdef debug_interp
  dmsg << "MPI_Allreduce ... " << "";
  #endif

  MPI_Allreduce(tagNd.data(), tmpL.data(), gnNo, cm_mod::mpint, MPI_MAX, cm.com());

  // tmpL was previously a list of element IDs but it now holds
  // processor IDs.
  //
  int nn = 0;
  int nbnd = 0;

  for (int a = 0; a < gnNo; a++) {
    if (tmpL(a) > 0) {
      nn = nn + 1;
    }
    if (tmpL(a) == 2*cm.np()) {
      nbnd = nbnd + 1;
    }
  }
  #ifdef debug_interp
  dmsg << "nn: " << nn;
  #endif

  // Assign tag for nodes that were interpolated in other processors
  //
  #ifdef debug_interp
  dmsg << "Assign tag for nodes ... " << "";
  #endif
  tagNd = 0;

  for (int a = 0; a < nNo; a++) {
    int Ac = gN(a);
    tagNd(Ac) = tmpL(Ac);
  }

  // Now use brute force to find left over non-interpolated nodes
  //
  #ifdef debug_interp
  dmsg << "Find left over non-interpolated nodes ... " << "";
  #endif
  for (int a = 0; a < nNo; a++) {
    int Ac = gN(a);

    if (tagNd(Ac) == 0) {
      Xp = 1.0;
      for (int i = 0; i < nsd; i++) {
        Xp(i) = tMsh.x(i,Ac);
      }

      int Ec = -1;
      find_n(com_mod, Xp, iM, Dg, masEList, Ec, Nsf, 0);

      if (Ec > -1) {
        gE(a) = Ec;
        tagNd(Ac) = cm.tF(cm_mod);
        for (int i = 0; i < Nsf.size(); i++) {
          gNsf(i,a) = Nsf(i);
        }
      }
    }
  }

  // Nodes belonging to other procs are reassigned 0
  #ifdef debug_interp
  dmsg << "Nodes in other procs set to 0 ..." << "";
  #endif
  for (int a = 0; a < nNo; a++) {
    int Ac = gN(a);
    if (tagNd(Ac) != cm.tF(cm_mod) && tagNd(Ac) != bTag) {
      gE(a) = 0;
      tagNd(Ac) = 0;
      for (int i = 0; i < gNsf.nrows(); i++) {
        gNsf(i,a) = 0.0;
      }
    }
  }

  // Now that all the elements have been found, data is interpolated
  // from the source to the target mesh
  //
  #ifdef debug_interp
  dmsg << "Interpolate from source to target mesh ... " << "";
  #endif
  Array<double> tmpX(lDof,nNo);

  for (int a = 0; a < nNo; a++) {
    int Ac = gN(a);

    if (tagNd(Ac) == cm.tF(cm_mod)) {
      int Ec = gE(a);
      for (int i = 0; i < gNsf.nrows(); i++) {
        Nsf(i) = gNsf(i,a);
      }

      for (int i = 0; i < eNoN; i++) {
        int Bc = msh[iM].IEN(i,Ec);
        Bc = msh[iM].lN(Bc);
        for (int j = 0; j < tmpX.nrows(); j++) {
          tmpX(j,a) = tmpX(j,a) + Nsf(i)*sD(j,Bc);
        }
      }
    }
  }

  // Since there is no direct mapping for face data, we use L2 norm
  // to find the nearest face node and copy its solution. This requires
  // face node/IEN structure to NOT be changed during remeshing.
  //
  tmpL.resize(nNo);
  bool flag = false; 

  for (int a = 0; a < nNo; a++) {
    int Ac = gN(a);

    if (srfNds(a) != 0) {      // srfNds is a bool (1|0) vector.
      flag = false; 

      for (int iFa = 0; iFa < msh[iM].nFa; iFa++) {
        auto& fa = msh[iM].fa[iFa];

        for (int b = 0; b <fa.nNo; b++) {
          int Bc = fa.gN(b);
          double dS = 0.0; 

          for (int i = 0; i < nsd; i++) {
            double sum = com_mod.x(i,Bc) + Dg(i,Bc) - tMsh.x(i,Ac);
            dS += sum * sum;
          }

          dS = sqrt(dS);

          if (dS < 1.E-12) {
            tmpL(a) = Bc;
            flag = true;
            break;
          }
        }

        if (flag) break;
      }

      if (flag) {
        int Bc = msh[iM].lN(tmpL(a));
        for (int i = 0; i < tmpX.nrows(); i++) { 
          tmpX(i,a) = sD(i,Bc);
        }
      } else { 
        tagNd(Ac) = 0;
      }
    }
  }

  // Map the tagged nodes and solution to local vector within a proc,
  // including boundary nodes. Since the boundary nodes can be overlapping
  // across different procs, these are repeated. But this will not cause
  // problem as the solution is simply overwritten depending on the face pointer.
  //
  #ifdef debug_interp
  dmsg << "Map the tagged nodes and solution ... " << "";
  #endif
  nn = 0;

  for (int a = 0; a < nNo; a++) {
    int Ac = gN(a);
    if (tagNd(Ac) > 0) nn = nn + 1;
  }
  tmpL.resize(nn);
  nn = 0;

  for (int a = 0; a < nNo; a++) {
    int Ac = gN(a);
    if (tagNd(Ac) > 0) {
      tmpL(nn) = a;
      nn = nn + 1;
    }
  }

  gNsf.resize(lDof,nn);

  for (int i = 0; i < nn; i++) {
    int a = tmpL(i);
    int Ac = gN(a);
    for (int j = 0; j < tmpX.nrows(); j++) { 
      gNsf(j,i) = tmpX(j,a);
    }
    tmpL(i) = Ac;
  }

  Vector<int> disp;
  if (cm.mas(cm_mod)) {
    disp.resize(cm.np());
  } 

  #ifdef debug_interp
  dmsg << "MPI_Gather ... " << "";
  #endif
  MPI_Gather(&nn, 1, cm_mod::mpint, disp.data(), 1, cm_mod::mpint, cm_mod.master, cm.com());

  Vector<int> sCount;
  Vector<double> gvec;

  if (cm.mas(cm_mod)) {
    int i = disp.sum();
    i = i * (1 + lDof);
    
    #ifdef debug_interp
    dmsg << "cm.mas ... " << "";
    dmsg << "disp(1): " << disp(0);
    dmsg << "i: " << i;
    #endif
    gvec.resize(i);
    sCount.resize(cm.np());
    for (int i = 0; i < cm.np(); i++) {
      sCount(i) = disp(i)*(1+lDof);
      //dmsg << "sCount(i): " << sCount(i);
    }
    disp(0) = 0;
    for (int i = 1; i < cm.np(); i++) {
      disp(i) = disp(i-1) + sCount(i-1);
      //dmsg << "disp(i): " << disp(i);
    }
  } 

  Vector<double> vec((lDof+1)*nn);
  int e = 0;

  for (int a = 0; a < nn; a++) {
    vec(e) = static_cast<double>(tmpL(a));
    e = e + 1;

    for (int b = 0; b < lDof; b++) {
      vec(e) = gNsf(b,a);
      e = e + 1;
    }
  }

  #ifdef debug_interp
  dmsg << "" << "";
  dmsg << "MPI_Gatherv ... " << "";
  dmsg << "nn: " << nn;
  dmsg << "(1+lDof)*nn: " << (1+lDof)*nn;
  dmsg << "vec.size(): " << vec.size();
  dmsg << "gvec.size(): " << gvec.size();
  dmsg << "sCount.size(): " << sCount.size();
  dmsg << "disp.size(): " << disp.size();
  dmsg << "cm_mod.master: " << cm_mod.master;
  int csize;
  MPI_Comm_size( cm.com(), &csize);
  dmsg << "csize: " << csize;
  #endif

  // Need to be careful here not to send null data. 
  if (disp.size() > -1) {
  MPI_Gatherv(vec.data(), (1+lDof)*nn, cm_mod::mpreal, gvec.data(), sCount.data(), disp.data(), 
      cm_mod::mpreal, cm_mod.master, cm.com());
  } else {
    gvec = vec;
  }

  #ifdef debug_interp
  dmsg << "done MPI_Gatherv " << "";
  dmsg << "" << "";
  #endif

  if (cm.mas(cm_mod)) {
    #ifdef debug_interp
    dmsg << "" << "";
    dmsg << "Set tgD ... " << "";
    #endif
    tgD = 0.0;
    nn = sCount.sum();
    i = 0;

    while (true) { 
      int Ac = round(gvec(i));
      i = i + 1;

      for (int b = 0; b < lDof; b++) {
        tgD(b,Ac) = gvec(i);
        i = i + 1;
      }
      if (i == nn) {
        break;
      }
    }
  }
}

/// @brief Reproduces 'SUBROUTINE INTMSHSRF(lM, lFa)'
//
void int_msh_srf(ComMod& com_mod, CmMod& cm_mod, mshType& lM, faceType& lFa)
{
  auto& cm = com_mod.cm;
  #define n_debug_int_msh_srf
  #ifdef debug_int_msh_srf 
  DebugMsg dmsg(__func__, cm.idcm());
  dmsg.banner();
  #endif

  auto& rmsh = com_mod.rmsh;

  int eNoNb = lM.fa[0].eNoN;
  lFa.eNoN = eNoNb;

  lFa.gnEl = 0;
  for (auto& face : lM.fa) {
    lFa.gnEl += face.gnEl;
  }

  lFa.nNo  = 0;
  lFa.name = "old_mesh_surface";
  lFa.nEl = lFa.gnEl;

  lFa.IEN.resize(lFa.eNoN,lFa.gnEl); 
  lFa.gE.resize(lFa.gnEl);

  if (eNoNb != lM.fa[lM.nFa-1].eNoN) {
    //wrn = "    Remeshing not formulated for different face types"
    return;
  }

  if (cm.mas(cm_mod)) {
    for (int iFa = 0; iFa < lM.nFa; iFa++) {
      //dmsg << "---------- iFa " << std::to_string(iFa+1) + " ----------";
      int eoff = 0;
      if (iFa > 0) {
        for (int i = 0; i < iFa; i++) {
          eoff += lM.fa[i].gnEl;
        }
      }
      //dmsg << "eoff: " << eoff;

      for (int e = eoff; e < eoff+lM.fa[iFa].gnEl; e++) {
        lFa.gE(e) = lM.fa[iFa].gebc(0,e-eoff);

        for (int i = 0; i < lFa.eNoN; i++) {
          lFa.IEN(i,e) = lM.fa[iFa].gebc(i+1,e-eoff);
        }
      }
    }

    read_msh_ns::calc_nbc(lM, lFa);

    Vector<int> incNd(lM.gnNo);

    for (int a = 0; a < lFa.nNo; a++) {
      int Ac = lFa.gN(a);
      incNd(Ac) = a;
    }

    for (int e = 0; e < lFa.gnEl; e++) {
      for (int a = 0; a < lFa.eNoN; a++) {
        int Ac = lFa.IEN(a,e);
        lFa.IEN(a,e) = incNd(Ac);
      }
    }
  }

  cm.bcast(cm_mod, &lFa.nNo);

  if (cm.slv(cm_mod)) {
    lFa.gN.resize(lFa.nNo);
   }
}

void remesher_3d(ComMod& com_mod, CmMod& cm_mod, int iM, faceType& lFa, mshType& lM)
{
  using namespace consts;

  #define n_debug_remesher_3d 
  #ifdef debug_remesher_3d
  auto& cm = com_mod.cm;
  DebugMsg dmsg(__func__, cm.idcm());
  dmsg.banner();
  dmsg << "iM: " << iM;
  #endif

  auto& rmsh = com_mod.rmsh;

  std::array<double,3> rparams = {
    rmsh.maxRadRatio,
    rmsh.minDihedAng,
    rmsh.maxEdgeSize(iM)
  };

  #ifdef debug_remesher_3d
  dmsg << "lFa.nNo: " << lFa.nNo; 
  dmsg << "lFa.nEl: " << lFa.nEl; 
  dmsg << "rmsh.maxEdgeSize(iM): " << rmsh.maxEdgeSize(iM);
  #endif

  int iOK = 0;

  if (rmsh.method == MeshGeneratorType::RMSH_TETGEN) {
     remesh3d_tetgen(lFa.nNo, lFa.nEl, lFa.x.data(), lFa.IEN.data(), rparams, &iOK);
  } else { 
     //err = "Unknown remesher choice."
  }

  std::string elem_file_name  = "new-vol-mesh-cpp.ele";
  std::ifstream new_elem_mesh;
  new_elem_mesh.open(elem_file_name);

  // Get the number of remeshed elements.
  //
  lM.gnEl = 0;
  int value;
  std::string line;

  while (std::getline(new_elem_mesh, line)) {
    lM.gnEl += 1;
  }

  new_elem_mesh.clear();
  new_elem_mesh.seekg(0);
  int min_e = 1000000;
  int max_e = -1000000;

  if (rmsh.method == MeshGeneratorType::RMSH_TETGEN) {
    lM.gnEl = lM.gnEl - 1;
  }
  #ifdef debug_remesher_3d
  dmsg << "Number of elements after remesh: " << lM.gnEl;
  #endif
  lM.gIEN.resize(lM.eNoN,lM.gnEl);

  // Read element connectivity.
  //
  bool first_list = true;
  int e = 0;

  while (std::getline(new_elem_mesh, line)) {
    if (line == "" || first_list) {
      first_list = false;
      continue;
    }
    std::istringstream line_input(line);
    std::vector<int> values;
    int i = 0;
    line_input >> value;

    while (line_input >> value) {
      lM.gIEN(i,e) = value;    
      if (value < min_e) min_e = value;
      if (value > max_e) max_e = value;
      i += 1;
    }

    e += 1;
  }

  // Get the number of remeshed nodes.
  //
  std::string node_file_name  = "new-vol-mesh-cpp.node";
  std::ifstream new_node_mesh;
  new_node_mesh.open(node_file_name);

  lM.gnNo = 0;
  while (std::getline(new_node_mesh, line)) {
    lM.gnNo += 1;
  }

  if (rmsh.method == MeshGeneratorType::RMSH_TETGEN) {
    lM.gnNo = lM.gnNo - 1;
  }
 
  #ifdef debug_remesher_3d
  dmsg << "Number of vertices after remesh: " << lM.gnNo;
  #endif
  new_node_mesh.clear();
  new_node_mesh.seekg(0);
  lM.x.resize(com_mod.nsd,lM.gnNo);

  // Read node coordinates.
  //
  first_list = true;
  int n = 0;

  while (std::getline(new_node_mesh, line)) {
    if (line == "" || first_list) {
      first_list = false;
      continue;
    }
    std::istringstream line_input(line);
    int i = 0;
    double value;
    line_input >> value;

    while (line_input >> value) {
      lM.x(i,n) = value;
      i += 1;
    }

    n += 1;
  }

  // Re-orient element connectivity.
  nn::select_ele(com_mod, lM);
}

/// @brief Reproduces Fortran 'SUBROUTINE REMESHRESTART(timeP)'
//
void remesh_restart(Simulation* simulation)
{
  using namespace consts;

  auto& com_mod = simulation->com_mod;
  auto& cm_mod = simulation->cm_mod;
  auto& cm = com_mod.cm;
  auto& chnl_mod = simulation->chnl_mod;

  #define n_debug_remesh_restart
  #ifdef debug_remesh_restart 
  DebugMsg dmsg(__func__, cm.idcm());
  dmsg.banner();
  #endif

  auto& stFileName = com_mod.stFileName;
  auto& rmsh = com_mod.rmsh;

  auto sTmp = stFileName + "_last.bin";
  #ifdef debug_remesh_restart 
  dmsg << "rmsh.rTS: " << rmsh.rTS;
  dmsg << "tDof: " << com_mod.tDof;
  #endif
  std::filesystem::remove(sTmp);

  // Create the file.
  //
  if (cm.mas(cm_mod)) {
    std::ofstream restart_file(sTmp, std::ios::out | std::ios::binary);
    restart_file.close();
  }

  // This call is to block all processors
  cm.bcast(cm_mod, &rmsh.rTS);

  // Write something.
  //
  auto fTmp = stFileName + "_" + std::to_string(rmsh.rTS) + ".bin";
  auto const recLn = com_mod.recLn;
  const bool dFlag = com_mod.dFlag;
  auto& timeP = com_mod.timeP;
 
  #ifdef debug_remesh_restart 
  dmsg << "dFlag: " << dFlag;
  #endif

  std::ofstream restart_file(fTmp, std::ios::out | std::ios::binary | std::ios::app);
  std::streampos write_pos = (cm.tF(cm_mod) - 1) * recLn;
  restart_file.seekp(write_pos);
  output::write_restart_header(com_mod, timeP, restart_file);

  auto& cplBC = com_mod.cplBC;
  restart_file.write((char*)cplBC.xn.data(), cplBC.xn.msize());
  restart_file.write((char*)rmsh.Y0.data(), rmsh.Y0.msize());
  restart_file.write((char*)rmsh.A0.data(), rmsh.A0.msize());

  if (dFlag) {
    restart_file.write((char*)rmsh.D0.data(), rmsh.D0.msize());
  }

  restart_file.close();

  if (cm.mas(cm_mod)) {
    std::string link_cmd = "ln -f " + fTmp + " " + sTmp;
    std::system(link_cmd.c_str());
  }

  auto& x = com_mod.x;
  auto& gtnNo = com_mod.gtnNo;
  const int nsd = com_mod.nsd;
  const int nMsh = com_mod.nMsh;
  const int tDof= com_mod.tDof;
  gtnNo = 0;
  int lDof = 3 * com_mod.tDof;
  mshType tMsh;
  tMsh.dname = "remesh_restart_tMsh";

  Array<double> gtX, gtD, gX, gD;
  Array<double> tempD, tempX, gnD;

  #ifdef debug_remesh_restart 
  dmsg << "rmsh.A0.nrows: " << rmsh.A0.nrows();
  dmsg << "        ncols: " << rmsh.A0.ncols();
  dmsg << "               " << " ";
  dmsg << "com_mod.x nrows: " << x.nrows();
  dmsg << "          ncols: " << x.ncols();
  #endif

  #ifdef debug_remesh_restart 
  dmsg << " " << " ";
  dmsg << "Remesh ... " << " ";
  #endif

  for (int iM = 0; iM < nMsh; iM++) {
    auto& msh = com_mod.msh[iM];
    #ifdef debug_remesh_restart 
    dmsg << "---------- iM " << iM;
    dmsg << "msh.gnNo: " << msh.gnNo;
    dmsg << "msh.nNo: " << msh.nNo;
    dmsg << "msh.gnEl: " << msh.gnEl;
    dmsg << "msh.nEl: " << msh.nEl;
    dmsg << "rmsh.flag: " << rmsh.flag[iM];
    #endif

    if (rmsh.flag[iM]) {
      if (rmsh.method == MeshGeneratorType::RMSH_TETGEN) {
         //std = " Remeshing <"//CLR(TRIM(msh(iM).name))// "> using <"//CLR("Tetgen")//"> library at time "// STR(rmsh.rTS)
      } else { 
         //err = "Unexpected behavior in Remesher"
      }

      tempX.resize(nsd,msh.nNo);
      gD.resize(lDof,msh.nNo);

      for (int a = 0; a < msh.nNo; a++) {
        int Ac = msh.gN(a);
        for (int i = 0; i < nsd; i++) {
          tempX(i,a) = x(i,Ac) + rmsh.D0(i+nsd+1, Ac);
        }
   
        for (int i = 0; i < tDof; i++) {
          gD(i,a) = rmsh.A0(i,Ac);
          gD(i+tDof,a) = rmsh.Y0(i,Ac);
          gD(2*tDof+i,a) = rmsh.D0(i,Ac);
        }
      }

      tMsh.nFa = 1;
      tMsh.fa.resize(tMsh.nFa);

      if (cm.mas(cm_mod)) {
        tMsh.gnNo = msh.gnNo;
        tMsh.gnEl = msh.gnEl;
        tMsh.eNoN = msh.eNoN;
      } else {
        tMsh.gnNo = 0;
        tMsh.gnEl = 0;
        tMsh.eNoN = 0;
      }

      gX.resize(nsd,tMsh.gnNo);
      gX = all_fun::global(com_mod, cm_mod, msh, tempX);

      if (cm.mas(cm_mod)) {
        tMsh.gIEN.resize(tMsh.eNoN,tMsh.gnEl);
        tMsh.gIEN = msh.gIEN;

        for (int e = 0; e < tMsh.gnEl; e++) {
          int Ec = msh.otnIEN(e);
          for (int a = 0; a < tMsh.eNoN; a++) {
            msh.gIEN(a,e) = tMsh.gIEN(a,Ec);
          }
        }
      }

      int_msh_srf(com_mod, cm_mod, msh, tMsh.fa[0]);

      if (cm.mas(cm_mod)) {
        tMsh.fa[0].x.resize(nsd,tMsh.fa[0].nNo);

        for (int a = 0; a < tMsh.fa[0].nNo; a++) {
          int Ac = tMsh.fa[0].gN(a);

          for (int i = 0; i < nsd; i++) {
            tMsh.fa[0].x(i,a) = gX(i,Ac);
          }
        }

        if (nsd == 2) {
          throw std::runtime_error("Remesher not yet developed for 2D objects.");
        } else {
          remesher_3d(com_mod, cm_mod, iM, tMsh.fa[0], tMsh);
        }

        gnD.resize(lDof,tMsh.gnNo);
      } 

      cm.bcast(cm_mod, &tMsh.gnNo);
      cm.bcast(cm_mod, &tMsh.eNoN);
      cm.bcast(cm_mod, &tMsh.gnEl);
      #ifdef debug_remesh_restart 
      dmsg << "  " << " ";
      dmsg << "Bcast " << " ...";
      dmsg << "tMsh.gnNo: " << tMsh.gnNo;
      dmsg << "tMsh.gnEl: " << tMsh.gnEl;
      #endif

      if (cm.slv(cm_mod)) {
        tMsh.x.resize(nsd,tMsh.gnNo);
        tMsh.gIEN.resize(tMsh.eNoN,tMsh.gnEl);
      }

      MPI_Bcast(tMsh.x.data(), nsd*tMsh.gnNo, cm_mod::mpreal, cm_mod.master, cm.com());

      MPI_Bcast(tMsh.gIEN.data(), tMsh.eNoN*tMsh.gnEl, cm_mod::mpint, cm_mod.master, cm.com());

      MPI_Bcast(tMsh.fa[0].IEN.data(), tMsh.fa[0].eNoN * tMsh.fa[0].nEl, cm_mod::mpint, cm_mod.master, cm.com());

      MPI_Bcast(tMsh.fa[0].gN.data(), tMsh.fa[0].nNo, cm_mod::mpint, cm_mod.master, cm.com());

      interp(com_mod, cm_mod, lDof, iM, tMsh, gD, gnD);

      if (cm.mas(cm_mod)) {
        msh.gnNo = tMsh.gnNo;
        int a = gtnNo + msh.gnNo;

        if (iM > 0) {
          tempX.resize(nsd,gtnNo); 
          tempD.resize(lDof,gtnNo);
          tempX = gtX;
          tempD = gtD;

          gtX.resize(nsd,a); 
          gtD.resize(lDof,a);

          for (int i = 0; i < lDof; i++) {
            for (int j = 0; j < gtnNo; j++) {
              gtD(i,j) = tempD(i,j);
            }
          }

          for (int i = 0; i < nsd; i++) {
            for (int j = 0; j < gtnNo; j++) {
              gtX(i,j) = tempX(i,j);
            }
          }
        } else {
          gtX.resize(nsd,a);
          gtD.resize(lDof,a);
        }

        msh.x.resize(nsd, msh.gnNo);

        for (int i = 0; i < lDof; i++) {
          for (int j = 0; j < tMsh.gnNo; j++) {
            gtD(i,j+gtnNo) = gnD(i,j);
          }
        }

        for (int i = 0; i < nsd; i++) {
          for (int j = 0; j < tMsh.gnNo; j++) {
            gtX(i,j+gtnNo) = tMsh.x(i,j) - gnD(i+2*tDof+nsd+1, j);
            msh.x(i,j) = gtX(i,j+gtnNo);
          }
        }

        gtnNo = a;

        set_face_ebc(com_mod, cm_mod, tMsh.fa[0], tMsh);

        msh.eNoN = tMsh.eNoN;
        msh.gnEl = tMsh.gnEl;
        #ifdef debug_remesh_restart 
        dmsg << "msh.eNoN: " << msh.eNoN;
        dmsg << "msh.gnEl: " << msh.gnEl;
        #endif

        msh.gIEN.resize(msh.eNoN,msh.gnEl);
        for (int i = 0; i < msh.eNoN; i++) {
          for (int j = 0; j < msh.gnEl; j++) {
            msh.gIEN(i,j) = tMsh.gIEN(i,j);
          }
        }

        dist_msh_srf(com_mod, chnl_mod, tMsh.fa[0], msh, 1);

        sTmp = chnl_mod.appPath + "/" + ".remesh_tmp.dir";
        fTmp = sTmp + "/" + msh.name +  "_" + std::to_string(rmsh.rTS) + ".vtu";
        vtk_xml::write_vtu(com_mod, msh, fTmp);
      } 

      MPI_Barrier(cm.com());

    } else {

      tempX.resize(nsd, msh.nNo);
      tempD.resize(lDof, msh.nNo);

      for (int a = 0; a < msh.nNo; a++) {
        int Ac = msh.gN(a);

        for (int i = 0; i < nsd; i++) {
          tempX(i,a) = x(i,Ac);
        }

        for (int i = 0; i < tDof; i++) {
          tempD(i,a) = rmsh.A0(i,Ac);
          tempD(i+tDof,a) = rmsh.Y0(i,Ac);
          tempD(2*tDof+i,a) = rmsh.D0(i,Ac);
        }
      }

      int a;

      if (cm.mas(cm_mod)) {
        a = msh.gnNo;
      } else {
        a = 0;
      }

      gX.resize(nsd,a); 
      gD.resize(lDof,a);

      gX = all_fun::global(com_mod, cm_mod, msh, tempX);
      gD = all_fun::global(com_mod, cm_mod, msh, tempD);

      if (cm.mas(cm_mod)) {
        tMsh.gIEN = msh.gIEN;

        for (int e = 0; e < msh.gnEl; e++) { 
          int Ec = msh.otnIEN(e);

          for (int i = 0; i < msh.eNoN; i++) { 
            msh.gIEN(i,e) = tMsh.gIEN(i,Ec);
          }
        }
      }

      tMsh.nFa = 1;
      tMsh.fa.resize(tMsh.nFa);

      int_msh_srf(com_mod, cm_mod, msh, tMsh.fa[0]); 

      if (cm.mas(cm_mod)) {
        tMsh.fa[0].x.resize(nsd, tMsh.fa[0].nNo);

        for (int i = 0; i < tMsh.fa[0].nNo; i++) {
          int Ac = tMsh.fa[0].gN(i);
          for (int j = 0; j < nsd; j++) {
            tMsh.fa[0].x(j,i) = gX(j,Ac);
          }
        }
        msh.x = gX;

        dist_msh_srf(com_mod, chnl_mod, tMsh.fa[0], msh, 2);

        a = gtnNo + a;

        if (iM > 0) {
          tempX = gtX;
          tempD = gtD;

          gtX.resize(nsd,a); 
          gtD.resize(lDof,a);

          for (int i = 0; i < nsd; i++) {
            for (int j = 0; j < gtnNo; j++) {
              gtX(i,j) = tempX(i,j);
            }
          }

          for (int i = 0; i < lDof; i++) {
            for (int j = 0; j < gtnNo; j++) {
              gtD(i,j) = tempD(i,j);
            }
          }

        } else {
          gtX.resize(nsd,a); 
          gtD.resize(lDof,a);
        }

        #ifdef debug_remesh_restart 
        dmsg << "Set gtX(i,j) " << " ...";
        dmsg << "gtX.nrows: " << gtX.nrows();
        dmsg << "gtX.ncols: " << gtX.ncols();
        dmsg << "gX.nrows: " << gX.nrows();
        dmsg << "gX.ncols: " << gX.ncols();
        dmsg << "a: " << a;
        dmsg << "gtnNo: " << gtnNo;
        #endif

        for (int i = 0; i < nsd; i++) {
          for (int j = 0; j < gX.ncols(); j++) {
            gtX(i,j+gtnNo) = gX(i,j);
          }
        }

        for (int i = 0; i < lDof; i++) {
          for (int j = 0; j < gD.ncols(); j++) {
            gtD(i,j+gtnNo) = gD(i,j);
          }
        }
        gtnNo = a;

      } 
    } // reMesh flag
  }

  cm.bcast(cm_mod, &gtnNo);
  #ifdef debug_remesh_restart 
  dmsg << "gtnNo: " << gtnNo;
  #endif

  if (cm.mas(cm_mod)) {
    com_mod.x.resize(nsd,gtnNo);
    rmsh.A0.resize(tDof,gtnNo);
    rmsh.Y0.resize(tDof,gtnNo);
    rmsh.D0.resize(tDof,gtnNo);
    x = gtX;

    for (int a = 0; a < gtnNo; a++) { 
      for (int i = 0; i < tDof; i++) { 
        rmsh.A0(i,a) = gtD(i,a);
        rmsh.Y0(i,a) = gtD(i+tDof,a);
        rmsh.D0(i,a) = gtD(i+2*tDof,a);
      }
    }

  } else {
    rmsh.A0.resize(0,0);
    rmsh.Y0.resize(0,0); 
    rmsh.D0.resize(0,0);
  }

  // We better free all of these just to be safe. 
  //
  for (int iM = 0; iM < nMsh; iM++) {
    auto& msh = com_mod.msh[iM];

    if (cm.mas(cm_mod)) {
      msh.eDist.clear();
      msh.eId.clear();
      msh.gN.clear();
      msh.gpN.clear();
      msh.IEN.clear();
      msh.otnIEN.clear();
      msh.INN.clear();
      msh.eIEN.clear();
      msh.sbc.clear();
      msh.iGC.clear();
      msh.nW.clear();
      msh.w.clear();
      msh.xi.clear();
      msh.xib.clear();
      msh.x.clear();
      msh.N.clear();
      msh.Nb.clear();
      msh.nV.clear();
      msh.fN.clear();
      msh.Nx.clear();
      msh.Nxx.clear();

      // Free all of the object member data.
      msh.nAdj.destroy();
      msh.eAdj.destroy();

      msh.fs.clear();
    } 
  }

  // [TODO:DaveP] i don't see where com_mod.msh is reallocated, just for master i guess..
  //
  if (cm.slv(cm_mod)) {
    com_mod.msh.clear();
  }

  // Free eq.
  //
  com_mod.eq.clear();
  com_mod.colPtr.clear();
  com_mod.dmnId.clear();
  com_mod.ltg.clear();
  com_mod.rowPtr.clear();
  com_mod.idMap.clear();
  com_mod.cmmBdry.clear();
  com_mod.iblank.clear();
  com_mod.Ao.clear();
  com_mod.An.clear();
  com_mod.Do.clear();
  com_mod.Dn.clear();
  com_mod.R.clear();
  com_mod.Val.clear();
  com_mod.Yo.clear();
  com_mod.Yn.clear();
  com_mod.Bf.clear();

  cplBC.nFa = 0;

  // Additional physics based variables to be deallocated
  com_mod.Ad.clear();
  com_mod.Rd.clear();
  com_mod.Kd.clear();
  com_mod.pS0.clear();
  com_mod.pSn.clear();
  com_mod.pSa.clear();
}

void set_face_ebc(ComMod& com_mod, CmMod& cm_mod, faceType& lFa, mshType& lM)
{
  const int nsd = com_mod.nsd;

  Vector<int> nAssocEl(lM.gnNo);

  for (int e = 0; e < lM.gnEl; e++) {
    for (int a = 0; a < lM.eNoN; a++) {
      int Ac = lM.gIEN(a,e);
      nAssocEl(Ac) = nAssocEl(Ac) + 1;
    }
  }

  int maxAssocEl = nAssocEl.max();
  Array<int> assocEl(maxAssocEl, lM.gnNo);
  nAssocEl = 0;

  for (int e = 0; e < lM.gnEl; e++) {
    for (int a = 0; a < lM.eNoN; a++) {
      int Ac = lM.gIEN(a,e);
      assocEl(nAssocEl(Ac),Ac) = e;
      nAssocEl(Ac) = nAssocEl(Ac) + 1;
    }
  }

  lFa.gE.resize(lFa.nEl);

  // 'bin(:,0) is an element ID
  // 'bin(:,1) is a counter.
  //
  Array<int> bin(maxAssocEl,2);

  for (int e = 0; e < lFa.gnEl; e++) {
    for (int j = 0; j < maxAssocEl; j++) {
      bin(j,0) = -1;
      bin(j,1) = 0;
    }

    for (int a = 0; a < lFa.eNoN; a++) {
      int Ac = lFa.IEN(a,e);

      if (a == 0) {
        for (int i = 0; i < nAssocEl(Ac); i++) {
          bin(i,0) = assocEl(i,Ac);
          bin(i,1) = 1;
        }
      } else {
        for (int i = 0; i < nAssocEl(Ac); i++) {
          for (int j = 0; j < maxAssocEl; j++) {
            if (bin(j,0) == -1) {
              bin(j,0) = assocEl(i,Ac);
              bin(j,1) = 1;
              break;
            } else if (bin(j,0) == assocEl(i,Ac)) {
              bin(j,1) = bin(j,1) + 1;
            }
          } 
        } 
      }
    } 

    for (int j = 0; j < maxAssocEl; j++) {
      if (bin(j,1) == lFa.eNoN) {
        lFa.gE(e) = bin(j,0);
        break;
      }
    } 
  } 

  Array<double> xl(nsd,lM.eNoN), v(nsd,lM.eNoN);

  // Check mesh quality and reset IEN if necessary
  //
  int e = lFa.gE(0);

  //read_msh_ns::check_tet_conn(lM);

  for (int i = 0; i < nsd; i++) {
    for (int j = 0; j < lM.eNoN; j++) {
      xl(i,j) = lM.x(i,lM.gIEN(j,e));
    }
  }
  //xl = lM.x(:,lM.gIEN(:,e))

  for (int i = 0; i < nsd; i++) {
    v(i,0) = xl(i,1) - xl(i,0);
    v(i,1) = xl(i,2) - xl(i,1);
    v(i,2) = xl(i,3) - xl(i,2);
  }

  Array<double> v12(3,2);
  for (int i = 0; i < nsd; i++) {
    v12(i,0) = v(i,0); 
    v12(i,1) = v(i,1); 
  }

  v.set_col(3, utils::cross(v12));

  double sum = 0.0;
  for (int i = 0; i < nsd; i++) {
    sum += v(i,2) * v(i,3);
  }
  int sn = utils::sign(sum);
  //sn = SGN(SUM(v(:,3)*v(:,4)))

  if (sn == 1) {
    int a=0, b=1;

    for (int e = 0; e < lM.gnEl; e++) {
      int Ac = lM.gIEN(a,e);
      lM.gIEN(a,e) = lM.gIEN(b,e);
      lM.gIEN(b,e) = Ac;
    }

  } else if (sn == 0) {
    throw std::runtime_error("[set_face_ebc] Surface element " + std::to_string(e+1) + " is distorted.");
  }

  for (int e = 0; e < lM.gnEl; e++) {
    for (int a = 0; a < lM.eNoN; a++) {
      int Ac = lM.gIEN(a,e);
      for (int i = 0; i < nsd; i++) {
        xl(i,a) = lM.x(i,Ac);
      }
    }

    double Jac = all_fun::jacobian(com_mod, nsd, lM.eNoN, xl, lM.Nx.rslice(0));

    if (Jac < 0.0) {
      throw std::runtime_error("[set_face_ebc] Remeshing didn't improve mesh quality.");
    }
  }
}


};
