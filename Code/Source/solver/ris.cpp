// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "ris.h"

#include "all_fun.h"
#include "lhsa.h"
#include "mat_fun.h"
#include "nn.h"
#include "utils.h"
#include "set_bc.h"

namespace ris {

/// @brief This subroutine computes the mean pressure and flux on the ris surface 
void ris_meanq(ComMod& com_mod, CmMod& cm_mod)
{
  #define n_debug_ris_meanq
  #ifdef debug_ris_meanq
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "risFlag: " << com_mod.risFlag;
  #endif

  using namespace consts;

  auto& cm = com_mod.cm;
  auto& eq = com_mod.eq;
  auto& RIS = com_mod.ris;
  auto& msh = com_mod.msh;

  const int nsd = com_mod.nsd;
  const int cEq = com_mod.cEq;

  auto& An = com_mod.An;
  auto& Ad = com_mod.Ad;
  auto& Dn = com_mod.Dn;
  auto& Yn = com_mod.Yn;

  Array<double> tmpV(maxNSD, com_mod.tnNo);

  // Get the number of projections from RIS
  int nPrj = RIS.nbrRIS;

  int iEq = 0;
  int m = 1;
  int s = eq[iEq].s + nsd;
  int e = s + m - 1;

  for (int j = 0; j < Yn.ncols(); j++) {
      tmpV(0,j) = Yn(s,j);
  }

  for (int iProj = 0; iProj < nPrj; iProj++){
    RIS.meanFl(iProj) = 0.0;
    for (int i = 0; i < 2; i++) { // Assuming two meshes
      RIS.meanP(iProj, i) = 0.0;
      int iM = RIS.lst(i,0,iProj);
      int iFa = RIS.lst(i,1,iProj);
      double tmp = msh[iM].fa[iFa].area;
      RIS.meanP(iProj,i) = all_fun::integ(com_mod, cm_mod, msh[iM].fa[iFa], tmpV, 0)/tmp;
    }
  }

  // For the velocity
  m = nsd; 
  s = eq[iEq].s;
  e = s + m - 1;

  for (int iProj = 0; iProj < nPrj; iProj++) {
    // tmpV[0:m,:] = Yn[s:e,:]; 
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < Yn.ncols(); j++) {
          tmpV(i,j) = Yn(s+i,j);
        }
    }
    int iM = RIS.lst(0,0,iProj);
    int iFa = RIS.lst(0,1,iProj);
    RIS.meanFl(iProj) = all_fun::integ(com_mod, cm_mod, msh[iM].fa[iFa], tmpV, 0, m-1);

    if (cm.mas(cm_mod)) {
      std::cout << "For RIS projection: " << iProj << std::endl;
      std::cout << "    The average pressure is: " << RIS.meanP(iProj,0) << ", " 
                << RIS.meanP(iProj,1) << std::endl;
      std::cout << "    The average flow is: " << RIS.meanFl(iProj) << std::endl;
    }
  }
}

/// @brief  Weak treatment of RIS resistance boundary conditions
void ris_resbc(ComMod& com_mod, const Array<double>& Yg, const Array<double>& Dg) 
{
  using namespace consts;
  #define n_debug_ris_resbc
  #ifdef debug_ris_resbc
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "risFlag: " << com_mod.risFlag;
  #endif

  // auto& eq = com_mod.eq[cEq];
  auto& eq = com_mod.eq;
  auto& RIS = com_mod.ris;
  auto& msh = com_mod.msh;

  const int nsd = com_mod.nsd;
  const int cEq = com_mod.cEq;
  const int cDmn = com_mod.cDmn;
  int nPrj = RIS.nbrRIS;

  bcType lBc;

  for (int iProj = 0; iProj < nPrj; iProj++) {
    if (!RIS.clsFlg[iProj]) {continue;}
    
    // Weak Dirichlet BC for fluid/FSI equations
    lBc.weakDir = true;
    lBc.tauB = RIS.Res[iProj];
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_Dir));
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_std));
    lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_flat));

    lBc.eDrn.resize(nsd);
    lBc.eDrn = 0;
    
    for (int i = 0; i < 2; i++) { // Always two meshes for one projection
      int iM = RIS.lst(i,0,iProj);
      int iFa = RIS.lst(i,1,iProj);
      lBc.gx.resize(msh[iM].fa[iFa].nNo);
      lBc.gx = 0.; //1.0;
      auto& cPhys = eq[cEq].dmn[cDmn].phys;

      double Yg_max, Yg_min, Yg_sum;
      double Dg_max, Dg_min, Dg_sum;
      std::array<fsType,2> fs;
      
      if (cPhys == EquationType::phys_fluid) {
        // Build the correct BC
        set_bc::set_bc_dir_wl(com_mod, lBc, msh[iM], msh[iM].fa[iFa], Yg, Dg);
      }
      lBc.gx.clear();
    }
    lBc.eDrn.clear();
  }
}


void setbc_ris(ComMod& com_mod, const bcType& lBc, const mshType& lM, const faceType& lFa, 
               const Array<double>& Yg, const Array<double>& Dg)
{
  // [HZ] looks not needed in the current implementation
}


/// @brief  This subroutine updates the resistance and activation flag for the 
/// closed and open configurations of the RIS surfaces 
void ris_updater(ComMod& com_mod, CmMod& cm_mod)
{
  #define n_debug_ris_updater
  #ifdef debug_ris_updater
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "risFlag: " << com_mod.risFlag;
  #endif

  // auto& eq = com_mod.eq[cEq];
  auto& cm = com_mod.cm;
  auto& eq = com_mod.eq;
  auto& RIS = com_mod.ris;

  int nPrj = RIS.nbrRIS;

  for (int iProj = 0; iProj < nPrj; iProj++) {
    // The valve is closed check if it should open
    if (RIS.clsFlg[iProj]) {
      // OPENING CONDITION: Check condition on the pressure difference 
      if (RIS.meanP(iProj,0) > RIS.meanP(iProj,1)) {
        RIS.clsFlg[iProj] = false;
        if (cm.mas(cm_mod)) {
          std::cout << "RIS Proj " << iProj << ": Going from close to open." << std::endl;
        }
        RIS.nbrIter(iProj) = 0;
        // I needed to update the state variables when the valve 
        // goes from close to open to prevent the valve goes back
        // to close at the next iteration. This is needed only for
        // close to open and cannot be used for open to close.
        com_mod.Ao = com_mod.An;
        com_mod.Yo = com_mod.Yn;
        if (com_mod.dFlag) {com_mod.Do = com_mod.Dn;}
        com_mod.cplBC.xo = com_mod.cplBC.xn;
      } 
    } else {
      // The valve is open, check if it should close.
      // CLOSING CONDITION: Check existence of a backflow
      if (RIS.meanFl(iProj) < 0.) {
        RIS.clsFlg[iProj] = true;
        if (cm.mas(cm_mod)) {
          std::cout << "RIS Proj " << iProj << ": Going from open to close." << std::endl;
        }
        RIS.nbrIter(iProj) = 0;
      }
    }
  }
}

/// @brief  This subroutine will check the valve status if it is admissible 
/// or not, if not admissible we recompute the iteration until it will be
void ris_status(ComMod& com_mod, CmMod& cm_mod)
{
  using namespace consts;

  #define n_debug_ris_status
  #ifdef debug_ris_status
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "risFlag: " << com_mod.risFlag;
  #endif

  // auto& eq = com_mod.eq[cEq];
  auto& cm = com_mod.cm;
  auto& eq = com_mod.eq;
  auto& RIS = com_mod.ris;

  int nPrj = RIS.nbrRIS;

  for (int iProj = 0; iProj < nPrj; iProj++) {
    RIS.nbrIter(iProj) += 1;
    RIS.status[iProj] = true;
    // If the valve is closed, check the pressure difference, 
    // if the pressure difference is negative the valve should be open
    // the status is then not admissible
    if (RIS.clsFlg[iProj]) {
      if (RIS.meanP(iProj,0) > RIS.meanP(iProj,1)) {
        if (cm.mas(cm_mod)) {
          std::cout << "RIS Proj " << iProj << ": **** Not admissible, it should be open ****" << std::endl;
        }
        RIS.status[iProj] = false;
      }
    } else {
      // If the valve is open, chech the flow, 
      // if the flow is negative the valve should be closed
      // the status is then not admissible
      if (RIS.meanFl(iProj) < 0.) {
        if (cm.mas(cm_mod)) {
          std::cout << "RIS Proj " << iProj << ": **** Not admissible, it should be closed ****" << std::endl;
        }
        RIS.status[iProj] = false;
      }
    }
  }
}

/// @brief This subroutine assembels the element stiffness matrix into the
/// global stiffness matrix (Val sparse matrix formatted as a vector)
void doassem_ris(ComMod& com_mod, const int d, const Vector<int>& eqN, 
     const Array3<double>& lK, const Array<double>& lR)
{ 
  using namespace consts;

  #define n_debug_doassem_ris
  #ifdef debug_doassem_ris
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "risFlag: " << com_mod.risFlag;
  #endif

  auto& cm = com_mod.cm;
  auto& eq = com_mod.eq;
  auto& RIS = com_mod.ris;

  const int nsd = com_mod.nsd;
  int nPrj = RIS.nbrRIS;
  int rowNadj = 0;
  double val_sum;

  for (int iProj = 0; iProj < nPrj; iProj++) {
    if (RIS.clsFlg[iProj]) {continue;}

    for (int a = 0; a < d; a++) {
      int rowN = eqN(a);
      if (rowN == -1) {continue;}

      std::array<int, 2> mapIdx;
      utils::find_loc(com_mod.grisMapList[iProj].map, rowN, mapIdx);

      if (mapIdx[0] == -1) {continue;}

      for (int jM = 0; jM < 2; jM++) {
        if (jM == mapIdx[0]) {continue;}
        if (com_mod.grisMapList[iProj].map(jM, mapIdx[1]) == -1) {continue;}
        rowNadj = com_mod.grisMapList[iProj].map(jM, mapIdx[1]);
      }

      if (rowNadj == -1) {continue;}

      for (int i = 0; i < com_mod.R.nrows(); i++) {
        com_mod.R(i,rowNadj) = com_mod.R(i,rowNadj) + lR(i,a);
      }

      for (int b = 0; b < d; b++) {
        int colN = eqN(b);
        std::array<int, 2> mapIdxC;
        utils::find_loc(com_mod.grisMapList[iProj].map, colN, mapIdxC);

        if (mapIdxC[0] != -1) {
          for (int jM = 0; jM < 2; jM++) {
            if (jM == mapIdxC[0]) {continue;}
            if (com_mod.grisMapList[iProj].map(jM, mapIdxC[1]) == -1) {continue;}
            colN = com_mod.grisMapList[iProj].map(jM, mapIdxC[1]);
          }
        } //else {
        //   continue;
        // }

        if (colN == -1) {continue;}
        int left = com_mod.rowPtr(rowNadj);
        int right = com_mod.rowPtr(rowNadj+1);
        int ptr = (right + left) / 2;

        while (colN != com_mod.colPtr(ptr)) {
          if (colN > com_mod.colPtr(ptr)) {
            left = ptr;
          } else {
            right = ptr;
          }
          ptr = (right + left) / 2;
        }

        for (int i = 0; i < com_mod.Val.nrows(); i++) {
          com_mod.Val(i,ptr) = com_mod.Val(i,ptr) + lK(i,a,b);
        }
      }
    }
  }
}

/// @brief This subroutine assembels the element stiffness matrix into the
/// global stiffness matrix (Val sparse matrix formatted as a vector) 
void doassem_velris(ComMod& com_mod, const int d, const Array<int>& eqN, 
     const Array3<double>& lK, const Array<double>& lR)
{
 // [HZ] looks not needed in the current implementation
}

void clean_r_ris(ComMod& com_mod)
{
  // [HZ] looks not needed in the current implementation
}

void setbcdir_ris(ComMod& com_mod, Array<double>& lA, Array<double>& lY, Array<double>& lD)
{
  // [HZ] looks not needed in the current implementation
}

/// RIS0D code
void ris0d_bc(ComMod& com_mod, CmMod& cm_mod, const Array<double>& Yg, const Array<double>& Dg) 
{
  using namespace consts;

  #define n_debug_ris0d_bc
  #ifdef debug_ris0d_bc
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "ris0DFlag: " << com_mod.ris0DFlag;
  #endif

  auto& cm = com_mod.cm;
  auto& eq = com_mod.eq;
  auto& msh = com_mod.msh;

  const int nsd = com_mod.nsd;
  const int cEq = com_mod.cEq;

  bcType lBc;
  
  for (int iBc = 0; iBc < eq[cEq].nBc; iBc++) {
    int iFa = eq[cEq].bc[iBc].iFa;
    int iM = eq[cEq].bc[iBc].iM;

    if (!utils::btest(eq[cEq].bc[iBc].bType, iBC_Ris0D)) {continue;}

    if (eq[cEq].bc[iBc].clsFlgRis == 1) {
      // Weak Dirichlet BC for fluid/FSI equations
      lBc.weakDir = true;
      lBc.tauB = eq[cEq].bc[iBc].resistance;
      lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_Dir));
      lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_std));
      lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_flat));

      lBc.eDrn.resize(nsd);
      lBc.eDrn = 0;

      // Apply bc Dir 
      lBc.gx.resize(msh[iM].fa[iFa].nNo);
      lBc.gx = 1.0;
      set_bc::set_bc_dir_wl(com_mod, lBc, msh[iM], msh[iM].fa[iFa], Yg, Dg);
      lBc.gx.clear();
      lBc.eDrn.clear();
    } else {
      // Apply Neu bc 
      set_bc::set_bc_neu_l(com_mod, cm_mod, eq[cEq].bc[iBc], msh[iM].fa[iFa], Yg, Dg);
    }

  }

}

void ris0d_status(ComMod& com_mod, CmMod& cm_mod)//, const Array<double>& Yg, const Array<double>& Dg) 
{
  using namespace consts;

  #define n_debug_status
  #ifdef debug_status
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "ris0DFlag: " << com_mod.ris0DFlag;
  #endif

  auto& cm = com_mod.cm;
  auto& eq = com_mod.eq;
  auto& msh = com_mod.msh;

  const int nsd = com_mod.nsd;
  const int cEq = com_mod.cEq;

  auto& An = com_mod.An;
  auto& Ad = com_mod.Ad;
  auto& Dn = com_mod.Dn;
  auto& Yn = com_mod.Yn;

  bcType lBc;
  faceType lFa;

  double meanP = 0.0;
  double meanFl = 0.0;
  double tmp, tmp_new;
  Vector<double> sA;
  Array<double> tmpV;

  for (int iBc = 0; iBc < eq[cEq].nBc; iBc++) {
    int iFa = eq[cEq].bc[iBc].iFa;
    int iM = eq[cEq].bc[iBc].iM;

    if (!utils::btest(eq[cEq].bc[iBc].bType, iBC_Ris0D)) {continue;}

    tmpV.resize(maxNSD, com_mod.tnNo);
    tmpV = 0.0;

    // Compute mean Q and pressure difference 
    int m = 1;
    int s = eq[cEq].s + nsd;
    // int e = s + m - 1;
  
    for (int j = 0; j < Yn.ncols(); j++) {
        tmpV(0,j) = Yn(s,j);
    }

    tmp = msh[iM].fa[iFa].area;
    sA.resize(com_mod.tnNo);
    sA = 1.0;
    lFa = msh[iM].fa[iFa];
    // such update may be not correct
    tmp_new = all_fun::integ(com_mod, cm_mod, lFa, sA);
    meanP = all_fun::integ(com_mod, cm_mod, msh[iM].fa[iFa], tmpV, 0, m-1)/tmp_new;

    // For the velocity
    m = nsd;
    s = eq[cEq].s;
    // e = s + m - 1;
    tmpV.resize(maxNSD,com_mod.tnNo);
    tmpV = 0.0;

    for (int i = 0; i < m; i++) {
      for (int j = 0; j < Yn.ncols(); j++) {
        tmpV(i,j) = Yn(s+i,j);
      }
    }

    meanFl = all_fun::integ(com_mod, cm_mod, msh[iM].fa[iFa], tmpV, 0, m-1);

    std::cout << "The average pressure is: " << meanP << std::endl;
    std::cout << "The pressure from 0D is: " << eq[cEq].bc[iBc].g << std::endl;
    std::cout << "The average flow is: " << meanFl << std::endl;

    com_mod.RisnbrIter = com_mod.RisnbrIter + 1;

    if (com_mod.RisnbrIter < 25 && com_mod.cTS > 0) {
      if (!cm.seq()) {cm.bcast(cm_mod, &com_mod.RisnbrIter);}
      return ;
    }

    // Update RES
    // Update the resistance - determine the configuration 
    // The valve is closed check if it should open
    if (eq[cEq].bc[iBc].clsFlgRis == 1) {
      // OPENING CONDITION: Check condition on the pressure difference 
      if (eq[cEq].bc[iBc].g < meanP) {
        eq[cEq].bc[iBc].clsFlgRis = 0;
        if (!cm.seq()) {
          cm.bcast(cm_mod, &eq[cEq].bc[iBc].clsFlgRis);
        }
        std::cout << "!!! -- Going from close to open " << std::endl;
        com_mod.RisnbrIter = 0;
      }
    } else {
      // The valve is open, check if it should close. 
      // CLOSING CONDITION: Check existence of a backflow
      if (meanFl < 0.) {
        eq[cEq].bc[iBc].clsFlgRis = 1;
        if (!cm.seq()) {
          cm.bcast(cm_mod, &eq[cEq].bc[iBc].clsFlgRis);
        }
        std::cout << "!!! -- Going from open to close " << std::endl;
        com_mod.RisnbrIter = 0;
      }
    }

    // Check for the status
    // If the valve is closed, check the pressure difference, 
    // if the pressure difference is negative the valve should be open
    // -> the status is then not admissible
    if (eq[cEq].bc[iBc].clsFlgRis == 1) {
      if (eq[cEq].bc[iBc].g < meanP) {
        std::cout << "** Not admissible, should be open **" << std::endl;
      }
    } else {
      // If the valve is open, chech the flow, 
      // if the flow is negative the valve should be closed
      // -> the status is then not admissible
      if (meanFl < 0.) {
        std::cout << "** Not admissible, should be closed **" << std::endl;
      }
    }

  }

  if (!cm.seq()) {
    cm.bcast(cm_mod, &com_mod.RisnbrIter);
  }

}

}