// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// The functions defined here replicate the Fortran functions defined in BAFINI.f.

#include "baf_ini.h"

#include "all_fun.h"
#include "consts.h"
#include "fs.h"
#include "nn.h"
#include "set_bc.h"
#include "utils.h"
#include "svZeroD_subroutines.h"

#include "fsils_api.hpp"
#include "fils_struct.hpp"

#include <numeric>
#include <set>
#include <math.h>

namespace baf_ini_ns {

/// @brief This routine initializes required structure for boundaries,
/// faces, those that interface with FSILS and cplBC.
///
/// Modifies:
/// \code {.cpp}
///  com_mod.cplBC.fa
///  com_mod.cplBC.xn
///
///  com_mod.cplBC.fa[i].RCR.Rp = bc.RCR.Rp;
///  com_mod.cplBC.fa[i].RCR.C  = bc.RCR.C;
///  com_mod.cplBC.fa[i].RCR.Rd = bc.RCR.Rd;
///  com_mod.cplBC.fa[i].RCR.Pd = bc.RCR.Pd;
///  com_mod.cplBC.fa[i].RCR.Xo = bc.RCR.Xo;
/// \endcode
///
/// Replicates 'SUBROUTINE BAFINI()' defined in BAFINIT.f
//
void baf_ini(Simulation* simulation)
{
  using namespace consts;
  using namespace fsi_linear_solver;

  auto& com_mod = simulation->com_mod;
  auto& cm_mod = simulation->cm_mod;
  auto& cm = com_mod.cm;
  const int nsd = com_mod.nsd;
  
  #define n_debug_baf_ini
  #ifdef debug_baf_ini
  DebugMsg dmsg(__func__, cm.idcm());
  dmsg.banner();
  #endif

  // Compute face normals and area
  //
  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    #ifdef debug_baf_ini
    dmsg << "iM: " << iM;
    #endif
    auto& msh = com_mod.msh[iM];
    for (int iFa = 0; iFa < msh.nFa; iFa++) {
      if (msh.lFib) {
        continue;
      }
      auto& face = msh.fa[iFa]; 
      face_ini(simulation, msh, face);
    }
    if (msh.lShl) {
      shl_ini(com_mod, cm_mod, com_mod.msh[iM]);
    }
  }

  // Initialize face BC profile
  //
  for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
    auto& eq = com_mod.eq[iEq];
    for (int iBc = 0; iBc < eq.nBc; iBc++) {
      auto& bc = eq.bc[iBc];
      int iFa = bc.iFa;
      int iM = bc.iM;
      bc_ini(com_mod, cm_mod, bc, com_mod.msh[iM].fa[iFa]);

      if (com_mod.msh[iM].lShl) {
        shl_bc_ini(com_mod, cm_mod, bc, com_mod.msh[iM].fa[iFa], com_mod.msh[iM]);
      }
    }
  }

  // cplBC faces are initialized here
  //
  int iEq = 0;
  com_mod.cplBC.fa.resize(com_mod.cplBC.nFa); 
  com_mod.cplBC.xn.resize(com_mod.cplBC.nX);
  
  // Assign cplBC internal variables
  if (com_mod.cplBC.coupled) {
    auto& eq = com_mod.eq[iEq];
    for (int iBc = 0; iBc < eq.nBc; iBc++) {
      auto& bc = eq.bc[iBc];
      int iFa = bc.iFa;
      int iM  = bc.iM;

      if (utils::btest(bc.bType,iBC_cpl) || utils::btest(bc.bType,iBC_RCR)) {
        int i = bc.cplBCptr;
        com_mod.cplBC.fa[i].name = com_mod.msh[iM].fa[iFa].name;
        com_mod.cplBC.fa[i].y = 0.0;

        if (utils::btest(bc.bType, iBC_Dir)) {
          com_mod.cplBC.fa[i].bGrp = CplBCType::cplBC_Dir;

        } else if (utils::btest(bc.bType, iBC_Neu)) {
          com_mod.cplBC.fa[i].bGrp = CplBCType::cplBC_Neu;
          // For implicit or semi-implicit (not explicit) Neumann 0D coupling scheme, 
          // set bType to resistance
          if (com_mod.cplBC.schm != CplBCType::cplBC_E) {
            bc.bType= utils::ibset(bc.bType, iBC_res);
          }

          // Copy RCR structure from bc() to cplBC()
          com_mod.cplBC.fa[i].RCR.Rp = bc.RCR.Rp;
          com_mod.cplBC.fa[i].RCR.C  = bc.RCR.C;
          com_mod.cplBC.fa[i].RCR.Rd = bc.RCR.Rd;
          com_mod.cplBC.fa[i].RCR.Pd = bc.RCR.Pd;
          com_mod.cplBC.fa[i].RCR.Xo = bc.RCR.Xo;
        } else { 
          throw std::runtime_error("Not a compatible cplBC_type");
        }
      }
    }

    if (!com_mod.stFileFlag) {
      set_bc::rcr_init(com_mod, cm_mod);
    }

    if (com_mod.cplBC.useGenBC) {
      set_bc::genBC_Integ_X(com_mod, cm_mod, "I");
    }

    if (com_mod.cplBC.useSvZeroD) {
      svZeroD::init_svZeroD(com_mod, cm_mod);
    }

    if (com_mod.cplBC.schm != CplBCType::cplBC_E) {
      set_bc::calc_der_cpl_bc(com_mod, cm_mod);
    }
  }

  // Setting up FSILS
  //
  int lsPtr = -1;
  #ifdef debug_baf_ini
  dmsg << "Setting up FSILS ... ";
  #endif

  for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
    auto& eq = com_mod.eq[iEq];
    for (int iBc = 0; iBc < eq.nBc; iBc++) {
      auto& bc = eq.bc[iBc];
      int iFa = bc.iFa;
      int iM = bc.iM;
      bc.lsPtr = 0;
      fsi_ls_ini(com_mod, cm_mod, bc, com_mod.msh[iM].fa[iFa], lsPtr);
    }
  }

  if (com_mod.mvMsh) {
    int i = 0;
    for (int a = 0; a < com_mod.tnNo; a++) {
      auto& eq = com_mod.eq[0];
      if (all_fun::is_domain(com_mod, eq, a, EquationType::phys_struct) || 
          all_fun::is_domain(com_mod, eq, a, EquationType::phys_ustruct) || 
          all_fun::is_domain(com_mod, eq, a, EquationType::phys_lElas)) {
        i = i + 1;
      }
    }

    Vector<int> gNodes(i);
    i = 0;

    for (int a = 0; a < com_mod.tnNo; a++) {
      auto& eq = com_mod.eq[0];
      if (all_fun::is_domain(com_mod, eq, a, EquationType::phys_struct) || 
          all_fun::is_domain(com_mod, eq, a, EquationType::phys_ustruct) || 
          all_fun::is_domain(com_mod, eq, a, EquationType::phys_lElas)) {
        gNodes(i) = a;
        i = i + 1;
      }
    }

    int lsPtr = com_mod.nFacesLS - 1;
    #ifdef debug_baf_ini
    dmsg << "lsPtr: " << lsPtr;
    dmsg << "i: " << i;
    #endif

    fsils_bc_create(com_mod.lhs, lsPtr, i, nsd, BcType::BC_TYPE_Dir, gNodes); 
  }

  if (com_mod.cmmInit) {
    #ifdef debug_baf_ini
    dmsg << "cmmInit ";
    dmsg << "cmmBdry size: " << com_mod.cmmBdry.size();
    #endif
    int i = std::accumulate(com_mod.cmmBdry.begin(), com_mod.cmmBdry.end(), 0);
    Vector<int> gNodes(i);

    i = 0;
    for (int a = 0; a < com_mod.tnNo; a++) {
      if (com_mod.cmmBdry[a] == 1) {
        gNodes(i) = a;
        i = i + 1;
      }
    }

    lsPtr = com_mod.nFacesLS - 1;
    fsils_bc_create(com_mod.lhs, lsPtr, i, nsd, BcType::BC_TYPE_Dir, gNodes); 
  }
}

//---------
// bc_ini
//---------
//
void bc_ini(const ComMod& com_mod, const CmMod& cm_mod, bcType& lBc, faceType& lFa)
{
  using namespace consts;
  using namespace utils;

  auto& cm = com_mod.cm;
  int nsd = com_mod.nsd;
  int tnNo = com_mod.tnNo;
 
  #define n_debug_bc_ini
  #ifdef debug_bc_ini
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  if (btest(lBc.bType, enum_int(BoundaryConditionType::bType_gen))) {
    if (btest(lBc.bType, enum_int(BoundaryConditionType::bType_Neu)) && lBc.gm.dof != 1) {
      throw std::runtime_error("Only for (int F=1 is accepted for Neu general BCs");
    }
    return;
  }

  if (btest(lBc.bType, enum_int(BoundaryConditionType::bType_Robin))) {
    if (!com_mod.dFlag) {
      throw std::runtime_error("Robin BC can be set for a displacement-based eqn only");
    }
  }

  int iM = lFa.iM;
  int iFa = lBc.iFa;
  lBc.gx.resize(lFa.nNo);
  //if (.NOT.ALLOCATED(lBc.gx)) ALLOCATE(lBc.gx(lFa.nNo))
  #ifdef debug_bc_ini
  dmsg << "iM: " << iM;
  dmsg << "iFa: " << iFa ;
  dmsg << "tnNo: " << tnNo ;

  dmsg << "cm.np(): " << cm.np();
  dmsg << "lBc.bType: " << lBc.bType;
  dmsg << "lFa.nNo: " << lFa.nNo;
  #endif

  Vector<double> s(tnNo);
  Vector<int> sCount(cm.np());
  Vector<int> disp(cm.np());

  // Just a constant value for Flat profile
  if (btest(lBc.bType, iBC_flat)) { 
    for (int a = 0; a < lFa.nNo; a++) {
      int Ac = lFa.gN(a);
      s(Ac) = 1.0;
    }

  // Here is the method that is used for imposing parabolic profile:
  // 1- Find the coordinate of the points on the boundary 2- find unit
  // vector from center to each of points on the boundary: ew
  // 3- maximize ew(i).e where e is the unit vector from current
  // point to the center 4- Use the point i as the diam here
  //
  } else if (btest(lBc.bType, iBC_para)) { 
    Vector<double> center(3);
    for (int i = 0; i < nsd; i++) {
      center(i) = all_fun::integ(com_mod, cm_mod, lFa, com_mod.x, i) / lFa.area;
    }

    // gNodes is one if a node located on the boundary (beside iFa)
    Vector<int> gNodes(tnNo); 
    Array<double> sVl(nsd,lFa.nNo);
    Array<double> sV(nsd,tnNo);

    for (int jFa = 0; jFa < com_mod.msh[iM].nFa; jFa++) {
      if (jFa == iFa) {
        continue;
      }

      for (int a = 0; a < com_mod.msh[iM].fa[jFa].nNo; a++) { 
        int Ac = com_mod.msh[iM].fa[jFa].gN(a);
        gNodes(Ac) = 1;
       }
     }

     // "j" is a counter for the number of nodes that are located on the
     // boundary of lFa and sVl contains the list of their coordinates
     //
     int j = 0;
     sVl = 0.0;

     for (int a = 0; a < lFa.nNo; a++) {
       int Ac = lFa.gN(a);
       if (gNodes(Ac) == 1) {
         sVl.set_col(j, com_mod.x.col(Ac));
         j = j + 1;
       }
     }

     // Getting the length data that is going to be received at each proc
     //
     if (!cm.seq()) {
        int i = j*nsd;
        MPI_Allgather(&i, 1, cm_mod::mpint, sCount.data(), 1, cm_mod::mpint, cm.com());
        disp(0) = 0;
        for (int i = 1; i < cm.np(); i++) {
          disp(i) = disp(i-1) + sCount(i-1);
        }

        MPI_Allgatherv(sVl.data(), j*nsd, cm_mod::mpreal, sV.data(), sCount.data(), disp.data(), cm_mod::mpreal, cm.com());

        j = sCount.sum() / nsd;
     } else {
       for (int i = 0; i < sV.nrows(); i++) {
         for (int k = 0; k < j; k++) {
           sV(i,k) = sVl(i,k);
         }
       }
     }

     if (cm.mas(cm_mod) && (j == 0)) {
       throw std::runtime_error("Face '" + lFa.name + "' has no perimeter.");
     }

     sVl.resize(nsd, j);

     for (int a = 0; a < j; a++) {
       for (int i = 0; i < sV.nrows(); i++) {
         sV(i,a) = sV(i,a) - center(i);
       }
       sVl.set_col(a, sV.col(a) / sqrt(norm(sV.col(a))));
     }

     // "s" is going to keep the ew.e value
     //
     for (int a = 0; a < lFa.nNo; a++) {
       int Ac = lFa.gN(a);
       auto nV = com_mod.x.col(Ac) - center;
       double maxN = norm(nV, sVl.col(0));
       int i = 0;
       for (int b = 1; b < j; b++) {
         double tmp = norm(nV, sVl.col(b));
         if (tmp > maxN) {
           maxN = tmp;
           i = b;
         }
       }
       s(Ac) = 1.0 - norm(nV) / norm(sV.col(i));
     }

  } else if (btest(lBc.bType, enum_int(BoundaryConditionType::bType_ud))) { 
    for (int a = 0; a < lFa.nNo; a++) {
      int Ac = lFa.gN(a);
      s(Ac) = lBc.gx(a);
    }
  }

  // Now correcting the inlet BC for the inlet ring
  //
  if (btest(lBc.bType, enum_int(BoundaryConditionType::bType_zp))) { 
    for (int jFa = 0; jFa < com_mod.msh[iM].nFa; jFa++) {
      if (jFa == iFa) {
        continue; 
      }
      for (int a = 0; a < com_mod.msh[iM].fa[jFa].nNo; a++) {
        int Ac = com_mod.msh[iM].fa[jFa].gN(a);
        s(Ac) = 0.0;
      }
    }
  }

  // Normalizing the profile for flux
  //
  double tmp = 1.0;
  if (btest(lBc.bType, enum_int(BoundaryConditionType::bType_flx))) { 
    tmp = all_fun::integ(com_mod, cm_mod, lFa, s);
    if (is_zero(tmp)) {
      tmp = 1.0;
      throw std::runtime_error("Face '" + lFa.name + "' used for a BC has no non-zero node.");
     }
  }

  for (int a = 0; a < lFa.nNo; a++) {
    int Ac = lFa.gN(a);
    lBc.gx(a) = s(Ac) / tmp;
  }
}

//----------
// face_ini
//----------
//
void face_ini(Simulation* simulation, mshType& lM, faceType& lFa)
{
  using namespace consts;
  auto& com_mod = simulation->com_mod;
  auto& cm = com_mod.cm;
  auto& cm_mod = simulation->cm_mod;
  int nsd = com_mod.nsd;

  #define n_debug_face_ini
  #ifdef debug_face_ini
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "lM.eType: " << lM.eType;
  dmsg << "lM.eNoN: " << lM.eNoN;
  dmsg << "lFa.eType: " << lFa.eType;
  dmsg << "lFa.eNoN: " << lFa.eNoN;
  dmsg << "lFa.nNo: " << lFa.nNo;
  #endif

  // Calculating face area
  //
  Vector<double> sA(com_mod.tnNo);
  sA = 1.0;
  double area = all_fun::integ(com_mod, cm_mod, lFa, sA);
  #ifdef debug_face_ini
  dmsg << "Face '" << lFa.name << "' area: " << area;
  #endif

  if (utils::is_zero(area)) {
     if (cm.mas(cm_mod)) {
      throw std::runtime_error("Face '" + lFa.name + "' has zero area.");
     }
  }
  lFa.area = area;

  // Compute face normals at nodes
  //
  lFa.nV.resize(nsd,lFa.nNo); 
  Array<double> sV(nsd,com_mod.tnNo);

  bool flag = false;

  if (std::set<ElementType>{ElementType::TRI6,ElementType::QUD8,ElementType::QUD9,
      ElementType::TET10,ElementType::HEX20, ElementType::HEX27}.count(lM.eType) != 0) {
    flag = true;
  }

  // For linear elements or NURBS, we simply project element normals to nodes
  //
  #ifdef debug_face_ini
  dmsg << "Flag: " << flag;
  #endif

  // Compute integral of normal vector over surface element
  if (!flag) {
    Vector<double> nV(nsd);
    for (int e = 0; e < lFa.nEl; e++) {

      if (lFa.eType == ElementType::NRB) {
        // [TODO:DaveP] not implemented. 
        // CALL NRBNNXB(lM, lFa, e);
      }

      for (int g = 0; g < lFa.nG; g++) {
        auto Nx = lFa.Nx.slice(g);
        nn::gnnb(com_mod, lFa, e, g, nsd, nsd-1, lFa.eNoN, Nx, nV);

        for (int a = 0; a < lFa.eNoN; a++) { 
          int Ac = lFa.IEN(a,e);
          for (int i = 0; i < sV.nrows(); i++) { 
            sV(i,Ac) = sV(i,Ac) + nV(i)*lFa.N(a,g)*lFa.w(g);
          }
        }
      }
    }

  // For higher order elements, use reduced order basis on mesh to project element normals. 
  // Lumping method is used to project to face corners. Normals at edge nodes are computed by 
  // simple interpolation from reduced basis. Standard lumping using higher order basis could 
  // lead to spurious errors.
  //
  } else {
    fsType fs;
    fs::set_thood_fs(fs, lFa.eType);
    fs::init_fs(fs, nsd, nsd-1);

    Array<double> xl(nsd,lM.eNoN); 
    Vector<int> ptr(lM.eNoN); 
    std::vector<bool> setIt(lM.eNoN);
    Array<double> xXi(nsd,nsd-1);
    Vector<double> nV(nsd); 

    for (int e = 0; e < lFa.nEl; e++) {
      int Ec = lFa.gE(e);
      std::fill(setIt.begin(), setIt.end(), true);

      for (int a = 0; a < lFa.eNoN; a++) {
        int Ac = lFa.IEN(a,e);
        int b = 0;

        for (int ib = 0; ib < lM.eNoN; ib++) {
          b = ib;
          if (setIt[ib]) {
            int Bc = lM.IEN(ib,Ec);
            if (Bc == Ac) {
              break;
            }
          }
        }

        if (b+1 > lM.eNoN) {
          throw std::runtime_error("Could not find matching face node on higher order mesh");
          //CALL STOPSIM()
        }

        ptr(a) = b;
        setIt[b] = false;
      }

      int a = lFa.eNoN;

      for (int b = 0; b < lM.eNoN; b++) {
        if (setIt[b]) {
          ptr(a) = b;
          a = a + 1;
        }
      }

      for (int a = 0; a < lM.eNoN; a++) {
        int Ac = lM.IEN(a,Ec);
        for (int i = 0; i < nsd; i++) {
          xl(i,a) = com_mod.x(i,Ac);
        }

        if (com_mod.mvMsh) {
          for (int i = 0; i < nsd; i++) {
            xl(i,a) = xl(i,a) + com_mod.Do(i+nsd+1,Ac);
          }
        }
      }

      for (int g = 0; g < fs.nG; g++) {
        xXi = 0.0;

        for (int a = 0; a < fs.eNoN; a++) {
          int b = ptr(a);

          for (int i = 0; i < nsd-1; i++) {
            for (int j = 0; j < xXi.nrows(); j++) {
              xXi(j,i) = xXi(j,i) + fs.Nx(i,a,g)*xl(j,b);
            }
          }
        }

        auto nV = utils::cross(xXi);

        int a = ptr(0);
        int b = ptr(lFa.eNoN);
        Vector<double> v(nsd);

        for (int i = 0; i < nsd; i++) {
          v(i)  = xl(i,a) - xl(i,b);
        }

        if (utils::norm(nV,v) < 0.0) {
          nV = -nV;
         }

        for (int a = 0; a < fs.eNoN; a++) {
          int Ac = lFa.IEN(a,e);
          for (int j = 0; j < sV.nrows(); j++) {
            sV(j,Ac) = sV(j,Ac) + fs.w(g)*fs.N(a,g)*nV(j);
          }
        }
      }

      int g;

      if (lFa.eNoN % 2 == 0) {
        g = lFa.eNoN;
      } else {
        g = lFa.eNoN - 1;
        int Ac = lFa.IEN(lFa.eNoN-1,e);

        for (int b = 0; b < fs.eNoN; b++) {
          int Bc = lFa.IEN(b,e);
          for (int i = 0; i < sV.nrows(); i++) {
            sV(i,Ac) = sV(i,Ac) + sV(i,Bc);
          }
        }

        for (int i = 0; i < sV.nrows(); i++) {
          sV(i,Ac) = sV(i,Ac) / static_cast<double>(fs.eNoN);
        }
      }

      for (int a = fs.eNoN; a < g; a++) {
        int b = a - fs.eNoN;
        int Ac = lFa.IEN(a,e);
        int Bc = lFa.IEN(b,e);

        for (int i = 0; i < sV.nrows(); i++) {
          nV(i) = sV(i,Bc);
        }

        if (b == fs.eNoN-1) {
          Bc = lFa.IEN(0,e);
        } else {
          Bc = lFa.IEN(b+1,e);
        }

        for (int i = 0; i < sV.nrows(); i++) {
          sV(i,Ac) = (nV(i) + sV(i,Bc)) * 0.5;
        }
      }
    }
  }

  all_fun::commu(com_mod, sV);
  flag = true;

  for (int a = 0; a < lFa.nNo; a++) {
    int Ac = lFa.gN(a);
    auto sV_col = sV.col(Ac);
    double sln = sqrt(utils::norm(sV_col));

    if (utils::is_zero(sln)) {
      if (flag) {
        throw std::runtime_error("Skipping normal calculation for node " + std::to_string(a) +  " in face '" +  lFa.name + "'.");
        flag = false;
      }
    }

    for (int i = 0; i < lFa.nV.nrows(); i++) {
      lFa.nV(i,a) = sV(i,Ac) / sln;
    }
  }
}

//------------
// fsi_ls_ini
//------------
//
// Modifies:
//   lBc.lsPtr 
//
// Replicates 'SUBROUTINE FSILSINI'.
//
void fsi_ls_ini(ComMod& com_mod, const CmMod& cm_mod, bcType& lBc, const faceType& lFa, int& lsPtr)
{
  using namespace consts;
  using namespace utils;
  using namespace fsi_linear_solver;

  auto& cm = com_mod.cm;
  int nsd = com_mod.nsd;
  int tnNo = com_mod.tnNo;

  #define n_debug_fsi_ls_ini
  #ifdef debug_fsi_ls_ini
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "lsPtr: " << lsPtr;
  #endif

  int iM = lFa.iM;
  int nNo = lFa.nNo;

  // [NOTE] 'nNo' can be zero so we must check for this.
  Array<double> sVl(nsd,nNo); 
  Array<double> sV(nsd,tnNo); 
  Vector<int> gNodes(nNo);

  // Copy mesh node id corresponding to face node id to gNodes
  for (int a= 0; a < nNo; a++) {
    gNodes(a) = lFa.gN(a);
  }

  if (btest(lBc.bType, iBC_Dir)) {
    if (lBc.weakDir) {
      lBc.lsPtr = -1;
    } else {
      lsPtr = lsPtr + 1;
      lBc.lsPtr = lsPtr;
      sVl = 0.0;
      bool eDrn = false;
      for (int i = 0; i < nsd; i++) {
        if (lBc.eDrn(i) != 0) {
          eDrn = true;
          break;
        }
      }

      if (eDrn) {
        sVl = 1.0;
        for (int i = 0; i < nsd; i++) {
          if ((lBc.eDrn(i) != 0) && (sVl.size() != 0)) {
            sVl.set_row(i, 0.0);
          }
        }
      }
      fsils_bc_create(com_mod.lhs, lsPtr, lFa.nNo, nsd, BcType::BC_TYPE_Dir, gNodes, sVl); 
    }

  } else if (btest(lBc.bType, iBC_Neu)) {
    // Compute integral of normal vector over the face (needed for resistance BC/0D-coupling)
    if (btest(lBc.bType, iBC_res)) {
      sV = 0.0;
      for (int e = 0; e < lFa.nEl; e++) {
        if (lFa.eType == ElementType::NRB) {
          // CALL NRBNNXB(msh(iM),lFa,e)
        }
        for (int g = 0; g < lFa.nG; g++) {
          Vector<double> n(nsd);
          auto Nx = lFa.Nx.slice(g);
          nn::gnnb(com_mod, lFa, e, g, nsd, nsd-1, lFa.eNoN, Nx, n);

          for (int a = 0; a < lFa.eNoN; a++) {
            int Ac = lFa.IEN(a,e);
            for (int i = 0; i < nsd; i++) {
              sV(i,Ac) = sV(i,Ac) + lFa.N(a,g)*lFa.w(g)*n(i);
            }
          }
        }
      }

      if (sVl.size() != 0) { 
        for (int a = 0; a < lFa.nNo; a++) {
          int Ac = lFa.gN(a);
          sVl.set_col(a, sV.col(Ac));
        }
      }

      lsPtr = lsPtr + 1;
      lBc.lsPtr = lsPtr;
      
      // Fills lhs.face(i) variables, including val is sVl exists
      fsils_bc_create(com_mod.lhs, lsPtr, lFa.nNo, nsd, BcType::BC_TYPE_Neu, gNodes, sVl); 
    } else {
      lBc.lsPtr = -1;
    }

  } else if (btest(lBc.bType, iBC_trac)) {
    lBc.lsPtr = -1;

  } else if (btest(lBc.bType, iBC_CMM)) {
    lsPtr = lsPtr + 1;
    lBc.lsPtr = lsPtr;

    nNo = 0;
    for (int a = 0; a < lFa.nNo; a++) {
      if (is_zero(lBc.gx(a))) {
        nNo = nNo + 1;
      }
    }

    sVl.resize(nsd,nNo); 
    gNodes.resize(nNo);
    bool eDrn = false; 

    for (int i = 0; i < nsd; i++) {
      if (lBc.eDrn(i) != 0) {
        eDrn = true;
        break;
      }
    }

    if (eDrn) {
      sVl = 1.0;
      for (int i = 0; i < nsd; i++) {
        if ((lBc.eDrn(i) != 0) && (sVl.size() != 0)) {
          sVl.set_row(i, 0.0);
        }
      }
    }

    nNo = 0;
    for (int a = 0; a < lFa.nNo; a++) {
      int Ac = lFa.gN(a);
      if (is_zero(lBc.gx(a))) {
        gNodes(nNo) = Ac;
        nNo = nNo + 1;
      }
    }

    fsils_bc_create(com_mod.lhs, lsPtr, nNo, nsd, BcType::BC_TYPE_Dir, gNodes, sVl); 
  } else {
    throw std::runtime_error("Unxpected bType in FSILSINI");
  }
}

//--------------
// set_shl_xien 
//--------------
// Compute shell extended IEN for triangular elements.
//
void set_shl_xien(Simulation* simulation, mshType& lM)
{
  using namespace consts;

  auto& com_mod = simulation->com_mod;
  int eNoN = lM.eNoN;
  int nEl = lM.nEl;

  std::array<std::array<int,3>,2> ep{1,2,0, 2,0,1};

  Vector<int> incN(eNoN);
  lM.eIEN.resize(eNoN,nEl); 
  lM.sbc.resize(eNoN,nEl);

  lM.eIEN = -1;
  lM.sbc = 0;

  for (int e = 0; e < nEl; e++) {
    for (int a = 0; a < eNoN; a++) {
      int Ac = lM.IEN(ep[0][a], e);
      int Bc = lM.IEN(ep[1][a], e);

      for (int f = 0; f < nEl; f++) {
        if (e == f) {
          continue; 
        }

        incN = 0;
        for (int b = 0; b < eNoN; b++) {
          if ((lM.IEN(b,f) == Ac) || (lM.IEN(b,f) == Bc)) {
            incN(b) = incN(b) + 1;
          }
        }

        if (incN.sum() == 2) {
          for (int b = 0; b < eNoN; b++) {
            if (incN(b) == 0) {
              lM.eIEN(a,e) = lM.IEN(b,f);
              break;
            }
          }
          break;
        }

        if (lM.eIEN(a,e) == -1) {
          lM.sbc(a,e) = utils::ibset(lM.sbc(a,e), enum_int(BoundaryConditionType::bType_free)); 
        }
      }
    }
  }
}

//------------
// shl_bc_ini 
//------------
// Initializing shell boundary condition variables.
//
// Reproduces 'SUBROUTINE SHLBCINI(lBc, lFa, lM)'.
//
void shl_bc_ini(const ComMod& com_mod, const CmMod& cm_mod, bcType& lBc, faceType& lFa, mshType& lM)
{
  using namespace consts;
  using namespace utils;

  auto& cm = com_mod.cm;
  int nsd = com_mod.nsd;
  int tnNo = com_mod.tnNo;

  int task_id = cm.idcm();
  std::string msg_prefix;

  if (lFa.eType == ElementType::NRB) {
    return; 
  }

  for (int e = 0; e < lFa.nEl; e++) {
    int Ec = lFa.gE(e);

    for (int a = 0; a < lM.eNoN; a++) {
      int Ac = lM.IEN(a,Ec);
      bool bFlag = false;

      for (int b = 0; b < lFa.eNoN; b++) {
        int Bc = lFa.IEN(b,e);
        if (Ac == Bc) {
          bFlag = true;
          break;
        }
      }

      if (!bFlag) {
        if (!btest(lM.sbc(a,Ec),iBC_free)) {
          throw std::runtime_error("BC detected on a non-boundary shell element. ");
        }

        lM.sbc(a,Ec) = ibclr(lM.sbc(a,Ec), iBC_free);

        if (btest(lBc.bType, iBC_free)) {
          lM.sbc(a,Ec) = ibset(lM.sbc(a,Ec), iBC_free);

        } else if (btest(lBc.bType, iBC_fix)) {
          lM.sbc(a,Ec) = ibset(lM.sbc(a,Ec), iBC_fix);

        } else if (btest(lBc.bType, iBC_hing)) {
          lM.sbc(a,Ec) = ibset(lM.sbc(a,Ec), iBC_hing);

        } else if (btest(lBc.bType, iBC_symm)) {
          lM.sbc(a,Ec) = ibset(lM.sbc(a,Ec), iBC_symm);
        }
        break; 
      }
    }
  }
}

//---------
// shl_ini
//---------
// Reproduces 'SUBROUTINE SHLINI(lM)' 
//
void shl_ini(const ComMod& com_mod, const CmMod& cm_mod, mshType& lM) 
{
  using namespace consts;
  using namespace utils;

  auto& cm = com_mod.cm;
  int nsd = com_mod.nsd;
  int tnNo = com_mod.tnNo;

  int nNo = lM.nNo;
  int nEl = lM.nEl;
  int eNoN = lM.eNoN;

  // Compute shell director (normal)
  Array<double> xl(nsd,eNoN); 
  Array<double> sV(nsd,tnNo); 
  lM.nV.resize(nsd,nNo);

  double area = 0.0;

  for (int e = 0; e < nEl; e++) {
    if (lM.eType == ElementType::NRB) {
      //CALL NRBNNX(lM, e)
    }

    for (int a = 0; a < eNoN; a++) {
      int Ac = lM.IEN(a,e);
      xl.set_col(a, com_mod.x.col(Ac));
    }

    for (int g = 0; g < lM.nG; g++) {
      Vector<double> nV(nsd);
      Array<double> tmpR(nsd,nsd-1);
      auto Nxi = lM.Nx.slice(g);
      nn::gnns(nsd, eNoN, Nxi, xl, nV, tmpR, tmpR);
      double Jac = sqrt(norm(nV));

      for (int a = 0; a < eNoN; a++) {
        int Ac = lM.IEN(a,e);
        for (int i = 0; i < nsd; i++) {
          sV(i,Ac) = sV(i,Ac) + lM.w(g)*lM.N(a,g)*nV(i);
        }
        area  = area  + lM.w(g)*lM.N(a,g)*Jac;
      }
    }
  }

  area = cm.reduce(cm_mod, area);
  all_fun::commu(com_mod, sV);
  bool flag = true;

  for (int a = 0; a < nNo; a++) {
    int Ac = lM.gN(a);
    double Jac = sqrt(norm(sV.col(Ac)));

    if (is_zero(Jac)) {
      if (flag) {
        throw std::runtime_error("Not a compatible cplBC_type");
        //wrn = "Skipping normal calculation of node "//a// " in mesh <"//TRIM(lM.name)//">"
        flag = false;
      }

      lM.nV.col(a) = 0.0;
      lM.nV(0,a) = 1.0;
      continue; 
    }
    lM.nV.set_col(a, sV.col(Ac) / Jac);
  }
}

};

