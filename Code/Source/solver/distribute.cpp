// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// The code here replicates the Fortran code in DISTRIBUTE.f.

#include "distribute.h"

#include "all_fun.h"
#include "ComMod.h"
#include "consts.h"
#include "nn.h"
#include "utils.h"

#include "CmMod.h"

#include "mpi.h"

#include <iostream>
#include <math.h>

extern "C" {

int split_(int *nElptr, int *eNoNptr, int *eNoNbptr, int *IEN, int *nPartsPtr, int *iElmdist, float *iWgt, int *part);

};
 
/// @brief Partition and distribute data across processors.
///
/// This function replicates the Fortran 'SUBROUTINE DISTRIBUTE' in DISTRIBUTE.f.
//
void distribute(Simulation* simulation)
{
  auto& cm_mod = simulation->cm_mod;
  auto& chnl_mod = simulation->chnl_mod;
  auto& com_mod = simulation->com_mod;
  auto& cm = com_mod.cm;

  #define n_debug_distribute
  #ifdef debug_distribute
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << " com_mod.resetSim: " << com_mod.resetSim;
  #endif

  if (!com_mod.resetSim) {
    cm.bcast(cm_mod, &chnl_mod.pClr);
    cm.bcast(cm_mod, chnl_mod.appPath);
    #ifdef debug_distribute
    dmsg << " chnl_mod.pClr: " << chnl_mod.pClr;
    dmsg << " chnl_mod.appPath: " << chnl_mod.appPath;
    #endif

    //[TODO:Davep] add this?
    //CALL cm%bcast(wrn%oTS)
    //wrn%oTF = wrn%oTS

    //  Constructing data structures one by one
    cm.bcast(cm_mod, &com_mod.nMsh);
    cm.bcast(cm_mod, &com_mod.nsd);
    cm.bcast(cm_mod, &com_mod.rmsh.isReqd);
  } 

  cm.bcast(cm_mod, &com_mod.gtnNo);

  if (cm.slv(cm_mod)) {
    com_mod.msh.resize(com_mod.nMsh);
    #ifdef debug_distribute
    dmsg << "slave process " << "";
    #endif
  } else {
    #ifdef debug_distribute
    dmsg << "master process " << "";
    #endif
  }

  #ifdef debug_distribute
  dmsg << "nMsh: " << com_mod.nMsh;
  dmsg << "gtnNo: " << com_mod.gtnNo;
  dmsg << "nsd: " << com_mod.nsd;
  #endif

  // wgt and wrk are the assigned portion of each mesh to the each
  // processor.
  //
  int nMsh = com_mod.nMsh;
  int num_proc = cm.np(); 
  Array<double> wgt(nMsh, num_proc); 
  Vector<double> wrk(nMsh); 

  // Here is rough estimation of how each mesh should be splited
  // between processors
  // 
  //wrk = REAL(msh%gnNo, KIND=RKIND)/REAL(gtnNo, KIND=RKIND)
  //
  #ifdef debug_distribute
  dmsg << "Rough estimation of how each mesh split ..." << " ";
  #endif

  for (int i = 0; i < com_mod.msh.size(); i++) {
    wrk[i] = static_cast<double>(com_mod.msh[i].gnNo) / com_mod.gtnNo;
    #ifdef debug_distribute
    dmsg << "---------- i " << i;
    dmsg << "msh[i].name: " << com_mod.msh[i].name;
    dmsg << "msh[i].gnNo: " << com_mod.msh[i].gnNo;
    dmsg << "com_mod.gtnNo: " << com_mod.gtnNo;
    dmsg << "msh[i].nNo: " << com_mod.msh[i].nNo;
    dmsg << "msh[i].msh.nEl: " << com_mod.msh[i].nEl;
    dmsg << "wrk[i]: "  << wrk[i];
    #endif
  }

  cm.bcast(cm_mod, wrk);

  int task_id = com_mod.cm.idcm();
  all_fun::split_jobs(task_id, nMsh, num_proc, wgt, wrk);

  #ifdef debug_distribute
  dmsg << " "  << " ";
  dmsg << "Split jobs .. "  << " ";
  dmsg << "wrk: "  << wrk;
  dmsg << "wgt: "  << wgt;
  #endif

  // First partitioning the meshes
  // gmtl:  gtnNo --> tnNo
  //
  #ifdef debug_distribute
  dmsg << "Partition meshes " << " ...";
  #endif
  com_mod.tnNo = 0;

  if (cm.seq()) {
    com_mod.tnNo = com_mod.gtnNo;
    com_mod.ltg.resize(com_mod.tnNo);

    for (int a = 0; a < com_mod.tnNo; a++) {
      com_mod.ltg[a] = a;
    }
  }

  Vector<int> gmtl(com_mod.gtnNo);
  Vector<float> iWgt(num_proc); 

  #ifdef debug_distribute
  dmsg << "          " << " ";
  dmsg << "wgt: " << wgt;
  #endif

  // Need the ris flag here
  cm.bcast(cm_mod, &com_mod.risFlag);
  // Communicating RIS info
  if (com_mod.risFlag) {
    dist_ris(com_mod, cm_mod, cm);
  }
  // Need the uris flag here
  cm.bcast(cm_mod, &com_mod.urisFlag);
  if (com_mod.urisFlag) {
    dist_uris(com_mod, cm_mod, cm);
  }


  int risProc = -1;
  int jM = 0;
  for (int iM = 0; iM < nMsh; iM++) {
    #ifdef debug_distribute
    dmsg << "          " << " ";
    dmsg << "Partitioning mesh: " << iM;
    #endif
    auto sum = wgt.sum_row(iM);
    for (int i = 0; i < num_proc; i++) {
      iWgt[i] = wgt(iM,i) / sum;
    }
    #ifdef debug_distribute
    dmsg << "iWgt: " << iWgt;
    #endif
    
    if (com_mod.risFlag) {
      for (int iProj = 0; iProj < com_mod.ris.nbrRIS; iProj++) {
        // Find the adjacent mesh
        if (com_mod.ris.lst(0,0,iProj) == iM) {
          jM = com_mod.ris.lst(1,0,iProj);
        } else if (com_mod.ris.lst(1,0,iProj) == iM) {
          jM = com_mod.ris.lst(0,0,iProj);
        } else {
          continue;
        }
        // Right now just assign all nodes to the same proc
        risProc = -1;
        for (int e = 0; e < com_mod.msh[jM].gnEl; e++) {
          if (com_mod.msh[jM].eRIS(e)) {
            risProc = com_mod.msh[jM].partRIS(e);
            break;
          }
        }
        if (risProc != -1) {
          for (int e = 0; e < com_mod.msh[iM].gnEl; e++) {
            if (com_mod.msh[iM].eRIS(e)) {
              com_mod.msh[iM].partRIS(e) = risProc;
            }
          }
        }
      }
    }
    part_msh(simulation, iM, com_mod.msh[iM], gmtl, num_proc, iWgt);
  }

  // Setting gtl pointer in case that it is needed and mapping IEN.
  //
  int tnNo = com_mod.tnNo;
  #ifdef debug_distribute
  dmsg << " " << " ";
  dmsg << "Setting gtl pointer " << " ...";
  dmsg << "tnNo: " << tnNo;
  #endif

  for (int iM = 0; iM < nMsh; iM++) {
    auto& msh = com_mod.msh[iM];
    msh.lN.resize(tnNo);
    #ifdef debug_distribute
    dmsg << "---------- iM " << iM;
    dmsg << "msh.gnNo: " << msh.gnNo;
    dmsg << "msh.nNo: " << msh.nNo;
    dmsg << "msh.nEl: " << msh.nEl;
    #endif

    for (int a = 0; a < msh.nNo; a++) {
      int Ac = msh.gN[a];
      msh.lN[Ac] = a;
    }

    for (int e = 0; e < msh.nEl; e++) {
      for (int a = 0; a < msh.eNoN; a++) { 
        int Ac = msh.IEN(a,e);
        msh.IEN(a,e) = msh.gN[Ac];
      }
    }
  }

  // Rearrange body force structure, if necessary
  //
  if (cm.seq()) {
    for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
      auto& eq = com_mod.eq[iEq];

      for (int iBf = 0; iBf < eq.nBf; iBf++) {
        auto& bf = eq.bf[iBf];

        if (bf.bx.size() != 0) {
          int i = bf.dof;
          int iM = bf.iM;
          auto& msh = com_mod.msh[iM];

          Array<double> tmpX(bf.bx);
          bf.bx.resize(i, msh.nNo);

          for (int a = 0; a < msh.nNo; a++) {
            int Ac = msh.gN[a];
            bf.bx.set_col(a, tmpX.col(Ac));
          } 

        } else if (bf.bm.defined()) {
          int i = bf.bm.dof;
          int a = bf.bm.nTP;
          int iM = bf.iM;
          auto tmpD = bf.bm.d;

          bf.bm.d.resize(i,com_mod.msh[iM].nNo,a);

          for (int i = 0; i < bf.bm.nTP; i++) {
            for (int a = 0; a < com_mod.msh[iM].nNo; a++) {
              int Ac = com_mod.msh[iM].gN[a];
              for (int j = 0; j < bf.bm.dof; j++) {
                bf.bm.d(j,a,i) = tmpD(j,Ac,i);
              }            
            }            
          }            
        }            
      }            
    }            
    return; 
  }            

  // Partitioning the faces
  //
  // tMs is a temporary variable to keep fa%gN of the old meshes.
  //
  #ifdef debug_distribute
  dmsg << " " << " ";
  dmsg << "Partitioning the faces" << " ...";
  #endif
  std::vector<mshType> tMs(nMsh);

  for (int iM = 0; iM < nMsh; iM++) {
    auto& msh = com_mod.msh[iM];
    tMs[iM].fa.resize(msh.nFa);

    for (int iFa = 0; iFa < msh.nFa; iFa++) {
      auto& face = msh.fa[iFa];
      part_face(simulation, msh, face, tMs[iM].fa[iFa], gmtl);
    }
  }

  #ifdef debug_distribute
  dmsg << " " << " ";
  dmsg << "Sending data read by master to slaves " << " ...";
  #endif


  if (!com_mod.resetSim) {
    cm.bcast(cm_mod, &com_mod.nsymd);
    cm.bcast(cm_mod, com_mod.stopTrigName);
    cm.bcast(cm_mod, com_mod.iniFilePath);
    cm.bcast(cm_mod, com_mod.stFileName);
    cm.bcast(cm_mod, &com_mod.stFileFlag);
    cm.bcast(cm_mod, &com_mod.stFileIncr);
    cm.bcast(cm_mod, &com_mod.stFileRepl);
    cm.bcast(cm_mod, &com_mod.saveIncr);

    cm.bcast(cm_mod, &com_mod.saveATS);
    cm.bcast(cm_mod, &com_mod.saveAve);
    cm.bcast(cm_mod, &com_mod.saveVTK);
    cm.bcast(cm_mod, &com_mod.bin2VTK);

    cm.bcast(cm_mod, &com_mod.mvMsh);

    cm.bcast(cm_mod, &com_mod.nITs);
    cm.bcast(cm_mod, &com_mod.nTS);
    cm.bcast(cm_mod, &com_mod.startTS);
    cm.bcast(cm_mod, &com_mod.nEq);
    cm.bcast(cm_mod, &com_mod.dt);
    cm.bcast(cm_mod, &com_mod.precompDt);

    cm.bcast(cm_mod, &com_mod.zeroAve);
    cm.bcast(cm_mod, &com_mod.cmmInit);
    cm.bcast(cm_mod, &com_mod.cmmVarWall);

    cm.bcast(cm_mod, &com_mod.shlEq);
    cm.bcast(cm_mod, &com_mod.pstEq);
    cm.bcast(cm_mod, &com_mod.sstEq);

    cm.bcast(cm_mod, &simulation->cep_mod.cepEq);
    cm.bcast(cm_mod, &com_mod.risFlag);
    cm.bcast(cm_mod, &com_mod.ris0DFlag);
    cm.bcast(cm_mod, &com_mod.urisFlag);
    cm.bcast(cm_mod, &com_mod.urisActFlag);
    cm.bcast(cm_mod, &com_mod.urisRes);
    cm.bcast(cm_mod, &com_mod.urisResClose);
    cm.bcast(cm_mod, &com_mod.usePrecomp);
    if (com_mod.rmsh.isReqd) {
      auto& rmsh = com_mod.rmsh;
      cm.bcast_enum(cm_mod, &rmsh.method);
      cm.bcast(cm_mod, &rmsh.freq);
      cm.bcast(cm_mod, &rmsh.cpVar);

      if (cm.slv(cm_mod)) {
        rmsh.maxEdgeSize.resize(com_mod.nMsh);
        rmsh.minDihedAng = 0.0;
        rmsh.maxRadRatio = 0.0;
      }

      cm.bcast(cm_mod, rmsh.maxEdgeSize);
    }

    cm.bcast(cm_mod, &com_mod.iCntct);

    if (com_mod.iCntct) {
      cm.bcast_enum(cm_mod, &com_mod.cntctM.cType);
      cm.bcast(cm_mod, &com_mod.cntctM.k);
      cm.bcast(cm_mod, &com_mod.cntctM.c);
      cm.bcast(cm_mod, &com_mod.cntctM.h);
      cm.bcast(cm_mod, &com_mod.cntctM.al);
    }

    cm.bcast(cm_mod, &com_mod.ibFlag);
    cm.bcast(cm_mod, &simulation->cep_mod.nXion);
  }  

  // Distributing X to processors
  //
  #ifdef debug_distribute
  dmsg << " " << " ";
  dmsg << "Distributing X to processors " << " ...";
  #endif
  Array<double> tmpX;

  if (cm.mas(cm_mod)) {
    tmpX.resize(com_mod.nsd, com_mod.gtnNo);
    tmpX = com_mod.x;
    com_mod.x.clear();
  } 

  com_mod.x.resize(com_mod.nsd, com_mod.tnNo);
  com_mod.x = all_fun::local(com_mod, cm_mod, cm, tmpX);

  // Distributing lM%dmnId if present to processors
  //
  bool flag = (com_mod.dmnId.size() != 0);
  cm.bcast(cm_mod, &flag);
  Vector<int> part;
  #ifdef debug_distribute
  dmsg << "dmnId.size(): " << com_mod.dmnId.size();
  dmsg << "flag: " << flag;
  #endif

  if (flag) {
  #ifdef debug_distribute
    dmsg << "Distributing dmnId " << " ... ";
  #endif
    if (cm.mas(cm_mod)) {
      part.resize(com_mod.gtnNo);
      part = com_mod.dmnId;
      com_mod.dmnId.clear();
    } else { 
      part.clear();
    }
    com_mod.dmnId.resize(com_mod.tnNo);
    com_mod.dmnId = all_fun::local(com_mod, cm_mod, cm, part);
    part.clear();
  }

  // Distribute prestress (pS0) to processors.
  //
  // pS0 is set in read_msh() and for CMM.
  //
  flag = (com_mod.pS0.size() != 0);
  cm.bcast(cm_mod, &flag);

  if (flag) {
    if (cm.mas(cm_mod)) {
      tmpX.resize(com_mod.nsymd, com_mod.gtnNo);
      tmpX = com_mod.pS0;
      com_mod.pS0.clear();
    } else { 
      tmpX.clear();
    }
    com_mod.pS0.resize(com_mod.nsymd, com_mod.tnNo);
    com_mod.pS0 = all_fun::local(com_mod, cm_mod, cm, tmpX);
    tmpX.clear();
  }

  //  Distribute initial flow quantities to processors
  //
  // Distribute Pinit (set in read_msh()).
  //
  flag = (com_mod.Pinit.size() != 0);
  cm.bcast(cm_mod, &flag);

  if (flag) {
    if (cm.mas(cm_mod)) {
      tmpX.resize(1, com_mod.gtnNo);
      for (int i = 0; i < tmpX.ncols(); i++) {
        tmpX(0,i) = com_mod.Pinit[i];
      }
      com_mod.Pinit.clear();
    } else { 
      tmpX.clear();
    }
    wgt.resize(1, com_mod.tnNo); 
    com_mod.Pinit.resize(com_mod.tnNo);
    wgt = all_fun::local(com_mod, cm_mod, cm, tmpX);
    for (int i = 0; i < wgt.ncols(); i++) {
      com_mod.Pinit[i] = wgt(0,i);
    }
    tmpX.clear(); 
    wgt.clear();
  }

  // Distribute Vinit
  //
  flag = (com_mod.Vinit.size() != 0);
  cm.bcast(cm_mod, &flag);

  if (flag) {
    if (cm.mas(cm_mod)) {
      tmpX.resize(com_mod.nsd, com_mod.gtnNo);
      tmpX = com_mod.Vinit;
      com_mod.Vinit.clear();
    } else {
      tmpX.clear();
    }
    com_mod.Vinit.resize(com_mod.nsd, com_mod.tnNo);
    com_mod.Vinit = all_fun::local(com_mod, cm_mod, cm, tmpX);
    tmpX.clear();
  }

  // Distribute Dinit
  //
  flag = (com_mod.Dinit.size() != 0);
  cm.bcast(cm_mod, &flag);

  if (flag) {
    if (cm.mas(cm_mod)) {
      tmpX.resize(com_mod.nsd, com_mod.gtnNo);
      tmpX = com_mod.Dinit;
      com_mod.Dinit.clear();
    } else {
      tmpX.clear();
    }
    com_mod.Dinit.resize(com_mod.nsd, com_mod.tnNo);
    com_mod.Dinit = all_fun::local(com_mod, cm_mod, cm, tmpX);
    tmpX.clear();
  }

  // And distributing eq to processors
  //
  if (cm.slv(cm_mod)) {
    com_mod.eq.resize(com_mod.nEq);
  }

  auto& cep_mod = simulation->cep_mod;
  for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
    auto& eq = com_mod.eq[iEq];
    dist_eq(com_mod, cm_mod, cm, tMs, gmtl, cep_mod, eq);
  }

  // For CMM initialization
  //
  flag = (com_mod.cmmBdry.size() != 0);
  cm.bcast(cm_mod, &flag);
  if (flag) {
    Vector<int> part;
    if (cm.mas(cm_mod)) {
      part.resize(com_mod.gtnNo);
      part = com_mod.cmmBdry;
    }
    com_mod.cmmBdry.resize(com_mod.tnNo);
    com_mod.cmmBdry = all_fun::local(com_mod, cm_mod, cm, part);
  }

  // For CMM variable wall properties
  //
  flag = (com_mod.varWallProps.size() != 0);
  cm.bcast(cm_mod, &flag);

  if (flag) {
    Array<double> tmpX;
    if (cm.mas(cm_mod)) {
      tmpX.resize(2, com_mod.gtnNo);
      tmpX = com_mod.varWallProps;
      com_mod.varWallProps.clear();
    }
    com_mod.varWallProps.resize(2, com_mod.tnNo);
    com_mod.varWallProps = all_fun::local(com_mod, cm_mod, cm, tmpX);
  }

  // Communicating cplBC data
  //
  auto& cplBC = com_mod.cplBC;
  cm.bcast(cm_mod, &cplBC.nFa);
  cm.bcast_enum(cm_mod, &cplBC.schm);
  cm.bcast(cm_mod, &cplBC.useGenBC);
  cm.bcast(cm_mod, &cplBC.useSvZeroD);

  if (cplBC.useGenBC) {   
    if (cm.slv(cm_mod)) {   
      cplBC.nX = 0;
      cplBC.xo.resize(cplBC.nX);
    }

  } else if (cplBC.useSvZeroD) {   
    if (cm.slv(cm_mod)) {   
      cplBC.nX = 0;
      cplBC.xo.resize(cplBC.nX);
    }

  } else { 
    cm.bcast(cm_mod, &cplBC.nX);
    if (cplBC.xo.size() == 0) {
       cplBC.xo.resize(cplBC.nX);
    }
    if (cplBC.nX != 0) {
      cm.bcast(cm_mod, cplBC.xo);
    }
  }

  cm.bcast(cm_mod, &cplBC.initRCR);

  if (com_mod.risFlag) {
    for (int iProj = 0; iProj < com_mod.ris.nbrRIS; iProj++) {
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < com_mod.risMapList[iProj].map.ncols(); j++) {
          int tmp = gmtl(com_mod.grisMapList[iProj].map(i,j));
          if (tmp != -1) {
            com_mod.risMapList[iProj].map(i,j) = com_mod.msh[i].lN(tmp);
          } else {
            com_mod.risMapList[iProj].map(i,j) = -1;
          }
          com_mod.grisMapList[iProj].map(i,j) = tmp;
        }
      }
    }
  }
}


void dist_bc(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, bcType& lBc, const std::vector<mshType>& tMs,
             const Vector<int>& gmtl)
{
  using namespace consts;

  #define n_debug_dist_bc
  #ifdef debug_dist_bc
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  int task_id = cm.idcm();
  bool is_slave = cm.slv(cm_mod);
  cm.bcast(cm_mod, &lBc.cplBCptr);
  cm.bcast(cm_mod, &lBc.bType);
  cm.bcast(cm_mod, &lBc.clsFlgRis);

  int nsd = com_mod.nsd;
  #ifdef debug_dist_bc
  dmsg << "nsd: " << nsd;
  dmsg << "lBc.bType: " << lBc.bType;
  dmsg << "is_slave: " << is_slave;
  #endif

  if (is_slave) {
    lBc.eDrn.resize(nsd); 
    lBc.h.resize(nsd);
  }

  cm.bcast(cm_mod, lBc.eDrn);

  cm.bcast(cm_mod, &lBc.iFa);
  cm.bcast(cm_mod, &lBc.iM);
  cm.bcast(cm_mod, &lBc.r);
  cm.bcast(cm_mod, &lBc.g);
  cm.bcast(cm_mod, lBc.h);
  cm.bcast(cm_mod, &lBc.weakDir);
  cm.bcast(cm_mod, lBc.tauB);
  cm.bcast(cm_mod, &lBc.flwP);


  if (utils::btest(lBc.bType, static_cast<int>(BoundaryConditionType::bType_RCR))) {
    cm.bcast(cm_mod, &lBc.RCR.Rp);
    cm.bcast(cm_mod, &lBc.RCR.C);
    cm.bcast(cm_mod, &lBc.RCR.Rd);
    cm.bcast(cm_mod, &lBc.RCR.Pd);
    cm.bcast(cm_mod, &lBc.RCR.Xo);
  }

  // Communicating time-dependent BC data
  //
  // lBc.gt is declare ALLOCATABLE in MOD.f but we don't
  // want to use pointers so use the define() method
  // to check if it has data define.
  //
  bool flag = lBc.gt.defined(); 
  cm.bcast(cm_mod, &flag);

  if (flag) {
    if (is_slave) {
      // [NOTE] This is allocated in ComMod.
      //lBc.gt = new fcType;
    }

    cm.bcast(cm_mod, &lBc.gt.lrmp);
    cm.bcast(cm_mod, &lBc.gt.d);
    cm.bcast(cm_mod, &lBc.gt.n);

    int j = lBc.gt.d;
    int i = lBc.gt.n;

    if (is_slave) { 
      lBc.gt.qi.resize(j);
      lBc.gt.qs.resize(j);
      lBc.gt.r.resize(j,i);
      lBc.gt.i.resize(j,i);
    }

    cm.bcast(cm_mod, &lBc.gt.ti);
    cm.bcast(cm_mod, &lBc.gt.T);
    cm.bcast(cm_mod, lBc.gt.qi);
    cm.bcast(cm_mod, lBc.gt.qs);
    cm.bcast(cm_mod, lBc.gt.r);
    cm.bcast(cm_mod, lBc.gt.i);
  }

  // Communicating moving BC data
  flag = lBc.gm.defined();
  cm.bcast(cm_mod, &flag);
  
  if (flag) {
    if (is_slave) {
      //lBc.gm = new MBType;
    }

    cm.bcast(cm_mod, &lBc.gm.period);

    // Communication the .t data
    cm.bcast(cm_mod, &lBc.gm.nTP);
    cm.bcast(cm_mod, &lBc.gm.dof);
    int nTp = lBc.gm.nTP;
    int iDof = lBc.gm.dof;

    if (is_slave) {
     lBc.gm.t.resize(nTp);
    }

    cm.bcast(cm_mod, lBc.gm.t);
    int nNo = tMs[lBc.iM].fa[lBc.iFa].nNo;
    int a = nTp*iDof*nNo;

    // Allocating the container and copying the nodes which belong to
    // this processor
    //
    Vector<double> tmp(a);

    if (!is_slave) {
      // Copy data row-wise to tmp.
      int n = 0;
      for (int k = 0; k < lBc.gm.d.nslices(); k++) {
        for (int j = 0; j < lBc.gm.d.ncols(); j++) {
          for (int i = 0; i < lBc.gm.d.nrows(); i++) {
            tmp[n] = lBc.gm.d(i,j,k);
            n += 1;
          }
        }
      }
    }

    cm.bcast(cm_mod, tmp);

    // This is the new number of nodes
    a = com_mod.msh[lBc.iM].fa[lBc.iFa].nNo;
    lBc.gm.d.resize(iDof, a, nTp);
    int b = 0;

    for (int a = 0; a < nNo; a++) {
      int Ac = tMs[lBc.iM].fa[lBc.iFa].gN[a];
      Ac = gmtl[Ac];
      if (Ac != -1) {
        for (int i = 0; i < nTp; i++) {
          int j = iDof * (i*nNo + a);
          for (int k = 0; k < iDof; k++) {
            lBc.gm.d(k,b,i) = tmp[k+j];
          }
        } 
        b = b + 1;
      }
    }
  }

  // Communicating profile data
  //
  flag = (lBc.gx.size() != 0);
  cm.bcast(cm_mod, &flag);

  if (flag) {
    int nNo = tMs[lBc.iM].fa[lBc.iFa].nNo;
    Vector<double> tmp(nNo);
    if (!is_slave) {
      tmp = lBc.gx;
      lBc.gx.clear();
    }

    cm.bcast(cm_mod, tmp);

    // This is the new number of nodes
    int a = com_mod.msh[lBc.iM].fa[lBc.iFa].nNo;
    lBc.gx.resize(a);
    int b = 0;
    for (int a = 0; a < nNo; a++) {
      int Ac = tMs[lBc.iM].fa[lBc.iFa].gN[a];
      Ac = gmtl[Ac];
      if (Ac != -1) {
        lBc.gx[b] = tmp[a];
        b = b + 1;
      }
    }
  }

  // Communicating Robin BC
  //
  bool has_robin_bc = lBc.robin_bc.is_initialized();
  cm.bcast(cm_mod, &has_robin_bc); // Master process broadcasts the flag to all processes

  if (has_robin_bc) {
    lBc.robin_bc.distribute(com_mod, cm_mod, cm, com_mod.msh[lBc.iM].fa[lBc.iFa]);
  }


  // Communicating and reordering master node data for 
  // undeforming Neumann BC faces.
  //
  if (utils::btest(lBc.bType, static_cast<int>(BoundaryConditionType::bType_undefNeu))) {
    cm.bcast(cm_mod, &lBc.masN);
    int iM = lBc.iM;
    int iFa = lBc.iFa;
    int nNo = com_mod.msh[iM].fa[iFa].nNo;
    bool flag = (nNo != 0);

    // Action performed only on the processes that share the face
    //
    if (flag) {
      int Ac = lBc.masN;
      int a = gmtl[Ac];

      // The process that owns the node is set to be master. For other
      // processes that share the face but do not own the node, we add
      // this node as a ghost master node. ltg pointer is reset.
      //
      if (a != 0) {
        lBc.masN = a;
      } else { 
        nNo = com_mod.tnNo;
        com_mod.tnNo = com_mod.tnNo + 1;

        // Remap ltg
        //
        Vector<int> tmpI(nNo);
        tmpI = com_mod.ltg;
        com_mod.ltg.resize(com_mod.tnNo);
        for (int i = 0;  i < nNo; i++) {
          com_mod.ltg[i] = tmpI[i];
        }
        com_mod.ltg[com_mod.tnNo] = Ac;   
        lBc.masN = com_mod.tnNo;

        // Add the ghost master node to the face data structure
        //
        nNo = com_mod.msh[iM].fa[iFa].nNo;
        com_mod.msh[iM].fa[iFa].nNo = nNo + 1;

        tmpI.resize(nNo);
        for (int i = 0;  i < nNo; i++) {
          tmpI[i] = com_mod.msh[iM].fa[iFa].gN[i];
        }
        com_mod.msh[iM].fa[iFa].gN.resize(com_mod.msh[iM].fa[iFa].nNo);

        for (int i = 0;  i < nNo; i++) {
          com_mod.msh[iM].fa[iFa].gN[i] = tmpI[i];
        }

        com_mod.msh[iM].fa[iFa].gN[nNo] = com_mod.tnNo;
      }

     // Zero out master node if not part of the face
     } else { 
       lBc.masN = 0;
    }
  } else { 
     lBc.masN = 0;
  }
}


//---------
// dist_bf
//---------
//
void dist_bf(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, bfType& lBf)
{
  using namespace consts;

  #define n_dist_bf
  #ifdef dist_bf
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "cm.seq(): " << cm.seq();
  #endif

  bool is_slave = cm.slv(cm_mod);

  cm.bcast(cm_mod, &lBf.bType);
  cm.bcast(cm_mod, &lBf.dof);
  cm.bcast(cm_mod, &lBf.iM);
  int iM = lBf.iM;

  // Steady value
  //
  bool flag = (lBf.b.size() != 0);
  cm.bcast(cm_mod, &flag);
  if (flag) {
    if (is_slave) { 
      lBf.b.resize(lBf.dof);
    }
    cm.bcast(cm_mod, lBf.b);
  }

  // Communicating spatially dependent BF
  //
  flag = (lBf.bx.size() != 0);
  cm.bcast(cm_mod, &flag);

  Array<double> tmpX;
  Array3<double> tmpD;

  if (flag) {
    if (!is_slave) {
      tmpX.resize(lBf.dof, com_mod.gtnNo);
      tmpX = lBf.bx;
      lBf.bx.clear();
    }

    lBf.bx.resize(lBf.dof, com_mod.tnNo);
    lBf.bx = all_fun::local(com_mod, cm_mod, cm, tmpX);
    tmpX.resize(lBf.dof, com_mod.tnNo);
    tmpX = lBf.bx;
    lBf.bx.clear();
    lBf.bx.resize(lBf.dof, com_mod.msh[iM].nNo);

    for (int a = 0; a < com_mod.msh[iM].nNo; a++) {
      int Ac = com_mod.msh[iM].gN[a];
      lBf.bx.set_col(a, tmpX.col(Ac));
    }
  }

  //  Communicating time-dependent BF data
  flag = lBf.bt.defined();
  cm.bcast(cm_mod, &flag);
  if (flag) {
    if (is_slave) {
      //ALLOCATE(lBf.bt)
    }
    cm.bcast(cm_mod, &lBf.bt.lrmp);
    cm.bcast(cm_mod, &lBf.bt.d);
    cm.bcast(cm_mod, &lBf.bt.n);

    if (is_slave) {
      int j = lBf.bt.d;
      int i = lBf.bt.n;
      lBf.bt.qi.resize(j);
      lBf.bt.qs.resize(j);
      lBf.bt.r.resize(j,i);
      lBf.bt.i.resize(j,i);
    }
    cm.bcast(cm_mod, &lBf.bt.ti);
    cm.bcast(cm_mod, &lBf.bt.T);
    cm.bcast(cm_mod, lBf.bt.qi);
    cm.bcast(cm_mod, lBf.bt.qs);
    cm.bcast(cm_mod, lBf.bt.r);
    cm.bcast(cm_mod, lBf.bt.i);
  }

  // Communicating moving BF data
  //
  flag = lBf.bm.defined();
  cm.bcast(cm_mod, &flag);

  if (flag) {
    if (is_slave) {
      //ALLOCATE(lBf.bm)
    }
    cm.bcast(cm_mod, &lBf.bm.period);
    cm.bcast(cm_mod, &lBf.bm.nTP);
    cm.bcast(cm_mod, &lBf.bm.dof);

    int nTP  = lBf.bm.nTP;
    int idof = lBf.bm.dof;

    if (is_slave) {
      lBf.bm.t.resize(nTP);
    }
    cm.bcast(cm_mod, lBf.bm.t);

    if (!is_slave) {
      tmpD.resize(lBf.bm.dof, com_mod.gtnNo, lBf.bm.nTP);
      tmpD = lBf.bm.d;
      lBf.bm.d.clear();
    }

    tmpX.resize(lBf.bm.dof, com_mod.tnNo);
    lBf.bm.d.resize(lBf.bm.dof, com_mod.msh[iM].nNo, lBf.bm.nTP);

    for (int i = 0; i < lBf.bm.nTP; i++) {
      // [NOTE] This is a big hack because tmpD can have 0 size
      // for a slave process.
      Array<double> tmpX_slice(lBf.bm.dof, com_mod.msh[iM].nNo);
      if (tmpD.size() != 0) {
        tmpX_slice = tmpD.slice(i);
      }
      tmpX = all_fun::local(com_mod, cm_mod, cm, tmpX_slice);

      for (int a = 0; a < com_mod.msh[iM].nNo; a++) {
        int Ac = com_mod.msh[iM].gN[a];
        for (int j = 0; j < lBf.bm.d.nrows(); j++) {
          lBf.bm.d(j,a,i) = tmpX(j,Ac);
        }
      }
    }
  }
}

//---------
// dist_ris
//---------
/// @brief This rountine distributes data structures for RIS
void dist_ris(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm)
{
  using namespace consts;

  #define n_debug_dist_ris
  #ifdef debug_dist_ris
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif
  int task_id = cm.idcm();
  auto& RIS = com_mod.ris;

  #ifdef debug_dist_ris
  dmsg << "risFlag: " << com_mod.risFlag;
  #endif

  Vector<int> dims(3);
  int lst_size = 0;
  int maplist_size = 0;
  Vector<int> nStks;
  int nProj = 0;

  // First find the RIS nodes between mesh pairs on the current processor
  if (cm.mas(cm_mod)) {
    nProj = RIS.nbrRIS;
    dims(0) = RIS.lst.nrows();
    dims(1) = RIS.lst.ncols();
    dims(2) = RIS.lst.nslices();
    lst_size = RIS.lst.size();
    nStks.resize(nProj);
    for (int iProj = 0; iProj < nProj; iProj++) {
      nStks(iProj) = com_mod.risMapList[iProj].map.ncols();
      maplist_size += 2*nStks(iProj);
    }
  }

  cm.bcast(cm_mod, &nProj);
  cm.bcast(cm_mod, dims);
  cm.bcast(cm_mod, &lst_size);
  if (cm.slv(cm_mod)) {nStks.resize(nProj);}
  cm.bcast(cm_mod, nStks);
  cm.bcast(cm_mod, &maplist_size);

  if (cm.slv(cm_mod)) {
    com_mod.risMapList.resize(nProj);
    com_mod.grisMapList.resize(nProj);
    for (int iProj = 0; iProj < nProj; iProj++) {
      com_mod.risMapList[iProj].map.resize(2, nStks(iProj));
      com_mod.grisMapList[iProj].map.resize(2, nStks(iProj));
    }
    RIS.nbrIter.resize(nProj);
    RIS.Res.resize(nProj);
    RIS.clsFlg.resize(nProj);
    RIS.meanP.resize(nProj,2);
    RIS.meanFl.resize(nProj);
    RIS.status.resize(nProj);
    RIS.lst.resize(dims(0),dims(1),dims(2));
  }
  dims.clear();

  Vector<int> flat_lst;
  Vector<int> flat_maplist;
  Vector<int> flat_gmaplist;

  if (cm.mas(cm_mod)) {
    flat_lst.resize(lst_size);
    int n = 0;
    for (int i = 0; i < RIS.lst.nrows(); i++) {
      for (int j = 0; j < RIS.lst.ncols(); j++) {
        for (int k = 0; k < RIS.lst.nslices(); k++) {
          flat_lst(n) = RIS.lst(i,j,k);
          n += 1;
        }
      }
    }
    n = 0;
    flat_maplist.resize(maplist_size);
    flat_gmaplist.resize(maplist_size);
    for (int iProj = 0; iProj < nProj; iProj++) {
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < nStks(iProj); j++) {
          flat_maplist(n) = com_mod.risMapList[iProj].map(i,j);
          flat_gmaplist(n) = com_mod.grisMapList[iProj].map(i,j);
          n += 1;
        }
      }
    }
  } else {
    flat_lst.resize(lst_size);
    flat_maplist.resize(maplist_size);
    flat_gmaplist.resize(maplist_size);
  }

  cm.bcast(cm_mod, flat_lst);
  cm.bcast(cm_mod, flat_maplist);
  cm.bcast(cm_mod, flat_gmaplist);
  cm.bcast(cm_mod, RIS.clsFlg);
  cm.bcast(cm_mod, RIS.nbrIter);
  cm.bcast(cm_mod, RIS.Res);
  cm.bcast(cm_mod, &RIS.nbrRIS);
  cm.bcast(cm_mod, RIS.meanP);
  cm.bcast(cm_mod, RIS.meanFl);
  cm.bcast(cm_mod, RIS.status);

  // Reshape the 1D array back to 3D after broadcasting
  if (cm.slv(cm_mod)) {
    int n = 0;
    for (int i = 0; i < RIS.lst.nrows(); i++) {
      for (int j = 0; j < RIS.lst.ncols(); j++) {
        for (int k = 0; k < RIS.lst.nslices(); k++) {
          RIS.lst(i,j,k) = flat_lst(n);
          n += 1;
        }
      }
    }
    n = 0;
    for (int iProj = 0; iProj < nProj; iProj++) {
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < nStks(iProj); j++) {
          com_mod.risMapList[iProj].map(i,j) = flat_maplist(n);
          com_mod.grisMapList[iProj].map(i,j) = flat_gmaplist(n);
          n += 1;
        }
      }
    }
  }
  flat_lst.clear();
  flat_maplist.clear();
  flat_gmaplist.clear();

  // // [HZ] broadcasting nestd vector causes a problem in parallel
  // for (int iProj = 0; iProj < nProj; iProj++) {
  //   cm.bcast(cm_mod, com_mod.risMapList[iProj].map);
  //   cm.bcast(cm_mod, com_mod.grisMapList[iProj].map);
  // }

  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    cm.bcast(cm_mod, &com_mod.msh[iM].gnEl);
    if (cm.slv(cm_mod)) {
      com_mod.msh[iM].eRIS.resize(com_mod.msh[iM].gnEl);
      com_mod.msh[iM].partRIS.resize(com_mod.msh[iM].gnEl);
    }
    cm.bcast(cm_mod, com_mod.msh[iM].eRIS);
    cm.bcast(cm_mod, com_mod.msh[iM].partRIS);
  }
  nStks.clear();
}


//---------
// dist_uris
//---------
/// @brief This routine distributes data structures for URIS
void dist_uris(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm) {
  using namespace consts;

  #define n_debug_dist_uris
  #ifdef debug_dist_uris
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  auto& uris = com_mod.uris;

  // We will copy the URIS valves to all processors 
  cm.bcast(cm_mod, &com_mod.nUris);
  if (cm.slv(cm_mod)) {
    uris.resize(com_mod.nUris);
    for (int iUris = 0; iUris < com_mod.nUris; iUris++) {
      uris[iUris].nrm.resize(com_mod.nsd);
    }
  }

  for (int iUris = 0; iUris < com_mod.nUris; iUris++) {
    cm.bcast(cm_mod, uris[iUris].name);
    cm.bcast(cm_mod, &uris[iUris].tnNo);
    cm.bcast(cm_mod, &uris[iUris].nFa);
    cm.bcast(cm_mod, &uris[iUris].sdf_default);
    cm.bcast(cm_mod, &uris[iUris].sdf_deps);
    cm.bcast(cm_mod, &uris[iUris].sdf_deps_close);
    cm.bcast(cm_mod, &uris[iUris].clsFlg);
    cm.bcast(cm_mod, &uris[iUris].cnt);
    cm.bcast(cm_mod, &uris[iUris].scF);
    cm.bcast(cm_mod, uris[iUris].nrm);
  }

  std::vector<Vector<int>> lM_gN_flat(com_mod.nUris);
  std::vector<Vector<int>> lM_lN_flat(com_mod.nUris);
  std::vector<Vector<int>> lM_IEN_flat(com_mod.nUris);
  Vector<int> lM_gN_flat_size(com_mod.nUris);
  lM_gN_flat_size = 0;
  Vector<int> lM_lN_flat_size(com_mod.nUris);
  lM_lN_flat_size = 0;
  Vector<int> lM_IEN_flat_size(com_mod.nUris);
  lM_IEN_flat_size = 0;
  int n0, n1, n2;

  if (cm.mas(cm_mod)) {
    for (int iUris = 0; iUris < com_mod.nUris; iUris++) {
      for (int iM = 0; iM < uris[iUris].nFa; iM++) {
        auto& lM = uris[iUris].msh[iM];
        lM_gN_flat_size(iUris) += lM.gnNo;
        lM_lN_flat_size(iUris) += uris[iUris].tnNo;
        lM_IEN_flat_size(iUris) += lM.eNoN*lM.gnEl;
      }
      // std::cout << "lM_gN_flat_size: " << lM_gN_flat_size(iUris) << " for uris: " << iUris << std::endl;
      // std::cout << "lM_lN_flat_size: " << lM_lN_flat_size(iUris) << " for uris: " << iUris << std::endl;
      // std::cout << "lM_IEN_flat_size: " << lM_IEN_flat_size(iUris) << " for uris: " << iUris << std::endl;
    }
  }

  cm.bcast(cm_mod, lM_gN_flat_size);
  cm.bcast(cm_mod, lM_lN_flat_size);
  cm.bcast(cm_mod, lM_IEN_flat_size);

  for (int iUris = 0; iUris < com_mod.nUris; iUris++) {
    lM_gN_flat[iUris].resize(lM_gN_flat_size(iUris));
    lM_lN_flat[iUris].resize(lM_lN_flat_size(iUris));
    lM_IEN_flat[iUris].resize(lM_IEN_flat_size(iUris));
  }

  if (cm.mas(cm_mod)) {
    for (int iUris = 0; iUris < com_mod.nUris; iUris++) {
      n0 = 0;
      n1 = 0;
      n2 = 0;
      for (int iM = 0; iM < uris[iUris].nFa; iM++) {
        auto& lM = uris[iUris].msh[iM];

        for (int i = 0; i < lM.gnNo; i++) {
          lM_gN_flat[iUris](n0) = lM.gN(i);
          n0 += 1;
        }

        for (int i = 0; i < uris[iUris].tnNo; i++) {
          lM_lN_flat[iUris](n1) = lM.lN(i);
          n1 += 1;
        }

        for (int i = 0; i < lM.eNoN; i++) {
          for (int j = 0; j < lM.gnEl; j++) {
            lM_IEN_flat[iUris](n2) = lM.IEN(i,j);
            n2 += 1;
          }
        }

      }
    }
  }

  for (int iUris = 0; iUris < com_mod.nUris; iUris++) {
    cm.bcast(cm_mod, lM_gN_flat[iUris]);
    cm.bcast(cm_mod, lM_lN_flat[iUris]);
    cm.bcast(cm_mod, lM_IEN_flat[iUris]);
  }

  for (int iUris = 0; iUris < com_mod.nUris; iUris++) {
    if (cm.slv(cm_mod)) {
      uris[iUris].msh.resize(uris[iUris].nFa);
    }
    for (int iM = 0; iM < uris[iUris].nFa; iM++) {
      dist_uris_msh(com_mod, cm_mod, cm, uris[iUris].msh[iM], iUris);
    }
  }

  if (cm.slv(cm_mod)) {
    for (int iUris = 0; iUris < com_mod.nUris; iUris++) {
      n0 = 0;
      n1 = 0;
      n2 = 0;
      for (int iM = 0; iM < uris[iUris].nFa; iM++) {
        auto& lM = uris[iUris].msh[iM];

        for (int i = 0; i < lM.gnNo; i++) {
          lM.gN(i) = lM_gN_flat[iUris](n0);
          n0 += 1;
        }

        for (int i = 0; i < uris[iUris].tnNo; i++) {
          lM.lN(i) = lM_lN_flat[iUris](n1);
          n1 += 1;
        }

        for (int i = 0; i < lM.eNoN; i++) {
          for (int j = 0; j < lM.gnEl; j++) {
            lM.IEN(i,j) = lM_IEN_flat[iUris](n2);
            n2 += 1;
          }
        }
      }
    }
  }

  Vector<int> OpenN(com_mod.nUris);
  Vector<int> CloseN(com_mod.nUris);
  std::vector<Vector<double>> DxOpenFlat(com_mod.nUris);
  std::vector<Vector<double>> DxCloseFlat(com_mod.nUris);
  Vector <int> DxOpenFlatSize(com_mod.nUris);
  Vector <int> DxCloseFlatSize(com_mod.nUris);
  int n;
  if (cm.mas(cm_mod)) {
    for (int iUris = 0; iUris < com_mod.nUris; iUris++) {
      OpenN(iUris) = uris[iUris].DxOpen.nrows();
      CloseN(iUris) = uris[iUris].DxClose.nrows();
      DxOpenFlatSize(iUris) = uris[iUris].DxOpen.size();
      DxCloseFlatSize(iUris) = uris[iUris].DxClose.size();
      DxOpenFlat[iUris].resize(DxOpenFlatSize(iUris));
      DxCloseFlat[iUris].resize(DxCloseFlatSize(iUris));
      n = 0;
      for (int i = 0; i < uris[iUris].DxOpen.nrows(); i++) {
        for (int j = 0; j < uris[iUris].DxOpen.ncols(); j++) {
          for (int k = 0; k < uris[iUris].DxOpen.nslices(); k++) {
            DxOpenFlat[iUris](n) = uris[iUris].DxOpen(i,j,k);
            n += 1;
          }
        }
      }
      n = 0;
      for (int i = 0; i < uris[iUris].DxClose.nrows(); i++) {
        for (int j = 0; j < uris[iUris].DxClose.ncols(); j++) {
          for (int k = 0; k < uris[iUris].DxClose.nslices(); k++) {
            DxCloseFlat[iUris](n) = uris[iUris].DxClose(i,j,k);
            n += 1;
          }
        }
      }
    }
  }

  cm.bcast(cm_mod, OpenN);
  cm.bcast(cm_mod, CloseN);
  cm.bcast(cm_mod, DxOpenFlatSize);
  cm.bcast(cm_mod, DxCloseFlatSize);

  if (cm.slv(cm_mod)) {
    for (int iUris = 0; iUris < com_mod.nUris; iUris++) {
      uris[iUris].DxOpen.resize(OpenN(iUris), com_mod.nsd, uris[iUris].tnNo);
      uris[iUris].DxClose.resize(CloseN(iUris), com_mod.nsd, uris[iUris].tnNo);
      uris[iUris].x.resize(com_mod.nsd, uris[iUris].tnNo);
      uris[iUris].Yd.resize(com_mod.nsd, uris[iUris].tnNo);
      DxOpenFlat[iUris].resize(DxOpenFlatSize(iUris));
      DxCloseFlat[iUris].resize(DxCloseFlatSize(iUris));
    }
  }

  for (int iUris = 0; iUris < com_mod.nUris; iUris++) {
    cm.bcast(cm_mod, uris[iUris].x);
    cm.bcast(cm_mod, uris[iUris].Yd);
    cm.bcast(cm_mod, DxOpenFlat[iUris]);
    cm.bcast(cm_mod, DxCloseFlat[iUris]);
  }

  if (cm.slv(cm_mod)) {
    for (int iUris = 0; iUris < com_mod.nUris; iUris++) {
      n = 0;
      for (int i = 0; i < uris[iUris].DxOpen.nrows(); i++) {
        for (int j = 0; j < uris[iUris].DxOpen.ncols(); j++) {
          for (int k = 0; k < uris[iUris].DxOpen.nslices(); k++) {
            uris[iUris].DxOpen(i,j,k) = DxOpenFlat[iUris](n);
            n += 1;
          }
        }
      }
      n = 0;
      for (int i = 0; i < uris[iUris].DxClose.nrows(); i++) {
        for (int j = 0; j < uris[iUris].DxClose.ncols(); j++) {
          for (int k = 0; k < uris[iUris].DxClose.nslices(); k++) {
            uris[iUris].DxClose(i,j,k) = DxCloseFlat[iUris](n);
            n += 1;
          }
        }
      }
    }
  }

}


/// @brief This rountine distributes data structures for URIS mesh
void dist_uris_msh(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, mshType& lM, const int iUris)
{
  using namespace consts;

  #define n_debug_dist_uris_msh
  #ifdef debug_dist_uris_msh
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  cm.bcast(cm_mod, lM.name);
  cm.bcast(cm_mod, &lM.lShpF);
  cm.bcast(cm_mod, &lM.lShl);
  cm.bcast(cm_mod, &lM.lFib);
  cm.bcast_enum(cm_mod, &lM.eType);
  cm.bcast(cm_mod, &lM.eNoN);
  cm.bcast(cm_mod, &lM.gnEl);
  cm.bcast(cm_mod, &lM.gnNo);
  cm.bcast(cm_mod, &lM.nFa);
  cm.bcast(cm_mod, &lM.nFs);
  cm.bcast(cm_mod, &lM.nFn);
  cm.bcast(cm_mod, &lM.nG);
  cm.bcast(cm_mod, &lM.vtkType);
  cm.bcast(cm_mod, &lM.scF);
  cm.bcast(cm_mod, &lM.dx);
  // cm.bcast(cm_mod, &lM.fib.nFn);
  // cm.bcast(cm_mod, &lM.fib.locNd);
  // cm.bcast(cm_mod, &lM.fib.locEl);
  // cm.bcast(cm_mod, &lM.fib.locGP);

  if (cm.slv(cm_mod)) { 
    lM.nNo = lM.gnNo;
    lM.nEl = lM.gnEl;
    lM.gN.resize(lM.nNo);
    lM.lN.resize(com_mod.uris[iUris].tnNo);
    lM.IEN.resize(lM.eNoN, lM.nEl);
    lM.fa.resize(lM.nFa);
    nn::select_ele(com_mod, lM);
    // if (lM.fib.locNd) {lM.fib.fN.resize(nsd*lM.fib.nFn, 1, lM.nNo);}
    // if (lM.fib.locEl) {lM.fib.fN.resize(nsd*lM.fib.nFn, 1, lM.nEl);}
    // if (lM.fib.locGP) {lM.fib.fN.resize(nsd*lM.fib.nFn, lM.nG, lM.nEl);}
  }
  
  // // [HZ] Broadcasting here causes problem, similar reason to the RIS map list
  // cm.bcast(cm_mod, lM.gN);
  // cm.bcast(cm_mod, lM.lN);
  // cm.bcast(cm_mod, lM.IEN);

  // if (lM.fib.locNd || lM.fib.locEl || lM.fib.locGP) {
  //   cm.bcast(com_mod, &lM.fib.fN)
  // }

  if (lM.eType == ElementType::NRB) {
    cm.bcast(cm_mod, &lM.nSl);
    int insd = com_mod.nsd;
    if (lM.lShl) {insd = com_mod.nsd - 1;}
    if (cm.slv(cm_mod)) {
      lM.nW.resize(lM.gnNo);
      lM.INN.resize(insd, lM.gnEl);
      lM.bs.resize(insd);
    }
    cm.bcast(cm_mod, lM.nW);
    cm.bcast(cm_mod, lM.INN);
    for (int i = 0; i < insd; i++) {
      cm.bcast(cm_mod, &lM.bs[i].n);
      cm.bcast(cm_mod, &lM.bs[i].nG);
      cm.bcast(cm_mod, &lM.bs[i].nEl);
      cm.bcast(cm_mod, &lM.bs[i].nSl);
      cm.bcast(cm_mod, &lM.bs[i].p);
      if (cm.slv(cm_mod)) {lM.bs[i].xi.resize(lM.bs[i].n);}
      cm.bcast(cm_mod, lM.bs[i].xi);
      lM.bs[i].nNo = lM.bs[i].n - lM.bs[i].p - 1;
    }
  }

}

void dist_eq(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, const std::vector<mshType>& tMs,
             const Vector<int>& gmtl, CepMod& cep_mod, eqType& lEq)
{
  using namespace consts;

  #define ndist_eq
  #ifdef dist_eq
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  // Distribute equation parameters
  //
  cm.bcast(cm_mod, &lEq.nOutput);
  cm.bcast(cm_mod, &lEq.coupled);
  cm.bcast(cm_mod, &lEq.maxItr);
  cm.bcast(cm_mod, &lEq.minItr);
  cm.bcast(cm_mod, &lEq.roInf);
  cm.bcast_enum(cm_mod, &lEq.phys);
  cm.bcast(cm_mod, &lEq.nDmn);
  cm.bcast(cm_mod, &lEq.nBc);
  cm.bcast(cm_mod, &lEq.nBf);
  cm.bcast(cm_mod, &lEq.tol);
  cm.bcast(cm_mod, &lEq.useTLS);
  cm.bcast(cm_mod, &lEq.assmTLS);

  #ifdef dist_eq
  dmsg << "lEq.nOutput: " << lEq.nOutput;
  dmsg << "lEq.nDmn: " << lEq.nDmn;
  dmsg << "lEq.phys: " << lEq.phys;
  #endif

  if (com_mod.ibFlag) {
    cm.bcast(cm_mod, &lEq.nDmnIB);
    cm.bcast(cm_mod, &lEq.nBcIB);
  }

  // Distribute linear solver settings
  //
  cm.bcast(cm_mod, &lEq.FSILS.foC);
  cm.bcast_enum(cm_mod, &lEq.FSILS.LS_type);
  cm.bcast(cm_mod, &lEq.FSILS.RI.relTol);
  cm.bcast(cm_mod, &lEq.FSILS.GM.relTol);
  cm.bcast(cm_mod, &lEq.FSILS.CG.relTol);
  cm.bcast(cm_mod, &lEq.FSILS.RI.absTol);
  cm.bcast(cm_mod, &lEq.FSILS.GM.absTol);
  cm.bcast(cm_mod, &lEq.FSILS.CG.absTol);
  cm.bcast(cm_mod, &lEq.FSILS.RI.mItr);
  cm.bcast(cm_mod, &lEq.FSILS.GM.mItr);
  cm.bcast(cm_mod, &lEq.FSILS.CG.mItr);
  cm.bcast(cm_mod, &lEq.FSILS.RI.sD);
  cm.bcast(cm_mod, &lEq.FSILS.GM.sD);
  cm.bcast(cm_mod, &lEq.FSILS.CG.sD);

  cm.bcast_enum(cm_mod, &lEq.ls.LS_type);

  cm.bcast_enum(cm_mod, &lEq.linear_algebra_type);
  cm.bcast_enum(cm_mod, &lEq.linear_algebra_preconditioner);
  cm.bcast_enum(cm_mod, &lEq.linear_algebra_assembly_type);

  cm.bcast(cm_mod, &lEq.ls.relTol);
  cm.bcast(cm_mod, &lEq.ls.absTol);
  cm.bcast(cm_mod, &lEq.ls.mItr);
  cm.bcast(cm_mod, &lEq.ls.sD);

  #ifdef dist_eq
  dmsg << "lEq.phys: " << lEq.phys;
  dmsg << "lEq.ls.relTol: " << lEq.ls.relTol;
  dmsg << "lEq.ls.absTol: " << lEq.ls.absTol;
  dmsg << "lEq.ls.mItr: " << lEq.ls.mItr;
  dmsg << "lEq.ls.sD: " << lEq.ls.sD;
  #endif

  // Distribute domain properties
  if (cm.slv(cm_mod)) {
    lEq.dmn.resize(lEq.nDmn);
  }

  for (int iDmn = 0; iDmn < lEq.nDmn; iDmn++) {
    auto& dmn = lEq.dmn[iDmn];
    cm.bcast_enum(cm_mod, &dmn.phys);
    cm.bcast(cm_mod, &dmn.Id);
    cm.bcast_prop(cm_mod, dmn.prop);

    if (dmn.phys == EquationType::phys_CEP) {
      auto& cep = dmn.cep;
      cm.bcast_enum(cm_mod, &cep.cepType);
      cm.bcast(cm_mod, &cep.nX);
      cm.bcast(cm_mod, &cep.nG);
      cm.bcast(cm_mod, &cep.nFn);
      cm.bcast(cm_mod, &cep.imyo);
      cm.bcast(cm_mod, &cep.dt);
      cm.bcast(cm_mod, &cep.Ksac);
      cm.bcast(cm_mod, &cep.Diso);

      if (cm.slv(cm_mod)) {
        cep.Dani.resize(cep.nFn);
      } 

      cm.bcast(cm_mod, cep.Dani);
      cm.bcast(cm_mod, &cep.Istim.Ts);
      cm.bcast(cm_mod, &cep.Istim.Td);
      cm.bcast(cm_mod, &cep.Istim.CL);
      cm.bcast(cm_mod, &cep.Istim.A);
      cm.bcast_enum(cm_mod, &cep.odes.tIntType);

      if (cep.odes.tIntType == TimeIntegratioType::CN2) {
        cm.bcast(cm_mod, &cep.odes.maxItr);
        cm.bcast(cm_mod, &cep.odes.absTol);
        cm.bcast(cm_mod, &cep.odes.relTol);
      }

      cm.bcast(cm_mod, &cep_mod.ttp.G_Na);
      cm.bcast(cm_mod, &cep_mod.ttp.G_CaL);
      cm.bcast(cm_mod, &cep_mod.ttp.G_Kr);
      cm.bcast(cm_mod, cep_mod.ttp.G_Ks);
      cm.bcast(cm_mod, cep_mod.ttp.G_to);

      cm.bcast(cm_mod, cep_mod.bo.tau_si);
      cm.bcast(cm_mod, cep_mod.bo.tau_fi);
    } 

    if ((dmn.phys == EquationType::phys_struct) || (dmn.phys == EquationType::phys_ustruct)) {
      dist_mat_consts(com_mod, cm_mod, cm, dmn.stM);
      dist_solid_visc_model(com_mod, cm_mod, cm, dmn.solid_visc);
    } 

    if ((dmn.phys == EquationType::phys_fluid) || 
        (dmn.phys == EquationType::phys_stokes) || 
        (dmn.phys == EquationType::phys_CMM && !com_mod.cmmInit)) {
      dist_fluid_visc_model(com_mod, cm_mod, cm, dmn.fluid_visc);
    }
  }

  // Distribute cardiac electromechanics parameters
  //
  cm.bcast(cm_mod, &cep_mod.cem.cpld);

  if (cep_mod.cem.cpld); {
    cm.bcast(cm_mod, &cep_mod.cem.aStress);
    cm.bcast(cm_mod, &cep_mod.cem.aStrain);
  } 

  if (com_mod.ibFlag) {
    if (cm.slv(cm_mod)) {
      lEq.dmnIB.resize(lEq.nDmnIB);
    }
    for (int iDmn = 0; iDmn < lEq.nDmnIB; iDmn++) {
      auto& dmnIB = lEq.dmnIB[iDmn];
      cm.bcast_enum(cm_mod, &dmnIB.phys);
      cm.bcast(cm_mod, &dmnIB.Id);
      cm.bcast_prop(cm_mod, dmnIB.prop);
      dist_mat_consts(com_mod, cm_mod, cm, dmnIB.stM);
    }
  } 

  // Distribute ECG leads parameters
  //
  cm.bcast(cm_mod, &cep_mod.ecgleads.num_leads);
  #ifdef dist_eq
  dmsg << "cep_mod.ecgleads.num_leads: " << cep_mod.ecgleads.num_leads;
  #endif

  if (cep_mod.ecgleads.num_leads != 0) {
    if (!cm.mas(cm_mod)) {
      cep_mod.ecgleads.x_coords.resize(cep_mod.ecgleads.num_leads);
      cep_mod.ecgleads.y_coords.resize(cep_mod.ecgleads.num_leads);
      cep_mod.ecgleads.z_coords.resize(cep_mod.ecgleads.num_leads);
      cep_mod.ecgleads.pseudo_ECG.resize(cep_mod.ecgleads.num_leads);
    }

    cm.bcast(cm_mod, cep_mod.ecgleads.x_coords);
    cm.bcast(cm_mod, cep_mod.ecgleads.y_coords);
    cm.bcast(cm_mod, cep_mod.ecgleads.z_coords);
    cm.bcast(cm_mod, cep_mod.ecgleads.pseudo_ECG);
  }

  // Distribute output parameters
  //
  if (cm.slv(cm_mod)) {
    lEq.output.resize(lEq.nOutput);
  }

  for (int iOut = 0; iOut < lEq.nOutput; iOut++) {
    auto& output = lEq.output[iOut];
    cm.bcast_enum(cm_mod, &output.grp);
    cm.bcast(cm_mod, &output.o);
    cm.bcast(cm_mod, &output.l);
    cm.bcast(cm_mod, output.name);

    cm.bcast(cm_mod, &output.options.boundary_integral);
    cm.bcast(cm_mod, &output.options.spatial);
    cm.bcast(cm_mod, &output.options.volume_integral);
  }

  // Distribute BC information
  //
  if (cm.slv(cm_mod)) {
    lEq.bc.resize(lEq.nBc);
  }

  for (int iBc = 0;  iBc < lEq.nBc; iBc++) {
    dist_bc(com_mod, cm_mod, cm, lEq.bc[iBc], tMs, gmtl);
  }

  if (com_mod.ibFlag) {
    if (cm.slv(cm_mod)) {
      lEq.bcIB.resize(lEq.nBcIB);
    }
    for (int iBc = 0; iBc < lEq.nBcIB; iBc++) {
      // DISTBCIB(lEq.bcIB(iBc);)
    }
  }

  // Distribute BF information
  //
  if (cm.slv(cm_mod)) {
    lEq.bf.resize(lEq.nBf);
  }

  for (int iBf = 0; iBf < lEq.nBf; iBf++) {
    dist_bf(com_mod, cm_mod, cm, lEq.bf[iBf]);
  }

} 


/// @brief Distribute material properties to all processors.
//
void dist_mat_consts(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, stModelType& lStM)
{
  using namespace consts;

  #define n_debug_dist_mat_consts
  #ifdef debug_dist_mat_consts
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  cm.bcast_enum(cm_mod, &lStM.volType);
  cm.bcast(cm_mod, &lStM.Kpen);
  cm.bcast_enum(cm_mod, &lStM.isoType);
  cm.bcast(cm_mod, &lStM.C01);
  cm.bcast(cm_mod, &lStM.C10);
  cm.bcast(cm_mod, &lStM.a);
  cm.bcast(cm_mod, &lStM.b);
  cm.bcast(cm_mod, &lStM.aff);
  cm.bcast(cm_mod, &lStM.bff);
  cm.bcast(cm_mod, &lStM.ass);
  cm.bcast(cm_mod, &lStM.bss);
  cm.bcast(cm_mod, &lStM.afs);
  cm.bcast(cm_mod, &lStM.bfs);
  cm.bcast(cm_mod, &lStM.kap);
  cm.bcast(cm_mod, &lStM.khs);
  cm.bcast(cm_mod, &lStM.a0);
  cm.bcast(cm_mod, &lStM.b1);
  cm.bcast(cm_mod, &lStM.b2);
  cm.bcast(cm_mod, &lStM.mu0);

  // Distribute fiber stress
  cm.bcast(cm_mod, &lStM.Tf.fType);

  if (utils::btest(lStM.Tf.fType, static_cast<int>(BoundaryConditionType::bType_std))) { 
    cm.bcast(cm_mod, &lStM.Tf.g);

  } else if (utils::btest(lStM.Tf.fType, static_cast<int>(BoundaryConditionType::bType_ustd))) {
    cm.bcast(cm_mod, &lStM.Tf.gt.lrmp);
    cm.bcast(cm_mod, &lStM.Tf.gt.d);
    cm.bcast(cm_mod, &lStM.Tf.gt.n);

    if (cm.slv(cm_mod)) {
      int j = lStM.Tf.gt.d;
      int i = lStM.Tf.gt.n;
      lStM.Tf.gt.qi.resize(j);
      lStM.Tf.gt.qs.resize(j);
      lStM.Tf.gt.r.resize(j,i);
      lStM.Tf.gt.i.resize(j,i);
   } 

   cm.bcast(cm_mod, &lStM.Tf.gt.ti);
   cm.bcast(cm_mod, &lStM.Tf.gt.T);

   cm.bcast(cm_mod, lStM.Tf.gt.qi, "lStM.Tf.gt.qi");
   cm.bcast(cm_mod, lStM.Tf.gt.qs, "lStM.Tf.gt.qs");
   cm.bcast(cm_mod, lStM.Tf.gt.r, "lStM.Tf.gt.r");
   cm.bcast(cm_mod, lStM.Tf.gt.i, "lStM.Tf.gt.i");
  }
  
  // Broadcast directional stress distribution parameters
  cm.bcast(cm_mod, &lStM.Tf.eta_f);
  cm.bcast(cm_mod, &lStM.Tf.eta_s);
  cm.bcast(cm_mod, &lStM.Tf.eta_n);

  // Distribute CANN parameter table
  if (lStM.isoType == ConstitutiveModelType::stArtificialNeuralNet) {
    cm.bcast(cm_mod, &lStM.paramTable.num_rows);

    if (lStM.paramTable.num_rows > 0) {
      if (cm.slv(cm_mod)) {
        lStM.paramTable.invariant_indices.resize(lStM.paramTable.num_rows);
        lStM.paramTable.activation_functions.resize(lStM.paramTable.num_rows, 3);
        lStM.paramTable.weights.resize(lStM.paramTable.num_rows, 3);
      }
      cm.bcast(cm_mod, lStM.paramTable.invariant_indices);
      cm.bcast(cm_mod, lStM.paramTable.activation_functions, "paramTable.act_func");
      cm.bcast(cm_mod, lStM.paramTable.weights, "paramTable.weights");
    }
}
  
}


void dist_fluid_visc_model(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, fluidViscModelType& lVis)
{
  using namespace consts;

  cm.bcast_enum(cm_mod, &lVis.viscType);
  cm.bcast(cm_mod, &lVis.mu_i);
  cm.bcast(cm_mod, &lVis.mu_o);
  cm.bcast(cm_mod, &lVis.lam);
  cm.bcast(cm_mod, &lVis.a);
  cm.bcast(cm_mod, &lVis.n);
}

void dist_solid_visc_model(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, solidViscModelType& lVis)
{
  using namespace consts;

  cm.bcast_enum(cm_mod, &lVis.viscType);
  cm.bcast(cm_mod, &lVis.mu);
}


void part_face(Simulation* simulation, mshType& lM, faceType& lFa, faceType& gFa, Vector<int>& gmtl)
{
  #define n_debug_part_face
  #ifdef debug_part_face
  DebugMsg dmsg(__func__, simulation->com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "lFa.name: " << lFa.name;
  #endif

  auto& cm_mod = simulation->cm_mod;
  auto& com_mod = simulation->com_mod;
  auto& chnl_mod = simulation->chnl_mod;
  auto& cm = com_mod.cm;
  int num_proc = cm.np();

  // Broadcasting the number of nodes and elements of to slaves and
  // populating gFa to all procs
  if (cm.mas(cm_mod)) {
    gFa.d = lFa.d;
    gFa.eNoN = lFa.eNoN;
    gFa.iM = lFa.iM;
    gFa.nEl = lFa.nEl;
    gFa.gnEl = lFa.gnEl;
    gFa.nNo = lFa.nNo;
    gFa.qmTRI3 = lFa.qmTRI3;

    if (com_mod.rmsh.isReqd) {
      gFa.gebc.resize(1+gFa.eNoN, gFa.gnEl);
    }
  } else { 
    if (com_mod.rmsh.isReqd) {
      //ALLOCATE(gFa%gebc(0,0))
    }
  }
  cm.bcast(cm_mod, &lFa.qmTRI3);

  cm.bcast(cm_mod, &gFa.d);
  cm.bcast(cm_mod, &gFa.eNoN);
  cm.bcast(cm_mod, &gFa.iM);
  cm.bcast(cm_mod, &gFa.nEl);
  cm.bcast(cm_mod, &gFa.gnEl);
  cm.bcast(cm_mod, &gFa.nNo);
  cm.bcast(cm_mod, &gFa.qmTRI3);

  #ifdef debug_part_face
  dmsg << "gFa.d: " << gFa.d;
  dmsg << "gFa.eNoN: " << gFa.eNoN;
  dmsg << "gFa.iM: " << gFa.iM;
  dmsg << "gFa.nEl: " << gFa.nEl;
  dmsg << "gFa.gnEl: " << gFa.gnEl;
  dmsg << "gFa.nNo: " << gFa.nNo;
  dmsg << "gFa.qmTRI3: " << gFa.qmTRI3;
  #endif

  // Set face properties for the input element type.
  nn::select_eleb(simulation, lM, gFa);

  int eNoNb = gFa.eNoN;
  int iM = gFa.iM;

  gFa.IEN.resize(eNoNb, gFa.nEl);
  gFa.gE.resize(gFa.nEl);
  gFa.gN.resize(gFa.nNo);

  Vector<int> ePtr(gFa.nEl);

  // [NOTE Not sure about destroying the 'lFa' passed in parameter.
  //
  if (cm.mas(cm_mod)) {
    gFa = lFa;
    lFa.destroy();
  }

  cm.bcast(cm_mod, gFa.name);

  lFa.name = gFa.name;
  lFa.d = gFa.d;
  lFa.eNoN = eNoNb;

  nn::select_eleb(simulation, lM, lFa);
  lFa.iM = iM;

  // Be careful with 'i', it seems to be the number of something 
  // and not a counter.
  //
  int i = gFa.nEl*(2+eNoNb) + gFa.nNo;
  Vector<int> part(i);

  // Broadcast data in gFa from master to all processes
  if (cm.mas(cm_mod)) {
    for (int e = 0; e < gFa.nEl; e++) {
      int Ec = gFa.gE[e];
      ePtr[e] = lM.otnIEN[Ec];
      int j = e * (2+eNoNb);
      part[j] = Ec;
      part[j+1] = ePtr[e];
      for (int k = 0; k < gFa.IEN.nrows(); k++) {
        part[k+j+2] = gFa.IEN(k,e);
      }
    }

    for (int a = 0; a < gFa.nNo; a++) {
      int j = gFa.nEl*(2+eNoNb) + a;
      part[j] = gFa.gN[a];
    }
  }

  cm.bcast(cm_mod, part);

  if (cm.slv(cm_mod)) {
    for (int e = 0; e < gFa.nEl; e++) {
      int j = e * (2+eNoNb);
      gFa.gE[e] = part[j];
      ePtr[e] = part[j+1];
      for (int i = 0; i < gFa.IEN.nrows(); i++) {
        gFa.IEN(i,e) = part[i+j+2];
      }
    }

    for (int a = 0; a < gFa.nNo; a++) {
      int j = gFa.nEl * (2+eNoNb) + a;
      gFa.gN[a] = part[j];
    }
  }

  part.clear();

  // Finding the number of lM%fas to allocate required space, also
  // maping global element number to processor element number
  //
  lFa.nEl = 0;
  int task_id = cm.idcm();

  for (int e = 0; e < gFa.nEl; e++) {
    int Ec = ePtr[e];
    gFa.gE[e] = Ec;
    if ((Ec < lM.eDist[task_id+1]) && (Ec >= lM.eDist[task_id])) {
      lFa.nEl = lFa.nEl + 1;
    }
  }

  lFa.gE.resize(lFa.nEl); 
  lFa.IEN.resize(eNoNb,lFa.nEl);
  lFa.nNo = 0;

  for (int a = 0; a < gFa.nNo; a++) {
    int Ac = gmtl[gFa.gN[a]];
    if (Ac != -1) {
      lFa.nNo = lFa.nNo + 1;
    }
  }

  lFa.gN.resize(lFa.nNo);

  int j = 0;
  for (int e = 0; e < gFa.nEl; e++) {
    int Ec = gFa.gE[e];
    if ((Ec < lM.eDist[task_id+1]) && (Ec >= lM.eDist[task_id])) {
      lFa.gE[j] = Ec - lM.eDist[task_id];
      for (int a = 0; a < eNoNb; a++) {
        lFa.IEN(a,j) = gmtl(gFa.IEN(a,e));
      }
      j = j + 1;
    }
  }

  // Analogously copying the nodes which belong to this processor
  j = 0;
  for (int a = 0; a < gFa.nNo; a++) {
    int Ac = gmtl(gFa.gN[a]);
    if (Ac != -1) {
      lFa.gN[j] = Ac;
      j = j + 1;
    }
  }

  lFa.gnEl = gFa.gnEl;

  if (com_mod.rmsh.isReqd) {
    if (cm.mas(cm_mod)) {
      lFa.gebc.resize(1+eNoNb, lFa.gnEl);
      for (int e = 0; e < gFa.gnEl; e++) {
        lFa.gebc(0,e) = gFa.gebc(0,e);

        // [NOTE] perhaps another way to do this using a range.
        // std::array<int,2> rows{1,eNoNb};
        // auto col = gFa.gebc.get_col(e, rows);
        // lFa.gebc.set_col(e, col, rows);

        for (int i = 1; i <= eNoNb; i++) {
          lFa.gebc(i,e) = gFa.gebc(i,e);
        }
      }
    } else { 
      lFa.gebc.clear();
    }
  }
}


/// @brief Reproduces the Fortran 'PARTMSH' subroutine.
/// Parameters for the part_msh function:
/// @param[in] simulation A pointer to the simulation object.
/// @param[in] iM The mesh index.
/// @param[in] lM The local mesh data.
/// @param[in] gmtl The global to local map.
/// @param[in] nP The number of processors.
/// @param[in] wgt The weights.
//
void part_msh(Simulation* simulation, int iM, mshType& lM, Vector<int>& gmtl, int nP, Vector<float>& wgt)
{
  auto& cm_mod = simulation->cm_mod;
  auto& com_mod = simulation->com_mod;
  auto& chnl_mod = simulation->chnl_mod;
  auto& cm = com_mod.cm;
  int num_proc = cm.np();
  int task_id = cm.idcm();

  #define n_dbg_part_msh
  #ifdef dbg_part_msh
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  // A single process simulation.
  //
  if (cm.seq()) {
    lM.nEl = lM.gnEl;
    lM.nNo = lM.gnNo;

    lM.IEN.resize(lM.eNoN, lM.nEl); 
    lM.eDist.resize(num_proc+1);

    lM.IEN = lM.gIEN;
    lM.eDist(0) = 0;
    lM.eDist(1) = lM.gnEl;

    lM.otnIEN.resize(lM.nEl);

    for (int e = 0; e < lM.nEl; e++) {
      lM.otnIEN[e] = e;
    }

    lM.iGC.resize(lM.nEl);
    return;
  }

  // Broadcast mesh data.
  //
  cm.bcast(cm_mod, &lM.lShpF);
  cm.bcast(cm_mod, &lM.lShl);
  cm.bcast(cm_mod, &lM.lFib);

  int eType = static_cast<int>(lM.eType);
  cm.bcast(cm_mod, &eType);

  cm.bcast(cm_mod, &lM.eNoN);
  cm.bcast(cm_mod, &lM.nFa);
  cm.bcast(cm_mod, &lM.nFs);
  cm.bcast(cm_mod, &lM.nG);
  cm.bcast(cm_mod, &lM.gnEl);
  cm.bcast(cm_mod, &lM.gnNo);
  cm.bcast(cm_mod, lM.name);
  cm.bcast(cm_mod, &lM.nFn);
  cm.bcast(cm_mod, &lM.scF);
  cm.bcast(cm_mod, &lM.qmTET4);

  // Number of fibers.
  int nFn = lM.nFn;

  // Set integration dimension.
  int nsd = com_mod.nsd;
  int insd = nsd;
  if (lM.lShl) {
    insd = nsd - 1;
  }
  if (lM.lFib) {
    insd = 1;
  }

  int eNoN = lM.eNoN;

  #ifdef dbg_part_msh
  dmsg << "lM.gnEl: " << lM.gnEl;
  dmsg << "lM.gnNo: " << lM.gnNo;
  #endif

  if (cm.slv(cm_mod)) {
    nn::select_ele(com_mod, lM);
    lM.gIEN.clear(); 
    lM.fa.resize(lM.nFa);
  }

  Vector<int> sCount(num_proc); 
  Vector<int> disp(num_proc); 

  // [NOTE] lM.eDist[] in Fortran starts from 0.
  lM.eDist.resize(num_proc+1);

  // And distributing bs for NURBS
  //
  if (lM.eType == consts::ElementType::NRB) { 
    /*
    IF (cm%slv()) ALLOCATE(lM%bs(insd))
     cm%bcast(lM%nSl)
    DO i=1, insd
       cm%bcast(lM%bs(i)%n)
       cm%bcast(lM%bs(i)%nG)
       cm%bcast(lM%bs(i)%nEl)
       cm%bcast(lM%bs(i)%nSl)
       cm%bcast(lM%bs(i)%p)
      IF (cm%slv()) ALLOCATE(lM%bs(i)%xi(lM%bs(i)%n))
       cm%bcast(lM%bs(i)%xi)
      lM%bs(i)%nNo = lM%bs(i)%n - lM%bs(i)%p - 1
    END DO

    a = lM%bs(2)%nEl
    IF (insd .EQ. 3) a = a*lM%bs(3)%nEl
    DO i=0, cm%np()
      lM%eDist(i) = a*NINT(SUM(wgt(1:i))*lM%bs(1)%nEl, KIND=IKIND)
      IF (lM%eDist(i) .GT. lM%gnEl) lM%eDist(i) = lM%gnEl
    END DO
    */

    //  A draft of splitting the mesh between processors
    //  lM%eDist(i) represents first element which belong to cm%id()=i
  } else { 
    for (int i = 0; i < lM.eDist.size(); i++) {
      double sum = 0.0;
      for (int j = 0; j < i; j++) {
        sum += wgt[j];
      }
      lM.eDist[i] = round(sum * lM.gnEl);
      if (lM.eDist[i] > lM.gnEl) {
        lM.eDist[i] = lM.gnEl;
      }
    }
  }

  lM.eDist[num_proc] = lM.gnEl;
  #ifdef dbg_part_msh
  dmsg << "wgt: " << wgt; 
  dmsg << "lM.eDist: " << lM.eDist;
  #endif

  for (int i = 0; i < num_proc; i++) { 
    disp[i] = lM.eDist[i] * eNoN;
    sCount[i] = lM.eDist[i+1] * eNoN - disp[i];
    #ifdef dbg_part_msh
    dmsg << ">>>> i: " << i;
    dmsg << "disp[i]: " << disp[i];
    dmsg << "sCount[i]: " << sCount[i];
    #endif
  }

  int nEl = lM.eDist(cm.id() + 1) - lM.eDist(cm.id());
  int idisp = lM.eDist(cm.id()) * sizeof(nEl);
  #ifdef dbg_part_msh
  dmsg << "cm.id(): " << cm.id();
  dmsg << "nEl: " << nEl;
  dmsg << "sizeof(nEl): " << sizeof(nEl);
  dmsg << "idisp: " << idisp;
  dmsg << "eNoN: " << eNoN;
  #endif

  Vector<int> part(nEl);

  std::string fTmp = chnl_mod.appPath + "partitioning_" + lM.name + ".bin";
  bool flag = false;
  FILE *fp = nullptr; 
  if (com_mod.rmsh.isReqd) { 
    fp = fopen(fTmp.c_str(), "r");
    if (FILE *fp = fopen(fTmp.c_str(), "r")) {
      flag = true;
      fclose(fp);
    }
  }
  //IF (rmsh%isReqd) INQUIRE(FILE=TRIM(fTmp), EXIST=flag)
  #ifdef dbg_part_msh
  dmsg << " " << " ";
  dmsg << "rmsh.isReqd: " << com_mod.rmsh.isReqd;
  dmsg << "fTmp: " << fTmp;
  dmsg << "flag: " << flag;
  dmsg << "com_mod.resetSim: " << com_mod.resetSim;
  dmsg << "fp: " << fp;
  if (fp) fclose(fp);
  #endif


  if (lM.eType == consts::ElementType::NRB) {
    part = cm.id();

  // [TODO:DaveP] Reading partition data does not seem to work.
  //
  } else if (false) { 
  //} else if (flag && !com_mod.resetSim) {
    #ifdef dbg_part_msh
    dmsg << " " << " ";
    dmsg << "Reading partition data from file " << fTmp;
    #endif
    MPI_File fid;
    MPI_File_open(cm.com(), fTmp.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fid);
    MPI_File_set_view(fid, idisp, cm_mod::mpint, cm_mod::mpint, "native", MPI_INFO_NULL);
    MPI_File_read(fid, part.data(), nEl, cm_mod::mpint, MPI_STATUS_IGNORE);
    MPI_File_close(&fid);
    #ifdef dbg_part_msh
    dmsg << "  nEl: " << nEl;
    dmsg << "  idisp: " << idisp;
    //if (cm.mas(cm_mod)) {
    //dmsg << "  part = " << part;
    //}
    dmsg << "---------- " << "---------- ";
    #endif

  // Scattering the lM.gIEN array to all processors.
  //
  } else { 
    #ifdef dbg_part_msh
    dmsg << " " << " ";
    dmsg << " " << " ";
    dmsg << "Scattering the lM%gIEN array to processors " << " ... ";
    #endif
    lM.IEN.resize(eNoN, nEl);


    // Send lM.gIEN array to all processor's lM.IEN[] array of siize nEl*eNoN.
    //
    #ifdef dbg_part_msh
    dmsg << "sCount: " << sCount;
    dmsg << "disp: " << disp;
    #endif
    MPI_Scatterv(lM.gIEN.data(), sCount.data(), disp.data(), cm_mod::mpint, lM.IEN.data(), 
        nEl*eNoN, cm_mod::mpint, cm_mod.master, cm.com());

    int eNoNb = consts::element_type_to_elem_nonb.at(lM.eType);
    #ifdef dbg_part_msh
    dmsg << "nEl: " << nEl;
    dmsg << "eNoNb: " << eNoNb;
    dmsg << "lM.IEN.size(): " << lM.IEN.size();
    #endif


    // The output of this process is "part" array which part(i) says
    // which processor element "i" belongs to
    // Doing partitioning, using ParMetis
    //
    auto edgecut = split_(&nEl, &eNoN, &eNoNb, lM.IEN.data(), &num_proc, lM.eDist.data(),  wgt.data(), part.data());
    #ifdef dbg_part_msh
    dmsg << "edgecut: " << edgecut;
    #endif

    if (edgecut == 0) {
      #ifdef dbg_part_msh
      dmsg << "ParMETIS failed to partition the mesh" << " ...";
      #endif
      part = cm.id();
    } else if (edgecut > 0) {
      #ifdef dbg_part_msh
      dmsg << "ParMETIS partitioned the mesh by cutting: " << std::to_string(edgecut) + " elements.";
      #endif
    } 

    lM.IEN.clear();

    if (com_mod.rmsh.isReqd) {
      #ifdef dbg_part_msh
      dmsg << "---------------------------" << "------ ";
      dmsg << "Writing partition data to file: " << fTmp;
      dmsg << "  nEl: " << nEl;
      dmsg << "  idisp: " << idisp;
      dmsg << "  idisp: " << idisp;
      //if (cm.mas(cm_mod)) {
      //dmsg << "  part = " << part;
      //}
      dmsg << "---------------------------" << "------ ";
      #endif
      MPI_File fid;
      MPI_File_open(cm.com(), fTmp.c_str(), MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, &fid);
      MPI_File_set_view(fid, idisp, cm_mod::mpint, cm_mod::mpint, "native", MPI_INFO_NULL);
      MPI_File_write(fid, part.data(), nEl, cm_mod::mpint, MPI_STATUS_IGNORE);
      MPI_File_close(&fid);
    }
  }

  for (int i = 0; i < num_proc; i++) { 
    disp[i] = lM.eDist[i];
    sCount[i] = lM.eDist[i+1] - disp[i];
  }

  // Gathering the parts inside master, part(e) is equal to the
  // cm%id() that the element e belong to.
  //
  Vector<int> gPart;
  if (cm.mas(cm_mod)) {
    gPart.resize(lM.gnEl);
  } else { 
    gPart.clear();
  } 

  // gpart is a global version of part in which processor p = gpart(e)
  // is the owner of element "e".
  //
  MPI_Gatherv(part.data(), nEl, cm_mod::mpint, gPart.data(), sCount.data(), disp.data(), 
      cm_mod::mpint, cm_mod.master, cm.com());

  part.clear();

  Array<int> tempIEN;
  Array<double> tmpFn;
  flag = false;
  bool fnFlag = false;

  #ifdef dbg_part_msh
  dmsg << "sCount: " << sCount;
  dmsg << "disp: " << disp;
  //dmsg << "gPart: " << gPart;
  dmsg << " " << " ";
  dmsg << "Making the lM%IEN array " << " ...";
  #endif

  if (cm.mas(cm_mod)) {
    sCount = 0;
    int eRisProc = -1;
    for (int e = 0; e < lM.gnEl; e++) {
      if (com_mod.risFlag) {
        // If we found that this element needs to be shared
        // since elements on the other mesh had been assigned 
        // a processor, then we change the gPart id
        if (lM.partRIS(e) != -1) {
          gPart[e] = lM.partRIS(e);
        } else if (lM.eRIS(e)) {
          // If this element is on a ris projection,
          // we record the processor id so that next mesh
          // will know that this element needs to be shared.
          // Ideally, we might want to handle the case where
          // elements next to the same projection are sent
          // to different processors. Yet, this doesn't often
          // happen and for simplicity, we
          // force the processor id to be the same for all 
          // elements next to the same projection.
          if (eRisProc == -1) {
            eRisProc = gPart[e];
          } else {
            gPart[e] = eRisProc;
          }
          lM.partRIS(e) = eRisProc;
        }
      }
      sCount[gPart[e]] = sCount[gPart[e]] + 1;
    }

    for (int i = 0; i < num_proc; i++) { 
      lM.eDist[i+1] = lM.eDist[i] + sCount[i];
    }

    #ifdef dbg_part_msh
    dmsg << "lM.eDist: " << lM.eDist;
    dmsg << "sCount: " << sCount;
    #endif

    tempIEN.resize(eNoN,lM.gnEl); 
    lM.otnIEN.resize(lM.gnEl);

    // Making the lM%IEN array in order, based on the cm%id() number in
    // master. lM%otnIEN maps old IEN order to new IEN order.
    //
    disp = 0;

    for (int e = 0; e < lM.gnEl; e++) { 
      int Ec = lM.eDist[gPart[e]];
      lM.eDist[gPart[e]] = Ec + 1;
      tempIEN.set_col(Ec, lM.gIEN.col(e));
      lM.otnIEN[e] = Ec;
    }

    lM.gIEN = tempIEN;
    lM.eDist[0] = 0;
    for (int i = 0; i < num_proc; i++) { 
      lM.eDist[i+1] = lM.eDist[i] + sCount[i];
    }
    #ifdef dbg_part_msh
    dmsg << "2 lM.eDist: " << lM.eDist;
    #endif

    // This it to distribute eId, if allocated
    //
    if (lM.eId.size() != 0) {
      flag = true;
      part.resize(lM.gnEl);
      for (int e = 0; e < lM.gnEl; e++) {
        int Ec = lM.otnIEN[e];
        part[Ec] = lM.eId[e];
      }
      lM.eId.clear();
    }

    // This it to distribute fN, if allocated
    if (lM.fN.size() != 0) {
      fnFlag = true;
      tmpFn.resize(nFn*nsd, lM.gnEl);
      for (int e = 0; e < lM.gnEl; e++) {
        int Ec = lM.otnIEN[e];
        tmpFn.set_col(Ec, lM.fN.col(e));
      }
      lM.fN.clear();
    }

  } else { 
    lM.otnIEN.clear();
  }

  gPart.clear();

  cm.bcast(cm_mod, &flag);
  cm.bcast(cm_mod, &fnFlag);
  cm.bcast(cm_mod, lM.eDist);
  if (com_mod.risFlag) {
    cm.bcast(cm_mod, lM.partRIS);
  }

  nEl = lM.eDist[cm.id()+1] - lM.eDist[cm.id()];
  #ifdef dbg_part_msh
  dmsg << "flag: " << flag;
  dmsg << "fnFlag: " << fnFlag;
  dmsg << "3 lM.eDist: " << lM.eDist;
  #endif
  lM.nEl = nEl;
  lM.IEN.resize(eNoN,nEl); 
  lM.iGC.resize(nEl);

  // Communicating eId, if neccessary.
  //
  if (flag) {
    lM.eId.resize(nEl);
    // [NOTE] Not sure what's going on here.
    if (part.size() == 0) {
      part.clear();
    }
    for (int i = 0; i < num_proc; i++) { 
      disp[i] = lM.eDist[i];
      sCount[i] = lM.eDist[i+1] - disp[i];
    }
    MPI_Scatterv(part.data(), sCount.data(), disp.data(), cm_mod::mpint, lM.eId.data(), nEl, 
        cm_mod::mpint, cm_mod.master, cm.com());
    part.clear();
  }

  // Communicating fN, if neccessary
  //
  if (fnFlag) { 
    #ifdef dbg_part_msh
    dmsg << "Communicating fN " << " ...";
    dmsg << "nFn: " << nFn;
    dmsg << "nsd: " << nsd;
    dmsg << "nEl: " << nEl;
    dmsg << "tmpFn.size(): " << tmpFn.size();
    #endif
    lM.fN.resize(nFn*nsd,nEl);
    if (tmpFn.size() == 0) {
      // ALLOCATE(tmpFn(0,0))
    }
    for (int i = 0; i < num_proc; i++) { 
      disp[i] = lM.eDist[i] * nFn * nsd;
      sCount[i] = lM.eDist[i+1] * nFn * nsd - disp[i];
    }
    MPI_Scatterv(tmpFn.data(), sCount.data(), disp.data(), cm_mod::mpreal, lM.fN.data(), nEl*nFn*nsd, 
        cm_mod::mpreal, cm_mod.master, cm.com());
    tmpFn.clear();
  }

  // Now scattering the sorted lM%IEN to all processors.
  //
  #ifdef dbg_part_msh
  dmsg << " " << " ";
  dmsg << "Now scattering the sorted lM%IEN to all processors " << " ...";
  #endif
  if (tempIEN.size() == 0) {
    tempIEN.clear();
  }
  for (int i = 0; i < num_proc; i++) {
    disp[i] = lM.eDist[i] * eNoN;
    sCount[i] = lM.eDist[i+1]*eNoN - disp[i];
  }

  MPI_Scatterv(tempIEN.data(), sCount.data(), disp.data(), cm_mod::mpint, lM.IEN.data(), 
      nEl*eNoN, cm_mod::mpint, cm_mod.master, cm.com());

  tempIEN.clear();

  // Constructing the initial global to local pointer
  // lM%IEN: eNoN,nEl --> gnNo
  // gtlPtr: gnNo     --> nNo
  // lM%IEN: eNoN,nEl --> nNo
  //
  #ifdef dbg_part_msh
  dmsg << " " << " ";
  dmsg << "Constructing the initial global to local pointer " << " ...";
  dmsg << "lM.gnNo: " << lM.gnNo;
  #endif

  Vector<int> gtlPtr(lM.gnNo);
  int nNo = 0;
  gtlPtr = -1;

  for (int e = 0; e < nEl; e++) {
    for (int a = 0; a < eNoN; a++) {
      int Ac = lM.IEN(a,e);
      if (gtlPtr[Ac] == -1) {
        gtlPtr[Ac] = nNo;
        nNo = nNo + 1;
      }
      lM.IEN(a,e) = gtlPtr[Ac];
    }
  }

  #ifdef dbg_part_msh
  dmsg << "nNo: " << nNo;
  #endif

  lM.nNo = nNo;
  if (cm.slv(cm_mod)) {
    lM.gN.resize(lM.gnNo);
  }
  cm.bcast(cm_mod, lM.gN);

  // lM%gN: gnNo --> gtnNo
  // part:  nNo  --> gtnNo
  part.resize(nNo);
  for (int Ac = 0; Ac < lM.gnNo; Ac++) {
    int a = gtlPtr[Ac];
    if (a != -1) {
      part[a] = lM.gN[Ac];
    }
  }

  #ifdef dbg_part_msh
  dmsg << "---------- partitioned mesh " << std::to_string(iM) + " ----------";
  dmsg << "lM.nNo: " << lM.nNo;
  dmsg << "lM.gnNo: " << lM.gnNo;
  dmsg << "lM.nEl: " << lM.nEl;
  dmsg << "lM.gnEl: " << lM.gnEl;
  #endif

  // mapping and converting other parameters.
  // I will use an upper bound for gPart as a container for ltg,
  // since there can be repeated nodes. gPart is just a temp variable.
  // gmtl:  gtnNo --> tnNo
  // gPart: tnNo  --> gtnNo
  // ltg:   tnNo  --> gtnNo
  // lM%gN: nNo   --> tnNo
  //
  #ifdef dbg_part_msh
  dmsg << "Mapping and converting other parameters " << " ... ";
  dmsg << "com_mod.tnNo: " << com_mod.tnNo;
  #endif

  lM.gN.clear();
  gPart.resize(com_mod.tnNo + nNo); 
  lM.gN.resize(nNo);
  gmtl = -1;

  for (int a = 0; a < com_mod.tnNo; a++) {
    int Ac = com_mod.ltg[a];
    gPart[a] = Ac;
    gmtl[Ac] = a;
  } 

  int tnNo = com_mod.tnNo;
  for (int a = 0; a < nNo; a++) {
    int Ac = part[a];
    if (gmtl[Ac] == -1) {
      gmtl[Ac] = tnNo;
      lM.gN[a] = tnNo;
      gPart(tnNo) = Ac;
      tnNo += 1;
    } else { 
      lM.gN[a] = gmtl[Ac];
    }
  }
  #ifdef dbg_part_msh
  dmsg << "tnNo: " << tnNo;
  #endif

  com_mod.ltg.resize(tnNo);

  for (int i = 0; i < tnNo; i++) {
    com_mod.ltg[i] = gPart(i);
  }

  gPart.clear();
  com_mod.tnNo = tnNo;

  // If neccessary communicate NURBS
  //
  // [TODO] Not implemented.
  //
  if (lM.eType == consts::ElementType::NRB) {
    /*
    ALLOCATE(tmpR(lM%gnNo))
    IF (cm%mas()) THEN
      tmpR = lM%nW
      DEALLOCATE(lM%nW)
    END IF
     cm%bcast(tmpR)
    ALLOCATE(lM%nW(lM%nNo))
    DO Ac=1, lM%gnNo
      a = gtlPtr(Ac)
      IF (a .NE. 0) THEN
        lM%nW(a) = tmpR(Ac)
      END IF
    END DO

    // Distributing INN, using tempIEN as tmp array
    IF (cm%mas()) THEN
      ALLOCATE(tempIEN(insd,lM%gnEl))
      DO e=1, lM%gnEl
        Ec = lM%otnIEN(e)
        tempIEN(:,Ec) = lM%INN(:,e)
      END DO
      DEALLOCATE(lM%INN)
    ELSE
      ALLOCATE(tempIEN(0,0))
    END IF
    DO i=1, cm%np()
      disp(i)   = lM%eDist(i-1)*insd
      sCount(i) = lM%eDist(i)*insd - disp(i)
    END DO
    ALLOCATE(lM%INN(insd,nEl))

    //  Now scattering the sorted lM%INN to all processors
     MPI_SCATTERV(tempIEN, sCount, disp, mpint, lM%INN,  nEl*insd, mpint, master, cm%com(), ierr)
    */
  }
  // If necessary, distribute precomputed state-variable data.
  //
  flag = (lM.Ys.size() != 0);
  cm.bcast(cm_mod, &flag);
  if (flag){
    #ifdef dbg_part_msh
    dmsg << "Distributing precomputed state-variable data " << " ...";
    #endif
    Array3<double> tmpYs;
    int nsYs = lM.Ys.nslices();
    if (cm.mas(cm_mod)) {
      tmpYs.resize(lM.Ys.nrows(), lM.Ys.ncols(), nsYs);
      tmpYs = lM.Ys;
      lM.Ys.clear();
    } else {
      tmpYs.clear();
    }
    lM.Ys.resize(com_mod.nsd, com_mod.tnNo, nsYs);
    lM.Ys = all_fun::local(com_mod, cm_mod, cm, tmpYs);
    tmpYs.clear();
  }
}