// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// The functions here reproduce the subroutines defined in svFSILS/LHS.f.

#include "lhs.h"
#include "CmMod.h"
#include "DebugMsg.h"

#include "mpi.h"

namespace fsi_linear_solver {

/// @brief Modifies:
///
///  lhs.foC 
///  lhs.gnNo 
///  lhs.nNo 
///  lhs.nnz 
///  lhs.commu 
///  lhs.nFaces 
///  lhs.mynNo 
///
///  lhs.colPtr
///  lhs.rowPtr
///  lhs.diagPtr
///  lhs.map
///  lhs.face
//
void fsils_lhs_create(FSILS_lhsType& lhs, FSILS_commuType& commu, int gnNo, int nNo, int nnz, Vector<int>& gNodes,  
       Vector<int> &rowPtr, Vector<int>& colPtr, int nFaces)
{
  #define n_debug_fsils_lhs_create
  #ifdef debug_fsils_lhs_create
  DebugMsg dmsg(__func__,  lhs.commu.task);
  dmsg.banner();
  #endif

  lhs.foC = true; 
  lhs.gnNo = gnNo;
  lhs.nNo = nNo;
  lhs.nnz = nnz;
  lhs.commu = commu;
  lhs.nFaces = nFaces;
  #ifdef debug_fsils_lhs_create
  dmsg << "gnNo: " << gnNo;
  dmsg << "nNo: " << nNo;
  dmsg << "nnz: " << nnz;
  dmsg << "nFaces: " << nFaces;
  #endif

  int nTasks = commu.nTasks;
  auto comm = commu.comm;
  auto tF = commu.tF;
  #ifdef debug_fsils_lhs_create
  dmsg << "nTasks: " << nTasks;
  dmsg << "tF: " << tF;
  #endif

  lhs.colPtr.resize(nnz); 
  lhs.rowPtr.resize(2,nNo); 
  lhs.diagPtr.resize(nNo);
  lhs.map.resize(nNo); 
  lhs.face.resize(nFaces);

  // For a sequential simulation. 
  //
  if (nTasks == 1) {

    for (int i = 0; i < nnz; i++) {
      lhs.colPtr(i) = colPtr(i);
    }

    for (int Ac = 0; Ac < nNo; Ac++) {
      int s = rowPtr(Ac);
      int e = rowPtr(Ac+1) - 1;

      for (int i = s; i <= e; i++) {
        int a = colPtr(i);
        if (Ac == a) {
          lhs.diagPtr(Ac) = i;
          break; 
        }
      }

      lhs.rowPtr(0,Ac) = s;
      lhs.rowPtr(1,Ac) = e;
      lhs.map(Ac) = Ac;
    }

    lhs.mynNo = nNo;
    return; 
  }

  // Get the number of nodes for this process and the max number of 
  // nodes for all processes.
  //
  // Combines values from all processes and distributes the result back to all processes.
  //
  int maxnNo;
  MPI_Allreduce(&nNo, &maxnNo, 1, cm_mod::mpint, MPI_MAX, comm);

  #ifdef debug_fsils_lhs_create
  dmsg << "nNo: " << nNo;
  dmsg << "maxnNo: " << maxnNo;
  #endif

  // Initialize some data structures.
  //
  Vector<int> part(maxnNo); 
  part = -1;

  for (int i = 0; i < nNo; i++) {
    part(i) = gNodes(i);
  }

  Vector<int> sCount(nTasks); 
  Vector<int> disp(nTasks); 

  for (int i = 0; i < nTasks; i++) {
     disp(i) = i*maxnNo;
     sCount(i) = maxnNo;
  }

  Array<int> aNodes(maxnNo,nTasks); 
  Vector<int> ltg(nNo);

  // Gather data from all processes and deliver it to all. Each process may contribute a different amount of data.
  MPI_Allgatherv(part.data(), maxnNo, cm_mod::mpint, aNodes.data(), sCount.data(), disp.data(), cm_mod::mpint, comm);

  Vector<int> gtlPtr(gnNo); 
  gtlPtr = -1;

  for (int a = 0; a < nNo; a++) {
     int Ac = gNodes(a);
     gtlPtr(Ac) = a;
  }

  // Including the nodes shared by processors with higher ID at the end,
  // and including the nodes shared by lower processors IDs at the front.
  // shnNo is counter for lower ID and mynNo is counter for higher ID
  //
  lhs.mynNo = nNo;
  lhs.shnNo = 0;
  #ifdef debug_fsils_lhs_create
  dmsg << "lhs.mynNo: " << lhs.mynNo;
  dmsg << "tF: " << tF;
  #endif

  for (int i = nTasks-1; i >= 0; i--) { 
    // Will include local nodes later
    if (i == tF) {
      continue; 
    }

    for (int a = 0; a < maxnNo; a++) { 
      // Global node number in processor i at location a
      int Ac = aNodes(a,i);
      // Exit if this is the last node
      if (Ac == -1) {
        break; 
      }

      // Corresponding local node in current processor.
      int ai = gtlPtr(Ac);
      if (ai != -1) {
        // If this node has not been included already
        if (aNodes(ai,tF) != -1) {
          // If the processor ID is lower, it is appended to the beginning
          if (i < tF) {
            ltg(lhs.shnNo) = Ac;
            lhs.shnNo = lhs.shnNo + 1;
          // If the processor ID is higher, it is appended to the end
          } else {
            ltg(lhs.mynNo-1) = Ac;
            lhs.mynNo = lhs.mynNo - 1;
          }
          aNodes(ai,tF) = -1;
        }
      }
    }
  }

  #ifdef debug_fsils_lhs_create
  dmsg << "lhs.shnNo: " << lhs.shnNo;
  dmsg << "lhs.mynNo: " << lhs.mynNo;
  #endif

  // Now including the local nodes that are left behind
  //
  int j = lhs.shnNo + 1;

  for (int a = 0; a < nNo; a++) {
    int Ac = aNodes(a,tF);
    //  If this node has not been included already
    if (Ac != -1) {
      ltg(j-1) = Ac;
      j = j + 1;
    }
  }

  if (j != lhs.mynNo+1) {
    throw std::runtime_error("FSILS: Unexpected behavior: j=" + std::to_string(j) + " lhs.mynNo: " + std::to_string(lhs.mynNo) + ".");
    MPI_Finalize();
  }

  // Having the new ltg pointer, map is constructed
  //
  gtlPtr = -1;
  for (int a = 0; a < nNo; a++) {
     int Ac = ltg(a);
     gtlPtr(Ac) = a;
  }

  for (int a = 0; a < nNo; a++) {
     int Ac = gNodes(a);
     lhs.map(a) = gtlPtr(Ac);
  }

  // Based on the new ordering of the nodes, rowPtr and colPtr are constructed
  //
  for (int a = 0; a < nNo; a++) {
    int Ac = lhs.map(a);
    lhs.rowPtr(0,Ac) = rowPtr(a);
    lhs.rowPtr(1,Ac) = rowPtr(a+1) - 1;
  }

  for (int i = 0; i < nnz; i++) {
    lhs.colPtr(i) = lhs.map(colPtr(i));
  }

  // diagPtr points to the diagonal entries of LHS
  for (int Ac = 0; Ac < nNo; Ac++) {
    for (int i = lhs.rowPtr(0,Ac); i <= lhs.rowPtr(1,Ac); i++) {
      int a = lhs.colPtr(i);
      if (Ac == a) {
        lhs.diagPtr(Ac) = i;
        break; 
      }
    }
  }

  // Constructing the communication data structure based on the ltg
  //
  for (int i = 0; i < nNo; i++) {
    part(i) = ltg(i);
  }

  MPI_Allgatherv(part.data(), maxnNo, cm_mod::mpint, aNodes.data(), sCount.data(), disp.data(), cm_mod::mpint, comm);

  // This variable keeps track of number of shared nodes
  disp = 0;
  lhs.nReq = 0;

  for (int i = 0; i < nTasks; i++) {
    if (i == tF) {
      continue; 
    }
    for (int a = 0; a <  maxnNo; a++) {
      // Global node number in processor i at location a
      int Ac = aNodes(a,i);
      // Exit if this is the last node
      if (Ac == -1) {
        break;
      }
      // Corresponding local node in current processor
      int ai = gtlPtr(Ac);
      if (ai != -1) {
        disp(i) = disp(i) + 1;
      }
    }

    if (disp(i) != 0) {
      lhs.nReq = lhs.nReq + 1;
    }
  }
  #ifdef debug_fsils_lhs_create
  dmsg << "lhs.nReq: " << lhs.nReq;
  #endif

  lhs.cS.resize(lhs.nReq);

  // Now that we know which processor is communicating to which, we can
  // setup the handles and structures
  //
  // lhs.cS[j].iP is the processor to communicate with.
  //
  #ifdef debug_fsils_lhs_create
  dmsg << "Setup the handles ...";
  #endif
  j = 0;
  for (int i = 0; i < nTasks; i++) {
    int a = disp(i);
    if (a != 0) {
      lhs.cS[j].iP = i;
      lhs.cS[j].n = a;
      lhs.cS[j].ptr.resize(a);
      j = j + 1;
    }
  }

  // Order of nodes in ptr is based on the node order in processor
  // with higher ID. ptr is calculated for tF+1:nTasks and will be
  // sent over.
  //
  #ifdef debug_fsils_lhs_create
  dmsg << "Order of nodes ...";
  #endif
  MPI_Status status;

  for (int i = 0; i < lhs.nReq; i++) {
    int iP = lhs.cS[i].iP;

    if (iP < tF) {
      MPI_Recv(lhs.cS[i].ptr.data(), lhs.cS[i].n, cm_mod::mpint, iP, 1,  comm, &status);

      for (int j = 0; j < lhs.cS[i].n; j++) { 
        lhs.cS[i].ptr[j] = gtlPtr(lhs.cS[i].ptr[j]);
      }
    } else {
      // This is a counter for the shared nodes
      j = 0;
      for (int a = 0; a < maxnNo; a++) {
        // Global node number in processor i at location a
        int Ac = aNodes(a,iP);
        // Exit if this is the last node
        if (Ac == -1) {
          break; 
        }
        // Corresponding local node in current processor
        int ai = gtlPtr(Ac);
        if (ai != -1) {
          // Just for now global node ID is used. Later on this will be changed
          // to make sure nodes corresponds to each other on both processors
          // and then will be transformed to local node IDs
          lhs.cS[i].ptr[j] = Ac;
          j = j + 1;
        }
      }

      MPI_Send(lhs.cS[i].ptr.data(), lhs.cS[i].n, cm_mod::mpint, iP, 1, comm);

      for (int j = 0; j < lhs.cS[i].n; j++) {
        lhs.cS[i].ptr[j] = gtlPtr(lhs.cS[i].ptr[j]);
      }
    }
  }
}

//----------------
// fsils_lhs_free
//----------------
//
void fsils_lhs_free(FSILS_lhsType& lhs)
{
  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    if (lhs.face[faIn].foC) {
      fsils_bc_free(lhs, faIn);
    }
  }

  for (int i = 0; i < lhs.nReq; i++) {
    //IF (ALLOCATED(lhs.cS(i).ptr)) DEALLOCATE(lhs.cS(i).ptr)
  }

  lhs.foC = false;
  lhs.gnNo   = 0;
  lhs.nNo    = 0;
  lhs.nnz    = 0;
  lhs.nFaces = 0;

  //IF (ALLOCATED(lhs.colPtr)) DEALLOCATE(lhs.colPtr)
  //IF (ALLOCATED(lhs.rowPtr)) DEALLOCATE(lhs.rowPtr)
  //IF (ALLOCATED(lhs.diagPtr)) DEALLOCATE(lhs.diagPtr)
  //IF (ALLOCATED(lhs.cS)) DEALLOCATE(lhs.cS)
  //IF (ALLOCATED(lhs.map)) DEALLOCATE(lhs.map)
  //IF (ALLOCATED(lhs.face)) DEALLOCATE(lhs.face)

}


};


