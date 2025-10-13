// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// The functions here reproduce the subroutines defined in svFSILS/INCOMMU.f.

// To syncronize the data on the boundaries between the processors.
//
// This a both way communication with three main part:
// 1 - rTmp {in master} = R          {from slave}
// 2 - R    {in master} = R + rTmp   {both from master}
// 3 - rTmp {in master} = R          {from master}
// 4 - R    {in slave}  = rTmp       {from master}

#include "fsils.hpp"
#include "CmMod.h"
#include "Array3.h"

#include "fsils_std.h"

namespace fsi_linear_solver {

void fsils_commus(const FSILS_lhsType& lhs, Vector<double>& R)
{
  if (lhs.commu.nTasks == 1) {
    return;
  }

  int nReq = lhs.nReq;
  if (lhs.cS.size() == 0) {
    return;
  }

  int nmax = std::max_element(lhs.cS.begin(), lhs.cS.end(), 
    [](const FSILS_cSType& a, const FSILS_cSType& b){return a.n < b.n;})->n; 

  Array<double> sB(nmax,nReq); 
  Array<double> rB(nmax,nReq); 
  std::vector<MPI_Request> rReq(nReq); 
  std::vector<MPI_Request> sReq(nReq);

  for (int i = 0; i < nReq; i++) {
    for (int j = 0; j < lhs.cS[i].n; j++) { 
      int k = lhs.cS[i].ptr(j);
      sB(j,i) = R(k);
    }
  }

  int mpi_tag = 1;

  for (int i = 0; i < nReq; i++) {
    auto rec_err = MPI_Irecv(rB.col_data(i), lhs.cS[i].n, mpreal, lhs.cS[i].iP, mpi_tag, lhs.commu.comm, &rReq[i]);
    auto send_err = MPI_Isend(sB.col_data(i), lhs.cS[i].n, mpreal, lhs.cS[i].iP, mpi_tag, lhs.commu.comm, &sReq[i]);
  }

  // Wait for the MPI receive to complete.
  //
  for (int i = 0; i < nReq; i++) {
    MPI_Status stat;
    auto err = MPI_Wait(&rReq[i], &stat);
  }

  for (int i = 0; i < nReq; i++) {
    for (int j = 0; j < lhs.cS[i].n; j++) { 
      int k = lhs.cS[i].ptr(j);
      R(k) = R(k) + rB(j,i);
    }
  }

  // Wait for the MPI send to complete.
  //
  for (int i = 0; i < nReq; i++) {
    MPI_Status stat;
    auto err = MPI_Wait(&sReq[i], &stat);
  }
}

/// @brief This a both way communication with three main part:
///
/// 1 - rTmp {in master} = R          {from slave}
/// 2 - R    {in master} = R + rTmp   {both from master}
/// 3 - rTmp {in master} = R          {from master}
/// 4 - R    {in slave}  = rTmp       {from master}
//
void fsils_commuv(const FSILS_lhsType& lhs, int dof, Array<double>& R)
{
  if (lhs.commu.nTasks == 1) {
    return;
  }

  if (lhs.cS.size() == 0) {
    return;
  }

  int nReq = lhs.nReq;
  int nmax = std::max_element(lhs.cS.begin(), lhs.cS.end(), 
      [](const FSILS_cSType& a, const FSILS_cSType& b){return a.n < b.n;})->n; 

  Array3<double> sB(dof,nmax,nReq); 
  Array3<double> rB(dof,nmax,nReq); 
  std::vector<MPI_Request> rReq(nReq); 
  std::vector<MPI_Request> sReq(nReq);

  int n = 1;

  for (int i = 0; i < nReq; i++) {
    for (int j = 0; j < lhs.cS[i].n; j++) { 
      int k = lhs.cS[i].ptr(j);
      for (int l = 0; l < dof; l++) { 
        sB(l,j,i) = R(l,k);
      }
    }
  }

  int mpi_tag = 1;

  for (int i = 0; i < nReq; i++) {
    auto rec_err = MPI_Irecv(rB.slice_data(i), lhs.cS[i].n*dof, mpreal, lhs.cS[i].iP, mpi_tag, lhs.commu.comm, &rReq[i]);
    auto send_err = MPI_Isend(sB.slice_data(i), lhs.cS[i].n*dof, mpreal, lhs.cS[i].iP, mpi_tag, lhs.commu.comm, &sReq[i]);
  }

  // Wait for the MPI receive to complete.
  //
  for (int i = 0; i < nReq; i++) {
    MPI_Status stat;
    auto err = MPI_Wait(&rReq[i], &stat);
  }

  for (int i = 0; i < nReq; i++) {
    for (int j = 0; j < lhs.cS[i].n; j++) { 
      int k = lhs.cS[i].ptr(j);
      for (int l = 0; l < dof; l++) { 
        R(l,k) = R(l,k) + rB(l,j,i);
      }
    }
  }

  // Wait for the MPI send to complete.
  //
  for (int i = 0; i < nReq; i++) {
    MPI_Status stat;
    auto err = MPI_Wait(&sReq[i], &stat);
  }
}

};


