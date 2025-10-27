// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// For calculating the norm of a scaler or vector based vector.
//
//
// Only the part of U which is owned by this processor is included in
// norm calculation, i.e. U(:,1:cS(tF)%ptr+cS(tF)%n-1)
// In order to have the correct answer it is needed that COMMU has
// been done before calling this function (or the ansesters of U
// are passed through COMMU)

#include "norm.h"

#include "CmMod.h"

#include "mpi.h"

#include <math.h>

namespace norm {

double fsi_ls_norms(const int nNo, FSILS_commuType& commu, const Vector<double>& U)
{
  double result = 0.0;

  for (int i = 0; i < nNo; i++) {
    result = result + U(i)*U(i);
  } 

  if (commu.nTasks != 1) {
    double tmp;
    MPI_Allreduce(&result, &tmp, 1, cm_mod::mpreal, MPI_SUM, commu.comm);
    result = tmp;
  }

  return sqrt(result);
}

double fsi_ls_normv(const int dof, const int nNo, FSILS_commuType& commu, const Array<double>& U)
{
  double result = 0.0;

  switch (dof) {
    case 1: {
      for (int i = 0; i < nNo; i++) {
        result = result + U(0,i)*U(0,i);
      }
    } break; 

    case 2: {
      for (int i = 0; i < nNo; i++) {
        result = result + U(0,i)*U(0,i) + U(1,i)*U(1,i);
      }
    } break; 

    case 3: {
      for (int i = 0; i < nNo; i++) {
        result = result + U(0,i)*U(0,i) + U(1,i)*U(1,i) + U(2,i)*U(2,i);
      }
    } break; 

    case 4: {
      for (int i = 0; i < nNo; i++) {
        result = result + U(0,i)*U(0,i) + U(1,i)*U(1,i) + U(2,i)*U(2,i) + U(3,i)*U(3,i);
      }
    } break; 

    default: { 
      for (int i = 0; i < nNo; i++) {
        for (int j = 0; j < U.nrows(); j++) {
          result = result + U(j,i)*U(j,i);
        }
      }
    } break; 
  } 

  if (commu.nTasks != 1) {
    double tmp;
    MPI_Allreduce(&result, &tmp, 1, cm_mod::mpreal, MPI_SUM, commu.comm);
    result = tmp;
  }

  return sqrt(result);
}

};
