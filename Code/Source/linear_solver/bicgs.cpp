// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "bicgs.h"

#include "fsils_api.hpp"

#include "add_bc_mul.h"
#include "bcast.h"
#include "dot.h"
#include "norm.h"
#include "omp_la.h"
#include "spar_mul.h"

#include "Array3.h"

#include <math.h>

namespace bicgs {

/// @brief Biconjugate-gradient algorithm, available for scaler and vectors.
void bicgsv (fsi_linear_solver::FSILS_lhsType& lhs, fsi_linear_solver::FSILS_subLsType& ls, const int dof, 
    const Array<double>& K, Array<double>& R)
{
  #define n_debug_bicgsv
  #ifdef debug_bicgsv
  DebugMsg dmsg(__func__,  lhs.commu.task);
  dmsg.banner();
  #endif

  using namespace fsi_linear_solver;

  int nNo = lhs.nNo;
  int mynNo = lhs.mynNo;
  #ifdef debug_bicgsv
  dmsg << "ls.mItr: " << ls.mItr;
  dmsg << "dof: " << dof;
  dmsg << "nNo: " << nNo;
  dmsg << "mynNo: " << mynNo;
  #endif

  Array<double> P(dof,nNo), Rh(dof,nNo), X(dof,nNo), V(dof,nNo),   
      S(dof,nNo), T(dof,nNo);

  ls.callD = fsi_linear_solver::fsils_cpu_t();
  ls.suc = false;
  double err = norm::fsi_ls_normv(dof, mynNo, lhs.commu, R);
  double errO = err;
  ls.iNorm = err;
  double eps = std::max(ls.absTol,ls.relTol*err);
  double rho = err*err;
  double beta = rho;
  X = 0.0;
  P = R;
  Rh = R;
  int i_itr = 1;
  #ifdef debug_bicgsv
  dmsg;
  dmsg << "err: " << err;
  dmsg << "eps: " << eps;
  #endif

  for (int i = 0; i < ls.mItr; i++) {
    #ifdef debug_bicgsv
    dmsg;
    dmsg << "----- i " << i+1 << " -----";
    dmsg << "err: " << err;
    dmsg << "eps: " << eps;
    #endif
    if (err < eps) { 
      ls.suc = true;
      break;
    }

    spar_mul::fsils_spar_mul_vv(lhs, lhs.rowPtr, lhs.colPtr, dof, K, P, V);
    double alpha = rho / dot::fsils_dot_v(dof, mynNo, lhs.commu, Rh, V);
    S = R - alpha*V;

    spar_mul::fsils_spar_mul_vv(lhs, lhs.rowPtr, lhs.colPtr, dof, K, S, T);
    double omega = norm::fsi_ls_normv(dof, mynNo, lhs.commu, T);
    omega = dot::fsils_dot_v(dof, mynNo, lhs.commu, T, S) / (omega * omega);

    X = X + alpha*P + omega*S;
    R = S - omega*T;

    errO = err;
    err =  norm::fsi_ls_normv(dof, mynNo, lhs.commu, R);
    double rhoO  = rho;
    rho = dot::fsils_dot_v(dof, mynNo, lhs.commu, R, Rh);
    beta = rho*alpha / (rhoO*omega);

    #ifdef debug_bicgsv
    dmsg << "alpha: " << alpha;
    dmsg << "omega: " << omega;
    dmsg << "rho: " << rho;
    dmsg << "beta: " << beta;
    #endif

    P = R + beta * (P - omega*V);
    i_itr += 1;
  } 

  R = X;
  ls.itr = i_itr - 1;
  ls.fNorm = err;
  ls.callD =  fsi_linear_solver::fsils_cpu_t() - ls.callD;
  #ifdef debug_bicgsv
  dmsg << "ls.itr: " << ls.itr;
  #endif

  if (errO < std::numeric_limits<double>::epsilon()) { 
     ls.dB = 0.0;
  } else { 
     ls.dB = 10.0 * log(err / errO);
  }

}

//--------
// bicgss
//--------
//
void bicgss(fsi_linear_solver::FSILS_lhsType& lhs, fsi_linear_solver::FSILS_subLsType& ls, const Vector<double>& K, Vector<double>& R)
{
  #define n_debug_bicgss
  #ifdef debug_bicgss
  DebugMsg dmsg(__func__,  lhs.commu.task);
  dmsg.banner();
  #endif

  using namespace fsi_linear_solver;

  int nNo = lhs.nNo;
  int mynNo = lhs.mynNo;
  #ifdef debug_bicgss
  dmsg << "ls.mItr: " << ls.mItr;
  dmsg << "dof: " << dof;
  dmsg << "nNo: " << nNo;
  dmsg << "mynNo: " << mynNo;
  #endif

  Vector<double> P(nNo), Rh(nNo), X(nNo), V(nNo), S(nNo), T(nNo);

  ls.callD = fsi_linear_solver::fsils_cpu_t();
  ls.suc = false;
  double err = norm::fsi_ls_norms(mynNo, lhs.commu, R);
  double errO = err;
  ls.iNorm = err;
  double eps = std::max(ls.absTol,ls.relTol*err);
  double rho = err*err;
  double beta = rho;
  X = 0.0;
  P = R;
  Rh = R;
  int i_itr = 1;
  #ifdef debug_bicgss
  dmsg;
  dmsg << "err: " << err;
  dmsg << "eps: " << eps;
  #endif

  for (int i = 0; i < ls.mItr; i++) {
    #ifdef debug_bicgss
    dmsg;
    dmsg << "----- i " << i+1 << " -----";
    dmsg << "err: " << err;
    dmsg << "eps: " << eps;
    #endif
    if (err < eps) { 
      ls.suc = true;
      break;
    }

    spar_mul::fsils_spar_mul_ss(lhs, lhs.rowPtr, lhs.colPtr, K, P, V);
    double alpha = rho / dot::fsils_dot_s(mynNo, lhs.commu, Rh, V);
    S = R - alpha*V;

    spar_mul::fsils_spar_mul_ss(lhs, lhs.rowPtr, lhs.colPtr, K, S, T);
    double omega = norm::fsi_ls_norms(mynNo, lhs.commu, T);
    omega = dot::fsils_dot_s(mynNo, lhs.commu, T, S) / (omega * omega);

    X = X + alpha*P + omega*S;
    R = S - omega*T;

    errO = err;
    err =  norm::fsi_ls_norms(mynNo, lhs.commu, R);
    double rhoO  = rho;
    rho = dot::fsils_dot_s(mynNo, lhs.commu, R, Rh);
    beta = rho*alpha / (rhoO*omega);

    #ifdef debug_bicgss
    dmsg << "alpha: " << alpha;
    dmsg << "omega: " << omega;
    dmsg << "rho: " << rho;
    dmsg << "beta: " << beta;
    #endif

    P = R + beta * (P - omega*V);
    i_itr += 1;
  } 

  R = X;
  ls.itr = i_itr - 1;
  ls.fNorm = err;
  ls.callD = fsi_linear_solver::fsils_cpu_t() - ls.callD;
  #ifdef debug_bicgss
  dmsg << "ls.itr: " << ls.itr;
  #endif

  if (errO < std::numeric_limits<double>::epsilon()) { 
     ls.dB = 0.0;
  } else { 
     ls.dB = 10.0 * log(err / errO);
  }
}

};



