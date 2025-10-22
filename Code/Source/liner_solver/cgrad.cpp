// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "cgrad.h"
#include "DebugMsg.h"

#include "fsils_api.hpp"
#include "add_bc_mul.h"
#include "dot.h"
#include "omp_la.h"
#include "norm.h"
#include "spar_mul.h"

#include <math.h>

namespace cgrad {

/// @brief Conjugate-gradient algorithm for scaler, vector and Schur
/// complement cases.
///
/// Reproduces 'SUBROUTINE CGRAD_SCHUR(lhs, ls, dof, D, G, L, R)'
//
void schur(FSILS_lhsType& lhs, FSILS_subLsType& ls, const int dof, const Array<double>& D, 
    const Array<double>& G, const Vector<double>& L, Vector<double>& R)
{
  #define n_debug_schur
  #ifdef debug_schur
  DebugMsg dmsg(__func__,  lhs.commu.task);
  dmsg.banner();
  #endif

  int nNo = lhs.nNo;
  int mynNo = lhs.mynNo;

  Vector<double> X(nNo), P(nNo), SP(nNo), DGP(nNo); 
  Array<double> GP(dof,nNo), unCondU(dof,nNo);

  double time = fsi_linear_solver::fsils_cpu_t();
  ls.suc = false;
  ls.iNorm = norm::fsi_ls_norms(mynNo, lhs.commu, R);
  double eps = pow(std::max(ls.absTol,ls.relTol*ls.iNorm),2.0);
  double errO = ls.iNorm*ls.iNorm;
  double err = errO;
  #ifdef debug_schur
  dmsg << "dof: " << dof;
  dmsg << "nNo: " << nNo;
  dmsg << "mynNo: " << mynNo;
  dmsg << "ls.iNorm: " << ls.iNorm;
  dmsg << "eps: " << eps;
  dmsg << "errO: " << errO;
  #endif

  X = 0.0;
  P = R;
  int last_i = 0;

  for (int i = 0; i < ls.mItr; i++) {
    #ifdef debug_schur
    dmsg;
    dmsg << "----- i " << i+1 << " -----";
    dmsg << "err: " << err;
    auto istr = "_" + std::to_string(i+1);
    #endif
    last_i = i;

    if (err < eps) {
      ls.suc = true;
      break;
    }

    errO = err;

    // GP = G * P
    spar_mul::fsils_spar_mul_sv(lhs, lhs.rowPtr, lhs.colPtr, dof, G, P, GP);

    for (auto& face : lhs.face) {
      if (face.coupledFlag) {
        auto unCondU = GP;
        add_bc_mul::add_bc_mul(lhs, BcopType::BCOP_TYPE_PRE, dof, unCondU, GP);
        break;
      }
    }

    // DGP = K * GP
    spar_mul::fsils_spar_mul_vs(lhs, lhs.rowPtr, lhs.colPtr, dof, D, GP, DGP);

    // SP = L * P
    spar_mul::fsils_spar_mul_ss(lhs, lhs.rowPtr, lhs.colPtr, L, P, SP);

    // SP = SP - DGP
    omp_la::omp_sum_s(nNo, -1.0, SP, DGP);

    double alpha = errO / dot::fsils_dot_s(mynNo, lhs.commu, P, SP);

    // X = X + alpha * P
    omp_la::omp_sum_s(nNo, alpha, X, P);

    // R = R - alpha * SP
    omp_la::omp_sum_s(nNo, -alpha, R, SP);

    err = norm::fsi_ls_norms(mynNo, lhs.commu, R);
    err = err * err;
    #ifdef debug_schur
    dmsg << "err: " << err;
    dmsg << "errO/err: " << errO/err;
    dmsg << "err/errO: " << err/errO;
    #endif

    // P = P + errO/err * R 
    double c1 = errO / err;
    omp_la::omp_sum_s(nNo, c1, P, R);

    // P = err/errO * P
    double c2 = err / errO;
    omp_la::omp_mul_s(nNo, c2, P);
  }

  R = X;
  ls.fNorm = sqrt(err);
  ls.callD = fsi_linear_solver::fsils_cpu_t() - time + ls.callD;
  ls.itr = ls.itr + last_i;
  #ifdef debug_schur
  dmsg << "errO: " << errO;
  dmsg << "ls.fNorm: " << ls.fNorm;
  dmsg << "ls.itr: " << ls.itr;
  #endif

  if (errO < std::numeric_limits<double>::epsilon()) {
    ls.dB = 0.0;
  } else {
    ls.dB = 5.0 * log(err/errO);
  }
}

//---------
// cgrad_v
//---------
//
void cgrad_v(FSILS_lhsType& lhs, FSILS_subLsType& ls, const int dof, const Array<double>& K, Array<double>& R)
{
  #define n_debug_cgrad_v 
  #ifdef debug_cgrad_v
  DebugMsg dmsg(__func__,  lhs.commu.task);
  dmsg.banner();
  double time = fsi_linear_solver::fsils_cpu_t();
  #endif

  int nNo = lhs.nNo;
  int mynNo = lhs.mynNo;
  #ifdef debug_cgrad_v
  dmsg << "nNo: " << nNo;
  dmsg << "mynNo: " << mynNo;
  dmsg << "ls.mItr: " << ls.mItr;
  #endif

  Array<double> P(dof,nNo), KP(dof,nNo), X(dof,nNo);

  ls.callD = fsi_linear_solver::fsils_cpu_t();
  ls.suc = false;
  ls.iNorm = norm::fsi_ls_normv(dof, mynNo, lhs.commu, R);
  double eps = pow(std::max(ls.absTol, ls.relTol* ls.iNorm), 2.0);

  double errO = ls.iNorm * ls.iNorm;
  double err  = errO;
  X = 0.0;
  P = R;
  int last_i = 0;
  #ifdef debug_cgrad_v
  dmsg << "ls.iNorm: " << ls.iNorm;
  dmsg << "eps: " << eps;
  dmsg << "err: " << eps;
  #endif

  for (int i = 0; i < ls.mItr; i++) {
    #ifdef debug_cgrad_v
    dmsg;
    dmsg << "----- i " << i+1 << " -----";
    dmsg << "err: " << err;
    auto istr = "_" + std::to_string(i+1);
    #endif
    last_i = i;

    if (err < eps) {
      ls.suc = true;
      break;
    }

    errO = err;

    spar_mul::fsils_spar_mul_vv(lhs, lhs.rowPtr, lhs.colPtr, dof, K, P, KP);

    double alpha = errO / dot::fsils_dot_v(dof, mynNo, lhs.commu, P, KP);
    omp_la::omp_sum_v(dof, nNo, alpha, X, P);
    omp_la::omp_sum_v(dof, nNo, -alpha, R, KP);

    err = norm::fsi_ls_normv(dof, mynNo, lhs.commu, R);
    err = err * err;

    omp_la::omp_sum_v(dof, nNo, errO/err, P, R);
    omp_la::omp_mul_v(dof, nNo, err/errO, P);
  }

  R = X;
  ls.itr = last_i;
  ls.fNorm = sqrt(err);
  ls.callD = fsi_linear_solver::fsils_cpu_t() - ls.callD;

  if (errO < std::numeric_limits<double>::epsilon()) {
    ls.dB = 0.0;
  } else {
    ls.dB = 5.0 * log(err/errO);
  }

  #ifdef debug_cgrad_v
  double exec_time = fsi_linear_solver::fsils_cpu_t() - time;
  dmsg << "Execution time: " << exec_time;
  dmsg << "Done";
  #endif
}

//---------
// cgrad_s
//---------
//
void cgrad_s(FSILS_lhsType& lhs, FSILS_subLsType& ls, const Vector<double>& K, Vector<double>& R)
{
  #define n_debug_cgrad_s 
  #ifdef debug_cgrad_s
  DebugMsg dmsg(__func__,  lhs.commu.task);
  dmsg.banner();
  double time = fsi_linear_solver::fsils_cpu_t();
  #endif

  int nNo = lhs.nNo;
  int mynNo = lhs.mynNo;
  #ifdef debug_cgrad_s
  dmsg << "nNo: " << nNo;
  dmsg << "mynNo: " << mynNo;
  dmsg << "ls.mItr: " << ls.mItr;
  #endif

  Vector<double> P(nNo), KP(nNo), X(nNo);

  ls.callD = fsi_linear_solver::fsils_cpu_t();
  ls.suc = false;
  ls.iNorm = norm::fsi_ls_norms(mynNo, lhs.commu, R);
  double eps = pow(std::max(ls.absTol, ls.relTol* ls.iNorm), 2.0);
  double errO = ls.iNorm * ls.iNorm;
  double err  = errO;
  X = 0.0;
  P = R;
  int last_i = 0;
  #ifdef debug_cgrad_s
  dmsg << "ls.iNorm: " << ls.iNorm;
  dmsg << "eps: " << eps;
  dmsg << "err: " << eps;
  #endif

  for (int i = 0; i < ls.mItr; i++) {
    #ifdef debug_cgrad_s
    dmsg;
    dmsg << "----- i " << i+1 << " -----";
    dmsg << "err: " << err;
    auto istr = "_" + std::to_string(i+1);
    #endif
    last_i = i;

    if (err < eps) {
      ls.suc = true;
      break;
    }

    errO = err;

    spar_mul::fsils_spar_mul_ss(lhs, lhs.rowPtr, lhs.colPtr, K, P, KP);

    double alpha = errO / dot::fsils_dot_s(mynNo, lhs.commu, P, KP);
    omp_la::omp_sum_s(nNo, alpha, X, P);
    omp_la::omp_sum_s(nNo, -alpha, R, KP);

    err = norm::fsi_ls_norms(mynNo, lhs.commu, R);
    err = err * err;

    omp_la::omp_sum_s(nNo, errO/err, P, R);
    omp_la::omp_mul_s(nNo, err/errO, P);
  }

  R = X;
  ls.itr = last_i;
  ls.fNorm = sqrt(err);
  ls.callD = fsi_linear_solver::fsils_cpu_t() - ls.callD;

  if (errO < std::numeric_limits<double>::epsilon()) {
    ls.dB = 0.0;
  } else {
    ls.dB = 5.0 * log(err/errO);
  }

  #ifdef debug_cgrad_s
  double exec_time = fsi_linear_solver::fsils_cpu_t() - time;
  dmsg << "Execution time: " << exec_time;
  dmsg << "Done";
  #endif
}

};


