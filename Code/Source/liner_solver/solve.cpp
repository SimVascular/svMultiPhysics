// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "lhs.h"
#include "CmMod.h"
#include "bicgs.h"
#include "cgrad.h"
#include "gmres.h"
#include "ns_solver.h"
#include "precond.h"

namespace fsi_linear_solver {

/// @brief In this routine, the appropriate LS algorithm is called and
/// the solution is returned.
/// Modifies: Val, Ri
///
/// Ri(dof,lhs.nNo): Residual
/// Val(dof*dof,lhs.nnz): LHS
///
/// Reproduces 'SUBROUTINE FSILS_SOLVE (lhs, ls, dof, Ri, Val, prec, incL, res)'.
//
void fsils_solve(FSILS_lhsType& lhs, FSILS_lsType& ls, const int dof, Array<double>& Ri, Array<double>& Val, 
    const consts::PreconditionerType prec, const Vector<int>& incL, const Vector<double>& res)
{
  using namespace consts;

  #define n_debug_fsils_solve
  #ifdef debug_fsils_solve
  DebugMsg dmsg(__func__,  lhs.commu.task);
  dmsg.banner();
  #endif

  const int nNo = lhs.nNo;
  const int nnz = lhs.nnz;
  const int nFaces = lhs.nFaces;
  #ifdef debug_fsils_solve
  dmsg << "nNo: " << nNo;
  dmsg << "nnz: " << nnz;
  dmsg << "nFaces: " << nFaces;
  dmsg << "dof: " << dof;
  dmsg << "ls.LS_type: " << ls.LS_type;
  #endif

  if (lhs.nFaces != 0) {
    for (auto& face : lhs.face) { 
      face.incFlag = true;
    }

    if (incL.size() != 0) {
      for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
        #ifdef debug_fsils_solve
        dmsg << "incL[" << faIn << "]: " << incL[faIn];
        #endif
        if (incL(faIn) == 0) {
          lhs.face[faIn].incFlag = false;
        }
      }
    }

    bool flag = false;
    for (auto& face : lhs.face) { 
      if (face.bGrp == BcType::BC_TYPE_Neu) {
        flag = true;
        break;
      }
    }

    if (res.size() == 0 && flag) {
      throw std::runtime_error("[fsils_solve] res is required for Neu surfaces");
    }

    for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
      auto& face = lhs.face[faIn];
      face.coupledFlag = false;
      if (!face.incFlag) {
        continue;
      }
      bool flag = (face.bGrp == BcType::BC_TYPE_Neu);
      if (flag && res(faIn) != 0.0) {
        face.res = res(faIn);
        face.coupledFlag = true;
      }
    }
  }

  Array<double> R(dof,nNo), Wr(dof,nNo), Wc(dof,nNo);

  for (int a = 0; a < nNo; a++) {
    for (int i = 0; i < dof; i++) {
      R(i,lhs.map(a)) = Ri(i,a);
    }
  }

  // Apply preconditioner.
  //
  // Modifies Val and R.
  //

  if (prec == PreconditionerType::PREC_FSILS) {
    precond::precond_diag(lhs, lhs.rowPtr, lhs.colPtr, lhs.diagPtr, dof, Val, R, Wc);
  } else if (prec == PreconditionerType::PREC_RCS) {
    precond::precond_rcs(lhs, lhs.rowPtr, lhs.colPtr, lhs.diagPtr, dof, Val, R, Wr, Wc);
  } else {
    //PRINT *, "This linear solver and preconditioner combination is not supported."
  }

  // Solve for 'R'.
  //
  switch (ls.LS_type) {
    case LinearSolverType::LS_TYPE_NS:
      ns_solver::ns_solver(lhs, ls, dof, Val, R);
    break;

    case LinearSolverType::LS_TYPE_GMRES:
      if (dof == 1) {
        auto Valv = Val.row(0);
        auto Rv = R.row(0);
        gmres::gmres_s(lhs, ls.RI, dof, Valv, Rv);
        Val.set_row(0,Valv);
        R.set_row(0,Rv);
      } else {
        gmres::gmres_v(lhs, ls.RI, dof, Val, R);
      }
    break;

    case LinearSolverType::LS_TYPE_CG:
      if (dof == 1) {
        auto Valv = Val.row(0);
        auto Rv = R.row(0);
        cgrad::cgrad_s(lhs, ls.RI, Valv, Rv);
        Val.set_row(0,Valv);
        R.set_row(0,Rv);
      } else {
        cgrad::cgrad_v(lhs, ls.RI, dof, Val, R);
      }
    break;

    case LinearSolverType::LS_TYPE_BICGS:
      if (dof == 1) {
        auto Valv = Val.row(0);
        auto Rv = R.row(0);
        bicgs::bicgss(lhs, ls.RI, Valv, Rv);
        Val.set_row(0,Valv);
        R.set_row(0,Rv);
      } else {
        bicgs::bicgsv(lhs, ls.RI, dof, Val, R);
      }
    break;

    default:
      throw std::runtime_error("FSILS: LS_type not defined");
  }

  // Element-wise multiplication.
  //
  for (int i = 0; i < Wc.size(); i++) {
    R(i) = Wc(i) * R(i);
  }

  for (int a = 0; a < nNo; a++) {
    for (int i = 0; i < R.nrows(); i++) {
      Ri(i,a) = R(i,lhs.map(a));
    }
  }
}

};
