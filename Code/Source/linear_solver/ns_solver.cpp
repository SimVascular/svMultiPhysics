// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// This routine is mainley intended for solving incompressible NS or
// FSI equations with a form of AU=R, in which A = [K D;-G L] and
// G = -D^t

#include "ns_solver.h"

#include "fsils_api.hpp"
#include "fils_struct.hpp"

#include "add_bc_mul.h"
#include "cgrad.h"
#include "dot.h"
#include "ge.h"
#include "gmres.h"
#include "norm.h"
#include "spar_mul.h"

#include "Array3.h"

#include <math.h>

namespace ns_solver {

/// @brief Modifies: lhs.face[].nS
//
void bc_pre(fsi_linear_solver::FSILS_lhsType& lhs, const int nsd, const int dof, const int nNo, const int mynNo) 
{
  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    auto& face = lhs.face[faIn];

    if (face.coupledFlag) {
      if (face.sharedFlag) {
        Array<double> v(nsd,nNo);

        for (int a = 0; a < face.nNo; a++) {
          int Ac = face.glob(a);
          for (int i = 0; i < nsd; i++) {
            v(i,Ac) = face.valM(i,a);
          }
        }

        face.nS = pow(norm::fsi_ls_normv(nsd, mynNo, lhs.commu, v),2.0);

      } else {
        face.nS = 0.0;
        for (int a = 0; a < face.nNo; a++) {
          int Ac = face.glob(a);
          for (int i = 0; i < nsd; i++) {
            face.nS = face.nS + pow(face.valM(i,a),2.0);
          }
        }
      }
    }
  }
}

/// @brief Store sections of the 'Val' into separate arrays: 'Gt', 'mK', etc.
///
/// Modifies: no globals
//
void depart(fsi_linear_solver::FSILS_lhsType& lhs, const int nsd, const int dof, const int nNo, const int nnz, 
    const Array<double>& Val, Array<double>& Gt, Array<double>& mK, Array<double>& mG, Array<double>& mD, Vector<double>& mL)
{
  Vector<double> tmp((nsd+1)*(nsd+1));

  if (nsd == 2) {
    for (int i = 0; i < nnz; i++) {
      auto tmp = Val.col(i);

      mK(0,i) = tmp(0);
      mK(1,i) = tmp(1);
      mK(2,i) = tmp(3);
      mK(3,i) = tmp(4);

      mG(0,i) = tmp(2);
      mG(1,i) = tmp(5);

      mD(0,i) = tmp(6);
      mD(1,i) = tmp(7);

      mL(i) = tmp(8);
    }

  } else if (nsd == 3) {

    for (int i = 0; i < nnz; i++) {
      auto tmp = Val.col(i);

      mK(0,i) = tmp(0);
      mK(1,i) = tmp(1);
      mK(2,i) = tmp(2);
      mK(3,i) = tmp(4);
      mK(4,i) = tmp(5);
      mK(5,i) = tmp(6);
      mK(6,i) = tmp(8);
      mK(7,i) = tmp(9);
      mK(8,i) = tmp(10);

      mG(0,i) = tmp(3);
      mG(1,i) = tmp(7);
      mG(2,i) = tmp(11);

      mD(0,i) = tmp(12);
      mD(1,i) = tmp(13);
      mD(2,i) = tmp(14);

      mL(i)   = tmp(15);
    }

  } else {
    //PRINT *, "FSILS: Not defined nsd for DEPART", nsd
    //STOP "FSILS: FATAL ERROR"
  }

  for (int i = 0; i < nNo; i++) {
    for (int j = lhs.rowPtr(0,i); j <= lhs.rowPtr(1,i); j++) {
      int k = lhs.colPtr(j);

      for (int l = lhs.rowPtr(0,k); l <= lhs.rowPtr(1,k); l++) {
        if (lhs.colPtr(l) == i) {
          for (int m = 0; m < Gt.nrows(); m++) {
            Gt(m,l) = -mG(m,j);
          }
          break;
        }
      }
    }
  }
}

/// @brief This routine is mainley intended for solving incompressible NS or
/// FSI equations with a form of AU=R, in which A = [K D;-G L] and
/// G = -D^t
///
/// Ri (dof, lhs.nNo )
//
void ns_solver(fsi_linear_solver::FSILS_lhsType& lhs, fsi_linear_solver::FSILS_lsType& ls, const int dof, const Array<double>& Val, Array<double>& Ri)
{
  using namespace consts;
  using namespace fsi_linear_solver;

  #define n_debug_ns_solver
  #ifdef debug_ns_solver
  DebugMsg dmsg(__func__,  lhs.commu.task);
  dmsg.banner();
  double time = fsi_linear_solver::fsils_cpu_t();
  #endif

  const int nNo = lhs.nNo;
  const int nnz = lhs.nnz;
  const int mynNo = lhs.mynNo;
  const int nsd = dof - 1;
  const int iB = ls.RI.mItr;
  const int nB = 2*iB;

  #ifdef debug_ns_solver
  dmsg << "dof: " << dof;
  dmsg << "nsd: " << nsd;
  dmsg << "nNo: " << nNo;
  dmsg << "nnz: " << nnz;
  dmsg << "mynNo: " << mynNo;
  dmsg << "iB: " << iB;
  dmsg << "nB: " << nB;
  Ri.write(msg_prefix+"Ri");
  #endif

  Vector<double> Rc(nNo), Rci(nNo), tmp(nB*nB+nB), tmpG(nB*nB+nB), B(nB), xB(nB), oldxB(nB); 
  Array<double> Rm(nsd,nNo), Rmi(nsd,nNo), A(nB,nB), P(nNo,iB), MP(nNo,nB);
  Array3<double> U(nsd,nNo,iB), MU(nsd,nNo,nB);

  for (int i = 0; i < nsd; i++) {
    for (int j = 0; j < Ri.ncols(); j++) {
      Rmi(i,j) = Ri(i,j);
    }
  }

  for (int i = 0; i < Ri.ncols(); i++) {
    Rci(i) = Ri(dof-1,i);
  }

  Rm = Rmi;
  Rc = Rci;

  double eps = sqrt(pow(norm::fsi_ls_normv(nsd,mynNo,lhs.commu,Rm),2.0) + 
                    pow(norm::fsi_ls_norms(mynNo,lhs.commu,Rc),2.0));

  #ifdef debug_ns_solver
  dmsg << "eps (Rm/Rc): " << eps;
  #endif

  ls.RI.iNorm = eps;
  ls.RI.fNorm = eps*eps;

  // Calling duration 
  ls.CG.callD = 0.0; 
  ls.GM.callD = 0.0;
  ls.RI.callD = fsi_linear_solver::fsils_cpu_t(); 

  ls.CG.itr   = 0;
  ls.GM.itr   = 0;
  ls.RI.suc   = false;
  eps = std::max(ls.RI.absTol, ls.RI.relTol*eps);
  #ifdef debug_ns_solver
  dmsg << "eps: " << eps;
  dmsg << "ls.RI.iNorm: " << ls.RI.iNorm;
  dmsg << "ls.RI.fNorm: " << ls.RI.fNorm;
  #endif

  Array<double> Gt(nsd,nnz), mK(nsd*nsd,nnz), mG(nsd,nnz), mD(nsd,nnz);
  Vector<double> mL(nnz); 

  // Store sections of the 'Val' array into separate arrays: 'Gt', 'mK', etc.
  //
  // Modfies: Gt, mK, mG, mD, and mL.
  //
  depart(lhs, nsd, dof, nNo, nnz, Val, Gt, mK, mG, mD, mL);

  // Computes lhs.face[].nS for each face.
  //
  bc_pre(lhs, nsd, dof, nNo, mynNo);

  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    auto& face = lhs.face[faIn];
    #ifdef debug_ns_solver
    dmsg << "faIn: " << faIn << "  face.nS: " << face.nS;
    #endif
   }

  #ifdef debug_ns_solver
  dmsg << "Loop i on ls.RI.mItr ... ";
  #endif
  int iBB{0}; 
  int i_count{0};
 
  // Note: iB and iBB appear to index into arrays.
  //
  for (int i = 0; i < ls.RI.mItr; i++) {
  //for (int i = 0; i < 1; i++) {
    #ifdef debug_ns_solver
    auto istr = "_" + std::to_string(i+1);
    dmsg << "---------- i " << i+1 << " ----------";
    #endif
   
    int iB = 2*i;
    iBB = 2*i + 1;
    ls.RI.dB = ls.RI.fNorm;
    i_count = i;
    #ifdef debug_ns_solver
    dmsg << "iB: " << iB;
    dmsg << "iBB: " << iBB;
    dmsg << "ls.RI.fNorm: " << ls.RI.fNorm;
    #endif

    // Solve for U = inv(mK) * Rm
    //
    auto U_slice = U.slice(i);
    gmres::gmres(lhs, ls.GM, nsd, mK, Rm, U_slice);
    U.set_slice(i, U_slice);

    // P = D*U
    //
    auto P_col = P.rcol(i);
    spar_mul::fsils_spar_mul_vs(lhs, lhs.rowPtr, lhs.colPtr, nsd, mD, U.rslice(i), P_col);
    //P.set_col(i, P_col);

    // P = Rc - P
    //
    P.set_col(i, Rc - P_col);

    // P = [L + G^t*G]^-1*P
    //
    P_col = P.rcol(i);
    cgrad::schur(lhs, ls.CG, nsd, Gt, mG, mL, P_col);
    //P.set_col(i, P_col);

    // MU1 = G*P
    //
    #ifdef debug_ns_solver
    dmsg << "i: " << i+1;
    dmsg << "iB: " << iB+1;
    #endif
    P_col = P.rcol(i);
    auto MU_iB = MU.rslice(iB);
    spar_mul::fsils_spar_mul_sv(lhs, lhs.rowPtr, lhs.colPtr, nsd, mG, P_col, MU_iB);
    //MU.set_slice(iB, MU_iB); 

    // MU2 = Rm - G*P
    //
    MU.set_slice(iBB, Rm - MU_iB);

    // U = inv(K) * [Rm - G*P]
    //
    lhs.debug_active = true;
    auto U_i = U.rslice(i);
    gmres::gmres(lhs, ls.GM, nsd, mK, MU.slice(iBB), U_i);
    //U.set_slice(i, U_i);

    // MU2 = K*U
    //
    auto MU_iBB = MU.rslice(iBB);
    spar_mul::fsils_spar_mul_vv(lhs, lhs.rowPtr, lhs.colPtr, nsd, mK, U.rslice(i), MU_iBB);
    //MU.set_slice(iBB, MU_iBB);

    add_bc_mul::add_bc_mul(lhs, BcopType::BCOP_TYPE_ADD, nsd, U.rslice(i), MU_iBB);
    //MU.set_slice(iBB, MU_iBB);

    // MP1 = L*P
    //
    auto MP_iB = MP.rcol(iB);
    spar_mul::fsils_spar_mul_ss(lhs, lhs.rowPtr, lhs.colPtr, mL, P.rcol(i), MP_iB);
    //MP.set_col(iB, MP_iB);

    // MP2 = D*U
    auto MP_iBB = MP.rcol(iBB);
    spar_mul::fsils_spar_mul_vs(lhs, lhs.rowPtr, lhs.colPtr, nsd, mD, U.rslice(i), MP_iBB);
    //MP.set_col(iBB, MP_iBB);

    int c = 0;

    for (int k = iB; k <= iBB; k++) {
      for (int j = 0; j <= k; j++) {
        tmp(c) = dot::fsils_nc_dot_v(nsd, mynNo, MU.slice(j), MU.slice(k))  +  
                 dot::fsils_nc_dot_s(mynNo, MP.col(j), MP.col(k));
        c = c + 1;
      }

      tmp(c) = dot::fsils_nc_dot_v(nsd, mynNo, MU.slice(k), Rmi)  +  
               dot::fsils_nc_dot_s(mynNo, MP.col(k), Rci);
      c = c + 1;
    }

    if (lhs.commu.nTasks > 1) {
      MPI_Allreduce(tmp.data(), tmpG.data(), c, cm_mod::mpreal, MPI_SUM, lhs.commu.comm);
      tmp = tmpG;
    }

    // Set arrays for Gauss elimination
    //
    c = 0;

    for (int k = iB; k <= iBB; k++) {
      for (int j = 0; j <= k; j++) {
        A(j,k) = tmp(c);
        A(k,j) = tmp(c);
        c = c + 1;
      }

      B(k) = tmp(c);
      c  = c + 1;
    }

    xB = B;

    // Perform Gauss elimination.
    //
    if (ge::ge(nB, iBB+1, A, xB)) {
      oldxB = xB;

    } else { 
      if (lhs.commu.masF) {
        throw std::runtime_error("FSILS: Singular matrix detected");
      }

      xB = oldxB;

      if (i > 0) {
        iB = iB - 2;
        iBB = iBB - 2;
      }
      break;
    }

    //dmsg << "Compute sum ... " ;
    double sum = 0.0;
    for (int i = 0; i <= iBB; i++) { 
      sum += xB(i) * B(i);
    }
    ls.RI.fNorm = pow(ls.RI.iNorm,2.0) - sum;
    #ifdef debug_ns_solver
    dmsg << "sum: " << sum;
    dmsg << "ls.RI.fNorm: " << ls.RI.fNorm;
    #endif

    if (ls.RI.fNorm < eps*eps) {
      ls.RI.suc = true;
      break;
    }

    Rm = Rmi - xB(0)*MU.slice(0);
    Rc = Rci - xB(0)*MP.col(0);

    for (int j = 1; j <= iBB; j++) {
      Rm = Rm - xB(j)*MU.slice(j);
      Rc = Rc - xB(j)*MP.col(j);
    }

  } // for i = 0; i < ls.RI.mItr

  if (i_count >= ls.RI.mItr) {
    ls.RI.itr = ls.RI.mItr;
  } else { 
    ls.RI.itr = i_count;
    Rc = Rci - xB(0)*MP.col(0);

    for (int j = 1; j <= iBB; j++) {
      Rc = Rc - xB(j)*MP.col(j);
    }
  }

  ls.Resc = static_cast<int>(100.0 * pow(norm::fsi_ls_norms(mynNo, lhs.commu, Rc),2.0) / ls.RI.fNorm);
  ls.Resm = 100 - ls.Resc;

  #ifdef debug_ns_solver
  dmsg << "ls.Resc: " << ls.Resc;
  dmsg << "ls.Resm: " << ls.Resm;
  dmsg << "ls.RI.itr: " << ls.RI.itr;
  #endif

  Rmi = xB(1) * U.slice(0);
  Rci = xB(0) * P.col(0);

  for (int i = 1; i <= ls.RI.itr; i++) {
    int iB = 2*i;
    int iBB = 2*i + 1;
    Rmi = Rmi + xB(iBB) * U.slice(i);
    Rci = Rci + xB(iB) * P.col(i);
  }

  // Set Calling duration.
  ls.RI.callD = fsi_linear_solver::fsils_cpu_t() - ls.RI.callD;

  ls.RI.dB = 5.0 * log(ls.RI.fNorm / ls.RI.dB);

  if (ls.Resc < 0.0 || ls.Resm < 0.0) {
    ls.Resc = 0;
    ls.Resm = 0;
    ls.RI.dB = 0;
    ls.RI.fNorm = 0.0;

    if (lhs.commu.masF) {
      throw std::runtime_error("FSILS: unexpected behavior in FSILS (likely due to the ill-conditioned LHS matrix)"); 
    }
  }

  ls.RI.fNorm = sqrt(ls.RI.fNorm);
  #ifdef debug_ns_solver
  dmsg << "ls.RI.callD: " << ls.RI.callD;
  dmsg << "ls.RI.dB: " << ls.RI.dB;
  dmsg << "ls.RI.fNorm: " << ls.RI.fNorm;
  #endif

  for (int i = 0; i < nsd; i++) {
    for (int j = 0; j < Rmi.ncols(); j++) {
      Ri(i, j) = Rmi(i, j);
    }
  }

  Ri.set_row(dof-1, Rci);

  if (lhs.commu.masF) {
    //CALL LOGFILE
  }

  #ifdef debug_ns_solver
  double exec_time = fsi_linear_solver::fsils_cpu_t() - time;
  dmsg << "Execution time: " << exec_time;
  dmsg << "Done";
  #endif
}

};

