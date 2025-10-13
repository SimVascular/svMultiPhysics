// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// In this routine, the appropriate LS algorithm is called and
// the solution is returned.

#include "precond.h"

#include "fsils_api.hpp"

#include <math.h>

namespace precond {

/// @brief Post-multipling Val by W: Val = Val*W
///
/// Modifies: Val
//
void pos_mul(const Array<int>& rowPtr, const Vector<int>& colPtr, const int nNo, const int nnz, const int dof, Array<double>& Val, const Array<double>& W)
{
  switch (dof) {
    case 1: {
      for (int Ac = 0; Ac < nNo; Ac++) { 
        for (int i = rowPtr(0,Ac); i <= rowPtr(1,Ac); i++) {
          int a = colPtr(i);
          Val(0,i) = Val(0,i)*W(0,a);
        }
      }
    } break; 

    case 2: {
      for (int Ac = 0; Ac < nNo; Ac++) { 
        for (int i = rowPtr(0,Ac); i <= rowPtr(1,Ac); i++) {
          int a = colPtr(i);
          for (int j = 0; j < 3; j += 2) {
            Val(j+0,i) = Val(j+0,i)*W(0,a);
            Val(j+1,i) = Val(j+1,i)*W(1,a);
          }
        }
      }
    } break; 

    case 3: {
      for (int Ac = 0; Ac < nNo; Ac++) { 
        for (int i = rowPtr(0,Ac); i <= rowPtr(1,Ac); i++) {
          int a = colPtr(i);
          for (int j = 0; j < 7; j += 3) {
            Val(j+0,i) = Val(j+0,i)*W(0,a);
            Val(j+1,i) = Val(j+1,i)*W(1,a);
            Val(j+2,i) = Val(j+2,i)*W(2,a);
          }
        }
      }
    } break; 

    case 4: {
      for (int Ac = 0; Ac < nNo; Ac++) { 
        for (int i = rowPtr(0,Ac); i <= rowPtr(1,Ac); i++) {
          int a = colPtr(i);
          for (int j = 0; j < 13; j += 4) {
            Val(j+0,i) = Val(j+0,i)*W(0,a);
            Val(j+1,i) = Val(j+1,i)*W(1,a);
            Val(j+2,i) = Val(j+2,i)*W(2,a);
            Val(j+3,i) = Val(j+3,i)*W(3,a);
          }
        }
      }
    } break; 

    default: {
      for (int Ac = 0; Ac < nNo; Ac++) { 
        for (int i = rowPtr(0,Ac); i <= rowPtr(1,Ac); i++) {
          int a = colPtr(i);
          for (int b = 0; b < dof; b++) {
            int j = dof*(dof-1) + b;
            for (int k = b; k < j; k += dof) {
              Val(k,i) = Val(k,i)*W(b,a);
            }
          }
        }
      }
    } break; 
  }
}

//--------------
// precond_diag 
//--------------
// Jacobi symmetic preconditioner, to precondition both LHS and RHS.
//
// Modifies: Val, R, W
//
// Reproduces Fortran 'PRECONDDIAG'.
//
void precond_diag(fsi_linear_solver::FSILS_lhsType& lhs, const Array<int>& rowPtr, const Vector<int>& colPtr, 
    const Vector<int>& diagPtr, const int dof, Array<double>& Val, Array<double>& R, Array<double>& W)
{
  #define n_debug_precond_diag
  #ifdef debug_precond_diag
  DebugMsg dmsg(__func__,  lhs.commu.task);
  dmsg.banner();
  #endif

  int nNo = lhs.nNo;
  #ifdef debug_precond_diag
  dmsg << "lhs.nFaces: " << lhs.nFaces;
  dmsg << "nNo: " << nNo;
  dmsg << "dof: " << dof;
  dmsg << "Val.nrows: " << Val.nrows_;
  dmsg << "Val.ncols: " << Val.ncols_;
  dmsg << "W.nrows: " << W.nrows_;
  dmsg << "W.ncols: " << W.ncols_;
  #endif

  // Calculating W: W = diag(K)
  //
  switch (dof) {
    case 1: {
      for (int Ac = 0; Ac < nNo; Ac++) {
        W(0,Ac) = Val(0,diagPtr(Ac));
      }
    } break;

    case 2: {
      for (int Ac = 0; Ac < nNo; Ac++) {
        int d = diagPtr(Ac);
        W(0,Ac) = Val(0,d);
        W(1,Ac) = Val(3,d);
      }
    } break;

    case 3: {
      for (int Ac = 0; Ac < nNo; Ac++) {
        int d = diagPtr(Ac);
        W(0,Ac) = Val(0,d);
        W(1,Ac) = Val(4,d);
        W(2,Ac) = Val(8,d);
      }
    } break;

    case 4: {
      for (int Ac = 0; Ac < nNo; Ac++) {
        int d = diagPtr(Ac);
        W(0,Ac) = Val(0,d);
        W(1,Ac) = Val(5,d);
        W(2,Ac) = Val(10,d);
        W(3,Ac) = Val(15,d);
      }
    } break;

    default: {
      for (int Ac = 0; Ac < nNo; Ac++) {
        int d = diagPtr(Ac);
        for (int i = 0; i < dof; i++) {
          W(i,Ac) = Val(i*dof-dof+i,d);
        }
      }
    } break;
  }

  fsils_commuv(lhs, dof, W);

  // Accounting for Dirichlet BC and inversing W = W^{-1/2}
  //
  for (int Ac = 0; Ac < nNo; Ac++) {
    for (int i = 0; i < dof; i++) {
      if (W(i,Ac) == 0.0) {
        W(i,Ac) = 1.0;
      }
    }
  }

  for (int i = 0; i < W.size(); i++) {
    W(i) = 1.0 / sqrt(fabs(W(i)));
  }

  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    auto& face = lhs.face[faIn];
    #ifdef debug_precond_diag
    dmsg << ">>> faIn: " << faIn;
    dmsg << "face.incFlag: " << face.incFlag;
    #endif

    if (!face.incFlag) {
      continue;
    }

    int n = std::min(face.dof,dof);

    if (face.bGrp == fsi_linear_solver::BcType::BC_TYPE_Dir) {
      for (int a = 0; a < face.nNo; a++) {
        int Ac = face.glob(a);
        for (int i = 0; i < n; i++) {
          W(i,Ac) = W(i,Ac) * face.val(i,a);
        }
      }
    }
  }

  // Pre-multipling K with W: K = W*K
  pre_mul(rowPtr, lhs.nNo, lhs.nnz, dof, Val, W);

  // Multipling R with W: R = W*R
  //
  // W ( dof, lhs.nNo )
  //
  // R ( dof, lhs.nNo )
  //
  // ELement-wise multiplication.
  //
  for (int i = 0; i < W.size(); i++) {
    R(i) = W(i) * R(i);
  }

  // Now post-multipling K by W: K = K*W
  pos_mul(rowPtr, colPtr, lhs.nNo, lhs.nnz, dof, Val, W);

  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    auto& face = lhs.face[faIn];

    if (face.coupledFlag) {
      for (int a = 0; a < face.nNo; a++) {
        int Ac = face.glob(a);
        for (int i = 0; i < std::min(face.dof,dof); i++) {
          face.valM(i,a) = face.val(i,a) * W(i,Ac);
        }
      }
    }
  }
}

//-------------
// precond_rcs
//-------------
// Row and column preconditioner, to precondition both LHS and RHS.
//
// Reproduces Fortran 'PRECONDRCS'.
//
void precond_rcs(fsi_linear_solver::FSILS_lhsType& lhs, const Array<int>& rowPtr, const Vector<int>& colPtr,
    const Vector<int>& diagPtr, const int dof, Array<double>& Val, Array<double>& R, Array<double>& W1, Array<double>& W2)
{
  const int nNo = lhs.nNo;
  int maxiter = 10;
  double tol = 2.0;
  int iter = 0;
  bool flag = true;
  W1 = 1.0;
  W2 = 1.0;

  //*****************************************************
  // Apply Dirichlet BC
  //*****************************************************
  //
  Array<double> Wr(dof,nNo), Wc(dof,nNo);
  Wr = 1.0;

  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    auto& face = lhs.face[faIn];
    if (!face.incFlag) {
      continue;
    }

    int n = std::min(face.dof,dof);

    if (face.bGrp == fsi_linear_solver::BcType::BC_TYPE_Dir) {
      for (int a = 0; a < face.nNo; a++) {
        int Ac = face.glob(a);
        for (int i = 0; i < n; i++) {
          Wr(i,Ac) = Wr(i,Ac) * face.val(i,a);
        }
      }
    }
  }

  fsils_commuv(lhs, dof, Wr);

  // For parallel case, val and Wr can be larger than 1 due to
  // the addition operator in FSILS_COMMUV. Hence need renormalization.
  //
  Wr = Wr - 0.5;
  Wr = Wr / abs(Wr);
  Wr = (Wr + abs(Wr)) * 0.5;

  // Kill the row and column corresponding to Dirichlet BC
  //
  // Modifies 'Val'.
  //
  pre_mul(rowPtr, lhs.nNo, lhs.nnz, dof, Val, Wr);

  R = Wr * R;

  pos_mul(rowPtr, colPtr, lhs.nNo, lhs.nnz, dof, Val, Wr);

  // Set diagonal term to one
  //
  switch (dof) {
    case 1:
      for (int Ac = 0; Ac < nNo; Ac++) {
        int d = diagPtr(Ac);
        Val(0,d) = Wr(0,Ac) * (Val(0,d) - 1.0) + 1.0;
      }
    break; 

    case 2:
      for (int Ac = 0; Ac < nNo; Ac++) {
        int d = diagPtr(Ac);
        Val(0,d) = Wr(0,Ac)*(Val(0,d)-1.0) + 1.0;
        Val(3,d) = Wr(1,Ac)*(Val(3,d)-1.0) + 1.0;
      }
    break; 

    case 3:
      for (int Ac = 0; Ac < nNo; Ac++) {
        int d = diagPtr(Ac);
        Val(0,d) = Wr(0,Ac)*(Val(0,d)-1.0) + 1.0;
        Val(4,d) = Wr(1,Ac)*(Val(4,d)-1.0) + 1.0;
        Val(8,d) = Wr(2,Ac)*(Val(8,d)-1.0) + 1.0;
      }
    break; 

    case 4:
      for (int Ac = 0; Ac < nNo; Ac++) {
        int d = diagPtr(Ac);
        Val(0 ,d) = Wr(0,Ac)*(Val(0 ,d)-1.0) + 1.0;
        Val(5 ,d) = Wr(1,Ac)*(Val(5 ,d)-1.0) + 1.0;
        Val(10,d) = Wr(2,Ac)*(Val(10,d)-1.0) + 1.0;
        Val(15,d) = Wr(3,Ac)*(Val(15,d)-1.0) + 1.0;
      }
    break; 

    default: 
      for (int Ac = 0; Ac < nNo; Ac++) {
        int d = diagPtr(Ac);
        for (int i = 0; i < dof; i++) {
          Val(i*dof+i,d) = Wr(i,Ac)*(Val(i*dof+i,d) - 1.0) + 1.0;
        }
      }
    break; 
  } 

  //*****************************************************
  // Row and column scaling
  //*****************************************************
  //
  // Define a lambda function for computing the maximum
  // absolute value of a list of values.
  //
  auto max_func = [](const double& a, const double& b) { return fabs(a) < fabs(b); };

  while (flag) {
    Wr = 0.0;
    Wc = 0.0;
    iter = iter + 1;

    if (iter >= maxiter) {
      std::cout << "[precond_rcs] Warning: maximum iteration number reached";
      flag = false; 
    }

    // Max norm along row and column
    //
    switch (dof) {
      case 1:
        for (int Ac = 0; Ac < nNo; Ac++) {
          int a = rowPtr(0,Ac);
          int b = rowPtr(1,Ac);
          auto values = Val.values({0,0}, {a,b});
          Wr(0,Ac) = fabs(*std::max_element(values.begin(), values.end(), max_func));

          for (int i = rowPtr(0,Ac); i <= rowPtr(1,Ac); i++) {
            a = colPtr(i);
            Wc(0,a) = std::max(fabs(Val(0,i)), Wc(0,a));
          }
        }
      break; 

      case 2:
        for (int Ac = 0; Ac < nNo; Ac++) {
          int a = rowPtr(0,Ac);
          int b = rowPtr(1,Ac);

          auto vals1 = Val.values({0,1}, {a,b});
          Wr(0,Ac) = fabs(*std::max_element(vals1.begin(), vals1.end(), max_func));

          auto vals2 = Val.values({2,3}, {a,b});
          Wr(1,Ac) = fabs(*std::max_element(vals2.begin(), vals2.end(), max_func));

          for (int i = rowPtr(0,Ac); i <= rowPtr(1,Ac); i++) {
            a = colPtr(i);
            auto vals1 = Val.values({0,2}, {i,i}, 2);
            Wc(0,a) = std::max(fabs(*std::max_element(vals1.begin(), vals1.end(), max_func)), Wc(0,a));

            auto vals2 = Val.values({1,3}, {i,i}, 2);
            Wc(1,a) = std::max(fabs(*std::max_element(vals2.begin(), vals2.end(), max_func)), Wc(1,a));
          }
        }
      break; 

      case 3:
        for (int Ac = 0; Ac < nNo; Ac++) {
          int a = rowPtr(0,Ac);
          int b = rowPtr(1,Ac);
          auto vals1 = Val.values({0,2}, {a,b});
          Wr(0,Ac) = fabs(*std::max_element(vals1.begin(), vals1.end(), max_func));

          auto vals2 = Val.values({3,5}, {a,b});
          Wr(1,Ac) = fabs(*std::max_element(vals2.begin(), vals2.end(), max_func));

          auto vals3 = Val.values({6,8}, {a,b});
          Wr(2,Ac) = fabs(*std::max_element(vals3.begin(), vals3.end(), max_func));

          for (int i = rowPtr(0,Ac); i <= rowPtr(1,Ac); i++) {
            a = colPtr(i);
            auto vals1 = Val.values({0,6}, {i,i}, 3);
            Wc(0,a) = std::max(fabs(*std::max_element(vals1.begin(), vals1.end(), max_func)), Wc(0,a));

            auto vals2 = Val.values({1,7}, {i,i}, 3);
            Wc(1,a) = std::max(fabs(*std::max_element(vals2.begin(), vals2.end(), max_func)), Wc(1,a));

            auto vals3 = Val.values({2,8}, {i,i}, 3);
            Wc(2,a) = std::max(fabs(*std::max_element(vals3.begin(), vals3.end(), max_func)), Wc(2,a));
          }
        }
      break; 

      case 4:
        for (int Ac = 0; Ac < nNo; Ac++) {
          int a = rowPtr(0,Ac);
          int b = rowPtr(1,Ac);

          auto vals1 = Val.values({0,3}, {a,b});
          Wr(0,Ac) = fabs(*std::max_element(vals1.begin(), vals1.end(), max_func));

          auto vals2 = Val.values({4,7}, {a,b});
          Wr(1,Ac) = fabs(*std::max_element(vals2.begin(), vals2.end(), max_func));

          auto vals3 = Val.values({8,11}, {a,b});
          Wr(2,Ac) = fabs(*std::max_element(vals3.begin(), vals3.end(), max_func));

          auto vals4 = Val.values({12,15}, {a,b});
          Wr(3,Ac) = fabs(*std::max_element(vals4.begin(), vals4.end(), max_func));

          for (int i = rowPtr(0,Ac); i <= rowPtr(1,Ac); i++) {
            a = colPtr(i);
            auto vals1 = Val.values({0,12}, {i,i}, 4);
            Wc(0,a) = std::max(fabs(*std::max_element(vals1.begin(), vals1.end(), max_func)), Wc(0,a));

            auto vals2 = Val.values({1,13}, {i,i}, 4);
            Wc(1,a) = std::max(fabs(*std::max_element(vals2.begin(), vals2.end(), max_func)), Wc(1,a));

            auto vals3 = Val.values({2,14}, {i,i}, 4);
            Wc(2,a) = std::max(fabs(*std::max_element(vals3.begin(), vals3.end(), max_func)), Wc(2,a));

            auto vals4 = Val.values({3,15}, {i,i}, 4);
            Wc(3,a) = std::max(fabs(*std::max_element(vals4.begin(), vals4.end(), max_func)), Wc(3,a));
          }
        }
      break; 

      default: 
        for (int Ac = 0; Ac < nNo; Ac++) {
          int a = rowPtr(0,Ac);
          int b = rowPtr(1,Ac);

          for (int i = 0; i < dof; i++) {
            int j = i*dof + 1;
            int k = (i+1)*dof - 1; 
            auto vals = Val.values({j,k}, {a,b});
            Wr(i,Ac) = fabs(*std::max_element(vals.begin(), vals.end(), max_func));
          }

          for (int i = rowPtr(0,Ac); i <= rowPtr(1,Ac); i++) {
            a = colPtr(i);
            for (int b = 0; b < dof; b++) { 
              int j = dof*(dof-1) + b;
              auto vals = Val.values({b,j}, {i,i}, dof);
              Wc(b,a) = std::max(fabs(*std::max_element(vals.begin(), vals.end(), max_func)), Wc(b,a));
            }
          }
        }
      break; 
    } 

    fsils_commuv(lhs, dof, Wr);
    fsils_commuv(lhs, dof, Wc);

    if ((max(abs(1.0 - Wr)) < tol) && (max(abs(1.0 - Wc)) < tol)) {
      flag = false;
    }

    Wr = 1.0 / sqrt(Wr);
    Wc = 1.0 / sqrt(Wc);

    pre_mul(rowPtr, lhs.nNo, lhs.nnz, dof, Val, Wr);

    pos_mul(rowPtr, colPtr, lhs.nNo, lhs.nnz, dof, Val, Wc);

    W1 = W1 * Wr;
    W2 = W2 * Wc;

    if (lhs.commu.nTasks > 1) {
      int iflag = flag;
      std::vector<int> gflag(lhs.commu.nTasks);
      MPI_Allgather(&iflag, 1, cm_mod::mplog, gflag.data(), 1, cm_mod::mplog, lhs.commu.comm);
      flag = std::find(gflag.begin(), gflag.end(), 1) != gflag.end();
    }
  } // while

  // Multipling R with Wr: R = Wr*R
  R = W1 * R;
}

//---------
// pre_mul
//---------
// Pre-multipling Val with W: Val = W*Val.
//
// Modifies: Val(dof*dof, nnz)
//
// W(dof,nNo)
//
void pre_mul(const Array<int>& rowPtr, const int nNo, const int nnz, const int dof, Array<double>& Val, const Array<double>& W)
{
  switch (dof) {
    case 1: {
      for (int Ac = 0; Ac < nNo; Ac++) {
        int a = rowPtr(0,Ac);
        int b = rowPtr(1,Ac);
        for (int j = a; j <= b; j++) {
          Val(0,j) = Val(0,j)*W(0,Ac);
        }
      }
    } break;
    
    case 2: {
      for (int Ac = 0; Ac < nNo; Ac++) {
        int a = rowPtr(0,Ac);
        int b = rowPtr(1,Ac);
        for (int i = 0; i < 2; i++) {
          for (int j = a; j <= b; j++) {
            Val(i+0,j) = Val(i+0,j)*W(0,Ac);
            Val(i+2,j) = Val(i+2,j)*W(1,Ac);
          }
        }
      }
    } break;

    case 3: {
      for (int Ac = 0; Ac < nNo; Ac++) {
        int a = rowPtr(0,Ac);
        int b = rowPtr(1,Ac);
        for (int i = 0; i < 3; i++) {
          for (int j = a; j <= b; j++) {
            Val(i+0,j) = Val(i+0,j)*W(0,Ac);
            Val(i+3,j) = Val(i+3,j)*W(1,Ac);
            Val(i+6,j) = Val(i+6,j)*W(2,Ac);
          }
        }
      }
    } break;

    case 4: {
      for (int Ac = 0; Ac < nNo; Ac++) {
        int a = rowPtr(0,Ac);
        int b = rowPtr(1,Ac);

        for (int i = 0; i < 4; i++) {
          for (int j = a; j <= b; j++) {
            Val(i+0,j) = Val(i+0,j)*W(0,Ac);
            Val(i+4,j) = Val(i+4,j)*W(1,Ac);
            Val(i+8,j) = Val(i+8,j)*W(2,Ac);
            Val(i+12,j) = Val(i+12,j)*W(3,Ac);
          }
        }
      }
    } break;

    // Fill rows of 'Val' with length 'dof'.
    //
    default: {
      for (int Ac = 0; Ac < nNo; Ac++) {
        int a = rowPtr(0,Ac);
        int b = rowPtr(1,Ac);

        for (int i = 0; i < dof; i++) {
          int j = i*dof;

          for (int m = j; m < j+dof; m++) {
            for (int n = a; n <= b; n++) {
              Val(m,n) = Val(m,n) * W(i,Ac);
            }
          }
        }
      }
    } break;
  }
}

};
