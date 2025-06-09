/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

//-------------------------------------------------------------------------
// Graduate [sic] minimum residual algorithm is implemented here for vector
// and scaler problems.
//-------------------------------------------------------------------------

#include "gmres.h"

#include "fsils_api.hpp"

#include "add_bc_mul.h"
#include "bcast.h"
#include "dot.h"
#include "norm.h"
#include "omp_la.h"
#include "spar_mul.h"

#include "Array3.h"

#include <math.h>

namespace gmres {

void bc_pre(fsi_linear_solver::FSILS_lhsType& lhs, fsi_linear_solver::FSILS_subLsType& ls, const int dof, 
    const int mynNo, const int nNo)
{
  int nsd = dof - 1;
  Array<double> v(nsd,nNo);

  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    auto &face = lhs.face[faIn];
    if (face.coupledFlag) {
      if (face.sharedFlag) {
        v = 0.0;

        for (int a = 0; a < face.nNo; a++) {
          int Ac = face.glob(a);
          for (int i = 0; i < nsd; i++) {
            v(i,Ac) = face.valM(i,a);
          }
        }

        face.nS = pow(norm::fsi_ls_normv(nsd, mynNo, lhs.commu, v), 2.0);

      } else { 
        face.nS = 0.0;
        for (int a = 0; a < face.nNo; a++) {
          int Ac = face.glob(a);
          for (int i = 0; i < nsd; i++) {
            face.nS = face.nS + pow(face.valM(i,a), 2.0);
          }
        }
      }
    }
  }
}

/// @brief Solver the system Val * X = R.
///
/// Reproduces the Fortran 'GMRES' subroutine.
//
void gmres(fsi_linear_solver::FSILS_lhsType& lhs, fsi_linear_solver::FSILS_subLsType& ls, const int dof, 
    const Array<double>& Val, const Array<double>& R, Array<double>& X)
{
  #define n_debug_gmres
  #ifdef debug_gmres
  DebugMsg dmsg(__func__,  lhs.commu.task);
  dmsg.banner();
  #endif

  using namespace fsi_linear_solver;

  int nNo = lhs.nNo;
  int mynNo = lhs.mynNo;
  #ifdef debug_gmres
  dmsg << "dof: " << dof;
  dmsg << "nNo: " << nNo;
  dmsg << "mynNo: " << mynNo;
  dmsg << "ls.sD: " << ls.sD;
  dmsg << "ls.mItr: " << ls.mItr;
  dmsg << "ls.absTol: " << ls.absTol;
  dmsg << "ls.relTol: " << ls.relTol;
  #endif

  Array<double> h(ls.sD+1,ls.sD); 
  Array3<double> u(dof,nNo,ls.sD+1); 
  Array<double> unCondU(dof,nNo);
  Vector<double> y(ls.sD), c(ls.sD), s(ls.sD), err(ls.sD+1);

  double time = fsi_linear_solver::fsils_cpu_t(); 
  ls.suc = false;
  double eps = 0.0;
  int last_i = 0;
  X = 0.0;

  for (int l = 0; l < ls.mItr; l++) {
    #ifdef debug_gmres
    dmsg;
    dmsg << "---------- l " << l+1 << " ----------";
    #endif

    if (l == 0) {
      u.set_slice(0, R);
    } else {
      auto u_slice = u.rslice(0);
      spar_mul::fsils_spar_mul_vv(lhs, lhs.rowPtr, lhs.colPtr, dof,  Val, X, u_slice);

      add_bc_mul::add_bc_mul(lhs, BcopType::BCOP_TYPE_ADD, dof, X, u_slice);

      ls.itr = ls.itr + 1;
      u.set_slice(0, R - u_slice);
    }

    for (auto& face : lhs.face) {
      if (face.coupledFlag) {
        auto u_slice = u.rslice(0);
        auto unCondU = u.rslice(0);
        add_bc_mul::add_bc_mul(lhs, BcopType::BCOP_TYPE_PRE, dof, unCondU, u_slice);
        break; 
      }
    }

    err[0] = norm::fsi_ls_normv(dof, mynNo, lhs.commu, u.rslice(0));
    #ifdef debug_gmres
    dmsg << "err(1): " << err[0];
    #endif

    if (l == 0) {
      eps = err[0];

      if (eps <= ls.absTol) {
        ls.callD = std::numeric_limits<double>::epsilon();
        ls.dB = 0.0;
        return; 
      }

      ls.iNorm = eps;
      ls.fNorm = eps;
      eps = std::max(ls.absTol, ls.relTol*eps);
    }
    #ifdef debug_gmres
    dmsg << "eps: " << eps;
    #endif

    ls.dB = ls.fNorm;
    auto u_slice = u.rslice(0);
    u_slice = u_slice / err(0);

    for (int i = 0; i < ls.sD; i++) {
      #ifdef debug_gmres
      dmsg;
      dmsg << "----- i " << i+1 << " -----";
      #endif
      last_i = i;
      auto u_slice = u.rslice(i);
      auto u_slice_1 = u.rslice(i+1);
      spar_mul::fsils_spar_mul_vv(lhs, lhs.rowPtr, lhs.colPtr, dof,  Val, u_slice, u_slice_1);

      add_bc_mul::add_bc_mul(lhs, BcopType::BCOP_TYPE_ADD, dof, u_slice, u_slice_1);

      ls.itr = ls.itr + 1;

      for (auto& face : lhs.face) {
        if (face.coupledFlag) {
          auto u_slice_1 = u.rslice(i+1);
          auto unCondU = u.rslice(i+1);
          add_bc_mul::add_bc_mul(lhs, BcopType::BCOP_TYPE_PRE, dof, unCondU, u_slice_1);
          break;
        }
      }

      for (int j = 0; j <= i+1; j++) {
        h(j,i) = dot::fsils_nc_dot_v(dof, mynNo, u.rslice(j), u.rslice(i+1));
      }

      // h_col is modofied here so don't use 'rcol() method'.
      auto h_col = h.col(i);
      bcast::fsils_bcast_v(i+2, h_col, lhs.commu);
      h.set_col(i, h_col);

      for (int j = 0; j <= i; j++) {
        auto u_slice_1 = u.rslice(i+1);
        omp_la::omp_sum_v(dof, nNo, -h(j,i), u_slice_1, u.rslice(j));
        h(i+1,i) = h(i+1,i) - h(j,i)*h(j,i);
      }

      h(i+1,i) = sqrt(fabs(h(i+1,i)));

      u_slice_1 = u.rslice(i+1);
      omp_la::omp_mul_v(dof, nNo, 1.0/h(i+1,i), u_slice_1);

      for (int j = 0; j <= i-1; j++) {
        double tmp = c(j)*h(j,i) + s(j)*h(j+1,i);
        h(j+1,i) = -s(j)*h(j,i) + c(j)*h(j+1,i);
        h(j,i) = tmp;
      }

      double tmp = sqrt(h(i,i)*h(i,i) + h(i+1,i)*h(i+1,i));
      c(i) = h(i,i) / tmp;
      s(i) = h(i+1,i) / tmp;
      h(i,i) = tmp;
      h(i+1,i) = 0.0;
      err(i+1) = -s(i)*err(i);
      err(i) = c(i)*err(i);
      #ifdef debug_gmres
      dmsg;
      dmsg << "tmp: " << tmp;
      dmsg << "err(i): " << err(i);
      dmsg << "err(i+1): " << err(i+1);
      dmsg << "eps: " << eps;
      #endif

      if (fabs(err(i+1)) < eps) {
        ls.suc = true;
        break;
      }
    } // for int i = 0; i < ls.sD

    if (last_i >= ls.sD) {
      last_i = ls.sD - 1;
    }

    for (int i = 0; i <= last_i; i++) {
      y(i) = err(i);
    }

    for (int j = last_i; j >= 0; j--) { 
      for (int k = j+1; k <= last_i; k++) {
        y(j) = y(j) - h(j,k)*y(k);
      }
      y(j) = y(j) / h(j,j);
    }

    for (int j = 0; j <= last_i; j++) {
      omp_la::omp_sum_v(dof, nNo, y(j), X, u.rslice(j));
    }

    ls.fNorm = fabs(err(last_i+1));
    if (ls.suc) {
      break;
    }

  } // for l = 0; l < ls.mItr

  ls.callD = fsi_linear_solver::fsils_cpu_t() - time + ls.callD;
  ls.dB  = 10.0 * log(ls.fNorm / ls.dB);

  #ifdef debug_gmres
  dmsg << "Done";
  #endif
}

//---------
// gmres_s
//---------
// Reproduces the Fortran 'GMRESS' subroutine.
//
void gmres_s(fsi_linear_solver::FSILS_lhsType& lhs, fsi_linear_solver::FSILS_subLsType& ls, const int dof,
    const Vector<double>& Val, Vector<double>& R)
{
  #define n_debug_gmres_s
  #ifdef debug_gmres_s
  DebugMsg dmsg(__func__,  lhs.commu.task);
  dmsg.banner();
  #endif

  using namespace fsi_linear_solver;

  bool flag = false;
  int nNo = lhs.nNo;
  int mynNo = lhs.mynNo;
  #ifdef debug_gmres_s
  dmsg << "dof: " << dof;
  dmsg << "nNo: " << nNo;
  dmsg << "mynNo: " << mynNo;
  dmsg << "ls.sD: " << ls.sD;
  dmsg << "ls.mItr: " << ls.mItr;
  dmsg << "ls.absTol: " << ls.absTol;
  dmsg << "ls.relTol: " << ls.relTol;
  #endif

  Array<double> h(ls.sD+1,ls.sD);
  Array<double> u(nNo,ls.sD+1);
  Vector<double> X(nNo), y(ls.sD), c(ls.sD), s(ls.sD), err(ls.sD+1);

  ls.callD = fsi_linear_solver::fsils_cpu_t();
  ls.suc = false;
  double eps = norm::fsi_ls_norms(mynNo, lhs.commu, R);
  ls.iNorm = eps;
  ls.fNorm = eps;
  eps = std::max(ls.absTol, ls.relTol*eps);
  ls.itr = 0;
  int last_i = 0;
  #ifdef debug_gmres_s
  dmsg << "ls.iNorm: " << ls.iNorm;
  dmsg << "eps: " << eps;
  #endif

  if (ls.iNorm <= ls.absTol) {
    ls.callD = std::numeric_limits<double>::epsilon();
    ls.dB = 0.0;
    ls.suc = true;
    return; 
  }

  for (int l = 0; l < ls.mItr; l++) {
    #ifdef debug_gmres_s
    dmsg;
    dmsg << "======== l " << l+1 << " ======== ";
    #endif
    ls.dB = ls.fNorm;
    ls.itr = ls.itr + 1;
    auto u_col = u.col(0);
    spar_mul::fsils_spar_mul_ss(lhs, lhs.rowPtr, lhs.colPtr, Val, X, u_col);
    u.set_col(0, R - u_col);

    err[0] = norm::fsi_ls_norms(mynNo, lhs.commu, u.col(0));
    u_col = u.col(0) / err[0];
    u.set_col(0, u_col);
    #ifdef debug_gmres_s
    dmsg << "err(1): " << err[0];
    #endif


    for (int i = 0; i < ls.sD; i++) {
      #ifdef debug_gmres_s
      dmsg;
      dmsg << "----- i " << i+1 << " ----- ";
      #endif
      ls.itr = ls.itr + 1;
      last_i = i;
      auto u_col = u.col(i);
      auto u_col_1 = u.col(i+1);
      spar_mul::fsils_spar_mul_ss(lhs, lhs.rowPtr, lhs.colPtr, Val, u_col, u_col_1);
      u.set_col(i+1, u_col_1);

      for (int j = 0; j <= i+1; j++) {
        h(j,i) = dot::fsils_nc_dot_s(mynNo, u.col(j), u.col(i+1));
      }

      auto h_col = h.col(i);
      bcast::fsils_bcast_v(i+2, h_col, lhs.commu);
      h.set_col(i, h_col);

      for (int j = 0; j <= i; j++) {
        auto u_col_1 = u.col(i+1);
        omp_la::omp_sum_s(nNo, -h(j,i), u_col_1, u.col(j));
        u.set_col(i+1, u_col_1);
        h(i+1,i) = h(i+1,i) - h(j,i)*h(j,i);
      }
      h(i+1,i) = sqrt(fabs(h(i+1,i)));

      u_col_1 = u.col(i+1);
      omp_la::omp_mul_s(nNo, 1.0/h(i+1,i), u_col_1);
      u.set_col(i+1, u_col_1);

      for (int j = 0; j <= i-1; j++) {
        double tmp = c(j)*h(j,i) + s(j)*h(j+1,i);
        h(j+1,i) = -s(j)*h(j,i) + c(j)*h(j+1,i);
        h(j,i) = tmp;
      }

      double tmp = sqrt(h(i,i)*h(i,i) + h(i+1,i)*h(i+1,i));
      c(i) = h(i,i) / tmp;
      s(i) = h(i+1,i) / tmp;
      h(i,i) = tmp;
      h(i+1,i) = 0.0;
      err(i+1) = -s(i)*err(i);
      err(i) = c(i)*err(i);
      #ifdef debug_gmres_s
      dmsg << "err(i+1): " << err(i+1);
      dmsg << "tmp: " << tmp;
      #endif

      if (fabs(err(i+1)) < eps) {
        ls.suc = true;
        break;
      }
    } // for int i = 0; i < ls.sD

    if (last_i >= ls.sD) {
      last_i = ls.sD - 1;
    }

    for (int i = 0; i <= last_i; i++) {
      y(i) = err(i);
    }

    for (int j = last_i; j >= 0; j--) { 
      for (int k = j+1; k <= last_i; k++) {
        y(j) = y(j) - h(j,k)*y(k);
      }
      y(j) = y(j) / h(j,j);
    }

    for (int j = 0; j <= last_i; j++) {
      omp_la::omp_sum_s(nNo, y(j), X, u.col(j));
    }

    ls.fNorm = fabs(err(last_i+1));
    if (ls.suc) {
      break;
    }
  }

  R = X;
  ls.callD = fsi_linear_solver::fsils_cpu_t() - ls.callD;
  ls.dB  = 10.0 * log(ls.fNorm / ls.dB);
}

//---------
// gmres_v
//---------
// Generalized minimum residual algorithm implemented vector problems.
//
// The Array3::rslice() method is used to create an Array object with 
// data directly referenced to the Array3 data. This eliminates the overhead 
// of copying data to and from an Array3 object.
//
// Reproduces the Fortran 'GMRESV' subroutine.
//
void gmres_v(fsi_linear_solver::FSILS_lhsType& lhs, fsi_linear_solver::FSILS_subLsType& ls, const int dof,
    const Array<double>& Val, Array<double>& R)
{
  using namespace fsi_linear_solver;

  #define n_debug_gmres_v
  #ifdef debug_gmres_v
  DebugMsg dmsg(__func__,  lhs.commu.task);
  dmsg.banner();
  double time = fsi_linear_solver::fsils_cpu_t();
  #endif

  bool flag = false;
  int nNo = lhs.nNo;
  int mynNo = lhs.mynNo;
  #ifdef debug_gmres_v
  dmsg << "dof: " << dof;
  dmsg << "nNo: " << nNo;
  dmsg << "mynNo: " << mynNo;
  dmsg << "ls.sD: " << ls.sD;
  dmsg << "ls.mItr: " << ls.mItr;
  dmsg << "ls.absTol: " << ls.absTol;
  dmsg << "ls.relTol: " << ls.relTol;
  #endif

  Array<double> h(ls.sD+1,ls.sD), X(dof,nNo);
  Array3<double> u(dof,nNo,ls.sD+1);
  Array<double> unCondU(dof,nNo);
  Vector<double> y(ls.sD), c(ls.sD), s(ls.sD), err(ls.sD+1);

  ls.callD = fsi_linear_solver::fsils_cpu_t();
  ls.suc = false;
  double eps = norm::fsi_ls_normv(dof, mynNo, lhs.commu, R);
  ls.iNorm = eps;
  ls.fNorm = eps;
  eps = std::max(ls.absTol, ls.relTol*eps);
  ls.itr = 0;
  int last_i = 0;
  #ifdef debug_gmres_v
  dmsg << "ls.iNorm: " << ls.iNorm;
  dmsg << "eps: " << eps;
  #endif

  bc_pre(lhs, ls, dof, mynNo, nNo);

  if (ls.iNorm <= ls.absTol) {
    ls.callD = std::numeric_limits<double>::epsilon();
    ls.dB = 0.0;
    ls.suc = true;
    return; 
  }

  for (int l = 0; l < ls.mItr; l++) {
    #ifdef debug_gmres_v
    dmsg << "===== l " << l+1 << " ===== ";
    #endif
    ls.dB = ls.fNorm;
    ls.itr = ls.itr + 1;
    auto u_slice = u.rslice(0);
    spar_mul::fsils_spar_mul_vv(lhs, lhs.rowPtr, lhs.colPtr, dof,  Val, X, u_slice);

    add_bc_mul::add_bc_mul(lhs, BcopType::BCOP_TYPE_ADD, dof, X, u_slice);

    u_slice = R - u_slice;

    for (auto& face : lhs.face) {
      if (face.coupledFlag && flag) {
        auto u_slice = u.rslice(0);
        auto unCondU = u.rslice(0);
        add_bc_mul::add_bc_mul(lhs, BcopType::BCOP_TYPE_PRE, dof, unCondU, u_slice);
        break;
      }
    }

    err[0] = norm::fsi_ls_normv(dof, mynNo, lhs.commu, u.rslice(0));
    u_slice = u.rslice(0) / err[0];
    #ifdef debug_gmres_v
    dmsg << "err(1): " << err[0];
    #endif

    for (int i = 0; i < ls.sD; i++) {
      #ifdef debug_gmres_v
      dmsg << "----- i " << i+1 << " ----- ";
      #endif
      ls.itr = ls.itr + 1;
      last_i = i;
      auto u_slice = u.rslice(i);
      auto u_slice_1 = u.rslice(i+1);
      spar_mul::fsils_spar_mul_vv(lhs, lhs.rowPtr, lhs.colPtr, dof,  Val, u_slice, u_slice_1);

      add_bc_mul::add_bc_mul(lhs, BcopType::BCOP_TYPE_ADD, dof, u_slice, u_slice_1);

      for (auto& face : lhs.face) {
        if (face.coupledFlag && flag) {
          auto u_slice_1 = u.rslice(i+1);
          auto unCondU = u.rslice(i+1);
          add_bc_mul::add_bc_mul(lhs, BcopType::BCOP_TYPE_PRE, dof, unCondU, u_slice_1);
          break;
        }
      }

      for (int j = 0; j <= i+1; j++) {
        h(j,i) = dot::fsils_nc_dot_v(dof, mynNo, u.rslice(j), u.rslice(i+1));
        #ifdef debug_gmres_v
        dmsg << "h(j,i): " << h(j,i);
        #endif
      }

      auto h_col = h.col(i);
      bcast::fsils_bcast_v(i+2, h_col, lhs.commu);
      h.set_col(i, h_col);

      for (int j = 0; j <= i; j++) {
        auto u_slice_1 = u.rslice(i+1);
        omp_la::omp_sum_v(dof, nNo, -h(j,i), u_slice_1, u.rslice(j));
        h(i+1,i) = h(i+1,i) - h(j,i)*h(j,i);
      }
      h(i+1,i) = sqrt(fabs(h(i+1,i)));

      u_slice_1 = u.rslice(i+1);
      omp_la::omp_mul_v(dof, nNo, 1.0/h(i+1,i), u_slice_1);

      for (int j = 0; j <= i-1; j++) {
        double tmp = c(j)*h(j,i) + s(j)*h(j+1,i);
        h(j+1,i) = -s(j)*h(j,i) + c(j)*h(j+1,i);
        h(j,i) = tmp;
      }

      double tmp = sqrt(h(i,i)*h(i,i) + h(i+1,i)*h(i+1,i));
      c(i) = h(i,i) / tmp;
      s(i) = h(i+1,i) / tmp;
      h(i,i) = tmp;
      h(i+1,i) = 0.0;
      err(i+1) = -s(i)*err(i);
      err(i) = c(i)*err(i);
      #ifdef debug_gmres_v
      dmsg << "err(i+1): " << err(i+1);
      dmsg << "tmp: " << tmp;
      #endif

      if (fabs(err(i+1)) < eps) {
        ls.suc = true;
        break;
      }
    } // for int i = 0; i < ls.sD

    if (last_i >= ls.sD) {
      last_i = ls.sD - 1;
    }

    for (int i = 0; i <= last_i; i++) {
      y(i) = err(i);
    }

    for (int j = last_i; j >= 0; j--) { 
      for (int k = j+1; k <= last_i; k++) {
        y(j) = y(j) - h(j,k)*y(k);
      }
      y(j) = y(j) / h(j,j);
    }

    for (int j = 0; j <= last_i; j++) {
      omp_la::omp_sum_v(dof, nNo, y(j), X, u.rslice(j));
    }

    ls.fNorm = fabs(err(last_i+1));
    if (ls.suc) {
      break;
    }
  }

  R = X;
  ls.callD = fsi_linear_solver::fsils_cpu_t() - ls.callD;
  ls.dB  = 10.0 * log(ls.fNorm / ls.dB);

  #ifdef debug_gmres_v
  double exec_time = fsi_linear_solver::fsils_cpu_t() - time;
  dmsg << "Execution time: " << exec_time;
  dmsg << "Done";
  #endif
}

};


