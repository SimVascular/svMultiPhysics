// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "add_bc_mul.h"

#include "CmMod.h"
#include "dot.h"

namespace add_bc_mul {

/// @brief The contribution of coupled BCs is added to the matrix-vector
/// product operation. Depending on the type of operation (adding the
/// contribution or computing the PC contribution) different
/// coefficients are used.
///
/// For reference, see 
/// Moghadam et al. 2013 eq. 27 (https://doi.org/10.1016/j.jcp.2012.07.035) and
/// Moghadam et al. 2013b (https://doi.org/10.1007/s00466-013-0868-1).
///
/// Reproduces code in ADDBCMUL.f.
/// @param lhs The left-hand side of the linear system. 0D resistance is stored in the face(i).res field.
/// @param op_Type The type of operation (addition or PC contribution)
/// @param dof The number of degrees of freedom.
/// @param X The input vector.
/// @param Y The current matrix-vector product (Y = K*X), to which we add K^BC * X = res * v * v^T * X.
/// The expression is slightly different if preconditioning.
void add_bc_mul(FSILS_lhsType& lhs, const BcopType op_Type, const int dof, const Array<double>& X, Array<double>& Y)
{
  Vector<double> coef(lhs.nFaces); 
  Array<double> v(dof,lhs.nNo);

  if (op_Type == BcopType::BCOP_TYPE_ADD) {
    for (int i = 0; i < lhs.nFaces; i++) {
      coef(i) = lhs.face[i].res;
    }
  } else if (op_Type == BcopType::BCOP_TYPE_PRE) {
    for (int i = 0; i < lhs.nFaces; i++) {
      coef(i) = -lhs.face[i].res / (1.0 + (lhs.face[i].res*lhs.face[i].nS));
    }
  } else { 
    //PRINT *, "FSILS: op_Type is not defined"
    //STOP "FSILS: FATAL ERROR"
  }

  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    auto& face = lhs.face[faIn];
    int nsd = std::min(face.dof, dof);

    if (face.coupledFlag) {
      bool isCapped = (face.cap_valM.size() > 0 && face.cap_glob.size() > 0);

      // If face is shared across procs: use global dot for S so all ranks use the same S
      if (face.sharedFlag) {
        v = 0.0;
        // Setting vector v = int{N_A n_i} dGamma (main face)
        for (int a = 0; a < face.nNo; a++) {
          int Ac = face.glob(a);
          for (int i = 0; i < nsd; i++) {
            v(i,Ac) = face.valM(i,a);
          }
        }
        // Global dot product for main face contribution to S
        double S = dot::fsils_dot_v(dof, lhs.mynNo, lhs.commu, v, X);
        // Cap contribution: local sum then reduce so S is consistent across ranks
        if (isCapped) {
          double local_cap_dot = 0.0;
          for (int a = 0; a < face.cap_valM.ncols(); a++) {
            int Ac = face.cap_glob(a);
            for (int i = 0; i < nsd; i++) {
              local_cap_dot += face.cap_valM(i, a) * X(i, Ac);
            }
          }
          double global_cap_dot = 0.0;
          MPI_Allreduce(&local_cap_dot, &global_cap_dot, 1, cm_mod::mpreal, MPI_SUM, lhs.commu.comm);
          S += global_cap_dot;
        }
        S = coef(faIn) * S;

        // Computing Y = Y + v * S
        for (int a = 0; a < face.nNo; a++) {
          int Ac = face.glob(a);
          for (int i = 0; i < nsd; i++) {
            Y(i,Ac) = Y(i,Ac) + v(i,Ac)*S;
          }
        }
      } else {
        // Face not shared: local S (and cap contribution) then local Y update
        double S = 0.0;
        for (int a = 0; a < face.nNo; a++) {
          int Ac = face.glob(a);
          for (int i = 0; i < nsd; i++) {
            S = S + face.valM(i,a)*X(i,Ac);
          }
        }
        if (isCapped) {
          for (int a = 0; a < face.cap_valM.ncols(); a++) {
            int Ac = face.cap_glob(a);
            for (int i = 0; i < nsd; i++) {
              S = S + face.cap_valM(i, a) * X(i, Ac);
            }
          }
        }
        S = coef(faIn) * S;

        for (int a = 0; a < face.nNo; a++) {
          int Ac = face.glob(a);
          for (int i = 0; i < nsd; i++) {
            Y(i,Ac) = Y(i,Ac) + face.valM(i,a)*S;
          }
        }
      }
    }
  }

}

};
