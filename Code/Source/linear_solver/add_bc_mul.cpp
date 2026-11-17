// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "add_bc_mul.h"

#include "bcast.h"

#include <stdexcept>

namespace add_bc_mul {

void add_bc_mul(FSILS_lhsType& lhs, const BcopType op_Type, const int dof, const Array<double>& X, Array<double>& Y)
{
  Vector<double> coef(lhs.nFaces);

  if (op_Type == BcopType::BCOP_TYPE_ADD) {
    for (int i = 0; i < lhs.nFaces; i++) {
      coef(i) = lhs.face[i].res;
    }
  } else if (op_Type == BcopType::BCOP_TYPE_PRE) {
    for (int i = 0; i < lhs.nFaces; i++) {
      coef(i) = -lhs.face[i].res / (1.0 + (lhs.face[i].res*lhs.face[i].nS));
    }
  } else { 
    throw std::runtime_error("[add_bc_mul] op_Type is not defined");
  }

  // Scalar product of the scaled flow rate derivative of each coupled face
  // with X, over the entries that this process stores
  Vector<double> S(lhs.nFaces);
  bool has_coupled_face = false;

  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    auto& face = lhs.face[faIn];
    int nsd = std::min(face.dof, dof);

    if (face.coupledFlag) {
      has_coupled_face = true;

      for (int a = 0; a < face.scaled_flowrate_derivative.ncols(); a++) {
        int Ac = face.flowrate_derivative_node_indices(a);
        for (int i = 0; i < nsd; i++) {
          S(faIn) = S(faIn) + face.scaled_flowrate_derivative(i, a) * X(i, Ac);
        }
      }
    }
  }

  // coupledFlag is the same on every process, so either all of them return
  // here or none does.
  if (!has_coupled_face) {
    return;
  }

  // A process stores the entries of the nodes it owns, so the scalar products
  // are summed over the processes.
  bcast::fsils_bcast_v(lhs.nFaces, S, lhs.commu);

  // Computing Y = Y + coef * scaled unit pressure load * S
  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    auto &face = lhs.face[faIn];
    int nsd = std::min(face.dof, dof);

    if (face.coupledFlag) {
      for (int a = 0; a < face.nNo; a++) {
        int Ac = face.glob(a);
        for (int i = 0; i < nsd; i++) {
          Y(i, Ac) = Y(i, Ac) + face.scaled_unit_pressure_load(i, a) *
                                    coef(faIn) * S(faIn);
        }
      }
    }
  }
}

};
