// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "add_bc_mul.h"

#include "bcast.h"

#include <algorithm>
#include <stdexcept>

namespace add_bc_mul {

namespace {

/// @brief Coefficient of the operator that a coupled Neumann face contributes.
///
/// @param[in] face The coupled face.
/// @param[in] op_Type Whether the contribution is to the matrix-vector product
///   or to the preconditioner.
///
/// @return face.res for the product, -face.res / (1 + face.res * face.nS) for
///   the preconditioner.
double coupling_coefficient(const FSILS_faceType &face,
                            const BcopType op_Type) {
  if (op_Type == BcopType::BCOP_TYPE_ADD) {
    return face.res;
  }

  if (op_Type == BcopType::BCOP_TYPE_PRE) {
    return -face.res / (1.0 + (face.res * face.nS));
  }

  throw std::runtime_error("[add_bc_mul] op_Type is not defined");
}

} // namespace

void add_bc_mul(const FSILS_lhsType &lhs, const BcopType op_Type,
                const Array<double> &X, Array<double> &Y) {
  // Scalar product of the scaled flow rate derivative of each coupled face
  // with X, over the entries that this process stores
  Vector<double> S(lhs.nFaces);
  bool has_coupled_face = false;

  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    const auto &face = lhs.face[faIn];
    const int num_components = std::min(face.dof, X.nrows());

    if (face.coupledFlag) {
      has_coupled_face = true;

      for (int a = 0; a < face.scaled_flowrate_derivative.ncols(); a++) {
        int Ac = face.flowrate_derivative_node_indices(a);
        for (int i = 0; i < num_components; i++) {
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
    const auto &face = lhs.face[faIn];
    const int num_components = std::min(face.dof, X.nrows());

    if (face.coupledFlag) {
      const double coef = coupling_coefficient(face, op_Type);

      for (int a = 0; a < face.nNo; a++) {
        int Ac = face.glob(a);
        for (int i = 0; i < num_components; i++) {
          Y(i, Ac) =
              Y(i, Ac) + face.scaled_unit_pressure_load(i, a) * coef * S(faIn);
        }
      }
    }
  }
}

void compute_coupled_face_norms(FSILS_lhsType &lhs, const int num_components) {
  // Sum of the squares of the scaled unit pressure load of each coupled face,
  // over the nodes that this process owns
  Vector<double> squared_norms(lhs.nFaces);
  bool has_coupled_face = false;

  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    const auto &face = lhs.face[faIn];

    if (face.coupledFlag) {
      has_coupled_face = true;

      for (int a = 0; a < face.nNo; a++) {
        // Skip the nodes that this process does not own.
        if (face.glob(a) >= lhs.mynNo) {
          continue;
        }

        for (int i = 0; i < num_components; i++) {
          const double value = face.scaled_unit_pressure_load(i, a);
          squared_norms(faIn) += value * value;
        }
      }
    }
  }

  // coupledFlag is the same on every process, so either all of them return
  // here or none does.
  if (!has_coupled_face) {
    return;
  }

  // Every node is owned by one process, so the sums are added up over the
  // processes.
  bcast::fsils_bcast_v(lhs.nFaces, squared_norms, lhs.commu);

  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    if (lhs.face[faIn].coupledFlag) {
      lhs.face[faIn].nS = squared_norms(faIn);
    }
  }
}
};
