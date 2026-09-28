// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file test_precond_rcs.cpp
 * @brief Matrix regression test for FSILS row and column scaling (RCS).
 *
 * @test RcsPreconditioner.PreservesConstrainedSystem verifies that
 * precond::precond_rcs preserves a Dirichlet-constrained system and its solution.
 *
 * @par Setup and coverage
 * A serial, three-node block tridiagonal matrix is stored in dense and FSILS
 * formats. The exact solution is zero at fixed components and `row + 1`
 * elsewhere; `b = A * exact` supplies the reference right-hand side.
 * The 75 cases combine `(dof, scale, boundary)`:
 * - Block sizes 1--5 exercise specialized and general scaling paths.
 * - Scales `1`, `1e-11`, `1e-15`, `1e-17`, and `1e-20` expose diagonal
 *   cancellation in the former update `(a - 1) + 1`.
 * - Boundary modes 0, 1, and 2 exclude the face, fix the first component of
 *   node 0, or fix all its components. Constrained rows retain only a unit
 *   diagonal, exposing omissions in row maxima; coupled blocks expose
 *   incomplete column scaling.
 *
 * @par Checks
 * Weights must be finite and positive. Removing the returned scaling must
 * recover the expected constrained operator and right-hand side. A dense solve
 * followed by right scaling must recover the exact solution and satisfy the
 * original free equations. Matrix, right-hand-side, and residual errors are
 * normalized by `scale` (unit diagonals use one), so erased small entries
 * cannot pass through an absolute tolerance.
 *
 * @par Reusing and running the test
 * Keep an independent dense reference and a known solution when adapting this
 * example. Update the FSILS block indices and diagonal pointers with the stencil,
 * and choose coefficient scales that expose the target failure.
 * Build with `-DENABLE_UNIT_TEST=ON`, then run beside the unit-test executable:
 * @code{.sh}
 * ./run_all_unit_tests --gtest_filter='Blocks/RcsPreconditioner.*'
 * @endcode
 */

#include "precond.h"
#include "FE/Math/DenseLinearAlgebra.h"

#include <gtest/gtest.h>

#include <cmath>
#include <tuple>
#include <vector>

namespace {

// Boundary modes: excluded face, first component fixed, all components fixed.
using RcsParameters = std::tuple<int, double, int>;
class RcsPreconditioner : public testing::TestWithParam<RcsParameters> {};

TEST_P(RcsPreconditioner, PreservesConstrainedSystem) {
  const auto [dof, scale, boundary] = GetParam();
  const int n = 3 * dof;
  const auto fixed = [=](int row) {
    return row < dof && boundary != 0 && (boundary == 2 || row == 0);
  };

  // Three-node tridiagonal stencil with coupled, positive-definite blocks.
  std::vector<double> original(n * n, 0.0), exact(n), original_rhs(n, 0.0);
  for (int row = 0; row < n; ++row) {
    exact[row] = fixed(row) ? 0.0 : row + 1.0;
    for (int col = 0; col < n; ++col) {
      const int a = row / dof;
      const int b = col / dof;
      if (std::abs(a - b) <= 1) {
        original[row * n + col] = scale * (a == b ? 2.0 : -1.0)
            * (row % dof == col % dof ? 4.0 : 0.25);
      }
    }
  }
  for (int row = 0; row < n; ++row) {
    for (int col = 0; col < n; ++col) {
      original_rhs[row] += original[row * n + col] * exact[col];
    }
  }

  fsi_linear_solver::FSILS_lhsType lhs{};
  lhs.nNo = lhs.gnNo = lhs.mynNo = 3;
  lhs.nnz = 7;
  lhs.commu.nTasks = 1;
  lhs.rowPtr.resize(2, 3);
  lhs.colPtr = Vector<int>{0, 1, 0, 1, 2, 1, 2};
  lhs.diagPtr = Vector<int>{0, 3, 6};
  const int starts[] = {0, 2, 5};
  const int ends[] = {1, 4, 6};
  for (int a = 0; a < 3; ++a) {
    lhs.rowPtr(0, a) = starts[a];
    lhs.rowPtr(1, a) = ends[a];
  }
  lhs.nFaces = 1;
  lhs.face.resize(1);
  auto& face = lhs.face[0];
  face.incFlag = boundary != 0;
  face.bGrp = fsi_linear_solver::BcType::BC_TYPE_Dir;
  face.nNo = 1;
  face.dof = dof;
  face.glob = Vector<int>{0};
  face.val.resize(dof, 1);
  for (int i = 0; i < dof; ++i) {
    face.val(i, 0) = (boundary == 2 || i == 0) ? 0.0 : 1.0;
  }

  Array<double> val(dof * dof, lhs.nnz), rhs(dof, 3);
  Array<double> left(dof, 3), right(dof, 3);
  for (int a = 0; a < 3; ++a) {
    for (int i = 0; i < dof; ++i) {
      rhs(i, a) = original_rhs[a * dof + i];
    }
    for (int k = starts[a]; k <= ends[a]; ++k) {
      for (int i = 0; i < dof; ++i) {
        for (int j = 0; j < dof; ++j) {
          val(i * dof + j, k) = original[(a * dof + i) * n + lhs.colPtr(k) * dof + j];
        }
      }
    }
  }

  precond::precond_rcs(lhs, lhs.rowPtr, lhs.colPtr, lhs.diagPtr,
                      dof, val, rhs, left, right);

  std::vector<double> matrix(n * n, 0.0), solution(n);
  for (int a = 0; a < 3; ++a) {
    for (int i = 0; i < dof; ++i) {
      ASSERT_TRUE(std::isfinite(left(i, a)));
      ASSERT_TRUE(std::isfinite(right(i, a)));
      ASSERT_GT(left(i, a), 0.0);
      ASSERT_GT(right(i, a), 0.0);
      const int row = a * dof + i;
      const double expected_rhs = fixed(row) ? 0.0 : original_rhs[row];
      ASSERT_NEAR(rhs(i, a) / left(i, a) / scale, expected_rhs / scale, 1e-12);
      solution[row] = rhs(i, a);
      for (int k = starts[a]; k <= ends[a]; ++k) {
        for (int j = 0; j < dof; ++j) {
          const int col = lhs.colPtr(k) * dof + j;
          const bool constrained = fixed(row) || fixed(col);
          const double expected = constrained ? (row == col ? 1.0 : 0.0)
                                               : original[row * n + col];
          const double restored = val(i * dof + j, k) / left(i, a) / right(j, lhs.colPtr(k));
          // Normalize by coefficient scale so erased small entries cannot pass.
          const double reference_scale = constrained && row == col ? 1.0 : scale;
          ASSERT_NEAR(restored / reference_scale, expected / reference_scale, 1e-12)
              << "row=" << row << ", col=" << col;
          matrix[row * n + col] = val(i * dof + j, k);
        }
      }
    }
  }

  svmp::FE::math::factor_dense_matrix(matrix, n, "RCS regression").solve_in_place(solution);
  for (int row = 0; row < n; ++row) {
    solution[row] *= right(row % dof, row / dof);
    EXPECT_NEAR(solution[row], exact[row], 1e-11);
  }
  for (int row = 0; row < n; ++row) {
    if (fixed(row)) continue;
    double residual = -original_rhs[row];
    for (int col = 0; col < n; ++col) {
      residual += original[row * n + col] * solution[col];
    }
    EXPECT_NEAR(residual / scale, 0.0, 1e-11);
  }
}

INSTANTIATE_TEST_SUITE_P(Blocks, RcsPreconditioner,
    testing::Combine(testing::Values(1, 2, 3, 4, 5), testing::Values(1.0, 1e-11, 1e-15, 1e-17, 1e-20),
                     testing::Values(0, 1, 2)));

} // namespace
