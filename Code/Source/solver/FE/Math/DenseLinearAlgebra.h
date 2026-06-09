// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_MATH_DENSELINEARALGEBRA_H
#define SVMP_FE_MATH_DENSELINEARALGEBRA_H

#include "Types.h"

#include <Eigen/Dense>

#include <cstddef>
#include <limits>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace svmp {
namespace FE {
namespace math {

// Dense solve, inverse, rank, and pseudo-inverse support for FE construction
// utilities, backed by Eigen. Matrices are row-major: matrix[row * cols + col].
[[nodiscard]] Real dense_matrix_max_abs(std::span<const Real> matrix) noexcept;

[[nodiscard]] Real dense_matrix_pivot_tolerance(std::size_t rows,
                                                std::size_t cols,
                                                Real max_abs,
                                                Real multiplier = Real(64)) noexcept;

[[nodiscard]] Real dense_matrix_singular_value_tolerance(std::size_t rows,
                                                         std::size_t cols,
                                                         Real largest_singular_value,
                                                         Real multiplier = Real(64)) noexcept;

struct DensePseudoInverseResult {
    std::vector<Real> inverse;
    std::size_t rank{0};
    Real tolerance{0};
    Real largest_singular_value{0};
    Real smallest_retained_singular_value{0};
};

struct DenseMatrixDiagnostics {
    std::size_t rank{0};
    Real tolerance{0};
    Real largest_singular_value{0};
    Real smallest_retained_singular_value{0};
    Real condition_estimate{std::numeric_limits<Real>::infinity()};
};

struct DenseInverseResult {
    std::vector<Real> inverse;
    DenseMatrixDiagnostics diagnostics;
    bool used_svd_fallback{false};
};

[[nodiscard]] Real dense_matrix_condition_fallback_threshold() noexcept;
[[nodiscard]] Real dense_matrix_condition_error_threshold() noexcept;

struct DenseLUSolver {
    using DenseMatrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

    std::size_t n{0};
    Eigen::PartialPivLU<DenseMatrix> lu;
    DenseMatrixDiagnostics diagnostics;
    Real pivot_tolerance{0};
    std::string label;

    [[nodiscard]] bool empty() const noexcept { return n == 0; }

    void solve_in_place(std::span<Real> rhs) const;
    void solve_in_place(std::span<Real> rhs, std::size_t rhs_count) const;
    [[nodiscard]] std::vector<Real> solve(std::span<const Real> rhs) const;
};

// Inverses and pseudo-inverses keep the same row-major convention for their
// returned dimensions.
[[nodiscard]] DenseMatrixDiagnostics dense_matrix_diagnostics(
    std::span<const Real> matrix,
    std::size_t rows,
    std::size_t cols,
    std::string_view label = "dense matrix");

[[nodiscard]] DenseLUSolver factor_dense_matrix(std::vector<Real> matrix,
                                                std::size_t n,
                                                std::string_view label = "dense matrix");

[[nodiscard]] std::vector<Real> invert_dense_matrix(std::vector<Real> matrix,
                                                    std::size_t n,
                                                    std::string_view label = "dense matrix");

[[nodiscard]] DenseInverseResult invert_dense_matrix_with_diagnostics(
    std::vector<Real> matrix,
    std::size_t n,
    std::string_view label = "dense matrix");

void validate_dense_inverse_diagnostics(
    const DenseInverseResult& result,
    std::size_t expected_rank,
    std::string_view label = "dense matrix",
    Real max_condition = dense_matrix_condition_error_threshold());

[[nodiscard]] std::size_t dense_matrix_rank(std::vector<Real> matrix,
                                            std::size_t rows,
                                            std::size_t cols);

[[nodiscard]] DensePseudoInverseResult rank_revealing_pseudo_inverse(
    std::span<const Real> matrix,
    std::size_t rows,
    std::size_t cols,
    std::string_view label = "dense matrix");

} // namespace math
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_MATH_DENSELINEARALGEBRA_H
