// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_MATH_DENSELINEARALGEBRA_H
#define SVMP_FE_MATH_DENSELINEARALGEBRA_H

#include "Types.h"

#include <cstddef>
#include <limits>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace svmp {
namespace FE {
namespace math {

// Dense solve, inverse, rank, and pseudo-inverse support for FE construction
// utilities. Matrices are row-major: matrix[row * cols + col].
[[nodiscard]] double dense_matrix_max_abs(std::span<const double> matrix) noexcept;

[[nodiscard]] double dense_matrix_pivot_tolerance(std::size_t rows,
                                                std::size_t cols,
                                                double max_abs,
                                                double multiplier = double(64)) noexcept;

[[nodiscard]] double dense_matrix_singular_value_tolerance(std::size_t rows,
                                                         std::size_t cols,
                                                         double largest_singular_value,
                                                         double multiplier = double(64)) noexcept;

struct DensePseudoInverseResult {
    std::vector<double> inverse;
    std::size_t rank{0};
    double tolerance{0};
    double largest_singular_value{0};
    double smallest_retained_singular_value{0};
};

struct DenseMatrixDiagnostics {
    std::size_t rank{0};
    double tolerance{0};
    double largest_singular_value{0};
    double smallest_retained_singular_value{0};
    double condition_estimate{std::numeric_limits<double>::infinity()};
};

struct DenseInverseResult {
    std::vector<double> inverse;
    DenseMatrixDiagnostics diagnostics;
    bool used_svd_fallback{false};
};

[[nodiscard]] double dense_matrix_condition_fallback_threshold() noexcept;
[[nodiscard]] double dense_matrix_condition_error_threshold() noexcept;

struct DenseLUSolver {
    struct Impl;

    DenseLUSolver();
    ~DenseLUSolver();
    DenseLUSolver(DenseLUSolver&&) noexcept;
    DenseLUSolver& operator=(DenseLUSolver&&) noexcept;
    DenseLUSolver(const DenseLUSolver&) = delete;
    DenseLUSolver& operator=(const DenseLUSolver&) = delete;

    std::size_t n{0};
    DenseMatrixDiagnostics diagnostics;
    double pivot_tolerance{0};
    double min_pivot{0};
    double max_pivot{0};
    std::string label;
    std::unique_ptr<Impl> impl;

    [[nodiscard]] bool empty() const noexcept { return n == 0; }

    void solve_in_place(std::span<double> rhs) const;
    void solve_in_place(std::span<double> rhs, std::size_t rhs_count) const;
    [[nodiscard]] std::vector<double> solve(std::span<const double> rhs) const;
};

// Inverses and pseudo-inverses keep the same row-major convention for their
// returned dimensions.
[[nodiscard]] DenseMatrixDiagnostics dense_matrix_diagnostics(
    std::span<const double> matrix,
    std::size_t rows,
    std::size_t cols,
    std::string_view label = "dense matrix");

[[nodiscard]] DenseLUSolver factor_dense_matrix(std::vector<double> matrix,
                                                std::size_t n,
                                                std::string_view label = "dense matrix");

[[nodiscard]] std::vector<double> invert_dense_matrix(std::vector<double> matrix,
                                                    std::size_t n,
                                                    std::string_view label = "dense matrix");

[[nodiscard]] DenseInverseResult invert_dense_matrix_with_diagnostics(
    std::vector<double> matrix,
    std::size_t n,
    std::string_view label = "dense matrix");

void validate_dense_inverse_diagnostics(
    const DenseInverseResult& result,
    std::size_t expected_rank,
    std::string_view label = "dense matrix",
    double max_condition = dense_matrix_condition_error_threshold());

[[nodiscard]] std::size_t dense_matrix_rank(std::vector<double> matrix,
                                            std::size_t rows,
                                            std::size_t cols);

[[nodiscard]] DensePseudoInverseResult rank_revealing_pseudo_inverse(
    std::span<const double> matrix,
    std::size_t rows,
    std::size_t cols,
    std::string_view label = "dense matrix");

} // namespace math
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_MATH_DENSELINEARALGEBRA_H
