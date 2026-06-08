// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "DenseLinearAlgebra.h"

#include "FEException.h"

#if defined(FE_HAS_EIGEN) && FE_HAS_EIGEN
#include <Eigen/Dense>
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

#define DENSE_LINALG_CHECK(condition, message) \
    ::svmp::FE::throw_if<::svmp::FE::FEException>(!(condition), SVMP_HERE, (message))

namespace svmp {
namespace FE {
namespace math {

namespace {

constexpr std::size_t kDenseSolveRhsBlock = 32u;

void materialize_inverse_from_solver(const DenseLUSolver& solver,
                                     std::vector<Real>& inverse) {
    const std::size_t n = solver.n;
    inverse.assign(n * n, Real(0));
    for (std::size_t diag = 0; diag < n; ++diag) {
        inverse[diag * n + diag] = Real(1);
    }
    solver.solve_in_place(std::span<Real>(inverse.data(), inverse.size()), n);
}

} // namespace

Real dense_matrix_max_abs(std::span<const Real> matrix) noexcept {
    Real max_abs = Real(0);
    for (const Real value : matrix) {
        max_abs = std::max(max_abs, std::abs(value));
    }
    return max_abs;
}

Real dense_matrix_pivot_tolerance(std::size_t rows,
                                  std::size_t cols,
                                  Real max_abs,
                                  Real multiplier) noexcept {
    const Real size_scale = static_cast<Real>(std::max<std::size_t>(rows, cols));
    const Real value_scale = std::max(Real(1), max_abs);
    return multiplier * std::numeric_limits<Real>::epsilon() *
           std::max(Real(1), size_scale) * value_scale;
}

Real dense_matrix_singular_value_tolerance(std::size_t rows,
                                           std::size_t cols,
                                           Real largest_singular_value,
                                           Real multiplier) noexcept {
    const Real size_scale = static_cast<Real>(std::max<std::size_t>(rows, cols));
    return multiplier * std::numeric_limits<Real>::epsilon() *
           std::max(Real(1), size_scale) *
           std::max(Real(1), largest_singular_value);
}

Real dense_matrix_condition_fallback_threshold() noexcept {
    return Real(1.0e12);
}

Real dense_matrix_condition_error_threshold() noexcept {
    return Real(1.0e14);
}

void DenseLUSolver::solve_in_place(std::span<Real> rhs) const {
    solve_in_place(rhs, 1u);
}

void DenseLUSolver::solve_in_place(std::span<Real> rhs,
                                   std::size_t rhs_count) const {
    DENSE_LINALG_CHECK(rhs_count > 0,
                             label + ": dense solve requires at least one right-hand side");
    DENSE_LINALG_CHECK(rhs.size() == n * rhs_count,
                             label + ": dense multi-RHS solve size mismatch");
    DENSE_LINALG_CHECK(lu.size() == n * n && pivots.size() == n,
                             label + ": dense solver is not factorized");

    for (std::size_t k = 0; k < n; ++k) {
        if (pivots[k] != k) {
            for (std::size_t block = 0; block < rhs_count; block += kDenseSolveRhsBlock) {
                const std::size_t end =
                    std::min(rhs_count, block + kDenseSolveRhsBlock);
                for (std::size_t r = block; r < end; ++r) {
                    std::swap(rhs[k * rhs_count + r],
                              rhs[pivots[k] * rhs_count + r]);
                }
            }
        }
    }

    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < row; ++col) {
            const Real factor = lu[row * n + col];
            for (std::size_t block = 0; block < rhs_count; block += kDenseSolveRhsBlock) {
                const std::size_t end =
                    std::min(rhs_count, block + kDenseSolveRhsBlock);
                for (std::size_t r = block; r < end; ++r) {
                    rhs[row * rhs_count + r] -= factor * rhs[col * rhs_count + r];
                }
            }
        }
    }

    for (std::size_t rev = 0; rev < n; ++rev) {
        const std::size_t row = n - 1u - rev;
        for (std::size_t col = row + 1u; col < n; ++col) {
            const Real factor = lu[row * n + col];
            for (std::size_t block = 0; block < rhs_count; block += kDenseSolveRhsBlock) {
                const std::size_t end =
                    std::min(rhs_count, block + kDenseSolveRhsBlock);
                for (std::size_t r = block; r < end; ++r) {
                    rhs[row * rhs_count + r] -= factor * rhs[col * rhs_count + r];
                }
            }
        }
        const Real pivot = lu[row * n + row];
        DENSE_LINALG_CHECK(
            std::abs(pivot) > pivot_tolerance,
            label + ": zero pivot during dense solve");
        for (std::size_t block = 0; block < rhs_count; block += kDenseSolveRhsBlock) {
            const std::size_t end =
                std::min(rhs_count, block + kDenseSolveRhsBlock);
            for (std::size_t r = block; r < end; ++r) {
                rhs[row * rhs_count + r] /= pivot;
            }
        }
    }
}

std::vector<Real> DenseLUSolver::solve(std::span<const Real> rhs) const {
    std::vector<Real> x(rhs.begin(), rhs.end());
    solve_in_place(std::span<Real>(x.data(), x.size()));
    return x;
}

DenseMatrixDiagnostics dense_matrix_diagnostics(
    std::span<const Real> matrix,
    std::size_t rows,
    std::size_t cols,
    std::string_view label) {
    DENSE_LINALG_CHECK(matrix.size() == rows * cols,
                             std::string(label) + ": diagnostic size mismatch");
    DENSE_LINALG_CHECK(rows > 0 && cols > 0,
                             std::string(label) + ": diagnostics require a nonempty matrix");

#if defined(FE_HAS_EIGEN) && FE_HAS_EIGEN
    using RowMajorMatrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using Matrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
    const Eigen::Map<const RowMajorMatrix> A(matrix.data(),
                                             static_cast<Eigen::Index>(rows),
                                             static_cast<Eigen::Index>(cols));
    const Matrix dense = A;
    Eigen::JacobiSVD<Matrix> svd(dense);

    DenseMatrixDiagnostics diagnostics;
    const auto& singular_values = svd.singularValues();
    diagnostics.largest_singular_value =
        (singular_values.size() > 0) ? singular_values[0] : Real(0);
    diagnostics.tolerance =
        dense_matrix_singular_value_tolerance(rows, cols,
                                              diagnostics.largest_singular_value);

    for (Eigen::Index i = 0; i < singular_values.size(); ++i) {
        const Real sigma = singular_values[i];
        if (sigma <= diagnostics.tolerance) {
            continue;
        }
        ++diagnostics.rank;
        diagnostics.smallest_retained_singular_value = sigma;
    }

    const std::size_t full_rank = std::min(rows, cols);
    if (diagnostics.rank == full_rank &&
        diagnostics.smallest_retained_singular_value > Real(0)) {
        diagnostics.condition_estimate =
            diagnostics.largest_singular_value /
            diagnostics.smallest_retained_singular_value;
    }
    return diagnostics;
#else
    DenseMatrixDiagnostics diagnostics;
    diagnostics.largest_singular_value = dense_matrix_max_abs(matrix);
    diagnostics.tolerance =
        dense_matrix_pivot_tolerance(rows, cols, diagnostics.largest_singular_value);
    diagnostics.rank =
        dense_matrix_rank(std::vector<Real>(matrix.begin(), matrix.end()), rows, cols);
    const std::size_t full_rank = std::min(rows, cols);
    if (diagnostics.rank == full_rank) {
        diagnostics.smallest_retained_singular_value = diagnostics.tolerance;
    }
    // Exact condition estimates require SVD diagnostics. In Eigen-disabled
    // builds this stays explicit instead of relying on a misleading estimate.
    diagnostics.condition_estimate = std::numeric_limits<Real>::infinity();
    return diagnostics;
#endif
}

DenseLUSolver factor_dense_matrix(std::vector<Real> matrix,
                                  std::size_t n,
                                  std::string_view label) {
    DENSE_LINALG_CHECK(matrix.size() == n * n,
                             std::string(label) + ": dense factorization size mismatch");

    DenseLUSolver solver;
    solver.n = n;
    solver.lu = std::move(matrix);
    solver.pivots.resize(n);
    const Real max_abs = dense_matrix_max_abs(solver.lu);
    solver.pivot_tolerance =
        dense_matrix_pivot_tolerance(n, n, max_abs);
    solver.label = std::string(label);

    Real max_pivot_abs = Real(0);
    Real min_pivot_abs = std::numeric_limits<Real>::infinity();
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot_row = col;
        Real pivot_abs = std::abs(solver.lu[col * n + col]);
        for (std::size_t row = col + 1; row < n; ++row) {
            const Real candidate = std::abs(solver.lu[row * n + col]);
            if (candidate > pivot_abs) {
                pivot_abs = candidate;
                pivot_row = row;
            }
        }

        DENSE_LINALG_CHECK(
            pivot_abs > solver.pivot_tolerance,
            solver.label + ": rank-deficient matrix (rank " +
                std::to_string(col) + " of " + std::to_string(n) +
                ", pivot below scale-aware tolerance " +
                std::to_string(solver.pivot_tolerance) + ")");

        solver.pivots[col] = pivot_row;
        if (pivot_row != col) {
            for (std::size_t j = 0; j < n; ++j) {
                std::swap(solver.lu[col * n + j], solver.lu[pivot_row * n + j]);
            }
        }

        const Real pivot = solver.lu[col * n + col];
        DENSE_LINALG_CHECK(
            std::abs(pivot) > solver.pivot_tolerance,
            solver.label + ": zero pivot after row exchange");
        const Real pivot_magnitude = std::abs(pivot);
        max_pivot_abs = std::max(max_pivot_abs, pivot_magnitude);
        min_pivot_abs = std::min(min_pivot_abs, pivot_magnitude);

        for (std::size_t row = col + 1; row < n; ++row) {
            const Real factor = solver.lu[row * n + col] / pivot;
            solver.lu[row * n + col] = factor;
            for (std::size_t j = col + 1; j < n; ++j) {
                solver.lu[row * n + j] -= factor * solver.lu[col * n + j];
            }
        }
    }

    solver.diagnostics.rank = n;
    solver.diagnostics.tolerance = solver.pivot_tolerance;
    solver.diagnostics.largest_singular_value = max_abs;
    solver.diagnostics.smallest_retained_singular_value =
        std::isfinite(min_pivot_abs) ? min_pivot_abs : Real(0);
    if (solver.diagnostics.smallest_retained_singular_value > Real(0)) {
        solver.diagnostics.condition_estimate =
            max_pivot_abs / solver.diagnostics.smallest_retained_singular_value;
    }
    return solver;
}

DenseInverseResult invert_dense_matrix_with_diagnostics(
    std::vector<Real> matrix,
    std::size_t n,
    std::string_view label) {
    DENSE_LINALG_CHECK(matrix.size() == n * n,
                             std::string(label) + ": dense inverse size mismatch");
    std::vector<Real> matrix_for_lu = matrix;
    const DenseLUSolver solver =
        factor_dense_matrix(std::move(matrix_for_lu), n, label);

    DenseInverseResult result;
    result.diagnostics =
        dense_matrix_diagnostics(std::span<const Real>(matrix.data(), matrix.size()),
                                 n, n, label);

#if defined(FE_HAS_EIGEN) && FE_HAS_EIGEN
    if (std::isfinite(solver.diagnostics.condition_estimate) &&
        std::isfinite(result.diagnostics.condition_estimate) &&
        result.diagnostics.condition_estimate > dense_matrix_condition_fallback_threshold()) {
        using RowMajorMatrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        using Matrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
        const Eigen::Map<const RowMajorMatrix> A(matrix.data(),
                                                 static_cast<Eigen::Index>(n),
                                                 static_cast<Eigen::Index>(n));
        const Matrix dense = A;
        Eigen::JacobiSVD<Matrix> svd(dense,
                                     Eigen::ComputeFullU | Eigen::ComputeFullV);
        Matrix sigma_inverse = Matrix::Zero(static_cast<Eigen::Index>(n),
                                            static_cast<Eigen::Index>(n));
        const auto& singular_values = svd.singularValues();
        for (Eigen::Index i = 0; i < singular_values.size(); ++i) {
            DENSE_LINALG_CHECK(
                singular_values[i] > solver.diagnostics.tolerance,
                std::string(label) + ": high-condition SVD fallback encountered a dropped singular value");
            sigma_inverse(i, i) = Real(1) / singular_values[i];
        }
        const Matrix inverse = svd.matrixV() * sigma_inverse * svd.matrixU().transpose();
        result.inverse.assign(n * n, Real(0));
        for (std::size_t row = 0; row < n; ++row) {
            for (std::size_t col = 0; col < n; ++col) {
                result.inverse[row * n + col] =
                    inverse(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col));
            }
        }
        result.used_svd_fallback = true;
        return result;
    }
#endif

    materialize_inverse_from_solver(solver, result.inverse);
    return result;
}

void validate_dense_inverse_diagnostics(
    const DenseInverseResult& result,
    std::size_t expected_rank,
    std::string_view label,
    Real max_condition) {
    DENSE_LINALG_CHECK(
        result.diagnostics.rank == expected_rank,
        std::string(label) + ": rank-deficient matrix (rank " +
            std::to_string(result.diagnostics.rank) + " of " +
            std::to_string(expected_rank) + ")");

    if (!std::isfinite(result.diagnostics.condition_estimate)) {
        return;
    }

    DENSE_LINALG_CHECK(
        result.diagnostics.condition_estimate <= max_condition,
        std::string(label) + ": condition estimate " +
            std::to_string(result.diagnostics.condition_estimate) +
            " exceeds supported threshold " + std::to_string(max_condition));
}

std::vector<Real> invert_dense_matrix(std::vector<Real> matrix,
                                      std::size_t n,
                                      std::string_view label) {
    const DenseLUSolver solver = factor_dense_matrix(std::move(matrix), n, label);
    std::vector<Real> inverse;
    materialize_inverse_from_solver(solver, inverse);
    return inverse;
}

std::size_t dense_matrix_rank(std::vector<Real> matrix,
                              std::size_t rows,
                              std::size_t cols) {
    DENSE_LINALG_CHECK(matrix.size() == rows * cols,
                             "dense_matrix_rank: size mismatch");
    const Real tolerance =
        dense_matrix_pivot_tolerance(rows, cols, dense_matrix_max_abs(matrix));

    std::size_t rank = 0;
    std::size_t pivot_row = 0;
    for (std::size_t col = 0; col < cols && pivot_row < rows; ++col) {
        std::size_t best_row = pivot_row;
        Real best_abs = std::abs(matrix[pivot_row * cols + col]);
        for (std::size_t row = pivot_row + 1; row < rows; ++row) {
            const Real candidate = std::abs(matrix[row * cols + col]);
            if (candidate > best_abs) {
                best_abs = candidate;
                best_row = row;
            }
        }
        if (best_abs <= tolerance) {
            continue;
        }

        if (best_row != pivot_row) {
            for (std::size_t c = col; c < cols; ++c) {
                std::swap(matrix[pivot_row * cols + c], matrix[best_row * cols + c]);
            }
        }

        const Real pivot = matrix[pivot_row * cols + col];
        for (std::size_t row = pivot_row + 1; row < rows; ++row) {
            const Real factor = matrix[row * cols + col] / pivot;
            if (std::abs(factor) <= tolerance) {
                matrix[row * cols + col] = Real(0);
                continue;
            }
            matrix[row * cols + col] = Real(0);
            for (std::size_t c = col + 1; c < cols; ++c) {
                matrix[row * cols + c] -= factor * matrix[pivot_row * cols + c];
            }
        }

        ++rank;
        ++pivot_row;
    }
    return rank;
}

DensePseudoInverseResult rank_revealing_pseudo_inverse(
    std::span<const Real> matrix,
    std::size_t rows,
    std::size_t cols,
    std::string_view label) {
    DENSE_LINALG_CHECK(matrix.size() == rows * cols,
                             std::string(label) + ": pseudo-inverse size mismatch");
    DENSE_LINALG_CHECK(rows > 0 && cols > 0,
                             std::string(label) + ": pseudo-inverse requires a nonempty matrix");

#if defined(FE_HAS_EIGEN) && FE_HAS_EIGEN
    using RowMajorMatrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using Matrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
    const Eigen::Map<const RowMajorMatrix> A(matrix.data(),
                                             static_cast<Eigen::Index>(rows),
                                             static_cast<Eigen::Index>(cols));
    const Matrix dense = A;
    Eigen::JacobiSVD<Matrix> svd(dense, Eigen::ComputeFullU | Eigen::ComputeFullV);

    DensePseudoInverseResult result;
    result.inverse.assign(cols * rows, Real(0));

    const auto& singular_values = svd.singularValues();
    result.largest_singular_value =
        (singular_values.size() > 0) ? singular_values[0] : Real(0);
    result.tolerance =
        dense_matrix_singular_value_tolerance(rows, cols, result.largest_singular_value);

    Matrix sigma_inverse = Matrix::Zero(static_cast<Eigen::Index>(cols),
                                        static_cast<Eigen::Index>(rows));
    for (Eigen::Index i = 0; i < singular_values.size(); ++i) {
        const Real sigma = singular_values[i];
        if (sigma <= result.tolerance) {
            continue;
        }
        sigma_inverse(i, i) = Real(1) / sigma;
        ++result.rank;
        result.smallest_retained_singular_value = sigma;
    }

    const Matrix pseudo_inverse =
        svd.matrixV() * sigma_inverse * svd.matrixU().transpose();
    for (std::size_t r = 0; r < cols; ++r) {
        for (std::size_t c = 0; c < rows; ++c) {
            result.inverse[r * rows + c] =
                pseudo_inverse(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(c));
        }
    }
    return result;
#else
    DENSE_LINALG_CHECK(
        false,
        std::string(label) +
            ": rank-revealing pseudo-inverse requires FE_ENABLE_EIGEN");
    return {};
#endif
}

} // namespace math
} // namespace FE
} // namespace svmp

#undef DENSE_LINALG_CHECK
