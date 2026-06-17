// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "DenseLinearAlgebra.h"

#include "FEException.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

namespace svmp {
namespace FE {
namespace math {

namespace {

using DenseMatrix = DenseLUSolver::DenseMatrix;
using RowMajorMatrix =
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using ConstRowMajorMap = Eigen::Map<const RowMajorMatrix>;

ConstRowMajorMap map_row_major(std::span<const Real> matrix,
                               std::size_t rows,
                               std::size_t cols) {
    return ConstRowMajorMap(matrix.data(),
                            static_cast<Eigen::Index>(rows),
                            static_cast<Eigen::Index>(cols));
}

void copy_to_row_major(const DenseMatrix& source, std::vector<Real>& dest) {
    const auto rows = static_cast<std::size_t>(source.rows());
    const auto cols = static_cast<std::size_t>(source.cols());
    dest.resize(rows * cols);
    Eigen::Map<RowMajorMatrix>(dest.data(), source.rows(), source.cols()) = source;
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
    ::svmp::FE::check_arg<FEException>(
        rhs_count > 0, SVMP_HERE,
        label + ": dense solve requires at least one right-hand side");
    ::svmp::FE::check_arg<FEException>(
        rhs.size() == n * rhs_count, SVMP_HERE,
        label + ": dense multi-RHS solve size mismatch");
    ::svmp::FE::check_arg<FEException>(
        lu.rows() == static_cast<Eigen::Index>(n), SVMP_HERE,
        label + ": dense solver is not factorized");
    if (n == 0) {
        return;
    }

    Eigen::Map<RowMajorMatrix> rhs_map(rhs.data(),
                                       static_cast<Eigen::Index>(n),
                                       static_cast<Eigen::Index>(rhs_count));
    // Evaluate into a temporary: lu.solve cannot alias its argument.
    const DenseMatrix solution = lu.solve(rhs_map);
    rhs_map = solution;
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
    ::svmp::FE::check_arg<FEException>(
        matrix.size() == rows * cols, SVMP_HERE,
        std::string(label) + ": diagnostic size mismatch");
    ::svmp::FE::check_arg<FEException>(
        rows > 0 && cols > 0, SVMP_HERE,
        std::string(label) + ": diagnostics require a nonempty matrix");

    const DenseMatrix dense = map_row_major(matrix, rows, cols);
    Eigen::JacobiSVD<DenseMatrix> svd(dense);

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
}

DenseLUSolver factor_dense_matrix(std::vector<Real> matrix,
                                  std::size_t n,
                                  std::string_view label) {
    ::svmp::FE::check_arg<FEException>(
        matrix.size() == n * n, SVMP_HERE,
        std::string(label) + ": dense factorization size mismatch");

    DenseLUSolver solver;
    solver.n = n;
    solver.label = std::string(label);
    const Real max_abs =
        dense_matrix_max_abs(std::span<const Real>(matrix.data(), matrix.size()));
    solver.pivot_tolerance = dense_matrix_pivot_tolerance(n, n, max_abs);

    solver.lu.compute(map_row_major(matrix, n, n));

    // Partial pivoting leaves the pivots on the diagonal of the packed LU
    // factor; a pivot below the scale-aware tolerance marks rank deficiency.
    Real max_pivot_abs = Real(0);
    Real min_pivot_abs = std::numeric_limits<Real>::infinity();
    const auto diagonal = solver.lu.matrixLU().diagonal();
    for (Eigen::Index col = 0; col < diagonal.size(); ++col) {
        const Real pivot_magnitude = std::abs(diagonal[col]);
        ::svmp::FE::check_arg<FEException>(
            pivot_magnitude > solver.pivot_tolerance, SVMP_HERE,
            solver.label + ": rank-deficient matrix (rank " +
                std::to_string(col) + " of " + std::to_string(n) +
                ", pivot below scale-aware tolerance " +
                std::to_string(solver.pivot_tolerance) + ")");
        max_pivot_abs = std::max(max_pivot_abs, pivot_magnitude);
        min_pivot_abs = std::min(min_pivot_abs, pivot_magnitude);
    }

    // PartialPivLU is not rank-revealing, so expose only what the pivots
    // legitimately convey: the factorization passed the pivot-tolerance check
    // above (full rank) and the pivot magnitudes.
    solver.diagnostics.rank = n;
    solver.diagnostics.tolerance = solver.pivot_tolerance;
    solver.max_pivot = max_pivot_abs;
    solver.min_pivot = std::isfinite(min_pivot_abs) ? min_pivot_abs : Real(0);
    return solver;
}

DenseInverseResult invert_dense_matrix_with_diagnostics(
    std::vector<Real> matrix,
    std::size_t n,
    std::string_view label) {
    ::svmp::FE::check_arg<FEException>(
        matrix.size() == n * n, SVMP_HERE,
        std::string(label) + ": dense inverse size mismatch");
    std::vector<Real> matrix_for_lu = matrix;
    const DenseLUSolver solver =
        factor_dense_matrix(std::move(matrix_for_lu), n, label);

    DenseInverseResult result;
    result.diagnostics =
        dense_matrix_diagnostics(std::span<const Real>(matrix.data(), matrix.size()),
                                 n, n, label);

    if (std::isfinite(result.diagnostics.condition_estimate) &&
        result.diagnostics.condition_estimate > dense_matrix_condition_fallback_threshold()) {
        const DenseMatrix dense = map_row_major(matrix, n, n);
        Eigen::JacobiSVD<DenseMatrix> svd(dense,
                                          Eigen::ComputeFullU | Eigen::ComputeFullV);
        DenseMatrix sigma_inverse = DenseMatrix::Zero(static_cast<Eigen::Index>(n),
                                                      static_cast<Eigen::Index>(n));
        const auto& singular_values = svd.singularValues();
        for (Eigen::Index i = 0; i < singular_values.size(); ++i) {
            ::svmp::FE::check_arg<FEException>(
                singular_values[i] > result.diagnostics.tolerance, SVMP_HERE,
                std::string(label) + ": high-condition SVD fallback encountered a dropped singular value");
            sigma_inverse(i, i) = Real(1) / singular_values[i];
        }
        const DenseMatrix inverse = svd.matrixV() * sigma_inverse * svd.matrixU().transpose();
        copy_to_row_major(inverse, result.inverse);
        result.used_svd_fallback = true;
        return result;
    }

    const DenseMatrix inverse = solver.lu.inverse();
    copy_to_row_major(inverse, result.inverse);
    return result;
}

void validate_dense_inverse_diagnostics(
    const DenseInverseResult& result,
    std::size_t expected_rank,
    std::string_view label,
    Real max_condition) {
    ::svmp::FE::check_arg<FEException>(
        result.diagnostics.rank == expected_rank, SVMP_HERE,
        std::string(label) + ": rank-deficient matrix (rank " +
            std::to_string(result.diagnostics.rank) + " of " +
            std::to_string(expected_rank) + ")");

    if (!std::isfinite(result.diagnostics.condition_estimate)) {
        return;
    }

    ::svmp::FE::check_arg<FEException>(
        result.diagnostics.condition_estimate <= max_condition, SVMP_HERE,
        std::string(label) + ": condition estimate " +
            std::to_string(result.diagnostics.condition_estimate) +
            " exceeds supported threshold " + std::to_string(max_condition));
}

std::vector<Real> invert_dense_matrix(std::vector<Real> matrix,
                                      std::size_t n,
                                      std::string_view label) {
    const DenseLUSolver solver = factor_dense_matrix(std::move(matrix), n, label);
    const DenseMatrix inverse = solver.lu.inverse();
    std::vector<Real> result;
    copy_to_row_major(inverse, result);
    return result;
}

std::size_t dense_matrix_rank(std::vector<Real> matrix,
                              std::size_t rows,
                              std::size_t cols) {
    ::svmp::FE::check_arg<FEException>(
        matrix.size() == rows * cols, SVMP_HERE,
        "dense_matrix_rank: size mismatch");

    const DenseMatrix dense =
        map_row_major(std::span<const Real>(matrix.data(), matrix.size()), rows, cols);
    Eigen::JacobiSVD<DenseMatrix> svd(dense);

    const auto& singular_values = svd.singularValues();
    const Real largest =
        (singular_values.size() > 0) ? singular_values[0] : Real(0);
    const Real tolerance =
        dense_matrix_singular_value_tolerance(rows, cols, largest);

    std::size_t rank = 0;
    for (Eigen::Index i = 0; i < singular_values.size(); ++i) {
        if (singular_values[i] > tolerance) {
            ++rank;
        }
    }
    return rank;
}

DensePseudoInverseResult rank_revealing_pseudo_inverse(
    std::span<const Real> matrix,
    std::size_t rows,
    std::size_t cols,
    std::string_view label) {
    ::svmp::FE::check_arg<FEException>(
        matrix.size() == rows * cols, SVMP_HERE,
        std::string(label) + ": pseudo-inverse size mismatch");
    ::svmp::FE::check_arg<FEException>(
        rows > 0 && cols > 0, SVMP_HERE,
        std::string(label) + ": pseudo-inverse requires a nonempty matrix");

    const DenseMatrix dense = map_row_major(matrix, rows, cols);
    Eigen::JacobiSVD<DenseMatrix> svd(dense, Eigen::ComputeFullU | Eigen::ComputeFullV);

    DensePseudoInverseResult result;

    const auto& singular_values = svd.singularValues();
    result.largest_singular_value =
        (singular_values.size() > 0) ? singular_values[0] : Real(0);
    result.tolerance =
        dense_matrix_singular_value_tolerance(rows, cols, result.largest_singular_value);

    DenseMatrix sigma_inverse = DenseMatrix::Zero(static_cast<Eigen::Index>(cols),
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

    const DenseMatrix pseudo_inverse =
        svd.matrixV() * sigma_inverse * svd.matrixU().transpose();
    copy_to_row_major(pseudo_inverse, result.inverse);
    return result;
}

} // namespace math
} // namespace FE
} // namespace svmp
