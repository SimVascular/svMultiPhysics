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

namespace svmp::FE::math {

namespace {

using DenseMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using RowMajorMatrix =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using ConstRowMajorMap = Eigen::Map<const RowMajorMatrix>;

ConstRowMajorMap map_row_major(std::span<const double> matrix,
                               std::size_t rows,
                               std::size_t cols) {
    return ConstRowMajorMap(matrix.data(),
                            static_cast<Eigen::Index>(rows),
                            static_cast<Eigen::Index>(cols));
}

void copy_to_row_major(const DenseMatrix& source, std::vector<double>& dest) {
    const auto rows = static_cast<std::size_t>(source.rows());
    const auto cols = static_cast<std::size_t>(source.cols());
    dest.resize(rows * cols);
    Eigen::Map<RowMajorMatrix>(dest.data(), source.rows(), source.cols()) = source;
}

} // namespace

struct DenseLUSolver::Impl {
    Eigen::PartialPivLU<DenseMatrix> lu;
};

DenseLUSolver::DenseLUSolver() = default;
DenseLUSolver::~DenseLUSolver() = default;
DenseLUSolver::DenseLUSolver(DenseLUSolver&&) noexcept = default;
DenseLUSolver& DenseLUSolver::operator=(DenseLUSolver&&) noexcept = default;

double dense_matrix_max_abs(std::span<const double> matrix) noexcept {
    double max_abs = double(0);
    for (const double value : matrix) {
        max_abs = std::max(max_abs, std::abs(value));
    }
    return max_abs;
}

double dense_matrix_pivot_tolerance(std::size_t rows,
                                  std::size_t cols,
                                  double max_abs,
                                  double multiplier) noexcept {
    const double size_scale = static_cast<double>(std::max<std::size_t>(rows, cols));
    const double value_scale = std::max(double(1), max_abs);
    return multiplier * std::numeric_limits<double>::epsilon() *
           std::max(double(1), size_scale) * value_scale;
}

double dense_matrix_singular_value_tolerance(std::size_t rows,
                                           std::size_t cols,
                                           double largest_singular_value,
                                           double multiplier) noexcept {
    const double size_scale = static_cast<double>(std::max<std::size_t>(rows, cols));
    return multiplier * std::numeric_limits<double>::epsilon() *
           std::max(double(1), size_scale) *
           std::max(double(1), largest_singular_value);
}

double dense_matrix_condition_fallback_threshold() noexcept {
    return double(1.0e12);
}

double dense_matrix_condition_error_threshold() noexcept {
    return double(1.0e14);
}

void DenseLUSolver::solve_in_place(std::span<double> rhs) const {
    solve_in_place(rhs, 1u);
}

void DenseLUSolver::solve_in_place(std::span<double> rhs,
                                   std::size_t rhs_count) const {
    ::svmp::check<FEException>(
        rhs_count > 0, error_message_label + ": dense solve requires at least one right-hand side");
    ::svmp::check<FEException>(
        rhs.size() == n * rhs_count, error_message_label + ": dense multi-RHS solve size mismatch");
    ::svmp::check<FEException>(
        impl && impl->lu.rows() == static_cast<Eigen::Index>(n), error_message_label + ": dense solver is not factorized");
    if (n == 0) {
        return;
    }

    Eigen::Map<RowMajorMatrix> rhs_map(rhs.data(),
                                       static_cast<Eigen::Index>(n),
                                       static_cast<Eigen::Index>(rhs_count));
    // Evaluate into a temporary: lu.solve cannot alias its argument.
    const DenseMatrix solution = impl->lu.solve(rhs_map);
    rhs_map = solution;
}

std::vector<double> DenseLUSolver::solve(std::span<const double> rhs) const {
    std::vector<double> x(rhs.begin(), rhs.end());
    solve_in_place(std::span<double>(x.data(), x.size()));
    return x;
}

DenseMatrixDiagnostics dense_matrix_diagnostics(
    std::span<const double> matrix,
    std::size_t rows,
    std::size_t cols,
    std::string_view error_message_label) {
    ::svmp::check<FEException>(
        matrix.size() == rows * cols, std::string(error_message_label) + ": diagnostic size mismatch");
    ::svmp::check<FEException>(
        rows > 0 && cols > 0, std::string(error_message_label) + ": diagnostics require a nonempty matrix");

    const DenseMatrix dense = map_row_major(matrix, rows, cols);
    Eigen::JacobiSVD<DenseMatrix> svd(dense);

    DenseMatrixDiagnostics diagnostics;
    const auto& singular_values = svd.singularValues();
    diagnostics.largest_singular_value =
        (singular_values.size() > 0) ? singular_values[0] : double(0);
    diagnostics.tolerance =
        dense_matrix_singular_value_tolerance(rows, cols,
                                              diagnostics.largest_singular_value);

    for (Eigen::Index i = 0; i < singular_values.size(); ++i) {
        const double sigma = singular_values[i];
        if (sigma <= diagnostics.tolerance) {
            continue;
        }
        ++diagnostics.rank;
        diagnostics.smallest_retained_singular_value = sigma;
    }

    const std::size_t full_rank = std::min(rows, cols);
    if (diagnostics.rank == full_rank &&
        diagnostics.smallest_retained_singular_value > double(0)) {
        diagnostics.condition_estimate =
            diagnostics.largest_singular_value /
            diagnostics.smallest_retained_singular_value;
    }
    return diagnostics;
}

DenseLUSolver factor_dense_matrix(std::vector<double> matrix,
                                  std::size_t n,
                                  std::string_view error_message_label) {
    ::svmp::check<FEException>(
        matrix.size() == n * n, std::string(error_message_label) + ": dense factorization size mismatch");

    DenseLUSolver solver;
    solver.n = n;
    solver.error_message_label = std::string(error_message_label);
    const double max_abs =
        dense_matrix_max_abs(std::span<const double>(matrix.data(), matrix.size()));
    solver.pivot_tolerance = dense_matrix_pivot_tolerance(n, n, max_abs);

    solver.impl = std::make_unique<DenseLUSolver::Impl>();
    solver.impl->lu.compute(map_row_major(matrix, n, n));

    // Partial pivoting leaves the pivots on the diagonal of the packed LU
    // factor; a pivot below the scale-aware tolerance marks rank deficiency.
    double max_pivot_abs = double(0);
    double min_pivot_abs = std::numeric_limits<double>::infinity();
    const auto diagonal = solver.impl->lu.matrixLU().diagonal();
    for (Eigen::Index col = 0; col < diagonal.size(); ++col) {
        const double pivot_magnitude = std::abs(diagonal[col]);
        ::svmp::check<FEException>(
            pivot_magnitude > solver.pivot_tolerance, solver.error_message_label + ": rank-deficient matrix (rank " +
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
    solver.min_pivot = std::isfinite(min_pivot_abs) ? min_pivot_abs : double(0);
    return solver;
}

DenseInverseResult invert_dense_matrix_with_diagnostics(
    std::vector<double> matrix,
    std::size_t n,
    std::string_view error_message_label) {
    ::svmp::check<FEException>(
        matrix.size() == n * n, std::string(error_message_label) + ": dense inverse size mismatch");
    std::vector<double> matrix_for_lu = matrix;
    const DenseLUSolver solver =
        factor_dense_matrix(std::move(matrix_for_lu), n, error_message_label);

    DenseInverseResult result;
    result.diagnostics =
        dense_matrix_diagnostics(std::span<const double>(matrix.data(), matrix.size()),
                                 n, n, error_message_label);

    if (std::isfinite(result.diagnostics.condition_estimate) &&
        result.diagnostics.condition_estimate > dense_matrix_condition_fallback_threshold()) {
        const DenseMatrix dense = map_row_major(matrix, n, n);
        Eigen::JacobiSVD<DenseMatrix> svd(dense,
                                          Eigen::ComputeFullU | Eigen::ComputeFullV);
        DenseMatrix sigma_inverse = DenseMatrix::Zero(static_cast<Eigen::Index>(n),
                                                      static_cast<Eigen::Index>(n));
        const auto& singular_values = svd.singularValues();
        for (Eigen::Index i = 0; i < singular_values.size(); ++i) {
            // Defensive: this branch runs only when condition_estimate is finite,
            // and dense_matrix_diagnostics leaves it infinite whenever it drops a
            // singular value (rank < full_rank). A sub-tolerance singular value
            // therefore cannot reach here in current code; the guard protects
            // against future refactors that derive the fallback condition differently.
            ::svmp::check<FEException>(
                singular_values[i] > result.diagnostics.tolerance, std::string(error_message_label) + ": high-condition SVD fallback encountered a dropped singular value");
            sigma_inverse(i, i) = double(1) / singular_values[i];
        }
        const DenseMatrix inverse = svd.matrixV() * sigma_inverse * svd.matrixU().transpose();
        copy_to_row_major(inverse, result.inverse);
        result.used_svd_fallback = true;
        return result;
    }

    const DenseMatrix inverse = solver.impl->lu.inverse();
    copy_to_row_major(inverse, result.inverse);
    return result;
}

void validate_dense_inverse_diagnostics(
    const DenseInverseResult& result,
    std::size_t expected_rank,
    std::string_view error_message_label,
    double max_condition) {
    ::svmp::check<FEException>(
        result.diagnostics.rank == expected_rank, std::string(error_message_label) + ": rank-deficient matrix (rank " +
            std::to_string(result.diagnostics.rank) + " of " +
            std::to_string(expected_rank) + ")");

    if (!std::isfinite(result.diagnostics.condition_estimate)) {
        return;
    }

    ::svmp::check<FEException>(
        result.diagnostics.condition_estimate <= max_condition, std::string(error_message_label) + ": condition estimate " +
            std::to_string(result.diagnostics.condition_estimate) +
            " exceeds supported threshold " + std::to_string(max_condition));
}

std::vector<double> invert_dense_matrix(std::vector<double> matrix,
                                      std::size_t n,
                                      std::string_view error_message_label) {
    const DenseLUSolver solver = factor_dense_matrix(std::move(matrix), n, error_message_label);
    const DenseMatrix inverse = solver.impl->lu.inverse();
    std::vector<double> result;
    copy_to_row_major(inverse, result);
    return result;
}

std::size_t dense_matrix_rank(std::vector<double> matrix,
                              std::size_t rows,
                              std::size_t cols) {
    ::svmp::check<FEException>(
        matrix.size() == rows * cols, "dense_matrix_rank: size mismatch");

    const DenseMatrix dense =
        map_row_major(std::span<const double>(matrix.data(), matrix.size()), rows, cols);
    Eigen::JacobiSVD<DenseMatrix> svd(dense);

    const auto& singular_values = svd.singularValues();
    const double largest =
        (singular_values.size() > 0) ? singular_values[0] : double(0);
    const double tolerance =
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
    std::span<const double> matrix,
    std::size_t rows,
    std::size_t cols,
    std::string_view error_message_label) {
    ::svmp::check<FEException>(
        matrix.size() == rows * cols, std::string(error_message_label) + ": pseudo-inverse size mismatch");
    ::svmp::check<FEException>(
        rows > 0 && cols > 0, std::string(error_message_label) + ": pseudo-inverse requires a nonempty matrix");

    const DenseMatrix dense = map_row_major(matrix, rows, cols);
    Eigen::JacobiSVD<DenseMatrix> svd(dense, Eigen::ComputeFullU | Eigen::ComputeFullV);

    DensePseudoInverseResult result;

    const auto& singular_values = svd.singularValues();
    result.largest_singular_value =
        (singular_values.size() > 0) ? singular_values[0] : double(0);
    result.tolerance =
        dense_matrix_singular_value_tolerance(rows, cols, result.largest_singular_value);

    DenseMatrix sigma_inverse = DenseMatrix::Zero(static_cast<Eigen::Index>(cols),
                                                  static_cast<Eigen::Index>(rows));
    for (Eigen::Index i = 0; i < singular_values.size(); ++i) {
        const double sigma = singular_values[i];
        if (sigma <= result.tolerance) {
            continue;
        }
        sigma_inverse(i, i) = double(1) / sigma;
        ++result.rank;
        result.smallest_retained_singular_value = sigma;
    }

    const DenseMatrix pseudo_inverse =
        svd.matrixV() * sigma_inverse * svd.matrixU().transpose();
    copy_to_row_major(pseudo_inverse, result.inverse);
    return result;
}

} // namespace svmp::FE::math
