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

namespace svmp::FE::math {

// Dense solve, inverse, rank, and pseudo-inverse support for FE construction
// utilities. Matrices are row-major: matrix[row * cols + col]. The
// error_message_label argument on the routines below is used only to prefix the
// diagnostic message of any exception they throw.

/**
 * @brief Largest absolute entry of a dense matrix.
 * @ingroup FE_Math
 * @param matrix Row-major matrix entries.
 * @return Maximum of |entry| over all entries, or 0 for an empty matrix.
 */
[[nodiscard]] double dense_matrix_max_abs(std::span<const double> matrix) noexcept;

/**
 * @brief Scale-aware pivot tolerance for dense factorization.
 * @ingroup FE_Math
 *
 * @details Proportional to machine epsilon scaled by the matrix size and
 * magnitude; pivots below it are treated as rank-deficient.
 *
 * @param rows Row count.
 * @param cols Column count.
 * @param max_abs Largest absolute matrix entry (see dense_matrix_max_abs()).
 * @param multiplier Safety factor applied to the epsilon-scaled tolerance.
 * @return Pivot magnitude threshold.
 */
[[nodiscard]] double dense_matrix_pivot_tolerance(std::size_t rows,
                                                std::size_t cols,
                                                double max_abs,
                                                double multiplier = double(64)) noexcept;

/**
 * @brief Scale-aware singular-value tolerance for rank decisions.
 * @ingroup FE_Math
 *
 * @details Singular values at or below the returned tolerance are treated as
 * zero when computing rank or a pseudo-inverse.
 *
 * @param rows Row count.
 * @param cols Column count.
 * @param largest_singular_value Largest singular value of the matrix.
 * @param multiplier Safety factor applied to the epsilon-scaled tolerance.
 * @return Singular-value threshold.
 */
[[nodiscard]] double dense_matrix_singular_value_tolerance(std::size_t rows,
                                                         std::size_t cols,
                                                         double largest_singular_value,
                                                         double multiplier = double(64)) noexcept;

/** @brief Result of a rank-revealing pseudo-inverse. @ingroup FE_Math */
struct DensePseudoInverseResult {
    std::vector<double> inverse;                 ///< Row-major pseudo-inverse.
    std::size_t rank{0};                         ///< Numerical rank at the chosen tolerance.
    double tolerance{0};                         ///< Singular-value tolerance used.
    double largest_singular_value{0};            ///< Largest singular value.
    double smallest_retained_singular_value{0};  ///< Smallest singular value kept.
};

/** @brief SVD-based conditioning and rank diagnostics for a dense matrix. @ingroup FE_Math */
struct DenseMatrixDiagnostics {
    std::size_t rank{0};                         ///< Numerical rank at @ref tolerance.
    double tolerance{0};                         ///< Singular-value tolerance used.
    double largest_singular_value{0};            ///< Largest singular value.
    double smallest_retained_singular_value{0};  ///< Smallest singular value kept.
    double condition_estimate{std::numeric_limits<double>::infinity()};  ///< Condition estimate; infinite when rank-deficient.
};

/** @brief A dense inverse together with its diagnostics. @ingroup FE_Math */
struct DenseInverseResult {
    std::vector<double> inverse;        ///< Row-major inverse.
    DenseMatrixDiagnostics diagnostics; ///< Conditioning/rank diagnostics of the input.
    bool used_svd_fallback{false};      ///< True when an SVD fallback was used for a high-condition matrix.
};

/**
 * @brief Condition estimate above which the inverse switches to an SVD fallback.
 * @ingroup FE_Math
 * @return The fallback condition-number threshold.
 */
[[nodiscard]] double dense_matrix_condition_fallback_threshold() noexcept;
/**
 * @brief Condition estimate above which validation rejects a dense inverse.
 * @ingroup FE_Math
 * @return The error condition-number threshold.
 */
[[nodiscard]] double dense_matrix_condition_error_threshold() noexcept;

/**
 * @brief LU factorization of a dense square matrix with a cached pivot summary.
 * @ingroup FE_Math
 *
 * @details Produced by factor_dense_matrix(); move-only because it owns the Eigen
 * factorization. @ref error_message_label prefixes the messages of exceptions
 * thrown by the solve methods.
 */
struct DenseLUSolver {
    struct Impl;

    DenseLUSolver();
    ~DenseLUSolver();
    DenseLUSolver(DenseLUSolver&&) noexcept;
    DenseLUSolver& operator=(DenseLUSolver&&) noexcept;
    DenseLUSolver(const DenseLUSolver&) = delete;
    DenseLUSolver& operator=(const DenseLUSolver&) = delete;

    std::size_t n{0};                    ///< Matrix dimension.
    DenseMatrixDiagnostics diagnostics;  ///< Pivot-derived diagnostics (rank, tolerance).
    double pivot_tolerance{0};           ///< Scale-aware pivot tolerance used.
    double min_pivot{0};                 ///< Smallest pivot magnitude.
    double max_pivot{0};                 ///< Largest pivot magnitude.
    std::string error_message_label;     ///< Prefix for solve-time exception messages.
    std::unique_ptr<Impl> impl;          ///< Eigen factorization (pimpl).

    /**
     * @brief Whether the factorization is empty (n == 0).
     * @return True when no matrix has been factored.
     */
    [[nodiscard]] bool empty() const noexcept { return n == 0; }

    /**
     * @brief Solve A x = rhs in place for a single right-hand side.
     * @param rhs On entry the right-hand side; on return the solution (size n).
     */
    void solve_in_place(std::span<double> rhs) const;
    /**
     * @brief Solve A X = RHS in place for several right-hand sides.
     * @param rhs Row-major block of size n * rhs_count (solutions on return).
     * @param rhs_count Number of right-hand sides.
     */
    void solve_in_place(std::span<double> rhs, std::size_t rhs_count) const;
    /**
     * @brief Solve A x = rhs and return the solution.
     * @param rhs Right-hand side of size n.
     * @return Solution vector of size n.
     */
    [[nodiscard]] std::vector<double> solve(std::span<const double> rhs) const;
};

// Inverses and pseudo-inverses keep the same row-major convention for their
// returned dimensions.

/**
 * @brief SVD-based rank and conditioning diagnostics for a dense matrix.
 * @ingroup FE_Math
 * @param matrix Row-major matrix of size rows * cols.
 * @param rows Row count.
 * @param cols Column count.
 * @param error_message_label Prefix for the message of any exception thrown.
 * @return Rank, tolerance, singular-value, and condition diagnostics.
 * @throws FEException If the matrix size is inconsistent or the matrix is empty.
 */
[[nodiscard]] DenseMatrixDiagnostics dense_matrix_diagnostics(
    std::span<const double> matrix,
    std::size_t rows,
    std::size_t cols,
    std::string_view error_message_label = "dense matrix");

/**
 * @brief LU-factor a dense square matrix.
 * @ingroup FE_Math
 * @param matrix Row-major n * n matrix (consumed).
 * @param n Matrix dimension.
 * @param error_message_label Prefix for the message of any exception thrown.
 * @return The factorization.
 * @throws FEException If the size is inconsistent or the matrix is rank-deficient.
 */
[[nodiscard]] DenseLUSolver factor_dense_matrix(std::vector<double> matrix,
                                                std::size_t n,
                                                std::string_view error_message_label = "dense matrix");

/**
 * @brief Invert a dense square matrix.
 * @ingroup FE_Math
 * @param matrix Row-major n * n matrix (consumed).
 * @param n Matrix dimension.
 * @param error_message_label Prefix for the message of any exception thrown.
 * @return Row-major inverse of size n * n.
 * @throws FEException If the size is inconsistent or the matrix is singular.
 */
[[nodiscard]] std::vector<double> invert_dense_matrix(std::vector<double> matrix,
                                                    std::size_t n,
                                                    std::string_view error_message_label = "dense matrix");

/**
 * @brief Invert a dense square matrix with diagnostics, using an SVD fallback for
 * high-condition matrices.
 * @ingroup FE_Math
 * @param matrix Row-major n * n matrix (consumed).
 * @param n Matrix dimension.
 * @param error_message_label Prefix for the message of any exception thrown.
 * @return Inverse plus diagnostics and whether the SVD fallback was used.
 * @throws FEException If the size is inconsistent or the matrix is rank-deficient.
 */
[[nodiscard]] DenseInverseResult invert_dense_matrix_with_diagnostics(
    std::vector<double> matrix,
    std::size_t n,
    std::string_view error_message_label = "dense matrix");

/**
 * @brief Validate that a dense inverse has full rank and acceptable conditioning.
 * @ingroup FE_Math
 * @param result Result from invert_dense_matrix_with_diagnostics().
 * @param expected_rank Required (full) rank.
 * @param error_message_label Prefix for the message of any exception thrown.
 * @param max_condition Largest acceptable condition estimate.
 * @throws FEException If the rank is below expected_rank or the condition exceeds max_condition.
 */
void validate_dense_inverse_diagnostics(
    const DenseInverseResult& result,
    std::size_t expected_rank,
    std::string_view error_message_label = "dense matrix",
    double max_condition = dense_matrix_condition_error_threshold());

/**
 * @brief Numerical rank of a dense matrix from its singular values.
 * @ingroup FE_Math
 * @param matrix Row-major matrix of size rows * cols (consumed).
 * @param rows Row count.
 * @param cols Column count.
 * @return Number of singular values above the scale-aware tolerance.
 * @throws FEException If the matrix size is inconsistent.
 */
[[nodiscard]] std::size_t dense_matrix_rank(std::vector<double> matrix,
                                            std::size_t rows,
                                            std::size_t cols);

/**
 * @brief Moore-Penrose pseudo-inverse via a rank-revealing SVD.
 * @ingroup FE_Math
 * @param matrix Row-major matrix of size rows * cols.
 * @param rows Row count.
 * @param cols Column count.
 * @param error_message_label Prefix for the message of any exception thrown.
 * @return Row-major pseudo-inverse (cols * rows) plus rank/tolerance diagnostics.
 * @throws FEException If the matrix size is inconsistent or the matrix is empty.
 */
[[nodiscard]] DensePseudoInverseResult rank_revealing_pseudo_inverse(
    std::span<const double> matrix,
    std::size_t rows,
    std::size_t cols,
    std::string_view error_message_label = "dense matrix");

} // namespace svmp::FE::math

#endif // SVMP_FE_MATH_DENSELINEARALGEBRA_H
