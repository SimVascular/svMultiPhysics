// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_COMMON_MATH_MATRIX_H
#define SVMP_COMMON_MATH_MATRIX_H

/**
 * @file Matrix.h
 * @brief Fixed-size matrix types for FE computations, backed by Eigen.
 *
 * The FE library standardizes on Eigen for linear algebra. These aliases give
 * element-level code a stable vocabulary type without re-exporting all of
 * Eigen. Storage is Eigen's default (column-major); element access through
 * operator()(row, col) is unchanged. Note that, unlike the previous in-house
 * implementation, Eigen types are NOT zero-initialized by default
 * construction; use Matrix::Zero() where a zeroed value is required.
 */

#include "Vector.h"

#include <Eigen/Core>

#include <cstddef>

/**
 * @defgroup SVMP_Math_Matrix Matrix
 * @ingroup SVMP_Math
 * @brief Fixed-size matrix type aliases.
 */

namespace svmp::math {

/**
 * @brief Fixed-size matrix for element-level computations
 * @ingroup SVMP_Math_Matrix
 * @tparam T Scalar type (float, double)
 * @tparam M Number of rows
 * @tparam N Number of columns
 */
template<typename T, std::size_t M, std::size_t N>
using Matrix = Eigen::Matrix<T, static_cast<int>(M), static_cast<int>(N)>;

} // namespace svmp::math

#endif // SVMP_COMMON_MATH_MATRIX_H
