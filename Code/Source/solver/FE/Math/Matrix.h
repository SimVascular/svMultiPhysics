// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_MATH_MATRIX_H
#define SVMP_FE_MATH_MATRIX_H

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

/// \defgroup FE_MatrixMath Matrix
/// \ingroup FE_Math
/// \brief Fixed-size matrix type aliases.

namespace svmp {
namespace FE {
namespace math {

/**
 * @brief Fixed-size matrix for element-level computations
 * @ingroup FE_MatrixMath
 * @tparam T Scalar type (float, double)
 * @tparam M Number of rows
 * @tparam N Number of columns
 */
template<typename T, std::size_t M, std::size_t N>
using Matrix = Eigen::Matrix<T, static_cast<int>(M), static_cast<int>(N)>;

// Type aliases for common matrix types
template<typename T> using Matrix2x2 = Matrix<T, 2, 2>;
template<typename T> using Matrix3x3 = Matrix<T, 3, 3>;
template<typename T> using Matrix4x4 = Matrix<T, 4, 4>;
template<typename T> using Matrix2x3 = Matrix<T, 2, 3>;
template<typename T> using Matrix3x2 = Matrix<T, 3, 2>;
template<typename T> using Matrix3x4 = Matrix<T, 3, 4>;
template<typename T> using Matrix4x3 = Matrix<T, 4, 3>;

// Double precision aliases
using Matrix2x2d = Matrix2x2<double>;
using Matrix3x3d = Matrix3x3<double>;
using Matrix4x4d = Matrix4x4<double>;

// Single precision aliases
using Matrix2x2f = Matrix2x2<float>;
using Matrix3x3f = Matrix3x3<float>;
using Matrix4x4f = Matrix4x4<float>;

} // namespace math
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_MATH_MATRIX_H
