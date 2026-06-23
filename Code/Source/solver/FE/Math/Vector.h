// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_MATH_VECTOR_H
#define SVMP_FE_MATH_VECTOR_H

/**
 * @file Vector.h
 * @brief Fixed-size vector types for FE computations, backed by Eigen.
 *
 * The FE library standardizes on Eigen for linear algebra. These aliases give
 * element-level code a stable vocabulary type without re-exporting all of
 * Eigen. Note that, unlike the previous in-house implementation, Eigen types
 * are NOT zero-initialized by default construction; use Vector::Zero() where a
 * zeroed value is required.
 */

#include <Eigen/Core>

#include <cstddef>

/**
 * @defgroup FE_Math Math
 * @ingroup FE
 * @brief Linear algebra vocabulary types and dense utilities for finite-element computations.
 *
 * @details The Math module defines the fixed-size vector and matrix types
 * used in element-level kernels (as aliases of Eigen types) and dense linear
 * algebra utilities used by basis construction and local transforms.
 *
 * @defgroup FE_VectorMath Vector
 * @ingroup FE_Math
 * @brief Fixed-size vector type aliases.
 */

namespace svmp {
namespace FE {
namespace math {

/**
 * @brief Fixed-size column vector for element-level computations
 * @ingroup FE_VectorMath
 * @tparam T Scalar type (float, double)
 * @tparam N Vector dimension
 */
template<typename T, std::size_t N>
using Vector = Eigen::Matrix<T, static_cast<int>(N), 1>;

// Type aliases for common vector types
template<typename T> using Vector2 = Vector<T, 2>;
template<typename T> using Vector3 = Vector<T, 3>;
template<typename T> using Vector4 = Vector<T, 4>;

// Double precision aliases
using Vector2d = Vector2<double>;
using Vector3d = Vector3<double>;
using Vector4d = Vector4<double>;

// Single precision aliases
using Vector2f = Vector2<float>;
using Vector3f = Vector3<float>;
using Vector4f = Vector4<float>;

// Integer aliases
using Vector2i = Vector2<int>;
using Vector3i = Vector3<int>;
using Vector4i = Vector4<int>;

} // namespace math
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_MATH_VECTOR_H
