// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_COMMON_MATH_VECTOR_H
#define SVMP_COMMON_MATH_VECTOR_H

/**
 * @file Vector.h
 * @brief Fixed-size vector types for FE computations, backed by Eigen.
 *
 * The FE library standardizes on Eigen for linear algebra. These aliases give
 * element-level code a stable vocabulary type without re-exporting all of
 * Eigen. Note that, unlike the previous in-house implementation, Eigen types
 * are NOT zero-initialized by default construction; use Vector::Zero() where a
 * zeroed value is required.
 *
 * This is a small, fixed-size (compile-time length) vector for element-level FE
 * kernels in namespace svmp::math. It is distinct from, and not a replacement
 * for, the legacy dynamically sized container in solver/Vector.h: the two differ
 * in namespace, size model (compile-time vs runtime), and memory management, and
 * coexist deliberately.
 */

#include <Eigen/Core>

#include <cstddef>

/**
 * @defgroup SVMP_Math Math
 * @brief Linear algebra vocabulary types and dense utilities for finite-element computations.
 *
 * @details The Math module defines the fixed-size vector and matrix types
 * used in element-level kernels (as aliases of Eigen types) and dense linear
 * algebra utilities used by basis construction and local transforms.
 *
 * @defgroup SVMP_Math_Vector Vector
 * @ingroup SVMP_Math
 * @brief Fixed-size vector type aliases.
 */

namespace svmp::math {

/**
 * @brief Fixed-size column vector for element-level computations
 * @ingroup SVMP_Math_Vector
 * @tparam T Scalar type (float, double)
 * @tparam N Vector dimension
 */
template<typename T, std::size_t N>
using Vector = Eigen::Matrix<T, static_cast<int>(N), 1>;

} // namespace svmp::math

#endif // SVMP_COMMON_MATH_VECTOR_H
