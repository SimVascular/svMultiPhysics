// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_MATH_VECTOR_H
#define SVMP_FE_MATH_VECTOR_H

/**
 * @file Vector.h
 * @brief Fixed-size vectors with expression templates for FE computations
 *
 * This header provides optimized fixed-size vector operations for element-level
 * computations. All operations use expression templates to eliminate temporaries
 * and are header-only for maximum inlining. Memory is aligned for SIMD operations.
 */

#include "VectorExpr.h"
#include "../Common/Types.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <type_traits>

namespace svmp {
namespace FE {
namespace math {

template<typename T>
inline constexpr T tolerance =
    std::is_floating_point_v<T> ? T(1000) * std::numeric_limits<T>::epsilon() : T(0);

template<typename T>
inline bool approx_zero(T value, T tol = tolerance<T>) {
    using std::abs;
    return abs(value) <= tol;
}

template<typename T>
inline bool approx_equal(T a, T b, T tol = tolerance<T>) {
    using std::abs;
    const T scale = std::max({abs(a), abs(b), T(1)});
    return abs(a - b) <= tol * scale;
}

/**
 * @brief Fixed-size vector for element-level computations
 * @tparam T Scalar type (float, double)
 * @tparam N Vector dimension
 *
 * This class provides small vector operations optimized for
 * compile-time known dimensions. Memory is aligned for SIMD operations.
 */
template<typename T, std::size_t N>
class Vector : public VectorExpr<Vector<T, N>> {
    static_assert(std::is_arithmetic_v<T>, "T must be an arithmetic type");
    static_assert(N > 0, "Vector dimension must be positive");

private:
    alignas(kFEFixedObjectAlignmentBytes) T data_[N];  // SIMD-friendly alignment

public:
    // Type definitions
    using value_type = T;
    using size_type = std::size_t;
    using reference = T&;
    using const_reference = const T&;
    using pointer = T*;
    using const_pointer = const T*;

    /**
     * @brief Default constructor - zero initializes all components
     */
    constexpr Vector() : data_{} {}

    /**
     * @brief Fill constructor - initializes all components with same value
     * @param value Value to fill vector with
     */
    constexpr explicit Vector(T value) {
        for (size_type i = 0; i < N; ++i) {
            data_[i] = value;
        }
    }

    /**
     * @brief Initializer list constructor
     * @param init List of values
     */
    constexpr Vector(std::initializer_list<T> init) : data_{} {
        auto it = init.begin();
        for (size_type i = 0; i < N && it != init.end(); ++i, ++it) {
            data_[i] = *it;
        }
    }

    /**
     * @brief Constructor from expression template
     * @tparam Expr Expression type
     * @param expr Vector expression to evaluate
     */
    template<typename Expr>
    Vector(const VectorExpr<Expr>& expr) {
        const auto& e = expr.derived();
        for (size_type i = 0; i < N; ++i) {
            data_[i] = e[i];
        }
    }

    /**
     * @brief Copy constructor
     */
    constexpr Vector(const Vector&) = default;

    /**
     * @brief Move constructor
     */
    constexpr Vector(Vector&&) noexcept = default;

    /**
     * @brief Copy assignment
     */
    Vector& operator=(const Vector&) = default;

    /**
     * @brief Move assignment
     */
    Vector& operator=(Vector&&) noexcept = default;

    /**
     * @brief Assignment from expression template
     * @tparam Expr Expression type
     * @param expr Vector expression to evaluate
     * @return Reference to this
     */
    template<typename Expr>
    Vector& operator=(const VectorExpr<Expr>& expr) {
        const auto& e = expr.derived();
        for (size_type i = 0; i < N; ++i) {
            data_[i] = e[i];
        }
        return *this;
    }

    /**
     * @brief Get vector size (compile-time constant)
     * @return Number of elements
     */
    static constexpr size_type size() { return N; }

    /**
     * @brief Element access (no bounds checking)
     * @param i Element index
     * @return Reference to element
     */
    constexpr T& operator[](size_type i) {
        return data_[i];
    }

    /**
     * @brief Element access (no bounds checking) - const version
     * @param i Element index
     * @return Const reference to element
     */
    constexpr const T& operator[](size_type i) const {
        return data_[i];
    }

    /**
     * @brief Element access with bounds checking
     * @param i Element index
     * @return Reference to element
     * @throws std::out_of_range if i >= N
     */
    T& at(size_type i) {
        if (i >= N) {
            throw std::out_of_range("Vector::at: index out of range");
        }
        return data_[i];
    }

    /**
     * @brief Element access with bounds checking - const version
     * @param i Element index
     * @return Const reference to element
     * @throws std::out_of_range if i >= N
     */
    const T& at(size_type i) const {
        if (i >= N) {
            throw std::out_of_range("Vector::at: index out of range");
        }
        return data_[i];
    }

    /**
     * @brief Access first element
     * @return Reference to first element
     */
    T& front() { return data_[0]; }
    const T& front() const { return data_[0]; }

    /**
     * @brief Access last element
     * @return Reference to last element
     */
    T& back() { return data_[N-1]; }
    const T& back() const { return data_[N-1]; }

    /**
     * @brief Get pointer to underlying data
     * @return Pointer to first element
     */
    T* data() { return data_; }
    const T* data() const { return data_; }

    /**
     * @brief Fill vector with value
     * @param value Value to fill with
     */
    void fill(T value) {
        for (size_type i = 0; i < N; ++i) {
            data_[i] = value;
        }
    }

    /**
     * @brief Set all components to zero
     */
    void set_zero() {
        fill(T{0});
    }

    // Arithmetic operators

    /**
     * @brief In-place addition
     * @param other Vector to add
     * @return Reference to this
     */
    Vector& operator+=(const Vector& other) {
        for (size_type i = 0; i < N; ++i) {
            data_[i] += other.data_[i];
        }
        return *this;
    }

    /**
     * @brief In-place subtraction
     * @param other Vector to subtract
     * @return Reference to this
     */
    Vector& operator-=(const Vector& other) {
        for (size_type i = 0; i < N; ++i) {
            data_[i] -= other.data_[i];
        }
        return *this;
    }

    /**
     * @brief In-place scalar multiplication
     * @param scalar Scalar to multiply by
     * @return Reference to this
     */
    Vector& operator*=(T scalar) {
        for (size_type i = 0; i < N; ++i) {
            data_[i] *= scalar;
        }
        return *this;
    }

    /**
     * @brief In-place scalar division
     * @param scalar Scalar to divide by
     * @return Reference to this
     */
    Vector& operator/=(T scalar) {
        const T inv = T(1) / scalar;
        return (*this) *= inv;
    }

    // Vector operations

    /**
     * @brief Compute dot product
     * @param other Other vector
     * @return Dot product
     */
    T dot(const Vector& other) const {
        T result = T(0);
        for (size_type i = 0; i < N; ++i) {
            result += data_[i] * other.data_[i];
        }
        return result;
    }

    /**
     * @brief Compute squared Euclidean norm
     * @return Squared norm
     */
    T norm_squared() const {
        return dot(*this);
    }

    /**
     * @brief Compute Euclidean norm
     * @return Norm
     */
    T norm() const {
        using std::sqrt;
        return sqrt(norm_squared());
    }

    /**
     * @brief Get normalized vector
     * @return Unit vector in same direction
     */
    Vector normalized() const {
        const T n = norm();
        if (approx_zero(n)) {
            return Vector();  // Return zero vector
        }
        return (*this) / n;
    }

    /**
     * @brief Normalize this vector in place
     * @return Reference to this
     */
    Vector& normalize() {
        const T n = norm();
        if (!approx_zero(n)) {
            (*this) /= n;
        }
        return *this;
    }

    /**
     * @brief Compute L1 norm (Manhattan norm)
     * @return Sum of absolute values
     */
    T norm_l1() const {
        T result = T(0);
        for (size_type i = 0; i < N; ++i) {
            using std::abs;
            result += abs(data_[i]);
        }
        return result;
    }

    /**
     * @brief Compute L-infinity norm (maximum norm)
     * @return Maximum absolute value
     */
    T norm_inf() const {
        T result = T(0);
        for (size_type i = 0; i < N; ++i) {
            using std::abs;
            result = std::max(result, abs(data_[i]));
        }
        return result;
    }

    /**
     * @brief Get minimum component
     * @return Minimum value
     */
    T min() const {
        T result = data_[0];
        for (size_type i = 1; i < N; ++i) {
            result = std::min(result, data_[i]);
        }
        return result;
    }

    /**
     * @brief Get maximum component
     * @return Maximum value
     */
    T max() const {
        T result = data_[0];
        for (size_type i = 1; i < N; ++i) {
            result = std::max(result, data_[i]);
        }
        return result;
    }

    /**
     * @brief Get sum of all components
     * @return Sum of components
     */
    T sum() const {
        T result = T(0);
        for (size_type i = 0; i < N; ++i) {
            result += data_[i];
        }
        return result;
    }

    /**
     * @brief Get product of all components
     * @return Product of components
     */
    T product() const {
        T result = data_[0];
        for (size_type i = 1; i < N; ++i) {
            result *= data_[i];
        }
        return result;
    }

    // Static factory functions

    /**
     * @brief Create zero vector
     * @return Vector with all components zero
     */
    static constexpr Vector zeros() {
        return Vector();
    }

    /**
     * @brief Create vector with all components one
     * @return Vector with all components one
     */
    static constexpr Vector ones() {
        return Vector(T(1));
    }

    /**
     * @brief Create unit vector along axis
     * @param axis Axis index (0-based)
     * @return Unit vector
     */
    static Vector unit(size_type axis) {
        Vector v;
        if (axis < N) {
            v[axis] = T(1);
        }
        return v;
    }

    /**
     * @brief Create basis vector (alias for unit)
     * @param i Axis index (0-based)
     * @return Basis vector
     */
    static Vector basis(size_type i) {
        return unit(i);
    }

    /**
     * @brief Create zero vector (alias for zeros)
     * @return Zero vector
     */
    static constexpr Vector zero() {
        return zeros();
    }

    /**
     * @brief Get index of minimum element
     * @return Index of minimum value
     */
    size_type min_index() const {
        size_type idx = 0;
        T min_val = data_[0];
        for (size_type i = 1; i < N; ++i) {
            if (data_[i] < min_val) {
                min_val = data_[i];
                idx = i;
            }
        }
        return idx;
    }

    /**
     * @brief Get index of maximum element
     * @return Index of maximum value
     */
    size_type max_index() const {
        size_type idx = 0;
        T max_val = data_[0];
        for (size_type i = 1; i < N; ++i) {
            if (data_[i] > max_val) {
                max_val = data_[i];
                idx = i;
            }
        }
        return idx;
    }

    /**
     * @brief Compute mean of all components
     * @return Average value
     */
    T mean() const {
        return sum() / static_cast<T>(N);
    }

    /**
     * @brief Cross product for 3D vectors
     * @param other Other vector
     * @return Cross product
     * @note Only available for 3D vectors
     */
    template<typename U = T>
    std::enable_if_t<N == 3, Vector<U, 3>> cross(const Vector<U, 3>& other) const {
        return Vector<U, 3>{
            data_[1] * other[2] - data_[2] * other[1],
            data_[2] * other[0] - data_[0] * other[2],
            data_[0] * other[1] - data_[1] * other[0]
        };
    }

    /**
     * @brief Check if vectors are approximately equal
     * @param other Other vector
     * @param tol Tolerance
     * @return true if equal within tolerance
     */
    bool approx_equal(const Vector& other, T tol = tolerance<T>) const {
        for (size_type i = 0; i < N; ++i) {
            using std::abs;
            if (abs(data_[i] - other.data_[i]) > tol) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Equality comparison
     * @param other Other vector
     * @return true if exactly equal
     */
    bool operator==(const Vector& other) const {
        for (size_type i = 0; i < N; ++i) {
            if (data_[i] != other.data_[i]) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Inequality comparison
     * @param other Other vector
     * @return true if not equal
     */
    bool operator!=(const Vector& other) const {
        return !(*this == other);
    }

    // Iterators
    T* begin() { return data_; }
    T* end() { return data_ + N; }
    const T* begin() const { return data_; }
    const T* end() const { return data_ + N; }
    const T* cbegin() const { return data_; }
    const T* cend() const { return data_ + N; }
};

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

/**
 * @brief 3D Cross product
 * @tparam T Scalar type
 * @param a First vector
 * @param b Second vector
 * @return Cross product a × b
 */
template<typename T>
inline Vector3<T> cross(const Vector3<T>& a, const Vector3<T>& b) {
    return Vector3<T>{
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    };
}

/**
 * @brief 2D Cross product (returns scalar - z component of 3D cross)
 * @tparam T Scalar type
 * @param a First vector
 * @param b Second vector
 * @return Scalar cross product
 */
template<typename T>
inline T cross(const Vector2<T>& a, const Vector2<T>& b) {
    return a[0] * b[1] - a[1] * b[0];
}

/**
 * @brief Triple scalar product (a · (b × c))
 * @tparam T Scalar type
 * @param a First vector
 * @param b Second vector
 * @param c Third vector
 * @return Scalar triple product
 */
template<typename T>
inline T triple_product(const Vector3<T>& a, const Vector3<T>& b, const Vector3<T>& c) {
    return a.dot(cross(b, c));
}

// Free functions for common operations

/**
 * @brief Compute dot product
 */
template<typename T, std::size_t N>
inline T dot(const Vector<T, N>& a, const Vector<T, N>& b) {
    return a.dot(b);
}

/**
 * @brief Compute Euclidean norm
 */
template<typename T, std::size_t N>
inline T norm(const Vector<T, N>& v) {
    return v.norm();
}

/**
 * @brief Compute squared Euclidean norm
 */
template<typename T, std::size_t N>
inline T norm_squared(const Vector<T, N>& v) {
    return v.norm_squared();
}

/**
 * @brief Get normalized vector
 */
template<typename T, std::size_t N>
inline Vector<T, N> normalize(const Vector<T, N>& v) {
    return v.normalized();
}

/**
 * @brief Component-wise absolute value
 */
template<typename T, std::size_t N>
inline Vector<T, N> abs(const Vector<T, N>& v) {
    Vector<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        using std::abs;
        result[i] = abs(v[i]);
    }
    return result;
}

/**
 * @brief Component-wise minimum
 */
template<typename T, std::size_t N>
inline Vector<T, N> min(const Vector<T, N>& a, const Vector<T, N>& b) {
    Vector<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = std::min(a[i], b[i]);
    }
    return result;
}

/**
 * @brief Component-wise maximum
 */
template<typename T, std::size_t N>
inline Vector<T, N> max(const Vector<T, N>& a, const Vector<T, N>& b) {
    Vector<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = std::max(a[i], b[i]);
    }
    return result;
}

/**
 * @brief Component-wise clamp
 */
template<typename T, std::size_t N>
inline Vector<T, N> clamp(const Vector<T, N>& v, const Vector<T, N>& min_v, const Vector<T, N>& max_v) {
    Vector<T, N> result;
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = std::clamp(v[i], min_v[i], max_v[i]);
    }
    return result;
}

/**
 * @brief Linear interpolation between vectors
 * @tparam T Scalar type
 * @tparam N Vector dimension
 * @param t Interpolation parameter [0, 1]
 * @param a Start vector (at t=0)
 * @param b End vector (at t=1)
 * @return Interpolated vector
 */
template<typename T, std::size_t N>
inline Vector<T, N> lerp(T t, const Vector<T, N>& a, const Vector<T, N>& b) {
    return a + t * (b - a);
}

/**
 * @brief Spherical linear interpolation (for unit vectors)
 * @tparam T Scalar type
 * @param t Interpolation parameter [0, 1]
 * @param a Start unit vector
 * @param b End unit vector
 * @return Interpolated unit vector
 */
template<typename T>
inline Vector3<T> slerp(T t, const Vector3<T>& a, const Vector3<T>& b) {
    T cos_angle = a.dot(b);

    // Handle numerical issues
    cos_angle = std::clamp(cos_angle, T(-1), T(1));

    // If vectors are nearly parallel, use linear interpolation
    if (cos_angle > T(0.9995)) {
        return normalize(lerp(t, a, b));
    }

    T angle = std::acos(cos_angle);
    T sin_angle = std::sin(angle);

    T t0 = std::sin((T(1) - t) * angle) / sin_angle;
    T t1 = std::sin(t * angle) / sin_angle;

    return t0 * a + t1 * b;
}

/**
 * @brief Reflect vector about normal
 * @tparam T Scalar type
 * @tparam N Vector dimension
 * @param v Incident vector
 * @param n Normal vector (should be unit)
 * @return Reflected vector
 */
template<typename T, std::size_t N>
inline Vector<T, N> reflect(const Vector<T, N>& v, const Vector<T, N>& n) {
    return v - T(2) * dot(v, n) * n;
}

/**
 * @brief Project vector onto another vector
 * @tparam T Scalar type
 * @tparam N Vector dimension
 * @param v Vector to project
 * @param onto Vector to project onto
 * @return Projection of v onto 'onto'
 */
template<typename T, std::size_t N>
inline Vector<T, N> project(const Vector<T, N>& v, const Vector<T, N>& onto) {
    T denom = onto.norm_squared();
    if (approx_zero(denom)) {
        return Vector<T, N>::zeros();
    }
    return (dot(v, onto) / denom) * onto;
}

/**
 * @brief Get perpendicular component of vector
 * @tparam T Scalar type
 * @tparam N Vector dimension
 * @param v Vector
 * @param direction Direction to remove
 * @return Component of v perpendicular to direction
 */
template<typename T, std::size_t N>
inline Vector<T, N> perpendicular(const Vector<T, N>& v, const Vector<T, N>& direction) {
    return v - project(v, direction);
}

/**
 * @brief Compute angle between two vectors
 * @tparam T Scalar type
 * @tparam N Vector dimension
 * @param a First vector
 * @param b Second vector
 * @return Angle in radians [0, π]
 */
template<typename T, std::size_t N>
inline T angle(const Vector<T, N>& a, const Vector<T, N>& b) {
    T cos_angle = dot(a, b) / (norm(a) * norm(b));
    cos_angle = std::clamp(cos_angle, T(-1), T(1));
    return std::acos(cos_angle);
}

/**
 * @brief Check if two vectors are approximately equal
 * @tparam T Scalar type
 * @tparam N Vector dimension
 * @param a First vector
 * @param b Second vector
 * @param tol Tolerance
 * @return true if vectors are equal within tolerance
 */
template<typename T, std::size_t N>
inline bool approx_equal(const Vector<T, N>& a, const Vector<T, N>& b, T tol = tolerance<T>) {
    for (std::size_t i = 0; i < N; ++i) {
        if (!approx_equal(a[i], b[i], tol)) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Stream output operator
 * @tparam T Scalar type
 * @tparam N Vector dimension
 * @param os Output stream
 * @param v Vector to output
 * @return Reference to output stream
 */
template<typename T, std::size_t N>
inline std::ostream& operator<<(std::ostream& os, const Vector<T, N>& v) {
    os << "[";
    for (std::size_t i = 0; i < N; ++i) {
        if (i > 0) os << ", ";
        os << v[i];
    }
    os << "]";
    return os;
}

} // namespace math
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_MATH_VECTOR_H
