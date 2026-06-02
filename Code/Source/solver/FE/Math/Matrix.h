#ifndef SVMP_FE_MATH_MATRIX_H
#define SVMP_FE_MATH_MATRIX_H

/**
 * @file Matrix.h
 * @brief Fixed-size matrices with expression templates and specializations for FE computations
 *
 * This header provides optimized fixed-size matrix operations for element-level
 * computations. Includes specialized analytical formulas for 2x2 and 3x3 matrices
 * (determinant, inverse using Cramer's rule) and Gauss elimination for larger matrices.
 * All operations use expression templates to eliminate temporaries.
 */

#include "MatrixExpr.h"
#include "Vector.h"
#include "MathConstants.h"
#include "../Common/Alignment.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <ostream>
#include <stdexcept>
#include <type_traits>

namespace svmp {
namespace FE {
namespace math {

/**
 * @brief Fixed-size matrix for element-level computations
 * @tparam T Scalar type (float, double)
 * @tparam M Number of rows
 * @tparam N Number of columns
 *
 * Storage is row-major for cache efficiency. Memory is aligned for SIMD operations.
 * Specializations exist for 2x2, 3x3, 4x4 matrices with analytical algorithms.
 */
template<typename T, std::size_t M, std::size_t N>
class Matrix : public MatrixExpr<Matrix<T, M, N>> {
    static_assert(std::is_arithmetic_v<T>, "T must be an arithmetic type");
    static_assert(M > 0 && N > 0, "Matrix dimensions must be positive");

private:
    alignas(kFEFixedObjectAlignmentBytes) T data_[M * N];  // Row-major, SIMD-friendly storage

    // Helper to compute linear index from (i,j)
    static constexpr std::size_t index(std::size_t i, std::size_t j) {
        return i * N + j;
    }

public:
    // Type definitions
    using value_type = T;
    using size_type = std::size_t;
    using reference = T&;
    using const_reference = const T&;
    using pointer = T*;
    using const_pointer = const T*;

    /**
     * @brief Default constructor - zero initializes all elements
     */
    constexpr Matrix() : data_{} {}

    /**
     * @brief Fill constructor - initializes all elements with same value
     * @param value Value to fill matrix with
     */
    constexpr explicit Matrix(T value) {
        for (size_type i = 0; i < M * N; ++i) {
            data_[i] = value;
        }
    }

    /**
     * @brief Initializer list constructor for row-wise initialization
     * @param init Nested initializer lists {{row0}, {row1}, ...}
     */
    constexpr Matrix(std::initializer_list<std::initializer_list<T>> init) : data_{} {
        size_type row = 0;
        for (auto row_init : init) {
            if (row >= M) break;
            size_type col = 0;
            for (auto val : row_init) {
                if (col >= N) break;
                (*this)(row, col) = val;
                ++col;
            }
            ++row;
        }
    }

    /**
     * @brief Constructor from expression template
     * @tparam Expr Expression type
     * @param expr Matrix expression to evaluate
     */
    template<typename Expr>
    Matrix(const MatrixExpr<Expr>& expr) {
        const auto& e = expr.derived();
        for (size_type i = 0; i < M; ++i) {
            for (size_type j = 0; j < N; ++j) {
                (*this)(i, j) = e(i, j);
            }
        }
    }

    /**
     * @brief Copy constructor
     */
    constexpr Matrix(const Matrix&) = default;

    /**
     * @brief Move constructor
     */
    constexpr Matrix(Matrix&&) noexcept = default;

    /**
     * @brief Copy assignment
     */
    Matrix& operator=(const Matrix&) = default;

    /**
     * @brief Move assignment
     */
    Matrix& operator=(Matrix&&) noexcept = default;

    /**
     * @brief Assignment from expression template
     * @tparam Expr Expression type
     * @param expr Matrix expression to evaluate
     * @return Reference to this
     */
    template<typename Expr>
    Matrix& operator=(const MatrixExpr<Expr>& expr) {
        const auto& e = expr.derived();
        for (size_type i = 0; i < M; ++i) {
            for (size_type j = 0; j < N; ++j) {
                (*this)(i, j) = e(i, j);
            }
        }
        return *this;
    }

    /**
     * @brief Get number of rows (compile-time constant)
     * @return Number of rows
     */
    static constexpr size_type rows() { return M; }

    /**
     * @brief Get number of columns (compile-time constant)
     * @return Number of columns
     */
    static constexpr size_type cols() { return N; }

    /**
     * @brief Get total number of elements
     * @return M * N
     */
    static constexpr size_type size() { return M * N; }

    /**
     * @brief Element access (no bounds checking)
     * @param i Row index
     * @param j Column index
     * @return Reference to element
     */
    constexpr T& operator()(size_type i, size_type j) {
        return data_[index(i, j)];
    }

    /**
     * @brief Element access (no bounds checking) - const version
     * @param i Row index
     * @param j Column index
     * @return Const reference to element
     */
    constexpr const T& operator()(size_type i, size_type j) const {
        return data_[index(i, j)];
    }

    /**
     * @brief Element access with bounds checking
     * @param i Row index
     * @param j Column index
     * @return Reference to element
     * @throws std::out_of_range if indices are out of bounds
     */
    T& at(size_type i, size_type j) {
        if (i >= M || j >= N) {
            throw std::out_of_range("Matrix::at: index out of range");
        }
        return (*this)(i, j);
    }

    /**
     * @brief Element access with bounds checking - const version
     * @param i Row index
     * @param j Column index
     * @return Const reference to element
     * @throws std::out_of_range if indices are out of bounds
     */
    const T& at(size_type i, size_type j) const {
        if (i >= M || j >= N) {
            throw std::out_of_range("Matrix::at: index out of range");
        }
        return (*this)(i, j);
    }

    /**
     * @brief Get row as vector
     * @param i Row index
     * @return Vector containing row elements
     */
    Vector<T, N> row(size_type i) const {
        Vector<T, N> result;
        for (size_type j = 0; j < N; ++j) {
            result[j] = (*this)(i, j);
        }
        return result;
    }

    /**
     * @brief Get column as vector
     * @param j Column index
     * @return Vector containing column elements
     */
    Vector<T, M> column(size_type j) const {
        Vector<T, M> result;
        for (size_type i = 0; i < M; ++i) {
            result[i] = (*this)(i, j);
        }
        return result;
    }

    /**
     * @brief Get column as vector (alias for column)
     * @param j Column index
     * @return Vector containing column elements
     */
    Vector<T, M> col(size_type j) const {
        return column(j);
    }

    /**
     * @brief Set row from vector
     * @param i Row index
     * @param v Vector of values
     */
    void set_row(size_type i, const Vector<T, N>& v) {
        for (size_type j = 0; j < N; ++j) {
            (*this)(i, j) = v[j];
        }
    }

    /**
     * @brief Set column from vector
     * @param j Column index
     * @param v Vector of values
     */
    void set_column(size_type j, const Vector<T, M>& v) {
        for (size_type i = 0; i < M; ++i) {
            (*this)(i, j) = v[i];
        }
    }

    /**
     * @brief Set column from vector (alias for set_column)
     * @param j Column index
     * @param v Vector of values
     */
    void set_col(size_type j, const Vector<T, M>& v) {
        set_column(j, v);
    }

    /**
     * @brief Get pointer to underlying data
     * @return Pointer to first element
     */
    T* data() { return data_; }
    const T* data() const { return data_; }

    /**
     * @brief Fill matrix with value
     * @param value Value to fill with
     */
    void fill(T value) {
        for (size_type i = 0; i < M * N; ++i) {
            data_[i] = value;
        }
    }

    /**
     * @brief Set all elements to zero
     */
    void set_zero() {
        fill(T{0});
    }

    // Arithmetic operators

    /**
     * @brief In-place addition
     * @param other Matrix to add
     * @return Reference to this
     */
    Matrix& operator+=(const Matrix& other) {
        for (size_type i = 0; i < M * N; ++i) {
            data_[i] += other.data_[i];
        }
        return *this;
    }

    /**
     * @brief In-place subtraction
     * @param other Matrix to subtract
     * @return Reference to this
     */
    Matrix& operator-=(const Matrix& other) {
        for (size_type i = 0; i < M * N; ++i) {
            data_[i] -= other.data_[i];
        }
        return *this;
    }

    /**
     * @brief In-place scalar multiplication
     * @param scalar Scalar to multiply by
     * @return Reference to this
     */
    Matrix& operator*=(T scalar) {
        for (size_type i = 0; i < M * N; ++i) {
            data_[i] *= scalar;
        }
        return *this;
    }

    /**
     * @brief In-place scalar division
     * @param scalar Scalar to divide by
     * @return Reference to this
     */
    Matrix& operator/=(T scalar) {
        const T inv = T(1) / scalar;
        return (*this) *= inv;
    }

    // Matrix operations

    /**
     * @brief Compute transpose
     * @return Transposed matrix
     */
    Matrix<T, N, M> transpose() const {
        Matrix<T, N, M> result;
        for (size_type i = 0; i < M; ++i) {
            for (size_type j = 0; j < N; ++j) {
                result(j, i) = (*this)(i, j);
            }
        }
        return result;
    }

    /**
     * @brief Compute trace (sum of diagonal elements)
     * @return Trace (only valid for square matrices)
     */
    template<std::size_t M2 = M, std::size_t N2 = N>
    std::enable_if_t<M2 == N2, T> trace() const {
        T result = T(0);
        for (size_type i = 0; i < M; ++i) {
            result += (*this)(i, i);
        }
        return result;
    }

    /**
     * @brief Compute Frobenius norm squared
     * @return Sum of squares of all elements
     */
    T frobenius_norm_squared() const {
        T result = T(0);
        for (size_type i = 0; i < M * N; ++i) {
            result += data_[i] * data_[i];
        }
        return result;
    }

    /**
     * @brief Compute Frobenius norm
     * @return Square root of sum of squares
     */
    T frobenius_norm() const {
        using std::sqrt;
        return sqrt(frobenius_norm_squared());
    }

    /**
     * @brief Compute infinity norm (maximum absolute row sum)
     * @return Infinity norm
     */
    T infinity_norm() const {
        T max_row_sum = T(0);
        for (size_type i = 0; i < M; ++i) {
            T row_sum = T(0);
            for (size_type j = 0; j < N; ++j) {
                using std::abs;
                row_sum += abs((*this)(i, j));
            }
            max_row_sum = std::max(max_row_sum, row_sum);
        }
        return max_row_sum;
    }

    /**
     * @brief Compute one norm (maximum absolute column sum)
     * @return One norm
     */
    T one_norm() const {
        T max_col_sum = T(0);
        for (size_type j = 0; j < N; ++j) {
            T col_sum = T(0);
            for (size_type i = 0; i < M; ++i) {
                using std::abs;
                col_sum += abs((*this)(i, j));
            }
            max_col_sum = std::max(max_col_sum, col_sum);
        }
        return max_col_sum;
    }

    /**
     * @brief Get minimum element
     * @return Minimum value
     */
    T min() const {
        return *std::min_element(data_, data_ + M * N);
    }

    /**
     * @brief Get maximum element
     * @return Maximum value
     */
    T max() const {
        return *std::max_element(data_, data_ + M * N);
    }

    /**
     * @brief Get sum of all elements
     * @return Sum of elements
     */
    T sum() const {
        T result = T(0);
        for (size_type i = 0; i < M * N; ++i) {
            result += data_[i];
        }
        return result;
    }

    // Static factory functions

    /**
     * @brief Create zero matrix
     * @return Matrix with all elements zero
     */
    static constexpr Matrix zeros() {
        return Matrix();
    }

    /**
     * @brief Create matrix with all elements one
     * @return Matrix with all elements one
     */
    static constexpr Matrix ones() {
        return Matrix(T(1));
    }

    /**
     * @brief Create identity matrix (only for square matrices)
     * @return Identity matrix
     */
    template<std::size_t M2 = M, std::size_t N2 = N>
    static std::enable_if_t<M2 == N2, Matrix> identity() {
        Matrix result;
        for (size_type i = 0; i < M; ++i) {
            result(i, i) = T(1);
        }
        return result;
    }

    /**
     * @brief Create diagonal matrix from vector (only for square matrices)
     * @param diag Vector of diagonal elements
     * @return Diagonal matrix
     */
    template<std::size_t M2 = M, std::size_t N2 = N>
    static std::enable_if_t<M2 == N2, Matrix> diagonal(const Vector<T, M>& diag) {
        Matrix result;
        for (size_type i = 0; i < M; ++i) {
            result(i, i) = diag[i];
        }
        return result;
    }

    /**
     * @brief Create zero matrix (static factory)
     * @return Zero matrix
     */
    static Matrix zero() {
        return zeros();
    }

    // Property checking methods

    /**
     * @brief Check if matrix is symmetric (only for square matrices)
     * @param tol Tolerance for comparison
     * @return true if symmetric
     */
    template<std::size_t M2 = M, std::size_t N2 = N>
    std::enable_if_t<M2 == N2, bool> is_symmetric(T tol = tolerance<T>) const {
        for (size_type i = 0; i < M; ++i) {
            for (size_type j = i + 1; j < N; ++j) {
                using std::abs;
                if (abs((*this)(i, j) - (*this)(j, i)) > tol) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * @brief Check if matrix is skew-symmetric (only for square matrices)
     * @param tol Tolerance for comparison
     * @return true if skew-symmetric
     */
    template<std::size_t M2 = M, std::size_t N2 = N>
    std::enable_if_t<M2 == N2, bool> is_skew_symmetric(T tol = tolerance<T>) const {
        for (size_type i = 0; i < M; ++i) {
            // Diagonal must be zero
            using std::abs;
            if (abs((*this)(i, i)) > tol) {
                return false;
            }
            for (size_type j = i + 1; j < N; ++j) {
                if (abs((*this)(i, j) + (*this)(j, i)) > tol) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * @brief Check if matrix is diagonal (only for square matrices)
     * @param tol Tolerance for comparison
     * @return true if diagonal
     */
    template<std::size_t M2 = M, std::size_t N2 = N>
    std::enable_if_t<M2 == N2, bool> is_diagonal(T tol = tolerance<T>) const {
        for (size_type i = 0; i < M; ++i) {
            for (size_type j = 0; j < N; ++j) {
                if (i != j) {
                    using std::abs;
                    if (abs((*this)(i, j)) > tol) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    // Determinant (general template, specialized for 2x2, 3x3)
    /**
     * @brief Compute determinant (only for square matrices)
     * @return Determinant value
     */
    template<std::size_t M2 = M, std::size_t N2 = N>
    std::enable_if_t<M2 == N2 && M2 != 2 && M2 != 3, T> determinant() const {
        // For 4x4 and larger, use LU decomposition
        return determinant_lu();
    }

    // Inverse (general template, specialized for 2x2, 3x3)
    /**
     * @brief Compute matrix inverse (only for square matrices)
     * @return Inverse matrix
     */
    template<std::size_t M2 = M, std::size_t N2 = N>
    std::enable_if_t<M2 == N2 && M2 != 2 && M2 != 3, Matrix> inverse() const {
        // For 4x4 and larger, use Gauss-Jordan elimination
        return inverse_gauss_jordan();
    }

private:
    // LU decomposition for determinant (4x4 and larger)
    T determinant_lu() const {
        Matrix<T, M, M> lu = *this;
        T det = T(1);

        for (size_type k = 0; k < M - 1; ++k) {
            // Find pivot
            size_type pivot = k;
            T max_val = std::abs(lu(k, k));
            for (size_type i = k + 1; i < M; ++i) {
                T val = std::abs(lu(i, k));
                if (val > max_val) {
                    max_val = val;
                    pivot = i;
                }
            }

            // Swap rows if needed
            if (pivot != k) {
                for (size_type j = 0; j < M; ++j) {
                    std::swap(lu(k, j), lu(pivot, j));
                }
                det = -det;  // Row swap changes sign
            }

            // Check for singularity
            if (approx_zero(lu(k, k))) {
                return T(0);
            }

            // Eliminate column
            for (size_type i = k + 1; i < M; ++i) {
                T factor = lu(i, k) / lu(k, k);
                for (size_type j = k + 1; j < M; ++j) {
                    lu(i, j) -= factor * lu(k, j);
                }
            }

            det *= lu(k, k);
        }
        det *= lu(M - 1, M - 1);

        return det;
    }

    // Gauss-Jordan elimination for inverse (4x4 and larger)
    Matrix inverse_gauss_jordan() const {
        Matrix<T, M, M> aug;  // Augmented matrix [A | I]
        Matrix<T, M, M> result = Matrix::identity();

        // Copy this matrix to augmented matrix
        for (size_type i = 0; i < M; ++i) {
            for (size_type j = 0; j < M; ++j) {
                aug(i, j) = (*this)(i, j);
            }
        }

        // Forward elimination with partial pivoting
        for (size_type k = 0; k < M; ++k) {
            // Find pivot
            size_type pivot = k;
            T max_val = std::abs(aug(k, k));
            for (size_type i = k + 1; i < M; ++i) {
                T val = std::abs(aug(i, k));
                if (val > max_val) {
                    max_val = val;
                    pivot = i;
                }
            }

            // Swap rows
            if (pivot != k) {
                for (size_type j = 0; j < M; ++j) {
                    std::swap(aug(k, j), aug(pivot, j));
                    std::swap(result(k, j), result(pivot, j));
                }
            }

            // Check for singularity
            if (approx_zero(aug(k, k))) {
                throw std::runtime_error("Matrix is singular");
            }

            // Scale pivot row
            T pivot_val = aug(k, k);
            for (size_type j = 0; j < M; ++j) {
                aug(k, j) /= pivot_val;
                result(k, j) /= pivot_val;
            }

            // Eliminate column
            for (size_type i = 0; i < M; ++i) {
                if (i != k) {
                    T factor = aug(i, k);
                    for (size_type j = 0; j < M; ++j) {
                        aug(i, j) -= factor * aug(k, j);
                        result(i, j) -= factor * result(k, j);
                    }
                }
            }
        }

        return result;
    }

    // Iterators
public:
    T* begin() { return data_; }
    T* end() { return data_ + M * N; }
    const T* begin() const { return data_; }
    const T* end() const { return data_ + M * N; }
    const T* cbegin() const { return data_; }
    const T* cend() const { return data_ + M * N; }
};

// Specialization for 2x2 determinant (analytical formula)
template<typename T>
inline T determinant_2x2(const Matrix<T, 2, 2>& m) {
    return m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0);
}

// Specialization for 2x2 inverse (Cramer's rule)
template<typename T>
inline Matrix<T, 2, 2> inverse_2x2(const Matrix<T, 2, 2>& m) {
    T det = determinant_2x2(m);
    if (approx_zero(det)) {
        throw std::runtime_error("Matrix is singular");
    }

    T inv_det = T(1) / det;
    return Matrix<T, 2, 2>{
        { m(1, 1) * inv_det, -m(0, 1) * inv_det},
        {-m(1, 0) * inv_det,  m(0, 0) * inv_det}
    };
}

// Specialization for 3x3 determinant (Sarrus rule)
template<typename T>
inline T determinant_3x3(const Matrix<T, 3, 3>& m) {
    return m(0, 0) * (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1))
         - m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0))
         + m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
}

// Specialization for 3x3 inverse (Cramer's rule / adjugate method)
template<typename T>
inline Matrix<T, 3, 3> inverse_3x3(const Matrix<T, 3, 3>& m) {
    T det = determinant_3x3(m);
    if (approx_zero(det)) {
        throw std::runtime_error("Matrix is singular");
    }

    T inv_det = T(1) / det;

    // Compute adjugate matrix (transpose of cofactor matrix)
    Matrix<T, 3, 3> adj;
    adj(0, 0) =  (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1));
    adj(0, 1) = -(m(0, 1) * m(2, 2) - m(0, 2) * m(2, 1));
    adj(0, 2) =  (m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1));

    adj(1, 0) = -(m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0));
    adj(1, 1) =  (m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0));
    adj(1, 2) = -(m(0, 0) * m(1, 2) - m(0, 2) * m(1, 0));

    adj(2, 0) =  (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
    adj(2, 1) = -(m(0, 0) * m(2, 1) - m(0, 1) * m(2, 0));
    adj(2, 2) =  (m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0));

    return adj * inv_det;
}

// Template specializations for 2x2 Matrix determinant and inverse
template<typename T>
class Matrix<T, 2, 2> : public MatrixExpr<Matrix<T, 2, 2>> {
    static constexpr std::size_t M = 2;
    static constexpr std::size_t N = 2;

private:
    alignas(kFEFixedObjectAlignmentBytes) T data_[4];

    static constexpr std::size_t index(std::size_t i, std::size_t j) {
        return i * 2 + j;
    }

public:
    using value_type = T;
    using size_type = std::size_t;

    // Include all the same constructors and methods as the general template
    constexpr Matrix() : data_{} {}
    constexpr explicit Matrix(T value) {
        for (size_type i = 0; i < 4; ++i) {
            data_[i] = value;
        }
    }
    constexpr Matrix(std::initializer_list<std::initializer_list<T>> init) : data_{} {
        size_type row = 0;
        for (auto row_init : init) {
            if (row >= 2) break;
            size_type col = 0;
            for (auto val : row_init) {
                if (col >= 2) break;
                (*this)(row, col) = val;
                ++col;
            }
            ++row;
        }
    }

    template<typename Expr>
    Matrix(const MatrixExpr<Expr>& expr) {
        const auto& e = expr.derived();
        for (size_type i = 0; i < 2; ++i) {
            for (size_type j = 0; j < 2; ++j) {
                (*this)(i, j) = e(i, j);
            }
        }
    }

    constexpr Matrix(const Matrix&) = default;
    constexpr Matrix(Matrix&&) noexcept = default;
    Matrix& operator=(const Matrix&) = default;
    Matrix& operator=(Matrix&&) noexcept = default;

    template<typename Expr>
    Matrix& operator=(const MatrixExpr<Expr>& expr) {
        const auto& e = expr.derived();
        for (size_type i = 0; i < 2; ++i) {
            for (size_type j = 0; j < 2; ++j) {
                (*this)(i, j) = e(i, j);
            }
        }
        return *this;
    }

    static constexpr size_type rows() { return 2; }
    static constexpr size_type cols() { return 2; }
    static constexpr size_type size() { return 4; }

    constexpr T& operator()(size_type i, size_type j) {
        return data_[index(i, j)];
    }
    constexpr const T& operator()(size_type i, size_type j) const {
        return data_[index(i, j)];
    }

    T* data() { return data_; }
    const T* data() const { return data_; }

    void fill(T value) {
        for (size_type i = 0; i < 4; ++i) {
            data_[i] = value;
        }
    }

    void set_zero() { fill(T{0}); }

    void set_row(size_type i, const Vector<T, 2>& v) {
        for (size_type j = 0; j < 2; ++j) {
            (*this)(i, j) = v[j];
        }
    }

    void set_column(size_type j, const Vector<T, 2>& v) {
        for (size_type i = 0; i < 2; ++i) {
            (*this)(i, j) = v[i];
        }
    }

    void set_col(size_type j, const Vector<T, 2>& v) {
        set_column(j, v);
    }

    Vector<T, 2> col(size_type j) const {
        return column(j);
    }

    static Matrix zero() {
        return zeros();
    }

    static Matrix diagonal(const Vector<T, 2>& diag) {
        Matrix result;
        result(0, 0) = diag[0];
        result(1, 1) = diag[1];
        return result;
    }

    bool is_symmetric(T tol = tolerance<T>) const {
        using std::abs;
        return abs((*this)(0, 1) - (*this)(1, 0)) <= tol;
    }

    bool is_skew_symmetric(T tol = tolerance<T>) const {
        using std::abs;
        // Diagonal must be zero
        if (abs((*this)(0, 0)) > tol || abs((*this)(1, 1)) > tol) {
            return false;
        }
        // Off-diagonal must be opposite
        return abs((*this)(0, 1) + (*this)(1, 0)) <= tol;
    }

    bool is_diagonal(T tol = tolerance<T>) const {
        using std::abs;
        return abs((*this)(0, 1)) <= tol && abs((*this)(1, 0)) <= tol;
    }

    T frobenius_norm() const {
        using std::sqrt;
        T sum = T(0);
        for (size_type i = 0; i < 4; ++i) {
            sum += data_[i] * data_[i];
        }
        return sqrt(sum);
    }

    T infinity_norm() const {
        using std::abs;
        T row0 = abs((*this)(0, 0)) + abs((*this)(0, 1));
        T row1 = abs((*this)(1, 0)) + abs((*this)(1, 1));
        return std::max(row0, row1);
    }

    T one_norm() const {
        using std::abs;
        T col0 = abs((*this)(0, 0)) + abs((*this)(1, 0));
        T col1 = abs((*this)(0, 1)) + abs((*this)(1, 1));
        return std::max(col0, col1);
    }

    Matrix& operator+=(const Matrix& other) {
        for (size_type i = 0; i < 4; ++i) {
            data_[i] += other.data_[i];
        }
        return *this;
    }

    Matrix& operator-=(const Matrix& other) {
        for (size_type i = 0; i < 4; ++i) {
            data_[i] -= other.data_[i];
        }
        return *this;
    }

    Matrix& operator*=(T scalar) {
        for (size_type i = 0; i < 4; ++i) {
            data_[i] *= scalar;
        }
        return *this;
    }

    Matrix& operator/=(T scalar) {
        const T inv = T(1) / scalar;
        return (*this) *= inv;
    }

    Matrix<T, 2, 2> transpose() const {
        return Matrix<T, 2, 2>{
            {(*this)(0, 0), (*this)(1, 0)},
            {(*this)(0, 1), (*this)(1, 1)}
        };
    }

    T trace() const {
        return (*this)(0, 0) + (*this)(1, 1);
    }

    static Matrix identity() {
        Matrix result;
        result(0, 0) = T(1);
        result(1, 1) = T(1);
        return result;
    }

    static Matrix zeros() {
        return Matrix();
    }

    static Matrix ones() {
        return Matrix(T(1));
    }

    // Specialized 2x2 determinant
    T determinant() const {
        return determinant_2x2(*this);
    }

    // Specialized 2x2 inverse
    Matrix inverse() const {
        return inverse_2x2(*this);
    }

    Vector<T, 2> row(size_type i) const {
        return Vector<T, 2>{(*this)(i, 0), (*this)(i, 1)};
    }

    Vector<T, 2> column(size_type j) const {
        return Vector<T, 2>{(*this)(0, j), (*this)(1, j)};
    }

    T* begin() { return data_; }
    T* end() { return data_ + 4; }
    const T* begin() const { return data_; }
    const T* end() const { return data_ + 4; }
};

// Template specialization for 3x3 Matrix
template<typename T>
class Matrix<T, 3, 3> : public MatrixExpr<Matrix<T, 3, 3>> {
    static constexpr std::size_t M = 3;
    static constexpr std::size_t N = 3;

private:
    alignas(kFEFixedObjectAlignmentBytes) T data_[9];

    static constexpr std::size_t index(std::size_t i, std::size_t j) {
        return i * 3 + j;
    }

public:
    using value_type = T;
    using size_type = std::size_t;

    constexpr Matrix() : data_{} {}
    constexpr explicit Matrix(T value) {
        for (size_type i = 0; i < 9; ++i) {
            data_[i] = value;
        }
    }
    constexpr Matrix(std::initializer_list<std::initializer_list<T>> init) : data_{} {
        size_type row = 0;
        for (auto row_init : init) {
            if (row >= 3) break;
            size_type col = 0;
            for (auto val : row_init) {
                if (col >= 3) break;
                (*this)(row, col) = val;
                ++col;
            }
            ++row;
        }
    }

    template<typename Expr>
    Matrix(const MatrixExpr<Expr>& expr) {
        const auto& e = expr.derived();
        for (size_type i = 0; i < 3; ++i) {
            for (size_type j = 0; j < 3; ++j) {
                (*this)(i, j) = e(i, j);
            }
        }
    }

    constexpr Matrix(const Matrix&) = default;
    constexpr Matrix(Matrix&&) noexcept = default;
    Matrix& operator=(const Matrix&) = default;
    Matrix& operator=(Matrix&&) noexcept = default;

    template<typename Expr>
    Matrix& operator=(const MatrixExpr<Expr>& expr) {
        const auto& e = expr.derived();
        for (size_type i = 0; i < 3; ++i) {
            for (size_type j = 0; j < 3; ++j) {
                (*this)(i, j) = e(i, j);
            }
        }
        return *this;
    }

    static constexpr size_type rows() { return 3; }
    static constexpr size_type cols() { return 3; }
    static constexpr size_type size() { return 9; }

    constexpr T& operator()(size_type i, size_type j) {
        return data_[index(i, j)];
    }
    constexpr const T& operator()(size_type i, size_type j) const {
        return data_[index(i, j)];
    }

    T* data() { return data_; }
    const T* data() const { return data_; }

    void fill(T value) {
        for (size_type i = 0; i < 9; ++i) {
            data_[i] = value;
        }
    }

    void set_zero() { fill(T{0}); }

    void set_row(size_type i, const Vector<T, 3>& v) {
        for (size_type j = 0; j < 3; ++j) {
            (*this)(i, j) = v[j];
        }
    }

    void set_column(size_type j, const Vector<T, 3>& v) {
        for (size_type i = 0; i < 3; ++i) {
            (*this)(i, j) = v[i];
        }
    }

    void set_col(size_type j, const Vector<T, 3>& v) {
        set_column(j, v);
    }

    Vector<T, 3> col(size_type j) const {
        return column(j);
    }

    static Matrix zero() {
        return zeros();
    }

    static Matrix diagonal(const Vector<T, 3>& diag) {
        Matrix result;
        result(0, 0) = diag[0];
        result(1, 1) = diag[1];
        result(2, 2) = diag[2];
        return result;
    }

    bool is_symmetric(T tol = tolerance<T>) const {
        using std::abs;
        for (size_type i = 0; i < 3; ++i) {
            for (size_type j = i + 1; j < 3; ++j) {
                if (abs((*this)(i, j) - (*this)(j, i)) > tol) {
                    return false;
                }
            }
        }
        return true;
    }

    bool is_skew_symmetric(T tol = tolerance<T>) const {
        using std::abs;
        // Diagonal must be zero
        for (size_type i = 0; i < 3; ++i) {
            if (abs((*this)(i, i)) > tol) {
                return false;
            }
        }
        // Off-diagonal must be opposite
        for (size_type i = 0; i < 3; ++i) {
            for (size_type j = i + 1; j < 3; ++j) {
                if (abs((*this)(i, j) + (*this)(j, i)) > tol) {
                    return false;
                }
            }
        }
        return true;
    }

    bool is_diagonal(T tol = tolerance<T>) const {
        using std::abs;
        for (size_type i = 0; i < 3; ++i) {
            for (size_type j = 0; j < 3; ++j) {
                if (i != j && abs((*this)(i, j)) > tol) {
                    return false;
                }
            }
        }
        return true;
    }

    T frobenius_norm() const {
        using std::sqrt;
        T sum = T(0);
        for (size_type i = 0; i < 9; ++i) {
            sum += data_[i] * data_[i];
        }
        return sqrt(sum);
    }

    T infinity_norm() const {
        using std::abs;
        T max_row_sum = T(0);
        for (size_type i = 0; i < 3; ++i) {
            T row_sum = T(0);
            for (size_type j = 0; j < 3; ++j) {
                row_sum += abs((*this)(i, j));
            }
            max_row_sum = std::max(max_row_sum, row_sum);
        }
        return max_row_sum;
    }

    T one_norm() const {
        using std::abs;
        T max_col_sum = T(0);
        for (size_type j = 0; j < 3; ++j) {
            T col_sum = T(0);
            for (size_type i = 0; i < 3; ++i) {
                col_sum += abs((*this)(i, j));
            }
            max_col_sum = std::max(max_col_sum, col_sum);
        }
        return max_col_sum;
    }

    Matrix& operator+=(const Matrix& other) {
        for (size_type i = 0; i < 9; ++i) {
            data_[i] += other.data_[i];
        }
        return *this;
    }

    Matrix& operator-=(const Matrix& other) {
        for (size_type i = 0; i < 9; ++i) {
            data_[i] -= other.data_[i];
        }
        return *this;
    }

    Matrix& operator*=(T scalar) {
        for (size_type i = 0; i < 9; ++i) {
            data_[i] *= scalar;
        }
        return *this;
    }

    Matrix& operator/=(T scalar) {
        const T inv = T(1) / scalar;
        return (*this) *= inv;
    }

    Matrix<T, 3, 3> transpose() const {
        Matrix<T, 3, 3> result;
        for (size_type i = 0; i < 3; ++i) {
            for (size_type j = 0; j < 3; ++j) {
                result(j, i) = (*this)(i, j);
            }
        }
        return result;
    }

    T trace() const {
        return (*this)(0, 0) + (*this)(1, 1) + (*this)(2, 2);
    }

    static Matrix identity() {
        Matrix result;
        result(0, 0) = T(1);
        result(1, 1) = T(1);
        result(2, 2) = T(1);
        return result;
    }

    static Matrix zeros() {
        return Matrix();
    }

    static Matrix ones() {
        return Matrix(T(1));
    }

    // Specialized 3x3 determinant
    T determinant() const {
        return determinant_3x3(*this);
    }

    // Specialized 3x3 inverse
    Matrix inverse() const {
        return inverse_3x3(*this);
    }

    Vector<T, 3> row(size_type i) const {
        return Vector<T, 3>{(*this)(i, 0), (*this)(i, 1), (*this)(i, 2)};
    }

    Vector<T, 3> column(size_type j) const {
        return Vector<T, 3>{(*this)(0, j), (*this)(1, j), (*this)(2, j)};
    }

    T* begin() { return data_; }
    T* end() { return data_ + 9; }
    const T* begin() const { return data_; }
    const T* end() const { return data_ + 9; }
};

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

// Matrix-vector multiplication
template<typename T, std::size_t M, std::size_t N>
inline Vector<T, M> operator*(const Matrix<T, M, N>& A, const Vector<T, N>& x) {
    Vector<T, M> result;
    for (std::size_t i = 0; i < M; ++i) {
        T sum = T(0);
        for (std::size_t j = 0; j < N; ++j) {
            sum += A(i, j) * x[j];
        }
        result[i] = sum;
    }
    return result;
}

// Vector-matrix multiplication (row vector * matrix)
template<typename T, std::size_t M, std::size_t N>
inline Vector<T, N> operator*(const Vector<T, M>& x, const Matrix<T, M, N>& A) {
    Vector<T, N> result;
    for (std::size_t j = 0; j < N; ++j) {
        T sum = T(0);
        for (std::size_t i = 0; i < M; ++i) {
            sum += x[i] * A(i, j);
        }
        result[j] = sum;
    }
    return result;
}

// Matrix-matrix multiplication
template<typename T, std::size_t M, std::size_t N, std::size_t P>
inline Matrix<T, M, P> operator*(const Matrix<T, M, N>& A, const Matrix<T, N, P>& B) {
    Matrix<T, M, P> result;
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t k = 0; k < N; ++k) {
            T a_ik = A(i, k);
            for (std::size_t j = 0; j < P; ++j) {
                result(i, j) += a_ik * B(k, j);
            }
        }
    }
    return result;
}

// Free functions

/**
 * @brief Compute matrix transpose
 */
template<typename T, std::size_t M, std::size_t N>
inline Matrix<T, N, M> transpose(const Matrix<T, M, N>& m) {
    return m.transpose();
}

/**
 * @brief Compute matrix trace
 */
template<typename T, std::size_t N>
inline T trace(const Matrix<T, N, N>& m) {
    return m.trace();
}

/**
 * @brief Compute matrix determinant
 */
template<typename T, std::size_t N>
inline T determinant(const Matrix<T, N, N>& m) {
    return m.determinant();
}

/**
 * @brief Compute matrix inverse
 */
template<typename T, std::size_t N>
inline Matrix<T, N, N> inverse(const Matrix<T, N, N>& m) {
    return m.inverse();
}

/**
 * @brief Compute Frobenius norm
 */
template<typename T, std::size_t M, std::size_t N>
inline T frobenius_norm(const Matrix<T, M, N>& m) {
    return m.frobenius_norm();
}

/**
 * @brief Component-wise absolute value
 */
template<typename T, std::size_t M, std::size_t N>
inline Matrix<T, M, N> abs(const Matrix<T, M, N>& m) {
    Matrix<T, M, N> result;
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            using std::abs;
            result(i, j) = abs(m(i, j));
        }
    }
    return result;
}

/**
 * @brief Component-wise minimum
 */
template<typename T, std::size_t M, std::size_t N>
inline Matrix<T, M, N> min(const Matrix<T, M, N>& a, const Matrix<T, M, N>& b) {
    Matrix<T, M, N> result;
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            result(i, j) = std::min(a(i, j), b(i, j));
        }
    }
    return result;
}

/**
 * @brief Component-wise maximum
 */
template<typename T, std::size_t M, std::size_t N>
inline Matrix<T, M, N> max(const Matrix<T, M, N>& a, const Matrix<T, M, N>& b) {
    Matrix<T, M, N> result;
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            result(i, j) = std::max(a(i, j), b(i, j));
        }
    }
    return result;
}

/**
 * @brief Outer product of two vectors
 */
template<typename T, std::size_t M, std::size_t N>
inline Matrix<T, M, N> outer_product(const Vector<T, M>& u, const Vector<T, N>& v) {
    Matrix<T, M, N> result;
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            result(i, j) = u[i] * v[j];
        }
    }
    return result;
}

/**
 * @brief Check if two matrices are approximately equal
 */
template<typename T, std::size_t M, std::size_t N>
inline bool approx_equal(const Matrix<T, M, N>& a, const Matrix<T, M, N>& b, T tol = tolerance<T>) {
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            if (!approx_equal(a(i, j), b(i, j), tol)) {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief Stream output operator for matrices
 * @tparam T Scalar type
 * @tparam M Number of rows
 * @tparam N Number of columns
 * @param os Output stream
 * @param m Matrix to output
 * @return Reference to output stream
 */
template<typename T, std::size_t M, std::size_t N>
inline std::ostream& operator<<(std::ostream& os, const Matrix<T, M, N>& m) {
    os << "[";
    for (std::size_t i = 0; i < M; ++i) {
        if (i > 0) os << "\n ";
        os << "[";
        for (std::size_t j = 0; j < N; ++j) {
            if (j > 0) os << ", ";
            os << m(i, j);
        }
        os << "]";
    }
    os << "]";
    return os;
}

} // namespace math
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_MATH_MATRIX_H
