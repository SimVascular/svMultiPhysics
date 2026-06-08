// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_MATH_MATRIX_EXPR_H
#define SVMP_FE_MATH_MATRIX_EXPR_H

/**
 * @file MatrixExpr.h
 * @brief Expression template infrastructure for lazy evaluation of matrix operations
 *
 * This header provides expression templates that enable compound matrix operations
 * without creating temporary objects. Operations are evaluated lazily at the point
 * of assignment, eliminating intermediate allocations and improving performance.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <type_traits>

#include "VectorExpr.h"

namespace svmp {
namespace FE {
namespace math {

/**
 * @brief Base class for all matrix expressions using CRTP
 * @tparam Derived The derived expression type
 *
 * This uses the Curiously Recurring Template Pattern (CRTP) to provide
 * static polymorphism for expression templates.
 */
template<typename Derived>
class MatrixExpr {
public:
    /**
     * @brief Get the derived expression
     * @return Reference to the derived type
     */
    const Derived& derived() const {
        return static_cast<const Derived&>(*this);
    }

    /**
     * @brief Get the derived expression (non-const)
     * @return Reference to the derived type
     */
    Derived& derived() {
        return static_cast<Derived&>(*this);
    }

    /**
     * @brief Access element by row and column indices
     * @param i Row index
     * @param j Column index
     * @return Value at (i,j)
     */
    auto operator()(std::size_t i, std::size_t j) const {
        return derived()(i, j);
    }

    /**
     * @brief Get number of rows
     * @return Number of rows
     */
    std::size_t rows() const {
        return derived().rows();
    }

    /**
     * @brief Get number of columns
     * @return Number of columns
     */
    std::size_t cols() const {
        return derived().cols();
    }
};

/**
 * @brief Binary expression for element-wise operations between two matrix expressions
 * @tparam LHS Left-hand side expression type
 * @tparam RHS Right-hand side expression type
 * @tparam Op Binary operation functor
 */
template<typename LHS, typename RHS, typename Op>
class MatrixBinaryExpr : public MatrixExpr<MatrixBinaryExpr<LHS, RHS, Op>> {
private:
    LHS lhs_;
    RHS rhs_;
    Op op_;

public:
    /**
     * @brief Construct binary expression
     * @param lhs Left operand
     * @param rhs Right operand
     * @param op Operation to apply
     */
    constexpr MatrixBinaryExpr(const LHS& lhs, const RHS& rhs, Op op = Op{})
        : lhs_(lhs), rhs_(rhs), op_(op) {}

    /**
     * @brief Access element at (i,j)
     * @param i Row index
     * @param j Column index
     * @return Result of operation on elements at (i,j)
     */
    constexpr auto operator()(std::size_t i, std::size_t j) const {
        return op_(lhs_(i, j), rhs_(i, j));
    }

    /**
     * @brief Get number of rows
     * @return Number of rows
     */
    constexpr std::size_t rows() const {
        return lhs_.rows();
    }

    /**
     * @brief Get number of columns
     * @return Number of columns
     */
    constexpr std::size_t cols() const {
        return lhs_.cols();
    }
};

/**
 * @brief Unary expression for element-wise operations on a single matrix expression
 * @tparam Expr Expression type
 * @tparam Op Unary operation functor
 */
template<typename Expr, typename Op>
class MatrixUnaryExpr : public MatrixExpr<MatrixUnaryExpr<Expr, Op>> {
private:
    Expr expr_;
    Op op_;

public:
    /**
     * @brief Construct unary expression
     * @param expr Operand expression
     * @param op Operation to apply
     */
    constexpr MatrixUnaryExpr(const Expr& expr, Op op = Op{})
        : expr_(expr), op_(op) {}

    /**
     * @brief Access element at (i,j)
     * @param i Row index
     * @param j Column index
     * @return Result of operation on element at (i,j)
     */
    constexpr auto operator()(std::size_t i, std::size_t j) const {
        return op_(expr_(i, j));
    }

    /**
     * @brief Get number of rows
     * @return Number of rows
     */
    constexpr std::size_t rows() const {
        return expr_.rows();
    }

    /**
     * @brief Get number of columns
     * @return Number of columns
     */
    constexpr std::size_t cols() const {
        return expr_.cols();
    }
};

/**
 * @brief Scalar multiplication expression
 * @tparam Expr Matrix expression type
 * @tparam Scalar Scalar type
 */
template<typename Expr, typename Scalar>
class MatrixScalarExpr : public MatrixExpr<MatrixScalarExpr<Expr, Scalar>> {
private:
    Expr expr_;
    Scalar scalar_;

public:
    /**
     * @brief Construct scalar multiplication expression
     * @param expr Matrix expression
     * @param scalar Scalar value
     */
    constexpr MatrixScalarExpr(const Expr& expr, Scalar scalar)
        : expr_(expr), scalar_(scalar) {}

    /**
     * @brief Access element at (i,j)
     * @param i Row index
     * @param j Column index
     * @return Element multiplied by scalar
     */
    constexpr auto operator()(std::size_t i, std::size_t j) const {
        return expr_(i, j) * scalar_;
    }

    /**
     * @brief Get number of rows
     * @return Number of rows
     */
    constexpr std::size_t rows() const {
        return expr_.rows();
    }

    /**
     * @brief Get number of columns
     * @return Number of columns
     */
    constexpr std::size_t cols() const {
        return expr_.cols();
    }
};

/**
 * @brief Scalar division expression
 * @tparam Expr Matrix expression type
 * @tparam Scalar Scalar type
 */
template<typename Expr, typename Scalar>
class MatrixScalarDivExpr : public MatrixExpr<MatrixScalarDivExpr<Expr, Scalar>> {
private:
    Expr expr_;
    Scalar scalar_;

public:
    /**
     * @brief Construct scalar division expression
     * @param expr Matrix expression
     * @param scalar Scalar divisor
     */
    constexpr MatrixScalarDivExpr(const Expr& expr, Scalar scalar)
        : expr_(expr), scalar_(scalar) {}

    /**
     * @brief Access element at (i,j)
     * @param i Row index
     * @param j Column index
     * @return Element divided by scalar
     */
    constexpr auto operator()(std::size_t i, std::size_t j) const {
        return expr_(i, j) / scalar_;
    }

    /**
     * @brief Get number of rows
     * @return Number of rows
     */
    constexpr std::size_t rows() const {
        return expr_.rows();
    }

    /**
     * @brief Get number of columns
     * @return Number of columns
     */
    constexpr std::size_t cols() const {
        return expr_.cols();
    }
};

/**
 * @brief Matrix multiplication expression (lazy evaluation)
 * @tparam LHS Left matrix expression type
 * @tparam RHS Right matrix expression type
 *
 * Computes matrix multiplication A*B lazily
 */
template<typename LHS, typename RHS>
class MatrixMulExpr : public MatrixExpr<MatrixMulExpr<LHS, RHS>> {
private:
    LHS lhs_;
    RHS rhs_;

public:
    /**
     * @brief Construct matrix multiplication expression
     * @param lhs Left matrix
     * @param rhs Right matrix
     */
    constexpr MatrixMulExpr(const LHS& lhs, const RHS& rhs)
        : lhs_(lhs), rhs_(rhs) {}

    /**
     * @brief Compute element at (i,j)
     * @param i Row index
     * @param j Column index
     * @return Dot product of row i of lhs and column j of rhs
     */
    constexpr auto operator()(std::size_t i, std::size_t j) const {
        using result_type = decltype(lhs_(0, 0) * rhs_(0, 0));
        result_type sum = result_type{0};
        const auto n = lhs_.cols();
        for (std::size_t k = 0; k < n; ++k) {
            sum += lhs_(i, k) * rhs_(k, j);
        }
        return sum;
    }

    /**
     * @brief Get number of rows (from left matrix)
     * @return Number of rows
     */
    constexpr std::size_t rows() const {
        return lhs_.rows();
    }

    /**
     * @brief Get number of columns (from right matrix)
     * @return Number of columns
     */
    constexpr std::size_t cols() const {
        return rhs_.cols();
    }
};

/**
 * @brief Transpose expression (lazy evaluation)
 * @tparam Expr Matrix expression type
 */
template<typename Expr>
class TransposeExpr : public MatrixExpr<TransposeExpr<Expr>> {
private:
    Expr expr_;

public:
    /**
     * @brief Construct transpose expression
     * @param expr Matrix expression to transpose
     */
    constexpr explicit TransposeExpr(const Expr& expr)
        : expr_(expr) {}

    /**
     * @brief Access transposed element
     * @param i Row index (becomes column in original)
     * @param j Column index (becomes row in original)
     * @return Element at (j,i) of original matrix
     */
    constexpr auto operator()(std::size_t i, std::size_t j) const {
        return expr_(j, i);
    }

    /**
     * @brief Get number of rows (columns of original)
     * @return Number of rows
     */
    constexpr std::size_t rows() const {
        return expr_.cols();
    }

    /**
     * @brief Get number of columns (rows of original)
     * @return Number of columns
     */
    constexpr std::size_t cols() const {
        return expr_.rows();
    }
};

/**
 * @brief Diagonal matrix expression (creates diagonal matrix from vector)
 * @tparam VecExpr Vector expression type
 */
template<typename VecExpr>
class DiagonalExpr : public MatrixExpr<DiagonalExpr<VecExpr>> {
private:
    VecExpr vec_;
    std::size_t n_;

public:
    /**
     * @brief Construct diagonal matrix from vector
     * @param vec Vector of diagonal elements
     * @param n Matrix dimension (default: vector size)
     */
    constexpr explicit DiagonalExpr(const VecExpr& vec, std::size_t n = 0)
        : vec_(vec), n_(n > 0 ? n : vec.size()) {}

    /**
     * @brief Access element
     * @param i Row index
     * @param j Column index
     * @return Diagonal element if i==j, zero otherwise
     */
    constexpr auto operator()(std::size_t i, std::size_t j) const {
        using result_type = decltype(vec_[0]);
        return (i == j && i < vec_.size()) ? vec_[i] : result_type{0};
    }

    /**
     * @brief Get number of rows
     * @return Number of rows
     */
    constexpr std::size_t rows() const {
        return n_;
    }

    /**
     * @brief Get number of columns
     * @return Number of columns
     */
    constexpr std::size_t cols() const {
        return n_;
    }
};

/**
 * @brief Addition operator for matrix expressions
 */
template<typename LHS, typename RHS,
         typename = std::enable_if_t<
             std::is_base_of_v<MatrixExpr<LHS>, LHS> &&
             std::is_base_of_v<MatrixExpr<RHS>, RHS>
         >>
constexpr auto operator+(const MatrixExpr<LHS>& lhs, const MatrixExpr<RHS>& rhs) {
    return MatrixBinaryExpr<LHS, RHS, detail::ops::Add>(
        lhs.derived(), rhs.derived(), detail::ops::Add{}
    );
}

/**
 * @brief Subtraction operator for matrix expressions
 */
template<typename LHS, typename RHS,
         typename = std::enable_if_t<
             std::is_base_of_v<MatrixExpr<LHS>, LHS> &&
             std::is_base_of_v<MatrixExpr<RHS>, RHS>
         >>
constexpr auto operator-(const MatrixExpr<LHS>& lhs, const MatrixExpr<RHS>& rhs) {
    return MatrixBinaryExpr<LHS, RHS, detail::ops::Sub>(
        lhs.derived(), rhs.derived(), detail::ops::Sub{}
    );
}

/**
 * @brief Matrix multiplication operator
 */
template<typename LHS, typename RHS,
         typename = std::enable_if_t<
             std::is_base_of_v<MatrixExpr<LHS>, LHS> &&
             std::is_base_of_v<MatrixExpr<RHS>, RHS>
         >>
constexpr auto operator*(const MatrixExpr<LHS>& lhs, const MatrixExpr<RHS>& rhs) {
    return MatrixMulExpr<LHS, RHS>(lhs.derived(), rhs.derived());
}

/**
 * @brief Element-wise multiplication (Hadamard product)
 */
template<typename LHS, typename RHS,
         typename = std::enable_if_t<
             std::is_base_of_v<MatrixExpr<LHS>, LHS> &&
             std::is_base_of_v<MatrixExpr<RHS>, RHS>
         >>
constexpr auto hadamard(const MatrixExpr<LHS>& lhs, const MatrixExpr<RHS>& rhs) {
    return MatrixBinaryExpr<LHS, RHS, detail::ops::Mul>(
        lhs.derived(), rhs.derived(), detail::ops::Mul{}
    );
}

/**
 * @brief Element-wise division
 */
template<typename LHS, typename RHS,
         typename = std::enable_if_t<
             std::is_base_of_v<MatrixExpr<LHS>, LHS> &&
             std::is_base_of_v<MatrixExpr<RHS>, RHS>
         >>
constexpr auto hadamard_div(const MatrixExpr<LHS>& lhs, const MatrixExpr<RHS>& rhs) {
    return MatrixBinaryExpr<LHS, RHS, detail::ops::Div>(
        lhs.derived(), rhs.derived(), detail::ops::Div{}
    );
}

/**
 * @brief Negation operator for matrix expressions
 */
template<typename Expr,
         typename = std::enable_if_t<
             std::is_base_of_v<MatrixExpr<Expr>, Expr>
         >>
constexpr auto operator-(const MatrixExpr<Expr>& expr) {
    return MatrixUnaryExpr<Expr, detail::ops::Negate>(
        expr.derived(), detail::ops::Negate{}
    );
}

/**
 * @brief Scalar multiplication operator (matrix * scalar)
 */
template<typename Expr, typename Scalar,
         typename = std::enable_if_t<
             std::is_base_of_v<MatrixExpr<Expr>, Expr> &&
             std::is_arithmetic_v<Scalar>
         >>
constexpr auto operator*(const MatrixExpr<Expr>& expr, Scalar scalar) {
    return MatrixScalarExpr<Expr, Scalar>(expr.derived(), scalar);
}

/**
 * @brief Scalar multiplication operator (scalar * matrix)
 */
template<typename Scalar, typename Expr,
         typename = std::enable_if_t<
             std::is_arithmetic_v<Scalar> &&
             std::is_base_of_v<MatrixExpr<Expr>, Expr>
         >>
constexpr auto operator*(Scalar scalar, const MatrixExpr<Expr>& expr) {
    return MatrixScalarExpr<Expr, Scalar>(expr.derived(), scalar);
}

/**
 * @brief Scalar division operator (matrix / scalar)
 */
template<typename Expr, typename Scalar,
         typename = std::enable_if_t<
             std::is_base_of_v<MatrixExpr<Expr>, Expr> &&
             std::is_arithmetic_v<Scalar>
         >>
constexpr auto operator/(const MatrixExpr<Expr>& expr, Scalar scalar) {
    return MatrixScalarDivExpr<Expr, Scalar>(expr.derived(), scalar);
}

/**
 * @brief Transpose function
 */
template<typename Expr,
         typename = std::enable_if_t<
             std::is_base_of_v<MatrixExpr<Expr>, Expr>
         >>
constexpr auto transpose(const MatrixExpr<Expr>& expr) {
    return TransposeExpr<Expr>(expr.derived());
}

/**
 * @brief Element-wise absolute value
 */
template<typename Expr,
         typename = std::enable_if_t<
             std::is_base_of_v<MatrixExpr<Expr>, Expr>
         >>
constexpr auto abs(const MatrixExpr<Expr>& expr) {
    return MatrixUnaryExpr<Expr, detail::ops::Abs>(expr.derived(), detail::ops::Abs{});
}

/**
 * @brief Element-wise square root
 */
template<typename Expr,
         typename = std::enable_if_t<
             std::is_base_of_v<MatrixExpr<Expr>, Expr>
         >>
constexpr auto sqrt(const MatrixExpr<Expr>& expr) {
    return MatrixUnaryExpr<Expr, detail::ops::Sqrt>(expr.derived(), detail::ops::Sqrt{});
}

/**
 * @brief Compute Frobenius norm squared of matrix expression
 * @tparam Expr Matrix expression type
 * @param expr Matrix expression
 * @return Square of the Frobenius norm
 */
template<typename Expr,
         typename = std::enable_if_t<
             std::is_base_of_v<MatrixExpr<Expr>, Expr>
         >>
constexpr auto frobenius_norm_squared(const MatrixExpr<Expr>& expr) {
    using result_type = decltype(expr.derived()(0, 0) * expr.derived()(0, 0));
    result_type sum = result_type{0};
    const auto m = expr.rows();
    const auto n = expr.cols();
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            auto val = expr.derived()(i, j);
            sum += val * val;
        }
    }
    return sum;
}

/**
 * @brief Compute Frobenius norm of matrix expression
 * @tparam Expr Matrix expression type
 * @param expr Matrix expression
 * @return Frobenius norm
 */
template<typename Expr,
         typename = std::enable_if_t<
             std::is_base_of_v<MatrixExpr<Expr>, Expr>
         >>
constexpr auto frobenius_norm(const MatrixExpr<Expr>& expr) {
    using std::sqrt;
    return sqrt(frobenius_norm_squared(expr));
}

/**
 * @brief Compute trace of square matrix expression
 * @tparam Expr Matrix expression type
 * @param expr Matrix expression
 * @return Sum of diagonal elements
 */
template<typename Expr,
         typename = std::enable_if_t<
             std::is_base_of_v<MatrixExpr<Expr>, Expr>
         >>
constexpr auto trace(const MatrixExpr<Expr>& expr) {
    using result_type = decltype(expr.derived()(0, 0));
    result_type sum = result_type{0};
    const auto n = std::min(expr.rows(), expr.cols());
    for (std::size_t i = 0; i < n; ++i) {
        sum += expr.derived()(i, i);
    }
    return sum;
}

} // namespace math
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_MATH_MATRIX_EXPR_H
