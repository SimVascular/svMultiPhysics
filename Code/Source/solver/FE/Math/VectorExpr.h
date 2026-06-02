#ifndef SVMP_FE_MATH_VECTOR_EXPR_H
#define SVMP_FE_MATH_VECTOR_EXPR_H

/**
 * @file VectorExpr.h
 * @brief Expression template infrastructure for lazy evaluation of vector operations
 *
 * This header provides expression templates that enable compound vector operations
 * without creating temporary objects. Operations are evaluated lazily at the point
 * of assignment, eliminating intermediate allocations and improving performance.
 */

#include <cstddef>
#include <type_traits>
#include <cmath>
#include "ExpressionOps.h"

namespace svmp {
namespace FE {
namespace math {

/**
 * @brief Base class for all vector expressions using CRTP
 * @tparam Derived The derived expression type
 *
 * This uses the Curiously Recurring Template Pattern (CRTP) to provide
 * static polymorphism for expression templates.
 */
template<typename Derived>
class VectorExpr {
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
     * @brief Access element by index
     * @param i Element index
     * @return Value at index i
     */
    auto operator[](std::size_t i) const {
        return derived()[i];
    }

    /**
     * @brief Get the size of the vector expression
     * @return Number of elements
     */
    std::size_t size() const {
        return derived().size();
    }
};

/**
 * @brief Binary expression for element-wise operations between two vector expressions
 * @tparam LHS Left-hand side expression type
 * @tparam RHS Right-hand side expression type
 * @tparam Op Binary operation functor
 */
template<typename LHS, typename RHS, typename Op>
class VectorBinaryExpr : public VectorExpr<VectorBinaryExpr<LHS, RHS, Op>> {
private:
    const LHS& lhs_;
    const RHS& rhs_;
    Op op_;

public:
    /**
     * @brief Construct binary expression
     * @param lhs Left operand
     * @param rhs Right operand
     * @param op Operation to apply
     */
    constexpr VectorBinaryExpr(const LHS& lhs, const RHS& rhs, Op op = Op{})
        : lhs_(lhs), rhs_(rhs), op_(op) {}

    /**
     * @brief Access element at index
     * @param i Element index
     * @return Result of operation on elements at index i
     */
    constexpr auto operator[](std::size_t i) const {
        return op_(lhs_[i], rhs_[i]);
    }

    /**
     * @brief Get size of expression (from left operand)
     * @return Number of elements
     */
    constexpr std::size_t size() const {
        return lhs_.size();
    }
};

/**
 * @brief Unary expression for element-wise operations on a single vector expression
 * @tparam Expr Expression type
 * @tparam Op Unary operation functor
 */
template<typename Expr, typename Op>
class VectorUnaryExpr : public VectorExpr<VectorUnaryExpr<Expr, Op>> {
private:
    const Expr& expr_;
    Op op_;

public:
    /**
     * @brief Construct unary expression
     * @param expr Operand expression
     * @param op Operation to apply
     */
    constexpr VectorUnaryExpr(const Expr& expr, Op op = Op{})
        : expr_(expr), op_(op) {}

    /**
     * @brief Access element at index
     * @param i Element index
     * @return Result of operation on element at index i
     */
    constexpr auto operator[](std::size_t i) const {
        return op_(expr_[i]);
    }

    /**
     * @brief Get size of expression
     * @return Number of elements
     */
    constexpr std::size_t size() const {
        return expr_.size();
    }
};

/**
 * @brief Scalar multiplication expression
 * @tparam Expr Vector expression type
 * @tparam Scalar Scalar type
 */
template<typename Expr, typename Scalar>
class VectorScalarExpr : public VectorExpr<VectorScalarExpr<Expr, Scalar>> {
private:
    const Expr& expr_;
    Scalar scalar_;

public:
    /**
     * @brief Construct scalar multiplication expression
     * @param expr Vector expression
     * @param scalar Scalar value
     */
    constexpr VectorScalarExpr(const Expr& expr, Scalar scalar)
        : expr_(expr), scalar_(scalar) {}

    /**
     * @brief Access element at index
     * @param i Element index
     * @return Element multiplied by scalar
     */
    constexpr auto operator[](std::size_t i) const {
        return expr_[i] * scalar_;
    }

    /**
     * @brief Get size of expression
     * @return Number of elements
     */
    constexpr std::size_t size() const {
        return expr_.size();
    }
};

/**
 * @brief Scalar division expression
 * @tparam Expr Vector expression type
 * @tparam Scalar Scalar type
 */
template<typename Expr, typename Scalar>
class VectorScalarDivExpr : public VectorExpr<VectorScalarDivExpr<Expr, Scalar>> {
private:
    const Expr& expr_;
    Scalar scalar_;

public:
    /**
     * @brief Construct scalar division expression
     * @param expr Vector expression
     * @param scalar Scalar divisor
     */
    constexpr VectorScalarDivExpr(const Expr& expr, Scalar scalar)
        : expr_(expr), scalar_(scalar) {}

    /**
     * @brief Access element at index
     * @param i Element index
     * @return Element divided by scalar
     */
    constexpr auto operator[](std::size_t i) const {
        return expr_[i] / scalar_;
    }

    /**
     * @brief Get size of expression
     * @return Number of elements
     */
    constexpr std::size_t size() const {
        return expr_.size();
    }
};

/**
 * @brief Addition operator for vector expressions
 */
template<typename LHS, typename RHS,
         typename = std::enable_if_t<
             std::is_base_of_v<VectorExpr<LHS>, LHS> &&
             std::is_base_of_v<VectorExpr<RHS>, RHS>
         >>
constexpr auto operator+(const VectorExpr<LHS>& lhs, const VectorExpr<RHS>& rhs) {
    return VectorBinaryExpr<LHS, RHS, detail::ops::Add>(
        lhs.derived(), rhs.derived(), detail::ops::Add{}
    );
}

/**
 * @brief Subtraction operator for vector expressions
 */
template<typename LHS, typename RHS,
         typename = std::enable_if_t<
             std::is_base_of_v<VectorExpr<LHS>, LHS> &&
             std::is_base_of_v<VectorExpr<RHS>, RHS>
         >>
constexpr auto operator-(const VectorExpr<LHS>& lhs, const VectorExpr<RHS>& rhs) {
    return VectorBinaryExpr<LHS, RHS, detail::ops::Sub>(
        lhs.derived(), rhs.derived(), detail::ops::Sub{}
    );
}

/**
 * @brief Element-wise multiplication operator for vector expressions
 */
template<typename LHS, typename RHS,
         typename = std::enable_if_t<
             std::is_base_of_v<VectorExpr<LHS>, LHS> &&
             std::is_base_of_v<VectorExpr<RHS>, RHS>
         >>
constexpr auto hadamard(const VectorExpr<LHS>& lhs, const VectorExpr<RHS>& rhs) {
    return VectorBinaryExpr<LHS, RHS, detail::ops::Mul>(
        lhs.derived(), rhs.derived(), detail::ops::Mul{}
    );
}

/**
 * @brief Element-wise division operator for vector expressions
 */
template<typename LHS, typename RHS,
         typename = std::enable_if_t<
             std::is_base_of_v<VectorExpr<LHS>, LHS> &&
             std::is_base_of_v<VectorExpr<RHS>, RHS>
         >>
constexpr auto hadamard_div(const VectorExpr<LHS>& lhs, const VectorExpr<RHS>& rhs) {
    return VectorBinaryExpr<LHS, RHS, detail::ops::Div>(
        lhs.derived(), rhs.derived(), detail::ops::Div{}
    );
}

/**
 * @brief Negation operator for vector expressions
 */
template<typename Expr,
         typename = std::enable_if_t<
             std::is_base_of_v<VectorExpr<Expr>, Expr>
         >>
constexpr auto operator-(const VectorExpr<Expr>& expr) {
    return VectorUnaryExpr<Expr, detail::ops::Negate>(
        expr.derived(), detail::ops::Negate{}
    );
}

/**
 * @brief Scalar multiplication operator (vector * scalar)
 */
template<typename Expr, typename Scalar,
         typename = std::enable_if_t<
             std::is_base_of_v<VectorExpr<Expr>, Expr> &&
             std::is_arithmetic_v<Scalar>
         >>
constexpr auto operator*(const VectorExpr<Expr>& expr, Scalar scalar) {
    return VectorScalarExpr<Expr, Scalar>(expr.derived(), scalar);
}

/**
 * @brief Scalar multiplication operator (scalar * vector)
 */
template<typename Scalar, typename Expr,
         typename = std::enable_if_t<
             std::is_arithmetic_v<Scalar> &&
             std::is_base_of_v<VectorExpr<Expr>, Expr>
         >>
constexpr auto operator*(Scalar scalar, const VectorExpr<Expr>& expr) {
    return VectorScalarExpr<Expr, Scalar>(expr.derived(), scalar);
}

/**
 * @brief Scalar division operator (vector / scalar)
 */
template<typename Expr, typename Scalar,
         typename = std::enable_if_t<
             std::is_base_of_v<VectorExpr<Expr>, Expr> &&
             std::is_arithmetic_v<Scalar>
         >>
constexpr auto operator/(const VectorExpr<Expr>& expr, Scalar scalar) {
    return VectorScalarDivExpr<Expr, Scalar>(expr.derived(), scalar);
}

/**
 * @brief Element-wise absolute value
 */
template<typename Expr,
         typename = std::enable_if_t<
             std::is_base_of_v<VectorExpr<Expr>, Expr>
         >>
constexpr auto abs(const VectorExpr<Expr>& expr) {
    return VectorUnaryExpr<Expr, detail::ops::Abs>(expr.derived(), detail::ops::Abs{});
}

/**
 * @brief Element-wise square root
 */
template<typename Expr,
         typename = std::enable_if_t<
             std::is_base_of_v<VectorExpr<Expr>, Expr>
         >>
constexpr auto sqrt(const VectorExpr<Expr>& expr) {
    return VectorUnaryExpr<Expr, detail::ops::Sqrt>(expr.derived(), detail::ops::Sqrt{});
}

/**
 * @brief Dot product for vector expressions
 * @tparam LHS Left vector expression type
 * @tparam RHS Right vector expression type
 * @param lhs Left operand
 * @param rhs Right operand
 * @return Dot product result
 */
template<typename LHS, typename RHS,
         typename = std::enable_if_t<
             std::is_base_of_v<VectorExpr<LHS>, LHS> &&
             std::is_base_of_v<VectorExpr<RHS>, RHS>
         >>
constexpr auto dot(const VectorExpr<LHS>& lhs, const VectorExpr<RHS>& rhs) {
    using result_type = decltype(lhs.derived()[0] * rhs.derived()[0]);
    result_type sum = result_type{0};
    const auto n = lhs.size();
    for (std::size_t i = 0; i < n; ++i) {
        sum += lhs.derived()[i] * rhs.derived()[i];
    }
    return sum;
}

/**
 * @brief Compute norm squared of vector expression
 * @tparam Expr Vector expression type
 * @param expr Vector expression
 * @return Square of the Euclidean norm
 */
template<typename Expr,
         typename = std::enable_if_t<
             std::is_base_of_v<VectorExpr<Expr>, Expr>
         >>
constexpr auto norm_squared(const VectorExpr<Expr>& expr) {
    return dot(expr, expr);
}

/**
 * @brief Compute norm of vector expression
 * @tparam Expr Vector expression type
 * @param expr Vector expression
 * @return Euclidean norm
 */
template<typename Expr,
         typename = std::enable_if_t<
             std::is_base_of_v<VectorExpr<Expr>, Expr>
         >>
constexpr auto norm(const VectorExpr<Expr>& expr) {
    using std::sqrt;
    return sqrt(norm_squared(expr));
}

/**
 * @brief Normalize vector expression
 * @tparam Expr Vector expression type
 * @param expr Vector expression
 * @return Normalized vector expression
 */
template<typename Expr,
         typename = std::enable_if_t<
             std::is_base_of_v<VectorExpr<Expr>, Expr>
         >>
constexpr auto normalize(const VectorExpr<Expr>& expr) {
    return expr / norm(expr);
}

} // namespace math
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_MATH_VECTOR_EXPR_H