#ifndef SVMP_FE_MATH_EXPRESSION_OPS_H
#define SVMP_FE_MATH_EXPRESSION_OPS_H

/**
 * @file ExpressionOps.h
 * @brief Common expression template operators for vector and matrix expressions
 *
 * This header provides shared operator functors used by both VectorExpr.h and
 * MatrixExpr.h to avoid code duplication and namespace conflicts. All operators
 * are defined in the detail::ops namespace for internal use by expression templates.
 */

#include <cmath>

namespace svmp {
namespace FE {
namespace math {
namespace detail {
namespace ops {

/**
 * @brief Addition operator functor
 */
struct Add {
    template<typename T1, typename T2>
    constexpr auto operator()(const T1& a, const T2& b) const {
        return a + b;
    }
};

/**
 * @brief Subtraction operator functor
 */
struct Sub {
    template<typename T1, typename T2>
    constexpr auto operator()(const T1& a, const T2& b) const {
        return a - b;
    }
};

/**
 * @brief Multiplication operator functor
 */
struct Mul {
    template<typename T1, typename T2>
    constexpr auto operator()(const T1& a, const T2& b) const {
        return a * b;
    }
};

/**
 * @brief Division operator functor
 */
struct Div {
    template<typename T1, typename T2>
    constexpr auto operator()(const T1& a, const T2& b) const {
        return a / b;
    }
};

/**
 * @brief Negation operator functor
 */
struct Negate {
    template<typename T>
    constexpr auto operator()(const T& a) const {
        return -a;
    }
};

/**
 * @brief Absolute value operator functor
 */
struct Abs {
    template<typename T>
    constexpr auto operator()(const T& a) const {
        using std::abs;
        return abs(a);
    }
};

/**
 * @brief Square root operator functor
 */
struct Sqrt {
    template<typename T>
    constexpr auto operator()(const T& a) const {
        using std::sqrt;
        return sqrt(a);
    }
};

} // namespace ops
} // namespace detail
} // namespace math
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_MATH_EXPRESSION_OPS_H
