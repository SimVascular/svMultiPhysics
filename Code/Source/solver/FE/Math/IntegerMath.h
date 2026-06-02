/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#ifndef SVMP_FE_MATH_INTEGERMATH_H
#define SVMP_FE_MATH_INTEGERMATH_H

#include "Types.h"

#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace svmp {
namespace FE {
namespace math {

[[nodiscard]] constexpr Real pow_int_nonnegative(Real base, int exponent) noexcept {
    Real result = Real(1);
    Real factor = base;
    int power = exponent;
    while (power > 0) {
        if ((power & 1) != 0) {
            result *= factor;
        }
        power >>= 1;
        if (power > 0) {
            factor *= factor;
        }
    }
    return result;
}

[[nodiscard]] constexpr Real pow_int(Real base, int exponent) noexcept {
    if (exponent < 0) {
        return Real(1) / pow_int_nonnegative(base, -exponent);
    }
    return pow_int_nonnegative(base, exponent);
}

[[nodiscard]] constexpr std::size_t binomial_size(int n, int k) {
    if (n < 0 || k < 0 || k > n) {
        return 0u;
    }
    if (k > n - k) {
        k = n - k;
    }

    std::size_t result = 1u;
    for (int i = 1; i <= k; ++i) {
        auto numerator = static_cast<std::size_t>(n - (k - i));
        auto denominator = static_cast<std::size_t>(i);

        const auto numerator_gcd = std::gcd(numerator, denominator);
        numerator /= numerator_gcd;
        denominator /= numerator_gcd;

        const auto result_gcd = std::gcd(result, denominator);
        result /= result_gcd;
        denominator /= result_gcd;
        if (denominator != 1u) {
            throw std::overflow_error(
                "binomial_size: failed to reduce exact binomial factor");
        }
        if (numerator != 0u &&
            result > std::numeric_limits<std::size_t>::max() / numerator) {
            throw std::overflow_error("binomial_size: result does not fit in size_t");
        }
        result *= numerator;
    }
    return result;
}

[[nodiscard]] constexpr Real binomial_real(int n, int k) noexcept {
    if (k < 0 || k > n) {
        return Real(0);
    }
    if (k > n - k) {
        k = n - k;
    }

    Real result = Real(1);
    for (int i = 1; i <= k; ++i) {
        result *= static_cast<Real>(n - (k - i));
        result /= static_cast<Real>(i);
    }
    return result;
}

} // namespace math
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_MATH_INTEGERMATH_H
