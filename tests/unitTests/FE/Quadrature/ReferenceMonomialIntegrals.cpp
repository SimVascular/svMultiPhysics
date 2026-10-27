// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/** @file ReferenceMonomialIntegrals.cpp
 * @brief Stable test-only moments on the reference line and triangle.
 */

#include "ReferenceMonomialIntegrals.h"

#include <algorithm>
#include <stdexcept>

namespace svmp::FE::quadrature::test {

long double line_monomial_integral(int power)
{
    if (power < 0) {
        throw std::invalid_argument("Line monomial exponent must be non-negative");
    }
    return power % 2 != 0 ? 0.0L : 2.0L / (static_cast<long double>(power) + 1.0L);
}

long double triangle_monomial_integral(int a, int b)
{
    if (a < 0 || b < 0) {
        throw std::invalid_argument("Triangle monomial exponents must be non-negative");
    }
    const int smaller = std::min(a, b);
    const long double larger = std::max(a, b);
    // a! b! / (a+b+2)! without separately forming large factorials.
    long double integral = 1.0L / ((larger + 1.0L) * (larger + 2.0L));
    for (int k = 0; k < smaller; ++k) {
        const long double numerator = static_cast<long double>(k) + 1.0L;
        integral *= numerator / (larger + 2.0L + numerator);
    }
    return integral;
}

} // namespace svmp::FE::quadrature::test
