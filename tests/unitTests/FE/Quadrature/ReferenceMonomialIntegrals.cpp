// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/** @file ReferenceMonomialIntegrals.cpp
 * @brief Stable test-only moments on the reference line and simplices.
 */

#include "ReferenceMonomialIntegrals.h"

#include "FE/Common/FEException.h"

#include <algorithm>
#include <array>

namespace svmp::FE::quadrature::test {

long double line_monomial_integral(int power)
{
    if (power < 0) {
        svmp::raise<InvalidArgumentException>("Line monomial exponent must be non-negative");
    }
    return power % 2 != 0 ? 0.0L : 2.0L / (static_cast<long double>(power) + 1.0L);
}

long double triangle_monomial_integral(int a, int b)
{
    if (a < 0 || b < 0) {
        svmp::raise<InvalidArgumentException>("Triangle monomial exponents must be non-negative");
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

long double tetrahedron_monomial_integral(int a, int b, int c)
{
    if (a < 0 || b < 0 || c < 0) {
        svmp::raise<InvalidArgumentException>("Tetrahedron monomial exponents must be non-negative");
    }
    std::array<int, 3> powers{a, b, c};
    std::sort(powers.begin(), powers.end());
    const long double largest = powers[2];
    // a! b! c! / (a+b+c+3)! without factorials or integer exponent sums.
    long double integral = 1.0L / ((largest + 1.0L) * (largest + 2.0L) * (largest + 3.0L));
    for (int k = 0; k < powers[1]; ++k) {
        const long double numerator = static_cast<long double>(k) + 1.0L;
        integral *= numerator / (largest + 3.0L + numerator);
    }
    const long double partial_degree = largest + static_cast<long double>(powers[1]);
    for (int k = 0; k < powers[0]; ++k) {
        const long double numerator = static_cast<long double>(k) + 1.0L;
        integral *= numerator / (partial_degree + 3.0L + numerator);
    }
    return integral;
}

} // namespace svmp::FE::quadrature::test
