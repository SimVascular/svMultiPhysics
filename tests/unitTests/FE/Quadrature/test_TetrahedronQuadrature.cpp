// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file test_TetrahedronQuadrature.cpp
 * @brief Tests for centroid and supplied-factor tetrahedron quadrature.
 */

#include <gtest/gtest.h>

#include "ReferenceMonomialIntegrals.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

namespace reference = svmp::FE::quadrature::test;

} // namespace

TEST(TetrahedronMonomialReferences, MatchesIndependentValues)
{
    EXPECT_EQ(reference::tetrahedron_monomial_integral(0, 0, 0), 1.0L / 6.0L);
    EXPECT_EQ(reference::tetrahedron_monomial_integral(1, 0, 0), 1.0L / 24.0L);
    EXPECT_EQ(reference::tetrahedron_monomial_integral(2, 0, 0), 1.0L / 60.0L);
    EXPECT_EQ(reference::tetrahedron_monomial_integral(1, 1, 0), 1.0L / 120.0L);
    EXPECT_NEAR(reference::tetrahedron_monomial_integral(1, 1, 1), 1.0L / 720.0L,
        std::numeric_limits<long double>::epsilon());

    // Exact integer factorials divided with 70-digit decimal precision.
    constexpr long double balanced = 3.607246336751302251048519164274536739744972043302599058518494788563773e-126L;
    constexpr long double axis = 6.030955689362359116875096495291029797745870001543924656476763933920025e-8L;
    constexpr long double mixed = 3.641858837811786551978253354892575639891623082005821937472743603847587e-40L;
    constexpr long double tolerance = 1024.0L * std::numeric_limits<long double>::epsilon();
    EXPECT_LE(std::abs(reference::tetrahedron_monomial_integral(84, 84, 85) / balanced - 1.0L), tolerance);
    EXPECT_LE(std::abs(reference::tetrahedron_monomial_integral(253, 0, 0) / axis - 1.0L), tolerance);
    EXPECT_LE(std::abs(reference::tetrahedron_monomial_integral(12, 28, 43) / mixed - 1.0L), tolerance);
}

TEST(TetrahedronMonomialReferences, IsPositiveFiniteAndPermutationInvariantThroughDegree253)
{
    for (int a = 0; a <= 253; ++a) {
        for (int b = a; b <= 253 - a; ++b) {
            for (int c = b; c <= 253 - a - b; ++c) {
                const long double expected = reference::tetrahedron_monomial_integral(a, b, c);
                ASSERT_TRUE(std::isfinite(expected));
                ASSERT_GT(expected, 0.0L);
                std::array<int, 3> powers{a, b, c};
                do {
                    EXPECT_EQ(reference::tetrahedron_monomial_integral(powers[0], powers[1], powers[2]), expected);
                } while (std::next_permutation(powers.begin(), powers.end()));
            }
        }
    }
}

TEST(TetrahedronMonomialReferences, RejectsEachNegativeExponent)
{
    for (const int negative : {-1, std::numeric_limits<int>::min()}) {
        EXPECT_THROW((void)reference::tetrahedron_monomial_integral(negative, 0, 0), std::invalid_argument);
        EXPECT_THROW((void)reference::tetrahedron_monomial_integral(0, negative, 0), std::invalid_argument);
        EXPECT_THROW((void)reference::tetrahedron_monomial_integral(0, 0, negative), std::invalid_argument);
    }
}
