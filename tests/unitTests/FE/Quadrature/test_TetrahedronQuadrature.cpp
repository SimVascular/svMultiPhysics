// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file test_TetrahedronQuadrature.cpp
 * @brief Tests for centroid and supplied-factor tetrahedron quadrature.
 */

#include <gtest/gtest.h>

#include "FE/Common/FEException.h"
#include "FE/Quadrature/GaussLobattoQuadrature.h"
#include "FE/Quadrature/GaussQuadrature.h"
#include "FE/Quadrature/TetrahedronQuadrature.h"
#include "ReferenceMonomialIntegrals.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

namespace reference = svmp::FE::quadrature::test;
using namespace svmp::FE::quadrature;
using svmp::FE::ConvergenceException;
using svmp::FE::InvalidArgumentException;

struct CompensatedSum {
    long double sum = 0.0L;
    long double correction = 0.0L;

    void add(long double term)
    {
        const long double next = sum + term;
        correction += std::abs(sum) >= std::abs(term)
            ? (sum - next) + term : (term - next) + sum;
        sum = next;
    }

    long double value() const { return sum + correction; }
};

long double direct_moment(const QuadratureRule& rule, int a, int b, int c)
{
    CompensatedSum sum;
    for (std::size_t q = 0; q < rule.num_points(); ++q) {
        sum.add(static_cast<long double>(rule.weight(q)) *
            std::pow(static_cast<long double>(rule.point(q)[0]), a) *
            std::pow(static_cast<long double>(rule.point(q)[1]), b) *
            std::pow(static_cast<long double>(rule.point(q)[2]), c));
    }
    return sum.value();
}

void expect_same_rule(const QuadratureRule& left, const QuadratureRule& right)
{
    ASSERT_EQ(left.num_points(), right.num_points());
    EXPECT_EQ(left.cell_family(), right.cell_family());
    EXPECT_EQ(left.dimension(), right.dimension());
    EXPECT_EQ(left.reference_cell_measure(), right.reference_cell_measure());
    EXPECT_EQ(left.polynomial_exactness(), right.polynomial_exactness());
    for (std::size_t q = 0; q < left.num_points(); ++q) {
        if (left.weight(q) != right.weight(q)) {
            FAIL() << "weight at sample " << q;
        }
        for (std::size_t axis = 0; axis < 3; ++axis) {
            if (left.point(q)[axis] != right.point(q)[axis]) {
                FAIL() << "coordinate at sample " << q << ", axis " << axis;
            }
        }
    }
}

QuadratureRule reversed(const QuadratureRule& rule)
{
    return QuadratureRule(rule.cell_family(), rule.polynomial_exactness(),
        {rule.points().rbegin(), rule.points().rend()},
        {rule.weights().rbegin(), rule.weights().rend()});
}

template<class Function>
void expect_invalid(Function&& function, std::string_view diagnostic)
{
    try {
        function();
        FAIL() << "Expected InvalidArgumentException containing " << diagnostic;
    } catch (const InvalidArgumentException& error) {
        EXPECT_NE(std::string_view(error.what()).find(diagnostic), std::string_view::npos)
            << error.what();
    }
}

static_assert(max_tetrahedron_exactness() == 253);
static_assert(noexcept(max_tetrahedron_exactness()));
static_assert(std::is_same_v<decltype(&make_tetrahedron_centroid_rule), QuadratureRule (*)()>);
static_assert(std::is_same_v<decltype(make_tetrahedron_rule(0, std::declval<const QuadratureRule&>())), QuadratureRule>);
static_assert(std::is_same_v<decltype(make_tetrahedron_rule(0,
    std::declval<const QuadratureRule&>(), std::declval<const QuadratureRule&>(),
    std::declval<const QuadratureRule&>())), QuadratureRule>);

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

TEST(TetrahedronCentroidImplementation, IntegratesConstantsAndLinearFunctions)
{
    const auto rule = make_tetrahedron_centroid_rule();
    expect_same_rule(rule, QuadratureRule(svmp::CellFamily::Tetra, 1, {{0.25, 0.25, 0.25}}, {1.0 / 6.0}));
    EXPECT_EQ(rule.dimension(), 3u);
    EXPECT_EQ(rule.reference_cell_measure(), 1.0 / 6.0);
    for (const auto powers : {std::array{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}) {
        EXPECT_NEAR(direct_moment(rule, powers[0], powers[1], powers[2]),
            reference::tetrahedron_monomial_integral(powers[0], powers[1], powers[2]), 2.0e-17L);
    }
    EXPECT_NEAR(direct_moment(rule, 2, 0, 0), 1.0L / 96.0L, 2.0e-18L);
    EXPECT_GT(std::abs(direct_moment(rule, 2, 0, 0) - reference::tetrahedron_monomial_integral(2, 0, 0)), 0.006L);
}

TEST(TetrahedronDuffyImplementation, UsesTheSuppliedFactorAtLowAndMaximumRequests)
{
    struct Case { int request; std::size_t count; int degree; };
    for (const auto fixture : {Case{0, 8, 1}, {1, 8, 1}, {2, 27, 3}, {3, 27, 3},
                              {4, 64, 5}, {5, 64, 5}, {252, 2097152, 253}, {253, 2097152, 253}}) {
        SCOPED_TRACE(fixture.request);
        const auto line = make_gauss_legendre_rule(fixture.request + 2);
        const auto rule = make_tetrahedron_rule(fixture.request, line);
        EXPECT_EQ(rule.num_points(), fixture.count);
        EXPECT_EQ(rule.polynomial_exactness(), fixture.degree);
        expect_same_rule(rule, make_tetrahedron_rule(fixture.request, line, line, line));
    }
    const auto high = make_gauss_legendre_rule(7);
    const auto expected = make_tetrahedron_rule(5, high);
    expect_same_rule(make_tetrahedron_rule(0, high), expected);
    expect_same_rule(make_tetrahedron_rule(1, high), expected);
    const QuadratureRule advertised(svmp::CellFamily::Line, 300,
        {high.points().begin(), high.points().end()}, {high.weights().begin(), high.weights().end()});
    EXPECT_EQ(make_tetrahedron_rule(0, advertised).polynomial_exactness(), 298);
}

TEST(TetrahedronDuffyImplementation, PreservesAnIndependentOrderedUnequalFixture)
{
    const QuadratureRule u(svmp::CellFamily::Line, 3,
        {{-1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}, {1.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0});
    const auto v = make_gauss_legendre_rule(1);
    const QuadratureRule w(svmp::CellFamily::Line, 1, {{-0.5, 0.0, 0.0}, {0.5, 0.0, 0.0}}, {1.0, 1.0});
    const QuadratureRule expected(svmp::CellFamily::Tetra, 0,
        {{0.0, 0.5, 0.125}, {0.0, 0.5, 0.375}, {0.5, 0.25, 0.0625},
         {0.5, 0.25, 0.1875}, {1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}},
        {1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0, 0.0, 0.0});
    expect_same_rule(make_tetrahedron_rule(0, u, v, w), expected);
    expect_same_rule(make_tetrahedron_rule(0, reversed(u), v, reversed(w)), reversed(expected));
}

TEST(TetrahedronDuffyImplementation, UsesDirectionalExactnessAndOwnsItsStorage)
{
    const std::array factors{make_gauss_legendre_rule(7), make_gauss_lobatto_rule(5), make_gauss_legendre_rule(3)};
    const auto snapshots = factors;
    struct Case { std::array<int, 3> order; int degree; };
    for (const auto fixture : {Case{{0, 1, 2}, 3}, {{0, 2, 1}, 2}, {{1, 0, 2}, 3},
                              {{1, 2, 0}, 2}, {{2, 0, 1}, 1}, {{2, 1, 0}, 1}}) {
        const auto rule = make_tetrahedron_rule(0,
            factors[fixture.order[0]], factors[fixture.order[1]], factors[fixture.order[2]]);
        EXPECT_EQ(rule.polynomial_exactness(), fixture.degree);
        EXPECT_EQ(rule.num_points(), 32u);
        EXPECT_NEAR(direct_moment(rule, 0, 0, 0), 1.0L / 6.0L, 2.0e-16L);
        for (const auto& factor : factors) {
            EXPECT_NE(rule.points().data(), factor.points().data());
            EXPECT_NE(rule.weights().data(), factor.weights().data());
        }
    }
    for (std::size_t i = 0; i < factors.size(); ++i) {
        expect_same_rule(factors[i], snapshots[i]);
    }
    const auto temporary = make_tetrahedron_rule(0,
        make_gauss_legendre_rule(3), make_gauss_lobatto_rule(1), make_gauss_legendre_rule(1));
    EXPECT_EQ(temporary.num_points(), 4u);
    EXPECT_NEAR(direct_moment(temporary, 1, 0, 0), 1.0L / 24.0L, 2.0e-17L);
}

TEST(TetrahedronDuffyValidation, RejectsExtremeRequestsBeforeInspectingFactors)
{
    auto empty = make_gauss_legendre_rule(3);
    const auto retained = std::move(empty);
    ASSERT_EQ(empty.num_points(), 0u);
    for (const int request : {std::numeric_limits<int>::min(), -1, 254, std::numeric_limits<int>::max()}) {
        expect_invalid([&] { (void)make_tetrahedron_rule(request, retained); }, "requested_exactness");
        expect_invalid([&] { (void)make_tetrahedron_rule(request, empty); }, "requested_exactness");
        expect_invalid([&] { (void)make_tetrahedron_rule(request, empty, empty, empty); }, "requested_exactness");
    }
}

TEST(TetrahedronDuffyValidation, ChecksEachFactorContractInEveryDirection)
{
    const auto qualified = make_gauss_legendre_rule(3);
    const auto reject = [&](const QuadratureRule& invalid, std::string_view diagnostic) {
        expect_invalid([&] { (void)make_tetrahedron_rule(0, invalid); }, diagnostic);
        expect_invalid([&] { (void)make_tetrahedron_rule(0, invalid, qualified, qualified); }, diagnostic);
        expect_invalid([&] { (void)make_tetrahedron_rule(0, qualified, invalid, qualified); }, diagnostic);
        expect_invalid([&] { (void)make_tetrahedron_rule(0, qualified, qualified, invalid); }, diagnostic);
    };
    reject(QuadratureRule(svmp::CellFamily::Triangle, 3, {{0.25, 0.25, 0.0}}, {0.5}), "Line");
    for (const double x : {-1.0001, 1.0001}) {
        reject(QuadratureRule(svmp::CellFamily::Line, 3, {{x, 0.0, 0.0}}, {2.0}), "[-1,1]");
    }
    reject(QuadratureRule(svmp::CellFamily::Line, 3, {{0.0, 0.0, 0.0}}, {3.0}), "normalization");
    auto empty = qualified;
    const auto retained = std::move(empty);
    ASSERT_EQ(empty.num_points(), 0u);
    reject(empty, "empty");
    expect_same_rule(retained, qualified);
    for (int role = 0; role < 3; ++role) {
        std::array factors{qualified, qualified, qualified};
        factors[role] = QuadratureRule(svmp::CellFamily::Line, 2 - role, {{0.0, 0.0, 0.0}}, {2.0});
        expect_invalid([&] { (void)make_tetrahedron_rule(1, factors[0], factors[1], factors[2]); }, "exactness");
    }
    const auto lobatto = make_gauss_lobatto_rule(253);
    for (const int request : {252, 253}) {
        expect_invalid([&] { (void)make_tetrahedron_rule(request, lobatto); }, "exactness");
    }
}

TEST(TetrahedronDuffyValidation, BoundsTheFullProductBeforeFactorChecks)
{
    constexpr std::size_t capacity = 2097152;
    const auto two = make_gauss_legendre_rule(3);
    const auto single = make_gauss_legendre_rule(1);
    for (int role = 0; role < 3; ++role) {
        const std::size_t count = role == 0 ? capacity : capacity / 2;
        std::vector<QuadPoint> points(count, QuadPoint::Zero());
        if (role == 0) {
            for (std::size_t i = 0; i < count; ++i) {
                points[i] = two.point(i % 2);
            }
        }
        const QuadratureRule many(svmp::CellFamily::Line, 3, std::move(points),
            std::vector<double>(count, 2.0 / count));
        std::array factors{role == 0 ? single : two, single, single};
        factors[role] = many;
        const auto rule = make_tetrahedron_rule(0, factors[0], factors[1], factors[2]);
        EXPECT_EQ(rule.num_points(), capacity);
        EXPECT_NEAR(direct_moment(rule, 0, 0, 0), 1.0L / 6.0L, 2.0e-16L);
        // Deliberately invalid metadata/measure: the product bound has precedence.
        factors[role] = QuadratureRule(svmp::CellFamily::Line, 0,
            std::vector<QuadPoint>(count + 1, QuadPoint::Zero()), std::vector<double>(count + 1, 0.0));
        expect_invalid([&] { (void)make_tetrahedron_rule(0, factors[0], factors[1], factors[2]); }, "2097152");
    }
    const QuadratureRule oversized(svmp::CellFamily::Line, 0,
        std::vector<QuadPoint>(129, QuadPoint::Zero()), std::vector<double>(129, 0.0));
    expect_invalid([&] { (void)make_tetrahedron_rule(0, oversized); }, "2097152");
}

TEST(TetrahedronDuffyValidation, PreservesSignedZeroDuplicateAndCollapsedSamples)
{
    const auto u = make_gauss_lobatto_rule(3);
    const QuadratureRule signed_rule(svmp::CellFamily::Line, 1,
        std::vector<QuadPoint>(3, QuadPoint::Zero()), {3.0, -1.0, 0.0});
    const auto rule = make_tetrahedron_rule(0, u, signed_rule, signed_rule);
    const std::array<QuadPoint, 3> points{QuadPoint{0.0, 0.5, 0.25}, {0.5, 0.25, 0.125}, {1.0, 0.0, 0.0}};
    const std::array<double, 9> weights{9.0 / 48.0, -3.0 / 48.0, 0.0,
        -3.0 / 48.0, 1.0 / 48.0, 0.0, 0.0, 0.0, 0.0};
    ASSERT_EQ(rule.num_points(), 27u);
    for (std::size_t q = 0; q < rule.num_points(); ++q) {
        for (std::size_t axis = 0; axis < 3; ++axis) {
            EXPECT_EQ(rule.point(q)[axis], points[q / 9][axis]);
        }
        EXPECT_DOUBLE_EQ(rule.weight(q), q < 18 ? weights[q % 9] : 0.0);
    }
}

TEST(TetrahedronDuffyValidation, HandlesFiniteWeightsWhoseAbsoluteSumOverflows)
{
    const double node = 1.0 / std::sqrt(3.0);
    for (const double magnitude : {1.0, 1.0e103, 1.0e200}) {
        const QuadratureRule factor(svmp::CellFamily::Line, 3,
            {{-node, 0.0, 0.0}, {-node, 0.0, 0.0}, {-node, 0.0, 0.0}, {node, 0.0, 0.0}},
            {magnitude, -magnitude, 1.0, 1.0});
        if (magnitude == 1.0e200) {
            EXPECT_THROW((void)make_tetrahedron_rule(1, factor), ConvergenceException);
            continue;
        }
        const auto rule = make_tetrahedron_rule(1, factor);
        CompensatedSum absolute;
        for (double weight : rule.weights()) {
            ASSERT_TRUE(std::isfinite(weight));
            absolute.add(std::abs(static_cast<long double>(weight)));
        }
        EXPECT_EQ(rule.weight(0), -rule.weight(1));
        EXPECT_EQ(rule.weight(0), -rule.weight(4));
        EXPECT_EQ(rule.weight(0), -rule.weight(16));
        if (magnitude == 1.0e103) {
            EXPECT_GT(absolute.value(), static_cast<long double>(std::numeric_limits<double>::max()));
        } else {
            EXPECT_NEAR(direct_moment(rule, 0, 0, 0), 1.0L / 6.0L, 2.0e-16L);
        }
    }
}

TEST(TetrahedronDuffyValidation, DetectsOutputMeasureMismatchWithoutCertifyingFactorMoments)
{
    const auto qualified = make_gauss_legendre_rule(3);
    const auto single = make_gauss_legendre_rule(1);
    const QuadratureRule false_u(svmp::CellFamily::Line, 2, {{0.0, 0.0, 0.0}}, {2.0});
    const QuadratureRule false_v(svmp::CellFamily::Line, 1, {{0.5, 0.0, 0.0}}, {2.0});
    EXPECT_NO_THROW((void)make_tetrahedron_rule(0, qualified, single, single));
    EXPECT_THROW((void)make_tetrahedron_rule(0, false_u, single, single), ConvergenceException);
    EXPECT_THROW((void)make_tetrahedron_rule(0, qualified, false_v, single), ConvergenceException);
    // The two directional errors compensate in the measure. Metadata remains trusted.
    const QuadratureRule compensating_v(svmp::CellFamily::Line, 1, {{-1.0 / 3.0, 0.0, 0.0}}, {2.0});
    const auto accepted = make_tetrahedron_rule(0, false_u, compensating_v, single);
    EXPECT_EQ(accepted.polynomial_exactness(), 0);
    EXPECT_NEAR(direct_moment(accepted, 0, 0, 0), 1.0L / 6.0L, 2.0e-17L);
}
