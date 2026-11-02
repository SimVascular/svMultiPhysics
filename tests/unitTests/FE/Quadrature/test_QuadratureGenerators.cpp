// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file test_QuadratureGenerators.cpp
 * @brief Tests for line and triangle quadrature generators.
 */

#include <gtest/gtest.h>

#include "FE/Basis/NodeOrderingConventions.h"
#include "FE/Common/FEException.h"
#include "FE/Quadrature/GaussLobattoNodes.h"
#include "FE/Quadrature/GaussLobattoQuadrature.h"
#include "FE/Quadrature/GaussQuadrature.h"
#include "FE/Quadrature/QuadratureRule.h"
#include "FE/Quadrature/TriangleQuadrature.h"
#include "ReferenceMonomialIntegrals.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

using namespace svmp::FE;
using namespace svmp::FE::quadrature;

namespace {

namespace reference = svmp::FE::quadrature::test;

// Exhaustive supported-domain sweeps observed structure, measure, and moment
// errors below 1.2e-15; these envelopes retain cross-toolchain headroom
// through degree 255.
constexpr double kStructureTolerance = 1.0e-12;
constexpr double kMomentTolerance = 2.0e-12;
// Computed fixtures follow the production 64*epsilon refinement scale.
constexpr double kFixtureTolerance =
    64.0 * std::numeric_limits<double>::epsilon();
// Basis and Quadrature use different correction tolerances, so terminal
// binary64 rounding can differ; this remains far below node spacing.
constexpr double kBasisConsistencyTolerance =
    128.0 * std::numeric_limits<double>::epsilon();

static_assert(max_gauss_legendre_exactness() == 255);
static_assert(noexcept(max_gauss_legendre_exactness()));
static_assert(
    std::is_same_v<decltype(&make_gauss_legendre_rule),
                   QuadratureRule (*)(int)>);
static_assert(max_gauss_lobatto_exactness() == 253);
static_assert(noexcept(max_gauss_lobatto_exactness()));
static_assert(
    std::is_same_v<decltype(&make_gauss_lobatto_rule),
                   QuadratureRule (*)(int)>);

enum class LineEndpointPolicy {
    Excluded,
    Included,
};

long double accumulate_line_moment(
    const QuadratureRule& rule,
    std::size_t power)
{
    long double sum = 0.0L;
    long double correction = 0.0L;

    for (std::size_t point_index = 0;
         point_index < rule.num_points();
         ++point_index) {
        const long double coordinate =
            static_cast<long double>(rule.point(point_index)[0]);
        const long double term =
            static_cast<long double>(rule.weight(point_index)) *
            std::pow(coordinate, static_cast<int>(power));
        const long double next_sum = sum + term;

        if (std::abs(sum) >= std::abs(term)) {
            correction += (sum - next_sum) + term;
        } else {
            correction += (term - next_sum) + sum;
        }
        sum = next_sum;
    }

    return sum + correction;
}

void expect_common_line_metadata(
    const QuadratureRule& rule,
    std::size_t expected_num_points,
    int expected_exactness)
{
    EXPECT_EQ(rule.cell_family(), svmp::CellFamily::Line);
    EXPECT_EQ(rule.dimension(), 1u);
    EXPECT_EQ(rule.reference_cell_measure(), 2.0);
    EXPECT_EQ(rule.polynomial_exactness(), expected_exactness);
    ASSERT_EQ(rule.num_points(), expected_num_points);
    ASSERT_EQ(rule.points().size(), expected_num_points);
    ASSERT_EQ(rule.weights().size(), expected_num_points);

    for (std::size_t point_index = 0;
         point_index < rule.num_points();
         ++point_index) {
        SCOPED_TRACE(::testing::Message() << "point index=" << point_index);
        EXPECT_EQ(rule.point(point_index)[1], 0.0);
        EXPECT_EQ(rule.point(point_index)[2], 0.0);
    }
}

void expect_line_rule_invariants(
    const QuadratureRule& rule,
    LineEndpointPolicy endpoint_policy)
{
    ASSERT_GT(rule.num_points(), 0u);
    ASSERT_EQ(rule.points().size(), rule.weights().size());

    for (std::size_t point_index = 0;
         point_index < rule.num_points();
         ++point_index) {
        SCOPED_TRACE(::testing::Message() << "point index=" << point_index);

        const double coordinate = rule.point(point_index)[0];
        const double weight = rule.weight(point_index);
        EXPECT_TRUE(std::isfinite(coordinate));
        EXPECT_TRUE(std::isfinite(weight));
        EXPECT_GT(weight, 0.0);

        if (endpoint_policy == LineEndpointPolicy::Included) {
            EXPECT_GE(coordinate, -1.0);
            EXPECT_LE(coordinate, 1.0);
            if (point_index > 0u &&
                point_index + 1u < rule.num_points()) {
                EXPECT_GT(coordinate, -1.0);
                EXPECT_LT(coordinate, 1.0);
            }
        } else {
            EXPECT_GT(coordinate, -1.0);
            EXPECT_LT(coordinate, 1.0);
        }

        if (point_index > 0u) {
            EXPECT_LT(
                rule.point(point_index - 1u)[0],
                coordinate);
        }

        const std::size_t mirror_index =
            rule.num_points() - 1u - point_index;
        EXPECT_NEAR(
            coordinate,
            -rule.point(mirror_index)[0],
            kStructureTolerance);
        EXPECT_NEAR(weight, rule.weight(mirror_index), kStructureTolerance);
    }

    if (endpoint_policy == LineEndpointPolicy::Included) {
        ASSERT_GE(rule.num_points(), 2u);
        EXPECT_EQ(rule.point(0)[0], -1.0);
        EXPECT_EQ(rule.point(rule.num_points() - 1u)[0], 1.0);
    }

    if (rule.num_points() % 2u == 1u) {
        EXPECT_EQ(rule.point(rule.num_points() / 2u)[0], 0.0);
    }

    const long double measure_error = std::abs(
        accumulate_line_moment(rule, 0u) -
        static_cast<long double>(rule.reference_cell_measure()));
    EXPECT_LE(measure_error, static_cast<long double>(kStructureTolerance));
}

void expect_advertised_line_exactness(const QuadratureRule& rule)
{
    ASSERT_GE(rule.polynomial_exactness(), 0);

    for (int power = 0;
         power <= rule.polynomial_exactness();
         ++power) {
        SCOPED_TRACE(::testing::Message() << "monomial power=" << power);
        const std::size_t nonnegative_power =
            static_cast<std::size_t>(power);
        const long double error = std::abs(
            accumulate_line_moment(rule, nonnegative_power) -
            reference::line_monomial_integral(power));
        EXPECT_LE(error, static_cast<long double>(kMomentTolerance));
    }
}

void expect_every_supported_line_rule(
    QuadratureRule (*generator)(int),
    std::size_t first_num_points,
    std::size_t last_num_points,
    int exactness_subtrahend,
    LineEndpointPolicy endpoint_policy)
{
    for (std::size_t num_points = first_num_points;
         num_points <= last_num_points;
         ++num_points) {
        SCOPED_TRACE(
            ::testing::Message() << "num_points=" << num_points);
        const int expected_exactness =
            static_cast<int>(2 * num_points) - exactness_subtrahend;
        // Both requests must select this minimum point count and report its
        // actual exactness, not just echo the requested degree.
        for (const int requested_exactness :
             {expected_exactness - 1, expected_exactness}) {
            SCOPED_TRACE(
                ::testing::Message() << "requested_exactness=" << requested_exactness);
            const QuadratureRule rule = generator(requested_exactness);

            expect_common_line_metadata(
                rule, num_points, expected_exactness);
            expect_line_rule_invariants(rule, endpoint_policy);
            expect_advertised_line_exactness(rule);
        }
    }
}

void expect_canonical_rule(
    const QuadratureRule& rule,
    int expected_exactness,
    std::span<const double> expected_points,
    std::span<const double> expected_weights)
{
    SCOPED_TRACE(::testing::Message() << "num_points=" << expected_points.size());
    ASSERT_EQ(rule.num_points(), expected_points.size());
    ASSERT_EQ(rule.num_points(), expected_weights.size());
    EXPECT_EQ(rule.polynomial_exactness(), expected_exactness);

    // The exhaustive sweeps cover metadata, invariants, and advertised moments.
    // These fixtures independently anchor the samples and the exactness limit.
    for (std::size_t point_index = 0; point_index < rule.num_points(); ++point_index) {
        SCOPED_TRACE(::testing::Message() << "point index=" << point_index);
        const double expected_coordinate = expected_points[point_index];
        if (expected_coordinate == -1.0 || expected_coordinate == 0.0 ||
            expected_coordinate == 1.0) {
            EXPECT_EQ(rule.point(point_index)[0], expected_coordinate);
        } else {
            EXPECT_NEAR(rule.point(point_index)[0], expected_coordinate, kFixtureTolerance);
        }
        EXPECT_NEAR(rule.weight(point_index), expected_weights[point_index], kFixtureTolerance);
    }

    const std::size_t first_unadvertised_even_power =
        static_cast<std::size_t>(expected_exactness + 1);
    const long double first_unadvertised_error = std::abs(
        accumulate_line_moment(rule, first_unadvertised_even_power) -
        reference::line_monomial_integral(expected_exactness + 1));
    EXPECT_GT(first_unadvertised_error, static_cast<long double>(kMomentTolerance));
}

template <typename Function>
void expect_invalid_argument_with_message(
    Function&& function,
    std::string_view expected_substring)
{
    ASSERT_FALSE(expected_substring.empty());

    try {
        std::forward<Function>(function)();
        FAIL() << "Expected InvalidArgumentException containing: "
               << expected_substring;
    } catch (const InvalidArgumentException& exception) {
        const std::string_view actual_message{exception.what()};
        EXPECT_NE(
            actual_message.find(expected_substring),
            std::string_view::npos)
            << "actual message: " << actual_message;
    } catch (const std::exception& exception) {
        FAIL() << "Received a different exception type: " << exception.what();
    } catch (...) {
        FAIL() << "Received an unknown exception type";
    }
}

} // namespace

/**
 * @test QuadratureGeneratorTestSupport.ExercisesSharedLineRuleChecks exercises
 * the analytic moments, metadata, invariants, and exception-message helpers
 * against explicit interior and endpoint-inclusive line rules.
 */
TEST(QuadratureGeneratorTestSupport, ExercisesSharedLineRuleChecks)
{
    const double abscissa = std::sqrt(3.0 / 5.0);
    const QuadratureRule interior_rule(
        svmp::CellFamily::Line,
        5,
        {{-abscissa, 0.0, 0.0},
         {0.0, 0.0, 0.0},
         {abscissa, 0.0, 0.0}},
        {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0});
    const QuadratureRule endpoint_rule(
        svmp::CellFamily::Line,
        1,
        {{-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}},
        {1.0, 1.0});

    EXPECT_EQ(reference::line_monomial_integral(0), 2.0L);
    EXPECT_EQ(reference::line_monomial_integral(1), 0.0L);
    EXPECT_EQ(
        reference::line_monomial_integral(2),
        2.0L / 3.0L);
    EXPECT_EQ(reference::line_monomial_integral(255), 0.0L);
    EXPECT_LE(
        std::abs(accumulate_line_moment(interior_rule, 255u)),
        static_cast<long double>(kMomentTolerance));

    expect_common_line_metadata(interior_rule, 3u, 5);
    expect_line_rule_invariants(
        interior_rule,
        LineEndpointPolicy::Excluded);
    expect_advertised_line_exactness(interior_rule);

    expect_common_line_metadata(endpoint_rule, 2u, 1);
    expect_line_rule_invariants(
        endpoint_rule,
        LineEndpointPolicy::Included);
    expect_advertised_line_exactness(endpoint_rule);
    expect_invalid_argument_with_message(
        [] {
            (void)QuadratureRule(
                svmp::CellFamily::Line, 1, {}, {});
        },
        "at least one point");
}

/**
 * @test GaussLegendreImplementation.GeneratesCanonicalLowOrderRules compares
 * requests 0, 2, and 4 with analytic one-, two-, and three-point rules, including
 * their weights and the actual exactness exceeding each request by one.
 */
TEST(GaussLegendreImplementation, GeneratesCanonicalLowOrderRules)
{
    expect_canonical_rule(
        make_gauss_legendre_rule(0), 1,
        std::array{0.0},
        std::array{2.0});

    const double two_point_abscissa = 1.0 / std::sqrt(3.0);
    expect_canonical_rule(
        make_gauss_legendre_rule(2), 3,
        std::array{-two_point_abscissa, two_point_abscissa},
        std::array{1.0, 1.0});

    const double three_point_abscissa = std::sqrt(3.0 / 5.0);
    expect_canonical_rule(
        make_gauss_legendre_rule(4), 5,
        std::array{
            -three_point_abscissa, 0.0, three_point_abscissa},
        std::array{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0});
}

/**
 * @test GaussLegendreImplementation.GeneratesEverySupportedRule checks every
 * request from 0 through 255 for minimum point count, actual exactness, ordered
 * symmetric interior samples, positive weights, and all advertised moments.
 */
TEST(GaussLegendreImplementation, GeneratesEverySupportedRule)
{
    expect_every_supported_line_rule(
        &make_gauss_legendre_rule,
        1,
        128,
        1,
        LineEndpointPolicy::Excluded);
}

/**
 * @test GaussLegendreImplementation.RejectsRequestsOutsideSupportedRange checks
 * that negative requests, the first unsupported degree, and integer extremes
 * raise InvalidArgumentException with the supported-range diagnostic.
 */
TEST(GaussLegendreImplementation, RejectsRequestsOutsideSupportedRange)
{
    constexpr std::string_view expected_message =
        "requested_exactness must be in [0, 255]";
    constexpr std::array invalid_exactness{
        std::numeric_limits<int>::min(),
        -1,
        256,
        std::numeric_limits<int>::max()};

    for (const int requested_exactness : invalid_exactness) {
        SCOPED_TRACE(
            ::testing::Message() << "requested_exactness=" << requested_exactness);
        expect_invalid_argument_with_message(
            [requested_exactness] {
                (void)make_gauss_legendre_rule(requested_exactness);
            },
            expected_message);
    }
}

/**
 * @test GaussLobattoImplementation.GeneratesCanonicalLowOrderRules compares
 * requests 0, 2, and 4 with analytic two-, three-, and four-point rules, including
 * exact endpoints, weights, and the actual exactness exceeding each request.
 */
TEST(GaussLobattoImplementation, GeneratesCanonicalLowOrderRules)
{
    expect_canonical_rule(
        make_gauss_lobatto_rule(0), 1,
        std::array{-1.0, 1.0},
        std::array{1.0, 1.0});
    expect_canonical_rule(
        make_gauss_lobatto_rule(2), 3,
        std::array{-1.0, 0.0, 1.0},
        std::array{1.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0});

    const double four_point_abscissa = 1.0 / std::sqrt(5.0);
    expect_canonical_rule(
        make_gauss_lobatto_rule(4), 5,
        std::array{
            -1.0, -four_point_abscissa, four_point_abscissa, 1.0},
        std::array{1.0 / 6.0, 5.0 / 6.0, 5.0 / 6.0, 1.0 / 6.0});
}

/**
 * @test GaussLobattoImplementation.GeneratesEverySupportedRule checks every
 * request from 0 through 253 for minimum point count, actual exactness, ordered
 * symmetric samples, exact endpoints, positive weights, and advertised moments.
 */
TEST(GaussLobattoImplementation, GeneratesEverySupportedRule)
{
    expect_every_supported_line_rule(
        &make_gauss_lobatto_rule,
        2,
        128,
        3,
        LineEndpointPolicy::Included);
}

/**
 * @test GaussLobattoImplementation.RejectsRequestsOutsideSupportedRange checks
 * that negative requests, the first unsupported degree, and integer extremes
 * raise InvalidArgumentException with the supported-range diagnostic.
 */
TEST(GaussLobattoImplementation, RejectsRequestsOutsideSupportedRange)
{
    constexpr std::string_view expected_message =
        "Gauss-Lobatto-Legendre generator: "
        "requested_exactness must be in [0, 253]";
    constexpr std::array invalid_exactness{
        std::numeric_limits<int>::min(),
        -1,
        254,
        std::numeric_limits<int>::max()};

    for (const int requested_exactness : invalid_exactness) {
        SCOPED_TRACE(
            ::testing::Message() << "requested_exactness=" << requested_exactness);
        expect_invalid_argument_with_message(
            [requested_exactness] {
                (void)make_gauss_lobatto_rule(requested_exactness);
            },
            expected_message);
    }
}

/**
 * @test GaussLobattoBasisConsistency.MatchesRepresentativeNodeDistributions
 * compares every supported Lobatto point count with Basis nodes, checking the
 * two callers' indexing, symmetry, endpoints, and exact centers.
 */
TEST(GaussLobattoBasisConsistency, MatchesRepresentativeNodeDistributions)
{
    for (std::size_t num_points = 2; num_points <= 128; ++num_points) {
        SCOPED_TRACE(
            ::testing::Message() << "num_points=" << num_points);
        const QuadratureRule rule =
            make_gauss_lobatto_rule(static_cast<int>(2 * num_points - 3));
        ASSERT_EQ(rule.num_points(), num_points);

        for (std::size_t point_index = 0;
             point_index < rule.num_points();
             ++point_index) {
            SCOPED_TRACE(
                ::testing::Message() << "point index=" << point_index);
            const double quadrature_coordinate =
                rule.point(point_index)[0];
            const double basis_coordinate =
                svmp::FE::basis::line_coord_pm_one(
                    static_cast<int>(point_index),
                    static_cast<int>(num_points - 1));

            if (point_index == 0u) {
                EXPECT_EQ(quadrature_coordinate, -1.0);
                EXPECT_EQ(basis_coordinate, -1.0);
            } else if (point_index + 1u == rule.num_points()) {
                EXPECT_EQ(quadrature_coordinate, 1.0);
                EXPECT_EQ(basis_coordinate, 1.0);
            } else if (num_points % 2 == 1 &&
                       point_index == rule.num_points() / 2u) {
                EXPECT_EQ(quadrature_coordinate, 0.0);
                EXPECT_EQ(basis_coordinate, 0.0);
            } else {
                EXPECT_NEAR(
                    quadrature_coordinate,
                    basis_coordinate,
                    kBasisConsistencyTolerance);
            }
        }
    }
}

TEST(GaussLobattoNodes, RefinesAnalyticNodesAndLegendreValues)
{
    const auto center = detail::refine_gauss_lobatto_node(3, 0, 1.0e-15);
    EXPECT_EQ(center.coordinate, 0.0);
    EXPECT_EQ(center.legendre_value, -0.5);

    const auto cubic = detail::refine_gauss_lobatto_node(4, 0, 1.0e-15);
    EXPECT_NEAR(cubic.coordinate, 1.0 / std::sqrt(5.0), kFixtureTolerance);
    EXPECT_NEAR(cubic.legendre_value, -1.0 / std::sqrt(5.0), kFixtureTolerance);

    const auto quartic = detail::refine_gauss_lobatto_node(5, 0, 1.0e-15);
    EXPECT_NEAR(quartic.coordinate, std::sqrt(3.0 / 7.0), kFixtureTolerance);
    EXPECT_NEAR(quartic.legendre_value, -3.0 / 7.0, kFixtureTolerance);
    const auto quartic_center = detail::refine_gauss_lobatto_node(5, 1, 1.0e-15);
    EXPECT_EQ(quartic_center.coordinate, 0.0);
    EXPECT_EQ(quartic_center.legendre_value, 3.0 / 8.0);
}

TEST(GaussLobattoNodes, RejectsInvalidRefinementRequests)
{
    for (const std::size_t num_points : {0u, 1u, 2u}) {
        EXPECT_THROW((void)detail::refine_gauss_lobatto_node(num_points, 0, 1.0e-15),
                     InvalidArgumentException);
    }
    EXPECT_THROW((void)detail::refine_gauss_lobatto_node(4, 1, 1.0e-15),
                 InvalidArgumentException);
    for (const double tolerance : {0.0, -1.0,
             std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::quiet_NaN()}) {
        EXPECT_THROW((void)detail::refine_gauss_lobatto_node(4, 0, tolerance),
                     InvalidArgumentException);
    }
}

namespace {

static_assert(max_triangle_exactness() == 254);
static_assert(noexcept(max_triangle_exactness()));
static_assert(std::is_same_v<decltype(&make_triangle_centroid_rule), QuadratureRule (*)()>);
static_assert(std::is_same_v<decltype(make_triangle_rule(0, std::declval<const QuadratureRule&>())), QuadratureRule>);
static_assert(std::is_same_v<decltype(make_triangle_rule(0,
    std::declval<const QuadratureRule&>(), std::declval<const QuadratureRule&>())), QuadratureRule>);

// Exhaustive qualification observed relative errors below 6e-14. This retains
// over 30-fold cross-toolchain headroom and stays meaningful for tiny moments.
constexpr long double kTriangleRelativeTolerance = 2.0e-12L;

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

// Group only samples whose actual returned x is exactly equal. For each group,
// sum actual w*y^b, then multiply by actual x^a. This is the full sample sum,
// algebraically regrouped: no assumption about row lengths, y, or weights.
// Entry [a][b] holds the moment of x^a y^b for a+b <= degree.
std::vector<std::vector<CompensatedSum>> triangle_moments(const QuadratureRule& rule, int degree)
{
    const std::size_t size = static_cast<std::size_t>(degree + 1);
    std::vector<std::vector<CompensatedSum>> moments(size, std::vector<CompensatedSum>(size));
    std::vector<CompensatedSum> row(size);
    for (std::size_t first = 0; first < rule.num_points();) {
        std::fill(row.begin(), row.end(), CompensatedSum{});
        const long double x = rule.point(first)[0];
        std::size_t end = first;
        do {
            const long double y = rule.point(end)[1];
            long double power = 1.0L;
            for (int b = 0; b <= degree; ++b) {
                row[b].add(static_cast<long double>(rule.weight(end)) * power);
                power *= y;
            }
            ++end;
        } while (end < rule.num_points() && rule.point(end)[0] == x);
        long double power = 1.0L;
        for (int a = 0; a <= degree; ++a) {
            for (int b = 0; b <= degree - a; ++b) {
                moments[a][b].add(power * row[b].value());
            }
            power *= x;
        }
        first = end;
    }
    return moments;
}

// A relative bound keeps tiny high-order moments from passing as zero.
void expect_triangle_moments(const QuadratureRule& rule, int degree)
{
    const auto moments = triangle_moments(rule, degree);
    for (int a = 0; a <= degree; ++a) {
        for (int b = 0; b <= degree - a; ++b) {
            const long double expected = reference::triangle_monomial_integral(a, b);
            const long double actual = moments[a][b].value();
            EXPECT_GT(actual, 0.0L) << "a=" << a << ", b=" << b;
            EXPECT_LE(std::abs(actual - expected) / expected, kTriangleRelativeTolerance)
                << "a=" << a << ", b=" << b;
        }
    }
}

void expect_same_ordered_rule(const QuadratureRule& left, const QuadratureRule& right)
{
    ASSERT_EQ(left.num_points(), right.num_points());
    EXPECT_EQ(left.cell_family(), right.cell_family());
    EXPECT_EQ(left.dimension(), right.dimension());
    EXPECT_EQ(left.reference_cell_measure(), right.reference_cell_measure());
    EXPECT_EQ(left.polynomial_exactness(), right.polynomial_exactness());
    for (std::size_t q = 0; q < left.num_points(); ++q) {
        for (std::size_t axis = 0; axis < 3; ++axis) {
            EXPECT_EQ(left.point(q)[axis], right.point(q)[axis]) << "point=" << q << ", axis=" << axis;
        }
        EXPECT_EQ(left.weight(q), right.weight(q)) << "point=" << q;
    }
}

void expect_every_triangle_rule(LineEndpointPolicy endpoint_policy)
{
    const bool lobatto = endpoint_policy == LineEndpointPolicy::Included;
    const int last_request = lobatto ? 252 : 254;
    int distinct_rules = 0;
    std::optional<QuadratureRule> previous;
    for (int p = 0; p <= last_request; ++p) {
        SCOPED_TRACE(::testing::Message() << "request=" << p << ", Lobatto=" << lobatto);
        const auto factor = lobatto ? make_gauss_lobatto_rule(p + 1)
                                    : make_gauss_legendre_rule(p + 1);
        const auto rule = make_triangle_rule(p, factor);
        expect_same_ordered_rule(rule, make_triangle_rule(p, factor, factor));
        const std::size_t n = static_cast<std::size_t>((p + 1) / 2 + (lobatto ? 2 : 1));
        const int degree = 2 * static_cast<int>(n) - (lobatto ? 4 : 2);
        ASSERT_EQ(factor.num_points(), n);
        ASSERT_EQ(rule.num_points(), n * n);
        ASSERT_EQ(rule.polynomial_exactness(), degree);
        EXPECT_GE(degree, p);
        EXPECT_EQ(rule.cell_family(), svmp::CellFamily::Triangle);
        EXPECT_EQ(rule.dimension(), 2u);
        EXPECT_EQ(rule.reference_cell_measure(), 0.5);
        for (std::size_t q = 0; q < rule.num_points(); ++q) {
            const double x = rule.point(q)[0];
            const double y = rule.point(q)[1];
            EXPECT_TRUE(std::isfinite(x) && std::isfinite(y) && std::isfinite(rule.weight(q)));
            EXPECT_EQ(rule.point(q)[2], 0.0);
            if (lobatto) {
                EXPECT_GE(x, 0.0);
                EXPECT_GE(y, 0.0);
                EXPECT_LE(x + y, 1.0 + kStructureTolerance);
                EXPECT_GE(rule.weight(q), 0.0);
            } else {
                EXPECT_GT(x, 0.0);
                EXPECT_GT(y, 0.0);
                EXPECT_LT(x + y, 1.0);
                EXPECT_GT(rule.weight(q), 0.0);
            }
        }
        if (lobatto) {
            EXPECT_EQ(rule.point(0)[0], 0.0);
            EXPECT_EQ(rule.point(0)[1], 0.0);
            EXPECT_EQ(rule.point(n - 1)[1], 1.0);
            for (std::size_t j = 0; j < n; ++j) {
                const auto q = (n - 1) * n + j;
                EXPECT_EQ(rule.point(q)[0], 1.0);
                EXPECT_EQ(rule.point(q)[1], 0.0);
                EXPECT_EQ(rule.weight(q), 0.0);
                EXPECT_EQ(rule.point(j * n)[1], 0.0);
            }
        }
        if (previous && previous->polynomial_exactness() == degree) {
            expect_same_ordered_rule(*previous, rule);
        } else {
            ++distinct_rules;
            expect_triangle_moments(rule, degree);
        }
        previous = rule;
    }
    EXPECT_EQ(previous->num_points(), 16384u);
    EXPECT_EQ(distinct_rules, lobatto ? 127 : 128);
}

} // namespace

TEST(QuadratureMonomialReferences, MatchesIndependentValuesAndRejectsNegativePowers)
{
    EXPECT_EQ(reference::line_monomial_integral(0), 2.0L);
    EXPECT_EQ(reference::line_monomial_integral(1), 0.0L);
    EXPECT_EQ(reference::line_monomial_integral(254), 2.0L / 255.0L);
    EXPECT_EQ(reference::triangle_monomial_integral(0, 0), 0.5L);
    EXPECT_EQ(reference::triangle_monomial_integral(1, 0), 1.0L / 6.0L);
    EXPECT_EQ(reference::triangle_monomial_integral(0, 1), 1.0L / 6.0L);
    EXPECT_EQ(reference::triangle_monomial_integral(2, 0), 1.0L / 12.0L);
    EXPECT_EQ(reference::triangle_monomial_integral(1, 1), 1.0L / 24.0L);
    // Independently evaluated with exact integer factorials and 65-digit division.
    // Headroom for at most 127 ratio multiplications, also when long double
    // has binary64 precision.
    constexpr long double reference_tolerance = 512.0L * std::numeric_limits<long double>::epsilon();
    constexpr long double balanced = 1.05804760028269023389812750598454900539e-80L;
    constexpr long double unbalanced = 6.87153782458757607635065222401753811610e-27L;
    EXPECT_LE(std::abs(reference::triangle_monomial_integral(127, 127) / balanced - 1.0L), reference_tolerance);
    EXPECT_LE(std::abs(reference::triangle_monomial_integral(32, 48) / unbalanced - 1.0L), reference_tolerance);
    for (int a = 0; a <= 254; ++a) {
        for (int b = 0; b <= 254 - a; ++b) {
            const long double value = reference::triangle_monomial_integral(a, b);
            EXPECT_TRUE(std::isfinite(value));
            ASSERT_GT(value, 0.0L);
            EXPECT_EQ(value, reference::triangle_monomial_integral(b, a));
        }
    }
    EXPECT_THROW(reference::line_monomial_integral(-1), std::invalid_argument);
    EXPECT_THROW(reference::line_monomial_integral(std::numeric_limits<int>::min()), std::invalid_argument);
    EXPECT_THROW(reference::triangle_monomial_integral(-1, 0), std::invalid_argument);
    EXPECT_THROW(reference::triangle_monomial_integral(0, -1), std::invalid_argument);
}

TEST(QuadratureGeneratorTestSupport, GroupsOnlyEqualActualCoordinates)
{
    const QuadratureRule samples(svmp::CellFamily::Triangle, 0,
        {{0.25, 0.125, 0.0}, {0.25, 0.5, 0.0}, {0.5, 0.25, 0.0},
         {0.25, 0.0, 0.0}, {0.0, 0.75, 0.0}}, {0.5, -0.125, 0.25, 0.0, 0.125});
    constexpr int degree = 6;
    const auto moments = triangle_moments(samples, degree);
    for (int a = 0; a <= degree; ++a) {
        for (int b = 0; b <= degree - a; ++b) {
            CompensatedSum direct;
            for (std::size_t q = 0; q < samples.num_points(); ++q) {
                direct.add(samples.weight(q) * std::pow(static_cast<long double>(samples.point(q)[0]), a)
                    * std::pow(static_cast<long double>(samples.point(q)[1]), b));
            }
            EXPECT_EQ(moments[a][b].value(), direct.value());
        }
    }
}

TEST(TriangleCentroidImplementation, IntegratesLinearButNotQuadraticMonomials)
{
    const auto rule = make_triangle_centroid_rule();
    expect_same_ordered_rule(rule, QuadratureRule(svmp::CellFamily::Triangle, 1,
        {{1.0 / 3.0, 1.0 / 3.0, 0.0}}, {0.5}));
    const auto moments = triangle_moments(rule, 2);
    EXPECT_EQ(moments[0][0].value(), 0.5L);
    EXPECT_NEAR(moments[0][1].value(), 1.0L / 6.0L, kFixtureTolerance);
    EXPECT_NEAR(moments[1][0].value(), 1.0L / 6.0L, kFixtureTolerance);
    EXPECT_GT(std::abs(moments[2][0].value() - 1.0L / 12.0L), 0.02L);
}

TEST(TriangleDuffyImplementation, PreservesIndependentOrderedFixtureAndInputOwnership)
{
    const QuadratureRule factor(svmp::CellFamily::Line, 1,
        {{-0.5, 0.0, 0.0}, {0.5, 0.0, 0.0}}, {1.0, 1.0});
    const auto snapshot = factor;
    const auto rule = make_triangle_rule(0, factor);
    const QuadratureRule expected(svmp::CellFamily::Triangle, 0,
        {{0.25, 0.1875, 0.0}, {0.25, 0.5625, 0.0}, {0.75, 0.0625, 0.0}, {0.75, 0.1875, 0.0}},
        {0.1875, 0.1875, 0.0625, 0.0625});
    expect_same_ordered_rule(rule, expected);
    expect_same_ordered_rule(factor, snapshot);
    expect_same_ordered_rule(rule, make_triangle_rule(0, factor));
    // The temporary factor is destroyed before inspecting the result.
    const auto owned = make_triangle_rule(0, QuadratureRule(factor));
    expect_same_ordered_rule(owned, expected);
    const QuadratureRule reversed(svmp::CellFamily::Line, 1,
        {{0.5, 0.0, 0.0}, {-0.5, 0.0, 0.0}}, {1.0, 1.0});
    expect_same_ordered_rule(make_triangle_rule(0, reversed),
        QuadratureRule(svmp::CellFamily::Triangle, 0,
            {expected.points().rbegin(), expected.points().rend()},
            {expected.weights().rbegin(), expected.weights().rend()}));
}

TEST(TriangleDuffyImplementation, UsesSuppliedFactorForLowRequestsAndReportsDeliveredDegree)
{
    expect_same_ordered_rule(make_triangle_rule(0, make_gauss_legendre_rule(1)),
        QuadratureRule(svmp::CellFamily::Triangle, 0, {{0.5, 0.25, 0.0}}, {0.5}));
    const auto factor = make_gauss_legendre_rule(7);
    const auto zero = make_triangle_rule(0, factor);
    EXPECT_EQ(zero.num_points(), 16u);
    EXPECT_EQ(zero.polynomial_exactness(), 6);
    expect_same_ordered_rule(zero, make_triangle_rule(1, factor));
    expect_same_ordered_rule(make_triangle_rule(0, make_gauss_lobatto_rule(1)),
        QuadratureRule(svmp::CellFamily::Triangle, 0,
            {{0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}},
            {0.25, 0.25, 0.0, 0.0}));
}

TEST(TriangleDuffyImplementation, PreservesUnequalOrderedFactorsAndOwnsTheirProduct)
{
    const QuadratureRule u(svmp::CellFamily::Line, 1,
        {{-0.5, 0.0, 0.0}, {0.5, 0.0, 0.0}}, {1.0, 1.0});
    // Unsorted endpoints, a duplicate, and signed/zero weights; only degree zero.
    const QuadratureRule v(svmp::CellFamily::Line, 0,
        {{1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}, {3.0, -1.0, 0.0});
    const auto u_snapshot = u;
    const auto v_snapshot = v;
    const QuadratureRule expected(svmp::CellFamily::Triangle, 0,
        {{0.25, 0.75, 0.0}, {0.25, 0.0, 0.0}, {0.25, 0.75, 0.0},
         {0.75, 0.25, 0.0}, {0.75, 0.0, 0.0}, {0.75, 0.25, 0.0}},
        {0.5625, -0.1875, 0.0, 0.1875, -0.0625, 0.0});
    expect_same_ordered_rule(make_triangle_rule(0, u, v), expected);
    expect_same_ordered_rule(u, u_snapshot);
    expect_same_ordered_rule(v, v_snapshot);
    const auto owned = make_triangle_rule(0, QuadratureRule(u), QuadratureRule(v));
    expect_same_ordered_rule(owned, expected);
}

TEST(TriangleDuffyImplementation, EnforcesDirectionalExactnessAndReportsTheLimitingFactor)
{
    const auto high = make_gauss_legendre_rule(7);
    const auto low = make_gauss_legendre_rule(5);
    const auto rule = make_triangle_rule(5, high, low);
    EXPECT_EQ(rule.num_points(), 12u);
    EXPECT_EQ(rule.polynomial_exactness(), 5);
    expect_same_ordered_rule(make_triangle_rule(0, high, low), rule);
    EXPECT_EQ(make_triangle_rule(0, low, high).polynomial_exactness(), 4);
    expect_invalid_argument_with_message([&] { (void)make_triangle_rule(5, low, high); }, "exactness");
    expect_invalid_argument_with_message([&] { (void)make_triangle_rule(6, high, low); }, "exactness");
    const QuadratureRule zero(svmp::CellFamily::Line, 0, {{0.0, 0.0, 0.0}}, {2.0});
    expect_invalid_argument_with_message([&] { (void)make_triangle_rule(0, zero, high); }, "exactness");
}

TEST(TriangleDuffyImplementation, BoundsTheProductForUnequalFactorSizes)
{
    const auto single = make_gauss_legendre_rule(1);
    const auto two = make_gauss_legendre_rule(3);
    // Non-dyadic weights also exercise accumulation error near the capacity.
    for (const std::size_t n : {10000u, 16383u, 16384u}) {
        SCOPED_TRACE(::testing::Message() << "factor points=" << n);
        const QuadratureRule many(svmp::CellFamily::Line, 1,
            std::vector<QuadPoint>(n, QuadPoint::Zero()), std::vector<double>(n, 2.0 / n));
        const QuadratureRule expected(svmp::CellFamily::Triangle, 0,
            std::vector<QuadPoint>(n, QuadPoint{0.5, 0.25, 0.0}), std::vector<double>(n, 0.5 / n));
        expect_same_ordered_rule(make_triangle_rule(0, single, many), expected);
        expect_same_ordered_rule(make_triangle_rule(0, many, single), expected);
        EXPECT_THROW((void)make_triangle_rule(0, many, two), InvalidArgumentException);
        EXPECT_THROW((void)make_triangle_rule(0, two, many), InvalidArgumentException);
    }
}

TEST(TriangleDuffyImplementation, PreservesSignedZeroAndDuplicateSamples)
{
    const QuadratureRule factor(svmp::CellFamily::Line, 1,
        {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}, {3.0, -1.0, 0.0});
    expect_same_ordered_rule(make_triangle_rule(0, factor),
        QuadratureRule(svmp::CellFamily::Triangle, 0,
            std::vector<QuadPoint>(9, QuadPoint{0.5, 0.25, 0.0}),
            {1.125, -0.375, 0.0, -0.375, 0.125, 0.0, 0.0, 0.0, 0.0}));
}

TEST(TriangleDuffyImplementation, AcceptsFiniteProductsWithLargeSignedWeights)
{
    const QuadratureRule factor(svmp::CellFamily::Line, 1,
        {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}, {3.0e154, -3.0e154, 2.0});
    // Each product is finite, but their absolute sum exceeds binary64 range.
    const auto rule = make_triangle_rule(0, factor);
    ASSERT_EQ(rule.num_points(), 9u);
    for (double weight : rule.weights()) {
        EXPECT_TRUE(std::isfinite(weight));
    }
    EXPECT_EQ(rule.weight(0), -rule.weight(1));
    EXPECT_EQ(rule.weight(3), -rule.weight(4));
    EXPECT_EQ(rule.weight(2), -rule.weight(5));
    EXPECT_EQ(rule.weight(6), -rule.weight(7));
    EXPECT_EQ(rule.weight(8), 0.5);
}

TEST(TriangleDuffyImplementation, RejectsInvalidRequestsBeforeDegreeArithmetic)
{
    const auto factor = make_gauss_legendre_rule(255);
    for (const int p : {std::numeric_limits<int>::min(), -1, 255, std::numeric_limits<int>::max()}) {
        expect_invalid_argument_with_message([&] { (void)make_triangle_rule(p, factor); }, "requested_exactness");
        expect_invalid_argument_with_message([&] { (void)make_triangle_rule(p, factor, factor); }, "requested_exactness");
    }
}

TEST(TriangleDuffyImplementation, RejectsInvalidFactorContracts)
{
    const auto qualified = make_gauss_legendre_rule(255);
    const auto reject = [&](const QuadratureRule& invalid, std::string_view condition) {
        expect_invalid_argument_with_message([&] { (void)make_triangle_rule(0, invalid); }, condition);
        expect_invalid_argument_with_message([&] { (void)make_triangle_rule(0, invalid, qualified); }, condition);
        expect_invalid_argument_with_message([&] { (void)make_triangle_rule(0, qualified, invalid); }, condition);
    };
    const auto insufficient = make_gauss_legendre_rule(1);
    const auto snapshot = insufficient;
    expect_invalid_argument_with_message([&] { (void)make_triangle_rule(1, insufficient); }, "exactness");
    expect_same_ordered_rule(insufficient, snapshot);
    const QuadratureRule degree_zero(svmp::CellFamily::Line, 0, {{0.0, 0.0, 0.0}}, {2.0});
    expect_invalid_argument_with_message([&] { (void)make_triangle_rule(0, degree_zero); }, "exactness");
    reject(make_triangle_centroid_rule(), "Line");
    for (const double x : {-1.0001, 1.0001}) {
        const QuadratureRule exterior(svmp::CellFamily::Line, 1, {{x, 0.0, 0.0}}, {2.0});
        reject(exterior, "[-1,1]");
    }
    const QuadratureRule wrong_measure(svmp::CellFamily::Line, 1, {{0.0, 0.0, 0.0}}, {3.0});
    reject(wrong_measure, "normalization");
    auto emptied = make_gauss_legendre_rule(1);
    const auto retained = std::move(emptied);
    ASSERT_EQ(emptied.num_points(), 0u);
    reject(emptied, "empty");
    EXPECT_EQ(retained.num_points(), 1u);
    const QuadratureRule oversized(svmp::CellFamily::Line, 255,
        std::vector<QuadPoint>(129, QuadPoint::Zero()), std::vector<double>(129, 2.0 / 129.0));
    reject(oversized, "16384");
    const auto largest_lobatto = make_gauss_lobatto_rule(253);
    expect_invalid_argument_with_message([&] { (void)make_triangle_rule(253, largest_lobatto); }, "exactness");
}

TEST(TriangleDuffyImplementation, ReportsNumericalFailureWithoutRepairingSamples)
{
    const QuadratureRule false_metadata(svmp::CellFamily::Line, 1, {{0.5, 0.0, 0.0}}, {2.0});
    EXPECT_THROW((void)make_triangle_rule(0, false_metadata), ConvergenceException);
    const QuadratureRule overflow(svmp::CellFamily::Line, 1,
        {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}, {1.0e200, -1.0e200, 2.0});
    EXPECT_THROW((void)make_triangle_rule(0, overflow), ConvergenceException);
}

TEST(TriangleDuffyQualification, IntegratesEveryGaussBackedRequest)
{
    expect_every_triangle_rule(LineEndpointPolicy::Excluded);
}

TEST(TriangleDuffyQualification, IntegratesEveryLobattoBackedRequest)
{
    expect_every_triangle_rule(LineEndpointPolicy::Included);
}

TEST(TriangleDuffyQualification, IntegratesIndependentGaussAndLobattoFactors)
{
    struct Pair {
        bool u_lobatto;
        bool v_lobatto;
        int maximum_degree;
        std::size_t degree_five_points;
    };
    for (const auto pair : {Pair{false, false, 254, 12}, Pair{false, true, 253, 16},
                           Pair{true, false, 252, 15}, Pair{true, true, 252, 20}}) {
        for (const int p : {0, 1, 2, 5, 31, 127, pair.maximum_degree}) {
            SCOPED_TRACE(::testing::Message() << "degree=" << p
                << ", u Lobatto=" << pair.u_lobatto << ", v Lobatto=" << pair.v_lobatto);
            const auto u = pair.u_lobatto ? make_gauss_lobatto_rule(p + 1) : make_gauss_legendre_rule(p + 1);
            const auto v = pair.v_lobatto ? make_gauss_lobatto_rule(p) : make_gauss_legendre_rule(p);
            const auto rule = make_triangle_rule(p, u, v);
            ASSERT_EQ(rule.polynomial_exactness(), p);
            EXPECT_EQ(rule.num_points(), u.num_points() * v.num_points());
            if (p == 5) {
                EXPECT_EQ(rule.num_points(), pair.degree_five_points);
            }
            expect_triangle_moments(rule, p);
            if (p == pair.maximum_degree) {
                EXPECT_THROW((void)make_triangle_rule(p + 1, u, v), InvalidArgumentException);
            }
        }
    }
}
