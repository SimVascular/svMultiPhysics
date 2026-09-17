// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file test_QuadratureGenerators.cpp
 * @brief Tests for exactness-requested one-dimensional quadrature generators.
 */

#include <gtest/gtest.h>

#include "FE/Basis/NodeOrderingConventions.h"
#include "FE/Common/FEException.h"
#include "FE/Quadrature/GaussLobattoQuadrature.h"
#include "FE/Quadrature/GaussQuadrature.h"
#include "FE/Quadrature/QuadratureRule.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <span>
#include <string_view>
#include <type_traits>
#include <utility>

using namespace svmp::FE;
using namespace svmp::FE::quadrature;

namespace {

// Exhaustive supported-domain sweeps observed structure, measure, and moment
// errors below 1.2e-15; these envelopes retain cross-toolchain headroom
// through degree 255.
constexpr double kStructureTolerance = 1.0e-12;
constexpr double kMomentTolerance = 2.0e-12;
// Computed fixtures follow the production 64*epsilon refinement scale.
constexpr double kFixtureTolerance =
    64.0 * std::numeric_limits<double>::epsilon();
// Basis and Quadrature independently refine computed interior GLL roots, so
// terminal binary64 rounding can differ; this remains far below node spacing.
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

long double analytic_line_monomial_integral(std::size_t power)
{
    if (power % 2u != 0u) {
        return 0.0L;
    }
    return 2.0L / (static_cast<long double>(power) + 1.0L);
}

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
            analytic_line_monomial_integral(nonnegative_power));
        EXPECT_LE(error, static_cast<long double>(kMomentTolerance));
    }
}

void expect_every_supported_line_rule(
    QuadratureRule (*generator)(int),
    int first_num_points,
    int last_num_points,
    int exactness_subtrahend,
    LineEndpointPolicy endpoint_policy)
{
    for (int num_points = first_num_points;
         num_points <= last_num_points;
         ++num_points) {
        SCOPED_TRACE(
            ::testing::Message() << "num_points=" << num_points);
        const int expected_exactness = 2 * num_points - exactness_subtrahend;
        // Both requests must select this minimum point count and report its
        // actual exactness, not just echo the requested degree.
        for (const int requested_exactness :
             {expected_exactness - 1, expected_exactness}) {
            SCOPED_TRACE(
                ::testing::Message() << "requested_exactness=" << requested_exactness);
            const QuadratureRule rule = generator(requested_exactness);

            expect_common_line_metadata(
                rule, static_cast<std::size_t>(num_points), expected_exactness);
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
        analytic_line_monomial_integral(first_unadvertised_even_power));
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

    EXPECT_EQ(analytic_line_monomial_integral(0u), 2.0L);
    EXPECT_EQ(analytic_line_monomial_integral(1u), 0.0L);
    EXPECT_EQ(
        analytic_line_monomial_integral(2u),
        2.0L / 3.0L);
    EXPECT_EQ(analytic_line_monomial_integral(255u), 0.0L);
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

TEST(GaussLegendreImplementation, GeneratesEverySupportedRule)
{
    expect_every_supported_line_rule(
        &make_gauss_legendre_rule,
        1,
        128,
        1,
        LineEndpointPolicy::Excluded);
}

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

TEST(GaussLobattoImplementation, GeneratesEverySupportedRule)
{
    expect_every_supported_line_rule(
        &make_gauss_lobatto_rule,
        2,
        128,
        3,
        LineEndpointPolicy::Included);
}

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

TEST(GaussLobattoBasisConsistency, MatchesRepresentativeNodeDistributions)
{
    constexpr std::array point_counts{2, 4, 65, 128};

    for (const int num_points : point_counts) {
        SCOPED_TRACE(
            ::testing::Message() << "num_points=" << num_points);
        const QuadratureRule rule =
            make_gauss_lobatto_rule(2 * num_points - 3);
        ASSERT_EQ(
            rule.num_points(),
            static_cast<std::size_t>(num_points));

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
                    num_points - 1);

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
