/**
 * @file test_QuadratureRules.cpp
 * @brief Unit tests for the core quadrature rule infrastructure.
 */

#include <gtest/gtest.h>

#include "FE/Common/FEException.h"
#include "FE/Quadrature/QuadratureRule.h"

#include <array>
#include <cmath>
#include <exception>
#include <limits>
#include <span>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace svmp::FE;
using namespace svmp::FE::quadrature;

namespace {

using ExpectedPoint = std::array<double, 3>;

QuadratureRule make_two_point_gauss_legendre_rule()
{
    const double abscissa = 1.0 / std::sqrt(3.0);
    return QuadratureRule(
        svmp::CellFamily::Line,
        3,
        {{-abscissa, 0.0, 0.0}, {abscissa, 0.0, 0.0}},
        {1.0, 1.0});
}

template <typename Function>
void expect_invalid_argument_with_message(
    Function&& function,
    const std::string& expected_message)
{
    try {
        std::forward<Function>(function)();
        FAIL() << "Expected InvalidArgumentException containing: "
               << expected_message;
    } catch (const InvalidArgumentException& exception) {
        const std::string actual_message = exception.what();
        EXPECT_NE(actual_message.find(expected_message), std::string::npos)
            << "actual message: " << actual_message;
    } catch (const std::exception& exception) {
        FAIL() << "Expected InvalidArgumentException, received: "
               << exception.what();
    } catch (...) {
        FAIL() << "Expected InvalidArgumentException, received an unknown exception";
    }
}

} // namespace

TEST(QuadPointContract, UsesFixedSizeFEVectorRepresentation)
{
    static_assert(std::is_same_v<QuadPoint, math::Vector<double, 3>>);
    static_assert(QuadPoint::RowsAtCompileTime == 3);
    static_assert(QuadPoint::ColsAtCompileTime == 1);

    const QuadPoint origin = QuadPoint::Zero();
    for (std::size_t component = 0; component < 3u; ++component) {
        EXPECT_DOUBLE_EQ(origin[component], 0.0);
    }

    const QuadPoint line_point{0.25, 0.0, 0.0};
    EXPECT_DOUBLE_EQ(line_point[0], 0.25);
    EXPECT_DOUBLE_EQ(line_point[1], 0.0);
    EXPECT_DOUBLE_EQ(line_point[2], 0.0);

    const QuadPoint surface_point{0.25, 0.5, 0.0};
    EXPECT_DOUBLE_EQ(surface_point[0], 0.25);
    EXPECT_DOUBLE_EQ(surface_point[1], 0.5);
    EXPECT_DOUBLE_EQ(surface_point[2], 0.0);

    QuadPoint mutable_point = QuadPoint::Zero();
    mutable_point[2] = 0.75;
    EXPECT_DOUBLE_EQ(mutable_point[2], 0.75);

    const std::vector<QuadPoint> points(2, QuadPoint::Zero());
    for (const auto& point : points) {
        for (std::size_t component = 0; component < 3u; ++component) {
            EXPECT_DOUBLE_EQ(point[component], 0.0);
        }
    }
}

TEST(QuadratureRuleValidation, AcceptsEverySupportedReferenceCell)
{
    struct Case {
        svmp::CellFamily family;
        std::size_t expected_dimension;
        double expected_measure;
        ExpectedPoint point;
    };

    const std::vector<Case> cases = {
        {svmp::CellFamily::Point, 0, 1.0, {0.0, 0.0, 0.0}},
        {svmp::CellFamily::Line, 1, 2.0, {0.0, 0.0, 0.0}},
        {svmp::CellFamily::Triangle, 2, 0.5, {0.25, 0.25, 0.0}},
        {svmp::CellFamily::Quad, 2, 4.0, {0.0, 0.0, 0.0}},
        {svmp::CellFamily::Tetra, 3, 1.0 / 6.0, {0.25, 0.25, 0.25}},
        {svmp::CellFamily::Hex, 3, 8.0, {0.0, 0.0, 0.0}},
        {svmp::CellFamily::Wedge, 3, 1.0, {0.25, 0.25, 0.0}},
    };

    for (const auto& c : cases) {
        const QuadratureRule rule(
            c.family,
            0,
            {{c.point[0], c.point[1], c.point[2]}},
            {c.expected_measure});
        EXPECT_EQ(rule.dimension(), c.expected_dimension);
        EXPECT_DOUBLE_EQ(
            rule.reference_cell_measure(),
            c.expected_measure);
    }
}

TEST(QuadratureRuleValidation, RejectsInvalidMetadata)
{
    expect_invalid_argument_with_message(
        [] {
            (void)QuadratureRule(
                svmp::CellFamily::Triangle, -1, {{0.0, 0.0, 0.0}}, {0.5});
        },
        "polynomial exactness must be non-negative");

    const std::array<svmp::CellFamily, 3> unsupported_families = {
        svmp::CellFamily::Pyramid,
        svmp::CellFamily::Polygon,
        svmp::CellFamily::Polyhedron,
    };
    for (const auto family : unsupported_families) {
        SCOPED_TRACE(static_cast<int>(family));
        expect_invalid_argument_with_message(
            [family] {
                (void)QuadratureRule(
                    family, 1, {{0.0, 0.0, 0.0}}, {1.0});
            },
            "unsupported reference-cell family");
    }

    expect_invalid_argument_with_message(
        [] {
            (void)QuadratureRule(
                static_cast<svmp::CellFamily>(255),
                1,
                {{0.0, 0.0, 0.0}},
                {1.0});
        },
        "unsupported reference-cell family");
}

TEST(QuadratureRuleValidation, RejectsMalformedStorageAndNonfiniteValues)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();

    expect_invalid_argument_with_message(
        [] { (void)QuadratureRule(svmp::CellFamily::Line, 1, {}, {}); },
        "at least one point");
    expect_invalid_argument_with_message(
        [] {
            (void)QuadratureRule(
                svmp::CellFamily::Line, 1, {{0.0, 0.0, 0.0}}, {});
        },
        "points/weights size mismatch");
    expect_invalid_argument_with_message(
        [nan] {
            (void)QuadratureRule(
                svmp::CellFamily::Line, 1, {{nan, 0.0, 0.0}}, {2.0});
        },
        "non-finite coordinate at point index 0");
    expect_invalid_argument_with_message(
        [inf] {
            (void)QuadratureRule(
                svmp::CellFamily::Line, 1, {{0.0, 0.0, 0.0}}, {inf});
        },
        "quadrature weight must be finite at point index 0");
    expect_invalid_argument_with_message(
        [nan] {
            (void)QuadratureRule(
                svmp::CellFamily::Line, 1, {{0.0, 0.0, 0.0}}, {nan});
        },
        "quadrature weight must be finite at point index 0");
}

TEST(QuadratureRuleValidation, RejectsAnyNonzeroInactiveCoordinate)
{
    constexpr double nonzero = std::numeric_limits<double>::epsilon();

    expect_invalid_argument_with_message(
        [nonzero] {
            (void)QuadratureRule(
                svmp::CellFamily::Point, 0, {{nonzero, 0.0, 0.0}}, {1.0});
        },
        "nonzero inactive coordinate at point index 0");
    expect_invalid_argument_with_message(
        [nonzero] {
            (void)QuadratureRule(
                svmp::CellFamily::Line, 1, {{0.0, -nonzero, 0.0}}, {2.0});
        },
        "nonzero inactive coordinate at point index 0");
}

TEST(QuadratureRuleValidation, AllowsZeroAndNegativeWeights)
{
    const QuadratureRule rule(
        svmp::CellFamily::Triangle,
        0,
        {{1.0 / 3.0, 1.0 / 3.0, 0.0},
         {0.2, 0.2, 0.0},
         {0.1, 0.1, 0.0}},
        {-0.25, 0.0, 0.75});
    EXPECT_LT(rule.weight(0), 0.0);
    EXPECT_DOUBLE_EQ(rule.weight(1), 0.0);
}

TEST(QuadratureRuleValidation, AllowsExteriorActiveCoordinates)
{
    const QuadratureRule rule(
        svmp::CellFamily::Line,
        0,
        {{2.0, 0.0, 0.0}},
        {2.0});
    EXPECT_DOUBLE_EQ(rule.point(0)[0], 2.0);
}

TEST(QuadratureRuleValidation, LeavesWeightNormalizationToGenerators)
{
    const QuadratureRule rule(
        svmp::CellFamily::Triangle,
        0,
        {{0.25, 0.25, 0.0}},
        {1.0});
    EXPECT_DOUBLE_EQ(rule.weight(0), 1.0);
    EXPECT_DOUBLE_EQ(rule.reference_cell_measure(), 0.5);
}

TEST(QuadratureRuleContract, SupportsValueSemanticsAndReadOnlyQueries)
{
    static_assert(!std::is_abstract_v<QuadratureRule>);
    static_assert(std::is_final_v<QuadratureRule>);
    static_assert(!std::has_virtual_destructor_v<QuadratureRule>);
    static_assert(std::is_copy_constructible_v<QuadratureRule>);
    static_assert(std::is_move_constructible_v<QuadratureRule>);
    static_assert(std::is_copy_assignable_v<QuadratureRule>);
    static_assert(std::is_move_assignable_v<QuadratureRule>);
    static_assert(
        std::is_same<decltype(std::declval<QuadratureRule&>().point(0)),
                     const QuadPoint&>::value,
        "A quadrature point must be exposed through a const reference");
    static_assert(
        std::is_same<decltype(std::declval<QuadratureRule&>().points()),
                     std::span<const QuadPoint>>::value,
        "Quadrature points must be exposed through a read-only view");
    static_assert(
        std::is_same<decltype(std::declval<QuadratureRule&>().weights()),
                     std::span<const double>>::value,
        "Quadrature weights must be exposed through a read-only view");

    const double a = 1.0 / std::sqrt(3.0);
    const QuadratureRule rule = make_two_point_gauss_legendre_rule();

    EXPECT_EQ(rule.cell_family(), svmp::CellFamily::Line);
    EXPECT_EQ(rule.dimension(), 1);
    EXPECT_EQ(rule.polynomial_exactness(), 3);
    EXPECT_DOUBLE_EQ(rule.reference_cell_measure(), 2.0);
    ASSERT_EQ(rule.num_points(), 2u);
    ASSERT_EQ(rule.points().size(), 2u);
    ASSERT_EQ(rule.weights().size(), 2u);
    EXPECT_EQ(rule.points().data(), &rule.point(0));
    EXPECT_DOUBLE_EQ(rule.point(0)[0], -a);
    EXPECT_DOUBLE_EQ(rule.point(1)[0], a);
    EXPECT_DOUBLE_EQ(rule.weight(0), 1.0);
    EXPECT_DOUBLE_EQ(rule.weight(1), 1.0);
}
