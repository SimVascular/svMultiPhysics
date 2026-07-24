/**
 * @file test_QuadratureRules.cpp
 * @brief Unit tests for the core quadrature rule infrastructure.
 */

#include <gtest/gtest.h>

#include "FE/Common/FEException.h"
#include "FE/Quadrature/QuadratureRule.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace svmp::FE;
using namespace svmp::FE::quadrature;

namespace {

class RuleProbe final : public QuadratureRule {
public:
    RuleProbe(svmp::CellFamily family,
              int polynomial_exactness,
              std::vector<QuadPoint> points,
              std::vector<double> weights)
        : QuadratureRule(
              family,
              RuleData{
                  polynomial_exactness,
                  std::move(points),
                  std::move(weights)})
    {
    }
};

constexpr double kTol = 1.0e-12;

using ExpectedPoint = std::array<double, 3>;

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

double weight_sum(const QuadratureRule& rule)
{
    double sum = 0.0;
    for (const double weight : rule.weights()) {
        sum += weight;
    }
    return sum;
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
        int expected_dimension;
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
        const RuleProbe rule(
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
            (void)RuleProbe(
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
                (void)RuleProbe(family, 1, {{0.0, 0.0, 0.0}}, {1.0});
            },
            "unsupported reference-cell family");
    }

    expect_invalid_argument_with_message(
        [] {
            (void)RuleProbe(
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
        [] { (void)RuleProbe(svmp::CellFamily::Line, 1, {}, {}); },
        "at least one sample");
    expect_invalid_argument_with_message(
        [] {
            (void)RuleProbe(
                svmp::CellFamily::Line, 1, {{0.0, 0.0, 0.0}}, {});
        },
        "points/weights size mismatch");
    expect_invalid_argument_with_message(
        [nan] {
            (void)RuleProbe(
                svmp::CellFamily::Line, 1, {{nan, 0.0, 0.0}}, {2.0});
        },
        "non-finite coordinate at sample 0");
    expect_invalid_argument_with_message(
        [inf] {
            (void)RuleProbe(
                svmp::CellFamily::Line, 1, {{0.0, 0.0, 0.0}}, {inf});
        },
        "quadrature weight must be finite at sample 0");
    expect_invalid_argument_with_message(
        [nan] {
            (void)RuleProbe(
                svmp::CellFamily::Line, 1, {{0.0, 0.0, 0.0}}, {nan});
        },
        "quadrature weight must be finite at sample 0");
}

TEST(QuadratureRuleValidation, RejectsInactiveCoordinatesAndOutOfCellPoints)
{
    expect_invalid_argument_with_message(
        [] {
            (void)RuleProbe(
                svmp::CellFamily::Point, 0, {{1.0e-4, 0.0, 0.0}}, {1.0});
        },
        "nonzero inactive coordinate at sample 0");
    expect_invalid_argument_with_message(
        [] {
            (void)RuleProbe(
                svmp::CellFamily::Line, 1, {{0.0, 1.0e-4, 0.0}}, {2.0});
        },
        "nonzero inactive coordinate at sample 0");
    expect_invalid_argument_with_message(
        [] {
            (void)RuleProbe(
                svmp::CellFamily::Triangle, 1, {{0.8, 0.3, 0.0}}, {0.5});
        },
        "outside the canonical reference cell at sample 0");
    EXPECT_THROW(
        (void)RuleProbe(svmp::CellFamily::Quad, 1, {{0.0, 1.1, 0.0}}, {4.0}),
        InvalidArgumentException);
    EXPECT_THROW(
        (void)RuleProbe(svmp::CellFamily::Tetra, 1, {{0.4, 0.4, 0.3}}, {1.0 / 6.0}),
        InvalidArgumentException);
    EXPECT_THROW(
        (void)RuleProbe(svmp::CellFamily::Hex, 1, {{0.0, 0.0, -1.1}}, {8.0}),
        InvalidArgumentException);
    EXPECT_THROW(
        (void)RuleProbe(svmp::CellFamily::Wedge, 1, {{0.6, 0.5, 0.0}}, {1.0}),
        InvalidArgumentException);
}

TEST(QuadratureRuleValidation, EnforcesReferenceCellMeasureButAllowsNegativeWeights)
{
    expect_invalid_argument_with_message(
        [] {
            (void)RuleProbe(
                svmp::CellFamily::Triangle,
                1,
                {{0.25, 0.25, 0.0}},
                {1.0});
        },
        "weights do not reproduce the reference-cell measure");

    const RuleProbe rule(
        svmp::CellFamily::Triangle,
        0,
        {{1.0 / 3.0, 1.0 / 3.0, 0.0}, {0.2, 0.2, 0.0}},
        {-0.25, 0.75});
    EXPECT_LT(rule.weight(0), 0.0);
    EXPECT_DOUBLE_EQ(
        weight_sum(rule),
        rule.reference_cell_measure());
}

TEST(QuadratureRuleValidation, RejectsIncorrectMeasureDespiteLargeCancellation)
{
    expect_invalid_argument_with_message(
        [] {
            (void)RuleProbe(
                svmp::CellFamily::Line,
                0,
                {{-0.5, 0.0, 0.0},
                 {0.0, 0.0, 0.0},
                 {0.5, 0.0, 0.0}},
                {1.0e20, 0.0, -1.0e20});
        },
        "weights do not reproduce the reference-cell measure");
}

TEST(QuadratureRuleValidation, RejectsIncorrectMeasureHiddenByCancellation)
{
    const double medium = std::ldexp(1.0, 93);
    const double large = std::ldexp(1.0, 233);
    const double residual = std::ldexp(1.0, 14);

    expect_invalid_argument_with_message(
        [medium, large, residual] {
            (void)RuleProbe(
                svmp::CellFamily::Line,
                0,
                std::vector<QuadPoint>(6u, QuadPoint::Zero()),
                {-medium, -large, residual, medium, large, 2.0});
        },
        "weights do not reproduce the reference-cell measure");
}

TEST(QuadratureRuleValidation, ExactMeasureValidationIsOrderIndependent)
{
    const double maximum = std::numeric_limits<double>::max();
    std::vector<double> valid_weights{
        -maximum, -maximum, 2.0, maximum, maximum};
    std::sort(valid_weights.begin(), valid_weights.end());

    std::size_t permutation_count = 0u;
    do {
        SCOPED_TRACE(permutation_count);
        EXPECT_NO_THROW(
            (void)RuleProbe(
                svmp::CellFamily::Line,
                0,
                std::vector<QuadPoint>(
                    valid_weights.size(),
                    QuadPoint::Zero()),
                valid_weights));
        ++permutation_count;
    } while (std::next_permutation(
        valid_weights.begin(),
        valid_weights.end()));
    EXPECT_EQ(permutation_count, 30u);

    const double excessive_residual = std::ldexp(1.0, -38);
    std::vector<double> invalid_weights{
        -maximum,
        -maximum,
        2.0,
        excessive_residual,
        maximum,
        maximum};
    std::sort(invalid_weights.begin(), invalid_weights.end());

    permutation_count = 0u;
    do {
        SCOPED_TRACE(permutation_count);
        expect_invalid_argument_with_message(
            [&invalid_weights] {
                (void)RuleProbe(
                    svmp::CellFamily::Line,
                    0,
                    std::vector<QuadPoint>(
                        invalid_weights.size(),
                        QuadPoint::Zero()),
                    invalid_weights);
            },
            "weights do not reproduce the reference-cell measure");
        ++permutation_count;
    } while (std::next_permutation(
        invalid_weights.begin(),
        invalid_weights.end()));
    EXPECT_EQ(permutation_count, 180u);
}

TEST(QuadratureRuleValidation, HandlesExtremeAndSubnormalCancellation)
{
    const double maximum = std::numeric_limits<double>::max();
    const double minimum_normal = std::numeric_limits<double>::min();
    const double minimum_subnormal =
        std::numeric_limits<double>::denorm_min();
    const double next_normal =
        std::nextafter(minimum_normal, std::numeric_limits<double>::infinity());

    EXPECT_NO_THROW(
        (void)RuleProbe(
            svmp::CellFamily::Line,
            0,
            std::vector<QuadPoint>(6u, QuadPoint::Zero()),
            {
                maximum,
                minimum_normal,
                minimum_subnormal,
                2.0,
                -next_normal,
                -maximum,
            }));

    expect_invalid_argument_with_message(
        [maximum] {
            (void)RuleProbe(
                svmp::CellFamily::Line,
                0,
                std::vector<QuadPoint>(2u, QuadPoint::Zero()),
                {maximum, maximum});
        },
        "weights do not reproduce the reference-cell measure");
}

TEST(QuadratureRuleValidation, AppliesToleranceToTheExactWeightSum)
{
    const double maximum = std::numeric_limits<double>::max();
    const double accepted_residual = std::ldexp(1.0, -42);
    const double rejected_residual = std::ldexp(1.0, -38);

    EXPECT_NO_THROW(
        (void)RuleProbe(
            svmp::CellFamily::Line,
            0,
            std::vector<QuadPoint>(4u, QuadPoint::Zero()),
            {maximum, accepted_residual, -maximum, 2.0}));

    expect_invalid_argument_with_message(
        [maximum, rejected_residual] {
            (void)RuleProbe(
                svmp::CellFamily::Line,
                0,
                std::vector<QuadPoint>(4u, QuadPoint::Zero()),
                {maximum, rejected_residual, -maximum, 2.0});
        },
        "weights do not reproduce the reference-cell measure");
}

TEST(QuadratureRuleValidation, AcceptsLargePositiveRule)
{
    // Repeating a non-dyadic weight exposes ordinary left-fold drift without
    // changing the mathematical normalization of the rule.
    constexpr std::size_t sample_count = 100000u;
    const double sample_weight = 2.0 / static_cast<double>(sample_count);
    std::vector<QuadPoint> points(sample_count, QuadPoint::Zero());
    std::vector<double> weights(sample_count, sample_weight);

    const RuleProbe rule(
        svmp::CellFamily::Line,
        0,
        std::move(points),
        std::move(weights));
    EXPECT_EQ(rule.num_points(), sample_count);
    EXPECT_DOUBLE_EQ(rule.weight(0), sample_weight);
}

TEST(QuadratureRuleValidation, AcceptsLargeSignedRule)
{
    constexpr std::size_t sample_count = 8192u;
    constexpr double cancelling_weight = 1.0 / 1048576.0;
    std::vector<QuadPoint> points(sample_count, QuadPoint::Zero());
    std::vector<double> weights(
        sample_count,
        2.0 / static_cast<double>(sample_count - 2u));
    weights[0] = -cancelling_weight;
    weights[1] = cancelling_weight;

    const RuleProbe rule(
        svmp::CellFamily::Line,
        0,
        std::move(points),
        std::move(weights));
    EXPECT_LT(rule.weight(0), 0.0);
    EXPECT_NEAR(
        weight_sum(rule),
        rule.reference_cell_measure(),
        kTol);
}

TEST(QuadratureRuleValidation, AppliesConstructionCoordinateTolerance)
{
    constexpr double accepted_offset = 1.0e-14;
    const RuleProbe rule(
        svmp::CellFamily::Line,
        0,
        {{1.0 + accepted_offset, accepted_offset, 0.0}},
        {2.0});
    EXPECT_DOUBLE_EQ(rule.point(0)[0], 1.0 + accepted_offset);

    constexpr double rejected_offset = 1.0e-10;
    expect_invalid_argument_with_message(
        [rejected_offset] {
            (void)RuleProbe(
                svmp::CellFamily::Line,
                0,
                {{1.0 + rejected_offset, 0.0, 0.0}},
                {2.0});
        },
        "outside the canonical reference cell");
}

TEST(QuadratureRuleValidation, AppliesConstructionWeightTolerance)
{
    constexpr double accepted_offset = 1.0e-14;
    const RuleProbe rule(
        svmp::CellFamily::Triangle,
        0,
        {{0.25, 0.25, 0.0}},
        {0.5 + accepted_offset});
    EXPECT_DOUBLE_EQ(rule.weight(0), 0.5 + accepted_offset);

    constexpr double rejected_offset = 1.0e-10;
    expect_invalid_argument_with_message(
        [rejected_offset] {
            (void)RuleProbe(
                svmp::CellFamily::Triangle,
                0,
                {{0.25, 0.25, 0.0}},
                {0.5 + rejected_offset});
        },
        "weights do not reproduce the reference-cell measure");
}

TEST(QuadratureRuleContract, PublishesOnlyACompleteImmutableQueryInterface)
{
    static_assert(std::is_abstract_v<QuadratureRule>);
    static_assert(!std::is_copy_constructible_v<QuadratureRule>);
    static_assert(!std::is_move_constructible_v<QuadratureRule>);
    static_assert(!std::is_copy_assignable_v<QuadratureRule>);
    static_assert(!std::is_move_assignable_v<QuadratureRule>);
    static_assert(!std::is_copy_constructible_v<RuleProbe>);
    static_assert(!std::is_move_constructible_v<RuleProbe>);
    static_assert(!std::is_copy_assignable_v<RuleProbe>);
    static_assert(!std::is_move_assignable_v<RuleProbe>);
    static_assert(
        std::is_same<decltype(std::declval<QuadratureRule&>().point(0)),
                     const QuadPoint&>::value,
        "A quadrature point must be exposed through an immutable reference");
    static_assert(
        std::is_same<decltype(std::declval<QuadratureRule&>().points()),
                     const std::vector<QuadPoint>&>::value,
        "Quadrature points must be exposed through an immutable view");
    static_assert(
        std::is_same<decltype(std::declval<QuadratureRule&>().weights()),
                     const std::vector<double>&>::value,
        "Quadrature weights must be exposed through an immutable view");

    const double a = 1.0 / std::sqrt(3.0);
    const RuleProbe rule(
        svmp::CellFamily::Line,
        3,
        {{-a, 0.0, 0.0}, {a, 0.0, 0.0}},
        {1.0, 1.0});

    EXPECT_EQ(rule.cell_family(), svmp::CellFamily::Line);
    EXPECT_EQ(rule.dimension(), 1);
    EXPECT_EQ(rule.polynomial_exactness(), 3);
    EXPECT_DOUBLE_EQ(rule.reference_cell_measure(), 2.0);
    ASSERT_EQ(rule.num_points(), 2u);
    ASSERT_EQ(rule.points().size(), 2u);
    ASSERT_EQ(rule.weights().size(), 2u);
    EXPECT_DOUBLE_EQ(rule.point(0)[0], -a);
    EXPECT_DOUBLE_EQ(rule.point(1)[0], a);
    EXPECT_DOUBLE_EQ(rule.weight(0), 1.0);
    EXPECT_DOUBLE_EQ(rule.weight(1), 1.0);
}
