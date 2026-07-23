/**
 * @file test_QuadratureRules.cpp
 * @brief Baseline solver data and core quadrature rule contract tests.
 */

#include <gtest/gtest.h>

#include "FE/Common/FEException.h"
#include "FE/Quadrature/QuadratureRule.h"
#include "nn.h"

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
constexpr double kExactnessTol = 3.0e-12;

using ExpectedPoint = std::array<double, 3>;

struct ExpectedSample {
    ExpectedPoint point;
    double weight;
};

using ExpectedSamples = std::vector<ExpectedSample>;

struct SolverRuleCase {
    const char* name;
    consts::ElementType solver_type;
    ElementType fe_type;
    int dimension;
    int requested_exactness;
    int advertised_exactness;
    double reference_measure;
    ExpectedSamples samples;
};

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

ExpectedSamples point_samples()
{
    return {{{0.0, 0.0, 0.0}, 1.0}};
}

ExpectedSamples line2_samples()
{
    const double a = 1.0 / std::sqrt(3.0);
    return {
        {{-a, 0.0, 0.0}, 1.0},
        {{ a, 0.0, 0.0}, 1.0},
    };
}

ExpectedSamples line3_samples()
{
    const double a = std::sqrt(0.6);
    return {
        {{-a, 0.0, 0.0}, 5.0 / 9.0},
        {{ a, 0.0, 0.0}, 5.0 / 9.0},
        {{0.0, 0.0, 0.0}, 8.0 / 9.0},
    };
}

ExpectedSamples triangle3_samples()
{
    const double s = 2.0 / 3.0;
    const double t = 1.0 / 6.0;
    return {
        {{t, t, 0.0}, 1.0 / 6.0},
        {{s, t, 0.0}, 1.0 / 6.0},
        {{t, s, 0.0}, 1.0 / 6.0},
    };
}

ExpectedSamples triangle6_samples()
{
    const double centroid = 0.333333333333333;
    const double a1 = 0.797426985353087;
    const double b1 = 0.101286507323456;
    const double a2 = 0.059715871789770;
    const double b2 = 0.470142064105115;
    const double w0 = 0.225000000000000 * 0.5;
    const double w1 = 0.125939180544827 * 0.5;
    const double w2 = 0.132394152788506 * 0.5;
    return {
        {{centroid, centroid, 0.0}, w0},
        {{a1, b1, 0.0}, w1},
        {{b1, a1, 0.0}, w1},
        {{b1, b1, 0.0}, w1},
        {{a2, b2, 0.0}, w2},
        {{b2, a2, 0.0}, w2},
        {{b2, b2, 0.0}, w2},
    };
}

ExpectedSamples quad4_samples()
{
    const double a = 1.0 / std::sqrt(3.0);
    return {
        {{-a, -a, 0.0}, 1.0},
        {{ a, -a, 0.0}, 1.0},
        {{ a,  a, 0.0}, 1.0},
        {{-a,  a, 0.0}, 1.0},
    };
}

ExpectedSamples quad9_samples(double a = std::sqrt(0.6))
{
    const double corner_weight = 25.0 / 81.0;
    const double edge_weight = 40.0 / 81.0;
    return {
        {{-a, -a, 0.0}, corner_weight},
        {{ a, -a, 0.0}, corner_weight},
        {{ a,  a, 0.0}, corner_weight},
        {{-a,  a, 0.0}, corner_weight},
        {{0.0, -a, 0.0}, edge_weight},
        {{ a, 0.0, 0.0}, edge_weight},
        {{0.0,  a, 0.0}, edge_weight},
        {{-a, 0.0, 0.0}, edge_weight},
        {{0.0, 0.0, 0.0}, 64.0 / 81.0},
    };
}

ExpectedSamples tetra4_samples()
{
    const double a = (5.0 + 3.0 * std::sqrt(5.0)) / 20.0;
    const double b = (1.0 - a) / 3.0;
    return {
        {{a, b, b}, 1.0 / 24.0},
        {{b, a, b}, 1.0 / 24.0},
        {{b, b, a}, 1.0 / 24.0},
        {{b, b, b}, 1.0 / 24.0},
    };
}

ExpectedSamples tetra10_samples()
{
    const double one_third = 0.3333333333333330;
    const double a1 = 0.0909090909090910;
    const double b1 = 0.7272727272727270;
    const double a2 = 0.0665501535736640;
    const double b2 = 0.4334498464263360;
    const double w0 = 0.0302836780970890;
    const double w1 = 0.0060267857142860;
    const double w2 = 0.0116452490860290;
    const double w3 = 0.0109491415613860;
    return {
        {{0.25, 0.25, 0.25}, w0},
        {{0.0, one_third, one_third}, w1},
        {{one_third, 0.0, one_third}, w1},
        {{one_third, one_third, 0.0}, w1},
        {{one_third, one_third, one_third}, w1},
        {{b1, a1, a1}, w2},
        {{a1, b1, a1}, w2},
        {{a1, a1, b1}, w2},
        {{a1, a1, a1}, w2},
        {{a2, a2, b2}, w3},
        {{a2, b2, a2}, w3},
        {{a2, b2, b2}, w3},
        {{b2, b2, a2}, w3},
        {{b2, a2, b2}, w3},
        {{b2, a2, a2}, w3},
    };
}

ExpectedSamples hex8_samples()
{
    const double a = 1.0 / std::sqrt(3.0);
    return {
        {{-a, -a, -a}, 1.0},
        {{ a, -a, -a}, 1.0},
        {{ a,  a, -a}, 1.0},
        {{-a,  a, -a}, 1.0},
        {{-a, -a,  a}, 1.0},
        {{ a, -a,  a}, 1.0},
        {{ a,  a,  a}, 1.0},
        {{-a,  a,  a}, 1.0},
    };
}

ExpectedSamples hex27_samples()
{
    const double a = std::sqrt(0.6);
    const double w0 = 125.0 / 729.0;
    const double w1 = 200.0 / 729.0;
    const double w2 = 320.0 / 729.0;
    return {
        {{-a, -a, -a}, w0}, {{ a, -a, -a}, w0},
        {{ a,  a, -a}, w0}, {{-a,  a, -a}, w0},
        {{-a, -a,  a}, w0}, {{ a, -a,  a}, w0},
        {{ a,  a,  a}, w0}, {{-a,  a,  a}, w0},

        {{0.0, -a, -a}, w1}, {{ a, 0.0, -a}, w1},
        {{0.0,  a, -a}, w1}, {{-a, 0.0, -a}, w1},
        {{0.0, -a,  a}, w1}, {{ a, 0.0,  a}, w1},
        {{0.0,  a,  a}, w1}, {{-a, 0.0,  a}, w1},
        {{-a, -a, 0.0}, w1}, {{ a, -a, 0.0}, w1},
        {{ a,  a, 0.0}, w1}, {{-a,  a, 0.0}, w1},

        {{-a, 0.0, 0.0}, w2}, {{ a, 0.0, 0.0}, w2},
        {{0.0, -a, 0.0}, w2}, {{0.0,  a, 0.0}, w2},
        {{0.0, 0.0, -a}, w2}, {{0.0, 0.0,  a}, w2},
        {{0.0, 0.0, 0.0}, 512.0 / 729.0},
    };
}

ExpectedSamples wedge6_samples()
{
    const double s = 2.0 / 3.0;
    const double t = 1.0 / 6.0;
    const double z = 1.0 / std::sqrt(3.0);
    return {
        {{s, t, -z}, 1.0 / 6.0},
        {{t, s, -z}, 1.0 / 6.0},
        {{t, t, -z}, 1.0 / 6.0},
        {{s, t,  z}, 1.0 / 6.0},
        {{t, s,  z}, 1.0 / 6.0},
        {{t, t,  z}, 1.0 / 6.0},
    };
}

const std::vector<SolverRuleCase>& standard_solver_cases()
{
    static const std::vector<SolverRuleCase> cases = {
        {"PNT", consts::ElementType::PNT, ElementType::Point1, 0, 1, 1, 1.0, point_samples()},
        {"LIN1", consts::ElementType::LIN1, ElementType::Line2, 1, 2, 3, 2.0, line2_samples()},
        {"LIN2", consts::ElementType::LIN2, ElementType::Line3, 1, 4, 5, 2.0, line3_samples()},
        {"TRI3", consts::ElementType::TRI3, ElementType::Triangle3, 2, 2, 2, 0.5, triangle3_samples()},
        {"TRI6", consts::ElementType::TRI6, ElementType::Triangle6, 2, 5, 5, 0.5, triangle6_samples()},
        {"QUD4", consts::ElementType::QUD4, ElementType::Quad4, 2, 2, 3, 4.0, quad4_samples()},
        {"QUD8", consts::ElementType::QUD8, ElementType::Quad8, 2, 4, 5, 4.0, quad9_samples()},
        {"QUD9", consts::ElementType::QUD9, ElementType::Quad9, 2, 4, 5, 4.0, quad9_samples()},
        {"TET4", consts::ElementType::TET4, ElementType::Tetra4, 3, 2, 2, 1.0 / 6.0, tetra4_samples()},
        {"TET10", consts::ElementType::TET10, ElementType::Tetra10, 3, 5, 5, 1.0 / 6.0, tetra10_samples()},
        {"HEX8", consts::ElementType::HEX8, ElementType::Hex8, 3, 2, 3, 8.0, hex8_samples()},
        {"HEX20", consts::ElementType::HEX20, ElementType::Hex20, 3, 4, 5, 8.0, hex27_samples()},
        {"HEX27", consts::ElementType::HEX27, ElementType::Hex27, 3, 4, 5, 8.0, hex27_samples()},
        {"WDG", consts::ElementType::WDG, ElementType::Wedge6, 3, 2, 2, 1.0, wedge6_samples()},
    };
    return cases;
}

std::vector<QuadPoint> rule_points(const ExpectedSamples& samples)
{
    std::vector<QuadPoint> points;
    points.reserve(samples.size());
    for (const auto& sample : samples) {
        points.push_back({
            sample.point[0],
            sample.point[1],
            sample.point[2],
        });
    }
    return points;
}

std::vector<double> rule_weights(const ExpectedSamples& samples)
{
    std::vector<double> weights;
    weights.reserve(samples.size());
    for (const auto& sample : samples) {
        weights.push_back(sample.weight);
    }
    return weights;
}

std::vector<QuadPoint> rule_points(const Array<double>& points, int dimension)
{
    std::vector<QuadPoint> result;
    result.reserve(static_cast<std::size_t>(points.ncols()));
    for (int q = 0; q < points.ncols(); ++q) {
        QuadPoint point = QuadPoint::Zero();
        for (int component = 0; component < dimension; ++component) {
            point[static_cast<std::size_t>(component)] = points(component, q);
        }
        result.push_back(point);
    }
    return result;
}

std::vector<double> rule_weights(const Vector<double>& weights)
{
    std::vector<double> result;
    result.reserve(static_cast<std::size_t>(weights.size()));
    for (int q = 0; q < weights.size(); ++q) {
        result.push_back(weights[q]);
    }
    return result;
}

double factorial(int n)
{
    double value = 1.0;
    for (int factor = 2; factor <= n; ++factor) {
        value *= static_cast<double>(factor);
    }
    return value;
}

double interval_moment(int power)
{
    return power % 2 == 0
               ? 2.0 / static_cast<double>(power + 1)
               : 0.0;
}

double simplex_moment(int px, int py, int pz, int dimension)
{
    const int total = px +
                      (dimension >= 2 ? py : 0) +
                      (dimension >= 3 ? pz : 0);
    double numerator = factorial(px);
    if (dimension >= 2) {
        numerator *= factorial(py);
    }
    if (dimension >= 3) {
        numerator *= factorial(pz);
    }
    return numerator / factorial(total + dimension);
}

double exact_reference_moment(
    svmp::CellFamily family,
    int px,
    int py,
    int pz)
{
    switch (family) {
        case svmp::CellFamily::Point:
            return px == 0 && py == 0 && pz == 0 ? 1.0 : 0.0;
        case svmp::CellFamily::Line:
            return interval_moment(px);
        case svmp::CellFamily::Triangle:
            return simplex_moment(px, py, 0, 2);
        case svmp::CellFamily::Quad:
            return interval_moment(px) * interval_moment(py);
        case svmp::CellFamily::Tetra:
            return simplex_moment(px, py, pz, 3);
        case svmp::CellFamily::Hex:
            return interval_moment(px) *
                   interval_moment(py) *
                   interval_moment(pz);
        case svmp::CellFamily::Wedge:
            return simplex_moment(px, py, 0, 2) * interval_moment(pz);
        default:
            svmp::raise<InvalidArgumentException>(
                "Unsupported reference cell in exactness test");
    }
}

double integer_power(double base, int exponent)
{
    double value = 1.0;
    for (int factor = 0; factor < exponent; ++factor) {
        value *= base;
    }
    return value;
}

double integrate_monomial(
    const QuadratureRule& rule,
    int px,
    int py,
    int pz)
{
    double result = 0.0;
    for (std::size_t q = 0; q < rule.num_points(); ++q) {
        const auto point = rule.point(q);
        result += rule.weight(q) *
                  integer_power(point[0], px) *
                  integer_power(point[1], py) *
                  integer_power(point[2], pz);
    }
    return result;
}

void expect_total_degree_exact(
    const QuadratureRule& rule,
    int degree,
    double tolerance = kExactnessTol)
{
    const int max_px = rule.dimension() >= 1 ? degree : 0;
    for (int px = 0; px <= max_px; ++px) {
        const int max_py = rule.dimension() >= 2 ? degree - px : 0;
        for (int py = 0; py <= max_py; ++py) {
            const int max_pz =
                rule.dimension() >= 3 ? degree - px - py : 0;
            for (int pz = 0; pz <= max_pz; ++pz) {
                EXPECT_NEAR(
                    integrate_monomial(rule, px, py, pz),
                    exact_reference_moment(
                        rule.cell_family(), px, py, pz),
                    tolerance)
                    << "powers=(" << px << ',' << py << ',' << pz << ')';
            }
        }
    }
}

void expect_samples_in_order(
    const QuadratureRule& rule,
    const ExpectedSamples& expected,
    double tolerance = kTol)
{
    ASSERT_EQ(rule.num_points(), expected.size());
    for (std::size_t q = 0; q < expected.size(); ++q) {
        EXPECT_NEAR(rule.weight(q), expected[q].weight, tolerance)
            << "sample=" << q;
        const auto point = rule.point(q);
        for (std::size_t component = 0; component < 3u; ++component) {
            EXPECT_NEAR(
                point[component],
                expected[q].point[component],
                tolerance)
                << "sample=" << q << " component=" << component;
        }
    }
}

RuleProbe make_rule(const SolverRuleCase& c)
{
    return RuleProbe(
        to_mesh_family(c.fe_type),
        c.advertised_exactness,
        rule_points(c.samples),
        rule_weights(c.samples));
}

void expect_legacy_arrays_in_order(const Vector<double>& weights,
                                   const Array<double>& points,
                                   int dimension,
                                   const ExpectedSamples& expected,
                                   double tol = kTol)
{
    ASSERT_EQ(weights.size(), static_cast<int>(expected.size()));
    ASSERT_EQ(points.nrows(), dimension);
    if (dimension > 0) {
        ASSERT_EQ(points.ncols(), static_cast<int>(expected.size()));
    }

    for (std::size_t q = 0; q < expected.size(); ++q) {
        EXPECT_NEAR(weights[static_cast<int>(q)], expected[q].weight, tol)
            << "sample=" << q;
        for (int d = 0; d < dimension; ++d) {
            EXPECT_NEAR(points(d, static_cast<int>(q)),
                        expected[q].point[static_cast<std::size_t>(d)], tol)
                << "sample=" << q << " component=" << d;
        }
    }
}

bool generic_context_supports(consts::ElementType type)
{
    switch (type) {
        case consts::ElementType::LIN1:
        case consts::ElementType::TRI3:
        case consts::ElementType::QUD4:
        case consts::ElementType::QUD9:
        case consts::ElementType::TET4:
        case consts::ElementType::HEX8:
        case consts::ElementType::HEX20:
        case consts::ElementType::HEX27:
        case consts::ElementType::WDG:
            return true;
        default:
            return false;
    }
}

bool volume_context_supports(consts::ElementType type)
{
    switch (type) {
        case consts::ElementType::LIN1:
        case consts::ElementType::LIN2:
        case consts::ElementType::TRI3:
        case consts::ElementType::TRI6:
        case consts::ElementType::QUD4:
        case consts::ElementType::QUD9:
        case consts::ElementType::TET4:
        case consts::ElementType::TET10:
        case consts::ElementType::HEX8:
        case consts::ElementType::HEX20:
        case consts::ElementType::HEX27:
        case consts::ElementType::WDG:
            return true;
        default:
            return false;
    }
}

bool face_context_supports(consts::ElementType type)
{
    switch (type) {
        case consts::ElementType::PNT:
        case consts::ElementType::LIN1:
        case consts::ElementType::LIN2:
        case consts::ElementType::TRI3:
        case consts::ElementType::TRI6:
        case consts::ElementType::QUD4:
        case consts::ElementType::QUD8:
        case consts::ElementType::QUD9:
            return true;
        default:
            return false;
    }
}

} // namespace

TEST(QuadraturePhase01Baseline, StandardSelectionTableHasOrderedPointWeightData)
{
    const auto& cases = standard_solver_cases();
    ASSERT_EQ(cases.size(), 14u);

    for (const auto& c : cases) {
        SCOPED_TRACE(c.name);
        EXPECT_FALSE(c.samples.empty());
        EXPECT_GT(c.requested_exactness, 0);
        EXPECT_GE(c.advertised_exactness, c.requested_exactness);
        EXPECT_GT(c.reference_measure, 0.0);

        double weight_sum = 0.0;
        for (const auto& sample : c.samples) {
            weight_sum += sample.weight;
        }
        EXPECT_NEAR(weight_sum, c.reference_measure, kTol);
    }
}

TEST(QuadraturePhase01Baseline, CanonicalFixturesSatisfyTheRuleAndExactnessContracts)
{
    for (const auto& c : standard_solver_cases()) {
        SCOPED_TRACE(c.name);
        const auto rule = make_rule(c);

        EXPECT_EQ(rule.cell_family(), to_mesh_family(c.fe_type));
        EXPECT_EQ(rule.dimension(), c.dimension);
        EXPECT_EQ(rule.polynomial_exactness(), c.advertised_exactness);
        EXPECT_EQ(rule.num_points(), c.samples.size());
        EXPECT_DOUBLE_EQ(rule.reference_measure(), c.reference_measure);
        EXPECT_TRUE(rule.is_structurally_valid());
        expect_samples_in_order(rule, c.samples);
        expect_total_degree_exact(rule, c.advertised_exactness);
    }
}

TEST(QuadraturePhase01Baseline, GenericFunctionSpaceRulesMatchOrderedLegacyData)
{
    for (const auto& c : standard_solver_cases()) {
        if (!generic_context_supports(c.solver_type) ||
            c.solver_type == consts::ElementType::QUD9) {
            // QUD9 is characterized separately because the legacy generic
            // table contains an out-of-domain sqrt(6) coordinate defect.
            continue;
        }

        SCOPED_TRACE(c.name);
        Vector<double> weights(static_cast<int>(c.samples.size()));
        Array<double> points(c.dimension, static_cast<int>(c.samples.size()));
        nn::get_gip(c.dimension, c.solver_type,
                    static_cast<int>(c.samples.size()), weights, points);
        expect_legacy_arrays_in_order(weights, points, c.dimension, c.samples);
    }
}

TEST(QuadraturePhase01Baseline, GenericQud9SqrtSixCoordinatesAreAnExplicitLegacyDefect)
{
    const auto legacy_defect = quad9_samples(std::sqrt(6.0));
    Vector<double> weights(static_cast<int>(legacy_defect.size()));
    Array<double> points(2, static_cast<int>(legacy_defect.size()));

    nn::get_gip(
        2,
        consts::ElementType::QUD9,
        static_cast<int>(legacy_defect.size()),
        weights,
        points);

    expect_legacy_arrays_in_order(weights, points, 2, legacy_defect);
    EXPECT_GT(std::abs(points(0, 0)), 1.0);

    double x_squared_moment = 0.0;
    for (int q = 0; q < weights.size(); ++q) {
        x_squared_moment += weights[q] * points(0, q) * points(0, q);
    }
    EXPECT_NEAR(x_squared_moment, 40.0 / 3.0, kExactnessTol);
    EXPECT_NEAR(
        exact_reference_moment(svmp::CellFamily::Quad, 2, 0, 0),
        4.0 / 3.0,
        kTol);

    EXPECT_THROW(
        (void)RuleProbe(
            svmp::CellFamily::Quad,
            5,
            rule_points(points, 2),
            rule_weights(weights)),
        InvalidArgumentException);
}

TEST(QuadraturePhase01Baseline, VolumeMeshRulesMatchOrderedLegacyData)
{
    for (const auto& c : standard_solver_cases()) {
        if (!volume_context_supports(c.solver_type)) {
            continue;
        }

        SCOPED_TRACE(c.name);
        mshType mesh;
        mesh.eType = c.solver_type;
        mesh.nG = static_cast<int>(c.samples.size());
        mesh.w = Vector<double>(mesh.nG);
        mesh.xi = Array<double>(c.dimension, mesh.nG);
        nn::get_gip(mesh);
        expect_legacy_arrays_in_order(mesh.w, mesh.xi, c.dimension, c.samples);
    }
}

TEST(QuadraturePhase01Baseline, Qud8VolumeSelectionHasNoLegacyPopulationEntry)
{
    const auto& qud8_case = standard_solver_cases()[6];
    ASSERT_EQ(qud8_case.solver_type, consts::ElementType::QUD8);

    mshType mesh;
    mesh.eType = qud8_case.solver_type;
    mesh.nG = static_cast<int>(qud8_case.samples.size());
    mesh.w = Vector<double>(mesh.nG);
    mesh.xi = Array<double>(qud8_case.dimension, mesh.nG);

    EXPECT_THROW(nn::get_gip(mesh), InvalidElementException);

    const auto canonical_rule = make_rule(qud8_case);
    expect_samples_in_order(canonical_rule, qud8_case.samples);
    expect_total_degree_exact(
        canonical_rule,
        qud8_case.advertised_exactness);
}

TEST(QuadraturePhase01Baseline, GenericFunctionSpaceMissingEntriesAreExplicit)
{
    for (const auto& c : standard_solver_cases()) {
        if (generic_context_supports(c.solver_type)) {
            continue;
        }

        SCOPED_TRACE(c.name);
        Vector<double> weights(static_cast<int>(c.samples.size()));
        Array<double> points(c.dimension, static_cast<int>(c.samples.size()));
        EXPECT_THROW(
            nn::get_gip(c.dimension, c.solver_type,
                        static_cast<int>(c.samples.size()), weights, points),
            InvalidElementException);
    }
}

TEST(QuadraturePhase01Baseline, FaceRulesMatchOrderedLegacyData)
{
    for (const auto& c : standard_solver_cases()) {
        if (!face_context_supports(c.solver_type)) {
            continue;
        }

        SCOPED_TRACE(c.name);
        faceType face;
        face.eType = c.solver_type;
        face.nG = static_cast<int>(c.samples.size());
        face.w = Vector<double>(face.nG);
        face.xi = Array<double>(c.dimension, face.nG);
        nn::get_gip(nullptr, face);
        expect_legacy_arrays_in_order(face.w, face.xi, c.dimension, c.samples);
    }
}

TEST(QuadraturePhase01Baseline, LegacyPositionModifiersAreOnlyDegreeOneExact)
{
    mshType tetra_mesh;
    tetra_mesh.eType = consts::ElementType::TET4;
    tetra_mesh.nG = 4;
    tetra_mesh.qmTET4 = 0.25;
    tetra_mesh.w = Vector<double>(tetra_mesh.nG);
    tetra_mesh.xi = Array<double>(3, tetra_mesh.nG);
    nn::get_gip(tetra_mesh);

    const RuleProbe tetra_rule(
        svmp::CellFamily::Tetra,
        1,
        rule_points(tetra_mesh.xi, 3),
        rule_weights(tetra_mesh.w));
    expect_total_degree_exact(tetra_rule, 1);
    EXPECT_NEAR(
        integrate_monomial(tetra_rule, 2, 0, 0),
        1.0 / 96.0,
        kExactnessTol);
    EXPECT_GT(
        std::abs(
            integrate_monomial(tetra_rule, 2, 0, 0) -
            exact_reference_moment(svmp::CellFamily::Tetra, 2, 0, 0)),
        1.0e-3);

    faceType triangle_face;
    triangle_face.eType = consts::ElementType::TRI3;
    triangle_face.nG = 3;
    triangle_face.qmTRI3 = 1.0 / 3.0;
    triangle_face.w = Vector<double>(triangle_face.nG);
    triangle_face.xi = Array<double>(2, triangle_face.nG);
    nn::get_gip(nullptr, triangle_face);

    const RuleProbe triangle_rule(
        svmp::CellFamily::Triangle,
        1,
        rule_points(triangle_face.xi, 2),
        rule_weights(triangle_face.w));
    expect_total_degree_exact(triangle_rule, 1);
    EXPECT_NEAR(
        integrate_monomial(triangle_rule, 2, 0, 0),
        1.0 / 18.0,
        kExactnessTol);
    EXPECT_GT(
        std::abs(
            integrate_monomial(triangle_rule, 2, 0, 0) -
            exact_reference_moment(svmp::CellFamily::Triangle, 2, 0, 0)),
        1.0e-3);
}

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
        EXPECT_DOUBLE_EQ(rule.reference_measure(), c.expected_measure);
        EXPECT_TRUE(rule.is_structurally_valid());
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
    expect_invalid_argument_with_message(
        [] {
            (void)RuleProbe(
                svmp::CellFamily::Line,
                1,
                {{0.0, 0.0, 0.0}, {0.5, 0.0, 0.0}},
                {std::numeric_limits<double>::max(),
                 std::numeric_limits<double>::max()});
        },
        "weight accumulation is not finite at sample 1");
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

TEST(QuadratureRuleValidation, EnforcesMeasureButAllowsNegativeWeights)
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
    EXPECT_TRUE(rule.is_structurally_valid());
    EXPECT_TRUE(rule.is_structurally_valid(0.0));
    EXPECT_LT(rule.weight(0), 0.0);
    EXPECT_DOUBLE_EQ(
        integrate_monomial(rule, 0, 0, 0),
        rule.reference_measure());
    EXPECT_FALSE(rule.is_structurally_valid(-1.0));
    EXPECT_FALSE(rule.is_structurally_valid(
        std::numeric_limits<double>::quiet_NaN()));
    EXPECT_FALSE(rule.is_structurally_valid(
        std::numeric_limits<double>::infinity()));
}

TEST(QuadratureRuleValidation, RejectsStoredOrderCancellation)
{
    expect_invalid_argument_with_message(
        [] {
            (void)RuleProbe(
                svmp::CellFamily::Line,
                0,
                {{-0.75, 0.0, 0.0},
                 {-0.25, 0.0, 0.0},
                 {0.25, 0.0, 0.0},
                 {0.75, 0.0, 0.0}},
                {1.0e16, 1.0, -1.0e16, 1.0});
        },
        "stored-order double-precision weight accumulation");
}

TEST(QuadratureRuleValidation, RejectsIllConditionedCancellationDespiteExactOrderedSum)
{
    expect_invalid_argument_with_message(
        [] {
            (void)RuleProbe(
                svmp::CellFamily::Line,
                0,
                {{-0.75, 0.0, 0.0},
                 {-0.25, 0.0, 0.0},
                 {0.25, 0.0, 0.0},
                 {0.75, 0.0, 0.0}},
                {1.0e16, -1.0e16, 1.0, 1.0});
        },
        "too ill-conditioned for stable double-precision integration");
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

TEST(QuadratureRuleValidation, AcceptsWellConditionedLargePositiveRule)
{
    constexpr std::size_t sample_count = 8192u;
    std::vector<QuadPoint> points(sample_count, QuadPoint::Zero());
    std::vector<double> weights(
        sample_count,
        2.0 / static_cast<double>(sample_count));

    const RuleProbe rule(
        svmp::CellFamily::Line,
        0,
        std::move(points),
        std::move(weights));
    EXPECT_TRUE(rule.is_structurally_valid());
    EXPECT_DOUBLE_EQ(
        integrate_monomial(rule, 0, 0, 0),
        rule.reference_measure());
}

TEST(QuadratureRuleValidation, AcceptsWellConditionedLargeSignedRule)
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
    EXPECT_TRUE(rule.is_structurally_valid());
    EXPECT_LT(rule.weight(0), 0.0);
    EXPECT_NEAR(
        integrate_monomial(rule, 0, 0, 0),
        rule.reference_measure(),
        QuadratureRule::default_validation_tolerance());
}

TEST(QuadratureRuleValidation, UsesDocumentedCoordinateTolerance)
{
    constexpr double offset = 0.5 * QuadratureRule::default_validation_tolerance();
    const RuleProbe rule(
        svmp::CellFamily::Line, 0, {{1.0 + offset, offset, 0.0}}, {2.0});
    EXPECT_TRUE(rule.is_structurally_valid());
    EXPECT_FALSE(rule.is_structurally_valid(offset / 2.0));
}

TEST(QuadratureRuleValidation, UsesDocumentedWeightTolerance)
{
    constexpr double offset =
        0.5 * QuadratureRule::default_validation_tolerance();
    const RuleProbe rule(
        svmp::CellFamily::Triangle,
        0,
        {{0.25, 0.25, 0.0}},
        {0.5 + offset});
    EXPECT_TRUE(rule.is_structurally_valid());
    EXPECT_FALSE(rule.is_structurally_valid(offset / 2.0));
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
    EXPECT_DOUBLE_EQ(rule.reference_measure(), 2.0);
    EXPECT_TRUE(rule.is_structurally_valid());
    ASSERT_EQ(rule.num_points(), 2u);
    ASSERT_EQ(rule.points().size(), 2u);
    ASSERT_EQ(rule.weights().size(), 2u);
    EXPECT_DOUBLE_EQ(rule.point(0)[0], -a);
    EXPECT_DOUBLE_EQ(rule.point(1)[0], a);
    EXPECT_DOUBLE_EQ(rule.weight(0), 1.0);
    EXPECT_DOUBLE_EQ(rule.weight(1), 1.0);
}
