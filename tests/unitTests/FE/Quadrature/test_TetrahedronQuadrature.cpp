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
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
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

const double gauss_node = 1.0 / std::sqrt(3.0);

QuadratureRule make_cancelling_factor(double magnitude, double node = -gauss_node)
{
    return QuadratureRule(svmp::CellFamily::Line, 3,
        {{node, 0.0, 0.0}, {node, 0.0, 0.0}, {-gauss_node, 0.0, 0.0}, {gauss_node, 0.0, 0.0}},
        {magnitude, -magnitude, 1.0, 1.0});
}

template<class Function>
void expect_invalid(Function&& function, std::string_view diagnostic)
{
    try {
        function();
        FAIL() << "Expected InvalidArgumentException containing " << diagnostic;
    } catch (const InvalidArgumentException& error) {
        EXPECT_EQ(error.status_code(), svmp::StatusCode::InvalidArgument);
        EXPECT_GT(error.context().line(), 0);
        EXPECT_NE(error.message().find(diagnostic), std::string_view::npos)
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

constexpr long double kMomentRelativeTolerance = 2.0e-12L;
using MomentTable = std::vector<std::vector<std::vector<CompensatedSum>>>;

// Regroup only equal actual x and (x,y) coordinates, including disjoint groups
// when samples are reordered. No factor data or reconstructed Duffy map enters
// the sum. Store just the moments with a+b+c <= degree.
MomentTable tetrahedron_moments(const QuadratureRule& rule, int degree)
{
    const std::size_t size = static_cast<std::size_t>(degree + 1);
    MomentTable moments(size);
    std::vector<std::vector<CompensatedSum>> plane(size);
    std::vector<CompensatedSum> row(size);
    for (std::size_t a = 0; a < size; ++a) {
        moments[a].resize(size - a);
        for (std::size_t b = 0; b < size - a; ++b) {
            moments[a][b].resize(size - a - b);
        }
        plane[a].resize(size - a);
    }
    for (std::size_t first = 0; first < rule.num_points();) {
        for (auto& row_of_plane : plane) {
            std::fill(row_of_plane.begin(), row_of_plane.end(), CompensatedSum{});
        }
        const long double x = rule.point(first)[0];
        do {
            std::fill(row.begin(), row.end(), CompensatedSum{});
            const long double y = rule.point(first)[1];
            do {
                const long double z = rule.point(first)[2];
                long double power = 1.0L;
                for (int c = 0; c <= degree; ++c) {
                    row[c].add(static_cast<long double>(rule.weight(first)) * power);
                    power *= z;
                }
                ++first;
            } while (first < rule.num_points() && rule.point(first)[0] == x && rule.point(first)[1] == y);
            long double power = 1.0L;
            for (int b = 0; b <= degree; ++b) {
                for (int c = 0; c <= degree - b; ++c) {
                    plane[b][c].add(power * row[c].value());
                }
                power *= y;
            }
        } while (first < rule.num_points() && rule.point(first)[0] == x);
        long double power = 1.0L;
        for (int a = 0; a <= degree; ++a) {
            for (int b = 0; b <= degree - a; ++b) {
                for (int c = 0; c <= degree - a - b; ++c) {
                    moments[a][b][c].add(power * plane[b][c].value());
                }
            }
            power *= x;
        }
    }
    return moments;
}

long double expect_all_moments(const QuadratureRule& rule, bool direct = false)
{
    const int degree = rule.polynomial_exactness();
    const auto moments = direct ? MomentTable{} : tetrahedron_moments(rule, degree);
    long double maximum_error = 0.0L;
    std::array<int, 3> worst{};
    std::size_t invalid = 0;
    for (int a = 0; a <= degree; ++a) {
        for (int b = 0; b <= degree - a; ++b) {
            for (int c = 0; c <= degree - a - b; ++c) {
                const long double expected = reference::tetrahedron_monomial_integral(a, b, c);
                const long double actual = direct ? direct_moment(rule, a, b, c) : moments[a][b][c].value();
                const long double error = std::isfinite(actual)
                    ? std::abs(actual - expected) / expected : std::numeric_limits<long double>::infinity();
                if (error > maximum_error) {
                    maximum_error = error;
                    worst = {a, b, c};
                }
                if (!(actual > 0.0L) || !(error <= kMomentRelativeTolerance)) {
                    ++invalid;
                }
            }
        }
    }
    EXPECT_EQ(invalid, 0u) << "degree " << degree << ", worst powers "
        << worst[0] << ',' << worst[1] << ',' << worst[2]
        << ", relative error " << maximum_error;
    return maximum_error;
}

void record_moment_error(long double maximum_error)
{
    std::ostringstream value;
    value << std::setprecision(21) << maximum_error;
    ::testing::Test::RecordProperty("max_relative_error", value.str());
    ::testing::Test::RecordProperty("long_double_digits", std::numeric_limits<long double>::digits);
    ::testing::Test::RecordProperty("long_double_max_exponent", std::numeric_limits<long double>::max_exponent);
    std::cout << "Maximum relative moment error: " << value.str() << std::endl;
}

QuadratureRule make_line(bool lobatto, int degree)
{
    return lobatto ? make_gauss_lobatto_rule(degree) : make_gauss_legendre_rule(degree);
}

int max_line_exactness(bool lobatto)
{
    return lobatto ? max_gauss_lobatto_exactness() : max_gauss_legendre_exactness();
}

void expect_structure(const QuadratureRule& rule, const QuadratureRule& u_rule,
    const QuadratureRule& v_rule, const QuadratureRule& w_rule, bool closed)
{
    const std::size_t count = rule.num_points();
    ASSERT_EQ(count, u_rule.num_points() * v_rule.num_points() * w_rule.num_points());
    EXPECT_EQ(rule.cell_family(), svmp::CellFamily::Tetra);
    EXPECT_EQ(rule.dimension(), 3u);
    EXPECT_EQ(rule.reference_cell_measure(), 1.0 / 6.0);
    EXPECT_EQ(rule.polynomial_exactness(), std::min({u_rule.polynomial_exactness() - 2,
        v_rule.polynomial_exactness() - 1, w_rule.polynomial_exactness()}));
    std::size_t first_nonfinite = count, first_exterior = count, first_bad_sign = count;
    std::size_t first_bad_order = count, first_bad_collapse = count;
    CompensatedSum measure;
    constexpr double tolerance = 8.0 * std::numeric_limits<double>::epsilon();
    for (std::size_t q = 0; q < count; ++q) {
        const double x = rule.point(q)[0], y = rule.point(q)[1], z = rule.point(q)[2];
        const double weight = rule.weight(q);
        if (!(std::isfinite(x) && std::isfinite(y) && std::isfinite(z) && std::isfinite(weight))) {
            first_nonfinite = std::min(first_nonfinite, q);
        }
        const bool inside = closed ? x >= 0.0 && y >= 0.0 && z >= 0.0 && x + y + z <= 1.0 + tolerance
                                   : x > 0.0 && y > 0.0 && z > 0.0 && x + y + z < 1.0;
        if (!inside) first_exterior = std::min(first_exterior, q);
        if (!(closed ? weight >= 0.0 : weight > 0.0)) first_bad_sign = std::min(first_bad_sign, q);
        const std::size_t i = q / (v_rule.num_points() * w_rule.num_points());
        const std::size_t j = (q / w_rule.num_points()) % v_rule.num_points();
        const std::size_t k = q % w_rule.num_points();
        const long double xi = u_rule.point(i)[0];
        const long double eta = v_rule.point(j)[0];
        const long double zeta = w_rule.point(k)[0];
        const long double expected_x = (1 + xi) / 2;
        const long double expected_y = (1 - xi) * (1 + eta) / 4;
        const long double expected_z = (1 - xi) * (1 - eta) * (1 + zeta) / 8;
        if (std::abs(x - expected_x) > tolerance * expected_x ||
            std::abs(y - expected_y) > tolerance * expected_y ||
            std::abs(z - expected_z) > tolerance * expected_z) {
            first_bad_order = std::min(first_bad_order, q);
        }
        if ((xi == 1 || eta == 1) && weight != 0.0) {
            first_bad_collapse = std::min(first_bad_collapse, q);
        }
        measure.add(weight);
    }
    EXPECT_EQ(first_nonfinite, count) << "first nonfinite sample";
    EXPECT_EQ(first_exterior, count) << "first exterior sample";
    EXPECT_EQ(first_bad_sign, count) << "first wrong weight sign";
    EXPECT_EQ(first_bad_order, count) << "first wrong ordered coordinate";
    EXPECT_EQ(first_bad_collapse, count) << "first nonzero collapsed weight";
    EXPECT_NEAR(measure.value(), 1.0L / 6.0L, kMomentRelativeTolerance / 6.0L);
}

// Evaluate P_(p-1)(2t-1) and P_p(2t-1) together for p >= 1.
// Advancing a block of independent samples allows the inner loop to vectorize.
template <std::size_t N>
auto shifted_legendre_pair(std::array<double, N> coordinates, int p)
{
    for (double& t : coordinates) {
        t = 2.0 * t - 1.0;
    }
    std::array<double, N> previous;
    previous.fill(1.0);
    auto current = coordinates;
    for (int k = 2; k <= p; ++k) {
        // P_k = s P_(k-1) + (k-1)/k [s P_(k-1) - P_(k-2)].
        const double ratio = static_cast<double>(k - 1) / k;
        for (std::size_t i = 0; i < N; ++i) {
            const double product = coordinates[i] * current[i];
            const double next = product + ratio * (product - previous[i]);
            previous[i] = current[i];
            current[i] = next;
        }
    }
    return std::pair{previous, current};
}

void expect_high_degree_legendre_moments(const QuadratureRule& rule, int p)
{
    // At odd p, symmetry can hide a short u factor in the top z moment.
    // The preceding z moment exposes that error. The x moment checks the outer factor independently.
    std::array<CompensatedSum, 3> sums{}, magnitudes{};
    std::optional<double> previous_x;
    double x_high = 0.0;
    constexpr std::size_t block_size = 128;
    for (std::size_t first = 0; first < rule.num_points(); first += block_size) {
        const std::size_t count = std::min(block_size, rule.num_points() - first);
        std::array<double, block_size> z{};
        for (std::size_t i = 0; i < count; ++i) {
            z[i] = rule.point(first + i)[2];
        }
        const auto [z_previous, z_high] = shifted_legendre_pair(z, p);
        for (std::size_t i = 0; i < count; ++i) {
            const double x = rule.point(first + i)[0];
            if (x != previous_x) {
                x_high = shifted_legendre_pair(std::array{x}, p).second[0];
                previous_x = x;
            }
            const std::array values{x_high, z_previous[i], z_high[i]};
            for (std::size_t j = 0; j < values.size(); ++j) {
                const long double term = static_cast<long double>(rule.weight(first + i)) * values[j];
                sums[j].add(term);
                magnitudes[j].add(std::abs(term));
            }
        }
    }

    // At fixed x or z = t the cross-section has area (1-t)^2/2, orthogonal
    // to every shifted Legendre polynomial of degree greater than two.
    // Allow for the O(p^2) endpoint sensitivity of the polynomial on double
    // coordinates, scaled by the absolute weighted terms in each zero moment.
    const long double roundoff = p * (p + 1) * std::numeric_limits<double>::epsilon();
    constexpr std::array labels{"P_p(2x-1)", "P_(p-1)(2z-1)", "P_p(2z-1)"};
    for (std::size_t i = 0; i < sums.size(); ++i) {
        EXPECT_LE(std::abs(sums[i].value()), roundoff * magnitudes[i].value()) << labels[i];
    }
}

void expect_representative_rules(bool lobatto)
{
    const int last_request = lobatto ? 251 : 253;
    long double maximum_error = 0.0L;
    for (const int p : {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 31, 63,
                       last_request - 1, last_request}) {
        SCOPED_TRACE(p);
        const auto line = make_line(lobatto, p + 2);
        const auto rule = make_tetrahedron_rule(p, line);
        const std::size_t n = static_cast<std::size_t>(p / 2 + (lobatto ? 3 : 2));
        EXPECT_EQ(rule.num_points(), n * n * n);
        EXPECT_EQ(rule.polynomial_exactness(), 2 * static_cast<int>(n) - (lobatto ? 5 : 3));
        expect_structure(rule, line, line, line, lobatto);
        expect_same_rule(rule, make_tetrahedron_rule(p, line, line, line));
        if (p <= 63) {
            maximum_error = std::max(maximum_error, expect_all_moments(rule, p <= 12));
        } else if (p == last_request) {
            expect_high_degree_legendre_moments(rule, p);
        }
    }
    record_moment_error(maximum_error);
}

// A set bit selects Lobatto: 4 = u, 2 = v, 1 = w.
int maximum_mixed_request(int combination)
{
    const int e_u = max_line_exactness(combination & 4);
    const int e_v = max_line_exactness(combination & 2);
    const int e_w = max_line_exactness(combination & 1);
    return std::min({max_tetrahedron_exactness(), e_u - 2, e_v - 1, e_w});
}

void expect_mixed_rules(bool exhaustive)
{
    long double maximum_error = 0.0L;
    for (int combination = 0; combination < 8; ++combination) {
        const std::vector requests = exhaustive ? std::vector{127, maximum_mixed_request(combination)}
                                                : std::vector{0, 1, 2, 5, 31};
        for (const int p : requests) {
            SCOPED_TRACE(::testing::Message() << "combination " << combination << ", request " << p);
            if (exhaustive) std::cout << "Mixed combination " << combination << ", request " << p << std::endl;
            const auto u = make_line(combination & 4, p + 2);
            const auto v = make_line(combination & 2, p + 1);
            const auto w = make_line(combination & 1, p);
            const auto rule = make_tetrahedron_rule(p, u, v, w);
            EXPECT_EQ(rule.polynomial_exactness(), p);
            expect_structure(rule, u, v, w, combination != 0);
            maximum_error = std::max(maximum_error, expect_all_moments(rule, p <= 5));
        }
    }
    record_moment_error(maximum_error);
}

void expect_every_rule(bool lobatto)
{
    const int last_request = lobatto ? 251 : 253;
    int distinct = 0;
    long double maximum_error = 0.0L;
    std::optional<QuadratureRule> previous;
    for (int p = 0; p <= last_request; ++p) {
        SCOPED_TRACE(p);
        const auto line = make_line(lobatto, p + 2);
        auto rule = make_tetrahedron_rule(p, line);
        const std::size_t n = static_cast<std::size_t>(p / 2 + (lobatto ? 3 : 2));
        const int degree = 2 * static_cast<int>(n) - (lobatto ? 5 : 3);
        ASSERT_EQ(rule.num_points(), n * n * n);
        ASSERT_EQ(rule.polynomial_exactness(), degree);
        expect_structure(rule, line, line, line, lobatto);
        expect_same_rule(rule, make_tetrahedron_rule(p, line, line, line));
        if (previous && previous->polynomial_exactness() == degree) {
            expect_same_rule(*previous, rule);
        } else {
            ++distinct;
            maximum_error = std::max(maximum_error, expect_all_moments(rule));
        }
        if (p % 32 == 0 || p == last_request) {
            std::cout << (lobatto ? "Lobatto" : "Gauss") << " qualified through request " << p << std::endl;
        }
        previous = std::move(rule);
    }
    EXPECT_EQ(distinct, lobatto ? 126 : 127);
    ASSERT_TRUE(previous.has_value());
    EXPECT_EQ(previous->num_points(), 2097152u);
    record_moment_error(maximum_error);
}

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
        expect_invalid([&] { reference::tetrahedron_monomial_integral(negative, 0, 0); },
            "Tetrahedron monomial exponents");
        EXPECT_THROW((void)reference::tetrahedron_monomial_integral(0, negative, 0), InvalidArgumentException);
        EXPECT_THROW((void)reference::tetrahedron_monomial_integral(0, 0, negative), InvalidArgumentException);
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
    const QuadratureRule reversed_middle(svmp::CellFamily::Tetra, 0,
        {{0.0, 1.0, 0.0}, {0.0, 0.0, 0.5}, {0.5, 0.5, 0.0},
         {0.5, 0.0, 0.25}, {1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}},
        {0.0, 1.0 / 12.0, 0.0, 1.0 / 12.0, 0.0, 0.0});
    expect_same_rule(make_tetrahedron_rule(0, u, reversed(make_gauss_lobatto_rule(1)), v), reversed_middle);
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
    constexpr std::size_t capacity = 128u * 128u * 128u;
    const auto pair = make_gauss_legendre_rule(3);
    const auto split_gauss_pair = [&pair](std::size_t count) {
        // Split each Gauss node's weight evenly among its copies.
        const std::size_t first_count = count / 2;
        std::vector<QuadPoint> points(count, pair.point(1));
        std::fill_n(points.begin(), first_count, pair.point(0));
        std::vector<double> weights(count, pair.weight(1) / (count - first_count));
        std::fill_n(weights.begin(), first_count, pair.weight(0) / first_count);
        return QuadratureRule(pair.cell_family(), pair.polynomial_exactness(), std::move(points), std::move(weights));
    };

    const std::array at_limit{split_gauss_pair(256), split_gauss_pair(128), split_gauss_pair(64)};
    for (std::size_t i = 0; i < at_limit.size(); ++i) {
        SCOPED_TRACE(i);
        const auto rule = make_tetrahedron_rule(0,
            at_limit[i], at_limit[(i + 1) % 3], at_limit[(i + 2) % 3]);
        EXPECT_EQ(rule.num_points(), capacity);
        EXPECT_NEAR(direct_moment(rule, 0, 0, 0), 1.0L / 6.0L, 2.0e-16L);
    }

    // The 9 x 43 x 5419 product exceeds the limit by one. Rotate the
    // unnormalized factor through all roles to check the order of validation.
    const std::array above_limit{split_gauss_pair(9), split_gauss_pair(43),
        QuadratureRule(svmp::CellFamily::Line, 0,
            std::vector<QuadPoint>(5419, QuadPoint::Zero()), std::vector<double>(5419, 0.0))};
    for (std::size_t i = 0; i < above_limit.size(); ++i) {
        SCOPED_TRACE(i);
        expect_invalid([&] {
            (void)make_tetrahedron_rule(0,
                above_limit[i], above_limit[(i + 1) % 3], above_limit[(i + 2) % 3]);
        }, "2097152");
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

TEST(TetrahedronDuffyValidation, PreservesCancellingFiniteWeights)
{
    const auto rule = make_tetrahedron_rule(1, make_cancelling_factor(1.0));
    ASSERT_EQ(rule.num_points(), 64u);
    for (double weight : rule.weights()) {
        ASSERT_TRUE(std::isfinite(weight));
    }
    EXPECT_EQ(rule.weight(0), -rule.weight(1));
    EXPECT_EQ(rule.weight(0), -rule.weight(4));
    EXPECT_EQ(rule.weight(0), -rule.weight(16));
    EXPECT_NEAR(direct_moment(rule, 0, 0, 0), 1.0L / 6.0L, 2.0e-16L);
}

TEST(TetrahedronDuffyValidation, PreservesInteriorSamplesNextToCollapsedFaces)
{
    const double near = std::nextafter(1.0, 0.0);
    const auto midpoint = make_gauss_legendre_rule(1);
    // The distance to the collapsed face in unit coordinates is exactly 2^-54.
    const auto u = make_tetrahedron_rule(0, make_cancelling_factor(1.0, near), midpoint, midpoint);
    EXPECT_EQ(u.point(0)[1], 0x1p-55);
    EXPECT_EQ(u.point(0)[2], 0x1p-56);
    EXPECT_EQ(u.weight(0), 0x1p-110);

    const QuadratureRule near_pair(svmp::CellFamily::Line, 1,
        {{-near, 0.0, 0.0}, {near, 0.0, 0.0}}, {1.0, 1.0});
    const auto v = make_tetrahedron_rule(0, make_cancelling_factor(1.0, -1.0), near_pair, midpoint);
    EXPECT_EQ(v.point(1)[2], 0x1p-55);
    EXPECT_EQ(v.weight(1), 0x1p-56);
}

TEST(TetrahedronDuffyValidation, IntegratesSignedFactorNearEndpoint)
{
    const double near = std::nextafter(1.0, 0.0);
    // The equally spaced final three nodes contribute zero to both line moments.
    const QuadratureRule factor(svmp::CellFamily::Line, 1,
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {near, 0.0, 0.0}, {std::nextafter(near, 0.0), 0.0, 0.0}},
        {2.0, 1024.0, -2048.0, 1024.0});
    const auto rule = make_tetrahedron_rule(0, make_gauss_legendre_rule(3), factor, make_gauss_legendre_rule(1));
    EXPECT_EQ(rule.weight(1), 0.0);
    EXPECT_LT(rule.weight(2), 0.0);
    EXPECT_NEAR(direct_moment(rule, 0, 0, 0), 1.0L / 6.0L, 2.0e-16L);
}

TEST(TetrahedronDuffyValidation, RecoversWeightsAfterPartialUnderflow)
{
    // At three zero nodes the Duffy weight is a*b*c/64. Every expected value is dyadic.
    const auto tiny = make_cancelling_factor(0x1p-540, 0.0);
    const auto large = make_cancelling_factor(0x1p+1020, 0.0);
    const auto rule = make_tetrahedron_rule(1, tiny, tiny, large);
    EXPECT_EQ(rule.weight(0), 0x1p-66);
    EXPECT_EQ(rule.weight(1), -0x1p-66);
    EXPECT_EQ(make_tetrahedron_rule(1, tiny, tiny, make_cancelling_factor(0x1p+12, 0.0)).weight(0),
        std::numeric_limits<double>::denorm_min());
    const auto underflow = make_tetrahedron_rule(1, tiny, tiny, make_cancelling_factor(0x1p+10, 0.0));
    EXPECT_EQ(underflow.weight(0), 0.0);
    EXPECT_EQ(underflow.weight(1), 0.0);
    EXPECT_TRUE(std::signbit(underflow.weight(1)));

    // Each of the three intermediate products can lose a bit before a large multiplier.
    const auto unit = make_cancelling_factor(1.0, 0.0);
    const auto small_u = make_cancelling_factor(3 * 0x1p-1022, 1.0 - 0x1p-24);
    EXPECT_EQ(make_tetrahedron_rule(1, small_u, large, unit).weight(0), 3 * 0x1p-56);
    const auto small_v = make_cancelling_factor(3 * 0x1p-1022, 1.0 - 0x1p-52);
    EXPECT_EQ(make_tetrahedron_rule(1, make_cancelling_factor(0x1p+1018, 0.0), small_v, unit).weight(0),
        3 * 0x1p-62);
    EXPECT_EQ(make_tetrahedron_rule(1, make_cancelling_factor(3 * 0x1p-534, 0.0),
        make_cancelling_factor(0x1p-535, 0.0), large).weight(0), 3 * 0x1p-55);
}

TEST(TetrahedronDuffyValidation, RejectsIntermediateOverflowEvenWhenFinalWeightsFit)
{
    const QuadratureRule w(svmp::CellFamily::Line, 1,
        std::vector<QuadPoint>(8, QuadPoint::Zero()), std::vector<double>(8, 0.25));
    // The partial product is 2^1025 at the first pair of nodes; each final weight would fit.
    EXPECT_THROW((void)make_tetrahedron_rule(0, make_cancelling_factor(0x1p+600, 0.0),
        make_cancelling_factor(0x1p+431, 0.0), w), ConvergenceException);
}

TEST(TetrahedronDuffyValidation, HandlesFiniteWeightsWhoseAbsoluteSumOverflows)
{
    const auto rule = make_tetrahedron_rule(1, make_cancelling_factor(1.0e103));
    ASSERT_EQ(rule.num_points(), 64u);
    // Scale before summing: long double may have binary64 exponent range.
    constexpr long double sum_scale = 1.0L / 64.0L;
    CompensatedSum absolute;
    for (double weight : rule.weights()) {
        ASSERT_TRUE(std::isfinite(weight));
        absolute.add(std::abs(static_cast<long double>(weight)) * sum_scale);
    }
    EXPECT_EQ(rule.weight(0), -rule.weight(1));
    EXPECT_EQ(rule.weight(0), -rule.weight(4));
    EXPECT_EQ(rule.weight(0), -rule.weight(16));
    EXPECT_GT(absolute.value(), static_cast<long double>(std::numeric_limits<double>::max()) * sum_scale);
}

TEST(TetrahedronDuffyValidation, RejectsNonfiniteProductWeights)
{
    try {
        (void)make_tetrahedron_rule(1, make_cancelling_factor(1.0e200));
        FAIL() << "Expected ConvergenceException for nonfinite weights";
    } catch (const ConvergenceException& error) {
        EXPECT_EQ(error.status_code(), svmp::StatusCode::InvalidState);
        EXPECT_GT(error.context().line(), 0);
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

TEST(TetrahedronMomentAccumulation, MatchesDirectSumsOnIrregularSignedAndReorderedSamples)
{
    const QuadratureRule irregular(svmp::CellFamily::Tetra, 7,
        {{0.25, 0.125, 0.25}, {0.25, 0.125, 0.25}, {0.25, 0.125, 0.5},
         {0.25, 0.5, 0.0625}, {0.5, 0.25, 0.125}, {0.5, 0.125, 0.25},
         {0.25, 0.125, 0.25}, {0.0, 0.0, 0.0}},
        {2.0, -1.0, 0.75, -0.25, 0.5, 0.0, -0.125, 0.125});
    for (const auto& rule : {irregular, reversed(irregular)}) {
        const auto moments = tetrahedron_moments(rule, 7);
        ASSERT_EQ(moments.size(), 8u);
        for (int a = 0; a <= 7; ++a) {
            ASSERT_EQ(moments[a].size(), 8u - a);
            for (int b = 0; b <= 7 - a; ++b) {
                ASSERT_EQ(moments[a][b].size(), 8u - a - b);
                for (int c = 0; c <= 7 - a - b; ++c) {
                    const long double direct = direct_moment(rule, a, b, c);
                    EXPECT_NEAR(moments[a][b][c].value(), direct,
                        256.0L * std::numeric_limits<long double>::epsilon() * std::max(1.0L, std::abs(direct)))
                        << "powers " << a << ',' << b << ',' << c;
                }
            }
        }
    }
}

TEST(TetrahedronDuffyQualification, RepresentativeGaussRules)
{
    expect_representative_rules(false);
}

TEST(TetrahedronDuffyQualification, RepresentativeLobattoRules)
{
    expect_representative_rules(true);
}

TEST(TetrahedronDuffyQualification, RepresentativeMixedFactors)
{
    expect_mixed_rules(false);
}

TEST(TetrahedronDuffyQualification, RejectsRequestsBeyondEachMixedDirectionalLimit)
{
    for (int combination = 0; combination < 8; ++combination) {
        SCOPED_TRACE(combination);
        const auto u = make_line(combination & 4, max_line_exactness(combination & 4));
        const auto v = make_line(combination & 2, max_line_exactness(combination & 2));
        const auto w = make_line(combination & 1, max_line_exactness(combination & 1));
        const int limit = maximum_mixed_request(combination);
        expect_invalid([&] { (void)make_tetrahedron_rule(limit + 1, u, v, w); },
            limit == max_tetrahedron_exactness() ? "requested_exactness" : "exactness");
    }
}

// Required Release qualification, separate from routine CTest:
// build/svMultiPhysics-build/bin/run_all_unit_tests --gtest_also_run_disabled_tests
//   --gtest_filter='TetrahedronDuffyQualification.DISABLED_*'
TEST(TetrahedronDuffyQualification, DISABLED_AllGaussMonomials)
{
    expect_every_rule(false);
}

TEST(TetrahedronDuffyQualification, DISABLED_AllLobattoMonomials)
{
    expect_every_rule(true);
}

TEST(TetrahedronDuffyQualification, DISABLED_MixedFactorMonomials)
{
    expect_mixed_rules(true);
}
