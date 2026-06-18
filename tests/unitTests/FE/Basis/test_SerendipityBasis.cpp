/**
 * @file test_SerendipityBasis.cpp
 * @brief Nodal-delta, partition-of-unity, and polynomial-reproduction tests for SerendipityBasis.
 */

#include <gtest/gtest.h>

#include "FE/Basis/LagrangeBasis.h"
#include "FE/Basis/NodeOrderingConventions.h"
#include "FE/Basis/SerendipityBasis.h"
#include "FE/Math/DenseLinearAlgebra.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

using namespace svmp::FE;
using namespace svmp::FE::basis;

namespace {

void expect_partition_of_unity(const SerendipityBasis& basis,
                               const math::Vector<Real, 3>& xi,
                               Real tolerance = Real(1e-10))
{
    std::vector<Real> values;
    std::vector<Gradient> gradients;
    basis.evaluate_values(xi, values);
    basis.evaluate_gradients(xi, gradients);

    Real value_sum = Real(0);
    Gradient gradient_sum = Gradient::Zero();
    for (std::size_t i = 0; i < values.size(); ++i) {
        value_sum += values[i];
        for (std::size_t component = 0; component < 3u; ++component) {
            gradient_sum[component] += gradients[i][component];
        }
    }

    EXPECT_NEAR(value_sum, Real(1), tolerance);
    for (int component = 0; component < basis.dimension(); ++component) {
        EXPECT_NEAR(gradient_sum[static_cast<std::size_t>(component)],
                    Real(0),
                    tolerance);
    }
}

void expect_nodal_delta(const SerendipityBasis& basis,
                        const std::vector<math::Vector<Real, 3>>& nodes,
                        Real tolerance)
{
    ASSERT_EQ(nodes.size(), basis.size());
    for (std::size_t node = 0; node < nodes.size(); ++node) {
        std::vector<Real> values;
        basis.evaluate_values(nodes[node], values);
        ASSERT_EQ(values.size(), basis.size());
        for (std::size_t dof = 0; dof < values.size(); ++dof) {
            EXPECT_NEAR(values[dof], dof == node ? Real(1) : Real(0), tolerance)
                << "node=" << node << " dof=" << dof;
        }
    }
}

std::vector<math::Vector<Real, 3>> reference_nodes(ElementType type,
                                                   std::size_t count)
{
    std::vector<math::Vector<Real, 3>> nodes;
    nodes.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        nodes.push_back(ReferenceNodeLayout::get_node_coords(type, i));
    }
    return nodes;
}

template<typename Function>
Real interpolate_nodal_function(const SerendipityBasis& basis,
                                const math::Vector<Real, 3>& xi,
                                Function&& nodal_function)
{
    std::vector<Real> values;
    basis.evaluate_values(xi, values);

    Real result = Real(0);
    const auto& nodes = basis.nodes();
    for (std::size_t i = 0; i < values.size(); ++i) {
        result += values[i] * nodal_function(nodes[i]);
    }
    return result;
}

int quad_serendipity_superlinear_degree_for_test(int ax, int ay) {
    return (ax > 1 ? ax : 0) + (ay > 1 ? ay : 0);
}

std::vector<std::array<int, 2>> quad_serendipity_exponents_for_test(int order) {
    std::vector<std::array<int, 2>> exponents;
    for (int ay = 0; ay <= order; ++ay) {
        for (int ax = 0; ax <= order; ++ax) {
            if (quad_serendipity_superlinear_degree_for_test(ax, ay) <= order) {
                exponents.push_back({ax, ay});
            }
        }
    }
    return exponents;
}

std::size_t expected_quad_serendipity_size(int order) {
    const auto p = static_cast<std::size_t>(order);
    const std::size_t boundary = 4u * p;
    if (order < 4) {
        return boundary;
    }
    const auto m = static_cast<std::size_t>(order - 4);
    return boundary + (m + 1u) * (m + 2u) / 2u;
}

Real integer_power_for_test(Real base, int exponent) {
    Real result = Real(1);
    for (int k = 0; k < exponent; ++k) {
        result *= base;
    }
    return result;
}

Real monomial_value_for_test(const math::Vector<Real, 3>& p,
                             const std::array<int, 2>& exponent) {
    return integer_power_for_test(p[0], exponent[0]) *
           integer_power_for_test(p[1], exponent[1]);
}

std::vector<Real> quadrilateral_vandermonde_for_test(
    const std::vector<math::Vector<Real, 3>>& nodes,
    const std::vector<std::array<int, 2>>& exponents)
{
    const std::size_t n = nodes.size();
    std::vector<Real> vandermonde(n * n, Real(0));
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            vandermonde[row * n + col] =
                monomial_value_for_test(nodes[row], exponents[col]);
        }
    }
    return vandermonde;
}

void expect_no_duplicate_nodes(const std::vector<math::Vector<Real, 3>>& nodes,
                               Real tolerance)
{
    for (std::size_t a = 0; a < nodes.size(); ++a) {
        for (std::size_t b = a + 1u; b < nodes.size(); ++b) {
            const Real dx = std::abs(nodes[a][0] - nodes[b][0]);
            const Real dy = std::abs(nodes[a][1] - nodes[b][1]);
            EXPECT_GT(std::max(dx, dy), tolerance)
                << "duplicate nodes " << a << " and " << b;
        }
    }
}

void expect_nodes_near(const std::vector<math::Vector<Real, 3>>& actual,
                       const std::vector<math::Vector<Real, 3>>& expected,
                       Real tolerance)
{
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i) {
        for (std::size_t d = 0; d < 3u; ++d) {
            EXPECT_NEAR(actual[i][d], expected[i][d], tolerance)
                << "node=" << i << " component=" << d;
        }
    }
}

// Every monomial here has superlinear degree at most three, so it lies in the
// order-three quadrilateral serendipity space.
Real cubic_serendipity_function(const math::Vector<Real, 3>& p) {
    const Real x = p[0];
    const Real y = p[1];
    return Real(1) + Real(2) * x - y + Real(3) * x * y +
           x * x * x - Real(2) * y * y * y +
           Real(0.5) * x * x * x * y - Real(0.25) * x * y * y * y;
}

Real bilinear_function(const math::Vector<Real, 3>& p) {
    return Real(2) - Real(3) * p[0] + Real(4) * p[1] + Real(0.5) * p[0] * p[1];
}

} // namespace

TEST(SerendipityBasis, Quad8IsNodalAndPartitionsUnity) {
    SerendipityBasis basis(ElementType::Quad8, 2);
    SerendipityBasis explicit_quad4_basis(ElementType::Quad4, 2);

    EXPECT_EQ(basis.size(), 8u);
    expect_nodes_near(basis.nodes(), explicit_quad4_basis.nodes(), Real(1e-14));
    expect_nodal_delta(basis, basis.nodes(), Real(1e-10));
    expect_partition_of_unity(basis, {Real(0.17), Real(-0.31), Real(0)});
}

TEST(SerendipityBasis, Hex20IsNodalAndPartitionsUnity) {
    SerendipityBasis basis(ElementType::Hex20, 2);

    EXPECT_EQ(basis.size(), 20u);
    expect_nodal_delta(basis,
                       reference_nodes(ElementType::Hex20, basis.size()),
                       Real(1e-10));
    expect_partition_of_unity(basis, {Real(0.2), Real(-0.1), Real(0.3)});
}

TEST(SerendipityBasis, Wedge15IsNodalAndPartitionsUnity) {
    SerendipityBasis basis(ElementType::Wedge15, 2);

    EXPECT_EQ(basis.size(), 15u);
    expect_nodal_delta(basis,
                       reference_nodes(ElementType::Wedge15, basis.size()),
                       Real(1e-9));
    expect_partition_of_unity(basis, {Real(0.2), Real(0.3), Real(0.1)});
}

TEST(SerendipityBasis, RejectsUnsupportedSerendipityAliases) {
    EXPECT_THROW(SerendipityBasis(ElementType::Quad9, 2), FEException);
    EXPECT_THROW(SerendipityBasis(ElementType::Pyramid13, 2), FEException);
    EXPECT_THROW(SerendipityBasis(ElementType::Pyramid14, 2), FEException);
    EXPECT_THROW(SerendipityBasis(ElementType::Quad8, 3), FEException);
    EXPECT_THROW(SerendipityBasis(ElementType::Quad8, 1), FEException);
}

TEST(SerendipityBasis, QuadrilateralOrderZeroNormalizesToLinear) {
    SerendipityBasis basis(ElementType::Quad4, 0);

    EXPECT_EQ(basis.order(), 1);
    EXPECT_EQ(basis.size(), 4u);
    expect_nodal_delta(basis, basis.nodes(), Real(1e-12));
}

// Explicit Quad4 serendipity orders run the documented monomial selection,
// boundary plus triangular interior node placement, and runtime Vandermonde
// inversion. Order four is the first order with an interior residual
// polynomial, so it is the first order that appends an interior node.
TEST(SerendipityBasis, QuadrilateralOrdersOneThreeFourAreNodalAndPartitionUnity) {
    const struct Case {
        int order;
        std::size_t size;
    } cases[] = {
        {1, 4u},
        {3, 12u},
        {4, 17u},
    };

    for (const auto& c : cases) {
        SerendipityBasis basis(ElementType::Quad4, c.order);
        EXPECT_EQ(basis.size(), c.size) << "order=" << c.order;
        EXPECT_EQ(basis.order(), c.order);
        EXPECT_EQ(basis.dimension(), 2);
        ASSERT_EQ(basis.nodes().size(), c.size);

        for (const auto& node : basis.nodes()) {
            EXPECT_LE(std::abs(node[0]), Real(1));
            EXPECT_LE(std::abs(node[1]), Real(1));
        }

        expect_nodal_delta(basis, basis.nodes(), Real(1e-9));
        expect_partition_of_unity(basis, {Real(0.17), Real(-0.31), Real(0)}, Real(1e-9));
        expect_partition_of_unity(basis, {Real(-0.45), Real(0.25), Real(0)}, Real(1e-9));
    }
}

TEST(SerendipityBasis, QuadrilateralNodesFollowDocumentedConstructionThroughOrderTen) {
    constexpr Real kTol = Real(1e-14);

    for (int order = 1; order <= 10; ++order) {
        SerendipityBasis basis(ElementType::Quad4, order);
        const auto& nodes = basis.nodes();
        const std::size_t expected_size = expected_quad_serendipity_size(order);
        const std::size_t boundary_count = static_cast<std::size_t>(4 * order);

        ASSERT_EQ(basis.size(), expected_size) << "order=" << order;
        ASSERT_EQ(nodes.size(), expected_size) << "order=" << order;
        EXPECT_EQ(quad_serendipity_exponents_for_test(order).size(),
                  expected_size) << "order=" << order;
        expect_no_duplicate_nodes(nodes, kTol);

        for (std::size_t i = 0; i < nodes.size(); ++i) {
            EXPECT_NEAR(nodes[i][2], Real(0), kTol) << "order=" << order
                                                    << " node=" << i;
            EXPECT_LE(std::abs(nodes[i][0]), Real(1)) << "order=" << order
                                                       << " node=" << i;
            EXPECT_LE(std::abs(nodes[i][1]), Real(1)) << "order=" << order
                                                       << " node=" << i;

            const bool on_boundary =
                std::abs(std::abs(nodes[i][0]) - Real(1)) <= kTol ||
                std::abs(std::abs(nodes[i][1]) - Real(1)) <= kTol;
            if (i < boundary_count) {
                EXPECT_TRUE(on_boundary) << "order=" << order << " node=" << i;
            } else {
                EXPECT_FALSE(on_boundary) << "order=" << order << " node=" << i;
                EXPECT_LT(std::abs(nodes[i][0]), Real(1)) << "order=" << order
                                                           << " node=" << i;
                EXPECT_LT(std::abs(nodes[i][1]), Real(1)) << "order=" << order
                                                           << " node=" << i;
            }
        }

        std::size_t index = boundary_count;
        if (order >= 4) {
            const int m = order - 4;
            const Real y_denominator = Real(m + 2);
            for (int row = 0; row <= m; ++row) {
                const int row_count = m + 1 - row;
                const Real expected_y =
                    Real(-1) + Real(2) * Real(row + 1) / y_denominator;
                const Real x_denominator = Real(row_count + 1);
                for (int col = 0; col < row_count; ++col) {
                    ASSERT_LT(index, nodes.size());
                    const Real expected_x =
                        Real(-1) + Real(2) * Real(col + 1) / x_denominator;
                    EXPECT_NEAR(nodes[index][0], expected_x, kTol)
                        << "order=" << order << " row=" << row << " col=" << col;
                    EXPECT_NEAR(nodes[index][1], expected_y, kTol)
                        << "order=" << order << " row=" << row << " col=" << col;
                    ++index;
                }
            }
        }
        EXPECT_EQ(index, nodes.size()) << "order=" << order;
    }
}

TEST(SerendipityBasis, QuadrilateralOrderOneReproducesBilinearFunctions) {
    SerendipityBasis basis(ElementType::Quad4, 1);

    const std::vector<math::Vector<Real, 3>> points = {
        {Real(0.25), Real(-0.4), Real(0)},
        {Real(-0.7), Real(0.6), Real(0)},
    };
    for (const auto& xi : points) {
        EXPECT_NEAR(interpolate_nodal_function(basis, xi, bilinear_function),
                    bilinear_function(xi),
                    Real(1e-12));
    }
}

TEST(SerendipityBasis, QuadrilateralOrderThreeReproducesSerendipityCubics) {
    SerendipityBasis basis(ElementType::Quad4, 3);

    const std::vector<math::Vector<Real, 3>> points = {
        {Real(0.25), Real(-0.4), Real(0)},
        {Real(-0.7), Real(0.6), Real(0)},
    };
    for (const auto& xi : points) {
        EXPECT_NEAR(interpolate_nodal_function(basis, xi, cubic_serendipity_function),
                    cubic_serendipity_function(xi),
                    Real(1e-11));
    }
}

TEST(SerendipityBasis, QuadrilateralOrdersReproduceEverySerendipityMonomial) {
    const std::vector<math::Vector<Real, 3>> points = {
        {Real(0.25), Real(-0.4), Real(0)},
        {Real(-0.7), Real(0.6), Real(0)},
        {Real(0.11), Real(0.23), Real(0)},
    };

    for (int order = 1; order <= 10; ++order) {
        SerendipityBasis basis(ElementType::Quad4, order);
        const auto exponents = quad_serendipity_exponents_for_test(order);
        ASSERT_EQ(exponents.size(), basis.size()) << "order=" << order;

        const Real tolerance = (order <= 7) ? Real(1e-10) : Real(2e-8);
        for (const auto& exponent : exponents) {
            for (const auto& xi : points) {
                const Real interpolated =
                    interpolate_nodal_function(
                        basis,
                        xi,
                        [&exponent](const math::Vector<Real, 3>& node) {
                            return monomial_value_for_test(node, exponent);
                        });
                const Real expected = monomial_value_for_test(xi, exponent);
                EXPECT_NEAR(interpolated, expected, tolerance)
                    << "order=" << order << " ax=" << exponent[0]
                    << " ay=" << exponent[1] << " xi=(" << xi[0] << ","
                    << xi[1] << ")";
            }
        }
    }
}

TEST(SerendipityBasis, QuadrilateralVandermondeHasFullRankThroughOrderTen) {
    for (int order = 1; order <= 10; ++order) {
        SerendipityBasis basis(ElementType::Quad4, order);
        const auto exponents = quad_serendipity_exponents_for_test(order);
        const auto vandermonde =
            quadrilateral_vandermonde_for_test(basis.nodes(), exponents);
        const std::size_t n = basis.size();

        ASSERT_EQ(exponents.size(), n) << "order=" << order;
        ASSERT_EQ(vandermonde.size(), n * n) << "order=" << order;
        EXPECT_EQ(math::dense_matrix_rank(vandermonde, n, n), n)
            << "order=" << order;
    }
}

// SerendipityBasis(Hex8, 1) is the only route to the hand-written trilinear
// corner evaluator (values, gradients, and Hessians); it must agree with the
// trilinear Lagrange basis on the same element.
TEST(SerendipityBasis, TrilinearHexMatchesLagrangeHex8) {
    SerendipityBasis serendipity(ElementType::Hex8, 1);
    LagrangeBasis lagrange(ElementType::Hex8, 1);

    EXPECT_EQ(serendipity.size(), 8u);
    EXPECT_EQ(serendipity.dimension(), 3);
    expect_nodal_delta(serendipity,
                       reference_nodes(ElementType::Hex8, serendipity.size()),
                       Real(1e-12));

    const std::vector<math::Vector<Real, 3>> points = {
        {Real(0.2), Real(-0.1), Real(0.3)},
        {Real(-0.35), Real(0.25), Real(-0.15)},
    };
    for (const auto& xi : points) {
        std::vector<Real> s_values;
        std::vector<Real> l_values;
        std::vector<Gradient> s_gradients;
        std::vector<Gradient> l_gradients;
        std::vector<Hessian> s_hessians;
        std::vector<Hessian> l_hessians;
        serendipity.evaluate_all(xi, s_values, s_gradients, s_hessians);
        lagrange.evaluate_all(xi, l_values, l_gradients, l_hessians);

        ASSERT_EQ(s_values.size(), l_values.size());
        for (std::size_t i = 0; i < s_values.size(); ++i) {
            EXPECT_NEAR(s_values[i], l_values[i], Real(1e-13));
            for (std::size_t d = 0; d < 3u; ++d) {
                EXPECT_NEAR(s_gradients[i][d], l_gradients[i][d], Real(1e-13));
                for (std::size_t e = 0; e < 3u; ++e) {
                    EXPECT_NEAR(s_hessians[i](d, e), l_hessians[i](d, e), Real(1e-13));
                }
            }
        }
    }
}
