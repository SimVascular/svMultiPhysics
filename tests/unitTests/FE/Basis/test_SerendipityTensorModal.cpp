/**
 * @file test_SerendipityTensorModal.cpp
 * @brief Tests for the migrated Serendipity basis subset.
 */

#include <gtest/gtest.h>

#include "FE/Basis/LagrangeBasis.h"
#include "FE/Basis/NodeOrderingConventions.h"
#include "FE/Basis/SerendipityBasis.h"

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

    EXPECT_EQ(basis.size(), 8u);
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
}

// Orders other than two run the generic quadrilateral path: serendipity
// monomial selection, boundary plus interior node placement, and a runtime
// Vandermonde inversion whose unisolvence is assumed rather than tabulated.
// Order four is the first order that selects an interior node.
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
