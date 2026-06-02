/**
 * @file test_SerendipityTensorModal.cpp
 * @brief Tests for the migrated Serendipity basis subset.
 */

#include <gtest/gtest.h>

#include "FE/Basis/NodeOrderingConventions.h"
#include "FE/Basis/SerendipityBasis.h"

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
    Gradient gradient_sum{};
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

TEST(SerendipityBasis, Pyramid13IsNodalAndPartitionsUnity) {
    SerendipityBasis basis(ElementType::Pyramid13, 2);

    EXPECT_EQ(basis.size(), 13u);
    expect_nodal_delta(basis,
                       reference_nodes(ElementType::Pyramid13, basis.size()),
                       Real(1e-8));
    expect_partition_of_unity(basis, {Real(0.1), Real(-0.2), Real(0.4)});
}

TEST(SerendipityBasis, RejectsUnsupportedSerendipityAliases) {
    EXPECT_THROW(SerendipityBasis(ElementType::Quad9, 2), FEException);
    EXPECT_THROW(SerendipityBasis(ElementType::Pyramid14, 2), FEException);
    EXPECT_THROW(SerendipityBasis(ElementType::Quad8, 3), FEException);
}

