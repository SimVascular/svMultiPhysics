/**
 * @file test_HigherOrderWedge.cpp
 * @brief Focused higher-order wedge checks for LagrangeBasis.
 */

#include <gtest/gtest.h>

#include "FE/Basis/LagrangeBasis.h"
#include "FE/Basis/NodeOrderingConventions.h"

#include <cmath>
#include <vector>

using namespace svmp::FE;
using namespace svmp::FE::basis;

namespace {

void expect_nodes_close(const std::vector<math::Vector<Real, 3>>& lhs,
                        const std::vector<math::Vector<Real, 3>>& rhs,
                        Real tol)
{
    ASSERT_EQ(lhs.size(), rhs.size());
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        EXPECT_NEAR(lhs[i][0], rhs[i][0], tol) << "node " << i;
        EXPECT_NEAR(lhs[i][1], rhs[i][1], tol) << "node " << i;
        EXPECT_NEAR(lhs[i][2], rhs[i][2], tol) << "node " << i;
    }
}

void expect_kronecker_at_nodes(const LagrangeBasis& basis, Real tol)
{
    const auto& nodes = basis.nodes();
    ASSERT_EQ(nodes.size(), basis.size());

    std::vector<Real> values;
    for (std::size_t node = 0; node < nodes.size(); ++node) {
        basis.evaluate_values(nodes[node], values);
        ASSERT_EQ(values.size(), basis.size());
        for (std::size_t i = 0; i < values.size(); ++i) {
            const Real expected = (i == node) ? Real(1) : Real(0);
            EXPECT_NEAR(values[i], expected, tol)
                << "node " << node << ", basis " << i;
        }
    }
}

void expect_partition_gradient_hessian_sums(const LagrangeBasis& basis,
                                            const std::vector<math::Vector<Real, 3>>& points,
                                            Real value_tol,
                                            Real derivative_tol)
{
    for (const auto& xi : points) {
        std::vector<Real> values;
        std::vector<Gradient> gradients;
        std::vector<Hessian> hessians;
        basis.evaluate_all(xi, values, gradients, hessians);

        Real value_sum = Real(0);
        Gradient gradient_sum = Gradient::Zero();
        Hessian hessian_sum = Hessian::Zero();
        for (std::size_t i = 0; i < values.size(); ++i) {
            value_sum += values[i];
            for (std::size_t d = 0; d < 3u; ++d) {
                gradient_sum[d] += gradients[i][d];
                for (std::size_t e = 0; e < 3u; ++e) {
                    hessian_sum(d, e) += hessians[i](d, e);
                }
            }
        }

        EXPECT_NEAR(value_sum, Real(1), value_tol);
        for (int d = 0; d < basis.dimension(); ++d) {
            EXPECT_NEAR(gradient_sum[static_cast<std::size_t>(d)], Real(0), derivative_tol);
            for (int e = 0; e < basis.dimension(); ++e) {
                EXPECT_NEAR(hessian_sum(static_cast<std::size_t>(d),
                                        static_cast<std::size_t>(e)),
                            Real(0),
                            derivative_tol);
            }
        }
    }
}

void expect_all_entries_finite(const LagrangeBasis& basis,
                               const math::Vector<Real, 3>& xi)
{
    std::vector<Real> values;
    std::vector<Gradient> gradients;
    std::vector<Hessian> hessians;
    basis.evaluate_all(xi, values, gradients, hessians);

    for (std::size_t i = 0; i < values.size(); ++i) {
        EXPECT_TRUE(std::isfinite(static_cast<double>(values[i]))) << "value " << i;
        for (std::size_t d = 0; d < 3u; ++d) {
            EXPECT_TRUE(std::isfinite(static_cast<double>(gradients[i][d])))
                << "gradient " << i << ", " << d;
            for (std::size_t e = 0; e < 3u; ++e) {
                EXPECT_TRUE(std::isfinite(static_cast<double>(hessians[i](d, e))))
                    << "hessian " << i << ", " << d << ", " << e;
            }
        }
    }
}

} // namespace

TEST(HigherOrderWedge, CompleteAliasMatchesGeneratedNodeLayout) {
    LagrangeBasis alias_basis(ElementType::Wedge18, 1);
    const auto generated =
        ReferenceNodeLayout::get_lagrange_node_coords(ElementType::Wedge6, 2);

    ASSERT_EQ(generated.size(), ReferenceNodeLayout::num_nodes(ElementType::Wedge18));
    EXPECT_EQ(alias_basis.element_type(), ElementType::Wedge6);
    EXPECT_EQ(alias_basis.order(), 2);
    expect_nodes_close(alias_basis.nodes(), generated, Real(1e-14));
}

TEST(HigherOrderWedge, OrderThreeIsNodalAndPartitionsUnity) {
    LagrangeBasis wedge(ElementType::Wedge6, 3);

    expect_kronecker_at_nodes(wedge, Real(2e-10));
    expect_partition_gradient_hessian_sums(
        wedge,
        {
            {Real(0.18), Real(0.22), Real(-0.2)},
            {Real(0.12), Real(0.16), Real(0.1)},
            {Real(0.25), Real(0.15), Real(0.45)},
        },
        Real(1e-12),
        Real(1e-9));
}

TEST(HigherOrderWedge, OrderFourEvaluationsRemainFinite) {
    LagrangeBasis wedge(ElementType::Wedge6, 4);

    expect_all_entries_finite(wedge, {Real(0.2), Real(0.1), Real(-0.6)});
    expect_all_entries_finite(wedge, {Real(0.05), Real(0.8), Real(0.3)});
}

// Finiteness alone cannot detect a wrong triangle-index or axis-index lookup;
// the Kronecker property validates the order-four node lattice and its inverse
// index mapping end to end.
TEST(HigherOrderWedge, OrderFourIsNodalAndPartitionsUnity) {
    LagrangeBasis wedge(ElementType::Wedge6, 4);

    // Order-4 wedge = triangle(order 4) x line(order 4) = 15 x 5 nodes.
    EXPECT_EQ(wedge.size(), 15u * 5u);
    expect_kronecker_at_nodes(wedge, Real(1e-9));
    expect_partition_gradient_hessian_sums(
        wedge,
        {
            {Real(0.18), Real(0.22), Real(-0.2)},
            {Real(0.25), Real(0.15), Real(0.45)},
        },
        Real(1e-12),
        Real(1e-7));
}
