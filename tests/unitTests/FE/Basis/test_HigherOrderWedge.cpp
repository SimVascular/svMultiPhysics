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

void expect_nodes_close(const std::vector<math::Vector<double, 3>>& lhs,
                        const std::vector<math::Vector<double, 3>>& rhs,
                        double tol)
{
    ASSERT_EQ(lhs.size(), rhs.size());
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        EXPECT_NEAR(lhs[i][0], rhs[i][0], tol) << "node " << i;
        EXPECT_NEAR(lhs[i][1], rhs[i][1], tol) << "node " << i;
        EXPECT_NEAR(lhs[i][2], rhs[i][2], tol) << "node " << i;
    }
}

void expect_kronecker_at_nodes(const LagrangeBasis& basis, double tol)
{
    const auto& nodes = basis.nodes();
    ASSERT_EQ(nodes.size(), basis.size());

    std::vector<double> values;
    for (std::size_t node = 0; node < nodes.size(); ++node) {
        basis.evaluate_values(nodes[node], values);
        ASSERT_EQ(values.size(), basis.size());
        for (std::size_t i = 0; i < values.size(); ++i) {
            const double expected = (i == node) ? double(1) : double(0);
            EXPECT_NEAR(values[i], expected, tol)
                << "node " << node << ", basis " << i;
        }
    }
}

void expect_partition_gradient_hessian_sums(const LagrangeBasis& basis,
                                            const std::vector<math::Vector<double, 3>>& points,
                                            double value_tol,
                                            double derivative_tol)
{
    for (const auto& xi : points) {
        std::vector<double> values;
        std::vector<Gradient> gradients;
        std::vector<Hessian> hessians;
        basis.evaluate_all(xi, values, gradients, hessians);

        double value_sum = double(0);
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

        EXPECT_NEAR(value_sum, double(1), value_tol);
        for (int d = 0; d < basis.dimension(); ++d) {
            EXPECT_NEAR(gradient_sum[static_cast<std::size_t>(d)], double(0), derivative_tol);
            for (int e = 0; e < basis.dimension(); ++e) {
                EXPECT_NEAR(hessian_sum(static_cast<std::size_t>(d),
                                        static_cast<std::size_t>(e)),
                            double(0),
                            derivative_tol);
            }
        }
    }
}

void expect_all_entries_finite(const LagrangeBasis& basis,
                               const math::Vector<double, 3>& xi)
{
    std::vector<double> values;
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
    LagrangeBasis alias_basis(ElementType::Wedge18, 2);
    const auto generated =
        ReferenceNodeLayout::get_lagrange_node_coords(ElementType::Wedge6, 2);

    ASSERT_EQ(generated.size(), ReferenceNodeLayout::num_nodes(ElementType::Wedge18));
    EXPECT_EQ(alias_basis.topology(), BasisTopology::Wedge);
    EXPECT_EQ(named_element_for(alias_basis.topology(), alias_basis.order(), alias_basis.basis_type()),
              ElementType::Wedge18);  // faithful round-trip
    EXPECT_EQ(alias_basis.order(), 2);
    expect_nodes_close(alias_basis.nodes(), generated, double(1e-14));
}

TEST(HigherOrderWedge, OrderThreeIsNodalAndPartitionsUnity) {
    LagrangeBasis wedge(BasisTopology::Wedge, 3);

    expect_kronecker_at_nodes(wedge, double(2e-10));
    expect_partition_gradient_hessian_sums(
        wedge,
        {
            {double(0.18), double(0.22), double(-0.2)},
            {double(0.12), double(0.16), double(0.1)},
            {double(0.25), double(0.15), double(0.45)},
        },
        double(1e-12),
        double(1e-9));
}

TEST(HigherOrderWedge, OrderFourEvaluationsRemainFinite) {
    LagrangeBasis wedge(BasisTopology::Wedge, 4);

    expect_all_entries_finite(wedge, {double(0.2), double(0.1), double(-0.6)});
    expect_all_entries_finite(wedge, {double(0.05), double(0.8), double(0.3)});
}

// Finiteness alone cannot detect a wrong triangle-index or axis-index lookup;
// the Kronecker property validates the order-four node lattice and its inverse
// index mapping end to end.
TEST(HigherOrderWedge, OrderFourIsNodalAndPartitionsUnity) {
    LagrangeBasis wedge(BasisTopology::Wedge, 4);

    // Order-4 wedge = triangle(order 4) x line(order 4) = 15 x 5 nodes.
    EXPECT_EQ(wedge.size(), 15u * 5u);
    expect_kronecker_at_nodes(wedge, double(1e-9));
    expect_partition_gradient_hessian_sums(
        wedge,
        {
            {double(0.18), double(0.22), double(-0.2)},
            {double(0.25), double(0.15), double(0.45)},
        },
        double(1e-12),
        double(1e-7));
}
