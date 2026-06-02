/**
 * @file test_HigherOrderWedgePyramid.cpp
 * @brief Focused higher-order wedge and pyramid checks for LagrangeBasis.
 */

#include <gtest/gtest.h>

#include "FE/Basis/LagrangeBasis.h"
#include "FE/Basis/NodeOrderingConventions.h"

#include <cmath>
#include <tuple>
#include <utility>
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
        Gradient gradient_sum{};
        Hessian hessian_sum{};
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

TEST(HigherOrderWedgePyramid, CompleteAliasesMatchGeneratedNodeLayouts) {
    const std::vector<std::tuple<ElementType, ElementType, int>> cases = {
        {ElementType::Wedge18, ElementType::Wedge6, 2},
        {ElementType::Pyramid14, ElementType::Pyramid5, 2},
    };

    for (const auto& [alias, canonical, order] : cases) {
        LagrangeBasis alias_basis(alias, order);
        const auto generated = ReferenceNodeLayout::get_lagrange_node_coords(canonical, order);
        ASSERT_EQ(generated.size(), ReferenceNodeLayout::num_nodes(alias));
        expect_nodes_close(alias_basis.nodes(), generated, Real(1e-14));

        for (std::size_t i = 0; i < generated.size(); ++i) {
            const auto public_node = ReferenceNodeLayout::get_node_coords(alias, i);
            EXPECT_NEAR(public_node[0], generated[i][0], Real(1e-14)) << "node " << i;
            EXPECT_NEAR(public_node[1], generated[i][1], Real(1e-14)) << "node " << i;
            EXPECT_NEAR(public_node[2], generated[i][2], Real(1e-14)) << "node " << i;
        }
    }
}

TEST(HigherOrderWedgePyramid, WedgeOrderThreeIsNodalAndPartitionsUnity) {
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

TEST(HigherOrderWedgePyramid, PyramidOrderThreeIsNodalAndPartitionsUnity) {
    LagrangeBasis pyramid(ElementType::Pyramid5, 3);

    expect_kronecker_at_nodes(pyramid, Real(5e-8));
    expect_partition_gradient_hessian_sums(
        pyramid,
        {
            {Real(0), Real(0), Real(0.2)},
            {Real(0.12), Real(-0.08), Real(0.24)},
            {Real(-0.08), Real(0.1), Real(0.55)},
        },
        Real(1e-11),
        Real(5e-7));
}

TEST(HigherOrderWedgePyramid, PyramidNearApexDerivativeQueriesRemainFinite) {
    const std::vector<std::pair<ElementType, int>> cases = {
        {ElementType::Pyramid5, 1},
        {ElementType::Pyramid14, 2},
        {ElementType::Pyramid5, 4},
    };

    for (const auto& [type, order] : cases) {
        LagrangeBasis basis(type, order);
        expect_all_entries_finite(basis, {Real(0.01), Real(-0.005), Real(0.92)});
        expect_all_entries_finite(basis, {Real(-0.004), Real(0.007), Real(0.98)});
    }
}
