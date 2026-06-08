/**
 * @file test_LagrangeBasis.cpp
 * @brief Unit tests for the reduced scalar Lagrange basis implementation.
 */

#include <gtest/gtest.h>

#include "FE/Basis/BasisExceptions.h"
#include "FE/Basis/BasisFactory.h"
#include "FE/Basis/LagrangeBasis.h"
#include "FE/Basis/NodeOrderingConventions.h"

#include <array>
#include <tuple>
#include <vector>

using namespace svmp::FE;
using namespace svmp::FE::basis;

namespace {

using Point = math::Vector<Real, 3>;

struct CanonicalCase {
    ElementType type;
    int order;
    std::size_t size;
    int dimension;
    std::vector<Point> points;
    Real derivative_tol;
};

const std::vector<CanonicalCase>& canonical_cases() {
    static const std::vector<CanonicalCase> cases = {
        {ElementType::Line2, 3, 4u, 1,
         {{Real(-0.35), Real(0), Real(0)}, {Real(0.2), Real(0), Real(0)}},
         Real(1e-11)},
        {ElementType::Triangle3, 3, 10u, 2,
         {{Real(0.15), Real(0.2), Real(0)}, {Real(0.25), Real(0.1), Real(0)}},
         Real(1e-9)},
        {ElementType::Quad4, 3, 16u, 2,
         {{Real(0.2), Real(-0.3), Real(0)}, {Real(-0.45), Real(0.25), Real(0)}},
         Real(1e-11)},
        {ElementType::Tetra4, 2, 10u, 3,
         {{Real(0.12), Real(0.18), Real(0.16)}, {Real(0.2), Real(0.1), Real(0.18)}},
         Real(1e-9)},
        {ElementType::Hex8, 2, 27u, 3,
         {{Real(0.1), Real(-0.2), Real(0.3)}, {Real(-0.35), Real(0.25), Real(-0.15)}},
         Real(1e-10)},
        {ElementType::Wedge6, 2, 18u, 3,
         {{Real(0.18), Real(0.22), Real(-0.2)}, {Real(0.12), Real(0.16), Real(0.1)}},
         Real(1e-9)},
    };
    return cases;
}

std::vector<Point> sample_points_for(ElementType type) {
    for (const auto& c : canonical_cases()) {
        if (c.type == type) {
            return c.points;
        }
    }
    return {};
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
            EXPECT_NEAR(values[i], i == node ? Real(1) : Real(0), tol)
                << "node=" << node << " basis=" << i;
        }
    }
}

void expect_partition_gradient_hessian_sums(const LagrangeBasis& basis,
                                            const std::vector<Point>& points,
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

        EXPECT_NEAR(value_sum, Real(1), Real(1e-12));
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

void expect_raw_sinks_match_vector_evaluation(const LagrangeBasis& basis,
                                              const Point& xi)
{
    std::vector<Real> values;
    std::vector<Gradient> gradients;
    std::vector<Hessian> hessians;
    basis.evaluate_all(xi, values, gradients, hessians);

    std::vector<Real> flat_values(basis.size());
    std::vector<Real> flat_gradients(basis.size() * 3u);
    std::vector<Real> flat_hessians(basis.size() * 9u);
    basis.evaluate_values_to(xi, flat_values.data());
    basis.evaluate_gradients_to(xi, flat_gradients.data());
    basis.evaluate_hessians_to(xi, flat_hessians.data());

    for (std::size_t i = 0; i < basis.size(); ++i) {
        EXPECT_NEAR(flat_values[i], values[i], Real(1e-14));
        for (std::size_t d = 0; d < 3u; ++d) {
            EXPECT_NEAR(flat_gradients[i * 3u + d], gradients[i][d], Real(1e-14));
            for (std::size_t e = 0; e < 3u; ++e) {
                EXPECT_NEAR(flat_hessians[i * 9u + d * 3u + e],
                            hessians[i](d, e),
                            Real(1e-14));
            }
        }
    }
}

void expect_nodes_close(const std::vector<Point>& lhs,
                        const std::vector<Point>& rhs,
                        Real tol)
{
    ASSERT_EQ(lhs.size(), rhs.size());
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        EXPECT_NEAR(lhs[i][0], rhs[i][0], tol) << "node=" << i;
        EXPECT_NEAR(lhs[i][1], rhs[i][1], tol) << "node=" << i;
        EXPECT_NEAR(lhs[i][2], rhs[i][2], tol) << "node=" << i;
    }
}

void expect_evaluations_match(const LagrangeBasis& lhs,
                              const LagrangeBasis& rhs,
                              const std::vector<Point>& points,
                              Real tol)
{
    ASSERT_EQ(lhs.size(), rhs.size());

    for (const auto& xi : points) {
        std::vector<Real> lhs_values;
        std::vector<Real> rhs_values;
        std::vector<Gradient> lhs_gradients;
        std::vector<Gradient> rhs_gradients;
        std::vector<Hessian> lhs_hessians;
        std::vector<Hessian> rhs_hessians;

        lhs.evaluate_all(xi, lhs_values, lhs_gradients, lhs_hessians);
        rhs.evaluate_all(xi, rhs_values, rhs_gradients, rhs_hessians);

        for (std::size_t i = 0; i < lhs.size(); ++i) {
            EXPECT_NEAR(lhs_values[i], rhs_values[i], tol);
            for (std::size_t d = 0; d < 3u; ++d) {
                EXPECT_NEAR(lhs_gradients[i][d], rhs_gradients[i][d], tol);
                for (std::size_t e = 0; e < 3u; ++e) {
                    EXPECT_NEAR(lhs_hessians[i](d, e), rhs_hessians[i](d, e), tol);
                }
            }
        }
    }
}

Real linear_function(const Point& p) {
    return Real(2) + Real(3) * p[0] - Real(4) * p[1] + Real(5) * p[2];
}

Gradient linear_gradient() {
    Gradient g{};
    g[0] = Real(3);
    g[1] = Real(-4);
    g[2] = Real(5);
    return g;
}

Real quadratic_function(const Point& p) {
    return Real(1) + Real(2) * p[0] - p[1] + Real(0.5) * p[2] +
           p[0] * p[0] + Real(0.75) * p[1] * p[1] - Real(0.25) * p[2] * p[2] +
           Real(0.2) * p[0] * p[1] - Real(0.3) * p[0] * p[2] +
           Real(0.4) * p[1] * p[2];
}

template<typename Function>
Real interpolate_value(const LagrangeBasis& basis,
                       const std::vector<Real>& values,
                       Function&& nodal_function)
{
    Real result = Real(0);
    const auto& nodes = basis.nodes();
    for (std::size_t i = 0; i < values.size(); ++i) {
        result += values[i] * nodal_function(nodes[i]);
    }
    return result;
}

} // namespace

TEST(LagrangeBasis, CanonicalTopologiesHaveExpectedSizesAndDimensions) {
    for (const auto& c : canonical_cases()) {
        LagrangeBasis basis(c.type, c.order);
        EXPECT_EQ(basis.basis_type(), BasisType::Lagrange);
        EXPECT_EQ(basis.element_type(), c.type);
        EXPECT_EQ(basis.order(), c.order);
        EXPECT_EQ(basis.size(), c.size);
        EXPECT_EQ(basis.dimension(), c.dimension);
    }
}

TEST(LagrangeBasis, CanonicalTopologiesAreNodalAndPartitionUnity) {
    for (const auto& c : canonical_cases()) {
        LagrangeBasis basis(c.type, c.order);
        expect_kronecker_at_nodes(basis, Real(2e-10));
        expect_partition_gradient_hessian_sums(basis, c.points, c.derivative_tol);
    }
}

TEST(LagrangeBasis, RawOutputSinksMatchVectorEvaluationAcrossTopologies) {
    for (const auto& c : canonical_cases()) {
        LagrangeBasis basis(c.type, c.order);
        expect_raw_sinks_match_vector_evaluation(basis, c.points.front());
    }
}

TEST(LagrangeBasis, CompleteAliasesNormalizeToCanonicalBases) {
    const std::vector<std::tuple<ElementType, ElementType, int>> aliases = {
        {ElementType::Line3, ElementType::Line2, 2},
        {ElementType::Triangle6, ElementType::Triangle3, 2},
        {ElementType::Quad9, ElementType::Quad4, 2},
        {ElementType::Tetra10, ElementType::Tetra4, 2},
        {ElementType::Hex27, ElementType::Hex8, 2},
        {ElementType::Wedge18, ElementType::Wedge6, 2},
    };

    for (const auto& [alias, canonical, order] : aliases) {
        LagrangeBasis alias_basis(alias, 1);
        LagrangeBasis canonical_basis(canonical, order);
        const auto generated = ReferenceNodeLayout::get_lagrange_node_coords(canonical, order);

        EXPECT_EQ(alias_basis.element_type(), canonical);
        EXPECT_EQ(alias_basis.order(), order);
        expect_nodes_close(alias_basis.nodes(), generated, Real(1e-14));
        expect_nodes_close(alias_basis.nodes(), canonical_basis.nodes(), Real(1e-14));
        expect_evaluations_match(alias_basis,
                                 canonical_basis,
                                 sample_points_for(canonical),
                                 Real(1e-12));
    }
}

TEST(LagrangeBasis, NodeOrderingMatchesPublicAliasLayouts) {
    const std::vector<std::tuple<ElementType, ElementType, int>> aliases = {
        {ElementType::Line2, ElementType::Line2, 1},
        {ElementType::Line3, ElementType::Line2, 2},
        {ElementType::Triangle3, ElementType::Triangle3, 1},
        {ElementType::Triangle6, ElementType::Triangle3, 2},
        {ElementType::Quad4, ElementType::Quad4, 1},
        {ElementType::Quad9, ElementType::Quad4, 2},
        {ElementType::Tetra4, ElementType::Tetra4, 1},
        {ElementType::Tetra10, ElementType::Tetra4, 2},
        {ElementType::Hex8, ElementType::Hex8, 1},
        {ElementType::Hex27, ElementType::Hex8, 2},
        {ElementType::Wedge6, ElementType::Wedge6, 1},
        {ElementType::Wedge18, ElementType::Wedge6, 2},
    };

    for (const auto& [alias, canonical, order] : aliases) {
        const auto generated = ReferenceNodeLayout::get_lagrange_node_coords(canonical, order);
        ASSERT_EQ(generated.size(), ReferenceNodeLayout::num_nodes(alias));

        for (std::size_t i = 0; i < generated.size(); ++i) {
            const auto public_node = ReferenceNodeLayout::get_node_coords(alias, i);
            EXPECT_NEAR(public_node[0], generated[i][0], Real(1e-14)) << "node=" << i;
            EXPECT_NEAR(public_node[1], generated[i][1], Real(1e-14)) << "node=" << i;
            EXPECT_NEAR(public_node[2], generated[i][2], Real(1e-14)) << "node=" << i;
        }
    }
}

TEST(LagrangeBasis, RemovedOrSerendipityFamiliesAreRejected) {
    const std::array<ElementType, 6> unsupported = {
        ElementType::Quad8,
        ElementType::Hex20,
        ElementType::Wedge15,
        ElementType::Pyramid5,
        ElementType::Pyramid13,
        ElementType::Pyramid14,
    };

    for (const auto type : unsupported) {
        EXPECT_THROW((void)LagrangeBasis(type, 2), BasisElementCompatibilityException)
            << "element=" << static_cast<int>(type);
    }
}

TEST(LagrangeBasis, LinearPolynomialReproductionAcrossLinearTopologies) {
    const std::vector<std::pair<ElementType, Point>> cases = {
        {ElementType::Line2, {Real(-0.2), Real(0), Real(0)}},
        {ElementType::Triangle3, {Real(0.2), Real(0.3), Real(0)}},
        {ElementType::Quad4, {Real(0.25), Real(-0.4), Real(0)}},
        {ElementType::Tetra4, {Real(0.1), Real(0.2), Real(0.3)}},
        {ElementType::Hex8, {Real(0.15), Real(-0.2), Real(0.25)}},
        {ElementType::Wedge6, {Real(0.2), Real(0.15), Real(-0.3)}},
    };
    const Gradient expected_gradient = linear_gradient();

    for (const auto& [type, point] : cases) {
        LagrangeBasis basis(type, 1);
        std::vector<Real> values;
        std::vector<Gradient> gradients;
        basis.evaluate_values(point, values);
        basis.evaluate_gradients(point, gradients);

        const Real interpolated =
            interpolate_value(basis, values, linear_function);
        EXPECT_NEAR(interpolated, linear_function(point), Real(1e-12));

        Gradient interpolated_gradient{};
        for (std::size_t i = 0; i < gradients.size(); ++i) {
            const Real nodal_value = linear_function(basis.nodes()[i]);
            for (int d = 0; d < basis.dimension(); ++d) {
                interpolated_gradient[static_cast<std::size_t>(d)] +=
                    nodal_value * gradients[i][static_cast<std::size_t>(d)];
            }
        }
        for (int d = 0; d < basis.dimension(); ++d) {
            EXPECT_NEAR(interpolated_gradient[static_cast<std::size_t>(d)],
                        expected_gradient[static_cast<std::size_t>(d)],
                        Real(1e-12));
        }
    }
}

TEST(LagrangeBasis, QuadraticPolynomialReproductionAcrossQuadraticAliases) {
    const std::vector<std::pair<ElementType, Point>> cases = {
        {ElementType::Line3, {Real(-0.2), Real(0), Real(0)}},
        {ElementType::Triangle6, {Real(0.2), Real(0.3), Real(0)}},
        {ElementType::Quad9, {Real(0.25), Real(-0.4), Real(0)}},
        {ElementType::Tetra10, {Real(0.1), Real(0.2), Real(0.3)}},
        {ElementType::Hex27, {Real(0.15), Real(-0.2), Real(0.25)}},
        {ElementType::Wedge18, {Real(0.2), Real(0.15), Real(-0.3)}},
    };

    for (const auto& [type, point] : cases) {
        LagrangeBasis basis(type, 1);
        std::vector<Real> values;
        basis.evaluate_values(point, values);

        const Real interpolated =
            interpolate_value(basis, values, quadratic_function);
        EXPECT_NEAR(interpolated, quadratic_function(point), Real(5e-12))
            << "element=" << static_cast<int>(type);
    }
}

TEST(LagrangeBasis, FactoryCreatesReducedScalarBasisFamilies) {
    auto lagrange =
        basis_factory::create(BasisRequest{ElementType::Hex27, BasisType::Lagrange, 1});
    ASSERT_NE(lagrange, nullptr);
    EXPECT_EQ(lagrange->basis_type(), BasisType::Lagrange);
    EXPECT_EQ(lagrange->element_type(), ElementType::Hex8);
    EXPECT_EQ(lagrange->order(), 2);

    auto serendipity =
        basis_factory::create(BasisRequest{ElementType::Quad8, BasisType::Serendipity, 2});
    ASSERT_NE(serendipity, nullptr);
    EXPECT_EQ(serendipity->basis_type(), BasisType::Serendipity);

    EXPECT_THROW((void)basis_factory::create(
                     BasisRequest{ElementType::Pyramid5, BasisType::Lagrange, 1}),
                 BasisElementCompatibilityException);
    EXPECT_THROW((void)basis_factory::create(
                     BasisRequest{ElementType::Pyramid13, BasisType::Serendipity, 2}),
                 BasisElementCompatibilityException);
}
