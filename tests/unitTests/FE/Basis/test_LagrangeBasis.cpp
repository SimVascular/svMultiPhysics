/**
 * @file test_LagrangeBasis.cpp
 * @brief Unit tests for the reduced scalar Lagrange basis implementation.
 */

#include <gtest/gtest.h>

#include "FE/Basis/BasisExceptions.h"
#include "FE/Basis/BasisFactory.h"
#include "FE/Basis/LagrangeBasis.h"
#include "FE/Basis/NodeOrderingConventions.h"

#include <algorithm>
#include <array>
#include <span>
#include <tuple>
#include <vector>

using namespace svmp::FE;
using namespace svmp::FE::basis;

namespace {

using Point = math::Vector<double, 3>;

// Tolerance convention (set from measured residuals, not guessed). Identities that
// share a generator -- node coordinates, the lattice forward-image, span vs vector
// evaluation, P0 constants -- are exact up to round-off and use 1e-14. The
// nodal/partition identities here -- the Kronecker delta at the nodes and the
// value/gradient/Hessian partition sums (sum_i N_i = 1, sum_i grad N_i = 0,
// sum_i Hess N_i = 0) -- are exact at every order (at a node an off-diagonal shape
// function has an exactly-zero factor, and the partition sums differentiate the
// constant 1), so one tolerance covers them all instead of an order-dependent
// ladder. Measured residuals stay near round-off (~7e-15 through the orders here).
constexpr double kPartitionTol = double(1e-12);

struct CanonicalCase {
    BasisTopology topology;
    ElementType representative;  // linear element for sample-point lookup and labeling
    int order;
    std::size_t size;
    int dimension;
    std::vector<Point> points;
};

const std::vector<CanonicalCase>& canonical_cases() {
    static const std::vector<CanonicalCase> cases = {
        {BasisTopology::Line, ElementType::Line2, 3, 4u, 1,
         {{double(-0.35), double(0), double(0)}, {double(0.2), double(0), double(0)}}},
        {BasisTopology::Triangle, ElementType::Triangle3, 3, 10u, 2,
         {{double(0.15), double(0.2), double(0)}, {double(0.25), double(0.1), double(0)}}},
        {BasisTopology::Quadrilateral, ElementType::Quad4, 3, 16u, 2,
         {{double(0.2), double(-0.3), double(0)}, {double(-0.45), double(0.25), double(0)}}},
        {BasisTopology::Tetrahedron, ElementType::Tetra4, 2, 10u, 3,
         {{double(0.12), double(0.18), double(0.16)}, {double(0.2), double(0.1), double(0.18)}}},
        {BasisTopology::Hexahedron, ElementType::Hex8, 2, 27u, 3,
         {{double(0.1), double(-0.2), double(0.3)}, {double(-0.35), double(0.25), double(-0.15)}}},
        {BasisTopology::Wedge, ElementType::Wedge6, 2, 18u, 3,
         {{double(0.18), double(0.22), double(-0.2)}, {double(0.12), double(0.16), double(0.1)}}},
    };
    return cases;
}

std::vector<Point> sample_points_for(ElementType representative) {
    for (const auto& c : canonical_cases()) {
        if (c.representative == representative) {
            return c.points;
        }
    }
    return {};
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
            EXPECT_NEAR(values[i], i == node ? double(1) : double(0), tol)
                << "node=" << node << " basis=" << i;
        }
    }
}

void expect_partition_gradient_hessian_sums(const LagrangeBasis& basis,
                                            const std::vector<Point>& points)
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

        EXPECT_NEAR(value_sum, double(1), kPartitionTol);
        for (int d = 0; d < basis.dimension(); ++d) {
            EXPECT_NEAR(gradient_sum[static_cast<std::size_t>(d)], double(0), kPartitionTol);
            for (int e = 0; e < basis.dimension(); ++e) {
                EXPECT_NEAR(hessian_sum(static_cast<std::size_t>(d),
                                        static_cast<std::size_t>(e)),
                            double(0),
                            kPartitionTol);
            }
        }
    }
}

void expect_span_sinks_match_vector_evaluation(const LagrangeBasis& basis,
                                               const Point& xi)
{
    std::vector<double> values;
    std::vector<Gradient> gradients;
    std::vector<Hessian> hessians;
    basis.evaluate_all(xi, values, gradients, hessians);

    std::vector<double> span_values(basis.size());
    std::vector<Gradient> span_gradients(basis.size());
    std::vector<Hessian> span_hessians(basis.size());
    basis.evaluate_values_to(xi, span_values);
    basis.evaluate_gradients_to(xi, span_gradients);
    basis.evaluate_hessians_to(xi, span_hessians);

    for (std::size_t i = 0; i < basis.size(); ++i) {
        EXPECT_NEAR(span_values[i], values[i], double(1e-14));
        for (std::size_t d = 0; d < 3u; ++d) {
            EXPECT_NEAR(span_gradients[i][d], gradients[i][d], double(1e-14));
            for (std::size_t e = 0; e < 3u; ++e) {
                EXPECT_NEAR(span_hessians[i](d, e),
                            hessians[i](d, e),
                            double(1e-14));
            }
        }
    }
}

void expect_nodes_close(const std::vector<Point>& lhs,
                        const std::vector<Point>& rhs,
                        double tol)
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
                              double tol)
{
    ASSERT_EQ(lhs.size(), rhs.size());

    for (const auto& xi : points) {
        std::vector<double> lhs_values;
        std::vector<double> rhs_values;
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

double linear_function(const Point& p) {
    return double(2) + double(3) * p[0] - double(4) * p[1] + double(5) * p[2];
}

Gradient linear_gradient() {
    Gradient g = Gradient::Zero();
    g[0] = double(3);
    g[1] = double(-4);
    g[2] = double(5);
    return g;
}

double quadratic_function(const Point& p) {
    return double(1) + double(2) * p[0] - p[1] + double(0.5) * p[2] +
           p[0] * p[0] + double(0.75) * p[1] * p[1] - double(0.25) * p[2] * p[2] +
           double(0.2) * p[0] * p[1] - double(0.3) * p[0] * p[2] +
           double(0.4) * p[1] * p[2];
}

// Total degree three, so it lies in both the P3 simplex space and the Q3
// tensor-product space.
double cubic_function(const Point& p) {
    return quadratic_function(p) +
           double(0.1) * p[0] * p[0] * p[0] -
           double(0.2) * p[1] * p[1] * p[1] +
           double(0.3) * p[2] * p[2] * p[2] +
           double(0.15) * p[0] * p[0] * p[1] -
           double(0.12) * p[0] * p[2] * p[2] +
           double(0.08) * p[0] * p[1] * p[2];
}

template<typename Function>
double interpolate_value(const LagrangeBasis& basis,
                       const std::vector<double>& values,
                       Function&& nodal_function)
{
    double result = double(0);
    const auto& nodes = basis.nodes();
    for (std::size_t i = 0; i < values.size(); ++i) {
        result += values[i] * nodal_function(nodes[i]);
    }
    return result;
}

} // namespace

TEST(LagrangeBasis, CanonicalTopologiesHaveExpectedSizesAndDimensions) {
    for (const auto& c : canonical_cases()) {
        LagrangeBasis basis(c.topology, c.order);
        EXPECT_EQ(basis.basis_type(), BasisType::Lagrange);
        EXPECT_EQ(basis.topology(), c.topology);
        EXPECT_EQ(basis.order(), c.order);
        EXPECT_EQ(basis.size(), c.size);
        EXPECT_EQ(basis.dimension(), c.dimension);
    }
}

TEST(LagrangeBasis, CanonicalTopologiesAreNodalAndPartitionUnity) {
    for (const auto& c : canonical_cases()) {
        LagrangeBasis basis(c.topology, c.order);
        expect_kronecker_at_nodes(basis, kPartitionTol);
        expect_partition_gradient_hessian_sums(basis, c.points);
    }
}

TEST(LagrangeBasis, SpanOutputSinksMatchVectorEvaluationAcrossTopologies) {
    for (const auto& c : canonical_cases()) {
        LagrangeBasis basis(c.topology, c.order);
        expect_span_sinks_match_vector_evaluation(basis, c.points.front());
    }
}

// A named quadratic alias is a fixed-order shorthand for the same topology at
// order 2: it builds the identical basis as the BasisTopology overload, reports
// that topology, and round-trips faithfully through named_element_for() (Hex27
// stays Hex27 rather than collapsing to a canonical linear type).
TEST(LagrangeBasis, CompleteAliasesMatchTopologyConstruction) {
    const std::vector<std::tuple<ElementType, BasisTopology, ElementType>> aliases = {
        {ElementType::Line3, BasisTopology::Line, ElementType::Line2},
        {ElementType::Triangle6, BasisTopology::Triangle, ElementType::Triangle3},
        {ElementType::Quad9, BasisTopology::Quadrilateral, ElementType::Quad4},
        {ElementType::Tetra10, BasisTopology::Tetrahedron, ElementType::Tetra4},
        {ElementType::Hex27, BasisTopology::Hexahedron, ElementType::Hex8},
        {ElementType::Wedge18, BasisTopology::Wedge, ElementType::Wedge6},
    };

    for (const auto& [alias, topo, representative] : aliases) {
        LagrangeBasis alias_basis(alias, 2);
        LagrangeBasis topo_basis(topo, 2);
        const auto generated = ReferenceNodeLayout::get_lagrange_node_coords(representative, 2);

        EXPECT_EQ(alias_basis.topology(), topo);
        EXPECT_EQ(named_element_for(alias_basis.topology(), alias_basis.order(), alias_basis.basis_type()),
                  alias);
        EXPECT_EQ(alias_basis.order(), 2);
        expect_nodes_close(alias_basis.nodes(), generated, double(1e-14));
        expect_nodes_close(alias_basis.nodes(), topo_basis.nodes(), double(1e-14));
        expect_evaluations_match(alias_basis,
                                 topo_basis,
                                 sample_points_for(representative),
                                 double(1e-12));
    }
}

// The arbitrary order a named alias rejects is exactly what the BasisTopology
// overload is for: a node-count-named element cannot carry a conflicting order,
// and an order with no named layout maps to named_element_for() == Unknown.
TEST(LagrangeBasis, ArbitraryOrderRequiresTopologyNotNamedAlias) {
    EXPECT_THROW((void)LagrangeBasis(ElementType::Hex27, 3), BasisConfigurationException);

    const LagrangeBasis basis(BasisTopology::Hexahedron, 5);
    EXPECT_EQ(basis.topology(), BasisTopology::Hexahedron);
    EXPECT_EQ(basis.order(), 5);
    EXPECT_EQ(basis.size(), 216u);  // (5 + 1)^3
    EXPECT_EQ(named_element_for(basis.topology(), basis.order(), basis.basis_type()),
              ElementType::Unknown);  // no named order-5 hex
}

// named_element_for() is the inverse of (topology, order): a named layout at
// orders 0-2, Unknown where no named element exists (order 0 on a volume
// topology, or any order >= 3). topology() + order() remain the authoritative
// identity; callers recover a named ElementType through this free helper.
TEST(LagrangeBasis, NamedElementForReflectsTopologyAndOrder) {
    EXPECT_EQ(named_element_for(BasisTopology::Point, 0, BasisType::Lagrange), ElementType::Point1);
    EXPECT_EQ(named_element_for(BasisTopology::Hexahedron, 1, BasisType::Lagrange), ElementType::Hex8);
    EXPECT_EQ(named_element_for(BasisTopology::Hexahedron, 2, BasisType::Lagrange), ElementType::Hex27);
    EXPECT_EQ(named_element_for(BasisTopology::Hexahedron, 3, BasisType::Lagrange), ElementType::Unknown);
    EXPECT_EQ(named_element_for(BasisTopology::Tetrahedron, 1, BasisType::Lagrange), ElementType::Tetra4);
    EXPECT_EQ(named_element_for(BasisTopology::Tetrahedron, 2, BasisType::Lagrange), ElementType::Tetra10);
    EXPECT_EQ(named_element_for(BasisTopology::Quadrilateral, 2, BasisType::Lagrange), ElementType::Quad9);
    // An order-0 P0 basis on a volume topology has no named element.
    EXPECT_EQ(named_element_for(BasisTopology::Hexahedron, 0, BasisType::Lagrange), ElementType::Unknown);
}

// The single-argument named overload infers the layout's baked-in order, so it
// builds the same basis as passing that order explicitly, and still rejects the
// non-Lagrange (serendipity/pyramid) layouts the two-argument overload rejects.
TEST(LagrangeBasis, SingleArgumentNamedOverloadInfersBakedOrder) {
    const std::vector<std::pair<ElementType, int>> named = {
        {ElementType::Point1, 0},
        {ElementType::Line2, 1},     {ElementType::Line3, 2},
        {ElementType::Triangle3, 1}, {ElementType::Triangle6, 2},
        {ElementType::Quad4, 1},     {ElementType::Quad9, 2},
        {ElementType::Tetra4, 1},    {ElementType::Tetra10, 2},
        {ElementType::Hex8, 1},      {ElementType::Hex27, 2},
        {ElementType::Wedge6, 1},    {ElementType::Wedge18, 2},
    };
    for (const auto& [type, baked_order] : named) {
        const LagrangeBasis inferred(type);
        const LagrangeBasis explicit_order(type, baked_order);
        EXPECT_EQ(inferred.order(), baked_order) << "type=" << static_cast<int>(type);
        EXPECT_EQ(inferred.topology(), explicit_order.topology()) << "type=" << static_cast<int>(type);
        EXPECT_EQ(inferred.size(), explicit_order.size()) << "type=" << static_cast<int>(type);
    }

    EXPECT_THROW((void)LagrangeBasis(ElementType::Quad8), BasisElementCompatibilityException);
    EXPECT_THROW((void)LagrangeBasis(ElementType::Pyramid5), BasisElementCompatibilityException);
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
            const auto public_node = ReferenceNodeLayout::node_coord_at(alias, i);
            EXPECT_NEAR(public_node[0], generated[i][0], double(1e-14)) << "node=" << i;
            EXPECT_NEAR(public_node[1], generated[i][1], double(1e-14)) << "node=" << i;
            EXPECT_NEAR(public_node[2], generated[i][2], double(1e-14)) << "node=" << i;
        }
    }
}

// The lattice emitted with each node must be the exact forward image of the
// coordinate: tensor axes invert through line_coord_pm_one, simplex axes through
// the [0, 1] equispaced map, and the wedge combines the two. This pins the
// integer-lattice contract that replaced the floating-point round-trip, so a
// generator that emitted a coordinate and a mismatched index would fail here.
TEST(LagrangeBasis, LatticeIsExactForwardImageOfCoordinates) {
    constexpr double kTol = double(1e-14);

    const std::vector<std::tuple<ElementType, int, int>> tensor_cases = {
        {ElementType::Line2, 1, 1}, {ElementType::Line2, 4, 1},
        {ElementType::Quad4, 1, 2}, {ElementType::Quad4, 4, 2},
        {ElementType::Hex8, 1, 3},  {ElementType::Hex8, 3, 3},
    };
    for (const auto& [type, order, dim] : tensor_cases) {
        const auto layout = ReferenceNodeLayout::get_lagrange_lattice(type, order);
        ASSERT_EQ(layout.coords.size(), layout.lattice.size())
            << "type=" << static_cast<int>(type);
        for (std::size_t n = 0; n < layout.coords.size(); ++n) {
            for (int d = 0; d < dim; ++d) {
                const auto sd = static_cast<std::size_t>(d);
                EXPECT_NEAR(layout.coords[n][sd],
                            line_coord_pm_one(layout.lattice[n][sd], order), kTol)
                    << "type=" << static_cast<int>(type) << " node=" << n << " axis=" << d;
            }
        }
    }

    const std::vector<std::tuple<ElementType, int, int>> simplex_cases = {
        {ElementType::Triangle3, 1, 2}, {ElementType::Triangle3, 4, 2},
        {ElementType::Tetra4, 1, 3},    {ElementType::Tetra4, 4, 3},
    };
    for (const auto& [type, order, dim] : simplex_cases) {
        const auto layout = ReferenceNodeLayout::get_lagrange_lattice(type, order);
        ASSERT_EQ(layout.coords.size(), layout.lattice.size())
            << "type=" << static_cast<int>(type);
        for (std::size_t n = 0; n < layout.coords.size(); ++n) {
            for (int d = 0; d < dim; ++d) {
                const auto sd = static_cast<std::size_t>(d);
                EXPECT_NEAR(layout.coords[n][sd],
                            static_cast<double>(layout.lattice[n][sd]) /
                                static_cast<double>(order),
                            kTol)
                    << "type=" << static_cast<int>(type) << " node=" << n << " axis=" << d;
            }
        }
    }

    for (const int order : {1, 2, 3, 4}) {
        const auto layout =
            ReferenceNodeLayout::get_lagrange_lattice(ElementType::Wedge6, order);
        ASSERT_EQ(layout.coords.size(), layout.lattice.size()) << "wedge order=" << order;
        for (std::size_t n = 0; n < layout.coords.size(); ++n) {
            // (x, y) are triangle [0, 1] indices; z inverts through line_coord_pm_one.
            EXPECT_NEAR(layout.coords[n][0],
                        static_cast<double>(layout.lattice[n][0]) / static_cast<double>(order),
                        kTol) << "wedge order=" << order << " node=" << n;
            EXPECT_NEAR(layout.coords[n][1],
                        static_cast<double>(layout.lattice[n][1]) / static_cast<double>(order),
                        kTol) << "wedge order=" << order << " node=" << n;
            EXPECT_NEAR(layout.coords[n][2],
                        line_coord_pm_one(layout.lattice[n][2], order), kTol)
                << "wedge order=" << order << " node=" << n;
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

// The polynomial-reproduction and higher-order-lattice tests here validate
// VALUES and derivative invariants (gradient/Hessian sums). The authoritative
// finite-difference checks of gradient and Hessian *values* live in
// test_BasisHessians.cpp (BasisGradients/BasisHessians suites), covering the
// canonical Lagrange topologies and the serendipity families.
TEST(LagrangeBasis, LinearPolynomialReproductionAcrossLinearTopologies) {
    const std::vector<std::pair<ElementType, Point>> cases = {
        {ElementType::Line2, {double(-0.2), double(0), double(0)}},
        {ElementType::Triangle3, {double(0.2), double(0.3), double(0)}},
        {ElementType::Quad4, {double(0.25), double(-0.4), double(0)}},
        {ElementType::Tetra4, {double(0.1), double(0.2), double(0.3)}},
        {ElementType::Hex8, {double(0.15), double(-0.2), double(0.25)}},
        {ElementType::Wedge6, {double(0.2), double(0.15), double(-0.3)}},
    };
    const Gradient expected_gradient = linear_gradient();

    for (const auto& [type, point] : cases) {
        LagrangeBasis basis(type, 1);
        std::vector<double> values;
        std::vector<Gradient> gradients;
        basis.evaluate_values(point, values);
        basis.evaluate_gradients(point, gradients);

        const double interpolated =
            interpolate_value(basis, values, linear_function);
        EXPECT_NEAR(interpolated, linear_function(point), double(1e-12));

        Gradient interpolated_gradient = Gradient::Zero();
        for (std::size_t i = 0; i < gradients.size(); ++i) {
            const double nodal_value = linear_function(basis.nodes()[i]);
            for (int d = 0; d < basis.dimension(); ++d) {
                interpolated_gradient[static_cast<std::size_t>(d)] +=
                    nodal_value * gradients[i][static_cast<std::size_t>(d)];
            }
        }
        for (int d = 0; d < basis.dimension(); ++d) {
            EXPECT_NEAR(interpolated_gradient[static_cast<std::size_t>(d)],
                        expected_gradient[static_cast<std::size_t>(d)],
                        double(1e-12));
        }
    }
}

TEST(LagrangeBasis, QuadraticPolynomialReproductionAcrossQuadraticAliases) {
    const std::vector<std::pair<ElementType, Point>> cases = {
        {ElementType::Line3, {double(-0.2), double(0), double(0)}},
        {ElementType::Triangle6, {double(0.2), double(0.3), double(0)}},
        {ElementType::Quad9, {double(0.25), double(-0.4), double(0)}},
        {ElementType::Tetra10, {double(0.1), double(0.2), double(0.3)}},
        {ElementType::Hex27, {double(0.15), double(-0.2), double(0.25)}},
        {ElementType::Wedge18, {double(0.2), double(0.15), double(-0.3)}},
    };

    for (const auto& [type, point] : cases) {
        LagrangeBasis basis(type, 2);
        std::vector<double> values;
        basis.evaluate_values(point, values);

        const double interpolated =
            interpolate_value(basis, values, quadratic_function);
        EXPECT_NEAR(interpolated, quadratic_function(point), double(5e-12))
            << "element=" << static_cast<int>(type);
    }
}

// Tetra order >= 3 activates the face-interior node loops, tetra order >= 4
// activates the volume-interior lattice, and hex order >= 3 activates the six
// orientation-specific face traversals in NodeOrderingConventions. None of
// those generation paths run at the orders covered elsewhere; the Kronecker
// test is what validates the node lattice together with the integer
// lattice-index mapping the basis builds from it (a duplicated or missing
// node makes the basis non-nodal here).
TEST(LagrangeBasis, HigherOrderLatticesAreNodalAndPartitionUnity) {
    const struct Case {
        BasisTopology topology;
        int order;
        std::size_t size;
        std::vector<Point> points;
    } cases[] = {
        {BasisTopology::Tetrahedron, 3, 20u,
         {{double(0.12), double(0.18), double(0.16)}, {double(0.3), double(0.2), double(0.25)}}},
        {BasisTopology::Tetrahedron, 4, 35u,
         {{double(0.12), double(0.18), double(0.16)}, {double(0.2), double(0.1), double(0.18)}}},
        {BasisTopology::Hexahedron, 3, 64u,
         {{double(0.1), double(-0.2), double(0.3)}, {double(-0.35), double(0.25), double(-0.15)}}},
    };

    for (const auto& c : cases) {
        LagrangeBasis basis(c.topology, c.order);
        EXPECT_EQ(basis.size(), c.size);
        expect_kronecker_at_nodes(basis, kPartitionTol);
        expect_partition_gradient_hessian_sums(basis, c.points);
    }
}

// The Kronecker test above proves the order-3 hex lattice is nodal, but a
// permuted-yet-consistent face ordering would also pass it. This pins the
// load-bearing external contract of the order>=3 face-interior emission: the
// six face-interior blocks appear in VTK face order (-X, +X, -Y, +Y, -Z, +Z)
// and lie on the correct face. (The within-face traversal is an internal
// convention and is not separately pinned here.)
TEST(LagrangeBasis, HigherOrderHexFaceInteriorFollowsVtkFaceOrder) {
    // Order-3 hex (64 nodes): 8 vertices + 24 edge nodes + 24 face-interior
    // (6 faces x (order-1)^2 = 4 each) + 8 volume. Face-interior block at [32, 56).
    const auto nodes = ReferenceNodeLayout::get_lagrange_node_coords(ElementType::Hex8, 3);
    ASSERT_EQ(nodes.size(), 64u);

    struct FaceBlock {
        std::size_t axis;  // constant axis: 0=x, 1=y, 2=z
        double value;        // constant coordinate on the face
    };
    const FaceBlock blocks[] = {
        {0u, double(-1)},  // -X
        {0u, double(1)},   // +X
        {1u, double(-1)},  // -Y
        {1u, double(1)},   // +Y
        {2u, double(-1)},  // -Z
        {2u, double(1)},   // +Z
    };

    constexpr std::size_t kFaceStart = 32u;
    constexpr std::size_t kPerFace = 4u;  // (order-1)^2 at order 3
    for (std::size_t f = 0; f < 6u; ++f) {
        for (std::size_t m = 0; m < kPerFace; ++m) {
            const auto& node = nodes[kFaceStart + f * kPerFace + m];
            EXPECT_NEAR(node[blocks[f].axis], blocks[f].value, double(1e-14))
                << "face block " << f << ", node " << m;
        }
    }
}

TEST(LagrangeBasis, CubicPolynomialReproductionAtOrderThree) {
    const std::vector<std::pair<BasisTopology, Point>> cases = {
        {BasisTopology::Tetrahedron, {double(0.15), double(0.2), double(0.25)}},
        {BasisTopology::Hexahedron, {double(0.15), double(-0.2), double(0.25)}},
    };

    for (const auto& [topo, point] : cases) {
        LagrangeBasis basis(topo, 3);
        std::vector<double> values;
        basis.evaluate_values(point, values);

        const double interpolated = interpolate_value(basis, values, cubic_function);
        EXPECT_NEAR(interpolated, cubic_function(point), double(1e-10))
            << "topology=" << static_cast<int>(topo);
    }
}

TEST(LagrangeBasis, PointTopologyEvaluatesConstantUnity) {
    LagrangeBasis basis(ElementType::Point1, 0);

    EXPECT_EQ(named_element_for(basis.topology(), basis.order(), basis.basis_type()),
              ElementType::Point1);
    EXPECT_EQ(basis.size(), 1u);
    EXPECT_EQ(basis.dimension(), 0);
    ASSERT_EQ(basis.nodes().size(), 1u);

    const Point xi{double(0.3), double(-0.4), double(0.1)};
    std::vector<double> values;
    std::vector<Gradient> gradients;
    std::vector<Hessian> hessians;
    basis.evaluate_all(xi, values, gradients, hessians);

    ASSERT_EQ(values.size(), 1u);
    EXPECT_EQ(values[0], double(1));
    for (std::size_t d = 0; d < 3u; ++d) {
        EXPECT_EQ(gradients[0][d], double(0));
        for (std::size_t e = 0; e < 3u; ++e) {
            EXPECT_EQ(hessians[0](d, e), double(0));
        }
    }

    double span_value = double(-1);
    Gradient span_gradient;
    span_gradient[0] = span_gradient[1] = span_gradient[2] = double(-1);
    Hessian span_hessian;
    for (std::size_t d = 0; d < 3u; ++d) {
        for (std::size_t e = 0; e < 3u; ++e) {
            span_hessian(d, e) = double(-1);
        }
    }
    basis.evaluate_values_to(xi, std::span<double>(&span_value, 1u));
    basis.evaluate_gradients_to(xi, std::span<Gradient>(&span_gradient, 1u));
    basis.evaluate_hessians_to(xi, std::span<Hessian>(&span_hessian, 1u));
    EXPECT_EQ(span_value, double(1));
    for (std::size_t d = 0; d < 3u; ++d) {
        EXPECT_EQ(span_gradient[d], double(0));
    }
    for (std::size_t d = 0; d < 3u; ++d) {
        for (std::size_t e = 0; e < 3u; ++e) {
            EXPECT_EQ(span_hessian(d, e), double(0));
        }
    }
}

// P0 bases back piecewise-constant fields (e.g. pressure in mixed elements);
// the order-zero branches in node generation and the simplex/tensor/wedge
// evaluators have no other coverage.
TEST(LagrangeBasis, OrderZeroBasesAreConstantUnity) {
    const std::array<std::pair<BasisTopology, ElementType>, 6> cases = {{
        {BasisTopology::Line, ElementType::Line2},
        {BasisTopology::Triangle, ElementType::Triangle3},
        {BasisTopology::Quadrilateral, ElementType::Quad4},
        {BasisTopology::Tetrahedron, ElementType::Tetra4},
        {BasisTopology::Hexahedron, ElementType::Hex8},
        {BasisTopology::Wedge, ElementType::Wedge6},
    }};

    for (const auto& [topo, representative] : cases) {
        LagrangeBasis basis(topo, 0);
        EXPECT_EQ(basis.order(), 0) << "topology=" << static_cast<int>(topo);
        EXPECT_EQ(basis.size(), 1u) << "topology=" << static_cast<int>(topo);

        for (const auto& xi : sample_points_for(representative)) {
            std::vector<double> values;
            std::vector<Gradient> gradients;
            std::vector<Hessian> hessians;
            basis.evaluate_all(xi, values, gradients, hessians);

            ASSERT_EQ(values.size(), 1u);
            EXPECT_NEAR(values[0], double(1), double(1e-14))
                << "topology=" << static_cast<int>(topo);
            for (std::size_t d = 0; d < 3u; ++d) {
                EXPECT_NEAR(gradients[0][d], double(0), double(1e-14));
                for (std::size_t e = 0; e < 3u; ++e) {
                    EXPECT_NEAR(hessians[0](d, e), double(0), double(1e-14));
                }
            }
        }
    }
}

// Pins the default basis selection for every supported element type. The
// solver adapter (nn.cpp) translates solver element names to ElementType and
// delegates the family/order choice to default_basis_request; a silent change
// here would change the discretization of every simulation using that element.
TEST(BasisFactoryDefaults, SelectionsArePinnedForAllSupportedElements) {
    struct Expected {
        ElementType type;
        BasisType family;
        int order;
        std::size_t size;
    };
    const std::vector<Expected> cases = {
        {ElementType::Point1,    BasisType::Lagrange,    0, 1u},
        {ElementType::Line2,     BasisType::Lagrange,    1, 2u},
        {ElementType::Line3,     BasisType::Lagrange,    2, 3u},
        {ElementType::Triangle3, BasisType::Lagrange,    1, 3u},
        {ElementType::Triangle6, BasisType::Lagrange,    2, 6u},
        {ElementType::Quad4,     BasisType::Lagrange,    1, 4u},
        {ElementType::Quad8,     BasisType::Serendipity, 2, 8u},
        {ElementType::Quad9,     BasisType::Lagrange,    2, 9u},
        {ElementType::Tetra4,    BasisType::Lagrange,    1, 4u},
        {ElementType::Tetra10,   BasisType::Lagrange,    2, 10u},
        {ElementType::Hex8,      BasisType::Lagrange,    1, 8u},
        {ElementType::Hex20,     BasisType::Serendipity, 2, 20u},
        {ElementType::Hex27,     BasisType::Lagrange,    2, 27u},
        {ElementType::Wedge6,    BasisType::Lagrange,    1, 6u},
        {ElementType::Wedge15,   BasisType::Serendipity, 2, 15u},
        {ElementType::Wedge18,   BasisType::Lagrange,    2, 18u},
    };

    for (const auto& expected : cases) {
        const auto request = basis_factory::default_basis_request(expected.type);
        EXPECT_EQ(request.element_type, expected.type)
            << "element=" << static_cast<int>(expected.type);
        EXPECT_EQ(request.basis_type, expected.family)
            << "element=" << static_cast<int>(expected.type);
        ASSERT_TRUE(request.order.has_value())
            << "element=" << static_cast<int>(expected.type);
        EXPECT_EQ(*request.order, expected.order)
            << "element=" << static_cast<int>(expected.type);

        auto basis = basis_factory::create_default_for(expected.type);
        ASSERT_NE(basis, nullptr);
        EXPECT_EQ(basis->basis_type(), expected.family)
            << "element=" << static_cast<int>(expected.type);
        EXPECT_EQ(basis->order(), expected.order)
            << "element=" << static_cast<int>(expected.type);
        EXPECT_EQ(basis->size(), expected.size)
            << "element=" << static_cast<int>(expected.type);
    }
}

TEST(BasisFactoryDefaults, RejectsElementsWithoutDefaultBasis) {
    EXPECT_THROW((void)basis_factory::default_basis_request(ElementType::Pyramid5),
                 BasisElementCompatibilityException);
    EXPECT_THROW((void)basis_factory::default_basis_request(ElementType::Pyramid13),
                 BasisElementCompatibilityException);
    EXPECT_THROW((void)basis_factory::create_default_for(ElementType::Unknown),
                 BasisElementCompatibilityException);
}

// The factory-creation paths below were one 90-line test; they are split by the
// path under test so an early ASSERT_NE in one cannot mask the others, and each
// name says what it covers.
TEST(LagrangeBasis, FactoryCreatesNamedLagrangeBasis) {
    auto lagrange =
        basis_factory::create(BasisRequest{ElementType::Hex27, BasisType::Lagrange, 2});
    ASSERT_NE(lagrange, nullptr);
    EXPECT_EQ(lagrange->basis_type(), BasisType::Lagrange);
    EXPECT_EQ(lagrange->topology(), BasisTopology::Hexahedron);
    EXPECT_EQ(named_element_for(lagrange->topology(), lagrange->order(), lagrange->basis_type()),
              ElementType::Hex27);
    EXPECT_EQ(lagrange->order(), 2);

    // The factory inherits the named-element order validation: Hex27 is order 2.
    EXPECT_THROW((void)basis_factory::create(
                     BasisRequest{ElementType::Hex27, BasisType::Lagrange, 1}),
                 BasisConfigurationException);
}

TEST(LagrangeBasis, FactoryCreatesArbitraryOrderLagrangeBasis) {
    BasisRequest arbitrary_lagrange;
    arbitrary_lagrange.basis_type = BasisType::Lagrange;
    arbitrary_lagrange.order = 5;
    arbitrary_lagrange.topology = BasisTopology::Hexahedron;
    auto high_order_lagrange = basis_factory::create(arbitrary_lagrange);
    ASSERT_NE(high_order_lagrange, nullptr);
    EXPECT_EQ(high_order_lagrange->basis_type(), BasisType::Lagrange);
    EXPECT_EQ(high_order_lagrange->topology(), BasisTopology::Hexahedron);
    EXPECT_EQ(named_element_for(high_order_lagrange->topology(), high_order_lagrange->order(),
                               high_order_lagrange->basis_type()),
              ElementType::Unknown);
    EXPECT_EQ(high_order_lagrange->order(), 5);
    EXPECT_EQ(high_order_lagrange->size(), 216u);
}

TEST(LagrangeBasis, FactoryCreatesNamedSerendipityBasis) {
    auto serendipity =
        basis_factory::create(BasisRequest{ElementType::Quad8, BasisType::Serendipity, 2});
    ASSERT_NE(serendipity, nullptr);
    EXPECT_EQ(serendipity->basis_type(), BasisType::Serendipity);
}

TEST(LagrangeBasis, FactoryCreatesArbitraryOrderQuadSerendipityBasis) {
    BasisRequest arbitrary_quad_serendipity;
    arbitrary_quad_serendipity.basis_type = BasisType::Serendipity;
    arbitrary_quad_serendipity.order = 4;
    arbitrary_quad_serendipity.topology = BasisTopology::Quadrilateral;
    auto high_order_quad_serendipity = basis_factory::create(arbitrary_quad_serendipity);
    ASSERT_NE(high_order_quad_serendipity, nullptr);
    EXPECT_EQ(high_order_quad_serendipity->basis_type(), BasisType::Serendipity);
    EXPECT_EQ(high_order_quad_serendipity->topology(), BasisTopology::Quadrilateral);
    EXPECT_EQ(named_element_for(high_order_quad_serendipity->topology(),
                               high_order_quad_serendipity->order(),
                               high_order_quad_serendipity->basis_type()),
              ElementType::Unknown);
    EXPECT_EQ(high_order_quad_serendipity->order(), 4);
    EXPECT_EQ(high_order_quad_serendipity->size(), 17u);
}

TEST(LagrangeBasis, FactoryCreatesArbitraryOrderHexSerendipityBasis) {
    BasisRequest arbitrary_hex_serendipity;
    arbitrary_hex_serendipity.basis_type = BasisType::Serendipity;
    arbitrary_hex_serendipity.order = 3;
    arbitrary_hex_serendipity.topology = BasisTopology::Hexahedron;
    auto high_order_hex_serendipity = basis_factory::create(arbitrary_hex_serendipity);
    ASSERT_NE(high_order_hex_serendipity, nullptr);
    EXPECT_EQ(high_order_hex_serendipity->basis_type(), BasisType::Serendipity);
    EXPECT_EQ(high_order_hex_serendipity->topology(), BasisTopology::Hexahedron);
    EXPECT_EQ(named_element_for(high_order_hex_serendipity->topology(),
                               high_order_hex_serendipity->order(),
                               high_order_hex_serendipity->basis_type()),
              ElementType::Unknown);
    EXPECT_EQ(high_order_hex_serendipity->order(), 3);
    EXPECT_EQ(high_order_hex_serendipity->size(), 32u);
}

TEST(LagrangeBasis, FactoryRejectsInvalidScalarRequests) {
    EXPECT_THROW((void)basis_factory::create(
                     BasisRequest{ElementType::Pyramid5, BasisType::Lagrange, 1}),
                 BasisElementCompatibilityException);
    EXPECT_THROW((void)basis_factory::create(
                     BasisRequest{ElementType::Pyramid13, BasisType::Serendipity, 2}),
                 BasisElementCompatibilityException);

    BasisRequest ambiguous;
    ambiguous.element_type = ElementType::Hex27;
    ambiguous.basis_type = BasisType::Lagrange;
    ambiguous.order = 2;
    ambiguous.topology = BasisTopology::Hexahedron;
    EXPECT_THROW((void)basis_factory::create(ambiguous), BasisConfigurationException);

    BasisRequest missing_target;
    missing_target.basis_type = BasisType::Lagrange;
    missing_target.order = 2;
    EXPECT_THROW((void)basis_factory::create(missing_target), BasisConfigurationException);

    BasisRequest unsupported_serendipity_topology;
    unsupported_serendipity_topology.basis_type = BasisType::Serendipity;
    unsupported_serendipity_topology.order = 2;
    unsupported_serendipity_topology.topology = BasisTopology::Wedge;
    EXPECT_THROW((void)basis_factory::create(unsupported_serendipity_topology),
                 BasisElementCompatibilityException);
}

// The shared 1D tensor-axis distribution (line_coord_pm_one) is Gauss-Lobatto-
// Legendre: endpoints exactly +/-1, symmetric about 0, strictly increasing, and
// coincident with the equispaced layout at the production orders 1 and 2.
TEST(LagrangeBasis, TensorAxisNodesAreGaussLobattoLegendre) {
    // Orders 1 and 2 coincide with equispaced (the production layouts).
    EXPECT_EQ(line_coord_pm_one(0, 1), double(-1));
    EXPECT_EQ(line_coord_pm_one(1, 1), double(1));
    EXPECT_EQ(line_coord_pm_one(0, 2), double(-1));
    EXPECT_EQ(line_coord_pm_one(1, 2), double(0));
    EXPECT_EQ(line_coord_pm_one(2, 2), double(1));

    for (int order = 1; order <= 12; ++order) {
        EXPECT_EQ(line_coord_pm_one(0, order), double(-1)) << "order=" << order;
        EXPECT_EQ(line_coord_pm_one(order, order), double(1)) << "order=" << order;

        double previous = line_coord_pm_one(0, order);
        for (int i = 1; i <= order; ++i) {
            const double x = line_coord_pm_one(i, order);
            EXPECT_GT(x, previous) << "order=" << order << " i=" << i;  // strictly increasing
            EXPECT_LE(x, double(1)) << "order=" << order << " i=" << i;
            EXPECT_GE(x, double(-1)) << "order=" << order << " i=" << i;
            EXPECT_NEAR(x, -line_coord_pm_one(order - i, order), double(1e-14))
                << "order=" << order << " i=" << i;  // symmetric about 0
            previous = x;
        }
    }

    // For order >= 3 GLL differs from equispaced: the interior nodes cluster toward
    // the endpoints. At order 4 the first interior node is -sqrt(3/7) ~ -0.6547,
    // past the equispaced position -0.5.
    EXPECT_LT(line_coord_pm_one(1, 4), double(-0.5) - double(1e-6));
}
