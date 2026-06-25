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
                               const math::Vector<double, 3>& xi,
                               double tolerance = double(1e-10))
{
    std::vector<double> values;
    std::vector<Gradient> gradients;
    basis.evaluate_values(xi, values);
    basis.evaluate_gradients(xi, gradients);

    double value_sum = double(0);
    Gradient gradient_sum = Gradient::Zero();
    for (std::size_t i = 0; i < values.size(); ++i) {
        value_sum += values[i];
        for (std::size_t component = 0; component < 3u; ++component) {
            gradient_sum[component] += gradients[i][component];
        }
    }

    EXPECT_NEAR(value_sum, double(1), tolerance);
    for (int component = 0; component < basis.dimension(); ++component) {
        EXPECT_NEAR(gradient_sum[static_cast<std::size_t>(component)],
                    double(0),
                    tolerance);
    }
}

void expect_nodal_delta(const SerendipityBasis& basis,
                        const std::vector<math::Vector<double, 3>>& nodes,
                        double tolerance)
{
    ASSERT_EQ(nodes.size(), basis.size());
    for (std::size_t node = 0; node < nodes.size(); ++node) {
        std::vector<double> values;
        basis.evaluate_values(nodes[node], values);
        ASSERT_EQ(values.size(), basis.size());
        for (std::size_t dof = 0; dof < values.size(); ++dof) {
            EXPECT_NEAR(values[dof], dof == node ? double(1) : double(0), tolerance)
                << "node=" << node << " dof=" << dof;
        }
    }
}

std::vector<math::Vector<double, 3>> reference_nodes(ElementType type,
                                                   std::size_t count)
{
    std::vector<math::Vector<double, 3>> nodes;
    nodes.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        nodes.push_back(ReferenceNodeLayout::node_coord_at(type, i));
    }
    return nodes;
}

template<typename Function>
double interpolate_nodal_function(const SerendipityBasis& basis,
                                const math::Vector<double, 3>& xi,
                                Function&& nodal_function)
{
    std::vector<double> values;
    basis.evaluate_values(xi, values);

    double result = double(0);
    const auto& nodes = basis.nodes();
    for (std::size_t i = 0; i < values.size(); ++i) {
        result += values[i] * nodal_function(nodes[i]);
    }
    return result;
}

// The _for_test helpers below intentionally re-derive the production monomial
// selection (superlinear-degree rule, exponent enumeration, and size formula)
// independently of SerendipityBasis, so the basis is checked against an external
// oracle rather than against its own code. If the production formula in
// SerendipityBasis.cpp is changed deliberately, update these copies to match; an
// accidental drift between the two is meant to surface here as a test failure.
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

double integer_power_for_test(double base, int exponent) {
    double result = double(1);
    for (int k = 0; k < exponent; ++k) {
        result *= base;
    }
    return result;
}

double monomial_value_for_test(const math::Vector<double, 3>& p,
                             const std::array<int, 2>& exponent) {
    return integer_power_for_test(p[0], exponent[0]) *
           integer_power_for_test(p[1], exponent[1]);
}

std::vector<double> quadrilateral_vandermonde_for_test(
    const std::vector<math::Vector<double, 3>>& nodes,
    const std::vector<std::array<int, 2>>& exponents)
{
    const std::size_t n = nodes.size();
    std::vector<double> vandermonde(n * n, double(0));
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            vandermonde[row * n + col] =
                monomial_value_for_test(nodes[row], exponents[col]);
        }
    }
    return vandermonde;
}

void expect_no_duplicate_nodes(const std::vector<math::Vector<double, 3>>& nodes,
                               double tolerance)
{
    for (std::size_t a = 0; a < nodes.size(); ++a) {
        for (std::size_t b = a + 1u; b < nodes.size(); ++b) {
            const double dx = std::abs(nodes[a][0] - nodes[b][0]);
            const double dy = std::abs(nodes[a][1] - nodes[b][1]);
            EXPECT_GT(std::max(dx, dy), tolerance)
                << "duplicate nodes " << a << " and " << b;
        }
    }
}

void expect_nodes_near(const std::vector<math::Vector<double, 3>>& actual,
                       const std::vector<math::Vector<double, 3>>& expected,
                       double tolerance)
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
double cubic_serendipity_function(const math::Vector<double, 3>& p) {
    const double x = p[0];
    const double y = p[1];
    return double(1) + double(2) * x - y + double(3) * x * y +
           x * x * x - double(2) * y * y * y +
           double(0.5) * x * x * x * y - double(0.25) * x * y * y * y;
}

double bilinear_function(const math::Vector<double, 3>& p) {
    return double(2) - double(3) * p[0] + double(4) * p[1] + double(0.5) * p[0] * p[1];
}

// --- 3D serendipity guard-test helpers (Hex20, Wedge15) --------------------
//
// Like the quadrilateral _for_test helpers above, these re-derive the monomial
// selection and reference-node placement from the mathematical definition of
// each serendipity space, independently of the production tables in
// SerendipityBasis.cpp

double monomial_value_3d_for_test(const math::Vector<double, 3>& p,
                                const std::array<int, 3>& exponent) {
    return integer_power_for_test(p[0], exponent[0]) *
           integer_power_for_test(p[1], exponent[1]) *
           integer_power_for_test(p[2], exponent[2]);
}

std::vector<double> vandermonde_3d_for_test(
    const std::vector<math::Vector<double, 3>>& nodes,
    const std::vector<std::array<int, 3>>& exponents) {
    const std::size_t n = nodes.size();
    std::vector<double> vandermonde(n * n, double(0));
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            vandermonde[row * n + col] =
                monomial_value_3d_for_test(nodes[row], exponents[col]);
        }
    }
    return vandermonde;
}

// Superlinear degree generalized to three axes (the quadrilateral rule extended
// to t). An exponent contributes only when it exceeds one.
int superlinear_degree_3d_for_test(int ax, int ay, int az) {
    return (ax > 1 ? ax : 0) + (ay > 1 ? ay : 0) + (az > 1 ? az : 0);
}

// Hex20 serendipity span: every (ax, ay, az) in {0,1,2}^3 with superlinear
// degree at most two.
std::vector<std::array<int, 3>> hex20_serendipity_exponents_for_test() {
    std::vector<std::array<int, 3>> exponents;
    for (int ax = 0; ax <= 2; ++ax) {
        for (int ay = 0; ay <= 2; ++ay) {
            for (int az = 0; az <= 2; ++az) {
                if (superlinear_degree_3d_for_test(ax, ay, az) <= 2) {
                    exponents.push_back({ax, ay, az});
                }
            }
        }
    }
    return exponents;
}

// Arbitrary-order hexahedral serendipity verification
std::vector<std::array<int, 3>> hex_serendipity_exponents_for_test(int order) {
    std::vector<std::array<int, 3>> exponents;
    for (int az = 0; az <= order; ++az) {
        for (int ay = 0; ay <= order; ++ay) {
            for (int ax = 0; ax <= order; ++ax) {
                if (superlinear_degree_3d_for_test(ax, ay, az) <= order) {
                    exponents.push_back({ax, ay, az});
                }
            }
        }
    }
    return exponents;
}

std::size_t quad_serendipity_interior_count_for_test(int order) {
    if (order < 4) {
        return 0u;
    }
    const auto m = static_cast<std::size_t>(order - 4);
    return (m + 1u) * (m + 2u) / 2u;
}

std::size_t hex_serendipity_volume_interior_count_for_test(int order) {
    if (order < 6) {
        return 0u;
    }
    const auto m = static_cast<std::size_t>(order - 6);
    return (m + 1u) * (m + 2u) * (m + 3u) / 6u;
}

// dim S_p from the node strata: 8 corners, 12 (p - 1) edge nodes, 6 q(p) face
// interiors, and the P_{p-6} volume residual.
std::size_t expected_hex_serendipity_size(int order) {
    const auto p = static_cast<std::size_t>(order);
    return 8u + 12u * (p - 1u) +
           6u * quad_serendipity_interior_count_for_test(order) +
           hex_serendipity_volume_interior_count_for_test(order);
}

// Wedge15 serendipity span: triangle monomials (ax, ay) with ax + ay <= 2,
// tensored with the through-axis. Linear triangle monomials (ax + ay <= 1) carry
// t-degree up to two; quadratic triangle monomials (ax + ay == 2) carry t-degree
// up to one.
std::vector<std::array<int, 3>> wedge15_serendipity_exponents_for_test() {
    std::vector<std::array<int, 3>> exponents;
    for (int ax = 0; ax <= 2; ++ax) {
        for (int ay = 0; ax + ay <= 2; ++ay) {
            const int triangle_degree = ax + ay;
            const int max_t = (triangle_degree <= 1) ? 2 : 1;
            for (int az = 0; az <= max_t; ++az) {
                exponents.push_back({ax, ay, az});
            }
        }
    }
    return exponents;
}

// Independent Hex20 reference layout: the eight cube corners followed by the
// twelve edge midpoints, in the corner/edge order the reference layout uses.
std::vector<math::Vector<double, 3>> hex20_reference_nodes_for_test() {
    std::vector<math::Vector<double, 3>> corners;
    corners.push_back({double(-1), double(-1), double(-1)});
    corners.push_back({double(1), double(-1), double(-1)});
    corners.push_back({double(1), double(1), double(-1)});
    corners.push_back({double(-1), double(1), double(-1)});
    corners.push_back({double(-1), double(-1), double(1)});
    corners.push_back({double(1), double(-1), double(1)});
    corners.push_back({double(1), double(1), double(1)});
    corners.push_back({double(-1), double(1), double(1)});

    const int edges[12][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0},
        {4, 5}, {5, 6}, {6, 7}, {7, 4},
        {0, 4}, {1, 5}, {2, 6}, {3, 7},
    };

    std::vector<math::Vector<double, 3>> nodes = corners;
    for (const auto& edge : edges) {
        const math::Vector<double, 3> midpoint =
            (corners[static_cast<std::size_t>(edge[0])] +
             corners[static_cast<std::size_t>(edge[1])]) * double(0.5);
        nodes.push_back(midpoint);
    }
    return nodes;
}

// Independent Wedge15 reference layout: the six prism corners followed by the
// nine edge midpoints, in reference-layout order.
std::vector<math::Vector<double, 3>> wedge15_reference_nodes_for_test() {
    std::vector<math::Vector<double, 3>> corners;
    corners.push_back({double(0), double(0), double(-1)});
    corners.push_back({double(1), double(0), double(-1)});
    corners.push_back({double(0), double(1), double(-1)});
    corners.push_back({double(0), double(0), double(1)});
    corners.push_back({double(1), double(0), double(1)});
    corners.push_back({double(0), double(1), double(1)});

    const int edges[9][2] = {
        {0, 1}, {1, 2}, {2, 0},
        {3, 4}, {4, 5}, {5, 3},
        {0, 3}, {1, 4}, {2, 5},
    };

    std::vector<math::Vector<double, 3>> nodes = corners;
    for (const auto& edge : edges) {
        const math::Vector<double, 3> midpoint =
            (corners[static_cast<std::size_t>(edge[0])] +
             corners[static_cast<std::size_t>(edge[1])]) * double(0.5);
        nodes.push_back(midpoint);
    }
    return nodes;
}

// Independent Quad8 reference layout: the four quad corners followed by the four
// edge midpoints, in the corner/edge order the reference layout uses (the VTK
// quad boundary traversal). Mirrors the Hex20/Wedge15 anchors above.
std::vector<math::Vector<double, 3>> quad8_reference_nodes_for_test() {
    std::vector<math::Vector<double, 3>> corners;
    corners.push_back({double(-1), double(-1), double(0)});
    corners.push_back({double(1), double(-1), double(0)});
    corners.push_back({double(1), double(1), double(0)});
    corners.push_back({double(-1), double(1), double(0)});

    const int edges[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

    std::vector<math::Vector<double, 3>> nodes = corners;
    for (const auto& edge : edges) {
        const math::Vector<double, 3> midpoint =
            (corners[static_cast<std::size_t>(edge[0])] +
             corners[static_cast<std::size_t>(edge[1])]) * double(0.5);
        nodes.push_back(midpoint);
    }
    return nodes;
}

// --- Conditioning oracles (Legendre Vandermonde + Lebesgue constant) ----------

double legendre_value_for_test(double x, int degree) {
    if (degree <= 0) {
        return double(1);
    }
    double p_km1 = double(1);
    double p_k = x;
    for (int k = 1; k < degree; ++k) {
        const double p_kp1 =
            ((double(2) * double(k) + double(1)) * x * p_k - double(k) * p_km1) /
            double(k + 1);
        p_km1 = p_k;
        p_k = p_kp1;
    }
    return p_k;
}

double legendre_mode_for_test(const math::Vector<double, 3>& p,
                              const std::array<int, 3>& mode) {
    return legendre_value_for_test(p[0], mode[0]) *
           legendre_value_for_test(p[1], mode[1]) *
           legendre_value_for_test(p[2], mode[2]);
}

double matrix_norm_inf_for_test(const std::vector<double>& matrix, std::size_t n) {
    double max_row = double(0);
    for (std::size_t row = 0; row < n; ++row) {
        double sum = double(0);
        for (std::size_t col = 0; col < n; ++col) {
            sum += std::abs(matrix[row * n + col]);
        }
        max_row = std::max(max_row, sum);
    }
    return max_row;
}

// Infinity-norm condition number of the Legendre generalized Vandermonde the
// production basis inverts, rebuilt from the basis nodes and the re-derived
// serendipity modes (an independent check that Source B is fixed).
double legendre_vandermonde_condition(const std::vector<math::Vector<double, 3>>& nodes,
                                      const std::vector<std::array<int, 3>>& modes) {
    const std::size_t n = nodes.size();
    std::vector<double> v(n * n, double(0));
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            v[row * n + col] = legendre_mode_for_test(nodes[row], modes[col]);
        }
    }
    const double norm_v = matrix_norm_inf_for_test(v, n);
    const auto inverse = math::invert_dense_matrix(v, n, "test Legendre Vandermonde");
    return norm_v * matrix_norm_inf_for_test(inverse, n);
}

std::vector<std::array<int, 3>> quad_serendipity_modes_3d_for_test(int order) {
    std::vector<std::array<int, 3>> modes;
    for (const auto& e : quad_serendipity_exponents_for_test(order)) {
        modes.push_back({e[0], e[1], 0});
    }
    return modes;
}

// Lebesgue constant of the nodal basis: the maximum over a dense reference-cell
// sample of sum_i |N_i(xi)|. Bounded and slowly growing for GLL nodes; it is the
// "are the shape functions good" metric (equispaced nodes make it blow up).
double serendipity_lebesgue_constant(const SerendipityBasis& basis, int samples) {
    const int dim = basis.dimension();
    const auto axis = [samples](int idx) {
        return double(-1) + double(2) * double(idx) / double(samples);
    };
    double max_sum = double(0);
    std::vector<double> values;
    for (int i = 0; i <= samples; ++i) {
        for (int j = 0; j <= samples; ++j) {
            const int kmax = (dim >= 3) ? samples : 0;
            for (int k = 0; k <= kmax; ++k) {
                const math::Vector<double, 3> xi{
                    axis(i), axis(j), dim >= 3 ? axis(k) : double(0)};
                basis.evaluate_values(xi, values);
                double sum = double(0);
                for (const double v : values) {
                    sum += std::abs(v);
                }
                max_sum = std::max(max_sum, sum);
            }
        }
    }
    return max_sum;
}

} // namespace

TEST(SerendipityBasis, Quad8IsNodalAndPartitionsUnity) {
    SerendipityBasis basis(ElementType::Quad8, 2);
    SerendipityBasis topology_quad_basis(BasisTopology::Quadrilateral, 2);

    EXPECT_EQ(basis.size(), 8u);
    // The named Quad8 and the arbitrary-order Quadrilateral path at order 2 now
    // share the single ReferenceNodeLayout serendipity generator, so this pins
    // that the named and topology overloads build the same object. The independent
    // node-coordinate oracle is Quad8ReferenceNodesMatchIndependentConstruction.
    expect_nodes_near(basis.nodes(), topology_quad_basis.nodes(), double(1e-14));
    expect_nodal_delta(basis, basis.nodes(), double(1e-10));
    expect_partition_of_unity(basis, {double(0.17), double(-0.31), double(0)});
}

// Quad8 takes its reference nodes from ReferenceNodeLayout -- the single public
// node-ordering source the solver adapter permutes against, the same source
// Hex20 and Wedge15 use.
TEST(SerendipityBasis, Quad8ReferenceNodesComeFromReferenceNodeLayout) {
    SerendipityBasis basis(ElementType::Quad8, 2);
    expect_nodes_near(basis.nodes(),
                      ReferenceNodeLayout::node_coords(ElementType::Quad8),
                      double(1e-14));
}

// Independent node-coordinate anchor for the Quad8 layout: the four corners
// followed by the four edge midpoints, breaking the loop where the basis and the
// reference layout are otherwise only checked against each other. Mirrors the
// Hex20/Wedge15 independent-construction anchors.
TEST(SerendipityBasis, Quad8ReferenceNodesMatchIndependentConstruction) {
    SerendipityBasis basis(ElementType::Quad8, 2);
    expect_nodes_near(basis.nodes(), quad8_reference_nodes_for_test(), double(1e-14));
}

TEST(SerendipityBasis, Hex20IsNodalAndPartitionsUnity) {
    SerendipityBasis basis(ElementType::Hex20, 2);

    EXPECT_EQ(basis.size(), 20u);
    expect_nodal_delta(basis,
                       reference_nodes(ElementType::Hex20, basis.size()),
                       double(1e-10));
    expect_partition_of_unity(basis, {double(0.2), double(-0.1), double(0.3)});
}

TEST(SerendipityBasis, Wedge15IsNodalAndPartitionsUnity) {
    SerendipityBasis basis(ElementType::Wedge15, 2);

    EXPECT_EQ(basis.size(), 15u);
    expect_nodal_delta(basis,
                       reference_nodes(ElementType::Wedge15, basis.size()),
                       double(1e-9));
    expect_partition_of_unity(basis, {double(0.2), double(0.3), double(0.1)});
}

TEST(SerendipityBasis, RejectsUnsupportedSerendipityAliases) {
    EXPECT_THROW(SerendipityBasis(ElementType::Quad9, 2), BasisElementCompatibilityException);
    EXPECT_THROW(SerendipityBasis(ElementType::Pyramid13, 2), BasisElementCompatibilityException);
    EXPECT_THROW(SerendipityBasis(ElementType::Pyramid14, 2), BasisElementCompatibilityException);
    EXPECT_THROW(SerendipityBasis(ElementType::Quad8, 3), BasisConfigurationException);
    EXPECT_THROW(SerendipityBasis(ElementType::Quad8, 1), BasisConfigurationException);
    // Quad4 is the linear Lagrange quad, not a named serendipity layout; arbitrary
    // quadrilateral serendipity is requested through BasisTopology::Quadrilateral.
    EXPECT_THROW(SerendipityBasis(ElementType::Quad4, 2), BasisElementCompatibilityException);
}

// Topology construction is the arbitrary-order entry point and exists only for
// the quadrilateral, the single serendipity family with a free order. Hex and
// wedge serendipity are fixed layouts requested through their named ElementType.
TEST(SerendipityBasis, TopologyConstructionSupportsQuadrilateralAndHexahedron) {
    EXPECT_NO_THROW((void)SerendipityBasis(BasisTopology::Quadrilateral, 3));
    EXPECT_NO_THROW((void)SerendipityBasis(BasisTopology::Hexahedron, 3));
    EXPECT_THROW(SerendipityBasis(BasisTopology::Wedge, 2),
                 BasisElementCompatibilityException);
    EXPECT_THROW(SerendipityBasis(BasisTopology::Triangle, 2),
                 BasisElementCompatibilityException);

    // Topology and named construction agree at the production order for both the
    // quadrilateral and the hexahedron.
    SerendipityBasis quad(BasisTopology::Quadrilateral, 2);
    EXPECT_EQ(quad.topology(), BasisTopology::Quadrilateral);
    EXPECT_EQ(quad.order(), 2);
    EXPECT_EQ(named_element_for(quad.topology(), quad.order(), quad.basis_type()),
              ElementType::Quad8);

    SerendipityBasis hex(BasisTopology::Hexahedron, 2);
    EXPECT_EQ(hex.topology(), BasisTopology::Hexahedron);
    EXPECT_EQ(hex.order(), 2);
    EXPECT_EQ(named_element_for(hex.topology(), hex.order(), hex.basis_type()),
              ElementType::Hex20);
}

TEST(SerendipityBasis, SingleArgumentNamedOverloadInfersFixedOrder) {
    const std::vector<std::pair<ElementType, int>> named = {
        {ElementType::Quad8, 2},
        {ElementType::Hex8, 1},
        {ElementType::Hex20, 2},
        {ElementType::Wedge15, 2},
    };
    for (const auto& [type, fixed_order] : named) {
        const SerendipityBasis inferred(type);
        const SerendipityBasis explicit_order(type, fixed_order);
        EXPECT_EQ(inferred.order(), fixed_order) << "type=" << static_cast<int>(type);
        EXPECT_EQ(inferred.topology(), explicit_order.topology()) << "type=" << static_cast<int>(type);
        EXPECT_EQ(inferred.size(), explicit_order.size()) << "type=" << static_cast<int>(type);
    }

    EXPECT_THROW((void)SerendipityBasis(ElementType::Quad4), BasisElementCompatibilityException);
    EXPECT_THROW((void)SerendipityBasis(ElementType::Tetra4), BasisElementCompatibilityException);
}

TEST(SerendipityBasis, QuadrilateralRejectsOrdersBelowOne) {
    // Serendipity bases require a positive polynomial order; orders <= 0 are
    // rejected rather than normalized up to the linear space.
    EXPECT_THROW(SerendipityBasis(BasisTopology::Quadrilateral, 0),
                 BasisConfigurationException);
    EXPECT_THROW(SerendipityBasis(BasisTopology::Quadrilateral, -1),
                 BasisConfigurationException);

    // Order 1 is the smallest valid quadrilateral serendipity (the bilinear Q1
    // space): four corner nodes and the nodal-interpolation property.
    SerendipityBasis basis(BasisTopology::Quadrilateral, 1);
    EXPECT_EQ(basis.order(), 1);
    EXPECT_EQ(basis.size(), 4u);
    expect_nodal_delta(basis, basis.nodes(), double(1e-12));
}

// Explicit quadrilateral-topology serendipity orders run the documented monomial
// selection, boundary plus triangular interior node placement, and runtime
// Vandermonde inversion. Order four is the first order with an interior residual
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
        SerendipityBasis basis(BasisTopology::Quadrilateral, c.order);
        EXPECT_EQ(basis.size(), c.size) << "order=" << c.order;
        EXPECT_EQ(basis.order(), c.order);
        EXPECT_EQ(basis.dimension(), 2);
        ASSERT_EQ(basis.nodes().size(), c.size);

        for (const auto& node : basis.nodes()) {
            EXPECT_LE(std::abs(node[0]), double(1));
            EXPECT_LE(std::abs(node[1]), double(1));
        }

        expect_nodal_delta(basis, basis.nodes(), double(1e-9));
        expect_partition_of_unity(basis, {double(0.17), double(-0.31), double(0)}, double(1e-9));
        expect_partition_of_unity(basis, {double(-0.45), double(0.25), double(0)}, double(1e-9));
    }
}

TEST(SerendipityBasis, QuadrilateralNodesFollowDocumentedConstructionThroughOrderTen) {
    constexpr double kTol = double(1e-14);

    for (int order = 1; order <= 10; ++order) {
        SerendipityBasis basis(BasisTopology::Quadrilateral, order);
        const auto& nodes = basis.nodes();
        const std::size_t expected_size = expected_quad_serendipity_size(order);
        const std::size_t boundary_count = static_cast<std::size_t>(4 * order);

        ASSERT_EQ(basis.size(), expected_size) << "order=" << order;
        ASSERT_EQ(nodes.size(), expected_size) << "order=" << order;
        EXPECT_EQ(quad_serendipity_exponents_for_test(order).size(),
                  expected_size) << "order=" << order;
        expect_no_duplicate_nodes(nodes, kTol);

        for (std::size_t i = 0; i < nodes.size(); ++i) {
            EXPECT_NEAR(nodes[i][2], double(0), kTol) << "order=" << order
                                                    << " node=" << i;
            EXPECT_LE(std::abs(nodes[i][0]), double(1)) << "order=" << order
                                                       << " node=" << i;
            EXPECT_LE(std::abs(nodes[i][1]), double(1)) << "order=" << order
                                                       << " node=" << i;

            const bool on_boundary =
                std::abs(std::abs(nodes[i][0]) - double(1)) <= kTol ||
                std::abs(std::abs(nodes[i][1]) - double(1)) <= kTol;
            if (i < boundary_count) {
                EXPECT_TRUE(on_boundary) << "order=" << order << " node=" << i;
            } else {
                EXPECT_FALSE(on_boundary) << "order=" << order << " node=" << i;
                EXPECT_LT(std::abs(nodes[i][0]), double(1)) << "order=" << order
                                                           << " node=" << i;
                EXPECT_LT(std::abs(nodes[i][1]), double(1)) << "order=" << order
                                                           << " node=" << i;
            }
        }

        std::size_t index = boundary_count;
        if (order >= 4) {
            // The interior staircase sits on Gauss-Lobatto-Legendre interior nodes:
            // row r at the (r+1)-th GLL node of order m+2, each row's columns at the
            // GLL interior of order row_count+1 (same line_coord_pm_one the basis
            // uses), re-derived here as an independent placement oracle.
            const int m = order - 4;
            for (int row = 0; row <= m; ++row) {
                const int row_count = m + 1 - row;
                const double expected_y = line_coord_pm_one(row + 1, m + 2);
                for (int col = 0; col < row_count; ++col) {
                    ASSERT_LT(index, nodes.size());
                    const double expected_x = line_coord_pm_one(col + 1, row_count + 1);
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
    SerendipityBasis basis(BasisTopology::Quadrilateral, 1);

    const std::vector<math::Vector<double, 3>> points = {
        {double(0.25), double(-0.4), double(0)},
        {double(-0.7), double(0.6), double(0)},
    };
    for (const auto& xi : points) {
        EXPECT_NEAR(interpolate_nodal_function(basis, xi, bilinear_function),
                    bilinear_function(xi),
                    double(1e-12));
    }
}

TEST(SerendipityBasis, QuadrilateralOrderThreeReproducesSerendipityCubics) {
    SerendipityBasis basis(BasisTopology::Quadrilateral, 3);

    const std::vector<math::Vector<double, 3>> points = {
        {double(0.25), double(-0.4), double(0)},
        {double(-0.7), double(0.6), double(0)},
    };
    for (const auto& xi : points) {
        EXPECT_NEAR(interpolate_nodal_function(basis, xi, cubic_serendipity_function),
                    cubic_serendipity_function(xi),
                    double(1e-11));
    }
}

TEST(SerendipityBasis, QuadrilateralOrdersReproduceEverySerendipityMonomial) {
    const std::vector<math::Vector<double, 3>> points = {
        {double(0.25), double(-0.4), double(0)},
        {double(-0.7), double(0.6), double(0)},
        {double(0.11), double(0.23), double(0)},
    };

    for (int order = 1; order <= 10; ++order) {
        SerendipityBasis basis(BasisTopology::Quadrilateral, order);
        const auto exponents = quad_serendipity_exponents_for_test(order);
        ASSERT_EQ(exponents.size(), basis.size()) << "order=" << order;

        // Uniformly tight across the whole range: GLL nodes and the Legendre modal
        // basis keep the reproduction accurate even at order 10 (the equispaced/
        // monomial construction needed 2e-8 here).
        const double tolerance = double(1e-10);
        for (const auto& exponent : exponents) {
            for (const auto& xi : points) {
                const double interpolated =
                    interpolate_nodal_function(
                        basis,
                        xi,
                        [&exponent](const math::Vector<double, 3>& node) {
                            return monomial_value_for_test(node, exponent);
                        });
                const double expected = monomial_value_for_test(xi, exponent);
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
        SerendipityBasis basis(BasisTopology::Quadrilateral, order);
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

// Hex8 serendipity is the order-1 instance of the generated hexahedral
// serendipity space (the eight multilinear monomials). It must still reproduce
// the trilinear Lagrange basis -- values, gradients, and Hessians -- which guards
// the generated order-1 coefficient table against the closed-form trilinear basis.
TEST(SerendipityBasis, TrilinearHexMatchesLagrangeHex8) {
    SerendipityBasis serendipity(ElementType::Hex8, 1);
    LagrangeBasis lagrange(ElementType::Hex8, 1);

    EXPECT_EQ(serendipity.size(), 8u);
    EXPECT_EQ(serendipity.dimension(), 3);
    expect_nodal_delta(serendipity,
                       reference_nodes(ElementType::Hex8, serendipity.size()),
                       double(1e-12));

    const std::vector<math::Vector<double, 3>> points = {
        {double(0.2), double(-0.1), double(0.3)},
        {double(-0.35), double(0.25), double(-0.15)},
    };
    for (const auto& xi : points) {
        std::vector<double> s_values;
        std::vector<double> l_values;
        std::vector<Gradient> s_gradients;
        std::vector<Gradient> l_gradients;
        std::vector<Hessian> s_hessians;
        std::vector<Hessian> l_hessians;
        serendipity.evaluate_all(xi, s_values, s_gradients, s_hessians);
        lagrange.evaluate_all(xi, l_values, l_gradients, l_hessians);

        ASSERT_EQ(s_values.size(), l_values.size());
        for (std::size_t i = 0; i < s_values.size(); ++i) {
            EXPECT_NEAR(s_values[i], l_values[i], double(1e-13));
            for (std::size_t d = 0; d < 3u; ++d) {
                EXPECT_NEAR(s_gradients[i][d], l_gradients[i][d], double(1e-13));
                for (std::size_t e = 0; e < 3u; ++e) {
                    EXPECT_NEAR(s_hessians[i](d, e), l_hessians[i](d, e), double(1e-13));
                }
            }
        }
    }
}

// The Hex20 and Wedge15 guard tests below pin the nodal coefficients away from
// the reference nodes. Each builds an independent Vandermonde from the
// re-derived monomial span and reference nodes, so a coefficient error that
// still vanishes at the nodes is caught.

TEST(SerendipityBasis, Hex20VandermondeHasFullRank) {
    SerendipityBasis basis(ElementType::Hex20, 2);
    const auto exponents = hex20_serendipity_exponents_for_test();
    const std::size_t n = basis.size();
    ASSERT_EQ(exponents.size(), n);
    const auto vandermonde = vandermonde_3d_for_test(basis.nodes(), exponents);
    ASSERT_EQ(vandermonde.size(), n * n);
    EXPECT_EQ(math::dense_matrix_rank(vandermonde, n, n), n);
}

TEST(SerendipityBasis, Wedge15VandermondeHasFullRank) {
    SerendipityBasis basis(ElementType::Wedge15, 2);
    const auto exponents = wedge15_serendipity_exponents_for_test();
    const std::size_t n = basis.size();
    ASSERT_EQ(exponents.size(), n);
    const auto vandermonde = vandermonde_3d_for_test(basis.nodes(), exponents);
    ASSERT_EQ(vandermonde.size(), n * n);
    EXPECT_EQ(math::dense_matrix_rank(vandermonde, n, n), n);
}

// V * C == I guard: independently invert the Vandermonde and confirm the basis
// evaluates to the same inverse-Vandermonde nodal functions, without reading the
// basis's internal coefficient table.
TEST(SerendipityBasis, Hex20MatchesIndependentlyInvertedVandermonde) {
    SerendipityBasis basis(ElementType::Hex20, 2);
    const auto exponents = hex20_serendipity_exponents_for_test();
    const std::size_t n = basis.size();
    ASSERT_EQ(exponents.size(), n);
    auto vandermonde = vandermonde_3d_for_test(basis.nodes(), exponents);
    const auto coefficients =
        math::invert_dense_matrix(std::move(vandermonde), n, "Hex20 test Vandermonde");

    const std::vector<math::Vector<double, 3>> points = {
        {double(0.2), double(-0.1), double(0.3)},
        {double(-0.35), double(0.25), double(-0.15)},
        {double(0.11), double(0.23), double(-0.42)},
    };
    for (const auto& xi : points) {
        std::vector<double> values;
        basis.evaluate_values(xi, values);
        ASSERT_EQ(values.size(), n);
        for (std::size_t i = 0; i < n; ++i) {
            double expected = double(0);
            for (std::size_t j = 0; j < n; ++j) {
                expected += coefficients[j * n + i] *
                            monomial_value_3d_for_test(xi, exponents[j]);
            }
            EXPECT_NEAR(values[i], expected, double(1e-10)) << "basis=" << i;
        }
    }
}

TEST(SerendipityBasis, Wedge15MatchesIndependentlyInvertedVandermonde) {
    SerendipityBasis basis(ElementType::Wedge15, 2);
    const auto exponents = wedge15_serendipity_exponents_for_test();
    const std::size_t n = basis.size();
    ASSERT_EQ(exponents.size(), n);
    auto vandermonde = vandermonde_3d_for_test(basis.nodes(), exponents);
    const auto coefficients =
        math::invert_dense_matrix(std::move(vandermonde), n, "Wedge15 test Vandermonde");

    const std::vector<math::Vector<double, 3>> points = {
        {double(0.2), double(0.3), double(0.1)},
        {double(0.25), double(0.25), double(-0.4)},
        {double(0.1), double(0.6), double(0.5)},
    };
    for (const auto& xi : points) {
        std::vector<double> values;
        basis.evaluate_values(xi, values);
        ASSERT_EQ(values.size(), n);
        for (std::size_t i = 0; i < n; ++i) {
            double expected = double(0);
            for (std::size_t j = 0; j < n; ++j) {
                expected += coefficients[j * n + i] *
                            monomial_value_3d_for_test(xi, exponents[j]);
            }
            EXPECT_NEAR(values[i], expected, double(1e-10)) << "basis=" << i;
        }
    }
}

// Non-nodal polynomial reproduction: the basis must reproduce every monomial in
// its span at interior points, not just interpolate at the nodes.
TEST(SerendipityBasis, Hex20ReproducesEverySerendipityMonomial) {
    SerendipityBasis basis(ElementType::Hex20, 2);
    const auto exponents = hex20_serendipity_exponents_for_test();
    ASSERT_EQ(exponents.size(), basis.size());

    const std::vector<math::Vector<double, 3>> points = {
        {double(0.2), double(-0.1), double(0.3)},
        {double(-0.35), double(0.25), double(-0.15)},
        {double(0.11), double(0.23), double(-0.42)},
    };
    for (const auto& exponent : exponents) {
        for (const auto& xi : points) {
            const double interpolated = interpolate_nodal_function(
                basis, xi,
                [&exponent](const math::Vector<double, 3>& node) {
                    return monomial_value_3d_for_test(node, exponent);
                });
            EXPECT_NEAR(interpolated, monomial_value_3d_for_test(xi, exponent),
                        double(1e-10))
                << "ax=" << exponent[0] << " ay=" << exponent[1]
                << " az=" << exponent[2];
        }
    }
}

TEST(SerendipityBasis, Wedge15ReproducesEverySerendipityMonomial) {
    SerendipityBasis basis(ElementType::Wedge15, 2);
    const auto exponents = wedge15_serendipity_exponents_for_test();
    ASSERT_EQ(exponents.size(), basis.size());

    const std::vector<math::Vector<double, 3>> points = {
        {double(0.2), double(0.3), double(0.1)},
        {double(0.25), double(0.25), double(-0.4)},
        {double(0.1), double(0.6), double(0.5)},
    };
    for (const auto& exponent : exponents) {
        for (const auto& xi : points) {
            const double interpolated = interpolate_nodal_function(
                basis, xi,
                [&exponent](const math::Vector<double, 3>& node) {
                    return monomial_value_3d_for_test(node, exponent);
                });
            EXPECT_NEAR(interpolated, monomial_value_3d_for_test(xi, exponent),
                        double(1e-10))
                << "ax=" << exponent[0] << " ay=" << exponent[1]
                << " az=" << exponent[2];
        }
    }
}

// Independent node-coordinate anchor: the reference nodes must be the cube/prism
// corners and edge midpoints, breaking the loop where the basis and its node
// table are otherwise only checked against each other.
TEST(SerendipityBasis, Hex20ReferenceNodesMatchIndependentConstruction) {
    SerendipityBasis basis(ElementType::Hex20, 2);
    expect_nodes_near(basis.nodes(), hex20_reference_nodes_for_test(), double(1e-14));
}

TEST(SerendipityBasis, Wedge15ReferenceNodesMatchIndependentConstruction) {
    SerendipityBasis basis(ElementType::Wedge15, 2);
    expect_nodes_near(basis.nodes(), wedge15_reference_nodes_for_test(), double(1e-14));
}

// --- Arbitrary-order hexahedral serendipity (BasisTopology::Hexahedron) -------

// dim S_p of the cube serendipity space for p = 1..6 (Hex8 = 8, Hex20 = 20),
// checked against both the re-derived monomial enumeration and the node-strata
// decomposition.
TEST(SerendipityBasis, HexahedralSerendipitySpaceHasExpectedDimensions) {
    const std::array<std::size_t, 7> expected = {0u, 8u, 20u, 32u, 50u, 74u, 105u};
    for (int order = 1; order <= 6; ++order) {
        const auto exponents = hex_serendipity_exponents_for_test(order);
        const auto p = static_cast<std::size_t>(order);
        EXPECT_EQ(exponents.size(), expected[p]) << "order=" << order;
        EXPECT_EQ(expected_hex_serendipity_size(order), expected[p]) << "order=" << order;
        for (const auto& e : exponents) {
            EXPECT_LE(superlinear_degree_3d_for_test(e[0], e[1], e[2]), order);
            for (int d = 0; d < 3; ++d) {
                EXPECT_GE(e[d], 0);
                EXPECT_LE(e[d], order);
            }
        }
    }

    // The order-2 hex serendipity span is exactly the Hex20 span (as a set).
    auto order_two = hex_serendipity_exponents_for_test(2);
    auto hex20 = hex20_serendipity_exponents_for_test();
    std::sort(order_two.begin(), order_two.end());
    std::sort(hex20.begin(), hex20.end());
    EXPECT_EQ(order_two, hex20);
}

// VTK conformance: the generated arbitrary-order layout reproduces the public
// Hex8 (order 1) and Hex20 (order 2) node ordering coordinate-for-coordinate.
TEST(SerendipityBasis, HexahedralTopologyNodesMatchPublicHex8AndHex20Layouts) {
    SerendipityBasis hex8(BasisTopology::Hexahedron, 1);
    EXPECT_EQ(hex8.size(), 8u);
    expect_nodes_near(hex8.nodes(),
                      ReferenceNodeLayout::node_coords(ElementType::Hex8),
                      double(1e-14));

    SerendipityBasis hex20(BasisTopology::Hexahedron, 2);
    EXPECT_EQ(hex20.size(), 20u);
    expect_nodes_near(hex20.nodes(),
                      ReferenceNodeLayout::node_coords(ElementType::Hex20),
                      double(1e-14));
}

TEST(SerendipityBasis, SkeletonMatchesCompleteLagrangePrefix) {
    for (int order = 1; order <= 8; ++order) {
        SerendipityBasis quad(BasisTopology::Quadrilateral, order);
        const auto quad_complete =
            ReferenceNodeLayout::get_lagrange_node_coords(ElementType::Quad4, order);
        const std::size_t quad_skeleton = static_cast<std::size_t>(4 * order);
        ASSERT_LE(quad_skeleton, quad.nodes().size()) << "quad order=" << order;
        ASSERT_LE(quad_skeleton, quad_complete.size()) << "quad order=" << order;
        for (std::size_t i = 0; i < quad_skeleton; ++i) {
            for (std::size_t d = 0; d < 3u; ++d) {
                EXPECT_EQ(quad.nodes()[i][d], quad_complete[i][d])
                    << "quad order=" << order << " node=" << i << " d=" << d;
            }
        }

        SerendipityBasis hex(BasisTopology::Hexahedron, order);
        const auto hex_complete =
            ReferenceNodeLayout::get_lagrange_node_coords(ElementType::Hex8, order);
        const std::size_t hex_skeleton =
            8u + 12u * static_cast<std::size_t>(order - 1);
        ASSERT_LE(hex_skeleton, hex.nodes().size()) << "hex order=" << order;
        ASSERT_LE(hex_skeleton, hex_complete.size()) << "hex order=" << order;
        for (std::size_t i = 0; i < hex_skeleton; ++i) {
            for (std::size_t d = 0; d < 3u; ++d) {
                EXPECT_EQ(hex.nodes()[i][d], hex_complete[i][d])
                    << "hex order=" << order << " node=" << i << " d=" << d;
            }
        }
    }
}

// The generated node set is unisolvent for the hex serendipity span at every
// supported order: the Vandermonde of the re-derived monomials at the generated
// nodes has full rank.
TEST(SerendipityBasis, HexahedralSerendipityVandermondeHasFullRankThroughOrderSix) {
    for (int order = 1; order <= 6; ++order) {
        SerendipityBasis basis(BasisTopology::Hexahedron, order);
        const auto exponents = hex_serendipity_exponents_for_test(order);
        const std::size_t n = basis.size();
        ASSERT_EQ(exponents.size(), n) << "order=" << order;
        const auto vandermonde = vandermonde_3d_for_test(basis.nodes(), exponents);
        ASSERT_EQ(vandermonde.size(), n * n) << "order=" << order;
        EXPECT_EQ(math::dense_matrix_rank(vandermonde, n, n), n) << "order=" << order;
    }
}

TEST(SerendipityBasis, HexahedralTopologyIsNodalAndPartitionsUnity) {
    const std::vector<math::Vector<double, 3>> points = {
        {double(0.2), double(-0.1), double(0.3)},
        {double(-0.35), double(0.25), double(-0.15)},
    };
    for (int order = 1; order <= 5; ++order) {
        SerendipityBasis basis(BasisTopology::Hexahedron, order);
        EXPECT_EQ(basis.dimension(), 3);
        EXPECT_EQ(basis.order(), order);
        EXPECT_EQ(basis.size(), expected_hex_serendipity_size(order)) << "order=" << order;
        ASSERT_EQ(basis.nodes().size(), basis.size());

        expect_nodal_delta(basis, basis.nodes(), double(1e-9));
        for (const auto& xi : points) {
            expect_partition_of_unity(basis, xi, double(1e-9));
        }
    }
}

// Non-nodal polynomial reproduction across orders: the basis reproduces every
// monomial in its span at interior points, pinning the production monomial space
// against the re-derived verification.
TEST(SerendipityBasis, HexahedralTopologyReproducesEverySerendipityMonomial) {
    const std::vector<math::Vector<double, 3>> points = {
        {double(0.2), double(-0.1), double(0.3)},
        {double(-0.35), double(0.25), double(-0.15)},
        {double(0.11), double(0.23), double(-0.42)},
    };
    for (int order = 1; order <= 5; ++order) {
        SerendipityBasis basis(BasisTopology::Hexahedron, order);
        const auto exponents = hex_serendipity_exponents_for_test(order);
        ASSERT_EQ(exponents.size(), basis.size()) << "order=" << order;
        for (const auto& exponent : exponents) {
            for (const auto& xi : points) {
                const double interpolated = interpolate_nodal_function(
                    basis, xi,
                    [&exponent](const math::Vector<double, 3>& node) {
                        return monomial_value_3d_for_test(node, exponent);
                    });
                EXPECT_NEAR(interpolated, monomial_value_3d_for_test(xi, exponent),
                            double(1e-9))
                    << "order=" << order << " ax=" << exponent[0]
                    << " ay=" << exponent[1] << " az=" << exponent[2];
            }
        }
    }
}


TEST(SerendipityBasis, NamedHexLayoutsMatchTopologyConstruction) {
    const struct Case { ElementType type; int order; } cases[] = {
        {ElementType::Hex8, 1},
        {ElementType::Hex20, 2},
    };
    const std::vector<math::Vector<double, 3>> points = {
        {double(0.2), double(-0.1), double(0.3)},
        {double(-0.35), double(0.25), double(-0.15)},
        {double(0.11), double(0.23), double(-0.42)},
    };
    for (const auto& c : cases) {
        SerendipityBasis named(c.type, c.order);
        SerendipityBasis topo(BasisTopology::Hexahedron, c.order);

        ASSERT_EQ(named.size(), topo.size());
        ASSERT_EQ(named.nodes().size(), topo.nodes().size());
        for (std::size_t i = 0; i < named.nodes().size(); ++i) {
            for (std::size_t d = 0; d < 3u; ++d) {
                EXPECT_EQ(named.nodes()[i][d], topo.nodes()[i][d])
                    << "node=" << i << " d=" << d;
            }
        }

        for (const auto& xi : points) {
            std::vector<double> nv, tv;
            std::vector<Gradient> ng, tg;
            std::vector<Hessian> nh, th;
            named.evaluate_all(xi, nv, ng, nh);
            topo.evaluate_all(xi, tv, tg, th);
            ASSERT_EQ(nv.size(), tv.size());
            for (std::size_t i = 0; i < nv.size(); ++i) {
                EXPECT_EQ(nv[i], tv[i]) << "value i=" << i;
                for (std::size_t d = 0; d < 3u; ++d) {
                    EXPECT_EQ(ng[i][d], tg[i][d]) << "grad i=" << i << " d=" << d;
                    for (std::size_t e = 0; e < 3u; ++e) {
                        EXPECT_EQ(nh[i](d, e), th[i](d, e))
                            << "hess i=" << i << " (" << d << "," << e << ")";
                    }
                }
            }
        }
    }
}

// Conditioning is a tested quantity, not a tolerance that quietly loosens. With
// the Legendre modal basis and Gauss-Lobatto-Legendre nodes, both the Vandermonde
// condition number and the Lebesgue constant stay small across the recommended
// range -- a logarithmic-style growth instead of the exponential blow-up of the
// previous equispaced/monomial construction (which lost ~8 digits by order 10).
TEST(SerendipityBasis, SerendipityStaysWellConditionedAcrossRecommendedRange) {
    for (int order = 1; order <= 10; ++order) {
        SerendipityBasis basis(BasisTopology::Quadrilateral, order);
        const double cond = legendre_vandermonde_condition(
            basis.nodes(), quad_serendipity_modes_3d_for_test(order));
        const double lebesgue = serendipity_lebesgue_constant(basis, 24);
        EXPECT_LT(cond, double(2.5e4)) << "quad order=" << order;
        EXPECT_LT(lebesgue, double(9e2)) << "quad order=" << order;
    }
    for (int order = 1; order <= 8; ++order) {
        SerendipityBasis basis(BasisTopology::Hexahedron, order);
        const double cond = legendre_vandermonde_condition(
            basis.nodes(), hex_serendipity_exponents_for_test(order));
        const double lebesgue = serendipity_lebesgue_constant(basis, 12);
        EXPECT_LT(cond, double(2e4)) << "hex order=" << order;
        EXPECT_LT(lebesgue, double(3.5e2)) << "hex order=" << order;
    }
}

// The condition-number guard is the numerical-soundness backstop: orders pushed
// far past the well-conditioned range throw rather than return shape functions
// whose coefficients have lost all precision. The recommended orders construct
// without complaint.
TEST(SerendipityBasis, RejectsOrdersBeyondTheWellConditionedRange) {
    EXPECT_NO_THROW((void)SerendipityBasis(BasisTopology::Quadrilateral, 10));
    EXPECT_NO_THROW((void)SerendipityBasis(BasisTopology::Hexahedron, 8));
    EXPECT_THROW((void)SerendipityBasis(BasisTopology::Quadrilateral, 20),
                 BasisConstructionException);
    EXPECT_THROW((void)SerendipityBasis(BasisTopology::Hexahedron, 16),
                 BasisConstructionException);
}
