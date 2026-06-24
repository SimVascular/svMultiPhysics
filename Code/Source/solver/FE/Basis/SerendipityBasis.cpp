// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "SerendipityBasis.h"
#include "NodeOrderingConventions.h"
#include "Math/DenseLinearAlgebra.h"

#include <algorithm>
#include <array>
#include <span>
#include <string>

namespace svmp {
namespace FE {
namespace basis {

namespace {
using Vec3 = math::Vector<double, 3>;

int quad_serendipity_superlinear_degree(int ax, int ay) {
    return (ax > 1 ? ax : 0) + (ay > 1 ? ay : 0);
}

inline double integer_power(double base, int exponent) {
    double result = double(1);
    for (int k = 0; k < exponent; ++k) {
        result *= base;
    }
    return result;
}

std::vector<std::array<int, 2>> quad_serendipity_exponents(int order) {
    std::vector<std::array<int, 2>> exponents;
    for (int ay = 0; ay <= order; ++ay) {
        for (int ax = 0; ax <= order; ++ax) {
            if (quad_serendipity_superlinear_degree(ax, ay) <= order) {
                exponents.push_back({ax, ay});
            }
        }
    }
    return exponents;
}

std::size_t quad_serendipity_interior_count(int order) {
    if (order < 4) {
        return 0u;
    }
    const auto m = static_cast<std::size_t>(order - 4);
    return (m + 1u) * (m + 2u) / 2u;
}

// Interior nodes are a triangular row set for P_m, m = order - 4. If a
// serendipity polynomial vanishes at the p + 1 boundary nodes on each edge,
// each edge restriction is identically zero and the polynomial factors as
// (1 - x^2)(1 - y^2) q with q in P_m. Row 0 has m + 1 distinct x-values; if q
// vanishes there, q(x, y_0) is the zero one-variable polynomial and
// q = (y - y_0) q_1 with q_1 in P_{m-1}. Repeating over the remaining rows
// proves q = 0, so the full quadrilateral serendipity Vandermonde is
// nonsingular for this node set.
void append_quad_serendipity_interior_nodes(std::vector<Vec3>& nodes, int order) {
    if (order < 4) {
        return;
    }

    const int m = order - 4;
    const double y_denominator = double(m + 2);
    for (int row = 0; row <= m; ++row) {
        const int row_count = m + 1 - row;
        const double y = double(-1) + double(2) * double(row + 1) / y_denominator;
        const double x_denominator = double(row_count + 1);
        for (int col = 0; col < row_count; ++col) {
            const double x = double(-1) + double(2) * double(col + 1) / x_denominator;
            nodes.push_back(Vec3{x, y, double(0)});
        }
    }
}

std::vector<Vec3> quad_serendipity_nodes(int order, std::size_t total_size) {
    if (order <= 0) {
        return {};
    }

    // The corner+edge skeleton is the leading prefix of the complete quadrilateral
    // Lagrange layout of the same order: 4 corners followed by 4(order-1) edge
    // nodes, in the same VTK boundary order. Source it from the single
    // ReferenceNodeLayout generator and drop that layout's interior, so the
    // reference-cell corner/edge geometry has one owner; only the reduced interior
    // appended below is serendipity-specific.
    std::vector<Vec3> nodes =
        ReferenceNodeLayout::get_lagrange_node_coords(ElementType::Quad4, order);
    const std::size_t boundary_count = static_cast<std::size_t>(4 * order);
    svmp::throw_if<BasisConstructionException>(
        boundary_count > nodes.size(), SVMP_HERE,
        "SerendipityBasis: quadrilateral skeleton exceeds the complete Lagrange layout");
    nodes.resize(boundary_count);

    svmp::throw_if<BasisConstructionException>(
        nodes.size() > total_size, SVMP_HERE,
        "SerendipityBasis: quadrilateral serendipity boundary nodes exceed requested size");

    const std::size_t interior_count = quad_serendipity_interior_count(order);
    svmp::throw_if<BasisConstructionException>(
        nodes.size() + interior_count != total_size, SVMP_HERE,
        "SerendipityBasis: quadrilateral serendipity monomial/node count mismatch");

    append_quad_serendipity_interior_nodes(nodes, order);
    return nodes;
}

int hex_serendipity_superlinear_degree(int ax, int ay, int az) {
    return (ax > 1 ? ax : 0) + (ay > 1 ? ay : 0) + (az > 1 ? az : 0);
}

// Hexahedral serendipity monomial space: every r^a s^b t^c whose superlinear
// degree is at most `order`. This is the three-axis generalization of
// quad_serendipity_exponents; at order 1 it is the eight multilinear monomials
// (the Hex8 space) and at order 2 it is the twenty-monomial Hex20 space. The
// enumeration order is internal -- evaluation sums over the monomials, so only
// the node order, not the monomial order, is observable to a caller.
std::vector<std::array<int, 3>> hex_serendipity_exponents(int order) {
    std::vector<std::array<int, 3>> exponents;
    for (int az = 0; az <= order; ++az) {
        for (int ay = 0; ay <= order; ++ay) {
            for (int ax = 0; ax <= order; ++ax) {
                if (hex_serendipity_superlinear_degree(ax, ay, az) <= order) {
                    exponents.push_back({ax, ay, az});
                }
            }
        }
    }
    return exponents;
}

// Volume-interior node count for hexahedral serendipity. Once the boundary trace
// is fixed, an interior serendipity function factors as the cube bubble
// (1 - r^2)(1 - s^2)(1 - t^2) times a quotient of total degree at most order - 6,
// so the interior space is P_{order-6} in three variables: empty until order 6,
// then dim P_{order-6} = (m+1)(m+2)(m+3)/6 with m = order - 6.
std::size_t hex_serendipity_volume_interior_count(int order) {
    if (order < 6) {
        return 0u;
    }
    const auto m = static_cast<std::size_t>(order - 6);
    return (m + 1u) * (m + 2u) * (m + 3u) / 6u;
}

// Append the face-interior nodes. The restriction of the order-`order` cube
// serendipity space to a face is the order-`order` quadrilateral serendipity
// space, so every face carries the same 2D quad-serendipity interior set,
// embedded into the face plane. Faces are visited in VTK face order
// (-X, +X, -Y, +Y, -Z, +Z); the in-plane (u, v) point maps to the two free axes
// of each face. Empty until order 4 (when the quad interior first appears).
void append_hex_serendipity_face_interior_nodes(std::vector<Vec3>& nodes, int order) {
    std::vector<Vec3> face_interior;  // (u, v, 0) interior points of one quad face
    append_quad_serendipity_interior_nodes(face_interior, order);
    if (face_interior.empty()) {
        return;
    }

    // Each face: the fixed axis (0 = r, 1 = s, 2 = t), its +/-1 value, and the two
    // in-plane axes that carry the 2D interior point (u, v).
    struct Face {
        int fixed_axis;
        double fixed_value;
        int u_axis;
        int v_axis;
    };
    static constexpr Face faces[6] = {
        {0, double(-1), 1, 2},  // -X: (s, t) in plane
        {0, double(1),  1, 2},  // +X
        {1, double(-1), 0, 2},  // -Y: (r, t) in plane
        {1, double(1),  0, 2},  // +Y
        {2, double(-1), 0, 1},  // -Z: (r, s) in plane
        {2, double(1),  0, 1},  // +Z
    };

    for (const auto& face : faces) {
        for (const auto& p : face_interior) {
            Vec3 node = Vec3::Zero();
            node[static_cast<std::size_t>(face.fixed_axis)] = face.fixed_value;
            node[static_cast<std::size_t>(face.u_axis)] = p[0];
            node[static_cast<std::size_t>(face.v_axis)] = p[1];
            nodes.push_back(node);
        }
    }
}

// Append the volume-interior nodes: a tetrahedral staircase unisolvent for the
// interior residual P_{order-6}. Each t-layer is a triangular staircase (the 2D
// construction reused) whose total degree decreases by one per layer, so the
// layers consume P_{order-6} by induction in t exactly as the quad interior
// consumes P_{order-4} by induction in s. Empty until order 6.
void append_hex_serendipity_volume_interior_nodes(std::vector<Vec3>& nodes, int order) {
    if (order < 6) {
        return;
    }
    const int m = order - 6;
    for (int layer = 0; layer <= m; ++layer) {
        const int tri_order = m - layer;
        const double t = double(-1) + double(2) * double(layer + 1) / double(m + 2);
        for (int row = 0; row <= tri_order; ++row) {
            const int row_count = tri_order + 1 - row;
            const double s = double(-1) + double(2) * double(row + 1) / double(tri_order + 2);
            for (int col = 0; col < row_count; ++col) {
                const double r = double(-1) + double(2) * double(col + 1) / double(row_count + 1);
                nodes.push_back(Vec3{r, s, t});
            }
        }
    }
}

// Generate the hexahedral serendipity reference nodes in the generalized
// right-hand-rule / VTK-consistent stratified order: 8 corners, then 12 edges in
// VTK quadratic-hex edge order, then the 6 face interiors in VTK face order, then
// the volume interior. The corner and edge strata are taken directly from the
// complete hexahedral Lagrange layout (generate_hex_nodes, via ReferenceNodeLayout),
// so they share that single generator's VTK ordering: at order 1 (corners only)
// and order 2 (corners plus edge midpoints) the layout is exactly the public
// Hex8 / Hex20 ordering, and for higher order the reduced face/volume sets are
// this module's own convention.
std::vector<Vec3> hex_serendipity_nodes(int order, std::size_t total_size) {
    std::vector<Vec3> nodes =
        ReferenceNodeLayout::get_lagrange_node_coords(ElementType::Hex8, order);
    const std::size_t skeleton_count =
        8u + 12u * static_cast<std::size_t>(order - 1);
    svmp::throw_if<BasisConstructionException>(
        skeleton_count > nodes.size(), SVMP_HERE,
        "SerendipityBasis: hexahedral skeleton exceeds the complete Lagrange layout");
    nodes.resize(skeleton_count);

    const std::size_t skeleton = nodes.size();
    append_hex_serendipity_face_interior_nodes(nodes, order);
    svmp::throw_if<BasisConstructionException>(
        nodes.size() - skeleton != 6u * quad_serendipity_interior_count(order), SVMP_HERE,
        "SerendipityBasis: hexahedral serendipity face-interior node count mismatch");

    const std::size_t before_volume = nodes.size();
    append_hex_serendipity_volume_interior_nodes(nodes, order);
    svmp::throw_if<BasisConstructionException>(
        nodes.size() - before_volume != hex_serendipity_volume_interior_count(order), SVMP_HERE,
        "SerendipityBasis: hexahedral serendipity volume-interior node count mismatch");

    svmp::throw_if<BasisConstructionException>(
        nodes.size() != total_size, SVMP_HERE,
        "SerendipityBasis: hexahedral serendipity node count does not match the monomial count");
    return nodes;
}

// Build the nodal coefficient table for a monomial-generated serendipity family:
// assemble V[node][monomial] = r^a s^b t^c at the public-order reference nodes and
// invert it. Because the nodes are in public order, the inverse is already in
// public basis order and needs no output permutation. The same routine serves the
// quadrilateral, Hex20, and Wedge15 spaces.
std::vector<double> build_inverse_vandermonde(
    std::span<const Vec3> nodes,
    std::span<const std::array<int, 3>> exponents,
    const std::string& label) {
    const std::size_t n = nodes.size();
    svmp::throw_if<BasisConstructionException>(
        n == 0 || exponents.size() != n, SVMP_HERE,
        "SerendipityBasis: invalid serendipity interpolation setup");

    std::vector<double> vandermonde(n * n, double(0));
    for (std::size_t row = 0; row < n; ++row) {
        const Vec3& p = nodes[row];
        for (std::size_t col = 0; col < n; ++col) {
            const auto& e = exponents[col];
            vandermonde[row * n + col] =
                integer_power(p[0], e[0]) * integer_power(p[1], e[1]) *
                integer_power(p[2], e[2]);
        }
    }

    return math::invert_dense_matrix(
        std::move(vandermonde), n,
        "SerendipityBasis interpolation matrix for " + label);
}

constexpr std::array<std::array<int, 3>, 15> kWedge15MonomialExponents = {{
    {{0, 0, 0}},
    {{0, 0, 1}},
    {{0, 0, 2}},
    {{0, 1, 0}},
    {{0, 1, 1}},
    {{0, 1, 2}},
    {{0, 2, 0}},
    {{0, 2, 1}},
    {{1, 0, 0}},
    {{1, 0, 1}},
    {{1, 0, 2}},
    {{1, 1, 0}},
    {{1, 1, 1}},
    {{2, 0, 0}},
    {{2, 0, 1}}
}};

// Value and first/second derivatives of the 1D monomial x^a. The derivative of
// a constant or linear term collapses to zero, so negative powers never arise.
struct MonomialAxis {
    double value;   ///< x^a
    double first;   ///< d/dx (x^a)     = a x^(a-1)
    double second;  ///< d^2/dx^2 (x^a) = a (a-1) x^(a-2)
};

inline MonomialAxis monomial_axis(double x, int exponent) {
    MonomialAxis axis;
    axis.value = integer_power(x, exponent);
    axis.first = (exponent > 0) ? double(exponent) * integer_power(x, exponent - 1) : double(0);
    axis.second = (exponent > 1)
                      ? double(exponent * (exponent - 1)) * integer_power(x, exponent - 2)
                      : double(0);
    return axis;
}

// Evaluate a nodal basis defined by a monomial coefficient table at one
// reference point. For each monomial j the routine forms x^a y^b z^c and the
// requested derivatives, then accumulates the coefficient-weighted sum into the
// output slots. `count` is both the number of monomials and the number of basis
// functions (the coefficient table is square). The table is in public basis
// order, so output slot i reads coefficient column i directly. Outputs are
// assumed pre-zeroed by the caller; an empty span skips that quantity.
template <typename ExponentFn, typename CoeffFn>
void eval_monomial_basis(double r, double s, double t,
                         std::size_t count,
                         ExponentFn&& exponent,
                         CoeffFn&& coeff,
                         std::span<double> values,
                         std::span<Gradient> gradients,
                         std::span<Hessian> hessians) {
    const bool want_values = !values.empty();
    const bool want_gradients = !gradients.empty();
    const bool want_hessians = !hessians.empty();

    for (std::size_t j = 0; j < count; ++j) {
        const std::array<int, 3> e = exponent(j);
        const MonomialAxis ax = monomial_axis(r, e[0]);
        const MonomialAxis ay = monomial_axis(s, e[1]);
        const MonomialAxis az = monomial_axis(t, e[2]);

        const double phi = ax.value * ay.value * az.value;

        double d_dr = double(0), d_ds = double(0), d_dt = double(0);
        if (want_gradients || want_hessians) {
            d_dr = ax.first * ay.value * az.value;
            d_ds = ax.value * ay.first * az.value;
            d_dt = ax.value * ay.value * az.first;
        }

        double d_drr = double(0), d_dss = double(0), d_dtt = double(0);
        double d_drs = double(0), d_drt = double(0), d_dst = double(0);
        if (want_hessians) {
            d_drr = ax.second * ay.value * az.value;
            d_dss = ax.value * ay.second * az.value;
            d_dtt = ax.value * ay.value * az.second;
            d_drs = ax.first * ay.first * az.value;
            d_drt = ax.first * ay.value * az.first;
            d_dst = ax.value * ay.first * az.first;
        }

        for (std::size_t slot = 0; slot < count; ++slot) {
            const double c = coeff(j, slot);
            if (want_values) {
                values[slot] += c * phi;
            }
            if (want_gradients) {
                Gradient& g = gradients[slot];
                g[0] += c * d_dr;
                g[1] += c * d_ds;
                g[2] += c * d_dt;
            }
            if (want_hessians) {
                Hessian& h = hessians[slot];
                h(0, 0) += c * d_drr;
                h(1, 1) += c * d_dss;
                h(2, 2) += c * d_dtt;
                h(0, 1) += c * d_drs;
                h(1, 0) += c * d_drs;
                h(0, 2) += c * d_drt;
                h(2, 0) += c * d_drt;
                h(1, 2) += c * d_dst;
                h(2, 1) += c * d_dst;
            }
        }
    }
}

struct NormalizedSerendipityRequest {
    BasisTopology topology;
    int dimension;
    int order;
};

// Validate a named serendipity element/order pairing and return its topology,
// reference dimension, and order. The named serendipity layouts (Quad8, Hex8,
// Hex20, Wedge15) are each pinned to a single polynomial order by their node
// count, so a mismatched explicit order is rejected. Arbitrary-order
// quadrilateral serendipity is not a named element: it is requested through the
// BasisTopology::Quadrilateral constructor.
NormalizedSerendipityRequest normalize_serendipity_request(ElementType type, int order) {
    // The named layouts carry an inferred fixed order (Hex8 -> 1; Quad8, Hex20,
    // and Wedge15 -> 2). The request must supply that exact order: it is never
    // floored or otherwise adjusted to fit, so order 0 and negative orders are
    // rejected rather than promoted to a valid layout.
    switch (type) {
        case ElementType::Quad8:
            svmp::throw_if<BasisConfigurationException>(order != 2, SVMP_HERE,
                "SerendipityBasis: Quad8 is the quadratic 8-node serendipity layout (order 2 only); "
                "use BasisTopology::Quadrilateral for higher-order quadrilateral serendipity");
            return {BasisTopology::Quadrilateral, 2, 2};
        case ElementType::Hex8:
            svmp::throw_if<BasisConfigurationException>(order != 1, SVMP_HERE,
                "SerendipityBasis: Hex8 is the trilinear 8-node basis (order 1 only); use Hex20 for quadratic serendipity");
            return {BasisTopology::Hexahedron, 3, 1};
        case ElementType::Hex20:
            svmp::throw_if<BasisConfigurationException>(order != 2, SVMP_HERE,
                "SerendipityBasis: Hex20 is the 20-node quadratic serendipity layout (order 2 only)");
            return {BasisTopology::Hexahedron, 3, 2};
        case ElementType::Wedge15:
            svmp::throw_if<BasisConfigurationException>(order != 2, SVMP_HERE,
                "SerendipityBasis: Wedge15 is the 15-node quadratic serendipity layout (order 2 only)");
            return {BasisTopology::Wedge, 3, 2};
        default:
            svmp::raise<BasisElementCompatibilityException>(SVMP_HERE,
                "SerendipityBasis named elements are Quad8, Hex8, Hex20, and Wedge15; "
                "use BasisTopology::Quadrilateral for arbitrary-order quadrilateral serendipity");
    }
}

} // namespace

SerendipityBasis::SerendipityBasis(BasisTopology topology, int order)
    : topology_(topology), dimension_(0), order_(0), size_(0) {
    const bool supported_topology = topology_ == BasisTopology::Quadrilateral ||
                                    topology_ == BasisTopology::Hexahedron;
    svmp::throw_if<BasisElementCompatibilityException>(
        !supported_topology, SVMP_HERE,
        "SerendipityBasis: arbitrary-order topology construction is supported for "
        "Quadrilateral and Hexahedron; use the named ElementType (Wedge15) for wedge serendipity");
    svmp::throw_if<BasisConfigurationException>(
        order < 1, SVMP_HERE,
        "SerendipityBasis: serendipity requires a polynomial order >= 1");
    dimension_ = topology_ == BasisTopology::Hexahedron ? 3 : 2;
    order_ = order;
    if (topology_ == BasisTopology::Hexahedron) {
        init_hexahedron(order_, /*nodes_from_reference_layout=*/false);
    } else {
        init_quadrilateral(order_, /*nodes_from_reference_layout=*/false);
    }
}

SerendipityBasis::SerendipityBasis(ElementType type, int order)
    : topology_(BasisTopology::Unknown), dimension_(0), order_(0), size_(0) {
    const NormalizedSerendipityRequest normalized = normalize_serendipity_request(type, order);
    topology_ = normalized.topology;
    dimension_ = normalized.dimension;
    order_ = normalized.order;

    switch (type) {
        case ElementType::Quad8:
            // Quad8 is the named quadratic layout; its nodes come from
            // ReferenceNodeLayout so the basis shares the single public Quad8
            // ordering (the same source Hex8/Hex20/Wedge15 use).
            init_quadrilateral(order_, /*nodes_from_reference_layout=*/true);
            return;
        case ElementType::Hex8:
            // Hex8 is the order-1 instance of the hexahedral serendipity space.
            init_hexahedron(1, /*nodes_from_reference_layout=*/true);
            return;
        case ElementType::Hex20:
            // Hex20 is the order-2 instance of the hexahedral serendipity space.
            init_hexahedron(2, /*nodes_from_reference_layout=*/true);
            return;
        case ElementType::Wedge15:
            init_fixed_named(type);
            return;
        default:
            // normalize_serendipity_request already rejected anything else.
            svmp::raise<BasisElementCompatibilityException>(SVMP_HERE,
                "SerendipityBasis: unsupported named serendipity element");
    }
}

// Build the quadrilateral serendipity monomial space, reference nodes, and nodal
// coefficient table for the given order. The coefficient table is the inverse
// Vandermonde of those monomials at the reference nodes; because the nodes are
// in public order, evaluation needs no output permutation. The named Quad8
// layout sources its nodes from ReferenceNodeLayout; the arbitrary-order
// topology path generates them.
void SerendipityBasis::init_quadrilateral(int order, bool nodes_from_reference_layout) {
    const auto quad_exponents = quad_serendipity_exponents(order);
    monomial_exponents_.clear();
    monomial_exponents_.reserve(quad_exponents.size());
    for (const auto& e : quad_exponents) {
        monomial_exponents_.push_back({e[0], e[1], 0});
    }
    size_ = monomial_exponents_.size();
    nodes_ = nodes_from_reference_layout
                 ? ReferenceNodeLayout::node_coords(ElementType::Quad8)
                 : quad_serendipity_nodes(order, size_);
    svmp::throw_if<BasisConstructionException>(
        nodes_.size() != size_, SVMP_HERE,
        "SerendipityBasis: quadrilateral serendipity setup produced inconsistent sizes");
    inv_vandermonde_ = build_inverse_vandermonde(
        nodes_, monomial_exponents_, "Quad order " + std::to_string(order));
}

// Build the hexahedral serendipity monomial space, reference nodes, and nodal
// coefficient table for the given order, mirroring init_quadrilateral. The
// arbitrary-order topology path generates its own VTK-consistent nodes; the named
// Hex8 (order 1) and Hex20 (order 2) layouts source their public-order nodes from
// ReferenceNodeLayout so the generated layout matches the public ordering exactly.
void SerendipityBasis::init_hexahedron(int order, bool nodes_from_reference_layout) {
    monomial_exponents_ = hex_serendipity_exponents(order);
    size_ = monomial_exponents_.size();
    if (nodes_from_reference_layout) {
        const ElementType named =
            (order == 1) ? ElementType::Hex8 : ElementType::Hex20;
        nodes_ = ReferenceNodeLayout::node_coords(named);
    } else {
        nodes_ = hex_serendipity_nodes(order, size_);
    }
    svmp::throw_if<BasisConstructionException>(
        nodes_.size() != size_, SVMP_HERE,
        "SerendipityBasis: hexahedral serendipity setup produced inconsistent sizes");
    inv_vandermonde_ = build_inverse_vandermonde(
        nodes_, monomial_exponents_, "Hex order " + std::to_string(order));
}

// Build the Wedge15 serendipity layout from its tabulated monomial space and
// public-order ReferenceNodeLayout nodes. Hexahedral serendipity (Hex8 and Hex20
// included) is generated by init_hexahedron, so the prism is the only named
// layout that still carries a fixed monomial table.
void SerendipityBasis::init_fixed_named(ElementType type) {
    svmp::throw_if<BasisConstructionException>(
        type != ElementType::Wedge15, SVMP_HERE,
        "SerendipityBasis: init_fixed_named builds only the Wedge15 layout");
    size_ = 15u;
    const std::span<const std::array<int, 3>> family_exponents(
        kWedge15MonomialExponents.data(), kWedge15MonomialExponents.size());
    nodes_ = ReferenceNodeLayout::node_coords(type);
    svmp::throw_if<BasisConstructionException>(
        nodes_.size() != size_, SVMP_HERE,
        "SerendipityBasis: Wedge15 layout node count does not match basis size");
    svmp::throw_if<BasisConstructionException>(
        family_exponents.size() != size_, SVMP_HERE,
        "SerendipityBasis: Wedge15 monomial count does not match basis size");
    monomial_exponents_.assign(family_exponents.begin(), family_exponents.end());
    inv_vandermonde_ = build_inverse_vandermonde(nodes_, monomial_exponents_, "Wedge15");
}

void SerendipityBasis::evaluate_all_to(const math::Vector<double, 3>& xi,
                                       std::span<double> values_out,
                                       std::span<Gradient> gradients_out,
                                       std::span<Hessian> hessians_out) const {
    require_requested_span_size(values_out, size_, "SerendipityBasis::evaluate_all_to values");
    require_requested_span_size(gradients_out, size_, "SerendipityBasis::evaluate_all_to gradients");
    require_requested_span_size(hessians_out, size_, "SerendipityBasis::evaluate_all_to hessians");

    if (values_out.empty() && gradients_out.empty() && hessians_out.empty()) {
        return;
    }

    if (!values_out.empty()) {
        std::fill(values_out.begin(), values_out.end(), double(0));
    }
    if (!gradients_out.empty()) {
        std::fill(gradients_out.begin(), gradients_out.end(), Gradient::Zero());
    }
    if (!hessians_out.empty()) {
        std::fill(hessians_out.begin(), hessians_out.end(), Hessian::Zero());
    }

    const double x = xi[0];
    const double y = xi[1];
    const double z = xi[2];

    // Every serendipity family evaluates through its generated coefficient table,
    // which is already in public basis order.
    svmp::throw_if<BasisEvaluationException>(
        monomial_exponents_.size() != size_ ||
            inv_vandermonde_.size() != size_ * size_,
        SVMP_HERE,
        "SerendipityBasis: interpolation tables are not initialized for evaluation");

    eval_monomial_basis(
        x, y, z, size_,
        [this](std::size_t j) { return monomial_exponents_[j]; },
        [this](std::size_t j, std::size_t i) {
            return inv_vandermonde_[j * size_ + i];
        },
        values_out, gradients_out, hessians_out);
}

void SerendipityBasis::evaluate_values(const math::Vector<double, 3>& xi,
                                       std::vector<double>& values) const {
    values.resize(size_);
    evaluate_values_to(xi, std::span<double>(values.data(), values.size()));
}

void SerendipityBasis::evaluate_gradients(const math::Vector<double, 3>& xi,
                                          std::vector<Gradient>& gradients) const {
    gradients.resize(size_);
    evaluate_gradients_to(xi, std::span<Gradient>(gradients.data(), gradients.size()));
}

void SerendipityBasis::evaluate_hessians(const math::Vector<double, 3>& xi,
                                         std::vector<Hessian>& hessians) const {
    hessians.resize(size_);
    evaluate_hessians_to(xi, std::span<Hessian>(hessians.data(), hessians.size()));
}

void SerendipityBasis::evaluate_all(const math::Vector<double, 3>& xi,
                                    std::vector<double>& values,
                                    std::vector<Gradient>& gradients,
                                    std::vector<Hessian>& hessians) const {
    values.resize(size_);
    gradients.resize(size_);
    hessians.resize(size_);
    evaluate_all_to(xi,
                    std::span<double>(values.data(), values.size()),
                    std::span<Gradient>(gradients.data(), gradients.size()),
                    std::span<Hessian>(hessians.data(), hessians.size()));
}

void SerendipityBasis::evaluate_values_to(const math::Vector<double, 3>& xi,
                                          std::span<double> values_out) const {
    require_span_size(values_out.size(), size_, "SerendipityBasis::evaluate_values_to");
    evaluate_all_to(xi, values_out, std::span<Gradient>{}, std::span<Hessian>{});
}

void SerendipityBasis::evaluate_gradients_to(const math::Vector<double, 3>& xi,
                                             std::span<Gradient> gradients_out) const {
    require_span_size(gradients_out.size(), size_, "SerendipityBasis::evaluate_gradients_to");
    evaluate_all_to(xi, std::span<double>{}, gradients_out, std::span<Hessian>{});
}

void SerendipityBasis::evaluate_hessians_to(const math::Vector<double, 3>& xi,
                                            std::span<Hessian> hessians_out) const {
    require_span_size(hessians_out.size(), size_, "SerendipityBasis::evaluate_hessians_to");
    evaluate_all_to(xi, std::span<double>{}, std::span<Gradient>{}, hessians_out);
}

} // namespace basis
} // namespace FE
} // namespace svmp
