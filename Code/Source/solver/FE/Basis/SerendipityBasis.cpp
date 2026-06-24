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

void evaluate_hex8_reference(double r,
                             double s,
                             double t,
                             std::span<double> values,
                             std::span<Gradient> gradients,
                             std::span<Hessian> hessians) {
    static constexpr int signs[8][3] = {
        {-1, -1, -1},
        { 1, -1, -1},
        { 1,  1, -1},
        {-1,  1, -1},
        {-1, -1,  1},
        { 1, -1,  1},
        { 1,  1,  1},
        {-1,  1,  1},
    };

    for (std::size_t i = 0; i < 8u; ++i) {
        const double a = double(signs[i][0]);
        const double b = double(signs[i][1]);
        const double c = double(signs[i][2]);
        const double ar = double(1) + a * r;
        const double bs = double(1) + b * s;
        const double ct = double(1) + c * t;

        if (!values.empty()) {
            values[i] = double(0.125) * ar * bs * ct;
        }
        if (!gradients.empty()) {
            Gradient& g = gradients[i];
            g[0] = double(0.125) * a * bs * ct;
            g[1] = double(0.125) * b * ar * ct;
            g[2] = double(0.125) * c * ar * bs;
        }
        if (!hessians.empty()) {
            Hessian& h = hessians[i];
            h(0, 0) = double(0);
            h(0, 1) = double(0.125) * a * b * ct;
            h(0, 2) = double(0.125) * a * c * bs;
            h(1, 0) = h(0, 1);
            h(1, 1) = double(0);
            h(1, 2) = double(0.125) * b * c * ar;
            h(2, 0) = h(0, 2);
            h(2, 1) = h(1, 2);
            h(2, 2) = double(0);
        }
    }
}

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
    std::vector<Vec3> nodes;
    if (order <= 0) {
        return nodes;
    }

    const double inv_order = double(1) / double(order);

    nodes.push_back(Vec3{double(-1), double(-1), double(0)});
    nodes.push_back(Vec3{double(1),  double(-1), double(0)});
    nodes.push_back(Vec3{double(1),  double(1),  double(0)});
    nodes.push_back(Vec3{double(-1), double(1),  double(0)});

    for (int i = 1; i < order; ++i) {
        nodes.push_back(Vec3{double(-1) + double(2 * i) * inv_order, double(-1), double(0)});
    }
    for (int i = 1; i < order; ++i) {
        nodes.push_back(Vec3{double(1), double(-1) + double(2 * i) * inv_order, double(0)});
    }
    for (int i = 1; i < order; ++i) {
        nodes.push_back(Vec3{double(1) - double(2 * i) * inv_order, double(1), double(0)});
    }
    for (int i = 1; i < order; ++i) {
        nodes.push_back(Vec3{double(-1), double(1) - double(2 * i) * inv_order, double(0)});
    }

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

constexpr std::array<std::array<int, 3>, 20> kHex20MonomialExponents = {{
    {{0, 0, 0}}, {{0, 0, 1}}, {{0, 0, 2}}, {{0, 1, 0}}, {{0, 1, 1}},
    {{0, 1, 2}}, {{0, 2, 0}}, {{0, 2, 1}}, {{1, 0, 0}}, {{1, 0, 1}},
    {{1, 0, 2}}, {{1, 1, 0}}, {{1, 1, 1}}, {{1, 1, 2}}, {{1, 2, 0}},
    {{1, 2, 1}}, {{2, 0, 0}}, {{2, 0, 1}}, {{2, 1, 0}}, {{2, 1, 1}}
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
    const int floored_order = std::max(order, 1);
    switch (type) {
        case ElementType::Quad8:
            svmp::throw_if<BasisConfigurationException>(floored_order != 2, SVMP_HERE,
                "SerendipityBasis: Quad8 is the quadratic 8-node serendipity layout (order 2 only); "
                "use BasisTopology::Quadrilateral for higher-order quadrilateral serendipity");
            return {BasisTopology::Quadrilateral, 2, 2};
        case ElementType::Hex8:
            svmp::throw_if<BasisConfigurationException>(floored_order != 1, SVMP_HERE,
                "SerendipityBasis: Hex8 is the trilinear 8-node basis (order 1 only); use Hex20 for quadratic serendipity");
            return {BasisTopology::Hexahedron, 3, 1};
        case ElementType::Hex20:
            svmp::throw_if<BasisConfigurationException>(floored_order != 2, SVMP_HERE,
                "SerendipityBasis: Hex20 is the 20-node quadratic serendipity layout (order 2 only)");
            return {BasisTopology::Hexahedron, 3, 2};
        case ElementType::Wedge15:
            svmp::throw_if<BasisConfigurationException>(floored_order != 2, SVMP_HERE,
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
    svmp::throw_if<BasisElementCompatibilityException>(
        topology_ != BasisTopology::Quadrilateral, SVMP_HERE,
        "SerendipityBasis: arbitrary-order topology construction is only supported for "
        "Quadrilateral; use the named ElementType (Hex8/Hex20/Wedge15) for hex/wedge serendipity");
    dimension_ = 2;
    order_ = std::max(order, 1);
    init_quadrilateral(order_, /*nodes_from_reference_layout=*/false);
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
            // ordering (the same source Hex20/Wedge15 use).
            init_quadrilateral(order_, /*nodes_from_reference_layout=*/true);
            return;
        case ElementType::Hex8:
            // Hex8 is the standard trilinear corner basis, evaluated directly
            // rather than through a generated coefficient table.
            size_ = 8u;
            nodes_ = ReferenceNodeLayout::node_coords(type);
            svmp::throw_if<BasisConstructionException>(
                nodes_.size() != size_, SVMP_HERE,
                "SerendipityBasis: Hex8 layout node count does not match basis size");
            return;
        case ElementType::Hex20:
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

// Build a fixed named volume serendipity layout (Hex20 or Wedge15). The nodal
// coefficient table is generated by inverting the Vandermonde built from the
// public-order ReferenceNodeLayout nodes, exactly like the quadrilateral, so no
// transcribed tables or output permutation are needed.
void SerendipityBasis::init_fixed_named(ElementType type) {
    std::span<const std::array<int, 3>> family_exponents;
    std::string label;
    if (type == ElementType::Hex20) {
        size_ = 20u;
        family_exponents = std::span<const std::array<int, 3>>(
            kHex20MonomialExponents.data(), kHex20MonomialExponents.size());
        label = "Hex20";
    } else {  // Wedge15
        size_ = 15u;
        family_exponents = std::span<const std::array<int, 3>>(
            kWedge15MonomialExponents.data(), kWedge15MonomialExponents.size());
        label = "Wedge15";
    }
    nodes_ = ReferenceNodeLayout::node_coords(type);
    svmp::throw_if<BasisConstructionException>(
        nodes_.size() != size_, SVMP_HERE,
        "SerendipityBasis: fixed serendipity layout node count does not match basis size");
    svmp::throw_if<BasisConstructionException>(
        family_exponents.size() != size_, SVMP_HERE,
        "SerendipityBasis: serendipity monomial count does not match basis size");
    monomial_exponents_.assign(family_exponents.begin(), family_exponents.end());
    inv_vandermonde_ = build_inverse_vandermonde(nodes_, monomial_exponents_, label);
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

    // Hex8 (Hexahedron at order 1) is the only serendipity basis evaluated
    // directly from the trilinear corner products rather than a coefficient table.
    if (topology_ == BasisTopology::Hexahedron && order_ == 1) {
        evaluate_hex8_reference(x, y, z, values_out, gradients_out, hessians_out);
        return;
    }

    // Quad, Hex20, and Wedge15 evaluate through their generated coefficient
    // table, which is already in public basis order.
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
