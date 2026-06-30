// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "SerendipityBasis.h"
#include "NodeOrderingConventions.h"
#include "Math/DenseLinearAlgebra.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <span>
#include <string>

namespace svmp::FE::basis {

namespace {
using Vec3 = math::Vector<double, 3>;

// Superlinear contribution of one axis exponent: degrees 0 and 1 are free (they
// do not raise the superlinear degree), every higher degree contributes its full
// value. Summed over the active axes this gives the serendipity superlinear
// degree that bounds the mode set.
int superlinear_term(int a) {
    return a > 1 ? a : 0;
}

inline double integer_power(double base, int exponent) {
    double result = double(1);
    for (int k = 0; k < exponent; ++k) {
        result *= base;
    }
    return result;
}

// Which 1D polynomial family the tensor modes are written in: the quadrilateral
// and hexahedral families use Legendre, the fixed Wedge15 layout (order 2) keeps
// the monomial form. Both span the same serendipity polynomial space; see the
// SerendipityBasis class documentation (Modal basis) for why Legendre is used to
// keep the Vandermonde well-conditioned.
enum class ModalAxisKind { Monomial, Legendre };

// Value and first/second derivative of every 1D mode phi_0..phi_{max_degree} at a
// fixed coordinate, indexed by per-axis degree.
struct AxisTable {
    std::vector<double> value;
    std::vector<double> first;
    std::vector<double> second;
};

// phi_k(x) = x^k and its derivatives. Matches the previous per-mode evaluation
// exactly, so the Wedge15 monomial path is unchanged.
void fill_monomial_table(double x, int max_degree, AxisTable& out) {
    const std::size_t n = static_cast<std::size_t>(max_degree) + 1u;
    out.value.assign(n, double(0));
    out.first.assign(n, double(0));
    out.second.assign(n, double(0));
    for (int k = 0; k <= max_degree; ++k) {
        const std::size_t kp = static_cast<std::size_t>(k);
        out.value[kp] = integer_power(x, k);
        out.first[kp] = (k > 0) ? double(k) * integer_power(x, k - 1) : double(0);
        out.second[kp] =
            (k > 1) ? double(k * (k - 1)) * integer_power(x, k - 2) : double(0);
    }
}

// phi_k(x) = P_k(x), the degree-k Legendre polynomial on [-1, 1], with first and
// second derivatives. Built from the three-term recurrences
//   (k+1) P_{k+1}   = (2k+1) x P_k          - k P_{k-1}
//   (k+1) P'_{k+1}  = (2k+1)(P_k + x P'_k)  - k P'_{k-1}
//   (k+1) P''_{k+1} = (2k+1)(2 P'_k + x P''_k) - k P''_{k-1}
// all regular at x = +/-1 (no division by 1 - x^2).
void fill_legendre_table(double x, int max_degree, AxisTable& out) {
    const std::size_t n = static_cast<std::size_t>(max_degree) + 1u;
    out.value.assign(n, double(0));
    out.first.assign(n, double(0));
    out.second.assign(n, double(0));
    out.value[0] = double(1);
    if (max_degree >= 1) {
        out.value[1] = x;
        out.first[1] = double(1);
    }
    for (int k = 1; k < max_degree; ++k) {
        const std::size_t kp = static_cast<std::size_t>(k);
        const double kk = static_cast<double>(k);
        const double two_k_plus_one = double(2) * kk + double(1);
        const double inv = double(1) / (kk + double(1));
        out.value[kp + 1] =
            (two_k_plus_one * x * out.value[kp] - kk * out.value[kp - 1]) * inv;
        out.first[kp + 1] =
            (two_k_plus_one * (out.value[kp] + x * out.first[kp]) -
             kk * out.first[kp - 1]) * inv;
        out.second[kp + 1] =
            (two_k_plus_one * (double(2) * out.first[kp] + x * out.second[kp]) -
             kk * out.second[kp - 1]) * inv;
    }
}

void fill_axis_table(ModalAxisKind kind, double x, int max_degree, AxisTable& out) {
    if (kind == ModalAxisKind::Legendre) {
        fill_legendre_table(x, max_degree, out);
    } else {
        fill_monomial_table(x, max_degree, out);
    }
}

// Maximum tolerated infinity-norm condition number of a serendipity interpolation
// (Vandermonde) matrix. Above this the inverse loses more than about half of
// double precision (~1/sqrt(epsilon)), so construction throws rather than return
// silently-degraded functions. With the Legendre modal basis and
// Gauss-Lobatto-Legendre nodes the condition number stays far below this across
// the recommended range (~1.7e4 at quadrilateral order 10, ~1.3e4 at hexahedral
// order 8); the bound is the numerical-soundness backstop for orders pushed well
// past it. The shape-function quality limit (Lebesgue constant) is the tighter,
// inherent constraint and is documented/tested separately.
constexpr double kSerendipityVandermondeMaxCond = double(1e8);

// Infinity norm (maximum absolute row sum) of a row-major n-by-n matrix.
double matrix_norm_inf(const std::vector<double>& matrix, std::size_t n) {
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

// Per-axis degree triples (ax, ay, az) of the serendipity mode space: every
// combination whose superlinear degree (the sum of superlinear_term over the axes;
// see the SerendipityBasis class documentation for the rule) is at most `order`.
// `active_axes` is 2 for the quadrilateral (az pinned to 0) and 3 for the
// hexahedron, so the quad space is exactly the hex space restricted to az = 0. The
// set is downward-closed (so it spans both the monomial and the tensor Legendre
// basis; see ModalAxisKind), and the resulting nodal basis is independent of how
// the set is ordered.
std::vector<std::array<int, 3>> serendipity_exponents(int order, int active_axes) {
    const int max_y = active_axes >= 2 ? order : 0;
    const int max_z = active_axes >= 3 ? order : 0;
    std::vector<std::array<int, 3>> exponents;
    for (int az = 0; az <= max_z; ++az) {
        for (int ay = 0; ay <= max_y; ++ay) {
            for (int ax = 0; ax <= order; ++ax) {
                if (superlinear_term(ax) + superlinear_term(ay) + superlinear_term(az) <= order) {
                    exponents.push_back({ax, ay, az});
                }
            }
        }
    }
    return exponents;
}

// Build the nodal coefficient table for a serendipity family: assemble the
// generalized Vandermonde V[node][mode] = phi_a(r) phi_b(s) phi_c(t) at the
// public-order reference nodes -- with phi the monomial or Legendre 1D modes per
// `kind` -- and invert it. Because the nodes are in public order, the inverse is
// already in public basis order and needs no output permutation. The same routine
// serves the quadrilateral, hexahedral, and Wedge15 spaces. `max_degree` bounds
// the per-axis mode degree (the family's order). Construction throws if the matrix
// is too ill-conditioned to trust (see kSerendipityVandermondeMaxCond).
std::vector<double> build_inverse_vandermonde(
    std::span<const Vec3> nodes,
    std::span<const std::array<int, 3>> exponents,
    const std::string& label,
    ModalAxisKind kind,
    int max_degree) {
    const std::size_t n = nodes.size();
    svmp::throw_if<BasisConstructionException>(
        n == 0 || exponents.size() != n, "SerendipityBasis: invalid serendipity interpolation setup");

    std::vector<double> vandermonde(n * n, double(0));
    AxisTable tx;
    AxisTable ty;
    AxisTable tz;
    for (std::size_t row = 0; row < n; ++row) {
        const Vec3& p = nodes[row];
        fill_axis_table(kind, p[0], max_degree, tx);
        fill_axis_table(kind, p[1], max_degree, ty);
        fill_axis_table(kind, p[2], max_degree, tz);
        for (std::size_t col = 0; col < n; ++col) {
            const auto& e = exponents[col];
            vandermonde[row * n + col] =
                tx.value[static_cast<std::size_t>(e[0])] *
                ty.value[static_cast<std::size_t>(e[1])] *
                tz.value[static_cast<std::size_t>(e[2])];
        }
    }

    // Condition-number backstop: the inverse is explicitly formed just above, so
    // this is the true infinity-norm condition number
    // cond_inf = ||V||_inf * ||V^{-1}||_inf, not an estimate. Reject orders where
    // the inverse can no longer be trusted rather than returning silently-degraded
    // shape functions; the comparison is negated so a non-finite value is rejected
    // too.
    const double norm_v = matrix_norm_inf(vandermonde, n);

    // invert_dense_matrix raises a generic FEException if the Vandermonde is
    // exactly singular (a rank-deficient pivot). For a serendipity family that
    // means the node set is not unisolvent at this order -- a construction failure
    // in basis terms -- so translate it to BasisConstructionException, presenting
    // the singular and the ill-conditioned cases (below) as one catchable type in
    // one vocabulary. The matrix was just built n-by-n from n modes, so a
    // size-mismatch FEException cannot originate here; rank deficiency is the only
    // FEException this call can raise. (Defensive: the supported node sets are
    // provably unisolvent, so this branch is not reachable for the shipped orders.)
    std::vector<double> inverse;
    try {
        inverse = math::invert_dense_matrix(
            std::move(vandermonde), n,
            "SerendipityBasis interpolation matrix for " + label);
    } catch (const FEException&) {
        svmp::raise<BasisConstructionException>("SerendipityBasis: " + label +
                " interpolation matrix is singular; the serendipity node set is not "
                "unisolvent at the requested order");
    }
    const double condition_number = norm_v * matrix_norm_inf(inverse, n);
    svmp::throw_if<BasisConstructionException>(
        !(condition_number <= kSerendipityVandermondeMaxCond), "SerendipityBasis: " + label +
            " interpolation matrix is too ill-conditioned (condition number ~ " +
            std::to_string(condition_number) +
            "); the requested order exceeds the well-conditioned range");
    return inverse;
}

// Wedge15 serendipity monomial space, as (x, y, z) exponent triples. The prism is
// the triangle cross-section (x, y) crossed with the through-axis (z): the 6
// triangle monomials x^a y^b with a + b <= 2 times the 3 line monomials z^c with
// c <= 2 form the complete 18-mode Wedge18 space. Wedge15 serendipity drops the 3
// superlinear modes -- a quadratic triangle monomial (a + b == 2) times z^2,
// i.e. (2,0,2), (1,1,2), (0,2,2) -- leaving 6*3 - 3 = 15. The set (not its order)
// fixes the space; the nodal basis is the inverse Vandermonde at the Wedge15 nodes.
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

struct NormalizedSerendipityRequest {
    BasisTopology topology;
    int dimension;
    int order;
};

int serendipity_named_order(ElementType type) {
    switch (type) {
        case ElementType::Hex8:
            return 1;
        case ElementType::Quad8:
        case ElementType::Hex20:
        case ElementType::Wedge15:
            return 2;
        default:
            return -1;
    }
}

// Validate a named serendipity element/order pairing and return its topology,
// reference dimension, and order. The named serendipity layouts (Quad8, Hex8,
// Hex20, Wedge15) are each pinned to a single polynomial order by their node
// count, so a mismatched explicit order is rejected. Arbitrary-order
// quadrilateral serendipity is not a named element: it is requested through the
// BasisTopology::Quadrilateral constructor.
NormalizedSerendipityRequest normalize_serendipity_request(ElementType type, int order) {
    // The request must supply the layout's fixed order (serendipity_named_order):
    // it is never floored or otherwise adjusted to fit, so order 0 and negative
    // orders are rejected rather than promoted to a valid layout.
    const int expected_order = serendipity_named_order(type);
    switch (type) {
        case ElementType::Quad8:
            svmp::throw_if<BasisConfigurationException>(order != expected_order, "SerendipityBasis: Quad8 is the quadratic 8-node serendipity layout (order 2 only); "
                "use BasisTopology::Quadrilateral for higher-order quadrilateral serendipity");
            return {BasisTopology::Quadrilateral, 2, expected_order};
        case ElementType::Hex8:
            svmp::throw_if<BasisConfigurationException>(order != expected_order, "SerendipityBasis: Hex8 is the trilinear 8-node basis (order 1 only); use Hex20 for quadratic serendipity");
            return {BasisTopology::Hexahedron, 3, expected_order};
        case ElementType::Hex20:
            svmp::throw_if<BasisConfigurationException>(order != expected_order, "SerendipityBasis: Hex20 is the 20-node quadratic serendipity layout (order 2 only)");
            return {BasisTopology::Hexahedron, 3, expected_order};
        case ElementType::Wedge15:
            svmp::throw_if<BasisConfigurationException>(order != expected_order, "SerendipityBasis: Wedge15 is the 15-node quadratic serendipity layout (order 2 only)");
            return {BasisTopology::Wedge, 3, expected_order};
        default:
            svmp::raise<BasisElementCompatibilityException>("SerendipityBasis named elements are Quad8, Hex8, Hex20, and Wedge15; "
                "use BasisTopology::Quadrilateral for arbitrary-order quadrilateral serendipity");
    }
}

} // namespace

SerendipityBasis::SerendipityBasis(BasisTopology topology, int order)
    : topology_(topology) {
    const bool supported_topology = topology_ == BasisTopology::Quadrilateral ||
                                    topology_ == BasisTopology::Hexahedron;
    svmp::throw_if<BasisElementCompatibilityException>(
        !supported_topology, "SerendipityBasis: arbitrary-order topology construction is supported for "
        "Quadrilateral and Hexahedron; use the named ElementType (Wedge15) for wedge serendipity");
    svmp::throw_if<BasisConfigurationException>(
        order < 1, "SerendipityBasis: serendipity requires a polynomial order >= 1");
    dimension_ = topology_ == BasisTopology::Hexahedron ? 3 : 2;
    order_ = order;
    if (topology_ == BasisTopology::Hexahedron) {
        init_hexahedron(order_);
    } else {
        init_quadrilateral(order_);
    }
}

SerendipityBasis::SerendipityBasis(ElementType type, int order) {
    const NormalizedSerendipityRequest normalized = normalize_serendipity_request(type, order);
    topology_ = normalized.topology;
    dimension_ = normalized.dimension;
    order_ = normalized.order;

    switch (type) {
        case ElementType::Quad8:
            // Quad8 is the order-2 instance of the quadrilateral serendipity
            // space; the named overload only pins the order.
            init_quadrilateral(order_);
            return;
        case ElementType::Hex8:
            // Hex8 is the order-1 instance of the hexahedral serendipity space.
            init_hexahedron(1);
            return;
        case ElementType::Hex20:
            // Hex20 is the order-2 instance of the hexahedral serendipity space.
            init_hexahedron(2);
            return;
        case ElementType::Wedge15:
            init_fixed_named(type);
            return;
        default:
            // normalize_serendipity_request already rejected anything else.
            svmp::raise<BasisElementCompatibilityException>("SerendipityBasis: unsupported named serendipity element");
    }
}

SerendipityBasis::SerendipityBasis(ElementType type)
    : SerendipityBasis(type, serendipity_named_order(type)) {}

// Build the quadrilateral serendipity mode set, reference nodes, and nodal
// coefficient table for the given order. The coefficient table is the inverse
// Vandermonde of tensor Legendre modes spanning the same polynomial space as the
// monomial degree triples; because the nodes are in public order, evaluation
// needs no output permutation. Reference nodes come from the single
// ReferenceNodeLayout serendipity generator for both the named Quad8 layout and
// the arbitrary-order path.
void SerendipityBasis::init_quadrilateral(int order) {
    mode_exponents_ = serendipity_exponents(order, /*active_axes=*/2);
    size_ = mode_exponents_.size();
    nodes_ = ReferenceNodeLayout::serendipity_node_coords(BasisTopology::Quadrilateral, order);
    svmp::throw_if<BasisConstructionException>(
        nodes_.size() != size_, "SerendipityBasis: quadrilateral serendipity setup produced inconsistent sizes");
    uses_legendre_modes_ = true;
    inv_vandermonde_ = build_inverse_vandermonde(
        nodes_, mode_exponents_, "Quad order " + std::to_string(order),
        ModalAxisKind::Legendre, order);
}

// Build the hexahedral serendipity mode set, reference nodes, and nodal
// coefficient table for the given order, mirroring init_quadrilateral. Reference
// nodes come from the single ReferenceNodeLayout serendipity generator; Hex8
// (order 1) and Hex20 (order 2) are its order-1/order-2 instances and match the
// public Hex8/Hex20 ordering exactly.
void SerendipityBasis::init_hexahedron(int order) {
    mode_exponents_ = serendipity_exponents(order, /*active_axes=*/3);
    size_ = mode_exponents_.size();
    nodes_ = ReferenceNodeLayout::serendipity_node_coords(BasisTopology::Hexahedron, order);
    svmp::throw_if<BasisConstructionException>(
        nodes_.size() != size_, "SerendipityBasis: hexahedral serendipity setup produced inconsistent sizes");
    uses_legendre_modes_ = true;
    inv_vandermonde_ = build_inverse_vandermonde(
        nodes_, mode_exponents_, "Hex order " + std::to_string(order),
        ModalAxisKind::Legendre, order);
}

// Build the Wedge15 serendipity layout from its tabulated monomial space and
// public-order ReferenceNodeLayout nodes. Hexahedral serendipity (Hex8 and Hex20
// included) is generated by init_hexahedron, so the prism is the only named
// layout that still carries a fixed monomial table.
void SerendipityBasis::init_fixed_named(ElementType type) {
    svmp::throw_if<BasisConstructionException>(
        type != ElementType::Wedge15, "SerendipityBasis: init_fixed_named builds only the Wedge15 layout");
    size_ = 15u;
    const std::span<const std::array<int, 3>> family_exponents(
        kWedge15MonomialExponents.data(), kWedge15MonomialExponents.size());
    nodes_ = ReferenceNodeLayout::node_coords(type);
    svmp::throw_if<BasisConstructionException>(
        nodes_.size() != size_, "SerendipityBasis: Wedge15 layout node count does not match basis size");
    svmp::throw_if<BasisConstructionException>(
        family_exponents.size() != size_, "SerendipityBasis: Wedge15 monomial count does not match basis size");
    mode_exponents_.assign(family_exponents.begin(), family_exponents.end());
    // Wedge15 is the fixed order-2 layout; its 15x15 system is trivially
    // well-conditioned, so it keeps the monomial modal basis.
    uses_legendre_modes_ = false;
    inv_vandermonde_ = build_inverse_vandermonde(
        nodes_, mode_exponents_, "Wedge15", ModalAxisKind::Monomial, order_);
}

void SerendipityBasis::evaluate_all_to(const math::Vector<double, 3>& xi,
                                       std::span<double> values_out,
                                       std::span<Gradient> gradients_out,
                                       std::span<Hessian> hessians_out) const {
    // Private sink: callers guarantee valid output spans -- the public *_to methods
    // validate their one output with require_span_size, and the vector evaluators
    // resize to size_. An empty span here means "skip that quantity".

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
        mode_exponents_.size() != size_ ||
            inv_vandermonde_.size() != size_ * size_,
        "SerendipityBasis: interpolation tables are not initialized for evaluation");

    // Build the per-axis modal tables once, then accumulate over the modes. The
    // mode family must match the one the coefficient table was built with.
    const ModalAxisKind kind =
        uses_legendre_modes_ ? ModalAxisKind::Legendre : ModalAxisKind::Monomial;
    AxisTable tx;
    AxisTable ty;
    AxisTable tz;
    fill_axis_table(kind, x, order_, tx);
    fill_axis_table(kind, y, order_, ty);
    fill_axis_table(kind, z, order_, tz);

    // Accumulate the nodal shape functions from the modal tables. For each mode j,
    // phi = phi_a(r) phi_b(s) phi_c(t) (and its derivatives) is weighted by the
    // inverse-Vandermonde coefficient for each basis slot; the table is already in
    // public basis order, so slot i reads column i directly. The spans were zeroed
    // above and an empty span is skipped.
    const bool want_values = !values_out.empty();
    const bool want_gradients = !gradients_out.empty();
    const bool want_hessians = !hessians_out.empty();

    for (std::size_t j = 0; j < size_; ++j) {
        const std::array<int, 3>& e = mode_exponents_[j];
        const std::size_t ex = static_cast<std::size_t>(e[0]);
        const std::size_t ey = static_cast<std::size_t>(e[1]);
        const std::size_t ez = static_cast<std::size_t>(e[2]);

        const double vx = tx.value[ex];
        const double vy = ty.value[ey];
        const double vz = tz.value[ez];
        const double phi = vx * vy * vz;

        double d_dr = double(0), d_ds = double(0), d_dt = double(0);
        if (want_gradients || want_hessians) {
            d_dr = tx.first[ex] * vy * vz;
            d_ds = vx * ty.first[ey] * vz;
            d_dt = vx * vy * tz.first[ez];
        }

        double d_drr = double(0), d_dss = double(0), d_dtt = double(0);
        double d_drs = double(0), d_drt = double(0), d_dst = double(0);
        if (want_hessians) {
            d_drr = tx.second[ex] * vy * vz;
            d_dss = vx * ty.second[ey] * vz;
            d_dtt = vx * vy * tz.second[ez];
            d_drs = tx.first[ex] * ty.first[ey] * vz;
            d_drt = tx.first[ex] * vy * tz.first[ez];
            d_dst = vx * ty.first[ey] * tz.first[ez];
        }

        for (std::size_t slot = 0; slot < size_; ++slot) {
            const double c = inv_vandermonde_[j * size_ + slot];
            if (want_values) {
                values_out[slot] += c * phi;
            }
            if (want_gradients) {
                Gradient& g = gradients_out[slot];
                g[0] += c * d_dr;
                g[1] += c * d_ds;
                g[2] += c * d_dt;
            }
            if (want_hessians) {
                Hessian& h = hessians_out[slot];
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

} // namespace svmp::FE::basis
