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

std::vector<double> invert_dense_matrix(std::vector<double> matrix, int n, const char* label) {
    return math::invert_dense_matrix(
        std::move(matrix),
        static_cast<std::size_t>(n),
        std::string("SerendipityBasis interpolation matrix for ") + label);
}

std::vector<double> quad_serendipity_inverse_vandermonde(
    std::span<const Vec3> nodes,
    std::span<const std::array<int, 2>> exponents,
    int order) {
    const int n = static_cast<int>(nodes.size());
    svmp::throw_if<BasisConstructionException>(
        n == 0 || exponents.size() != nodes.size(), SVMP_HERE,
        "SerendipityBasis: invalid quadrilateral serendipity interpolation setup");

    std::vector<double> vandermonde(static_cast<std::size_t>(n * n), double(0));
    auto idx = [n](int row, int col) -> std::size_t {
        return static_cast<std::size_t>(row * n + col);
    };

    for (int row = 0; row < n; ++row) {
        const double x = nodes[static_cast<std::size_t>(row)][0];
        const double y = nodes[static_cast<std::size_t>(row)][1];
        for (int col = 0; col < n; ++col) {
            const auto [ax, ay] = exponents[static_cast<std::size_t>(col)];
            vandermonde[idx(row, col)] = integer_power(x, ax) * integer_power(y, ay);
        }
    }

    // Quadrilateral serendipity bases are generated from the requested
    // monomial space, so a small dense inverse produces the nodal coefficient
    // table at construction time. Hex20 and Wedge15 use fixed tables because
    // only their quadratic layouts are supported here.
    const std::string label = "Quad order " + std::to_string(order);
    return invert_dense_matrix(std::move(vandermonde), n, label.c_str());
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

// Coefficients for the quadratic Wedge15 nodal serendipity basis. Rows are
// monomials in kWedge15MonomialExponents order; columns are basis functions in
// public Wedge15 node order. The table is the inverse of
// V[node][monomial] = r^a s^b t^c evaluated at ReferenceNodeLayout Wedge15
// nodes, so V * kWedge15Coefficients is the identity.
constexpr std::array<std::array<double, 15>, 15> kWedge15Coefficients = {{
    {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}},
    {{-0.5, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
    {{0.5, -0, -0, 0.5, -0, -0, -0, -0, -0, -0, -0, -0, -1, -0, -0}},
    {{-1, 0, -1, -1, 0, -1, 0, 0, 2, 0, 0, 2, -1, 0, 1}},
    {{1.5, 0, 0.5, -1.5, 0, -0.5, 0, 0, -2, 0, 0, 2, 0, 0, 0}},
    {{-0.5, -0, 0.5, -0.5, -0, 0.5, -0, -0, -0, -0, -0, -0, 1, -0, -1}},
    {{1, 0, 1, 1, 0, 1, 0, 0, -2, 0, 0, -2, 0, 0, 0}},
    {{-1, 0, -1, 1, 0, 1, 0, 0, 2, 0, 0, -2, 0, 0, 0}},
    {{-1, -1, 0, -1, -1, 0, 2, 0, 0, 2, 0, 0, -1, 1, 0}},
    {{1.5, 0.5, 0, -1.5, -0.5, 0, -2, 0, 0, 2, 0, 0, 0, 0, 0}},
    {{-0.5, 0.5, -0, -0.5, 0.5, -0, -0, -0, -0, -0, -0, -0, 1, -1, -0}},
    {{2, 0, -0, 2, 0, -0, -2, 2, -2, -2, 2, -2, -0, -0, -0}},
    {{-2, 0, 0, 2, 0, 0, 2, -2, 2, -2, 2, -2, 0, 0, 0}},
    {{1, 1, -0, 1, 1, -0, -2, -0, -0, -2, -0, -0, -0, -0, -0}},
    {{-1, -1, -0, 1, 1, -0, 2, -0, -0, -2, -0, -0, -0, -0, -0}}
}};

constexpr std::array<std::array<int, 3>, 20> kHex20MonomialExponents = {{
    {{0, 0, 0}}, {{0, 0, 1}}, {{0, 0, 2}}, {{0, 1, 0}}, {{0, 1, 1}},
    {{0, 1, 2}}, {{0, 2, 0}}, {{0, 2, 1}}, {{1, 0, 0}}, {{1, 0, 1}},
    {{1, 0, 2}}, {{1, 1, 0}}, {{1, 1, 1}}, {{1, 1, 2}}, {{1, 2, 0}},
    {{1, 2, 1}}, {{2, 0, 0}}, {{2, 0, 1}}, {{2, 1, 0}}, {{2, 1, 1}}
}};

// Coefficients for the quadratic Hex20 nodal serendipity basis. Rows are
// monomials in kHex20MonomialExponents order; columns are basis functions in
// the internal Hex20 coefficient-table order. The table is the inverse of
// V[node][monomial] = r^a s^b t^c evaluated at the corresponding Hex20
// reference nodes, so V * kHex20Coefficients is the identity. Evaluation
// remaps public output slots through ReferenceNodeLayout::mesh_to_basis_ordering.
constexpr std::array<std::array<double, 20>, 20> kHex20Coefficients = {{
    {{-0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25}},
    {{0.125, 0.125, 0.125, 0.125, -0.125, -0.125, -0.125, -0.125, -0.25, 0.25, -0.25, 0.25, -0.25, -0.25, 0.25, 0.25, 0, 0, 0, 0}},
    {{0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0, 0, 0, 0, 0, 0, 0, 0, -0.25, -0.25, -0.25, -0.25}},
    {{0.125, 0.125, -0.125, -0.125, 0.125, 0.125, -0.125, -0.125, -0.25, -0.25, 0.25, 0.25, 0, 0, 0, 0, -0.25, -0.25, 0.25, 0.25}},
    {{0, 0, 0, 0, 0, 0, 0, 0, 0.25, -0.25, -0.25, 0.25, 0, 0, 0, 0, 0, 0, 0, 0}},
    {{-0.125, -0.125, 0.125, 0.125, -0.125, -0.125, 0.125, 0.125, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, 0.25, -0.25, -0.25}},
    {{0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0, 0, 0, 0, -0.25, -0.25, -0.25, -0.25, 0, 0, 0, 0}},
    {{-0.125, -0.125, -0.125, -0.125, 0.125, 0.125, 0.125, 0.125, 0, 0, 0, 0, 0.25, 0.25, -0.25, -0.25, 0, 0, 0, 0}},
    {{0.125, -0.125, -0.125, 0.125, 0.125, -0.125, -0.125, 0.125, 0, 0, 0, 0, -0.25, 0.25, -0.25, 0.25, -0.25, 0.25, -0.25, 0.25}},
    {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, -0.25, -0.25, 0.25, 0, 0, 0, 0}},
    {{-0.125, 0.125, 0.125, -0.125, -0.125, 0.125, 0.125, -0.125, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, -0.25, 0.25, -0.25}},
    {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, -0.25, -0.25, 0.25}},
    {{-0.125, 0.125, -0.125, 0.125, 0.125, -0.125, 0.125, -0.125, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
    {{0.125, -0.125, 0.125, -0.125, 0.125, -0.125, 0.125, -0.125, 0, 0, 0, 0, 0, 0, 0, 0, -0.25, 0.25, 0.25, -0.25}},
    {{-0.125, 0.125, 0.125, -0.125, -0.125, 0.125, 0.125, -0.125, 0, 0, 0, 0, 0.25, -0.25, 0.25, -0.25, 0, 0, 0, 0}},
    {{0.125, -0.125, -0.125, 0.125, -0.125, 0.125, 0.125, -0.125, 0, 0, 0, 0, -0.25, 0.25, 0.25, -0.25, 0, 0, 0, 0}},
    {{0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, -0.25, -0.25, -0.25, -0.25, 0, 0, 0, 0, 0, 0, 0, 0}},
    {{-0.125, -0.125, -0.125, -0.125, 0.125, 0.125, 0.125, 0.125, 0.25, -0.25, 0.25, -0.25, 0, 0, 0, 0, 0, 0, 0, 0}},
    {{-0.125, -0.125, 0.125, 0.125, -0.125, -0.125, 0.125, 0.125, 0.25, 0.25, -0.25, -0.25, 0, 0, 0, 0, 0, 0, 0, 0}},
    {{0.125, 0.125, -0.125, -0.125, -0.125, -0.125, 0.125, 0.125, -0.25, 0.25, 0.25, -0.25, 0, 0, 0, 0, 0, 0, 0, 0}}
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
// functions (the coefficient table is square). Outputs are assumed pre-zeroed by
// the caller; an empty span skips that quantity.
//
// `table_to_output_order` maps each output slot to the basis column of the
// coefficient table. An empty span means the table is already in output (public
// node) order, i.e. the identity permutation: Hex20 supplies a real permutation
// because its table is authored in an internal node order, while Wedge15 and the
// quadrilateral serendipity tables are authored directly in public order.
template <typename ExponentFn, typename CoeffFn>
void eval_monomial_basis(double r, double s, double t,
                         std::size_t count,
                         ExponentFn&& exponent,
                         CoeffFn&& coeff,
                         std::span<const std::size_t> table_to_output_order,
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
            const std::size_t basis_index =
                table_to_output_order.empty() ? slot : table_to_output_order[slot];
            const double c = coeff(j, basis_index);
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

} // namespace

SerendipityBasis::SerendipityBasis(ElementType type, int order)
    : element_type_(type), dimension_(0), order_(order), size_(0) {
    if (type == ElementType::Quad4 || type == ElementType::Quad8) {
        dimension_ = 2;
        if (order_ < 1) {
            order_ = 1;
        }
        svmp::throw_if<BasisConfigurationException>(
            type == ElementType::Quad8 && order_ != 2, SVMP_HERE,
            "SerendipityBasis: Quad8 is only valid for quadratic order 2; use Quad4 for higher-order quadrilateral serendipity");
        quad_monomial_exponents_ = quad_serendipity_exponents(order_);
        size_ = quad_monomial_exponents_.size();
        nodes_ = quad_serendipity_nodes(order_, size_);
        svmp::throw_if<BasisConstructionException>(
            nodes_.size() != size_, SVMP_HERE,
            "SerendipityBasis: quadrilateral serendipity setup produced inconsistent sizes");
        quad_inv_vandermonde_ = quad_serendipity_inverse_vandermonde(nodes_, quad_monomial_exponents_, order_);
    } else if (type == ElementType::Hex8 || type == ElementType::Hex20) {
        dimension_ = 3;
        if (order_ < 1) {
            order_ = 1;
        }
        svmp::throw_if<BasisConfigurationException>(
            type == ElementType::Hex8 && order_ != 1, SVMP_HERE,
            "SerendipityBasis: Hex8 is the trilinear 8-node basis (order 1 only); use Hex20 for quadratic serendipity");
        svmp::throw_if<BasisConfigurationException>(
            type == ElementType::Hex20 && order_ != 2, SVMP_HERE,
            "SerendipityBasis: Hex20 is the 20-node quadratic serendipity layout (order 2 only)");
        size_ = (type == ElementType::Hex8) ? 8u : 20u;
    } else if (type == ElementType::Wedge15) {
        dimension_ = 3;
        if (order_ < 2) {
            order_ = 2;
        }
        if (order_ == 2) {
            size_ = 15;
        } else {
            svmp::raise<BasisConfigurationException>(SVMP_HERE,
                "SerendipityBasis supports up to quadratic on wedge15");
        }
    } else {
        svmp::raise<BasisElementCompatibilityException>(SVMP_HERE,
            "SerendipityBasis supports Quad4/Quad8, Hex8/Hex20, and Wedge15 elements");
    }

    if (nodes_.empty()) {
        nodes_.reserve(size_);
        for (std::size_t i = 0; i < size_; ++i) {
            nodes_.push_back(ReferenceNodeLayout::get_node_coords(element_type_, i));
        }
    }
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

    if (dimension_ == 2) {
        svmp::throw_if<BasisEvaluationException>(
            quad_monomial_exponents_.size() != size_ ||
                quad_inv_vandermonde_.size() != size_ * size_,
            SVMP_HERE,
            "SerendipityBasis: quadrilateral interpolation tables are not initialized for value evaluation");

        // Quadrilateral serendipity monomials are planar; the through-axis
        // exponent is zero, so all z derivatives vanish. The inverse Vandermonde
        // is already in public node order (identity output ordering).
        eval_monomial_basis(
            x, y, z, size_,
            [this](std::size_t j) {
                const auto& e = quad_monomial_exponents_[j];
                return std::array<int, 3>{e[0], e[1], 0};
            },
            [this](std::size_t j, std::size_t i) {
                return quad_inv_vandermonde_[j * size_ + i];
            },
            std::span<const std::size_t>{},
            values_out, gradients_out, hessians_out);
        return;
    }

    if (dimension_ == 3 && order_ == 1) {
        evaluate_hex8_reference(x, y, z, values_out, gradients_out, hessians_out);
        return;
    }

    if (element_type_ == ElementType::Hex20) {
        // The Hex20 coefficient table is authored in an internal node order, so
        // results are remapped into the public node layout through mesh_to_basis.
        const auto mesh_to_basis = ReferenceNodeLayout::mesh_to_basis_ordering(element_type_);
        svmp::throw_if<BasisEvaluationException>(mesh_to_basis.size() != size_, SVMP_HERE,
                                               "Hex20 mesh-to-basis ordering is not registered");
        eval_monomial_basis(
            x, y, z, size_,
            [](std::size_t j) { return kHex20MonomialExponents[j]; },
            [](std::size_t j, std::size_t i) { return kHex20Coefficients[j][i]; },
            mesh_to_basis,
            values_out, gradients_out, hessians_out);
        return;
    }

    if (element_type_ == ElementType::Wedge15) {
        // The Wedge15 coefficient table is authored directly in public node
        // order, so no output reordering is applied (identity permutation).
        eval_monomial_basis(
            x, y, z, size_,
            [](std::size_t j) { return kWedge15MonomialExponents[j]; },
            [](std::size_t j, std::size_t i) { return kWedge15Coefficients[j][i]; },
            std::span<const std::size_t>{},
            values_out, gradients_out, hessians_out);
        return;
    }

    svmp::raise<BasisEvaluationException>(SVMP_HERE,
        "SerendipityBasis::evaluate_all_to: unsupported serendipity configuration");
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
