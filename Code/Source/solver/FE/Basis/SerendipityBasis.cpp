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

namespace svmp {
namespace FE {
namespace basis {

namespace {
using Vec3 = math::Vector<Real, 3>;

void evaluate_hex8_reference(Real r,
                             Real s,
                             Real t,
                             std::span<Real> values,
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
        const Real a = Real(signs[i][0]);
        const Real b = Real(signs[i][1]);
        const Real c = Real(signs[i][2]);
        const Real ar = Real(1) + a * r;
        const Real bs = Real(1) + b * s;
        const Real ct = Real(1) + c * t;

        if (!values.empty()) {
            values[i] = Real(0.125) * ar * bs * ct;
        }
        if (!gradients.empty()) {
            Gradient& g = gradients[i];
            g[0] = Real(0.125) * a * bs * ct;
            g[1] = Real(0.125) * b * ar * ct;
            g[2] = Real(0.125) * c * ar * bs;
        }
        if (!hessians.empty()) {
            Hessian& h = hessians[i];
            h(0, 0) = Real(0);
            h(0, 1) = Real(0.125) * a * b * ct;
            h(0, 2) = Real(0.125) * a * c * bs;
            h(1, 0) = h(0, 1);
            h(1, 1) = Real(0);
            h(1, 2) = Real(0.125) * b * c * ar;
            h(2, 0) = h(0, 2);
            h(2, 1) = h(1, 2);
            h(2, 2) = Real(0);
        }
    }
}

int quad_serendipity_superlinear_degree(int ax, int ay) {
    return (ax > 1 ? ax : 0) + (ay > 1 ? ay : 0);
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

std::vector<Vec3> quad_serendipity_nodes(int order, std::size_t total_size) {
    std::vector<Vec3> nodes;
    if (order <= 0) {
        return nodes;
    }

    const Real inv_order = Real(1) / Real(order);

    nodes.push_back(Vec3{Real(-1), Real(-1), Real(0)});
    nodes.push_back(Vec3{Real(1),  Real(-1), Real(0)});
    nodes.push_back(Vec3{Real(1),  Real(1),  Real(0)});
    nodes.push_back(Vec3{Real(-1), Real(1),  Real(0)});

    for (int i = 1; i < order; ++i) {
        nodes.push_back(Vec3{Real(-1) + Real(2 * i) * inv_order, Real(-1), Real(0)});
    }
    for (int i = 1; i < order; ++i) {
        nodes.push_back(Vec3{Real(1), Real(-1) + Real(2 * i) * inv_order, Real(0)});
    }
    for (int i = 1; i < order; ++i) {
        nodes.push_back(Vec3{Real(1) - Real(2 * i) * inv_order, Real(1), Real(0)});
    }
    for (int i = 1; i < order; ++i) {
        nodes.push_back(Vec3{Real(-1), Real(1) - Real(2 * i) * inv_order, Real(0)});
    }

    FE::throw_if<BasisConstructionException>(
        nodes.size() > total_size, SVMP_HERE,
        "SerendipityBasis: quadrilateral serendipity boundary nodes exceed requested size");

    const std::size_t interior_count = total_size - nodes.size();
    if (interior_count == 0u) {
        return nodes;
    }

    std::vector<Vec3> interior_candidates;
    interior_candidates.reserve(static_cast<std::size_t>((order - 1) * (order - 1)));
    for (int j = 1; j < order; ++j) {
        for (int i = 1; i < order; ++i) {
            interior_candidates.push_back(
                Vec3{Real(-1) + Real(2 * i) * inv_order,
                     Real(-1) + Real(2 * j) * inv_order,
                     Real(0)});
        }
    }

    std::sort(interior_candidates.begin(), interior_candidates.end(),
              [](const Vec3& a, const Vec3& b) {
                  const Real a_linf = std::max(std::abs(a[0]), std::abs(a[1]));
                  const Real b_linf = std::max(std::abs(b[0]), std::abs(b[1]));
                  if (a_linf != b_linf) {
                      return a_linf < b_linf;
                  }

                  const Real a_l1 = std::abs(a[0]) + std::abs(a[1]);
                  const Real b_l1 = std::abs(b[0]) + std::abs(b[1]);
                  if (a_l1 != b_l1) {
                      return a_l1 < b_l1;
                  }

                  if (a[1] != b[1]) {
                      return a[1] < b[1];
                  }
                  return a[0] < b[0];
              });

    FE::throw_if<BasisConstructionException>(
        interior_count > interior_candidates.size(), SVMP_HERE,
        "SerendipityBasis: insufficient quadrilateral interior nodes for requested serendipity order");

    nodes.insert(nodes.end(),
                 interior_candidates.begin(),
                 interior_candidates.begin() + static_cast<std::ptrdiff_t>(interior_count));
    return nodes;
}

std::vector<Real> invert_dense_matrix(std::vector<Real> matrix, int n, const char* label) {
    return math::invert_dense_matrix(
        std::move(matrix),
        static_cast<std::size_t>(n),
        std::string("SerendipityBasis interpolation matrix for ") + label);
}

std::vector<Real> quad_serendipity_inverse_vandermonde(
    std::span<const Vec3> nodes,
    std::span<const std::array<int, 2>> exponents,
    int order) {
    const int n = static_cast<int>(nodes.size());
    FE::throw_if<BasisConstructionException>(
        n == 0 || exponents.size() != nodes.size(), SVMP_HERE,
        "SerendipityBasis: invalid quadrilateral serendipity interpolation setup");

    std::vector<Real> vandermonde(static_cast<std::size_t>(n * n), Real(0));
    auto idx = [n](int row, int col) -> std::size_t {
        return static_cast<std::size_t>(row * n + col);
    };

    for (int row = 0; row < n; ++row) {
        const Real x = nodes[static_cast<std::size_t>(row)][0];
        const Real y = nodes[static_cast<std::size_t>(row)][1];
        for (int col = 0; col < n; ++col) {
            const auto [ax, ay] = exponents[static_cast<std::size_t>(col)];
            vandermonde[idx(row, col)] = std::pow(x, ax) * std::pow(y, ay);
        }
    }

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

constexpr std::array<std::array<Real, 15>, 15> kWedge15Coefficients = {{
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

constexpr std::array<std::array<Real, 20>, 20> kHex20Coefficients = {{
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
    Real value;   ///< x^a
    Real first;   ///< d/dx (x^a)     = a x^(a-1)
    Real second;  ///< d^2/dx^2 (x^a) = a (a-1) x^(a-2)
};

inline Real integer_power(Real base, int exponent) {
    Real result = Real(1);
    for (int k = 0; k < exponent; ++k) {
        result *= base;
    }
    return result;
}

inline MonomialAxis monomial_axis(Real x, int exponent) {
    MonomialAxis axis;
    axis.value = integer_power(x, exponent);
    axis.first = (exponent > 0) ? Real(exponent) * integer_power(x, exponent - 1) : Real(0);
    axis.second = (exponent > 1)
                      ? Real(exponent * (exponent - 1)) * integer_power(x, exponent - 2)
                      : Real(0);
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
void eval_monomial_basis(Real r, Real s, Real t,
                         std::size_t count,
                         ExponentFn&& exponent,
                         CoeffFn&& coeff,
                         std::span<const std::size_t> table_to_output_order,
                         std::span<Real> values,
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

        const Real phi = ax.value * ay.value * az.value;

        Real d_dr = Real(0), d_ds = Real(0), d_dt = Real(0);
        if (want_gradients || want_hessians) {
            d_dr = ax.first * ay.value * az.value;
            d_ds = ax.value * ay.first * az.value;
            d_dt = ax.value * ay.value * az.first;
        }

        Real d_drr = Real(0), d_dss = Real(0), d_dtt = Real(0);
        Real d_drs = Real(0), d_drt = Real(0), d_dst = Real(0);
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
            const Real c = coeff(j, basis_index);
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
        FE::throw_if<BasisConfigurationException>(
            type == ElementType::Quad8 && order_ != 2, SVMP_HERE,
            "SerendipityBasis: Quad8 is only valid for quadratic order 2; use Quad4 for higher-order quadrilateral serendipity");
        quad_monomial_exponents_ = quad_serendipity_exponents(order_);
        size_ = quad_monomial_exponents_.size();
        nodes_ = quad_serendipity_nodes(order_, size_);
        FE::throw_if<BasisConstructionException>(
            nodes_.size() != size_, SVMP_HERE,
            "SerendipityBasis: quadrilateral serendipity setup produced inconsistent sizes");
        quad_inv_vandermonde_ = quad_serendipity_inverse_vandermonde(nodes_, quad_monomial_exponents_, order_);
    } else if (type == ElementType::Hex8 || type == ElementType::Hex20) {
        dimension_ = 3;
        if (order_ < 1) order_ = 1;
        if (order_ == 1) {
            size_ = 8;
        } else if (order_ == 2) {
            size_ = 20;
        } else {
            FE::raise<BasisConfigurationException>(SVMP_HERE,
                "SerendipityBasis supports up to quadratic on hexahedra");
        }
    } else if (type == ElementType::Wedge15) {
        dimension_ = 3;
        if (order_ < 2) {
            order_ = 2;
        }
        if (order_ == 2) {
            size_ = 15;
        } else {
            FE::raise<BasisConfigurationException>(SVMP_HERE,
                "SerendipityBasis supports up to quadratic on wedge15");
        }
    } else {
        FE::raise<BasisElementCompatibilityException>(SVMP_HERE,
            "SerendipityBasis supports Quad4/Quad8, Hex8/Hex20, and Wedge15 elements");
    }

    if (nodes_.empty()) {
        nodes_.reserve(size_);
        for (std::size_t i = 0; i < size_; ++i) {
            nodes_.push_back(ReferenceNodeLayout::get_node_coords(element_type_, i));
        }
    }
}

void SerendipityBasis::evaluate_all_to(const math::Vector<Real, 3>& xi,
                                       std::span<Real> values_out,
                                       std::span<Gradient> gradients_out,
                                       std::span<Hessian> hessians_out) const {
    require_requested_span_size(values_out, size_, "SerendipityBasis::evaluate_all_to values");
    require_requested_span_size(gradients_out, size_, "SerendipityBasis::evaluate_all_to gradients");
    require_requested_span_size(hessians_out, size_, "SerendipityBasis::evaluate_all_to hessians");

    if (values_out.empty() && gradients_out.empty() && hessians_out.empty()) {
        return;
    }

    if (!values_out.empty()) {
        std::fill(values_out.begin(), values_out.end(), Real(0));
    }
    if (!gradients_out.empty()) {
        std::fill(gradients_out.begin(), gradients_out.end(), Gradient::Zero());
    }
    if (!hessians_out.empty()) {
        std::fill(hessians_out.begin(), hessians_out.end(), Hessian::Zero());
    }

    const Real x = xi[0];
    const Real y = xi[1];
    const Real z = xi[2];

    if (dimension_ == 2) {
        FE::throw_if<BasisEvaluationException>(
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
        FE::throw_if<BasisEvaluationException>(mesh_to_basis.size() != size_, SVMP_HERE,
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

    FE::raise<BasisEvaluationException>(SVMP_HERE,
        "SerendipityBasis::evaluate_all_to: unsupported serendipity configuration");
}

void SerendipityBasis::evaluate_values(const math::Vector<Real, 3>& xi,
                                       std::vector<Real>& values) const {
    values.resize(size_);
    evaluate_values_to(xi, std::span<Real>(values.data(), values.size()));
}

void SerendipityBasis::evaluate_gradients(const math::Vector<Real, 3>& xi,
                                          std::vector<Gradient>& gradients) const {
    gradients.resize(size_);
    evaluate_gradients_to(xi, std::span<Gradient>(gradients.data(), gradients.size()));
}

void SerendipityBasis::evaluate_hessians(const math::Vector<Real, 3>& xi,
                                         std::vector<Hessian>& hessians) const {
    hessians.resize(size_);
    evaluate_hessians_to(xi, std::span<Hessian>(hessians.data(), hessians.size()));
}

void SerendipityBasis::evaluate_all(const math::Vector<Real, 3>& xi,
                                    std::vector<Real>& values,
                                    std::vector<Gradient>& gradients,
                                    std::vector<Hessian>& hessians) const {
    values.resize(size_);
    gradients.resize(size_);
    hessians.resize(size_);
    evaluate_all_to(xi,
                    std::span<Real>(values.data(), values.size()),
                    std::span<Gradient>(gradients.data(), gradients.size()),
                    std::span<Hessian>(hessians.data(), hessians.size()));
}

void SerendipityBasis::evaluate_values_to(const math::Vector<Real, 3>& xi,
                                          std::span<Real> values_out) const {
    require_span_size(values_out.size(), size_, "SerendipityBasis::evaluate_values_to");
    evaluate_all_to(xi, values_out, std::span<Gradient>{}, std::span<Hessian>{});
}

void SerendipityBasis::evaluate_gradients_to(const math::Vector<Real, 3>& xi,
                                             std::span<Gradient> gradients_out) const {
    require_span_size(gradients_out.size(), size_, "SerendipityBasis::evaluate_gradients_to");
    evaluate_all_to(xi, std::span<Real>{}, gradients_out, std::span<Hessian>{});
}

void SerendipityBasis::evaluate_hessians_to(const math::Vector<Real, 3>& xi,
                                            std::span<Hessian> hessians_out) const {
    require_span_size(hessians_out.size(), size_, "SerendipityBasis::evaluate_hessians_to");
    evaluate_all_to(xi, std::span<Real>{}, std::span<Gradient>{}, hessians_out);
}

} // namespace basis
} // namespace FE
} // namespace svmp
