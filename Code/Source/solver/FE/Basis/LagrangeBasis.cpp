// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "LagrangeBasis.h"
#include "NodeOrderingConventions.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <span>
#include <string>

namespace svmp {
namespace FE {
namespace basis {

namespace {

using Vec3 = math::Vector<double, 3>;

struct AxisEval {
    std::vector<double> value;
    std::vector<double> first;
    std::vector<double> second;
};

struct SimplexEval {
    std::vector<double> value;
    std::vector<Gradient> gradient;
    std::vector<Hessian> hessian;
};

struct NormalizedLagrangeRequest {
    ElementType element_type;
    int order;
};

// Validate and return the supported basis topology for a Lagrange element type.
BasisTopology supported_lagrange_topology(ElementType type) {
    const BasisTopology top = topology(type);
    svmp::throw_if<BasisElementCompatibilityException>(top == BasisTopology::Unknown, SVMP_HERE,
                                                     "LagrangeBasis: unsupported element type");
    return top;
}

// Normalize named higher-order element requests to base Lagrange topologies.
//
// This function only adds the LagrangeBasis routing policy: serendipity
// layouts and pyramids are rejected, and a named quadratic alias
// (Line3, Triangle6, Quad9, Tetra10, Hex27, Wedge18) is floored
// to at least quadratic order. The floor only raises the
// order: a higher requested order on an alias is honored, so
// LagrangeBasis(Hex27, 5) yields an order-5 basis on the Hex8 topology rather
// than rejecting the over-specified order.
NormalizedLagrangeRequest normalize_lagrange_request(ElementType element_type, int order) {
    switch (element_type) {
        case ElementType::Quad8:
            svmp::raise<BasisElementCompatibilityException>(SVMP_HERE,
                "LagrangeBasis: Quad8 is serendipity; use SerendipityBasis");
        case ElementType::Hex20:
            svmp::raise<BasisElementCompatibilityException>(SVMP_HERE,
                "LagrangeBasis: Hex20 is serendipity; use SerendipityBasis");
        case ElementType::Wedge15:
            svmp::raise<BasisElementCompatibilityException>(SVMP_HERE,
                "LagrangeBasis: Wedge15 is serendipity; use SerendipityBasis");
        case ElementType::Pyramid5:
        case ElementType::Pyramid13:
        case ElementType::Pyramid14:
            svmp::raise<BasisElementCompatibilityException>(SVMP_HERE,
                "LagrangeBasis: pyramid support is not within the current solver basis scope");
        default:
            break;
    }

    const ElementType canonical = canonical_lagrange_type(element_type);
    const bool is_quadratic_alias = canonical != element_type;
    const int normalized_order = is_quadratic_alias ? std::max(order, 2) : order;
    return {canonical, normalized_order};
}

// Convert a coordinate on [-1, 1] to an equispaced axis node index.
std::size_t axis_index_pm_one(double coord, int order) {
    if (order <= 0) {
        return 0u;
    }
    const double scaled = (coord + double(1)) * double(order) / double(2);
    const long long rounded = std::llround(scaled);
    svmp::throw_if<BasisConstructionException>(
        rounded < 0 || rounded > static_cast<long long>(order) ||
            !detail::basis_nearly_equal(scaled, static_cast<double>(rounded)),
        SVMP_HERE,
        "LagrangeBasis: tensor-product node coordinate is off the equispaced lattice");
    return static_cast<std::size_t>(rounded);
}

// Convert a simplex barycentric coordinate to a lattice index.
int simplex_lattice_index(double value, int order) {
    if (order <= 0) {
        return 0;
    }
    const double scaled = value * double(order);
    const long long rounded = std::llround(scaled);
    svmp::throw_if<BasisConstructionException>(
        rounded < 0 || rounded > static_cast<long long>(order) ||
            !detail::basis_nearly_equal(scaled, static_cast<double>(rounded)),
        SVMP_HERE,
        "LagrangeBasis: simplex node coordinate is off the lattice");
    return static_cast<int>(rounded);
}

// Compute simplex interpolation exponents from a reference node.
LagrangeBasis::SimplexExponent simplex_exponent_from_point(const Vec3& p,
                                                           BasisTopology top,
                                                           int order) {
    LagrangeBasis::SimplexExponent e{0, 0, 0, 0};
    if (order <= 0) {
        return e;
    }
    if (top == BasisTopology::Triangle) {
        e[1] = simplex_lattice_index(p[0], order);
        e[2] = simplex_lattice_index(p[1], order);
        e[0] = order - e[1] - e[2];
    } else {
        e[1] = simplex_lattice_index(p[0], order);
        e[2] = simplex_lattice_index(p[1], order);
        e[3] = simplex_lattice_index(p[2], order);
        e[0] = order - e[1] - e[2] - e[3];
    }
    // e[0] is order minus the other exponents, so the exponents sum to order by
    // construction; a negative e[0] means the node coordinates are off-lattice.
    svmp::throw_if<BasisConstructionException>(
        e[0] < 0, SVMP_HERE,
        "LagrangeBasis: simplex node coordinate yields a negative implied exponent");
    return e;
}

// Sentinel node index meaning "skip nothing" in product_excluding below.
constexpr std::size_t kNoSkip = std::numeric_limits<std::size_t>::max();

// Evaluate 1D Lagrange polynomials and derivatives at a point. `level` selects
// how many derivative orders to compute: 0 for values only, 1 to also fill the
// first derivative, and 2 to also fill the second. The output arrays stay sized
// at n regardless of level so the tensor-product assembly can index them
// unconditionally; only the higher-order computation loops are skipped.
void evaluate_1d_lagrange(double x, const std::vector<double>& nodes, AxisEval& out,
                          int level) {
    const std::size_t n = nodes.size();
    out.value.assign(n, double(0));
    out.first.assign(n, double(0));
    out.second.assign(n, double(0));

    if (n == 1u) {
        out.value[0] = double(1);
        return;
    }

    for (std::size_t i = 0; i < n; ++i) {
        // Product of (x - nodes[j]) over all j except i and the listed skips.
        // Each derivative order drops one additional factor from the product.
        const auto product_excluding = [&](std::size_t skip1 = kNoSkip,
                                           std::size_t skip2 = kNoSkip) {
            double product = double(1);
            for (std::size_t j = 0; j < n; ++j) {
                if (j != i && j != skip1 && j != skip2) {
                    product *= x - nodes[j];
                }
            }
            return product;
        };

        double denom = double(1);
        for (std::size_t j = 0; j < n; ++j) {
            if (j != i) {
                denom *= nodes[i] - nodes[j];
            }
        }

        out.value[i] = product_excluding() / denom;

        if (level >= 1) {
            double first = double(0);
            for (std::size_t m = 0; m < n; ++m) {
                if (m != i) {
                    first += product_excluding(m);
                }
            }
            out.first[i] = first / denom;
        }

        if (level >= 2) {
            double second = double(0);
            for (std::size_t m = 0; m < n; ++m) {
                if (m == i) {
                    continue;
                }
                for (std::size_t l = 0; l < n; ++l) {
                    if (l != i && l != m) {
                        second += product_excluding(m, l);
                    }
                }
            }
            out.second[i] = second / denom;
        }
    }
}

// Evaluate one barycentric polynomial factor and its derivatives. `level`
// selects how far the recurrence runs: 0 for the value only, 1 to also produce
// the first derivative, and 2 to also produce the second.
std::array<double, 3> simplex_factor(int alpha, double lambda, int order, int level) {
    double value = double(1);
    double first = double(0);
    double second = double(0);

    for (int m = 0; m < alpha; ++m) {
        const double factor = double(order) * lambda - double(m);
        const double inv = double(1) / double(m + 1);
        const double old_value = value;
        const double old_first = first;
        const double old_second = second;
        value = old_value * factor * inv;
        if (level >= 1) {
            first = (old_first * factor + old_value * double(order)) * inv;
        }
        if (level >= 2) {
            second = (old_second * factor + double(2) * old_first * double(order)) * inv;
        }
    }

    return {value, first, second};
}

// Evaluate simplex Lagrange basis functions and the requested derivatives.
// Gradients and Hessians are assembled only when asked for; `out.gradient` and
// `out.hessian` are left empty otherwise so a values-only request neither
// allocates those buffers nor runs the derivative loops.
void evaluate_simplex(const Vec3& xi,
                      BasisTopology top,
                      int order,
                      const std::vector<LagrangeBasis::SimplexExponent>& exponents,
                      SimplexEval& out,
                      bool want_gradient,
                      bool want_hessian) {
    const std::size_t n = exponents.size();
    out.value.assign(n, double(0));
    out.gradient.assign(want_gradient ? n : std::size_t{0}, Gradient::Zero());
    out.hessian.assign(want_hessian ? n : std::size_t{0}, Hessian::Zero());

    if (n == 1u && order == 0) {
        out.value[0] = double(1);
        return;
    }

    // A Hessian factor also needs the first-derivative recurrence, so the
    // per-factor work runs to the highest requested order.
    const int factor_level = want_hessian ? 2 : (want_gradient ? 1 : 0);

    const std::size_t bary_count = top == BasisTopology::Triangle ? 3u : 4u;
    std::array<double, 4> lambda{double(0), double(0), double(0), double(0)};
    std::array<Gradient, 4> lambda_grad;
    lambda_grad.fill(Gradient::Zero());

    lambda[1] = xi[0];
    lambda[2] = xi[1];
    lambda_grad[1][0] = double(1);
    lambda_grad[2][1] = double(1);
    if (top == BasisTopology::Triangle) {
        lambda[0] = double(1) - xi[0] - xi[1];
        lambda_grad[0][0] = double(-1);
        lambda_grad[0][1] = double(-1);
    } else {
        lambda[3] = xi[2];
        lambda[0] = double(1) - xi[0] - xi[1] - xi[2];
        lambda_grad[0][0] = double(-1);
        lambda_grad[0][1] = double(-1);
        lambda_grad[0][2] = double(-1);
        lambda_grad[3][2] = double(1);
    }

    for (std::size_t i = 0; i < n; ++i) {
        std::array<std::array<double, 3>, 4> f{};
        for (std::size_t a = 0; a < bary_count; ++a) {
            f[a] = simplex_factor(exponents[i][a], lambda[a], order, factor_level);
        }

        double value = double(1);
        for (std::size_t a = 0; a < bary_count; ++a) {
            value *= f[a][0];
        }
        out.value[i] = value;

        if (want_gradient) {
            for (std::size_t a = 0; a < bary_count; ++a) {
                double product = f[a][1];
                for (std::size_t b = 0; b < bary_count; ++b) {
                    if (b != a) {
                        product *= f[b][0];
                    }
                }
                for (std::size_t c = 0; c < 3u; ++c) {
                    out.gradient[i][c] += product * lambda_grad[a][c];
                }
            }
        }

        if (want_hessian) {
            for (std::size_t a = 0; a < bary_count; ++a) {
                for (std::size_t b = 0; b < bary_count; ++b) {
                    double product = (a == b) ? f[a][2] : f[a][1] * f[b][1];
                    for (std::size_t k = 0; k < bary_count; ++k) {
                        if (k != a && k != b) {
                            product *= f[k][0];
                        }
                    }
                    for (std::size_t r = 0; r < 3u; ++r) {
                        for (std::size_t c = 0; c < 3u; ++c) {
                            out.hessian[i](r, c) +=
                                product * lambda_grad[a][r] * lambda_grad[b][c];
                        }
                    }
                }
            }
        }
    }
}

} // namespace

LagrangeBasis::LagrangeBasis(ElementType type, int order)
    : element_type_(type), order_(order) {
    const auto normalized = normalize_lagrange_request(element_type_, order_);
    element_type_ = normalized.element_type;
    order_ = normalized.order;
    svmp::throw_if<BasisConfigurationException>(order_ < 0, SVMP_HERE,
                                              "LagrangeBasis requires non-negative polynomial order");

    topology_ = supported_lagrange_topology(element_type_);
    dimension_ = reference_dimension(element_type_);
    init_nodes();
}

// Initialize equispaced 1D interpolation nodes for tensor-product axes.
void LagrangeBasis::init_equispaced_1d_nodes() {
    nodes_1d_.resize(static_cast<std::size_t>(order_ + 1));
    for (int i = 0; i <= order_; ++i) {
        nodes_1d_[static_cast<std::size_t>(i)] =
            line_coord_pm_one(i, order_);
    }
}

// Initialize reference nodes and topology-specific lookup data.
void LagrangeBasis::init_nodes() {
    nodes_.clear();
    nodes_1d_.clear();
    tensor_indices_.clear();
    simplex_exponents_.clear();
    wedge_indices_.clear();

    switch (topology_) {
        case BasisTopology::Point:
            build_point_nodes();
            return;
        case BasisTopology::Line:
        case BasisTopology::Quadrilateral:
        case BasisTopology::Hexahedron:
            build_tensor_product_nodes();
            return;
        case BasisTopology::Triangle:
        case BasisTopology::Tetrahedron:
            build_simplex_nodes();
            return;
        case BasisTopology::Wedge:
            build_wedge_nodes();
            return;
        default:
            break;
    }

    svmp::raise<BasisElementCompatibilityException>(SVMP_HERE,
        "Unsupported element type in LagrangeBasis::init_nodes");
}

// Build the single reference node for a point basis.
void LagrangeBasis::build_point_nodes() {
    nodes_.push_back(Vec3{double(0), double(0), double(0)});
}

// Build nodes and axis indices for tensor-product elements.
void LagrangeBasis::build_tensor_product_nodes() {
    init_equispaced_1d_nodes();
    nodes_ = ReferenceNodeLayout::get_lagrange_node_coords(element_type_, order_);
    tensor_indices_.reserve(nodes_.size());
    for (const auto& node : nodes_) {
        TensorNodeIndex idx{0u, 0u, 0u};
        idx[0] = axis_index_pm_one(node[0], order_);
        if (dimension_ >= 2) {
            idx[1] = axis_index_pm_one(node[1], order_);
        }
        if (dimension_ >= 3) {
            idx[2] = axis_index_pm_one(node[2], order_);
        }
        tensor_indices_.push_back(idx);
    }
}

// Build nodes and barycentric exponents for simplex elements.
void LagrangeBasis::build_simplex_nodes() {
    nodes_ = ReferenceNodeLayout::get_lagrange_node_coords(element_type_, order_);
    simplex_exponents_.reserve(nodes_.size());
    for (const auto& node : nodes_) {
        simplex_exponents_.push_back(simplex_exponent_from_point(node, topology_, order_));
    }
}

// Build nodes and mixed triangle-axis lookup data for wedge elements.
void LagrangeBasis::build_wedge_nodes() {
    init_equispaced_1d_nodes();
    nodes_ = ReferenceNodeLayout::get_lagrange_node_coords(element_type_, order_);
    const auto tri_nodes =
        ReferenceNodeLayout::get_lagrange_node_coords(ElementType::Triangle3, order_);
    simplex_exponents_.reserve(tri_nodes.size());
    for (const auto& tri_node : tri_nodes) {
        simplex_exponents_.push_back(
            simplex_exponent_from_point(tri_node, BasisTopology::Triangle, order_));
    }

    wedge_indices_.reserve(nodes_.size());
    for (const auto& node : nodes_) {
        const auto tri_exp =
            simplex_exponent_from_point(node, BasisTopology::Triangle, order_);
        auto it = std::find(simplex_exponents_.begin(), simplex_exponents_.end(), tri_exp);
        svmp::throw_if<BasisConstructionException>(it == simplex_exponents_.end(), SVMP_HERE,
                                                 "LagrangeBasis: wedge node triangle index lookup failed");
        const std::size_t tri_index =
            static_cast<std::size_t>(std::distance(simplex_exponents_.begin(), it));
        wedge_indices_.push_back({tri_index, axis_index_pm_one(node[2], order_)});
    }
}

// Evaluate the constant point basis.
void LagrangeBasis::evaluate_point_to(std::span<double> values_out,
                                      std::span<Gradient> gradients_out,
                                      std::span<Hessian> hessians_out) const {
    if (!values_out.empty()) {
        values_out[0] = double(1);
    }
    if (!gradients_out.empty()) {
        gradients_out[0] = Gradient::Zero();
    }
    if (!hessians_out.empty()) {
        hessians_out[0] = Hessian::Zero();
    }
}

// Evaluate line, quadrilateral, and hexahedron bases as axis-polynomial products.
void LagrangeBasis::evaluate_tensor_product_to(const Vec3& xi,
                                               std::span<double> values_out,
                                               std::span<Gradient> gradients_out,
                                               std::span<Hessian> hessians_out) const {
    const int level = !hessians_out.empty() ? 2 : (!gradients_out.empty() ? 1 : 0);

    AxisEval ax;
    AxisEval ay;
    AxisEval az;
    evaluate_1d_lagrange(xi[0], nodes_1d_, ax, level);
    if (dimension_ >= 2) {
        evaluate_1d_lagrange(xi[1], nodes_1d_, ay, level);
    }
    if (dimension_ >= 3) {
        evaluate_1d_lagrange(xi[2], nodes_1d_, az, level);
    }

    for (std::size_t node = 0; node < tensor_indices_.size(); ++node) {
        const auto& idx = tensor_indices_[node];
        const double vx = ax.value[idx[0]];
        const double dx = ax.first[idx[0]];
        const double d2x = ax.second[idx[0]];
        const double vy = dimension_ >= 2 ? ay.value[idx[1]] : double(1);
        const double dy = dimension_ >= 2 ? ay.first[idx[1]] : double(0);
        const double d2y = dimension_ >= 2 ? ay.second[idx[1]] : double(0);
        const double vz = dimension_ >= 3 ? az.value[idx[2]] : double(1);
        const double dz = dimension_ >= 3 ? az.first[idx[2]] : double(0);
        const double d2z = dimension_ >= 3 ? az.second[idx[2]] : double(0);

        if (!values_out.empty()) {
            values_out[node] = vx * vy * vz;
        }
        if (!gradients_out.empty()) {
            Gradient& g = gradients_out[node];
            g[0] = dx * vy * vz;
            g[1] = vx * dy * vz;
            g[2] = vx * vy * dz;
        }
        if (!hessians_out.empty()) {
            Hessian& h = hessians_out[node];
            h(0, 0) = d2x * vy * vz;
            h(0, 1) = dx * dy * vz;
            h(0, 2) = dx * vy * dz;
            h(1, 0) = h(0, 1);
            h(1, 1) = vx * d2y * vz;
            h(1, 2) = vx * dy * dz;
            h(2, 0) = h(0, 2);
            h(2, 1) = h(1, 2);
            h(2, 2) = vx * vy * d2z;
        }
    }
}

// Evaluate triangle and tetrahedron bases from barycentric factors.
void LagrangeBasis::evaluate_simplex_to(const Vec3& xi,
                                        std::span<double> values_out,
                                        std::span<Gradient> gradients_out,
                                        std::span<Hessian> hessians_out) const {
    const bool want_values = !values_out.empty();
    const bool want_gradients = !gradients_out.empty();
    const bool want_hessians = !hessians_out.empty();

    SimplexEval simplex;
    evaluate_simplex(xi, topology_, order_, simplex_exponents_, simplex,
                     want_gradients, want_hessians);
    for (std::size_t i = 0; i < simplex.value.size(); ++i) {
        if (want_values) {
            values_out[i] = simplex.value[i];
        }
        if (want_gradients) {
            gradients_out[i] = simplex.gradient[i];
        }
        if (want_hessians) {
            hessians_out[i] = simplex.hessian[i];
        }
    }
}

// Evaluate wedge bases as triangle/through-axis products.
void LagrangeBasis::evaluate_wedge_to(const Vec3& xi,
                                      std::span<double> values_out,
                                      std::span<Gradient> gradients_out,
                                      std::span<Hessian> hessians_out) const {
    const bool want_values = !values_out.empty();
    const bool want_gradients = !gradients_out.empty();
    const bool want_hessians = !hessians_out.empty();

    // The wedge gradient pairs the triangle gradient with the through-axis value,
    // and the wedge Hessian reuses the triangle gradient for its mixed terms, so
    // the triangle factor must supply gradients whenever the wedge needs either
    // gradients or Hessians.
    const bool want_tri_gradient = want_gradients || want_hessians;
    const int z_level = want_hessians ? 2 : (want_gradients ? 1 : 0);

    SimplexEval tri;
    AxisEval z_axis;
    evaluate_simplex(xi, BasisTopology::Triangle, order_, simplex_exponents_, tri,
                     want_tri_gradient, want_hessians);
    evaluate_1d_lagrange(xi[2], nodes_1d_, z_axis, z_level);

    for (std::size_t node = 0; node < wedge_indices_.size(); ++node) {
        const auto [tri_idx, z_idx] = wedge_indices_[node];
        const double tv = tri.value[tri_idx];
        const double zv = z_axis.value[z_idx];

        if (want_values) {
            values_out[node] = tv * zv;
        }
        if (want_gradients) {
            const double dz = z_axis.first[z_idx];
            Gradient& g = gradients_out[node];
            g[0] = tri.gradient[tri_idx][0] * zv;
            g[1] = tri.gradient[tri_idx][1] * zv;
            g[2] = tv * dz;
        }
        if (want_hessians) {
            const double dz = z_axis.first[z_idx];
            const double d2z = z_axis.second[z_idx];
            Hessian& h = hessians_out[node];
            const Hessian& th = tri.hessian[tri_idx];
            const Gradient& tg = tri.gradient[tri_idx];
            h(0, 0) = th(0, 0) * zv;
            h(0, 1) = th(0, 1) * zv;
            h(0, 2) = tg[0] * dz;
            h(1, 0) = h(0, 1);
            h(1, 1) = th(1, 1) * zv;
            h(1, 2) = tg[1] * dz;
            h(2, 0) = h(0, 2);
            h(2, 1) = h(1, 2);
            h(2, 2) = tv * d2z;
        }
    }
}

// Evaluate requested basis quantities into caller-provided spans.
void LagrangeBasis::evaluate_all_to(const Vec3& xi,
                                    std::span<double> values_out,
                                    std::span<Gradient> gradients_out,
                                    std::span<Hessian> hessians_out) const {
    require_requested_span_size(values_out, size(), "LagrangeBasis::evaluate_all_to values");
    require_requested_span_size(gradients_out, size(), "LagrangeBasis::evaluate_all_to gradients");
    require_requested_span_size(hessians_out, size(), "LagrangeBasis::evaluate_all_to hessians");

    if (values_out.empty() && gradients_out.empty() && hessians_out.empty()) {
        return;
    }

    switch (topology_) {
        case BasisTopology::Point:
            evaluate_point_to(values_out, gradients_out, hessians_out);
            return;
        case BasisTopology::Line:
        case BasisTopology::Quadrilateral:
        case BasisTopology::Hexahedron:
            evaluate_tensor_product_to(xi, values_out, gradients_out, hessians_out);
            return;
        case BasisTopology::Triangle:
        case BasisTopology::Tetrahedron:
            evaluate_simplex_to(xi, values_out, gradients_out, hessians_out);
            return;
        case BasisTopology::Wedge:
            evaluate_wedge_to(xi, values_out, gradients_out, hessians_out);
            return;
        default:
            break;
    }

    svmp::raise<BasisEvaluationException>(SVMP_HERE,
        "Unsupported element in LagrangeBasis evaluation");
}

void LagrangeBasis::evaluate_values(const Vec3& xi,
                                    std::vector<double>& values) const {
    values.resize(size());
    evaluate_values_to(xi, std::span<double>(values.data(), values.size()));
}

void LagrangeBasis::evaluate_gradients(const Vec3& xi,
                                       std::vector<Gradient>& gradients) const {
    gradients.resize(size());
    evaluate_gradients_to(xi, std::span<Gradient>(gradients.data(), gradients.size()));
}

void LagrangeBasis::evaluate_hessians(const Vec3& xi,
                                      std::vector<Hessian>& hessians) const {
    hessians.resize(size());
    evaluate_hessians_to(xi, std::span<Hessian>(hessians.data(), hessians.size()));
}

void LagrangeBasis::evaluate_all(const Vec3& xi,
                                 std::vector<double>& values,
                                 std::vector<Gradient>& gradients,
                                 std::vector<Hessian>& hessians) const {
    values.resize(size());
    gradients.resize(size());
    hessians.resize(size());
    evaluate_all_to(xi,
                    std::span<double>(values.data(), values.size()),
                    std::span<Gradient>(gradients.data(), gradients.size()),
                    std::span<Hessian>(hessians.data(), hessians.size()));
}

void LagrangeBasis::evaluate_values_to(const Vec3& xi,
                                       std::span<double> values_out) const {
    require_span_size(values_out.size(), size(), "LagrangeBasis::evaluate_values_to");
    evaluate_all_to(xi, values_out, std::span<Gradient>{}, std::span<Hessian>{});
}

void LagrangeBasis::evaluate_gradients_to(const Vec3& xi,
                                          std::span<Gradient> gradients_out) const {
    require_span_size(gradients_out.size(), size(), "LagrangeBasis::evaluate_gradients_to");
    evaluate_all_to(xi, std::span<double>{}, gradients_out, std::span<Hessian>{});
}

void LagrangeBasis::evaluate_hessians_to(const Vec3& xi,
                                         std::span<Hessian> hessians_out) const {
    require_span_size(hessians_out.size(), size(), "LagrangeBasis::evaluate_hessians_to");
    evaluate_all_to(xi, std::span<double>{}, std::span<Gradient>{}, hessians_out);
}

} // namespace basis
} // namespace FE
} // namespace svmp
