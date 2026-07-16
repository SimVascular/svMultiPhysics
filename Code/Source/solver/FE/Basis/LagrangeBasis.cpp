// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "LagrangeBasis.h"
#include "NodeOrderingConventions.h"

#include <algorithm>
#include <array>
#include <span>
#include <string>

namespace svmp::FE::basis {

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

// Validate a named-element Lagrange request and return its reference topology.
//
// Serendipity layouts and pyramids are rejected. The requested order must equal
// the order baked into the named layout (0 for Point1, 1 for the linear
// elements, 2 for the complete-quadratic aliases Line3/Triangle6/Quad9/Tetra10/
// Hex27/Wedge18). A named element therefore cannot carry a conflicting order;
// arbitrary orders are requested through the BasisTopology constructor, so a
// reader never has to read a polynomial order out of a node-count name such as
// Hex8 or Tetra10.
BasisTopology validated_lagrange_topology(ElementType element_type, int order) {
    switch (element_type) {
        case ElementType::Quad8:
            svmp::raise<BasisElementCompatibilityException>("LagrangeBasis: Quad8 is serendipity; use SerendipityBasis");
        case ElementType::Hex20:
            svmp::raise<BasisElementCompatibilityException>("LagrangeBasis: Hex20 is serendipity; use SerendipityBasis");
        case ElementType::Wedge15:
            svmp::raise<BasisElementCompatibilityException>("LagrangeBasis: Wedge15 is serendipity; use SerendipityBasis");
        case ElementType::Pyramid5:
        case ElementType::Pyramid13:
        case ElementType::Pyramid14:
            svmp::raise<BasisElementCompatibilityException>("LagrangeBasis: pyramid support is not within the current solver basis scope");
        default:
            break;
    }

    const BasisTopology top = topology(element_type);
    svmp::throw_if<BasisElementCompatibilityException>(top == BasisTopology::Unknown, "LagrangeBasis: unsupported element type");

    const int baked_order = named_lagrange_order(element_type);
    svmp::throw_if<BasisConfigurationException>(order != baked_order, "LagrangeBasis: a named element layout has a fixed polynomial order; request the matching "
        "BasisTopology with an explicit order to choose a different order");
    return top;
}

// Convert an integer lattice index (i, j[, k]) into the barycentric exponent
// tuple (order - i - j - k, i, j, k). The lattice already carries the exact
// coordinate indices, the accessor's structural invariants guarantee
//  i + j + k <= order, hence e[0] >= 0.
LagrangeBasis::SimplexExponent simplex_exponent_from_lattice(const std::array<int, 3>& idx,
                                                            BasisTopology top,
                                                            int order) {
    LagrangeBasis::SimplexExponent e{0, 0, 0, 0};
    e[1] = idx[0];
    e[2] = idx[1];
    if (top == BasisTopology::Tetrahedron) {
        e[3] = idx[2];
        e[0] = order - idx[0] - idx[1] - idx[2];
    } else {
        e[0] = order - idx[0] - idx[1];
    }
    return e;
}

// Evaluate 1D Lagrange polynomials and their derivatives at a point in the
// barycentric form l_i(x) = w_i * prod_{j!=i}(x - x_j), where the weights
// w_i = 1 / prod_{j!=i}(x_i - x_j) depend only on the fixed node set and are
// precomputed once by the caller. For each i the numerator and its first and
// second derivatives are built by a single product-rule accumulation over the
// remaining nodes.
void evaluate_1d_lagrange(double x,
                          const std::vector<double>& nodes,
                          const std::vector<double>& weights,
                          AxisEval& out,
                          int level) {
    const std::size_t n = nodes.size();
    out.value.assign(n, double(0));
    out.first.assign(n, double(0));
    out.second.assign(n, double(0));

    for (std::size_t i = 0; i < n; ++i) {
        double value = double(1);
        double first = double(0);
        double second = double(0);
        for (std::size_t j = 0; j < n; ++j) {
            if (j == i) {
                continue;
            }
            const double f = x - nodes[j];
            if (level >= 2) {
                second = second * f + double(2) * first;
            }
            if (level >= 1) {
                first = first * f + value;
            }
            value = value * f;
        }

        const double w = weights[i];
        out.value[i] = value * w;
        if (level >= 1) {
            out.first[i] = first * w;
        }
        if (level >= 2) {
            out.second[i] = second * w;
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

LagrangeBasis::LagrangeBasis(BasisTopology topology, int order)
    : topology_(topology), order_(order) {
    svmp::throw_if<BasisElementCompatibilityException>(topology_ == BasisTopology::Unknown, "LagrangeBasis: unknown reference topology");
    svmp::throw_if<BasisConfigurationException>(order_ < 0, "LagrangeBasis requires non-negative polynomial order");
    svmp::throw_if<BasisConfigurationException>(
        topology_ == BasisTopology::Point && order_ != 0, "LagrangeBasis: Point topology supports order 0 only");
    dimension_ = topology_dimension(topology_);
    init_nodes();
}

LagrangeBasis::LagrangeBasis(ElementType type, int order)
    : LagrangeBasis(validated_lagrange_topology(type, order), order) {}

LagrangeBasis::LagrangeBasis(ElementType type)
    : LagrangeBasis(type, named_lagrange_order(type)) {}

// Initialize the 1D tensor-axis interpolation nodes (Gauss-Lobatto-Legendre, via
// line_coord_pm_one) and their barycentric weights for tensor-product axes.
void LagrangeBasis::init_tensor_axis_nodes() {
    const std::size_t n = static_cast<std::size_t>(order_ + 1);
    nodes_1d_.resize(n);
    for (int i = 0; i <= order_; ++i) {
        nodes_1d_[static_cast<std::size_t>(i)] =
            line_coord_pm_one(i, order_);
    }

    // Barycentric weights w_i = 1 / prod_{j!=i}(x_i - x_j); the nodes are
    // distinct so every denominator is nonzero.
    nodes_1d_weights_.assign(n, double(1));
    for (std::size_t i = 0; i < n; ++i) {
        double denom = double(1);
        for (std::size_t j = 0; j < n; ++j) {
            if (j != i) {
                denom *= nodes_1d_[i] - nodes_1d_[j];
            }
        }
        nodes_1d_weights_[i] = double(1) / denom;
    }
}

// Initialize reference nodes and topology-specific lookup data.
void LagrangeBasis::init_nodes() {
    nodes_.clear();
    nodes_1d_.clear();
    nodes_1d_weights_.clear();
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

    svmp::raise<BasisElementCompatibilityException>("Unsupported element type in LagrangeBasis::init_nodes");
}

// Build the single reference node for a point basis.
void LagrangeBasis::build_point_nodes() {
    nodes_.push_back(Vec3{double(0), double(0), double(0)});
}

// Build nodes and axis indices for tensor-product elements.
void LagrangeBasis::build_tensor_product_nodes() {
    init_tensor_axis_nodes();
    const auto layout =
        ReferenceNodeLayout::get_lagrange_lattice(lagrange_topology_representative(topology_), order_);
    nodes_ = layout.coords;
    tensor_indices_.reserve(layout.lattice.size());
    for (const auto& idx : layout.lattice) {
        // The lattice already holds the per-axis node index 0..order along each
        // axis (unused axes are zero; the coordinate itself is the GLL node for
        // that index), so no coordinate-to-index inversion is needed.
        tensor_indices_.push_back(TensorNodeIndex{
            static_cast<std::size_t>(idx[0]),
            static_cast<std::size_t>(idx[1]),
            static_cast<std::size_t>(idx[2])});
    }
}

// Build nodes and barycentric exponents for simplex elements.
void LagrangeBasis::build_simplex_nodes() {
    const auto layout =
        ReferenceNodeLayout::get_lagrange_lattice(lagrange_topology_representative(topology_), order_);
    nodes_ = layout.coords;
    simplex_exponents_.reserve(layout.lattice.size());
    for (const auto& idx : layout.lattice) {
        simplex_exponents_.push_back(simplex_exponent_from_lattice(idx, topology_, order_));
    }
}

// Build nodes and mixed triangle-axis lookup data for wedge elements.
void LagrangeBasis::build_wedge_nodes() {
    init_tensor_axis_nodes();
    const auto layout =
        ReferenceNodeLayout::get_lagrange_lattice(lagrange_topology_representative(topology_), order_);
    nodes_ = layout.coords;

    const auto tri_layout =
        ReferenceNodeLayout::get_lagrange_lattice(ElementType::Triangle3, order_);
    simplex_exponents_.reserve(tri_layout.lattice.size());
    for (const auto& idx : tri_layout.lattice) {
        simplex_exponents_.push_back(
            simplex_exponent_from_lattice(idx, BasisTopology::Triangle, order_));
    }

    // Map a triangle cross-section lattice (i, j) to its triangle-node ordinal
    // through the flat key i * (order + 1) + j, so each wedge node's triangle
    // index is an exact integer lookup.
    const int stride = order_ + 1;
    std::vector<int> tri_ordinal_for_key(static_cast<std::size_t>(stride * stride), -1);
    for (std::size_t t = 0; t < tri_layout.lattice.size(); ++t) {
        const auto& idx = tri_layout.lattice[t];
        tri_ordinal_for_key[static_cast<std::size_t>(idx[0] * stride + idx[1])] =
            static_cast<int>(t);
    }

    wedge_indices_.reserve(layout.lattice.size());
    for (const auto& idx : layout.lattice) {
        const int tri_ordinal =
            tri_ordinal_for_key[static_cast<std::size_t>(idx[0] * stride + idx[1])];
        svmp::throw_if<BasisConstructionException>(tri_ordinal < 0, "LagrangeBasis: wedge node triangle index lookup failed");
        wedge_indices_.push_back({static_cast<std::size_t>(tri_ordinal),
                                  static_cast<std::size_t>(idx[2])});
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
    evaluate_1d_lagrange(xi[0], nodes_1d_, nodes_1d_weights_, ax, level);
    if (dimension_ >= 2) {
        evaluate_1d_lagrange(xi[1], nodes_1d_, nodes_1d_weights_, ay, level);
    }
    if (dimension_ >= 3) {
        evaluate_1d_lagrange(xi[2], nodes_1d_, nodes_1d_weights_, az, level);
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
    evaluate_1d_lagrange(xi[2], nodes_1d_, nodes_1d_weights_, z_axis, z_level);

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
    // Private sink: callers guarantee valid output spans -- the public *_to methods
    // validate their one output with require_span_size, and the vector evaluators
    // resize to size(). An empty span here means "skip that quantity".

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

    svmp::raise<BasisEvaluationException>("Unsupported element in LagrangeBasis evaluation");
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

} // namespace svmp::FE::basis
