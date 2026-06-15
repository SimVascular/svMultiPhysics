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

using Vec3 = math::Vector<Real, 3>;

// Return the equispaced 1D reference coordinate in [-1, 1].
inline constexpr Real equispaced_pm_one_coord(int i, int order) {
    if (order <= 0) {
        return Real(0);
    }
    return Real(-1) + Real(2) * static_cast<Real>(i) / static_cast<Real>(order);
}

struct AxisEval {
    std::vector<Real> value;
    std::vector<Real> first;
    std::vector<Real> second;
};

struct SimplexEval {
    std::vector<Real> value;
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
    FE::throw_if<BasisElementCompatibilityException>(top == BasisTopology::Unknown, SVMP_HERE,
                                                     "LagrangeBasis: unsupported element type");
    return top;
}

// Normalize named higher-order element requests to base Lagrange topologies.
NormalizedLagrangeRequest normalize_lagrange_request(ElementType element_type, int order) {
    switch (element_type) {
        case ElementType::Line3:
            return {ElementType::Line2, std::max(order, 2)};
        case ElementType::Triangle6:
            return {ElementType::Triangle3, std::max(order, 2)};
        case ElementType::Quad9:
            return {ElementType::Quad4, std::max(order, 2)};
        case ElementType::Tetra10:
            return {ElementType::Tetra4, std::max(order, 2)};
        case ElementType::Hex27:
            return {ElementType::Hex8, std::max(order, 2)};
        case ElementType::Wedge18:
            return {ElementType::Wedge6, std::max(order, 2)};
        case ElementType::Quad8:
            FE::raise<BasisElementCompatibilityException>(SVMP_HERE,
                "LagrangeBasis: Quad8 is serendipity; use SerendipityBasis");
        case ElementType::Hex20:
            FE::raise<BasisElementCompatibilityException>(SVMP_HERE,
                "LagrangeBasis: Hex20 is serendipity; use SerendipityBasis");
        case ElementType::Wedge15:
            FE::raise<BasisElementCompatibilityException>(SVMP_HERE,
                "LagrangeBasis: Wedge15 is serendipity; use SerendipityBasis");
        case ElementType::Pyramid5:
        case ElementType::Pyramid13:
        case ElementType::Pyramid14:
            FE::raise<BasisElementCompatibilityException>(SVMP_HERE,
                "LagrangeBasis: pyramid support is not within the current solver basis scope");
        default:
            return {element_type, order};
    }
}

// Convert a coordinate on [-1, 1] to an equispaced axis node index.
std::size_t axis_index_pm_one(Real coord, int order) {
    if (order <= 0) {
        return 0u;
    }
    const Real scaled = (coord + Real(1)) * Real(order) / Real(2);
    return static_cast<std::size_t>(std::llround(scaled));
}

// Convert a simplex barycentric coordinate to a lattice index.
int simplex_lattice_index(Real value, int order) {
    if (order <= 0) {
        return 0;
    }
    return static_cast<int>(std::llround(value * Real(order)));
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
    return e;
}

// Sentinel node index meaning "skip nothing" in product_excluding below.
constexpr std::size_t kNoSkip = std::numeric_limits<std::size_t>::max();

// Evaluate 1D Lagrange polynomials and derivatives at a point.
void evaluate_1d_lagrange(Real x, const std::vector<Real>& nodes, AxisEval& out) {
    const std::size_t n = nodes.size();
    out.value.assign(n, Real(0));
    out.first.assign(n, Real(0));
    out.second.assign(n, Real(0));

    if (n == 1u) {
        out.value[0] = Real(1);
        return;
    }

    for (std::size_t i = 0; i < n; ++i) {
        // Product of (x - nodes[j]) over all j except i and the listed skips.
        // Each derivative order drops one additional factor from the product.
        const auto product_excluding = [&](std::size_t skip1 = kNoSkip,
                                           std::size_t skip2 = kNoSkip) {
            Real product = Real(1);
            for (std::size_t j = 0; j < n; ++j) {
                if (j != i && j != skip1 && j != skip2) {
                    product *= x - nodes[j];
                }
            }
            return product;
        };

        Real denom = Real(1);
        for (std::size_t j = 0; j < n; ++j) {
            if (j != i) {
                denom *= nodes[i] - nodes[j];
            }
        }

        out.value[i] = product_excluding() / denom;

        Real first = Real(0);
        for (std::size_t m = 0; m < n; ++m) {
            if (m != i) {
                first += product_excluding(m);
            }
        }
        out.first[i] = first / denom;

        Real second = Real(0);
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

// Evaluate one barycentric polynomial factor and derivatives.
std::array<Real, 3> simplex_factor(int alpha, Real lambda, int order) {
    Real value = Real(1);
    Real first = Real(0);
    Real second = Real(0);

    for (int m = 0; m < alpha; ++m) {
        const Real factor = Real(order) * lambda - Real(m);
        const Real inv = Real(1) / Real(m + 1);
        const Real old_value = value;
        const Real old_first = first;
        const Real old_second = second;
        value = old_value * factor * inv;
        first = (old_first * factor + old_value * Real(order)) * inv;
        second = (old_second * factor + Real(2) * old_first * Real(order)) * inv;
    }

    return {value, first, second};
}

// Evaluate simplex Lagrange basis functions and derivatives.
void evaluate_simplex(const Vec3& xi,
                      BasisTopology top,
                      int order,
                      const std::vector<LagrangeBasis::SimplexExponent>& exponents,
                      SimplexEval& out) {
    const std::size_t n = exponents.size();
    out.value.assign(n, Real(0));
    out.gradient.assign(n, Gradient::Zero());
    out.hessian.assign(n, Hessian::Zero());

    if (n == 1u && order == 0) {
        out.value[0] = Real(1);
        return;
    }

    const std::size_t bary_count = top == BasisTopology::Triangle ? 3u : 4u;
    std::array<Real, 4> lambda{Real(0), Real(0), Real(0), Real(0)};
    std::array<Gradient, 4> lambda_grad;
    lambda_grad.fill(Gradient::Zero());

    lambda[1] = xi[0];
    lambda[2] = xi[1];
    lambda_grad[1][0] = Real(1);
    lambda_grad[2][1] = Real(1);
    if (top == BasisTopology::Triangle) {
        lambda[0] = Real(1) - xi[0] - xi[1];
        lambda_grad[0][0] = Real(-1);
        lambda_grad[0][1] = Real(-1);
    } else {
        lambda[3] = xi[2];
        lambda[0] = Real(1) - xi[0] - xi[1] - xi[2];
        lambda_grad[0][0] = Real(-1);
        lambda_grad[0][1] = Real(-1);
        lambda_grad[0][2] = Real(-1);
        lambda_grad[3][2] = Real(1);
    }

    for (std::size_t i = 0; i < n; ++i) {
        std::array<std::array<Real, 3>, 4> f{};
        for (std::size_t a = 0; a < bary_count; ++a) {
            f[a] = simplex_factor(exponents[i][a], lambda[a], order);
        }

        Real value = Real(1);
        for (std::size_t a = 0; a < bary_count; ++a) {
            value *= f[a][0];
        }
        out.value[i] = value;

        for (std::size_t a = 0; a < bary_count; ++a) {
            Real product = f[a][1];
            for (std::size_t b = 0; b < bary_count; ++b) {
                if (b != a) {
                    product *= f[b][0];
                }
            }
            for (std::size_t c = 0; c < 3u; ++c) {
                out.gradient[i][c] += product * lambda_grad[a][c];
            }
        }

        for (std::size_t a = 0; a < bary_count; ++a) {
            for (std::size_t b = 0; b < bary_count; ++b) {
                Real product = (a == b) ? f[a][2] : f[a][1] * f[b][1];
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

void require_output_span_size(std::size_t actual,
                              std::size_t expected,
                              const char* label) {
    FE::throw_if<BasisEvaluationException>(actual < expected, SVMP_HERE,
        std::string(label) + ": output span is smaller than basis size");
}

template<typename T>
void require_requested_span_size(std::span<T> output,
                                 std::size_t expected,
                                 const char* label) {
    if (!output.empty()) {
        require_output_span_size(output.size(), expected, label);
    }
}

} // namespace

LagrangeBasis::LagrangeBasis(ElementType type, int order)
    : element_type_(type), order_(order) {
    const auto normalized = normalize_lagrange_request(element_type_, order_);
    element_type_ = normalized.element_type;
    order_ = normalized.order;
    FE::throw_if<BasisConfigurationException>(order_ < 0, SVMP_HERE,
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
            equispaced_pm_one_coord(i, order_);
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

    FE::raise<BasisElementCompatibilityException>(SVMP_HERE,
        "Unsupported element type in LagrangeBasis::init_nodes");
}

// Build the single reference node for a point basis.
void LagrangeBasis::build_point_nodes() {
    nodes_.push_back(Vec3{Real(0), Real(0), Real(0)});
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
        FE::throw_if<BasisConstructionException>(it == simplex_exponents_.end(), SVMP_HERE,
                                                 "LagrangeBasis: wedge node triangle index lookup failed");
        const std::size_t tri_index =
            static_cast<std::size_t>(std::distance(simplex_exponents_.begin(), it));
        wedge_indices_.push_back({tri_index, axis_index_pm_one(node[2], order_)});
    }
}

// Evaluate the constant point basis.
void LagrangeBasis::evaluate_point_to(std::span<Real> values_out,
                                      std::span<Gradient> gradients_out,
                                      std::span<Hessian> hessians_out) const {
    if (!values_out.empty()) {
        values_out[0] = Real(1);
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
                                               std::span<Real> values_out,
                                               std::span<Gradient> gradients_out,
                                               std::span<Hessian> hessians_out) const {
    AxisEval ax;
    AxisEval ay;
    AxisEval az;
    evaluate_1d_lagrange(xi[0], nodes_1d_, ax);
    if (dimension_ >= 2) {
        evaluate_1d_lagrange(xi[1], nodes_1d_, ay);
    }
    if (dimension_ >= 3) {
        evaluate_1d_lagrange(xi[2], nodes_1d_, az);
    }

    for (std::size_t node = 0; node < tensor_indices_.size(); ++node) {
        const auto& idx = tensor_indices_[node];
        const Real vx = ax.value[idx[0]];
        const Real dx = ax.first[idx[0]];
        const Real d2x = ax.second[idx[0]];
        const Real vy = dimension_ >= 2 ? ay.value[idx[1]] : Real(1);
        const Real dy = dimension_ >= 2 ? ay.first[idx[1]] : Real(0);
        const Real d2y = dimension_ >= 2 ? ay.second[idx[1]] : Real(0);
        const Real vz = dimension_ >= 3 ? az.value[idx[2]] : Real(1);
        const Real dz = dimension_ >= 3 ? az.first[idx[2]] : Real(0);
        const Real d2z = dimension_ >= 3 ? az.second[idx[2]] : Real(0);

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
                                        std::span<Real> values_out,
                                        std::span<Gradient> gradients_out,
                                        std::span<Hessian> hessians_out) const {
    SimplexEval simplex;
    evaluate_simplex(xi, topology_, order_, simplex_exponents_, simplex);
    for (std::size_t i = 0; i < simplex.value.size(); ++i) {
        if (!values_out.empty()) {
            values_out[i] = simplex.value[i];
        }
        if (!gradients_out.empty()) {
            gradients_out[i] = simplex.gradient[i];
        }
        if (!hessians_out.empty()) {
            hessians_out[i] = simplex.hessian[i];
        }
    }
}

// Evaluate wedge bases as triangle/through-axis products.
void LagrangeBasis::evaluate_wedge_to(const Vec3& xi,
                                      std::span<Real> values_out,
                                      std::span<Gradient> gradients_out,
                                      std::span<Hessian> hessians_out) const {
    SimplexEval tri;
    AxisEval z_axis;
    evaluate_simplex(xi, BasisTopology::Triangle, order_, simplex_exponents_, tri);
    evaluate_1d_lagrange(xi[2], nodes_1d_, z_axis);

    for (std::size_t node = 0; node < wedge_indices_.size(); ++node) {
        const auto [tri_idx, z_idx] = wedge_indices_[node];
        const Real tv = tri.value[tri_idx];
        const Real zv = z_axis.value[z_idx];
        const Real dz = z_axis.first[z_idx];
        const Real d2z = z_axis.second[z_idx];

        if (!values_out.empty()) {
            values_out[node] = tv * zv;
        }
        if (!gradients_out.empty()) {
            Gradient& g = gradients_out[node];
            g[0] = tri.gradient[tri_idx][0] * zv;
            g[1] = tri.gradient[tri_idx][1] * zv;
            g[2] = tv * dz;
        }
        if (!hessians_out.empty()) {
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
                                    std::span<Real> values_out,
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

    FE::raise<BasisEvaluationException>(SVMP_HERE,
        "Unsupported element in LagrangeBasis evaluation");
}

void LagrangeBasis::evaluate_values(const Vec3& xi,
                                    std::vector<Real>& values) const {
    values.resize(size());
    evaluate_values_to(xi, std::span<Real>(values.data(), values.size()));
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
                                 std::vector<Real>& values,
                                 std::vector<Gradient>& gradients,
                                 std::vector<Hessian>& hessians) const {
    values.resize(size());
    gradients.resize(size());
    hessians.resize(size());
    evaluate_all_to(xi,
                    std::span<Real>(values.data(), values.size()),
                    std::span<Gradient>(gradients.data(), gradients.size()),
                    std::span<Hessian>(hessians.data(), hessians.size()));
}

void LagrangeBasis::evaluate_values_to(const Vec3& xi,
                                       std::span<Real> values_out) const {
    require_output_span_size(values_out.size(), size(), "LagrangeBasis::evaluate_values_to");
    evaluate_all_to(xi, values_out, std::span<Gradient>{}, std::span<Hessian>{});
}

void LagrangeBasis::evaluate_gradients_to(const Vec3& xi,
                                          std::span<Gradient> gradients_out) const {
    require_output_span_size(gradients_out.size(), size(), "LagrangeBasis::evaluate_gradients_to");
    evaluate_all_to(xi, std::span<Real>{}, gradients_out, std::span<Hessian>{});
}

void LagrangeBasis::evaluate_hessians_to(const Vec3& xi,
                                         std::span<Hessian> hessians_out) const {
    require_output_span_size(hessians_out.size(), size(), "LagrangeBasis::evaluate_hessians_to");
    evaluate_all_to(xi, std::span<Real>{}, std::span<Gradient>{}, hessians_out);
}

} // namespace basis
} // namespace FE
} // namespace svmp
