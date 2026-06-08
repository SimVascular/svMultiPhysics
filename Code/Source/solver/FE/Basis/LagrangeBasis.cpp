/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#include "LagrangeBasis.h"
#include "NodeOrderingConventions.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace svmp {
namespace FE {
namespace basis {

namespace {

using Vec3 = math::Vector<Real, 3>;

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

BasisTopology supported_lagrange_topology(ElementType type) {
    const BasisTopology top = topology(type);
    if (top == BasisTopology::Unknown) {
        throw BasisElementCompatibilityException("LagrangeBasis: unsupported element type",
                                                __FILE__, __LINE__, __func__);
    }
    return top;
}

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
            throw BasisElementCompatibilityException(
                "LagrangeBasis: Quad8 is serendipity; use SerendipityBasis",
                __FILE__, __LINE__, __func__);
        case ElementType::Hex20:
            throw BasisElementCompatibilityException(
                "LagrangeBasis: Hex20 is serendipity; use SerendipityBasis",
                __FILE__, __LINE__, __func__);
        case ElementType::Wedge15:
            throw BasisElementCompatibilityException(
                "LagrangeBasis: Wedge15 is serendipity; use SerendipityBasis",
                __FILE__, __LINE__, __func__);
        case ElementType::Pyramid5:
        case ElementType::Pyramid13:
        case ElementType::Pyramid14:
            throw BasisElementCompatibilityException(
                "LagrangeBasis: pyramid support has been removed from the current solver basis scope",
                __FILE__, __LINE__, __func__);
        default:
            return {element_type, order};
    }
}

std::size_t axis_index_pm_one(Real coord, int order) {
    if (order <= 0) {
        return 0u;
    }
    const Real scaled = (coord + Real(1)) * Real(order) / Real(2);
    return static_cast<std::size_t>(std::llround(scaled));
}

int simplex_lattice_index(Real value, int order) {
    if (order <= 0) {
        return 0;
    }
    return static_cast<int>(std::llround(value * Real(order)));
}

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
        Real denom = Real(1);
        for (std::size_t j = 0; j < n; ++j) {
            if (j != i) {
                denom *= nodes[i] - nodes[j];
            }
        }

        Real value = Real(1);
        for (std::size_t j = 0; j < n; ++j) {
            if (j != i) {
                value *= x - nodes[j];
            }
        }
        out.value[i] = value / denom;

        Real first = Real(0);
        for (std::size_t m = 0; m < n; ++m) {
            if (m == i) {
                continue;
            }
            Real product = Real(1);
            for (std::size_t j = 0; j < n; ++j) {
                if (j != i && j != m) {
                    product *= x - nodes[j];
                }
            }
            first += product;
        }
        out.first[i] = first / denom;

        Real second = Real(0);
        for (std::size_t m = 0; m < n; ++m) {
            if (m == i) {
                continue;
            }
            for (std::size_t l = 0; l < n; ++l) {
                if (l == i || l == m) {
                    continue;
                }
                Real product = Real(1);
                for (std::size_t j = 0; j < n; ++j) {
                    if (j != i && j != m && j != l) {
                        product *= x - nodes[j];
                    }
                }
                second += product;
            }
        }
        out.second[i] = second / denom;
    }
}

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

void evaluate_simplex(const Vec3& xi,
                      BasisTopology top,
                      int order,
                      const std::vector<LagrangeBasis::SimplexExponent>& exponents,
                      SimplexEval& out) {
    const std::size_t n = exponents.size();
    out.value.assign(n, Real(0));
    out.gradient.assign(n, Gradient{});
    out.hessian.assign(n, Hessian{});

    if (n == 1u && order == 0) {
        out.value[0] = Real(1);
        return;
    }

    const int bary_count = top == BasisTopology::Triangle ? 3 : 4;
    std::array<Real, 4> lambda{Real(0), Real(0), Real(0), Real(0)};
    std::array<Gradient, 4> lambda_grad{};

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
        for (int a = 0; a < bary_count; ++a) {
            f[static_cast<std::size_t>(a)] =
                simplex_factor(exponents[i][static_cast<std::size_t>(a)],
                               lambda[static_cast<std::size_t>(a)],
                               order);
        }

        Real value = Real(1);
        for (int a = 0; a < bary_count; ++a) {
            value *= f[static_cast<std::size_t>(a)][0];
        }
        out.value[i] = value;

        for (int a = 0; a < bary_count; ++a) {
            Real product = f[static_cast<std::size_t>(a)][1];
            for (int b = 0; b < bary_count; ++b) {
                if (b != a) {
                    product *= f[static_cast<std::size_t>(b)][0];
                }
            }
            for (std::size_t c = 0; c < 3u; ++c) {
                out.gradient[i][c] += product * lambda_grad[static_cast<std::size_t>(a)][c];
            }
        }

        for (int a = 0; a < bary_count; ++a) {
            for (int b = 0; b < bary_count; ++b) {
                Real product = (a == b)
                    ? f[static_cast<std::size_t>(a)][2]
                    : f[static_cast<std::size_t>(a)][1] *
                      f[static_cast<std::size_t>(b)][1];
                for (int c = 0; c < bary_count; ++c) {
                    if (c != a && c != b) {
                        product *= f[static_cast<std::size_t>(c)][0];
                    }
                }
                for (std::size_t r = 0; r < 3u; ++r) {
                    for (std::size_t c = 0; c < 3u; ++c) {
                        out.hessian[i](r, c) +=
                            product *
                            lambda_grad[static_cast<std::size_t>(a)][r] *
                            lambda_grad[static_cast<std::size_t>(b)][c];
                    }
                }
            }
        }
    }
}

void store_gradient(const Gradient& gradient, Real* dst) {
    dst[0] = gradient[0];
    dst[1] = gradient[1];
    dst[2] = gradient[2];
}

} // namespace

void prewarm_lagrange_basis_scratch(int max_order, std::size_t max_qpts) {
    const auto n = static_cast<std::size_t>(std::max(0, max_order) + 1);
    prewarm_basis_function_scratch(std::max(n * n * n, max_qpts));
}

LagrangeBasis::LagrangeBasis(ElementType type, int order)
    : element_type_(type), order_(order) {
    const auto normalized = normalize_lagrange_request(element_type_, order_);
    element_type_ = normalized.element_type;
    order_ = normalized.order;
    if (order_ < 0) {
        throw BasisConfigurationException("LagrangeBasis requires non-negative polynomial order",
                                          __FILE__, __LINE__, __func__);
    }

    topology_ = supported_lagrange_topology(element_type_);
    dimension_ = reference_dimension(element_type_);
    init_nodes();
}

void LagrangeBasis::init_equispaced_1d_nodes() {
    nodes_1d_.resize(static_cast<std::size_t>(order_ + 1));
    for (int i = 0; i <= order_; ++i) {
        nodes_1d_[static_cast<std::size_t>(i)] =
            equispaced_pm_one_coord(i, order_);
    }
}

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
            build_tensor_product_nodes(1);
            return;
        case BasisTopology::Quadrilateral:
            build_tensor_product_nodes(2);
            return;
        case BasisTopology::Hexahedron:
            build_tensor_product_nodes(3);
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

    throw BasisElementCompatibilityException("Unsupported element type in LagrangeBasis::init_nodes",
                                             __FILE__, __LINE__, __func__);
}

void LagrangeBasis::build_point_nodes() {
    nodes_.push_back(Vec3{Real(0), Real(0), Real(0)});
}

void LagrangeBasis::build_tensor_product_nodes(int dimensions) {
    init_equispaced_1d_nodes();
    nodes_ = ReferenceNodeLayout::get_lagrange_node_coords(element_type_, order_);
    tensor_indices_.reserve(nodes_.size());
    for (const auto& node : nodes_) {
        TensorNodeIndex idx{0u, 0u, 0u};
        idx[0] = axis_index_pm_one(node[0], order_);
        if (dimensions >= 2) {
            idx[1] = axis_index_pm_one(node[1], order_);
        }
        if (dimensions >= 3) {
            idx[2] = axis_index_pm_one(node[2], order_);
        }
        tensor_indices_.push_back(idx);
    }
}

void LagrangeBasis::build_simplex_nodes() {
    nodes_ = ReferenceNodeLayout::get_lagrange_node_coords(element_type_, order_);
    simplex_exponents_.reserve(nodes_.size());
    for (const auto& node : nodes_) {
        simplex_exponents_.push_back(simplex_exponent_from_point(node, topology_, order_));
    }
}

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
        if (it == simplex_exponents_.end()) {
            throw BasisConstructionException("LagrangeBasis: wedge node triangle index lookup failed",
                                             __FILE__, __LINE__, __func__);
        }
        const std::size_t tri_index =
            static_cast<std::size_t>(std::distance(simplex_exponents_.begin(), it));
        wedge_indices_.push_back({tri_index, axis_index_pm_one(node[2], order_)});
    }
}

void LagrangeBasis::evaluate_all_to(const Vec3& xi,
                                    Real* SVMP_RESTRICT values_out,
                                    Real* SVMP_RESTRICT gradients_out,
                                    Real* SVMP_RESTRICT hessians_out) const {
    if (topology_ == BasisTopology::Point) {
        if (values_out) {
            values_out[0] = Real(1);
        }
        if (gradients_out) {
            gradients_out[0] = gradients_out[1] = gradients_out[2] = Real(0);
        }
        if (hessians_out) {
            std::fill_n(hessians_out, 9u, Real(0));
        }
        return;
    }

    if (topology_ == BasisTopology::Line ||
        topology_ == BasisTopology::Quadrilateral ||
        topology_ == BasisTopology::Hexahedron) {
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

            if (values_out) {
                values_out[node] = vx * vy * vz;
            }
            if (gradients_out) {
                Real* g = gradients_out + node * 3u;
                g[0] = dx * vy * vz;
                g[1] = vx * dy * vz;
                g[2] = vx * vy * dz;
            }
            if (hessians_out) {
                Real* h = hessians_out + node * 9u;
                h[0] = d2x * vy * vz;
                h[1] = dx * dy * vz;
                h[2] = dx * vy * dz;
                h[3] = h[1];
                h[4] = vx * d2y * vz;
                h[5] = vx * dy * dz;
                h[6] = h[2];
                h[7] = h[5];
                h[8] = vx * vy * d2z;
            }
        }
        return;
    }

    if (topology_ == BasisTopology::Triangle || topology_ == BasisTopology::Tetrahedron) {
        SimplexEval simplex;
        evaluate_simplex(xi, topology_, order_, simplex_exponents_, simplex);
        for (std::size_t i = 0; i < simplex.value.size(); ++i) {
            if (values_out) {
                values_out[i] = simplex.value[i];
            }
            if (gradients_out) {
                store_gradient(simplex.gradient[i], gradients_out + i * 3u);
            }
            if (hessians_out) {
                store_hessian(simplex.hessian[i], hessians_out + i * 9u);
            }
        }
        return;
    }

    if (topology_ == BasisTopology::Wedge) {
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

            if (values_out) {
                values_out[node] = tv * zv;
            }
            if (gradients_out) {
                Real* g = gradients_out + node * 3u;
                g[0] = tri.gradient[tri_idx][0] * zv;
                g[1] = tri.gradient[tri_idx][1] * zv;
                g[2] = tv * dz;
            }
            if (hessians_out) {
                Real* h = hessians_out + node * 9u;
                const Hessian& th = tri.hessian[tri_idx];
                const Gradient& tg = tri.gradient[tri_idx];
                h[0] = th(0, 0) * zv;
                h[1] = th(0, 1) * zv;
                h[2] = tg[0] * dz;
                h[3] = h[1];
                h[4] = th(1, 1) * zv;
                h[5] = tg[1] * dz;
                h[6] = h[2];
                h[7] = h[5];
                h[8] = tv * d2z;
            }
        }
        return;
    }

    throw BasisEvaluationException("Unsupported element in LagrangeBasis evaluation",
                                   __FILE__, __LINE__, __func__);
}

void LagrangeBasis::evaluate_values(const Vec3& xi,
                                    std::vector<Real>& values) const {
    values.resize(size());
    evaluate_values_to(xi, values.data());
}

void LagrangeBasis::evaluate_gradients(const Vec3& xi,
                                       std::vector<Gradient>& gradients) const {
    gradients.resize(size());
    std::vector<Real> flat(size() * 3u, Real(0));
    evaluate_gradients_to(xi, flat.data());
    for (std::size_t i = 0; i < size(); ++i) {
        gradients[i][0] = flat[i * 3u + 0u];
        gradients[i][1] = flat[i * 3u + 1u];
        gradients[i][2] = flat[i * 3u + 2u];
    }
}

void LagrangeBasis::evaluate_hessians(const Vec3& xi,
                                      std::vector<Hessian>& hessians) const {
    hessians.resize(size());
    std::vector<Real> flat(size() * 9u, Real(0));
    evaluate_hessians_to(xi, flat.data());
    for (std::size_t i = 0; i < size(); ++i) {
        hessians[i] = load_hessian(flat.data() + i * 9u);
    }
}

void LagrangeBasis::evaluate_all(const Vec3& xi,
                                 std::vector<Real>& values,
                                 std::vector<Gradient>& gradients,
                                 std::vector<Hessian>& hessians) const {
    values.resize(size());
    gradients.resize(size());
    hessians.resize(size());
    std::vector<Real> flat_g(size() * 3u, Real(0));
    std::vector<Real> flat_h(size() * 9u, Real(0));
    evaluate_all_to(xi, values.data(), flat_g.data(), flat_h.data());
    for (std::size_t i = 0; i < size(); ++i) {
        gradients[i][0] = flat_g[i * 3u + 0u];
        gradients[i][1] = flat_g[i * 3u + 1u];
        gradients[i][2] = flat_g[i * 3u + 2u];
        hessians[i] = load_hessian(flat_h.data() + i * 9u);
    }
}

void LagrangeBasis::evaluate_values_to(const Vec3& xi,
                                       Real* SVMP_RESTRICT values_out) const {
    evaluate_all_to(xi, values_out, nullptr, nullptr);
}

void LagrangeBasis::evaluate_gradients_to(const Vec3& xi,
                                          Real* SVMP_RESTRICT gradients_out) const {
    evaluate_all_to(xi, nullptr, gradients_out, nullptr);
}

void LagrangeBasis::evaluate_hessians_to(const Vec3& xi,
                                         Real* SVMP_RESTRICT hessians_out) const {
    evaluate_all_to(xi, nullptr, nullptr, hessians_out);
}

} // namespace basis
} // namespace FE
} // namespace svmp
