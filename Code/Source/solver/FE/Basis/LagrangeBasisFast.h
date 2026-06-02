/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#ifndef SVMP_FE_BASIS_LAGRANGEBASISFAST_H
#define SVMP_FE_BASIS_LAGRANGEBASISFAST_H

/**
 * @file LagrangeBasisFast.h
 * @brief Header-only zero-overhead specializations of the Lagrange basis
 *
 * Provides templated static methods for the common nodal Lagrange families
 * with compile-time-known polynomial order. Callers that know their basis
 * type and order at compile time use these directly — there is no virtual
 * dispatch, no std::vector allocation, no scratch lookup, and no topology
 * switch. The output buffers are stack-allocated std::array, sized at
 * compile time. The compiler fully unrolls and constant-folds.
 *
 * These specializations are an alternative entry point to the runtime path
 * provided by `LagrangeBasis`. The runtime path remains the canonical API
 * for generic callers; these specializations serve hot loops that know the
 * element type.
 *
 * Node orderings match `ReferenceNodeLayout::get_lagrange_node_coords(...)` (VTK).
 */

#include "Types.h"
#include "Math/Vector.h"
#include "Math/Matrix.h"
#include <array>
#include <cstddef>

namespace svmp {
namespace FE {
namespace basis {

using Gradient = math::Vector<Real, 3>;
using Hessian  = math::Matrix<Real, 3, 3>;

namespace detail {

constexpr Gradient scaled_gradient(const Gradient& gradient, Real scale) {
    return Gradient{scale * gradient[0], scale * gradient[1], scale * gradient[2]};
}

constexpr Gradient p2_edge_gradient(Real left,
                                    const Gradient& left_gradient,
                                    Real right,
                                    const Gradient& right_gradient) {
    return Gradient{
        Real(4) * (left_gradient[0] * right + right_gradient[0] * left),
        Real(4) * (left_gradient[1] * right + right_gradient[1] * left),
        Real(4) * (left_gradient[2] * right + right_gradient[2] * left),
    };
}

constexpr Hessian p2_vertex_hessian(const Gradient& gradient) {
    Hessian hessian{};
    for (std::size_t row = 0; row < 3u; ++row) {
        for (std::size_t col = 0; col < 3u; ++col) {
            hessian(row, col) = Real(4) * gradient[row] * gradient[col];
        }
    }
    return hessian;
}

constexpr Hessian p2_edge_hessian(const Gradient& left_gradient,
                                  const Gradient& right_gradient) {
    Hessian hessian{};
    for (std::size_t row = 0; row < 3u; ++row) {
        for (std::size_t col = 0; col < 3u; ++col) {
            hessian(row, col) = Real(4) * (
                left_gradient[row] * right_gradient[col] +
                right_gradient[row] * left_gradient[col]);
        }
    }
    return hessian;
}

constexpr std::size_t public_axis_index(int lattice, int order) noexcept {
    return lattice == 0 ? 0u :
           lattice == order ? 1u :
           static_cast<std::size_t>(lattice + 1);
}

template<int Order>
constexpr Real public_axis_coord(std::size_t public_index) noexcept {
    const int lattice = public_index == 0u ? 0 :
                        public_index == 1u ? Order :
                        static_cast<int>(public_index) - 1;
    return Real(-1) + Real(2) * static_cast<Real>(lattice) / static_cast<Real>(Order);
}

template<int Order>
constexpr std::array<Real, Order + 1> make_public_axis_nodes() {
    std::array<Real, Order + 1> nodes{};
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        nodes[i] = public_axis_coord<Order>(i);
    }
    return nodes;
}

template<int Order>
constexpr std::array<Real, Order + 1> make_public_axis_inverse_denominators() {
    constexpr auto nodes = make_public_axis_nodes<Order>();
    std::array<Real, Order + 1> inv_denominators{};
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        Real denominator = Real(1);
        for (std::size_t j = 0; j < nodes.size(); ++j) {
            if (j != i) {
                denominator *= nodes[i] - nodes[j];
            }
        }
        inv_denominators[i] = Real(1) / denominator;
    }
    return inv_denominators;
}

template<int Order, bool NeedFirst, bool NeedSecond>
void fill_axis_lagrange(Real x,
                        std::array<Real, Order + 1>& values,
                        std::array<Real, Order + 1>* first,
                        std::array<Real, Order + 1>* second) {
    constexpr auto nodes = make_public_axis_nodes<Order>();
    constexpr auto inv_denominators = make_public_axis_inverse_denominators<Order>();
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        Real product = Real(1);
        for (std::size_t j = 0; j < nodes.size(); ++j) {
            if (j != i) {
                product *= x - nodes[j];
            }
        }
        values[i] = product * inv_denominators[i];

        if constexpr (NeedFirst) {
            Real derivative = Real(0);
            for (std::size_t m = 0; m < nodes.size(); ++m) {
                if (m == i) {
                    continue;
                }
                Real term = Real(1);
                for (std::size_t j = 0; j < nodes.size(); ++j) {
                    if (j != i && j != m) {
                        term *= x - nodes[j];
                    }
                }
                derivative += term;
            }
            (*first)[i] = derivative * inv_denominators[i];
        }

        if constexpr (NeedSecond) {
            Real curvature = Real(0);
            for (std::size_t m = 0; m < nodes.size(); ++m) {
                if (m == i) {
                    continue;
                }
                for (std::size_t l = 0; l < nodes.size(); ++l) {
                    if (l == i || l == m) {
                        continue;
                    }
                    Real term = Real(1);
                    for (std::size_t j = 0; j < nodes.size(); ++j) {
                        if (j != i && j != m && j != l) {
                            term *= x - nodes[j];
                        }
                    }
                    curvature += term;
                }
            }
            (*second)[i] = curvature * inv_denominators[i];
        }
    }
}

template<int Order>
void fill_axis_values(Real x, std::array<Real, Order + 1>& values) {
    fill_axis_lagrange<Order, false, false>(x, values, nullptr, nullptr);
}

template<int Order>
void fill_axis_values_first(Real x,
                            std::array<Real, Order + 1>& values,
                            std::array<Real, Order + 1>& first) {
    fill_axis_lagrange<Order, true, false>(x, values, &first, nullptr);
}

template<int Order>
void fill_axis_values_first_second(Real x,
                                   std::array<Real, Order + 1>& values,
                                   std::array<Real, Order + 1>& first,
                                   std::array<Real, Order + 1>& second) {
    fill_axis_lagrange<Order, true, true>(x, values, &first, &second);
}

template<int Order>
constexpr std::array<std::array<std::size_t, 2>, (Order + 1) * (Order + 1)>
make_quad_tensor_node_axes() {
    std::array<std::array<std::size_t, 2>, (Order + 1) * (Order + 1)> axes{};
    std::size_t n = 0;

    axes[n++] = {{0u, 0u}};
    axes[n++] = {{1u, 0u}};
    axes[n++] = {{1u, 1u}};
    axes[n++] = {{0u, 1u}};

    for (int i = 1; i < Order; ++i) {
        axes[n++] = {{public_axis_index(i, Order), 0u}};
    }
    for (int j = 1; j < Order; ++j) {
        axes[n++] = {{1u, public_axis_index(j, Order)}};
    }
    for (int i = Order - 1; i >= 1; --i) {
        axes[n++] = {{public_axis_index(i, Order), 1u}};
    }
    for (int j = Order - 1; j >= 1; --j) {
        axes[n++] = {{0u, public_axis_index(j, Order)}};
    }

    for (int j = 1; j < Order; ++j) {
        for (int i = 1; i < Order; ++i) {
            axes[n++] = {{public_axis_index(i, Order), public_axis_index(j, Order)}};
        }
    }

    return axes;
}

template<int Order>
constexpr std::array<std::array<std::size_t, 3>, (Order + 1) * (Order + 1) * (Order + 1)>
make_hex_tensor_node_axes() {
    std::array<std::array<std::size_t, 3>, (Order + 1) * (Order + 1) * (Order + 1)> axes{};
    std::size_t n = 0;

    axes[n++] = {{0u, 0u, 0u}};
    axes[n++] = {{1u, 0u, 0u}};
    axes[n++] = {{1u, 1u, 0u}};
    axes[n++] = {{0u, 1u, 0u}};
    axes[n++] = {{0u, 0u, 1u}};
    axes[n++] = {{1u, 0u, 1u}};
    axes[n++] = {{1u, 1u, 1u}};
    axes[n++] = {{0u, 1u, 1u}};

    for (int i = 1; i < Order; ++i) {
        axes[n++] = {{public_axis_index(i, Order), 0u, 0u}};
    }
    for (int j = 1; j < Order; ++j) {
        axes[n++] = {{1u, public_axis_index(j, Order), 0u}};
    }
    for (int i = Order - 1; i >= 1; --i) {
        axes[n++] = {{public_axis_index(i, Order), 1u, 0u}};
    }
    for (int j = Order - 1; j >= 1; --j) {
        axes[n++] = {{0u, public_axis_index(j, Order), 0u}};
    }
    for (int i = 1; i < Order; ++i) {
        axes[n++] = {{public_axis_index(i, Order), 0u, 1u}};
    }
    for (int j = 1; j < Order; ++j) {
        axes[n++] = {{1u, public_axis_index(j, Order), 1u}};
    }
    for (int i = Order - 1; i >= 1; --i) {
        axes[n++] = {{public_axis_index(i, Order), 1u, 1u}};
    }
    for (int j = Order - 1; j >= 1; --j) {
        axes[n++] = {{0u, public_axis_index(j, Order), 1u}};
    }
    for (int k = 1; k < Order; ++k) {
        axes[n++] = {{0u, 0u, public_axis_index(k, Order)}};
    }
    for (int k = 1; k < Order; ++k) {
        axes[n++] = {{1u, 0u, public_axis_index(k, Order)}};
    }
    for (int k = 1; k < Order; ++k) {
        axes[n++] = {{1u, 1u, public_axis_index(k, Order)}};
    }
    for (int k = 1; k < Order; ++k) {
        axes[n++] = {{0u, 1u, public_axis_index(k, Order)}};
    }

    for (int j = 1; j < Order; ++j) {
        for (int i = 1; i < Order; ++i) {
            axes[n++] = {{public_axis_index(i, Order), public_axis_index(j, Order), 0u}};
        }
    }
    for (int j = 1; j < Order; ++j) {
        for (int i = 1; i < Order; ++i) {
            axes[n++] = {{public_axis_index(i, Order), public_axis_index(j, Order), 1u}};
        }
    }
    for (int k = 1; k < Order; ++k) {
        for (int i = 1; i < Order; ++i) {
            axes[n++] = {{public_axis_index(i, Order), 0u, public_axis_index(k, Order)}};
        }
    }
    for (int k = 1; k < Order; ++k) {
        for (int j = 1; j < Order; ++j) {
            axes[n++] = {{1u, public_axis_index(j, Order), public_axis_index(k, Order)}};
        }
    }
    for (int k = 1; k < Order; ++k) {
        for (int i = Order - 1; i >= 1; --i) {
            axes[n++] = {{public_axis_index(i, Order), 1u, public_axis_index(k, Order)}};
        }
    }
    for (int k = 1; k < Order; ++k) {
        for (int j = Order - 1; j >= 1; --j) {
            axes[n++] = {{0u, public_axis_index(j, Order), public_axis_index(k, Order)}};
        }
    }

    for (int k = 1; k < Order; ++k) {
        for (int j = 1; j < Order; ++j) {
            for (int i = 1; i < Order; ++i) {
                axes[n++] = {{public_axis_index(i, Order),
                              public_axis_index(j, Order),
                              public_axis_index(k, Order)}};
            }
        }
    }

    return axes;
}

template<int Order>
constexpr std::array<std::array<std::size_t, 3>, (Order + 1) * (Order + 2) / 2>
make_triangle_simplex_exponents() {
    std::array<std::array<std::size_t, 3>, (Order + 1) * (Order + 2) / 2> exponents{};
    std::size_t n = 0;

    exponents[n++] = {{static_cast<std::size_t>(Order), 0u, 0u}};
    exponents[n++] = {{0u, static_cast<std::size_t>(Order), 0u}};
    exponents[n++] = {{0u, 0u, static_cast<std::size_t>(Order)}};

    for (int m = 1; m < Order; ++m) {
        exponents[n++] = {{static_cast<std::size_t>(Order - m), static_cast<std::size_t>(m), 0u}};
    }
    for (int m = 1; m < Order; ++m) {
        exponents[n++] = {{0u, static_cast<std::size_t>(Order - m), static_cast<std::size_t>(m)}};
    }
    for (int m = 1; m < Order; ++m) {
        exponents[n++] = {{static_cast<std::size_t>(m), 0u, static_cast<std::size_t>(Order - m)}};
    }

    for (int c = 1; c <= Order - 2; ++c) {
        for (int b = 1; b <= Order - c - 1; ++b) {
            const int a = Order - b - c;
            exponents[n++] = {{static_cast<std::size_t>(a),
                               static_cast<std::size_t>(b),
                               static_cast<std::size_t>(c)}};
        }
    }

    return exponents;
}

template<int Order>
constexpr std::array<std::array<std::size_t, 4>, (Order + 1) * (Order + 2) * (Order + 3) / 6>
make_tetrahedron_simplex_exponents() {
    std::array<std::array<std::size_t, 4>, (Order + 1) * (Order + 2) * (Order + 3) / 6> exponents{};
    std::size_t n = 0;

    exponents[n++] = {{static_cast<std::size_t>(Order), 0u, 0u, 0u}};
    exponents[n++] = {{0u, static_cast<std::size_t>(Order), 0u, 0u}};
    exponents[n++] = {{0u, 0u, static_cast<std::size_t>(Order), 0u}};
    exponents[n++] = {{0u, 0u, 0u, static_cast<std::size_t>(Order)}};

    constexpr int edges[6][2] = {
        {0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}
    };
    for (const auto& edge : edges) {
        for (int m = 1; m < Order; ++m) {
            std::array<std::size_t, 4> e{};
            e[static_cast<std::size_t>(edge[0])] = static_cast<std::size_t>(Order - m);
            e[static_cast<std::size_t>(edge[1])] = static_cast<std::size_t>(m);
            exponents[n++] = e;
        }
    }

    constexpr int faces[4][3] = {
        {0, 1, 2},
        {0, 1, 3},
        {1, 2, 3},
        {0, 2, 3},
    };
    for (const auto& face : faces) {
        for (int c = 1; c <= Order - 2; ++c) {
            for (int b = 1; b <= Order - c - 1; ++b) {
                const int a = Order - b - c;
                std::array<std::size_t, 4> e{};
                e[static_cast<std::size_t>(face[0])] = static_cast<std::size_t>(a);
                e[static_cast<std::size_t>(face[1])] = static_cast<std::size_t>(b);
                e[static_cast<std::size_t>(face[2])] = static_cast<std::size_t>(c);
                exponents[n++] = e;
            }
        }
    }

    for (int l = 1; l <= Order - 3; ++l) {
        for (int k = 1; k <= Order - l - 2; ++k) {
            for (int j = 1; j <= Order - l - k - 1; ++j) {
                const int i = Order - j - k - l;
                exponents[n++] = {{static_cast<std::size_t>(i),
                                   static_cast<std::size_t>(j),
                                   static_cast<std::size_t>(k),
                                   static_cast<std::size_t>(l)}};
            }
        }
    }

    return exponents;
}

template<int Order, bool NeedFirst, bool NeedSecond>
void fill_simplex_factor_sequence(Real lambda,
                                  std::array<Real, Order + 1>& phi,
                                  std::array<Real, Order + 1>* dphi,
                                  std::array<Real, Order + 1>* d2phi) {
    phi[0] = Real(1);
    if constexpr (NeedFirst) {
        (*dphi)[0] = Real(0);
    }
    if constexpr (NeedSecond) {
        (*d2phi)[0] = Real(0);
    }

    const Real t = static_cast<Real>(Order) * lambda;
    constexpr Real dt_dlambda = static_cast<Real>(Order);
    Real dphi_dt_prev = Real(0);
    Real d2phi_dt2_prev = Real(0);

    for (int a = 1; a <= Order; ++a) {
        const std::size_t au = static_cast<std::size_t>(a);
        const Real inv_a = Real(1) / static_cast<Real>(a);
        const Real s = (t - static_cast<Real>(a - 1)) * inv_a;
        phi[au] = s * phi[au - 1];

        if constexpr (NeedFirst) {
            const Real dphi_dt = inv_a * phi[au - 1] + s * dphi_dt_prev;
            (*dphi)[au] = dt_dlambda * dphi_dt;

            if constexpr (NeedSecond) {
                const Real d2phi_dt2 = Real(2) * inv_a * dphi_dt_prev + s * d2phi_dt2_prev;
                (*d2phi)[au] = dt_dlambda * dt_dlambda * d2phi_dt2;
                d2phi_dt2_prev = d2phi_dt2;
            }

            dphi_dt_prev = dphi_dt;
        }
    }
}

template<int Order>
void fill_simplex_factor_values(Real lambda, std::array<Real, Order + 1>& phi) {
    fill_simplex_factor_sequence<Order, false, false>(lambda, phi, nullptr, nullptr);
}

template<int Order>
void fill_simplex_factor_values_first(Real lambda,
                                      std::array<Real, Order + 1>& phi,
                                      std::array<Real, Order + 1>& dphi) {
    fill_simplex_factor_sequence<Order, true, false>(lambda, phi, &dphi, nullptr);
}

template<int Order>
void fill_simplex_factor_values_first_second(Real lambda,
                                             std::array<Real, Order + 1>& phi,
                                             std::array<Real, Order + 1>& dphi,
                                             std::array<Real, Order + 1>& d2phi) {
    fill_simplex_factor_sequence<Order, true, true>(lambda, phi, &dphi, &d2phi);
}

} // namespace detail

// ---------------------------------------------------------------------------
// LagrangeLineFast<Order>
// ---------------------------------------------------------------------------
template<int Order>
struct LagrangeLineFast;

template<>
struct LagrangeLineFast<1> {
    static constexpr int n_dofs = 2;

    static constexpr void evaluate(const math::Vector<Real, 3>& xi, std::array<Real, n_dofs>& out) {
        out[0] = (Real(1) - xi[0]) * Real(0.5);
        out[1] = (Real(1) + xi[0]) * Real(0.5);
    }

    static constexpr void evaluate_gradients(const math::Vector<Real, 3>& /*xi*/,
                                             std::array<Gradient, n_dofs>& out) {
        out[0] = Gradient{Real(-0.5), Real(0), Real(0)};
        out[1] = Gradient{Real( 0.5), Real(0), Real(0)};
    }

    static constexpr void evaluate_hessians(const math::Vector<Real, 3>& /*xi*/,
                                            std::array<Hessian, n_dofs>& out) {
        out[0] = Hessian{};
        out[1] = Hessian{};
    }
};

template<>
struct LagrangeLineFast<2> {
    static constexpr int n_dofs = 3;

    static constexpr void evaluate(const math::Vector<Real, 3>& xi, std::array<Real, n_dofs>& out) {
        const Real x = xi[0];
        out[0] = x * (x - Real(1)) * Real(0.5);
        out[1] = x * (x + Real(1)) * Real(0.5);
        out[2] = (Real(1) - x) * (Real(1) + x);
    }

    static constexpr void evaluate_gradients(const math::Vector<Real, 3>& xi,
                                             std::array<Gradient, n_dofs>& out) {
        const Real x = xi[0];
        out[0] = Gradient{x - Real(0.5), Real(0), Real(0)};
        out[1] = Gradient{x + Real(0.5), Real(0), Real(0)};
        out[2] = Gradient{Real(-2) * x, Real(0), Real(0)};
    }

    static constexpr void evaluate_hessians(const math::Vector<Real, 3>& /*xi*/,
                                            std::array<Hessian, n_dofs>& out) {
        out[0] = Hessian{};
        out[1] = Hessian{};
        out[2] = Hessian{};
        out[0](0, 0) = Real(1);
        out[1](0, 0) = Real(1);
        out[2](0, 0) = Real(-2);
    }
};

template<>
struct LagrangeLineFast<3> {
    static constexpr int n_dofs = 4;

    static void evaluate(const math::Vector<Real, 3>& xi, std::array<Real, n_dofs>& out) {
        detail::fill_axis_values<3>(xi[0], out);
    }

    static void evaluate_gradients(const math::Vector<Real, 3>& xi,
                                   std::array<Gradient, n_dofs>& out) {
        std::array<Real, n_dofs> values{};
        std::array<Real, n_dofs> first{};
        detail::fill_axis_values_first<3>(xi[0], values, first);
        for (std::size_t i = 0; i < first.size(); ++i) {
            out[i] = Gradient{first[i], Real(0), Real(0)};
        }
    }

    static void evaluate_hessians(const math::Vector<Real, 3>& xi,
                                  std::array<Hessian, n_dofs>& out) {
        std::array<Real, n_dofs> values{};
        std::array<Real, n_dofs> first{};
        std::array<Real, n_dofs> second{};
        detail::fill_axis_values_first_second<3>(xi[0], values, first, second);
        for (std::size_t i = 0; i < second.size(); ++i) {
            Hessian H{};
            H(0, 0) = second[i];
            out[i] = H;
        }
    }
};

// ---------------------------------------------------------------------------
// LagrangeQuadFast<Order>
// ---------------------------------------------------------------------------
template<int Order>
struct LagrangeQuadFast;

template<>
struct LagrangeQuadFast<1> {
    static constexpr int n_dofs = 4;

    // VTK Quad4 corner ordering: (-,-), (+,-), (+,+), (-,+).
    static constexpr void evaluate(const math::Vector<Real, 3>& xi, std::array<Real, n_dofs>& out) {
        const Real lx = (Real(1) - xi[0]) * Real(0.5);
        const Real ly = (Real(1) - xi[1]) * Real(0.5);
        const Real ux = (Real(1) + xi[0]) * Real(0.5);
        const Real uy = (Real(1) + xi[1]) * Real(0.5);
        out[0] = lx * ly;
        out[1] = ux * ly;
        out[2] = ux * uy;
        out[3] = lx * uy;
    }

    static constexpr void evaluate_gradients(const math::Vector<Real, 3>& xi,
                                             std::array<Gradient, n_dofs>& out) {
        const Real lx = (Real(1) - xi[0]) * Real(0.5);
        const Real ly = (Real(1) - xi[1]) * Real(0.5);
        const Real ux = (Real(1) + xi[0]) * Real(0.5);
        const Real uy = (Real(1) + xi[1]) * Real(0.5);
        out[0] = Gradient{Real(-0.5) * ly, Real(-0.5) * lx, Real(0)};
        out[1] = Gradient{Real( 0.5) * ly, Real(-0.5) * ux, Real(0)};
        out[2] = Gradient{Real( 0.5) * uy, Real( 0.5) * ux, Real(0)};
        out[3] = Gradient{Real(-0.5) * uy, Real( 0.5) * lx, Real(0)};
    }

    static constexpr void evaluate_hessians(const math::Vector<Real, 3>& /*xi*/,
                                            std::array<Hessian, n_dofs>& out) {
        out[0] = Hessian{};
        out[1] = Hessian{};
        out[2] = Hessian{};
        out[3] = Hessian{};
        constexpr Real qrt = Real(0.25);
        out[0](0, 1) = qrt;  out[0](1, 0) = qrt;
        out[1](0, 1) = -qrt; out[1](1, 0) = -qrt;
        out[2](0, 1) = qrt;  out[2](1, 0) = qrt;
        out[3](0, 1) = -qrt; out[3](1, 0) = -qrt;
    }
};

template<>
struct LagrangeQuadFast<2> {
    static constexpr int n_dofs = 9;

    static constexpr std::array<std::array<std::size_t, 2>, n_dofs> node_axes = {{
        {{0u, 0u}}, {{1u, 0u}}, {{1u, 1u}}, {{0u, 1u}},
        {{2u, 0u}}, {{1u, 2u}}, {{2u, 1u}}, {{0u, 2u}},
        {{2u, 2u}},
    }};

    static void evaluate(const math::Vector<Real, 3>& xi, std::array<Real, n_dofs>& out) {
        std::array<Real, LagrangeLineFast<2>::n_dofs> lx{};
        std::array<Real, LagrangeLineFast<2>::n_dofs> ly{};
        LagrangeLineFast<2>::evaluate({xi[0], Real(0), Real(0)}, lx);
        LagrangeLineFast<2>::evaluate({xi[1], Real(0), Real(0)}, ly);
        for (std::size_t n = 0; n < node_axes.size(); ++n) {
            out[n] = lx[node_axes[n][0]] * ly[node_axes[n][1]];
        }
    }

    static void evaluate_gradients(const math::Vector<Real, 3>& xi,
                                   std::array<Gradient, n_dofs>& out) {
        std::array<Real, LagrangeLineFast<2>::n_dofs> lx{};
        std::array<Real, LagrangeLineFast<2>::n_dofs> ly{};
        std::array<Gradient, LagrangeLineFast<2>::n_dofs> gx{};
        std::array<Gradient, LagrangeLineFast<2>::n_dofs> gy{};
        LagrangeLineFast<2>::evaluate({xi[0], Real(0), Real(0)}, lx);
        LagrangeLineFast<2>::evaluate({xi[1], Real(0), Real(0)}, ly);
        LagrangeLineFast<2>::evaluate_gradients({xi[0], Real(0), Real(0)}, gx);
        LagrangeLineFast<2>::evaluate_gradients({xi[1], Real(0), Real(0)}, gy);
        for (std::size_t n = 0; n < node_axes.size(); ++n) {
            const auto i = node_axes[n][0];
            const auto j = node_axes[n][1];
            out[n] = Gradient{gx[i][0] * ly[j], lx[i] * gy[j][0], Real(0)};
        }
    }

    static void evaluate_hessians(const math::Vector<Real, 3>& xi,
                                  std::array<Hessian, n_dofs>& out) {
        std::array<Real, LagrangeLineFast<2>::n_dofs> lx{};
        std::array<Real, LagrangeLineFast<2>::n_dofs> ly{};
        std::array<Gradient, LagrangeLineFast<2>::n_dofs> gx{};
        std::array<Gradient, LagrangeLineFast<2>::n_dofs> gy{};
        std::array<Hessian, LagrangeLineFast<2>::n_dofs> hx{};
        std::array<Hessian, LagrangeLineFast<2>::n_dofs> hy{};
        LagrangeLineFast<2>::evaluate({xi[0], Real(0), Real(0)}, lx);
        LagrangeLineFast<2>::evaluate({xi[1], Real(0), Real(0)}, ly);
        LagrangeLineFast<2>::evaluate_gradients({xi[0], Real(0), Real(0)}, gx);
        LagrangeLineFast<2>::evaluate_gradients({xi[1], Real(0), Real(0)}, gy);
        LagrangeLineFast<2>::evaluate_hessians({xi[0], Real(0), Real(0)}, hx);
        LagrangeLineFast<2>::evaluate_hessians({xi[1], Real(0), Real(0)}, hy);
        for (std::size_t n = 0; n < node_axes.size(); ++n) {
            const auto i = node_axes[n][0];
            const auto j = node_axes[n][1];
            Hessian H{};
            H(0, 0) = hx[i](0, 0) * ly[j];
            H(1, 1) = lx[i] * hy[j](0, 0);
            H(0, 1) = gx[i][0] * gy[j][0];
            H(1, 0) = H(0, 1);
            out[n] = H;
        }
    }
};

template<>
struct LagrangeQuadFast<3> {
    static constexpr int n_dofs = 16;

    static constexpr std::array<std::array<std::size_t, 2>, n_dofs> node_axes =
        detail::make_quad_tensor_node_axes<3>();

    static void evaluate(const math::Vector<Real, 3>& xi, std::array<Real, n_dofs>& out) {
        std::array<Real, LagrangeLineFast<3>::n_dofs> lx{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> ly{};
        detail::fill_axis_values<3>(xi[0], lx);
        detail::fill_axis_values<3>(xi[1], ly);
        for (std::size_t n = 0; n < node_axes.size(); ++n) {
            out[n] = lx[node_axes[n][0]] * ly[node_axes[n][1]];
        }
    }

    static void evaluate_gradients(const math::Vector<Real, 3>& xi,
                                   std::array<Gradient, n_dofs>& out) {
        std::array<Real, LagrangeLineFast<3>::n_dofs> lx{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> ly{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> gx{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> gy{};
        detail::fill_axis_values_first<3>(xi[0], lx, gx);
        detail::fill_axis_values_first<3>(xi[1], ly, gy);
        for (std::size_t n = 0; n < node_axes.size(); ++n) {
            const auto i = node_axes[n][0];
            const auto j = node_axes[n][1];
            out[n] = Gradient{gx[i] * ly[j], lx[i] * gy[j], Real(0)};
        }
    }

    static void evaluate_hessians(const math::Vector<Real, 3>& xi,
                                  std::array<Hessian, n_dofs>& out) {
        std::array<Real, LagrangeLineFast<3>::n_dofs> lx{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> ly{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> gx{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> gy{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> hx{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> hy{};
        detail::fill_axis_values_first_second<3>(xi[0], lx, gx, hx);
        detail::fill_axis_values_first_second<3>(xi[1], ly, gy, hy);
        for (std::size_t n = 0; n < node_axes.size(); ++n) {
            const auto i = node_axes[n][0];
            const auto j = node_axes[n][1];
            Hessian H{};
            H(0, 0) = hx[i] * ly[j];
            H(1, 1) = lx[i] * hy[j];
            H(0, 1) = gx[i] * gy[j];
            H(1, 0) = H(0, 1);
            out[n] = H;
        }
    }
};

// ---------------------------------------------------------------------------
// LagrangeHexFast<Order>
// ---------------------------------------------------------------------------
template<int Order>
struct LagrangeHexFast;

template<>
struct LagrangeHexFast<1> {
    static constexpr int n_dofs = 8;

    // VTK Hex8 corner ordering: (-,-,-), (+,-,-), (+,+,-), (-,+,-),
    //                           (-,-,+), (+,-,+), (+,+,+), (-,+,+).
    static constexpr void evaluate(const math::Vector<Real, 3>& xi, std::array<Real, n_dofs>& out) {
        const Real lx = (Real(1) - xi[0]) * Real(0.5);
        const Real ly = (Real(1) - xi[1]) * Real(0.5);
        const Real lz = (Real(1) - xi[2]) * Real(0.5);
        const Real ux = (Real(1) + xi[0]) * Real(0.5);
        const Real uy = (Real(1) + xi[1]) * Real(0.5);
        const Real uz = (Real(1) + xi[2]) * Real(0.5);
        // Precompute z-plane partial products (sum factorization).
        const Real lxly = lx * ly;
        const Real uxly = ux * ly;
        const Real uxuy = ux * uy;
        const Real lxuy = lx * uy;
        out[0] = lxly * lz;
        out[1] = uxly * lz;
        out[2] = uxuy * lz;
        out[3] = lxuy * lz;
        out[4] = lxly * uz;
        out[5] = uxly * uz;
        out[6] = uxuy * uz;
        out[7] = lxuy * uz;
    }

    static constexpr void evaluate_gradients(const math::Vector<Real, 3>& xi,
                                             std::array<Gradient, n_dofs>& out) {
        const Real lx = (Real(1) - xi[0]) * Real(0.5);
        const Real ly = (Real(1) - xi[1]) * Real(0.5);
        const Real lz = (Real(1) - xi[2]) * Real(0.5);
        const Real ux = (Real(1) + xi[0]) * Real(0.5);
        const Real uy = (Real(1) + xi[1]) * Real(0.5);
        const Real uz = (Real(1) + xi[2]) * Real(0.5);
        // dL_0(x)/dx = -0.5, dL_1(x)/dx = +0.5 along each axis.
        out[0] = Gradient{Real(-0.5) * ly * lz, Real(-0.5) * lx * lz, Real(-0.5) * lx * ly};
        out[1] = Gradient{Real( 0.5) * ly * lz, Real(-0.5) * ux * lz, Real(-0.5) * ux * ly};
        out[2] = Gradient{Real( 0.5) * uy * lz, Real( 0.5) * ux * lz, Real(-0.5) * ux * uy};
        out[3] = Gradient{Real(-0.5) * uy * lz, Real( 0.5) * lx * lz, Real(-0.5) * lx * uy};
        out[4] = Gradient{Real(-0.5) * ly * uz, Real(-0.5) * lx * uz, Real( 0.5) * lx * ly};
        out[5] = Gradient{Real( 0.5) * ly * uz, Real(-0.5) * ux * uz, Real( 0.5) * ux * ly};
        out[6] = Gradient{Real( 0.5) * uy * uz, Real( 0.5) * ux * uz, Real( 0.5) * ux * uy};
        out[7] = Gradient{Real(-0.5) * uy * uz, Real( 0.5) * lx * uz, Real( 0.5) * lx * uy};
    }

    static constexpr void evaluate_hessians(const math::Vector<Real, 3>& xi,
                                            std::array<Hessian, n_dofs>& out) {
        const Real lx = (Real(1) - xi[0]) * Real(0.5);
        const Real ly = (Real(1) - xi[1]) * Real(0.5);
        const Real lz = (Real(1) - xi[2]) * Real(0.5);
        const Real ux = (Real(1) + xi[0]) * Real(0.5);
        const Real uy = (Real(1) + xi[1]) * Real(0.5);
        const Real uz = (Real(1) + xi[2]) * Real(0.5);
        const Real ax[8] = {lx, ux, ux, lx, lx, ux, ux, lx};
        const Real ay[8] = {ly, ly, uy, uy, ly, ly, uy, uy};
        const Real az[8] = {lz, lz, lz, lz, uz, uz, uz, uz};
        const int sx[8] = {-1, 1, 1, -1, -1, 1, 1, -1};
        const int sy[8] = {-1, -1, 1, 1, -1, -1, 1, 1};
        const int sz[8] = {-1, -1, -1, -1, 1, 1, 1, 1};
        constexpr Real qrt = Real(0.25);
        for (std::size_t n = 0; n < static_cast<std::size_t>(n_dofs); ++n) {
            out[n] = Hessian{};
            out[n](0, 1) = static_cast<Real>(sx[n] * sy[n]) * qrt * az[n];
            out[n](1, 0) = out[n](0, 1);
            out[n](0, 2) = static_cast<Real>(sx[n] * sz[n]) * qrt * ay[n];
            out[n](2, 0) = out[n](0, 2);
            out[n](1, 2) = static_cast<Real>(sy[n] * sz[n]) * qrt * ax[n];
            out[n](2, 1) = out[n](1, 2);
        }
    }
};

template<>
struct LagrangeHexFast<2> {
    static constexpr int n_dofs = 27;

    static constexpr std::array<std::array<std::size_t, 3>, n_dofs> node_axes = {{
        {{0u, 0u, 0u}}, {{1u, 0u, 0u}}, {{1u, 1u, 0u}}, {{0u, 1u, 0u}},
        {{0u, 0u, 1u}}, {{1u, 0u, 1u}}, {{1u, 1u, 1u}}, {{0u, 1u, 1u}},
        {{2u, 0u, 0u}}, {{1u, 2u, 0u}}, {{2u, 1u, 0u}}, {{0u, 2u, 0u}},
        {{2u, 0u, 1u}}, {{1u, 2u, 1u}}, {{2u, 1u, 1u}}, {{0u, 2u, 1u}},
        {{0u, 0u, 2u}}, {{1u, 0u, 2u}}, {{1u, 1u, 2u}}, {{0u, 1u, 2u}},
        {{2u, 2u, 0u}}, {{2u, 2u, 1u}}, {{2u, 0u, 2u}}, {{1u, 2u, 2u}},
        {{2u, 1u, 2u}}, {{0u, 2u, 2u}}, {{2u, 2u, 2u}},
    }};

    static void evaluate(const math::Vector<Real, 3>& xi, std::array<Real, n_dofs>& out) {
        std::array<Real, LagrangeLineFast<2>::n_dofs> lx{};
        std::array<Real, LagrangeLineFast<2>::n_dofs> ly{};
        std::array<Real, LagrangeLineFast<2>::n_dofs> lz{};
        LagrangeLineFast<2>::evaluate({xi[0], Real(0), Real(0)}, lx);
        LagrangeLineFast<2>::evaluate({xi[1], Real(0), Real(0)}, ly);
        LagrangeLineFast<2>::evaluate({xi[2], Real(0), Real(0)}, lz);
        for (std::size_t n = 0; n < node_axes.size(); ++n) {
            out[n] = lx[node_axes[n][0]] * ly[node_axes[n][1]] * lz[node_axes[n][2]];
        }
    }

    static void evaluate_gradients(const math::Vector<Real, 3>& xi,
                                   std::array<Gradient, n_dofs>& out) {
        std::array<Real, LagrangeLineFast<2>::n_dofs> lx{};
        std::array<Real, LagrangeLineFast<2>::n_dofs> ly{};
        std::array<Real, LagrangeLineFast<2>::n_dofs> lz{};
        std::array<Gradient, LagrangeLineFast<2>::n_dofs> gx{};
        std::array<Gradient, LagrangeLineFast<2>::n_dofs> gy{};
        std::array<Gradient, LagrangeLineFast<2>::n_dofs> gz{};
        LagrangeLineFast<2>::evaluate({xi[0], Real(0), Real(0)}, lx);
        LagrangeLineFast<2>::evaluate({xi[1], Real(0), Real(0)}, ly);
        LagrangeLineFast<2>::evaluate({xi[2], Real(0), Real(0)}, lz);
        LagrangeLineFast<2>::evaluate_gradients({xi[0], Real(0), Real(0)}, gx);
        LagrangeLineFast<2>::evaluate_gradients({xi[1], Real(0), Real(0)}, gy);
        LagrangeLineFast<2>::evaluate_gradients({xi[2], Real(0), Real(0)}, gz);
        for (std::size_t n = 0; n < node_axes.size(); ++n) {
            const auto i = node_axes[n][0];
            const auto j = node_axes[n][1];
            const auto k = node_axes[n][2];
            out[n] = Gradient{
                gx[i][0] * ly[j] * lz[k],
                lx[i] * gy[j][0] * lz[k],
                lx[i] * ly[j] * gz[k][0],
            };
        }
    }

    static void evaluate_hessians(const math::Vector<Real, 3>& xi,
                                  std::array<Hessian, n_dofs>& out) {
        std::array<Real, LagrangeLineFast<2>::n_dofs> lx{};
        std::array<Real, LagrangeLineFast<2>::n_dofs> ly{};
        std::array<Real, LagrangeLineFast<2>::n_dofs> lz{};
        std::array<Gradient, LagrangeLineFast<2>::n_dofs> gx{};
        std::array<Gradient, LagrangeLineFast<2>::n_dofs> gy{};
        std::array<Gradient, LagrangeLineFast<2>::n_dofs> gz{};
        std::array<Hessian, LagrangeLineFast<2>::n_dofs> hx{};
        std::array<Hessian, LagrangeLineFast<2>::n_dofs> hy{};
        std::array<Hessian, LagrangeLineFast<2>::n_dofs> hz{};
        LagrangeLineFast<2>::evaluate({xi[0], Real(0), Real(0)}, lx);
        LagrangeLineFast<2>::evaluate({xi[1], Real(0), Real(0)}, ly);
        LagrangeLineFast<2>::evaluate({xi[2], Real(0), Real(0)}, lz);
        LagrangeLineFast<2>::evaluate_gradients({xi[0], Real(0), Real(0)}, gx);
        LagrangeLineFast<2>::evaluate_gradients({xi[1], Real(0), Real(0)}, gy);
        LagrangeLineFast<2>::evaluate_gradients({xi[2], Real(0), Real(0)}, gz);
        LagrangeLineFast<2>::evaluate_hessians({xi[0], Real(0), Real(0)}, hx);
        LagrangeLineFast<2>::evaluate_hessians({xi[1], Real(0), Real(0)}, hy);
        LagrangeLineFast<2>::evaluate_hessians({xi[2], Real(0), Real(0)}, hz);
        for (std::size_t n = 0; n < node_axes.size(); ++n) {
            const auto i = node_axes[n][0];
            const auto j = node_axes[n][1];
            const auto k = node_axes[n][2];
            Hessian H{};
            H(0, 0) = hx[i](0, 0) * ly[j] * lz[k];
            H(1, 1) = lx[i] * hy[j](0, 0) * lz[k];
            H(2, 2) = lx[i] * ly[j] * hz[k](0, 0);
            H(0, 1) = gx[i][0] * gy[j][0] * lz[k];
            H(1, 0) = H(0, 1);
            H(0, 2) = gx[i][0] * ly[j] * gz[k][0];
            H(2, 0) = H(0, 2);
            H(1, 2) = lx[i] * gy[j][0] * gz[k][0];
            H(2, 1) = H(1, 2);
            out[n] = H;
        }
    }
};

template<>
struct LagrangeHexFast<3> {
    static constexpr int n_dofs = 64;

    static constexpr std::array<std::array<std::size_t, 3>, n_dofs> node_axes =
        detail::make_hex_tensor_node_axes<3>();

    static void evaluate(const math::Vector<Real, 3>& xi, std::array<Real, n_dofs>& out) {
        std::array<Real, LagrangeLineFast<3>::n_dofs> lx{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> ly{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> lz{};
        detail::fill_axis_values<3>(xi[0], lx);
        detail::fill_axis_values<3>(xi[1], ly);
        detail::fill_axis_values<3>(xi[2], lz);
        for (std::size_t n = 0; n < node_axes.size(); ++n) {
            out[n] = lx[node_axes[n][0]] * ly[node_axes[n][1]] * lz[node_axes[n][2]];
        }
    }

    static void evaluate_gradients(const math::Vector<Real, 3>& xi,
                                   std::array<Gradient, n_dofs>& out) {
        std::array<Real, LagrangeLineFast<3>::n_dofs> lx{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> ly{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> lz{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> gx{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> gy{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> gz{};
        detail::fill_axis_values_first<3>(xi[0], lx, gx);
        detail::fill_axis_values_first<3>(xi[1], ly, gy);
        detail::fill_axis_values_first<3>(xi[2], lz, gz);
        for (std::size_t n = 0; n < node_axes.size(); ++n) {
            const auto i = node_axes[n][0];
            const auto j = node_axes[n][1];
            const auto k = node_axes[n][2];
            out[n] = Gradient{
                gx[i] * ly[j] * lz[k],
                lx[i] * gy[j] * lz[k],
                lx[i] * ly[j] * gz[k],
            };
        }
    }

    static void evaluate_hessians(const math::Vector<Real, 3>& xi,
                                  std::array<Hessian, n_dofs>& out) {
        std::array<Real, LagrangeLineFast<3>::n_dofs> lx{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> ly{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> lz{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> gx{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> gy{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> gz{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> hx{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> hy{};
        std::array<Real, LagrangeLineFast<3>::n_dofs> hz{};
        detail::fill_axis_values_first_second<3>(xi[0], lx, gx, hx);
        detail::fill_axis_values_first_second<3>(xi[1], ly, gy, hy);
        detail::fill_axis_values_first_second<3>(xi[2], lz, gz, hz);
        for (std::size_t n = 0; n < node_axes.size(); ++n) {
            const auto i = node_axes[n][0];
            const auto j = node_axes[n][1];
            const auto k = node_axes[n][2];
            Hessian H{};
            H(0, 0) = hx[i] * ly[j] * lz[k];
            H(1, 1) = lx[i] * hy[j] * lz[k];
            H(2, 2) = lx[i] * ly[j] * hz[k];
            H(0, 1) = gx[i] * gy[j] * lz[k];
            H(1, 0) = H(0, 1);
            H(0, 2) = gx[i] * ly[j] * gz[k];
            H(2, 0) = H(0, 2);
            H(1, 2) = lx[i] * gy[j] * gz[k];
            H(2, 1) = H(1, 2);
            out[n] = H;
        }
    }
};

// ---------------------------------------------------------------------------
// LagrangeTriFast<Order>
// ---------------------------------------------------------------------------
template<int Order>
struct LagrangeTriFast;

template<>
struct LagrangeTriFast<1> {
    static constexpr int n_dofs = 3;

    static constexpr void evaluate(const math::Vector<Real, 3>& xi, std::array<Real, n_dofs>& out) {
        out[0] = Real(1) - xi[0] - xi[1];
        out[1] = xi[0];
        out[2] = xi[1];
    }

    static constexpr void evaluate_gradients(const math::Vector<Real, 3>& /*xi*/,
                                             std::array<Gradient, n_dofs>& out) {
        out[0] = Gradient{Real(-1), Real(-1), Real(0)};
        out[1] = Gradient{Real( 1), Real( 0), Real(0)};
        out[2] = Gradient{Real( 0), Real( 1), Real(0)};
    }

    static constexpr void evaluate_hessians(const math::Vector<Real, 3>& /*xi*/,
                                            std::array<Hessian, n_dofs>& out) {
        out[0] = Hessian{};
        out[1] = Hessian{};
        out[2] = Hessian{};
    }
};

template<>
struct LagrangeTriFast<2> {
    static constexpr int n_dofs = 6;

    static constexpr void evaluate(const math::Vector<Real, 3>& xi, std::array<Real, n_dofs>& out) {
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l0 = Real(1) - l1 - l2;

        out[0] = l0 * (Real(2) * l0 - Real(1));
        out[1] = l1 * (Real(2) * l1 - Real(1));
        out[2] = l2 * (Real(2) * l2 - Real(1));
        out[3] = Real(4) * l0 * l1;
        out[4] = Real(4) * l1 * l2;
        out[5] = Real(4) * l0 * l2;
    }

    static constexpr void evaluate_gradients(const math::Vector<Real, 3>& xi,
                                             std::array<Gradient, n_dofs>& out) {
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l0 = Real(1) - l1 - l2;
        constexpr Gradient g0{Real(-1), Real(-1), Real(0)};
        constexpr Gradient g1{Real( 1), Real( 0), Real(0)};
        constexpr Gradient g2{Real( 0), Real( 1), Real(0)};

        out[0] = detail::scaled_gradient(g0, Real(4) * l0 - Real(1));
        out[1] = detail::scaled_gradient(g1, Real(4) * l1 - Real(1));
        out[2] = detail::scaled_gradient(g2, Real(4) * l2 - Real(1));
        out[3] = detail::p2_edge_gradient(l0, g0, l1, g1);
        out[4] = detail::p2_edge_gradient(l1, g1, l2, g2);
        out[5] = detail::p2_edge_gradient(l0, g0, l2, g2);
    }

    static constexpr void evaluate_hessians(const math::Vector<Real, 3>& /*xi*/,
                                            std::array<Hessian, n_dofs>& out) {
        constexpr Gradient g0{Real(-1), Real(-1), Real(0)};
        constexpr Gradient g1{Real( 1), Real( 0), Real(0)};
        constexpr Gradient g2{Real( 0), Real( 1), Real(0)};

        out[0] = detail::p2_vertex_hessian(g0);
        out[1] = detail::p2_vertex_hessian(g1);
        out[2] = detail::p2_vertex_hessian(g2);
        out[3] = detail::p2_edge_hessian(g0, g1);
        out[4] = detail::p2_edge_hessian(g1, g2);
        out[5] = detail::p2_edge_hessian(g0, g2);
    }
};

template<>
struct LagrangeTriFast<3> {
    static constexpr int n_dofs = 10;

    static constexpr std::array<std::array<std::size_t, 3>, n_dofs> exponents =
        detail::make_triangle_simplex_exponents<3>();

    static void evaluate(const math::Vector<Real, 3>& xi, std::array<Real, n_dofs>& out) {
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l0 = Real(1) - l1 - l2;
        std::array<Real, 4> phi0{};
        std::array<Real, 4> phi1{};
        std::array<Real, 4> phi2{};
        detail::fill_simplex_factor_values<3>(l0, phi0);
        detail::fill_simplex_factor_values<3>(l1, phi1);
        detail::fill_simplex_factor_values<3>(l2, phi2);

        for (std::size_t n = 0; n < exponents.size(); ++n) {
            const auto& e = exponents[n];
            out[n] = phi0[e[0]] * phi1[e[1]] * phi2[e[2]];
        }
    }

    static void evaluate_gradients(const math::Vector<Real, 3>& xi,
                                   std::array<Gradient, n_dofs>& out) {
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l0 = Real(1) - l1 - l2;
        std::array<Real, 4> phi0{};
        std::array<Real, 4> phi1{};
        std::array<Real, 4> phi2{};
        std::array<Real, 4> dphi0{};
        std::array<Real, 4> dphi1{};
        std::array<Real, 4> dphi2{};
        detail::fill_simplex_factor_values_first<3>(l0, phi0, dphi0);
        detail::fill_simplex_factor_values_first<3>(l1, phi1, dphi1);
        detail::fill_simplex_factor_values_first<3>(l2, phi2, dphi2);

        for (std::size_t n = 0; n < exponents.size(); ++n) {
            const auto& e = exponents[n];
            const Real v0 = phi0[e[0]];
            const Real v1 = phi1[e[1]];
            const Real v2 = phi2[e[2]];
            const Real dl0 = dphi0[e[0]] * v1 * v2;
            const Real dl1 = v0 * dphi1[e[1]] * v2;
            const Real dl2 = v0 * v1 * dphi2[e[2]];
            out[n] = Gradient{dl1 - dl0, dl2 - dl0, Real(0)};
        }
    }

    static void evaluate_hessians(const math::Vector<Real, 3>& xi,
                                  std::array<Hessian, n_dofs>& out) {
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l0 = Real(1) - l1 - l2;
        std::array<Real, 4> phi0{};
        std::array<Real, 4> phi1{};
        std::array<Real, 4> phi2{};
        std::array<Real, 4> dphi0{};
        std::array<Real, 4> dphi1{};
        std::array<Real, 4> dphi2{};
        std::array<Real, 4> d2phi0{};
        std::array<Real, 4> d2phi1{};
        std::array<Real, 4> d2phi2{};
        detail::fill_simplex_factor_values_first_second<3>(l0, phi0, dphi0, d2phi0);
        detail::fill_simplex_factor_values_first_second<3>(l1, phi1, dphi1, d2phi1);
        detail::fill_simplex_factor_values_first_second<3>(l2, phi2, dphi2, d2phi2);

        for (std::size_t n = 0; n < exponents.size(); ++n) {
            const auto& e = exponents[n];
            const Real v0 = phi0[e[0]];
            const Real v1 = phi1[e[1]];
            const Real v2 = phi2[e[2]];
            const Real D0 = dphi0[e[0]];
            const Real D1 = dphi1[e[1]];
            const Real D2 = dphi2[e[2]];
            const Real H00 = d2phi0[e[0]] * v1 * v2;
            const Real H11 = v0 * d2phi1[e[1]] * v2;
            const Real H22 = v0 * v1 * d2phi2[e[2]];
            const Real H01 = D0 * D1 * v2;
            const Real H02 = D0 * v1 * D2;
            const Real H12 = v0 * D1 * D2;

            Hessian H{};
            H(0, 0) = H00 - Real(2) * H01 + H11;
            H(1, 1) = H00 - Real(2) * H02 + H22;
            H(0, 1) = H00 - H01 - H02 + H12;
            H(1, 0) = H(0, 1);
            out[n] = H;
        }
    }
};

// ---------------------------------------------------------------------------
// LagrangeTetFast<Order>
// ---------------------------------------------------------------------------
template<int Order>
struct LagrangeTetFast;

template<>
struct LagrangeTetFast<1> {
    static constexpr int n_dofs = 4;

    static constexpr void evaluate(const math::Vector<Real, 3>& xi, std::array<Real, n_dofs>& out) {
        out[0] = Real(1) - xi[0] - xi[1] - xi[2];
        out[1] = xi[0];
        out[2] = xi[1];
        out[3] = xi[2];
    }

    static constexpr void evaluate_gradients(const math::Vector<Real, 3>& /*xi*/,
                                             std::array<Gradient, n_dofs>& out) {
        out[0] = Gradient{Real(-1), Real(-1), Real(-1)};
        out[1] = Gradient{Real( 1), Real( 0), Real( 0)};
        out[2] = Gradient{Real( 0), Real( 1), Real( 0)};
        out[3] = Gradient{Real( 0), Real( 0), Real( 1)};
    }

    static constexpr void evaluate_hessians(const math::Vector<Real, 3>& /*xi*/,
                                            std::array<Hessian, n_dofs>& out) {
        out[0] = Hessian{};
        out[1] = Hessian{};
        out[2] = Hessian{};
        out[3] = Hessian{};
    }
};

template<>
struct LagrangeTetFast<2> {
    static constexpr int n_dofs = 10;

    static constexpr void evaluate(const math::Vector<Real, 3>& xi, std::array<Real, n_dofs>& out) {
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l3 = xi[2];
        const Real l0 = Real(1) - l1 - l2 - l3;

        out[0] = l0 * (Real(2) * l0 - Real(1));
        out[1] = l1 * (Real(2) * l1 - Real(1));
        out[2] = l2 * (Real(2) * l2 - Real(1));
        out[3] = l3 * (Real(2) * l3 - Real(1));
        out[4] = Real(4) * l0 * l1;
        out[5] = Real(4) * l1 * l2;
        out[6] = Real(4) * l0 * l2;
        out[7] = Real(4) * l0 * l3;
        out[8] = Real(4) * l1 * l3;
        out[9] = Real(4) * l2 * l3;
    }

    static constexpr void evaluate_gradients(const math::Vector<Real, 3>& xi,
                                             std::array<Gradient, n_dofs>& out) {
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l3 = xi[2];
        const Real l0 = Real(1) - l1 - l2 - l3;
        constexpr Gradient g0{Real(-1), Real(-1), Real(-1)};
        constexpr Gradient g1{Real( 1), Real( 0), Real( 0)};
        constexpr Gradient g2{Real( 0), Real( 1), Real( 0)};
        constexpr Gradient g3{Real( 0), Real( 0), Real( 1)};

        out[0] = detail::scaled_gradient(g0, Real(4) * l0 - Real(1));
        out[1] = detail::scaled_gradient(g1, Real(4) * l1 - Real(1));
        out[2] = detail::scaled_gradient(g2, Real(4) * l2 - Real(1));
        out[3] = detail::scaled_gradient(g3, Real(4) * l3 - Real(1));
        out[4] = detail::p2_edge_gradient(l0, g0, l1, g1);
        out[5] = detail::p2_edge_gradient(l1, g1, l2, g2);
        out[6] = detail::p2_edge_gradient(l0, g0, l2, g2);
        out[7] = detail::p2_edge_gradient(l0, g0, l3, g3);
        out[8] = detail::p2_edge_gradient(l1, g1, l3, g3);
        out[9] = detail::p2_edge_gradient(l2, g2, l3, g3);
    }

    static constexpr void evaluate_hessians(const math::Vector<Real, 3>& /*xi*/,
                                            std::array<Hessian, n_dofs>& out) {
        constexpr Gradient g0{Real(-1), Real(-1), Real(-1)};
        constexpr Gradient g1{Real( 1), Real( 0), Real( 0)};
        constexpr Gradient g2{Real( 0), Real( 1), Real( 0)};
        constexpr Gradient g3{Real( 0), Real( 0), Real( 1)};

        out[0] = detail::p2_vertex_hessian(g0);
        out[1] = detail::p2_vertex_hessian(g1);
        out[2] = detail::p2_vertex_hessian(g2);
        out[3] = detail::p2_vertex_hessian(g3);
        out[4] = detail::p2_edge_hessian(g0, g1);
        out[5] = detail::p2_edge_hessian(g1, g2);
        out[6] = detail::p2_edge_hessian(g0, g2);
        out[7] = detail::p2_edge_hessian(g0, g3);
        out[8] = detail::p2_edge_hessian(g1, g3);
        out[9] = detail::p2_edge_hessian(g2, g3);
    }
};

template<>
struct LagrangeTetFast<3> {
    static constexpr int n_dofs = 20;

    static constexpr std::array<std::array<std::size_t, 4>, n_dofs> exponents =
        detail::make_tetrahedron_simplex_exponents<3>();

    static void evaluate(const math::Vector<Real, 3>& xi, std::array<Real, n_dofs>& out) {
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l3 = xi[2];
        const Real l0 = Real(1) - l1 - l2 - l3;
        std::array<Real, 4> phi0{};
        std::array<Real, 4> phi1{};
        std::array<Real, 4> phi2{};
        std::array<Real, 4> phi3{};
        detail::fill_simplex_factor_values<3>(l0, phi0);
        detail::fill_simplex_factor_values<3>(l1, phi1);
        detail::fill_simplex_factor_values<3>(l2, phi2);
        detail::fill_simplex_factor_values<3>(l3, phi3);

        for (std::size_t n = 0; n < exponents.size(); ++n) {
            const auto& e = exponents[n];
            out[n] = phi0[e[0]] * phi1[e[1]] * phi2[e[2]] * phi3[e[3]];
        }
    }

    static void evaluate_gradients(const math::Vector<Real, 3>& xi,
                                   std::array<Gradient, n_dofs>& out) {
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l3 = xi[2];
        const Real l0 = Real(1) - l1 - l2 - l3;
        std::array<Real, 4> phi0{};
        std::array<Real, 4> phi1{};
        std::array<Real, 4> phi2{};
        std::array<Real, 4> phi3{};
        std::array<Real, 4> dphi0{};
        std::array<Real, 4> dphi1{};
        std::array<Real, 4> dphi2{};
        std::array<Real, 4> dphi3{};
        detail::fill_simplex_factor_values_first<3>(l0, phi0, dphi0);
        detail::fill_simplex_factor_values_first<3>(l1, phi1, dphi1);
        detail::fill_simplex_factor_values_first<3>(l2, phi2, dphi2);
        detail::fill_simplex_factor_values_first<3>(l3, phi3, dphi3);

        for (std::size_t n = 0; n < exponents.size(); ++n) {
            const auto& e = exponents[n];
            const Real v0 = phi0[e[0]];
            const Real v1 = phi1[e[1]];
            const Real v2 = phi2[e[2]];
            const Real v3 = phi3[e[3]];
            const Real dl0 = dphi0[e[0]] * v1 * v2 * v3;
            const Real dl1 = v0 * dphi1[e[1]] * v2 * v3;
            const Real dl2 = v0 * v1 * dphi2[e[2]] * v3;
            const Real dl3 = v0 * v1 * v2 * dphi3[e[3]];
            out[n] = Gradient{dl1 - dl0, dl2 - dl0, dl3 - dl0};
        }
    }

    static void evaluate_hessians(const math::Vector<Real, 3>& xi,
                                  std::array<Hessian, n_dofs>& out) {
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l3 = xi[2];
        const Real l0 = Real(1) - l1 - l2 - l3;
        std::array<Real, 4> phi0{};
        std::array<Real, 4> phi1{};
        std::array<Real, 4> phi2{};
        std::array<Real, 4> phi3{};
        std::array<Real, 4> dphi0{};
        std::array<Real, 4> dphi1{};
        std::array<Real, 4> dphi2{};
        std::array<Real, 4> dphi3{};
        std::array<Real, 4> d2phi0{};
        std::array<Real, 4> d2phi1{};
        std::array<Real, 4> d2phi2{};
        std::array<Real, 4> d2phi3{};
        detail::fill_simplex_factor_values_first_second<3>(l0, phi0, dphi0, d2phi0);
        detail::fill_simplex_factor_values_first_second<3>(l1, phi1, dphi1, d2phi1);
        detail::fill_simplex_factor_values_first_second<3>(l2, phi2, dphi2, d2phi2);
        detail::fill_simplex_factor_values_first_second<3>(l3, phi3, dphi3, d2phi3);

        for (std::size_t n = 0; n < exponents.size(); ++n) {
            const auto& e = exponents[n];
            const Real v0 = phi0[e[0]];
            const Real v1 = phi1[e[1]];
            const Real v2 = phi2[e[2]];
            const Real v3 = phi3[e[3]];
            const Real D0 = dphi0[e[0]];
            const Real D1 = dphi1[e[1]];
            const Real D2 = dphi2[e[2]];
            const Real D3 = dphi3[e[3]];

            const Real H00 = d2phi0[e[0]] * v1 * v2 * v3;
            const Real H11 = v0 * d2phi1[e[1]] * v2 * v3;
            const Real H22 = v0 * v1 * d2phi2[e[2]] * v3;
            const Real H33 = v0 * v1 * v2 * d2phi3[e[3]];
            const Real H01 = D0 * D1 * v2 * v3;
            const Real H02 = D0 * v1 * D2 * v3;
            const Real H03 = D0 * v1 * v2 * D3;
            const Real H12 = v0 * D1 * D2 * v3;
            const Real H13 = v0 * D1 * v2 * D3;
            const Real H23 = v0 * v1 * D2 * D3;

            Hessian H{};
            H(0, 0) = H00 - Real(2) * H01 + H11;
            H(1, 1) = H00 - Real(2) * H02 + H22;
            H(2, 2) = H00 - Real(2) * H03 + H33;
            H(0, 1) = H00 - H01 - H02 + H12;
            H(1, 0) = H(0, 1);
            H(0, 2) = H00 - H01 - H03 + H13;
            H(2, 0) = H(0, 2);
            H(1, 2) = H00 - H02 - H03 + H23;
            H(2, 1) = H(1, 2);
            out[n] = H;
        }
    }
};

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_LAGRANGEBASISFAST_H
