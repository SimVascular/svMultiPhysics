#include "LagrangeBasisSimplex.h"

#include <array>
#include <vector>

namespace svmp {
namespace FE {
namespace basis {
namespace detail {

// Falling-factorial (equispaced barycentric) Lagrange factors for simplex nodes.
//
// For a fixed polynomial order p and barycentric coordinate lambda in [0, 1],
// define
//   phi_a(lambda) = product_{m=0}^{a-1} (p * lambda - m) / (a - m), a = 0..p
// Then for a multi-index (i0, i1, ..., id) with sum i_k = p, the simplex
// Lagrange basis function is product_k phi_{i_k}(lambda_k), nodal on the
// barycentric lattice.
//
// Output buffers must each be sized to at least p+1 entries; the function
// writes every output slot (no pre-zero required by the caller).
template <bool NeedFirst, bool NeedSecond>
void simplex_lagrange_factor_sequence_impl(int p,
                                           Real lambda,
                                           Real* phi,
                                           Real* dphi,
                                           Real* d2phi) {
    static_assert(!NeedSecond || NeedFirst,
                  "second derivative factors require first-derivative recurrence state");

    phi[0] = Real(1);
    if constexpr (NeedFirst) {
        dphi[0] = Real(0);
    }
    if constexpr (NeedSecond) {
        d2phi[0] = Real(0);
    }
    if (p == 0) {
        return;
    }

    const Real t = static_cast<Real>(p) * lambda;
    const Real dt_dlambda = static_cast<Real>(p);

    Real dphi_dt_prev = Real(0);
    Real d2phi_dt2_prev = Real(0);

    for (int a = 1; a <= p; ++a) {
        const std::size_t au = static_cast<std::size_t>(a);
        const Real inv_a = Real(1) / static_cast<Real>(a);
        const Real s = (t - static_cast<Real>(a - 1)) * inv_a;

        phi[au] = s * phi[au - 1];

        if constexpr (NeedFirst) {
            const Real dphi_dt_old = dphi_dt_prev;
            const Real dphi_dt = inv_a * phi[au - 1] + s * dphi_dt_old;
            dphi[au] = dt_dlambda * dphi_dt;

            if constexpr (NeedSecond) {
                const Real d2phi_dt2 = Real(2) * inv_a * dphi_dt_old + s * d2phi_dt2_prev;
                d2phi[au] = dt_dlambda * dt_dlambda * d2phi_dt2;
                d2phi_dt2_prev = d2phi_dt2;
            }

            dphi_dt_prev = dphi_dt;
        }
    }
}

void simplex_lagrange_factor_sequence(int p,
                                      Real lambda,
                                      Real* phi,
                                      Real* dphi,
                                      Real* d2phi) {
    if (d2phi != nullptr) {
        simplex_lagrange_factor_sequence_impl<true, true>(p, lambda, phi, dphi, d2phi);
    } else if (dphi != nullptr) {
        simplex_lagrange_factor_sequence_impl<true, false>(p, lambda, phi, dphi, nullptr);
    } else {
        simplex_lagrange_factor_sequence_impl<false, false>(p, lambda, phi, nullptr, nullptr);
    }
}

constexpr int kFixedSimplexAxisOrder = 12;
constexpr std::size_t kFixedSimplexAxisSize =
    static_cast<std::size_t>(kFixedSimplexAxisOrder + 1);
constexpr std::size_t kFixedSimplexBatchEntries = 512;

template <int Order>
inline void simplex_lagrange_factor_values_product(Real lambda,
                                                   Real* SVMP_RESTRICT values) {
    static_assert(Order >= 0, "simplex order must be non-negative");
    values[0] = Real(1);
    const Real t = static_cast<Real>(Order) * lambda;
    for (int a = 1; a <= Order; ++a) {
        const Real inv_a = Real(1) / static_cast<Real>(a);
        values[a] = values[a - 1] * (t - static_cast<Real>(a - 1)) * inv_a;
    }
}

template <int Order>
void evaluate_triangle_simplex_values_q4(
    const std::vector<std::array<int, 4>>& simplex_exponents,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    static_assert(Order >= 4 && Order <= 8, "specialized simplex path covers orders 4..8");

    Real phi0[4][Order + 1];
    Real phi1[4][Order + 1];
    Real phi2[4][Order + 1];

    for (std::size_t q = 0; q < 4u; ++q) {
        const auto& xi = points[q];
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l0 = Real(1) - l1 - l2;
        simplex_lagrange_factor_values_product<Order>(l0, phi0[q]);
        simplex_lagrange_factor_values_product<Order>(l1, phi1[q]);
        simplex_lagrange_factor_values_product<Order>(l2, phi2[q]);
    }

    const std::size_t num_nodes = simplex_exponents.size();
    for (std::size_t node = 0; node < num_nodes; ++node) {
        const auto& e = simplex_exponents[node];
        const std::size_t i0 = static_cast<std::size_t>(e[0]);
        const std::size_t i1 = static_cast<std::size_t>(e[1]);
        const std::size_t i2 = static_cast<std::size_t>(e[2]);
        Real* SVMP_RESTRICT row = values_out + node * output_stride;
        row[0] = phi0[0][i0] * phi1[0][i1] * phi2[0][i2];
        row[1] = phi0[1][i0] * phi1[1][i1] * phi2[1][i2];
        row[2] = phi0[2][i0] * phi1[2][i1] * phi2[2][i2];
        row[3] = phi0[3][i0] * phi1[3][i1] * phi2[3][i2];
    }
}

bool try_evaluate_triangle_simplex_values_q4(
    const std::vector<std::array<int, 4>>& simplex_exponents,
    int order,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    switch (order) {
    case 4:
        evaluate_triangle_simplex_values_q4<4>(
            simplex_exponents, points, output_stride, values_out);
        return true;
    case 5:
        evaluate_triangle_simplex_values_q4<5>(
            simplex_exponents, points, output_stride, values_out);
        return true;
    case 6:
        evaluate_triangle_simplex_values_q4<6>(
            simplex_exponents, points, output_stride, values_out);
        return true;
    case 7:
        evaluate_triangle_simplex_values_q4<7>(
            simplex_exponents, points, output_stride, values_out);
        return true;
    case 8:
        evaluate_triangle_simplex_values_q4<8>(
            simplex_exponents, points, output_stride, values_out);
        return true;
    default:
        return false;
    }
}

template <int Order>
void evaluate_tetrahedron_simplex_values_q4(
    const std::vector<std::array<int, 4>>& simplex_exponents,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    static_assert(Order >= 4 && Order <= 8, "specialized simplex path covers orders 4..8");

    Real phi0[4][Order + 1];
    Real phi1[4][Order + 1];
    Real phi2[4][Order + 1];
    Real phi3[4][Order + 1];

    for (std::size_t q = 0; q < 4u; ++q) {
        const auto& xi = points[q];
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l3 = xi[2];
        const Real l0 = Real(1) - l1 - l2 - l3;
        simplex_lagrange_factor_values_product<Order>(l0, phi0[q]);
        simplex_lagrange_factor_values_product<Order>(l1, phi1[q]);
        simplex_lagrange_factor_values_product<Order>(l2, phi2[q]);
        simplex_lagrange_factor_values_product<Order>(l3, phi3[q]);
    }

    const std::size_t num_nodes = simplex_exponents.size();
    for (std::size_t node = 0; node < num_nodes; ++node) {
        const auto& e = simplex_exponents[node];
        const std::size_t i0 = static_cast<std::size_t>(e[0]);
        const std::size_t i1 = static_cast<std::size_t>(e[1]);
        const std::size_t i2 = static_cast<std::size_t>(e[2]);
        const std::size_t i3 = static_cast<std::size_t>(e[3]);
        Real* SVMP_RESTRICT row = values_out + node * output_stride;
        row[0] = phi0[0][i0] * phi1[0][i1] * phi2[0][i2] * phi3[0][i3];
        row[1] = phi0[1][i0] * phi1[1][i1] * phi2[1][i2] * phi3[1][i3];
        row[2] = phi0[2][i0] * phi1[2][i1] * phi2[2][i2] * phi3[2][i3];
        row[3] = phi0[3][i0] * phi1[3][i1] * phi2[3][i2] * phi3[3][i3];
    }
}

bool try_evaluate_tetrahedron_simplex_values_q4(
    const std::vector<std::array<int, 4>>& simplex_exponents,
    int order,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    switch (order) {
    case 4:
        evaluate_tetrahedron_simplex_values_q4<4>(
            simplex_exponents, points, output_stride, values_out);
        return true;
    case 5:
        evaluate_tetrahedron_simplex_values_q4<5>(
            simplex_exponents, points, output_stride, values_out);
        return true;
    case 6:
        evaluate_tetrahedron_simplex_values_q4<6>(
            simplex_exponents, points, output_stride, values_out);
        return true;
    case 7:
        evaluate_tetrahedron_simplex_values_q4<7>(
            simplex_exponents, points, output_stride, values_out);
        return true;
    case 8:
        evaluate_tetrahedron_simplex_values_q4<8>(
            simplex_exponents, points, output_stride, values_out);
        return true;
    default:
        return false;
    }
}

template <int Order>
void evaluate_tetrahedron_simplex_gradients_q4(
    const std::vector<std::array<int, 4>>& simplex_exponents,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT gradients_out) {
    static_assert(Order >= 3 && Order <= 8,
                  "specialized tetrahedron gradient path covers orders 3..8");

    Real phi0[4][Order + 1];
    Real phi1[4][Order + 1];
    Real phi2[4][Order + 1];
    Real phi3[4][Order + 1];
    Real dphi0[4][Order + 1];
    Real dphi1[4][Order + 1];
    Real dphi2[4][Order + 1];
    Real dphi3[4][Order + 1];

    for (std::size_t q = 0; q < 4u; ++q) {
        const auto& xi = points[q];
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l3 = xi[2];
        const Real l0 = Real(1) - l1 - l2 - l3;
        simplex_lagrange_factor_sequence_impl<true, false>(
            Order, l0, phi0[q], dphi0[q], nullptr);
        simplex_lagrange_factor_sequence_impl<true, false>(
            Order, l1, phi1[q], dphi1[q], nullptr);
        simplex_lagrange_factor_sequence_impl<true, false>(
            Order, l2, phi2[q], dphi2[q], nullptr);
        simplex_lagrange_factor_sequence_impl<true, false>(
            Order, l3, phi3[q], dphi3[q], nullptr);
    }

    const std::size_t num_nodes = simplex_exponents.size();
    for (std::size_t node = 0; node < num_nodes; ++node) {
        const auto& e = simplex_exponents[node];
        const std::size_t i0 = static_cast<std::size_t>(e[0]);
        const std::size_t i1 = static_cast<std::size_t>(e[1]);
        const std::size_t i2 = static_cast<std::size_t>(e[2]);
        const std::size_t i3 = static_cast<std::size_t>(e[3]);
        Real gx[4];
        Real gy[4];
        Real gz[4];

        for (std::size_t q = 0; q < 4u; ++q) {
            const Real v0 = phi0[q][i0];
            const Real v1 = phi1[q][i1];
            const Real v2 = phi2[q][i2];
            const Real v3 = phi3[q][i3];
            const Real D0 = dphi0[q][i0];
            const Real D1 = dphi1[q][i1];
            const Real D2 = dphi2[q][i2];
            const Real D3 = dphi3[q][i3];
            const Real v23 = v2 * v3;
            const Real v01 = v0 * v1;
            const Real dl0 = D0 * v1 * v23;
            gx[q] = v0 * D1 * v23 - dl0;
            gy[q] = v01 * D2 * v3 - dl0;
            gz[q] = v01 * v2 * D3 - dl0;
        }

        Real* SVMP_RESTRICT g = gradients_out + node * 3u * output_stride;
        g[0u] = gx[0];
        g[1u] = gx[1];
        g[2u] = gx[2];
        g[3u] = gx[3];
        g[output_stride + 0u] = gy[0];
        g[output_stride + 1u] = gy[1];
        g[output_stride + 2u] = gy[2];
        g[output_stride + 3u] = gy[3];
        g[2u * output_stride + 0u] = gz[0];
        g[2u * output_stride + 1u] = gz[1];
        g[2u * output_stride + 2u] = gz[2];
        g[2u * output_stride + 3u] = gz[3];
    }
}

template <int Order>
void evaluate_triangle_simplex_gradients_q4(
    const std::vector<std::array<int, 4>>& simplex_exponents,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT gradients_out) {
    static_assert((Order == 2) || (Order >= 4 && Order <= 8),
                  "specialized simplex path covers order 2 and orders 4..8");

    Real phi0[4][Order + 1];
    Real phi1[4][Order + 1];
    Real phi2[4][Order + 1];
    Real dphi0[4][Order + 1];
    Real dphi1[4][Order + 1];
    Real dphi2[4][Order + 1];

    for (std::size_t q = 0; q < 4u; ++q) {
        const auto& xi = points[q];
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l0 = Real(1) - l1 - l2;
        simplex_lagrange_factor_sequence_impl<true, false>(
            Order, l0, phi0[q], dphi0[q], nullptr);
        simplex_lagrange_factor_sequence_impl<true, false>(
            Order, l1, phi1[q], dphi1[q], nullptr);
        simplex_lagrange_factor_sequence_impl<true, false>(
            Order, l2, phi2[q], dphi2[q], nullptr);
    }

    const std::size_t num_nodes = simplex_exponents.size();
    for (std::size_t node = 0; node < num_nodes; ++node) {
        const auto& e = simplex_exponents[node];
        const std::size_t i0 = static_cast<std::size_t>(e[0]);
        const std::size_t i1 = static_cast<std::size_t>(e[1]);
        const std::size_t i2 = static_cast<std::size_t>(e[2]);
        Real* SVMP_RESTRICT g = gradients_out + node * 3u * output_stride;

        for (std::size_t q = 0; q < 4u; ++q) {
            const Real v0 = phi0[q][i0];
            const Real v1 = phi1[q][i1];
            const Real v2 = phi2[q][i2];
            const Real D0 = dphi0[q][i0];
            const Real D1 = dphi1[q][i1];
            const Real D2 = dphi2[q][i2];
            const Real dl0 = D0 * v1 * v2;
            g[0u * output_stride + q] = v0 * D1 * v2 - dl0;
            g[1u * output_stride + q] = v0 * v1 * D2 - dl0;
            g[2u * output_stride + q] = Real(0);
        }
    }
}

bool try_evaluate_triangle_simplex_gradients_q4(
    const std::vector<std::array<int, 4>>& simplex_exponents,
    int order,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT gradients_out) {
    switch (order) {
    case 2:
        evaluate_triangle_simplex_gradients_q4<2>(
            simplex_exponents, points, output_stride, gradients_out);
        return true;
    case 4:
        evaluate_triangle_simplex_gradients_q4<4>(
            simplex_exponents, points, output_stride, gradients_out);
        return true;
    case 5:
        evaluate_triangle_simplex_gradients_q4<5>(
            simplex_exponents, points, output_stride, gradients_out);
        return true;
    case 6:
        evaluate_triangle_simplex_gradients_q4<6>(
            simplex_exponents, points, output_stride, gradients_out);
        return true;
    case 7:
        evaluate_triangle_simplex_gradients_q4<7>(
            simplex_exponents, points, output_stride, gradients_out);
        return true;
    case 8:
        evaluate_triangle_simplex_gradients_q4<8>(
            simplex_exponents, points, output_stride, gradients_out);
        return true;
    default:
        return false;
    }
}

template <int Order>
void evaluate_triangle_simplex_hessian_outputs_q4(
    const std::vector<std::array<int, 4>>& simplex_exponents,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    static_assert(Order >= 2 && Order <= 8, "specialized simplex path covers orders 2..8");

    Real phi0[4][Order + 1];
    Real phi1[4][Order + 1];
    Real phi2[4][Order + 1];
    Real dphi0[4][Order + 1];
    Real dphi1[4][Order + 1];
    Real dphi2[4][Order + 1];
    Real d2phi0[4][Order + 1];
    Real d2phi1[4][Order + 1];
    Real d2phi2[4][Order + 1];

    for (std::size_t q = 0; q < 4u; ++q) {
        const auto& xi = points[q];
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l0 = Real(1) - l1 - l2;
        simplex_lagrange_factor_sequence_impl<true, true>(
            Order, l0, phi0[q], dphi0[q], d2phi0[q]);
        simplex_lagrange_factor_sequence_impl<true, true>(
            Order, l1, phi1[q], dphi1[q], d2phi1[q]);
        simplex_lagrange_factor_sequence_impl<true, true>(
            Order, l2, phi2[q], dphi2[q], d2phi2[q]);
    }

    const std::size_t num_nodes = simplex_exponents.size();
    for (std::size_t node = 0; node < num_nodes; ++node) {
        const auto& e = simplex_exponents[node];
        const std::size_t i0 = static_cast<std::size_t>(e[0]);
        const std::size_t i1 = static_cast<std::size_t>(e[1]);
        const std::size_t i2 = static_cast<std::size_t>(e[2]);
        Real* SVMP_RESTRICT value_row = values_out ? values_out + node * output_stride : nullptr;
        Real* SVMP_RESTRICT g = gradients_out ? gradients_out + node * 3u * output_stride : nullptr;
        Real* SVMP_RESTRICT H = hessians_out + node * 9u * output_stride;
        H[2u * output_stride + 0u] = Real(0);
        H[2u * output_stride + 1u] = Real(0);
        H[2u * output_stride + 2u] = Real(0);
        H[2u * output_stride + 3u] = Real(0);
        H[5u * output_stride + 0u] = Real(0);
        H[5u * output_stride + 1u] = Real(0);
        H[5u * output_stride + 2u] = Real(0);
        H[5u * output_stride + 3u] = Real(0);
        H[6u * output_stride + 0u] = Real(0);
        H[6u * output_stride + 1u] = Real(0);
        H[6u * output_stride + 2u] = Real(0);
        H[6u * output_stride + 3u] = Real(0);
        H[7u * output_stride + 0u] = Real(0);
        H[7u * output_stride + 1u] = Real(0);
        H[7u * output_stride + 2u] = Real(0);
        H[7u * output_stride + 3u] = Real(0);
        H[8u * output_stride + 0u] = Real(0);
        H[8u * output_stride + 1u] = Real(0);
        H[8u * output_stride + 2u] = Real(0);
        H[8u * output_stride + 3u] = Real(0);

        for (std::size_t q = 0; q < 4u; ++q) {
            const Real v0 = phi0[q][i0];
            const Real v1 = phi1[q][i1];
            const Real v2 = phi2[q][i2];
            if (value_row != nullptr) {
                value_row[q] = v0 * v1 * v2;
            }

            const Real D0 = dphi0[q][i0];
            const Real D1 = dphi1[q][i1];
            const Real D2 = dphi2[q][i2];
            if (g != nullptr) {
                const Real dl0 = D0 * v1 * v2;
                g[0u * output_stride + q] = v0 * D1 * v2 - dl0;
                g[1u * output_stride + q] = v0 * v1 * D2 - dl0;
                g[2u * output_stride + q] = Real(0);
            }

            const Real DD0 = d2phi0[q][i0];
            const Real DD1 = d2phi1[q][i1];
            const Real DD2 = d2phi2[q][i2];
            const Real H00 = DD0 * v1 * v2;
            const Real H11 = v0 * DD1 * v2;
            const Real H22 = v0 * v1 * DD2;
            const Real H01 = D0 * D1 * v2;
            const Real H02 = D0 * v1 * D2;
            const Real H12 = v0 * D1 * D2;
            const Real h01 = H00 - H01 - H02 + H12;
            H[0u * output_stride + q] = H00 - Real(2) * H01 + H11;
            H[1u * output_stride + q] = h01;
            H[3u * output_stride + q] = h01;
            H[4u * output_stride + q] = H00 - Real(2) * H02 + H22;
        }
    }
}

bool try_evaluate_triangle_simplex_hessian_outputs_q4(
    const std::vector<std::array<int, 4>>& simplex_exponents,
    int order,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    switch (order) {
    case 2:
        evaluate_triangle_simplex_hessian_outputs_q4<2>(
            simplex_exponents, points, output_stride, values_out, gradients_out, hessians_out);
        return true;
    case 3:
        evaluate_triangle_simplex_hessian_outputs_q4<3>(
            simplex_exponents, points, output_stride, values_out, gradients_out, hessians_out);
        return true;
    case 4:
        evaluate_triangle_simplex_hessian_outputs_q4<4>(
            simplex_exponents, points, output_stride, values_out, gradients_out, hessians_out);
        return true;
    case 5:
        evaluate_triangle_simplex_hessian_outputs_q4<5>(
            simplex_exponents, points, output_stride, values_out, gradients_out, hessians_out);
        return true;
    case 6:
        evaluate_triangle_simplex_hessian_outputs_q4<6>(
            simplex_exponents, points, output_stride, values_out, gradients_out, hessians_out);
        return true;
    case 7:
        evaluate_triangle_simplex_hessian_outputs_q4<7>(
            simplex_exponents, points, output_stride, values_out, gradients_out, hessians_out);
        return true;
    case 8:
        evaluate_triangle_simplex_hessian_outputs_q4<8>(
            simplex_exponents, points, output_stride, values_out, gradients_out, hessians_out);
        return true;
    default:
        return false;
    }
}

template <int Order, std::size_t Q>
inline void write_tetrahedron_simplex_hessian_q4(
    const Real (&phi0)[4][Order + 1],
    const Real (&phi1)[4][Order + 1],
    const Real (&phi2)[4][Order + 1],
    const Real (&phi3)[4][Order + 1],
    const Real (&dphi0)[4][Order + 1],
    const Real (&dphi1)[4][Order + 1],
    const Real (&dphi2)[4][Order + 1],
    const Real (&dphi3)[4][Order + 1],
    const Real (&d2phi0)[4][Order + 1],
    const Real (&d2phi1)[4][Order + 1],
    const Real (&d2phi2)[4][Order + 1],
    const Real (&d2phi3)[4][Order + 1],
    std::size_t i0,
    std::size_t i1,
    std::size_t i2,
    std::size_t i3,
    std::size_t output_stride,
    Real* SVMP_RESTRICT H) {
    const Real v0 = phi0[Q][i0];
    const Real v1 = phi1[Q][i1];
    const Real v2 = phi2[Q][i2];
    const Real v3 = phi3[Q][i3];
    const Real D0 = dphi0[Q][i0];
    const Real D1 = dphi1[Q][i1];
    const Real D2 = dphi2[Q][i2];
    const Real D3 = dphi3[Q][i3];
    const Real DD0 = d2phi0[Q][i0];
    const Real DD1 = d2phi1[Q][i1];
    const Real DD2 = d2phi2[Q][i2];
    const Real DD3 = d2phi3[Q][i3];
    const Real H00 = DD0 * v1 * v2 * v3;
    const Real H11 = v0 * DD1 * v2 * v3;
    const Real H22 = v0 * v1 * DD2 * v3;
    const Real H33 = v0 * v1 * v2 * DD3;
    const Real H01 = D0 * D1 * v2 * v3;
    const Real H02 = D0 * v1 * D2 * v3;
    const Real H03 = D0 * v1 * v2 * D3;
    const Real H12 = v0 * D1 * D2 * v3;
    const Real H13 = v0 * D1 * v2 * D3;
    const Real H23 = v0 * v1 * D2 * D3;
    const Real h01 = H00 - H01 - H02 + H12;
    const Real h02 = H00 - H01 - H03 + H13;
    const Real h12 = H00 - H02 - H03 + H23;
    H[0u * output_stride + Q] = H00 - Real(2) * H01 + H11;
    H[1u * output_stride + Q] = h01;
    H[2u * output_stride + Q] = h02;
    H[3u * output_stride + Q] = h01;
    H[4u * output_stride + Q] = H00 - Real(2) * H02 + H22;
    H[5u * output_stride + Q] = h12;
    H[6u * output_stride + Q] = h02;
    H[7u * output_stride + Q] = h12;
    H[8u * output_stride + Q] = H00 - Real(2) * H03 + H33;
}

template <int Order, std::size_t Q>
inline void write_tetrahedron_simplex_hessian_stride4_q(
    const Real (&phi0)[4][Order + 1],
    const Real (&phi1)[4][Order + 1],
    const Real (&phi2)[4][Order + 1],
    const Real (&phi3)[4][Order + 1],
    const Real (&dphi0)[4][Order + 1],
    const Real (&dphi1)[4][Order + 1],
    const Real (&dphi2)[4][Order + 1],
    const Real (&dphi3)[4][Order + 1],
    const Real (&d2phi0)[4][Order + 1],
    const Real (&d2phi1)[4][Order + 1],
    const Real (&d2phi2)[4][Order + 1],
    const Real (&d2phi3)[4][Order + 1],
    std::size_t i0,
    std::size_t i1,
    std::size_t i2,
    std::size_t i3,
    Real* SVMP_RESTRICT H) {
    const Real v0 = phi0[Q][i0];
    const Real v1 = phi1[Q][i1];
    const Real v2 = phi2[Q][i2];
    const Real v3 = phi3[Q][i3];
    const Real D0 = dphi0[Q][i0];
    const Real D1 = dphi1[Q][i1];
    const Real D2 = dphi2[Q][i2];
    const Real D3 = dphi3[Q][i3];
    const Real DD0 = d2phi0[Q][i0];
    const Real DD1 = d2phi1[Q][i1];
    const Real DD2 = d2phi2[Q][i2];
    const Real DD3 = d2phi3[Q][i3];
    const Real v12 = v1 * v2;
    const Real v13 = v1 * v3;
    const Real v23 = v2 * v3;
    const Real v123 = v1 * v23;
    const Real v023 = v0 * v23;
    const Real v013 = v0 * v13;
    const Real v012 = v0 * v12;
    const Real H00 = DD0 * v123;
    const Real H11 = DD1 * v023;
    const Real H22 = DD2 * v013;
    const Real H33 = DD3 * v012;
    const Real H01 = D0 * D1 * v23;
    const Real H02 = D0 * D2 * v13;
    const Real H03 = D0 * D3 * v12;
    const Real H12 = D1 * D2 * v0 * v3;
    const Real H13 = D1 * D3 * v0 * v2;
    const Real H23 = D2 * D3 * v0 * v1;
    const Real h01 = H00 - H01 - H02 + H12;
    const Real h02 = H00 - H01 - H03 + H13;
    const Real h12 = H00 - H02 - H03 + H23;
    H[Q] = H00 - Real(2) * H01 + H11;
    H[4u + Q] = h01;
    H[8u + Q] = h02;
    H[12u + Q] = h01;
    H[16u + Q] = H00 - Real(2) * H02 + H22;
    H[20u + Q] = h12;
    H[24u + Q] = h02;
    H[28u + Q] = h12;
    H[32u + Q] = H00 - Real(2) * H03 + H33;
}

template <int Order, std::size_t Q>
inline void write_tetrahedron_simplex_all_stride4_q(
    const Real (&phi0)[4][Order + 1],
    const Real (&phi1)[4][Order + 1],
    const Real (&phi2)[4][Order + 1],
    const Real (&phi3)[4][Order + 1],
    const Real (&dphi0)[4][Order + 1],
    const Real (&dphi1)[4][Order + 1],
    const Real (&dphi2)[4][Order + 1],
    const Real (&dphi3)[4][Order + 1],
    const Real (&d2phi0)[4][Order + 1],
    const Real (&d2phi1)[4][Order + 1],
    const Real (&d2phi2)[4][Order + 1],
    const Real (&d2phi3)[4][Order + 1],
    std::size_t i0,
    std::size_t i1,
    std::size_t i2,
    std::size_t i3,
    Real* SVMP_RESTRICT value_row,
    Real* SVMP_RESTRICT g,
    Real* SVMP_RESTRICT H) {
    const Real v0 = phi0[Q][i0];
    const Real v1 = phi1[Q][i1];
    const Real v2 = phi2[Q][i2];
    const Real v3 = phi3[Q][i3];
    const Real D0 = dphi0[Q][i0];
    const Real D1 = dphi1[Q][i1];
    const Real D2 = dphi2[Q][i2];
    const Real D3 = dphi3[Q][i3];
    const Real DD0 = d2phi0[Q][i0];
    const Real DD1 = d2phi1[Q][i1];
    const Real DD2 = d2phi2[Q][i2];
    const Real DD3 = d2phi3[Q][i3];
    const Real v12 = v1 * v2;
    const Real v13 = v1 * v3;
    const Real v23 = v2 * v3;
    const Real v123 = v1 * v23;
    const Real v023 = v0 * v23;
    const Real v013 = v0 * v13;
    const Real v012 = v0 * v12;
    const Real dl0 = D0 * v123;
    const Real H00 = DD0 * v123;
    const Real H11 = DD1 * v023;
    const Real H22 = DD2 * v013;
    const Real H33 = DD3 * v012;
    const Real H01 = D0 * D1 * v23;
    const Real H02 = D0 * D2 * v13;
    const Real H03 = D0 * D3 * v12;
    const Real H12 = D1 * D2 * v0 * v3;
    const Real H13 = D1 * D3 * v0 * v2;
    const Real H23 = D2 * D3 * v0 * v1;
    const Real h01 = H00 - H01 - H02 + H12;
    const Real h02 = H00 - H01 - H03 + H13;
    const Real h12 = H00 - H02 - H03 + H23;

    value_row[Q] = v0 * v123;
    g[Q] = D1 * v023 - dl0;
    g[4u + Q] = D2 * v013 - dl0;
    g[8u + Q] = D3 * v012 - dl0;
    H[Q] = H00 - Real(2) * H01 + H11;
    H[4u + Q] = h01;
    H[8u + Q] = h02;
    H[12u + Q] = h01;
    H[16u + Q] = H00 - Real(2) * H02 + H22;
    H[20u + Q] = h12;
    H[24u + Q] = h02;
    H[28u + Q] = h12;
    H[32u + Q] = H00 - Real(2) * H03 + H33;
}

template <int Order>
void evaluate_tetrahedron_simplex_hessian_outputs_q4(
    const std::vector<std::array<int, 4>>& simplex_exponents,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    static_assert(Order >= 2 && Order <= 8, "specialized simplex path covers orders 2..8");

    Real phi0[4][Order + 1];
    Real phi1[4][Order + 1];
    Real phi2[4][Order + 1];
    Real phi3[4][Order + 1];
    Real dphi0[4][Order + 1];
    Real dphi1[4][Order + 1];
    Real dphi2[4][Order + 1];
    Real dphi3[4][Order + 1];
    Real d2phi0[4][Order + 1];
    Real d2phi1[4][Order + 1];
    Real d2phi2[4][Order + 1];
    Real d2phi3[4][Order + 1];

    for (std::size_t q = 0; q < 4u; ++q) {
        const auto& xi = points[q];
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l3 = xi[2];
        const Real l0 = Real(1) - l1 - l2 - l3;
        simplex_lagrange_factor_sequence_impl<true, true>(
            Order, l0, phi0[q], dphi0[q], d2phi0[q]);
        simplex_lagrange_factor_sequence_impl<true, true>(
            Order, l1, phi1[q], dphi1[q], d2phi1[q]);
        simplex_lagrange_factor_sequence_impl<true, true>(
            Order, l2, phi2[q], dphi2[q], d2phi2[q]);
        simplex_lagrange_factor_sequence_impl<true, true>(
            Order, l3, phi3[q], dphi3[q], d2phi3[q]);
    }

    const std::size_t num_nodes = simplex_exponents.size();
    if (values_out == nullptr && gradients_out == nullptr) {
        if (output_stride == 4u) {
            for (std::size_t node = 0; node < num_nodes; ++node) {
                const auto& e = simplex_exponents[node];
                const std::size_t i0 = static_cast<std::size_t>(e[0]);
                const std::size_t i1 = static_cast<std::size_t>(e[1]);
                const std::size_t i2 = static_cast<std::size_t>(e[2]);
                const std::size_t i3 = static_cast<std::size_t>(e[3]);
                Real* SVMP_RESTRICT H = hessians_out + node * 36u;
                write_tetrahedron_simplex_hessian_stride4_q<Order, 0>(
                    phi0, phi1, phi2, phi3, dphi0, dphi1, dphi2, dphi3,
                    d2phi0, d2phi1, d2phi2, d2phi3, i0, i1, i2, i3, H);
                write_tetrahedron_simplex_hessian_stride4_q<Order, 1>(
                    phi0, phi1, phi2, phi3, dphi0, dphi1, dphi2, dphi3,
                    d2phi0, d2phi1, d2phi2, d2phi3, i0, i1, i2, i3, H);
                write_tetrahedron_simplex_hessian_stride4_q<Order, 2>(
                    phi0, phi1, phi2, phi3, dphi0, dphi1, dphi2, dphi3,
                    d2phi0, d2phi1, d2phi2, d2phi3, i0, i1, i2, i3, H);
                write_tetrahedron_simplex_hessian_stride4_q<Order, 3>(
                    phi0, phi1, phi2, phi3, dphi0, dphi1, dphi2, dphi3,
                    d2phi0, d2phi1, d2phi2, d2phi3, i0, i1, i2, i3, H);
            }
        } else {
            for (std::size_t node = 0; node < num_nodes; ++node) {
                const auto& e = simplex_exponents[node];
                const std::size_t i0 = static_cast<std::size_t>(e[0]);
                const std::size_t i1 = static_cast<std::size_t>(e[1]);
                const std::size_t i2 = static_cast<std::size_t>(e[2]);
                const std::size_t i3 = static_cast<std::size_t>(e[3]);
                Real* SVMP_RESTRICT H = hessians_out + node * 9u * output_stride;
                write_tetrahedron_simplex_hessian_q4<Order, 0>(
                    phi0, phi1, phi2, phi3, dphi0, dphi1, dphi2, dphi3,
                    d2phi0, d2phi1, d2phi2, d2phi3, i0, i1, i2, i3, output_stride, H);
                write_tetrahedron_simplex_hessian_q4<Order, 1>(
                    phi0, phi1, phi2, phi3, dphi0, dphi1, dphi2, dphi3,
                    d2phi0, d2phi1, d2phi2, d2phi3, i0, i1, i2, i3, output_stride, H);
                write_tetrahedron_simplex_hessian_q4<Order, 2>(
                    phi0, phi1, phi2, phi3, dphi0, dphi1, dphi2, dphi3,
                    d2phi0, d2phi1, d2phi2, d2phi3, i0, i1, i2, i3, output_stride, H);
                write_tetrahedron_simplex_hessian_q4<Order, 3>(
                    phi0, phi1, phi2, phi3, dphi0, dphi1, dphi2, dphi3,
                    d2phi0, d2phi1, d2phi2, d2phi3, i0, i1, i2, i3, output_stride, H);
            }
        }
        return;
    }

    if (values_out != nullptr && gradients_out != nullptr) {
        if (output_stride == 4u) {
            for (std::size_t node = 0; node < num_nodes; ++node) {
                const auto& e = simplex_exponents[node];
                const std::size_t i0 = static_cast<std::size_t>(e[0]);
                const std::size_t i1 = static_cast<std::size_t>(e[1]);
                const std::size_t i2 = static_cast<std::size_t>(e[2]);
                const std::size_t i3 = static_cast<std::size_t>(e[3]);
                Real* SVMP_RESTRICT value_row = values_out + node * output_stride;
                Real* SVMP_RESTRICT g = gradients_out + node * 3u * output_stride;
                Real* SVMP_RESTRICT H = hessians_out + node * 9u * output_stride;
                write_tetrahedron_simplex_all_stride4_q<Order, 0>(
                    phi0, phi1, phi2, phi3, dphi0, dphi1, dphi2, dphi3,
                    d2phi0, d2phi1, d2phi2, d2phi3, i0, i1, i2, i3, value_row, g, H);
                write_tetrahedron_simplex_all_stride4_q<Order, 1>(
                    phi0, phi1, phi2, phi3, dphi0, dphi1, dphi2, dphi3,
                    d2phi0, d2phi1, d2phi2, d2phi3, i0, i1, i2, i3, value_row, g, H);
                write_tetrahedron_simplex_all_stride4_q<Order, 2>(
                    phi0, phi1, phi2, phi3, dphi0, dphi1, dphi2, dphi3,
                    d2phi0, d2phi1, d2phi2, d2phi3, i0, i1, i2, i3, value_row, g, H);
                write_tetrahedron_simplex_all_stride4_q<Order, 3>(
                    phi0, phi1, phi2, phi3, dphi0, dphi1, dphi2, dphi3,
                    d2phi0, d2phi1, d2phi2, d2phi3, i0, i1, i2, i3, value_row, g, H);
            }
            return;
        }

        for (std::size_t node = 0; node < num_nodes; ++node) {
            const auto& e = simplex_exponents[node];
            const std::size_t i0 = static_cast<std::size_t>(e[0]);
            const std::size_t i1 = static_cast<std::size_t>(e[1]);
            const std::size_t i2 = static_cast<std::size_t>(e[2]);
            const std::size_t i3 = static_cast<std::size_t>(e[3]);
            Real* SVMP_RESTRICT value_row = values_out + node * output_stride;
            Real* SVMP_RESTRICT g = gradients_out + node * 3u * output_stride;
            Real* SVMP_RESTRICT H = hessians_out + node * 9u * output_stride;

            for (std::size_t q = 0; q < 4u; ++q) {
                const Real v0 = phi0[q][i0];
                const Real v1 = phi1[q][i1];
                const Real v2 = phi2[q][i2];
                const Real v3 = phi3[q][i3];
                const Real D0 = dphi0[q][i0];
                const Real D1 = dphi1[q][i1];
                const Real D2 = dphi2[q][i2];
                const Real D3 = dphi3[q][i3];
                const Real DD0 = d2phi0[q][i0];
                const Real DD1 = d2phi1[q][i1];
                const Real DD2 = d2phi2[q][i2];
                const Real DD3 = d2phi3[q][i3];
                const Real v12 = v1 * v2;
                const Real v13 = v1 * v3;
                const Real v23 = v2 * v3;
                const Real v123 = v1 * v23;
                const Real v023 = v0 * v23;
                const Real v013 = v0 * v13;
                const Real v012 = v0 * v12;
                const Real dl0 = D0 * v123;
                const Real H00 = DD0 * v123;
                const Real H11 = DD1 * v023;
                const Real H22 = DD2 * v013;
                const Real H33 = DD3 * v012;
                const Real H01 = D0 * D1 * v23;
                const Real H02 = D0 * D2 * v13;
                const Real H03 = D0 * D3 * v12;
                const Real H12 = D1 * D2 * v0 * v3;
                const Real H13 = D1 * D3 * v0 * v2;
                const Real H23 = D2 * D3 * v0 * v1;
                const Real h01 = H00 - H01 - H02 + H12;
                const Real h02 = H00 - H01 - H03 + H13;
                const Real h12 = H00 - H02 - H03 + H23;

                value_row[q] = v0 * v123;
                g[0u * output_stride + q] = D1 * v023 - dl0;
                g[1u * output_stride + q] = D2 * v013 - dl0;
                g[2u * output_stride + q] = D3 * v012 - dl0;
                H[0u * output_stride + q] = H00 - Real(2) * H01 + H11;
                H[1u * output_stride + q] = h01;
                H[2u * output_stride + q] = h02;
                H[3u * output_stride + q] = h01;
                H[4u * output_stride + q] = H00 - Real(2) * H02 + H22;
                H[5u * output_stride + q] = h12;
                H[6u * output_stride + q] = h02;
                H[7u * output_stride + q] = h12;
                H[8u * output_stride + q] = H00 - Real(2) * H03 + H33;
            }
        }
        return;
    }

    for (std::size_t node = 0; node < num_nodes; ++node) {
        const auto& e = simplex_exponents[node];
        const std::size_t i0 = static_cast<std::size_t>(e[0]);
        const std::size_t i1 = static_cast<std::size_t>(e[1]);
        const std::size_t i2 = static_cast<std::size_t>(e[2]);
        const std::size_t i3 = static_cast<std::size_t>(e[3]);
        Real* SVMP_RESTRICT value_row = values_out ? values_out + node * output_stride : nullptr;
        Real* SVMP_RESTRICT g = gradients_out ? gradients_out + node * 3u * output_stride : nullptr;
        Real* SVMP_RESTRICT H = hessians_out + node * 9u * output_stride;

        for (std::size_t q = 0; q < 4u; ++q) {
            const Real v0 = phi0[q][i0];
            const Real v1 = phi1[q][i1];
            const Real v2 = phi2[q][i2];
            const Real v3 = phi3[q][i3];
            if (value_row != nullptr) {
                value_row[q] = v0 * v1 * v2 * v3;
            }

            const Real D0 = dphi0[q][i0];
            const Real D1 = dphi1[q][i1];
            const Real D2 = dphi2[q][i2];
            const Real D3 = dphi3[q][i3];
            if (g != nullptr) {
                const Real dl0 = D0 * v1 * v2 * v3;
                g[0u * output_stride + q] = v0 * D1 * v2 * v3 - dl0;
                g[1u * output_stride + q] = v0 * v1 * D2 * v3 - dl0;
                g[2u * output_stride + q] = v0 * v1 * v2 * D3 - dl0;
            }

            const Real DD0 = d2phi0[q][i0];
            const Real DD1 = d2phi1[q][i1];
            const Real DD2 = d2phi2[q][i2];
            const Real DD3 = d2phi3[q][i3];
            const Real H00 = DD0 * v1 * v2 * v3;
            const Real H11 = v0 * DD1 * v2 * v3;
            const Real H22 = v0 * v1 * DD2 * v3;
            const Real H33 = v0 * v1 * v2 * DD3;
            const Real H01 = D0 * D1 * v2 * v3;
            const Real H02 = D0 * v1 * D2 * v3;
            const Real H03 = D0 * v1 * v2 * D3;
            const Real H12 = v0 * D1 * D2 * v3;
            const Real H13 = v0 * D1 * v2 * D3;
            const Real H23 = v0 * v1 * D2 * D3;
            const Real h01 = H00 - H01 - H02 + H12;
            const Real h02 = H00 - H01 - H03 + H13;
            const Real h12 = H00 - H02 - H03 + H23;
            H[0u * output_stride + q] = H00 - Real(2) * H01 + H11;
            H[1u * output_stride + q] = h01;
            H[2u * output_stride + q] = h02;
            H[3u * output_stride + q] = h01;
            H[4u * output_stride + q] = H00 - Real(2) * H02 + H22;
            H[5u * output_stride + q] = h12;
            H[6u * output_stride + q] = h02;
            H[7u * output_stride + q] = h12;
            H[8u * output_stride + q] = H00 - Real(2) * H03 + H33;
        }
    }
}

bool try_evaluate_tetrahedron_simplex_hessian_outputs_q4(
    const std::vector<std::array<int, 4>>& simplex_exponents,
    int order,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    switch (order) {
    case 2:
        evaluate_tetrahedron_simplex_hessian_outputs_q4<2>(
            simplex_exponents, points, output_stride, values_out, gradients_out, hessians_out);
        return true;
    case 3:
        evaluate_tetrahedron_simplex_hessian_outputs_q4<3>(
            simplex_exponents, points, output_stride, values_out, gradients_out, hessians_out);
        return true;
    case 4:
        evaluate_tetrahedron_simplex_hessian_outputs_q4<4>(
            simplex_exponents, points, output_stride, values_out, gradients_out, hessians_out);
        return true;
    case 5:
        evaluate_tetrahedron_simplex_hessian_outputs_q4<5>(
            simplex_exponents, points, output_stride, values_out, gradients_out, hessians_out);
        return true;
    case 6:
        evaluate_tetrahedron_simplex_hessian_outputs_q4<6>(
            simplex_exponents, points, output_stride, values_out, gradients_out, hessians_out);
        return true;
    case 7:
        evaluate_tetrahedron_simplex_hessian_outputs_q4<7>(
            simplex_exponents, points, output_stride, values_out, gradients_out, hessians_out);
        return true;
    case 8:
        evaluate_tetrahedron_simplex_hessian_outputs_q4<8>(
            simplex_exponents, points, output_stride, values_out, gradients_out, hessians_out);
        return true;
    default:
        return false;
    }
}

// Per-thread scratch space for simplex factor sequences. Common low orders use
// fixed storage; higher orders fall back to dynamic vectors.
struct SimplexAxisScratch {
    std::size_t size{0};
    std::array<Real, kFixedSimplexAxisSize> phi_fixed{};
    std::array<Real, kFixedSimplexAxisSize> dphi_fixed{};
    std::array<Real, kFixedSimplexAxisSize> d2phi_fixed{};
    std::vector<Real> phi_dynamic;
    std::vector<Real> dphi_dynamic;
    std::vector<Real> d2phi_dynamic;

    void reserveFor(std::size_t n) {
        size = n;
        if (n <= kFixedSimplexAxisSize) {
            return;
        }
        if (phi_dynamic.size() < n) phi_dynamic.resize(n);
        if (dphi_dynamic.size() < n) dphi_dynamic.resize(n);
        if (d2phi_dynamic.size() < n) d2phi_dynamic.resize(n);
    }

    Real* phi() noexcept {
        return size <= kFixedSimplexAxisSize ? phi_fixed.data() : phi_dynamic.data();
    }

    Real* dphi() noexcept {
        return size <= kFixedSimplexAxisSize ? dphi_fixed.data() : dphi_dynamic.data();
    }

    Real* d2phi() noexcept {
        return size <= kFixedSimplexAxisSize ? d2phi_fixed.data() : d2phi_dynamic.data();
    }

    const Real* phi() const noexcept {
        return size <= kFixedSimplexAxisSize ? phi_fixed.data() : phi_dynamic.data();
    }

    const Real* dphi() const noexcept {
        return size <= kFixedSimplexAxisSize ? dphi_fixed.data() : dphi_dynamic.data();
    }

    const Real* d2phi() const noexcept {
        return size <= kFixedSimplexAxisSize ? d2phi_fixed.data() : d2phi_dynamic.data();
    }
};

SimplexAxisScratch& simplex_axis_scratch_slot(int slot) {
    thread_local SimplexAxisScratch s[4];
    return s[slot];
}

struct SimplexVectorSink {
    std::vector<Real>* values;
    std::vector<Gradient>* gradients;
    std::vector<Hessian>* hessians;

    bool wants_values() const noexcept { return values != nullptr; }
    bool wants_gradients() const noexcept { return gradients != nullptr; }
    bool wants_hessians() const noexcept { return hessians != nullptr; }

    void prepare(std::size_t n_nodes) const {
        if (values)    values->resize(n_nodes);
        if (gradients) gradients->resize(n_nodes);
        if (hessians)  hessians->resize(n_nodes);
    }

    void write_value(std::size_t n, Real value) const {
        (*values)[n] = value;
    }

    void write_gradient(std::size_t n, Real x, Real y, Real z) const {
        auto& gradient = (*gradients)[n];
        gradient[0] = x;
        gradient[1] = y;
        gradient[2] = z;
    }

    void write_hessian(std::size_t n,
                       Real xx,
                       Real yy,
                       Real zz,
                       Real xy,
                       Real xz,
                       Real yz) const {
        Hessian hessian{};
        hessian(0, 0) = xx;
        hessian(1, 1) = yy;
        hessian(2, 2) = zz;
        hessian(0, 1) = xy; hessian(1, 0) = xy;
        hessian(0, 2) = xz; hessian(2, 0) = xz;
        hessian(1, 2) = yz; hessian(2, 1) = yz;
        (*hessians)[n] = hessian;
    }
};

struct SimplexRawSink {
    Real* values;
    Real* gradients;
    Real* hessians;

    bool wants_values() const noexcept { return values != nullptr; }
    bool wants_gradients() const noexcept { return gradients != nullptr; }
    bool wants_hessians() const noexcept { return hessians != nullptr; }

    void prepare(std::size_t) const {}

    void write_value(std::size_t n, Real value) const {
        values[n] = value;
    }

    void write_gradient(std::size_t n, Real x, Real y, Real z) const {
        Real* gradient = gradients + n * 3u;
        gradient[0] = x;
        gradient[1] = y;
        gradient[2] = z;
    }

    void write_hessian(std::size_t n,
                       Real xx,
                       Real yy,
                       Real zz,
                       Real xy,
                       Real xz,
                       Real yz) const {
        Real* hessian = hessians + n * 9u;
        hessian[0] = xx;
        hessian[1] = xy;
        hessian[2] = xz;
        hessian[3] = xy;
        hessian[4] = yy;
        hessian[5] = yz;
        hessian[6] = xz;
        hessian[7] = yz;
        hessian[8] = zz;
    }
};

template <typename Sink>
void evaluate_triangle_simplex_basis_impl(const std::vector<std::array<int, 4>>& simplex_exponents,
                                          int order,
                                          const math::Vector<Real, 3>& xi,
                                          const Sink& sink) {
    const Real l1 = xi[0];
    const Real l2 = xi[1];
    const Real l0 = Real(1) - l1 - l2;

    const std::size_t n = static_cast<std::size_t>(order + 1);
    SimplexAxisScratch& s0 = simplex_axis_scratch_slot(0);
    SimplexAxisScratch& s1 = simplex_axis_scratch_slot(1);
    SimplexAxisScratch& s2 = simplex_axis_scratch_slot(2);
    s0.reserveFor(n);
    s1.reserveFor(n);
    s2.reserveFor(n);

    const std::size_t num_nodes = simplex_exponents.size();
    sink.prepare(num_nodes);
    const bool need_values = sink.wants_values();
    const bool need_gradients = sink.wants_gradients();
    const bool need_hessians = sink.wants_hessians();
    Real* d0_out = (need_gradients || need_hessians) ? s0.dphi() : nullptr;
    Real* d1_out = (need_gradients || need_hessians) ? s1.dphi() : nullptr;
    Real* d2_out = (need_gradients || need_hessians) ? s2.dphi() : nullptr;
    Real* d20_out = need_hessians ? s0.d2phi() : nullptr;
    Real* d21_out = need_hessians ? s1.d2phi() : nullptr;
    Real* d22_out = need_hessians ? s2.d2phi() : nullptr;

    simplex_lagrange_factor_sequence(order, l0, s0.phi(), d0_out, d20_out);
    simplex_lagrange_factor_sequence(order, l1, s1.phi(), d1_out, d21_out);
    simplex_lagrange_factor_sequence(order, l2, s2.phi(), d2_out, d22_out);
    const Real* phi0 = s0.phi();
    const Real* phi1 = s1.phi();
    const Real* phi2 = s2.phi();
    const Real* dphi0 = s0.dphi();
    const Real* dphi1 = s1.dphi();
    const Real* dphi2 = s2.dphi();
    const Real* d2phi0 = s0.d2phi();
    const Real* d2phi1 = s1.d2phi();
    const Real* d2phi2 = s2.d2phi();

    for (std::size_t n_idx = 0; n_idx < num_nodes; ++n_idx) {
        const auto& e = simplex_exponents[n_idx];
        const std::size_t i0 = static_cast<std::size_t>(e[0]);
        const std::size_t i1 = static_cast<std::size_t>(e[1]);
        const std::size_t i2 = static_cast<std::size_t>(e[2]);

        const Real v0 = phi0[i0];
        const Real v1 = phi1[i1];
        const Real v2 = phi2[i2];
        if (need_values) {
            sink.write_value(n_idx, v0 * v1 * v2);
        }
        if (!need_gradients && !need_hessians) {
            continue;
        }

        const Real D0 = dphi0[i0];
        const Real D1 = dphi1[i1];
        const Real D2 = dphi2[i2];

        if (need_gradients) {
            const Real dl0 = D0 * v1 * v2;
            const Real dl1 = v0 * D1 * v2;
            const Real dl2 = v0 * v1 * D2;
            sink.write_gradient(n_idx, dl1 - dl0, dl2 - dl0, Real(0));
        }

        if (need_hessians) {
            const Real DD0 = d2phi0[i0];
            const Real DD1 = d2phi1[i1];
            const Real DD2 = d2phi2[i2];

            const Real H00 = DD0 * v1 * v2;
            const Real H11 = v0 * DD1 * v2;
            const Real H22 = v0 * v1 * DD2;
            const Real H01 = D0 * D1 * v2;
            const Real H02 = D0 * v1 * D2;
            const Real H12 = v0 * D1 * D2;

            sink.write_hessian(n_idx,
                               H00 - Real(2) * H01 + H11,
                               H00 - Real(2) * H02 + H22,
                               Real(0),
                               H00 - H01 - H02 + H12,
                               Real(0),
                               Real(0));
        }
    }
}

void evaluate_triangle_simplex_basis(const std::vector<std::array<int, 4>>& simplex_exponents,
                                     int order,
                                     const math::Vector<Real, 3>& xi,
                                     std::vector<Real>* values,
                                     std::vector<Gradient>* gradients,
                                     std::vector<Hessian>* hessians) {
    const SimplexVectorSink sink{values, gradients, hessians};
    evaluate_triangle_simplex_basis_impl(simplex_exponents, order, xi, sink);
}

void evaluate_triangle_simplex_basis_to(const std::vector<std::array<int, 4>>& simplex_exponents,
                                        int order,
                                        const math::Vector<Real, 3>& xi,
                                        Real* SVMP_RESTRICT values_out,
                                        Real* SVMP_RESTRICT gradients_out,
                                        Real* SVMP_RESTRICT hessians_out) {
    const SimplexRawSink sink{values_out, gradients_out, hessians_out};
    evaluate_triangle_simplex_basis_impl(simplex_exponents, order, xi, sink);
}

void evaluate_triangle_simplex_basis_strided(
    const std::vector<std::array<int, 4>>& simplex_exponents,
    int order,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    const std::size_t num_nodes = simplex_exponents.size();
    if (points.empty() || num_nodes == 0u) {
        return;
    }

    const std::size_t sequence_size = static_cast<std::size_t>(order + 1);
    const std::size_t num_qpts = points.size();
    const bool need_gradients = gradients_out != nullptr;
    const bool need_hessians = hessians_out != nullptr;
    if (num_qpts == 4u &&
        values_out != nullptr &&
        !need_gradients &&
        !need_hessians &&
        try_evaluate_triangle_simplex_values_q4(
            simplex_exponents, order, points, output_stride, values_out)) {
        return;
    }
    if (num_qpts == 4u &&
        values_out == nullptr &&
        need_gradients &&
        !need_hessians &&
        try_evaluate_triangle_simplex_gradients_q4(
            simplex_exponents, order, points, output_stride, gradients_out)) {
        return;
    }
    if (num_qpts == 4u &&
        need_hessians &&
        try_evaluate_triangle_simplex_hessian_outputs_q4(
            simplex_exponents, order, points, output_stride,
            values_out, gradients_out, hessians_out)) {
        return;
    }
    const std::size_t batch_entries = sequence_size * num_qpts;
    if (batch_entries <= kFixedSimplexBatchEntries) {
        if (values_out != nullptr && gradients_out == nullptr && hessians_out == nullptr) {
            std::array<Real, kFixedSimplexBatchEntries> phi0_batch;
            std::array<Real, kFixedSimplexBatchEntries> phi1_batch;
            std::array<Real, kFixedSimplexBatchEntries> phi2_batch;

            for (std::size_t q = 0; q < num_qpts; ++q) {
                const auto& xi = points[q];
                const Real l1 = xi[0];
                const Real l2 = xi[1];
                const Real l0 = Real(1) - l1 - l2;
                const std::size_t offset = q * sequence_size;
                simplex_lagrange_factor_sequence(
                    order, l0, phi0_batch.data() + offset, nullptr, nullptr);
                simplex_lagrange_factor_sequence(
                    order, l1, phi1_batch.data() + offset, nullptr, nullptr);
                simplex_lagrange_factor_sequence(
                    order, l2, phi2_batch.data() + offset, nullptr, nullptr);
            }

            for (std::size_t node = 0; node < num_nodes; ++node) {
                const auto& e = simplex_exponents[node];
                const std::size_t i0 = static_cast<std::size_t>(e[0]);
                const std::size_t i1 = static_cast<std::size_t>(e[1]);
                const std::size_t i2 = static_cast<std::size_t>(e[2]);
                Real* value_row = values_out + node * output_stride;

                for (std::size_t q = 0; q < num_qpts; ++q) {
                    const std::size_t offset = q * sequence_size;
                    value_row[q] =
                        phi0_batch[offset + i0] *
                        phi1_batch[offset + i1] *
                        phi2_batch[offset + i2];
                }
            }
            return;
        }

        if (values_out == nullptr && gradients_out != nullptr && hessians_out == nullptr) {
            std::array<Real, kFixedSimplexBatchEntries> phi0_batch;
            std::array<Real, kFixedSimplexBatchEntries> phi1_batch;
            std::array<Real, kFixedSimplexBatchEntries> phi2_batch;
            std::array<Real, kFixedSimplexBatchEntries> dphi0_batch;
            std::array<Real, kFixedSimplexBatchEntries> dphi1_batch;
            std::array<Real, kFixedSimplexBatchEntries> dphi2_batch;

            for (std::size_t q = 0; q < num_qpts; ++q) {
                const auto& xi = points[q];
                const Real l1 = xi[0];
                const Real l2 = xi[1];
                const Real l0 = Real(1) - l1 - l2;
                const std::size_t offset = q * sequence_size;
                simplex_lagrange_factor_sequence(
                    order, l0, phi0_batch.data() + offset, dphi0_batch.data() + offset, nullptr);
                simplex_lagrange_factor_sequence(
                    order, l1, phi1_batch.data() + offset, dphi1_batch.data() + offset, nullptr);
                simplex_lagrange_factor_sequence(
                    order, l2, phi2_batch.data() + offset, dphi2_batch.data() + offset, nullptr);
            }

            for (std::size_t node = 0; node < num_nodes; ++node) {
                const auto& e = simplex_exponents[node];
                const std::size_t i0 = static_cast<std::size_t>(e[0]);
                const std::size_t i1 = static_cast<std::size_t>(e[1]);
                const std::size_t i2 = static_cast<std::size_t>(e[2]);
                Real* g = gradients_out + node * 3u * output_stride;

                for (std::size_t q = 0; q < num_qpts; ++q) {
                    const std::size_t offset = q * sequence_size;
                    const Real v0 = phi0_batch[offset + i0];
                    const Real v1 = phi1_batch[offset + i1];
                    const Real v2 = phi2_batch[offset + i2];
                    const Real D0 = dphi0_batch[offset + i0];
                    const Real D1 = dphi1_batch[offset + i1];
                    const Real D2 = dphi2_batch[offset + i2];
                    const Real dl0 = D0 * v1 * v2;
                    g[0u * output_stride + q] = v0 * D1 * v2 - dl0;
                    g[1u * output_stride + q] = v0 * v1 * D2 - dl0;
                    g[2u * output_stride + q] = Real(0);
                }
            }
            return;
        }

        if (order >= 4 &&
            values_out == nullptr &&
            gradients_out == nullptr &&
            hessians_out != nullptr) {
            std::array<Real, kFixedSimplexBatchEntries> phi0_batch;
            std::array<Real, kFixedSimplexBatchEntries> phi1_batch;
            std::array<Real, kFixedSimplexBatchEntries> phi2_batch;
            std::array<Real, kFixedSimplexBatchEntries> dphi0_batch;
            std::array<Real, kFixedSimplexBatchEntries> dphi1_batch;
            std::array<Real, kFixedSimplexBatchEntries> dphi2_batch;
            std::array<Real, kFixedSimplexBatchEntries> d2phi0_batch;
            std::array<Real, kFixedSimplexBatchEntries> d2phi1_batch;
            std::array<Real, kFixedSimplexBatchEntries> d2phi2_batch;

            for (std::size_t q = 0; q < num_qpts; ++q) {
                const auto& xi = points[q];
                const Real l1 = xi[0];
                const Real l2 = xi[1];
                const Real l0 = Real(1) - l1 - l2;
                const std::size_t offset = q * sequence_size;
                simplex_lagrange_factor_sequence(
                    order, l0, phi0_batch.data() + offset,
                    dphi0_batch.data() + offset, d2phi0_batch.data() + offset);
                simplex_lagrange_factor_sequence(
                    order, l1, phi1_batch.data() + offset,
                    dphi1_batch.data() + offset, d2phi1_batch.data() + offset);
                simplex_lagrange_factor_sequence(
                    order, l2, phi2_batch.data() + offset,
                    dphi2_batch.data() + offset, d2phi2_batch.data() + offset);
            }

            for (std::size_t node = 0; node < num_nodes; ++node) {
                const auto& e = simplex_exponents[node];
                const std::size_t i0 = static_cast<std::size_t>(e[0]);
                const std::size_t i1 = static_cast<std::size_t>(e[1]);
                const std::size_t i2 = static_cast<std::size_t>(e[2]);
                Real* H = hessians_out + node * 9u * output_stride;

                for (std::size_t q = 0; q < num_qpts; ++q) {
                    const std::size_t offset = q * sequence_size;
                    const Real v0 = phi0_batch[offset + i0];
                    const Real v1 = phi1_batch[offset + i1];
                    const Real v2 = phi2_batch[offset + i2];
                    const Real D0 = dphi0_batch[offset + i0];
                    const Real D1 = dphi1_batch[offset + i1];
                    const Real D2 = dphi2_batch[offset + i2];
                    const Real DD0 = d2phi0_batch[offset + i0];
                    const Real DD1 = d2phi1_batch[offset + i1];
                    const Real DD2 = d2phi2_batch[offset + i2];
                    const Real H00 = DD0 * v1 * v2;
                    const Real H11 = v0 * DD1 * v2;
                    const Real H22 = v0 * v1 * DD2;
                    const Real H01 = D0 * D1 * v2;
                    const Real H02 = D0 * v1 * D2;
                    const Real H12 = v0 * D1 * D2;
                    const Real h01 = H00 - H01 - H02 + H12;

                    H[0u * output_stride + q] = H00 - Real(2) * H01 + H11;
                    H[1u * output_stride + q] = h01;
                    H[2u * output_stride + q] = Real(0);
                    H[3u * output_stride + q] = h01;
                    H[4u * output_stride + q] = H00 - Real(2) * H02 + H22;
                    H[5u * output_stride + q] = Real(0);
                    H[6u * output_stride + q] = Real(0);
                    H[7u * output_stride + q] = Real(0);
                    H[8u * output_stride + q] = Real(0);
                }
            }
            return;
        }

        std::array<Real, kFixedSimplexBatchEntries> phi0_batch;
        std::array<Real, kFixedSimplexBatchEntries> phi1_batch;
        std::array<Real, kFixedSimplexBatchEntries> phi2_batch;
        std::array<Real, kFixedSimplexBatchEntries> dphi0_batch;
        std::array<Real, kFixedSimplexBatchEntries> dphi1_batch;
        std::array<Real, kFixedSimplexBatchEntries> dphi2_batch;
        std::array<Real, kFixedSimplexBatchEntries> d2phi0_batch;
        std::array<Real, kFixedSimplexBatchEntries> d2phi1_batch;
        std::array<Real, kFixedSimplexBatchEntries> d2phi2_batch;

        for (std::size_t q = 0; q < num_qpts; ++q) {
            const auto& xi = points[q];
            const Real l1 = xi[0];
            const Real l2 = xi[1];
            const Real l0 = Real(1) - l1 - l2;
            const std::size_t offset = q * sequence_size;
            Real* d0_out = (need_gradients || need_hessians) ? dphi0_batch.data() + offset : nullptr;
            Real* d1_out = (need_gradients || need_hessians) ? dphi1_batch.data() + offset : nullptr;
            Real* d2_out = (need_gradients || need_hessians) ? dphi2_batch.data() + offset : nullptr;
            Real* d20_out = need_hessians ? d2phi0_batch.data() + offset : nullptr;
            Real* d21_out = need_hessians ? d2phi1_batch.data() + offset : nullptr;
            Real* d22_out = need_hessians ? d2phi2_batch.data() + offset : nullptr;
            simplex_lagrange_factor_sequence(order, l0, phi0_batch.data() + offset, d0_out, d20_out);
            simplex_lagrange_factor_sequence(order, l1, phi1_batch.data() + offset, d1_out, d21_out);
            simplex_lagrange_factor_sequence(order, l2, phi2_batch.data() + offset, d2_out, d22_out);
        }

        for (std::size_t node = 0; node < num_nodes; ++node) {
            const auto& e = simplex_exponents[node];
            const std::size_t i0 = static_cast<std::size_t>(e[0]);
            const std::size_t i1 = static_cast<std::size_t>(e[1]);
            const std::size_t i2 = static_cast<std::size_t>(e[2]);
            Real* value_row = values_out ? values_out + node * output_stride : nullptr;
            Real* g = gradients_out ? gradients_out + node * 3u * output_stride : nullptr;
            Real* H = hessians_out ? hessians_out + node * 9u * output_stride : nullptr;

            for (std::size_t q = 0; q < num_qpts; ++q) {
                const std::size_t offset = q * sequence_size;
                const Real v0 = phi0_batch[offset + i0];
                const Real v1 = phi1_batch[offset + i1];
                const Real v2 = phi2_batch[offset + i2];
                if (value_row != nullptr) {
                    value_row[q] = v0 * v1 * v2;
                }
                if (!need_gradients && !need_hessians) {
                    continue;
                }

                const Real D0 = dphi0_batch[offset + i0];
                const Real D1 = dphi1_batch[offset + i1];
                const Real D2 = dphi2_batch[offset + i2];

                if (gradients_out != nullptr) {
                    const Real dl0 = D0 * v1 * v2;
                    const Real dl1 = v0 * D1 * v2;
                    const Real dl2 = v0 * v1 * D2;
                    g[0u * output_stride + q] = dl1 - dl0;
                    g[1u * output_stride + q] = dl2 - dl0;
                    g[2u * output_stride + q] = Real(0);
                }

                if (hessians_out != nullptr) {
                    const Real DD0 = d2phi0_batch[offset + i0];
                    const Real DD1 = d2phi1_batch[offset + i1];
                    const Real DD2 = d2phi2_batch[offset + i2];
                    const Real H00 = DD0 * v1 * v2;
                    const Real H11 = v0 * DD1 * v2;
                    const Real H22 = v0 * v1 * DD2;
                    const Real H01 = D0 * D1 * v2;
                    const Real H02 = D0 * v1 * D2;
                    const Real H12 = v0 * D1 * D2;
                    const Real h01 = H00 - H01 - H02 + H12;
                    H[0u * output_stride + q] = H00 - Real(2) * H01 + H11;
                    H[1u * output_stride + q] = h01;
                    H[2u * output_stride + q] = Real(0);
                    H[3u * output_stride + q] = h01;
                    H[4u * output_stride + q] = H00 - Real(2) * H02 + H22;
                    H[5u * output_stride + q] = Real(0);
                    H[6u * output_stride + q] = Real(0);
                    H[7u * output_stride + q] = Real(0);
                    H[8u * output_stride + q] = Real(0);
                }
            }
        }
        return;
    }

    SimplexAxisScratch& s0 = simplex_axis_scratch_slot(0);
    SimplexAxisScratch& s1 = simplex_axis_scratch_slot(1);
    SimplexAxisScratch& s2 = simplex_axis_scratch_slot(2);
    s0.reserveFor(sequence_size);
    s1.reserveFor(sequence_size);
    s2.reserveFor(sequence_size);

    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l0 = Real(1) - l1 - l2;

        Real* d0_out = (need_gradients || need_hessians) ? s0.dphi() : nullptr;
        Real* d1_out = (need_gradients || need_hessians) ? s1.dphi() : nullptr;
        Real* d2_out = (need_gradients || need_hessians) ? s2.dphi() : nullptr;
        Real* d20_out = need_hessians ? s0.d2phi() : nullptr;
        Real* d21_out = need_hessians ? s1.d2phi() : nullptr;
        Real* d22_out = need_hessians ? s2.d2phi() : nullptr;

        simplex_lagrange_factor_sequence(order, l0, s0.phi(), d0_out, d20_out);
        simplex_lagrange_factor_sequence(order, l1, s1.phi(), d1_out, d21_out);
        simplex_lagrange_factor_sequence(order, l2, s2.phi(), d2_out, d22_out);
        const Real* phi0 = s0.phi();
        const Real* phi1 = s1.phi();
        const Real* phi2 = s2.phi();
        const Real* dphi0 = s0.dphi();
        const Real* dphi1 = s1.dphi();
        const Real* dphi2 = s2.dphi();
        const Real* d2phi0 = s0.d2phi();
        const Real* d2phi1 = s1.d2phi();
        const Real* d2phi2 = s2.d2phi();

        for (std::size_t node = 0; node < num_nodes; ++node) {
            const auto& e = simplex_exponents[node];
            const std::size_t i0 = static_cast<std::size_t>(e[0]);
            const std::size_t i1 = static_cast<std::size_t>(e[1]);
            const std::size_t i2 = static_cast<std::size_t>(e[2]);

            const Real v0 = phi0[i0];
            const Real v1 = phi1[i1];
            const Real v2 = phi2[i2];
            const Real value = v0 * v1 * v2;
            if (values_out != nullptr) {
                values_out[node * output_stride + q] = value;
            }
            if (!need_gradients && !need_hessians) {
                continue;
            }

            const Real D0 = dphi0[i0];
            const Real D1 = dphi1[i1];
            const Real D2 = dphi2[i2];

            if (gradients_out != nullptr) {
                const Real dl0 = D0 * v1 * v2;
                const Real dl1 = v0 * D1 * v2;
                const Real dl2 = v0 * v1 * D2;
                Real* g = gradients_out + node * 3u * output_stride;
                g[0u * output_stride + q] = dl1 - dl0;
                g[1u * output_stride + q] = dl2 - dl0;
                g[2u * output_stride + q] = Real(0);
            }

            if (hessians_out != nullptr) {
                const Real DD0 = d2phi0[i0];
                const Real DD1 = d2phi1[i1];
                const Real DD2 = d2phi2[i2];

                const Real H00 = DD0 * v1 * v2;
                const Real H11 = v0 * DD1 * v2;
                const Real H22 = v0 * v1 * DD2;
                const Real H01 = D0 * D1 * v2;
                const Real H02 = D0 * v1 * D2;
                const Real H12 = v0 * D1 * D2;

                Real* H = hessians_out + node * 9u * output_stride;
                const Real h01 = H00 - H01 - H02 + H12;
                H[0u * output_stride + q] = H00 - Real(2) * H01 + H11;
                H[1u * output_stride + q] = h01;
                H[2u * output_stride + q] = Real(0);
                H[3u * output_stride + q] = h01;
                H[4u * output_stride + q] = H00 - Real(2) * H02 + H22;
                H[5u * output_stride + q] = Real(0);
                H[6u * output_stride + q] = Real(0);
                H[7u * output_stride + q] = Real(0);
                H[8u * output_stride + q] = Real(0);
            }
        }
    }
}

void evaluate_triangle_simplex_basis_wedge_components_strided(
    const std::vector<std::array<int, 4>>& simplex_exponents,
    int order,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_xy_out,
    Real* SVMP_RESTRICT hessians_xx_xy_yy_out) {
    const std::size_t num_nodes = simplex_exponents.size();
    if (points.empty() || num_nodes == 0u) {
        return;
    }

    const std::size_t sequence_size = static_cast<std::size_t>(order + 1);
    const std::size_t num_qpts = points.size();
    const bool need_gradients = gradients_xy_out != nullptr;
    const bool need_hessians = hessians_xx_xy_yy_out != nullptr;
    const std::size_t batch_entries = sequence_size * num_qpts;

    if (batch_entries <= kFixedSimplexBatchEntries) {
        if (values_out != nullptr &&
            gradients_xy_out != nullptr &&
            hessians_xx_xy_yy_out == nullptr) {
            std::array<Real, kFixedSimplexBatchEntries> phi0_batch;
            std::array<Real, kFixedSimplexBatchEntries> phi1_batch;
            std::array<Real, kFixedSimplexBatchEntries> phi2_batch;
            std::array<Real, kFixedSimplexBatchEntries> dphi0_batch;
            std::array<Real, kFixedSimplexBatchEntries> dphi1_batch;
            std::array<Real, kFixedSimplexBatchEntries> dphi2_batch;

            for (std::size_t q = 0; q < num_qpts; ++q) {
                const auto& xi = points[q];
                const Real l1 = xi[0];
                const Real l2 = xi[1];
                const Real l0 = Real(1) - l1 - l2;
                const std::size_t offset = q * sequence_size;
                simplex_lagrange_factor_sequence(
                    order, l0, phi0_batch.data() + offset, dphi0_batch.data() + offset, nullptr);
                simplex_lagrange_factor_sequence(
                    order, l1, phi1_batch.data() + offset, dphi1_batch.data() + offset, nullptr);
                simplex_lagrange_factor_sequence(
                    order, l2, phi2_batch.data() + offset, dphi2_batch.data() + offset, nullptr);
            }

            for (std::size_t node = 0; node < num_nodes; ++node) {
                const auto& e = simplex_exponents[node];
                const std::size_t i0 = static_cast<std::size_t>(e[0]);
                const std::size_t i1 = static_cast<std::size_t>(e[1]);
                const std::size_t i2 = static_cast<std::size_t>(e[2]);
                Real* value_row = values_out + node * output_stride;
                Real* g = gradients_xy_out + node * 2u * output_stride;

                for (std::size_t q = 0; q < num_qpts; ++q) {
                    const std::size_t offset = q * sequence_size;
                    const Real v0 = phi0_batch[offset + i0];
                    const Real v1 = phi1_batch[offset + i1];
                    const Real v2 = phi2_batch[offset + i2];
                    const Real D0 = dphi0_batch[offset + i0];
                    const Real D1 = dphi1_batch[offset + i1];
                    const Real D2 = dphi2_batch[offset + i2];
                    const Real dl0 = D0 * v1 * v2;
                    value_row[q] = v0 * v1 * v2;
                    g[0u * output_stride + q] = v0 * D1 * v2 - dl0;
                    g[1u * output_stride + q] = v0 * v1 * D2 - dl0;
                }
            }
            return;
        }

        if (values_out != nullptr &&
            gradients_xy_out != nullptr &&
            hessians_xx_xy_yy_out != nullptr) {
            std::array<Real, kFixedSimplexBatchEntries> phi0_batch;
            std::array<Real, kFixedSimplexBatchEntries> phi1_batch;
            std::array<Real, kFixedSimplexBatchEntries> phi2_batch;
            std::array<Real, kFixedSimplexBatchEntries> dphi0_batch;
            std::array<Real, kFixedSimplexBatchEntries> dphi1_batch;
            std::array<Real, kFixedSimplexBatchEntries> dphi2_batch;
            std::array<Real, kFixedSimplexBatchEntries> d2phi0_batch;
            std::array<Real, kFixedSimplexBatchEntries> d2phi1_batch;
            std::array<Real, kFixedSimplexBatchEntries> d2phi2_batch;

            for (std::size_t q = 0; q < num_qpts; ++q) {
                const auto& xi = points[q];
                const Real l1 = xi[0];
                const Real l2 = xi[1];
                const Real l0 = Real(1) - l1 - l2;
                const std::size_t offset = q * sequence_size;
                simplex_lagrange_factor_sequence_impl<true, true>(
                    order, l0, phi0_batch.data() + offset,
                    dphi0_batch.data() + offset, d2phi0_batch.data() + offset);
                simplex_lagrange_factor_sequence_impl<true, true>(
                    order, l1, phi1_batch.data() + offset,
                    dphi1_batch.data() + offset, d2phi1_batch.data() + offset);
                simplex_lagrange_factor_sequence_impl<true, true>(
                    order, l2, phi2_batch.data() + offset,
                    dphi2_batch.data() + offset, d2phi2_batch.data() + offset);
            }

            for (std::size_t node = 0; node < num_nodes; ++node) {
                const auto& e = simplex_exponents[node];
                const std::size_t i0 = static_cast<std::size_t>(e[0]);
                const std::size_t i1 = static_cast<std::size_t>(e[1]);
                const std::size_t i2 = static_cast<std::size_t>(e[2]);
                Real* SVMP_RESTRICT value_row = values_out + node * output_stride;
                Real* SVMP_RESTRICT g = gradients_xy_out + node * 2u * output_stride;
                Real* SVMP_RESTRICT H = hessians_xx_xy_yy_out + node * 3u * output_stride;

                for (std::size_t q = 0; q < num_qpts; ++q) {
                    const std::size_t offset = q * sequence_size;
                    const Real v0 = phi0_batch[offset + i0];
                    const Real v1 = phi1_batch[offset + i1];
                    const Real v2 = phi2_batch[offset + i2];
                    const Real D0 = dphi0_batch[offset + i0];
                    const Real D1 = dphi1_batch[offset + i1];
                    const Real D2 = dphi2_batch[offset + i2];
                    const Real dl0 = D0 * v1 * v2;
                    const Real dl1 = v0 * D1 * v2;
                    const Real dl2 = v0 * v1 * D2;
                    const Real DD0 = d2phi0_batch[offset + i0];
                    const Real DD1 = d2phi1_batch[offset + i1];
                    const Real DD2 = d2phi2_batch[offset + i2];
                    const Real H00 = DD0 * v1 * v2;
                    const Real H11 = v0 * DD1 * v2;
                    const Real H22 = v0 * v1 * DD2;
                    const Real H01 = D0 * D1 * v2;
                    const Real H02 = D0 * v1 * D2;
                    const Real H12 = v0 * D1 * D2;

                    value_row[q] = v0 * v1 * v2;
                    g[0u * output_stride + q] = dl1 - dl0;
                    g[1u * output_stride + q] = dl2 - dl0;
                    H[0u * output_stride + q] = H00 - Real(2) * H01 + H11;
                    H[1u * output_stride + q] = H00 - H01 - H02 + H12;
                    H[2u * output_stride + q] = H00 - Real(2) * H02 + H22;
                }
            }
            return;
        }

        std::array<Real, kFixedSimplexBatchEntries> phi0_batch;
        std::array<Real, kFixedSimplexBatchEntries> phi1_batch;
        std::array<Real, kFixedSimplexBatchEntries> phi2_batch;
        std::array<Real, kFixedSimplexBatchEntries> dphi0_batch;
        std::array<Real, kFixedSimplexBatchEntries> dphi1_batch;
        std::array<Real, kFixedSimplexBatchEntries> dphi2_batch;
        std::array<Real, kFixedSimplexBatchEntries> d2phi0_batch;
        std::array<Real, kFixedSimplexBatchEntries> d2phi1_batch;
        std::array<Real, kFixedSimplexBatchEntries> d2phi2_batch;

        for (std::size_t q = 0; q < num_qpts; ++q) {
            const auto& xi = points[q];
            const Real l1 = xi[0];
            const Real l2 = xi[1];
            const Real l0 = Real(1) - l1 - l2;
            const std::size_t offset = q * sequence_size;
            Real* d0_out = (need_gradients || need_hessians) ? dphi0_batch.data() + offset : nullptr;
            Real* d1_out = (need_gradients || need_hessians) ? dphi1_batch.data() + offset : nullptr;
            Real* d2_out = (need_gradients || need_hessians) ? dphi2_batch.data() + offset : nullptr;
            Real* d20_out = need_hessians ? d2phi0_batch.data() + offset : nullptr;
            Real* d21_out = need_hessians ? d2phi1_batch.data() + offset : nullptr;
            Real* d22_out = need_hessians ? d2phi2_batch.data() + offset : nullptr;
            simplex_lagrange_factor_sequence(order, l0, phi0_batch.data() + offset, d0_out, d20_out);
            simplex_lagrange_factor_sequence(order, l1, phi1_batch.data() + offset, d1_out, d21_out);
            simplex_lagrange_factor_sequence(order, l2, phi2_batch.data() + offset, d2_out, d22_out);
        }

        for (std::size_t node = 0; node < num_nodes; ++node) {
            const auto& e = simplex_exponents[node];
            const std::size_t i0 = static_cast<std::size_t>(e[0]);
            const std::size_t i1 = static_cast<std::size_t>(e[1]);
            const std::size_t i2 = static_cast<std::size_t>(e[2]);
            Real* value_row = values_out ? values_out + node * output_stride : nullptr;
            Real* g = gradients_xy_out ? gradients_xy_out + node * 2u * output_stride : nullptr;
            Real* H = hessians_xx_xy_yy_out ? hessians_xx_xy_yy_out + node * 3u * output_stride : nullptr;

            for (std::size_t q = 0; q < num_qpts; ++q) {
                const std::size_t offset = q * sequence_size;
                const Real v0 = phi0_batch[offset + i0];
                const Real v1 = phi1_batch[offset + i1];
                const Real v2 = phi2_batch[offset + i2];
                if (value_row != nullptr) {
                    value_row[q] = v0 * v1 * v2;
                }
                if (!need_gradients && !need_hessians) {
                    continue;
                }

                const Real D0 = dphi0_batch[offset + i0];
                const Real D1 = dphi1_batch[offset + i1];
                const Real D2 = dphi2_batch[offset + i2];
                const Real dl0 = D0 * v1 * v2;
                const Real dl1 = v0 * D1 * v2;
                const Real dl2 = v0 * v1 * D2;

                if (gradients_xy_out != nullptr) {
                    g[0u * output_stride + q] = dl1 - dl0;
                    g[1u * output_stride + q] = dl2 - dl0;
                }

                if (hessians_xx_xy_yy_out != nullptr) {
                    const Real DD0 = d2phi0_batch[offset + i0];
                    const Real DD1 = d2phi1_batch[offset + i1];
                    const Real DD2 = d2phi2_batch[offset + i2];
                    const Real H00 = DD0 * v1 * v2;
                    const Real H11 = v0 * DD1 * v2;
                    const Real H22 = v0 * v1 * DD2;
                    const Real H01 = D0 * D1 * v2;
                    const Real H02 = D0 * v1 * D2;
                    const Real H12 = v0 * D1 * D2;
                    H[0u * output_stride + q] = H00 - Real(2) * H01 + H11;
                    H[1u * output_stride + q] = H00 - H01 - H02 + H12;
                    H[2u * output_stride + q] = H00 - Real(2) * H02 + H22;
                }
            }
        }
        return;
    }

    SimplexAxisScratch& s0 = simplex_axis_scratch_slot(0);
    SimplexAxisScratch& s1 = simplex_axis_scratch_slot(1);
    SimplexAxisScratch& s2 = simplex_axis_scratch_slot(2);
    s0.reserveFor(sequence_size);
    s1.reserveFor(sequence_size);
    s2.reserveFor(sequence_size);

    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l0 = Real(1) - l1 - l2;

        Real* d0_out = (need_gradients || need_hessians) ? s0.dphi() : nullptr;
        Real* d1_out = (need_gradients || need_hessians) ? s1.dphi() : nullptr;
        Real* d2_out = (need_gradients || need_hessians) ? s2.dphi() : nullptr;
        Real* d20_out = need_hessians ? s0.d2phi() : nullptr;
        Real* d21_out = need_hessians ? s1.d2phi() : nullptr;
        Real* d22_out = need_hessians ? s2.d2phi() : nullptr;
        simplex_lagrange_factor_sequence(order, l0, s0.phi(), d0_out, d20_out);
        simplex_lagrange_factor_sequence(order, l1, s1.phi(), d1_out, d21_out);
        simplex_lagrange_factor_sequence(order, l2, s2.phi(), d2_out, d22_out);

        const Real* phi0 = s0.phi();
        const Real* phi1 = s1.phi();
        const Real* phi2 = s2.phi();
        const Real* dphi0 = s0.dphi();
        const Real* dphi1 = s1.dphi();
        const Real* dphi2 = s2.dphi();
        const Real* d2phi0 = s0.d2phi();
        const Real* d2phi1 = s1.d2phi();
        const Real* d2phi2 = s2.d2phi();

        for (std::size_t node = 0; node < num_nodes; ++node) {
            const auto& e = simplex_exponents[node];
            const std::size_t i0 = static_cast<std::size_t>(e[0]);
            const std::size_t i1 = static_cast<std::size_t>(e[1]);
            const std::size_t i2 = static_cast<std::size_t>(e[2]);
            const Real v0 = phi0[i0];
            const Real v1 = phi1[i1];
            const Real v2 = phi2[i2];

            if (values_out != nullptr) {
                values_out[node * output_stride + q] = v0 * v1 * v2;
            }
            if (!need_gradients && !need_hessians) {
                continue;
            }

            const Real D0 = dphi0[i0];
            const Real D1 = dphi1[i1];
            const Real D2 = dphi2[i2];
            const Real dl0 = D0 * v1 * v2;
            const Real dl1 = v0 * D1 * v2;
            const Real dl2 = v0 * v1 * D2;

            if (gradients_xy_out != nullptr) {
                Real* g = gradients_xy_out + node * 2u * output_stride;
                g[0u * output_stride + q] = dl1 - dl0;
                g[1u * output_stride + q] = dl2 - dl0;
            }

            if (hessians_xx_xy_yy_out != nullptr) {
                const Real DD0 = d2phi0[i0];
                const Real DD1 = d2phi1[i1];
                const Real DD2 = d2phi2[i2];
                const Real H00 = DD0 * v1 * v2;
                const Real H11 = v0 * DD1 * v2;
                const Real H22 = v0 * v1 * DD2;
                const Real H01 = D0 * D1 * v2;
                const Real H02 = D0 * v1 * D2;
                const Real H12 = v0 * D1 * D2;
                Real* H = hessians_xx_xy_yy_out + node * 3u * output_stride;
                H[0u * output_stride + q] = H00 - Real(2) * H01 + H11;
                H[1u * output_stride + q] = H00 - H01 - H02 + H12;
                H[2u * output_stride + q] = H00 - Real(2) * H02 + H22;
            }
        }
    }
}

template <typename Sink>
void evaluate_tetrahedron_simplex_basis_impl(const std::vector<std::array<int, 4>>& simplex_exponents,
                                             int order,
                                             const math::Vector<Real, 3>& xi,
                                             const Sink& sink) {
    const Real l1 = xi[0];
    const Real l2 = xi[1];
    const Real l3 = xi[2];
    const Real l0 = Real(1) - l1 - l2 - l3;

    const std::size_t n = static_cast<std::size_t>(order + 1);
    SimplexAxisScratch& s0 = simplex_axis_scratch_slot(0);
    SimplexAxisScratch& s1 = simplex_axis_scratch_slot(1);
    SimplexAxisScratch& s2 = simplex_axis_scratch_slot(2);
    SimplexAxisScratch& s3 = simplex_axis_scratch_slot(3);
    s0.reserveFor(n);
    s1.reserveFor(n);
    s2.reserveFor(n);
    s3.reserveFor(n);

    const std::size_t num_nodes = simplex_exponents.size();
    sink.prepare(num_nodes);
    const bool need_values = sink.wants_values();
    const bool need_gradients = sink.wants_gradients();
    const bool need_hessians = sink.wants_hessians();
    Real* d0_out = (need_gradients || need_hessians) ? s0.dphi() : nullptr;
    Real* d1_out = (need_gradients || need_hessians) ? s1.dphi() : nullptr;
    Real* d2_out = (need_gradients || need_hessians) ? s2.dphi() : nullptr;
    Real* d3_out = (need_gradients || need_hessians) ? s3.dphi() : nullptr;
    Real* d20_out = need_hessians ? s0.d2phi() : nullptr;
    Real* d21_out = need_hessians ? s1.d2phi() : nullptr;
    Real* d22_out = need_hessians ? s2.d2phi() : nullptr;
    Real* d23_out = need_hessians ? s3.d2phi() : nullptr;

    simplex_lagrange_factor_sequence(order, l0, s0.phi(), d0_out, d20_out);
    simplex_lagrange_factor_sequence(order, l1, s1.phi(), d1_out, d21_out);
    simplex_lagrange_factor_sequence(order, l2, s2.phi(), d2_out, d22_out);
    simplex_lagrange_factor_sequence(order, l3, s3.phi(), d3_out, d23_out);
    const Real* phi0 = s0.phi();
    const Real* phi1 = s1.phi();
    const Real* phi2 = s2.phi();
    const Real* phi3 = s3.phi();
    const Real* dphi0 = s0.dphi();
    const Real* dphi1 = s1.dphi();
    const Real* dphi2 = s2.dphi();
    const Real* dphi3 = s3.dphi();
    const Real* d2phi0 = s0.d2phi();
    const Real* d2phi1 = s1.d2phi();
    const Real* d2phi2 = s2.d2phi();
    const Real* d2phi3 = s3.d2phi();

    for (std::size_t n_idx = 0; n_idx < num_nodes; ++n_idx) {
        const auto& e = simplex_exponents[n_idx];
        const std::size_t i0 = static_cast<std::size_t>(e[0]);
        const std::size_t i1 = static_cast<std::size_t>(e[1]);
        const std::size_t i2 = static_cast<std::size_t>(e[2]);
        const std::size_t i3 = static_cast<std::size_t>(e[3]);

        const Real v0 = phi0[i0];
        const Real v1 = phi1[i1];
        const Real v2 = phi2[i2];
        const Real v3 = phi3[i3];
        if (need_values) {
            sink.write_value(n_idx, v0 * v1 * v2 * v3);
        }
        if (!need_gradients && !need_hessians) {
            continue;
        }

        const Real D0 = dphi0[i0];
        const Real D1 = dphi1[i1];
        const Real D2 = dphi2[i2];
        const Real D3 = dphi3[i3];

        if (need_gradients) {
            const Real dl0 = D0 * v1 * v2 * v3;
            const Real dl1 = v0 * D1 * v2 * v3;
            const Real dl2 = v0 * v1 * D2 * v3;
            const Real dl3 = v0 * v1 * v2 * D3;
            sink.write_gradient(n_idx, dl1 - dl0, dl2 - dl0, dl3 - dl0);
        }

        if (need_hessians) {
            const Real DD0 = d2phi0[i0];
            const Real DD1 = d2phi1[i1];
            const Real DD2 = d2phi2[i2];
            const Real DD3 = d2phi3[i3];

            const Real H00 = DD0 * v1 * v2 * v3;
            const Real H11 = v0 * DD1 * v2 * v3;
            const Real H22 = v0 * v1 * DD2 * v3;
            const Real H33 = v0 * v1 * v2 * DD3;

            const Real H01 = D0 * D1 * v2 * v3;
            const Real H02 = D0 * v1 * D2 * v3;
            const Real H03 = D0 * v1 * v2 * D3;
            const Real H12 = v0 * D1 * D2 * v3;
            const Real H13 = v0 * D1 * v2 * D3;
            const Real H23 = v0 * v1 * D2 * D3;

            sink.write_hessian(n_idx,
                               H00 - Real(2) * H01 + H11,
                               H00 - Real(2) * H02 + H22,
                               H00 - Real(2) * H03 + H33,
                               H00 - H01 - H02 + H12,
                               H00 - H01 - H03 + H13,
                               H00 - H02 - H03 + H23);
        }
    }
}

void evaluate_tetrahedron_simplex_basis(const std::vector<std::array<int, 4>>& simplex_exponents,
                                        int order,
                                        const math::Vector<Real, 3>& xi,
                                        std::vector<Real>* values,
                                        std::vector<Gradient>* gradients,
                                        std::vector<Hessian>* hessians) {
    const SimplexVectorSink sink{values, gradients, hessians};
    evaluate_tetrahedron_simplex_basis_impl(simplex_exponents, order, xi, sink);
}

void evaluate_tetrahedron_simplex_basis_to(const std::vector<std::array<int, 4>>& simplex_exponents,
                                           int order,
                                           const math::Vector<Real, 3>& xi,
                                           Real* SVMP_RESTRICT values_out,
                                           Real* SVMP_RESTRICT gradients_out,
                                           Real* SVMP_RESTRICT hessians_out) {
    const SimplexRawSink sink{values_out, gradients_out, hessians_out};
    evaluate_tetrahedron_simplex_basis_impl(simplex_exponents, order, xi, sink);
}

void evaluate_tetrahedron_simplex_basis_strided(
    const std::vector<std::array<int, 4>>& simplex_exponents,
    int order,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    const std::size_t num_nodes = simplex_exponents.size();
    if (points.empty() || num_nodes == 0u) {
        return;
    }

    const std::size_t sequence_size = static_cast<std::size_t>(order + 1);
    const std::size_t num_qpts = points.size();
    const bool need_gradients = gradients_out != nullptr;
    const bool need_hessians = hessians_out != nullptr;
    if (num_qpts == 4u &&
        values_out != nullptr &&
        !need_gradients &&
        !need_hessians &&
        try_evaluate_tetrahedron_simplex_values_q4(
            simplex_exponents, order, points, output_stride, values_out)) {
        return;
    }
    if (num_qpts == 4u &&
        values_out == nullptr &&
        need_gradients &&
        !need_hessians) {
        switch (order) {
        case 3:
            evaluate_tetrahedron_simplex_gradients_q4<3>(
                simplex_exponents, points, output_stride, gradients_out);
            return;
        case 4:
            evaluate_tetrahedron_simplex_gradients_q4<4>(
                simplex_exponents, points, output_stride, gradients_out);
            return;
        case 5:
            evaluate_tetrahedron_simplex_gradients_q4<5>(
                simplex_exponents, points, output_stride, gradients_out);
            return;
        case 6:
            evaluate_tetrahedron_simplex_gradients_q4<6>(
                simplex_exponents, points, output_stride, gradients_out);
            return;
        case 7:
            evaluate_tetrahedron_simplex_gradients_q4<7>(
                simplex_exponents, points, output_stride, gradients_out);
            return;
        case 8:
            evaluate_tetrahedron_simplex_gradients_q4<8>(
                simplex_exponents, points, output_stride, gradients_out);
            return;
        default:
            break;
        }
    }
    if (num_qpts == 4u &&
        need_hessians &&
        try_evaluate_tetrahedron_simplex_hessian_outputs_q4(
            simplex_exponents, order, points, output_stride,
            values_out, gradients_out, hessians_out)) {
        return;
    }
    const std::size_t batch_entries = sequence_size * num_qpts;
    if (batch_entries <= kFixedSimplexBatchEntries) {
        if (values_out != nullptr && gradients_out == nullptr && hessians_out == nullptr) {
            std::array<Real, kFixedSimplexBatchEntries> phi0_batch;
            std::array<Real, kFixedSimplexBatchEntries> phi1_batch;
            std::array<Real, kFixedSimplexBatchEntries> phi2_batch;
            std::array<Real, kFixedSimplexBatchEntries> phi3_batch;

            for (std::size_t q = 0; q < num_qpts; ++q) {
                const auto& xi = points[q];
                const Real l1 = xi[0];
                const Real l2 = xi[1];
                const Real l3 = xi[2];
                const Real l0 = Real(1) - l1 - l2 - l3;
                const std::size_t offset = q * sequence_size;
                simplex_lagrange_factor_sequence(
                    order, l0, phi0_batch.data() + offset, nullptr, nullptr);
                simplex_lagrange_factor_sequence(
                    order, l1, phi1_batch.data() + offset, nullptr, nullptr);
                simplex_lagrange_factor_sequence(
                    order, l2, phi2_batch.data() + offset, nullptr, nullptr);
                simplex_lagrange_factor_sequence(
                    order, l3, phi3_batch.data() + offset, nullptr, nullptr);
            }

            for (std::size_t node = 0; node < num_nodes; ++node) {
                const auto& e = simplex_exponents[node];
                const std::size_t i0 = static_cast<std::size_t>(e[0]);
                const std::size_t i1 = static_cast<std::size_t>(e[1]);
                const std::size_t i2 = static_cast<std::size_t>(e[2]);
                const std::size_t i3 = static_cast<std::size_t>(e[3]);
                Real* value_row = values_out + node * output_stride;

                for (std::size_t q = 0; q < num_qpts; ++q) {
                    const std::size_t offset = q * sequence_size;
                    value_row[q] =
                        phi0_batch[offset + i0] *
                        phi1_batch[offset + i1] *
                        phi2_batch[offset + i2] *
                        phi3_batch[offset + i3];
                }
            }
            return;
        }

        if (values_out == nullptr && gradients_out != nullptr && hessians_out == nullptr) {
            std::array<Real, kFixedSimplexBatchEntries> phi0_batch;
            std::array<Real, kFixedSimplexBatchEntries> phi1_batch;
            std::array<Real, kFixedSimplexBatchEntries> phi2_batch;
            std::array<Real, kFixedSimplexBatchEntries> phi3_batch;
            std::array<Real, kFixedSimplexBatchEntries> dphi0_batch;
            std::array<Real, kFixedSimplexBatchEntries> dphi1_batch;
            std::array<Real, kFixedSimplexBatchEntries> dphi2_batch;
            std::array<Real, kFixedSimplexBatchEntries> dphi3_batch;

            for (std::size_t q = 0; q < num_qpts; ++q) {
                const auto& xi = points[q];
                const Real l1 = xi[0];
                const Real l2 = xi[1];
                const Real l3 = xi[2];
                const Real l0 = Real(1) - l1 - l2 - l3;
                const std::size_t offset = q * sequence_size;
                simplex_lagrange_factor_sequence(
                    order, l0, phi0_batch.data() + offset, dphi0_batch.data() + offset, nullptr);
                simplex_lagrange_factor_sequence(
                    order, l1, phi1_batch.data() + offset, dphi1_batch.data() + offset, nullptr);
                simplex_lagrange_factor_sequence(
                    order, l2, phi2_batch.data() + offset, dphi2_batch.data() + offset, nullptr);
                simplex_lagrange_factor_sequence(
                    order, l3, phi3_batch.data() + offset, dphi3_batch.data() + offset, nullptr);
            }

            for (std::size_t node = 0; node < num_nodes; ++node) {
                const auto& e = simplex_exponents[node];
                const std::size_t i0 = static_cast<std::size_t>(e[0]);
                const std::size_t i1 = static_cast<std::size_t>(e[1]);
                const std::size_t i2 = static_cast<std::size_t>(e[2]);
                const std::size_t i3 = static_cast<std::size_t>(e[3]);
                Real* g = gradients_out + node * 3u * output_stride;

                for (std::size_t q = 0; q < num_qpts; ++q) {
                    const std::size_t offset = q * sequence_size;
                    const Real v0 = phi0_batch[offset + i0];
                    const Real v1 = phi1_batch[offset + i1];
                    const Real v2 = phi2_batch[offset + i2];
                    const Real v3 = phi3_batch[offset + i3];
                    const Real D0 = dphi0_batch[offset + i0];
                    const Real D1 = dphi1_batch[offset + i1];
                    const Real D2 = dphi2_batch[offset + i2];
                    const Real D3 = dphi3_batch[offset + i3];
                    const Real v23 = v2 * v3;
                    const Real dl0 = D0 * v1 * v23;
                    g[0u * output_stride + q] = v0 * D1 * v23 - dl0;
                    g[1u * output_stride + q] = v0 * v1 * D2 * v3 - dl0;
                    g[2u * output_stride + q] = v0 * v1 * v2 * D3 - dl0;
                }
            }
            return;
        }

        std::array<Real, kFixedSimplexBatchEntries> phi0_batch;
        std::array<Real, kFixedSimplexBatchEntries> phi1_batch;
        std::array<Real, kFixedSimplexBatchEntries> phi2_batch;
        std::array<Real, kFixedSimplexBatchEntries> phi3_batch;
        std::array<Real, kFixedSimplexBatchEntries> dphi0_batch;
        std::array<Real, kFixedSimplexBatchEntries> dphi1_batch;
        std::array<Real, kFixedSimplexBatchEntries> dphi2_batch;
        std::array<Real, kFixedSimplexBatchEntries> dphi3_batch;
        std::array<Real, kFixedSimplexBatchEntries> d2phi0_batch;
        std::array<Real, kFixedSimplexBatchEntries> d2phi1_batch;
        std::array<Real, kFixedSimplexBatchEntries> d2phi2_batch;
        std::array<Real, kFixedSimplexBatchEntries> d2phi3_batch;

        for (std::size_t q = 0; q < num_qpts; ++q) {
            const auto& xi = points[q];
            const Real l1 = xi[0];
            const Real l2 = xi[1];
            const Real l3 = xi[2];
            const Real l0 = Real(1) - l1 - l2 - l3;
            const std::size_t offset = q * sequence_size;
            Real* d0_out = (need_gradients || need_hessians) ? dphi0_batch.data() + offset : nullptr;
            Real* d1_out = (need_gradients || need_hessians) ? dphi1_batch.data() + offset : nullptr;
            Real* d2_out = (need_gradients || need_hessians) ? dphi2_batch.data() + offset : nullptr;
            Real* d3_out = (need_gradients || need_hessians) ? dphi3_batch.data() + offset : nullptr;
            Real* d20_out = need_hessians ? d2phi0_batch.data() + offset : nullptr;
            Real* d21_out = need_hessians ? d2phi1_batch.data() + offset : nullptr;
            Real* d22_out = need_hessians ? d2phi2_batch.data() + offset : nullptr;
            Real* d23_out = need_hessians ? d2phi3_batch.data() + offset : nullptr;
            simplex_lagrange_factor_sequence(order, l0, phi0_batch.data() + offset, d0_out, d20_out);
            simplex_lagrange_factor_sequence(order, l1, phi1_batch.data() + offset, d1_out, d21_out);
            simplex_lagrange_factor_sequence(order, l2, phi2_batch.data() + offset, d2_out, d22_out);
            simplex_lagrange_factor_sequence(order, l3, phi3_batch.data() + offset, d3_out, d23_out);
        }

        for (std::size_t node = 0; node < num_nodes; ++node) {
            const auto& e = simplex_exponents[node];
            const std::size_t i0 = static_cast<std::size_t>(e[0]);
            const std::size_t i1 = static_cast<std::size_t>(e[1]);
            const std::size_t i2 = static_cast<std::size_t>(e[2]);
            const std::size_t i3 = static_cast<std::size_t>(e[3]);
            Real* value_row = values_out ? values_out + node * output_stride : nullptr;
            Real* g = gradients_out ? gradients_out + node * 3u * output_stride : nullptr;
            Real* H = hessians_out ? hessians_out + node * 9u * output_stride : nullptr;

            for (std::size_t q = 0; q < num_qpts; ++q) {
                const std::size_t offset = q * sequence_size;
                const Real v0 = phi0_batch[offset + i0];
                const Real v1 = phi1_batch[offset + i1];
                const Real v2 = phi2_batch[offset + i2];
                const Real v3 = phi3_batch[offset + i3];
                if (value_row != nullptr) {
                    value_row[q] = v0 * v1 * v2 * v3;
                }
                if (!need_gradients && !need_hessians) {
                    continue;
                }

                const Real D0 = dphi0_batch[offset + i0];
                const Real D1 = dphi1_batch[offset + i1];
                const Real D2 = dphi2_batch[offset + i2];
                const Real D3 = dphi3_batch[offset + i3];

                if (gradients_out != nullptr) {
                    const Real dl0 = D0 * v1 * v2 * v3;
                    const Real dl1 = v0 * D1 * v2 * v3;
                    const Real dl2 = v0 * v1 * D2 * v3;
                    const Real dl3 = v0 * v1 * v2 * D3;
                    g[0u * output_stride + q] = dl1 - dl0;
                    g[1u * output_stride + q] = dl2 - dl0;
                    g[2u * output_stride + q] = dl3 - dl0;
                }

                if (hessians_out != nullptr) {
                    const Real DD0 = d2phi0_batch[offset + i0];
                    const Real DD1 = d2phi1_batch[offset + i1];
                    const Real DD2 = d2phi2_batch[offset + i2];
                    const Real DD3 = d2phi3_batch[offset + i3];
                    const Real H00 = DD0 * v1 * v2 * v3;
                    const Real H11 = v0 * DD1 * v2 * v3;
                    const Real H22 = v0 * v1 * DD2 * v3;
                    const Real H33 = v0 * v1 * v2 * DD3;
                    const Real H01 = D0 * D1 * v2 * v3;
                    const Real H02 = D0 * v1 * D2 * v3;
                    const Real H03 = D0 * v1 * v2 * D3;
                    const Real H12 = v0 * D1 * D2 * v3;
                    const Real H13 = v0 * D1 * v2 * D3;
                    const Real H23 = v0 * v1 * D2 * D3;
                    const Real h01 = H00 - H01 - H02 + H12;
                    const Real h02 = H00 - H01 - H03 + H13;
                    const Real h12 = H00 - H02 - H03 + H23;
                    H[0u * output_stride + q] = H00 - Real(2) * H01 + H11;
                    H[1u * output_stride + q] = h01;
                    H[2u * output_stride + q] = h02;
                    H[3u * output_stride + q] = h01;
                    H[4u * output_stride + q] = H00 - Real(2) * H02 + H22;
                    H[5u * output_stride + q] = h12;
                    H[6u * output_stride + q] = h02;
                    H[7u * output_stride + q] = h12;
                    H[8u * output_stride + q] = H00 - Real(2) * H03 + H33;
                }
            }
        }
        return;
    }

    SimplexAxisScratch& s0 = simplex_axis_scratch_slot(0);
    SimplexAxisScratch& s1 = simplex_axis_scratch_slot(1);
    SimplexAxisScratch& s2 = simplex_axis_scratch_slot(2);
    SimplexAxisScratch& s3 = simplex_axis_scratch_slot(3);
    s0.reserveFor(sequence_size);
    s1.reserveFor(sequence_size);
    s2.reserveFor(sequence_size);
    s3.reserveFor(sequence_size);

    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l3 = xi[2];
        const Real l0 = Real(1) - l1 - l2 - l3;

        Real* d0_out = (need_gradients || need_hessians) ? s0.dphi() : nullptr;
        Real* d1_out = (need_gradients || need_hessians) ? s1.dphi() : nullptr;
        Real* d2_out = (need_gradients || need_hessians) ? s2.dphi() : nullptr;
        Real* d3_out = (need_gradients || need_hessians) ? s3.dphi() : nullptr;
        Real* d20_out = need_hessians ? s0.d2phi() : nullptr;
        Real* d21_out = need_hessians ? s1.d2phi() : nullptr;
        Real* d22_out = need_hessians ? s2.d2phi() : nullptr;
        Real* d23_out = need_hessians ? s3.d2phi() : nullptr;

        simplex_lagrange_factor_sequence(order, l0, s0.phi(), d0_out, d20_out);
        simplex_lagrange_factor_sequence(order, l1, s1.phi(), d1_out, d21_out);
        simplex_lagrange_factor_sequence(order, l2, s2.phi(), d2_out, d22_out);
        simplex_lagrange_factor_sequence(order, l3, s3.phi(), d3_out, d23_out);
        const Real* phi0 = s0.phi();
        const Real* phi1 = s1.phi();
        const Real* phi2 = s2.phi();
        const Real* phi3 = s3.phi();
        const Real* dphi0 = s0.dphi();
        const Real* dphi1 = s1.dphi();
        const Real* dphi2 = s2.dphi();
        const Real* dphi3 = s3.dphi();
        const Real* d2phi0 = s0.d2phi();
        const Real* d2phi1 = s1.d2phi();
        const Real* d2phi2 = s2.d2phi();
        const Real* d2phi3 = s3.d2phi();

        for (std::size_t node = 0; node < num_nodes; ++node) {
            const auto& e = simplex_exponents[node];
            const std::size_t i0 = static_cast<std::size_t>(e[0]);
            const std::size_t i1 = static_cast<std::size_t>(e[1]);
            const std::size_t i2 = static_cast<std::size_t>(e[2]);
            const std::size_t i3 = static_cast<std::size_t>(e[3]);

            const Real v0 = phi0[i0];
            const Real v1 = phi1[i1];
            const Real v2 = phi2[i2];
            const Real v3 = phi3[i3];
            if (values_out != nullptr) {
                values_out[node * output_stride + q] = v0 * v1 * v2 * v3;
            }
            if (!need_gradients && !need_hessians) {
                continue;
            }

            const Real D0 = dphi0[i0];
            const Real D1 = dphi1[i1];
            const Real D2 = dphi2[i2];
            const Real D3 = dphi3[i3];

            if (gradients_out != nullptr) {
                const Real dl0 = D0 * v1 * v2 * v3;
                const Real dl1 = v0 * D1 * v2 * v3;
                const Real dl2 = v0 * v1 * D2 * v3;
                const Real dl3 = v0 * v1 * v2 * D3;
                Real* g = gradients_out + node * 3u * output_stride;
                g[0u * output_stride + q] = dl1 - dl0;
                g[1u * output_stride + q] = dl2 - dl0;
                g[2u * output_stride + q] = dl3 - dl0;
            }

            if (hessians_out != nullptr) {
                const Real DD0 = d2phi0[i0];
                const Real DD1 = d2phi1[i1];
                const Real DD2 = d2phi2[i2];
                const Real DD3 = d2phi3[i3];

                const Real H00 = DD0 * v1 * v2 * v3;
                const Real H11 = v0 * DD1 * v2 * v3;
                const Real H22 = v0 * v1 * DD2 * v3;
                const Real H33 = v0 * v1 * v2 * DD3;

                const Real H01 = D0 * D1 * v2 * v3;
                const Real H02 = D0 * v1 * D2 * v3;
                const Real H03 = D0 * v1 * v2 * D3;
                const Real H12 = v0 * D1 * D2 * v3;
                const Real H13 = v0 * D1 * v2 * D3;
                const Real H23 = v0 * v1 * D2 * D3;

                const Real h01 = H00 - H01 - H02 + H12;
                const Real h02 = H00 - H01 - H03 + H13;
                const Real h12 = H00 - H02 - H03 + H23;

                Real* H = hessians_out + node * 9u * output_stride;
                H[0u * output_stride + q] = H00 - Real(2) * H01 + H11;
                H[1u * output_stride + q] = h01;
                H[2u * output_stride + q] = h02;
                H[3u * output_stride + q] = h01;
                H[4u * output_stride + q] = H00 - Real(2) * H02 + H22;
                H[5u * output_stride + q] = h12;
                H[6u * output_stride + q] = h02;
                H[7u * output_stride + q] = h12;
                H[8u * output_stride + q] = H00 - Real(2) * H03 + H33;
            }
        }
    }
}

} // namespace detail
} // namespace basis
} // namespace FE
} // namespace svmp
