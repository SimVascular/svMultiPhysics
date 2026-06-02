/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#include "LagrangeBasis.h"
#include "BasisTraits.h"
#include "BasisTolerance.h"
#include "LagrangeBasisFast.h"
#include "NodeOrderingConventions.h"
#include "LagrangeBasisPyramid.h"
#include "LagrangeBasisSimplex.h"
#include "LagrangeBasisUtility.h"
#include <algorithm>
#include <cmath>
#include <unordered_map>

namespace svmp {
namespace FE {
namespace basis {

namespace {

using LagrangeTopology = BasisTopology;

#if defined(_MSC_VER)
#define SVMP_LAGRANGE_NOINLINE __declspec(noinline)
#define SVMP_LAGRANGE_ALIGN64
#elif defined(__GNUC__) || defined(__clang__)
#define SVMP_LAGRANGE_NOINLINE __attribute__((noinline))
#define SVMP_LAGRANGE_ALIGN64 __attribute__((aligned(64)))
#else
#define SVMP_LAGRANGE_NOINLINE
#define SVMP_LAGRANGE_ALIGN64
#endif

#ifndef FE_ALWAYS_INLINE
#if defined(_MSC_VER)
#define FE_ALWAYS_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define FE_ALWAYS_INLINE __attribute__((always_inline)) inline
#else
#define FE_ALWAYS_INLINE inline
#endif
#endif

SVMP_LAGRANGE_NOINLINE void evaluate_triangle_order1_gradients_strided(
    std::size_t num_qpts,
    std::size_t output_stride,
    Real* SVMP_RESTRICT gradients_out);

struct LagrangeTopologyTraits {
    LagrangeTopology topology;
    int dimension;
};

struct SimplexExponentHash {
    std::size_t operator()(const std::array<int, 4>& exponents) const noexcept {
        std::size_t seed = 0x9e3779b97f4a7c15ull;
        for (const int exponent : exponents) {
            const auto value = static_cast<std::size_t>(exponent);
            seed ^= value + 0x9e3779b97f4a7c15ull + (seed << 6u) + (seed >> 2u);
        }
        return seed;
    }
};

template<typename T, std::size_t N>
void assign_array(std::vector<T>& out, const std::array<T, N>& values) {
    out.assign(values.begin(), values.end());
}

bool coordinate_matches_expected(Real coord, Real expected) noexcept {
    return detail::basis_nearly_equal(coord, expected);
}

template<typename FastBasis>
void evaluate_fast_outputs(const math::Vector<Real, 3>& xi,
                           std::vector<Real>* values,
                           std::vector<Gradient>* gradients,
                           std::vector<Hessian>* hessians) {
    if (values != nullptr) {
        std::array<Real, FastBasis::n_dofs> fast_values{};
        FastBasis::evaluate(xi, fast_values);
        assign_array(*values, fast_values);
    }
    if (gradients != nullptr) {
        std::array<Gradient, FastBasis::n_dofs> fast_gradients{};
        FastBasis::evaluate_gradients(xi, fast_gradients);
        assign_array(*gradients, fast_gradients);
    }
    if (hessians != nullptr) {
        std::array<Hessian, FastBasis::n_dofs> fast_hessians{};
        FastBasis::evaluate_hessians(xi, fast_hessians);
        assign_array(*hessians, fast_hessians);
    }
}

template<typename FastBasis>
void evaluate_fast_outputs_to(const math::Vector<Real, 3>& xi,
                              Real* SVMP_RESTRICT values_out,
                              Real* SVMP_RESTRICT gradients_out,
                              Real* SVMP_RESTRICT hessians_out) {
    if (values_out != nullptr) {
        std::array<Real, FastBasis::n_dofs> fast_values{};
        FastBasis::evaluate(xi, fast_values);
        for (std::size_t i = 0; i < fast_values.size(); ++i) {
            values_out[i] = fast_values[i];
        }
    }
    if (gradients_out != nullptr) {
        std::array<Gradient, FastBasis::n_dofs> fast_gradients{};
        FastBasis::evaluate_gradients(xi, fast_gradients);
        for (std::size_t i = 0; i < fast_gradients.size(); ++i) {
            gradients_out[i * 3u + 0u] = fast_gradients[i][0];
            gradients_out[i * 3u + 1u] = fast_gradients[i][1];
            gradients_out[i * 3u + 2u] = fast_gradients[i][2];
        }
    }
    if (hessians_out != nullptr) {
        std::array<Hessian, FastBasis::n_dofs> fast_hessians{};
        FastBasis::evaluate_hessians(xi, fast_hessians);
        for (std::size_t i = 0; i < fast_hessians.size(); ++i) {
            store_hessian(fast_hessians[i], hessians_out + i * 9u);
        }
    }
}

template<typename FastBasis>
void evaluate_fast_outputs_strided(const std::vector<math::Vector<Real, 3>>& points,
                                   std::size_t output_stride,
                                   Real* SVMP_RESTRICT values_out,
                                   Real* SVMP_RESTRICT gradients_out,
                                   Real* SVMP_RESTRICT hessians_out) {
    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        if (values_out != nullptr) {
            std::array<Real, FastBasis::n_dofs> fast_values{};
            FastBasis::evaluate(xi, fast_values);
            for (std::size_t i = 0; i < fast_values.size(); ++i) {
                values_out[i * output_stride + q] = fast_values[i];
            }
        }
        if (gradients_out != nullptr) {
            std::array<Gradient, FastBasis::n_dofs> fast_gradients{};
            FastBasis::evaluate_gradients(xi, fast_gradients);
            for (std::size_t i = 0; i < fast_gradients.size(); ++i) {
                Real* g = gradients_out + i * 3u * output_stride;
                g[0u * output_stride + q] = fast_gradients[i][0];
                g[1u * output_stride + q] = fast_gradients[i][1];
                g[2u * output_stride + q] = fast_gradients[i][2];
            }
        }
        if (hessians_out != nullptr) {
            std::array<Hessian, FastBasis::n_dofs> fast_hessians{};
            FastBasis::evaluate_hessians(xi, fast_hessians);
            for (std::size_t i = 0; i < fast_hessians.size(); ++i) {
                const Hessian& hessian = fast_hessians[i];
                Real* H = hessians_out + i * 9u * output_stride;
                H[0u * output_stride + q] = hessian(0, 0);
                H[1u * output_stride + q] = hessian(0, 1);
                H[2u * output_stride + q] = hessian(0, 2);
                H[3u * output_stride + q] = hessian(1, 0);
                H[4u * output_stride + q] = hessian(1, 1);
                H[5u * output_stride + q] = hessian(1, 2);
                H[6u * output_stride + q] = hessian(2, 0);
                H[7u * output_stride + q] = hessian(2, 1);
                H[8u * output_stride + q] = hessian(2, 2);
            }
        }
    }
}

template<int Order>
bool evaluate_fixed_lagrange_fast_order(LagrangeTopology topology,
                                        const math::Vector<Real, 3>& xi,
                                        std::vector<Real>* values,
                                        std::vector<Gradient>* gradients,
                                        std::vector<Hessian>* hessians) {
    switch (topology) {
        case LagrangeTopology::Line:
            evaluate_fast_outputs<LagrangeLineFast<Order>>(xi, values, gradients, hessians);
            return true;
        case LagrangeTopology::Quadrilateral:
            evaluate_fast_outputs<LagrangeQuadFast<Order>>(xi, values, gradients, hessians);
            return true;
        case LagrangeTopology::Hexahedron:
            evaluate_fast_outputs<LagrangeHexFast<Order>>(xi, values, gradients, hessians);
            return true;
        case LagrangeTopology::Triangle:
            evaluate_fast_outputs<LagrangeTriFast<Order>>(xi, values, gradients, hessians);
            return true;
        case LagrangeTopology::Tetrahedron:
            evaluate_fast_outputs<LagrangeTetFast<Order>>(xi, values, gradients, hessians);
            return true;
        default:
            return false;
    }
}

template<int Order>
bool evaluate_fixed_lagrange_fast_to_order(LagrangeTopology topology,
                                           const math::Vector<Real, 3>& xi,
                                           Real* SVMP_RESTRICT values_out,
                                           Real* SVMP_RESTRICT gradients_out,
                                           Real* SVMP_RESTRICT hessians_out) {
    switch (topology) {
        case LagrangeTopology::Line:
            evaluate_fast_outputs_to<LagrangeLineFast<Order>>(xi, values_out, gradients_out, hessians_out);
            return true;
        case LagrangeTopology::Quadrilateral:
            evaluate_fast_outputs_to<LagrangeQuadFast<Order>>(xi, values_out, gradients_out, hessians_out);
            return true;
        case LagrangeTopology::Hexahedron:
            evaluate_fast_outputs_to<LagrangeHexFast<Order>>(xi, values_out, gradients_out, hessians_out);
            return true;
        case LagrangeTopology::Triangle:
            evaluate_fast_outputs_to<LagrangeTriFast<Order>>(xi, values_out, gradients_out, hessians_out);
            return true;
        case LagrangeTopology::Tetrahedron:
            evaluate_fast_outputs_to<LagrangeTetFast<Order>>(xi, values_out, gradients_out, hessians_out);
            return true;
        default:
            return false;
    }
}

template<int Order>
bool evaluate_fixed_lagrange_fast_strided_order(
    LagrangeTopology topology,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    switch (topology) {
        case LagrangeTopology::Line:
            evaluate_fast_outputs_strided<LagrangeLineFast<Order>>(
                points, output_stride, values_out, gradients_out, hessians_out);
            return true;
        case LagrangeTopology::Quadrilateral:
            evaluate_fast_outputs_strided<LagrangeQuadFast<Order>>(
                points, output_stride, values_out, gradients_out, hessians_out);
            return true;
        case LagrangeTopology::Hexahedron:
            evaluate_fast_outputs_strided<LagrangeHexFast<Order>>(
                points, output_stride, values_out, gradients_out, hessians_out);
            return true;
        case LagrangeTopology::Triangle:
            evaluate_fast_outputs_strided<LagrangeTriFast<Order>>(
                points, output_stride, values_out, gradients_out, hessians_out);
            return true;
        case LagrangeTopology::Tetrahedron:
            evaluate_fast_outputs_strided<LagrangeTetFast<Order>>(
                points, output_stride, values_out, gradients_out, hessians_out);
            return true;
        default:
            return false;
    }
}

void evaluate_triangle_order3_values_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    Real* row0 = values_out + 0u * output_stride;
    Real* row1 = values_out + 1u * output_stride;
    Real* row2 = values_out + 2u * output_stride;
    Real* row3 = values_out + 3u * output_stride;
    Real* row4 = values_out + 4u * output_stride;
    Real* row5 = values_out + 5u * output_stride;
    Real* row6 = values_out + 6u * output_stride;
    Real* row7 = values_out + 7u * output_stride;
    Real* row8 = values_out + 8u * output_stride;
    Real* row9 = values_out + 9u * output_stride;

    if (points.size() == 4u && output_stride == 4u) {
        Real p10[4];
        Real p11[4];
        Real p12[4];
        Real p20[4];
        Real p21[4];
        Real p22[4];
        Real p30[4];
        Real p31[4];
        Real p32[4];

        for (std::size_t q = 0; q < 4u; ++q) {
            const auto& xi = points[q];
            const Real l1 = xi[0];
            const Real l2 = xi[1];
            const Real l0 = Real(1) - l1 - l2;

            p10[q] = Real(3) * l0;
            p11[q] = Real(3) * l1;
            p12[q] = Real(3) * l2;
            p20[q] = Real(0.5) * p10[q] * (p10[q] - Real(1));
            p21[q] = Real(0.5) * p11[q] * (p11[q] - Real(1));
            p22[q] = Real(0.5) * p12[q] * (p12[q] - Real(1));
            p30[q] = (p10[q] * (p10[q] - Real(1)) * (p10[q] - Real(2))) / Real(6);
            p31[q] = (p11[q] * (p11[q] - Real(1)) * (p11[q] - Real(2))) / Real(6);
            p32[q] = (p12[q] * (p12[q] - Real(1)) * (p12[q] - Real(2))) / Real(6);
        }

        row0[0] = p30[0]; row0[1] = p30[1]; row0[2] = p30[2]; row0[3] = p30[3];
        row1[0] = p31[0]; row1[1] = p31[1]; row1[2] = p31[2]; row1[3] = p31[3];
        row2[0] = p32[0]; row2[1] = p32[1]; row2[2] = p32[2]; row2[3] = p32[3];
        row3[0] = p20[0] * p11[0];
        row3[1] = p20[1] * p11[1];
        row3[2] = p20[2] * p11[2];
        row3[3] = p20[3] * p11[3];
        row4[0] = p10[0] * p21[0];
        row4[1] = p10[1] * p21[1];
        row4[2] = p10[2] * p21[2];
        row4[3] = p10[3] * p21[3];
        row5[0] = p21[0] * p12[0];
        row5[1] = p21[1] * p12[1];
        row5[2] = p21[2] * p12[2];
        row5[3] = p21[3] * p12[3];
        row6[0] = p11[0] * p22[0];
        row6[1] = p11[1] * p22[1];
        row6[2] = p11[2] * p22[2];
        row6[3] = p11[3] * p22[3];
        row7[0] = p10[0] * p22[0];
        row7[1] = p10[1] * p22[1];
        row7[2] = p10[2] * p22[2];
        row7[3] = p10[3] * p22[3];
        row8[0] = p20[0] * p12[0];
        row8[1] = p20[1] * p12[1];
        row8[2] = p20[2] * p12[2];
        row8[3] = p20[3] * p12[3];
        row9[0] = p10[0] * p11[0] * p12[0];
        row9[1] = p10[1] * p11[1] * p12[1];
        row9[2] = p10[2] * p11[2] * p12[2];
        row9[3] = p10[3] * p11[3] * p12[3];
        return;
    }

    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l0 = Real(1) - l1 - l2;

        const Real p10 = Real(3) * l0;
        const Real p11 = Real(3) * l1;
        const Real p12 = Real(3) * l2;
        const Real p20 = Real(0.5) * p10 * (p10 - Real(1));
        const Real p21 = Real(0.5) * p11 * (p11 - Real(1));
        const Real p22 = Real(0.5) * p12 * (p12 - Real(1));
        const Real p30 = (p10 * (p10 - Real(1)) * (p10 - Real(2))) / Real(6);
        const Real p31 = (p11 * (p11 - Real(1)) * (p11 - Real(2))) / Real(6);
        const Real p32 = (p12 * (p12 - Real(1)) * (p12 - Real(2))) / Real(6);

        row0[q] = p30;
        row1[q] = p31;
        row2[q] = p32;
        row3[q] = p20 * p11;
        row4[q] = p10 * p21;
        row5[q] = p21 * p12;
        row6[q] = p11 * p22;
        row7[q] = p10 * p22;
        row8[q] = p20 * p12;
        row9[q] = p10 * p11 * p12;
    }
}

void evaluate_triangle_order2_values_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    Real* row0 = values_out + 0u * output_stride;
    Real* row1 = values_out + 1u * output_stride;
    Real* row2 = values_out + 2u * output_stride;
    Real* row3 = values_out + 3u * output_stride;
    Real* row4 = values_out + 4u * output_stride;
    Real* row5 = values_out + 5u * output_stride;

    if (points.size() == 4u && output_stride == 4u) {
        Real l0[4];
        Real l1[4];
        Real l2[4];
        for (std::size_t q = 0; q < 4u; ++q) {
            const auto& xi = points[q];
            l1[q] = xi[0];
            l2[q] = xi[1];
            l0[q] = Real(1) - l1[q] - l2[q];
        }

        row0[0] = l0[0] * (Real(2) * l0[0] - Real(1));
        row0[1] = l0[1] * (Real(2) * l0[1] - Real(1));
        row0[2] = l0[2] * (Real(2) * l0[2] - Real(1));
        row0[3] = l0[3] * (Real(2) * l0[3] - Real(1));
        row1[0] = l1[0] * (Real(2) * l1[0] - Real(1));
        row1[1] = l1[1] * (Real(2) * l1[1] - Real(1));
        row1[2] = l1[2] * (Real(2) * l1[2] - Real(1));
        row1[3] = l1[3] * (Real(2) * l1[3] - Real(1));
        row2[0] = l2[0] * (Real(2) * l2[0] - Real(1));
        row2[1] = l2[1] * (Real(2) * l2[1] - Real(1));
        row2[2] = l2[2] * (Real(2) * l2[2] - Real(1));
        row2[3] = l2[3] * (Real(2) * l2[3] - Real(1));
        row3[0] = Real(4) * l0[0] * l1[0];
        row3[1] = Real(4) * l0[1] * l1[1];
        row3[2] = Real(4) * l0[2] * l1[2];
        row3[3] = Real(4) * l0[3] * l1[3];
        row4[0] = Real(4) * l1[0] * l2[0];
        row4[1] = Real(4) * l1[1] * l2[1];
        row4[2] = Real(4) * l1[2] * l2[2];
        row4[3] = Real(4) * l1[3] * l2[3];
        row5[0] = Real(4) * l0[0] * l2[0];
        row5[1] = Real(4) * l0[1] * l2[1];
        row5[2] = Real(4) * l0[2] * l2[2];
        row5[3] = Real(4) * l0[3] * l2[3];
        return;
    }

    auto write_q = [&](std::size_t q) {
        const auto& xi = points[q];
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l0 = Real(1) - l1 - l2;
        row0[q] = l0 * (Real(2) * l0 - Real(1));
        row1[q] = l1 * (Real(2) * l1 - Real(1));
        row2[q] = l2 * (Real(2) * l2 - Real(1));
        row3[q] = Real(4) * l0 * l1;
        row4[q] = Real(4) * l1 * l2;
        row5[q] = Real(4) * l0 * l2;
    };

    if (points.size() == 4u) {
        write_q(0u);
        write_q(1u);
        write_q(2u);
        write_q(3u);
        return;
    }

    for (std::size_t q = 0; q < points.size(); ++q) {
        write_q(q);
    }
}

void evaluate_triangle_order1_values_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    Real* row0 = values_out + 0u * output_stride;
    Real* row1 = values_out + 1u * output_stride;
    Real* row2 = values_out + 2u * output_stride;

    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        row0[q] = Real(1) - xi[0] - xi[1];
        row1[q] = xi[0];
        row2[q] = xi[1];
    }
}

void evaluate_triangle_order2_gradients_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT gradients_out) {
    Real* row0 = gradients_out + 0u * 3u * output_stride;
    Real* row1 = gradients_out + 1u * 3u * output_stride;
    Real* row2 = gradients_out + 2u * 3u * output_stride;
    Real* row3 = gradients_out + 3u * 3u * output_stride;
    Real* row4 = gradients_out + 4u * 3u * output_stride;
    Real* row5 = gradients_out + 5u * 3u * output_stride;

    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l0 = Real(1) - l1 - l2;
        const Real g0 = Real(1) - Real(4) * l0;
        row0[0u * output_stride + q] = g0;
        row0[1u * output_stride + q] = g0;
        row0[2u * output_stride + q] = Real(0);
        row1[0u * output_stride + q] = Real(4) * l1 - Real(1);
        row1[1u * output_stride + q] = Real(0);
        row1[2u * output_stride + q] = Real(0);
        row2[0u * output_stride + q] = Real(0);
        row2[1u * output_stride + q] = Real(4) * l2 - Real(1);
        row2[2u * output_stride + q] = Real(0);
        row3[0u * output_stride + q] = Real(4) * (l0 - l1);
        row3[1u * output_stride + q] = Real(-4) * l1;
        row3[2u * output_stride + q] = Real(0);
        row4[0u * output_stride + q] = Real(4) * l2;
        row4[1u * output_stride + q] = Real(4) * l1;
        row4[2u * output_stride + q] = Real(0);
        row5[0u * output_stride + q] = Real(-4) * l2;
        row5[1u * output_stride + q] = Real(4) * (l0 - l2);
        row5[2u * output_stride + q] = Real(0);
    }
}

inline void write_constant_hessian_q4(Real* SVMP_RESTRICT row,
                                      std::size_t output_stride,
                                      Real h00,
                                      Real h01,
                                      Real h02,
                                      Real h10,
                                      Real h11,
                                      Real h12,
                                      Real h20,
                                      Real h21,
                                      Real h22) {
    Real* c0 = row + 0u * output_stride;
    Real* c1 = row + 1u * output_stride;
    Real* c2 = row + 2u * output_stride;
    Real* c3 = row + 3u * output_stride;
    Real* c4 = row + 4u * output_stride;
    Real* c5 = row + 5u * output_stride;
    Real* c6 = row + 6u * output_stride;
    Real* c7 = row + 7u * output_stride;
    Real* c8 = row + 8u * output_stride;

    c0[0] = h00; c0[1] = h00; c0[2] = h00; c0[3] = h00;
    c1[0] = h01; c1[1] = h01; c1[2] = h01; c1[3] = h01;
    c2[0] = h02; c2[1] = h02; c2[2] = h02; c2[3] = h02;
    c3[0] = h10; c3[1] = h10; c3[2] = h10; c3[3] = h10;
    c4[0] = h11; c4[1] = h11; c4[2] = h11; c4[3] = h11;
    c5[0] = h12; c5[1] = h12; c5[2] = h12; c5[3] = h12;
    c6[0] = h20; c6[1] = h20; c6[2] = h20; c6[3] = h20;
    c7[0] = h21; c7[1] = h21; c7[2] = h21; c7[3] = h21;
    c8[0] = h22; c8[1] = h22; c8[2] = h22; c8[3] = h22;
}

SVMP_LAGRANGE_NOINLINE void evaluate_triangle_order2_hessians_q4(
    std::size_t output_stride,
    Real* SVMP_RESTRICT hessians_out) {
    write_constant_hessian_q4(hessians_out + 0u * 9u * output_stride,
                              output_stride,
                              Real(4), Real(4), Real(0),
                              Real(4), Real(4), Real(0),
                              Real(0), Real(0), Real(0));
    write_constant_hessian_q4(hessians_out + 1u * 9u * output_stride,
                              output_stride,
                              Real(4), Real(0), Real(0),
                              Real(0), Real(0), Real(0),
                              Real(0), Real(0), Real(0));
    write_constant_hessian_q4(hessians_out + 2u * 9u * output_stride,
                              output_stride,
                              Real(0), Real(0), Real(0),
                              Real(0), Real(4), Real(0),
                              Real(0), Real(0), Real(0));
    write_constant_hessian_q4(hessians_out + 3u * 9u * output_stride,
                              output_stride,
                              Real(-8), Real(-4), Real(0),
                              Real(-4), Real(0), Real(0),
                              Real(0), Real(0), Real(0));
    write_constant_hessian_q4(hessians_out + 4u * 9u * output_stride,
                              output_stride,
                              Real(0), Real(4), Real(0),
                              Real(4), Real(0), Real(0),
                              Real(0), Real(0), Real(0));
    write_constant_hessian_q4(hessians_out + 5u * 9u * output_stride,
                              output_stride,
                              Real(0), Real(-4), Real(0),
                              Real(-4), Real(-8), Real(0),
                              Real(0), Real(0), Real(0));
}

SVMP_LAGRANGE_NOINLINE void evaluate_tet_order2_hessians_q4(
    std::size_t output_stride,
    Real* SVMP_RESTRICT hessians_out) {
    write_constant_hessian_q4(hessians_out + 0u * 9u * output_stride,
                              output_stride,
                              Real(4), Real(4), Real(4),
                              Real(4), Real(4), Real(4),
                              Real(4), Real(4), Real(4));
    write_constant_hessian_q4(hessians_out + 1u * 9u * output_stride,
                              output_stride,
                              Real(4), Real(0), Real(0),
                              Real(0), Real(0), Real(0),
                              Real(0), Real(0), Real(0));
    write_constant_hessian_q4(hessians_out + 2u * 9u * output_stride,
                              output_stride,
                              Real(0), Real(0), Real(0),
                              Real(0), Real(4), Real(0),
                              Real(0), Real(0), Real(0));
    write_constant_hessian_q4(hessians_out + 3u * 9u * output_stride,
                              output_stride,
                              Real(0), Real(0), Real(0),
                              Real(0), Real(0), Real(0),
                              Real(0), Real(0), Real(4));
    write_constant_hessian_q4(hessians_out + 4u * 9u * output_stride,
                              output_stride,
                              Real(-8), Real(-4), Real(-4),
                              Real(-4), Real(0), Real(0),
                              Real(-4), Real(0), Real(0));
    write_constant_hessian_q4(hessians_out + 5u * 9u * output_stride,
                              output_stride,
                              Real(0), Real(4), Real(0),
                              Real(4), Real(0), Real(0),
                              Real(0), Real(0), Real(0));
    write_constant_hessian_q4(hessians_out + 6u * 9u * output_stride,
                              output_stride,
                              Real(0), Real(-4), Real(0),
                              Real(-4), Real(-8), Real(-4),
                              Real(0), Real(-4), Real(0));
    write_constant_hessian_q4(hessians_out + 7u * 9u * output_stride,
                              output_stride,
                              Real(0), Real(0), Real(-4),
                              Real(0), Real(0), Real(-4),
                              Real(-4), Real(-4), Real(-8));
    write_constant_hessian_q4(hessians_out + 8u * 9u * output_stride,
                              output_stride,
                              Real(0), Real(0), Real(4),
                              Real(0), Real(0), Real(0),
                              Real(4), Real(0), Real(0));
    write_constant_hessian_q4(hessians_out + 9u * 9u * output_stride,
                              output_stride,
                              Real(0), Real(0), Real(0),
                              Real(0), Real(0), Real(4),
                              Real(0), Real(4), Real(0));
}

void evaluate_tet_order1_values_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    Real* row0 = values_out + 0u * output_stride;
    Real* row1 = values_out + 1u * output_stride;
    Real* row2 = values_out + 2u * output_stride;
    Real* row3 = values_out + 3u * output_stride;

    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        row0[q] = Real(1) - xi[0] - xi[1] - xi[2];
        row1[q] = xi[0];
        row2[q] = xi[1];
        row3[q] = xi[2];
    }
}

void evaluate_tet_order1_gradients_strided(
    std::size_t num_qpts,
    std::size_t output_stride,
    Real* SVMP_RESTRICT gradients_out) {
    Real* row0 = gradients_out + 0u * 3u * output_stride;
    Real* row1 = gradients_out + 1u * 3u * output_stride;
    Real* row2 = gradients_out + 2u * 3u * output_stride;
    Real* row3 = gradients_out + 3u * 3u * output_stride;

    for (std::size_t q = 0; q < num_qpts; ++q) {
        row0[0u * output_stride + q] = Real(-1);
        row0[1u * output_stride + q] = Real(-1);
        row0[2u * output_stride + q] = Real(-1);
        row1[0u * output_stride + q] = Real(1);
        row1[1u * output_stride + q] = Real(0);
        row1[2u * output_stride + q] = Real(0);
        row2[0u * output_stride + q] = Real(0);
        row2[1u * output_stride + q] = Real(1);
        row2[2u * output_stride + q] = Real(0);
        row3[0u * output_stride + q] = Real(0);
        row3[1u * output_stride + q] = Real(0);
        row3[2u * output_stride + q] = Real(1);
    }
}

void evaluate_zero_hessians_strided(
    std::size_t num_nodes,
    std::size_t num_qpts,
    std::size_t output_stride,
    Real* SVMP_RESTRICT hessians_out) {
    if (num_qpts == 4u) {
        for (std::size_t node = 0; node < num_nodes; ++node) {
            write_constant_hessian_q4(hessians_out + node * 9u * output_stride,
                                      output_stride,
                                      Real(0), Real(0), Real(0),
                                      Real(0), Real(0), Real(0),
                                      Real(0), Real(0), Real(0));
        }
        return;
    }

    for (std::size_t node = 0; node < num_nodes; ++node) {
        Real* row = hessians_out + node * 9u * output_stride;
        for (std::size_t q = 0; q < num_qpts; ++q) {
            row[0u * output_stride + q] = Real(0);
            row[1u * output_stride + q] = Real(0);
            row[2u * output_stride + q] = Real(0);
            row[3u * output_stride + q] = Real(0);
            row[4u * output_stride + q] = Real(0);
            row[5u * output_stride + q] = Real(0);
            row[6u * output_stride + q] = Real(0);
            row[7u * output_stride + q] = Real(0);
            row[8u * output_stride + q] = Real(0);
        }
    }
}

void evaluate_tet_order2_values_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    Real* row0 = values_out + 0u * output_stride;
    Real* row1 = values_out + 1u * output_stride;
    Real* row2 = values_out + 2u * output_stride;
    Real* row3 = values_out + 3u * output_stride;
    Real* row4 = values_out + 4u * output_stride;
    Real* row5 = values_out + 5u * output_stride;
    Real* row6 = values_out + 6u * output_stride;
    Real* row7 = values_out + 7u * output_stride;
    Real* row8 = values_out + 8u * output_stride;
    Real* row9 = values_out + 9u * output_stride;

    if (points.size() == 4u && output_stride == 4u) {
        Real l0[4];
        Real l1[4];
        Real l2[4];
        Real l3[4];
        for (std::size_t q = 0; q < 4u; ++q) {
            const auto& xi = points[q];
            l1[q] = xi[0];
            l2[q] = xi[1];
            l3[q] = xi[2];
            l0[q] = Real(1) - l1[q] - l2[q] - l3[q];
        }

        row0[0] = l0[0] * (Real(2) * l0[0] - Real(1));
        row0[1] = l0[1] * (Real(2) * l0[1] - Real(1));
        row0[2] = l0[2] * (Real(2) * l0[2] - Real(1));
        row0[3] = l0[3] * (Real(2) * l0[3] - Real(1));
        row1[0] = l1[0] * (Real(2) * l1[0] - Real(1));
        row1[1] = l1[1] * (Real(2) * l1[1] - Real(1));
        row1[2] = l1[2] * (Real(2) * l1[2] - Real(1));
        row1[3] = l1[3] * (Real(2) * l1[3] - Real(1));
        row2[0] = l2[0] * (Real(2) * l2[0] - Real(1));
        row2[1] = l2[1] * (Real(2) * l2[1] - Real(1));
        row2[2] = l2[2] * (Real(2) * l2[2] - Real(1));
        row2[3] = l2[3] * (Real(2) * l2[3] - Real(1));
        row3[0] = l3[0] * (Real(2) * l3[0] - Real(1));
        row3[1] = l3[1] * (Real(2) * l3[1] - Real(1));
        row3[2] = l3[2] * (Real(2) * l3[2] - Real(1));
        row3[3] = l3[3] * (Real(2) * l3[3] - Real(1));
        row4[0] = Real(4) * l0[0] * l1[0];
        row4[1] = Real(4) * l0[1] * l1[1];
        row4[2] = Real(4) * l0[2] * l1[2];
        row4[3] = Real(4) * l0[3] * l1[3];
        row5[0] = Real(4) * l1[0] * l2[0];
        row5[1] = Real(4) * l1[1] * l2[1];
        row5[2] = Real(4) * l1[2] * l2[2];
        row5[3] = Real(4) * l1[3] * l2[3];
        row6[0] = Real(4) * l0[0] * l2[0];
        row6[1] = Real(4) * l0[1] * l2[1];
        row6[2] = Real(4) * l0[2] * l2[2];
        row6[3] = Real(4) * l0[3] * l2[3];
        row7[0] = Real(4) * l0[0] * l3[0];
        row7[1] = Real(4) * l0[1] * l3[1];
        row7[2] = Real(4) * l0[2] * l3[2];
        row7[3] = Real(4) * l0[3] * l3[3];
        row8[0] = Real(4) * l1[0] * l3[0];
        row8[1] = Real(4) * l1[1] * l3[1];
        row8[2] = Real(4) * l1[2] * l3[2];
        row8[3] = Real(4) * l1[3] * l3[3];
        row9[0] = Real(4) * l2[0] * l3[0];
        row9[1] = Real(4) * l2[1] * l3[1];
        row9[2] = Real(4) * l2[2] * l3[2];
        row9[3] = Real(4) * l2[3] * l3[3];
        return;
    }

    auto write_q = [&](std::size_t q) {
        const auto& xi = points[q];
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l3 = xi[2];
        const Real l0 = Real(1) - l1 - l2 - l3;
        row0[q] = l0 * (Real(2) * l0 - Real(1));
        row1[q] = l1 * (Real(2) * l1 - Real(1));
        row2[q] = l2 * (Real(2) * l2 - Real(1));
        row3[q] = l3 * (Real(2) * l3 - Real(1));
        row4[q] = Real(4) * l0 * l1;
        row5[q] = Real(4) * l1 * l2;
        row6[q] = Real(4) * l0 * l2;
        row7[q] = Real(4) * l0 * l3;
        row8[q] = Real(4) * l1 * l3;
        row9[q] = Real(4) * l2 * l3;
    };

    if (points.size() == 4u) {
        write_q(0u);
        write_q(1u);
        write_q(2u);
        write_q(3u);
        return;
    }

    for (std::size_t q = 0; q < points.size(); ++q) {
        write_q(q);
    }
}

inline void write_tet_order2_gradient_q(Real* SVMP_RESTRICT row,
                                        std::size_t output_stride,
                                        std::size_t q,
                                        Real gx,
                                        Real gy,
                                        Real gz) {
    row[0u * output_stride + q] = gx;
    row[1u * output_stride + q] = gy;
    row[2u * output_stride + q] = gz;
}

SVMP_LAGRANGE_NOINLINE void evaluate_tet_order2_gradients_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT gradients_out) {
    Real* row0 = gradients_out + 0u * 3u * output_stride;
    Real* row1 = gradients_out + 1u * 3u * output_stride;
    Real* row2 = gradients_out + 2u * 3u * output_stride;
    Real* row3 = gradients_out + 3u * 3u * output_stride;
    Real* row4 = gradients_out + 4u * 3u * output_stride;
    Real* row5 = gradients_out + 5u * 3u * output_stride;
    Real* row6 = gradients_out + 6u * 3u * output_stride;
    Real* row7 = gradients_out + 7u * 3u * output_stride;
    Real* row8 = gradients_out + 8u * 3u * output_stride;
    Real* row9 = gradients_out + 9u * 3u * output_stride;

    auto write_q = [&](std::size_t q) {
        const auto& xi = points[q];
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l3 = xi[2];
        const Real l0 = Real(1) - l1 - l2 - l3;
        const Real four = Real(4);
        const Real g0 = Real(1) - four * l0;

        write_tet_order2_gradient_q(row0, output_stride, q, g0, g0, g0);
        write_tet_order2_gradient_q(row1, output_stride, q, four * l1 - Real(1), Real(0), Real(0));
        write_tet_order2_gradient_q(row2, output_stride, q, Real(0), four * l2 - Real(1), Real(0));
        write_tet_order2_gradient_q(row3, output_stride, q, Real(0), Real(0), four * l3 - Real(1));
        write_tet_order2_gradient_q(row4, output_stride, q, four * (l0 - l1), -four * l1, -four * l1);
        write_tet_order2_gradient_q(row5, output_stride, q, four * l2, four * l1, Real(0));
        write_tet_order2_gradient_q(row6, output_stride, q, -four * l2, four * (l0 - l2), -four * l2);
        write_tet_order2_gradient_q(row7, output_stride, q, -four * l3, -four * l3, four * (l0 - l3));
        write_tet_order2_gradient_q(row8, output_stride, q, four * l3, Real(0), four * l1);
        write_tet_order2_gradient_q(row9, output_stride, q, Real(0), four * l3, four * l2);
    };

    if (points.size() == 4u) {
        write_q(0u);
        write_q(1u);
        write_q(2u);
        write_q(3u);
        return;
    }

    for (std::size_t q = 0; q < points.size(); ++q) {
        write_q(q);
    }
}

inline void fill_simplex_order3_factor_values(Real lambda, Real* SVMP_RESTRICT phi) {
    const Real t = Real(3) * lambda;
    phi[0] = Real(1);
    phi[1] = t;
    phi[2] = phi[1] * (t - Real(1)) * Real(0.5);
    phi[3] = phi[2] * (t - Real(2)) / Real(3);
}

void evaluate_tet_order3_values_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    Real* row0 = values_out + 0u * output_stride;
    Real* row1 = values_out + 1u * output_stride;
    Real* row2 = values_out + 2u * output_stride;
    Real* row3 = values_out + 3u * output_stride;
    Real* row4 = values_out + 4u * output_stride;
    Real* row5 = values_out + 5u * output_stride;
    Real* row6 = values_out + 6u * output_stride;
    Real* row7 = values_out + 7u * output_stride;
    Real* row8 = values_out + 8u * output_stride;
    Real* row9 = values_out + 9u * output_stride;
    Real* row10 = values_out + 10u * output_stride;
    Real* row11 = values_out + 11u * output_stride;
    Real* row12 = values_out + 12u * output_stride;
    Real* row13 = values_out + 13u * output_stride;
    Real* row14 = values_out + 14u * output_stride;
    Real* row15 = values_out + 15u * output_stride;
    Real* row16 = values_out + 16u * output_stride;
    Real* row17 = values_out + 17u * output_stride;
    Real* row18 = values_out + 18u * output_stride;
    Real* row19 = values_out + 19u * output_stride;

    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l3 = xi[2];
        const Real l0 = Real(1) - l1 - l2 - l3;
        Real p0[4];
        Real p1[4];
        Real p2[4];
        Real p3[4];
        fill_simplex_order3_factor_values(l0, p0);
        fill_simplex_order3_factor_values(l1, p1);
        fill_simplex_order3_factor_values(l2, p2);
        fill_simplex_order3_factor_values(l3, p3);

        row0[q] = p0[3];
        row1[q] = p1[3];
        row2[q] = p2[3];
        row3[q] = p3[3];
        row4[q] = p0[2] * p1[1];
        row5[q] = p0[1] * p1[2];
        row6[q] = p1[2] * p2[1];
        row7[q] = p1[1] * p2[2];
        row8[q] = p0[1] * p2[2];
        row9[q] = p0[2] * p2[1];
        row10[q] = p0[2] * p3[1];
        row11[q] = p0[1] * p3[2];
        row12[q] = p1[2] * p3[1];
        row13[q] = p1[1] * p3[2];
        row14[q] = p2[2] * p3[1];
        row15[q] = p2[1] * p3[2];
        row16[q] = p0[1] * p1[1] * p2[1];
        row17[q] = p0[1] * p1[1] * p3[1];
        row18[q] = p1[1] * p2[1] * p3[1];
        row19[q] = p0[1] * p2[1] * p3[1];
    }
}

void evaluate_triangle_order3_gradients_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT gradients_out) {
    Real* rows[10] = {
        gradients_out + 0u * 3u * output_stride,
        gradients_out + 1u * 3u * output_stride,
        gradients_out + 2u * 3u * output_stride,
        gradients_out + 3u * 3u * output_stride,
        gradients_out + 4u * 3u * output_stride,
        gradients_out + 5u * 3u * output_stride,
        gradients_out + 6u * 3u * output_stride,
        gradients_out + 7u * 3u * output_stride,
        gradients_out + 8u * 3u * output_stride,
        gradients_out + 9u * 3u * output_stride,
    };

    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l0 = Real(1) - l1 - l2;

        const Real p10 = Real(3) * l0;
        const Real p11 = Real(3) * l1;
        const Real p12 = Real(3) * l2;
        const Real p20 = Real(0.5) * p10 * (p10 - Real(1));
        const Real p21 = Real(0.5) * p11 * (p11 - Real(1));
        const Real p22 = Real(0.5) * p12 * (p12 - Real(1));
        const Real d10 = Real(3);
        const Real d11 = Real(3);
        const Real d12 = Real(3);
        const Real d20 = Real(3) * p10 - Real(1.5);
        const Real d21 = Real(3) * p11 - Real(1.5);
        const Real d22 = Real(3) * p12 - Real(1.5);
        const Real d30 = Real(1.5) * p10 * p10 - Real(3) * p10 + Real(1);
        const Real d31 = Real(1.5) * p11 * p11 - Real(3) * p11 + Real(1);
        const Real d32 = Real(1.5) * p12 * p12 - Real(3) * p12 + Real(1);

        const Real dl0[10] = {
            d30,
            Real(0),
            Real(0),
            d20 * p11,
            d10 * p21,
            Real(0),
            Real(0),
            d10 * p22,
            d20 * p12,
            d10 * p11 * p12,
        };
        const Real dl1[10] = {
            Real(0),
            d31,
            Real(0),
            p20 * d11,
            p10 * d21,
            d21 * p12,
            d11 * p22,
            Real(0),
            Real(0),
            p10 * d11 * p12,
        };
        const Real dl2[10] = {
            Real(0),
            Real(0),
            d32,
            Real(0),
            Real(0),
            p21 * d12,
            p11 * d22,
            p10 * d22,
            p20 * d12,
            p10 * p11 * d12,
        };

        for (std::size_t node = 0; node < 10u; ++node) {
            Real* g = rows[node];
            g[0u * output_stride + q] = dl1[node] - dl0[node];
            g[1u * output_stride + q] = dl2[node] - dl0[node];
            g[2u * output_stride + q] = Real(0);
        }
    }
}

void evaluate_hex_order1_values_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    Real* row0 = values_out + 0u * output_stride;
    Real* row1 = values_out + 1u * output_stride;
    Real* row2 = values_out + 2u * output_stride;
    Real* row3 = values_out + 3u * output_stride;
    Real* row4 = values_out + 4u * output_stride;
    Real* row5 = values_out + 5u * output_stride;
    Real* row6 = values_out + 6u * output_stride;
    Real* row7 = values_out + 7u * output_stride;

    const auto write_q = [&](std::size_t q) {
        const auto& xi = points[q];
        const Real lx = (Real(1) - xi[0]) * Real(0.5);
        const Real ly = (Real(1) - xi[1]) * Real(0.5);
        const Real lz = (Real(1) - xi[2]) * Real(0.5);
        const Real ux = (Real(1) + xi[0]) * Real(0.5);
        const Real uy = (Real(1) + xi[1]) * Real(0.5);
        const Real uz = (Real(1) + xi[2]) * Real(0.5);
        const Real lxly = lx * ly;
        const Real uxly = ux * ly;
        const Real uxuy = ux * uy;
        const Real lxuy = lx * uy;
        row0[q] = lxly * lz;
        row1[q] = uxly * lz;
        row2[q] = uxuy * lz;
        row3[q] = lxuy * lz;
        row4[q] = lxly * uz;
        row5[q] = uxly * uz;
        row6[q] = uxuy * uz;
        row7[q] = lxuy * uz;
    };
    if (points.size() == 4u) {
        write_q(0u);
        write_q(1u);
        write_q(2u);
        write_q(3u);
        return;
    }

    for (std::size_t q = 0; q < points.size(); ++q) {
        write_q(q);
    }
}

template <bool NeedValues, bool NeedGradients, bool NeedHessians>
void evaluate_hex_order1_outputs_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    constexpr Real half = Real(0.5);
    constexpr std::array<Real, 8> dx{{-half, half, half, -half, -half, half, half, -half}};
    constexpr std::array<Real, 8> dy{{-half, -half, half, half, -half, -half, half, half}};
    constexpr std::array<Real, 8> dz{{-half, -half, -half, -half, half, half, half, half}};

    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        const Real lx = (Real(1) - xi[0]) * half;
        const Real ly = (Real(1) - xi[1]) * half;
        const Real lz = (Real(1) - xi[2]) * half;
        const Real ux = (Real(1) + xi[0]) * half;
        const Real uy = (Real(1) + xi[1]) * half;
        const Real uz = (Real(1) + xi[2]) * half;
        const Real xval[8] = {lx, ux, ux, lx, lx, ux, ux, lx};
        const Real yval[8] = {ly, ly, uy, uy, ly, ly, uy, uy};
        const Real zval[8] = {lz, lz, lz, lz, uz, uz, uz, uz};

        for (std::size_t node = 0; node < 8u; ++node) {
            if constexpr (NeedValues) {
                values_out[node * output_stride + q] =
                    xval[node] * yval[node] * zval[node];
            }
            if constexpr (NeedGradients) {
                Real* SVMP_RESTRICT g = gradients_out + node * 3u * output_stride;
                g[0u * output_stride + q] = dx[node] * yval[node] * zval[node];
                g[1u * output_stride + q] = xval[node] * dy[node] * zval[node];
                g[2u * output_stride + q] = xval[node] * yval[node] * dz[node];
            }
            if constexpr (NeedHessians) {
                Real* SVMP_RESTRICT H = hessians_out + node * 9u * output_stride;
                const Real hxy = dx[node] * dy[node] * zval[node];
                const Real hxz = dx[node] * yval[node] * dz[node];
                const Real hyz = xval[node] * dy[node] * dz[node];
                H[0u * output_stride + q] = Real(0);
                H[1u * output_stride + q] = hxy;
                H[2u * output_stride + q] = hxz;
                H[3u * output_stride + q] = hxy;
                H[4u * output_stride + q] = Real(0);
                H[5u * output_stride + q] = hyz;
                H[6u * output_stride + q] = hxz;
                H[7u * output_stride + q] = hyz;
                H[8u * output_stride + q] = Real(0);
            }
        }
    }
}

void evaluate_quad_order1_values_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    Real* row0 = values_out + 0u * output_stride;
    Real* row1 = values_out + 1u * output_stride;
    Real* row2 = values_out + 2u * output_stride;
    Real* row3 = values_out + 3u * output_stride;

    if (points.size() == 4u && output_stride == 4u) {
        Real lx[4];
        Real ux[4];
        Real ly[4];
        Real uy[4];
        for (std::size_t q = 0; q < 4u; ++q) {
            const auto& xi = points[q];
            lx[q] = (Real(1) - xi[0]) * Real(0.5);
            ux[q] = (Real(1) + xi[0]) * Real(0.5);
            ly[q] = (Real(1) - xi[1]) * Real(0.5);
            uy[q] = (Real(1) + xi[1]) * Real(0.5);
        }
        row0[0] = lx[0] * ly[0];
        row0[1] = lx[1] * ly[1];
        row0[2] = lx[2] * ly[2];
        row0[3] = lx[3] * ly[3];
        row1[0] = ux[0] * ly[0];
        row1[1] = ux[1] * ly[1];
        row1[2] = ux[2] * ly[2];
        row1[3] = ux[3] * ly[3];
        row2[0] = ux[0] * uy[0];
        row2[1] = ux[1] * uy[1];
        row2[2] = ux[2] * uy[2];
        row2[3] = ux[3] * uy[3];
        row3[0] = lx[0] * uy[0];
        row3[1] = lx[1] * uy[1];
        row3[2] = lx[2] * uy[2];
        row3[3] = lx[3] * uy[3];
        return;
    }

    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        const Real lx = (Real(1) - xi[0]) * Real(0.5);
        const Real ly = (Real(1) - xi[1]) * Real(0.5);
        const Real ux = (Real(1) + xi[0]) * Real(0.5);
        const Real uy = (Real(1) + xi[1]) * Real(0.5);
        row0[q] = lx * ly;
        row1[q] = ux * ly;
        row2[q] = ux * uy;
        row3[q] = lx * uy;
    }
}

void evaluate_quad_order1_gradients_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT gradients_out) {
    Real* row0 = gradients_out + 0u * 3u * output_stride;
    Real* row1 = gradients_out + 1u * 3u * output_stride;
    Real* row2 = gradients_out + 2u * 3u * output_stride;
    Real* row3 = gradients_out + 3u * 3u * output_stride;

    if (points.size() == 4u) {
        Real lx[4];
        Real ly[4];
        Real ux[4];
        Real uy[4];
        for (std::size_t q = 0; q < 4u; ++q) {
            const auto& xi = points[q];
            lx[q] = (Real(1) - xi[0]) * Real(0.5);
            ly[q] = (Real(1) - xi[1]) * Real(0.5);
            ux[q] = (Real(1) + xi[0]) * Real(0.5);
            uy[q] = (Real(1) + xi[1]) * Real(0.5);
        }

        auto write_component = [](Real* SVMP_RESTRICT row,
                                  Real a0,
                                  Real a1,
                                  Real a2,
                                  Real a3) {
            row[0] = a0;
            row[1] = a1;
            row[2] = a2;
            row[3] = a3;
        };

        write_component(row0, Real(-0.5) * ly[0], Real(-0.5) * ly[1],
                        Real(-0.5) * ly[2], Real(-0.5) * ly[3]);
        write_component(row0 + output_stride, Real(-0.5) * lx[0], Real(-0.5) * lx[1],
                        Real(-0.5) * lx[2], Real(-0.5) * lx[3]);
        write_component(row0 + 2u * output_stride, Real(0), Real(0), Real(0), Real(0));

        write_component(row1, Real(0.5) * ly[0], Real(0.5) * ly[1],
                        Real(0.5) * ly[2], Real(0.5) * ly[3]);
        write_component(row1 + output_stride, Real(-0.5) * ux[0], Real(-0.5) * ux[1],
                        Real(-0.5) * ux[2], Real(-0.5) * ux[3]);
        write_component(row1 + 2u * output_stride, Real(0), Real(0), Real(0), Real(0));

        write_component(row2, Real(0.5) * uy[0], Real(0.5) * uy[1],
                        Real(0.5) * uy[2], Real(0.5) * uy[3]);
        write_component(row2 + output_stride, Real(0.5) * ux[0], Real(0.5) * ux[1],
                        Real(0.5) * ux[2], Real(0.5) * ux[3]);
        write_component(row2 + 2u * output_stride, Real(0), Real(0), Real(0), Real(0));

        write_component(row3, Real(-0.5) * uy[0], Real(-0.5) * uy[1],
                        Real(-0.5) * uy[2], Real(-0.5) * uy[3]);
        write_component(row3 + output_stride, Real(0.5) * lx[0], Real(0.5) * lx[1],
                        Real(0.5) * lx[2], Real(0.5) * lx[3]);
        write_component(row3 + 2u * output_stride, Real(0), Real(0), Real(0), Real(0));
        return;
    }

    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        const Real lx = (Real(1) - xi[0]) * Real(0.5);
        const Real ly = (Real(1) - xi[1]) * Real(0.5);
        const Real ux = (Real(1) + xi[0]) * Real(0.5);
        const Real uy = (Real(1) + xi[1]) * Real(0.5);
        row0[0u * output_stride + q] = Real(-0.5) * ly;
        row0[1u * output_stride + q] = Real(-0.5) * lx;
        row0[2u * output_stride + q] = Real(0);
        row1[0u * output_stride + q] = Real( 0.5) * ly;
        row1[1u * output_stride + q] = Real(-0.5) * ux;
        row1[2u * output_stride + q] = Real(0);
        row2[0u * output_stride + q] = Real( 0.5) * uy;
        row2[1u * output_stride + q] = Real( 0.5) * ux;
        row2[2u * output_stride + q] = Real(0);
        row3[0u * output_stride + q] = Real(-0.5) * uy;
        row3[1u * output_stride + q] = Real( 0.5) * lx;
        row3[2u * output_stride + q] = Real(0);
    }
}

inline void write_quad_order1_hessian_q(
    Real* SVMP_RESTRICT row,
    std::size_t output_stride,
    std::size_t q,
    Real xy) {
    row[0u * output_stride + q] = Real(0);
    row[1u * output_stride + q] = xy;
    row[2u * output_stride + q] = Real(0);
    row[3u * output_stride + q] = xy;
    row[4u * output_stride + q] = Real(0);
    row[5u * output_stride + q] = Real(0);
    row[6u * output_stride + q] = Real(0);
    row[7u * output_stride + q] = Real(0);
    row[8u * output_stride + q] = Real(0);
}

void evaluate_quad_order1_hessians_strided(
    std::size_t num_qpts,
    std::size_t output_stride,
    Real* SVMP_RESTRICT hessians_out) {
    Real* row0 = hessians_out + 0u * 9u * output_stride;
    Real* row1 = hessians_out + 1u * 9u * output_stride;
    Real* row2 = hessians_out + 2u * 9u * output_stride;
    Real* row3 = hessians_out + 3u * 9u * output_stride;

    constexpr Real positive = Real(0.25);
    constexpr Real negative = Real(-0.25);
    for (std::size_t q = 0; q < num_qpts; ++q) {
        write_quad_order1_hessian_q(row0, output_stride, q, positive);
        write_quad_order1_hessian_q(row1, output_stride, q, negative);
        write_quad_order1_hessian_q(row2, output_stride, q, positive);
        write_quad_order1_hessian_q(row3, output_stride, q, negative);
    }
}

template <std::size_t Q>
inline void write_quad_order1_all_q4(
    std::size_t output_stride,
    std::size_t i,
    std::size_t j,
    const Real lx[4][2],
    const Real ly[4][2],
    Real* SVMP_RESTRICT value_row,
    Real* SVMP_RESTRICT grad_row,
    Real* SVMP_RESTRICT hess_row) {
    const Real xv = lx[Q][i];
    const Real yv = ly[Q][j];
    const Real xd = (i == 0u) ? Real(-0.5) : Real(0.5);
    const Real yd = (j == 0u) ? Real(-0.5) : Real(0.5);
    const Real hxy = xd * yd;

    value_row[Q] = xv * yv;
    grad_row[0u * output_stride + Q] = xd * yv;
    grad_row[1u * output_stride + Q] = xv * yd;
    grad_row[2u * output_stride + Q] = Real(0);
    hess_row[0u * output_stride + Q] = Real(0);
    hess_row[4u * output_stride + Q] = Real(0);
    hess_row[8u * output_stride + Q] = Real(0);
    hess_row[1u * output_stride + Q] = hxy;
    hess_row[3u * output_stride + Q] = hxy;
    hess_row[2u * output_stride + Q] = Real(0);
    hess_row[6u * output_stride + Q] = Real(0);
    hess_row[5u * output_stride + Q] = Real(0);
    hess_row[7u * output_stride + Q] = Real(0);
}

void evaluate_quad_order1_all_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    constexpr auto node_axes = detail::make_quad_tensor_node_axes<1>();

    Real lx[4][2];
    Real ly[4][2];
    for (std::size_t q = 0; q < 4u; ++q) {
        const auto& xi = points[q];
        lx[q][0] = (Real(1) - xi[0]) * Real(0.5);
        lx[q][1] = (Real(1) + xi[0]) * Real(0.5);
        ly[q][0] = (Real(1) - xi[1]) * Real(0.5);
        ly[q][1] = (Real(1) + xi[1]) * Real(0.5);
    }

    for (std::size_t node = 0; node < node_axes.size(); ++node) {
        const auto& axes = node_axes[node];
        const std::size_t i = axes[0];
        const std::size_t j = axes[1];
        Real* value_row = values_out + node * output_stride;
        Real* grad_row = gradients_out + node * 3u * output_stride;
        Real* hess_row = hessians_out + node * 9u * output_stride;
        write_quad_order1_all_q4<0u>(
            output_stride, i, j, lx, ly, value_row, grad_row, hess_row);
        write_quad_order1_all_q4<1u>(
            output_stride, i, j, lx, ly, value_row, grad_row, hess_row);
        write_quad_order1_all_q4<2u>(
            output_stride, i, j, lx, ly, value_row, grad_row, hess_row);
        write_quad_order1_all_q4<3u>(
            output_stride, i, j, lx, ly, value_row, grad_row, hess_row);
    }
}

void evaluate_quad_order2_values_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    Real* row0 = values_out + 0u * output_stride;
    Real* row1 = values_out + 1u * output_stride;
    Real* row2 = values_out + 2u * output_stride;
    Real* row3 = values_out + 3u * output_stride;
    Real* row4 = values_out + 4u * output_stride;
    Real* row5 = values_out + 5u * output_stride;
    Real* row6 = values_out + 6u * output_stride;
    Real* row7 = values_out + 7u * output_stride;
    Real* row8 = values_out + 8u * output_stride;

    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        const Real x = xi[0];
        const Real y = xi[1];
        const Real x0 = x * (x - Real(1)) * Real(0.5);
        const Real x1 = x * (x + Real(1)) * Real(0.5);
        const Real x2 = Real(1) - x * x;
        const Real y0 = y * (y - Real(1)) * Real(0.5);
        const Real y1 = y * (y + Real(1)) * Real(0.5);
        const Real y2 = Real(1) - y * y;

        row0[q] = x0 * y0;
        row1[q] = x1 * y0;
        row2[q] = x1 * y1;
        row3[q] = x0 * y1;
        row4[q] = x2 * y0;
        row5[q] = x1 * y2;
        row6[q] = x2 * y1;
        row7[q] = x0 * y2;
        row8[q] = x2 * y2;
    }
}

inline void write_quad_order2_gradient_q(
    Real* SVMP_RESTRICT row,
    std::size_t output_stride,
    std::size_t q,
    Real dx,
    Real dy) {
    row[0u * output_stride + q] = dx;
    row[1u * output_stride + q] = dy;
    row[2u * output_stride + q] = Real(0);
}

void evaluate_quad_order2_gradients_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT gradients_out) {
    if (points.size() == 4u) {
        Real xv[4][3];
        Real yv[4][3];
        Real xd[4][3];
        Real yd[4][3];
        for (std::size_t q = 0; q < 4u; ++q) {
            const auto& xi = points[q];
            const Real x = xi[0];
            const Real y = xi[1];
            xv[q][0] = x * (x - Real(1)) * Real(0.5);
            xv[q][1] = x * (x + Real(1)) * Real(0.5);
            xv[q][2] = Real(1) - x * x;
            yv[q][0] = y * (y - Real(1)) * Real(0.5);
            yv[q][1] = y * (y + Real(1)) * Real(0.5);
            yv[q][2] = Real(1) - y * y;
            xd[q][0] = x - Real(0.5);
            xd[q][1] = x + Real(0.5);
            xd[q][2] = Real(-2) * x;
            yd[q][0] = y - Real(0.5);
            yd[q][1] = y + Real(0.5);
            yd[q][2] = Real(-2) * y;
        }

        auto write_node = [&](std::size_t node, std::size_t i, std::size_t j) {
            Real* SVMP_RESTRICT row = gradients_out + node * 3u * output_stride;
            row[0u] = xd[0][i] * yv[0][j];
            row[1u] = xd[1][i] * yv[1][j];
            row[2u] = xd[2][i] * yv[2][j];
            row[3u] = xd[3][i] * yv[3][j];
            row[output_stride + 0u] = xv[0][i] * yd[0][j];
            row[output_stride + 1u] = xv[1][i] * yd[1][j];
            row[output_stride + 2u] = xv[2][i] * yd[2][j];
            row[output_stride + 3u] = xv[3][i] * yd[3][j];
            row[2u * output_stride + 0u] = Real(0);
            row[2u * output_stride + 1u] = Real(0);
            row[2u * output_stride + 2u] = Real(0);
            row[2u * output_stride + 3u] = Real(0);
        };

        write_node(0u, 0u, 0u);
        write_node(1u, 1u, 0u);
        write_node(2u, 1u, 1u);
        write_node(3u, 0u, 1u);
        write_node(4u, 2u, 0u);
        write_node(5u, 1u, 2u);
        write_node(6u, 2u, 1u);
        write_node(7u, 0u, 2u);
        write_node(8u, 2u, 2u);
        return;
    }

    Real* row0 = gradients_out + 0u * 3u * output_stride;
    Real* row1 = gradients_out + 1u * 3u * output_stride;
    Real* row2 = gradients_out + 2u * 3u * output_stride;
    Real* row3 = gradients_out + 3u * 3u * output_stride;
    Real* row4 = gradients_out + 4u * 3u * output_stride;
    Real* row5 = gradients_out + 5u * 3u * output_stride;
    Real* row6 = gradients_out + 6u * 3u * output_stride;
    Real* row7 = gradients_out + 7u * 3u * output_stride;
    Real* row8 = gradients_out + 8u * 3u * output_stride;

    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        const Real x = xi[0];
        const Real y = xi[1];
        const Real x0 = x * (x - Real(1)) * Real(0.5);
        const Real x1 = x * (x + Real(1)) * Real(0.5);
        const Real x2 = Real(1) - x * x;
        const Real y0 = y * (y - Real(1)) * Real(0.5);
        const Real y1 = y * (y + Real(1)) * Real(0.5);
        const Real y2 = Real(1) - y * y;
        const Real dx0 = x - Real(0.5);
        const Real dx1 = x + Real(0.5);
        const Real dx2 = Real(-2) * x;
        const Real dy0 = y - Real(0.5);
        const Real dy1 = y + Real(0.5);
        const Real dy2 = Real(-2) * y;

        write_quad_order2_gradient_q(row0, output_stride, q, dx0 * y0, x0 * dy0);
        write_quad_order2_gradient_q(row1, output_stride, q, dx1 * y0, x1 * dy0);
        write_quad_order2_gradient_q(row2, output_stride, q, dx1 * y1, x1 * dy1);
        write_quad_order2_gradient_q(row3, output_stride, q, dx0 * y1, x0 * dy1);
        write_quad_order2_gradient_q(row4, output_stride, q, dx2 * y0, x2 * dy0);
        write_quad_order2_gradient_q(row5, output_stride, q, dx1 * y2, x1 * dy2);
        write_quad_order2_gradient_q(row6, output_stride, q, dx2 * y1, x2 * dy1);
        write_quad_order2_gradient_q(row7, output_stride, q, dx0 * y2, x0 * dy2);
        write_quad_order2_gradient_q(row8, output_stride, q, dx2 * y2, x2 * dy2);
    }
}

inline void write_quad_order2_hessian_q(
    Real* SVMP_RESTRICT row,
    std::size_t output_stride,
    std::size_t q,
    Real hxx,
    Real hxy,
    Real hyy) {
    row[0u * output_stride + q] = hxx;
    row[1u * output_stride + q] = hxy;
    row[2u * output_stride + q] = Real(0);
    row[3u * output_stride + q] = hxy;
    row[4u * output_stride + q] = hyy;
    row[5u * output_stride + q] = Real(0);
    row[6u * output_stride + q] = Real(0);
    row[7u * output_stride + q] = Real(0);
    row[8u * output_stride + q] = Real(0);
}

void evaluate_quad_order2_hessians_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT hessians_out) {
    Real* row0 = hessians_out + 0u * 9u * output_stride;
    Real* row1 = hessians_out + 1u * 9u * output_stride;
    Real* row2 = hessians_out + 2u * 9u * output_stride;
    Real* row3 = hessians_out + 3u * 9u * output_stride;
    Real* row4 = hessians_out + 4u * 9u * output_stride;
    Real* row5 = hessians_out + 5u * 9u * output_stride;
    Real* row6 = hessians_out + 6u * 9u * output_stride;
    Real* row7 = hessians_out + 7u * 9u * output_stride;
    Real* row8 = hessians_out + 8u * 9u * output_stride;

    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        const Real x = xi[0];
        const Real y = xi[1];
        const Real x0 = x * (x - Real(1)) * Real(0.5);
        const Real x1 = x * (x + Real(1)) * Real(0.5);
        const Real x2 = Real(1) - x * x;
        const Real y0 = y * (y - Real(1)) * Real(0.5);
        const Real y1 = y * (y + Real(1)) * Real(0.5);
        const Real y2 = Real(1) - y * y;
        const Real dx0 = x - Real(0.5);
        const Real dx1 = x + Real(0.5);
        const Real dx2 = Real(-2) * x;
        const Real dy0 = y - Real(0.5);
        const Real dy1 = y + Real(0.5);
        const Real dy2 = Real(-2) * y;

        write_quad_order2_hessian_q(row0, output_stride, q, y0, dx0 * dy0, x0);
        write_quad_order2_hessian_q(row1, output_stride, q, y0, dx1 * dy0, x1);
        write_quad_order2_hessian_q(row2, output_stride, q, y1, dx1 * dy1, x1);
        write_quad_order2_hessian_q(row3, output_stride, q, y1, dx0 * dy1, x0);
        write_quad_order2_hessian_q(row4, output_stride, q, Real(-2) * y0, dx2 * dy0, x2);
        write_quad_order2_hessian_q(row5, output_stride, q, y2, dx1 * dy2, Real(-2) * x1);
        write_quad_order2_hessian_q(row6, output_stride, q, Real(-2) * y1, dx2 * dy1, x2);
        write_quad_order2_hessian_q(row7, output_stride, q, y2, dx0 * dy2, Real(-2) * x0);
        write_quad_order2_hessian_q(row8, output_stride, q, Real(-2) * y2, dx2 * dy2, Real(-2) * x2);
    }
}

inline void fill_order3_axis_values(Real x, Real* SVMP_RESTRICT values) {
    const Real x2 = x * x;
    values[0] = Real(-9.0 / 16.0) * (x - Real(1)) * (x2 - Real(1.0 / 9.0));
    values[1] = Real( 9.0 / 16.0) * (x + Real(1)) * (x2 - Real(1.0 / 9.0));
    values[2] = Real(27.0 / 16.0) * (x2 - Real(1)) * (x - Real(1.0 / 3.0));
    values[3] = Real(-27.0 / 16.0) * (x2 - Real(1)) * (x + Real(1.0 / 3.0));
}

inline void fill_order3_axis_value_scalars(Real x,
                                           Real& v0,
                                           Real& v1,
                                           Real& v2,
                                           Real& v3) {
    const Real x2 = x * x;
    v0 = Real(-9.0 / 16.0) * (x - Real(1)) * (x2 - Real(1.0 / 9.0));
    v1 = Real( 9.0 / 16.0) * (x + Real(1)) * (x2 - Real(1.0 / 9.0));
    v2 = Real(27.0 / 16.0) * (x2 - Real(1)) * (x - Real(1.0 / 3.0));
    v3 = Real(-27.0 / 16.0) * (x2 - Real(1)) * (x + Real(1.0 / 3.0));
}

void evaluate_line_order1_values_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    Real* row0 = values_out + 0u * output_stride;
    Real* row1 = values_out + 1u * output_stride;
    for (std::size_t q = 0; q < 4u; ++q) {
        const Real x = points[q][0];
        row0[q] = (Real(1) - x) * Real(0.5);
        row1[q] = (Real(1) + x) * Real(0.5);
    }
}

void evaluate_line_order2_values_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    Real* row0 = values_out + 0u * output_stride;
    Real* row1 = values_out + 1u * output_stride;
    Real* row2 = values_out + 2u * output_stride;
    for (std::size_t q = 0; q < 4u; ++q) {
        const Real x = points[q][0];
        row0[q] = x * (x - Real(1)) * Real(0.5);
        row1[q] = x * (x + Real(1)) * Real(0.5);
        row2[q] = Real(1) - x * x;
    }
}

void evaluate_line_order3_values_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    Real* row0 = values_out + 0u * output_stride;
    Real* row1 = values_out + 1u * output_stride;
    Real* row2 = values_out + 2u * output_stride;
    Real* row3 = values_out + 3u * output_stride;
    for (std::size_t q = 0; q < 4u; ++q) {
        Real values[4];
        fill_order3_axis_values(points[q][0], values);
        row0[q] = values[0];
        row1[q] = values[1];
        row2[q] = values[2];
        row3[q] = values[3];
    }
}

inline void fill_order3_axis_values_first(Real x,
                                          Real* SVMP_RESTRICT values,
                                          Real* SVMP_RESTRICT first);

inline void fill_order3_axis_values_first_second(Real x,
                                                 Real* SVMP_RESTRICT values,
                                                 Real* SVMP_RESTRICT first,
                                                 Real* SVMP_RESTRICT second);

inline void write_line_gradient_q4_row(Real* SVMP_RESTRICT row,
                                       std::size_t output_stride,
                                       Real g0,
                                       Real g1,
                                       Real g2,
                                       Real g3) {
    row[0] = g0;
    row[1] = g1;
    row[2] = g2;
    row[3] = g3;
    row[output_stride + 0u] = Real(0);
    row[output_stride + 1u] = Real(0);
    row[output_stride + 2u] = Real(0);
    row[output_stride + 3u] = Real(0);
    row[2u * output_stride + 0u] = Real(0);
    row[2u * output_stride + 1u] = Real(0);
    row[2u * output_stride + 2u] = Real(0);
    row[2u * output_stride + 3u] = Real(0);
}

inline void write_line_hessian_q4_row(Real* SVMP_RESTRICT row,
                                      std::size_t output_stride,
                                      Real h0,
                                      Real h1,
                                      Real h2,
                                      Real h3) {
    row[0] = h0;
    row[1] = h1;
    row[2] = h2;
    row[3] = h3;
    for (std::size_t component = 1u; component < 9u; ++component) {
        Real* slot = row + component * output_stride;
        slot[0] = Real(0);
        slot[1] = Real(0);
        slot[2] = Real(0);
        slot[3] = Real(0);
    }
}

SVMP_LAGRANGE_NOINLINE void evaluate_line_order1_gradients_q4(
    std::size_t output_stride,
    Real* SVMP_RESTRICT gradients_out) {
    write_line_gradient_q4_row(gradients_out + 0u * 3u * output_stride,
                               output_stride,
                               Real(-0.5), Real(-0.5), Real(-0.5), Real(-0.5));
    write_line_gradient_q4_row(gradients_out + 1u * 3u * output_stride,
                               output_stride,
                               Real(0.5), Real(0.5), Real(0.5), Real(0.5));
}

SVMP_LAGRANGE_NOINLINE void evaluate_line_order1_hessians_q4(
    std::size_t output_stride,
    Real* SVMP_RESTRICT hessians_out) {
    write_line_hessian_q4_row(hessians_out + 0u * 9u * output_stride,
                              output_stride, Real(0), Real(0), Real(0), Real(0));
    write_line_hessian_q4_row(hessians_out + 1u * 9u * output_stride,
                              output_stride, Real(0), Real(0), Real(0), Real(0));
}

SVMP_LAGRANGE_NOINLINE void evaluate_line_order1_all_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    evaluate_line_order1_values_q4(points, output_stride, values_out);
    evaluate_line_order1_gradients_q4(output_stride, gradients_out);
    evaluate_line_order1_hessians_q4(output_stride, hessians_out);
}

SVMP_LAGRANGE_NOINLINE void evaluate_line_order2_gradients_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT gradients_out) {
    const Real x0 = points[0][0];
    const Real x1 = points[1][0];
    const Real x2 = points[2][0];
    const Real x3 = points[3][0];
    write_line_gradient_q4_row(gradients_out + 0u * 3u * output_stride,
                               output_stride,
                               x0 - Real(0.5), x1 - Real(0.5),
                               x2 - Real(0.5), x3 - Real(0.5));
    write_line_gradient_q4_row(gradients_out + 1u * 3u * output_stride,
                               output_stride,
                               x0 + Real(0.5), x1 + Real(0.5),
                               x2 + Real(0.5), x3 + Real(0.5));
    write_line_gradient_q4_row(gradients_out + 2u * 3u * output_stride,
                               output_stride,
                               Real(-2) * x0, Real(-2) * x1,
                               Real(-2) * x2, Real(-2) * x3);
}

SVMP_LAGRANGE_NOINLINE void evaluate_line_order2_hessians_q4(
    std::size_t output_stride,
    Real* SVMP_RESTRICT hessians_out) {
    write_line_hessian_q4_row(hessians_out + 0u * 9u * output_stride,
                              output_stride, Real(1), Real(1), Real(1), Real(1));
    write_line_hessian_q4_row(hessians_out + 1u * 9u * output_stride,
                              output_stride, Real(1), Real(1), Real(1), Real(1));
    write_line_hessian_q4_row(hessians_out + 2u * 9u * output_stride,
                              output_stride, Real(-2), Real(-2), Real(-2), Real(-2));
}

SVMP_LAGRANGE_NOINLINE void evaluate_line_order2_all_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    evaluate_line_order2_values_q4(points, output_stride, values_out);
    evaluate_line_order2_gradients_q4(points, output_stride, gradients_out);
    evaluate_line_order2_hessians_q4(output_stride, hessians_out);
}

SVMP_LAGRANGE_NOINLINE void evaluate_line_order3_gradients_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT gradients_out) {
    Real first[4][4];
    for (std::size_t q = 0; q < 4u; ++q) {
        Real values[4];
        fill_order3_axis_values_first(points[q][0], values, first[q]);
    }
    for (std::size_t node = 0; node < 4u; ++node) {
        write_line_gradient_q4_row(gradients_out + node * 3u * output_stride,
                                   output_stride,
                                   first[0][node], first[1][node],
                                   first[2][node], first[3][node]);
    }
}

SVMP_LAGRANGE_NOINLINE void evaluate_line_order3_hessians_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT hessians_out) {
    Real second[4][4];
    for (std::size_t q = 0; q < 4u; ++q) {
        Real values[4];
        Real first[4];
        fill_order3_axis_values_first_second(points[q][0], values, first, second[q]);
    }
    for (std::size_t node = 0; node < 4u; ++node) {
        write_line_hessian_q4_row(hessians_out + node * 9u * output_stride,
                                  output_stride,
                                  second[0][node], second[1][node],
                                  second[2][node], second[3][node]);
    }
}

SVMP_LAGRANGE_NOINLINE void evaluate_line_order3_all_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    Real values[4][4];
    Real first[4][4];
    Real second[4][4];
    for (std::size_t q = 0; q < 4u; ++q) {
        fill_order3_axis_values_first_second(points[q][0], values[q], first[q], second[q]);
    }
    for (std::size_t node = 0; node < 4u; ++node) {
        Real* value_row = values_out + node * output_stride;
        value_row[0] = values[0][node];
        value_row[1] = values[1][node];
        value_row[2] = values[2][node];
        value_row[3] = values[3][node];
        write_line_gradient_q4_row(gradients_out + node * 3u * output_stride,
                                   output_stride,
                                   first[0][node], first[1][node],
                                   first[2][node], first[3][node]);
        write_line_hessian_q4_row(hessians_out + node * 9u * output_stride,
                                  output_stride,
                                  second[0][node], second[1][node],
                                  second[2][node], second[3][node]);
    }
}

inline void fill_order3_axis_values_first(Real x,
                                          Real* SVMP_RESTRICT values,
                                          Real* SVMP_RESTRICT first) {
    fill_order3_axis_values(x, values);
    const Real x2 = x * x;
    first[0] = Real(-9.0 / 16.0) * (Real(3) * x2 - Real(2) * x - Real(1.0 / 9.0));
    first[1] = Real( 9.0 / 16.0) * (Real(3) * x2 + Real(2) * x - Real(1.0 / 9.0));
    first[2] = Real(27.0 / 16.0) * (Real(3) * x2 - Real(2.0 / 3.0) * x - Real(1));
    first[3] = Real(-27.0 / 16.0) * (Real(3) * x2 + Real(2.0 / 3.0) * x - Real(1));
}

inline void fill_order3_axis_values_first_second(Real x,
                                                 Real* SVMP_RESTRICT values,
                                                 Real* SVMP_RESTRICT first,
                                                 Real* SVMP_RESTRICT second) {
    fill_order3_axis_values_first(x, values, first);
    second[0] = Real(-9.0 / 16.0) * (Real(6) * x - Real(2));
    second[1] = Real( 9.0 / 16.0) * (Real(6) * x + Real(2));
    second[2] = Real(27.0 / 16.0) * (Real(6) * x - Real(2.0 / 3.0));
    second[3] = Real(-27.0 / 16.0) * (Real(6) * x + Real(2.0 / 3.0));
}

inline void write_quad_order3_value_row_q4(Real* SVMP_RESTRICT row,
                                           const Real lx[4][4],
                                           const Real ly[4][4],
                                           std::size_t i,
                                           std::size_t j) {
    row[0] = lx[0][i] * ly[0][j];
    row[1] = lx[1][i] * ly[1][j];
    row[2] = lx[2][i] * ly[2][j];
    row[3] = lx[3][i] * ly[3][j];
}

void evaluate_quad_order3_values_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    if (output_stride == 4u) {
        Real* row0 = values_out + 0u * 4u;
        Real* row1 = values_out + 1u * 4u;
        Real* row2 = values_out + 2u * 4u;
        Real* row3 = values_out + 3u * 4u;
        Real* row4 = values_out + 4u * 4u;
        Real* row5 = values_out + 5u * 4u;
        Real* row6 = values_out + 6u * 4u;
        Real* row7 = values_out + 7u * 4u;
        Real* row8 = values_out + 8u * 4u;
        Real* row9 = values_out + 9u * 4u;
        Real* row10 = values_out + 10u * 4u;
        Real* row11 = values_out + 11u * 4u;
        Real* row12 = values_out + 12u * 4u;
        Real* row13 = values_out + 13u * 4u;
        Real* row14 = values_out + 14u * 4u;
        Real* row15 = values_out + 15u * 4u;

        auto write_q = [&](std::size_t q) {
            const auto& xi = points[q];
            Real x0;
            Real x1;
            Real x2;
            Real x3;
            Real y0;
            Real y1;
            Real y2;
            Real y3;
            fill_order3_axis_value_scalars(xi[0], x0, x1, x2, x3);
            fill_order3_axis_value_scalars(xi[1], y0, y1, y2, y3);
            row0[q] = x0 * y0;
            row1[q] = x1 * y0;
            row2[q] = x1 * y1;
            row3[q] = x0 * y1;
            row4[q] = x2 * y0;
            row5[q] = x3 * y0;
            row6[q] = x1 * y2;
            row7[q] = x1 * y3;
            row8[q] = x3 * y1;
            row9[q] = x2 * y1;
            row10[q] = x0 * y3;
            row11[q] = x0 * y2;
            row12[q] = x2 * y2;
            row13[q] = x3 * y2;
            row14[q] = x2 * y3;
            row15[q] = x3 * y3;
        };

        write_q(0u);
        write_q(1u);
        write_q(2u);
        write_q(3u);
        return;
    }

    Real lx[4][4];
    Real ly[4][4];
    for (std::size_t q = 0; q < 4u; ++q) {
        const auto& xi = points[q];
        fill_order3_axis_values(xi[0], lx[q]);
        fill_order3_axis_values(xi[1], ly[q]);
    }

    write_quad_order3_value_row_q4(values_out + 0u * output_stride, lx, ly, 0u, 0u);
    write_quad_order3_value_row_q4(values_out + 1u * output_stride, lx, ly, 1u, 0u);
    write_quad_order3_value_row_q4(values_out + 2u * output_stride, lx, ly, 1u, 1u);
    write_quad_order3_value_row_q4(values_out + 3u * output_stride, lx, ly, 0u, 1u);
    write_quad_order3_value_row_q4(values_out + 4u * output_stride, lx, ly, 2u, 0u);
    write_quad_order3_value_row_q4(values_out + 5u * output_stride, lx, ly, 3u, 0u);
    write_quad_order3_value_row_q4(values_out + 6u * output_stride, lx, ly, 1u, 2u);
    write_quad_order3_value_row_q4(values_out + 7u * output_stride, lx, ly, 1u, 3u);
    write_quad_order3_value_row_q4(values_out + 8u * output_stride, lx, ly, 3u, 1u);
    write_quad_order3_value_row_q4(values_out + 9u * output_stride, lx, ly, 2u, 1u);
    write_quad_order3_value_row_q4(values_out + 10u * output_stride, lx, ly, 0u, 3u);
    write_quad_order3_value_row_q4(values_out + 11u * output_stride, lx, ly, 0u, 2u);
    write_quad_order3_value_row_q4(values_out + 12u * output_stride, lx, ly, 2u, 2u);
    write_quad_order3_value_row_q4(values_out + 13u * output_stride, lx, ly, 3u, 2u);
    write_quad_order3_value_row_q4(values_out + 14u * output_stride, lx, ly, 2u, 3u);
    write_quad_order3_value_row_q4(values_out + 15u * output_stride, lx, ly, 3u, 3u);
}

void evaluate_quad_order3_values_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    if (points.size() == 4u) {
        evaluate_quad_order3_values_q4(points, output_stride, values_out);
        return;
    }

    Real* row0 = values_out + 0u * output_stride;
    Real* row1 = values_out + 1u * output_stride;
    Real* row2 = values_out + 2u * output_stride;
    Real* row3 = values_out + 3u * output_stride;
    Real* row4 = values_out + 4u * output_stride;
    Real* row5 = values_out + 5u * output_stride;
    Real* row6 = values_out + 6u * output_stride;
    Real* row7 = values_out + 7u * output_stride;
    Real* row8 = values_out + 8u * output_stride;
    Real* row9 = values_out + 9u * output_stride;
    Real* row10 = values_out + 10u * output_stride;
    Real* row11 = values_out + 11u * output_stride;
    Real* row12 = values_out + 12u * output_stride;
    Real* row13 = values_out + 13u * output_stride;
    Real* row14 = values_out + 14u * output_stride;
    Real* row15 = values_out + 15u * output_stride;

    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        Real lx[4];
        Real ly[4];
        fill_order3_axis_values(xi[0], lx);
        fill_order3_axis_values(xi[1], ly);
        row0[q] = lx[0] * ly[0];
        row1[q] = lx[1] * ly[0];
        row2[q] = lx[1] * ly[1];
        row3[q] = lx[0] * ly[1];
        row4[q] = lx[2] * ly[0];
        row5[q] = lx[3] * ly[0];
        row6[q] = lx[1] * ly[2];
        row7[q] = lx[1] * ly[3];
        row8[q] = lx[3] * ly[1];
        row9[q] = lx[2] * ly[1];
        row10[q] = lx[0] * ly[3];
        row11[q] = lx[0] * ly[2];
        row12[q] = lx[2] * ly[2];
        row13[q] = lx[3] * ly[2];
        row14[q] = lx[2] * ly[3];
        row15[q] = lx[3] * ly[3];
    }
}

template <std::size_t N>
inline void write_quad_gradient_row_q4(
    Real* SVMP_RESTRICT row,
    std::size_t output_stride,
    const Real (&lx)[4][N],
    const Real (&ly)[4][N],
    const Real (&dx)[4][N],
    const Real (&dy)[4][N],
    std::size_t i,
    std::size_t j) {
    row[0u] = dx[0][i] * ly[0][j];
    row[1u] = dx[1][i] * ly[1][j];
    row[2u] = dx[2][i] * ly[2][j];
    row[3u] = dx[3][i] * ly[3][j];
    row[output_stride + 0u] = lx[0][i] * dy[0][j];
    row[output_stride + 1u] = lx[1][i] * dy[1][j];
    row[output_stride + 2u] = lx[2][i] * dy[2][j];
    row[output_stride + 3u] = lx[3][i] * dy[3][j];
    row[2u * output_stride + 0u] = Real(0);
    row[2u * output_stride + 1u] = Real(0);
    row[2u * output_stride + 2u] = Real(0);
    row[2u * output_stride + 3u] = Real(0);
}

inline void fill_order4_axis_values_first(Real x,
                                          Real* SVMP_RESTRICT values,
                                          Real* SVMP_RESTRICT first) {
    const Real r = (x + Real(1)) * Real(2);
    const Real r2 = r * r;
    const Real r3 = r2 * r;
    const Real f0 = r;
    const Real f1 = r - Real(1);
    const Real f2 = r - Real(2);
    const Real f3 = r - Real(3);
    const Real f4 = r - Real(4);
    const Real f01 = f0 * f1;
    const Real f12 = f1 * f2;
    const Real f23 = f2 * f3;
    const Real f34 = f3 * f4;

    values[0] = (f12 * f34) / Real(24);
    values[1] = (f01 * f23) / Real(24);
    values[2] = -(f0 * f2 * f34) / Real(6);
    values[3] = (f01 * f34) / Real(4);
    values[4] = -(f01 * f2 * f4) / Real(6);

    first[0] = (Real(4) * r3 - Real(30) * r2 + Real(70) * r - Real(50)) / Real(12);
    first[1] = (Real(4) * r3 - Real(18) * r2 + Real(22) * r - Real(6)) / Real(12);
    first[2] = (-Real(4) * r3 + Real(27) * r2 - Real(52) * r + Real(24)) / Real(3);
    first[3] = Real(2) * r3 - Real(12) * r2 + Real(19) * r - Real(6);
    first[4] = (-Real(4) * r3 + Real(21) * r2 - Real(28) * r + Real(8)) / Real(3);
}

SVMP_LAGRANGE_NOINLINE void evaluate_quad_order3_gradients_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT gradients_out) {
    constexpr auto node_axes = detail::make_quad_tensor_node_axes<3>();

    Real lx[4][4];
    Real ly[4][4];
    Real dx[4][4];
    Real dy[4][4];
    for (std::size_t q = 0; q < 4u; ++q) {
        const auto& xi = points[q];
        fill_order3_axis_values_first(xi[0], lx[q], dx[q]);
        fill_order3_axis_values_first(xi[1], ly[q], dy[q]);
    }

    for (std::size_t node = 0; node < node_axes.size(); ++node) {
        const auto& axes = node_axes[node];
        write_quad_gradient_row_q4(
            gradients_out + node * 3u * output_stride,
            output_stride,
            lx,
            ly,
            dx,
            dy,
            axes[0],
            axes[1]);
    }
}

SVMP_LAGRANGE_NOINLINE void evaluate_quad_order4_gradients_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT gradients_out) {
    constexpr auto node_axes = detail::make_quad_tensor_node_axes<4>();

    Real lx[4][5];
    Real ly[4][5];
    Real dx[4][5];
    Real dy[4][5];
    for (std::size_t q = 0; q < 4u; ++q) {
        const auto& xi = points[q];
        fill_order4_axis_values_first(xi[0], lx[q], dx[q]);
        fill_order4_axis_values_first(xi[1], ly[q], dy[q]);
    }

    for (std::size_t node = 0; node < node_axes.size(); ++node) {
        const auto& axes = node_axes[node];
        write_quad_gradient_row_q4(
            gradients_out + node * 3u * output_stride,
            output_stride,
            lx,
            ly,
            dx,
            dy,
            axes[0],
            axes[1]);
    }
}

void evaluate_quad_order3_gradients_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT gradients_out) {
    if (points.size() == 4u) {
        evaluate_quad_order3_gradients_q4(points, output_stride, gradients_out);
        return;
    }

    Real* row0 = gradients_out + 0u * 3u * output_stride;
    Real* row1 = gradients_out + 1u * 3u * output_stride;
    Real* row2 = gradients_out + 2u * 3u * output_stride;
    Real* row3 = gradients_out + 3u * 3u * output_stride;
    Real* row4 = gradients_out + 4u * 3u * output_stride;
    Real* row5 = gradients_out + 5u * 3u * output_stride;
    Real* row6 = gradients_out + 6u * 3u * output_stride;
    Real* row7 = gradients_out + 7u * 3u * output_stride;
    Real* row8 = gradients_out + 8u * 3u * output_stride;
    Real* row9 = gradients_out + 9u * 3u * output_stride;
    Real* row10 = gradients_out + 10u * 3u * output_stride;
    Real* row11 = gradients_out + 11u * 3u * output_stride;
    Real* row12 = gradients_out + 12u * 3u * output_stride;
    Real* row13 = gradients_out + 13u * 3u * output_stride;
    Real* row14 = gradients_out + 14u * 3u * output_stride;
    Real* row15 = gradients_out + 15u * 3u * output_stride;

    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        Real lx[4];
        Real ly[4];
        Real dx[4];
        Real dy[4];
        fill_order3_axis_values_first(xi[0], lx, dx);
        fill_order3_axis_values_first(xi[1], ly, dy);
        write_quad_order2_gradient_q(row0, output_stride, q, dx[0] * ly[0], lx[0] * dy[0]);
        write_quad_order2_gradient_q(row1, output_stride, q, dx[1] * ly[0], lx[1] * dy[0]);
        write_quad_order2_gradient_q(row2, output_stride, q, dx[1] * ly[1], lx[1] * dy[1]);
        write_quad_order2_gradient_q(row3, output_stride, q, dx[0] * ly[1], lx[0] * dy[1]);
        write_quad_order2_gradient_q(row4, output_stride, q, dx[2] * ly[0], lx[2] * dy[0]);
        write_quad_order2_gradient_q(row5, output_stride, q, dx[3] * ly[0], lx[3] * dy[0]);
        write_quad_order2_gradient_q(row6, output_stride, q, dx[1] * ly[2], lx[1] * dy[2]);
        write_quad_order2_gradient_q(row7, output_stride, q, dx[1] * ly[3], lx[1] * dy[3]);
        write_quad_order2_gradient_q(row8, output_stride, q, dx[3] * ly[1], lx[3] * dy[1]);
        write_quad_order2_gradient_q(row9, output_stride, q, dx[2] * ly[1], lx[2] * dy[1]);
        write_quad_order2_gradient_q(row10, output_stride, q, dx[0] * ly[3], lx[0] * dy[3]);
        write_quad_order2_gradient_q(row11, output_stride, q, dx[0] * ly[2], lx[0] * dy[2]);
        write_quad_order2_gradient_q(row12, output_stride, q, dx[2] * ly[2], lx[2] * dy[2]);
        write_quad_order2_gradient_q(row13, output_stride, q, dx[3] * ly[2], lx[3] * dy[2]);
        write_quad_order2_gradient_q(row14, output_stride, q, dx[2] * ly[3], lx[2] * dy[3]);
        write_quad_order2_gradient_q(row15, output_stride, q, dx[3] * ly[3], lx[3] * dy[3]);
    }
}

void evaluate_quad_order3_hessians_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT hessians_out) {
    Real* row0 = hessians_out + 0u * 9u * output_stride;
    Real* row1 = hessians_out + 1u * 9u * output_stride;
    Real* row2 = hessians_out + 2u * 9u * output_stride;
    Real* row3 = hessians_out + 3u * 9u * output_stride;
    Real* row4 = hessians_out + 4u * 9u * output_stride;
    Real* row5 = hessians_out + 5u * 9u * output_stride;
    Real* row6 = hessians_out + 6u * 9u * output_stride;
    Real* row7 = hessians_out + 7u * 9u * output_stride;
    Real* row8 = hessians_out + 8u * 9u * output_stride;
    Real* row9 = hessians_out + 9u * 9u * output_stride;
    Real* row10 = hessians_out + 10u * 9u * output_stride;
    Real* row11 = hessians_out + 11u * 9u * output_stride;
    Real* row12 = hessians_out + 12u * 9u * output_stride;
    Real* row13 = hessians_out + 13u * 9u * output_stride;
    Real* row14 = hessians_out + 14u * 9u * output_stride;
    Real* row15 = hessians_out + 15u * 9u * output_stride;

    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        Real lx[4];
        Real ly[4];
        Real dx[4];
        Real dy[4];
        Real hx[4];
        Real hy[4];
        fill_order3_axis_values_first_second(xi[0], lx, dx, hx);
        fill_order3_axis_values_first_second(xi[1], ly, dy, hy);
        write_quad_order2_hessian_q(row0, output_stride, q, hx[0] * ly[0], dx[0] * dy[0], lx[0] * hy[0]);
        write_quad_order2_hessian_q(row1, output_stride, q, hx[1] * ly[0], dx[1] * dy[0], lx[1] * hy[0]);
        write_quad_order2_hessian_q(row2, output_stride, q, hx[1] * ly[1], dx[1] * dy[1], lx[1] * hy[1]);
        write_quad_order2_hessian_q(row3, output_stride, q, hx[0] * ly[1], dx[0] * dy[1], lx[0] * hy[1]);
        write_quad_order2_hessian_q(row4, output_stride, q, hx[2] * ly[0], dx[2] * dy[0], lx[2] * hy[0]);
        write_quad_order2_hessian_q(row5, output_stride, q, hx[3] * ly[0], dx[3] * dy[0], lx[3] * hy[0]);
        write_quad_order2_hessian_q(row6, output_stride, q, hx[1] * ly[2], dx[1] * dy[2], lx[1] * hy[2]);
        write_quad_order2_hessian_q(row7, output_stride, q, hx[1] * ly[3], dx[1] * dy[3], lx[1] * hy[3]);
        write_quad_order2_hessian_q(row8, output_stride, q, hx[3] * ly[1], dx[3] * dy[1], lx[3] * hy[1]);
        write_quad_order2_hessian_q(row9, output_stride, q, hx[2] * ly[1], dx[2] * dy[1], lx[2] * hy[1]);
        write_quad_order2_hessian_q(row10, output_stride, q, hx[0] * ly[3], dx[0] * dy[3], lx[0] * hy[3]);
        write_quad_order2_hessian_q(row11, output_stride, q, hx[0] * ly[2], dx[0] * dy[2], lx[0] * hy[2]);
        write_quad_order2_hessian_q(row12, output_stride, q, hx[2] * ly[2], dx[2] * dy[2], lx[2] * hy[2]);
        write_quad_order2_hessian_q(row13, output_stride, q, hx[3] * ly[2], dx[3] * dy[2], lx[3] * hy[2]);
        write_quad_order2_hessian_q(row14, output_stride, q, hx[2] * ly[3], dx[2] * dy[3], lx[2] * hy[3]);
        write_quad_order2_hessian_q(row15, output_stride, q, hx[3] * ly[3], dx[3] * dy[3], lx[3] * hy[3]);
    }
}

template <std::size_t Q>
inline void write_quad_order3_all_q4(
    std::size_t output_stride,
    std::size_t i,
    std::size_t j,
    const Real lx[4][4],
    const Real ly[4][4],
    const Real dx[4][4],
    const Real dy[4][4],
    const Real hx[4][4],
    const Real hy[4][4],
    Real* SVMP_RESTRICT value_row,
    Real* SVMP_RESTRICT grad_row,
    Real* SVMP_RESTRICT hess_row) {
    const Real xv = lx[Q][i];
    const Real yv = ly[Q][j];
    const Real xd = dx[Q][i];
    const Real yd = dy[Q][j];
    const Real hxy = xd * yd;

    value_row[Q] = xv * yv;
    grad_row[0u * output_stride + Q] = xd * yv;
    grad_row[1u * output_stride + Q] = xv * yd;
    grad_row[2u * output_stride + Q] = Real(0);
    hess_row[0u * output_stride + Q] = hx[Q][i] * yv;
    hess_row[4u * output_stride + Q] = xv * hy[Q][j];
    hess_row[8u * output_stride + Q] = Real(0);
    hess_row[1u * output_stride + Q] = hxy;
    hess_row[3u * output_stride + Q] = hxy;
    hess_row[2u * output_stride + Q] = Real(0);
    hess_row[6u * output_stride + Q] = Real(0);
    hess_row[5u * output_stride + Q] = Real(0);
    hess_row[7u * output_stride + Q] = Real(0);
}

void evaluate_quad_order3_all_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    constexpr auto node_axes = detail::make_quad_tensor_node_axes<3>();

    Real lx[4][4];
    Real ly[4][4];
    Real dx[4][4];
    Real dy[4][4];
    Real hx[4][4];
    Real hy[4][4];
    for (std::size_t q = 0; q < 4u; ++q) {
        const auto& xi = points[q];
        fill_order3_axis_values_first_second(xi[0], lx[q], dx[q], hx[q]);
        fill_order3_axis_values_first_second(xi[1], ly[q], dy[q], hy[q]);
    }

    for (std::size_t node = 0; node < node_axes.size(); ++node) {
        const auto& axes = node_axes[node];
        const std::size_t i = axes[0];
        const std::size_t j = axes[1];
        Real* value_row = values_out + node * output_stride;
        Real* grad_row = gradients_out + node * 3u * output_stride;
        Real* hess_row = hessians_out + node * 9u * output_stride;
        write_quad_order3_all_q4<0u>(
            output_stride, i, j, lx, ly, dx, dy, hx, hy, value_row, grad_row, hess_row);
        write_quad_order3_all_q4<1u>(
            output_stride, i, j, lx, ly, dx, dy, hx, hy, value_row, grad_row, hess_row);
        write_quad_order3_all_q4<2u>(
            output_stride, i, j, lx, ly, dx, dy, hx, hy, value_row, grad_row, hess_row);
        write_quad_order3_all_q4<3u>(
            output_stride, i, j, lx, ly, dx, dy, hx, hy, value_row, grad_row, hess_row);
    }
}

void evaluate_hex_order3_values_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    constexpr auto node_axes = detail::make_hex_tensor_node_axes<3>();

    Real lx[4][4];
    Real ly[4][4];
    Real lz[4][4];
    for (std::size_t q = 0; q < 4u; ++q) {
        const auto& xi = points[q];
        fill_order3_axis_values(xi[0], lx[q]);
        fill_order3_axis_values(xi[1], ly[q]);
        fill_order3_axis_values(xi[2], lz[q]);
    }

    for (std::size_t node = 0; node < node_axes.size(); ++node) {
        const auto& axes = node_axes[node];
        const std::size_t i = axes[0];
        const std::size_t j = axes[1];
        const std::size_t k = axes[2];
        Real* row = values_out + node * output_stride;
        row[0] = lx[0][i] * ly[0][j] * lz[0][k];
        row[1] = lx[1][i] * ly[1][j] * lz[1][k];
        row[2] = lx[2][i] * ly[2][j] * lz[2][k];
        row[3] = lx[3][i] * ly[3][j] * lz[3][k];
    }
}

void evaluate_hex_order3_gradients_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT gradients_out) {
    constexpr auto node_axes = detail::make_hex_tensor_node_axes<3>();

    Real lx[4][4];
    Real ly[4][4];
    Real lz[4][4];
    Real dx[4][4];
    Real dy[4][4];
    Real dz[4][4];
    for (std::size_t q = 0; q < 4u; ++q) {
        const auto& xi = points[q];
        fill_order3_axis_values_first(xi[0], lx[q], dx[q]);
        fill_order3_axis_values_first(xi[1], ly[q], dy[q]);
        fill_order3_axis_values_first(xi[2], lz[q], dz[q]);
    }

    for (std::size_t node = 0; node < node_axes.size(); ++node) {
        const auto& axes = node_axes[node];
        const std::size_t i = axes[0];
        const std::size_t j = axes[1];
        const std::size_t k = axes[2];
        Real* row = gradients_out + node * 3u * output_stride;
        row[0] = dx[0][i] * ly[0][j] * lz[0][k];
        row[1] = dx[1][i] * ly[1][j] * lz[1][k];
        row[2] = dx[2][i] * ly[2][j] * lz[2][k];
        row[3] = dx[3][i] * ly[3][j] * lz[3][k];
        row[output_stride + 0u] = lx[0][i] * dy[0][j] * lz[0][k];
        row[output_stride + 1u] = lx[1][i] * dy[1][j] * lz[1][k];
        row[output_stride + 2u] = lx[2][i] * dy[2][j] * lz[2][k];
        row[output_stride + 3u] = lx[3][i] * dy[3][j] * lz[3][k];
        row[2u * output_stride + 0u] = lx[0][i] * ly[0][j] * dz[0][k];
        row[2u * output_stride + 1u] = lx[1][i] * ly[1][j] * dz[1][k];
        row[2u * output_stride + 2u] = lx[2][i] * ly[2][j] * dz[2][k];
        row[2u * output_stride + 3u] = lx[3][i] * ly[3][j] * dz[3][k];
    }
}

template <std::size_t Q, bool WriteValue, bool WriteGradient>
inline void write_hex_order3_q4_hessian_outputs(
    std::size_t output_stride,
    std::size_t i,
    std::size_t j,
    std::size_t k,
    const Real lx[4][4],
    const Real ly[4][4],
    const Real lz[4][4],
    const Real dx[4][4],
    const Real dy[4][4],
    const Real dz[4][4],
    const Real hx[4][4],
    const Real hy[4][4],
    const Real hz[4][4],
    Real* SVMP_RESTRICT value_row,
    Real* SVMP_RESTRICT grad_row,
    Real* SVMP_RESTRICT hess_row) {
    const Real xv = lx[Q][i];
    const Real yv = ly[Q][j];
    const Real zv = lz[Q][k];
    const Real yz = yv * zv;

    if constexpr (WriteValue) {
        value_row[Q] = xv * yz;
    }

    const Real xd = dx[Q][i];
    const Real yd = dy[Q][j];
    const Real zd = dz[Q][k];
    const Real yd_z = yd * zv;
    const Real yv_zd = yv * zd;

    if constexpr (WriteGradient) {
        grad_row[0u * output_stride + Q] = xd * yz;
        grad_row[1u * output_stride + Q] = xv * yd_z;
        grad_row[2u * output_stride + Q] = xv * yv_zd;
    }

    const Real hxy = xd * yd_z;
    const Real hxz = xd * yv_zd;
    const Real hyz = xv * yd * zd;
    hess_row[0u * output_stride + Q] = hx[Q][i] * yz;
    hess_row[4u * output_stride + Q] = xv * hy[Q][j] * zv;
    hess_row[8u * output_stride + Q] = xv * yv * hz[Q][k];
    hess_row[1u * output_stride + Q] = hxy;
    hess_row[3u * output_stride + Q] = hxy;
    hess_row[2u * output_stride + Q] = hxz;
    hess_row[6u * output_stride + Q] = hxz;
    hess_row[5u * output_stride + Q] = hyz;
    hess_row[7u * output_stride + Q] = hyz;
}

template <bool WriteValue, bool WriteGradient>
void evaluate_hex_order3_q4_hessian_outputs(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    constexpr auto node_axes = detail::make_hex_tensor_node_axes<3>();

    Real lx[4][4];
    Real ly[4][4];
    Real lz[4][4];
    Real dx[4][4];
    Real dy[4][4];
    Real dz[4][4];
    Real hx[4][4];
    Real hy[4][4];
    Real hz[4][4];
    for (std::size_t q = 0; q < 4u; ++q) {
        const auto& xi = points[q];
        fill_order3_axis_values_first_second(xi[0], lx[q], dx[q], hx[q]);
        fill_order3_axis_values_first_second(xi[1], ly[q], dy[q], hy[q]);
        fill_order3_axis_values_first_second(xi[2], lz[q], dz[q], hz[q]);
    }

    for (std::size_t node = 0; node < node_axes.size(); ++node) {
        const auto& axes = node_axes[node];
        const std::size_t i = axes[0];
        const std::size_t j = axes[1];
        const std::size_t k = axes[2];
        Real* value_row = values_out ? values_out + node * output_stride : nullptr;
        Real* grad_row = gradients_out ? gradients_out + node * 3u * output_stride : nullptr;
        Real* hess_row = hessians_out + node * 9u * output_stride;
        write_hex_order3_q4_hessian_outputs<0u, WriteValue, WriteGradient>(
            output_stride, i, j, k, lx, ly, lz, dx, dy, dz, hx, hy, hz,
            value_row, grad_row, hess_row);
        write_hex_order3_q4_hessian_outputs<1u, WriteValue, WriteGradient>(
            output_stride, i, j, k, lx, ly, lz, dx, dy, dz, hx, hy, hz,
            value_row, grad_row, hess_row);
        write_hex_order3_q4_hessian_outputs<2u, WriteValue, WriteGradient>(
            output_stride, i, j, k, lx, ly, lz, dx, dy, dz, hx, hy, hz,
            value_row, grad_row, hess_row);
        write_hex_order3_q4_hessian_outputs<3u, WriteValue, WriteGradient>(
            output_stride, i, j, k, lx, ly, lz, dx, dy, dz, hx, hy, hz,
            value_row, grad_row, hess_row);
    }
}

void evaluate_hex_order3_hessians_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT hessians_out) {
    evaluate_hex_order3_q4_hessian_outputs<false, false>(
        points, output_stride, nullptr, nullptr, hessians_out);
}

void evaluate_hex_order3_all_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    evaluate_hex_order3_q4_hessian_outputs<true, true>(
        points, output_stride, values_out, gradients_out, hessians_out);
}

void evaluate_hex_order2_values_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    Real* row0 = values_out + 0u * output_stride;
    Real* row1 = values_out + 1u * output_stride;
    Real* row2 = values_out + 2u * output_stride;
    Real* row3 = values_out + 3u * output_stride;
    Real* row4 = values_out + 4u * output_stride;
    Real* row5 = values_out + 5u * output_stride;
    Real* row6 = values_out + 6u * output_stride;
    Real* row7 = values_out + 7u * output_stride;
    Real* row8 = values_out + 8u * output_stride;
    Real* row9 = values_out + 9u * output_stride;
    Real* row10 = values_out + 10u * output_stride;
    Real* row11 = values_out + 11u * output_stride;
    Real* row12 = values_out + 12u * output_stride;
    Real* row13 = values_out + 13u * output_stride;
    Real* row14 = values_out + 14u * output_stride;
    Real* row15 = values_out + 15u * output_stride;
    Real* row16 = values_out + 16u * output_stride;
    Real* row17 = values_out + 17u * output_stride;
    Real* row18 = values_out + 18u * output_stride;
    Real* row19 = values_out + 19u * output_stride;
    Real* row20 = values_out + 20u * output_stride;
    Real* row21 = values_out + 21u * output_stride;
    Real* row22 = values_out + 22u * output_stride;
    Real* row23 = values_out + 23u * output_stride;
    Real* row24 = values_out + 24u * output_stride;
    Real* row25 = values_out + 25u * output_stride;
    Real* row26 = values_out + 26u * output_stride;

    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        const Real x = xi[0];
        const Real y = xi[1];
        const Real z = xi[2];
        const Real x0 = x * (x - Real(1)) * Real(0.5);
        const Real x1 = x * (x + Real(1)) * Real(0.5);
        const Real x2 = Real(1) - x * x;
        const Real y0 = y * (y - Real(1)) * Real(0.5);
        const Real y1 = y * (y + Real(1)) * Real(0.5);
        const Real y2 = Real(1) - y * y;
        const Real z0 = z * (z - Real(1)) * Real(0.5);
        const Real z1 = z * (z + Real(1)) * Real(0.5);
        const Real z2 = Real(1) - z * z;
        const Real x0y0 = x0 * y0;
        const Real x1y0 = x1 * y0;
        const Real x1y1 = x1 * y1;
        const Real x0y1 = x0 * y1;
        const Real x2y0 = x2 * y0;
        const Real x1y2 = x1 * y2;
        const Real x2y1 = x2 * y1;
        const Real x0y2 = x0 * y2;
        const Real x2y2 = x2 * y2;

        row0[q] = x0y0 * z0;
        row1[q] = x1y0 * z0;
        row2[q] = x1y1 * z0;
        row3[q] = x0y1 * z0;
        row4[q] = x0y0 * z1;
        row5[q] = x1y0 * z1;
        row6[q] = x1y1 * z1;
        row7[q] = x0y1 * z1;
        row8[q] = x2y0 * z0;
        row9[q] = x1y2 * z0;
        row10[q] = x2y1 * z0;
        row11[q] = x0y2 * z0;
        row12[q] = x2y0 * z1;
        row13[q] = x1y2 * z1;
        row14[q] = x2y1 * z1;
        row15[q] = x0y2 * z1;
        row16[q] = x0y0 * z2;
        row17[q] = x1y0 * z2;
        row18[q] = x1y1 * z2;
        row19[q] = x0y1 * z2;
        row20[q] = x2y2 * z0;
        row21[q] = x2y2 * z1;
        row22[q] = x2y0 * z2;
        row23[q] = x1y2 * z2;
        row24[q] = x2y1 * z2;
        row25[q] = x0y2 * z2;
        row26[q] = x2y2 * z2;
    }
}

inline void fill_order2_axis_values_first(Real x,
                                          Real* SVMP_RESTRICT values,
                                          Real* SVMP_RESTRICT first) {
    values[0] = x * (x - Real(1)) * Real(0.5);
    values[1] = x * (x + Real(1)) * Real(0.5);
    values[2] = Real(1) - x * x;
    first[0] = x - Real(0.5);
    first[1] = x + Real(0.5);
    first[2] = Real(-2) * x;
}

inline void fill_order2_axis_values_first_second(Real x,
                                                 Real* SVMP_RESTRICT values,
                                                 Real* SVMP_RESTRICT first,
                                                 Real* SVMP_RESTRICT second) {
    fill_order2_axis_values_first(x, values, first);
    second[0] = Real(1);
    second[1] = Real(1);
    second[2] = Real(-2);
}

template <std::size_t Q>
inline void write_hex_order2_hessian_q4(
    std::size_t output_stride,
    std::size_t i,
    std::size_t j,
    std::size_t k,
    const Real lx[4][3],
    const Real ly[4][3],
    const Real lz[4][3],
    const Real dx[4][3],
    const Real dy[4][3],
    const Real dz[4][3],
    const Real hx[4][3],
    const Real hy[4][3],
    const Real hz[4][3],
    Real* SVMP_RESTRICT hess_row) {
    const Real xv = lx[Q][i];
    const Real yv = ly[Q][j];
    const Real zv = lz[Q][k];
    const Real yz = yv * zv;
    const Real xd = dx[Q][i];
    const Real yd = dy[Q][j];
    const Real zd = dz[Q][k];
    const Real yd_z = yd * zv;
    const Real yv_zd = yv * zd;
    const Real hxy = xd * yd_z;
    const Real hxz = xd * yv_zd;
    const Real hyz = xv * yd * zd;
    hess_row[0u * output_stride + Q] = hx[Q][i] * yz;
    hess_row[4u * output_stride + Q] = xv * hy[Q][j] * zv;
    hess_row[8u * output_stride + Q] = xv * yv * hz[Q][k];
    hess_row[1u * output_stride + Q] = hxy;
    hess_row[3u * output_stride + Q] = hxy;
    hess_row[2u * output_stride + Q] = hxz;
    hess_row[6u * output_stride + Q] = hxz;
    hess_row[5u * output_stride + Q] = hyz;
    hess_row[7u * output_stride + Q] = hyz;
}

void evaluate_hex_order2_hessians_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT hessians_out) {
    constexpr auto node_axes = detail::make_hex_tensor_node_axes<2>();

    Real lx[4][3];
    Real ly[4][3];
    Real lz[4][3];
    Real dx[4][3];
    Real dy[4][3];
    Real dz[4][3];
    Real hx[4][3];
    Real hy[4][3];
    Real hz[4][3];
    for (std::size_t q = 0; q < 4u; ++q) {
        const auto& xi = points[q];
        fill_order2_axis_values_first_second(xi[0], lx[q], dx[q], hx[q]);
        fill_order2_axis_values_first_second(xi[1], ly[q], dy[q], hy[q]);
        fill_order2_axis_values_first_second(xi[2], lz[q], dz[q], hz[q]);
    }

    for (std::size_t node = 0; node < node_axes.size(); ++node) {
        const auto& axes = node_axes[node];
        const std::size_t i = axes[0];
        const std::size_t j = axes[1];
        const std::size_t k = axes[2];
        Real* hess_row = hessians_out + node * 9u * output_stride;
        write_hex_order2_hessian_q4<0u>(
            output_stride, i, j, k, lx, ly, lz, dx, dy, dz, hx, hy, hz, hess_row);
        write_hex_order2_hessian_q4<1u>(
            output_stride, i, j, k, lx, ly, lz, dx, dy, dz, hx, hy, hz, hess_row);
        write_hex_order2_hessian_q4<2u>(
            output_stride, i, j, k, lx, ly, lz, dx, dy, dz, hx, hy, hz, hess_row);
        write_hex_order2_hessian_q4<3u>(
            output_stride, i, j, k, lx, ly, lz, dx, dy, dz, hx, hy, hz, hess_row);
    }
}

template <std::size_t Q>
inline void write_quad_order2_all_q4(
    std::size_t output_stride,
    std::size_t i,
    std::size_t j,
    const Real lx[4][3],
    const Real ly[4][3],
    const Real dx[4][3],
    const Real dy[4][3],
    const Real hx[4][3],
    const Real hy[4][3],
    Real* SVMP_RESTRICT value_row,
    Real* SVMP_RESTRICT grad_row,
    Real* SVMP_RESTRICT hess_row) {
    const Real xv = lx[Q][i];
    const Real yv = ly[Q][j];
    const Real xd = dx[Q][i];
    const Real yd = dy[Q][j];
    const Real hxy = xd * yd;

    value_row[Q] = xv * yv;
    grad_row[0u * output_stride + Q] = xd * yv;
    grad_row[1u * output_stride + Q] = xv * yd;
    grad_row[2u * output_stride + Q] = Real(0);
    hess_row[0u * output_stride + Q] = hx[Q][i] * yv;
    hess_row[4u * output_stride + Q] = xv * hy[Q][j];
    hess_row[8u * output_stride + Q] = Real(0);
    hess_row[1u * output_stride + Q] = hxy;
    hess_row[3u * output_stride + Q] = hxy;
    hess_row[2u * output_stride + Q] = Real(0);
    hess_row[6u * output_stride + Q] = Real(0);
    hess_row[5u * output_stride + Q] = Real(0);
    hess_row[7u * output_stride + Q] = Real(0);
}

void evaluate_quad_order2_all_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    constexpr auto node_axes = detail::make_quad_tensor_node_axes<2>();

    Real lx[4][3];
    Real ly[4][3];
    Real dx[4][3];
    Real dy[4][3];
    Real hx[4][3];
    Real hy[4][3];
    for (std::size_t q = 0; q < 4u; ++q) {
        const auto& xi = points[q];
        fill_order2_axis_values_first_second(xi[0], lx[q], dx[q], hx[q]);
        fill_order2_axis_values_first_second(xi[1], ly[q], dy[q], hy[q]);
    }

    for (std::size_t node = 0; node < node_axes.size(); ++node) {
        const auto& axes = node_axes[node];
        const std::size_t i = axes[0];
        const std::size_t j = axes[1];
        Real* value_row = values_out + node * output_stride;
        Real* grad_row = gradients_out + node * 3u * output_stride;
        Real* hess_row = hessians_out + node * 9u * output_stride;
        write_quad_order2_all_q4<0u>(
            output_stride, i, j, lx, ly, dx, dy, hx, hy, value_row, grad_row, hess_row);
        write_quad_order2_all_q4<1u>(
            output_stride, i, j, lx, ly, dx, dy, hx, hy, value_row, grad_row, hess_row);
        write_quad_order2_all_q4<2u>(
            output_stride, i, j, lx, ly, dx, dy, hx, hy, value_row, grad_row, hess_row);
        write_quad_order2_all_q4<3u>(
            output_stride, i, j, lx, ly, dx, dy, hx, hy, value_row, grad_row, hess_row);
    }
}

void evaluate_hex_order2_gradients_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT gradients_out) {
    constexpr std::array<std::array<std::size_t, 3>, 27> node_axes = {{
        {{0u, 0u, 0u}}, {{1u, 0u, 0u}}, {{1u, 1u, 0u}}, {{0u, 1u, 0u}},
        {{0u, 0u, 1u}}, {{1u, 0u, 1u}}, {{1u, 1u, 1u}}, {{0u, 1u, 1u}},
        {{2u, 0u, 0u}}, {{1u, 2u, 0u}}, {{2u, 1u, 0u}}, {{0u, 2u, 0u}},
        {{2u, 0u, 1u}}, {{1u, 2u, 1u}}, {{2u, 1u, 1u}}, {{0u, 2u, 1u}},
        {{0u, 0u, 2u}}, {{1u, 0u, 2u}}, {{1u, 1u, 2u}}, {{0u, 1u, 2u}},
        {{2u, 2u, 0u}}, {{2u, 2u, 1u}}, {{2u, 0u, 2u}}, {{1u, 2u, 2u}},
        {{2u, 1u, 2u}}, {{0u, 2u, 2u}}, {{2u, 2u, 2u}},
    }};

    Real lx[4][3];
    Real ly[4][3];
    Real lz[4][3];
    Real dx[4][3];
    Real dy[4][3];
    Real dz[4][3];
    for (std::size_t q = 0; q < 4u; ++q) {
        const auto& xi = points[q];
        fill_order2_axis_values_first(xi[0], lx[q], dx[q]);
        fill_order2_axis_values_first(xi[1], ly[q], dy[q]);
        fill_order2_axis_values_first(xi[2], lz[q], dz[q]);
    }

    for (std::size_t node = 0; node < node_axes.size(); ++node) {
        const auto& axes = node_axes[node];
        const std::size_t i = axes[0];
        const std::size_t j = axes[1];
        const std::size_t k = axes[2];
        Real* row = gradients_out + node * 3u * output_stride;
        row[0] = dx[0][i] * ly[0][j] * lz[0][k];
        row[1] = dx[1][i] * ly[1][j] * lz[1][k];
        row[2] = dx[2][i] * ly[2][j] * lz[2][k];
        row[3] = dx[3][i] * ly[3][j] * lz[3][k];
        row[output_stride + 0u] = lx[0][i] * dy[0][j] * lz[0][k];
        row[output_stride + 1u] = lx[1][i] * dy[1][j] * lz[1][k];
        row[output_stride + 2u] = lx[2][i] * dy[2][j] * lz[2][k];
        row[output_stride + 3u] = lx[3][i] * dy[3][j] * lz[3][k];
        row[2u * output_stride + 0u] = lx[0][i] * ly[0][j] * dz[0][k];
        row[2u * output_stride + 1u] = lx[1][i] * ly[1][j] * dz[1][k];
        row[2u * output_stride + 2u] = lx[2][i] * ly[2][j] * dz[2][k];
        row[2u * output_stride + 3u] = lx[3][i] * ly[3][j] * dz[3][k];
    }
}

template<typename FastBasis>
void evaluate_constant_fast_hessians_strided(
    std::size_t num_qpts,
    std::size_t output_stride,
    Real* SVMP_RESTRICT hessians_out) {
    std::array<Hessian, FastBasis::n_dofs> fast_hessians{};
    FastBasis::evaluate_hessians(math::Vector<Real, 3>{}, fast_hessians);
    for (std::size_t i = 0; i < fast_hessians.size(); ++i) {
        const Hessian& hessian = fast_hessians[i];
        Real* H = hessians_out + i * 9u * output_stride;
        const Real h00 = hessian(0, 0);
        const Real h01 = hessian(0, 1);
        const Real h02 = hessian(0, 2);
        const Real h10 = hessian(1, 0);
        const Real h11 = hessian(1, 1);
        const Real h12 = hessian(1, 2);
        const Real h20 = hessian(2, 0);
        const Real h21 = hessian(2, 1);
        const Real h22 = hessian(2, 2);
        for (std::size_t q = 0; q < num_qpts; ++q) {
            H[0u * output_stride + q] = h00;
            H[1u * output_stride + q] = h01;
            H[2u * output_stride + q] = h02;
            H[3u * output_stride + q] = h10;
            H[4u * output_stride + q] = h11;
            H[5u * output_stride + q] = h12;
            H[6u * output_stride + q] = h20;
            H[7u * output_stride + q] = h21;
            H[8u * output_stride + q] = h22;
        }
    }
}

template<typename FastBasis>
void evaluate_fast_outputs_with_constant_hessians_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        if (values_out != nullptr) {
            std::array<Real, FastBasis::n_dofs> fast_values{};
            FastBasis::evaluate(xi, fast_values);
            for (std::size_t i = 0; i < fast_values.size(); ++i) {
                values_out[i * output_stride + q] = fast_values[i];
            }
        }
        if (gradients_out != nullptr) {
            std::array<Gradient, FastBasis::n_dofs> fast_gradients{};
            FastBasis::evaluate_gradients(xi, fast_gradients);
            for (std::size_t i = 0; i < fast_gradients.size(); ++i) {
                Real* g = gradients_out + i * 3u * output_stride;
                g[0u * output_stride + q] = fast_gradients[i][0];
                g[1u * output_stride + q] = fast_gradients[i][1];
                g[2u * output_stride + q] = fast_gradients[i][2];
            }
        }
    }
    evaluate_constant_fast_hessians_strided<FastBasis>(
        points.size(), output_stride, hessians_out);
}

template<int Order>
void evaluate_wedge_fast_outputs_strided(
    const std::vector<std::array<std::size_t, 2>>& wedge_indices,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    static_assert(Order >= 1 && Order <= 2,
                  "wedge fast outputs rely on low-order public triangle ordering");
    using TriFast = LagrangeTriFast<Order>;
    constexpr std::size_t axis_size = static_cast<std::size_t>(Order + 1);
    const bool need_grad = gradients_out != nullptr;
    const bool need_hess = hessians_out != nullptr;

    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        std::array<Real, TriFast::n_dofs> tri_values{};
        std::array<Gradient, TriFast::n_dofs> tri_gradients{};
        std::array<Hessian, TriFast::n_dofs> tri_hessians{};
        std::array<Real, axis_size> z_values{};
        std::array<Real, axis_size> z_first{};
        std::array<Real, axis_size> z_second{};

        TriFast::evaluate(xi, tri_values);
        if (need_grad || need_hess) {
            TriFast::evaluate_gradients(xi, tri_gradients);
        }
        if (need_hess) {
            TriFast::evaluate_hessians(xi, tri_hessians);
            detail::fill_axis_values_first_second<Order>(xi[2], z_values, z_first, z_second);
        } else if (need_grad) {
            detail::fill_axis_values_first<Order>(xi[2], z_values, z_first);
        } else {
            detail::fill_axis_values<Order>(xi[2], z_values);
        }

        for (std::size_t node = 0; node < wedge_indices.size(); ++node) {
            const auto& index = wedge_indices[node];
            const std::size_t tri = index[0];
            const std::size_t z = index[1];
            const Real tri_v = tri_values[tri];
            const Real zv = z_values[z];

            if (values_out != nullptr) {
                values_out[node * output_stride + q] = tri_v * zv;
            }

            if (gradients_out != nullptr) {
                Real* g = gradients_out + node * 3u * output_stride;
                const Gradient& tri_g = tri_gradients[tri];
                g[0u * output_stride + q] = tri_g[0] * zv;
                g[1u * output_stride + q] = tri_g[1] * zv;
                g[2u * output_stride + q] = tri_v * z_first[z];
            }

            if (hessians_out != nullptr) {
                Real* H = hessians_out + node * 9u * output_stride;
                const Gradient& tri_g = tri_gradients[tri];
                const Hessian& tri_H = tri_hessians[tri];
                const Real zd = z_first[z];
                const Real hxz = tri_g[0] * zd;
                const Real hxy = tri_H(0, 1) * zv;
                const Real hyz = tri_g[1] * zd;
                H[0u * output_stride + q] = tri_H(0, 0) * zv;
                H[1u * output_stride + q] = hxy;
                H[2u * output_stride + q] = hxz;
                H[3u * output_stride + q] = hxy;
                H[4u * output_stride + q] = tri_H(1, 1) * zv;
                H[5u * output_stride + q] = hyz;
                H[6u * output_stride + q] = hxz;
                H[7u * output_stride + q] = hyz;
                H[8u * output_stride + q] = tri_v * z_second[z];
            }
        }
    }
}

template <int Order>
inline void fill_triangle_simplex_product_factors(Real lambda, Real* SVMP_RESTRICT factors) {
    const Real t = static_cast<Real>(Order) * lambda;
    factors[0] = Real(1);
    for (int a = 1; a <= Order; ++a) {
        factors[a] =
            factors[a - 1] *
            (t - static_cast<Real>(a - 1)) /
            static_cast<Real>(a);
    }
}

template <int Order>
SVMP_LAGRANGE_NOINLINE SVMP_LAGRANGE_ALIGN64 bool evaluate_wedge_values_product_q4(
    const std::vector<std::array<int, 4>>& simplex_exponents,
    const std::vector<std::array<std::size_t, 2>>& wedge_indices,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    constexpr std::size_t tri_count =
        static_cast<std::size_t>((Order + 1) * (Order + 2) / 2);
    if (simplex_exponents.size() != tri_count || points.size() != 4u) {
        return false;
    }

    Real tri_values[4][tri_count];
    std::array<Real, Order + 1> z_values[4];
    for (std::size_t q = 0; q < 4u; ++q) {
        const auto& xi = points[q];
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l0 = Real(1) - l1 - l2;
        Real f0[Order + 1];
        Real f1[Order + 1];
        Real f2[Order + 1];
        fill_triangle_simplex_product_factors<Order>(l0, f0);
        fill_triangle_simplex_product_factors<Order>(l1, f1);
        fill_triangle_simplex_product_factors<Order>(l2, f2);
        detail::fill_axis_values<Order>(xi[2], z_values[q]);

        for (std::size_t tri = 0; tri < tri_count; ++tri) {
            const auto& e = simplex_exponents[tri];
            tri_values[q][tri] =
                f0[static_cast<std::size_t>(e[0])] *
                f1[static_cast<std::size_t>(e[1])] *
                f2[static_cast<std::size_t>(e[2])];
        }
    }

    for (std::size_t node = 0; node < wedge_indices.size(); ++node) {
        const auto& index = wedge_indices[node];
        const std::size_t tri = index[0];
        const std::size_t z = index[1];
        Real* SVMP_RESTRICT row = values_out + node * output_stride;
        row[0] = tri_values[0][tri] * z_values[0][z];
        row[1] = tri_values[1][tri] * z_values[1][z];
        row[2] = tri_values[2][tri] * z_values[2][z];
        row[3] = tri_values[3][tri] * z_values[3][z];
    }
    return true;
}

bool try_evaluate_wedge_values_product_q4(
    const std::vector<std::array<int, 4>>& simplex_exponents,
    const std::vector<std::array<std::size_t, 2>>& wedge_indices,
    int order,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    switch (order) {
        case 4:
            return evaluate_wedge_values_product_q4<4>(
                simplex_exponents, wedge_indices, points, output_stride, values_out);
        case 5:
            return evaluate_wedge_values_product_q4<5>(
                simplex_exponents, wedge_indices, points, output_stride, values_out);
        case 6:
            return evaluate_wedge_values_product_q4<6>(
                simplex_exponents, wedge_indices, points, output_stride, values_out);
        case 7:
            return evaluate_wedge_values_product_q4<7>(
                simplex_exponents, wedge_indices, points, output_stride, values_out);
        case 8:
            return evaluate_wedge_values_product_q4<8>(
                simplex_exponents, wedge_indices, points, output_stride, values_out);
        default:
            return false;
    }
}

void evaluate_wedge_order1_values_q4(
    const std::vector<std::array<std::size_t, 2>>& wedge_indices,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    Real tri[4][3];
    Real axis[4][2];
    for (std::size_t q = 0; q < 4u; ++q) {
        const auto& xi = points[q];
        tri[q][0] = Real(1) - xi[0] - xi[1];
        tri[q][1] = xi[0];
        tri[q][2] = xi[1];
        axis[q][0] = (Real(1) - xi[2]) * Real(0.5);
        axis[q][1] = (Real(1) + xi[2]) * Real(0.5);
    }

    for (std::size_t node = 0; node < wedge_indices.size(); ++node) {
        const auto& index = wedge_indices[node];
        const std::size_t tri_node = index[0];
        const std::size_t axis_node = index[1];
        Real* row = values_out + node * output_stride;
        row[0] = tri[0][tri_node] * axis[0][axis_node];
        row[1] = tri[1][tri_node] * axis[1][axis_node];
        row[2] = tri[2][tri_node] * axis[2][axis_node];
        row[3] = tri[3][tri_node] * axis[3][axis_node];
    }
}

bool evaluate_wedge_fast_strided(
    int order,
    const std::vector<std::array<std::size_t, 2>>& wedge_indices,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    if (order == 3) {
        return false;
    }
    if (order == 1 &&
        points.size() == 4u &&
        values_out != nullptr &&
        gradients_out == nullptr &&
        hessians_out == nullptr) {
        evaluate_wedge_order1_values_q4(wedge_indices, points, output_stride, values_out);
        return true;
    }

    switch (order) {
        case 1:
            evaluate_wedge_fast_outputs_strided<1>(
                wedge_indices, points, output_stride, values_out, gradients_out, hessians_out);
            return true;
        case 2:
            evaluate_wedge_fast_outputs_strided<2>(
                wedge_indices, points, output_stride, values_out, gradients_out, hessians_out);
            return true;
        default:
            return false;
    }
}

bool evaluate_fixed_lagrange_fast(LagrangeTopology topology,
                                  int order,
                                  const math::Vector<Real, 3>& xi,
                                  std::vector<Real>* values,
                                  std::vector<Gradient>* gradients,
                                  std::vector<Hessian>* hessians) {
    switch (order) {
        case 1:
            return evaluate_fixed_lagrange_fast_order<1>(
                topology, xi, values, gradients, hessians);
        case 2:
            return evaluate_fixed_lagrange_fast_order<2>(
                topology, xi, values, gradients, hessians);
        case 3:
            return evaluate_fixed_lagrange_fast_order<3>(
                topology, xi, values, gradients, hessians);
        default:
            return false;
    }
}

bool evaluate_fixed_lagrange_fast_strided(LagrangeTopology topology,
                                          int order,
                                          const std::vector<math::Vector<Real, 3>>& points,
                                          std::size_t output_stride,
                                          Real* SVMP_RESTRICT values_out,
                                          Real* SVMP_RESTRICT gradients_out,
                                          Real* SVMP_RESTRICT hessians_out) {
    if (topology == LagrangeTopology::Line &&
        points.size() == 4u) {
        const bool values_only =
            values_out != nullptr && gradients_out == nullptr && hessians_out == nullptr;
        const bool gradients_only =
            values_out == nullptr && gradients_out != nullptr && hessians_out == nullptr;
        const bool hessians_only =
            values_out == nullptr && gradients_out == nullptr && hessians_out != nullptr;
        const bool all_outputs =
            values_out != nullptr && gradients_out != nullptr && hessians_out != nullptr;
        if (values_only) {
            if (order == 1) {
                evaluate_line_order1_values_q4(points, output_stride, values_out);
                return true;
            }
            if (order == 2) {
                evaluate_line_order2_values_q4(points, output_stride, values_out);
                return true;
            }
            if (order == 3) {
                evaluate_line_order3_values_q4(points, output_stride, values_out);
                return true;
            }
        }
        if (order == 1) {
            if (gradients_only) {
                evaluate_line_order1_gradients_q4(output_stride, gradients_out);
                return true;
            }
            if (hessians_only) {
                evaluate_line_order1_hessians_q4(output_stride, hessians_out);
                return true;
            }
            if (all_outputs) {
                evaluate_line_order1_all_q4(
                    points, output_stride, values_out, gradients_out, hessians_out);
                return true;
            }
        }
        if (order == 2) {
            if (gradients_only) {
                evaluate_line_order2_gradients_q4(points, output_stride, gradients_out);
                return true;
            }
            if (hessians_only) {
                evaluate_line_order2_hessians_q4(output_stride, hessians_out);
                return true;
            }
            if (all_outputs) {
                evaluate_line_order2_all_q4(
                    points, output_stride, values_out, gradients_out, hessians_out);
                return true;
            }
        }
        if (order == 3) {
            if (gradients_only) {
                evaluate_line_order3_gradients_q4(points, output_stride, gradients_out);
                return true;
            }
            if (hessians_only) {
                evaluate_line_order3_hessians_q4(points, output_stride, hessians_out);
                return true;
            }
            if (all_outputs) {
                evaluate_line_order3_all_q4(
                    points, output_stride, values_out, gradients_out, hessians_out);
                return true;
            }
        }
    }

    if (topology == LagrangeTopology::Tetrahedron &&
        order == 3 &&
        (gradients_out != nullptr || hessians_out != nullptr)) {
        return false;
    }
    if (topology == LagrangeTopology::Triangle &&
        order == 3 &&
        hessians_out != nullptr) {
        return false;
    }
    if (topology == LagrangeTopology::Triangle &&
        order == 1 &&
        values_out != nullptr &&
        gradients_out == nullptr &&
        hessians_out == nullptr) {
        evaluate_triangle_order1_values_strided(points, output_stride, values_out);
        return true;
    }
    if (topology == LagrangeTopology::Triangle &&
        order == 1 &&
        values_out == nullptr &&
        gradients_out != nullptr &&
        hessians_out == nullptr) {
        evaluate_triangle_order1_gradients_strided(points.size(), output_stride, gradients_out);
        return true;
    }
    if (topology == LagrangeTopology::Tetrahedron &&
        order == 1 &&
        values_out != nullptr &&
        gradients_out == nullptr &&
        hessians_out == nullptr) {
        evaluate_tet_order1_values_strided(points, output_stride, values_out);
        return true;
    }
    if (topology == LagrangeTopology::Tetrahedron &&
        order == 1 &&
        values_out == nullptr &&
        gradients_out != nullptr &&
        hessians_out == nullptr) {
        evaluate_tet_order1_gradients_strided(points.size(), output_stride, gradients_out);
        return true;
    }
    if (topology == LagrangeTopology::Triangle &&
        order == 1 &&
        values_out == nullptr &&
        gradients_out == nullptr &&
        hessians_out != nullptr) {
        evaluate_zero_hessians_strided(3u, points.size(), output_stride, hessians_out);
        return true;
    }
    if (topology == LagrangeTopology::Tetrahedron &&
        order == 1 &&
        values_out == nullptr &&
        gradients_out == nullptr &&
        hessians_out != nullptr) {
        evaluate_zero_hessians_strided(4u, points.size(), output_stride, hessians_out);
        return true;
    }
    if (topology == LagrangeTopology::Triangle &&
        order == 1 &&
        values_out != nullptr &&
        gradients_out != nullptr &&
        hessians_out != nullptr) {
        evaluate_triangle_order1_values_strided(points, output_stride, values_out);
        evaluate_triangle_order1_gradients_strided(points.size(), output_stride, gradients_out);
        evaluate_zero_hessians_strided(3u, points.size(), output_stride, hessians_out);
        return true;
    }
    if (topology == LagrangeTopology::Tetrahedron &&
        order == 1 &&
        values_out != nullptr &&
        gradients_out != nullptr &&
        hessians_out != nullptr) {
        evaluate_tet_order1_values_strided(points, output_stride, values_out);
        evaluate_tet_order1_gradients_strided(points.size(), output_stride, gradients_out);
        evaluate_zero_hessians_strided(4u, points.size(), output_stride, hessians_out);
        return true;
    }
    if (topology == LagrangeTopology::Triangle &&
        order == 2 &&
        values_out != nullptr &&
        gradients_out == nullptr &&
        hessians_out == nullptr) {
        evaluate_triangle_order2_values_strided(points, output_stride, values_out);
        return true;
    }
    if (topology == LagrangeTopology::Triangle &&
        order == 2 &&
        values_out == nullptr &&
        gradients_out != nullptr &&
        hessians_out == nullptr) {
        evaluate_triangle_order2_gradients_strided(points, output_stride, gradients_out);
        return true;
    }
    if (topology == LagrangeTopology::Triangle &&
        order == 2 &&
        points.size() == 4u &&
        values_out == nullptr &&
        gradients_out == nullptr &&
        hessians_out != nullptr) {
        evaluate_triangle_order2_hessians_q4(output_stride, hessians_out);
        return true;
    }
    if (topology == LagrangeTopology::Triangle &&
        order == 2 &&
        points.size() == 4u &&
        values_out != nullptr &&
        gradients_out != nullptr &&
        hessians_out != nullptr) {
        evaluate_triangle_order2_values_strided(points, output_stride, values_out);
        evaluate_triangle_order2_gradients_strided(points, output_stride, gradients_out);
        evaluate_triangle_order2_hessians_q4(output_stride, hessians_out);
        return true;
    }
    if (topology == LagrangeTopology::Tetrahedron &&
        order == 2 &&
        values_out != nullptr &&
        gradients_out == nullptr &&
        hessians_out == nullptr) {
        evaluate_tet_order2_values_strided(points, output_stride, values_out);
        return true;
    }
    if (topology == LagrangeTopology::Tetrahedron &&
        order == 2 &&
        values_out == nullptr &&
        gradients_out != nullptr &&
        hessians_out == nullptr) {
        evaluate_tet_order2_gradients_strided(points, output_stride, gradients_out);
        return true;
    }
    if (topology == LagrangeTopology::Tetrahedron &&
        order == 2 &&
        points.size() == 4u &&
        values_out == nullptr &&
        gradients_out == nullptr &&
        hessians_out != nullptr) {
        evaluate_tet_order2_hessians_q4(output_stride, hessians_out);
        return true;
    }
    if (topology == LagrangeTopology::Tetrahedron &&
        order == 2 &&
        points.size() == 4u &&
        values_out != nullptr &&
        gradients_out != nullptr &&
        hessians_out != nullptr) {
        evaluate_tet_order2_values_strided(points, output_stride, values_out);
        evaluate_tet_order2_gradients_strided(points, output_stride, gradients_out);
        evaluate_tet_order2_hessians_q4(output_stride, hessians_out);
        return true;
    }
    if (topology == LagrangeTopology::Tetrahedron &&
        order == 3 &&
        values_out != nullptr &&
        gradients_out == nullptr &&
        hessians_out == nullptr) {
        evaluate_tet_order3_values_strided(points, output_stride, values_out);
        return true;
    }
    if (topology == LagrangeTopology::Triangle &&
        order == 3 &&
        values_out != nullptr &&
        gradients_out == nullptr &&
        hessians_out == nullptr) {
        evaluate_triangle_order3_values_strided(points, output_stride, values_out);
        return true;
    }
    if (topology == LagrangeTopology::Triangle &&
        order == 3 &&
        values_out == nullptr &&
        gradients_out != nullptr &&
        hessians_out == nullptr) {
        evaluate_triangle_order3_gradients_strided(points, output_stride, gradients_out);
        return true;
    }
    if (topology == LagrangeTopology::Hexahedron &&
        order == 1 &&
        values_out != nullptr &&
        gradients_out == nullptr &&
        hessians_out == nullptr) {
        evaluate_hex_order1_values_strided(points, output_stride, values_out);
        return true;
    }
    if (topology == LagrangeTopology::Hexahedron &&
        order == 1 &&
        values_out == nullptr &&
        gradients_out != nullptr &&
        hessians_out == nullptr) {
        evaluate_hex_order1_outputs_strided<false, true, false>(
            points, output_stride, values_out, gradients_out, hessians_out);
        return true;
    }
    if (topology == LagrangeTopology::Hexahedron &&
        order == 1 &&
        values_out == nullptr &&
        gradients_out == nullptr &&
        hessians_out != nullptr) {
        evaluate_hex_order1_outputs_strided<false, false, true>(
            points, output_stride, values_out, gradients_out, hessians_out);
        return true;
    }
    if (topology == LagrangeTopology::Hexahedron &&
        order == 1 &&
        values_out != nullptr &&
        gradients_out != nullptr &&
        hessians_out != nullptr) {
        evaluate_hex_order1_outputs_strided<true, true, true>(
            points, output_stride, values_out, gradients_out, hessians_out);
        return true;
    }
    if (topology == LagrangeTopology::Quadrilateral &&
        order == 1 &&
        values_out != nullptr &&
        gradients_out == nullptr &&
        hessians_out == nullptr) {
        evaluate_quad_order1_values_strided(points, output_stride, values_out);
        return true;
    }
    if (topology == LagrangeTopology::Quadrilateral &&
        order == 1 &&
        values_out == nullptr &&
        gradients_out != nullptr &&
        hessians_out == nullptr) {
        evaluate_quad_order1_gradients_strided(points, output_stride, gradients_out);
        return true;
    }
    if (topology == LagrangeTopology::Quadrilateral &&
        order == 1 &&
        values_out == nullptr &&
        gradients_out == nullptr &&
        hessians_out != nullptr) {
        evaluate_quad_order1_hessians_strided(points.size(), output_stride, hessians_out);
        return true;
    }
    if (topology == LagrangeTopology::Quadrilateral &&
        order == 1 &&
        points.size() == 4u &&
        values_out != nullptr &&
        gradients_out != nullptr &&
        hessians_out != nullptr) {
        evaluate_quad_order1_all_q4(points, output_stride, values_out, gradients_out, hessians_out);
        return true;
    }
    if (topology == LagrangeTopology::Quadrilateral &&
        order == 2 &&
        values_out != nullptr &&
        gradients_out == nullptr &&
        hessians_out == nullptr) {
        evaluate_quad_order2_values_strided(points, output_stride, values_out);
        return true;
    }
    if (topology == LagrangeTopology::Quadrilateral &&
        order == 2 &&
        values_out == nullptr &&
        gradients_out != nullptr &&
        hessians_out == nullptr) {
        evaluate_quad_order2_gradients_strided(points, output_stride, gradients_out);
        return true;
    }
    if (topology == LagrangeTopology::Quadrilateral &&
        order == 2 &&
        values_out == nullptr &&
        gradients_out == nullptr &&
        hessians_out != nullptr) {
        evaluate_quad_order2_hessians_strided(points, output_stride, hessians_out);
        return true;
    }
    if (topology == LagrangeTopology::Quadrilateral &&
        order == 2 &&
        points.size() == 4u &&
        values_out != nullptr &&
        gradients_out != nullptr &&
        hessians_out != nullptr) {
        evaluate_quad_order2_all_q4(points, output_stride, values_out, gradients_out, hessians_out);
        return true;
    }
    if (topology == LagrangeTopology::Quadrilateral &&
        order == 3 &&
        values_out != nullptr &&
        gradients_out == nullptr &&
        hessians_out == nullptr) {
        evaluate_quad_order3_values_strided(points, output_stride, values_out);
        return true;
    }
    if (topology == LagrangeTopology::Quadrilateral &&
        order == 3 &&
        values_out == nullptr &&
        gradients_out != nullptr &&
        hessians_out == nullptr) {
        evaluate_quad_order3_gradients_strided(points, output_stride, gradients_out);
        return true;
    }
    if (topology == LagrangeTopology::Quadrilateral &&
        order == 3 &&
        values_out == nullptr &&
        gradients_out == nullptr &&
        hessians_out != nullptr) {
        evaluate_quad_order3_hessians_strided(points, output_stride, hessians_out);
        return true;
    }
    if (topology == LagrangeTopology::Quadrilateral &&
        order == 3 &&
        points.size() == 4u &&
        values_out != nullptr &&
        gradients_out != nullptr &&
        hessians_out != nullptr) {
        evaluate_quad_order3_all_q4(points, output_stride, values_out, gradients_out, hessians_out);
        return true;
    }
    if (topology == LagrangeTopology::Hexahedron &&
        order == 2 &&
        values_out != nullptr &&
        gradients_out == nullptr &&
        hessians_out == nullptr) {
        evaluate_hex_order2_values_strided(points, output_stride, values_out);
        return true;
    }
    if (topology == LagrangeTopology::Hexahedron &&
        order == 2 &&
        points.size() == 4u &&
        values_out == nullptr &&
        gradients_out != nullptr &&
        hessians_out == nullptr) {
        evaluate_hex_order2_gradients_q4(points, output_stride, gradients_out);
        return true;
    }
    if (topology == LagrangeTopology::Hexahedron &&
        order == 2 &&
        points.size() == 4u &&
        values_out == nullptr &&
        gradients_out == nullptr &&
        hessians_out != nullptr) {
        evaluate_hex_order2_hessians_q4(points, output_stride, hessians_out);
        return true;
    }
    if (topology == LagrangeTopology::Hexahedron &&
        order == 2 &&
        points.size() == 4u &&
        values_out != nullptr &&
        gradients_out != nullptr &&
        hessians_out != nullptr) {
        evaluate_hex_order2_values_strided(points, output_stride, values_out);
        evaluate_hex_order2_gradients_q4(points, output_stride, gradients_out);
        evaluate_hex_order2_hessians_q4(points, output_stride, hessians_out);
        return true;
    }
    if (topology == LagrangeTopology::Hexahedron &&
        order == 3 &&
        points.size() == 4u &&
        output_stride == 4u &&
        hessians_out != nullptr) {
        return false;
    }
    if (topology == LagrangeTopology::Hexahedron &&
        order == 3 &&
        points.size() == 4u &&
        values_out != nullptr &&
        gradients_out == nullptr &&
        hessians_out == nullptr) {
        evaluate_hex_order3_values_q4(points, output_stride, values_out);
        return true;
    }
    if (topology == LagrangeTopology::Hexahedron &&
        order == 3 &&
        points.size() == 4u &&
        values_out == nullptr &&
        gradients_out != nullptr &&
        hessians_out == nullptr) {
        evaluate_hex_order3_gradients_q4(points, output_stride, gradients_out);
        return true;
    }
    if (topology == LagrangeTopology::Hexahedron &&
        order == 3 &&
        points.size() == 4u &&
        values_out == nullptr &&
        gradients_out == nullptr &&
        hessians_out != nullptr) {
        evaluate_hex_order3_hessians_q4(points, output_stride, hessians_out);
        return true;
    }
    if (topology == LagrangeTopology::Hexahedron &&
        order == 3 &&
        points.size() == 4u &&
        values_out != nullptr &&
        gradients_out != nullptr &&
        hessians_out != nullptr) {
        evaluate_hex_order3_all_q4(points, output_stride, values_out, gradients_out, hessians_out);
        return true;
    }
    if (hessians_out != nullptr && order > 1 &&
        (topology == LagrangeTopology::Quadrilateral ||
         topology == LagrangeTopology::Hexahedron)) {
        return false;
    }
    if (hessians_out != nullptr) {
        const bool hessians_only = values_out == nullptr && gradients_out == nullptr;
        if (order == 1) {
            if (topology == LagrangeTopology::Triangle && hessians_only) {
                evaluate_fast_outputs_with_constant_hessians_strided<LagrangeTriFast<1>>(
                    points, output_stride, values_out, gradients_out, hessians_out);
                return true;
            }
            if (topology == LagrangeTopology::Tetrahedron) {
                evaluate_fast_outputs_with_constant_hessians_strided<LagrangeTetFast<1>>(
                    points, output_stride, values_out, gradients_out, hessians_out);
                return true;
            }
        } else if (order == 2) {
            if (topology == LagrangeTopology::Triangle && hessians_only) {
                evaluate_fast_outputs_with_constant_hessians_strided<LagrangeTriFast<2>>(
                    points, output_stride, values_out, gradients_out, hessians_out);
                return true;
            }
            if (topology == LagrangeTopology::Tetrahedron) {
                evaluate_fast_outputs_with_constant_hessians_strided<LagrangeTetFast<2>>(
                    points, output_stride, values_out, gradients_out, hessians_out);
                return true;
            }
        }
    }

    switch (order) {
        case 1:
            return evaluate_fixed_lagrange_fast_strided_order<1>(
                topology, points, output_stride, values_out, gradients_out, hessians_out);
        case 2:
            return evaluate_fixed_lagrange_fast_strided_order<2>(
                topology, points, output_stride, values_out, gradients_out, hessians_out);
        case 3:
            return evaluate_fixed_lagrange_fast_strided_order<3>(
                topology, points, output_stride, values_out, gradients_out, hessians_out);
        default:
            return false;
    }
}

bool evaluate_fixed_lagrange_fast_to(LagrangeTopology topology,
                                     int order,
                                     const math::Vector<Real, 3>& xi,
                                     Real* SVMP_RESTRICT values_out,
                                     Real* SVMP_RESTRICT gradients_out,
                                     Real* SVMP_RESTRICT hessians_out) {
    switch (order) {
        case 1:
            return evaluate_fixed_lagrange_fast_to_order<1>(
                topology, xi, values_out, gradients_out, hessians_out);
        case 2:
            return evaluate_fixed_lagrange_fast_to_order<2>(
                topology, xi, values_out, gradients_out, hessians_out);
        case 3:
            return evaluate_fixed_lagrange_fast_to_order<3>(
                topology, xi, values_out, gradients_out, hessians_out);
        default:
            return false;
    }
}

template<std::size_t N>
struct AxisMonomialCoefficientTable {
    std::array<Real, N * N> values{};
    std::array<Real, N * (N > 1 ? N - 1 : 0)> first{};
    std::array<Real, N * (N > 2 ? N - 2 : 0)> second{};
};

template<std::size_t N>
constexpr AxisMonomialCoefficientTable<N> make_axis_monomial_coefficient_table() {
    AxisMonomialCoefficientTable<N> table{};
    std::array<Real, N> nodes{};
    constexpr int order = static_cast<int>(N) - 1;
    for (std::size_t i = 0; i < N; ++i) {
        nodes[i] = detail::equispaced_pm_one_coord(static_cast<int>(i), order);
    }

    for (std::size_t i = 0; i < N; ++i) {
        std::array<Real, N> coeffs{};
        std::array<Real, N> next{};
        coeffs[0] = Real(1);
        std::size_t degree = 0;
        for (std::size_t j = 0; j < N; ++j) {
            if (j == i) {
                continue;
            }
            next = {};
            for (std::size_t k = 0; k <= degree; ++k) {
                next[k] -= nodes[j] * coeffs[k];
                next[k + 1] += coeffs[k];
            }
            coeffs = next;
            ++degree;
        }

        Real denominator = Real(1);
        for (std::size_t j = 0; j < N; ++j) {
            if (j != i) {
                denominator *= nodes[i] - nodes[j];
            }
        }
        const Real inv_denominator = Real(1) / denominator;
        for (std::size_t k = 0; k < N; ++k) {
            table.values[i * N + k] = coeffs[k] * inv_denominator;
        }
        if constexpr (N >= 2) {
            for (std::size_t k = 1; k < N; ++k) {
                table.first[i * (N - 1) + (k - 1)] =
                    static_cast<Real>(k) * table.values[i * N + k];
            }
        }
        if constexpr (N >= 3) {
            for (std::size_t k = 2; k < N; ++k) {
                table.second[i * (N - 2) + (k - 2)] =
                    static_cast<Real>(k * (k - 1)) * table.values[i * N + k];
            }
        }
    }

    return table;
}

template<std::size_t N>
void assign_axis_coefficient_table(const AxisMonomialCoefficientTable<N>& table,
                                   std::vector<Real>& values,
                                   std::vector<Real>& first,
                                   std::vector<Real>& second) {
    assign_array(values, table.values);
    assign_array(first, table.first);
    assign_array(second, table.second);
}

bool assign_precomputed_axis_coefficients(int n_axis,
                                          std::vector<Real>& values,
                                          std::vector<Real>& first,
                                          std::vector<Real>& second) {
    static constexpr auto kAxisCoefficients1 = make_axis_monomial_coefficient_table<1>();
    static constexpr auto kAxisCoefficients2 = make_axis_monomial_coefficient_table<2>();
    static constexpr auto kAxisCoefficients3 = make_axis_monomial_coefficient_table<3>();
    static constexpr auto kAxisCoefficients4 = make_axis_monomial_coefficient_table<4>();
    static constexpr auto kAxisCoefficients5 = make_axis_monomial_coefficient_table<5>();

    switch (n_axis) {
        case 1:
            assign_axis_coefficient_table(kAxisCoefficients1, values, first, second);
            return true;
        case 2:
            assign_axis_coefficient_table(kAxisCoefficients2, values, first, second);
            return true;
        case 3:
            assign_axis_coefficient_table(kAxisCoefficients3, values, first, second);
            return true;
        case 4:
            assign_axis_coefficient_table(kAxisCoefficients4, values, first, second);
            return true;
        case 5:
            assign_axis_coefficient_table(kAxisCoefficients5, values, first, second);
            return true;
        default:
            return false;
    }
}

LagrangeTopologyTraits lagrange_topology_traits(ElementType type) {
    const auto topo = topology(type);
    if (topo != LagrangeTopology::Unknown) {
        return {topo, reference_dimension(type)};
    }

    throw BasisElementCompatibilityException("Unsupported element type for LagrangeBasis",
                                             __FILE__, __LINE__, __func__);
}

std::size_t lattice_index_pm_one(Real coord, int order, const char* context) {
    if (order <= 0) {
        if (!coordinate_matches_expected(coord, Real(0))) {
            throw BasisNodeOrderingException(context, __FILE__, __LINE__, __func__);
        }
        return 0;
    }

    const Real scaled = (coord + Real(1)) * static_cast<Real>(order) / Real(2);
    const long idx = std::lround(scaled);
    if (idx < 0 || idx > order ||
        !coordinate_matches_expected(
            coord,
            detail::equispaced_pm_one_coord(static_cast<int>(idx), order))) {
        throw BasisNodeOrderingException(context, __FILE__, __LINE__, __func__);
    }
    return static_cast<std::size_t>(idx);
}

int simplex_lattice_index(Real coord, int order, const char* context) {
    if (order <= 0) {
        if (!coordinate_matches_expected(coord, Real(0)) &&
            !coordinate_matches_expected(coord, Real(0.25)) &&
            !coordinate_matches_expected(coord, Real(1) / Real(3))) {
            throw BasisNodeOrderingException(context, __FILE__, __LINE__, __func__);
        }
        return 0;
    }

    const Real scaled = coord * static_cast<Real>(order);
    const long idx = std::lround(scaled);
    const Real reconstructed = static_cast<Real>(idx) / static_cast<Real>(order);
    if (idx < 0 || idx > order || !coordinate_matches_expected(coord, reconstructed)) {
        throw BasisNodeOrderingException(context, __FILE__, __LINE__, __func__);
    }
    return static_cast<int>(idx);
}

std::array<int, 4> triangle_exponents_from_public_node(const math::Vector<Real, 3>& node,
                                                       int order) {
    if (order == 0) {
        return {0, 0, 0, 0};
    }

    const int j = simplex_lattice_index(node[0], order,
                                        "LagrangeBasis: invalid triangle node coordinate for public ordering");
    const int k = simplex_lattice_index(node[1], order,
                                        "LagrangeBasis: invalid triangle node coordinate for public ordering");
    const int i = order - j - k;
    if (i < 0) {
        throw BasisNodeOrderingException("LagrangeBasis: invalid triangle barycentric coordinates for public ordering",
                                         __FILE__, __LINE__, __func__);
    }
    return {i, j, k, 0};
}

std::array<int, 4> tetrahedron_exponents_from_public_node(const math::Vector<Real, 3>& node,
                                                          int order) {
    if (order == 0) {
        return {0, 0, 0, 0};
    }

    const int j = simplex_lattice_index(node[0], order,
                                        "LagrangeBasis: invalid tetrahedron node x-coordinate for public ordering");
    const int k = simplex_lattice_index(node[1], order,
                                        "LagrangeBasis: invalid tetrahedron node y-coordinate for public ordering");
    const int l = simplex_lattice_index(node[2], order,
                                        "LagrangeBasis: invalid tetrahedron node z-coordinate for public ordering");
    const int i = order - j - k - l;
    if (i < 0) {
        throw BasisNodeOrderingException("LagrangeBasis: invalid tetrahedron barycentric coordinates for public ordering",
                                         __FILE__, __LINE__, __func__);
    }
    return {i, j, k, l};
}

struct NormalizedLagrangeRequest {
    ElementType element_type;
    int order;
};

// Non-owning view of the per-axis 1D Lagrange basis evaluations
// (values, first derivative, second derivative), each of length `size`.
struct AxisBasisEvaluations {
    const Real* values;
    const Real* first;
    const Real* second;
    std::size_t size;
};

AxisBasisEvaluations constant_axis_basis() {
    static const Real kOne[1]  = {Real(1)};
    static const Real kZero[1] = {Real(0)};
    return AxisBasisEvaluations{kOne, kZero, kZero, 1};
}

// Horner-form evaluator for the precomputed 1D Lagrange basis.
//
// Inputs are precomputed monomial coefficients of L_i(x), L_i'(x), L_i''(x)
// (built once at LagrangeBasis construction). Evaluation is purely
// multiply-add on the coefficients — no divisions and no node-position
// lookups in the hot path. Templated on N for compile-time loop unrolling
// and FMA-friendly straight-line code on the common Hex/Quad/Line orders.
//
// Layout:
//   v_coeffs:  N * N entries; row i holds [c_i0, c_i1, ..., c_i(N-1)]
//              such that L_i(x) = sum_k c_ik * x^k
//   d_coeffs:  N * (N-1) entries; row i holds derivative coefficients of L_i'(x)
//   d2_coeffs: N * (N-2) entries; row i holds coefficients of L_i''(x)
//              (only valid when N >= 3)
template<int N>
inline void evaluate_1d_horner_impl(const Real* v_coeffs,
                                    const Real* d_coeffs,
                                    const Real* d2_coeffs,
                                    Real xi,
                                    Real* values, Real* first, Real* second) {
    if constexpr (N == 1) {
        values[0] = v_coeffs[0];
        if (first)  first[0]  = Real(0);
        if (second) second[0] = Real(0);
        return;
    } else {
        // Values: degree N-1 polynomials.
        for (int i = 0; i < N; ++i) {
            const Real* c = v_coeffs + i * N;
            Real r = c[N - 1];
            for (int k = N - 1; k > 0; --k) {
                r = r * xi + c[k - 1];
            }
            values[i] = r;
        }

        if (!first && !second) return;

        if (first) {
            // First derivatives: degree N-2 polynomials (per row of d_coeffs).
            for (int i = 0; i < N; ++i) {
                const Real* c = d_coeffs + i * (N - 1);
                Real r = c[N - 2];
                for (int k = N - 2; k > 0; --k) {
                    r = r * xi + c[k - 1];
                }
                first[i] = r;
            }
        }

        if (!second) return;

        if constexpr (N <= 2) {
            for (int i = 0; i < N; ++i) second[i] = Real(0);
        } else {
            // Second derivatives: degree N-3 polynomials (per row of d2_coeffs).
            for (int i = 0; i < N; ++i) {
                const Real* c = d2_coeffs + i * (N - 2);
                Real r = c[N - 3];
                for (int k = N - 3; k > 0; --k) {
                    r = r * xi + c[k - 1];
                }
                second[i] = r;
            }
        }
    }
}

void fill_equispaced_barycentric_weights(int n_axis, Real* weights) {
    const int order = n_axis - 1;
    Real weight = (order % 2 == 0) ? Real(1) : Real(-1);
    Real max_abs = Real(0);
    for (int i = 0; i < n_axis; ++i) {
        weights[i] = weight;
        max_abs = std::max(max_abs, std::abs(weight));
        if (i < order) {
            weight *= -static_cast<Real>(order - i) / static_cast<Real>(i + 1);
        }
    }

    if (max_abs > Real(0)) {
        const Real inv_scale = Real(1) / max_abs;
        for (int i = 0; i < n_axis; ++i) {
            weights[i] *= inv_scale;
        }
    }
}

bool coordinate_matches_axis_node(Real xi, Real node) {
    return coordinate_matches_expected(xi, node);
}

struct CompensatedSum {
    Real sum{Real(0)};
    Real compensation{Real(0)};

    void add(Real value) noexcept {
        const Real y = value - compensation;
        const Real t = sum + y;
        compensation = (t - sum) - y;
        sum = t;
    }
};

void distribute_residual_by_abs(int n_axis, Real* values, Real residual) {
    if (values == nullptr || n_axis <= 0 || residual == Real(0)) {
        return;
    }

    CompensatedSum abs_sum;
    int largest_index = 0;
    Real largest_abs = Real(0);
    for (int i = 0; i < n_axis; ++i) {
        const Real magnitude = std::abs(values[i]);
        abs_sum.add(magnitude);
        if (magnitude > largest_abs) {
            largest_abs = magnitude;
            largest_index = i;
        }
    }

    if (abs_sum.sum <= Real(0)) {
        values[0] += residual;
        return;
    }

    const Real inv_abs_sum = Real(1) / abs_sum.sum;
    CompensatedSum applied;
    for (int i = 0; i < n_axis; ++i) {
        const Real correction = residual * std::abs(values[i]) * inv_abs_sum;
        values[i] += correction;
        applied.add(correction);
    }
    values[largest_index] += residual - applied.sum;
}

void evaluate_1d_barycentric_runtime(int n_axis,
                                     Real xi,
                                     const Real* weights,
                                     Real* values,
                                     Real* first,
                                     Real* second) {
    const int order = n_axis - 1;
    BASIS_CHECK_EVAL(weights != nullptr,
                     "LagrangeBasis: missing cached barycentric weights for runtime axis evaluation");

    int node_index = -1;
    for (int i = 0; i < n_axis; ++i) {
        const Real node = detail::equispaced_pm_one_coord(i, order);
        if (coordinate_matches_axis_node(xi, node)) {
            node_index = i;
            break;
        }
    }

    if (node_index >= 0) {
        std::fill(values, values + n_axis, Real(0));
        values[node_index] = Real(1);
        if (!first && !second) {
            return;
        }

        const Real xk = detail::equispaced_pm_one_coord(node_index, order);
        const Real wk = weights[static_cast<std::size_t>(node_index)];
        Real reciprocal_sum = Real(0);
        if (second) {
            for (int m = 0; m < n_axis; ++m) {
                if (m == node_index) {
                    continue;
                }
                const Real xm = detail::equispaced_pm_one_coord(m, order);
                reciprocal_sum += Real(1) / (xk - xm);
            }
        }

        Real first_diagonal = Real(0);
        Real second_diagonal = Real(0);
        if (first) {
            std::fill(first, first + n_axis, Real(0));
        }
        if (second) {
            std::fill(second, second + n_axis, Real(0));
        }

        for (int j = 0; j < n_axis; ++j) {
            if (j == node_index) {
                continue;
            }
            const Real xj = detail::equispaced_pm_one_coord(j, order);
            const Real distance = xk - xj;
            const Real offdiag_first = weights[static_cast<std::size_t>(j)] / (wk * distance);
            first_diagonal -= offdiag_first;
            if (first) {
                first[j] = offdiag_first;
            }
            if (second) {
                const Real offdiag_second =
                    Real(2) * offdiag_first * (reciprocal_sum - Real(1) / distance);
                second[j] = offdiag_second;
                second_diagonal -= offdiag_second;
            }
        }
        if (first) {
            first[node_index] = first_diagonal;
        }
        if (second) {
            second[node_index] = second_diagonal;
        }
        return;
    }

    Real sum0 = Real(0);
    Real sum1 = Real(0);
    Real sum2 = Real(0);
    for (int i = 0; i < n_axis; ++i) {
        const Real node = detail::equispaced_pm_one_coord(i, order);
        const Real inv_distance = Real(1) / (xi - node);
        const Real weighted = weights[static_cast<std::size_t>(i)] * inv_distance;
        sum0 += weighted;
        sum1 += weighted * inv_distance;
        sum2 += weighted * inv_distance * inv_distance;
    }

    const Real inv_sum0 = Real(1) / sum0;
    const Real first_ratio = sum1 * inv_sum0;
    const Real second_ratio = sum2 * inv_sum0;
    const Real first_ratio_sq = first_ratio * first_ratio;

    CompensatedSum value_sum;
    CompensatedSum first_sum;
    CompensatedSum second_sum;
    for (int i = 0; i < n_axis; ++i) {
        const Real node = detail::equispaced_pm_one_coord(i, order);
        const Real inv_distance = Real(1) / (xi - node);
        const Real value = weights[static_cast<std::size_t>(i)] * inv_distance * inv_sum0;
        values[i] = value;
        value_sum.add(value);
        if (first || second) {
            const Real derivative_factor = first_ratio - inv_distance;
            if (first) {
                first[i] = value * derivative_factor;
                first_sum.add(first[i]);
            }
            if (second) {
                second[i] = value * (derivative_factor * derivative_factor +
                                     inv_distance * inv_distance -
                                     Real(2) * second_ratio +
                                     first_ratio_sq);
                second_sum.add(second[i]);
            }
        }
    }

    distribute_residual_by_abs(n_axis, values, Real(1) - value_sum.sum);
    if (first) {
        distribute_residual_by_abs(n_axis, first, -first_sum.sum);
    }
    if (second) {
        distribute_residual_by_abs(n_axis, second, -second_sum.sum);
    }
}

// 1D Lagrange-basis evaluator. Writes n_axis entries to each non-null output
// buffer. Dispatches to compile-time Horner specializations for sizes 1..9
// (orders 0..8 — the Lagrange performance sweep) and uses barycentric
// evaluation above that threshold to avoid high-order monomial conditioning
// issues.
void evaluate_1d_basis_to(const Real* v_coeffs,
                          const Real* d_coeffs,
                          const Real* d2_coeffs,
                          const Real* barycentric_weights,
                          int n_axis, Real xi,
                          Real* values, Real* first, Real* second) {
    switch (n_axis) {
        case 1: evaluate_1d_horner_impl<1>(v_coeffs, d_coeffs, d2_coeffs, xi, values, first, second); return;
        case 2: evaluate_1d_horner_impl<2>(v_coeffs, d_coeffs, d2_coeffs, xi, values, first, second); return;
        case 3: evaluate_1d_horner_impl<3>(v_coeffs, d_coeffs, d2_coeffs, xi, values, first, second); return;
        case 4: evaluate_1d_horner_impl<4>(v_coeffs, d_coeffs, d2_coeffs, xi, values, first, second); return;
        case 5: evaluate_1d_horner_impl<5>(v_coeffs, d_coeffs, d2_coeffs, xi, values, first, second); return;
        case 6: evaluate_1d_horner_impl<6>(v_coeffs, d_coeffs, d2_coeffs, xi, values, first, second); return;
        case 7: evaluate_1d_horner_impl<7>(v_coeffs, d_coeffs, d2_coeffs, xi, values, first, second); return;
        case 8: evaluate_1d_horner_impl<8>(v_coeffs, d_coeffs, d2_coeffs, xi, values, first, second); return;
        case 9: evaluate_1d_horner_impl<9>(v_coeffs, d_coeffs, d2_coeffs, xi, values, first, second); return;
        default:
            evaluate_1d_barycentric_runtime(n_axis, xi, barycentric_weights, values, first, second);
            return;
    }
}

// Selects which derivative passes are computed by the 1D evaluator.
enum class AxisDeriv {
    ValuesOnly,           // skip first and second
    ValuesAndFirst,       // for gradients
    ValuesAndFirstAndSecond, // for hessians or fused evaluate_all
};

// Per-axis storage (values, first derivative, second derivative). Backed by
// per-thread scratch that grows lazily; subsequent calls reuse capacity with no
// reallocation.
struct AxisScratch {
    std::vector<Real> values;
    std::vector<Real> first;
    std::vector<Real> second;

    void reserveFor(std::size_t n) {
        if (values.size() < n) values.resize(n);
        if (first.size() < n) first.resize(n);
        if (second.size() < n) second.resize(n);
    }
};

struct AxisBatchScratch {
    std::vector<Real> values;
    std::vector<Real> first;
    std::vector<Real> second;

    void resizeFor(std::size_t count, AxisDeriv level) {
        if (values.size() < count) values.resize(count);
        if (level != AxisDeriv::ValuesOnly && first.size() < count) first.resize(count);
        if (level == AxisDeriv::ValuesAndFirstAndSecond && second.size() < count) second.resize(count);
    }
};

template<int Order, bool NeedFirst, bool NeedSecond>
inline void fill_simplex_factor_sequence_fixed(Real lambda,
                                               Real* SVMP_RESTRICT phi,
                                               Real* SVMP_RESTRICT dphi,
                                               Real* SVMP_RESTRICT d2phi) {
    static_assert(!NeedSecond || NeedFirst,
                  "second derivative factors require first-derivative recurrence state");
    phi[0] = Real(1);
    if constexpr (NeedFirst) {
        dphi[0] = Real(0);
    }
    if constexpr (NeedSecond) {
        d2phi[0] = Real(0);
    }

    const Real t = static_cast<Real>(Order) * lambda;
    const Real dt_dlambda = static_cast<Real>(Order);
    Real dphi_dt_prev = Real(0);
    Real d2phi_dt2_prev = Real(0);
    for (int a = 1; a <= Order; ++a) {
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

template<int Order, bool NeedSecond>
inline void fill_triangle_factors_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    Real (&phi0)[4][Order + 1],
    Real (&phi1)[4][Order + 1],
    Real (&phi2)[4][Order + 1],
    Real (&dphi0)[4][Order + 1],
    Real (&dphi1)[4][Order + 1],
    Real (&dphi2)[4][Order + 1],
    Real (&d2phi0)[4][Order + 1],
    Real (&d2phi1)[4][Order + 1],
    Real (&d2phi2)[4][Order + 1]) {
    for (std::size_t q = 0; q < 4u; ++q) {
        const auto& xi = points[q];
        const Real l1 = xi[0];
        const Real l2 = xi[1];
        const Real l0 = Real(1) - l1 - l2;
        if constexpr (NeedSecond) {
            fill_simplex_factor_sequence_fixed<Order, true, true>(
                l0, phi0[q], dphi0[q], d2phi0[q]);
            fill_simplex_factor_sequence_fixed<Order, true, true>(
                l1, phi1[q], dphi1[q], d2phi1[q]);
            fill_simplex_factor_sequence_fixed<Order, true, true>(
                l2, phi2[q], dphi2[q], d2phi2[q]);
        } else {
            fill_simplex_factor_sequence_fixed<Order, true, false>(
                l0, phi0[q], dphi0[q], nullptr);
            fill_simplex_factor_sequence_fixed<Order, true, false>(
                l1, phi1[q], dphi1[q], nullptr);
            fill_simplex_factor_sequence_fixed<Order, true, false>(
                l2, phi2[q], dphi2[q], nullptr);
        }
    }
}

template<std::size_t Q>
inline void write_wedge_gradient_strided_q(std::size_t tri_stride,
                                           std::size_t axis_stride,
                                           std::size_t tri,
                                           std::size_t z,
                                           std::size_t output_stride,
                                           const Real* SVMP_RESTRICT tri_values,
                                           const Real* SVMP_RESTRICT tri_g,
                                           const AxisBatchScratch& axis_batch,
                                           Real* SVMP_RESTRICT g) {
    const std::size_t tri_q = tri * tri_stride + Q;
    const std::size_t z_q = Q * axis_stride + z;
    const Real tri_v = tri_values[tri_q];
    const Real zv = axis_batch.values[z_q];
    g[0u * output_stride + Q] = tri_g[0u * tri_stride + Q] * zv;
    g[1u * output_stride + Q] = tri_g[1u * tri_stride + Q] * zv;
    g[2u * output_stride + Q] = tri_v * axis_batch.first[z_q];
}

template<std::size_t Q>
inline void write_wedge_gradient_stride4_q(std::size_t tri_stride,
                                           std::size_t axis_stride,
                                           std::size_t tri,
                                           std::size_t z,
                                           const Real* SVMP_RESTRICT tri_values,
                                           const Real* SVMP_RESTRICT tri_g,
                                           const AxisBatchScratch& axis_batch,
                                           Real* SVMP_RESTRICT g) {
    const std::size_t tri_q = tri * tri_stride + Q;
    const std::size_t z_q = Q * axis_stride + z;
    const Real tri_v = tri_values[tri_q];
    const Real zv = axis_batch.values[z_q];
    g[Q] = tri_g[0u * tri_stride + Q] * zv;
    g[4u + Q] = tri_g[1u * tri_stride + Q] * zv;
    g[8u + Q] = tri_v * axis_batch.first[z_q];
}

template<std::size_t Q>
inline void write_wedge_hessian_strided_q(std::size_t tri_stride,
                                          std::size_t axis_stride,
                                          std::size_t tri,
                                          std::size_t z,
                                          std::size_t output_stride,
                                          const Real* SVMP_RESTRICT tri_values,
                                          const Real* SVMP_RESTRICT tri_g,
                                          const Real* SVMP_RESTRICT tri_H,
                                          const AxisBatchScratch& axis_batch,
                                          Real* SVMP_RESTRICT H) {
    const std::size_t tri_q = tri * tri_stride + Q;
    const std::size_t z_q = Q * axis_stride + z;
    const Real tri_v = tri_values[tri_q];
    const Real zv = axis_batch.values[z_q];
    const Real zd = axis_batch.first[z_q];
    const Real tri_gx = tri_g[0u * tri_stride + Q];
    const Real tri_gy = tri_g[1u * tri_stride + Q];
    const Real tri_hxx = tri_H[0u * tri_stride + Q];
    const Real tri_hxy = tri_H[1u * tri_stride + Q];
    const Real tri_hyy = tri_H[2u * tri_stride + Q];
    const Real hxz = tri_gx * zd;
    const Real hxy = tri_hxy * zv;
    const Real hyz = tri_gy * zd;

    H[0u * output_stride + Q] = tri_hxx * zv;
    H[1u * output_stride + Q] = hxy;
    H[2u * output_stride + Q] = hxz;
    H[3u * output_stride + Q] = hxy;
    H[4u * output_stride + Q] = tri_hyy * zv;
    H[5u * output_stride + Q] = hyz;
    H[6u * output_stride + Q] = hxz;
    H[7u * output_stride + Q] = hyz;
    H[8u * output_stride + Q] = tri_v * axis_batch.second[z_q];
}

template<std::size_t Q>
inline void write_wedge_hessian_stride4_q(std::size_t tri_stride,
                                          std::size_t axis_stride,
                                          std::size_t tri,
                                          std::size_t z,
                                          const Real* SVMP_RESTRICT tri_values,
                                          const Real* SVMP_RESTRICT tri_g,
                                          const Real* SVMP_RESTRICT tri_H,
                                          const AxisBatchScratch& axis_batch,
                                          Real* SVMP_RESTRICT H) {
    const std::size_t tri_q = tri * tri_stride + Q;
    const std::size_t z_q = Q * axis_stride + z;
    const Real tri_v = tri_values[tri_q];
    const Real zv = axis_batch.values[z_q];
    const Real zd = axis_batch.first[z_q];
    const Real tri_gx = tri_g[0u * tri_stride + Q];
    const Real tri_gy = tri_g[1u * tri_stride + Q];
    const Real tri_hxx = tri_H[0u * tri_stride + Q];
    const Real tri_hxy = tri_H[1u * tri_stride + Q];
    const Real tri_hyy = tri_H[2u * tri_stride + Q];
    const Real hxz = tri_gx * zd;
    const Real hxy = tri_hxy * zv;
    const Real hyz = tri_gy * zd;

    H[Q] = tri_hxx * zv;
    H[4u + Q] = hxy;
    H[8u + Q] = hxz;
    H[12u + Q] = hxy;
    H[16u + Q] = tri_hyy * zv;
    H[20u + Q] = hyz;
    H[24u + Q] = hxz;
    H[28u + Q] = hyz;
    H[32u + Q] = tri_v * axis_batch.second[z_q];
}

template<std::size_t Q>
inline void write_wedge_all_strided_q(std::size_t tri_stride,
                                      std::size_t axis_stride,
                                      std::size_t tri,
                                      std::size_t z,
                                      std::size_t output_stride,
                                      const Real* SVMP_RESTRICT tri_values,
                                      const Real* SVMP_RESTRICT tri_g,
                                      const Real* SVMP_RESTRICT tri_H,
                                      const AxisBatchScratch& axis_batch,
                                      Real* SVMP_RESTRICT value_row,
                                      Real* SVMP_RESTRICT g,
                                      Real* SVMP_RESTRICT H) {
    const std::size_t tri_q = tri * tri_stride + Q;
    const std::size_t z_q = Q * axis_stride + z;
    const Real tri_v = tri_values[tri_q];
    const Real zv = axis_batch.values[z_q];
    const Real zd = axis_batch.first[z_q];
    const Real tri_gx = tri_g[0u * tri_stride + Q];
    const Real tri_gy = tri_g[1u * tri_stride + Q];
    const Real tri_hxx = tri_H[0u * tri_stride + Q];
    const Real tri_hxy = tri_H[1u * tri_stride + Q];
    const Real tri_hyy = tri_H[2u * tri_stride + Q];
    const Real hxz = tri_gx * zd;
    const Real hxy = tri_hxy * zv;
    const Real hyz = tri_gy * zd;

    value_row[Q] = tri_v * zv;
    g[0u * output_stride + Q] = tri_gx * zv;
    g[1u * output_stride + Q] = tri_gy * zv;
    g[2u * output_stride + Q] = tri_v * zd;
    H[0u * output_stride + Q] = tri_hxx * zv;
    H[1u * output_stride + Q] = hxy;
    H[2u * output_stride + Q] = hxz;
    H[3u * output_stride + Q] = hxy;
    H[4u * output_stride + Q] = tri_hyy * zv;
    H[5u * output_stride + Q] = hyz;
    H[6u * output_stride + Q] = hxz;
    H[7u * output_stride + Q] = hyz;
    H[8u * output_stride + Q] = tri_v * axis_batch.second[z_q];
}

template<std::size_t Q>
inline void write_wedge_all_stride4_q(std::size_t tri_stride,
                                      std::size_t axis_stride,
                                      std::size_t tri,
                                      std::size_t z,
                                      const Real* SVMP_RESTRICT tri_values,
                                      const Real* SVMP_RESTRICT tri_g,
                                      const Real* SVMP_RESTRICT tri_H,
                                      const AxisBatchScratch& axis_batch,
                                      Real* SVMP_RESTRICT value_row,
                                      Real* SVMP_RESTRICT g,
                                      Real* SVMP_RESTRICT H) {
    const std::size_t tri_q = tri * tri_stride + Q;
    const std::size_t z_q = Q * axis_stride + z;
    const Real tri_v = tri_values[tri_q];
    const Real zv = axis_batch.values[z_q];
    const Real zd = axis_batch.first[z_q];
    const Real tri_gx = tri_g[0u * tri_stride + Q];
    const Real tri_gy = tri_g[1u * tri_stride + Q];
    const Real tri_hxx = tri_H[0u * tri_stride + Q];
    const Real tri_hxy = tri_H[1u * tri_stride + Q];
    const Real tri_hyy = tri_H[2u * tri_stride + Q];
    const Real hxz = tri_gx * zd;
    const Real hxy = tri_hxy * zv;
    const Real hyz = tri_gy * zd;

    value_row[Q] = tri_v * zv;
    g[Q] = tri_gx * zv;
    g[4u + Q] = tri_gy * zv;
    g[8u + Q] = tri_v * zd;
    H[Q] = tri_hxx * zv;
    H[4u + Q] = hxy;
    H[8u + Q] = hxz;
    H[12u + Q] = hxy;
    H[16u + Q] = tri_hyy * zv;
    H[20u + Q] = hyz;
    H[24u + Q] = hxz;
    H[28u + Q] = hyz;
    H[32u + Q] = tri_v * axis_batch.second[z_q];
}

template<int Order, bool NeedHess>
bool evaluate_wedge_fused_stride4_q4(
    const std::vector<std::array<int, 4>>& simplex_exponents,
    const std::vector<std::size_t>& wedge_node_by_tri_z,
    const std::vector<math::Vector<Real, 3>>& points,
    const AxisBatchScratch& axis_batch,
    int n_axis,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    static_assert(Order >= 3 && Order <= 8, "fused wedge q4 path covers orders 3..8");
    const std::size_t tri_count = simplex_exponents.size();
    const std::size_t z_count = static_cast<std::size_t>(n_axis);
    if (points.size() != 4u ||
        z_count != static_cast<std::size_t>(Order + 1) ||
        wedge_node_by_tri_z.size() != tri_count * z_count) {
        return false;
    }

    Real phi0[4][Order + 1];
    Real phi1[4][Order + 1];
    Real phi2[4][Order + 1];
    Real dphi0[4][Order + 1];
    Real dphi1[4][Order + 1];
    Real dphi2[4][Order + 1];
    Real d2phi0[4][Order + 1];
    Real d2phi1[4][Order + 1];
    Real d2phi2[4][Order + 1];
    fill_triangle_factors_q4<Order, NeedHess>(
        points, phi0, phi1, phi2, dphi0, dphi1, dphi2, d2phi0, d2phi1, d2phi2);

    for (std::size_t tri = 0; tri < tri_count; ++tri) {
        const auto& e = simplex_exponents[tri];
        const std::size_t i0 = static_cast<std::size_t>(e[0]);
        const std::size_t i1 = static_cast<std::size_t>(e[1]);
        const std::size_t i2 = static_cast<std::size_t>(e[2]);

        Real tri_v[4];
        Real tri_gx[4];
        Real tri_gy[4];
        Real tri_hxx[4];
        Real tri_hxy[4];
        Real tri_hyy[4];
        for (std::size_t q = 0; q < 4u; ++q) {
            const Real v0 = phi0[q][i0];
            const Real v1 = phi1[q][i1];
            const Real v2 = phi2[q][i2];
            const Real D0 = dphi0[q][i0];
            const Real D1 = dphi1[q][i1];
            const Real D2 = dphi2[q][i2];
            const Real dl0 = D0 * v1 * v2;
            tri_v[q] = v0 * v1 * v2;
            tri_gx[q] = v0 * D1 * v2 - dl0;
            tri_gy[q] = v0 * v1 * D2 - dl0;

            if constexpr (NeedHess) {
                const Real DD0 = d2phi0[q][i0];
                const Real DD1 = d2phi1[q][i1];
                const Real DD2 = d2phi2[q][i2];
                const Real H00 = DD0 * v1 * v2;
                const Real H11 = v0 * DD1 * v2;
                const Real H22 = v0 * v1 * DD2;
                const Real H01 = D0 * D1 * v2;
                const Real H02 = D0 * v1 * D2;
                const Real H12 = v0 * D1 * D2;
                tri_hxx[q] = H00 - Real(2) * H01 + H11;
                tri_hxy[q] = H00 - H01 - H02 + H12;
                tri_hyy[q] = H00 - Real(2) * H02 + H22;
            }
        }

        for (std::size_t z = 0; z < z_count; ++z) {
            const std::size_t node = wedge_node_by_tri_z[tri * z_count + z];
            Real* SVMP_RESTRICT value_row =
                values_out != nullptr ? values_out + node * 4u : nullptr;
            Real* SVMP_RESTRICT g =
                gradients_out != nullptr ? gradients_out + node * 12u : nullptr;
            Real* SVMP_RESTRICT H =
                hessians_out != nullptr ? hessians_out + node * 36u : nullptr;

            const Real z0 = axis_batch.values[z];
            const Real z1 = axis_batch.values[z_count + z];
            const Real z2 = axis_batch.values[2u * z_count + z];
            const Real z3 = axis_batch.values[3u * z_count + z];
            const Real dz0 = axis_batch.first[z];
            const Real dz1 = axis_batch.first[z_count + z];
            const Real dz2 = axis_batch.first[2u * z_count + z];
            const Real dz3 = axis_batch.first[3u * z_count + z];

            if (value_row != nullptr) {
                value_row[0] = tri_v[0] * z0;
                value_row[1] = tri_v[1] * z1;
                value_row[2] = tri_v[2] * z2;
                value_row[3] = tri_v[3] * z3;
            }
            if (g != nullptr) {
                g[0] = tri_gx[0] * z0;
                g[1] = tri_gx[1] * z1;
                g[2] = tri_gx[2] * z2;
                g[3] = tri_gx[3] * z3;
                g[4] = tri_gy[0] * z0;
                g[5] = tri_gy[1] * z1;
                g[6] = tri_gy[2] * z2;
                g[7] = tri_gy[3] * z3;
                g[8] = tri_v[0] * dz0;
                g[9] = tri_v[1] * dz1;
                g[10] = tri_v[2] * dz2;
                g[11] = tri_v[3] * dz3;
            }
            if constexpr (NeedHess) {
                if (H != nullptr) {
                    const Real d2z0 = axis_batch.second[z];
                    const Real d2z1 = axis_batch.second[z_count + z];
                    const Real d2z2 = axis_batch.second[2u * z_count + z];
                    const Real d2z3 = axis_batch.second[3u * z_count + z];
                    const Real hxz0 = tri_gx[0] * dz0;
                    const Real hxz1 = tri_gx[1] * dz1;
                    const Real hxz2 = tri_gx[2] * dz2;
                    const Real hxz3 = tri_gx[3] * dz3;
                    const Real hyz0 = tri_gy[0] * dz0;
                    const Real hyz1 = tri_gy[1] * dz1;
                    const Real hyz2 = tri_gy[2] * dz2;
                    const Real hyz3 = tri_gy[3] * dz3;
                    H[0] = tri_hxx[0] * z0;
                    H[1] = tri_hxx[1] * z1;
                    H[2] = tri_hxx[2] * z2;
                    H[3] = tri_hxx[3] * z3;
                    H[4] = tri_hxy[0] * z0;
                    H[5] = tri_hxy[1] * z1;
                    H[6] = tri_hxy[2] * z2;
                    H[7] = tri_hxy[3] * z3;
                    H[8] = hxz0;
                    H[9] = hxz1;
                    H[10] = hxz2;
                    H[11] = hxz3;
                    H[12] = H[4];
                    H[13] = H[5];
                    H[14] = H[6];
                    H[15] = H[7];
                    H[16] = tri_hyy[0] * z0;
                    H[17] = tri_hyy[1] * z1;
                    H[18] = tri_hyy[2] * z2;
                    H[19] = tri_hyy[3] * z3;
                    H[20] = hyz0;
                    H[21] = hyz1;
                    H[22] = hyz2;
                    H[23] = hyz3;
                    H[24] = hxz0;
                    H[25] = hxz1;
                    H[26] = hxz2;
                    H[27] = hxz3;
                    H[28] = hyz0;
                    H[29] = hyz1;
                    H[30] = hyz2;
                    H[31] = hyz3;
                    H[32] = tri_v[0] * d2z0;
                    H[33] = tri_v[1] * d2z1;
                    H[34] = tri_v[2] * d2z2;
                    H[35] = tri_v[3] * d2z3;
                }
            }
        }
    }
    return true;
}

template<bool NeedHess>
bool try_evaluate_wedge_fused_stride4_q4(
    const std::vector<std::array<int, 4>>& simplex_exponents,
    const std::vector<std::size_t>& wedge_node_by_tri_z,
    int order,
    const std::vector<math::Vector<Real, 3>>& points,
    const AxisBatchScratch& axis_batch,
    int n_axis,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    switch (order) {
        case 3:
            return evaluate_wedge_fused_stride4_q4<3, NeedHess>(
                simplex_exponents, wedge_node_by_tri_z, points, axis_batch, n_axis,
                values_out, gradients_out, hessians_out);
        case 4:
            return evaluate_wedge_fused_stride4_q4<4, NeedHess>(
                simplex_exponents, wedge_node_by_tri_z, points, axis_batch, n_axis,
                values_out, gradients_out, hessians_out);
        case 5:
            return evaluate_wedge_fused_stride4_q4<5, NeedHess>(
                simplex_exponents, wedge_node_by_tri_z, points, axis_batch, n_axis,
                values_out, gradients_out, hessians_out);
        case 6:
            return evaluate_wedge_fused_stride4_q4<6, NeedHess>(
                simplex_exponents, wedge_node_by_tri_z, points, axis_batch, n_axis,
                values_out, gradients_out, hessians_out);
        case 7:
            return evaluate_wedge_fused_stride4_q4<7, NeedHess>(
                simplex_exponents, wedge_node_by_tri_z, points, axis_batch, n_axis,
                values_out, gradients_out, hessians_out);
        case 8:
            return evaluate_wedge_fused_stride4_q4<8, NeedHess>(
                simplex_exponents, wedge_node_by_tri_z, points, axis_batch, n_axis,
                values_out, gradients_out, hessians_out);
        default:
            return false;
    }
}

struct TensorProductTableScratch {
    std::vector<Real> vv;
    std::vector<Real> dv;
    std::vector<Real> vd;
    std::vector<Real> d2v;
    std::vector<Real> vd2;
    std::vector<Real> dd;

    void resizeFor(std::size_t count) {
        if (vv.size() < count) vv.resize(count);
        if (dv.size() < count) dv.resize(count);
        if (vd.size() < count) vd.resize(count);
        if (d2v.size() < count) d2v.resize(count);
        if (vd2.size() < count) vd2.resize(count);
        if (dd.size() < count) dd.resize(count);
    }
};

// Caller-provided scratch buffers used by tensor-product evaluation. Three
// independent axes plus reusable simplex/wedge intermediates.
struct LagrangeEvaluateScratch {
    AxisScratch axis_x;
    AxisScratch axis_y;
    AxisScratch axis_z;
    AxisBatchScratch axis_x_batch;
    AxisBatchScratch axis_y_batch;
    AxisBatchScratch axis_z_batch;
    TensorProductTableScratch tensor_tables;

    std::vector<Real> tri_values;
    std::vector<Gradient> tri_gradients;
    std::vector<Hessian> tri_hessians;
    std::vector<Real> tri_gradient_components;
    std::vector<Real> tri_hessian_components;
    std::vector<Real> wedge_tri_values_batch;
    std::vector<Real> wedge_tri_gradient_batch;
    std::vector<Real> wedge_tri_hessian_batch;

    std::vector<Real> strided_values_tmp;
    std::vector<Real> strided_gradients_tmp;
    std::vector<Real> strided_hessians_tmp;

    void prewarm(int max_order, std::size_t max_qpts) {
        const int clamped_order = std::max(max_order, 0);
        const std::size_t axis_size = static_cast<std::size_t>(clamped_order) + 1u;
        const std::size_t axis_batch_size = axis_size * max_qpts;
        const std::size_t tensor_table_size =
            axis_size * axis_size * std::max<std::size_t>(max_qpts, 1u);
        const std::size_t tensor_dofs = tensor_table_size * axis_size;
        const std::size_t tri_count = axis_size * (axis_size + 1u) / 2u;

        axis_x.reserveFor(axis_size);
        axis_y.reserveFor(axis_size);
        axis_z.reserveFor(axis_size);
        axis_x_batch.resizeFor(axis_batch_size, AxisDeriv::ValuesAndFirstAndSecond);
        axis_y_batch.resizeFor(axis_batch_size, AxisDeriv::ValuesAndFirstAndSecond);
        axis_z_batch.resizeFor(axis_batch_size, AxisDeriv::ValuesAndFirstAndSecond);
        tensor_tables.resizeFor(tensor_table_size);
        tri_values.reserve(tri_count);
        tri_gradients.reserve(tri_count);
        tri_hessians.reserve(tri_count);
        tri_gradient_components.reserve(tri_count * 3u);
        tri_hessian_components.reserve(tri_count * 9u);
        wedge_tri_values_batch.reserve(tri_count * max_qpts);
        wedge_tri_gradient_batch.reserve(tri_count * 3u * max_qpts);
        wedge_tri_hessian_batch.reserve(tri_count * 9u * max_qpts);
        strided_values_tmp.reserve(tensor_dofs);
        strided_gradients_tmp.reserve(tensor_dofs * 3u);
        strided_hessians_tmp.reserve(tensor_dofs * 9u);
    }
};

LagrangeEvaluateScratch& evaluate_scratch() {
    // Scratch is intentionally thread-local: assembly and benchmark callers run
    // evaluation on persistent worker threads, so capacity is reused by thread.
    static thread_local LagrangeEvaluateScratch s;
    return s;
}

// Fill axis scratch and return a non-owning view. Uncomputed slots still have
// valid pointers to scratch storage (they may hold stale data) — callers must
// only read the slots they requested via `level`. Common low orders use
// precomputed Horner coefficients; high orders use barycentric axis evaluation.
AxisBasisEvaluations fill_axis_scratch(AxisScratch& s,
                                       const Real* v_coeffs,
                                       const Real* d_coeffs,
                                       const Real* d2_coeffs,
                                       const Real* barycentric_weights,
                                       int n_axis, Real xi,
                                       AxisDeriv level) {
    const std::size_t n = static_cast<std::size_t>(n_axis);
    s.reserveFor(n);
    Real* first  = (level == AxisDeriv::ValuesOnly) ? nullptr : s.first.data();
    Real* second = (level == AxisDeriv::ValuesAndFirstAndSecond) ? s.second.data() : nullptr;
    evaluate_1d_basis_to(v_coeffs, d_coeffs, d2_coeffs, barycentric_weights,
                         n_axis, xi, s.values.data(), first, second);
    return AxisBasisEvaluations{s.values.data(), s.first.data(), s.second.data(), n};
}

void fill_axis_batch(AxisBatchScratch& scratch,
                     const std::vector<math::Vector<Real, 3>>& points,
                     std::size_t component,
                     const Real* v_coeffs,
                     const Real* d_coeffs,
                     const Real* d2_coeffs,
                     const Real* barycentric_weights,
                     int n_axis,
                     AxisDeriv level) {
    const std::size_t count = points.size() * static_cast<std::size_t>(n_axis);
    scratch.resizeFor(count, level);
    Real* first = (level == AxisDeriv::ValuesOnly) ? nullptr : scratch.first.data();
    Real* second = (level == AxisDeriv::ValuesAndFirstAndSecond) ? scratch.second.data() : nullptr;
    const std::size_t axis_stride = static_cast<std::size_t>(n_axis);
    for (std::size_t q = 0; q < points.size(); ++q) {
        evaluate_1d_basis_to(v_coeffs, d_coeffs, d2_coeffs, barycentric_weights, n_axis,
                             points[q][component],
                             scratch.values.data() + q * axis_stride,
                             first ? first + q * axis_stride : nullptr,
                             second ? second + q * axis_stride : nullptr);
    }
}

// Maximum yz-table footprint that fits comfortably on the stack for the
// Lagrange performance sweep. Order-8 hex q=4 needs 4*(9x9) entries per table.
// Higher orders fall back to thread_local heap buffers.
inline constexpr std::size_t kMaxStackYZ = 384;

struct TensorProductVectorSink {
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

    void write_gradient(std::size_t n, Real dx, Real dy, Real dz) const {
        auto& g = (*gradients)[n];
        g[0] = dx;
        g[1] = dy;
        g[2] = dz;
    }

    void write_hessian(std::size_t n,
                       Real xx,
                       Real yy,
                       Real zz,
                       Real xy,
                       Real xz,
                       Real yz) const {
        (*hessians)[n] = make_symmetric_hessian(xx, yy, zz, xy, xz, yz);
    }
};

struct TensorProductRawSink {
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

    void write_gradient(std::size_t n, Real dx, Real dy, Real dz) const {
        Real* g = gradients + n * 3u;
        g[0] = dx;
        g[1] = dy;
        g[2] = dz;
    }

    void write_hessian(std::size_t n,
                       Real xx,
                       Real yy,
                       Real zz,
                       Real xy,
                       Real xz,
                       Real yz) const {
        Real* H = hessians + n * 9u;
        H[0] = xx;
        H[4] = yy;
        H[8] = zz;
        H[1] = xy; H[3] = xy;
        H[2] = xz; H[6] = xz;
        H[5] = yz; H[7] = yz;
    }
};

// Fused sum-factorized tensor-product evaluator.
//
// Precomputes one to six (ny x nz)-shaped tables of partial products
// `M_xy[j*nz + k]` so that the inner per-node loop performs at most one
// multiplication per output instead of two. With all three output buffers
// supplied, this is the fused values + gradients + hessians path that shares
// every per-axis evaluation.
//
// Per-node multiply count (vs. the unfactored variants):
//   values only       : 1  (was 2)
//   gradients only    : 3  (was 6)
//   hessians only     : 6  (was 12)
//   all three         : 10 (was 20)
//
// Dimensional scope: works uniformly for Line/Quadrilateral/Hexahedron with
// the unused axes' size folded to 1 via constant_axis_basis().
template <typename Sink>
void evaluate_tensor_product_factorized_impl(
    const std::vector<std::array<std::size_t, 3>>& tensor_indices,
    const AxisBasisEvaluations& x_axis,
    const AxisBasisEvaluations& y_axis,
    const AxisBasisEvaluations& z_axis,
    const Sink& sink) {
    const std::size_t ny = y_axis.size;
    const std::size_t nz = z_axis.size;
    const std::size_t nyz = ny * nz;
    const bool need_values = sink.wants_values();
    const bool need_grad = sink.wants_gradients();
    const bool need_hess = sink.wants_hessians();

    Real Mvv_stack[kMaxStackYZ];
    Real Mdv_stack[kMaxStackYZ];
    Real Mvd_stack[kMaxStackYZ];
    Real Md2v_stack[kMaxStackYZ];
    Real Mvd2_stack[kMaxStackYZ];
    Real Mdd_stack[kMaxStackYZ];

    Real* Mvv;
    Real* Mdv;
    Real* Mvd;
    Real* Md2v;
    Real* Mvd2;
    Real* Mdd;
    if (nyz <= kMaxStackYZ) {
        Mvv = Mvv_stack;
        Mdv = Mdv_stack;
        Mvd = Mvd_stack;
        Md2v = Md2v_stack;
        Mvd2 = Mvd2_stack;
        Mdd = Mdd_stack;
    } else {
        auto& tables = evaluate_scratch().tensor_tables;
        tables.resizeFor(nyz);
        Mvv = tables.vv.data();
        Mdv = tables.dv.data();
        Mvd = tables.vd.data();
        Md2v = tables.d2v.data();
        Mvd2 = tables.vd2.data();
        Mdd = tables.dd.data();
    }

    // M_vv is required by every output (values, ∂ξ, ∂ξ²).
    for (std::size_t j = 0; j < ny; ++j) {
        const Real yv = y_axis.values[j];
        for (std::size_t k = 0; k < nz; ++k) {
            Mvv[j * nz + k] = yv * z_axis.values[k];
        }
    }

    if (need_grad || need_hess) {
        for (std::size_t j = 0; j < ny; ++j) {
            const Real yv = y_axis.values[j];
            const Real yd = y_axis.first[j];
            for (std::size_t k = 0; k < nz; ++k) {
                Mdv[j * nz + k] = yd * z_axis.values[k];
                Mvd[j * nz + k] = yv * z_axis.first[k];
            }
        }
    }

    if (need_hess) {
        for (std::size_t j = 0; j < ny; ++j) {
            const Real yv = y_axis.values[j];
            const Real yd = y_axis.first[j];
            const Real yd2 = y_axis.second[j];
            for (std::size_t k = 0; k < nz; ++k) {
                Md2v[j * nz + k] = yd2 * z_axis.values[k];
                Mvd2[j * nz + k] = yv  * z_axis.second[k];
                Mdd[j * nz + k]  = yd  * z_axis.first[k];
            }
        }
    }

    const std::size_t n_nodes = tensor_indices.size();
    sink.prepare(n_nodes);

    for (std::size_t n = 0; n < n_nodes; ++n) {
        const auto& idx = tensor_indices[n];
        const std::size_t i = idx[0];
        const std::size_t jk = idx[1] * nz + idx[2];

        const Real Lx = x_axis.values[i];

        if (need_values) {
            sink.write_value(n, Lx * Mvv[jk]);
        }

        if (need_grad) {
            const Real dLx = x_axis.first[i];
            sink.write_gradient(n,
                                dLx * Mvv[jk],
                                Lx  * Mdv[jk],
                                Lx  * Mvd[jk]);
        }

        if (need_hess) {
            const Real dLx  = x_axis.first[i];
            const Real d2Lx = x_axis.second[i];
            sink.write_hessian(n,
                               d2Lx * Mvv[jk],
                               Lx   * Md2v[jk],
                               Lx   * Mvd2[jk],
                               dLx  * Mdv[jk],
                               dLx  * Mvd[jk],
                               Lx   * Mdd[jk]);
        }
    }
}

void evaluate_tensor_product_factorized(
    const std::vector<std::array<std::size_t, 3>>& tensor_indices,
    const AxisBasisEvaluations& x_axis,
    const AxisBasisEvaluations& y_axis,
    const AxisBasisEvaluations& z_axis,
    std::vector<Real>* values_out,
    std::vector<Gradient>* gradients_out,
    std::vector<Hessian>* hessians_out) {
    const TensorProductVectorSink sink{values_out, gradients_out, hessians_out};
    evaluate_tensor_product_factorized_impl(tensor_indices, x_axis, y_axis, z_axis, sink);
}

void evaluate_tensor_product_factorized_to(
    const std::vector<std::array<std::size_t, 3>>& tensor_indices,
    const AxisBasisEvaluations& x_axis,
    const AxisBasisEvaluations& y_axis,
    const AxisBasisEvaluations& z_axis,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    const TensorProductRawSink sink{values_out, gradients_out, hessians_out};
    evaluate_tensor_product_factorized_impl(tensor_indices, x_axis, y_axis, z_axis, sink);
}

template <std::size_t Q>
inline void write_tensor_product_value_strided_q(
    std::size_t axis_stride,
    std::size_t nyz,
    std::size_t i,
    std::size_t jk,
    const AxisBatchScratch& x_batch,
    const Real* SVMP_RESTRICT Mvv,
    Real* SVMP_RESTRICT value_row) {
    const std::size_t q_axis = Q * axis_stride;
    const std::size_t slot = Q * nyz + jk;
    value_row[Q] = x_batch.values[q_axis + i] * Mvv[slot];
}

template <std::size_t Q>
inline void write_tensor_product_hessian_strided_q(
    std::size_t axis_stride,
    std::size_t nyz,
    std::size_t i,
    std::size_t jk,
    std::size_t output_stride,
    const AxisBatchScratch& x_batch,
    const Real* SVMP_RESTRICT Mvv,
    const Real* SVMP_RESTRICT Mdv,
    const Real* SVMP_RESTRICT Mvd,
    const Real* SVMP_RESTRICT Md2v,
    const Real* SVMP_RESTRICT Mvd2,
    const Real* SVMP_RESTRICT Mdd,
    Real* SVMP_RESTRICT hess_row) {
    const std::size_t q_axis = Q * axis_stride;
    const std::size_t slot = Q * nyz + jk;
    const Real xv = x_batch.values[q_axis + i];
    const Real xd = x_batch.first[q_axis + i];
    const Real x2 = x_batch.second[q_axis + i];
    const Real hxy = xd * Mdv[slot];
    const Real hxz = xd * Mvd[slot];
    const Real hyz = xv * Mdd[slot];
    hess_row[0u * output_stride + Q] = x2 * Mvv[slot];
    hess_row[4u * output_stride + Q] = xv * Md2v[slot];
    hess_row[8u * output_stride + Q] = xv * Mvd2[slot];
    hess_row[1u * output_stride + Q] = hxy;
    hess_row[3u * output_stride + Q] = hxy;
    hess_row[2u * output_stride + Q] = hxz;
    hess_row[6u * output_stride + Q] = hxz;
    hess_row[5u * output_stride + Q] = hyz;
    hess_row[7u * output_stride + Q] = hyz;
}

template <std::size_t Q>
inline void write_tensor_product_hessian_stride4_q(
    std::size_t axis_stride,
    std::size_t nyz,
    std::size_t i,
    std::size_t jk,
    const AxisBatchScratch& x_batch,
    const Real* SVMP_RESTRICT Mvv,
    const Real* SVMP_RESTRICT Mdv,
    const Real* SVMP_RESTRICT Mvd,
    const Real* SVMP_RESTRICT Md2v,
    const Real* SVMP_RESTRICT Mvd2,
    const Real* SVMP_RESTRICT Mdd,
    Real* SVMP_RESTRICT hess_row) {
    const std::size_t q_axis = Q * axis_stride;
    const std::size_t slot = Q * nyz + jk;
    const Real xv = x_batch.values[q_axis + i];
    const Real xd = x_batch.first[q_axis + i];
    const Real x2 = x_batch.second[q_axis + i];
    const Real hxy = xd * Mdv[slot];
    const Real hxz = xd * Mvd[slot];
    const Real hyz = xv * Mdd[slot];
    hess_row[Q] = x2 * Mvv[slot];
    hess_row[16u + Q] = xv * Md2v[slot];
    hess_row[32u + Q] = xv * Mvd2[slot];
    hess_row[4u + Q] = hxy;
    hess_row[12u + Q] = hxy;
    hess_row[8u + Q] = hxz;
    hess_row[24u + Q] = hxz;
    hess_row[20u + Q] = hyz;
    hess_row[28u + Q] = hyz;
}

template <std::size_t Q>
inline void write_tensor_product_gradient_strided_q(
    std::size_t axis_stride,
    std::size_t nyz,
    std::size_t i,
    std::size_t jk,
    std::size_t output_stride,
    const AxisBatchScratch& x_batch,
    const Real* SVMP_RESTRICT Mvv,
    const Real* SVMP_RESTRICT Mdv,
    const Real* SVMP_RESTRICT Mvd,
    Real* SVMP_RESTRICT grad_row) {
    const std::size_t q_axis = Q * axis_stride;
    const std::size_t slot = Q * nyz + jk;
    const Real xv = x_batch.values[q_axis + i];
    const Real xd = x_batch.first[q_axis + i];
    grad_row[0u * output_stride + Q] = xd * Mvv[slot];
    grad_row[1u * output_stride + Q] = xv * Mdv[slot];
    grad_row[2u * output_stride + Q] = xv * Mvd[slot];
}

template <std::size_t Q>
inline void write_tensor_product_gradient_stride4_q(
    std::size_t axis_stride,
    std::size_t nyz,
    std::size_t i,
    std::size_t jk,
    const AxisBatchScratch& x_batch,
    const Real* SVMP_RESTRICT Mvv,
    const Real* SVMP_RESTRICT Mdv,
    const Real* SVMP_RESTRICT Mvd,
    Real* SVMP_RESTRICT grad_row) {
    const std::size_t q_axis = Q * axis_stride;
    const std::size_t slot = Q * nyz + jk;
    const Real xv = x_batch.values[q_axis + i];
    const Real xd = x_batch.first[q_axis + i];
    grad_row[Q] = xd * Mvv[slot];
    grad_row[4u + Q] = xv * Mdv[slot];
    grad_row[8u + Q] = xv * Mvd[slot];
}

template <std::size_t Q>
inline void write_tensor_product_all_strided_q(
    std::size_t axis_stride,
    std::size_t nyz,
    std::size_t i,
    std::size_t jk,
    std::size_t output_stride,
    const AxisBatchScratch& x_batch,
    const Real* SVMP_RESTRICT Mvv,
    const Real* SVMP_RESTRICT Mdv,
    const Real* SVMP_RESTRICT Mvd,
    const Real* SVMP_RESTRICT Md2v,
    const Real* SVMP_RESTRICT Mvd2,
    const Real* SVMP_RESTRICT Mdd,
    Real* SVMP_RESTRICT value_row,
    Real* SVMP_RESTRICT grad_row,
    Real* SVMP_RESTRICT hess_row) {
    const std::size_t q_axis = Q * axis_stride;
    const std::size_t slot = Q * nyz + jk;
    const Real xv = x_batch.values[q_axis + i];
    const Real xd = x_batch.first[q_axis + i];
    value_row[Q] = xv * Mvv[slot];
    grad_row[0u * output_stride + Q] = xd * Mvv[slot];
    grad_row[1u * output_stride + Q] = xv * Mdv[slot];
    grad_row[2u * output_stride + Q] = xv * Mvd[slot];

    const Real x2 = x_batch.second[q_axis + i];
    const Real hxy = xd * Mdv[slot];
    const Real hxz = xd * Mvd[slot];
    const Real hyz = xv * Mdd[slot];
    hess_row[0u * output_stride + Q] = x2 * Mvv[slot];
    hess_row[4u * output_stride + Q] = xv * Md2v[slot];
    hess_row[8u * output_stride + Q] = xv * Mvd2[slot];
    hess_row[1u * output_stride + Q] = hxy;
    hess_row[3u * output_stride + Q] = hxy;
    hess_row[2u * output_stride + Q] = hxz;
    hess_row[6u * output_stride + Q] = hxz;
    hess_row[5u * output_stride + Q] = hyz;
    hess_row[7u * output_stride + Q] = hyz;
}

template <std::size_t Q>
inline void write_tensor_product_all_stride4_q(
    std::size_t axis_stride,
    std::size_t nyz,
    std::size_t i,
    std::size_t jk,
    const AxisBatchScratch& x_batch,
    const Real* SVMP_RESTRICT Mvv,
    const Real* SVMP_RESTRICT Mdv,
    const Real* SVMP_RESTRICT Mvd,
    const Real* SVMP_RESTRICT Md2v,
    const Real* SVMP_RESTRICT Mvd2,
    const Real* SVMP_RESTRICT Mdd,
    Real* SVMP_RESTRICT value_row,
    Real* SVMP_RESTRICT grad_row,
    Real* SVMP_RESTRICT hess_row) {
    const std::size_t q_axis = Q * axis_stride;
    const std::size_t slot = Q * nyz + jk;
    const Real xv = x_batch.values[q_axis + i];
    const Real xd = x_batch.first[q_axis + i];
    value_row[Q] = xv * Mvv[slot];
    grad_row[Q] = xd * Mvv[slot];
    grad_row[4u + Q] = xv * Mdv[slot];
    grad_row[8u + Q] = xv * Mvd[slot];

    const Real x2 = x_batch.second[q_axis + i];
    const Real hxy = xd * Mdv[slot];
    const Real hxz = xd * Mvd[slot];
    const Real hyz = xv * Mdd[slot];
    hess_row[Q] = x2 * Mvv[slot];
    hess_row[16u + Q] = xv * Md2v[slot];
    hess_row[32u + Q] = xv * Mvd2[slot];
    hess_row[4u + Q] = hxy;
    hess_row[12u + Q] = hxy;
    hess_row[8u + Q] = hxz;
    hess_row[24u + Q] = hxz;
    hess_row[20u + Q] = hyz;
    hess_row[28u + Q] = hyz;
}

SVMP_LAGRANGE_NOINLINE SVMP_LAGRANGE_ALIGN64 bool
evaluate_tensor_product_values_stride4_q4_transposed(
    const std::vector<std::array<std::size_t, 3>>& tensor_indices,
    std::size_t axis_stride,
    const AxisBatchScratch& x_batch,
    const AxisBatchScratch& y_batch,
    const AxisBatchScratch& z_batch,
    Real* SVMP_RESTRICT values_out) {
    const std::size_t nyz = axis_stride * axis_stride;
    const std::size_t table_count = 4u * nyz;
    if (table_count > kMaxStackYZ || values_out == nullptr) {
        return false;
    }

    Real Mvv_stack[kMaxStackYZ];
    for (std::size_t j = 0; j < axis_stride; ++j) {
        const Real yv0 = y_batch.values[j];
        const Real yv1 = y_batch.values[axis_stride + j];
        const Real yv2 = y_batch.values[2u * axis_stride + j];
        const Real yv3 = y_batch.values[3u * axis_stride + j];
        for (std::size_t k = 0; k < axis_stride; ++k) {
            const std::size_t base = (j * axis_stride + k) * 4u;
            Mvv_stack[base + 0u] = yv0 * z_batch.values[k];
            Mvv_stack[base + 1u] = yv1 * z_batch.values[axis_stride + k];
            Mvv_stack[base + 2u] = yv2 * z_batch.values[2u * axis_stride + k];
            Mvv_stack[base + 3u] = yv3 * z_batch.values[3u * axis_stride + k];
        }
    }

    for (std::size_t node = 0; node < tensor_indices.size(); ++node) {
        const auto& idx = tensor_indices[node];
        const std::size_t i = idx[0];
        const std::size_t jk = (idx[1] * axis_stride + idx[2]) * 4u;
        Real* SVMP_RESTRICT value_row = values_out + node * 4u;
        value_row[0u] = x_batch.values[i] * Mvv_stack[jk + 0u];
        value_row[1u] = x_batch.values[axis_stride + i] * Mvv_stack[jk + 1u];
        value_row[2u] = x_batch.values[2u * axis_stride + i] * Mvv_stack[jk + 2u];
        value_row[3u] = x_batch.values[3u * axis_stride + i] * Mvv_stack[jk + 3u];
    }

    return true;
}

SVMP_LAGRANGE_NOINLINE SVMP_LAGRANGE_ALIGN64 bool
evaluate_tensor_product_gradients_stride4_q4_transposed(
    const std::vector<std::array<std::size_t, 3>>& tensor_indices,
    std::size_t axis_stride,
    const AxisBatchScratch& x_batch,
    const AxisBatchScratch& y_batch,
    const AxisBatchScratch& z_batch,
    Real* SVMP_RESTRICT gradients_out) {
    const std::size_t nyz = axis_stride * axis_stride;
    const std::size_t table_count = 4u * nyz;
    if (table_count > kMaxStackYZ || gradients_out == nullptr) {
        return false;
    }

    Real Mvv_stack[kMaxStackYZ];
    Real Mdv_stack[kMaxStackYZ];
    Real Mvd_stack[kMaxStackYZ];
    for (std::size_t j = 0; j < axis_stride; ++j) {
        const Real yv0 = y_batch.values[j];
        const Real yv1 = y_batch.values[axis_stride + j];
        const Real yv2 = y_batch.values[2u * axis_stride + j];
        const Real yv3 = y_batch.values[3u * axis_stride + j];
        const Real yd0 = y_batch.first[j];
        const Real yd1 = y_batch.first[axis_stride + j];
        const Real yd2 = y_batch.first[2u * axis_stride + j];
        const Real yd3 = y_batch.first[3u * axis_stride + j];
        for (std::size_t k = 0; k < axis_stride; ++k) {
            const std::size_t base = (j * axis_stride + k) * 4u;
            const Real zv0 = z_batch.values[k];
            const Real zv1 = z_batch.values[axis_stride + k];
            const Real zv2 = z_batch.values[2u * axis_stride + k];
            const Real zv3 = z_batch.values[3u * axis_stride + k];
            const Real zd0 = z_batch.first[k];
            const Real zd1 = z_batch.first[axis_stride + k];
            const Real zd2 = z_batch.first[2u * axis_stride + k];
            const Real zd3 = z_batch.first[3u * axis_stride + k];

            Mvv_stack[base + 0u] = yv0 * zv0;
            Mvv_stack[base + 1u] = yv1 * zv1;
            Mvv_stack[base + 2u] = yv2 * zv2;
            Mvv_stack[base + 3u] = yv3 * zv3;
            Mdv_stack[base + 0u] = yd0 * zv0;
            Mdv_stack[base + 1u] = yd1 * zv1;
            Mdv_stack[base + 2u] = yd2 * zv2;
            Mdv_stack[base + 3u] = yd3 * zv3;
            Mvd_stack[base + 0u] = yv0 * zd0;
            Mvd_stack[base + 1u] = yv1 * zd1;
            Mvd_stack[base + 2u] = yv2 * zd2;
            Mvd_stack[base + 3u] = yv3 * zd3;
        }
    }

    for (std::size_t node = 0; node < tensor_indices.size(); ++node) {
        const auto& idx = tensor_indices[node];
        const std::size_t i = idx[0];
        const std::size_t jk = (idx[1] * axis_stride + idx[2]) * 4u;

        const Real xv0 = x_batch.values[i];
        const Real xv1 = x_batch.values[axis_stride + i];
        const Real xv2 = x_batch.values[2u * axis_stride + i];
        const Real xv3 = x_batch.values[3u * axis_stride + i];
        const Real xd0 = x_batch.first[i];
        const Real xd1 = x_batch.first[axis_stride + i];
        const Real xd2 = x_batch.first[2u * axis_stride + i];
        const Real xd3 = x_batch.first[3u * axis_stride + i];

        Real* SVMP_RESTRICT grad_row = gradients_out + node * 12u;
        grad_row[0u] = xd0 * Mvv_stack[jk + 0u];
        grad_row[1u] = xd1 * Mvv_stack[jk + 1u];
        grad_row[2u] = xd2 * Mvv_stack[jk + 2u];
        grad_row[3u] = xd3 * Mvv_stack[jk + 3u];
        grad_row[4u] = xv0 * Mdv_stack[jk + 0u];
        grad_row[5u] = xv1 * Mdv_stack[jk + 1u];
        grad_row[6u] = xv2 * Mdv_stack[jk + 2u];
        grad_row[7u] = xv3 * Mdv_stack[jk + 3u];
        grad_row[8u] = xv0 * Mvd_stack[jk + 0u];
        grad_row[9u] = xv1 * Mvd_stack[jk + 1u];
        grad_row[10u] = xv2 * Mvd_stack[jk + 2u];
        grad_row[11u] = xv3 * Mvd_stack[jk + 3u];
    }

    return true;
}

template<bool NeedAllOutputs>
SVMP_LAGRANGE_NOINLINE SVMP_LAGRANGE_ALIGN64 bool
evaluate_tensor_product_second_stride4_q4_transposed(
    const std::vector<std::array<std::size_t, 3>>& tensor_indices,
    std::size_t axis_stride,
    const AxisBatchScratch& x_batch,
    const AxisBatchScratch& y_batch,
    const AxisBatchScratch& z_batch,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    const std::size_t nyz = axis_stride * axis_stride;
    const std::size_t table_count = 4u * nyz;
    if (table_count > kMaxStackYZ || hessians_out == nullptr) {
        return false;
    }
    if constexpr (NeedAllOutputs) {
        if (values_out == nullptr || gradients_out == nullptr) {
            return false;
        }
    }

    Real Mvv_stack[kMaxStackYZ];
    Real Mdv_stack[kMaxStackYZ];
    Real Mvd_stack[kMaxStackYZ];
    Real Md2v_stack[kMaxStackYZ];
    Real Mvd2_stack[kMaxStackYZ];
    Real Mdd_stack[kMaxStackYZ];

    for (std::size_t j = 0; j < axis_stride; ++j) {
        const Real yv0 = y_batch.values[j];
        const Real yv1 = y_batch.values[axis_stride + j];
        const Real yv2 = y_batch.values[2u * axis_stride + j];
        const Real yv3 = y_batch.values[3u * axis_stride + j];
        const Real yd0 = y_batch.first[j];
        const Real yd1 = y_batch.first[axis_stride + j];
        const Real yd2 = y_batch.first[2u * axis_stride + j];
        const Real yd3 = y_batch.first[3u * axis_stride + j];
        const Real y20 = y_batch.second[j];
        const Real y21 = y_batch.second[axis_stride + j];
        const Real y22 = y_batch.second[2u * axis_stride + j];
        const Real y23 = y_batch.second[3u * axis_stride + j];

        for (std::size_t k = 0; k < axis_stride; ++k) {
            const std::size_t base = (j * axis_stride + k) * 4u;
            const Real zv0 = z_batch.values[k];
            const Real zv1 = z_batch.values[axis_stride + k];
            const Real zv2 = z_batch.values[2u * axis_stride + k];
            const Real zv3 = z_batch.values[3u * axis_stride + k];
            const Real zd0 = z_batch.first[k];
            const Real zd1 = z_batch.first[axis_stride + k];
            const Real zd2 = z_batch.first[2u * axis_stride + k];
            const Real zd3 = z_batch.first[3u * axis_stride + k];
            const Real z20 = z_batch.second[k];
            const Real z21 = z_batch.second[axis_stride + k];
            const Real z22 = z_batch.second[2u * axis_stride + k];
            const Real z23 = z_batch.second[3u * axis_stride + k];

            Mvv_stack[base + 0u] = yv0 * zv0;
            Mvv_stack[base + 1u] = yv1 * zv1;
            Mvv_stack[base + 2u] = yv2 * zv2;
            Mvv_stack[base + 3u] = yv3 * zv3;
            Mdv_stack[base + 0u] = yd0 * zv0;
            Mdv_stack[base + 1u] = yd1 * zv1;
            Mdv_stack[base + 2u] = yd2 * zv2;
            Mdv_stack[base + 3u] = yd3 * zv3;
            Mvd_stack[base + 0u] = yv0 * zd0;
            Mvd_stack[base + 1u] = yv1 * zd1;
            Mvd_stack[base + 2u] = yv2 * zd2;
            Mvd_stack[base + 3u] = yv3 * zd3;
            Md2v_stack[base + 0u] = y20 * zv0;
            Md2v_stack[base + 1u] = y21 * zv1;
            Md2v_stack[base + 2u] = y22 * zv2;
            Md2v_stack[base + 3u] = y23 * zv3;
            Mvd2_stack[base + 0u] = yv0 * z20;
            Mvd2_stack[base + 1u] = yv1 * z21;
            Mvd2_stack[base + 2u] = yv2 * z22;
            Mvd2_stack[base + 3u] = yv3 * z23;
            Mdd_stack[base + 0u] = yd0 * zd0;
            Mdd_stack[base + 1u] = yd1 * zd1;
            Mdd_stack[base + 2u] = yd2 * zd2;
            Mdd_stack[base + 3u] = yd3 * zd3;
        }
    }

    for (std::size_t node = 0; node < tensor_indices.size(); ++node) {
        const auto& idx = tensor_indices[node];
        const std::size_t i = idx[0];
        const std::size_t jk = (idx[1] * axis_stride + idx[2]) * 4u;

        const Real xv0 = x_batch.values[i];
        const Real xv1 = x_batch.values[axis_stride + i];
        const Real xv2 = x_batch.values[2u * axis_stride + i];
        const Real xv3 = x_batch.values[3u * axis_stride + i];
        const Real xd0 = x_batch.first[i];
        const Real xd1 = x_batch.first[axis_stride + i];
        const Real xd2 = x_batch.first[2u * axis_stride + i];
        const Real xd3 = x_batch.first[3u * axis_stride + i];
        const Real x20 = x_batch.second[i];
        const Real x21 = x_batch.second[axis_stride + i];
        const Real x22 = x_batch.second[2u * axis_stride + i];
        const Real x23 = x_batch.second[3u * axis_stride + i];

        const Real mvv0 = Mvv_stack[jk + 0u];
        const Real mvv1 = Mvv_stack[jk + 1u];
        const Real mvv2 = Mvv_stack[jk + 2u];
        const Real mvv3 = Mvv_stack[jk + 3u];
        const Real mdv0 = Mdv_stack[jk + 0u];
        const Real mdv1 = Mdv_stack[jk + 1u];
        const Real mdv2 = Mdv_stack[jk + 2u];
        const Real mdv3 = Mdv_stack[jk + 3u];
        const Real mvd0 = Mvd_stack[jk + 0u];
        const Real mvd1 = Mvd_stack[jk + 1u];
        const Real mvd2 = Mvd_stack[jk + 2u];
        const Real mvd3 = Mvd_stack[jk + 3u];
        const Real md2v0 = Md2v_stack[jk + 0u];
        const Real md2v1 = Md2v_stack[jk + 1u];
        const Real md2v2 = Md2v_stack[jk + 2u];
        const Real md2v3 = Md2v_stack[jk + 3u];
        const Real mvd20 = Mvd2_stack[jk + 0u];
        const Real mvd21 = Mvd2_stack[jk + 1u];
        const Real mvd22 = Mvd2_stack[jk + 2u];
        const Real mvd23 = Mvd2_stack[jk + 3u];
        const Real mdd0 = Mdd_stack[jk + 0u];
        const Real mdd1 = Mdd_stack[jk + 1u];
        const Real mdd2 = Mdd_stack[jk + 2u];
        const Real mdd3 = Mdd_stack[jk + 3u];

        if constexpr (NeedAllOutputs) {
            Real* SVMP_RESTRICT value_row = values_out + node * 4u;
            value_row[0u] = xv0 * mvv0;
            value_row[1u] = xv1 * mvv1;
            value_row[2u] = xv2 * mvv2;
            value_row[3u] = xv3 * mvv3;

            Real* SVMP_RESTRICT grad_row = gradients_out + node * 12u;
            grad_row[0u] = xd0 * mvv0;
            grad_row[1u] = xd1 * mvv1;
            grad_row[2u] = xd2 * mvv2;
            grad_row[3u] = xd3 * mvv3;
            grad_row[4u] = xv0 * mdv0;
            grad_row[5u] = xv1 * mdv1;
            grad_row[6u] = xv2 * mdv2;
            grad_row[7u] = xv3 * mdv3;
            grad_row[8u] = xv0 * mvd0;
            grad_row[9u] = xv1 * mvd1;
            grad_row[10u] = xv2 * mvd2;
            grad_row[11u] = xv3 * mvd3;
        }

        const Real hxy0 = xd0 * mdv0;
        const Real hxy1 = xd1 * mdv1;
        const Real hxy2 = xd2 * mdv2;
        const Real hxy3 = xd3 * mdv3;
        const Real hxz0 = xd0 * mvd0;
        const Real hxz1 = xd1 * mvd1;
        const Real hxz2 = xd2 * mvd2;
        const Real hxz3 = xd3 * mvd3;
        const Real hyz0 = xv0 * mdd0;
        const Real hyz1 = xv1 * mdd1;
        const Real hyz2 = xv2 * mdd2;
        const Real hyz3 = xv3 * mdd3;

        Real* SVMP_RESTRICT hess_row = hessians_out + node * 36u;
        hess_row[0u] = x20 * mvv0;
        hess_row[1u] = x21 * mvv1;
        hess_row[2u] = x22 * mvv2;
        hess_row[3u] = x23 * mvv3;
        hess_row[4u] = hxy0;
        hess_row[5u] = hxy1;
        hess_row[6u] = hxy2;
        hess_row[7u] = hxy3;
        hess_row[8u] = hxz0;
        hess_row[9u] = hxz1;
        hess_row[10u] = hxz2;
        hess_row[11u] = hxz3;
        hess_row[12u] = hxy0;
        hess_row[13u] = hxy1;
        hess_row[14u] = hxy2;
        hess_row[15u] = hxy3;
        hess_row[16u] = xv0 * md2v0;
        hess_row[17u] = xv1 * md2v1;
        hess_row[18u] = xv2 * md2v2;
        hess_row[19u] = xv3 * md2v3;
        hess_row[20u] = hyz0;
        hess_row[21u] = hyz1;
        hess_row[22u] = hyz2;
        hess_row[23u] = hyz3;
        hess_row[24u] = hxz0;
        hess_row[25u] = hxz1;
        hess_row[26u] = hxz2;
        hess_row[27u] = hxz3;
        hess_row[28u] = hyz0;
        hess_row[29u] = hyz1;
        hess_row[30u] = hyz2;
        hess_row[31u] = hyz3;
        hess_row[32u] = xv0 * mvd20;
        hess_row[33u] = xv1 * mvd21;
        hess_row[34u] = xv2 * mvd22;
        hess_row[35u] = xv3 * mvd23;
    }

    return true;
}

template<int N>
constexpr std::size_t line_public_axis_index(std::size_t node) noexcept {
    return node == 0u ? 0u : (node == 1u ? static_cast<std::size_t>(N - 1) : node - 1u);
}

template<int N>
constexpr std::array<Real, N> make_line_axis_inv_denoms() noexcept {
    std::array<Real, N> inv_denoms{};
    for (int i = 0; i < N; ++i) {
        Real denom = Real(1);
        for (int j = 0; j < N; ++j) {
            if (j != i) {
                denom *= static_cast<Real>(i - j);
            }
        }
        inv_denoms[static_cast<std::size_t>(i)] = Real(1) / denom;
    }
    return inv_denoms;
}

template<int N>
void fill_line_values_product(Real x, Real* SVMP_RESTRICT values) {
    static constexpr auto inv_denoms = make_line_axis_inv_denoms<N>();
    const Real p = static_cast<Real>(N - 1);
    const Real r = (x + Real(1)) * p * Real(0.5);
    Real prefix[N];
    Real suffix[N];
    prefix[0] = Real(1);
    for (int i = 1; i < N; ++i) {
        prefix[i] = prefix[i - 1] * (r - static_cast<Real>(i - 1));
    }
    suffix[N - 1] = Real(1);
    for (int i = N - 2; i >= 0; --i) {
        suffix[i] = suffix[i + 1] * (r - static_cast<Real>(i + 1));
    }
    for (int i = 0; i < N; ++i) {
        const std::size_t slot = static_cast<std::size_t>(i);
        values[slot] = prefix[i] * suffix[i] * inv_denoms[slot];
    }
}

template<int N>
void fill_line_values_product_derivatives(Real x,
                                          Real* SVMP_RESTRICT values,
                                          Real* SVMP_RESTRICT first,
                                          Real* SVMP_RESTRICT second) {
    static constexpr auto inv_denoms = make_line_axis_inv_denoms<N>();
    const Real p = static_cast<Real>(N - 1);
    const Real drdx = p * Real(0.5);
    const Real d2rdx2 = drdx * drdx;
    const Real r = (x + Real(1)) * drdx;

    Real prefix[N + 1];
    Real prefix_d1[N + 1];
    Real prefix_d2[N + 1];
    Real suffix[N + 1];
    Real suffix_d1[N + 1];
    Real suffix_d2[N + 1];

    const bool need_second = second != nullptr;

    prefix[0] = Real(1);
    prefix_d1[0] = Real(0);
    if (need_second) {
        prefix_d2[0] = Real(0);
    }
    for (int i = 0; i < N; ++i) {
        const Real factor = r - static_cast<Real>(i);
        prefix[i + 1] = prefix[i] * factor;
        prefix_d1[i + 1] = prefix_d1[i] * factor + prefix[i];
        if (need_second) {
            prefix_d2[i + 1] = prefix_d2[i] * factor + Real(2) * prefix_d1[i];
        }
    }

    suffix[N] = Real(1);
    suffix_d1[N] = Real(0);
    if (need_second) {
        suffix_d2[N] = Real(0);
    }
    for (int i = N - 1; i >= 0; --i) {
        const Real factor = r - static_cast<Real>(i);
        suffix[i] = suffix[i + 1] * factor;
        suffix_d1[i] = suffix_d1[i + 1] * factor + suffix[i + 1];
        if (need_second) {
            suffix_d2[i] = suffix_d2[i + 1] * factor + Real(2) * suffix_d1[i + 1];
        }
    }

    for (int i = 0; i < N; ++i) {
        const std::size_t slot = static_cast<std::size_t>(i);
        const Real inv = inv_denoms[slot];
        const Real pre = prefix[i];
        const Real suf = suffix[i + 1];
        const Real pre_d1 = prefix_d1[i];
        const Real suf_d1 = suffix_d1[i + 1];
        values[slot] = pre * suf * inv;
        if (first != nullptr) {
            first[slot] = (pre_d1 * suf + pre * suf_d1) * inv * drdx;
        }
        if (second != nullptr) {
            const Real d2 =
                prefix_d2[i] * suf +
                Real(2) * pre_d1 * suf_d1 +
                pre * suffix_d2[i + 1];
            second[slot] = d2 * inv * d2rdx2;
        }
    }
}

template<int N>
void fill_axis_batch_product_q4(
    AxisBatchScratch& scratch,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t component,
    AxisDeriv level) {
    constexpr std::size_t axis_stride = static_cast<std::size_t>(N);
    scratch.resizeFor(4u * axis_stride, level);
    for (std::size_t q = 0; q < 4u; ++q) {
        Real* values = scratch.values.data() + q * axis_stride;
        if (level == AxisDeriv::ValuesOnly) {
            fill_line_values_product<N>(points[q][component], values);
        } else {
            Real* first = scratch.first.data() + q * axis_stride;
            Real* second = level == AxisDeriv::ValuesAndFirstAndSecond
                ? scratch.second.data() + q * axis_stride
                : nullptr;
            fill_line_values_product_derivatives<N>(
                points[q][component], values, first, second);
        }
    }
}

bool try_fill_axis_batch_product_q4(
    AxisBatchScratch& scratch,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t component,
    int n_axis,
    AxisDeriv level) {
    switch (n_axis) {
        case 5:
            fill_axis_batch_product_q4<5>(scratch, points, component, level);
            return true;
        case 6:
            fill_axis_batch_product_q4<6>(scratch, points, component, level);
            return true;
        case 7:
            fill_axis_batch_product_q4<7>(scratch, points, component, level);
            return true;
        case 8:
            fill_axis_batch_product_q4<8>(scratch, points, component, level);
            return true;
        case 9:
            fill_axis_batch_product_q4<9>(scratch, points, component, level);
            return true;
        default:
            return false;
    }
}

template<int N>
SVMP_LAGRANGE_NOINLINE SVMP_LAGRANGE_ALIGN64 void evaluate_line_values_product_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    Real q0[N];
    Real q1[N];
    Real q2[N];
    Real q3[N];
    fill_line_values_product<N>(points[0][0], q0);
    fill_line_values_product<N>(points[1][0], q1);
    fill_line_values_product<N>(points[2][0], q2);
    fill_line_values_product<N>(points[3][0], q3);

    for (std::size_t node = 0; node < static_cast<std::size_t>(N); ++node) {
        const std::size_t i = line_public_axis_index<N>(node);
        Real* row = values_out + node * output_stride;
        row[0] = q0[i];
        row[1] = q1[i];
        row[2] = q2[i];
        row[3] = q3[i];
    }
}

FE_ALWAYS_INLINE void write_line_order4_values_q(
    Real x,
    std::size_t q,
    Real* SVMP_RESTRICT row0,
    Real* SVMP_RESTRICT row1,
    Real* SVMP_RESTRICT row2,
    Real* SVMP_RESTRICT row3,
    Real* SVMP_RESTRICT row4) {
    const Real r = (x + Real(1)) * Real(2);
    const Real f0 = r;
    const Real f1 = r - Real(1);
    const Real f2 = r - Real(2);
    const Real f3 = r - Real(3);
    const Real f4 = r - Real(4);
    const Real f01 = f0 * f1;
    const Real f12 = f1 * f2;
    const Real f23 = f2 * f3;
    const Real f34 = f3 * f4;
    const Real v0 = (f12 * f34) / Real(24);
    const Real v1 = -(f0 * f2 * f34) / Real(6);
    const Real v2 = (f01 * f34) / Real(4);
    const Real v3 = -(f01 * f2 * f4) / Real(6);
    const Real v4 = (f01 * f23) / Real(24);
    row0[q] = v0;
    row1[q] = v4;
    row2[q] = v1;
    row3[q] = v2;
    row4[q] = v3;
}

SVMP_LAGRANGE_NOINLINE SVMP_LAGRANGE_ALIGN64 void evaluate_line_order4_values_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    Real* row0 = values_out + 0u * output_stride;
    Real* row1 = values_out + 1u * output_stride;
    Real* row2 = values_out + 2u * output_stride;
    Real* row3 = values_out + 3u * output_stride;
    Real* row4 = values_out + 4u * output_stride;
    write_line_order4_values_q(points[0][0], 0u, row0, row1, row2, row3, row4);
    write_line_order4_values_q(points[1][0], 1u, row0, row1, row2, row3, row4);
    write_line_order4_values_q(points[2][0], 2u, row0, row1, row2, row3, row4);
    write_line_order4_values_q(points[3][0], 3u, row0, row1, row2, row3, row4);
}

SVMP_LAGRANGE_NOINLINE void evaluate_triangle_order1_gradients_strided(
    std::size_t num_qpts,
    std::size_t output_stride,
    Real* SVMP_RESTRICT gradients_out) {
    Real* SVMP_RESTRICT row0 = gradients_out + 0u * 3u * output_stride;
    Real* SVMP_RESTRICT row1 = gradients_out + 1u * 3u * output_stride;
    Real* SVMP_RESTRICT row2 = gradients_out + 2u * 3u * output_stride;

    for (std::size_t q = 0; q < num_qpts; ++q) {
        row0[0u * output_stride + q] = Real(-1);
        row0[1u * output_stride + q] = Real(-1);
        row0[2u * output_stride + q] = Real(0);
        row1[0u * output_stride + q] = Real(1);
        row1[1u * output_stride + q] = Real(0);
        row1[2u * output_stride + q] = Real(0);
        row2[0u * output_stride + q] = Real(0);
        row2[1u * output_stride + q] = Real(1);
        row2[2u * output_stride + q] = Real(0);
    }
}

template<int N>
SVMP_LAGRANGE_NOINLINE void evaluate_line_hessians_product_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT hessians_out) {
    Real values[4][N];
    Real second[4][N];
    for (std::size_t q = 0; q < 4u; ++q) {
        fill_line_values_product_derivatives<N>(
            points[q][0], values[q], nullptr, second[q]);
    }
    for (std::size_t node = 0; node < static_cast<std::size_t>(N); ++node) {
        const std::size_t i = line_public_axis_index<N>(node);
        write_line_hessian_q4_row(hessians_out + node * 9u * output_stride,
                                  output_stride,
                                  second[0][i], second[1][i],
                                  second[2][i], second[3][i]);
    }
}

template<int N>
SVMP_LAGRANGE_NOINLINE void evaluate_line_all_product_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    Real values[4][N];
    Real first[4][N];
    Real second[4][N];
    for (std::size_t q = 0; q < 4u; ++q) {
        fill_line_values_product_derivatives<N>(
            points[q][0], values[q], first[q], second[q]);
    }
    for (std::size_t node = 0; node < static_cast<std::size_t>(N); ++node) {
        const std::size_t i = line_public_axis_index<N>(node);
        Real* value_row = values_out + node * output_stride;
        value_row[0] = values[0][i];
        value_row[1] = values[1][i];
        value_row[2] = values[2][i];
        value_row[3] = values[3][i];
        write_line_gradient_q4_row(gradients_out + node * 3u * output_stride,
                                   output_stride,
                                   first[0][i], first[1][i],
                                   first[2][i], first[3][i]);
        write_line_hessian_q4_row(hessians_out + node * 9u * output_stride,
                                  output_stride,
                                  second[0][i], second[1][i],
                                  second[2][i], second[3][i]);
    }
}

inline void write_quad_product_value_row_q4(
    Real* SVMP_RESTRICT row,
    const Real* SVMP_RESTRICT x0,
    const Real* SVMP_RESTRICT x1,
    const Real* SVMP_RESTRICT x2,
    const Real* SVMP_RESTRICT x3,
    const Real* SVMP_RESTRICT y0,
    const Real* SVMP_RESTRICT y1,
    const Real* SVMP_RESTRICT y2,
    const Real* SVMP_RESTRICT y3,
    std::size_t i,
    std::size_t j) {
    row[0] = x0[i] * y0[j];
    row[1] = x1[i] * y1[j];
    row[2] = x2[i] * y2[j];
    row[3] = x3[i] * y3[j];
}

template<int N>
void evaluate_quad_values_product_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out) {
    Real x0[N];
    Real x1[N];
    Real x2[N];
    Real x3[N];
    Real y0[N];
    Real y1[N];
    Real y2[N];
    Real y3[N];
    fill_line_values_product<N>(points[0][0], x0);
    fill_line_values_product<N>(points[1][0], x1);
    fill_line_values_product<N>(points[2][0], x2);
    fill_line_values_product<N>(points[3][0], x3);
    fill_line_values_product<N>(points[0][1], y0);
    fill_line_values_product<N>(points[1][1], y1);
    fill_line_values_product<N>(points[2][1], y2);
    fill_line_values_product<N>(points[3][1], y3);

    constexpr std::size_t p = static_cast<std::size_t>(N - 1);
    std::size_t node = 0u;
    write_quad_product_value_row_q4(values_out + node++ * output_stride,
                                    x0, x1, x2, x3, y0, y1, y2, y3, 0u, 0u);
    write_quad_product_value_row_q4(values_out + node++ * output_stride,
                                    x0, x1, x2, x3, y0, y1, y2, y3, p, 0u);
    write_quad_product_value_row_q4(values_out + node++ * output_stride,
                                    x0, x1, x2, x3, y0, y1, y2, y3, p, p);
    write_quad_product_value_row_q4(values_out + node++ * output_stride,
                                    x0, x1, x2, x3, y0, y1, y2, y3, 0u, p);

    for (std::size_t i = 1u; i < p; ++i) {
        write_quad_product_value_row_q4(values_out + node++ * output_stride,
                                        x0, x1, x2, x3, y0, y1, y2, y3, i, 0u);
    }
    for (std::size_t j = 1u; j < p; ++j) {
        write_quad_product_value_row_q4(values_out + node++ * output_stride,
                                        x0, x1, x2, x3, y0, y1, y2, y3, p, j);
    }
    for (std::size_t i = p - 1u; i > 0u; --i) {
        write_quad_product_value_row_q4(values_out + node++ * output_stride,
                                        x0, x1, x2, x3, y0, y1, y2, y3, i, p);
    }
    for (std::size_t j = p - 1u; j > 0u; --j) {
        write_quad_product_value_row_q4(values_out + node++ * output_stride,
                                        x0, x1, x2, x3, y0, y1, y2, y3, 0u, j);
    }
    for (std::size_t j = 1u; j < p; ++j) {
        for (std::size_t i = 1u; i < p; ++i) {
            write_quad_product_value_row_q4(values_out + node++ * output_stride,
                                            x0, x1, x2, x3, y0, y1, y2, y3, i, j);
        }
    }
}

template<int N>
void evaluate_quad_derivatives_product_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    const bool need_grad = gradients_out != nullptr;
    const bool need_hess = hessians_out != nullptr;
    Real xv[4][N];
    Real xd[4][N];
    Real x2[4][N];
    Real yv[4][N];
    Real yd[4][N];
    Real y2[4][N];

    for (std::size_t q = 0; q < 4u; ++q) {
        fill_line_values_product_derivatives<N>(
            points[q][0], xv[q], (need_grad || need_hess) ? xd[q] : nullptr,
            need_hess ? x2[q] : nullptr);
        fill_line_values_product_derivatives<N>(
            points[q][1], yv[q], (need_grad || need_hess) ? yd[q] : nullptr,
            need_hess ? y2[q] : nullptr);
    }

    constexpr std::size_t p = static_cast<std::size_t>(N - 1);
    std::size_t node = 0u;
    auto write_node = [&](std::size_t i, std::size_t j) {
        Real* value_row = values_out != nullptr ? values_out + node * output_stride : nullptr;
        Real* grad_row = gradients_out != nullptr ? gradients_out + node * 3u * output_stride : nullptr;
        Real* hess_row = hessians_out != nullptr ? hessians_out + node * 9u * output_stride : nullptr;
        if (grad_row != nullptr) {
            grad_row[2u * output_stride + 0u] = Real(0);
            grad_row[2u * output_stride + 1u] = Real(0);
            grad_row[2u * output_stride + 2u] = Real(0);
            grad_row[2u * output_stride + 3u] = Real(0);
        }
        if (hess_row != nullptr) {
            hess_row[2u * output_stride + 0u] = Real(0);
            hess_row[2u * output_stride + 1u] = Real(0);
            hess_row[2u * output_stride + 2u] = Real(0);
            hess_row[2u * output_stride + 3u] = Real(0);
            hess_row[5u * output_stride + 0u] = Real(0);
            hess_row[5u * output_stride + 1u] = Real(0);
            hess_row[5u * output_stride + 2u] = Real(0);
            hess_row[5u * output_stride + 3u] = Real(0);
            hess_row[6u * output_stride + 0u] = Real(0);
            hess_row[6u * output_stride + 1u] = Real(0);
            hess_row[6u * output_stride + 2u] = Real(0);
            hess_row[6u * output_stride + 3u] = Real(0);
            hess_row[7u * output_stride + 0u] = Real(0);
            hess_row[7u * output_stride + 1u] = Real(0);
            hess_row[7u * output_stride + 2u] = Real(0);
            hess_row[7u * output_stride + 3u] = Real(0);
            hess_row[8u * output_stride + 0u] = Real(0);
            hess_row[8u * output_stride + 1u] = Real(0);
            hess_row[8u * output_stride + 2u] = Real(0);
            hess_row[8u * output_stride + 3u] = Real(0);
        }
        for (std::size_t q = 0; q < 4u; ++q) {
            const Real x_value = xv[q][i];
            const Real y_value = yv[q][j];
            if (value_row != nullptr) {
                value_row[q] = x_value * y_value;
            }
            if (grad_row != nullptr) {
                grad_row[0u * output_stride + q] = xd[q][i] * y_value;
                grad_row[1u * output_stride + q] = x_value * yd[q][j];
            }
            if (hess_row != nullptr) {
                const Real hxy = xd[q][i] * yd[q][j];
                hess_row[0u * output_stride + q] = x2[q][i] * y_value;
                hess_row[1u * output_stride + q] = hxy;
                hess_row[3u * output_stride + q] = hxy;
                hess_row[4u * output_stride + q] = x_value * y2[q][j];
            }
        }
        ++node;
    };

    write_node(0u, 0u);
    write_node(p, 0u);
    write_node(p, p);
    write_node(0u, p);
    for (std::size_t i = 1u; i < p; ++i) {
        write_node(i, 0u);
    }
    for (std::size_t j = 1u; j < p; ++j) {
        write_node(p, j);
    }
    for (std::size_t i = p - 1u; i > 0u; --i) {
        write_node(i, p);
    }
    for (std::size_t j = p - 1u; j > 0u; --j) {
        write_node(0u, j);
    }
    for (std::size_t j = 1u; j < p; ++j) {
        for (std::size_t i = 1u; i < p; ++i) {
            write_node(i, j);
        }
    }
}

SVMP_LAGRANGE_NOINLINE SVMP_LAGRANGE_ALIGN64 void evaluate_quad_order8_gradients_product_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT gradients_out) {
    constexpr int N = 9;
    constexpr std::size_t p = 8u;
    Real xv[4][N];
    Real xd[4][N];
    Real yv[4][N];
    Real yd[4][N];

    for (std::size_t q = 0; q < 4u; ++q) {
        fill_line_values_product_derivatives<N>(points[q][0], xv[q], xd[q], nullptr);
        fill_line_values_product_derivatives<N>(points[q][1], yv[q], yd[q], nullptr);
    }

    std::size_t node = 0u;
    auto write_node = [&](std::size_t i, std::size_t j) {
        Real* SVMP_RESTRICT row = gradients_out + node * 3u * output_stride;
        row[0u] = xd[0][i] * yv[0][j];
        row[1u] = xd[1][i] * yv[1][j];
        row[2u] = xd[2][i] * yv[2][j];
        row[3u] = xd[3][i] * yv[3][j];
        row[output_stride + 0u] = xv[0][i] * yd[0][j];
        row[output_stride + 1u] = xv[1][i] * yd[1][j];
        row[output_stride + 2u] = xv[2][i] * yd[2][j];
        row[output_stride + 3u] = xv[3][i] * yd[3][j];
        row[2u * output_stride + 0u] = Real(0);
        row[2u * output_stride + 1u] = Real(0);
        row[2u * output_stride + 2u] = Real(0);
        row[2u * output_stride + 3u] = Real(0);
        ++node;
    };

    write_node(0u, 0u);
    write_node(p, 0u);
    write_node(p, p);
    write_node(0u, p);
    for (std::size_t i = 1u; i < p; ++i) {
        write_node(i, 0u);
    }
    for (std::size_t j = 1u; j < p; ++j) {
        write_node(p, j);
    }
    for (std::size_t i = p - 1u; i > 0u; --i) {
        write_node(i, p);
    }
    for (std::size_t j = p - 1u; j > 0u; --j) {
        write_node(0u, j);
    }
    for (std::size_t j = 1u; j < p; ++j) {
        for (std::size_t i = 1u; i < p; ++i) {
            write_node(i, j);
        }
    }
}

template<int N>
void evaluate_line_gradients_horner_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    const Real* SVMP_RESTRICT d_coeffs,
    Real* SVMP_RESTRICT gradients_out) {
    const Real x0 = points[0][0];
    const Real x1 = points[1][0];
    const Real x2 = points[2][0];
    const Real x3 = points[3][0];

    for (std::size_t node = 0; node < static_cast<std::size_t>(N); ++node) {
        const std::size_t i = line_public_axis_index<N>(node);
        const Real* c = d_coeffs + i * static_cast<std::size_t>(N - 1);
        Real r0 = c[N - 2];
        Real r1 = c[N - 2];
        Real r2 = c[N - 2];
        Real r3 = c[N - 2];
        for (int k = N - 2; k > 0; --k) {
            const Real ck = c[k - 1];
            r0 = r0 * x0 + ck;
            r1 = r1 * x1 + ck;
            r2 = r2 * x2 + ck;
            r3 = r3 * x3 + ck;
        }
        Real* row = gradients_out + node * 3u * output_stride;
        row[0] = r0;
        row[1] = r1;
        row[2] = r2;
        row[3] = r3;
        row[output_stride + 0u] = Real(0);
        row[output_stride + 1u] = Real(0);
        row[output_stride + 2u] = Real(0);
        row[output_stride + 3u] = Real(0);
        row[2u * output_stride + 0u] = Real(0);
        row[2u * output_stride + 1u] = Real(0);
        row[2u * output_stride + 2u] = Real(0);
        row[2u * output_stride + 3u] = Real(0);
    }
}

bool try_evaluate_line_values_horner_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    const Real* SVMP_RESTRICT v_coeffs,
    int n_axis,
    Real* SVMP_RESTRICT values_out) {
    (void)v_coeffs;
    switch (n_axis) {
        case 5:
            evaluate_line_order4_values_q4(points, output_stride, values_out);
            return true;
        case 6:
            evaluate_line_values_product_q4<6>(points, output_stride, values_out);
            return true;
        case 7:
            evaluate_line_values_product_q4<7>(points, output_stride, values_out);
            return true;
        case 8:
            evaluate_line_values_product_q4<8>(points, output_stride, values_out);
            return true;
        case 9:
            evaluate_line_values_product_q4<9>(points, output_stride, values_out);
            return true;
        default:
            return false;
    }
}

bool try_evaluate_line_gradients_horner_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    const Real* SVMP_RESTRICT d_coeffs,
    int n_axis,
    Real* SVMP_RESTRICT gradients_out) {
    switch (n_axis) {
        case 5:
            evaluate_line_gradients_horner_q4<5>(points, output_stride, d_coeffs, gradients_out);
            return true;
        case 6:
            evaluate_line_gradients_horner_q4<6>(points, output_stride, d_coeffs, gradients_out);
            return true;
        case 7:
            evaluate_line_gradients_horner_q4<7>(points, output_stride, d_coeffs, gradients_out);
            return true;
        case 8:
            evaluate_line_gradients_horner_q4<8>(points, output_stride, d_coeffs, gradients_out);
            return true;
        case 9:
            evaluate_line_gradients_horner_q4<9>(points, output_stride, d_coeffs, gradients_out);
            return true;
        default:
            return false;
    }
}

SVMP_LAGRANGE_NOINLINE bool try_evaluate_line_hessians_product_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    int n_axis,
    Real* SVMP_RESTRICT hessians_out) {
    switch (n_axis) {
        case 5:
            evaluate_line_hessians_product_q4<5>(points, output_stride, hessians_out);
            return true;
        case 6:
            evaluate_line_hessians_product_q4<6>(points, output_stride, hessians_out);
            return true;
        case 7:
            evaluate_line_hessians_product_q4<7>(points, output_stride, hessians_out);
            return true;
        case 8:
            evaluate_line_hessians_product_q4<8>(points, output_stride, hessians_out);
            return true;
        case 9:
            evaluate_line_hessians_product_q4<9>(points, output_stride, hessians_out);
            return true;
        default:
            return false;
    }
}

SVMP_LAGRANGE_NOINLINE bool try_evaluate_line_all_product_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    int n_axis,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    switch (n_axis) {
        case 5:
            evaluate_line_all_product_q4<5>(
                points, output_stride, values_out, gradients_out, hessians_out);
            return true;
        case 6:
            evaluate_line_all_product_q4<6>(
                points, output_stride, values_out, gradients_out, hessians_out);
            return true;
        case 7:
            evaluate_line_all_product_q4<7>(
                points, output_stride, values_out, gradients_out, hessians_out);
            return true;
        case 8:
            evaluate_line_all_product_q4<8>(
                points, output_stride, values_out, gradients_out, hessians_out);
            return true;
        case 9:
            evaluate_line_all_product_q4<9>(
                points, output_stride, values_out, gradients_out, hessians_out);
            return true;
        default:
            return false;
    }
}

SVMP_LAGRANGE_NOINLINE bool try_evaluate_quad_values_product_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    int n_axis,
    Real* SVMP_RESTRICT values_out) {
    switch (n_axis) {
        case 5:
            evaluate_quad_values_product_q4<5>(points, output_stride, values_out);
            return true;
        case 6:
            evaluate_quad_values_product_q4<6>(points, output_stride, values_out);
            return true;
        case 7:
            evaluate_quad_values_product_q4<7>(points, output_stride, values_out);
            return true;
        case 8:
            evaluate_quad_values_product_q4<8>(points, output_stride, values_out);
            return true;
        case 9:
            evaluate_quad_values_product_q4<9>(points, output_stride, values_out);
            return true;
        default:
            return false;
    }
}

SVMP_LAGRANGE_NOINLINE bool try_evaluate_quad_derivatives_product_q4(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    int n_axis,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    switch (n_axis) {
        case 5:
            evaluate_quad_derivatives_product_q4<5>(
                points, output_stride, values_out, gradients_out, hessians_out);
            return true;
        case 6:
            evaluate_quad_derivatives_product_q4<6>(
                points, output_stride, values_out, gradients_out, hessians_out);
            return true;
        case 7:
            evaluate_quad_derivatives_product_q4<7>(
                points, output_stride, values_out, gradients_out, hessians_out);
            return true;
        case 8:
            evaluate_quad_derivatives_product_q4<8>(
                points, output_stride, values_out, gradients_out, hessians_out);
            return true;
        case 9:
            evaluate_quad_derivatives_product_q4<9>(
                points, output_stride, values_out, gradients_out, hessians_out);
            return true;
        default:
            return false;
    }
}

void evaluate_tensor_product_points_strided(
    LagrangeTopology topology,
    const std::vector<std::array<std::size_t, 3>>& tensor_indices,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    const Real* v_coeffs,
    const Real* d_coeffs,
    const Real* d2_coeffs,
    const Real* barycentric_weights,
    int n_axis,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    const std::size_t num_qpts = points.size();
    if (num_qpts == 0 || tensor_indices.empty()) {
        return;
    }

    const bool need_grad = gradients_out != nullptr;
    const bool need_hess = hessians_out != nullptr;
    const bool values_only = values_out != nullptr && !need_grad && !need_hess;
    const bool gradients_only = values_out == nullptr && need_grad && !need_hess;
    const bool hessians_only = values_out == nullptr && gradients_out == nullptr && need_hess;
    const bool all_outputs = values_out != nullptr && need_grad && need_hess;
    const AxisDeriv level = need_hess
        ? AxisDeriv::ValuesAndFirstAndSecond
        : (need_grad ? AxisDeriv::ValuesAndFirst : AxisDeriv::ValuesOnly);

    if (topology == LagrangeTopology::Line && num_qpts == 4u) {
        if (values_only &&
            try_evaluate_line_values_horner_q4(
                points, output_stride, v_coeffs, n_axis, values_out)) {
            return;
        }
        if (gradients_only &&
            try_evaluate_line_gradients_horner_q4(
                points, output_stride, d_coeffs, n_axis, gradients_out)) {
            return;
        }
        if (hessians_only &&
            try_evaluate_line_hessians_product_q4(
                points, output_stride, n_axis, hessians_out)) {
            return;
        }
        if (all_outputs &&
            try_evaluate_line_all_product_q4(
                points, output_stride, n_axis, values_out, gradients_out, hessians_out)) {
            return;
        }
    }
    if (topology == LagrangeTopology::Quadrilateral &&
        values_only &&
        num_qpts == 4u &&
        try_evaluate_quad_values_product_q4(points, output_stride, n_axis, values_out)) {
        return;
    }
    if (topology == LagrangeTopology::Quadrilateral &&
        gradients_only &&
        num_qpts == 4u &&
        n_axis == 5) {
        evaluate_quad_order4_gradients_q4(points, output_stride, gradients_out);
        return;
    }
    if (topology == LagrangeTopology::Quadrilateral &&
        gradients_only &&
        num_qpts == 4u &&
        n_axis == 9) {
        evaluate_quad_order8_gradients_product_q4(points, output_stride, gradients_out);
        return;
    }
    if (topology == LagrangeTopology::Quadrilateral &&
        (gradients_only || hessians_only || all_outputs) &&
        num_qpts == 4u &&
        try_evaluate_quad_derivatives_product_q4(
            points, output_stride, n_axis, values_out, gradients_out, hessians_out)) {
        return;
    }

    auto& scratch = evaluate_scratch();
    AxisBatchScratch& x_batch = scratch.axis_x_batch;
    AxisBatchScratch& y_batch = scratch.axis_y_batch;
    AxisBatchScratch& z_batch = scratch.axis_z_batch;

    const bool has_y = topology != LagrangeTopology::Line;
    const bool has_z = topology == LagrangeTopology::Hexahedron;
    const std::size_t axis_stride = static_cast<std::size_t>(n_axis);
    const bool use_product_axis_batch =
        has_z &&
        gradients_only &&
        num_qpts == 4u &&
        n_axis >= 5 &&
        n_axis <= 9;
    auto fill_tensor_axis_batch = [&](AxisBatchScratch& batch, std::size_t component) {
        if (use_product_axis_batch &&
            try_fill_axis_batch_product_q4(batch, points, component, n_axis, level)) {
            return;
        }
        fill_axis_batch(batch, points, component, v_coeffs, d_coeffs, d2_coeffs,
                        barycentric_weights, n_axis, level);
    };

    fill_tensor_axis_batch(x_batch, 0u);
    if (!has_y) {
        if (values_only) {
            if (num_qpts == 4u) {
                for (std::size_t node = 0; node < tensor_indices.size(); ++node) {
                    const std::size_t i = tensor_indices[node][0];
                    Real* value_row = values_out + node * output_stride;
                    value_row[0] = x_batch.values[i];
                    value_row[1] = x_batch.values[axis_stride + i];
                    value_row[2] = x_batch.values[2u * axis_stride + i];
                    value_row[3] = x_batch.values[3u * axis_stride + i];
                }
                return;
            }
            for (std::size_t node = 0; node < tensor_indices.size(); ++node) {
                const std::size_t i = tensor_indices[node][0];
                Real* value_row = values_out + node * output_stride;

                for (std::size_t q = 0; q < num_qpts; ++q) {
                    value_row[q] = x_batch.values[q * axis_stride + i];
                }
            }
            return;
        }

        if (gradients_only) {
            for (std::size_t node = 0; node < tensor_indices.size(); ++node) {
                const std::size_t i = tensor_indices[node][0];
                Real* grad_row = gradients_out + node * 3u * output_stride;

                for (std::size_t q = 0; q < num_qpts; ++q) {
                    grad_row[0u * output_stride + q] =
                        x_batch.first[q * axis_stride + i];
                    grad_row[1u * output_stride + q] = Real(0);
                    grad_row[2u * output_stride + q] = Real(0);
                }
            }
            return;
        }

        for (std::size_t node = 0; node < tensor_indices.size(); ++node) {
            const std::size_t i = tensor_indices[node][0];
            Real* value_row = values_out ? values_out + node * output_stride : nullptr;
            Real* grad_row = gradients_out ? gradients_out + node * 3u * output_stride : nullptr;
            Real* hess_row = hessians_out ? hessians_out + node * 9u * output_stride : nullptr;

            for (std::size_t q = 0; q < num_qpts; ++q) {
                const std::size_t q_axis = q * axis_stride + i;
                if (value_row != nullptr) {
                    value_row[q] = x_batch.values[q_axis];
                }
                if (need_grad) {
                    grad_row[0u * output_stride + q] = x_batch.first[q_axis];
                    grad_row[1u * output_stride + q] = Real(0);
                    grad_row[2u * output_stride + q] = Real(0);
                }
                if (need_hess) {
                    hess_row[0u * output_stride + q] = x_batch.second[q_axis];
                    hess_row[1u * output_stride + q] = Real(0);
                    hess_row[2u * output_stride + q] = Real(0);
                    hess_row[3u * output_stride + q] = Real(0);
                    hess_row[4u * output_stride + q] = Real(0);
                    hess_row[5u * output_stride + q] = Real(0);
                    hess_row[6u * output_stride + q] = Real(0);
                    hess_row[7u * output_stride + q] = Real(0);
                    hess_row[8u * output_stride + q] = Real(0);
                }
            }
        }
        return;
    }
    const bool use_tensor_tables =
        has_z ||
        (axis_stride == 8u && !(need_hess && values_out == nullptr && gradients_out == nullptr));
    if (use_tensor_tables) {
        fill_tensor_axis_batch(y_batch, 1u);
    } else if (has_y) {
        fill_tensor_axis_batch(y_batch, 1u);
    }
    if (has_z) {
        fill_tensor_axis_batch(z_batch, 2u);
    }

    if (use_tensor_tables) {
        const std::size_t ny = axis_stride;
        const std::size_t nz = has_z ? axis_stride : 1u;
        const std::size_t nyz = ny * nz;
        const std::size_t table_count = num_qpts * nyz;

        if (has_z && num_qpts == 4u && output_stride == 4u) {
            if (values_only &&
                evaluate_tensor_product_values_stride4_q4_transposed(
                    tensor_indices, axis_stride, x_batch, y_batch, z_batch, values_out)) {
                return;
            }
            if (gradients_only &&
                evaluate_tensor_product_gradients_stride4_q4_transposed(
                    tensor_indices, axis_stride, x_batch, y_batch, z_batch, gradients_out)) {
                return;
            }
            if (hessians_only &&
                evaluate_tensor_product_second_stride4_q4_transposed<false>(
                    tensor_indices, axis_stride, x_batch, y_batch, z_batch,
                    nullptr, nullptr, hessians_out)) {
                return;
            }
            if (all_outputs &&
                evaluate_tensor_product_second_stride4_q4_transposed<true>(
                    tensor_indices, axis_stride, x_batch, y_batch, z_batch,
                    values_out, gradients_out, hessians_out)) {
                return;
            }
        }

        Real Mvv_stack[kMaxStackYZ];
        Real Mdv_stack[kMaxStackYZ];
        Real Mvd_stack[kMaxStackYZ];
        Real Md2v_stack[kMaxStackYZ];
        Real Mvd2_stack[kMaxStackYZ];
        Real Mdd_stack[kMaxStackYZ];

        Real* Mvv;
        Real* Mdv;
        Real* Mvd;
        Real* Md2v;
        Real* Mvd2;
        Real* Mdd;
        if (table_count <= kMaxStackYZ) {
            Mvv = Mvv_stack;
            Mdv = Mdv_stack;
            Mvd = Mvd_stack;
            Md2v = Md2v_stack;
            Mvd2 = Mvd2_stack;
            Mdd = Mdd_stack;
        } else {
            auto& tables = scratch.tensor_tables;
            tables.resizeFor(table_count);
            Mvv = tables.vv.data();
            Mdv = tables.dv.data();
            Mvd = tables.vd.data();
            Md2v = tables.d2v.data();
            Mvd2 = tables.vd2.data();
            Mdd = tables.dd.data();
        }

        for (std::size_t q = 0; q < num_qpts; ++q) {
            const std::size_t q_axis = q * axis_stride;
            const std::size_t q_table = q * nyz;
            for (std::size_t j = 0; j < ny; ++j) {
                const Real yv = y_batch.values[q_axis + j];
                const Real yd = (need_grad || need_hess) ? y_batch.first[q_axis + j] : Real(0);
                const Real y2 = need_hess ? y_batch.second[q_axis + j] : Real(0);
                for (std::size_t k = 0; k < nz; ++k) {
                    const std::size_t slot = q_table + j * nz + k;
                    const Real zv = has_z ? z_batch.values[q_axis + k] : Real(1);
                    Mvv[slot] = yv * zv;
                    if (need_grad || need_hess) {
                        const Real zd = has_z ? z_batch.first[q_axis + k] : Real(0);
                        Mdv[slot] = yd * zv;
                        Mvd[slot] = yv * zd;
                    }
                    if (need_hess) {
                        const Real zd = has_z ? z_batch.first[q_axis + k] : Real(0);
                        const Real z2 = has_z ? z_batch.second[q_axis + k] : Real(0);
                        Md2v[slot] = y2 * zv;
                        Mvd2[slot] = yv * z2;
                        Mdd[slot] = yd * zd;
                    }
                }
            }
        }

        if (values_only) {
            if (has_z && num_qpts == 4u) {
                for (std::size_t node = 0; node < tensor_indices.size(); ++node) {
                    const auto& idx = tensor_indices[node];
                    const std::size_t i = idx[0];
                    const std::size_t jk = idx[1] * nz + idx[2];
                    Real* value_row = values_out + node * output_stride;

                    write_tensor_product_value_strided_q<0>(
                        axis_stride, nyz, i, jk, x_batch, Mvv, value_row);
                    write_tensor_product_value_strided_q<1>(
                        axis_stride, nyz, i, jk, x_batch, Mvv, value_row);
                    write_tensor_product_value_strided_q<2>(
                        axis_stride, nyz, i, jk, x_batch, Mvv, value_row);
                    write_tensor_product_value_strided_q<3>(
                        axis_stride, nyz, i, jk, x_batch, Mvv, value_row);
                }
                return;
            }
            for (std::size_t node = 0; node < tensor_indices.size(); ++node) {
                const auto& idx = tensor_indices[node];
                const std::size_t i = idx[0];
                const std::size_t jk = idx[1] * nz + idx[2];
                Real* value_row = values_out + node * output_stride;

                for (std::size_t q = 0; q < num_qpts; ++q) {
                    const std::size_t q_axis = q * axis_stride;
                    const std::size_t slot = q * nyz + jk;
                    value_row[q] = x_batch.values[q_axis + i] * Mvv[slot];
                }
            }
            return;
        }

        if (gradients_only) {
            if (has_z && num_qpts == 4u) {
                if (output_stride == 4u) {
                    for (std::size_t node = 0; node < tensor_indices.size(); ++node) {
                        const auto& idx = tensor_indices[node];
                        const std::size_t i = idx[0];
                        const std::size_t jk = idx[1] * nz + idx[2];
                        Real* grad_row = gradients_out + node * 3u * output_stride;

                        write_tensor_product_gradient_stride4_q<0>(
                            axis_stride, nyz, i, jk, x_batch, Mvv, Mdv, Mvd, grad_row);
                        write_tensor_product_gradient_stride4_q<1>(
                            axis_stride, nyz, i, jk, x_batch, Mvv, Mdv, Mvd, grad_row);
                        write_tensor_product_gradient_stride4_q<2>(
                            axis_stride, nyz, i, jk, x_batch, Mvv, Mdv, Mvd, grad_row);
                        write_tensor_product_gradient_stride4_q<3>(
                            axis_stride, nyz, i, jk, x_batch, Mvv, Mdv, Mvd, grad_row);
                    }
                } else {
                    for (std::size_t node = 0; node < tensor_indices.size(); ++node) {
                        const auto& idx = tensor_indices[node];
                        const std::size_t i = idx[0];
                        const std::size_t jk = idx[1] * nz + idx[2];
                        Real* grad_row = gradients_out + node * 3u * output_stride;

                        write_tensor_product_gradient_strided_q<0>(
                            axis_stride, nyz, i, jk, output_stride, x_batch,
                            Mvv, Mdv, Mvd, grad_row);
                        write_tensor_product_gradient_strided_q<1>(
                            axis_stride, nyz, i, jk, output_stride, x_batch,
                            Mvv, Mdv, Mvd, grad_row);
                        write_tensor_product_gradient_strided_q<2>(
                            axis_stride, nyz, i, jk, output_stride, x_batch,
                            Mvv, Mdv, Mvd, grad_row);
                        write_tensor_product_gradient_strided_q<3>(
                            axis_stride, nyz, i, jk, output_stride, x_batch,
                            Mvv, Mdv, Mvd, grad_row);
                    }
                }
                return;
            }

            for (std::size_t node = 0; node < tensor_indices.size(); ++node) {
                const auto& idx = tensor_indices[node];
                const std::size_t i = idx[0];
                const std::size_t jk = idx[1] * nz + idx[2];
                Real* grad_row = gradients_out + node * 3u * output_stride;

                for (std::size_t q = 0; q < num_qpts; ++q) {
                    const std::size_t q_axis = q * axis_stride;
                    const std::size_t slot = q * nyz + jk;
                    const Real xv = x_batch.values[q_axis + i];
                    const Real xd = x_batch.first[q_axis + i];
                    grad_row[0u * output_stride + q] = xd * Mvv[slot];
                    grad_row[1u * output_stride + q] = xv * Mdv[slot];
                    grad_row[2u * output_stride + q] = xv * Mvd[slot];
                }
            }
            return;
        }

        if (has_z && num_qpts == 4u && hessians_only) {
            if (output_stride == 4u) {
                for (std::size_t node = 0; node < tensor_indices.size(); ++node) {
                    const auto& idx = tensor_indices[node];
                    const std::size_t i = idx[0];
                    const std::size_t jk = idx[1] * nz + idx[2];
                    Real* hess_row = hessians_out + node * 9u * output_stride;

                    write_tensor_product_hessian_stride4_q<0>(
                        axis_stride, nyz, i, jk, x_batch,
                        Mvv, Mdv, Mvd, Md2v, Mvd2, Mdd, hess_row);
                    write_tensor_product_hessian_stride4_q<1>(
                        axis_stride, nyz, i, jk, x_batch,
                        Mvv, Mdv, Mvd, Md2v, Mvd2, Mdd, hess_row);
                    write_tensor_product_hessian_stride4_q<2>(
                        axis_stride, nyz, i, jk, x_batch,
                        Mvv, Mdv, Mvd, Md2v, Mvd2, Mdd, hess_row);
                    write_tensor_product_hessian_stride4_q<3>(
                        axis_stride, nyz, i, jk, x_batch,
                        Mvv, Mdv, Mvd, Md2v, Mvd2, Mdd, hess_row);
                }
            } else {
                for (std::size_t node = 0; node < tensor_indices.size(); ++node) {
                    const auto& idx = tensor_indices[node];
                    const std::size_t i = idx[0];
                    const std::size_t jk = idx[1] * nz + idx[2];
                    Real* hess_row = hessians_out + node * 9u * output_stride;

                    write_tensor_product_hessian_strided_q<0>(
                        axis_stride, nyz, i, jk, output_stride, x_batch,
                        Mvv, Mdv, Mvd, Md2v, Mvd2, Mdd, hess_row);
                    write_tensor_product_hessian_strided_q<1>(
                        axis_stride, nyz, i, jk, output_stride, x_batch,
                        Mvv, Mdv, Mvd, Md2v, Mvd2, Mdd, hess_row);
                    write_tensor_product_hessian_strided_q<2>(
                        axis_stride, nyz, i, jk, output_stride, x_batch,
                        Mvv, Mdv, Mvd, Md2v, Mvd2, Mdd, hess_row);
                    write_tensor_product_hessian_strided_q<3>(
                        axis_stride, nyz, i, jk, output_stride, x_batch,
                        Mvv, Mdv, Mvd, Md2v, Mvd2, Mdd, hess_row);
                }
            }
            return;
        }

        if (has_z && num_qpts == 4u && all_outputs) {
            if (output_stride == 4u) {
                for (std::size_t node = 0; node < tensor_indices.size(); ++node) {
                    const auto& idx = tensor_indices[node];
                    const std::size_t i = idx[0];
                    const std::size_t jk = idx[1] * nz + idx[2];
                    Real* value_row = values_out + node * output_stride;
                    Real* grad_row = gradients_out + node * 3u * output_stride;
                    Real* hess_row = hessians_out + node * 9u * output_stride;

                    write_tensor_product_all_stride4_q<0>(
                        axis_stride, nyz, i, jk, x_batch,
                        Mvv, Mdv, Mvd, Md2v, Mvd2, Mdd, value_row, grad_row, hess_row);
                    write_tensor_product_all_stride4_q<1>(
                        axis_stride, nyz, i, jk, x_batch,
                        Mvv, Mdv, Mvd, Md2v, Mvd2, Mdd, value_row, grad_row, hess_row);
                    write_tensor_product_all_stride4_q<2>(
                        axis_stride, nyz, i, jk, x_batch,
                        Mvv, Mdv, Mvd, Md2v, Mvd2, Mdd, value_row, grad_row, hess_row);
                    write_tensor_product_all_stride4_q<3>(
                        axis_stride, nyz, i, jk, x_batch,
                        Mvv, Mdv, Mvd, Md2v, Mvd2, Mdd, value_row, grad_row, hess_row);
                }
            } else {
                for (std::size_t node = 0; node < tensor_indices.size(); ++node) {
                    const auto& idx = tensor_indices[node];
                    const std::size_t i = idx[0];
                    const std::size_t jk = idx[1] * nz + idx[2];
                    Real* value_row = values_out + node * output_stride;
                    Real* grad_row = gradients_out + node * 3u * output_stride;
                    Real* hess_row = hessians_out + node * 9u * output_stride;

                    write_tensor_product_all_strided_q<0>(
                        axis_stride, nyz, i, jk, output_stride, x_batch,
                        Mvv, Mdv, Mvd, Md2v, Mvd2, Mdd, value_row, grad_row, hess_row);
                    write_tensor_product_all_strided_q<1>(
                        axis_stride, nyz, i, jk, output_stride, x_batch,
                        Mvv, Mdv, Mvd, Md2v, Mvd2, Mdd, value_row, grad_row, hess_row);
                    write_tensor_product_all_strided_q<2>(
                        axis_stride, nyz, i, jk, output_stride, x_batch,
                        Mvv, Mdv, Mvd, Md2v, Mvd2, Mdd, value_row, grad_row, hess_row);
                    write_tensor_product_all_strided_q<3>(
                        axis_stride, nyz, i, jk, output_stride, x_batch,
                        Mvv, Mdv, Mvd, Md2v, Mvd2, Mdd, value_row, grad_row, hess_row);
                }
            }
            return;
        }

        for (std::size_t node = 0; node < tensor_indices.size(); ++node) {
            const auto& idx = tensor_indices[node];
            const std::size_t i = idx[0];
            const std::size_t jk = idx[1] * nz + idx[2];

            Real* value_row = values_out ? values_out + node * output_stride : nullptr;
            Real* grad_row = gradients_out ? gradients_out + node * 3u * output_stride : nullptr;
            Real* hess_row = hessians_out ? hessians_out + node * 9u * output_stride : nullptr;

            for (std::size_t q = 0; q < num_qpts; ++q) {
                const std::size_t q_axis = q * axis_stride;
                const std::size_t slot = q * nyz + jk;
                const Real xv = x_batch.values[q_axis + i];

                if (value_row != nullptr) {
                    value_row[q] = xv * Mvv[slot];
                }

                if (need_grad) {
                    const Real xd = x_batch.first[q_axis + i];
                    grad_row[0u * output_stride + q] = xd * Mvv[slot];
                    grad_row[1u * output_stride + q] = xv * Mdv[slot];
                    grad_row[2u * output_stride + q] = xv * Mvd[slot];
                }

                if (need_hess) {
                    const Real xd = x_batch.first[q_axis + i];
                    const Real x2 = x_batch.second[q_axis + i];
                    const Real hxy = xd * Mdv[slot];
                    const Real hxz = xd * Mvd[slot];
                    const Real hyz = xv * Mdd[slot];
                    hess_row[0u * output_stride + q] = x2 * Mvv[slot];
                    hess_row[4u * output_stride + q] = xv * Md2v[slot];
                    hess_row[8u * output_stride + q] = xv * Mvd2[slot];
                    hess_row[1u * output_stride + q] = hxy;
                    hess_row[3u * output_stride + q] = hxy;
                    hess_row[2u * output_stride + q] = hxz;
                    hess_row[6u * output_stride + q] = hxz;
                    hess_row[5u * output_stride + q] = hyz;
                    hess_row[7u * output_stride + q] = hyz;
                }
            }
        }
        return;
    }

    for (std::size_t node = 0; node < tensor_indices.size(); ++node) {
        const auto& idx = tensor_indices[node];
        const std::size_t i = idx[0];
        const std::size_t j = idx[1];

        Real* value_row = values_out ? values_out + node * output_stride : nullptr;
        Real* grad_row = gradients_out ? gradients_out + node * 3u * output_stride : nullptr;
        Real* hess_row = hessians_out ? hessians_out + node * 9u * output_stride : nullptr;

        for (std::size_t q = 0; q < num_qpts; ++q) {
            const std::size_t q_axis = q * axis_stride;
            const Real xv = x_batch.values[q_axis + i];
            const Real yv = y_batch.values[q_axis + j];

            if (value_row != nullptr) {
                value_row[q] = xv * yv;
            }

            if (need_grad) {
                const Real xd = x_batch.first[q_axis + i];
                const Real yd = y_batch.first[q_axis + j];
                grad_row[0u * output_stride + q] = xd * yv;
                grad_row[1u * output_stride + q] = xv * yd;
                grad_row[2u * output_stride + q] = Real(0);
            }

            if (need_hess) {
                const Real xd = x_batch.first[q_axis + i];
                const Real yd = y_batch.first[q_axis + j];
                const Real x2 = x_batch.second[q_axis + i];
                const Real y2 = y_batch.second[q_axis + j];
                const Real hxy = xd * yd;

                hess_row[0u * output_stride + q] = x2 * yv;
                hess_row[4u * output_stride + q] = xv * y2;
                hess_row[8u * output_stride + q] = Real(0);
                hess_row[1u * output_stride + q] = hxy;
                hess_row[3u * output_stride + q] = hxy;
                hess_row[2u * output_stride + q] = Real(0);
                hess_row[6u * output_stride + q] = Real(0);
                hess_row[5u * output_stride + q] = Real(0);
                hess_row[7u * output_stride + q] = Real(0);
            }
        }
    }
}

void evaluate_wedge_points_strided(
    const std::vector<std::array<int, 4>>& simplex_exponents,
    const std::vector<std::array<std::size_t, 2>>& wedge_indices,
    const std::vector<std::size_t>& wedge_node_by_tri_z,
    int order,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    const Real* v_coeffs,
    const Real* d_coeffs,
    const Real* d2_coeffs,
    const Real* barycentric_weights,
    int n_axis,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) {
    if (points.empty() || wedge_indices.empty()) {
        return;
    }

    const bool want_values = values_out != nullptr;
    const bool need_grad = gradients_out != nullptr;
    const bool need_hess = hessians_out != nullptr;
    const bool values_only = want_values && !need_grad && !need_hess;
    const bool gradients_only = !want_values && need_grad && !need_hess;
    const bool hessians_only = !want_values && !need_grad && need_hess;
    const bool all_outputs = want_values && need_grad && need_hess;
    const bool use_batched_wedge =
        (values_only && order <= 3) ||
        (gradients_only && order >= 2) ||
        (hessians_only && order >= 3) ||
        (all_outputs && order >= 3);
    if (values_only &&
        order >= 4 &&
        order <= 8 &&
        try_evaluate_wedge_values_product_q4(
            simplex_exponents, wedge_indices, order, points, output_stride, values_out)) {
        return;
    }
    const AxisDeriv level = need_hess
        ? AxisDeriv::ValuesAndFirstAndSecond
        : (need_grad ? AxisDeriv::ValuesAndFirst : AxisDeriv::ValuesOnly);

    LagrangeEvaluateScratch& scratch = evaluate_scratch();
    const std::size_t tri_count = simplex_exponents.size();
    if (use_batched_wedge) {
        const std::size_t num_qpts = points.size();
        const std::size_t tri_stride = num_qpts;
        if (num_qpts == 4u &&
            output_stride == 4u &&
            (gradients_only || hessians_only || all_outputs) &&
            order >= 3 &&
            order <= 8 &&
            wedge_node_by_tri_z.size() == tri_count * static_cast<std::size_t>(n_axis)) {
            const bool use_product_axis_batch =
                gradients_only &&
                n_axis >= 5 &&
                n_axis <= 9;
            if (!use_product_axis_batch ||
                !try_fill_axis_batch_product_q4(
                    scratch.axis_z_batch, points, 2u, n_axis, level)) {
                fill_axis_batch(scratch.axis_z_batch,
                                points,
                                2u,
                                v_coeffs,
                                d_coeffs,
                                d2_coeffs,
                                barycentric_weights,
                                n_axis,
                                level);
            }
            if (need_hess) {
                if (try_evaluate_wedge_fused_stride4_q4<true>(
                        simplex_exponents, wedge_node_by_tri_z, order, points,
                        scratch.axis_z_batch, n_axis, values_out, gradients_out, hessians_out)) {
                    return;
                }
            } else if (try_evaluate_wedge_fused_stride4_q4<false>(
                           simplex_exponents, wedge_node_by_tri_z, order, points,
                           scratch.axis_z_batch, n_axis, values_out, gradients_out, hessians_out)) {
                return;
            }
        }

        const std::size_t tri_values_size = tri_count * tri_stride;
        scratch.wedge_tri_values_batch.resize(tri_values_size);
        if (need_grad || need_hess) {
            scratch.wedge_tri_gradient_batch.resize(tri_count * 2u * tri_stride);
        }
        if (need_hess) {
            scratch.wedge_tri_hessian_batch.resize(tri_count * 3u * tri_stride);
        }

        detail::evaluate_triangle_simplex_basis_wedge_components_strided(
            simplex_exponents,
            order,
            points,
            tri_stride,
            scratch.wedge_tri_values_batch.data(),
            (need_grad || need_hess) ? scratch.wedge_tri_gradient_batch.data() : nullptr,
            need_hess ? scratch.wedge_tri_hessian_batch.data() : nullptr);

        const bool use_product_axis_batch =
            gradients_only &&
            points.size() == 4u &&
            n_axis >= 5 &&
            n_axis <= 9;
        if (!use_product_axis_batch ||
            !try_fill_axis_batch_product_q4(
                scratch.axis_z_batch, points, 2u, n_axis, level)) {
            fill_axis_batch(scratch.axis_z_batch,
                            points,
                            2u,
                            v_coeffs,
                            d_coeffs,
                            d2_coeffs,
                            barycentric_weights,
                            n_axis,
                            level);
        }

        const std::size_t axis_stride = static_cast<std::size_t>(n_axis);
        if (all_outputs) {
            if (num_qpts == 4u) {
                if (output_stride == 4u) {
                    for (std::size_t node = 0; node < wedge_indices.size(); ++node) {
                        const auto& index = wedge_indices[node];
                        const std::size_t tri = index[0];
                        const std::size_t z = index[1];
                        Real* value_row = values_out + node * output_stride;
                        Real* g = gradients_out + node * 3u * output_stride;
                        Real* H = hessians_out + node * 9u * output_stride;
                        const Real* tri_g = scratch.wedge_tri_gradient_batch.data() + tri * 2u * tri_stride;
                        const Real* tri_H = scratch.wedge_tri_hessian_batch.data() + tri * 3u * tri_stride;
                        const Real* tri_values = scratch.wedge_tri_values_batch.data();

                        write_wedge_all_stride4_q<0>(
                            tri_stride, axis_stride, tri, z,
                            tri_values, tri_g, tri_H, scratch.axis_z_batch, value_row, g, H);
                        write_wedge_all_stride4_q<1>(
                            tri_stride, axis_stride, tri, z,
                            tri_values, tri_g, tri_H, scratch.axis_z_batch, value_row, g, H);
                        write_wedge_all_stride4_q<2>(
                            tri_stride, axis_stride, tri, z,
                            tri_values, tri_g, tri_H, scratch.axis_z_batch, value_row, g, H);
                        write_wedge_all_stride4_q<3>(
                            tri_stride, axis_stride, tri, z,
                            tri_values, tri_g, tri_H, scratch.axis_z_batch, value_row, g, H);
                    }
                } else {
                    for (std::size_t node = 0; node < wedge_indices.size(); ++node) {
                        const auto& index = wedge_indices[node];
                        const std::size_t tri = index[0];
                        const std::size_t z = index[1];
                        Real* value_row = values_out + node * output_stride;
                        Real* g = gradients_out + node * 3u * output_stride;
                        Real* H = hessians_out + node * 9u * output_stride;
                        const Real* tri_g = scratch.wedge_tri_gradient_batch.data() + tri * 2u * tri_stride;
                        const Real* tri_H = scratch.wedge_tri_hessian_batch.data() + tri * 3u * tri_stride;
                        const Real* tri_values = scratch.wedge_tri_values_batch.data();

                        write_wedge_all_strided_q<0>(
                            tri_stride, axis_stride, tri, z, output_stride,
                            tri_values, tri_g, tri_H, scratch.axis_z_batch, value_row, g, H);
                        write_wedge_all_strided_q<1>(
                            tri_stride, axis_stride, tri, z, output_stride,
                            tri_values, tri_g, tri_H, scratch.axis_z_batch, value_row, g, H);
                        write_wedge_all_strided_q<2>(
                            tri_stride, axis_stride, tri, z, output_stride,
                            tri_values, tri_g, tri_H, scratch.axis_z_batch, value_row, g, H);
                        write_wedge_all_strided_q<3>(
                            tri_stride, axis_stride, tri, z, output_stride,
                            tri_values, tri_g, tri_H, scratch.axis_z_batch, value_row, g, H);
                    }
                }
                return;
            }

            for (std::size_t node = 0; node < wedge_indices.size(); ++node) {
                const auto& index = wedge_indices[node];
                const std::size_t tri = index[0];
                const std::size_t z = index[1];
                Real* value_row = values_out + node * output_stride;
                Real* g = gradients_out + node * 3u * output_stride;
                Real* H = hessians_out + node * 9u * output_stride;
                const Real* tri_g = scratch.wedge_tri_gradient_batch.data() + tri * 2u * tri_stride;
                const Real* tri_H = scratch.wedge_tri_hessian_batch.data() + tri * 3u * tri_stride;

                for (std::size_t q = 0; q < num_qpts; ++q) {
                    const std::size_t tri_q = tri * tri_stride + q;
                    const std::size_t z_q = q * axis_stride + z;
                    const Real tri_v = scratch.wedge_tri_values_batch[tri_q];
                    const Real zv = scratch.axis_z_batch.values[z_q];
                    const Real zd = scratch.axis_z_batch.first[z_q];
                    const Real tri_gx = tri_g[0u * tri_stride + q];
                    const Real tri_gy = tri_g[1u * tri_stride + q];
                    const Real tri_hxx = tri_H[0u * tri_stride + q];
                    const Real tri_hxy = tri_H[1u * tri_stride + q];
                    const Real tri_hyy = tri_H[2u * tri_stride + q];
                    const Real hxz = tri_gx * zd;
                    const Real hxy = tri_hxy * zv;
                    const Real hyz = tri_gy * zd;

                    value_row[q] = tri_v * zv;
                    g[0u * output_stride + q] = tri_gx * zv;
                    g[1u * output_stride + q] = tri_gy * zv;
                    g[2u * output_stride + q] = tri_v * zd;
                    H[0u * output_stride + q] = tri_hxx * zv;
                    H[1u * output_stride + q] = hxy;
                    H[2u * output_stride + q] = hxz;
                    H[3u * output_stride + q] = hxy;
                    H[4u * output_stride + q] = tri_hyy * zv;
                    H[5u * output_stride + q] = hyz;
                    H[6u * output_stride + q] = hxz;
                    H[7u * output_stride + q] = hyz;
                    H[8u * output_stride + q] = tri_v * scratch.axis_z_batch.second[z_q];
                }
            }
            return;
        }

        if (hessians_only) {
            if (num_qpts == 4u) {
                if (output_stride == 4u) {
                    for (std::size_t node = 0; node < wedge_indices.size(); ++node) {
                        const auto& index = wedge_indices[node];
                        const std::size_t tri = index[0];
                        const std::size_t z = index[1];
                        Real* H = hessians_out + node * 9u * output_stride;
                        const Real* tri_g = scratch.wedge_tri_gradient_batch.data() + tri * 2u * tri_stride;
                        const Real* tri_H = scratch.wedge_tri_hessian_batch.data() + tri * 3u * tri_stride;
                        const Real* tri_values = scratch.wedge_tri_values_batch.data();

                        write_wedge_hessian_stride4_q<0>(
                            tri_stride, axis_stride, tri, z,
                            tri_values, tri_g, tri_H, scratch.axis_z_batch, H);
                        write_wedge_hessian_stride4_q<1>(
                            tri_stride, axis_stride, tri, z,
                            tri_values, tri_g, tri_H, scratch.axis_z_batch, H);
                        write_wedge_hessian_stride4_q<2>(
                            tri_stride, axis_stride, tri, z,
                            tri_values, tri_g, tri_H, scratch.axis_z_batch, H);
                        write_wedge_hessian_stride4_q<3>(
                            tri_stride, axis_stride, tri, z,
                            tri_values, tri_g, tri_H, scratch.axis_z_batch, H);
                    }
                } else {
                    for (std::size_t node = 0; node < wedge_indices.size(); ++node) {
                        const auto& index = wedge_indices[node];
                        const std::size_t tri = index[0];
                        const std::size_t z = index[1];
                        Real* H = hessians_out + node * 9u * output_stride;
                        const Real* tri_g = scratch.wedge_tri_gradient_batch.data() + tri * 2u * tri_stride;
                        const Real* tri_H = scratch.wedge_tri_hessian_batch.data() + tri * 3u * tri_stride;
                        const Real* tri_values = scratch.wedge_tri_values_batch.data();

                        write_wedge_hessian_strided_q<0>(
                            tri_stride, axis_stride, tri, z, output_stride,
                            tri_values, tri_g, tri_H, scratch.axis_z_batch, H);
                        write_wedge_hessian_strided_q<1>(
                            tri_stride, axis_stride, tri, z, output_stride,
                            tri_values, tri_g, tri_H, scratch.axis_z_batch, H);
                        write_wedge_hessian_strided_q<2>(
                            tri_stride, axis_stride, tri, z, output_stride,
                            tri_values, tri_g, tri_H, scratch.axis_z_batch, H);
                        write_wedge_hessian_strided_q<3>(
                            tri_stride, axis_stride, tri, z, output_stride,
                            tri_values, tri_g, tri_H, scratch.axis_z_batch, H);
                    }
                }
                return;
            }

            for (std::size_t node = 0; node < wedge_indices.size(); ++node) {
                const auto& index = wedge_indices[node];
                const std::size_t tri = index[0];
                const std::size_t z = index[1];
                Real* H = hessians_out + node * 9u * output_stride;
                const Real* tri_g = scratch.wedge_tri_gradient_batch.data() + tri * 2u * tri_stride;
                const Real* tri_H = scratch.wedge_tri_hessian_batch.data() + tri * 3u * tri_stride;

                for (std::size_t q = 0; q < num_qpts; ++q) {
                    const std::size_t tri_q = tri * tri_stride + q;
                    const std::size_t z_q = q * axis_stride + z;
                    const Real tri_v = scratch.wedge_tri_values_batch[tri_q];
                    const Real zv = scratch.axis_z_batch.values[z_q];
                    const Real zd = scratch.axis_z_batch.first[z_q];
                    const Real tri_gx = tri_g[0u * tri_stride + q];
                    const Real tri_gy = tri_g[1u * tri_stride + q];
                    const Real tri_hxx = tri_H[0u * tri_stride + q];
                    const Real tri_hxy = tri_H[1u * tri_stride + q];
                    const Real tri_hyy = tri_H[2u * tri_stride + q];
                    const Real hxz = tri_gx * zd;
                    const Real hxy = tri_hxy * zv;
                    const Real hyz = tri_gy * zd;

                    H[0u * output_stride + q] = tri_hxx * zv;
                    H[1u * output_stride + q] = hxy;
                    H[2u * output_stride + q] = hxz;
                    H[3u * output_stride + q] = hxy;
                    H[4u * output_stride + q] = tri_hyy * zv;
                    H[5u * output_stride + q] = hyz;
                    H[6u * output_stride + q] = hxz;
                    H[7u * output_stride + q] = hyz;
                    H[8u * output_stride + q] = tri_v * scratch.axis_z_batch.second[z_q];
                }
            }
            return;
        }

        if (gradients_only) {
            if (num_qpts == 4u) {
                if (output_stride == 4u) {
                    for (std::size_t node = 0; node < wedge_indices.size(); ++node) {
                        const auto& index = wedge_indices[node];
                        const std::size_t tri = index[0];
                        const std::size_t z = index[1];
                        Real* g = gradients_out + node * 3u * output_stride;
                        const Real* tri_g = scratch.wedge_tri_gradient_batch.data() + tri * 2u * tri_stride;
                        const Real* tri_values = scratch.wedge_tri_values_batch.data();

                        write_wedge_gradient_stride4_q<0>(
                            tri_stride, axis_stride, tri, z,
                            tri_values, tri_g, scratch.axis_z_batch, g);
                        write_wedge_gradient_stride4_q<1>(
                            tri_stride, axis_stride, tri, z,
                            tri_values, tri_g, scratch.axis_z_batch, g);
                        write_wedge_gradient_stride4_q<2>(
                            tri_stride, axis_stride, tri, z,
                            tri_values, tri_g, scratch.axis_z_batch, g);
                        write_wedge_gradient_stride4_q<3>(
                            tri_stride, axis_stride, tri, z,
                            tri_values, tri_g, scratch.axis_z_batch, g);
                    }
                } else {
                    for (std::size_t node = 0; node < wedge_indices.size(); ++node) {
                        const auto& index = wedge_indices[node];
                        const std::size_t tri = index[0];
                        const std::size_t z = index[1];
                        Real* g = gradients_out + node * 3u * output_stride;
                        const Real* tri_g = scratch.wedge_tri_gradient_batch.data() + tri * 2u * tri_stride;
                        const Real* tri_values = scratch.wedge_tri_values_batch.data();

                        write_wedge_gradient_strided_q<0>(
                            tri_stride, axis_stride, tri, z, output_stride,
                            tri_values, tri_g, scratch.axis_z_batch, g);
                        write_wedge_gradient_strided_q<1>(
                            tri_stride, axis_stride, tri, z, output_stride,
                            tri_values, tri_g, scratch.axis_z_batch, g);
                        write_wedge_gradient_strided_q<2>(
                            tri_stride, axis_stride, tri, z, output_stride,
                            tri_values, tri_g, scratch.axis_z_batch, g);
                        write_wedge_gradient_strided_q<3>(
                            tri_stride, axis_stride, tri, z, output_stride,
                            tri_values, tri_g, scratch.axis_z_batch, g);
                    }
                }
                return;
            }

            for (std::size_t node = 0; node < wedge_indices.size(); ++node) {
                const auto& index = wedge_indices[node];
                const std::size_t tri = index[0];
                const std::size_t z = index[1];
                Real* g = gradients_out + node * 3u * output_stride;
                const Real* tri_g = scratch.wedge_tri_gradient_batch.data() + tri * 2u * tri_stride;

                for (std::size_t q = 0; q < num_qpts; ++q) {
                    const std::size_t tri_q = tri * tri_stride + q;
                    const std::size_t z_q = q * axis_stride + z;
                    const Real tri_v = scratch.wedge_tri_values_batch[tri_q];
                    const Real zv = scratch.axis_z_batch.values[z_q];
                    g[0u * output_stride + q] = tri_g[0u * tri_stride + q] * zv;
                    g[1u * output_stride + q] = tri_g[1u * tri_stride + q] * zv;
                    g[2u * output_stride + q] = tri_v * scratch.axis_z_batch.first[z_q];
                }
            }
            return;
        }

        for (std::size_t node = 0; node < wedge_indices.size(); ++node) {
            const auto& index = wedge_indices[node];
            const std::size_t tri = index[0];
            const std::size_t z = index[1];
            Real* value_row = values_out ? values_out + node * output_stride : nullptr;
            Real* g = gradients_out ? gradients_out + node * 3u * output_stride : nullptr;
            Real* H = hessians_out ? hessians_out + node * 9u * output_stride : nullptr;

            for (std::size_t q = 0; q < num_qpts; ++q) {
                const std::size_t tri_q = tri * tri_stride + q;
                const std::size_t z_q = q * axis_stride + z;
                const Real tri_v = scratch.wedge_tri_values_batch[tri_q];
                const Real zv = scratch.axis_z_batch.values[z_q];
                if (values_out != nullptr) {
                    value_row[q] = tri_v * zv;
                }

                if (need_grad) {
                    const Real* tri_g = scratch.wedge_tri_gradient_batch.data() + tri * 2u * tri_stride;
                    g[0u * output_stride + q] = tri_g[0u * tri_stride + q] * zv;
                    g[1u * output_stride + q] = tri_g[1u * tri_stride + q] * zv;
                    g[2u * output_stride + q] = tri_v * scratch.axis_z_batch.first[z_q];
                }

                if (need_hess) {
                    const Real* tri_g = scratch.wedge_tri_gradient_batch.data() + tri * 2u * tri_stride;
                    const Real* tri_H = scratch.wedge_tri_hessian_batch.data() + tri * 3u * tri_stride;
                    const Real zd = scratch.axis_z_batch.first[z_q];
                    const Real hxz = tri_g[0u * tri_stride + q] * zd;
                    const Real hxy = tri_H[1u * tri_stride + q] * zv;
                    const Real hyz = tri_g[1u * tri_stride + q] * zd;
                    H[0u * output_stride + q] = tri_H[0u * tri_stride + q] * zv;
                    H[1u * output_stride + q] = hxy;
                    H[2u * output_stride + q] = hxz;
                    H[3u * output_stride + q] = hxy;
                    H[4u * output_stride + q] = tri_H[2u * tri_stride + q] * zv;
                    H[5u * output_stride + q] = hyz;
                    H[6u * output_stride + q] = hxz;
                    H[7u * output_stride + q] = hyz;
                    H[8u * output_stride + q] = tri_v * scratch.axis_z_batch.second[z_q];
                }
            }
        }

        return;
    }

    scratch.tri_values.resize(tri_count);
    if (need_grad || need_hess) {
        scratch.tri_gradient_components.resize(tri_count * 3u);
    }
    if (need_hess) {
        scratch.tri_hessian_components.resize(tri_count * 9u);
    }

    for (std::size_t q = 0; q < points.size(); ++q) {
        const auto& xi = points[q];
        const AxisBasisEvaluations z_axis =
            fill_axis_scratch(scratch.axis_z,
                              v_coeffs,
                              d_coeffs,
                              d2_coeffs,
                              barycentric_weights,
                              n_axis,
                              xi[2],
                              level);
        detail::evaluate_triangle_simplex_basis_to(
            simplex_exponents,
            order,
            xi,
            scratch.tri_values.data(),
            (need_grad || need_hess) ? scratch.tri_gradient_components.data() : nullptr,
            need_hess ? scratch.tri_hessian_components.data() : nullptr);

        for (std::size_t node = 0; node < wedge_indices.size(); ++node) {
            const auto& index = wedge_indices[node];
            const std::size_t tri = index[0];
            const std::size_t z = index[1];
            const Real tri_v = scratch.tri_values[tri];
            const Real zv = z_axis.values[z];

            if (values_out != nullptr) {
                values_out[node * output_stride + q] = tri_v * zv;
            }

            if (need_grad) {
                const Real* tri_g = scratch.tri_gradient_components.data() + tri * 3u;
                Real* g = gradients_out + node * 3u * output_stride;
                g[0u * output_stride + q] = tri_g[0] * zv;
                g[1u * output_stride + q] = tri_g[1] * zv;
                g[2u * output_stride + q] = tri_v * z_axis.first[z];
            }

            if (need_hess) {
                const Real* tri_g = scratch.tri_gradient_components.data() + tri * 3u;
                const Real* tri_H = scratch.tri_hessian_components.data() + tri * 9u;
                const Real zd = z_axis.first[z];
                const Real hxz = tri_g[0] * zd;
                const Real hxy = tri_H[1] * zv;
                const Real hyz = tri_g[1] * zd;
                Real* H = hessians_out + node * 9u * output_stride;
                H[0u * output_stride + q] = tri_H[0] * zv;
                H[1u * output_stride + q] = hxy;
                H[2u * output_stride + q] = hxz;
                H[3u * output_stride + q] = hxy;
                H[4u * output_stride + q] = tri_H[4] * zv;
                H[5u * output_stride + q] = hyz;
                H[6u * output_stride + q] = hxz;
                H[7u * output_stride + q] = hyz;
                H[8u * output_stride + q] = tri_v * z_axis.second[z];
            }
        }
    }
}

NormalizedLagrangeRequest normalize_lagrange_request(ElementType element_type, int order) {
    switch (element_type) {
        case ElementType::Line3:
            return {ElementType::Line2, std::max(order, 2)};
        case ElementType::Triangle6:
            return {ElementType::Triangle3, std::max(order, 2)};
        case ElementType::Quad9:
            return {ElementType::Quad4, std::max(order, 2)};
        case ElementType::Quad8:
            throw BasisElementCompatibilityException(
                "Quad8 is a serendipity element; use SerendipityBasis for Quad8",
                __FILE__, __LINE__, __func__);
        case ElementType::Tetra10:
            return {ElementType::Tetra4, std::max(order, 2)};
        case ElementType::Hex27:
            return {ElementType::Hex8, std::max(order, 2)};
        case ElementType::Hex20:
            throw BasisElementCompatibilityException(
                "Hex20 is a serendipity element; use SerendipityBasis for Hex20",
                __FILE__, __LINE__, __func__);
        case ElementType::Wedge18:
            return {ElementType::Wedge6, std::max(order, 2)};
        case ElementType::Wedge15:
            throw BasisElementCompatibilityException(
                "Wedge15 is a serendipity element; use SerendipityBasis for Wedge15",
                __FILE__, __LINE__, __func__);
        case ElementType::Pyramid13:
            throw BasisElementCompatibilityException(
                "Pyramid13 is a serendipity variant; use SerendipityBasis (Pyramid13) or the complete-family Lagrange path via LagrangeBasis (Pyramid5, order >= 2)",
                __FILE__, __LINE__, __func__);
        case ElementType::Pyramid14:
            return {ElementType::Pyramid5, std::max(order, 2)};
        default:
            return {element_type, order};
    }
}

} // namespace

void prewarm_lagrange_basis_scratch(int max_order, std::size_t max_qpts) {
    evaluate_scratch().prewarm(max_order, max_qpts);
}

LagrangeBasis::LagrangeBasis(ElementType type, int order)
    : element_type_(type), dimension_(0), order_(order) {
    const NormalizedLagrangeRequest normalized = normalize_lagrange_request(element_type_, order_);
    element_type_ = normalized.element_type;
    order_ = normalized.order;

    if (order_ < 0) {
        throw BasisConfigurationException("LagrangeBasis requires non-negative polynomial order",
                                          __FILE__, __LINE__, __func__);
    }

    dimension_ = lagrange_topology_traits(element_type_).dimension;

    init_nodes();
    init_evaluation_dispatch();
}

void LagrangeBasis::init_nodes() {
    nodes_.clear();
    nodes_1d_.clear();
    tensor_indices_.clear();
    simplex_exponents_.clear();
    wedge_indices_.clear();
    wedge_node_by_tri_z_.clear();
    axis_v_coeffs_.clear();
    axis_d_coeffs_.clear();
    axis_d2_coeffs_.clear();
    axis_barycentric_weights_.clear();
    const auto topology = lagrange_topology_traits(element_type_).topology;
    topology_id_ = static_cast<int>(topology);
    switch (topology) {
        case LagrangeTopology::Point:
            build_point_nodes();
            return;
        case LagrangeTopology::Line:
            build_tensor_product_nodes(1);
            compute_axis_monomial_coefficients();
            return;
        case LagrangeTopology::Quadrilateral:
            build_tensor_product_nodes(2);
            compute_axis_monomial_coefficients();
            return;
        case LagrangeTopology::Hexahedron:
            build_tensor_product_nodes(3);
            compute_axis_monomial_coefficients();
            return;
        case LagrangeTopology::Triangle:
        case LagrangeTopology::Tetrahedron:
            build_simplex_nodes();
            return;
        case LagrangeTopology::Wedge:
            build_wedge_nodes();
            compute_axis_monomial_coefficients();
            return;
        case LagrangeTopology::Pyramid:
            build_pyramid_nodes();
            return;
        case LagrangeTopology::Unknown:
            break;
    }

    throw BasisElementCompatibilityException("Unsupported element type in LagrangeBasis::init_nodes",
                                             __FILE__, __LINE__, __func__);
}

void LagrangeBasis::init_evaluation_dispatch() {
    const auto topology = static_cast<LagrangeTopology>(topology_id_);
    switch (topology) {
        case LagrangeTopology::Point:
            vector_evaluation_dispatch_ = &LagrangeBasis::evaluate_point_vectors;
            return;
        case LagrangeTopology::Line:
        case LagrangeTopology::Quadrilateral:
        case LagrangeTopology::Hexahedron:
            vector_evaluation_dispatch_ = &LagrangeBasis::evaluate_tensor_product_vectors;
            return;
        case LagrangeTopology::Triangle:
            vector_evaluation_dispatch_ = &LagrangeBasis::evaluate_triangle_vectors;
            return;
        case LagrangeTopology::Tetrahedron:
            vector_evaluation_dispatch_ = &LagrangeBasis::evaluate_tetrahedron_vectors;
            return;
        case LagrangeTopology::Wedge:
            vector_evaluation_dispatch_ = &LagrangeBasis::evaluate_wedge_vectors;
            return;
        case LagrangeTopology::Pyramid:
            vector_evaluation_dispatch_ = &LagrangeBasis::evaluate_pyramid_vectors;
            return;
        case LagrangeTopology::Unknown:
            break;
    }
    vector_evaluation_dispatch_ = &LagrangeBasis::evaluate_unsupported_vectors;
}

void LagrangeBasis::compute_axis_monomial_coefficients() {
    const int N = static_cast<int>(nodes_1d_.size());
    if (N == 0) return;

    axis_barycentric_weights_.resize(static_cast<std::size_t>(N));
    fill_equispaced_barycentric_weights(N, axis_barycentric_weights_.data());

    if (assign_precomputed_axis_coefficients(N, axis_v_coeffs_, axis_d_coeffs_, axis_d2_coeffs_)) {
        return;
    }

    axis_v_coeffs_.assign(static_cast<std::size_t>(N) * static_cast<std::size_t>(N), Real(0));
    if (N >= 2) {
        axis_d_coeffs_.assign(static_cast<std::size_t>(N) * static_cast<std::size_t>(N - 1), Real(0));
    }
    if (N >= 3) {
        axis_d2_coeffs_.assign(static_cast<std::size_t>(N) * static_cast<std::size_t>(N - 2), Real(0));
    }

    if (N == 1) {
        axis_v_coeffs_[0] = Real(1);
        return;
    }

    // For each L_i, compute monomial coefficients of P_i(x) = prod_{j != i} (x - x_j),
    // then divide by w_i = prod_{j != i} (x_i - x_j).
    std::vector<Real> coeffs;
    coeffs.reserve(static_cast<std::size_t>(N));
    for (int i = 0; i < N; ++i) {
        coeffs.assign(1, Real(1));  // start with constant polynomial 1
        for (int j = 0; j < N; ++j) {
            if (j == i) continue;
            // Multiply (x - x_j) into coeffs (in-place via temp).
            std::vector<Real> next(coeffs.size() + 1, Real(0));
            for (std::size_t k = 0; k < coeffs.size(); ++k) {
                next[k]     -= nodes_1d_[static_cast<std::size_t>(j)] * coeffs[k];
                next[k + 1] += coeffs[k];
            }
            coeffs.swap(next);
        }
        // Divide by w_i.
        Real denom = Real(1);
        for (int j = 0; j < N; ++j) {
            if (j == i) continue;
            denom *= (nodes_1d_[static_cast<std::size_t>(i)] - nodes_1d_[static_cast<std::size_t>(j)]);
        }
        const Real inv_denom = Real(1) / denom;
        for (int k = 0; k < N; ++k) {
            axis_v_coeffs_[static_cast<std::size_t>(i) * static_cast<std::size_t>(N) + static_cast<std::size_t>(k)]
                = coeffs[static_cast<std::size_t>(k)] * inv_denom;
        }

        // First derivative coefficients: d/dx (sum_k c_ik * x^k) = sum_{k>=1} k*c_ik * x^(k-1).
        if (N >= 2) {
            for (int k = 1; k < N; ++k) {
                axis_d_coeffs_[static_cast<std::size_t>(i) * static_cast<std::size_t>(N - 1)
                              + static_cast<std::size_t>(k - 1)]
                    = static_cast<Real>(k)
                      * axis_v_coeffs_[static_cast<std::size_t>(i) * static_cast<std::size_t>(N)
                                       + static_cast<std::size_t>(k)];
            }
        }

        // Second derivative coefficients: d^2/dx^2 = sum_{k>=2} k*(k-1)*c_ik * x^(k-2).
        if (N >= 3) {
            for (int k = 2; k < N; ++k) {
                axis_d2_coeffs_[static_cast<std::size_t>(i) * static_cast<std::size_t>(N - 2)
                              + static_cast<std::size_t>(k - 2)]
                    = static_cast<Real>(k * (k - 1))
                      * axis_v_coeffs_[static_cast<std::size_t>(i) * static_cast<std::size_t>(N)
                                       + static_cast<std::size_t>(k)];
            }
        }
    }
}

void LagrangeBasis::build_point_nodes() {
    nodes_.push_back(math::Vector<Real, 3>{Real(0), Real(0), Real(0)});
}

void LagrangeBasis::init_equispaced_1d_nodes() {
    nodes_1d_.clear();
    for (int i = 0; i <= std::max(order_, 0); ++i) {
        nodes_1d_.push_back(detail::equispaced_pm_one_coord(i, order_));
    }
}

void LagrangeBasis::build_tensor_product_nodes(int dimensions) {
    init_equispaced_1d_nodes();

    if (dimensions < 1 || dimensions > 3) {
        throw BasisConfigurationException("LagrangeBasis::build_tensor_product_nodes requires dimension 1, 2, or 3",
                                          __FILE__, __LINE__, __func__);
    }

    nodes_ = ReferenceNodeLayout::get_lagrange_node_coords(element_type_, order_);
    tensor_indices_.resize(nodes_.size(), TensorNodeIndex{0u, 0u, 0u});
    for (std::size_t n = 0; n < nodes_.size(); ++n) {
        tensor_indices_[n][0] = lattice_index_pm_one(
            nodes_[n][0], order_,
            "LagrangeBasis: invalid tensor-product x-coordinate in public node ordering");
        if (dimensions >= 2) {
            tensor_indices_[n][1] = lattice_index_pm_one(
                nodes_[n][1], order_,
                "LagrangeBasis: invalid tensor-product y-coordinate in public node ordering");
        }
        if (dimensions == 3) {
            tensor_indices_[n][2] = lattice_index_pm_one(
                nodes_[n][2], order_,
                "LagrangeBasis: invalid tensor-product z-coordinate in public node ordering");
        }
    }
}

void LagrangeBasis::build_simplex_nodes() {
    nodes_ = ReferenceNodeLayout::get_lagrange_node_coords(element_type_, order_);
    const auto topology = static_cast<LagrangeTopology>(topology_id_);
    simplex_exponents_.clear();
    simplex_exponents_.reserve(nodes_.size());
    for (const auto& node : nodes_) {
        switch (topology) {
            case LagrangeTopology::Triangle:
                simplex_exponents_.push_back(triangle_exponents_from_public_node(node, order_));
                break;
            case LagrangeTopology::Tetrahedron:
                simplex_exponents_.push_back(tetrahedron_exponents_from_public_node(node, order_));
                break;
            default:
                throw BasisElementCompatibilityException("LagrangeBasis::build_simplex_nodes requires simplex topology",
                                                         __FILE__, __LINE__, __func__);
        }
    }
}

void LagrangeBasis::build_wedge_nodes() {
    init_equispaced_1d_nodes();
    const auto triangle_nodes = ReferenceNodeLayout::get_lagrange_node_coords(ElementType::Triangle3, order_);
    simplex_exponents_.clear();
    simplex_exponents_.reserve(triangle_nodes.size());
    std::unordered_map<std::array<int, 4>, std::size_t, SimplexExponentHash> triangle_index_by_exponent;
    triangle_index_by_exponent.reserve(triangle_nodes.size());
    for (std::size_t tri = 0; tri < triangle_nodes.size(); ++tri) {
        const auto exponents = triangle_exponents_from_public_node(triangle_nodes[tri], order_);
        simplex_exponents_.push_back(exponents);
        const auto inserted = triangle_index_by_exponent.emplace(exponents, tri);
        if (!inserted.second) {
            throw BasisNodeOrderingException("LagrangeBasis: duplicate wedge triangle descriptor",
                                             __FILE__, __LINE__, __func__);
        }
    }

    nodes_ = ReferenceNodeLayout::get_lagrange_node_coords(element_type_, order_);
    wedge_indices_.clear();
    wedge_indices_.reserve(nodes_.size());
    const std::size_t z_count = static_cast<std::size_t>(order_ + 1);
    const std::size_t missing_node = nodes_.size();
    wedge_node_by_tri_z_.assign(triangle_nodes.size() * z_count, missing_node);
    for (std::size_t node_index = 0; node_index < nodes_.size(); ++node_index) {
        const auto& node = nodes_[node_index];
        const auto exponents = triangle_exponents_from_public_node(node, order_);
        const auto found = triangle_index_by_exponent.find(exponents);
        if (found == triangle_index_by_exponent.end()) {
            throw BasisNodeOrderingException("LagrangeBasis: failed to resolve wedge triangle descriptor in public ordering",
                                             __FILE__, __LINE__, __func__);
        }
        const std::size_t tri = found->second;
        const std::size_t z =
            lattice_index_pm_one(node[2], order_,
                                 "LagrangeBasis: invalid wedge z-coordinate in public node ordering");
        wedge_indices_.push_back(WedgeNodeIndex{tri, z});
        wedge_node_by_tri_z_[tri * z_count + z] = node_index;
    }
    for (std::size_t entry = 0; entry < wedge_node_by_tri_z_.size(); ++entry) {
        if (wedge_node_by_tri_z_[entry] == missing_node) {
            throw BasisNodeOrderingException("LagrangeBasis: incomplete wedge tensor-product node map",
                                             __FILE__, __LINE__, __func__);
        }
    }
}

void LagrangeBasis::build_pyramid_nodes() {
    nodes_ = detail::lagrange_pyramid::nodes(order_);
}

void LagrangeBasis::evaluate_point_vectors(const math::Vector<Real, 3>&,
                                           std::vector<Real>* values,
                                           std::vector<Gradient>* gradients,
                                           std::vector<Hessian>* hessians) const {
    if (values != nullptr) {
        values->resize(1u);
        (*values)[0] = Real(1);
    }
    if (gradients != nullptr) {
        gradients->resize(1u);
        (*gradients)[0] = Gradient{};
    }
    if (hessians != nullptr) {
        hessians->resize(1u);
        (*hessians)[0] = Hessian{};
    }
}

void LagrangeBasis::evaluate_tensor_product_vectors(const math::Vector<Real, 3>& xi,
                                                    std::vector<Real>* values,
                                                    std::vector<Gradient>* gradients,
                                                    std::vector<Hessian>* hessians) const {
    const auto topology = static_cast<LagrangeTopology>(topology_id_);
    if (evaluate_fixed_lagrange_fast(topology, order_, xi, values, gradients, hessians)) {
        return;
    }

    const int n_axis = static_cast<int>(nodes_1d_.size());
    const Real* vc = axis_v_coeffs_.data();
    const Real* dc = axis_d_coeffs_.data();
    const Real* d2c = axis_d2_coeffs_.data();
    const Real* bw = axis_barycentric_weights_.data();
    const AxisDeriv level = hessians != nullptr ? AxisDeriv::ValuesAndFirstAndSecond
                           : gradients != nullptr ? AxisDeriv::ValuesAndFirst
                                                  : AxisDeriv::ValuesOnly;

    LagrangeEvaluateScratch& scratch = evaluate_scratch();
    const AxisBasisEvaluations x_axis =
        fill_axis_scratch(scratch.axis_x, vc, dc, d2c, bw, n_axis, xi[0], level);
    AxisBasisEvaluations y_axis = constant_axis_basis();
    AxisBasisEvaluations z_axis = constant_axis_basis();

    if (topology != LagrangeTopology::Line) {
        y_axis = fill_axis_scratch(scratch.axis_y, vc, dc, d2c, bw, n_axis, xi[1], level);
    }
    if (topology == LagrangeTopology::Hexahedron) {
        z_axis = fill_axis_scratch(scratch.axis_z, vc, dc, d2c, bw, n_axis, xi[2], level);
    }

    evaluate_tensor_product_factorized(tensor_indices_, x_axis, y_axis, z_axis,
                                       values, gradients, hessians);
}

void LagrangeBasis::evaluate_triangle_vectors(const math::Vector<Real, 3>& xi,
                                              std::vector<Real>* values,
                                              std::vector<Gradient>* gradients,
                                              std::vector<Hessian>* hessians) const {
    const auto topology = static_cast<LagrangeTopology>(topology_id_);
    if (evaluate_fixed_lagrange_fast(topology, order_, xi, values, gradients, hessians)) {
        return;
    }
    detail::evaluate_triangle_simplex_basis(simplex_exponents_, order_, xi,
                                            values, gradients, hessians);
}

void LagrangeBasis::evaluate_tetrahedron_vectors(const math::Vector<Real, 3>& xi,
                                                 std::vector<Real>* values,
                                                 std::vector<Gradient>* gradients,
                                                 std::vector<Hessian>* hessians) const {
    const auto topology = static_cast<LagrangeTopology>(topology_id_);
    if (evaluate_fixed_lagrange_fast(topology, order_, xi, values, gradients, hessians)) {
        return;
    }
    detail::evaluate_tetrahedron_simplex_basis(simplex_exponents_, order_, xi,
                                               values, gradients, hessians);
}

void LagrangeBasis::evaluate_wedge_vectors(const math::Vector<Real, 3>& xi,
                                           std::vector<Real>* values,
                                           std::vector<Gradient>* gradients,
                                           std::vector<Hessian>* hessians) const {
    const int n_axis = static_cast<int>(nodes_1d_.size());
    const Real* vc = axis_v_coeffs_.data();
    const Real* dc = axis_d_coeffs_.data();
    const Real* d2c = axis_d2_coeffs_.data();
    const Real* bw = axis_barycentric_weights_.data();
    const AxisDeriv level = hessians != nullptr ? AxisDeriv::ValuesAndFirstAndSecond
                           : gradients != nullptr ? AxisDeriv::ValuesAndFirst
                                                  : AxisDeriv::ValuesOnly;

    LagrangeEvaluateScratch& scratch = evaluate_scratch();
    const AxisBasisEvaluations z_axis =
        fill_axis_scratch(scratch.axis_z, vc, dc, d2c, bw, n_axis, xi[2], level);

    if (hessians != nullptr) {
        detail::evaluate_triangle_simplex_basis(
            simplex_exponents_, order_, xi,
            &scratch.tri_values, &scratch.tri_gradients, &scratch.tri_hessians);
    } else if (gradients != nullptr) {
        detail::evaluate_triangle_simplex_basis(
            simplex_exponents_, order_, xi,
            &scratch.tri_values, &scratch.tri_gradients, nullptr);
    } else {
        detail::evaluate_triangle_simplex_basis(
            simplex_exponents_, order_, xi,
            &scratch.tri_values, nullptr, nullptr);
    }

    const std::size_t n_nodes = wedge_indices_.size();
    if (values != nullptr) {
        values->resize(n_nodes);
    }
    if (gradients != nullptr) {
        gradients->resize(n_nodes);
    }
    if (hessians != nullptr) {
        hessians->resize(n_nodes);
    }

    for (std::size_t n = 0; n < n_nodes; ++n) {
        const auto& index = wedge_indices_[n];
        const std::size_t tri_idx = index[0];
        const std::size_t z_idx = index[1];
        const Real zv = z_axis.values[z_idx];
        const Real tri_v = scratch.tri_values[tri_idx];

        if (values != nullptr) {
            (*values)[n] = tri_v * zv;
        }
        if (gradients != nullptr) {
            const Real zd = z_axis.first[z_idx];
            (*gradients)[n][0] = scratch.tri_gradients[tri_idx][0] * zv;
            (*gradients)[n][1] = scratch.tri_gradients[tri_idx][1] * zv;
            (*gradients)[n][2] = tri_v * zd;
        }
        if (hessians != nullptr) {
            const Real zd = z_axis.first[z_idx];
            const Real zd2 = z_axis.second[z_idx];
            Hessian H{};
            H(0, 0) = scratch.tri_hessians[tri_idx](0, 0) * zv;
            H(1, 1) = scratch.tri_hessians[tri_idx](1, 1) * zv;
            H(0, 1) = scratch.tri_hessians[tri_idx](0, 1) * zv;
            H(1, 0) = H(0, 1);
            H(2, 2) = tri_v * zd2;
            H(0, 2) = scratch.tri_gradients[tri_idx][0] * zd;
            H(2, 0) = H(0, 2);
            H(1, 2) = scratch.tri_gradients[tri_idx][1] * zd;
            H(2, 1) = H(1, 2);
            (*hessians)[n] = H;
        }
    }
}

void LagrangeBasis::evaluate_pyramid_vectors(const math::Vector<Real, 3>& xi,
                                             std::vector<Real>* values,
                                             std::vector<Gradient>* gradients,
                                             std::vector<Hessian>* hessians) const {
    if (values != nullptr && gradients != nullptr && hessians != nullptr) {
        detail::lagrange_pyramid::evaluate_all(order_, xi, *values, *gradients, *hessians);
        return;
    }
    if (values != nullptr) {
        detail::lagrange_pyramid::evaluate_values(order_, xi, *values);
    }
    if (gradients != nullptr) {
        detail::lagrange_pyramid::evaluate_gradients(order_, xi, *gradients);
    }
    if (hessians != nullptr) {
        detail::lagrange_pyramid::evaluate_hessians(order_, xi, *hessians);
    }
}

void LagrangeBasis::evaluate_unsupported_vectors(const math::Vector<Real, 3>&,
                                                 std::vector<Real>*,
                                                 std::vector<Gradient>*,
                                                 std::vector<Hessian>*) const {
    throw BasisEvaluationException("Unsupported element in LagrangeBasis vector evaluation",
                                   __FILE__, __LINE__, __func__);
}

void LagrangeBasis::evaluate_values(const math::Vector<Real, 3>& xi,
                                    std::vector<Real>& values) const {
    (this->*vector_evaluation_dispatch_)(xi, &values, nullptr, nullptr);
}

void LagrangeBasis::evaluate_gradients(const math::Vector<Real, 3>& xi,
                                       std::vector<Gradient>& gradients) const {
    (this->*vector_evaluation_dispatch_)(xi, nullptr, &gradients, nullptr);
}

void LagrangeBasis::evaluate_hessians(const math::Vector<Real, 3>& xi,
                                      std::vector<Hessian>& hessians) const {
    (this->*vector_evaluation_dispatch_)(xi, nullptr, nullptr, &hessians);
}

void LagrangeBasis::evaluate_all(const math::Vector<Real, 3>& xi,
                                 std::vector<Real>& values,
                                 std::vector<Gradient>& gradients,
                                 std::vector<Hessian>& hessians) const {
    (this->*vector_evaluation_dispatch_)(xi, &values, &gradients, &hessians);
}

void LagrangeBasis::evaluate_values_to(const math::Vector<Real, 3>& xi,
                                       Real* SVMP_RESTRICT values_out) const {
    const auto topology = static_cast<LagrangeTopology>(topology_id_);
    if (evaluate_fixed_lagrange_fast_to(topology, order_, xi, values_out, nullptr, nullptr)) {
        return;
    }

    const int n_axis = static_cast<int>(nodes_1d_.size());
    const Real* vc = axis_v_coeffs_.data();
    const Real* dc = axis_d_coeffs_.data();
    const Real* d2c = axis_d2_coeffs_.data();
    const Real* bw = axis_barycentric_weights_.data();
    switch (topology) {
        case LagrangeTopology::Point:
            values_out[0] = Real(1);
            return;
        case LagrangeTopology::Line:
        case LagrangeTopology::Quadrilateral:
        case LagrangeTopology::Hexahedron: {
            LagrangeEvaluateScratch& scratch = evaluate_scratch();
            const AxisBasisEvaluations x_axis =
                fill_axis_scratch(scratch.axis_x, vc, dc, d2c, bw, n_axis, xi[0], AxisDeriv::ValuesOnly);
            AxisBasisEvaluations y_axis = constant_axis_basis();
            AxisBasisEvaluations z_axis = constant_axis_basis();
            if (topology != LagrangeTopology::Line) {
                y_axis = fill_axis_scratch(scratch.axis_y, vc, dc, d2c, bw, n_axis, xi[1], AxisDeriv::ValuesOnly);
            }
            if (topology == LagrangeTopology::Hexahedron) {
                z_axis = fill_axis_scratch(scratch.axis_z, vc, dc, d2c, bw, n_axis, xi[2], AxisDeriv::ValuesOnly);
            }
            evaluate_tensor_product_factorized_to(tensor_indices_, x_axis, y_axis, z_axis,
                                                  values_out, nullptr, nullptr);
            return;
        }
        case LagrangeTopology::Triangle:
            detail::evaluate_triangle_simplex_basis_to(simplex_exponents_, order_, xi,
                                                       values_out, nullptr, nullptr);
            return;
        case LagrangeTopology::Tetrahedron:
            detail::evaluate_tetrahedron_simplex_basis_to(simplex_exponents_, order_, xi,
                                                          values_out, nullptr, nullptr);
            return;
        case LagrangeTopology::Wedge: {
            LagrangeEvaluateScratch& scratch = evaluate_scratch();
            const AxisBasisEvaluations z_axis =
                fill_axis_scratch(scratch.axis_z, vc, dc, d2c, bw, n_axis, xi[2], AxisDeriv::ValuesOnly);
            scratch.tri_values.resize(simplex_exponents_.size());
            detail::evaluate_triangle_simplex_basis_to(simplex_exponents_, order_, xi,
                                                       scratch.tri_values.data(), nullptr, nullptr);
            for (std::size_t n = 0; n < wedge_indices_.size(); ++n) {
                const auto& index = wedge_indices_[n];
                values_out[n] = scratch.tri_values[index[0]] * z_axis.values[index[1]];
            }
            return;
        }
        case LagrangeTopology::Pyramid: {
            detail::lagrange_pyramid::evaluate_values_to(order_, xi, values_out);
            return;
        }
        case LagrangeTopology::Unknown:
            break;
    }

    throw BasisEvaluationException("Unsupported element in evaluate_values_to",
                                   __FILE__, __LINE__, __func__);
}

void LagrangeBasis::evaluate_gradients_to(const math::Vector<Real, 3>& xi,
                                          Real* SVMP_RESTRICT gradients_out) const {
    const auto topology = static_cast<LagrangeTopology>(topology_id_);
    if (evaluate_fixed_lagrange_fast_to(topology, order_, xi, nullptr, gradients_out, nullptr)) {
        return;
    }

    const int n_axis = static_cast<int>(nodes_1d_.size());
    const Real* vc = axis_v_coeffs_.data();
    const Real* dc = axis_d_coeffs_.data();
    const Real* d2c = axis_d2_coeffs_.data();
    const Real* bw = axis_barycentric_weights_.data();
    switch (topology) {
        case LagrangeTopology::Point:
            gradients_out[0] = Real(0);
            gradients_out[1] = Real(0);
            gradients_out[2] = Real(0);
            return;
        case LagrangeTopology::Line:
        case LagrangeTopology::Quadrilateral:
        case LagrangeTopology::Hexahedron: {
            LagrangeEvaluateScratch& scratch = evaluate_scratch();
            const AxisBasisEvaluations x_axis =
                fill_axis_scratch(scratch.axis_x, vc, dc, d2c, bw, n_axis, xi[0], AxisDeriv::ValuesAndFirst);
            AxisBasisEvaluations y_axis = constant_axis_basis();
            AxisBasisEvaluations z_axis = constant_axis_basis();
            if (topology != LagrangeTopology::Line) {
                y_axis = fill_axis_scratch(scratch.axis_y, vc, dc, d2c, bw, n_axis, xi[1], AxisDeriv::ValuesAndFirst);
            }
            if (topology == LagrangeTopology::Hexahedron) {
                z_axis = fill_axis_scratch(scratch.axis_z, vc, dc, d2c, bw, n_axis, xi[2], AxisDeriv::ValuesAndFirst);
            }
            evaluate_tensor_product_factorized_to(tensor_indices_, x_axis, y_axis, z_axis,
                                                  nullptr, gradients_out, nullptr);
            return;
        }
        case LagrangeTopology::Triangle:
            detail::evaluate_triangle_simplex_basis_to(simplex_exponents_, order_, xi,
                                                       nullptr, gradients_out, nullptr);
            return;
        case LagrangeTopology::Tetrahedron:
            detail::evaluate_tetrahedron_simplex_basis_to(simplex_exponents_, order_, xi,
                                                          nullptr, gradients_out, nullptr);
            return;
        case LagrangeTopology::Wedge: {
            LagrangeEvaluateScratch& scratch = evaluate_scratch();
            const AxisBasisEvaluations z_axis =
                fill_axis_scratch(scratch.axis_z, vc, dc, d2c, bw, n_axis, xi[2], AxisDeriv::ValuesAndFirst);
            const std::size_t tri_count = simplex_exponents_.size();
            scratch.tri_values.resize(tri_count);
            scratch.tri_gradient_components.resize(tri_count * 3u);
            detail::evaluate_triangle_simplex_basis_to(simplex_exponents_, order_, xi,
                                                       scratch.tri_values.data(),
                                                       scratch.tri_gradient_components.data(),
                                                       nullptr);
            for (std::size_t n = 0; n < wedge_indices_.size(); ++n) {
                const auto& index = wedge_indices_[n];
                const std::size_t tri = index[0];
                const std::size_t z = index[1];
                const Real* tri_g = scratch.tri_gradient_components.data() + tri * 3u;
                Real* g = gradients_out + n * 3u;
                g[0] = tri_g[0] * z_axis.values[z];
                g[1] = tri_g[1] * z_axis.values[z];
                g[2] = scratch.tri_values[tri] * z_axis.first[z];
            }
            return;
        }
        case LagrangeTopology::Pyramid: {
            detail::lagrange_pyramid::evaluate_gradients_to(order_, xi, gradients_out);
            return;
        }
        case LagrangeTopology::Unknown:
            break;
    }

    throw BasisEvaluationException("Unsupported element in evaluate_gradients_to",
                                   __FILE__, __LINE__, __func__);
}

void LagrangeBasis::evaluate_hessians_to(const math::Vector<Real, 3>& xi,
                                         Real* SVMP_RESTRICT hessians_out) const {
    const auto topology = static_cast<LagrangeTopology>(topology_id_);
    if (evaluate_fixed_lagrange_fast_to(topology, order_, xi, nullptr, nullptr, hessians_out)) {
        return;
    }

    const int n_axis = static_cast<int>(nodes_1d_.size());
    const Real* vc = axis_v_coeffs_.data();
    const Real* dc = axis_d_coeffs_.data();
    const Real* d2c = axis_d2_coeffs_.data();
    const Real* bw = axis_barycentric_weights_.data();
    switch (topology) {
        case LagrangeTopology::Point:
            for (std::size_t i = 0; i < 9; ++i) {
                hessians_out[i] = Real(0);
            }
            return;
        case LagrangeTopology::Line:
        case LagrangeTopology::Quadrilateral:
        case LagrangeTopology::Hexahedron: {
            LagrangeEvaluateScratch& scratch = evaluate_scratch();
            const AxisBasisEvaluations x_axis =
                fill_axis_scratch(scratch.axis_x, vc, dc, d2c, bw, n_axis, xi[0], AxisDeriv::ValuesAndFirstAndSecond);
            AxisBasisEvaluations y_axis = constant_axis_basis();
            AxisBasisEvaluations z_axis = constant_axis_basis();
            if (topology != LagrangeTopology::Line) {
                y_axis = fill_axis_scratch(scratch.axis_y, vc, dc, d2c, bw, n_axis, xi[1], AxisDeriv::ValuesAndFirstAndSecond);
            }
            if (topology == LagrangeTopology::Hexahedron) {
                z_axis = fill_axis_scratch(scratch.axis_z, vc, dc, d2c, bw, n_axis, xi[2], AxisDeriv::ValuesAndFirstAndSecond);
            }
            evaluate_tensor_product_factorized_to(tensor_indices_, x_axis, y_axis, z_axis,
                                                  nullptr, nullptr, hessians_out);
            return;
        }
        case LagrangeTopology::Triangle:
            detail::evaluate_triangle_simplex_basis_to(simplex_exponents_, order_, xi,
                                                       nullptr, nullptr, hessians_out);
            return;
        case LagrangeTopology::Tetrahedron:
            detail::evaluate_tetrahedron_simplex_basis_to(simplex_exponents_, order_, xi,
                                                          nullptr, nullptr, hessians_out);
            return;
        case LagrangeTopology::Wedge: {
            LagrangeEvaluateScratch& scratch = evaluate_scratch();
            const AxisBasisEvaluations z_axis =
                fill_axis_scratch(scratch.axis_z, vc, dc, d2c, bw, n_axis, xi[2], AxisDeriv::ValuesAndFirstAndSecond);
            const std::size_t tri_count = simplex_exponents_.size();
            scratch.tri_values.resize(tri_count);
            scratch.tri_gradient_components.resize(tri_count * 3u);
            scratch.tri_hessian_components.resize(tri_count * 9u);
            detail::evaluate_triangle_simplex_basis_to(simplex_exponents_, order_, xi,
                                                       scratch.tri_values.data(),
                                                       scratch.tri_gradient_components.data(),
                                                       scratch.tri_hessian_components.data());
            for (std::size_t n = 0; n < wedge_indices_.size(); ++n) {
                const auto& index = wedge_indices_[n];
                const std::size_t tri = index[0];
                const std::size_t z = index[1];
                const Real zv = z_axis.values[z];
                const Real zd = z_axis.first[z];
                const Real zd2 = z_axis.second[z];
                const Real* tri_g = scratch.tri_gradient_components.data() + tri * 3u;
                const Real* tri_H = scratch.tri_hessian_components.data() + tri * 9u;
                const Real hxy = tri_H[1] * zv;
                const Real hxz = tri_g[0] * zd;
                const Real hyz = tri_g[1] * zd;
                Real* H = hessians_out + n * 9u;
                H[0] = tri_H[0] * zv;
                H[4] = tri_H[4] * zv;
                H[1] = hxy;
                H[3] = hxy;
                H[8] = scratch.tri_values[tri] * zd2;
                H[2] = hxz;
                H[6] = hxz;
                H[5] = hyz;
                H[7] = hyz;
            }
            return;
        }
        case LagrangeTopology::Pyramid: {
            detail::lagrange_pyramid::evaluate_hessians_to(order_, xi, hessians_out);
            return;
        }
        case LagrangeTopology::Unknown:
            break;
    }

    throw BasisEvaluationException("Unsupported element in evaluate_hessians_to",
                                   __FILE__, __LINE__, __func__);
}

void LagrangeBasis::evaluate_all_to(const math::Vector<Real, 3>& xi,
                                    Real* SVMP_RESTRICT values_out,
                                    Real* SVMP_RESTRICT gradients_out,
                                    Real* SVMP_RESTRICT hessians_out) const {
    const auto topology = static_cast<LagrangeTopology>(topology_id_);
    if (evaluate_fixed_lagrange_fast_to(topology, order_, xi, values_out, gradients_out, hessians_out)) {
        return;
    }

    const int n_axis = static_cast<int>(nodes_1d_.size());
    const Real* vc = axis_v_coeffs_.data();
    const Real* dc = axis_d_coeffs_.data();
    const Real* d2c = axis_d2_coeffs_.data();
    const Real* bw = axis_barycentric_weights_.data();
    switch (topology) {
        case LagrangeTopology::Point:
            values_out[0] = Real(1);
            gradients_out[0] = Real(0);
            gradients_out[1] = Real(0);
            gradients_out[2] = Real(0);
            for (std::size_t i = 0; i < 9; ++i) {
                hessians_out[i] = Real(0);
            }
            return;
        case LagrangeTopology::Line:
        case LagrangeTopology::Quadrilateral:
        case LagrangeTopology::Hexahedron: {
            LagrangeEvaluateScratch& scratch = evaluate_scratch();
            const AxisBasisEvaluations x_axis =
                fill_axis_scratch(scratch.axis_x, vc, dc, d2c, bw, n_axis, xi[0], AxisDeriv::ValuesAndFirstAndSecond);
            AxisBasisEvaluations y_axis = constant_axis_basis();
            AxisBasisEvaluations z_axis = constant_axis_basis();
            if (topology != LagrangeTopology::Line) {
                y_axis = fill_axis_scratch(scratch.axis_y, vc, dc, d2c, bw, n_axis, xi[1], AxisDeriv::ValuesAndFirstAndSecond);
            }
            if (topology == LagrangeTopology::Hexahedron) {
                z_axis = fill_axis_scratch(scratch.axis_z, vc, dc, d2c, bw, n_axis, xi[2], AxisDeriv::ValuesAndFirstAndSecond);
            }
            evaluate_tensor_product_factorized_to(tensor_indices_, x_axis, y_axis, z_axis,
                                                  values_out, gradients_out, hessians_out);
            return;
        }
        case LagrangeTopology::Triangle:
            detail::evaluate_triangle_simplex_basis_to(simplex_exponents_, order_, xi,
                                                       values_out, gradients_out, hessians_out);
            return;
        case LagrangeTopology::Tetrahedron:
            detail::evaluate_tetrahedron_simplex_basis_to(simplex_exponents_, order_, xi,
                                                          values_out, gradients_out, hessians_out);
            return;
        case LagrangeTopology::Wedge: {
            LagrangeEvaluateScratch& scratch = evaluate_scratch();
            const AxisBasisEvaluations z_axis =
                fill_axis_scratch(scratch.axis_z, vc, dc, d2c, bw, n_axis, xi[2], AxisDeriv::ValuesAndFirstAndSecond);
            const std::size_t tri_count = simplex_exponents_.size();
            scratch.tri_values.resize(tri_count);
            scratch.tri_gradient_components.resize(tri_count * 3u);
            scratch.tri_hessian_components.resize(tri_count * 9u);
            detail::evaluate_triangle_simplex_basis_to(simplex_exponents_, order_, xi,
                                                       scratch.tri_values.data(),
                                                       scratch.tri_gradient_components.data(),
                                                       scratch.tri_hessian_components.data());
            for (std::size_t n = 0; n < wedge_indices_.size(); ++n) {
                const auto& index = wedge_indices_[n];
                const std::size_t tri = index[0];
                const std::size_t z = index[1];
                const Real zv = z_axis.values[z];
                const Real zd = z_axis.first[z];
                const Real zd2 = z_axis.second[z];
                const Real tri_v = scratch.tri_values[tri];
                const Real* tri_g = scratch.tri_gradient_components.data() + tri * 3u;
                const Real* tri_H = scratch.tri_hessian_components.data() + tri * 9u;
                const Real hxy = tri_H[1] * zv;
                const Real hxz = tri_g[0] * zd;
                const Real hyz = tri_g[1] * zd;

                values_out[n] = tri_v * zv;

                Real* g = gradients_out + n * 3u;
                g[0] = tri_g[0] * zv;
                g[1] = tri_g[1] * zv;
                g[2] = tri_v * zd;

                Real* H = hessians_out + n * 9u;
                H[0] = tri_H[0] * zv;
                H[4] = tri_H[4] * zv;
                H[1] = hxy;
                H[3] = hxy;
                H[8] = tri_v * zd2;
                H[2] = hxz;
                H[6] = hxz;
                H[5] = hyz;
                H[7] = hyz;
            }
            return;
        }
        case LagrangeTopology::Pyramid: {
            detail::lagrange_pyramid::evaluate_all_to(
                order_, xi, values_out, gradients_out, hessians_out);
            return;
        }
        case LagrangeTopology::Unknown:
            break;
    }

    throw BasisEvaluationException("Unsupported element in evaluate_all_to",
                                   __FILE__, __LINE__, __func__);
}

void LagrangeBasis::evaluate_at_quadrature_points(
    const std::vector<math::Vector<Real, 3>>& points,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) const {
    evaluate_at_quadrature_points_strided(points, points.size(), values_out, gradients_out, hessians_out);
}

void LagrangeBasis::evaluate_at_quadrature_points_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) const {
    const std::size_t num_qpts = points.size();
    const std::size_t num_dofs = size();
    if (output_stride < num_qpts) {
        throw BasisConfigurationException("LagrangeBasis strided evaluation requires output_stride >= points.size()",
                                          __FILE__, __LINE__, __func__);
    }
    if (values_out == nullptr && gradients_out == nullptr && hessians_out == nullptr) {
        return;
    }

    const auto topology = static_cast<LagrangeTopology>(topology_id_);
    if (evaluate_fixed_lagrange_fast_strided(topology,
                                             order_,
                                             points,
                                             output_stride,
                                             values_out,
                                             gradients_out,
                                             hessians_out)) {
        return;
    }

    if (topology == LagrangeTopology::Line ||
        topology == LagrangeTopology::Quadrilateral ||
        topology == LagrangeTopology::Hexahedron) {
        evaluate_tensor_product_points_strided(topology,
                                               tensor_indices_,
                                               points,
                                               output_stride,
                                               axis_v_coeffs_.data(),
                                               axis_d_coeffs_.data(),
                                               axis_d2_coeffs_.data(),
                                               axis_barycentric_weights_.data(),
                                               static_cast<int>(nodes_1d_.size()),
                                               values_out,
                                               gradients_out,
                                               hessians_out);
        return;
    }

    if (topology == LagrangeTopology::Triangle) {
        detail::evaluate_triangle_simplex_basis_strided(
            simplex_exponents_, order_, points, output_stride, values_out, gradients_out, hessians_out);
        return;
    }

    if (topology == LagrangeTopology::Tetrahedron) {
        detail::evaluate_tetrahedron_simplex_basis_strided(
            simplex_exponents_, order_, points, output_stride, values_out, gradients_out, hessians_out);
        return;
    }

    if (topology == LagrangeTopology::Wedge &&
        evaluate_wedge_fast_strided(order_,
                                    wedge_indices_,
                                    points,
                                    output_stride,
                                    values_out,
                                    gradients_out,
                                    hessians_out)) {
        return;
    }

    const bool wedge_scalar_hessian_fallback =
        topology == LagrangeTopology::Wedge &&
        values_out == nullptr &&
        gradients_out == nullptr &&
        hessians_out != nullptr &&
        order_ <= 2;
    if (topology == LagrangeTopology::Wedge && !wedge_scalar_hessian_fallback) {
        evaluate_wedge_points_strided(simplex_exponents_,
                                      wedge_indices_,
                                      wedge_node_by_tri_z_,
                                      order_,
                                      points,
                                      output_stride,
                                      axis_v_coeffs_.data(),
                                      axis_d_coeffs_.data(),
                                      axis_d2_coeffs_.data(),
                                      axis_barycentric_weights_.data(),
                                      static_cast<int>(nodes_1d_.size()),
                                      values_out,
                                      gradients_out,
                                      hessians_out);
        return;
    }

    if (topology == LagrangeTopology::Pyramid) {
        detail::lagrange_pyramid::evaluate_at_quadrature_points_strided(
            order_, points, output_stride, values_out, gradients_out, hessians_out);
        return;
    }

    auto& scratch = evaluate_scratch();
    auto& v_tmp = scratch.strided_values_tmp;
    auto& g_tmp = scratch.strided_gradients_tmp;
    auto& h_tmp = scratch.strided_hessians_tmp;

    if (values_out)    v_tmp.resize(num_dofs);
    if (gradients_out) g_tmp.resize(num_dofs * 3u);
    if (hessians_out)  h_tmp.resize(num_dofs * 9u);

    for (std::size_t q = 0; q < num_qpts; ++q) {
        if (values_out && gradients_out && hessians_out) {
            evaluate_all_to(points[q], v_tmp.data(), g_tmp.data(), h_tmp.data());
        } else {
            if (values_out)    evaluate_values_to(points[q], v_tmp.data());
            if (gradients_out) evaluate_gradients_to(points[q], g_tmp.data());
            if (hessians_out)  evaluate_hessians_to(points[q], h_tmp.data());
        }

        if (values_out) {
            for (std::size_t d = 0; d < num_dofs; ++d) {
                values_out[d * output_stride + q] = v_tmp[d];
            }
        }
        if (gradients_out) {
            for (std::size_t d = 0; d < num_dofs; ++d) {
                gradients_out[(d * 3u + 0u) * output_stride + q] = g_tmp[d * 3u + 0u];
                gradients_out[(d * 3u + 1u) * output_stride + q] = g_tmp[d * 3u + 1u];
                gradients_out[(d * 3u + 2u) * output_stride + q] = g_tmp[d * 3u + 2u];
            }
        }
        if (hessians_out) {
            for (std::size_t d = 0; d < num_dofs; ++d) {
                scatter_hessian_components_strided(
                    h_tmp.data() + d * 9u,
                    hessians_out + d * 9u * output_stride,
                    output_stride,
                    q);
            }
        }
    }
}

} // namespace basis
} // namespace FE
} // namespace svmp
