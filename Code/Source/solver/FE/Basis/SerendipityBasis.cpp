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

static const int hex20_monomial_exponents[20][3] = {
    {0, 0, 0}, {0, 0, 1}, {0, 0, 2}, {0, 1, 0}, {0, 1, 1},
    {0, 1, 2}, {0, 2, 0}, {0, 2, 1}, {1, 0, 0}, {1, 0, 1},
    {1, 0, 2}, {1, 1, 0}, {1, 1, 1}, {1, 1, 2}, {1, 2, 0},
    {1, 2, 1}, {2, 0, 0}, {2, 0, 1}, {2, 1, 0}, {2, 1, 1}
};

static const Real hex20_coeffs[20][20] = {
    {-0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25},
    {0.125, 0.125, 0.125, 0.125, -0.125, -0.125, -0.125, -0.125, -0.25, 0.25, -0.25, 0.25, -0.25, -0.25, 0.25, 0.25, 0, 0, 0, 0},
    {0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0, 0, 0, 0, 0, 0, 0, 0, -0.25, -0.25, -0.25, -0.25},
    {0.125, 0.125, -0.125, -0.125, 0.125, 0.125, -0.125, -0.125, -0.25, -0.25, 0.25, 0.25, 0, 0, 0, 0, -0.25, -0.25, 0.25, 0.25},
    {0, 0, 0, 0, 0, 0, 0, 0, 0.25, -0.25, -0.25, 0.25, 0, 0, 0, 0, 0, 0, 0, 0},
    {-0.125, -0.125, 0.125, 0.125, -0.125, -0.125, 0.125, 0.125, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, 0.25, -0.25, -0.25},
    {0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0, 0, 0, 0, -0.25, -0.25, -0.25, -0.25, 0, 0, 0, 0},
    {-0.125, -0.125, -0.125, -0.125, 0.125, 0.125, 0.125, 0.125, 0, 0, 0, 0, 0.25, 0.25, -0.25, -0.25, 0, 0, 0, 0},
    {0.125, -0.125, -0.125, 0.125, 0.125, -0.125, -0.125, 0.125, 0, 0, 0, 0, -0.25, 0.25, -0.25, 0.25, -0.25, 0.25, -0.25, 0.25},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, -0.25, -0.25, 0.25, 0, 0, 0, 0},
    {-0.125, 0.125, 0.125, -0.125, -0.125, 0.125, 0.125, -0.125, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, -0.25, 0.25, -0.25},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, -0.25, -0.25, 0.25},
    {-0.125, 0.125, -0.125, 0.125, 0.125, -0.125, 0.125, -0.125, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0.125, -0.125, 0.125, -0.125, 0.125, -0.125, 0.125, -0.125, 0, 0, 0, 0, 0, 0, 0, 0, -0.25, 0.25, 0.25, -0.25},
    {-0.125, 0.125, 0.125, -0.125, -0.125, 0.125, 0.125, -0.125, 0, 0, 0, 0, 0.25, -0.25, 0.25, -0.25, 0, 0, 0, 0},
    {0.125, -0.125, -0.125, 0.125, -0.125, 0.125, 0.125, -0.125, 0, 0, 0, 0, -0.25, 0.25, 0.25, -0.25, 0, 0, 0, 0},
    {0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, -0.25, -0.25, -0.25, -0.25, 0, 0, 0, 0, 0, 0, 0, 0},
    {-0.125, -0.125, -0.125, -0.125, 0.125, 0.125, 0.125, 0.125, 0.25, -0.25, 0.25, -0.25, 0, 0, 0, 0, 0, 0, 0, 0},
    {-0.125, -0.125, 0.125, 0.125, -0.125, -0.125, 0.125, 0.125, 0.25, 0.25, -0.25, -0.25, 0, 0, 0, 0, 0, 0, 0, 0},
    {0.125, 0.125, -0.125, -0.125, -0.125, -0.125, 0.125, 0.125, -0.25, 0.25, 0.25, -0.25, 0, 0, 0, 0, 0, 0, 0, 0}
};

inline std::array<Real, 3> quadratic_powers(Real x) {
    return {Real(1), x, x * x};
}

void eval_hex20_internal(Real r, Real s, Real t, std::span<Real> internal_vals) {
    const auto rp = quadratic_powers(r);
    const auto sp = quadratic_powers(s);
    const auto tp = quadratic_powers(t);
    Real phi[20];
    for (int j = 0; j < 20; ++j) {
        const int a = hex20_monomial_exponents[j][0];
        const int b = hex20_monomial_exponents[j][1];
        const int c = hex20_monomial_exponents[j][2];
        phi[j] = rp[static_cast<std::size_t>(a)] *
                 sp[static_cast<std::size_t>(b)] *
                 tp[static_cast<std::size_t>(c)];
    }
    for (int i = 0; i < 20; ++i) {
        Real v = Real(0);
        for (int j = 0; j < 20; ++j) {
            v += hex20_coeffs[j][i] * phi[j];
        }
        internal_vals[i] = v;
    }
}

void eval_hex20_grad_internal(Real r, Real s, Real t, std::span<Gradient> internal_grads) {
    const auto rp = quadratic_powers(r);
    const auto sp = quadratic_powers(s);
    const auto tp = quadratic_powers(t);
    Real dphi_dr[20], dphi_ds[20], dphi_dt[20];
    for (int j = 0; j < 20; ++j) {
        const int a = hex20_monomial_exponents[j][0];
        const int b = hex20_monomial_exponents[j][1];
        const int c = hex20_monomial_exponents[j][2];

        dphi_dr[j] = (a > 0) ? Real(a) * rp[static_cast<std::size_t>(a - 1)] *
                                    sp[static_cast<std::size_t>(b)] *
                                    tp[static_cast<std::size_t>(c)]
                              : Real(0);
        dphi_ds[j] = (b > 0) ? rp[static_cast<std::size_t>(a)] *
                                    Real(b) * sp[static_cast<std::size_t>(b - 1)] *
                                    tp[static_cast<std::size_t>(c)]
                              : Real(0);
        dphi_dt[j] = (c > 0) ? rp[static_cast<std::size_t>(a)] *
                                    sp[static_cast<std::size_t>(b)] *
                                    Real(c) * tp[static_cast<std::size_t>(c - 1)]
                              : Real(0);
    }

    for (int i = 0; i < 20; ++i) {
        Real gr = Real(0), gs = Real(0), gt = Real(0);
        for (int j = 0; j < 20; ++j) {
            gr += hex20_coeffs[j][i] * dphi_dr[j];
            gs += hex20_coeffs[j][i] * dphi_ds[j];
            gt += hex20_coeffs[j][i] * dphi_dt[j];
        }
        internal_grads[i][0] = gr;
        internal_grads[i][1] = gs;
        internal_grads[i][2] = gt;
    }
}

void eval_hex20_hess_internal(Real r, Real s, Real t, std::span<Hessian> internal_hessians) {
    const auto rp = quadratic_powers(r);
    const auto sp = quadratic_powers(s);
    const auto tp = quadratic_powers(t);
    Real d2phi_drr[20], d2phi_dss[20], d2phi_dtt[20];
    Real d2phi_drs[20], d2phi_drt[20], d2phi_dst[20];
    for (int j = 0; j < 20; ++j) {
        const int a = hex20_monomial_exponents[j][0];
        const int b = hex20_monomial_exponents[j][1];
        const int c = hex20_monomial_exponents[j][2];

        d2phi_drr[j] = (a > 1) ? Real(a * (a - 1)) *
                                      rp[static_cast<std::size_t>(a - 2)] *
                                      sp[static_cast<std::size_t>(b)] *
                                      tp[static_cast<std::size_t>(c)]
                                : Real(0);
        d2phi_dss[j] = (b > 1) ? rp[static_cast<std::size_t>(a)] *
                                      Real(b * (b - 1)) *
                                      sp[static_cast<std::size_t>(b - 2)] *
                                      tp[static_cast<std::size_t>(c)]
                                : Real(0);
        d2phi_dtt[j] = (c > 1) ? rp[static_cast<std::size_t>(a)] *
                                      sp[static_cast<std::size_t>(b)] *
                                      Real(c * (c - 1)) *
                                      tp[static_cast<std::size_t>(c - 2)]
                                : Real(0);
        d2phi_drs[j] = (a > 0 && b > 0) ? Real(a * b) *
                                              rp[static_cast<std::size_t>(a - 1)] *
                                              sp[static_cast<std::size_t>(b - 1)] *
                                              tp[static_cast<std::size_t>(c)]
                                        : Real(0);
        d2phi_drt[j] = (a > 0 && c > 0) ? Real(a * c) *
                                              rp[static_cast<std::size_t>(a - 1)] *
                                              sp[static_cast<std::size_t>(b)] *
                                              tp[static_cast<std::size_t>(c - 1)]
                                        : Real(0);
        d2phi_dst[j] = (b > 0 && c > 0) ? rp[static_cast<std::size_t>(a)] *
                                              Real(b * c) *
                                              sp[static_cast<std::size_t>(b - 1)] *
                                              tp[static_cast<std::size_t>(c - 1)]
                                        : Real(0);
    }

    for (int i = 0; i < 20; ++i) {
        Hessian H = Hessian::Zero();
        for (int j = 0; j < 20; ++j) {
            H(0, 0) += hex20_coeffs[j][i] * d2phi_drr[j];
            H(1, 1) += hex20_coeffs[j][i] * d2phi_dss[j];
            H(2, 2) += hex20_coeffs[j][i] * d2phi_dtt[j];
            H(0, 1) += hex20_coeffs[j][i] * d2phi_drs[j];
            H(0, 2) += hex20_coeffs[j][i] * d2phi_drt[j];
            H(1, 2) += hex20_coeffs[j][i] * d2phi_dst[j];
        }
        H(1, 0) = H(0, 1);
        H(2, 0) = H(0, 2);
        H(2, 1) = H(1, 2);
        internal_hessians[i] = H;
    }
}

void eval_wedge15_polynomial(Real r,
                             Real s,
                             Real t,
                             std::span<Real> values,
                             std::span<Gradient> gradients,
                             std::span<Hessian> hessians) {
    Real phi[15]{};
    Real dr[15]{};
    Real ds[15]{};
    Real dt[15]{};
    Real drr[15]{};
    Real dss[15]{};
    Real dtt[15]{};
    Real drs[15]{};
    Real drt[15]{};
    Real dst[15]{};

    const auto rp = quadratic_powers(r);
    const auto sp = quadratic_powers(s);
    const auto tp = quadratic_powers(t);

    for (int j = 0; j < 15; ++j) {
        const auto& exponent = kWedge15MonomialExponents[static_cast<std::size_t>(j)];
        const int a = exponent[0];
        const int b = exponent[1];
        const int c = exponent[2];
        const auto ar = static_cast<std::size_t>(a);
        const auto bs = static_cast<std::size_t>(b);
        const auto ct = static_cast<std::size_t>(c);

        const Real ra = rp[ar];
        const Real sb = sp[bs];
        const Real tc = tp[ct];

        if (!values.empty()) {
            phi[j] = ra * sb * tc;
        }
        if (!gradients.empty()) {
            dr[j] = (a > 0) ? Real(a) * rp[ar - 1u] * sb * tc : Real(0);
            ds[j] = (b > 0) ? ra * Real(b) * sp[bs - 1u] * tc : Real(0);
            dt[j] = (c > 0) ? ra * sb * Real(c) * tp[ct - 1u] : Real(0);
        }
        if (!hessians.empty()) {
            drr[j] = (a > 1) ? Real(a * (a - 1)) * rp[ar - 2u] * sb * tc : Real(0);
            dss[j] = (b > 1) ? ra * Real(b * (b - 1)) * sp[bs - 2u] * tc : Real(0);
            dtt[j] = (c > 1) ? ra * sb * Real(c * (c - 1)) * tp[ct - 2u] : Real(0);
            drs[j] = (a > 0 && b > 0) ? Real(a * b) * rp[ar - 1u] * sp[bs - 1u] * tc : Real(0);
            drt[j] = (a > 0 && c > 0) ? Real(a * c) * rp[ar - 1u] * sb * tp[ct - 1u] : Real(0);
            dst[j] = (b > 0 && c > 0) ? ra * Real(b * c) * sp[bs - 1u] * tp[ct - 1u] : Real(0);
        }
    }

    for (int i = 0; i < 15; ++i) {
        Real value = Real(0);
        Real gr = Real(0);
        Real gs = Real(0);
        Real gt = Real(0);
        Hessian H = Hessian::Zero();
        for (int j = 0; j < 15; ++j) {
            const Real coefficient =
                kWedge15Coefficients[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)];
            if (!values.empty()) {
                value += coefficient * phi[j];
            }
            if (!gradients.empty()) {
                gr += coefficient * dr[j];
                gs += coefficient * ds[j];
                gt += coefficient * dt[j];
            }
            if (!hessians.empty()) {
                H(0, 0) += coefficient * drr[j];
                H(1, 1) += coefficient * dss[j];
                H(2, 2) += coefficient * dtt[j];
                H(0, 1) += coefficient * drs[j];
                H(0, 2) += coefficient * drt[j];
                H(1, 2) += coefficient * dst[j];
            }
        }

        const std::size_t index = static_cast<std::size_t>(i);
        if (!values.empty()) {
            values[index] = value;
        }
        if (!gradients.empty()) {
            gradients[index][0] = gr;
            gradients[index][1] = gs;
            gradients[index][2] = gt;
        }
        if (!hessians.empty()) {
            H(1, 0) = H(0, 1);
            H(2, 0) = H(0, 2);
            H(2, 1) = H(1, 2);
            hessians[index] = H;
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

SerendipityBasis::SerendipityBasis(ElementType type, int order, bool geometry_mode)
    : element_type_(type), dimension_(0), order_(order), size_(0), geometry_mode_(geometry_mode) {
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

        for (std::size_t j = 0; j < size_; ++j) {
            const auto [ax, ay] = quad_monomial_exponents_[j];
            const Real value = std::pow(x, ax) * std::pow(y, ay);
            const Real dx =
                (ax > 0) ? Real(ax) * std::pow(x, ax - 1) * std::pow(y, ay) : Real(0);
            const Real dy =
                (ay > 0) ? std::pow(x, ax) * Real(ay) * std::pow(y, ay - 1) : Real(0);
            const Real dxx =
                (ax > 1) ? Real(ax * (ax - 1)) * std::pow(x, ax - 2) * std::pow(y, ay)
                         : Real(0);
            const Real dxy =
                (ax > 0 && ay > 0)
                    ? Real(ax * ay) * std::pow(x, ax - 1) * std::pow(y, ay - 1)
                    : Real(0);
            const Real dyy =
                (ay > 1) ? Real(ay * (ay - 1)) * std::pow(x, ax) * std::pow(y, ay - 2)
                         : Real(0);

            for (std::size_t i = 0; i < size_; ++i) {
                const Real coeff = quad_inv_vandermonde_[j * size_ + i];
                if (!values_out.empty()) {
                    values_out[i] += value * coeff;
                }
                if (!gradients_out.empty()) {
                    Gradient& g = gradients_out[i];
                    g[0] += dx * coeff;
                    g[1] += dy * coeff;
                }
                if (!hessians_out.empty()) {
                    Hessian& h = hessians_out[i];
                    h(0, 0) += dxx * coeff;
                    h(0, 1) += dxy * coeff;
                    h(1, 0) += dxy * coeff;
                    h(1, 1) += dyy * coeff;
                }
            }
        }
        return;
    }

    if (dimension_ == 3 && order_ == 1) {
        evaluate_hex8_reference(x, y, z, values_out, gradients_out, hessians_out);
        return;
    }

    if (geometry_mode_ && element_type_ == ElementType::Hex20) {
        evaluate_hex8_reference(x, y, z, values_out, gradients_out, hessians_out);
        return;
    }

    if (element_type_ == ElementType::Hex20) {
        const auto mesh_to_basis = ReferenceNodeLayout::mesh_to_basis_ordering(element_type_);
        FE::throw_if<BasisEvaluationException>(mesh_to_basis.size() != size_, SVMP_HERE,
                                               "Hex20 mesh-to-basis ordering is not registered");

        if (!values_out.empty()) {
            std::array<Real, 20u> internal_vals{};
            eval_hex20_internal(x, y, z, internal_vals);
            for (std::size_t i = 0; i < 20u; ++i) {
                values_out[i] = internal_vals[mesh_to_basis[i]];
            }
        }
        if (!gradients_out.empty()) {
            std::array<Gradient, 20u> internal_grads{};
            eval_hex20_grad_internal(x, y, z, internal_grads);
            for (std::size_t i = 0; i < 20u; ++i) {
                gradients_out[i] = internal_grads[mesh_to_basis[i]];
            }
        }
        if (!hessians_out.empty()) {
            std::array<Hessian, 20u> internal_hessians{};
            eval_hex20_hess_internal(x, y, z, internal_hessians);
            for (std::size_t i = 0; i < 20u; ++i) {
                hessians_out[i] = internal_hessians[mesh_to_basis[i]];
            }
        }
        return;
    }

    if (element_type_ == ElementType::Wedge15) {
        eval_wedge15_polynomial(x,
                                 y,
                                 z,
                                 values_out,
                                 gradients_out,
                                 hessians_out);
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
    require_output_span_size(values_out.size(), size_, "SerendipityBasis::evaluate_values_to");
    evaluate_all_to(xi, values_out, std::span<Gradient>{}, std::span<Hessian>{});
}

void SerendipityBasis::evaluate_gradients_to(const math::Vector<Real, 3>& xi,
                                             std::span<Gradient> gradients_out) const {
    require_output_span_size(gradients_out.size(), size_, "SerendipityBasis::evaluate_gradients_to");
    evaluate_all_to(xi, std::span<Real>{}, gradients_out, std::span<Hessian>{});
}

void SerendipityBasis::evaluate_hessians_to(const math::Vector<Real, 3>& xi,
                                            std::span<Hessian> hessians_out) const {
    require_output_span_size(hessians_out.size(), size_, "SerendipityBasis::evaluate_hessians_to");
    evaluate_all_to(xi, std::span<Real>{}, std::span<Gradient>{}, hessians_out);
}

} // namespace basis
} // namespace FE
} // namespace svmp
