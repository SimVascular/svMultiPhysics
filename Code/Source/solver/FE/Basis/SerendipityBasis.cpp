/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#include "SerendipityBasis.h"
#include "LagrangeBasis.h"
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

    if (nodes.size() > total_size) {
        throw BasisConstructionException(
            "SerendipityBasis: quadrilateral serendipity boundary nodes exceed requested size",
            __FILE__, __LINE__, __func__);
    }

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

    if (interior_count > interior_candidates.size()) {
        throw BasisConstructionException(
            "SerendipityBasis: insufficient quadrilateral interior nodes for requested serendipity order",
            __FILE__, __LINE__, __func__);
    }

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
    if (n == 0 || exponents.size() != nodes.size()) {
        throw BasisConstructionException(
            "SerendipityBasis: invalid quadrilateral serendipity interpolation setup",
            __FILE__, __LINE__, __func__);
    }

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

void eval_hex20_internal(Real r, Real s, Real t, Real* internal_vals) {
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

void eval_hex20_grad_internal(Real r, Real s, Real t, Gradient* internal_grads) {
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

void eval_hex20_hess_internal(Real r, Real s, Real t, Hessian* internal_hessians) {
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
        Hessian H{};
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
                             Real* values,
                             Gradient* gradients,
                             Hessian* hessians) {
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

        if (values) {
            phi[j] = ra * sb * tc;
        }
        if (gradients) {
            dr[j] = (a > 0) ? Real(a) * rp[ar - 1u] * sb * tc : Real(0);
            ds[j] = (b > 0) ? ra * Real(b) * sp[bs - 1u] * tc : Real(0);
            dt[j] = (c > 0) ? ra * sb * Real(c) * tp[ct - 1u] : Real(0);
        }
        if (hessians) {
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
        Hessian H{};
        for (int j = 0; j < 15; ++j) {
            const Real coefficient =
                kWedge15Coefficients[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)];
            if (values) {
                value += coefficient * phi[j];
            }
            if (gradients) {
                gr += coefficient * dr[j];
                gs += coefficient * ds[j];
                gt += coefficient * dt[j];
            }
            if (hessians) {
                H(0, 0) += coefficient * drr[j];
                H(1, 1) += coefficient * dss[j];
                H(2, 2) += coefficient * dtt[j];
                H(0, 1) += coefficient * drs[j];
                H(0, 2) += coefficient * drt[j];
                H(1, 2) += coefficient * dst[j];
            }
        }

        const std::size_t index = static_cast<std::size_t>(i);
        if (values) {
            values[index] = value;
        }
        if (gradients) {
            gradients[index][0] = gr;
            gradients[index][1] = gs;
            gradients[index][2] = gt;
        }
        if (hessians) {
            H(1, 0) = H(0, 1);
            H(2, 0) = H(0, 2);
            H(2, 1) = H(1, 2);
            hessians[index] = H;
        }
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
        if (type == ElementType::Quad8 && order_ != 2) {
            throw BasisConfigurationException(
                "SerendipityBasis: Quad8 is only valid for quadratic order 2; use Quad4 for higher-order quadrilateral serendipity",
                __FILE__, __LINE__, __func__);
        }
        quad_monomial_exponents_ = quad_serendipity_exponents(order_);
        size_ = quad_monomial_exponents_.size();
        nodes_ = quad_serendipity_nodes(order_, size_);
        if (nodes_.size() != size_) {
            throw BasisConstructionException(
                "SerendipityBasis: quadrilateral serendipity setup produced inconsistent sizes",
                __FILE__, __LINE__, __func__);
        }
        quad_inv_vandermonde_ = quad_serendipity_inverse_vandermonde(nodes_, quad_monomial_exponents_, order_);
    } else if (type == ElementType::Hex8 || type == ElementType::Hex20) {
        dimension_ = 3;
        if (order_ < 1) order_ = 1;
        if (order_ == 1) {
            size_ = 8;
        } else if (order_ == 2) {
            size_ = 20;
        } else {
            throw BasisConfigurationException(
                "SerendipityBasis supports up to quadratic on hexahedra",
                __FILE__, __LINE__, __func__);
        }
    } else if (type == ElementType::Wedge15) {
        dimension_ = 3;
        if (order_ < 2) {
            order_ = 2;
        }
        if (order_ == 2) {
            size_ = 15;
        } else {
            throw BasisConfigurationException(
                "SerendipityBasis supports up to quadratic on wedge15",
                __FILE__, __LINE__, __func__);
        }
    } else {
        throw BasisElementCompatibilityException("SerendipityBasis supports Quad4/Quad8, Hex8/Hex20, and Wedge15 elements",
                                                 __FILE__, __LINE__, __func__);
    }

    if (nodes_.empty()) {
        nodes_.reserve(size_);
        for (std::size_t i = 0; i < size_; ++i) {
            nodes_.push_back(ReferenceNodeLayout::get_node_coords(element_type_, i));
        }
    }
}

void SerendipityBasis::evaluate_values(const math::Vector<Real, 3>& xi,
                                       std::vector<Real>& values) const {
    values.assign(size_, Real(0));
    const Real x = xi[0];
    const Real y = xi[1];
    const Real z = xi[2];

    if (dimension_ == 2) {
        if (quad_monomial_exponents_.size() != size_ ||
            quad_inv_vandermonde_.size() != size_ * size_) {
            throw BasisEvaluationException(
                "SerendipityBasis: quadrilateral interpolation tables are not initialized for value evaluation",
                __FILE__, __LINE__, __func__);
        }

        std::vector<Real> monomials(size_, Real(0));
        for (std::size_t j = 0; j < size_; ++j) {
            const auto [ax, ay] = quad_monomial_exponents_[j];
            monomials[j] = std::pow(x, ax) * std::pow(y, ay);
        }

        for (std::size_t i = 0; i < size_; ++i) {
            Real value = Real(0);
            for (std::size_t j = 0; j < size_; ++j) {
                value += monomials[j] * quad_inv_vandermonde_[j * size_ + i];
            }
            values[i] = value;
        }
        return;
    }

    if (dimension_ == 3 && order_ == 1) {
        // Hex8 trilinear shape functions
        const Real r = x;
        const Real s = y;
        const Real t = z;
        values[0] = Real(0.125) * (Real(1) - r) * (Real(1) - s) * (Real(1) - t);
        values[1] = Real(0.125) * (Real(1) + r) * (Real(1) - s) * (Real(1) - t);
        values[2] = Real(0.125) * (Real(1) + r) * (Real(1) + s) * (Real(1) - t);
        values[3] = Real(0.125) * (Real(1) - r) * (Real(1) + s) * (Real(1) - t);
        values[4] = Real(0.125) * (Real(1) - r) * (Real(1) - s) * (Real(1) + t);
        values[5] = Real(0.125) * (Real(1) + r) * (Real(1) - s) * (Real(1) + t);
        values[6] = Real(0.125) * (Real(1) + r) * (Real(1) + s) * (Real(1) + t);
        values[7] = Real(0.125) * (Real(1) - r) * (Real(1) + s) * (Real(1) + t);
        return;
    }

    const Real r = x;
    const Real s = y;
    const Real t = z;

    if (geometry_mode_ && element_type_ == ElementType::Hex20) {
        // Hex20 geometry mode: use trilinear Hex8 shape functions on corners, edges zero.
        values[0] = Real(0.125) * (Real(1) - r) * (Real(1) - s) * (Real(1) - t);
        values[1] = Real(0.125) * (Real(1) + r) * (Real(1) - s) * (Real(1) - t);
        values[2] = Real(0.125) * (Real(1) + r) * (Real(1) + s) * (Real(1) - t);
        values[3] = Real(0.125) * (Real(1) - r) * (Real(1) + s) * (Real(1) - t);
        values[4] = Real(0.125) * (Real(1) - r) * (Real(1) - s) * (Real(1) + t);
        values[5] = Real(0.125) * (Real(1) + r) * (Real(1) - s) * (Real(1) + t);
        values[6] = Real(0.125) * (Real(1) + r) * (Real(1) + s) * (Real(1) + t);
        values[7] = Real(0.125) * (Real(1) - r) * (Real(1) + s) * (Real(1) + t);
        for (std::size_t i = 8; i < 20; ++i) {
            values[i] = Real(0);
        }
        return;
    }

    if (element_type_ == ElementType::Hex20) {
        Real internal_vals[20];
        eval_hex20_internal(r, s, t, internal_vals);
        const auto mesh_to_basis = ReferenceNodeLayout::mesh_to_basis_ordering(element_type_);
        BASIS_CHECK_EVAL(mesh_to_basis.size() == size_,
                         "Hex20 mesh-to-basis ordering is not registered");
        for (std::size_t i = 0; i < 20; ++i) {
            values[i] = internal_vals[mesh_to_basis[i]];
        }
        return;
    }

    if (element_type_ == ElementType::Wedge15) {
        eval_wedge15_polynomial(r, s, t, values.data(), nullptr, nullptr);
        return;
    }

}

void SerendipityBasis::evaluate_gradients(const math::Vector<Real, 3>& xi,
                                          std::vector<Gradient>& gradients) const {
    gradients.assign(size_, Gradient{});

    const Real x = xi[0];
    const Real y = xi[1];
    const Real z = xi[2];

    if (dimension_ == 2) {
        if (quad_monomial_exponents_.size() != size_ ||
            quad_inv_vandermonde_.size() != size_ * size_) {
            throw BasisEvaluationException(
                "SerendipityBasis: quadrilateral interpolation tables are not initialized for gradient evaluation",
                __FILE__, __LINE__, __func__);
        }

        std::vector<Real> dmon_dx(size_, Real(0));
        std::vector<Real> dmon_dy(size_, Real(0));
        for (std::size_t j = 0; j < size_; ++j) {
            const auto [ax, ay] = quad_monomial_exponents_[j];
            dmon_dx[j] =
                (ax > 0) ? Real(ax) * std::pow(x, ax - 1) * std::pow(y, ay) : Real(0);
            dmon_dy[j] =
                (ay > 0) ? std::pow(x, ax) * Real(ay) * std::pow(y, ay - 1) : Real(0);
        }

        for (std::size_t i = 0; i < size_; ++i) {
            Real gx = Real(0);
            Real gy = Real(0);
            for (std::size_t j = 0; j < size_; ++j) {
                const Real coeff = quad_inv_vandermonde_[j * size_ + i];
                gx += dmon_dx[j] * coeff;
                gy += dmon_dy[j] * coeff;
            }
            gradients[i][0] = gx;
            gradients[i][1] = gy;
        }
        return;
    }

    // 3D linear hex (Hex8)
    if (dimension_ == 3 && order_ == 1) {
        const Real r = x, s = y, t = z;
        gradients[0][0] = -Real(0.125) * (Real(1) - s) * (Real(1) - t);
        gradients[0][1] = -Real(0.125) * (Real(1) - r) * (Real(1) - t);
        gradients[0][2] = -Real(0.125) * (Real(1) - r) * (Real(1) - s);

        gradients[1][0] =  Real(0.125) * (Real(1) - s) * (Real(1) - t);
        gradients[1][1] = -Real(0.125) * (Real(1) + r) * (Real(1) - t);
        gradients[1][2] = -Real(0.125) * (Real(1) + r) * (Real(1) - s);

        gradients[2][0] =  Real(0.125) * (Real(1) + s) * (Real(1) - t);
        gradients[2][1] =  Real(0.125) * (Real(1) + r) * (Real(1) - t);
        gradients[2][2] = -Real(0.125) * (Real(1) + r) * (Real(1) + s);

        gradients[3][0] = -Real(0.125) * (Real(1) + s) * (Real(1) - t);
        gradients[3][1] =  Real(0.125) * (Real(1) - r) * (Real(1) - t);
        gradients[3][2] = -Real(0.125) * (Real(1) - r) * (Real(1) + s);

        gradients[4][0] = -Real(0.125) * (Real(1) - s) * (Real(1) + t);
        gradients[4][1] = -Real(0.125) * (Real(1) - r) * (Real(1) + t);
        gradients[4][2] =  Real(0.125) * (Real(1) - r) * (Real(1) - s);

        gradients[5][0] =  Real(0.125) * (Real(1) - s) * (Real(1) + t);
        gradients[5][1] = -Real(0.125) * (Real(1) + r) * (Real(1) + t);
        gradients[5][2] =  Real(0.125) * (Real(1) + r) * (Real(1) - s);

        gradients[6][0] =  Real(0.125) * (Real(1) + s) * (Real(1) + t);
        gradients[6][1] =  Real(0.125) * (Real(1) + r) * (Real(1) + t);
        gradients[6][2] =  Real(0.125) * (Real(1) + r) * (Real(1) + s);

        gradients[7][0] = -Real(0.125) * (Real(1) + s) * (Real(1) + t);
        gradients[7][1] =  Real(0.125) * (Real(1) - r) * (Real(1) + t);
        gradients[7][2] =  Real(0.125) * (Real(1) - r) * (Real(1) + s);
        return;
    }

    // Hex20 geometry mode: use Hex8 gradients
    if (dimension_ == 3 && order_ == 2 && geometry_mode_ &&
        (element_type_ == ElementType::Hex20 || element_type_ == ElementType::Quad8)) {
        const Real r = x, s = y, t = z;
        gradients[0][0] = -Real(0.125) * (Real(1) - s) * (Real(1) - t);
        gradients[0][1] = -Real(0.125) * (Real(1) - r) * (Real(1) - t);
        gradients[0][2] = -Real(0.125) * (Real(1) - r) * (Real(1) - s);

        gradients[1][0] =  Real(0.125) * (Real(1) - s) * (Real(1) - t);
        gradients[1][1] = -Real(0.125) * (Real(1) + r) * (Real(1) - t);
        gradients[1][2] = -Real(0.125) * (Real(1) + r) * (Real(1) - s);

        gradients[2][0] =  Real(0.125) * (Real(1) + s) * (Real(1) - t);
        gradients[2][1] =  Real(0.125) * (Real(1) + r) * (Real(1) - t);
        gradients[2][2] = -Real(0.125) * (Real(1) + r) * (Real(1) + s);

        gradients[3][0] = -Real(0.125) * (Real(1) + s) * (Real(1) - t);
        gradients[3][1] =  Real(0.125) * (Real(1) - r) * (Real(1) - t);
        gradients[3][2] = -Real(0.125) * (Real(1) - r) * (Real(1) + s);

        gradients[4][0] = -Real(0.125) * (Real(1) - s) * (Real(1) + t);
        gradients[4][1] = -Real(0.125) * (Real(1) - r) * (Real(1) + t);
        gradients[4][2] =  Real(0.125) * (Real(1) - r) * (Real(1) - s);

        gradients[5][0] =  Real(0.125) * (Real(1) - s) * (Real(1) + t);
        gradients[5][1] = -Real(0.125) * (Real(1) + r) * (Real(1) + t);
        gradients[5][2] =  Real(0.125) * (Real(1) + r) * (Real(1) - s);

        gradients[6][0] =  Real(0.125) * (Real(1) + s) * (Real(1) + t);
        gradients[6][1] =  Real(0.125) * (Real(1) + r) * (Real(1) + t);
        gradients[6][2] =  Real(0.125) * (Real(1) + r) * (Real(1) + s);

        gradients[7][0] = -Real(0.125) * (Real(1) + s) * (Real(1) + t);
        gradients[7][1] =  Real(0.125) * (Real(1) - r) * (Real(1) + t);
        gradients[7][2] =  Real(0.125) * (Real(1) - r) * (Real(1) + s);
        // Edge-node gradients remain zero
        return;
    }

    // Hex20 analytical gradients using monomial differentiation
    if (element_type_ == ElementType::Hex20 && order_ == 2) {
        const Real r = x, s = y, t = z;
        Gradient internal_grads[20];
        eval_hex20_grad_internal(r, s, t, internal_grads);
        const auto mesh_to_basis = ReferenceNodeLayout::mesh_to_basis_ordering(element_type_);
        BASIS_CHECK_EVAL(mesh_to_basis.size() == size_,
                         "Hex20 mesh-to-basis ordering is not registered");
        for (std::size_t i = 0; i < 20; ++i) {
            gradients[i] = internal_grads[mesh_to_basis[i]];
        }
        return;
    }

    // Wedge15 analytical gradients using monomial differentiation
    if (element_type_ == ElementType::Wedge15 && order_ == 2) {
        eval_wedge15_polynomial(x, y, z, nullptr, gradients.data(), nullptr);
        return;
    }

    throw BasisEvaluationException("SerendipityBasis::evaluate_gradients: unsupported serendipity configuration",
                                   __FILE__, __LINE__, __func__);
}

void SerendipityBasis::evaluate_hessians(const math::Vector<Real, 3>& xi,
                                         std::vector<Hessian>& hessians) const {
    hessians.assign(size_, Hessian{});
    const Real x = xi[0];
    const Real y = xi[1];
    const Real z = xi[2];

    if (dimension_ == 2) {
        if (quad_monomial_exponents_.size() != size_ ||
            quad_inv_vandermonde_.size() != size_ * size_) {
            throw BasisEvaluationException(
                "SerendipityBasis: quadrilateral interpolation tables are not initialized for Hessian evaluation",
                __FILE__, __LINE__, __func__);
        }

        std::vector<Real> dxx(size_, Real(0));
        std::vector<Real> dxy(size_, Real(0));
        std::vector<Real> dyy(size_, Real(0));
        for (std::size_t j = 0; j < size_; ++j) {
            const auto [ax, ay] = quad_monomial_exponents_[j];
            dxx[j] = (ax > 1)
                         ? Real(ax * (ax - 1)) * std::pow(x, ax - 2) * std::pow(y, ay)
                         : Real(0);
            dxy[j] = (ax > 0 && ay > 0)
                         ? Real(ax * ay) * std::pow(x, ax - 1) * std::pow(y, ay - 1)
                         : Real(0);
            dyy[j] = (ay > 1)
                         ? Real(ay * (ay - 1)) * std::pow(x, ax) * std::pow(y, ay - 2)
                         : Real(0);
        }

        for (std::size_t i = 0; i < size_; ++i) {
            for (std::size_t j = 0; j < size_; ++j) {
                const Real coeff = quad_inv_vandermonde_[j * size_ + i];
                hessians[i](0, 0) += dxx[j] * coeff;
                hessians[i](0, 1) += dxy[j] * coeff;
                hessians[i](1, 1) += dyy[j] * coeff;
            }
            hessians[i](1, 0) = hessians[i](0, 1);
        }
        return;
    }

    if (element_type_ == ElementType::Hex8 && order_ == 1) {
        static const LagrangeBasis parent(ElementType::Hex8, 1);
        parent.evaluate_hessians(xi, hessians);
        return;
    }

    if (geometry_mode_ && element_type_ == ElementType::Hex20) {
        static const LagrangeBasis parent(ElementType::Hex8, 1);
        std::array<Real, 8u * 9u> parent_hessians{};
        parent.evaluate_hessians_to(xi, parent_hessians.data());
        for (std::size_t i = 0; i < 8; ++i) {
            for (std::size_t r = 0; r < 3; ++r) {
                for (std::size_t c = 0; c < 3; ++c) {
                    hessians[i](r, c) = parent_hessians[i * 9u + r * 3u + c];
                }
            }
        }
        return;
    }

    if (element_type_ == ElementType::Hex20 && order_ == 2) {
        Hessian internal_hessians[20];
        eval_hex20_hess_internal(x, y, z, internal_hessians);
        const auto mesh_to_basis = ReferenceNodeLayout::mesh_to_basis_ordering(element_type_);
        BASIS_CHECK_EVAL(mesh_to_basis.size() == size_,
                         "Hex20 mesh-to-basis ordering is not registered");
        for (std::size_t i = 0; i < 20; ++i) {
            hessians[i] = internal_hessians[mesh_to_basis[i]];
        }
        return;
    }

    if (element_type_ == ElementType::Wedge15 && order_ == 2) {
        eval_wedge15_polynomial(x, y, z, nullptr, nullptr, hessians.data());
        return;
    }

    throw BasisEvaluationException("SerendipityBasis::evaluate_hessians: unsupported serendipity configuration",
                                   __FILE__, __LINE__, __func__);
}

} // namespace basis
} // namespace FE
} // namespace svmp
