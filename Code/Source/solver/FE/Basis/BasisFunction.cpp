// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "BasisFunction.h"

#include <algorithm>

namespace svmp {
namespace FE {
namespace basis {

namespace {

struct BasisFunctionScratch {
    std::vector<Real> values;
    std::vector<Gradient> gradients;
    std::vector<Hessian> hessians;
};

BasisFunctionScratch& scratch() {
    static thread_local BasisFunctionScratch data;
    return data;
}

} // namespace

void BasisFunction::evaluate_gradients(const math::Vector<Real, 3>& xi,
                                       std::vector<Gradient>& gradients) const {
    (void)xi;
    (void)gradients;
    throw BasisEvaluationException("Analytic gradient evaluation is not implemented for this basis",
                                   __FILE__, __LINE__, __func__);
}

void BasisFunction::evaluate_hessians(const math::Vector<Real, 3>& xi,
                                      std::vector<Hessian>& hessians) const {
    (void)xi;
    (void)hessians;
    throw BasisEvaluationException("Analytic Hessian evaluation is not implemented for this basis",
                                   __FILE__, __LINE__, __func__);
}

void BasisFunction::evaluate_all(const math::Vector<Real, 3>& xi,
                                 std::vector<Real>& values,
                                 std::vector<Gradient>& gradients,
                                 std::vector<Hessian>& hessians) const {
    evaluate_values(xi, values);
    evaluate_gradients(xi, gradients);
    evaluate_hessians(xi, hessians);
}

void BasisFunction::evaluate_values_to(const math::Vector<Real, 3>& xi,
                                       Real* SVMP_RESTRICT values_out) const {
    auto& tmp = scratch().values;
    tmp.resize(size());
    evaluate_values(xi, tmp);
    std::copy_n(tmp.data(), tmp.size(), values_out);
}

void BasisFunction::evaluate_gradients_to(const math::Vector<Real, 3>& xi,
                                          Real* SVMP_RESTRICT gradients_out) const {
    auto& tmp = scratch().gradients;
    tmp.resize(size());
    evaluate_gradients(xi, tmp);
    for (std::size_t i = 0; i < tmp.size(); ++i) {
        gradients_out[i * 3u + 0u] = tmp[i][0];
        gradients_out[i * 3u + 1u] = tmp[i][1];
        gradients_out[i * 3u + 2u] = tmp[i][2];
    }
}

void BasisFunction::evaluate_hessians_to(const math::Vector<Real, 3>& xi,
                                         Real* SVMP_RESTRICT hessians_out) const {
    auto& tmp = scratch().hessians;
    tmp.resize(size());
    evaluate_hessians(xi, tmp);
    for (std::size_t i = 0; i < tmp.size(); ++i) {
        store_hessian(tmp[i], hessians_out + i * 9u);
    }
}

void BasisFunction::numerical_gradient(const math::Vector<Real, 3>& xi,
                                       std::vector<Gradient>& gradients,
                                       Real eps) const {
    std::vector<Real> base;
    evaluate_values(xi, base);
    gradients.assign(base.size(), Gradient{});

    for (int d = 0; d < dimension(); ++d) {
        math::Vector<Real, 3> forward = xi;
        math::Vector<Real, 3> backward = xi;
        const auto idx = static_cast<std::size_t>(d);
        forward[idx] += eps;
        backward[idx] -= eps;

        std::vector<Real> fwd;
        std::vector<Real> bwd;
        evaluate_values(forward, fwd);
        evaluate_values(backward, bwd);

        for (std::size_t i = 0; i < base.size(); ++i) {
            gradients[i][idx] = (fwd[i] - bwd[i]) / (Real(2) * eps);
        }
    }
}

void BasisFunction::numerical_hessian(const math::Vector<Real, 3>& xi,
                                      std::vector<Hessian>& hessians,
                                      Real eps) const {
    std::vector<Gradient> base_grad;
    evaluate_gradients(xi, base_grad);
    hessians.assign(base_grad.size(), Hessian{});

    for (int d = 0; d < dimension(); ++d) {
        math::Vector<Real, 3> forward = xi;
        math::Vector<Real, 3> backward = xi;
        const auto col = static_cast<std::size_t>(d);
        forward[col] += eps;
        backward[col] -= eps;

        std::vector<Gradient> g_forward;
        std::vector<Gradient> g_backward;
        evaluate_gradients(forward, g_forward);
        evaluate_gradients(backward, g_backward);

        for (std::size_t i = 0; i < base_grad.size(); ++i) {
            for (int k = 0; k < dimension(); ++k) {
                const auto row = static_cast<std::size_t>(k);
                hessians[i](row, col) =
                    (g_forward[i][row] - g_backward[i][row]) / (Real(2) * eps);
            }
        }
    }
}

} // namespace basis
} // namespace FE
} // namespace svmp
