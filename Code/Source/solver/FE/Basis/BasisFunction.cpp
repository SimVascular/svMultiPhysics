// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "BasisFunction.h"

#include <algorithm>
#include <string>

namespace svmp::FE::basis {

void require_span_size(std::size_t actual,
                       std::size_t expected,
                       const char* label) {
    svmp::throw_if<BasisEvaluationException>(actual < expected, std::string(label) + ": output span is smaller than basis size");
}

const std::vector<math::Vector<double, 3>>& BasisFunction::nodes() const noexcept {
    // Default for bases that do not expose interpolation nodes; nodal families
    // (LagrangeBasis, SerendipityBasis) override this to return their layout.
    static const std::vector<math::Vector<double, 3>> kNoNodes;
    return kNoNodes;
}

// Vector-output overloads: size the container and forward to the matching span
// primitive. Defined once here so concrete families implement only the span
// primitives below.
void BasisFunction::evaluate_values(const math::Vector<double, 3>& xi,
                                    std::vector<double>& values) const {
    values.resize(size());
    evaluate_values_to(xi, std::span<double>(values.data(), values.size()));
}

void BasisFunction::evaluate_gradients(const math::Vector<double, 3>& xi,
                                       std::vector<Gradient>& gradients) const {
    gradients.resize(size());
    evaluate_gradients_to(xi, std::span<Gradient>(gradients.data(), gradients.size()));
}

void BasisFunction::evaluate_hessians(const math::Vector<double, 3>& xi,
                                      std::vector<Hessian>& hessians) const {
    hessians.resize(size());
    evaluate_hessians_to(xi, std::span<Hessian>(hessians.data(), hessians.size()));
}

void BasisFunction::evaluate_all(const math::Vector<double, 3>& xi,
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

// The gradient/Hessian span primitives default to reporting "not implemented"; a
// family supplies analytical derivatives by overriding them. evaluate_values_to
// has no base definition: every basis must provide values.
void BasisFunction::evaluate_gradients_to(const math::Vector<double, 3>& xi,
                                          std::span<Gradient> gradients_out) const {
    (void)xi;
    (void)gradients_out;
    svmp::raise<BasisEvaluationException>("Analytic gradient evaluation is not implemented for this basis");
}

void BasisFunction::evaluate_hessians_to(const math::Vector<double, 3>& xi,
                                         std::span<Hessian> hessians_out) const {
    (void)xi;
    (void)hessians_out;
    svmp::raise<BasisEvaluationException>("Analytic Hessian evaluation is not implemented for this basis");
}

// Combined evaluator default: forward each requested (non-empty) quantity to its
// single-quantity span primitive. Families override this to share per-point setup
// across the requested quantities.
void BasisFunction::evaluate_all_to(const math::Vector<double, 3>& xi,
                                    std::span<double> values_out,
                                    std::span<Gradient> gradients_out,
                                    std::span<Hessian> hessians_out) const {
    if (!values_out.empty()) {
        evaluate_values_to(xi, values_out);
    }
    if (!gradients_out.empty()) {
        evaluate_gradients_to(xi, gradients_out);
    }
    if (!hessians_out.empty()) {
        evaluate_hessians_to(xi, hessians_out);
    }
}

void BasisFunction::numerical_gradient(const math::Vector<double, 3>& xi,
                                       std::vector<Gradient>& gradients,
                                       double eps) const {
    std::vector<double> base;
    evaluate_values(xi, base);
    gradients.assign(base.size(), Gradient::Zero());

    for (int d = 0; d < dimension(); ++d) {
        math::Vector<double, 3> forward = xi;
        math::Vector<double, 3> backward = xi;
        const auto idx = static_cast<std::size_t>(d);
        forward[idx] += eps;
        backward[idx] -= eps;

        std::vector<double> fwd;
        std::vector<double> bwd;
        evaluate_values(forward, fwd);
        evaluate_values(backward, bwd);

        for (std::size_t i = 0; i < base.size(); ++i) {
            gradients[i][idx] = (fwd[i] - bwd[i]) / (double(2) * eps);
        }
    }
}

void BasisFunction::numerical_hessian(const math::Vector<double, 3>& xi,
                                      std::vector<Hessian>& hessians,
                                      double eps) const {
    std::vector<Gradient> base_grad;
    evaluate_gradients(xi, base_grad);
    hessians.assign(base_grad.size(), Hessian::Zero());

    for (int d = 0; d < dimension(); ++d) {
        math::Vector<double, 3> forward = xi;
        math::Vector<double, 3> backward = xi;
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
                    (g_forward[i][row] - g_backward[i][row]) / (double(2) * eps);
            }
        }
    }
}

} // namespace svmp::FE::basis
