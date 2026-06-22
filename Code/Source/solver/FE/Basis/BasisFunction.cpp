// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "BasisFunction.h"

#include <algorithm>
#include <string>

namespace svmp {
namespace FE {
namespace basis {

void require_span_size(std::size_t actual,
                       std::size_t expected,
                       const char* label) {
    svmp::throw_if<BasisEvaluationException>(actual < expected, SVMP_HERE,
        std::string(label) + ": output span is smaller than basis size");
}

const std::vector<math::Vector<double, 3>>& BasisFunction::nodes() const noexcept {
    // Default for bases that do not expose interpolation nodes; nodal families
    // (LagrangeBasis, SerendipityBasis) override this to return their layout.
    static const std::vector<math::Vector<double, 3>> kNoNodes;
    return kNoNodes;
}

void BasisFunction::evaluate_gradients(const math::Vector<double, 3>& xi,
                                       std::vector<Gradient>& gradients) const {
    (void)xi;
    (void)gradients;
    svmp::raise<BasisEvaluationException>(SVMP_HERE,
        "Analytic gradient evaluation is not implemented for this basis");
}

void BasisFunction::evaluate_hessians(const math::Vector<double, 3>& xi,
                                      std::vector<Hessian>& hessians) const {
    (void)xi;
    (void)hessians;
    svmp::raise<BasisEvaluationException>(SVMP_HERE,
        "Analytic Hessian evaluation is not implemented for this basis");
}

void BasisFunction::evaluate_all(const math::Vector<double, 3>& xi,
                                 std::vector<double>& values,
                                 std::vector<Gradient>& gradients,
                                 std::vector<Hessian>& hessians) const {
    evaluate_values(xi, values);
    evaluate_gradients(xi, gradients);
    evaluate_hessians(xi, hessians);
}

void BasisFunction::evaluate_values_to(const math::Vector<double, 3>& xi,
                                       std::span<double> values_out) const {
    require_span_size(values_out.size(), size(), "BasisFunction::evaluate_values_to");
    std::vector<double> tmp(size());
    evaluate_values(xi, tmp);
    std::copy_n(tmp.begin(), tmp.size(), values_out.begin());
}

void BasisFunction::evaluate_gradients_to(const math::Vector<double, 3>& xi,
                                          std::span<Gradient> gradients_out) const {
    require_span_size(gradients_out.size(), size(), "BasisFunction::evaluate_gradients_to");
    std::vector<Gradient> tmp(size());
    evaluate_gradients(xi, tmp);
    std::copy_n(tmp.begin(), tmp.size(), gradients_out.begin());
}

void BasisFunction::evaluate_hessians_to(const math::Vector<double, 3>& xi,
                                         std::span<Hessian> hessians_out) const {
    require_span_size(hessians_out.size(), size(), "BasisFunction::evaluate_hessians_to");
    std::vector<Hessian> tmp(size());
    evaluate_hessians(xi, tmp);
    std::copy_n(tmp.begin(), tmp.size(), hessians_out.begin());
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

} // namespace basis
} // namespace FE
} // namespace svmp
