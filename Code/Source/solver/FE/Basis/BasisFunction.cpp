/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#include "BasisFunction.h"
#include "VectorBasisEvaluationHelpers.h"
#include <algorithm>
#include <iomanip>
#include <limits>
#include <sstream>

namespace svmp {
namespace FE {
namespace basis {

namespace {

struct BasisFunctionScratch {
    std::vector<Real> scalar_values;
    std::vector<Gradient> scalar_gradients;
    std::vector<Hessian> scalar_hessians;
    std::vector<math::Vector<Real, 3>> vector_values;
    std::vector<VectorJacobian> vector_jacobians;
    std::vector<math::Vector<Real, 3>> vector_curls;
    std::vector<Real> vector_divergences;

    void prewarm(std::size_t max_size) {
        scalar_values.reserve(max_size);
        scalar_gradients.reserve(max_size);
        scalar_hessians.reserve(max_size);
        vector_values.reserve(max_size);
        vector_jacobians.reserve(max_size);
        vector_curls.reserve(max_size);
        vector_divergences.reserve(max_size);
    }
};

BasisFunctionScratch& basis_function_scratch() {
    // Scratch is intentionally thread-local: production assembly uses a
    // persistent worker-thread team, so buffers stay warm on each worker.
    static thread_local BasisFunctionScratch scratch;
    return scratch;
}

void mix_identity_hash_word(std::uint64_t word,
                            std::uint64_t& hash_a,
                            std::uint64_t& hash_b) noexcept {
    hash_a ^= word + 0x9e3779b97f4a7c15ULL + (hash_a << 6u) + (hash_a >> 2u);
    hash_b ^= (word + 0xbf58476d1ce4e5b9ULL) + (hash_b << 7u) + (hash_b >> 3u);
}

} // namespace

BasisIdentityFingerprint
compute_basis_identity_fingerprint(std::span<const std::uint64_t> words) noexcept {
    BasisIdentityFingerprint fingerprint{0x243f6a8885a308d3ULL,
                                         0x13198a2e03707344ULL};
    mix_identity_hash_word(static_cast<std::uint64_t>(words.size()),
                           fingerprint.hash_a,
                           fingerprint.hash_b);
    for (const auto word : words) {
        mix_identity_hash_word(word, fingerprint.hash_a, fingerprint.hash_b);
    }
    return fingerprint;
}

std::string BasisFunction::cache_identity() const {
    std::ostringstream oss;
    oss << "basis=" << static_cast<int>(basis_type())
        << "|elem=" << static_cast<int>(element_type())
        << "|dim=" << dimension()
        << "|order=" << order()
        << "|size=" << size()
        << "|vector=" << is_vector_valued();
    return oss.str();
}

bool BasisFunction::cache_identity_words(std::vector<std::uint64_t>& words) const {
    (void)words;
    return false;
}

bool BasisFunction::cache_identity_fingerprint(std::uint64_t& hash_a,
                                               std::uint64_t& hash_b) const {
    (void)hash_a;
    (void)hash_b;
    return false;
}

void prewarm_basis_function_scratch(std::size_t max_size,
                                    std::size_t max_qpts) {
    (void)max_qpts;
    basis_function_scratch().prewarm(max_size);
}

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
    auto& tmp = basis_function_scratch().scalar_values;
    tmp.resize(size());
    evaluate_values(xi, tmp);
    std::copy_n(tmp.data(), tmp.size(), values_out);
}

void BasisFunction::evaluate_gradients_to(const math::Vector<Real, 3>& xi,
                                          Real* SVMP_RESTRICT gradients_out) const {
    auto& tmp = basis_function_scratch().scalar_gradients;
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
    auto& tmp = basis_function_scratch().scalar_hessians;
    tmp.resize(size());
    evaluate_hessians(xi, tmp);
    for (std::size_t i = 0; i < tmp.size(); ++i) {
        store_hessian(tmp[i], hessians_out + i * 9u);
    }
}

void BasisFunction::evaluate_at_quadrature_points(
    const std::vector<math::Vector<Real, 3>>& points,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) const {
    evaluate_at_quadrature_points_strided(
        points, points.size(), values_out, gradients_out, hessians_out);
}

void BasisFunction::evaluate_at_quadrature_points_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) const {
    const std::size_t num_qpts = points.size();
    const std::size_t num_dofs = size();
    if (output_stride < num_qpts) {
        throw BasisConfigurationException(
            "BasisFunction strided evaluation requires output_stride >= points.size()",
            __FILE__, __LINE__, __func__);
    }

    auto& scratch = basis_function_scratch();
    auto& v_tmp = scratch.scalar_values;
    auto& g_tmp = scratch.scalar_gradients;
    auto& h_tmp = scratch.scalar_hessians;
    if (values_out) v_tmp.resize(num_dofs);
    if (gradients_out) g_tmp.resize(num_dofs);
    if (hessians_out) h_tmp.resize(num_dofs);

    for (std::size_t q = 0; q < num_qpts; ++q) {
        if (values_out && gradients_out && hessians_out) {
            evaluate_all(points[q], v_tmp, g_tmp, h_tmp);
        } else {
            if (values_out) evaluate_values(points[q], v_tmp);
            if (gradients_out) evaluate_gradients(points[q], g_tmp);
            if (hessians_out) evaluate_hessians(points[q], h_tmp);
        }

        if (values_out) {
            for (std::size_t dof = 0; dof < num_dofs; ++dof) {
                values_out[dof * output_stride + q] = v_tmp[dof];
            }
        }
        if (gradients_out) {
            for (std::size_t dof = 0; dof < num_dofs; ++dof) {
                for (std::size_t component = 0; component < 3u; ++component) {
                    gradients_out[(dof * 3u + component) * output_stride + q] =
                        g_tmp[dof][component];
                }
            }
        }
        if (hessians_out) {
            for (std::size_t dof = 0; dof < num_dofs; ++dof) {
                store_hessian_strided(
                    h_tmp[dof], hessians_out + dof * 9u * output_stride, output_stride, q);
            }
        }
    }
}

void BasisFunction::fill_scalar_cache_entry(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out) const {
    evaluate_at_quadrature_points_strided(
        points, output_stride, values_out, gradients_out, hessians_out);
}

void BasisFunction::evaluate_vector_at_quadrature_points(
    const std::vector<math::Vector<Real, 3>>& points,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT jacobians_out,
    Real* SVMP_RESTRICT curls_out,
    Real* SVMP_RESTRICT divergence_out) const {
    evaluate_vector_at_quadrature_points_strided(
        points, points.size(), values_out, jacobians_out, curls_out, divergence_out);
}

void BasisFunction::evaluate_vector_at_quadrature_points_strided(
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT jacobians_out,
    Real* SVMP_RESTRICT curls_out,
    Real* SVMP_RESTRICT divergence_out) const {
    const std::size_t num_qpts = points.size();
    const std::size_t num_dofs = size();
    detail::vector_common::validate_vector_strided_outputs(
        num_qpts, output_stride, "BasisFunction");

    auto& scratch = basis_function_scratch();
    auto& v_tmp = scratch.vector_values;
    auto& j_tmp = scratch.vector_jacobians;
    auto& c_tmp = scratch.vector_curls;
    auto& d_tmp = scratch.vector_divergences;
    if (values_out) v_tmp.resize(num_dofs);
    if (jacobians_out) j_tmp.resize(num_dofs);
    if (curls_out) c_tmp.resize(num_dofs);
    if (divergence_out) d_tmp.resize(num_dofs);

    for (std::size_t q = 0; q < num_qpts; ++q) {
        if (values_out) {
            evaluate_vector_values(points[q], v_tmp);
            detail::vector_common::write_vector_values_strided(
                v_tmp, num_dofs, output_stride, q, values_out);
        }

        if (jacobians_out) {
            evaluate_vector_jacobians(points[q], j_tmp);
            detail::vector_common::write_vector_jacobians_strided(
                j_tmp, num_dofs, output_stride, q, jacobians_out);
        }

        if (curls_out) {
            evaluate_curl(points[q], c_tmp);
            detail::vector_common::write_vector_curl_strided(
                c_tmp, num_dofs, output_stride, q, curls_out);
        }

        if (divergence_out) {
            evaluate_divergence(points[q], d_tmp);
            detail::vector_common::write_vector_divergence_strided(
                d_tmp, num_dofs, output_stride, q, divergence_out);
        }
    }
}

void BasisFunction::evaluate_vector_values(
    const math::Vector<Real, 3>&,
    std::vector<math::Vector<Real, 3>>&) const {
    throw BasisEvaluationException("Vector-valued evaluation requested on scalar basis",
                                   __FILE__, __LINE__, __func__);
}

void BasisFunction::evaluate_vector_jacobians(
    const math::Vector<Real, 3>&,
    std::vector<VectorJacobian>&) const {
    throw BasisEvaluationException("Vector-basis Jacobian evaluation requested on scalar basis",
                                   __FILE__, __LINE__, __func__);
}

void BasisFunction::evaluate_divergence(
    const math::Vector<Real, 3>&,
    std::vector<Real>&) const {
    throw BasisEvaluationException("Divergence requested on scalar basis",
                                   __FILE__, __LINE__, __func__);
}

void BasisFunction::evaluate_curl(
    const math::Vector<Real, 3>&,
    std::vector<math::Vector<Real, 3>>&) const {
    throw BasisEvaluationException("Curl requested on scalar basis",
                                   __FILE__, __LINE__, __func__);
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
        const std::size_t idx = static_cast<std::size_t>(d);
        forward[idx] += eps;
        backward[idx] -= eps;

        std::vector<Real> fwd, bwd;
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
        const std::size_t col = static_cast<std::size_t>(d);
        forward[col] += eps;
        backward[col] -= eps;

        std::vector<Gradient> g_forward, g_backward;
        evaluate_gradients(forward, g_forward);
        evaluate_gradients(backward, g_backward);

        for (std::size_t i = 0; i < base_grad.size(); ++i) {
            for (int k = 0; k < dimension(); ++k) {
                const std::size_t row = static_cast<std::size_t>(k);
                hessians[i](row, col) = (g_forward[i][row] - g_backward[i][row]) / (Real(2) * eps);
            }
        }
    }
}

} // namespace basis
} // namespace FE
} // namespace svmp
