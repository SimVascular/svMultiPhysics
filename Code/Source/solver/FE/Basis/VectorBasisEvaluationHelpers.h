/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#ifndef SVMP_FE_BASIS_VECTORBASISEVALUATIONHELPERS_H
#define SVMP_FE_BASIS_VECTORBASISEVALUATIONHELPERS_H

#include "VectorBasis.h"
#include "Basis/BasisTraits.h"

#include <algorithm>
#include <array>
#include <limits>
#include <string>
#include <vector>

namespace svmp {
namespace FE {
namespace basis {
namespace detail {
namespace vector_common {

using Vec3 = math::Vector<Real, 3>;

struct VectorBasisScratch {
    std::vector<Real> px;
    std::vector<Real> py;
    std::vector<Real> pz;
    std::vector<Real> batched_px;
    std::vector<Real> batched_py;
    std::vector<Real> batched_pz;
    std::vector<Real> candidate_values;
    std::vector<Real> candidate_dx;
    std::vector<Real> candidate_dy;
    std::vector<Real> candidate_dz;
    std::vector<Real> modal_values_batched;
    std::vector<Real> modal_jacobians_batched;
    std::vector<Real> modal_curls_batched;
    std::vector<Real> modal_divergence_batched;
    std::vector<Vec3> vector_values;
    std::vector<VectorJacobian> vector_jacobians;
    std::vector<Real> scalars;
    std::vector<Vec3> api_values;
    std::vector<VectorJacobian> api_jacobians;
    std::vector<Vec3> api_curl;
    std::vector<Real> api_divergence;

    void prewarm(std::size_t max_size, std::size_t max_qpts) {
        const std::size_t batched_size = max_size * std::max<std::size_t>(max_qpts, 1u);
        px.reserve(max_size);
        py.reserve(max_size);
        pz.reserve(max_size);
        batched_px.reserve(batched_size);
        batched_py.reserve(batched_size);
        batched_pz.reserve(batched_size);
        candidate_values.reserve(max_size);
        candidate_dx.reserve(max_size);
        candidate_dy.reserve(max_size);
        candidate_dz.reserve(max_size);
        modal_values_batched.reserve(batched_size * 3u);
        modal_jacobians_batched.reserve(batched_size * 9u);
        modal_curls_batched.reserve(batched_size * 3u);
        modal_divergence_batched.reserve(batched_size);
        vector_values.reserve(max_size);
        vector_jacobians.reserve(max_size);
        scalars.reserve(max_size);
        api_values.reserve(max_size);
        api_jacobians.reserve(max_size);
        api_curl.reserve(max_size);
        api_divergence.reserve(max_size);
    }
};

VectorBasisScratch& vector_basis_scratch();
void prewarm_vector_basis_scratch(std::size_t max_size, std::size_t max_qpts = 0);

void fill_powers(Real x, int max_p, std::vector<Real>& out);
void fill_power_tables(const Vec3& xi,
                       const std::array<int, 3>& limits,
                       VectorBasisScratch& scratch);
void fill_batched_power_tables(const std::vector<Vec3>& points,
                               const std::array<int, 3>& limits,
                               VectorBasisScratch& scratch);
void validate_vector_strided_outputs(std::size_t num_qpts,
                                     std::size_t output_stride,
                                     const char* family_name);
void zero_active_strided_rows(Real* output,
                              std::size_t rows,
                              std::size_t output_stride,
                              std::size_t num_qpts);
SparseModalCoefficientMatrix build_sparse_modal_coefficients(
    const std::vector<Real>& dense_coefficients,
    std::size_t rows,
    std::size_t cols);
Vec3 curl_from_jacobian(const VectorJacobian& J) noexcept;
Real divergence_from_jacobian(const VectorJacobian& J) noexcept;

inline Real batched_power_product(const std::vector<Real>& px,
                                  const std::vector<Real>& py,
                                  const std::vector<Real>& pz,
                                  std::size_t stride,
                                  int px_pow,
                                  int py_pow,
                                  int pz_pow,
                                  std::size_t q) noexcept {
    return px[static_cast<std::size_t>(px_pow) * stride + q] *
           py[static_cast<std::size_t>(py_pow) * stride + q] *
           pz[static_cast<std::size_t>(pz_pow) * stride + q];
}

inline Real batched_component_partial(const std::vector<Real>& px,
                                      const std::vector<Real>& py,
                                      const std::vector<Real>& pz,
                                      std::size_t stride,
                                      int px_pow,
                                      int py_pow,
                                      int pz_pow,
                                      int derivative_axis,
                                      std::size_t q) noexcept {
    if (derivative_axis == 0) {
        if (px_pow == 0) {
            return Real(0);
        }
        return Real(px_pow) *
               px[static_cast<std::size_t>(px_pow - 1) * stride + q] *
               py[static_cast<std::size_t>(py_pow) * stride + q] *
               pz[static_cast<std::size_t>(pz_pow) * stride + q];
    }
    if (derivative_axis == 1) {
        if (py_pow == 0) {
            return Real(0);
        }
        return Real(py_pow) *
               px[static_cast<std::size_t>(px_pow) * stride + q] *
               py[static_cast<std::size_t>(py_pow - 1) * stride + q] *
               pz[static_cast<std::size_t>(pz_pow) * stride + q];
    }
    if (pz_pow == 0) {
        return Real(0);
    }
    return Real(pz_pow) *
           px[static_cast<std::size_t>(px_pow) * stride + q] *
           py[static_cast<std::size_t>(py_pow) * stride + q] *
           pz[static_cast<std::size_t>(pz_pow - 1) * stride + q];
}

inline Vec3 curl_from_component_gradient(int component,
                                         Real dphidx,
                                         Real dphidy,
                                         Real dphidz) noexcept {
    if (component == 0) {
        return Vec3{Real(0), dphidz, -dphidy};
    }
    if (component == 1) {
        return Vec3{-dphidz, Real(0), dphidx};
    }
    return Vec3{dphidy, -dphidx, Real(0)};
}

inline void axpy_qpoints(Real* target,
                         const Real* source,
                         Real coefficient,
                         std::size_t num_qpts) noexcept {
    for (std::size_t q = 0; q < num_qpts; ++q) {
        target[q] += coefficient * source[q];
    }
}

void write_vector_values_strided(const std::vector<Vec3>& values,
                                 std::size_t num_dofs,
                                 std::size_t output_stride,
                                 std::size_t q,
                                 Real* SVMP_RESTRICT values_out);
void write_vector_jacobians_strided(const std::vector<VectorJacobian>& jacobians,
                                    std::size_t num_dofs,
                                    std::size_t output_stride,
                                    std::size_t q,
                                    Real* SVMP_RESTRICT jacobians_out);
void write_vector_curl_strided(const std::vector<Vec3>& curl,
                               std::size_t num_dofs,
                               std::size_t output_stride,
                               std::size_t q,
                               Real* SVMP_RESTRICT curls_out);
void write_vector_divergence_strided(const std::vector<Real>& divergence,
                                     std::size_t num_dofs,
                                     std::size_t output_stride,
                                     std::size_t q,
                                     Real* SVMP_RESTRICT divergence_out);
void write_curl_and_divergence_from_jacobians_strided(
    const std::vector<VectorJacobian>& jacobians,
    std::size_t num_dofs,
    std::size_t output_stride,
    std::size_t q,
    Real* SVMP_RESTRICT curls_out,
    Real* SVMP_RESTRICT divergence_out);

template <typename BasisLike>
void evaluate_vector_public_api_strided(
    const BasisLike& basis,
    const std::vector<Vec3>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT jacobians_out,
    Real* SVMP_RESTRICT curls_out,
    Real* SVMP_RESTRICT divergence_out,
    bool use_direct_curl,
    bool use_direct_divergence,
    const char* family_name) {
    const std::size_t num_qpts = points.size();
    const std::size_t num_dofs = basis.size();
    validate_vector_strided_outputs(num_qpts, output_stride, family_name);

    auto& scratch = vector_basis_scratch();
    for (std::size_t q = 0; q < num_qpts; ++q) {
        if (values_out != nullptr) {
            basis.evaluate_vector_values(points[q], scratch.api_values);
            write_vector_values_strided(
                scratch.api_values, num_dofs, output_stride, q, values_out);
        }

        const bool needs_jacobians =
            jacobians_out != nullptr ||
            (curls_out != nullptr && !use_direct_curl) ||
            (divergence_out != nullptr && !use_direct_divergence);

        if (needs_jacobians) {
            basis.evaluate_vector_jacobians(points[q], scratch.api_jacobians);
            write_vector_jacobians_strided(
                scratch.api_jacobians, num_dofs, output_stride, q, jacobians_out);
            write_curl_and_divergence_from_jacobians_strided(
                scratch.api_jacobians,
                num_dofs,
                output_stride,
                q,
                curls_out,
                divergence_out);
            continue;
        }

        if (curls_out != nullptr) {
            basis.evaluate_curl(points[q], scratch.api_curl);
            write_vector_curl_strided(
                scratch.api_curl, num_dofs, output_stride, q, curls_out);
        }
        if (divergence_out != nullptr) {
            basis.evaluate_divergence(points[q], scratch.api_divergence);
            write_vector_divergence_strided(
                scratch.api_divergence, num_dofs, output_stride, q, divergence_out);
        }
    }
}

Vec3 lerp(const Vec3& a, const Vec3& b, Real s);
Vec3 bilinear(const std::array<Vec3, 4>& v, Real u, Real w);
Vec3 bilinear_du(const std::array<Vec3, 4>& v, Real u, Real w);
Vec3 bilinear_dw(const std::array<Vec3, 4>& v, Real u, Real w);
Vec3 cross3(const Vec3& a, const Vec3& b);
Vec3 normalize3(const Vec3& v);

template <typename ModalPolynomials>
std::array<int, 3> modal_power_limits(const ModalPolynomials& monomials) {
    std::array<int, 3> limits{{0, 0, 0}};
    for (const auto& poly : monomials) {
        for (int t = 0; t < poly.num_terms; ++t) {
            const auto& m = poly.terms[static_cast<std::size_t>(t)];
            limits[0] = std::max(limits[0], m.px);
            limits[1] = std::max(limits[1], m.py);
            limits[2] = std::max(limits[2], m.pz);
        }
    }
    return limits;
}

std::array<int, 3> component_monomial_power_limits(
    const std::vector<std::array<int, 4>>& candidates);
std::size_t triangle_poly_dim(std::size_t k);
std::size_t tetra_poly_dim(std::size_t k);
std::size_t rt_wedge_size(int order);
std::size_t rt_pyramid_size(int order);
std::size_t nd_wedge_size(int order);
std::size_t nd_pyramid_size(int order);
void ensure_supported_hybrid_vector_order(ElementType type,
                                          int order,
                                          const char* family_name);
std::vector<std::array<int, 4>> make_component_monomial_candidates(int max_total_degree);
std::vector<std::array<int, 4>> make_rt_extra_monomial_candidates(ElementType type,
                                                                  int order);
Real eval_transformed_rt_monomial_scalar(const std::array<int, 4>& mono,
                                         const std::vector<Real>& px,
                                         const std::vector<Real>& py,
                                         const std::vector<Real>& pz);
Real eval_transformed_rt_monomial_divergence(const std::array<int, 4>& mono,
                                             const std::vector<Real>& px,
                                             const std::vector<Real>& py,
                                             const std::vector<Real>& pz);

void add_component_monomial_jacobian(VectorJacobian& J,
                                     int component,
                                     int px_pow,
                                     int py_pow,
                                     int pz_pow,
                                     Real coefficient,
                                     const std::vector<Real>& px,
                                     const std::vector<Real>& py,
                                     const std::vector<Real>& pz);
VectorJacobian eval_transformed_component_monomial_jacobian(
    const std::array<int, 4>& mono,
    const std::vector<Real>& px,
    const std::vector<Real>& py,
    const std::vector<Real>& pz);
void add_component_monomial_curl(Vec3& curl,
                                 int component,
                                 int px_pow,
                                 int py_pow,
                                 int pz_pow,
                                 Real coefficient,
                                 const std::vector<Real>& px,
                                 const std::vector<Real>& py,
                                 const std::vector<Real>& pz);

template <typename ModalPolynomials>
void evaluate_nodal_modal_vector_values_with_limits(const ModalPolynomials& monomials,
                                                    const SparseModalCoefficientMatrix& sparse_coeffs,
                                                    std::size_t n,
                                                    const Vec3& xi,
                                                    const std::array<int, 3>& power_limits,
                                                    std::vector<Vec3>& values) {
    values.assign(n, Vec3{});

    auto& scratch = vector_basis_scratch();
    fill_power_tables(xi, power_limits, scratch);
    const auto& px = scratch.px;
    const auto& py = scratch.py;
    const auto& pz = scratch.pz;

    auto& modal_vals = scratch.vector_values;
    modal_vals.assign(n, Vec3{});
    for (std::size_t p = 0; p < n; ++p) {
        const auto& poly = monomials[p];
        auto& v = modal_vals[p];
        for (int t = 0; t < poly.num_terms; ++t) {
            const auto& m = poly.terms[static_cast<std::size_t>(t)];
            const Real mv =
                px[static_cast<std::size_t>(m.px)] *
                py[static_cast<std::size_t>(m.py)] *
                pz[static_cast<std::size_t>(m.pz)];
            v[static_cast<std::size_t>(m.component)] += m.coefficient * mv;
        }
    }

    BASIS_CHECK_CONSTRUCTION(sparse_coeffs.rows == n &&
                     sparse_coeffs.cols == n &&
                     sparse_coeffs.row_offsets.size() == n + 1u,
                 "evaluate_nodal_modal_vector_values: sparse coefficient size mismatch");
    BASIS_CHECK_CONSTRUCTION(sparse_coeffs.dofs.size() == sparse_coeffs.coefficients.size(),
                 "evaluate_nodal_modal_vector_values: sparse coefficient entry mismatch");
    for (std::size_t p = 0; p < n; ++p) {
        const Vec3& mv = modal_vals[p];
        const std::size_t row_begin = sparse_coeffs.row_offsets[p];
        const std::size_t row_end = sparse_coeffs.row_offsets[p + 1u];
        for (std::size_t entry = row_begin; entry < row_end; ++entry) {
            const std::size_t dof = sparse_coeffs.dofs[entry];
            const Real c = sparse_coeffs.coefficients[entry];
            values[dof][0] += c * mv[0];
            values[dof][1] += c * mv[1];
            values[dof][2] += c * mv[2];
        }
    }
}

template <typename ModalPolynomials>
void evaluate_nodal_modal_vector_jacobians_with_limits(const ModalPolynomials& monomials,
                                                       const SparseModalCoefficientMatrix& sparse_coeffs,
                                                       std::size_t n,
                                                       const Vec3& xi,
                                                       const std::array<int, 3>& power_limits,
                                                       std::vector<VectorJacobian>& jacobians) {
    jacobians.assign(n, VectorJacobian{});

    auto& scratch = vector_basis_scratch();
    fill_power_tables(xi, power_limits, scratch);
    const auto& px = scratch.px;
    const auto& py = scratch.py;
    const auto& pz = scratch.pz;

    auto& modal_jacs = scratch.vector_jacobians;
    modal_jacs.assign(n, VectorJacobian{});
    for (std::size_t p = 0; p < n; ++p) {
        const auto& poly = monomials[p];
        auto& J = modal_jacs[p];
        for (int t = 0; t < poly.num_terms; ++t) {
            const auto& m = poly.terms[static_cast<std::size_t>(t)];
            add_component_monomial_jacobian(J, m.component, m.px, m.py, m.pz, m.coefficient, px, py, pz);
        }
    }

    BASIS_CHECK_CONSTRUCTION(sparse_coeffs.rows == n &&
                     sparse_coeffs.cols == n &&
                     sparse_coeffs.row_offsets.size() == n + 1u,
                 "evaluate_nodal_modal_vector_jacobians: sparse coefficient size mismatch");
    BASIS_CHECK_CONSTRUCTION(sparse_coeffs.dofs.size() == sparse_coeffs.coefficients.size(),
                 "evaluate_nodal_modal_vector_jacobians: sparse coefficient entry mismatch");
    for (std::size_t p = 0; p < n; ++p) {
        const auto& Jp = modal_jacs[p];
        const std::size_t row_begin = sparse_coeffs.row_offsets[p];
        const std::size_t row_end = sparse_coeffs.row_offsets[p + 1u];
        for (std::size_t entry = row_begin; entry < row_end; ++entry) {
            const std::size_t dof = sparse_coeffs.dofs[entry];
            const Real c = sparse_coeffs.coefficients[entry];
            for (std::size_t r = 0; r < 3; ++r) {
                for (std::size_t col = 0; col < 3; ++col) {
                    jacobians[dof](r, col) += c * Jp(r, col);
                }
            }
        }
    }
}

template <typename ModalPolynomials>
void evaluate_nodal_modal_vector_curl_with_limits(const ModalPolynomials& monomials,
                                                  const SparseModalCoefficientMatrix& sparse_coeffs,
                                                  std::size_t n,
                                                  const Vec3& xi,
                                                  const std::array<int, 3>& power_limits,
                                                  std::vector<Vec3>& curl) {
    curl.assign(n, Vec3{});

    auto& scratch = vector_basis_scratch();
    fill_power_tables(xi, power_limits, scratch);
    const auto& px = scratch.px;
    const auto& py = scratch.py;
    const auto& pz = scratch.pz;

    auto& modal_curl = scratch.vector_values;
    modal_curl.assign(n, Vec3{});
    for (std::size_t p = 0; p < n; ++p) {
        const auto& poly = monomials[p];
        auto& c = modal_curl[p];
        for (int t = 0; t < poly.num_terms; ++t) {
            const auto& m = poly.terms[static_cast<std::size_t>(t)];
            add_component_monomial_curl(c, m.component, m.px, m.py, m.pz, m.coefficient, px, py, pz);
        }
    }

    BASIS_CHECK_CONSTRUCTION(sparse_coeffs.rows == n &&
                     sparse_coeffs.cols == n &&
                     sparse_coeffs.row_offsets.size() == n + 1u,
                 "evaluate_nodal_modal_vector_curl: sparse coefficient size mismatch");
    BASIS_CHECK_CONSTRUCTION(sparse_coeffs.dofs.size() == sparse_coeffs.coefficients.size(),
                 "evaluate_nodal_modal_vector_curl: sparse coefficient entry mismatch");
    for (std::size_t p = 0; p < n; ++p) {
        const Vec3& cm = modal_curl[p];
        const std::size_t row_begin = sparse_coeffs.row_offsets[p];
        const std::size_t row_end = sparse_coeffs.row_offsets[p + 1u];
        for (std::size_t entry = row_begin; entry < row_end; ++entry) {
            const std::size_t dof = sparse_coeffs.dofs[entry];
            const Real c = sparse_coeffs.coefficients[entry];
            curl[dof][0] += c * cm[0];
            curl[dof][1] += c * cm[1];
            curl[dof][2] += c * cm[2];
        }
    }
}

template <typename ModalPolynomials>
void evaluate_nodal_modal_divergence_with_limits(const ModalPolynomials& monomials,
                                                 const SparseModalCoefficientMatrix& sparse_coeffs,
                                                 std::size_t n,
                                                 const Vec3& xi,
                                                 const std::array<int, 3>& power_limits,
                                                 std::vector<Real>& divergence) {
    divergence.assign(n, Real(0));

    auto& scratch = vector_basis_scratch();
    fill_power_tables(xi, power_limits, scratch);
    const auto& px = scratch.px;
    const auto& py = scratch.py;
    const auto& pz = scratch.pz;

    auto& modal_divergence = scratch.scalars;
    modal_divergence.assign(n, Real(0));
    for (std::size_t p = 0; p < n; ++p) {
        const auto& poly = monomials[p];
        Real div = Real(0);
        for (int t = 0; t < poly.num_terms; ++t) {
            const auto& m = poly.terms[static_cast<std::size_t>(t)];
            if (m.component == 0 && m.px > 0) {
                div += m.coefficient * Real(m.px) *
                       px[static_cast<std::size_t>(m.px - 1)] *
                       py[static_cast<std::size_t>(m.py)] *
                       pz[static_cast<std::size_t>(m.pz)];
            } else if (m.component == 1 && m.py > 0) {
                div += m.coefficient * Real(m.py) *
                       px[static_cast<std::size_t>(m.px)] *
                       py[static_cast<std::size_t>(m.py - 1)] *
                       pz[static_cast<std::size_t>(m.pz)];
            } else if (m.component == 2 && m.pz > 0) {
                div += m.coefficient * Real(m.pz) *
                       px[static_cast<std::size_t>(m.px)] *
                       py[static_cast<std::size_t>(m.py)] *
                       pz[static_cast<std::size_t>(m.pz - 1)];
            }
        }
        modal_divergence[p] = div;
    }

    BASIS_CHECK_CONSTRUCTION(sparse_coeffs.rows == n &&
                     sparse_coeffs.cols == n &&
                     sparse_coeffs.row_offsets.size() == n + 1u,
                 "evaluate_nodal_modal_divergence: sparse coefficient size mismatch");
    BASIS_CHECK_CONSTRUCTION(sparse_coeffs.dofs.size() == sparse_coeffs.coefficients.size(),
                 "evaluate_nodal_modal_divergence: sparse coefficient entry mismatch");
    for (std::size_t p = 0; p < n; ++p) {
        const Real div = modal_divergence[p];
        if (div == Real(0)) {
            continue;
        }
        const std::size_t row_begin = sparse_coeffs.row_offsets[p];
        const std::size_t row_end = sparse_coeffs.row_offsets[p + 1u];
        for (std::size_t entry = row_begin; entry < row_end; ++entry) {
            divergence[sparse_coeffs.dofs[entry]] +=
                sparse_coeffs.coefficients[entry] * div;
        }
    }
}

template <typename ModalPolynomials>
void evaluate_nodal_modal_vector_strided_with_limits(
    const ModalPolynomials& monomials,
    const SparseModalCoefficientMatrix& sparse_coeffs,
    std::size_t n,
    const std::vector<Vec3>& points,
    std::size_t output_stride,
    const std::array<int, 3>& power_limits,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT jacobians_out,
    Real* SVMP_RESTRICT curls_out,
    Real* SVMP_RESTRICT divergence_out,
    const char* family_name) {
    const std::size_t num_qpts = points.size();
    validate_vector_strided_outputs(num_qpts, output_stride, family_name);
    BASIS_CHECK_CONSTRUCTION(sparse_coeffs.rows == n &&
                     sparse_coeffs.cols == n &&
                     sparse_coeffs.row_offsets.size() == n + 1u,
                 "evaluate_nodal_modal_vector_strided: sparse coefficient size mismatch");
    BASIS_CHECK_CONSTRUCTION(sparse_coeffs.dofs.size() == sparse_coeffs.coefficients.size(),
                 "evaluate_nodal_modal_vector_strided: sparse coefficient entry mismatch");

    auto& scratch = vector_basis_scratch();
    const bool need_values = values_out != nullptr;
    const bool need_jacobians = jacobians_out != nullptr;
    const bool need_curls = curls_out != nullptr;
    const bool need_divergence = divergence_out != nullptr;

    if (need_values) {
        zero_active_strided_rows(values_out, n * 3u, output_stride, num_qpts);
    }
    if (need_jacobians) {
        zero_active_strided_rows(jacobians_out, n * 9u, output_stride, num_qpts);
    }
    if (need_curls) {
        zero_active_strided_rows(curls_out, n * 3u, output_stride, num_qpts);
    }
    if (need_divergence) {
        zero_active_strided_rows(divergence_out, n, output_stride, num_qpts);
    }
    if (num_qpts == 0 || n == 0) {
        return;
    }

    fill_batched_power_tables(points, power_limits, scratch);
    const auto& px = scratch.batched_px;
    const auto& py = scratch.batched_py;
    const auto& pz = scratch.batched_pz;
    const std::size_t power_stride = num_qpts;
    const bool need_modal_gradient = need_jacobians || need_curls || need_divergence;

    auto& modal_values = scratch.modal_values_batched;
    auto& modal_jacobians = scratch.modal_jacobians_batched;
    auto& modal_curls = scratch.modal_curls_batched;
    auto& modal_divergence = scratch.modal_divergence_batched;

    for (std::size_t p = 0; p < n; ++p) {
        if (need_values) {
            modal_values.assign(3u * num_qpts, Real(0));
        }
        if (need_jacobians) {
            modal_jacobians.assign(9u * num_qpts, Real(0));
        }
        if (need_curls) {
            modal_curls.assign(3u * num_qpts, Real(0));
        }
        if (need_divergence) {
            modal_divergence.assign(num_qpts, Real(0));
        }

        const auto& poly = monomials[p];
        for (int term_index = 0; term_index < poly.num_terms; ++term_index) {
            const auto& term = poly.terms[static_cast<std::size_t>(term_index)];
            const std::size_t component = static_cast<std::size_t>(term.component);
            Real* modal_value_row = need_values
                ? modal_values.data() + component * num_qpts
                : nullptr;
            Real* modal_jacobian_row = need_jacobians
                ? modal_jacobians.data() + component * 3u * num_qpts
                : nullptr;
            Real* modal_curl_rows = need_curls ? modal_curls.data() : nullptr;
            Real* modal_divergence_row =
                need_divergence ? modal_divergence.data() : nullptr;

            if (need_values) {
                for (std::size_t q = 0; q < num_qpts; ++q) {
                    modal_value_row[q] +=
                        term.coefficient *
                        batched_power_product(px,
                                              py,
                                              pz,
                                              power_stride,
                                              term.px,
                                              term.py,
                                              term.pz,
                                              q);
                }
            }

            if (need_modal_gradient) {
                for (std::size_t q = 0; q < num_qpts; ++q) {
                    const Real dphidx =
                        term.coefficient *
                        batched_component_partial(px,
                                                  py,
                                                  pz,
                                                  power_stride,
                                                  term.px,
                                                  term.py,
                                                  term.pz,
                                                  0,
                                                  q);
                    const Real dphidy =
                        term.coefficient *
                        batched_component_partial(px,
                                                  py,
                                                  pz,
                                                  power_stride,
                                                  term.px,
                                                  term.py,
                                                  term.pz,
                                                  1,
                                                  q);
                    const Real dphidz =
                        term.coefficient *
                        batched_component_partial(px,
                                                  py,
                                                  pz,
                                                  power_stride,
                                                  term.px,
                                                  term.py,
                                                  term.pz,
                                                  2,
                                                  q);

                    if (need_jacobians) {
                        modal_jacobian_row[q] += dphidx;
                        modal_jacobian_row[num_qpts + q] += dphidy;
                        modal_jacobian_row[2u * num_qpts + q] += dphidz;
                    }
                    if (need_curls) {
                        const Vec3 curl =
                            curl_from_component_gradient(term.component,
                                                         dphidx,
                                                         dphidy,
                                                         dphidz);
                        modal_curl_rows[q] += curl[0];
                        modal_curl_rows[num_qpts + q] += curl[1];
                        modal_curl_rows[2u * num_qpts + q] += curl[2];
                    }
                    if (need_divergence) {
                        const Real div = term.component == 0 ? dphidx
                                       : term.component == 1 ? dphidy
                                                            : dphidz;
                        modal_divergence_row[q] += div;
                    }
                }
            }
        }

        const std::size_t row_begin = sparse_coeffs.row_offsets[p];
        const std::size_t row_end = sparse_coeffs.row_offsets[p + 1u];
        for (std::size_t entry = row_begin; entry < row_end; ++entry) {
            const std::size_t dof = sparse_coeffs.dofs[entry];
            const Real c = sparse_coeffs.coefficients[entry];
            if (need_values) {
                for (std::size_t component = 0; component < 3u; ++component) {
                    axpy_qpoints(values_out + (dof * 3u + component) * output_stride,
                                 modal_values.data() + component * num_qpts,
                                 c,
                                 num_qpts);
                }
            }
            if (need_jacobians) {
                for (std::size_t row = 0; row < 3u; ++row) {
                    for (std::size_t col = 0; col < 3u; ++col) {
                        axpy_qpoints(jacobians_out +
                                         (dof * 9u + row * 3u + col) * output_stride,
                                     modal_jacobians.data() +
                                         (row * 3u + col) * num_qpts,
                                     c,
                                     num_qpts);
                    }
                }
            }
            if (need_curls) {
                for (std::size_t component = 0; component < 3u; ++component) {
                    axpy_qpoints(curls_out + (dof * 3u + component) * output_stride,
                                 modal_curls.data() + component * num_qpts,
                                 c,
                                 num_qpts);
                }
            }
            if (need_divergence) {
                axpy_qpoints(divergence_out + dof * output_stride,
                             modal_divergence.data(),
                             c,
                             num_qpts);
            }
        }
    }
}

std::vector<std::array<int, 4>> make_nd_extra_monomial_candidates(ElementType type,
                                                                  int order);
Real eval_transformed_nd_monomial_scalar(const std::array<int, 4>& mono,
                                         const std::vector<Real>& px,
                                         const std::vector<Real>& py,
                                         const std::vector<Real>& pz);
Vec3 eval_transformed_nd_monomial_curl(const std::array<int, 4>& mono,
                                       const std::vector<Real>& px,
                                       const std::vector<Real>& py,
                                       const std::vector<Real>& pz);


} // namespace vector_common
} // namespace detail
} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_VECTORBASISEVALUATIONHELPERS_H
