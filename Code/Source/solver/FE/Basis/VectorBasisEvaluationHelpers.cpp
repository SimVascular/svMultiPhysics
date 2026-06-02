/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#include "VectorBasisEvaluationHelpers.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace svmp {
namespace FE {
namespace basis {
namespace detail {
namespace vector_common {

VectorBasisScratch& vector_basis_scratch() {
    // Scratch is intentionally thread-local: production assembly uses a
    // persistent worker-thread team, so buffers stay warm on each worker.
    static thread_local VectorBasisScratch scratch;
    return scratch;
}

void prewarm_vector_basis_scratch(std::size_t max_size, std::size_t max_qpts) {
    vector_basis_scratch().prewarm(max_size, max_qpts);
}

void fill_powers(Real x, int max_p, std::vector<Real>& out) {
    BASIS_CHECK_CONSTRUCTION(max_p >= 0, "powers: negative max_p");
    out.assign(static_cast<std::size_t>(max_p + 1), Real(1));
    for (int i = 1; i <= max_p; ++i) {
        out[static_cast<std::size_t>(i)] =
            out[static_cast<std::size_t>(i - 1)] * x;
    }
}

void fill_power_tables(const Vec3& xi,
                       const std::array<int, 3>& limits,
                       VectorBasisScratch& scratch) {
    fill_powers(xi[0], limits[0], scratch.px);
    fill_powers(xi[1], limits[1], scratch.py);
    fill_powers(xi[2], limits[2], scratch.pz);
}

namespace {

constexpr Real kSparseCoefficientRelativeTolerance =
    Real(256) * std::numeric_limits<Real>::epsilon();

void fill_batched_axis_powers(const std::vector<Vec3>& points,
                              std::size_t axis,
                              int max_power,
                              std::vector<Real>& out) {
    BASIS_CHECK_CONSTRUCTION(max_power >= 0, "batched powers: negative max_p");
    const std::size_t num_qpts = points.size();
    out.assign(static_cast<std::size_t>(max_power + 1) * num_qpts, Real(1));
    if (num_qpts == 0 || max_power == 0) {
        return;
    }

    Real* first_power = out.data() + num_qpts;
    for (std::size_t q = 0; q < num_qpts; ++q) {
        first_power[q] = points[q][axis];
    }
    for (int power = 2; power <= max_power; ++power) {
        const Real* previous =
            out.data() + static_cast<std::size_t>(power - 1) * num_qpts;
        Real* current = out.data() + static_cast<std::size_t>(power) * num_qpts;
        for (std::size_t q = 0; q < num_qpts; ++q) {
            current[q] = previous[q] * points[q][axis];
        }
    }
}

} // namespace

void fill_batched_power_tables(const std::vector<Vec3>& points,
                               const std::array<int, 3>& limits,
                               VectorBasisScratch& scratch) {
    fill_batched_axis_powers(points, 0u, limits[0], scratch.batched_px);
    fill_batched_axis_powers(points, 1u, limits[1], scratch.batched_py);
    fill_batched_axis_powers(points, 2u, limits[2], scratch.batched_pz);
}

void validate_vector_strided_outputs(std::size_t num_qpts,
                                     std::size_t output_stride,
                                     const char* family_name) {
    if (output_stride < num_qpts) {
        throw BasisConfigurationException(
            std::string(family_name) +
                " strided vector evaluation requires output_stride >= points.size()",
            __FILE__, __LINE__, __func__);
    }
}

void zero_active_strided_rows(Real* output,
                              std::size_t rows,
                              std::size_t output_stride,
                              std::size_t num_qpts) {
    for (std::size_t row = 0; row < rows; ++row) {
        std::fill_n(output + row * output_stride, num_qpts, Real(0));
    }
}

SparseModalCoefficientMatrix build_sparse_modal_coefficients(
    const std::vector<Real>& dense_coefficients,
    std::size_t rows,
    std::size_t cols) {
    BASIS_CHECK_CONSTRUCTION(dense_coefficients.size() == rows * cols,
                 "build_sparse_modal_coefficients: dense coefficient size mismatch");

    SparseModalCoefficientMatrix sparse;
    sparse.rows = rows;
    sparse.cols = cols;
    sparse.row_offsets.reserve(rows + 1u);
    sparse.row_offsets.push_back(0u);

    Real max_abs = Real(0);
    for (const Real coefficient : dense_coefficients) {
        max_abs = std::max(max_abs, std::abs(coefficient));
    }
    const Real prune_threshold = kSparseCoefficientRelativeTolerance * max_abs;

    for (std::size_t row = 0; row < rows; ++row) {
        const Real* dense_row = dense_coefficients.data() + row * cols;
        for (std::size_t col = 0; col < cols; ++col) {
            const Real coefficient = dense_row[col];
            if (std::abs(coefficient) > prune_threshold) {
                sparse.dofs.push_back(col);
                sparse.coefficients.push_back(coefficient);
            }
        }
        sparse.row_offsets.push_back(sparse.dofs.size());
    }

    return sparse;
}

Vec3 curl_from_jacobian(const VectorJacobian& J) noexcept {
    return Vec3{J(2u, 1u) - J(1u, 2u),
                J(0u, 2u) - J(2u, 0u),
                J(1u, 0u) - J(0u, 1u)};
}

Real divergence_from_jacobian(const VectorJacobian& J) noexcept {
    return J(0u, 0u) + J(1u, 1u) + J(2u, 2u);
}

void write_vector_values_strided(const std::vector<Vec3>& values,
                                 std::size_t num_dofs,
                                 std::size_t output_stride,
                                 std::size_t q,
                                 Real* SVMP_RESTRICT values_out) {
    if (values_out == nullptr) {
        return;
    }
    BASIS_CHECK_CONSTRUCTION(values.size() == num_dofs,
                 "vector value evaluation returned the wrong number of DOFs");
    for (std::size_t dof = 0; dof < num_dofs; ++dof) {
        for (std::size_t component = 0; component < 3u; ++component) {
            values_out[(dof * 3u + component) * output_stride + q] =
                values[dof][component];
        }
    }
}

void write_vector_jacobians_strided(const std::vector<VectorJacobian>& jacobians,
                                    std::size_t num_dofs,
                                    std::size_t output_stride,
                                    std::size_t q,
                                    Real* SVMP_RESTRICT jacobians_out) {
    if (jacobians_out == nullptr) {
        return;
    }
    BASIS_CHECK_CONSTRUCTION(jacobians.size() == num_dofs,
                 "vector Jacobian evaluation returned the wrong number of DOFs");
    for (std::size_t dof = 0; dof < num_dofs; ++dof) {
        const auto& J = jacobians[dof];
        for (std::size_t component = 0; component < 3u; ++component) {
            for (std::size_t derivative = 0; derivative < 3u; ++derivative) {
                jacobians_out[(dof * 9u + component * 3u + derivative) *
                                  output_stride + q] = J(component, derivative);
            }
        }
    }
}

void write_vector_curl_strided(const std::vector<Vec3>& curl,
                               std::size_t num_dofs,
                               std::size_t output_stride,
                               std::size_t q,
                               Real* SVMP_RESTRICT curls_out) {
    if (curls_out == nullptr) {
        return;
    }
    BASIS_CHECK_CONSTRUCTION(curl.size() == num_dofs,
                 "vector curl evaluation returned the wrong number of DOFs");
    for (std::size_t dof = 0; dof < num_dofs; ++dof) {
        for (std::size_t component = 0; component < 3u; ++component) {
            curls_out[(dof * 3u + component) * output_stride + q] =
                curl[dof][component];
        }
    }
}

void write_vector_divergence_strided(const std::vector<Real>& divergence,
                                     std::size_t num_dofs,
                                     std::size_t output_stride,
                                     std::size_t q,
                                     Real* SVMP_RESTRICT divergence_out) {
    if (divergence_out == nullptr) {
        return;
    }
    BASIS_CHECK_CONSTRUCTION(divergence.size() == num_dofs,
                 "vector divergence evaluation returned the wrong number of DOFs");
    for (std::size_t dof = 0; dof < num_dofs; ++dof) {
        divergence_out[dof * output_stride + q] = divergence[dof];
    }
}

void write_curl_and_divergence_from_jacobians_strided(
    const std::vector<VectorJacobian>& jacobians,
    std::size_t num_dofs,
    std::size_t output_stride,
    std::size_t q,
    Real* SVMP_RESTRICT curls_out,
    Real* SVMP_RESTRICT divergence_out) {
    BASIS_CHECK_CONSTRUCTION(jacobians.size() == num_dofs,
                 "vector Jacobian evaluation returned the wrong number of DOFs");
    for (std::size_t dof = 0; dof < num_dofs; ++dof) {
        const auto& J = jacobians[dof];
        if (curls_out != nullptr) {
            const Vec3 curl = curl_from_jacobian(J);
            for (std::size_t component = 0; component < 3u; ++component) {
                curls_out[(dof * 3u + component) * output_stride + q] =
                    curl[component];
            }
        }
        if (divergence_out != nullptr) {
            divergence_out[dof * output_stride + q] = divergence_from_jacobian(J);
        }
    }
}

Vec3 lerp(const Vec3& a, const Vec3& b, Real s) {
    const Real t = (s + Real(1)) * Real(0.5);
    return a * (Real(1) - t) + b * t;
}

Vec3 bilinear(const std::array<Vec3, 4>& v, Real u, Real w) {
    const Real N0 = Real(0.25) * (Real(1) - u) * (Real(1) - w);
    const Real N1 = Real(0.25) * (Real(1) + u) * (Real(1) - w);
    const Real N2 = Real(0.25) * (Real(1) + u) * (Real(1) + w);
    const Real N3 = Real(0.25) * (Real(1) - u) * (Real(1) + w);
    return v[0] * N0 + v[1] * N1 + v[2] * N2 + v[3] * N3;
}

Vec3 bilinear_du(const std::array<Vec3, 4>& v, Real u, Real w) {
    (void)u;
    const Real dN0 = -Real(0.25) * (Real(1) - w);
    const Real dN1 =  Real(0.25) * (Real(1) - w);
    const Real dN2 =  Real(0.25) * (Real(1) + w);
    const Real dN3 = -Real(0.25) * (Real(1) + w);
    return v[0] * dN0 + v[1] * dN1 + v[2] * dN2 + v[3] * dN3;
}

Vec3 bilinear_dw(const std::array<Vec3, 4>& v, Real u, Real w) {
    (void)w;
    const Real dN0 = -Real(0.25) * (Real(1) - u);
    const Real dN1 = -Real(0.25) * (Real(1) + u);
    const Real dN2 =  Real(0.25) * (Real(1) + u);
    const Real dN3 =  Real(0.25) * (Real(1) - u);
    return v[0] * dN0 + v[1] * dN1 + v[2] * dN2 + v[3] * dN3;
}

Vec3 cross3(const Vec3& a, const Vec3& b) {
    return Vec3{a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0]};
}

Vec3 normalize3(const Vec3& v) {
    const Real n = v.norm();
    BASIS_CHECK_CONSTRUCTION(n > std::numeric_limits<Real>::epsilon(),
                 "normalize3: zero-length vector");
    return v / n;
}

std::array<int, 3> component_monomial_power_limits(
    const std::vector<std::array<int, 4>>& candidates) {
    std::array<int, 3> limits{{0, 0, 0}};
    for (const auto& mono : candidates) {
        limits[0] = std::max(limits[0], mono[1]);
        limits[1] = std::max(limits[1], mono[2]);
        limits[2] = std::max(limits[2], mono[3]);
    }
    return limits;
}

std::size_t triangle_poly_dim(std::size_t k) {
    return (k + 1u) * (k + 2u) / 2u;
}

std::size_t tetra_poly_dim(std::size_t k) {
    return (k + 1u) * (k + 2u) * (k + 3u) / 6u;
}

std::size_t rt_wedge_size(int order) {
    const std::size_t k = static_cast<std::size_t>(order);
    const std::size_t face_dofs =
        2u * triangle_poly_dim(k) + 3u * (k + 1u) * (k + 1u);
    const std::size_t interior_dofs =
        (k >= 1u) ? (3u * k * (k + 1u) * (k + 1u) / 2u) : 0u;
    return face_dofs + interior_dofs;
}

std::size_t rt_pyramid_size(int order) {
    const std::size_t k = static_cast<std::size_t>(order);
    const std::size_t face_dofs = (k + 1u) * (k + 1u) + 4u * triangle_poly_dim(k);
    const std::size_t interior_dofs = (k >= 1u) ? (3u * k * k * k) : 0u;
    return face_dofs + interior_dofs;
}

std::size_t nd_wedge_size(int order) {
    const std::size_t k = static_cast<std::size_t>(order);
    const std::size_t edge_dofs = 9u * (k + 1u);
    const std::size_t face_dofs = (k >= 1u) ? (8u * k * (k + 1u)) : 0u;
    const std::size_t interior_dofs =
        (k >= 2u) ? (3u * k * (k - 1u) * (k + 1u) / 2u) : 0u;
    return edge_dofs + face_dofs + interior_dofs;
}

std::size_t nd_pyramid_size(int order) {
    const std::size_t k = static_cast<std::size_t>(order);
    const std::size_t edge_dofs = 8u * (k + 1u);
    const std::size_t face_dofs = (k >= 1u) ? (6u * k * (k + 1u)) : 0u;
    const std::size_t interior_dofs =
        (k >= 2u) ? (k * (k - 1u) * (k + 1u) / 2u) : 0u;
    return edge_dofs + face_dofs + interior_dofs;
}

void ensure_supported_hybrid_vector_order(ElementType type,
                                          int order,
                                          const char* family_name) {
    (void)type;
    (void)order;
    (void)family_name;
}

std::vector<std::array<int, 4>> make_component_monomial_candidates(
    int max_total_degree) {
    BASIS_CHECK_CONSTRUCTION(max_total_degree >= 0,
                 "make_component_monomial_candidates: negative total degree");

    std::vector<std::array<int, 4>> candidates;
    for (int component = 0; component < 3; ++component) {
        for (int total = 0; total <= max_total_degree; ++total) {
            for (int pz = 0; pz <= total; ++pz) {
                for (int py = 0; py <= total - pz; ++py) {
                    const int px = total - py - pz;
                    candidates.push_back({component, px, py, pz});
                }
            }
        }
    }
    return candidates;
}

std::vector<std::array<int, 4>> make_rt_extra_monomial_candidates(ElementType type,
                                                                  int order) {
    if (order >= 3) {
        return make_component_monomial_candidates(3 * order);
    }

    std::vector<std::array<int, 4>> candidates;
    if (!is_pyramid(type) || order != 2) {
        return candidates;
    }

    for (int component = 0; component < 3; ++component) {
        for (int pz = 0; pz <= 2; ++pz) {
            for (int py = 0; py <= 2 - pz; ++py) {
                for (int px = 0; px <= 2 - py - pz; ++px) {
                    candidates.push_back({component, px, py, pz});
                }
            }
        }
    }
    return candidates;
}

Real eval_transformed_rt_monomial_scalar(const std::array<int, 4>& mono,
                                         const std::vector<Real>& px,
                                         const std::vector<Real>& py,
                                         const std::vector<Real>& pz) {
    return px[static_cast<std::size_t>(mono[1])] *
           py[static_cast<std::size_t>(mono[2])] *
           pz[static_cast<std::size_t>(mono[3])];
}

Real eval_transformed_rt_monomial_divergence(const std::array<int, 4>& mono,
                                             const std::vector<Real>& px,
                                             const std::vector<Real>& py,
                                             const std::vector<Real>& pz) {
    const int component = mono[0];
    const int px_pow = mono[1];
    const int py_pow = mono[2];
    const int pz_pow = mono[3];

    if (component == 0) {
        if (px_pow == 0) {
            return Real(0);
        }
        return Real(px_pow) *
               px[static_cast<std::size_t>(px_pow - 1)] *
               py[static_cast<std::size_t>(py_pow)] *
               pz[static_cast<std::size_t>(pz_pow)];
    }
    if (component == 1) {
        if (py_pow == 0) {
            return Real(0);
        }
        return Real(py_pow) *
               px[static_cast<std::size_t>(px_pow)] *
               py[static_cast<std::size_t>(py_pow - 1)] *
               pz[static_cast<std::size_t>(pz_pow)];
    }
    if (pz_pow == 0) {
        return Real(0);
    }
    return Real(pz_pow) *
           px[static_cast<std::size_t>(px_pow)] *
           py[static_cast<std::size_t>(py_pow)] *
           pz[static_cast<std::size_t>(pz_pow - 1)];
}

void add_component_monomial_jacobian(VectorJacobian& J,
                                     int component,
                                     int px_pow,
                                     int py_pow,
                                     int pz_pow,
                                     Real coefficient,
                                     const std::vector<Real>& px,
                                     const std::vector<Real>& py,
                                     const std::vector<Real>& pz) {
    const auto comp = static_cast<std::size_t>(component);
    if (px_pow > 0) {
        J(comp, 0) += coefficient * Real(px_pow) *
                      px[static_cast<std::size_t>(px_pow - 1)] *
                      py[static_cast<std::size_t>(py_pow)] *
                      pz[static_cast<std::size_t>(pz_pow)];
    }
    if (py_pow > 0) {
        J(comp, 1) += coefficient * Real(py_pow) *
                      px[static_cast<std::size_t>(px_pow)] *
                      py[static_cast<std::size_t>(py_pow - 1)] *
                      pz[static_cast<std::size_t>(pz_pow)];
    }
    if (pz_pow > 0) {
        J(comp, 2) += coefficient * Real(pz_pow) *
                      px[static_cast<std::size_t>(px_pow)] *
                      py[static_cast<std::size_t>(py_pow)] *
                      pz[static_cast<std::size_t>(pz_pow - 1)];
    }
}

VectorJacobian eval_transformed_component_monomial_jacobian(
    const std::array<int, 4>& mono,
    const std::vector<Real>& px,
    const std::vector<Real>& py,
    const std::vector<Real>& pz) {
    VectorJacobian J{};
    add_component_monomial_jacobian(
        J, mono[0], mono[1], mono[2], mono[3], Real(1), px, py, pz);
    return J;
}

void add_component_monomial_curl(Vec3& curl,
                                 int component,
                                 int px_pow,
                                 int py_pow,
                                 int pz_pow,
                                 Real coefficient,
                                 const std::vector<Real>& px,
                                 const std::vector<Real>& py,
                                 const std::vector<Real>& pz) {
    const Real dphidx = (px_pow == 0)
        ? Real(0)
        : coefficient * Real(px_pow) *
              px[static_cast<std::size_t>(px_pow - 1)] *
              py[static_cast<std::size_t>(py_pow)] *
              pz[static_cast<std::size_t>(pz_pow)];
    const Real dphidy = (py_pow == 0)
        ? Real(0)
        : coefficient * Real(py_pow) *
              px[static_cast<std::size_t>(px_pow)] *
              py[static_cast<std::size_t>(py_pow - 1)] *
              pz[static_cast<std::size_t>(pz_pow)];
    const Real dphidz = (pz_pow == 0)
        ? Real(0)
        : coefficient * Real(pz_pow) *
              px[static_cast<std::size_t>(px_pow)] *
              py[static_cast<std::size_t>(py_pow)] *
              pz[static_cast<std::size_t>(pz_pow - 1)];

    if (component == 0) {
        curl[1] += dphidz;
        curl[2] -= dphidy;
    } else if (component == 1) {
        curl[0] -= dphidz;
        curl[2] += dphidx;
    } else {
        curl[0] += dphidy;
        curl[1] -= dphidx;
    }
}

std::vector<std::array<int, 4>> make_nd_extra_monomial_candidates(ElementType,
                                                                  int order) {
    if (order >= 3) {
        return make_component_monomial_candidates(3 * order);
    }

    std::vector<std::array<int, 4>> candidates;
    const int max_total_degree = (order == 1) ? 4 : 5;
    for (int component = 0; component < 3; ++component) {
        for (int total = 0; total <= max_total_degree; ++total) {
            for (int pz = 0; pz <= total; ++pz) {
                for (int py = 0; py <= total - pz; ++py) {
                    const int px = total - py - pz;
                    candidates.push_back({component, px, py, pz});
                }
            }
        }
    }
    return candidates;
}

Real eval_transformed_nd_monomial_scalar(const std::array<int, 4>& mono,
                                         const std::vector<Real>& px,
                                         const std::vector<Real>& py,
                                         const std::vector<Real>& pz) {
    return px[static_cast<std::size_t>(mono[1])] *
           py[static_cast<std::size_t>(mono[2])] *
           pz[static_cast<std::size_t>(mono[3])];
}

Vec3 eval_transformed_nd_monomial_curl(const std::array<int, 4>& mono,
                                       const std::vector<Real>& px,
                                       const std::vector<Real>& py,
                                       const std::vector<Real>& pz) {
    const int component = mono[0];
    const int px_pow = mono[1];
    const int py_pow = mono[2];
    const int pz_pow = mono[3];

    const Real dphidx = (px_pow == 0)
        ? Real(0)
        : Real(px_pow) *
              px[static_cast<std::size_t>(px_pow - 1)] *
              py[static_cast<std::size_t>(py_pow)] *
              pz[static_cast<std::size_t>(pz_pow)];
    const Real dphidy = (py_pow == 0)
        ? Real(0)
        : Real(py_pow) *
              px[static_cast<std::size_t>(px_pow)] *
              py[static_cast<std::size_t>(py_pow - 1)] *
              pz[static_cast<std::size_t>(pz_pow)];
    const Real dphidz = (pz_pow == 0)
        ? Real(0)
        : Real(pz_pow) *
              px[static_cast<std::size_t>(px_pow)] *
              py[static_cast<std::size_t>(py_pow)] *
              pz[static_cast<std::size_t>(pz_pow - 1)];

    if (component == 0) {
        return Vec3{Real(0), dphidz, -dphidy};
    }
    if (component == 1) {
        return Vec3{-dphidz, Real(0), dphidx};
    }
    return Vec3{dphidy, -dphidx, Real(0)};
}

} // namespace vector_common
} // namespace detail
} // namespace basis
} // namespace FE
} // namespace svmp
