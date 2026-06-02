/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#ifndef SVMP_FE_BASIS_BASISFUNCTION_H
#define SVMP_FE_BASIS_BASISFUNCTION_H

/**
 * @file BasisFunction.h
 * @brief Abstract interface for basis function evaluation on reference elements
 *
 * The Basis module operates purely on reference elements and is independent of
 * mesh-specific data structures. Implementations may leverage Math and
 * Quadrature utilities but must not read mesh connectivity or geometry.
 */

#include "Types.h"
#include "BasisExceptions.h"
#include "Math/Vector.h"
#include "Math/Matrix.h"
#include <algorithm>
#include <cstdint>
#include <functional>
#include <span>
#include <string>
#include <vector>

namespace svmp {
namespace FE {
namespace basis {

using Gradient = math::Vector<Real, 3>;
using Hessian  = math::Matrix<Real, 3, 3>;
using VectorJacobian = math::Matrix<Real, 3, 3>;

struct BasisIdentityFingerprint {
    std::uint64_t hash_a{0};
    std::uint64_t hash_b{0};
};

[[nodiscard]] BasisIdentityFingerprint
compute_basis_identity_fingerprint(std::span<const std::uint64_t> words) noexcept;

void prewarm_basis_function_scratch(std::size_t max_size,
                                    std::size_t max_qpts = 0);

[[nodiscard]] inline Hessian make_symmetric_hessian(Real xx,
                                                    Real yy,
                                                    Real zz,
                                                    Real xy,
                                                    Real xz,
                                                    Real yz) {
    Hessian hessian{};
    hessian(0, 0) = xx;
    hessian(1, 1) = yy;
    hessian(2, 2) = zz;
    hessian(0, 1) = xy;
    hessian(1, 0) = xy;
    hessian(0, 2) = xz;
    hessian(2, 0) = xz;
    hessian(1, 2) = yz;
    hessian(2, 1) = yz;
    return hessian;
}

// Raw Hessian buffers use row-major 3x3 blocks:
// dst[row * 3 + col] = H(row, col).
inline void store_hessian(const Hessian& hessian, Real* dst) noexcept {
    dst[0u] = hessian(0u, 0u);
    dst[1u] = hessian(0u, 1u);
    dst[2u] = hessian(0u, 2u);
    dst[3u] = hessian(1u, 0u);
    dst[4u] = hessian(1u, 1u);
    dst[5u] = hessian(1u, 2u);
    dst[6u] = hessian(2u, 0u);
    dst[7u] = hessian(2u, 1u);
    dst[8u] = hessian(2u, 2u);
}

inline void store_hessian_strided(const Hessian& hessian,
                                  Real* dst,
                                  std::size_t stride,
                                  std::size_t offset) noexcept {
    dst[0u * stride + offset] = hessian(0u, 0u);
    dst[1u * stride + offset] = hessian(0u, 1u);
    dst[2u * stride + offset] = hessian(0u, 2u);
    dst[3u * stride + offset] = hessian(1u, 0u);
    dst[4u * stride + offset] = hessian(1u, 1u);
    dst[5u * stride + offset] = hessian(1u, 2u);
    dst[6u * stride + offset] = hessian(2u, 0u);
    dst[7u * stride + offset] = hessian(2u, 1u);
    dst[8u * stride + offset] = hessian(2u, 2u);
}

inline void scatter_hessian_components_strided(const Real* src,
                                               Real* dst,
                                               std::size_t stride,
                                               std::size_t offset) noexcept {
    dst[0u * stride + offset] = src[0u];
    dst[1u * stride + offset] = src[1u];
    dst[2u * stride + offset] = src[2u];
    dst[3u * stride + offset] = src[3u];
    dst[4u * stride + offset] = src[4u];
    dst[5u * stride + offset] = src[5u];
    dst[6u * stride + offset] = src[6u];
    dst[7u * stride + offset] = src[7u];
    dst[8u * stride + offset] = src[8u];
}

[[nodiscard]] inline Hessian load_hessian(const Real* src) noexcept {
    Hessian hessian{};
    hessian(0u, 0u) = src[0u];
    hessian(0u, 1u) = src[1u];
    hessian(0u, 2u) = src[2u];
    hessian(1u, 0u) = src[3u];
    hessian(1u, 1u) = src[4u];
    hessian(1u, 2u) = src[5u];
    hessian(2u, 0u) = src[6u];
    hessian(2u, 1u) = src[7u];
    hessian(2u, 2u) = src[8u];
    return hessian;
}

inline void add_scaled_hessian(Hessian& target,
                               const Hessian& source,
                               Real scale) noexcept {
    target(0u, 0u) += scale * source(0u, 0u);
    target(0u, 1u) += scale * source(0u, 1u);
    target(0u, 2u) += scale * source(0u, 2u);
    target(1u, 0u) += scale * source(1u, 0u);
    target(1u, 1u) += scale * source(1u, 1u);
    target(1u, 2u) += scale * source(1u, 2u);
    target(2u, 0u) += scale * source(2u, 0u);
    target(2u, 1u) += scale * source(2u, 1u);
    target(2u, 2u) += scale * source(2u, 2u);
}

/**
 * @brief Base interface for scalar and vector-valued basis families
 *
 * All basis implementations operate in reference space. Physical mappings are
 * handled by the Geometry module. Derivatives are returned with unused
 * components set to zero for lower dimensional elements.
 */
class BasisFunction {
public:
    virtual ~BasisFunction() = default;

    /// Basis family identifier
    virtual BasisType basis_type() const noexcept = 0;

    /// Underlying element type on the reference domain
    virtual ElementType element_type() const noexcept = 0;

    /// Reference dimensionality (1, 2, or 3)
    virtual int dimension() const noexcept = 0;

    /// Polynomial order (modal/nodal definition dependent)
    virtual int order() const noexcept = 0;

    /// Number of basis functions (scalar or vector-valued)
    virtual std::size_t size() const noexcept = 0;

    /**
     * @brief Whether BasisCache can key this basis from common structural fields.
     *
     * Return true only when basis_type/element_type/dimension/order/size and
     * vector-valued status fully determine evaluation behavior. Parameterized
     * bases such as splines and custom user bases should keep the default false
     * so BasisCache includes cache_identity() in the key.
     */
    virtual bool cache_identity_is_structural() const noexcept { return false; }

    /// Whether the basis is vector-valued (H(div)/H(curl))
    virtual bool is_vector_valued() const noexcept { return false; }

    /// Whether vector-valued basis Jacobians are available.
    virtual bool supports_vector_jacobians() const noexcept { return false; }

    /// Whether vector-valued basis curls are available.
    virtual bool supports_curl() const noexcept { return false; }

    /// Whether vector-valued basis divergences are available.
    virtual bool supports_divergence() const noexcept { return false; }

    /**
     * @brief Stable semantic identity used by BasisCache
     *
     * Derived classes should override this when evaluation depends on
     * additional state beyond basis family / element / order metadata.
     */
    virtual std::string cache_identity() const;

    /**
     * @brief Optional exact structured identity payload for BasisCache keys.
     *
     * Parameterized bases may append stable integer/bit-pattern words and
     * return true to let BasisCache avoid using cache_identity() as the exact
     * key payload. The human-readable cache_identity() remains available for
     * diagnostics and for custom bases that do not implement this path.
     */
    virtual bool cache_identity_words(std::vector<std::uint64_t>& words) const;

    /**
     * @brief Optional cached fingerprint for structured identity words.
     *
     * Implementations that precompute cache_identity_words() may also cache the
     * corresponding fingerprint. BasisCache still retains exact identity words
     * for equality after hash matches.
     */
    virtual bool cache_identity_fingerprint(std::uint64_t& hash_a,
                                            std::uint64_t& hash_b) const;

    /**
     * @brief Evaluate scalar basis values at a reference point
     * @param xi Reference coordinates (unused entries are ignored)
     * @param[out] values Output array resized to size()
     */
    virtual void evaluate_values(const math::Vector<Real, 3>& xi,
                                 std::vector<Real>& values) const = 0;

    /**
     * @brief Evaluate gradients of scalar basis functions
     *
     * Production bases must override this with analytic derivatives.
     * Use numerical_gradient explicitly in tests or diagnostics when a finite
     * difference approximation is intended.
     */
    virtual void evaluate_gradients(const math::Vector<Real, 3>& xi,
                                    std::vector<Gradient>& gradients) const;

    /**
     * @brief Evaluate Hessians of scalar basis functions
     *
     * Production bases must override this with analytic second derivatives.
     * Use numerical_hessian explicitly in tests or diagnostics when a finite
     * difference approximation is intended.
     */
    virtual void evaluate_hessians(const math::Vector<Real, 3>& xi,
                                   std::vector<Hessian>& hessians) const;

    /**
     * @brief Fused evaluation of values, gradients, and Hessians at one point
     *
     * Default implementation calls evaluate_values, evaluate_gradients, and
     * evaluate_hessians in sequence. Bases that share intermediate
     * computations (e.g., LagrangeBasis sharing per-axis 1D evaluations)
     * should override this to avoid redundant work.
     */
    virtual void evaluate_all(const math::Vector<Real, 3>& xi,
                              std::vector<Real>& values,
                              std::vector<Gradient>& gradients,
                              std::vector<Hessian>& hessians) const;

    /**
     * @brief Fill SoA buffers with basis evaluations at all quadrature points
     *
     * Outputs are written directly to caller-provided strided buffers in
     * DOF-major SoA layout — no scratch+transpose required by the caller.
     * Pass `nullptr` for any output that is not needed.
     *
     *   values_out:    size num_dofs * num_qpts; element [d * num_qpts + q]
     *   gradients_out: size num_dofs * 3 * num_qpts; element [(d*3 + c) * num_qpts + q]
     *   hessians_out:  size num_dofs * 9 * num_qpts; element [(d*9 + r*3 + c) * num_qpts + q]
     *
     * Non-null output ranges must not overlap each other. Implementations may
     * fill requested quantities in any order that is efficient for the basis.
     *
     * Default implementation calls evaluate_all (or evaluate_values/gradients/
     * hessians as appropriate) per QP, materializing into temp buffers then
     * scatter-writing to the output. Performance-sensitive bases must override
     * this path so batched assembly does not fall back to Q virtual point
     * evaluations. Unit coverage keeps an explicit list of hot bases that are
     * expected to provide a direct strided implementation.
     */
    virtual void evaluate_at_quadrature_points(
        const std::vector<math::Vector<Real, 3>>& points,
        Real* SVMP_RESTRICT values_out,
        Real* SVMP_RESTRICT gradients_out,
        Real* SVMP_RESTRICT hessians_out) const;

    /**
     * @brief Fill strided SoA buffers with basis evaluations at quadrature points
     *
     * Same component layout as evaluate_at_quadrature_points, but each
     * dof/component row advances by `output_stride` rather than `points.size()`.
     * This lets padded SIMD cache storage be filled directly. Non-null output
     * ranges have the same non-overlap requirement.
     */
    virtual void evaluate_at_quadrature_points_strided(
        const std::vector<math::Vector<Real, 3>>& points,
        std::size_t output_stride,
        Real* SVMP_RESTRICT values_out,
        Real* SVMP_RESTRICT gradients_out,
        Real* SVMP_RESTRICT hessians_out) const;

    /**
     * @brief Fill zero-initialized scalar cache storage.
     *
     * BasisCache allocates and zero-initializes its scalar SoA buffers before
     * calling this hook. The default implementation overwrites all requested
     * entries through the public strided evaluator. Sparse-support bases may
     * override this and write only active entries, relying on the caller's
     * zero-initialization for inactive DOFs and unused derivative components.
     */
    virtual void fill_scalar_cache_entry(
        const std::vector<math::Vector<Real, 3>>& points,
        std::size_t output_stride,
        Real* SVMP_RESTRICT values_out,
        Real* SVMP_RESTRICT gradients_out,
        Real* SVMP_RESTRICT hessians_out) const;

    /**
     * @brief Fill SoA buffers with vector-basis evaluations at all quadrature points
     *
     * Outputs are written in DOF-major SoA layout. Pass `nullptr` for any
     * quantity that is not needed.
     *
     *   values_out:     size num_dofs * 3 * num_qpts; element [(d*3 + c) * num_qpts + q]
     *   jacobians_out:  size num_dofs * 9 * num_qpts; element [(d*9 + c*3 + r) * num_qpts + q]
     *   curls_out:      size num_dofs * 3 * num_qpts; element [(d*3 + c) * num_qpts + q]
     *   divergence_out: size num_dofs * num_qpts; element [d * num_qpts + q]
     *
     * Non-null output ranges must not overlap each other. Implementations may
     * fill requested quantities in any order that is efficient for the basis.
     */
    virtual void evaluate_vector_at_quadrature_points(
        const std::vector<math::Vector<Real, 3>>& points,
        Real* SVMP_RESTRICT values_out,
        Real* SVMP_RESTRICT jacobians_out,
        Real* SVMP_RESTRICT curls_out,
        Real* SVMP_RESTRICT divergence_out) const;

    /**
     * @brief Fill strided SoA buffers with vector-basis evaluations
     *
     * Same component layout as evaluate_vector_at_quadrature_points, but each
     * dof/component row advances by `output_stride` rather than `points.size()`.
     * Non-null output ranges have the same non-overlap requirement.
     *
     * The base fallback loops over quadrature points through virtual point
     * evaluation. H(div)/H(curl) bases used in assembly should override this
     * method directly, and tests track the current hot vector families.
     */
    virtual void evaluate_vector_at_quadrature_points_strided(
        const std::vector<math::Vector<Real, 3>>& points,
        std::size_t output_stride,
        Real* SVMP_RESTRICT values_out,
        Real* SVMP_RESTRICT jacobians_out,
        Real* SVMP_RESTRICT curls_out,
        Real* SVMP_RESTRICT divergence_out) const;

    /**
     * @brief Evaluate scalar basis values into a caller-provided raw buffer
     *
     * Caller is responsible for providing a buffer of at least size() Real
     * entries. This avoids the per-call std::vector::resize() cost of the
     * vector-output overload. Default implementation forwards through a temp
     * vector; bases should override for direct write.
     */
    virtual void evaluate_values_to(const math::Vector<Real, 3>& xi,
                                    Real* SVMP_RESTRICT values_out) const;

    /**
     * @brief Evaluate gradients into a flat caller-provided buffer
     *
     * Layout: gradients_out[i * 3 + c] = component c of gradient of basis i.
     * Caller provides a buffer of size() * 3 Real entries.
     */
    virtual void evaluate_gradients_to(const math::Vector<Real, 3>& xi,
                                       Real* SVMP_RESTRICT gradients_out) const;

    /**
     * @brief Evaluate Hessians into a flat caller-provided buffer
     *
     * Layout: hessians_out[i * 9 + r * 3 + c] = H_i(r, c).
     */
    virtual void evaluate_hessians_to(const math::Vector<Real, 3>& xi,
                                      Real* SVMP_RESTRICT hessians_out) const;

    /**
     * @brief Evaluate vector-valued basis functions (H(div)/H(curl))
     *
     * Default implementation throws; vector bases must override.
     */
    virtual void evaluate_vector_values(const math::Vector<Real, 3>& xi,
                                        std::vector<math::Vector<Real, 3>>& values) const;

    /**
     * @brief Evaluate reference-space Jacobians of vector-valued basis functions
     *
     * The returned matrix for basis function `i` has entries
     * `jacobians[i](component, derivative_direction) = d phi_i_component / d xi_direction`.
     * Unused rows/columns are zero-filled for lower-dimensional elements.
     */
    virtual void evaluate_vector_jacobians(const math::Vector<Real, 3>& xi,
                                           std::vector<VectorJacobian>& jacobians) const;

    /// Evaluate divergence of vector-valued basis functions (if applicable)
    virtual void evaluate_divergence(const math::Vector<Real, 3>& xi,
                                     std::vector<Real>& divergence) const;

    /// Evaluate curl of vector-valued basis functions (if applicable)
    virtual void evaluate_curl(const math::Vector<Real, 3>& xi,
                               std::vector<math::Vector<Real, 3>>& curl) const;

protected:
    /// Finite-difference helper for gradients of scalar bases
    void numerical_gradient(const math::Vector<Real, 3>& xi,
                            std::vector<Gradient>& gradients,
                            Real eps = Real(1e-6)) const;

    /// Finite-difference helper for Hessians of scalar bases
    void numerical_hessian(const math::Vector<Real, 3>& xi,
                           std::vector<Hessian>& hessians,
                           Real eps = Real(1e-5)) const;
};

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_BASISFUNCTION_H
