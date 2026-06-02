/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#ifndef SVMP_FE_BASIS_LAGRANGEBASIS_H
#define SVMP_FE_BASIS_LAGRANGEBASIS_H

/**
 * @file LagrangeBasis.h
 * @brief Nodal Lagrange polynomial basis on reference elements
 */

#include "BasisFunction.h"
#include <array>
#include <cstddef>

namespace svmp {
namespace FE {
namespace basis {

void prewarm_lagrange_basis_scratch(int max_order, std::size_t max_qpts = 0);

/**
 * @brief Complete nodal H1 Lagrange basis on canonical reference topologies
 *
 * Supports arbitrary polynomial order on the canonical complete families:
 * `Line2`, `Triangle3`, `Quad4`, `Tetra4`, `Hex8`, `Wedge6`, and `Pyramid5`.
 * Low-order complete-family aliases (`Line3`, `Triangle6`, `Quad9`,
 * `Tetra10`, `Hex27`, `Wedge18`, `Pyramid14`) normalize to their canonical
 * topology plus order. Serendipity variants remain intentionally excluded.
 *
 * Node locations are generated on canonical reference elements using
 * equispaced coordinates on tensor-product elements, barycentric grids on
 * simplices, tensorized triangle-line grids on wedges, and a rational nodal
 * pyramid construction on `Pyramid5`.
 *
 * The evaluator is numerically stabilized for those nodes, but the
 * interpolation problem itself remains the equispaced Lagrange problem. For
 * high-order interpolation, especially order >= 4, prefer `SpectralBasis`
 * (GLL / Warp & Blend nodes) unless exact equispaced nodal placement is part
 * of the requested discretization.
 *
 * For the rational pyramid family, basis values remain exact at the apex.
 * Gradients and Hessians are analytic on the supported interior reference
 * domain, but the exact-apex nodal derivative limit is not unique and those
 * derivative queries throw at the exact apex.
 */
class LagrangeBasis : public BasisFunction {
public:
    LagrangeBasis(ElementType type, int order);

    BasisType basis_type() const noexcept override { return BasisType::Lagrange; }
    ElementType element_type() const noexcept override { return element_type_; }
    int dimension() const noexcept override { return dimension_; }
    int order() const noexcept override { return order_; }
    std::size_t size() const noexcept override { return nodes_.size(); }
    bool cache_identity_is_structural() const noexcept override { return true; }

    const std::vector<math::Vector<Real, 3>>& nodes() const noexcept { return nodes_; }

    void evaluate_values(const math::Vector<Real, 3>& xi,
                         std::vector<Real>& values) const final;
    void evaluate_gradients(const math::Vector<Real, 3>& xi,
                            std::vector<Gradient>& gradients) const final;
    void evaluate_hessians(const math::Vector<Real, 3>& xi,
                           std::vector<Hessian>& hessians) const final;
    void evaluate_all(const math::Vector<Real, 3>& xi,
                      std::vector<Real>& values,
                      std::vector<Gradient>& gradients,
                      std::vector<Hessian>& hessians) const final;

    void evaluate_at_quadrature_points(
        const std::vector<math::Vector<Real, 3>>& points,
        Real* SVMP_RESTRICT values_out,
        Real* SVMP_RESTRICT gradients_out,
        Real* SVMP_RESTRICT hessians_out) const final;
    void evaluate_at_quadrature_points_strided(
        const std::vector<math::Vector<Real, 3>>& points,
        std::size_t output_stride,
        Real* SVMP_RESTRICT values_out,
        Real* SVMP_RESTRICT gradients_out,
        Real* SVMP_RESTRICT hessians_out) const final;

    // Raw-pointer output API. Caller must pre-size buffers to size().
    void evaluate_values_to(const math::Vector<Real, 3>& xi, Real* SVMP_RESTRICT values_out) const final;
    void evaluate_gradients_to(const math::Vector<Real, 3>& xi, Real* SVMP_RESTRICT gradients_out) const final;
    void evaluate_hessians_to(const math::Vector<Real, 3>& xi, Real* SVMP_RESTRICT hessians_out) const final;

private:
    using TensorNodeIndex = std::array<std::size_t, 3>;
    using WedgeNodeIndex = std::array<std::size_t, 2>;
    using VectorEvaluationDispatch = void (LagrangeBasis::*)(
        const math::Vector<Real, 3>&,
        std::vector<Real>*,
        std::vector<Gradient>*,
        std::vector<Hessian>*) const;

    // Cached topology encoded as int because the topology enum lives in
    // the .cpp anon namespace. Set once in init_nodes.
    int topology_id_ = 0;

    ElementType element_type_;
    int dimension_;
    int order_;

    std::vector<Real> nodes_1d_;
    std::vector<math::Vector<Real, 3>> nodes_;
    std::vector<TensorNodeIndex> tensor_indices_;
    std::vector<std::array<int, 4>> simplex_exponents_;
    std::vector<WedgeNodeIndex> wedge_indices_;
    std::vector<std::size_t> wedge_node_by_tri_z_;

    // Precomputed Horner-form coefficients of the 1D Lagrange basis.
    // Layout per axis (n_axis = nodes_1d_.size() = order_+1):
    //   axis_v_coeffs_[i * n_axis + k] = coeff of x^k in L_i(x), 0 <= i,k < n_axis
    //   axis_d_coeffs_[i * (n_axis - 1) + k] = coeff of x^k in L_i'(x)
    //   axis_d2_coeffs_[i * (n_axis - 2) + k] = coeff of x^k in L_i''(x)  (only if n_axis >= 3)
    // Populated by build_tensor_product_nodes / build_wedge_nodes.
    std::vector<Real> axis_v_coeffs_;
    std::vector<Real> axis_d_coeffs_;
    std::vector<Real> axis_d2_coeffs_;
    std::vector<Real> axis_barycentric_weights_;
    VectorEvaluationDispatch vector_evaluation_dispatch_{nullptr};

    void init_nodes();
    void init_evaluation_dispatch();
    void build_point_nodes();
    void build_tensor_product_nodes(int dimensions);
    void build_simplex_nodes();
    void build_wedge_nodes();
    void build_pyramid_nodes();
    void init_equispaced_1d_nodes();
    void compute_axis_monomial_coefficients();
    void evaluate_point_vectors(const math::Vector<Real, 3>& xi,
                                std::vector<Real>* values,
                                std::vector<Gradient>* gradients,
                                std::vector<Hessian>* hessians) const;
    void evaluate_tensor_product_vectors(const math::Vector<Real, 3>& xi,
                                         std::vector<Real>* values,
                                         std::vector<Gradient>* gradients,
                                         std::vector<Hessian>* hessians) const;
    void evaluate_triangle_vectors(const math::Vector<Real, 3>& xi,
                                   std::vector<Real>* values,
                                   std::vector<Gradient>* gradients,
                                   std::vector<Hessian>* hessians) const;
    void evaluate_tetrahedron_vectors(const math::Vector<Real, 3>& xi,
                                      std::vector<Real>* values,
                                      std::vector<Gradient>* gradients,
                                      std::vector<Hessian>* hessians) const;
    void evaluate_wedge_vectors(const math::Vector<Real, 3>& xi,
                                std::vector<Real>* values,
                                std::vector<Gradient>* gradients,
                                std::vector<Hessian>* hessians) const;
    void evaluate_pyramid_vectors(const math::Vector<Real, 3>& xi,
                                  std::vector<Real>* values,
                                  std::vector<Gradient>* gradients,
                                  std::vector<Hessian>* hessians) const;
    void evaluate_unsupported_vectors(const math::Vector<Real, 3>& xi,
                                      std::vector<Real>* values,
                                      std::vector<Gradient>* gradients,
                                      std::vector<Hessian>* hessians) const;
    void evaluate_all_to(const math::Vector<Real, 3>& xi,
                         Real* SVMP_RESTRICT values_out,
                         Real* SVMP_RESTRICT gradients_out,
                         Real* SVMP_RESTRICT hessians_out) const;
};

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_LAGRANGEBASIS_H
