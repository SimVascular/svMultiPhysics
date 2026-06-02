/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#ifndef SVMP_FE_BASIS_VECTORBASIS_H
#define SVMP_FE_BASIS_VECTORBASIS_H

/**
 * @file VectorBasis.h
 * @brief Vector-valued bases for H(div) and H(curl) conforming spaces
 */

#include "BasisFunction.h"
#include "VectorBasisModalPolynomial.h"
#include <array>
#include <cstddef>
#include <vector>

namespace svmp {
namespace FE {
namespace basis {

/**
 * @brief DOF entity type for vector-valued basis functions
 */
enum class DofEntity {
    Vertex,   ///< DOF associated with a vertex
    Edge,     ///< DOF associated with an edge (tangential moments for H(curl))
    Face,     ///< DOF associated with a face (normal moments for H(div), tangential for H(curl))
    Interior  ///< DOF associated with element interior
};

/**
 * @brief DOF association metadata for a single DOF
 */
struct DofAssociation {
    DofEntity entity_type{DofEntity::Interior};
    int entity_id{-1};      ///< Local index of the entity (edge/face/vertex)
    int moment_index{0};    ///< Index within the entity's moment space
};

struct SparseModalCoefficientMatrix {
    std::size_t rows{0};
    std::size_t cols{0};
    std::vector<std::size_t> row_offsets;
    std::vector<std::size_t> dofs;
    std::vector<Real> coefficients;
};

class VectorBasisFunction : public BasisFunction {
public:
    bool is_vector_valued() const noexcept override { return true; }
    bool supports_vector_jacobians() const noexcept override { return true; }
    void evaluate_values(const math::Vector<Real, 3>&,
                         std::vector<Real>&) const override {
        throw BasisEvaluationException("Vector basis uses evaluate_vector_values",
                                       __FILE__, __LINE__, __func__);
    }

    void evaluate_vector_at_quadrature_points_strided(
        const std::vector<math::Vector<Real, 3>>& points,
        std::size_t output_stride,
        Real* SVMP_RESTRICT values_out,
        Real* SVMP_RESTRICT jacobians_out,
        Real* SVMP_RESTRICT curls_out,
        Real* SVMP_RESTRICT divergence_out) const override;

    /**
     * @brief Get DOF association metadata for all basis functions
     *
     * Returns a vector of size(), where each entry describes which
     * geometric entity (vertex/edge/face/interior) the corresponding
     * DOF is associated with. This is essential for orientation-aware
     * assembly of H(div) and H(curl) spaces.
     */
    virtual std::vector<DofAssociation> dof_associations() const {
        // Default: all interior DOFs (subclasses should override)
        std::vector<DofAssociation> result(size());
        for (std::size_t i = 0; i < size(); ++i) {
            result[i].entity_type = DofEntity::Interior;
            result[i].entity_id = 0;
            result[i].moment_index = static_cast<int>(i);
        }
        return result;
    }
};

/**
 * @brief Raviart-Thomas H(div) basis on supported element families
 */
class RaviartThomasBasis : public VectorBasisFunction {
public:
    RaviartThomasBasis(ElementType type, int order = 0);

    BasisType basis_type() const noexcept override { return BasisType::RaviartThomas; }
    ElementType element_type() const noexcept override { return element_type_; }
    int dimension() const noexcept override { return dimension_; }
    int order() const noexcept override { return order_; }
    std::size_t size() const noexcept override { return size_; }
    bool cache_identity_is_structural() const noexcept override { return true; }

    void evaluate_vector_values(const math::Vector<Real, 3>& xi,
                                std::vector<math::Vector<Real, 3>>& values) const override;
    void evaluate_vector_jacobians(const math::Vector<Real, 3>& xi,
                                   std::vector<VectorJacobian>& jacobians) const override;
    void evaluate_divergence(const math::Vector<Real, 3>& xi,
                             std::vector<Real>& divergence) const override;
    bool supports_divergence() const noexcept override { return true; }
    void evaluate_vector_at_quadrature_points_strided(
        const std::vector<math::Vector<Real, 3>>& points,
        std::size_t output_stride,
        Real* SVMP_RESTRICT values_out,
        Real* SVMP_RESTRICT jacobians_out,
        Real* SVMP_RESTRICT curls_out,
        Real* SVMP_RESTRICT divergence_out) const override;

    /// Get DOF associations (face/edge DOFs for 2D, face DOFs for 3D H(div))
    std::vector<DofAssociation> dof_associations() const override;

private:
    using ModalTerm = VectorBasisModalTerm;
    using ModalPolynomial = VectorBasisModalPolynomial;
    using SeedJacobianEvaluator = void (*)(
        const math::Vector<Real, 3>&,
        std::vector<VectorJacobian>&);

    ElementType element_type_;
    int dimension_;
    int order_;
    std::size_t size_{0};

    bool nodal_generated_{false};
    bool use_transformed_direct_seed_{false};  ///< True for wedge/pyramid RT(k=1,2) transformed from direct seed functions
    std::vector<int> transformed_seed_indices_;
    std::vector<std::array<int, 4>> transformed_monomial_candidates_; ///< {component, px, py, pz}
    std::vector<ModalPolynomial> monomials_;
    std::array<int, 3> modal_power_limits_{{0, 0, 0}};
    std::array<int, 3> transformed_power_limits_{{0, 0, 0}};
    SeedJacobianEvaluator transformed_seed_jacobian_evaluator_{nullptr};
    // Sparse coefficients for nodal basis in modal monomial basis:
    //   phi_j = sum_p c(p,j) * modal_p.
    // Rows index modal functions; entries target nodal DOFs.
    SparseModalCoefficientMatrix modal_sparse_coeffs_;
    SparseModalCoefficientMatrix transformed_sparse_coeffs_;
};

/**
 * @brief First-kind Nedelec H(curl) basis on supported element families
 */
class NedelecBasis : public VectorBasisFunction {
public:
    NedelecBasis(ElementType type, int order = 0);

    BasisType basis_type() const noexcept override { return BasisType::Nedelec; }
    ElementType element_type() const noexcept override { return element_type_; }
    int dimension() const noexcept override { return dimension_; }
    int order() const noexcept override { return order_; }
    std::size_t size() const noexcept override { return size_; }
    bool cache_identity_is_structural() const noexcept override { return true; }

    void evaluate_vector_values(const math::Vector<Real, 3>& xi,
                                std::vector<math::Vector<Real, 3>>& values) const override;
    void evaluate_vector_jacobians(const math::Vector<Real, 3>& xi,
                                   std::vector<VectorJacobian>& jacobians) const override;
    void evaluate_curl(const math::Vector<Real, 3>& xi,
                       std::vector<math::Vector<Real, 3>>& curl) const override;
    bool supports_curl() const noexcept override { return true; }
    void evaluate_vector_at_quadrature_points_strided(
        const std::vector<math::Vector<Real, 3>>& points,
        std::size_t output_stride,
        Real* SVMP_RESTRICT values_out,
        Real* SVMP_RESTRICT jacobians_out,
        Real* SVMP_RESTRICT curls_out,
        Real* SVMP_RESTRICT divergence_out) const override;

    /// Get DOF associations (edge DOFs for H(curl), face DOFs for 3D interior)
    std::vector<DofAssociation> dof_associations() const override;

private:
    using ModalTerm = VectorBasisModalTerm;
    using ModalPolynomial = VectorBasisModalPolynomial;
    using SeedJacobianEvaluator = void (*)(
        const math::Vector<Real, 3>&,
        std::vector<VectorJacobian>&);

    ElementType element_type_;
    int dimension_;
    int order_;
    std::size_t size_{0};

    bool nodal_generated_{false};
    bool use_transformed_direct_seed_{false};  ///< True for wedge/pyramid ND(k=1,2) transformed from direct seed/candidate functions
    std::vector<std::array<int, 4>> transformed_monomial_candidates_; ///< {component, px, py, pz}
    std::vector<ModalPolynomial> monomials_;
    SparseModalCoefficientMatrix modal_sparse_coeffs_;
    SparseModalCoefficientMatrix transformed_sparse_coeffs_;
    std::array<int, 3> modal_power_limits_{{0, 0, 0}};
    std::array<int, 3> transformed_power_limits_{{0, 0, 0}};
    SeedJacobianEvaluator transformed_seed_jacobian_evaluator_{nullptr};
};

/**
 * @brief Brezzi-Douglas-Marini basis (simple linear variant)
 */
class BDMBasis : public VectorBasisFunction {
public:
    BDMBasis(ElementType type, int order = 1);

    BasisType basis_type() const noexcept override { return BasisType::BDM; }
    ElementType element_type() const noexcept override { return element_type_; }
    int dimension() const noexcept override { return dimension_; }
    int order() const noexcept override { return order_; }
    std::size_t size() const noexcept override { return size_; }
    bool cache_identity_is_structural() const noexcept override { return true; }

    void evaluate_vector_values(const math::Vector<Real, 3>& xi,
                                std::vector<math::Vector<Real, 3>>& values) const override;
    void evaluate_vector_jacobians(const math::Vector<Real, 3>& xi,
                                   std::vector<VectorJacobian>& jacobians) const override;
    void evaluate_divergence(const math::Vector<Real, 3>& xi,
                             std::vector<Real>& divergence) const override;
    bool supports_divergence() const noexcept override { return true; }
    void evaluate_vector_at_quadrature_points_strided(
        const std::vector<math::Vector<Real, 3>>& points,
        std::size_t output_stride,
        Real* SVMP_RESTRICT values_out,
        Real* SVMP_RESTRICT jacobians_out,
        Real* SVMP_RESTRICT curls_out,
        Real* SVMP_RESTRICT divergence_out) const override;

    /// Get DOF associations (face/edge DOFs for H(div))
    std::vector<DofAssociation> dof_associations() const override;

private:
    using ModalTerm = VectorBasisModalTerm;
    using ModalPolynomial = VectorBasisModalPolynomial;

    ElementType element_type_;
    int dimension_;
    int order_;
    std::size_t size_{0};
    bool nodal_generated_{false};
    std::vector<ModalPolynomial> monomials_;
    SparseModalCoefficientMatrix modal_sparse_coeffs_;
    std::array<int, 3> modal_power_limits_{{0, 0, 0}};
};

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_VECTORBASIS_H
