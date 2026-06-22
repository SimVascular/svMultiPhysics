// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_BASIS_LAGRANGEBASIS_H
#define SVMP_FE_BASIS_LAGRANGEBASIS_H

#include "BasisFunction.h"
#include "BasisTraits.h"

#include <array>
#include <cstddef>
#include <span>

namespace svmp {
namespace FE {
namespace basis {

/// \defgroup FE_LagrangeBasis LagrangeBasis
/// \ingroup FE_Basis
/// \brief Construction and evaluation API for nodal Lagrange finite-element bases.
///
/// \details This group documents the complete nodal Lagrange basis evaluator
/// used by the FE library. The implementation covers tensor-product,
/// simplex, and wedge reference topologies with exact analytical first and
/// second derivatives in reference coordinates.
/// @{

/// \brief Nodal Lagrange basis on supported reference finite elements.
///
/// \details LagrangeBasis represents the nodal interpolation basis associated
/// with an equispaced reference-node lattice. It supports point, line,
/// quadrilateral, hexahedron, triangle, tetrahedron, and wedge reference
/// elements. Named complete quadratic elements such as Line3, Triangle6,
/// Quad9, Tetra10, Hex27, and Wedge18 are normalized to their canonical
/// linear topology plus effective order 2.
///
/// Tensor-product elements use the one-dimensional nodal polynomials
/// \f[
///   l_i(x) = \prod_{j \ne i} \frac{x - x_j}{x_i - x_j}
/// \f]
/// on equispaced coordinates in \f$[-1, 1]\f$. Multi-dimensional basis
/// functions are products of the active axis polynomials, for example
/// \f$N_{ijk}(r,s,t) = l_i(r)l_j(s)l_k(t)\f$ on a hexahedron.
///
/// Simplex elements use barycentric coordinates and integer lattice
/// exponents. For a node with exponent tuple \f$\alpha\f$, where
/// \f$\sum_a \alpha_a = p\f$, the basis is assembled from scaled
/// falling-factorial factors,
/// \f[
///   N_\alpha(\lambda) =
///   \prod_a \prod_{m=0}^{\alpha_a-1}
///   \frac{p\lambda_a - m}{m + 1}.
/// \f]
/// Gradients and Hessians are evaluated analytically by differentiating these
/// factors and applying the barycentric-coordinate chain rule.
///
/// Wedge elements are treated as a tensor product between a triangle simplex
/// basis and a one-dimensional through-axis basis:
/// \f$N_{a k}(r,s,t) = T_a(r,s)l_k(t)\f$.
///
/// The vector-returning evaluators are convenient API wrappers. The `*_to`
/// methods write to caller-provided spans and are intended for assembly paths
/// that avoid temporary allocations.
class LagrangeBasis final : public BasisFunction {
public:
    /// \brief Axis-index tuple for tensor-product reference nodes.
    using TensorNodeIndex = std::array<std::size_t, 3>;

    /// \brief Barycentric exponent tuple for simplex reference nodes.
    using SimplexExponent = std::array<int, 4>;

    /// \brief Triangle-node and axis-node tuple for wedge reference nodes.
    using WedgeNodeIndex = std::array<std::size_t, 2>;

    /// \brief Construct a Lagrange basis for an element type and polynomial order.
    ///
    /// \details The constructor normalizes complete higher-order aliases to the
    /// canonical topology and effective polynomial order, builds the reference
    /// node coordinates, and precomputes topology-specific lookup data used by
    /// evaluation. Tensor-product bases store per-axis node indices, simplex
    /// bases store barycentric exponent tuples, and wedge bases store the
    /// triangle-node/axis-node decomposition.
    ///
    /// \param type Element type used to determine topology and reference-node layout.
    /// \param order Requested polynomial order.
    /// \throws BasisConfigurationException If the effective order is negative.
    /// \throws BasisElementCompatibilityException If the element type is unsupported.
    LagrangeBasis(ElementType type, int order);

    /// \copydoc BasisFunction::basis_type()
    BasisType basis_type() const noexcept final { return BasisType::Lagrange; }

    /// \copydoc BasisFunction::element_type()
    ElementType element_type() const noexcept final { return element_type_; }

    /// \copydoc BasisFunction::dimension()
    int dimension() const noexcept final { return dimension_; }

    /// \copydoc BasisFunction::order()
    int order() const noexcept final { return order_; }

    /// \copydoc BasisFunction::size()
    std::size_t size() const noexcept final { return nodes_.size(); }

    /// \brief Return the reference interpolation nodes in basis ordering.
    ///
    /// \details The returned node order matches the basis-function order used
    /// by all evaluators. Coordinates are reference-element coordinates:
    /// tensor-product axes use \f$[-1,1]\f$, triangles and tetrahedra use the
    /// repository's simplex reference coordinates, and wedges combine triangle
    /// reference coordinates with a \f$[-1,1]\f$ through-axis coordinate.
    ///
    /// \return Reference node coordinates, one per basis function.
    const std::vector<math::Vector<double, 3>>& nodes() const noexcept final { return nodes_; }

    /// \brief Evaluate Lagrange basis function values at a reference coordinate.
    ///
    /// \details Values satisfy the nodal interpolation property
    /// \f$N_i(x_j)=\delta_{ij}\f$ at the basis nodes. Tensor-product values are
    /// products of one-dimensional Lagrange polynomials. Simplex values are
    /// products of barycentric falling-factorial factors. Wedge values are
    /// products of triangle simplex values and through-axis Lagrange values.
    ///
    /// \param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
    /// \param values Receives one value per basis function.
    void evaluate_values(const math::Vector<double, 3>& xi,
                         std::vector<double>& values) const final;

    /// \brief Evaluate analytical Lagrange basis gradients at a reference coordinate.
    ///
    /// \details Gradients are derivatives with respect to reference
    /// coordinates, not physical coordinates. Tensor-product gradients apply
    /// the product rule to the active axis polynomials. Simplex gradients
    /// differentiate the barycentric factors and multiply by the constant
    /// gradients of the barycentric coordinates. Wedge gradients combine the
    /// triangle gradient in the first two components with the through-axis
    /// derivative in the third component.
    ///
    /// \param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
    /// \param gradients Receives one three-component gradient per basis function.
    void evaluate_gradients(const math::Vector<double, 3>& xi,
                            std::vector<Gradient>& gradients) const final;

    /// \brief Evaluate analytical Lagrange basis Hessians at a reference coordinate.
    ///
    /// \details Hessians are second derivatives in reference coordinates and
    /// are stored as 3-by-3 matrices. Tensor-product Hessians contain pure
    /// second axis derivatives on the diagonal and mixed product-rule terms
    /// off diagonal. Simplex Hessians are assembled from first and second
    /// derivatives of the barycentric factors. Wedge Hessians contain triangle
    /// Hessian terms, through-axis second derivatives, and mixed
    /// triangle/through-axis derivative products.
    ///
    /// \param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
    /// \param hessians Receives one 3-by-3 Hessian per basis function.
    void evaluate_hessians(const math::Vector<double, 3>& xi,
                           std::vector<Hessian>& hessians) const final;

    /// \brief Evaluate Lagrange values, gradients, and Hessians together.
    ///
    /// \details This is the allocation-friendly vector API for callers that
    /// need all basis quantities at the same quadrature point. The underlying
    /// evaluator computes only topology-local polynomial data once and then
    /// fills all requested outputs.
    ///
    /// \param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
    /// \param values Receives one value per basis function.
    /// \param gradients Receives one three-component gradient per basis function.
    /// \param hessians Receives one 3-by-3 Hessian per basis function.
    void evaluate_all(const math::Vector<double, 3>& xi,
                      std::vector<double>& values,
                      std::vector<Gradient>& gradients,
                      std::vector<Hessian>& hessians) const final;

    /// \brief Evaluate Lagrange basis values into caller-provided storage.
    ///
    /// \details This is the low-allocation API intended for element assembly
    /// loops. The span is filled in basis-node order and no vector resizing is
    /// performed.
    ///
    /// \param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
    /// \param values_out Output span with at least size() entries.
    void evaluate_values_to(const math::Vector<double, 3>& xi,
                            std::span<double> values_out) const final;

    /// \brief Evaluate Lagrange basis gradients into caller-provided storage.
    ///
    /// \details Gradients are written in basis-node order with one
    /// three-component gradient per node.
    ///
    /// \param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
    /// \param gradients_out Output span with at least size() entries.
    void evaluate_gradients_to(const math::Vector<double, 3>& xi,
                               std::span<Gradient> gradients_out) const final;

    /// \brief Evaluate Lagrange basis Hessians into caller-provided storage.
    ///
    /// \details Hessians are written in basis-node order with one 3-by-3
    /// Hessian per node.
    ///
    /// \param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
    /// \param hessians_out Output span with at least size() entries.
    void evaluate_hessians_to(const math::Vector<double, 3>& xi,
                              std::span<Hessian> hessians_out) const final;

private:
    ElementType element_type_;
    BasisTopology topology_{BasisTopology::Unknown};
    int dimension_{0};
    int order_{0};

    std::vector<double> nodes_1d_;
    std::vector<math::Vector<double, 3>> nodes_;
    std::vector<TensorNodeIndex> tensor_indices_;
    std::vector<SimplexExponent> simplex_exponents_;
    std::vector<WedgeNodeIndex> wedge_indices_;

    void init_nodes();
    void build_point_nodes();
    void build_tensor_product_nodes();
    void build_simplex_nodes();
    void build_wedge_nodes();
    void init_equispaced_1d_nodes();

    void evaluate_all_to(const math::Vector<double, 3>& xi,
                         std::span<double> values_out,
                         std::span<Gradient> gradients_out,
                         std::span<Hessian> hessians_out) const;
    void evaluate_point_to(std::span<double> values_out,
                           std::span<Gradient> gradients_out,
                           std::span<Hessian> hessians_out) const;
    void evaluate_tensor_product_to(const math::Vector<double, 3>& xi,
                                    std::span<double> values_out,
                                    std::span<Gradient> gradients_out,
                                    std::span<Hessian> hessians_out) const;
    void evaluate_simplex_to(const math::Vector<double, 3>& xi,
                             std::span<double> values_out,
                             std::span<Gradient> gradients_out,
                             std::span<Hessian> hessians_out) const;
    void evaluate_wedge_to(const math::Vector<double, 3>& xi,
                           std::span<double> values_out,
                           std::span<Gradient> gradients_out,
                           std::span<Hessian> hessians_out) const;
};

/// @}

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_LAGRANGEBASIS_H
