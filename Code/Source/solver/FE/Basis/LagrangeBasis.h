// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_BASIS_LAGRANGEBASIS_H
#define SVMP_FE_BASIS_LAGRANGEBASIS_H

#include "BasisFunction.h"
#include "BasisTraits.h"

#include <array>
#include <cstddef>
#include <span>

namespace svmp::FE::basis {

/**
 * @defgroup FE_LagrangeBasis LagrangeBasis
 * @ingroup FE_Basis
 * @brief Construction and evaluation API for nodal Lagrange finite-element bases.
 *
 * @details This group documents the complete nodal Lagrange basis evaluator
 * used by the FE library. The implementation covers tensor-product,
 * simplex, and wedge reference topologies with exact analytical first and
 * second derivatives in reference coordinates.
 * @{
 */

/**
 * @brief Nodal Lagrange basis on supported reference finite elements.
 *
 * @details LagrangeBasis represents the complete (full-degree) nodal
 * interpolation basis on a reference topology. It supports point, line,
 * quadrilateral, hexahedron, triangle, tetrahedron, and wedge reference
 * topologies. The primary constructor takes a BasisTopology and an explicit
 * polynomial order, so an arbitrary order carries no node-count assumption
 * (an order-2 hexahedron is BasisTopology::Hexahedron with order 2). A named
 * ElementType such as Line3, Quad9, Tetra10, or Hex27 is a fixed-order
 * shorthand: it maps to the same (topology, order) pair and the requested order
 * must equal the order baked into that layout (1 for the linear elements, 2 for
 * the complete-quadratic aliases, 0 for Point1).
 *
 * ## Reference-node distribution
 *
 * The interpolation nodes are not a single distribution across topologies; each
 * family uses the node set its evaluator is built for:
 * - **Tensor-product (line, quadrilateral, hexahedron):** the shared
 *   Gauss-Lobatto-Legendre (GLL) tensor-axis nodes -- see line_coord_pm_one for
 *   the distribution and its conditioning -- not an equispaced layout.
 * - **Simplex (triangle, tetrahedron):** the equispaced barycentric lattice
 *   (each barycentric coordinate at @f$i/p@f$). The closed-form evaluator below
 *   is specific to this equispaced lattice.
 * - **Wedge:** the tensor product of an equispaced triangle cross-section with a
 *   GLL through-axis.
 *
 * Because GLL coincides with the equispaced layout at orders 1 and 2
 * (line_coord_pm_one), the linear and quadratic tensor elements -- Line2/Line3,
 * Quad4/Quad9, Hex8/Hex27, and the wedge through-axis -- are built on equispaced
 * nodes, and the GLL/equispaced distinction appears only for order >= 3.
 *
 * ## Evaluation
 *
 * Tensor-product elements use the one-dimensional nodal polynomials
 * @f[
 *   l_i(x) = \prod_{j \ne i} \frac{x - x_j}{x_i - x_j}
 * @f]
 * on the per-axis GLL coordinates @f$x_j \in [-1, 1]@f$ (the barycentric-weight
 * form, valid for any distinct node set). Multi-dimensional basis functions are
 * products of the active axis polynomials, for example
 * @f$N_{ijk}(r,s,t) = l_i(r)l_j(s)l_k(t)@f$ on a hexahedron.
 *
 * Simplex elements use barycentric coordinates and integer lattice
 * exponents on the equispaced lattice. For a node with exponent tuple
 * @f$\alpha@f$, where
 * @f$\sum_a \alpha_a = p@f$, the basis is assembled from scaled
 * falling-factorial factors,
 * @f[
 *   N_\alpha(\lambda) =
 *   \prod_a \prod_{m=0}^{\alpha_a-1}
 *   \frac{p\lambda_a - m}{m + 1}.
 * @f]
 * Gradients and Hessians are evaluated analytically by differentiating these
 * factors and applying the barycentric-coordinate chain rule.
 *
 * Wedge elements are treated as a tensor product between a triangle simplex
 * basis and a one-dimensional through-axis basis:
 * @f$N_{a k}(r,s,t) = T_a(r,s)l_k(t)@f$.
 *
 * ## Conditioning and the supported order range
 *
 * Interpolation conditioning is governed by the node distribution and so differs
 * by topology:
 * - **Tensor-product topologies stay well-conditioned at high order.** GLL nodes
 *   have a logarithmic Lebesgue constant, so line/quadrilateral/hexahedron bases
 *   remain trustworthy well beyond the production orders.
 * - **Simplex topologies degrade at high order.** The equispaced barycentric
 *   lattice has a Lebesgue constant that grows roughly exponentially with order
 *   (the Runge phenomenon), so triangle and tetrahedron bases are reliable
 *   through low orders but become increasingly ill-conditioned beyond them. The
 *   wedge inherits this through its equispaced triangle cross-section.
 *
 * The vector-returning evaluators are convenient API wrappers. The `*_to`
 * methods write to caller-provided spans and are intended for assembly paths
 * that avoid temporary allocations.
 */
class LagrangeBasis final : public BasisFunction {
public:
    /** @brief Axis-index tuple for tensor-product reference nodes. */
    using TensorNodeIndex = std::array<std::size_t, 3>;

    /** @brief Barycentric exponent tuple for simplex reference nodes. */
    using SimplexExponent = std::array<int, 4>;

    /** @brief Triangle-node and axis-node tuple for wedge reference nodes. */
    using WedgeNodeIndex = std::array<std::size_t, 2>;

    /**
     * @brief Construct a Lagrange basis on a reference topology at a polynomial order.
     *
     * @details This is the primary, arbitrary-order entry point: a BasisTopology
     * carries no node-count assumption, so any supported order is requested
     * explicitly (e.g. an order-5 hexahedron is BasisTopology::Hexahedron with
     * order 5). The constructor builds the reference node coordinates and the
     * topology-specific lookup data used by evaluation. Tensor-product bases
     * store per-axis node indices, simplex bases store barycentric exponent
     * tuples, and wedge bases store the triangle-node/axis-node decomposition.
     *
     * Reference nodes follow the per-topology distribution described in the class
     * documentation (Reference-node distribution). Unlike SerendipityBasis, this
     * constructor does not reject ill-conditioned high-order simplex/wedge requests
     * (where the equispaced barycentric lattice degrades); that choice is the
     * caller's.
     *
     * @param topology Reference topology; Point through the volume topologies.
     * @param order Polynomial order; must be non-negative. Point is order 0.
     * @throws BasisConfigurationException If the order is negative, or if Point
     *         is requested with a nonzero order.
     * @throws BasisElementCompatibilityException If the topology is Unknown.
     */
    LagrangeBasis(BasisTopology topology, int order);

    /**
     * @brief Construct a Lagrange basis from a named element layout.
     *
     * @details Convenience overload for a named mesh element. The order is baked
     * into the layout (0 for Point1, 1 for the linear elements, 2 for the
     * complete-quadratic aliases such as Hex27/Tetra10) and the requested
     * @p order must match it; arbitrary orders must be requested through the
     * BasisTopology overload. Serendipity and pyramid layouts are rejected.
     *
     * @param type Named element type used to determine topology and baked-in order.
     * @param order Requested order; must equal the element's baked-in order.
     * @throws BasisConfigurationException If @p order does not match the element's baked-in order.
     * @throws BasisElementCompatibilityException If the element type is unsupported.
     */
    LagrangeBasis(ElementType type, int order);

    /**
     * @brief Construct a Lagrange basis from a named element layout at its baked-in order.
     *
     * @details Single-argument convenience overload: the polynomial order is the
     * one baked into the layout (0 for Point1, 1 for the linear elements, 2 for
     * the complete-quadratic aliases such as Hex27/Tetra10), so the caller does
     * not repeat it. Equivalent to LagrangeBasis(type, <baked-in order>).
     * Serendipity and pyramid layouts are rejected, as for the two-argument
     * overload.
     *
     * @param type Named element type; determines both topology and order.
     * @throws BasisElementCompatibilityException If the element type is unsupported.
     */
    explicit LagrangeBasis(ElementType type);

    /** @copydoc BasisFunction::basis_type() */
    BasisType basis_type() const noexcept final { return BasisType::Lagrange; }

    /** @copydoc BasisFunction::topology() */
    BasisTopology topology() const noexcept final { return topology_; }

    /** @copydoc BasisFunction::dimension() */
    int dimension() const noexcept final { return dimension_; }

    /** @copydoc BasisFunction::order() */
    int order() const noexcept final { return order_; }

    /** @copydoc BasisFunction::size() */
    std::size_t size() const noexcept final { return nodes_.size(); }

    /**
     * @brief Return the reference interpolation nodes in basis ordering.
     *
     * @details The returned node order matches the basis-function order used by
     * all evaluators; the coordinates follow the per-topology distribution
     * described in the class documentation (Reference-node distribution).
     *
     * @return Reference node coordinates, one per basis function.
     */
    const std::vector<math::Vector<double, 3>>& nodes() const noexcept final { return nodes_; }

    /**
     * @brief Evaluate Lagrange basis values into caller-provided storage.
     *
     * @details This is the low-allocation API intended for element assembly
     * loops. The span is filled in basis-node order and no vector resizing is
     * performed.
     *
     * @param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
     * @param values_out Output span with at least size() entries.
     */
    void evaluate_values_to(const math::Vector<double, 3>& xi,
                            std::span<double> values_out) const final;

    /**
     * @brief Evaluate Lagrange basis gradients into caller-provided storage.
     *
     * @details Gradients are written in basis-node order with one
     * three-component gradient per node.
     *
     * @param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
     * @param gradients_out Output span with at least size() entries.
     */
    void evaluate_gradients_to(const math::Vector<double, 3>& xi,
                               std::span<Gradient> gradients_out) const final;

    /**
     * @brief Evaluate Lagrange basis Hessians into caller-provided storage.
     *
     * @details Hessians are written in basis-node order with one 3-by-3
     * Hessian per node.
     *
     * @param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
     * @param hessians_out Output span with at least size() entries.
     */
    void evaluate_hessians_to(const math::Vector<double, 3>& xi,
                              std::span<Hessian> hessians_out) const final;

private:
    BasisTopology topology_{BasisTopology::Unknown};
    int dimension_{0};
    int order_{0};

    // Topology-specific construction data. nodes_ (the reference nodes in basis
    // order) is populated for every topology and backs size(); each remaining
    // vector is filled only for the topologies that use it and stays empty
    // otherwise:
    //   line/quad/hex : nodes_1d_, nodes_1d_weights_, tensor_indices_
    //   triangle/tetra: simplex_exponents_
    //   wedge         : nodes_1d_, nodes_1d_weights_, wedge_indices_, and
    //                   simplex_exponents_ (the triangle cross-section exponents)
    //   point         : nodes_ only
    std::vector<double> nodes_1d_;
    std::vector<double> nodes_1d_weights_;
    std::vector<math::Vector<double, 3>> nodes_;
    std::vector<TensorNodeIndex> tensor_indices_;
    std::vector<SimplexExponent> simplex_exponents_;
    std::vector<WedgeNodeIndex> wedge_indices_;

    void init_nodes();
    void build_point_nodes();
    void build_tensor_product_nodes();
    void build_simplex_nodes();
    void build_wedge_nodes();
    void init_tensor_axis_nodes();

    void evaluate_all_to(const math::Vector<double, 3>& xi,
                         std::span<double> values_out,
                         std::span<Gradient> gradients_out,
                         std::span<Hessian> hessians_out) const override;
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

/** @} */

} // namespace svmp::FE::basis

#endif // SVMP_FE_BASIS_LAGRANGEBASIS_H
