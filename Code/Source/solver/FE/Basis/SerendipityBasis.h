// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_BASIS_SERENDIPITYBASIS_H
#define SVMP_FE_BASIS_SERENDIPITYBASIS_H

/**
 * @file SerendipityBasis.h
 * @brief Reduced-degree-of-freedom serendipity bases
 */

#include "BasisFunction.h"

#include <array>
#include <span>

namespace svmp::FE::basis {

/**
 * @defgroup FE_SerendipityBasis SerendipityBasis
 * @ingroup FE_Basis
 * @brief Construction and evaluation API for reduced serendipity finite-element bases.
 *
 * @details This group documents reduced degree-of-freedom basis families that
 * preserve nodal interpolation on supported element boundaries while omitting
 * selected interior tensor-product modes. These bases are used for standard
 * serendipity elements and geometry-mode mappings that intentionally use a
 * lower-order interpolation space.
 * @{
 */

/**
 * @brief Reduced-degree-of-freedom serendipity basis on supported reference elements.
 *
 * @details SerendipityBasis implements nodal bases for the quadrilateral and
 * hexahedral serendipity families at arbitrary order, plus the Wedge15 prism
 * layout. Compared with a complete tensor-product Lagrange basis of the same
 * nominal order, a serendipity basis removes selected interior modes while
 * retaining nodal interpolation on the supported node layout. The named layouts
 * Quad8, Hex8, and Hex20 are the fixed-order instances of these families
 * (quadrilateral order 2, hexahedron orders 1 and 2).
 *
 * The quadrilateral serendipity polynomial space is described by monomials
 * @f$x^{a_x}y^{a_y}@f$ whose superlinear degree is at most the requested
 * order. The implementation evaluates this space through tensor Legendre
 * modes, which span the same polynomial space but give a better-conditioned
 * Vandermonde. The superlinear degree is
 * @f[
 *   sldeg(x^{a_x}y^{a_y}) =
 *   \begin{cases} a_x, & a_x > 1 \\ 0, & a_x \le 1 \end{cases}
 *   +
 *   \begin{cases} a_y, & a_y > 1 \\ 0, & a_y \le 1 \end{cases}.
 * @f]
 * The nodal basis is recovered by inverting the generalized Vandermonde
 * interpolation matrix at the selected reference nodes. Values, gradients, and
 * Hessians are then evaluated by differentiating the modal vector and applying
 * the inverse Vandermonde coefficients.
 * For order @f$p \ge 1@f$, this space has @f$4p@f$ boundary modes for
 * @f$p \le 3@f$ and
 * @f[
 *   4p + \frac{(p - 3)(p - 2)}{2}
 * @f]
 * modes for @f$p \ge 4@f$.
 *
 * The quadrilateral node set is unisolvent by construction. If
 * @f$s(x,y)@f$ in this space vanishes at the @f$p + 1@f$ distinct nodes on
 * every edge, each edge restriction is a degree-@f$p@f$ one-variable
 * polynomial with @f$p + 1@f$ roots, so all edge restrictions vanish. Thus
 * @f$s@f$ is divisible by the boundary bubble
 * @f$(1 - x^2)(1 - y^2)@f$, and the quotient lies in
 * @f$P_{p-4}@f$ (with no quotient for @f$p < 4@f$). For @f$p \ge 4@f$, the
 * interior nodes form triangular rows for @f$P_{p-4}@f$: the first row has
 * @f$m + 1@f$ distinct @f$x@f$ values, the next row has @f$m@f$, and so on
 * for @f$m = p - 4@f$. A total-degree polynomial that vanishes on those rows
 * is zero by induction over rows, because each vanished row factors out one
 * linear term in @f$y@f$. The interpolation Vandermonde is therefore
 * nonsingular for the implemented quadrilateral serendipity space.
 *
 * Hexahedral serendipity generalizes the same construction to the cube. The
 * polynomial space is described by every monomial
 * @f$r^{a_r}s^{a_s}t^{a_t}@f$ whose superlinear degree (the three-axis form of
 * the rule above) is at most @f$p@f$, and the nodal basis is again the inverse
 * Vandermonde at the reference nodes. Those nodes are
 * distributed by boundary stratum: 8 corners, @f$12(p-1)@f$ edge nodes,
 * @f$6\,q(p)@f$ face-interior nodes -- each face carries the 2D quadrilateral
 * serendipity interior, since the trace of the cube space on a face is the
 * square space -- and a volume interior that is empty until @f$p \ge 6@f$.
 * Unisolvence follows the same factorization: a function vanishing on every
 * boundary node vanishes on each face by the quadrilateral result above, hence
 * is divisible by the cube bubble @f$(1 - r^2)(1 - s^2)(1 - t^2)@f$ with quotient
 * in @f$P_{p-6}@f$; the volume-interior nodes form a tetrahedral staircase that
 * is unisolvent for @f$P_{p-6}@f$ by induction over @f$t@f$-layers, so the cube
 * Vandermonde is nonsingular.
 *
 * `SerendipityBasis(BasisTopology::Quadrilateral, p)` and
 * `SerendipityBasis(BasisTopology::Hexahedron, p)` are the arbitrary-order entry
 * points (@f$p \ge 1@f$; orders below one are rejected). Reference nodes for both
 * the arbitrary-order and the named paths come from the single
 * ReferenceNodeLayout serendipity generator, in a VTK-consistent stratified
 * order; for @f$p \ge 3@f$ the interior ordering is an implementation convention
 * rather than a public layout. The named fixed layouts -- `ElementType::Quad8`
 * (order 2), `Hex8` (order 1), and `Hex20` (order 2) -- are the same construction
 * at those orders; the named overload only pins the order, so the named and
 * topology constructions produce identical objects and share the single public
 * node ordering the solver permutes against (order 1 and order 2 reuse the VTK
 * corner/edge ordering exactly). Wedge serendipity remains a single fixed layout
 * (Wedge15), constructed only from its named ElementType. Solver-default basis
 * selection is separate: `basis_factory` maps the complete Quad4 layout to the
 * default linear Lagrange basis and maps Quad8/Hex20 to serendipity unless a
 * caller explicitly requests a different supported basis.
 *
 * Every supported family -- quadrilateral, hexahedral, and Wedge15 -- is built by
 * inverting the generalized Vandermonde of its mode space at the public-order
 * reference nodes. Quadrilateral and hexahedral bases use tensor Legendre modes;
 * the fixed Wedge15 table uses monomial modes. Values, gradients, and Hessians
 * are evaluated by differentiating the matching mode vector and applying the
 * inverse-Vandermonde coefficients. Because the tables are generated in public
 * node order, evaluation needs no output reordering, and there is no hand-written
 * special case -- the Hex8 basis is the order-1 instance of the generated
 * hexahedral space, not a separate trilinear evaluator.
 *
 * ## Conditioning and the well-conditioned order range
 *
 * High-order nodal interpolation is governed by two conditioning factors, both
 * addressed so that arbitrary orders produce trustworthy shape functions:
 * - **Node distribution.** The quadrilateral and hexahedral families place their
 *   nodes on the shared Gauss-Lobatto-Legendre (GLL) distribution -- edges, faces,
 *   and the interior staircase all use the GLL 1D nodes (line_coord_pm_one), whose
 *   logarithmic Lebesgue constant keeps high-order interpolation well-conditioned.
 *   The named production layouts are unaffected, since GLL coincides with the
 *   equispaced layout at orders 1 and 2 (so Quad8/Hex8/Hex20 keep their exact
 *   public coordinates); the layout is this module's own convention only for
 *   order >= 3.
 * - **Modal basis.** The quadrilateral and hexahedral Vandermondes are assembled
 *   in a tensor **Legendre** basis rather than raw monomials. The serendipity
 *   exponent set is downward-closed, so the Legendre and monomial spans are
 *   identical (the change of basis is triangular) -- the nodal shape functions are
 *   unchanged -- but the Legendre Vandermonde is far better conditioned. (The
 *   fixed Wedge15 layout, order 2, keeps the monomial form; it is trivially
 *   well-conditioned.)
 */
class SerendipityBasis final : public BasisFunction {
public:
    /**
     * @brief Construct an arbitrary-order quadrilateral or hexahedral serendipity basis.
     *
     * @details This is the arbitrary-order entry point for the serendipity
     * families with a free order: the quadrilateral and the hexahedron. The
     * topology carries no node-count assumption; the serendipity polynomial
     * space, reference nodes (generated here in VTK-consistent stratified order),
     * and nodal coefficient table are built from the requested order (which must
     * be @f$p \ge 1@f$). Wedge serendipity is a single fixed layout and is not
     * constructed this way -- use the named ElementType overload (Wedge15).
     *
     * @param topology Must be BasisTopology::Quadrilateral or BasisTopology::Hexahedron.
     * @param order Polynomial order @f$p \ge 1@f$; orders below 1 are rejected.
     * @throws BasisConfigurationException If @p order is less than 1.
     * @throws BasisElementCompatibilityException If @p topology is not Quadrilateral or Hexahedron.
     */
    SerendipityBasis(BasisTopology topology, int order);

    /**
     * @brief Construct a serendipity basis from a named element layout.
     *
     * @details Convenience overload for the named, fixed serendipity layouts.
     * Each layout is the fixed-order instance of its family, built through the
     * same generated construction as the arbitrary-order path and taking its
     * nodes from ReferenceNodeLayout: Quad8 is the quadrilateral at order 2, Hex8
     * and Hex20 are the hexahedron at orders 1 and 2, and Wedge15 is the prism
     * layout. Each layout carries an inferred fixed order (Hex8 to 1; Quad8,
     * Hex20, and Wedge15 to 2); the requested @p order must equal that inferred
     * order and is never adjusted to fit, so a mismatched request (including
     * order 0 or negative) is rejected. Arbitrary-order quadrilateral and
     * hexahedral serendipity is requested through the BasisTopology overload.
     *
     * @param type Named serendipity element type (Quad8, Hex8, Hex20, or Wedge15).
     * @param order Requested order; must equal the layout's inferred fixed order
     *        (1 for Hex8; 2 for Quad8, Hex20, and Wedge15).
     * @throws BasisConfigurationException If @p order does not match the layout's inferred order.
     * @throws BasisElementCompatibilityException If the element type is unsupported.
     */
    SerendipityBasis(ElementType type, int order);

    /**
     * @brief Construct a serendipity basis from a named layout at its fixed order.
     *
     * @details Single-argument convenience overload for the named serendipity
     * layouts: the order is the one fixed by the layout (1 for Hex8; 2 for Quad8,
     * Hex20, and Wedge15), so the caller does not repeat it. Equivalent to
     * SerendipityBasis(type, <fixed order>).
     *
     * @param type Named serendipity element type (Quad8, Hex8, Hex20, or Wedge15).
     * @throws BasisElementCompatibilityException If the element type is unsupported.
     */
    explicit SerendipityBasis(ElementType type);

    /** @copydoc BasisFunction::basis_type() */
    BasisType basis_type() const noexcept final { return BasisType::Serendipity; }

    /** @copydoc BasisFunction::topology() */
    BasisTopology topology() const noexcept final { return topology_; }

    /** @copydoc BasisFunction::dimension() */
    int dimension() const noexcept final { return dimension_; }

    /** @copydoc BasisFunction::order() */
    int order() const noexcept final { return order_; }

    /** @copydoc BasisFunction::size() */
    std::size_t size() const noexcept final { return size_; }

    /**
     * @brief Return the reference interpolation nodes in basis ordering.
     *
     * @details Node coordinates are the points at which the serendipity basis
     * satisfies the nodal interpolation property. All families take their nodes
     * from ReferenceNodeLayout, the public node-ordering source the solver adapter
     * permutes against: the fixed Wedge15 layout and the quadrilateral/hexahedral
     * families (named or arbitrary-order) alike, in VTK-consistent stratified
     * order -- corners and edges first (matching the public Quad8/Hex8/Hex20
     * ordering at the named orders), then the face and volume interior points
     * needed to make the reduced polynomial space unisolvent. For @f$p \ge 3@f$
     * that interior ordering is an implementation convention; callers should pair
     * it with basis values from the same object rather than assume an external
     * mesh ordering contract beyond the supported named production layouts.
     *
     * @return Reference node coordinates, one per basis function.
     */
    const std::vector<math::Vector<double, 3>>& nodes() const noexcept final { return nodes_; }

    /**
     * @brief Evaluate serendipity basis values into caller-provided storage.
     * @param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
     * @param values_out Output span with at least size() entries.
     */
    void evaluate_values_to(const math::Vector<double, 3>& xi,
                            std::span<double> values_out) const final;

    /**
     * @brief Evaluate serendipity basis gradients into caller-provided storage.
     * @param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
     * @param gradients_out Output span with at least size() entries.
     */
    void evaluate_gradients_to(const math::Vector<double, 3>& xi,
                               std::span<Gradient> gradients_out) const final;

    /**
     * @brief Evaluate serendipity basis Hessians into caller-provided storage.
     * @param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
     * @param hessians_out Output span with at least size() entries.
     */
    void evaluate_hessians_to(const math::Vector<double, 3>& xi,
                              std::span<Hessian> hessians_out) const final;

private:
    BasisTopology topology_{BasisTopology::Unknown};
    int dimension_{0};
    int order_{0};
    std::size_t size_{0};
    std::vector<math::Vector<double, 3>> nodes_;
    // Per-axis degrees (a, b, c) of the tensor modes spanning the family's
    // polynomial space. Interpreted as monomial powers r^a s^b t^c or, when
    // uses_legendre_modes_ is set, as tensor Legendre degrees P_a(r) P_b(s) P_c(t)
    // (the same space; see ModalAxisKind in SerendipityBasis.cpp).
    std::vector<std::array<int, 3>> mode_exponents_;
    // Row-major inverse (generalized) Vandermonde, indexed as [mode, basis].
    std::vector<double> inv_vandermonde_;
    // Whether the tensor modes are Legendre polynomials (quadrilateral/hexahedral
    // families) or plain monomials (the fixed Wedge15 layout). Evaluation must use
    // the same family the coefficient table was built with.
    bool uses_legendre_modes_{false};

    // Build the quadrilateral serendipity mode set, nodes, and Legendre
    // coefficient table for the given order. (Details at the definition.)
    void init_quadrilateral(int order);
    // Build the hexahedral serendipity mode set, nodes, and Legendre coefficient
    // table for the given order; Hex8/Hex20 are its order-1/order-2 instances.
    void init_hexahedron(int order);
    // Build the fixed Wedge15 layout from its tabulated monomial mode space.
    void init_fixed_named(ElementType type);

    void evaluate_all_to(const math::Vector<double, 3>& xi,
                         std::span<double> values_out,
                         std::span<Gradient> gradients_out,
                         std::span<Hessian> hessians_out) const override;
};

/** @} */

} // namespace svmp::FE::basis

#endif // SVMP_FE_BASIS_SERENDIPITYBASIS_H
