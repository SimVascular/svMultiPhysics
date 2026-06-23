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

namespace svmp {
namespace FE {
namespace basis {

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
 * @details SerendipityBasis implements nodal bases for Quad4/Quad8,
 * Hex8/Hex20, and Wedge15. Compared with a complete tensor-product Lagrange
 * basis of the same nominal order, a serendipity basis removes selected
 * interior modes while retaining nodal interpolation on the supported node
 * layout.
 *
 * Quadrilateral serendipity bases are built from monomials
 * @f$x^{a_x}y^{a_y}@f$ whose superlinear degree is at most the requested
 * order. In this implementation the superlinear degree is
 * @f[
 *   sldeg(x^{a_x}y^{a_y}) =
 *   \begin{cases} a_x, & a_x > 1 \\ 0, & a_x \le 1 \end{cases}
 *   +
 *   \begin{cases} a_y, & a_y > 1 \\ 0, & a_y \le 1 \end{cases}.
 * @f]
 * The nodal basis is recovered by inverting the Vandermonde interpolation
 * matrix at the selected reference nodes. Values, gradients, and Hessians are
 * then evaluated by differentiating the monomial vector and applying the
 * inverse Vandermonde coefficients.
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
 * `SerendipityBasis(ElementType::Quad4, p)` supports explicit
 * arbitrary-order quadrilateral serendipity requests for @f$p \ge 1@f$
 * (requests below one are normalized to one). `ElementType::Quad8` remains
 * the standard quadratic eight-node layout and is valid only with order 2.
 * Solver-default basis selection remains separate: `basis_factory` maps the
 * complete Quad4 layout to the default linear Lagrange basis and maps Quad8 to
 * quadratic serendipity unless a caller explicitly requests a different
 * supported basis.
 *
 * Hex8 uses the standard trilinear corner basis
 * @f$(1 \pm r)(1 \pm s)(1 \pm t)/8@f$. Hex20 and Wedge15 use tabulated
 * polynomial coefficient tables over monomial bases; analytical gradients and
 * Hessians are obtained by differentiating those monomials. Hex20 evaluation
 * is reordered through ReferenceNodeLayout so the output matches the public
 * basis ordering.
 */
class SerendipityBasis final : public BasisFunction {
public:
    /**
     * @brief Construct a serendipity basis for an element type and polynomial order.
     *
     * @details The constructor selects the topology-specific interpolation
     * space, computes the reference node coordinates, and initializes any
     * coefficient tables needed for evaluation. Quadrilateral bases build and
     * invert a Vandermonde matrix for the selected serendipity monomials.
     * Hex20 and Wedge15 use fixed coefficient tables. For hexahedra, only
     * linear Hex8 and quadratic Hex20 serendipity spaces are supported. For
     * wedges, only quadratic Wedge15 is supported. Quad4 supports explicit
     * quadrilateral serendipity requests of any order @f$p \ge 1@f$; Quad8 is
     * restricted to order 2.
     *
     * @param type Element type used to determine topology and reference-node layout.
     * @param order Requested polynomial order.
     * @throws BasisConfigurationException If the requested order is invalid.
     * @throws BasisElementCompatibilityException If the element type is unsupported.
     */
    SerendipityBasis(ElementType type, int order);

    /** @copydoc BasisFunction::basis_type() */
    BasisType basis_type() const noexcept final { return BasisType::Serendipity; }

    /** @copydoc BasisFunction::element_type() */
    ElementType element_type() const noexcept final { return element_type_; }

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
     * satisfies the nodal interpolation property. Quadrilateral nodes are
     * placed first on the boundary and then, for higher order requests, at the
     * selected interior points needed to make the reduced monomial space
     * unisolvent. Hexahedral and wedge nodes are taken from
     * ReferenceNodeLayout. For high-order Quad4 serendipity, the deterministic
     * interior row ordering is an implementation convention; callers should
     * pair it with basis values from the same object rather than assume an
     * external mesh ordering contract beyond the supported Quad4/Quad8
     * production layouts.
     *
     * @return Reference node coordinates, one per basis function.
     */
    const std::vector<math::Vector<double, 3>>& nodes() const noexcept final { return nodes_; }

    /**
     * @brief Evaluate serendipity basis function values at a reference coordinate.
     *
     * @details For quadrilateral bases, this evaluates the serendipity
     * monomial vector and multiplies by the inverse Vandermonde matrix to
     * obtain nodal shape-function values. For Hex8, values are the standard
     * trilinear corner products. For Hex20 and Wedge15, values are evaluated
     * from the stored polynomial coefficient tables.
     *
     * @param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
     * @param values Receives one value per basis function.
     */
    void evaluate_values(const math::Vector<double, 3>& xi,
                         std::vector<double>& values) const final;

    /**
     * @brief Evaluate analytical serendipity basis gradients at a reference coordinate.
     *
     * @details Gradients are derivatives with respect to reference
     * coordinates. Quadrilateral gradients differentiate the monomial vector
     * before applying the inverse Vandermonde coefficients. Hex8 gradients are
     * direct derivatives of the trilinear corner products. Hex20 and Wedge15
     * gradients are computed by differentiating the tabulated monomial
     * expansions.
     *
     * @param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
     * @param gradients Receives one three-component gradient per basis function.
     */
    void evaluate_gradients(const math::Vector<double, 3>& xi,
                            std::vector<Gradient>& gradients) const final;

    /**
     * @brief Evaluate analytical serendipity basis Hessians at a reference coordinate.
     *
     * @details Hessians are second derivatives in reference coordinates and
     * are stored as 3-by-3 matrices. Quadrilateral Hessians use second
     * derivatives of the monomial vector and inverse Vandermonde coefficients.
     * Hex8 Hessians are computed directly from the trilinear corner products.
     * Hex20 and Wedge15 Hessians are computed by differentiating their
     * polynomial coefficient tables twice.
     *
     * @param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
     * @param hessians Receives one 3-by-3 Hessian per basis function.
     */
    void evaluate_hessians(const math::Vector<double, 3>& xi,
                           std::vector<Hessian>& hessians) const final;

    /**
     * @brief Evaluate serendipity values, gradients, and Hessians together.
     *
     * @details This vector API is backed by the same span-based evaluator as
     * the assembly-oriented `*_to` methods, so topology-specific polynomial
     * setup can be shared for a quadrature point.
     *
     * @param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
     * @param values Receives one value per basis function.
     * @param gradients Receives one three-component gradient per basis function.
     * @param hessians Receives one 3-by-3 Hessian per basis function.
     */
    void evaluate_all(const math::Vector<double, 3>& xi,
                      std::vector<double>& values,
                      std::vector<Gradient>& gradients,
                      std::vector<Hessian>& hessians) const final;

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
    ElementType element_type_;
    int dimension_;
    int order_;
    std::size_t size_;
    std::vector<math::Vector<double, 3>> nodes_;
    std::vector<std::array<int, 2>> quad_monomial_exponents_;
    // Row-major inverse Vandermonde, indexed as [monomial, basis].
    std::vector<double> quad_inv_vandermonde_;

    void evaluate_all_to(const math::Vector<double, 3>& xi,
                         std::span<double> values_out,
                         std::span<Gradient> gradients_out,
                         std::span<Hessian> hessians_out) const;
};

/** @} */

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_SERENDIPITYBASIS_H
