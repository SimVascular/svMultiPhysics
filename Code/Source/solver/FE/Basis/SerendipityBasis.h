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

namespace svmp {
namespace FE {
namespace basis {

/// \defgroup FE_SerendipityBasis SerendipityBasis
/// \ingroup FE_Basis
/// \brief Construction and evaluation API for reduced serendipity finite-element bases.
///
/// \details This group documents reduced degree-of-freedom basis families that
/// preserve nodal interpolation on supported element boundaries while omitting
/// selected interior tensor-product modes. These bases are used for standard
/// serendipity elements and geometry-mode mappings that intentionally use a
/// lower-order interpolation space.
/// @{

/// \brief Reduced-degree-of-freedom serendipity basis on supported reference elements.
///
/// \details SerendipityBasis implements nodal bases for Quad4/Quad8,
/// Hex8/Hex20, and Wedge15. Compared with a complete tensor-product Lagrange
/// basis of the same nominal order, a serendipity basis removes selected
/// interior modes while retaining nodal interpolation on the supported node
/// layout.
///
/// Quadrilateral serendipity bases are built from monomials
/// \f$x^{a_x}y^{a_y}\f$ whose superlinear degree is at most the requested
/// order. In this implementation the superlinear degree is
/// \f[
///   sldeg(x^{a_x}y^{a_y}) =
///   \begin{cases} a_x, & a_x > 1 \\ 0, & a_x \le 1 \end{cases}
///   +
///   \begin{cases} a_y, & a_y > 1 \\ 0, & a_y \le 1 \end{cases}.
/// \f]
/// The nodal basis is recovered by inverting the Vandermonde interpolation
/// matrix at the selected reference nodes. Values, gradients, and Hessians are
/// then evaluated by differentiating the monomial vector and applying the
/// inverse Vandermonde coefficients.
///
/// Hex8 uses the standard trilinear corner basis
/// \f$(1 \pm r)(1 \pm s)(1 \pm t)/8\f$. Hex20 and Wedge15 use tabulated
/// polynomial coefficient tables over monomial bases; analytical gradients and
/// Hessians are obtained by differentiating those monomials. Hex20 evaluation
/// is reordered through ReferenceNodeLayout so the output matches the public
/// basis ordering.
///
/// When `geometry_mode` is enabled for Hex20, the basis uses the trilinear
/// Hex8 corner functions for geometry mapping and assigns zero contribution to
/// the quadratic edge nodes. This preserves the public Hex20 node count while
/// intentionally reducing the geometry interpolation order.
class SerendipityBasis : public BasisFunction {
public:
    /// \brief Construct a serendipity basis for an element type and polynomial order.
    ///
    /// \details The constructor selects the topology-specific interpolation
    /// space, computes the reference node coordinates, and initializes any
    /// coefficient tables needed for evaluation. Quadrilateral bases build and
    /// invert a Vandermonde matrix for the selected serendipity monomials.
    /// Hex20 and Wedge15 use fixed coefficient tables. For hexahedra, only
    /// linear Hex8 and quadratic Hex20 serendipity spaces are supported. For
    /// wedges, only quadratic Wedge15 is supported.
    ///
    /// \param type Element type used to determine topology and reference-node layout.
    /// \param order Requested polynomial order.
    /// \param geometry_mode When true, allow reduced geometry-mapping behavior for supported elements.
    /// \throws BasisConfigurationException If the requested order or mode is invalid.
    /// \throws BasisElementCompatibilityException If the element type is unsupported.
    SerendipityBasis(ElementType type, int order, bool geometry_mode = false);

    /// \copydoc BasisFunction::basis_type()
    BasisType basis_type() const noexcept override { return BasisType::Serendipity; }

    /// \copydoc BasisFunction::element_type()
    ElementType element_type() const noexcept override { return element_type_; }

    /// \copydoc BasisFunction::dimension()
    int dimension() const noexcept override { return dimension_; }

    /// \copydoc BasisFunction::order()
    int order() const noexcept override { return order_; }

    /// \copydoc BasisFunction::size()
    std::size_t size() const noexcept override { return size_; }

    /// \brief Return the reference interpolation nodes in basis ordering.
    ///
    /// \details Node coordinates are the points at which the serendipity basis
    /// satisfies the nodal interpolation property. Quadrilateral nodes are
    /// placed first on the boundary and then, for higher order requests, at the
    /// selected interior points needed to make the reduced monomial space
    /// unisolvent. Hexahedral and wedge nodes are taken from
    /// ReferenceNodeLayout.
    ///
    /// \return Reference node coordinates, one per basis function.
    const std::vector<math::Vector<Real, 3>>& nodes() const noexcept { return nodes_; }

    /// \brief Evaluate serendipity basis function values at a reference coordinate.
    ///
    /// \details For quadrilateral bases, this evaluates the serendipity
    /// monomial vector and multiplies by the inverse Vandermonde matrix to
    /// obtain nodal shape-function values. For Hex8, values are the standard
    /// trilinear corner products. For Hex20 and Wedge15, values are evaluated
    /// from the stored polynomial coefficient tables. In Hex20 geometry mode,
    /// only the first eight corner values are nonzero and they match Hex8.
    ///
    /// \param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
    /// \param values Receives one value per basis function.
    void evaluate_values(const math::Vector<Real, 3>& xi,
                         std::vector<Real>& values) const override;

    /// \brief Evaluate analytical serendipity basis gradients at a reference coordinate.
    ///
    /// \details Gradients are derivatives with respect to reference
    /// coordinates. Quadrilateral gradients differentiate the monomial vector
    /// before applying the inverse Vandermonde coefficients. Hex8 gradients are
    /// direct derivatives of the trilinear corner products. Hex20 and Wedge15
    /// gradients are computed by differentiating the tabulated monomial
    /// expansions. In Hex20 geometry mode, edge-node gradients are zero and the
    /// corner gradients match Hex8.
    ///
    /// \param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
    /// \param gradients Receives one three-component gradient per basis function.
    void evaluate_gradients(const math::Vector<Real, 3>& xi,
                            std::vector<Gradient>& gradients) const override;

    /// \brief Evaluate analytical serendipity basis Hessians at a reference coordinate.
    ///
    /// \details Hessians are second derivatives in reference coordinates and
    /// are stored as 3-by-3 matrices. Quadrilateral Hessians use second
    /// derivatives of the monomial vector and inverse Vandermonde coefficients.
    /// Hex8 Hessians are delegated to the linear Lagrange Hex8 basis. Hex20 and
    /// Wedge15 Hessians are computed by differentiating their polynomial
    /// coefficient tables twice. In Hex20 geometry mode, only the corner
    /// Hessians from the Hex8 geometry mapping are populated.
    ///
    /// \param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
    /// \param hessians Receives one 3-by-3 Hessian per basis function.
    void evaluate_hessians(const math::Vector<Real, 3>& xi,
                           std::vector<Hessian>& hessians) const override;

private:
    ElementType element_type_;
    int dimension_;
    int order_;
    std::size_t size_;
    std::vector<math::Vector<Real, 3>> nodes_;
    std::vector<std::array<int, 2>> quad_monomial_exponents_;
    // Row-major inverse Vandermonde, indexed as [monomial, basis].
    std::vector<Real> quad_inv_vandermonde_;

    // When true, this basis is used purely for geometry mapping and may use
    // reduced polynomial order (e.g., Hex20 geometry as Hex8).
    bool geometry_mode_;
};

/// @}

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_SERENDIPITYBASIS_H
