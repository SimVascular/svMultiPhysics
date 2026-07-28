// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_QUADRATURE_RULE_H
#define SVMP_FE_QUADRATURE_RULE_H

/**
 * @file QuadratureRule.h
 * @brief Validated reference-space quadrature rule value type.
 * @ingroup FE_Quadrature
 *
 * This header defines an owning quadrature-rule value with read-only access to
 * its metadata, points, and weights.
 */

/**
 * @defgroup FE_Quadrature Quadrature
 * @ingroup FE
 * @brief Validated integration rules on canonical finite-element reference cells.
 *
 * @details
 * ## Scope
 *
 * The Quadrature module owns ordered reference coordinates and weights used to
 * approximate an unweighted integral on a canonical reference cell:
 * @f[
 *   \int_{\hat K} f(\hat x)\,d\hat x
 *   \approx \sum_q w_q f(\hat x_q).
 * @f]
 * The reference-cell measure satisfies
 * @f[
 *   |\hat K| = \int_{\hat K} 1\,d\hat x = \sum_q w_q.
 * @f]
 * A rule identifies its reference-cell family, reports its dimension and
 * declared polynomial exactness, and keeps every point paired with its
 * corresponding weight.
 *
 * The module does not choose the exactness required by an equation term, apply
 * reduced-integration policy, select a basis, own mesh storage, embed or orient
 * a face in its parent element, or map a reference integral into physical
 * space. Those operations require solver, basis, mesh, and geometry context and
 * remain at the caller's integration boundary.
 *
 * ## Public API
 *
 * @ref svmp::FE::quadrature::QuadPoint "QuadPoint" and the const query surface
 * of @ref svmp::FE::quadrature::QuadratureRule "QuadratureRule" form the public
 * API. Integration consumers read rule metadata, points, and weights; they do
 * not derive new rules or modify rule storage.
 *
 * Reference-cell metadata, point-containment checks, exact weight summation,
 * concrete generators, caches, and rule-selection facilities are module
 * implementation details.
 *
 * ## Constructing rules
 *
 * Construct a QuadratureRule from its cell family, declared polynomial
 * exactness, ordered points, and paired weights. Rule-generating functions can
 * compute those values and return the resulting rule. A generator must
 * advertise only exactness established through analytic moment tests.
 *
 * ## Rule contract
 *
 * Construction is the sole validity boundary: a rule is complete and
 * structurally valid when its constructor returns, and consumers do not perform
 * a separate revalidation step. The constructor rejects unsupported cells,
 * negative exactness, empty or mismatched storage, non-finite coordinates or
 * weights, points outside the declared reference cell, and weights whose sum
 * does not equal the canonical reference-cell measure within the scaled
 * measure tolerance. The sum of the stored binary64 weights is evaluated
 * exactly and independently of their order.
 *
 * Structural validation does not require unique points or nonzero, positive
 * individual weights. It verifies metadata, containment, finiteness, and the
 * reference-cell measure; it does not prove higher-order polynomial moments. A
 * polynomial exactness of @f$p@f$ guarantees every polynomial of total degree
 * at most @f$p@f$. A rule can integrate selected higher-degree polynomials
 * without increasing that common guarantee. Rule generators are responsible
 * for establishing their advertised exactness with analytic moment tests.
 *
 * Points use one fixed-size three-component representation. Generators initialize
 * all three coordinates explicitly because the Eigen-backed vector is not
 * zero-initialized by default. Only the first dimension() components
 * are active, and every inactive component is zero within the coordinate
 * tolerance. The supported canonical domains are:
 *
 * | Cell family | Canonical reference domain | Reference-cell measure |
 * | ----------- | -------------------------- | ------------- |
 * | Point | @f$(0,0,0)@f$ | @f$1@f$ |
 * | Line | @f$[-1,1]@f$ | @f$2@f$ |
 * | Triangle | @f$\xi,\eta\geq0;\ \xi+\eta\leq1@f$ | @f$1/2@f$ |
 * | Quad | @f$[-1,1]^2@f$ | @f$4@f$ |
 * | Tetra | @f$\xi,\eta,\zeta\geq0;\ \xi+\eta+\zeta\leq1@f$ | @f$1/6@f$ |
 * | Hex | @f$[-1,1]^3@f$ | @f$8@f$ |
 * | Wedge | reference triangle @f$\times[-1,1]@f$ | @f$1@f$ |
 *
 * Pyramid, Polygon, Polyhedron, and unknown cell families are intentionally
 * unsupported.
 *
 * ## Value semantics and lifetime
 *
 * Rules can be copied, moved, and assigned. Each value owns its points and
 * weights. References returned by points() and weights() remain valid until
 * that rule is assigned or destroyed. Concurrent const access is safe while
 * no thread assigns to the same rule.
 */

#include "FE/Common/Types.h"
#include "FE/Math/Vector.h"

#include <cstddef>
#include <vector>

namespace svmp::FE::quadrature {

/** @addtogroup FE_Quadrature
 * @{
 */

/**
 * @brief Three-component coordinate used for every reference quadrature point.
 *
 * Only the first QuadratureRule::dimension() components are active. Generators
 * explicitly zero remaining components, giving point, line, surface, and volume
 * rules a uniform representation directly compatible with FE math consumers.
 */
using QuadPoint = math::Vector<double, 3>;

/**
 * @brief Validated quadrature rule on a canonical reference cell.
 *
 * Construct the rule from its family, declared exactness, ordered points, and
 * paired weights. Construction validates the supplied data. Public queries do
 * not permit editing the owned point and weight storage.
 */
class QuadratureRule final {
public:
    /**
     * @brief Construct and validate one complete rule.
     *
     * Dimension and reference-cell measure are derived from @p family; callers
     * cannot supply redundant topology metadata.
     *
     * @param family Supported canonical reference-cell family.
     * @param polynomial_exactness Declared total-degree polynomial exactness.
     * @param points Ordered canonical reference coordinates.
     * @param weights Weights paired with @p points in the same order.
     * @note Duplicate points and zero or negative weights remain admissible
     * when every other rule invariant is satisfied.
     * @throws InvalidArgumentException If the family is unsupported, exactness
     * is negative, storage is empty or mismatched, a value is non-finite, a point
     * is outside the reference cell, or the weights do not reproduce the
     * reference-cell measure within the scaled measure tolerance.
     */
    explicit QuadratureRule(
        svmp::CellFamily family,
        int polynomial_exactness,
        std::vector<QuadPoint> points,
        std::vector<double> weights);

    /**
     * @brief Return the number of ordered point/weight pairs.
     * @return Quadrature point count.
     */
    std::size_t num_points() const noexcept { return points_.size(); }

    /**
     * @brief Return the total-degree polynomial exactness declared by the rule.
     *
     * A value @f$p@f$ guarantees exact integration of every polynomial with
     * total degree at most @f$p@f$. The rule can also integrate selected
     * higher-degree polynomials. Structural validation does not independently
     * prove this guarantee; rule generators establish it through analytic
     * moment tests.
     *
     * @return Declared total degree that a conforming rule integrates exactly.
     */
    int polynomial_exactness() const noexcept { return polynomial_exactness_; }

    /**
     * @brief Return the dimension of the canonical reference cell.
     *
     * This is also the number of active components in each QuadPoint.
     *
     * @return Reference dimension, from zero for Point through three for volume cells.
     */
    int dimension() const noexcept;

    /**
     * @brief Return the canonical reference-cell family.
     * @return Reference topology integrated by this rule.
     */
    svmp::CellFamily cell_family() const noexcept { return cell_family_; }

    /**
     * @brief Return one reference coordinate without bounds checking.
     * @param i Point index in the half-open range `[0, num_points())`.
     * @return Const reference to the indexed quadrature point, valid for
     * the lifetime of this rule.
     * @pre @p i is less than num_points().
     */
    const QuadPoint& point(std::size_t i) const noexcept { return points_[i]; }

    /**
     * @brief Return one reference weight without bounds checking.
     * @param i Weight index in the half-open range `[0, num_points())`.
     * @return Weight paired with point(@p i).
     * @pre @p i is less than num_points().
     */
    double weight(std::size_t i) const noexcept { return weights_[i]; }

    /**
     * @brief Return all reference coordinates in integration order.
     * @return Read-only point storage, valid for the lifetime of this rule.
     */
    const std::vector<QuadPoint>& points() const noexcept { return points_; }

    /**
     * @brief Return all reference weights in point order.
     * @return Read-only weight storage, valid for the lifetime of this rule.
     */
    const std::vector<double>& weights() const noexcept { return weights_; }

    /**
     * @brief Return the measure of the canonical reference cell.
     *
     * This is the integral of the constant function one. All supported rules
     * are unweighted rules on complete canonical reference cells, so the
     * constructor derives this value from cell_family() and it equals the
     * geometric measure of that cell.
     * @f[
     *   |\hat K| = \int_{\hat K} 1\,d\hat x = \sum_q w_q.
     * @f]
     *
     * @return Geometric measure of the canonical reference cell.
     */
    double reference_cell_measure() const noexcept { return reference_cell_measure_; }

private:
    svmp::CellFamily cell_family_;          ///< Canonical reference topology.
    int polynomial_exactness_;              ///< Exactness declared by the generator.
    double reference_cell_measure_;         ///< Canonical reference-cell measure.
    std::vector<QuadPoint> points_;          ///< Ordered reference coordinates.
    std::vector<double> weights_;            ///< Weights paired with points_.
};

/** @} */

} // namespace svmp::FE::quadrature

#endif // SVMP_FE_QUADRATURE_RULE_H
