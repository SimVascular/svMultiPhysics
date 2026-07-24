// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_QUADRATURE_RULE_H
#define SVMP_FE_QUADRATURE_RULE_H

/**
 * @file QuadratureRule.h
 * @brief Immutable reference-space quadrature rule contract.
 * @ingroup FE_Quadrature
 *
 * This header defines the consumer-facing representation of a finite-element
 * quadrature rule. Rule construction and validation are implemented separately
 * so consumers depend only on the stable query interface.
 */

/**
 * @defgroup FE_Quadrature Quadrature
 * @ingroup FE
 * @brief Immutable integration rules on canonical finite-element reference cells.
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
 * The zeroth moment is the analytic normalization
 * @f[
 *   M_0 = \int_{\hat K} 1\,d\hat x = \sum_q w_q.
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
 * API. Integration consumers read rule metadata and ordered samples; they do
 * not derive new rules or modify rule storage.
 *
 * Complete-data construction, reference-cell traits, point-containment checks,
 * exact weight summation, concrete generators, caches, and rule-selection
 * facilities are module implementation details.
 *
 * ## Rule-provider contract
 *
 * Concrete rule providers are the only supported subclasses. A provider builds
 * one complete RuleData payload before invoking the protected constructor and
 * exposes no mutation afterward. It must advertise only exactness established
 * for every rule it supplies through analytic moment tests. Derivation is a
 * provider extension seam, not an integration-consumer customization point.
 * The protected constructor is the only provider extension point; public
 * metadata and sample queries are fixed, nonvirtual operations on validated
 * base-class state.
 *
 * ## Rule contract
 *
 * Construction is the sole validity boundary: a rule is complete and
 * structurally valid when its constructor returns, and consumers do not perform
 * a separate revalidation step. The constructor rejects unsupported cells,
 * negative exactness, empty or mismatched storage, non-finite coordinates or
 * weights, points outside the declared reference cell, and weights whose sum
 * does not equal the canonical rule's zeroth moment within the scaled moment
 * tolerance. The sum of the stored binary64 weights is evaluated exactly and
 * independently of their order.
 *
 * Structural validation does not require unique points or nonzero, positive
 * individual weights. It verifies metadata, containment, finiteness, and the
 * zeroth moment; it does not prove higher-order polynomial moments. A
 * polynomial exactness of @f$p@f$ guarantees every polynomial of total degree
 * at most @f$p@f$. A rule can integrate selected higher-degree polynomials
 * without increasing that common guarantee. Every concrete provider is
 * responsible for establishing its advertised exactness with analytic moment
 * tests.
 *
 * Points use one fixed-size three-component representation. Providers initialize
 * all three coordinates explicitly because the Eigen-backed vector is not
 * zero-initialized by default. Only the first dimension() components
 * are active, and every inactive component is zero within the coordinate
 * tolerance. The supported canonical domains are:
 *
 * | Cell family | Canonical reference domain | Zeroth moment |
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
 * ## Ownership and lifetime
 *
 * Rule objects are non-copyable and non-movable. They are intended to be built
 * once, retained through a const owning handle, and shared across integrations.
 * References returned by points() and weights() remain valid for the lifetime
 * of the rule. Consumers should retain that rule rather than copying its point
 * count or samples into parallel authoritative state. Concurrent const access
 * is safe while an owning handle keeps the rule alive.
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
 * Only the first QuadratureRule::dimension() components are active. Providers
 * explicitly zero remaining components, giving point, line, surface, and volume
 * rules a uniform representation directly compatible with FE math consumers.
 */
using QuadPoint = math::Vector<double, 3>;

/**
 * @brief Immutable consumer interface for a quadrature rule on a reference cell.
 *
 * Concrete rule providers initialize the base with one complete RuleData
 * payload. The payload is validated before construction returns, and the object
 * exposes no mutation or assignment path afterward. General solver consumers
 * use only the public const query interface. Providers supply data through the
 * protected constructor; they do not override the fixed metadata or sample
 * queries. Successful construction establishes the rule invariant; no public
 * revalidation operation is required or exposed. Concurrent const access is
 * safe while an owning handle keeps the rule alive.
 */
class QuadratureRule {
public:
    /** @brief Destroy a quadrature rule through the abstract interface. */
    virtual ~QuadratureRule() = 0;

    /** @brief Rule objects cannot be copied through the abstract interface. */
    QuadratureRule(const QuadratureRule&) = delete;

    /** @brief Rule objects cannot be moved through the abstract interface. */
    QuadratureRule(QuadratureRule&&) = delete;

    /** @brief Rule objects cannot be replaced through base assignment. */
    QuadratureRule& operator=(const QuadratureRule&) = delete;

    /** @brief Rule objects cannot be replaced through base move assignment. */
    QuadratureRule& operator=(QuadratureRule&&) = delete;

    /**
     * @brief Return the number of ordered point/weight pairs.
     * @return Quadrature sample count.
     */
    std::size_t num_points() const noexcept { return points_.size(); }

    /**
     * @brief Return the total-degree polynomial exactness declared by the rule.
     *
     * A value @f$p@f$ guarantees exact integration of every polynomial with
     * total degree at most @f$p@f$. The rule can also integrate selected
     * higher-degree polynomials. Structural validation does not independently
     * prove this guarantee; concrete providers establish it through analytic
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
    int dimension() const noexcept { return dimension_; }

    /**
     * @brief Return the canonical reference-cell family.
     * @return Reference topology integrated by this rule.
     */
    svmp::CellFamily cell_family() const noexcept { return cell_family_; }

    /**
     * @brief Return one reference coordinate without bounds checking.
     * @param i Point index in the half-open range `[0, num_points())`.
     * @return Immutable reference to the indexed quadrature point, valid for
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
     * @return Immutable point storage, valid for the lifetime of this rule.
     */
    const std::vector<QuadPoint>& points() const noexcept { return points_; }

    /**
     * @brief Return all reference weights in point order.
     * @return Immutable weight storage, valid for the lifetime of this rule.
     */
    const std::vector<double>& weights() const noexcept { return weights_; }

    /**
     * @brief Return the rule's zeroth moment.
     *
     * This is the integral of the constant function one. All supported rules
     * are unweighted rules on complete canonical reference cells, so the
     * constructor derives this value from cell_family() and it equals the
     * geometric measure of that cell.
     * @f[
     *   M_0 = \int_{\hat K} 1\,d\hat x = \sum_q w_q.
     * @f]
     *
     * @return Expected sum of the quadrature weights.
     */
    double zeroth_moment() const noexcept { return zeroth_moment_; }

protected:
    /**
     * @brief Complete construction payload used by concrete rule providers.
     *
     * Integration consumers do not construct this payload. A provider computes
     * all three fields before invoking the protected base constructor and must
     * substantiate its declared exactness with analytic moment tests.
     */
    struct RuleData {
        int polynomial_exactness{-1};    ///< Declared total-degree polynomial exactness.
        std::vector<QuadPoint> points;   ///< Ordered canonical reference coordinates.
        std::vector<double> weights;     ///< Weights paired with points in the same order.
    };

    /**
     * @brief Construct and validate one complete immutable rule.
     *
     * Dimension and zeroth moment are derived from @p family; callers cannot
     * supply redundant topology metadata.
     *
     * @param family Supported canonical reference-cell family.
     * @param data Complete exactness, point, and weight payload.
     * @note Duplicate points and zero or negative weights remain admissible
     * when every other rule invariant is satisfied.
     * @throws InvalidArgumentException If the family is unsupported, exactness
     * is negative, storage is empty or mismatched, a value is non-finite, a point
     * is outside the reference cell, the weights do not reproduce its zeroth
     * moment within the scaled moment tolerance.
     */
    explicit QuadratureRule(svmp::CellFamily family, RuleData data);

private:
    /** @brief Fully checked state used by the delegating constructor. */
    struct ValidatedState {
        svmp::CellFamily cell_family;
        int dimension;
        int polynomial_exactness;
        double zeroth_moment;
        std::vector<QuadPoint> points;
        std::vector<double> weights;
    };

    /** @brief Validate a local payload before initializing immutable members. */
    static ValidatedState validate(svmp::CellFamily family, RuleData data);

    /** @brief Initialize members from state already checked by validate(). */
    explicit QuadratureRule(ValidatedState state);

    const svmp::CellFamily cell_family_;     ///< Canonical reference topology.
    const int dimension_;                    ///< Number of active coordinate components.
    const int polynomial_exactness_;         ///< Exactness declared by the concrete generator.
    const double zeroth_moment_;              ///< Canonical reference-cell measure.
    const std::vector<QuadPoint> points_;      ///< Ordered immutable reference coordinates.
    const std::vector<double> weights_;        ///< Immutable weights paired with points_.
};

/** @} */

} // namespace svmp::FE::quadrature

#endif // SVMP_FE_QUADRATURE_RULE_H
