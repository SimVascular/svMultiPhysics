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
 * approximate an integral on a canonical reference cell:
 * @f[
 *   \int_{\hat K} f(\hat x)\,d\hat x
 *   \approx \sum_q w_q f(\hat x_q).
 * @f]
 * A rule identifies its reference-cell family, reports its intrinsic dimension
 * and declared polynomial exactness, and keeps every point paired with its
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
 * compensated weight summation, concrete generators, caches, and rule-selection
 * facilities are module implementation details.
 *
 * ## Rule-provider contract
 *
 * Concrete rule providers are the only supported subclasses. A provider builds
 * one complete RuleData payload before invoking the protected constructor and
 * exposes no mutation afterward. It must advertise only exactness established
 * for every rule it supplies through analytic moment tests. Derivation is a
 * provider extension seam, not an integration-consumer customization point.
 *
 * ## Rule contract
 *
 * A rule is complete and structurally valid when construction returns. The
 * constructor rejects unsupported cells, negative exactness, empty or mismatched
 * storage, non-finite coordinates or weights, points outside the declared
 * reference cell, and weights whose zeroth moment does not equal the reference
 * measure in both compensated arithmetic and ordinary stored-order double
 * accumulation. A condition estimate based on the absolute weight sum rejects
 * signed rules whose cancellation is too sensitive for double precision.
 * Negative individual weights remain valid because some quadrature families
 * require them.
 *
 * polynomial_exactness() is declared by the concrete construction algorithm.
 * Structural validation verifies metadata, containment, finiteness, and the
 * zeroth moment; it does not prove higher-order polynomial moments. Every
 * concrete rule provider is responsible for making its declaration a guarantee
 * by establishing the advertised exactness with analytic moment tests.
 *
 * Points use one zero-initialized three-component representation. Constructors
 * that receive fewer than three coordinates zero-fill the omitted components.
 * Only the first dimension() components are active, and every inactive
 * component is zero within the construction tolerance. The supported canonical
 * domains are:
 *
 * | Cell family | Canonical reference domain | Measure |
 * | ----------- | -------------------------- | ------- |
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
 * count or samples into parallel authoritative state.
 */

#include "FE/Common/Types.h"
#include "Math/Vector.h"

#include <cstddef>
#include <vector>

namespace svmp::FE::quadrature {

/** @addtogroup FE_Quadrature
 * @{
 */

/**
 * @brief Zero-initialized three-component reference coordinate.
 *
 * Only the first QuadratureRule::dimension() components are active. Remaining
 * components are zero, giving point, line, surface, and volume rules a uniform
 * representation.
 */
class QuadPoint {
public:
    /** @brief Fixed-size FE vector used as the coordinate storage type. */
    using CoordinateVector = math::Vector<double, 3>;

    /** @brief Construct the origin, with all three components initialized to zero. */
    QuadPoint() : coordinates_(CoordinateVector::Zero()) {}

    /**
     * @brief Construct a line coordinate and zero-fill the remaining components.
     * @param x First reference coordinate.
     */
    explicit QuadPoint(double x)
        : coordinates_(x, 0.0, 0.0)
    {
    }

    /**
     * @brief Construct a surface coordinate and zero-fill the third component.
     * @param x First reference coordinate.
     * @param y Second reference coordinate.
     */
    QuadPoint(double x, double y)
        : coordinates_(x, y, 0.0)
    {
    }

    /**
     * @brief Construct a three-component coordinate.
     * @param x First reference coordinate.
     * @param y Second reference coordinate.
     * @param z Third reference coordinate.
     */
    QuadPoint(double x, double y, double z)
        : coordinates_(x, y, z)
    {
    }

    /**
     * @brief Construct from an initialized fixed-size FE vector.
     * @param coordinates Three reference-coordinate components.
     */
    explicit QuadPoint(const CoordinateVector& coordinates)
        : coordinates_(coordinates)
    {
    }

    /**
     * @brief Return the zero coordinate.
     * @return Three-component coordinate initialized to the origin.
     */
    static QuadPoint Zero() { return {}; }

    /**
     * @brief Access one coordinate component without bounds checking.
     * @param component Component index in the half-open range `[0, 3)`.
     * @return Mutable component reference for rule construction.
     */
    double& operator[](std::size_t component) noexcept
    {
        return coordinates_[static_cast<Eigen::Index>(component)];
    }

    /**
     * @brief Access one coordinate component without bounds checking.
     * @param component Component index in the half-open range `[0, 3)`.
     * @return Immutable component reference.
     */
    const double& operator[](std::size_t component) const noexcept
    {
        return coordinates_[static_cast<Eigen::Index>(component)];
    }

    /**
     * @brief Return mutable fixed-size FE coordinate storage for rule generation.
     * @return Mutable three-component coordinate storage.
     */
    CoordinateVector& coordinates() noexcept { return coordinates_; }

    /**
     * @brief Return immutable fixed-size FE coordinate storage.
     * @return Immutable three-component coordinate storage.
     */
    const CoordinateVector& coordinates() const noexcept { return coordinates_; }

private:
    CoordinateVector coordinates_;
};

/**
 * @brief Immutable consumer interface for a quadrature rule on a reference cell.
 *
 * Concrete rule providers initialize the base with one complete RuleData
 * payload. The payload is validated before construction returns, and the object
 * exposes no mutation or assignment path afterward. General solver consumers
 * use only the public const query interface.
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
     * Structural validation does not independently prove this degree.
     * Concrete rule providers establish it through analytic moment tests.
     *
     * @return Declared total degree that a conforming rule integrates exactly.
     */
    int polynomial_exactness() const noexcept { return polynomial_exactness_; }

    /**
     * @brief Return the intrinsic dimension of the reference integration domain.
     * @return Active coordinate count, from zero for Point through three for volume cells.
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
     * @brief Return the default tolerance used during structural validation.
     *
     * Coordinate bounds and inactive components use this as an absolute
     * tolerance. Compensated and ordinary stored-order zeroth-moment checks
     * scale it by the larger of one and the reference-cell measure. For signed
     * rules, a separate cancellation-sensitivity bound based on double epsilon
     * and the absolute weight sum always uses this default tolerance, so callers
     * cannot waive the minimum double-precision stability requirement.
     *
     * @return Default absolute validation tolerance.
     */
    static constexpr double default_validation_tolerance() noexcept { return 1.0e-12; }

    /**
     * @brief Recheck the rule's structural invariants using a caller tolerance.
     *
     * This diagnostic verifies metadata, storage, finite values, point
     * containment, and the zeroth moment. It does not prove
     * polynomial_exactness(). Construction already performs this check with
     * default_validation_tolerance().
     *
     * @param tol Finite, non-negative absolute coordinate and scaled measure tolerance.
     * @return True when the stored rule satisfies every structural invariant.
     */
    bool is_structurally_valid(
        double tol = default_validation_tolerance()) const noexcept;

    /**
     * @brief Return the measure of the canonical reference cell.
     * @return Reference length, area, volume, or unit point measure.
     */
    double reference_measure() const noexcept { return reference_measure_; }

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
     * Dimension and reference measure are derived from @p family; callers cannot
     * supply redundant topology metadata.
     *
     * @param family Supported canonical reference-cell family.
     * @param data Complete exactness, point, and weight payload.
     * @throws InvalidArgumentException If the family is unsupported, exactness
     * is negative, storage is empty or mismatched, a value is non-finite, a point
     * is outside the reference cell, the weights do not reproduce its measure in
     * compensated and stored-order arithmetic, or their cancellation is too
     * ill-conditioned for stable double-precision integration.
     */
    explicit QuadratureRule(svmp::CellFamily family, RuleData data);

private:
    /** @brief Fully checked state used by the delegating constructor. */
    struct ValidatedState {
        svmp::CellFamily cell_family;
        int dimension;
        int polynomial_exactness;
        double reference_measure;
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
    const double reference_measure_;         ///< Canonical cell measure.
    const std::vector<QuadPoint> points_;     ///< Ordered immutable reference coordinates.
    const std::vector<double> weights_;       ///< Immutable weights paired with points_.
};

/** @} */

} // namespace svmp::FE::quadrature

#endif // SVMP_FE_QUADRATURE_RULE_H
