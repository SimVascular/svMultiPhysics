// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_QUADRATURE_RULE_H
#define SVMP_FE_QUADRATURE_RULE_H

/**
 * @file QuadratureRule.h
 * @brief Validated reference-space quadrature rule value type.
 * @ingroup FE_Quadrature
 */

/**
 * @defgroup FE_Quadrature Quadrature
 * @ingroup FE
 * @brief Integration rules on canonical reference cells.
 *
 * @details
 * A QuadratureRule owns ordered reference coordinates and weights for
 * @f[
 *   \int_{\hat K} f(\hat x)\,d\hat x
 *   \approx \sum_q w_q f(\hat x_q).
 * @f]
 * Supported families are Point, Line, Triangle, Quad, Tetra, Hex, and Wedge.
 * The family determines the reference dimension and cell measure. Generating
 * code is responsible for establishing declared polynomial exactness through
 * analytic moment tests.
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
 * Only the first QuadratureRule::dimension() components are active; remaining
 * components must be zero.
 */
using QuadPoint = math::Vector<double, 3>;

/**
 * @brief Owning quadrature rule on a canonical reference cell.
 *
 * Construction requires:
 *
 * - a supported cell family and non-negative polynomial exactness;
 * - at least one point and the same number of points and weights;
 * - finite coordinates and weights, with inactive coordinates equal to zero
 *   within the coordinate tolerance; and
 * - weights whose sum matches the reference-cell measure within the scaled
 *   measure tolerance.
 *
 * Points may be duplicate or outside the reference cell, and weights may be
 * zero or negative. Construction does not verify polynomial exactness.
 */
class QuadratureRule final {
public:
    /**
     * @brief Construct a rule from complete point and weight data.
     * @param family Reference-cell family; also determines dimension and measure.
     * @param polynomial_exactness Declared total-degree polynomial exactness.
     * @param points Ordered reference coordinates.
     * @param weights Weights paired with @p points.
     * @throws InvalidArgumentException If a construction requirement is violated.
     */
    explicit QuadratureRule(
        svmp::CellFamily family,
        int polynomial_exactness,
        std::vector<QuadPoint> points,
        std::vector<double> weights);

    /** @brief Return the number of point/weight pairs. */
    std::size_t num_points() const noexcept { return points_.size(); }

    /** @brief Return the declared total-degree polynomial exactness. */
    int polynomial_exactness() const noexcept { return polynomial_exactness_; }

    /**
     * @brief Return the reference dimension and active QuadPoint component count.
     */
    int dimension() const noexcept;

    /** @brief Return the canonical reference-cell family. */
    svmp::CellFamily cell_family() const noexcept { return cell_family_; }

    /**
     * @brief Return point @p i without bounds checking.
     * @pre @p i is less than num_points().
     */
    const QuadPoint& point(std::size_t i) const noexcept { return points_[i]; }

    /**
     * @brief Return the weight paired with point @p i without bounds checking.
     * @pre @p i is less than num_points().
     */
    double weight(std::size_t i) const noexcept { return weights_[i]; }

    /** @brief Return all points in integration order. */
    const std::vector<QuadPoint>& points() const noexcept { return points_; }

    /** @brief Return all weights in point order. */
    const std::vector<double>& weights() const noexcept { return weights_; }

    /** @brief Return the reference-cell measure derived from cell_family(). */
    double reference_cell_measure() const noexcept;

private:
    svmp::CellFamily cell_family_;          ///< Canonical reference topology.
    int polynomial_exactness_;              ///< Exactness declared by the generator.
    std::vector<QuadPoint> points_;          ///< Ordered reference coordinates.
    std::vector<double> weights_;            ///< Weights paired with points_.
};

/** @} */

} // namespace svmp::FE::quadrature

#endif // SVMP_FE_QUADRATURE_RULE_H
