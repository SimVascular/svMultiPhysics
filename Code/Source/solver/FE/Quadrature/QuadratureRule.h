// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_QUADRATURE_RULE_H
#define SVMP_FE_QUADRATURE_RULE_H

/**
 * @file QuadratureRule.h
 * @brief Value type for reference-space quadrature rules.
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
 * Supported families use these canonical reference cells:
 *
 * | Family   | Reference cell                            | Measure   |
 * |----------|-------------------------------------------|-----------|
 * | Point    | @f$(0,0,0)@f$                             | @f$1@f$   |
 * | Line     | @f$[-1,1]@f$                              | @f$2@f$   |
 * | Triangle | @f$\{(x,y):x,y\geq0,\ x+y\leq1\}@f$       | @f$1/2@f$ |
 * | Quad     | @f$[-1,1]^2@f$                            | @f$4@f$   |
 * | Tetra    | @f$\{(x,y,z):x,y,z\geq0,\ x+y+z\leq1\}@f$ | @f$1/6@f$ |
 * | Hex      | @f$[-1,1]^3@f$                            | @f$8@f$   |
 * | Wedge    | unit triangle @f$\times[-1,1]@f$          | @f$1@f$   |
 *
 * The family therefore determines both reference dimension and cell measure.
 * Quadrature points are not required to lie inside the reference cell.
 * Generating code is responsible for verifying weight normalization and
 * declared polynomial exactness through analytic moment tests.
 */

#include "FE/Common/Types.h"
#include "FE/Math/Vector.h"

#include <cstddef>
#include <span>
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
 * @brief Owning value type for quadrature rules on canonical reference cells.
 *
 * Construction requires:
 *
 * - a supported cell family and non-negative polynomial exactness;
 * - a vector of points, with at least one element;
 * - for all points, coordinates must be finite;
 * - for all points, coordinates beyond the reference dimension of the cell
 *   family must be equal to zero; and
 * - a vector of finite weights, with as many elements as the points.
 *
 * The constructor checks each requirement and throws InvalidArgumentException
 * when one is violated.
 *
 * Points may be duplicate or outside the reference cell, and weights may be
 * zero or negative. Construction does not verify weight normalization or
 * polynomial exactness.
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
    std::size_t dimension() const;

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

    /** @brief Return a read-only view of all points in integration order. */
    std::span<const QuadPoint> points() const noexcept { return points_; }

    /** @brief Return a read-only view of all weights in point order. */
    std::span<const double> weights() const noexcept { return weights_; }

    /** @brief Return the reference-cell measure derived from cell_family(). */
    double reference_cell_measure() const;

private:
    svmp::CellFamily cell_family_;          ///< Canonical reference topology.
    int polynomial_exactness_;              ///< Exactness declared by the generator.
    std::vector<QuadPoint> points_;          ///< Ordered reference coordinates.
    std::vector<double> weights_;            ///< Weights paired with points_.
};

/** @} */

} // namespace svmp::FE::quadrature

#endif // SVMP_FE_QUADRATURE_RULE_H
