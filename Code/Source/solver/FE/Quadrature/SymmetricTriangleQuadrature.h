// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_SYMMETRIC_TRIANGLE_QUADRATURE_H
#define SVMP_FE_SYMMETRIC_TRIANGLE_QUADRATURE_H

/**
 * @file SymmetricTriangleQuadrature.h
 * @brief Contained Dunavant rules on the unit right triangle.
 * @ingroup FE_Quadrature
 */

#include "FE/Quadrature/QuadratureRule.h"

namespace svmp::FE::quadrature {

/** @addtogroup FE_Quadrature
 * @{
 */

/**
 * @brief Return the largest supported symmetric triangle exactness.
 * @return 19.
 */
[[nodiscard]] constexpr int max_symmetric_triangle_exactness() noexcept
{
    return 19;
}

/**
 * @brief Return the smallest Dunavant rule meeting the requested degree.
 *
 * @details **Reference:** D. A. Dunavant (1985), High degree efficient symmetrical
 * Gaussian quadrature rules for the triangle. https://doi.org/10.1002/nme.1620210612
 *
 * The reference tabulates rules of degrees 1 to 20. Rules 11, 15, 16, 18, and 20
 * place some points outside the triangle, where the integrand may not be defined,
 * so only degrees 1–10, 12, 13, 14, 17, and 19 are implemented. A request for any
 * other degree returns the next higher implemented one (request 0 returns the
 * degree-1 centroid rule); polynomial_exactness() reports the delivered degree.
 *
 * The rules are symmetric: points come in groups of one (the centroid), three, or
 * six that share a weight, each group being the permutations of one point's
 * barycentric coordinates @f$(1-x-y,x,y)@f$. Groups follow the reference's order.
 *
 * Points lie in the unit triangle @f$x,y\geq0,\ x+y\leq1@f$, with third coordinate
 * exactly zero. Weights sum to the area @f$1/2@f$ within @f$5\cdot10^{-15}@f$, the
 * rounding of the reference's 15-digit tables. The degree-3 and degree-7 rules have
 * a negative weight at the centroid, as in the reference. No clipping or
 * renormalization is applied.
 *
 * @param[in] requested_exactness Minimum total degree in [0, 19], inclusive.
 * @return Owning Triangle rule reporting its delivered polynomial exactness.
 * @throws InvalidArgumentException If the request is outside [0, 19], checked
 *         before table selection or allocation.
 */
[[nodiscard]] QuadratureRule make_symmetric_triangle_rule(int requested_exactness);

/** @} */

} // namespace svmp::FE::quadrature

#endif // SVMP_FE_SYMMETRIC_TRIANGLE_QUADRATURE_H
