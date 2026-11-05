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
 * @brief Return the smallest retained Dunavant rule meeting the requested degree.
 * @details The retained degrees are 1–10, 12, 13, 14, 17, and 19. Requests zero
 * and one return the centroid. Gaps advance to the next retained degree;
 * Dunavant tables 11, 15, 16, 18, and 20 are excluded for exterior points.
 *
 * Points lie in @f$x,y\geq0,\ x+y\leq1@f$, with third coordinate exactly zero.
 * Weights sum to area @f$1/2@f$ within tabulation precision. Degrees 3 and 7
 * retain their negative centroid weights; no clipping or renormalization is applied.
 *
 * Coefficients and orbit order follow Dunavant (1985), doi:10.1002/nme.1620210612.
 * Each orbit has equal weights. A centroid emits @f$(1/3,1/3)@f$; a three-point
 * orbit with repeated barycentric coordinate @f$a@f$ emits @f$(a,c),(a,a),(c,a)@f$,
 * @f$c=1-2a@f$; a six-point orbit with first two barycentric coordinates
 * @f$a,b@f$ emits @f$(b,c),(c,b),(a,c),(c,a),(a,b),(b,a)@f$, @f$c=1-a-b@f$.
 * @param requested_exactness Minimum total degree in [0, 19], inclusive.
 * @return Owning Triangle rule reporting its delivered polynomial exactness.
 * @throws InvalidArgumentException If the request is outside [0, 19], checked
 *         before table selection or allocation.
 */
[[nodiscard]] QuadratureRule make_symmetric_triangle_rule(int requested_exactness);

/** @} */

} // namespace svmp::FE::quadrature

#endif // SVMP_FE_SYMMETRIC_TRIANGLE_QUADRATURE_H
