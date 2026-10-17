// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_GAUSS_QUADRATURE_H
#define SVMP_FE_GAUSS_QUADRATURE_H

/**
 * @file GaussQuadrature.h
 * @brief Bounded Gauss-Legendre quadrature generation on the reference line.
 * @ingroup FE_Quadrature
 */

#include "FE/Quadrature/QuadratureRule.h"

namespace svmp::FE::quadrature {

/** @addtogroup FE_Quadrature
 * @{
 */

/**
 * @brief Return the largest supported Gauss-Legendre exactness request.
 * @details Degree 255 corresponds to the internal 128-point project support
 * bound, which limits generator work and downstream product-rule growth; it
 * is not a mathematical or convergence limit.
 * @return The inclusive requested-exactness limit, 255.
 */
[[nodiscard]] constexpr int max_gauss_legendre_exactness() noexcept
{
    return 255;
}

/**
 * @brief Generate a Gauss-Legendre rule on @f$[-1,1]@f$ with at least the
 *        requested exactness.
 *
 * @details For @f$d@f$ = @p requested_exactness, integer division gives the
 * minimum point count @f$n=\lfloor d/2\rfloor+1@f$. The returned metadata reports
 * the actual polynomial exactness @f$2n-1@f$, which exceeds even requests by one.
 * Degree zero produces one point with exactness one. Points are the roots of
 * the Legendre polynomial @f$P_n@f$ of degree @f$n@f$, in strictly increasing
 * order inside @f$(-1,1)@f$; neither endpoint is included. Weights are positive
 * and aligned with their points.
 *
 * @see [NIST DLMF: Legendre polynomials](https://dlmf.nist.gov/18.3)
 *
 * @param requested_exactness Minimum polynomial degree to integrate exactly;
 *        must be in [0, 255], inclusive (see max_gauss_legendre_exactness()).
 * @return A complete QuadratureRule value for CellFamily::Line.
 * @throws InvalidArgumentException If @p requested_exactness is outside the
 *         supported range; checked before point-count conversion or allocation.
 * @throws ConvergenceException If root refinement or final numerical
 *         validation fails.
 */
[[nodiscard]] QuadratureRule
make_gauss_legendre_rule(int requested_exactness);

/** @} */

} // namespace svmp::FE::quadrature

#endif // SVMP_FE_GAUSS_QUADRATURE_H
