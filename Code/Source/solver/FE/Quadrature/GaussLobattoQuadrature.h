// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_GAUSS_LOBATTO_QUADRATURE_H
#define SVMP_FE_GAUSS_LOBATTO_QUADRATURE_H

/**
 * @file GaussLobattoQuadrature.h
 * @brief Bounded Gauss-Lobatto-Legendre quadrature generation on the reference line.
 * @ingroup FE_Quadrature
 */

#include "FE/Quadrature/QuadratureRule.h"

namespace svmp::FE::quadrature {

/** @addtogroup FE_Quadrature
 * @{
 */

/**
 * @brief Return the largest supported Gauss-Lobatto-Legendre exactness request.
 * @details Degree 253 corresponds to the internal 128-point project support
 * bound, including both endpoints, which limits generator work and downstream
 * product-rule growth; it is not a mathematical or convergence limit.
 * @return The inclusive requested-exactness limit, 253.
 */
[[nodiscard]] constexpr int max_gauss_lobatto_exactness() noexcept
{
    return 253;
}

/**
 * @brief Generate a Gauss-Lobatto-Legendre rule on @f$[-1,1]@f$ with at least
 *        the requested exactness.
 *
 * @details For @f$d@f$ = @p requested_exactness, integer division gives the
 * minimum point count @f$n=\lfloor d/2\rfloor+2@f$. The returned metadata reports
 * the actual polynomial exactness @f$2n-3@f$, which exceeds even requests by one.
 * Degree zero produces two points with exactness one. The first and last points
 * are exactly @f$-1@f$ and @f$+1@f$. When present, the @f$n-2@f$ interior points
 * are the roots of @f$P'_{n-1}@f$, the derivative of the Legendre polynomial
 * of degree @f$n-1@f$. Points are strictly increasing, and weights are positive
 * and aligned with their points.
 *
 * @see [NIST DLMF: Legendre polynomials](https://dlmf.nist.gov/18.3)
 *
 * @param requested_exactness Minimum polynomial degree to integrate exactly;
 *        must be in [0, 253], inclusive (see max_gauss_lobatto_exactness()).
 * @return A complete QuadratureRule value for CellFamily::Line.
 * @throws InvalidArgumentException If @p requested_exactness is outside the
 *         supported range; checked before point-count conversion or allocation.
 * @throws ConvergenceException If root refinement or final numerical
 *         validation fails.
 */
[[nodiscard]] QuadratureRule
make_gauss_lobatto_rule(int requested_exactness);

/** @} */

} // namespace svmp::FE::quadrature

#endif // SVMP_FE_GAUSS_LOBATTO_QUADRATURE_H
