// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_TRIANGLE_QUADRATURE_H
#define SVMP_FE_TRIANGLE_QUADRATURE_H

/**
 * @file TriangleQuadrature.h
 * @brief Centroid and bounded Duffy quadrature on the unit right triangle.
 * @ingroup FE_Quadrature
 */

#include "FE/Quadrature/QuadratureRule.h"

namespace svmp::FE::quadrature {

/** @addtogroup FE_Quadrature
 * @{
 */

/**
 * @brief Return the project bound on requested triangle exactness.
 * @return 254, permitting a 128-point Gauss factor.
 */
[[nodiscard]] constexpr int max_triangle_exactness() noexcept
{
    return 254;
}

/**
 * @brief Return the degree-one centroid rule: point @f$(1/3,1/3,0)@f$, weight @f$1/2@f$.
 * @return Owning CellFamily::Triangle rule for constant or linear integration.
 */
[[nodiscard]] QuadratureRule make_triangle_centroid_rule();

/**
 * @brief Compose two canonical line rules into a Duffy triangle rule.
 *
 * @details Duffy maps the unit square onto @f$x,y\geq0,\ x+y\leq1@f$ by
 * @f$(x,y)=(u,(1-u)v)@f$, collapsing the edge @f$u=1@f$ to vertex @f$(1,0)@f$.
 * For supplied points @f$\xi_i,\eta_j@f$ and weights @f$a_i,b_j@f$,
 * @f[
 * u_i=(\xi_i+1)/2,\quad v_j=(\eta_j+1)/2,\quad
 * \mathbf{x}_{ij}=(u_i,(1-u_i)v_j,0),\quad W_{ij}=a_i b_j(1-u_i)/4.
 * @f]
 * Weights include the Jacobian @f$1-u@f$ and both interval scalings.
 * Including the Jacobian, @f$x^a y^b@f$ becomes @f$u^a(1-u)^{b+1}v^b@f$;
 * total degree @f$p@f$ therefore requires @f$e_u\geq p+1@f$ and @f$e_v\geq p@f$.
 *
 * Factors may come from different line generators and have different point
 * counts; both are used as supplied. Select make_triangle_centroid_rule()
 * for the centroid. All samples are kept in input order (i outer, j inner),
 * including duplicates and collapsed points. Signed/zero weights are allowed.
 * Duffy rules are not permutation-symmetric triangle tables.
 *
 * Each rule requires @f$|\sum w-M|\leq\tau\max(M,\sum|w|)@f$, with @f$M=2@f$
 * for each line, @f$M=1/2@f$ for the triangle, and @f$\tau=256\epsilon@f$ for
 * binary64 machine epsilon @f$\epsilon@f$. Mapped points lie in the triangle up
 * to rounding; no clipping or weight rescaling is performed.
 *
 * @param requested_exactness Minimum total degree @f$p@f$; must be in [0, 254],
 *        inclusive (see max_triangle_exactness()).
 * @param u_rule Nonempty canonical Line rule on @f$[-1,1]@f$ with trusted
 *        @f$e_u\geq p+1@f$.
 * @param v_rule Nonempty canonical Line rule on @f$[-1,1]@f$ with trusted
 *        @f$e_v\geq p@f$.
 * @return Owning Triangle rule with exactness @f$\min(e_u-1,e_v)@f$ and
 *         @f$n_u n_v@f$ samples.
 * @throws InvalidArgumentException If the request or factors violate these
 *         contracts, or @f$n_u n_v@f$ exceeds 16384; checked before product
 *         allocation.
 * @throws ConvergenceException If the mapped weights do not reproduce measure
 *         @f$1/2@f$, for example after an overflowed product or when @p u_rule
 *         does not integrate linear functions exactly.
 */
[[nodiscard]] QuadratureRule make_triangle_rule(
    int requested_exactness, const QuadratureRule& u_rule, const QuadratureRule& v_rule);

/**
 * @brief Reuse one line rule in both Duffy directions.
 * @details Equivalent to make_triangle_rule(requested_exactness, line_rule, line_rule),
 * with the same validation and sample-preservation contracts.
 * @param requested_exactness Minimum total degree @f$p@f$; must be in [0, 254],
 *        inclusive (see max_triangle_exactness()).
 * @param line_rule Canonical Line rule with @f$1\leq n\leq128@f$ points and
 *        trusted exactness @f$e\geq p+1@f$.
 * @return Owning Triangle rule with exactness @f$e-1@f$ and @f$n^2@f$ samples.
 */
[[nodiscard]] QuadratureRule make_triangle_rule(int requested_exactness, const QuadratureRule& line_rule);

/** @} */

} // namespace svmp::FE::quadrature

#endif // SVMP_FE_TRIANGLE_QUADRATURE_H
