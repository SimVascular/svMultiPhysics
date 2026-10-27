// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_TRIANGLE_QUADRATURE_H
#define SVMP_FE_TRIANGLE_QUADRATURE_H

/** @file TriangleQuadrature.h
 * @brief Centroid and bounded Duffy quadrature on the unit right triangle.
 * @ingroup FE_Quadrature
 */

#include "FE/Quadrature/QuadratureRule.h"

namespace svmp::FE::quadrature {

/** @addtogroup FE_Quadrature
 * @{
 */

/** @brief Return the project bound on requested triangle exactness.
 * @return 254, permitting a 128-point Gauss factor.
 */
[[nodiscard]] constexpr int max_triangle_exactness() noexcept { return 254; }

/** @brief Return the degree-one centroid rule: point (1/3,1/3,0), weight 1/2.
 * @return Owning CellFamily::Triangle rule for constant or linear integration.
 */
[[nodiscard]] QuadratureRule make_triangle_centroid_rule();

/** @brief Compose two canonical line rules into a Duffy triangle rule.
 * @details Duffy maps the unit square onto x,y >= 0, x+y <= 1 by
 * (x,y)=(u,(1-u)v), collapsing the edge u=1 to vertex (1,0).
 * For supplied points xi_i, eta_j and weights a_i, b_j,
 * @f[
 * u_i=(\xi_i+1)/2,\quad v_j=(\eta_j+1)/2,\quad
 * \mathbf{x}_{ij}=(u_i,(1-u_i)v_j,0),\quad W_{ij}=a_i b_j(1-u_i)/4.
 * @f]
 * Weights include the Jacobian (1-u) and both interval scalings.
 * Including the Jacobian, x^a y^b becomes u^a(1-u)^(b+1)v^b;
 * total degree p therefore requires e_u >= p+1 and e_v >= p.
 *
 * Factors may come from different line generators and have different point
 * counts; both are used as supplied. Select make_triangle_centroid_rule()
 * for the centroid. All samples are kept in input order (i outer, j inner),
 * including duplicates and collapsed points. Signed/zero weights are allowed.
 * Duffy rules are not permutation-symmetric triangle tables.
 *
 * Each rule requires abs(sum(w)-M) <= tau*max(M,sum(abs(w))), with M=2
 * for each line, M=1/2 for the triangle, and tau=256*epsilon(double).
 * Mapped points lie in the triangle up to rounding; no clipping or weight
 * rescaling is performed.
 *
 * @param requested_exactness Minimum total degree p in [0,254].
 * @param u_rule Nonempty canonical Line rule on [-1,1] with trusted e_u >= p+1.
 * @param v_rule Nonempty canonical Line rule on [-1,1] with trusted e_v >= p.
 * @return Owning Triangle rule with exactness min(e_u-1,e_v) and n_u*n_v samples.
 * @throws InvalidArgumentException If the request or factors violate these
 * contracts, or n_u*n_v exceeds 16384; checked before product allocation.
 * @throws ConvergenceException If mapped data fail numerical validation.
 */
[[nodiscard]] QuadratureRule make_triangle_rule(
    int requested_exactness, const QuadratureRule& u_rule, const QuadratureRule& v_rule);

/** @brief Reuse one line rule in both Duffy directions.
 * @details Equivalent to make_triangle_rule(requested_exactness, line_rule, line_rule),
 * with the same validation and sample-preservation contracts.
 * @param requested_exactness Minimum total degree p in [0,254].
 * @param line_rule Canonical Line rule with n=1..128 and trusted exactness e >= p+1.
 * @return Owning Triangle rule with exactness e-1 and n*n samples.
 */
[[nodiscard]] QuadratureRule make_triangle_rule(int requested_exactness, const QuadratureRule& line_rule);

/** @} */
} // namespace svmp::FE::quadrature

#endif // SVMP_FE_TRIANGLE_QUADRATURE_H
