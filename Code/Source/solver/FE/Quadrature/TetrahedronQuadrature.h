// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_TETRAHEDRON_QUADRATURE_H
#define SVMP_FE_TETRAHEDRON_QUADRATURE_H

/**
 * @file TetrahedronQuadrature.h
 * @brief Centroid and bounded Duffy quadrature on the unit right tetrahedron.
 * @ingroup FE_Quadrature
 */

#include "FE/Quadrature/QuadratureRule.h"

namespace svmp::FE::quadrature {

/** @addtogroup FE_Quadrature
 * @{
 */

/**
 * @brief Return the project bound on requested tetrahedron exactness.
 * @return 253, permitting a 128-point Gauss factor.
 */
[[nodiscard]] constexpr int max_tetrahedron_exactness() noexcept
{
    return 253;
}

/**
 * @brief Return the degree-one centroid: point @f$(1/4,1/4,1/4)@f$, weight @f$1/6@f$.
 * @return Owning CellFamily::Tetra rule for constant or linear integration.
 */
[[nodiscard]] QuadratureRule make_tetrahedron_centroid_rule();

/**
 * @brief Compose canonical line rules into a Duffy tetrahedron rule.
 *
 * @details For input nodes @f$\xi_i,\eta_j,\zeta_k@f$ and weights @f$a_i,b_j,c_k@f$,
 * @f[
 * u_i=(\xi_i+1)/2,\quad v_j=(\eta_j+1)/2,\quad w_k=(\zeta_k+1)/2,
 * @f]
 * @f[
 * \mathbf{x}_{ijk}=(u_i,(1-u_i)v_j,(1-u_i)(1-v_j)w_k),\qquad
 * W_{ijk}=a_i b_j c_k(1-u_i)^2(1-v_j)/8.
 * @f]
 * The map covers @f$x,y,z\geq0,\ x+y+z\leq1@f$ with Jacobian
 * @f$(1-u)^2(1-v)@f$. It collapses @f$u=1@f$ to vertex @f$(1,0,0)@f$
 * and @f$v=1@f$ to edge @f$(u,1-u,0)@f$.
 *
 * Factors must be nonempty canonical Line rules on @f$[-1,1]@f$; generators,
 * counts, order, and exactness may differ. All samples are kept in input
 * order @f$q=(i n_v+j)n_w+k@f$, including duplicates, collapsed points, and
 * signed/zero weights. Even requests 0 and 1 use the supplied factors;
 * select make_tetrahedron_centroid_rule() explicitly. Permutation symmetry
 * is not guaranteed. Distances to collapsed faces are computed from the
 * original line nodes. The returned points may lie on a face to rounding
 * even when the Jacobian is positive. Weights include the Jacobian without
 * renormalization.
 *
 * Each rule must satisfy @f$|\sum w-M|\leq256\epsilon\max(M,\sum|w|)@f$,
 * where @f$\epsilon@f$ is binary64 machine epsilon, @f$M=2@f$ for inputs,
 * and @f$M=1/6@f$ for output. Exactness metadata is trusted; this check
 * cannot certify factor moments because directional errors can cancel.
 *
 * @param requested_exactness Minimum total degree @f$p\in[0,253]@f$
 *        (see max_tetrahedron_exactness()).
 * @param u_rule Line factor with @f$e_u\geq p+2@f$.
 * @param v_rule Line factor with @f$e_v\geq p+1@f$.
 * @param w_rule Line factor with @f$e_w\geq p@f$.
 * @return Owning Tetra rule with @f$n_u n_v n_w@f$ samples and exactness
 *         @f$\min(e_u-2,e_v-1,e_w)@f$, which may exceed 253.
 * @throws InvalidArgumentException For invalid requests or factors, or
 *         @f$n_u n_v n_w>2097152@f$. Request, nonemptiness, and capacity
 *         checks precede forming directional requirements @f$p+2@f$ and
 *         @f$p+1@f$, multiplying counts, and allocating output.
 * @throws ConvergenceException For nonfinite generated weights (including
 *         intermediate product overflow) or failed output normalization.
 */
[[nodiscard]] QuadratureRule make_tetrahedron_rule(
    int requested_exactness, const QuadratureRule& u_rule,
    const QuadratureRule& v_rule, const QuadratureRule& w_rule);

/**
 * @brief Reuse one line rule in all three Duffy directions.
 * @details Equivalent to make_tetrahedron_rule(requested_exactness, line_rule,
 * line_rule, line_rule); all its contracts apply. For validated @f$p@f$,
 * make_gauss_legendre_rule(p + 2) gives the smallest repeated Gauss factor:
 * @f$n=\lfloor p/2\rfloor+2@f$ points and tetrahedron degree @f$2n-3@f$.
 * @param requested_exactness Minimum total degree @f$p\in[0,253]@f$.
 * @param line_rule Line factor with @f$1\leq n\leq128@f$ and @f$e\geq p+2@f$.
 * @return Owning Tetra rule with @f$n^3@f$ samples and exactness @f$e-2@f$.
 */
[[nodiscard]] QuadratureRule make_tetrahedron_rule(
    int requested_exactness, const QuadratureRule& line_rule);

/** @} */

} // namespace svmp::FE::quadrature

#endif // SVMP_FE_TETRAHEDRON_QUADRATURE_H
