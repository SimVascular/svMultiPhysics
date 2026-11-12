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
 * @brief Compose three canonical line rules into a Duffy tetrahedron rule.
 *
 * @details Duffy maps the unit cube onto @f$x,y,z\geq0,\ x+y+z\leq1@f$.
 * For supplied points @f$\xi_i,\eta_j,\zeta_k@f$ and weights @f$a_i,b_j,c_k@f$,
 * @f[
 * u_i=(\xi_i+1)/2,\quad v_j=(\eta_j+1)/2,\quad w_k=(\zeta_k+1)/2,
 * @f]
 * @f[
 * \mathbf{x}_{ijk}=(u_i,(1-u_i)v_j,(1-u_i)(1-v_j)w_k),\qquad
 * W_{ijk}=a_i b_j c_k(1-u_i)^2(1-v_j)/8.
 * @f]
 * Weights include the Jacobian and all three interval scalings. Including the
 * Jacobian, @f$x^a y^b z^c@f$ becomes
 * @f$u^a(1-u)^{b+c+2}v^b(1-v)^{c+1}w^c@f$; total degree @f$p@f$ requires
 * @f$e_u\geq p+2@f$, @f$e_v\geq p+1@f$, and @f$e_w\geq p@f$.
 *
 * Factors may differ in generator, count, order, and exactness. All samples
 * are used in input order, @f$q=(i n_v+j)n_w+k@f$, including duplicates,
 * collapsed points, and signed/zero weights. Requests zero and one still use
 * the supplied factors; select make_tetrahedron_centroid_rule() explicitly.
 * Duffy rules are not permutation-symmetric tetrahedron tables.
 *
 * Each rule requires @f$|\sum w-M|\leq\tau\max(M,\sum|w|)@f$, with @f$M=2@f$
 * for each line, @f$M=1/6@f$ for the tetrahedron, and @f$\tau=256\epsilon@f$
 * for binary64 machine epsilon @f$\epsilon@f$. Factor exactness metadata is
 * trusted: this measure check does not certify individual factor moments,
 * since directional errors can compensate. Mapped points lie in the
 * tetrahedron up to rounding; no clipping or weight rescaling is performed.
 *
 * @param requested_exactness Minimum total degree @f$p@f$ in [0, 253],
 *        inclusive (see max_tetrahedron_exactness()).
 * @param u_rule Nonempty canonical Line rule on @f$[-1,1]@f$ with trusted
 *        @f$e_u\geq p+2@f$.
 * @param v_rule Nonempty canonical Line rule on @f$[-1,1]@f$ with trusted
 *        @f$e_v\geq p+1@f$.
 * @param w_rule Nonempty canonical Line rule on @f$[-1,1]@f$ with trusted
 *        @f$e_w\geq p@f$.
 * @return Owning Tetra rule with exactness @f$\min(e_u-2,e_v-1,e_w)@f$ and
 *         @f$n_u n_v n_w@f$ samples; the request limit does not clamp metadata.
 * @throws InvalidArgumentException If the request or factors violate these
 *         contracts, or @f$n_u n_v n_w@f$ exceeds 2097152. The request,
 *         nonempty inputs, and capacity are checked before degree conversion,
 *         product multiplication, and output allocation.
 * @throws ConvergenceException If generated weights are nonfinite, including
 *         after intermediate product overflow, or do not reproduce measure
 *         @f$1/6@f$ within the stated tolerance.
 */
[[nodiscard]] QuadratureRule make_tetrahedron_rule(
    int requested_exactness, const QuadratureRule& u_rule,
    const QuadratureRule& v_rule, const QuadratureRule& w_rule);

/**
 * @brief Reuse one line rule in all three Duffy directions.
 * @details Equivalent to make_tetrahedron_rule(requested_exactness, line_rule,
 * line_rule, line_rule), with the same validation and sample-preservation
 * contracts. Selecting make_gauss_legendre_rule(p + 2) for a validated request
 * gives @f$n=\lfloor p/2\rfloor+2@f$ points per direction and degree @f$2n-3@f$.
 * @param requested_exactness Minimum total degree @f$p@f$ in [0, 253], inclusive.
 * @param line_rule Canonical Line rule with @f$1\leq n\leq128@f$ points and
 *        trusted exactness @f$e\geq p+2@f$.
 * @return Owning Tetra rule with exactness @f$e-2@f$ and @f$n^3@f$ samples.
 */
[[nodiscard]] QuadratureRule make_tetrahedron_rule(
    int requested_exactness, const QuadratureRule& line_rule);

/** @} */

} // namespace svmp::FE::quadrature

#endif // SVMP_FE_TETRAHEDRON_QUADRATURE_H
