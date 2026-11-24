// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file TriangleQuadrature.cpp
 * @brief Validated composition of supplied line quadrature into triangle rules.
 * @ingroup FE_Quadrature
 */

#include "FE/Quadrature/TriangleQuadrature.h"

#include "FE/Common/FEException.h"
#include "FE/Quadrature/DuffyQuadrature.h"

#include <algorithm>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace svmp::FE::quadrature {
namespace {

constexpr std::size_t kMaximumProductPoints = 128u * 128u;

} // namespace

QuadratureRule make_triangle_centroid_rule()
{
    return QuadratureRule(svmp::CellFamily::Triangle, 1, {{1.0 / 3.0, 1.0 / 3.0, 0.0}}, {0.5});
}

QuadratureRule make_triangle_rule(
    int requested_exactness, const QuadratureRule& u_rule, const QuadratureRule& v_rule)
{
    svmp::check<InvalidArgumentException>(
        requested_exactness >= 0 && requested_exactness <= max_triangle_exactness(),
        "Triangle Duffy generator: requested_exactness must be in [0, " +
            std::to_string(max_triangle_exactness()) + ']');
    const std::size_t n_u = u_rule.num_points();
    const std::size_t n_v = v_rule.num_points();
    svmp::check<InvalidArgumentException>(
        n_u != 0 && n_v != 0, "Triangle Duffy generator: empty input rule");
    // Bound the product and both validation sums before multiplication or allocation.
    svmp::check<InvalidArgumentException>(
        n_u <= kMaximumProductPoints / n_v,
        "Triangle Duffy generator: product exceeds " + std::to_string(kMaximumProductPoints) + " samples");
    detail::validate_line_rule(u_rule, requested_exactness + 1, kMaximumProductPoints,
        "Triangle Duffy generator: u_rule");
    detail::validate_line_rule(v_rule, requested_exactness, kMaximumProductPoints,
        "Triangle Duffy generator: v_rule");

    // Checked interval nodes and the Duffy map establish containment up to rounding.
    std::vector<QuadPoint> points(n_u * n_v, QuadPoint::Zero());
    std::vector<double> weights(points.size());
    for (std::size_t i = 0; i < n_u; ++i) {
        const double xi = u_rule.point(i)[0];
        const double u = (1.0 + xi) / 2.0;
        const double remaining = detail::unit_interval_complement(xi);
        // Two interval factors of 1/2 and the Duffy Jacobian (1-u).
        const double u_weight = u_rule.weight(i);
        const double u_scale = remaining / 4.0;
        const double outer_weight = u_weight * u_scale;
        const bool scale_outer = detail::needs_scaled_product(outer_weight, u_weight, u_scale);
        for (std::size_t j = 0; j < n_v; ++j) {
            const std::size_t q = i * n_v + j;
            const double v = (v_rule.point(j)[0] + 1.0) / 2.0;
            points[q] = QuadPoint{u, remaining * v, 0.0};
            const double v_weight = v_rule.weight(j);
            weights[q] = scale_outer
                ? detail::scaled_product({u_weight, v_weight, remaining}, -2)
                : outer_weight * v_weight;
        }
    }
    // Given both factor measures, this checks the u_rule first moment and
    // rejects overflowed products.
    svmp::check<ConvergenceException>(detail::reproduces_measure(weights, 0.5, kMaximumProductPoints),
        "Triangle Duffy generator: output weights do not reproduce measure 1/2");

    // x^a y^b maps to u^a (1-u)^(b+1) v^b, so u_rule needs one degree more than it delivers.
    const int polynomial_exactness = std::min(u_rule.polynomial_exactness() - 1, v_rule.polynomial_exactness());
    return QuadratureRule(
        svmp::CellFamily::Triangle, polynomial_exactness, std::move(points), std::move(weights));
}

QuadratureRule make_triangle_rule(int requested_exactness, const QuadratureRule& line_rule)
{
    return make_triangle_rule(requested_exactness, line_rule, line_rule);
}

} // namespace svmp::FE::quadrature
