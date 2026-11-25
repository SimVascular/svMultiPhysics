// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file TetrahedronQuadrature.cpp
 * @brief Validated composition of supplied line quadrature into tetrahedron rules.
 * @ingroup FE_Quadrature
 */

#include "FE/Quadrature/TetrahedronQuadrature.h"

#include "FE/Common/FEException.h"
#include "FE/Quadrature/DuffyQuadrature.h"

#include <algorithm>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace svmp::FE::quadrature {
namespace {

constexpr std::size_t kMaximumProductPoints = 128u * 128u * 128u;

} // namespace

QuadratureRule make_tetrahedron_centroid_rule()
{
    return QuadratureRule(svmp::CellFamily::Tetra, 1, {{0.25, 0.25, 0.25}}, {1.0 / 6.0});
}

QuadratureRule make_tetrahedron_rule(
    int requested_exactness, const QuadratureRule& u_rule,
    const QuadratureRule& v_rule, const QuadratureRule& w_rule)
{
    svmp::check<InvalidArgumentException>(
        requested_exactness >= 0 && requested_exactness <= max_tetrahedron_exactness(),
        "Tetrahedron Duffy generator: requested_exactness must be in [0, " +
            std::to_string(max_tetrahedron_exactness()) + ']');
    const std::size_t n_u = u_rule.num_points();
    const std::size_t n_v = v_rule.num_points();
    const std::size_t n_w = w_rule.num_points();
    svmp::check<InvalidArgumentException>(n_u != 0 && n_v != 0 && n_w != 0,
        "Tetrahedron Duffy generator: empty input rule");
    // Division bounds the full product and each validation sum before multiplication.
    svmp::check<InvalidArgumentException>(
        n_u <= kMaximumProductPoints / n_v / n_w,
        "Tetrahedron Duffy generator: product exceeds " +
            std::to_string(kMaximumProductPoints) + " samples");
    detail::validate_line_rule(u_rule, requested_exactness + 2, kMaximumProductPoints,
        "Tetrahedron Duffy generator: u_rule");
    detail::validate_line_rule(v_rule, requested_exactness + 1, kMaximumProductPoints,
        "Tetrahedron Duffy generator: v_rule");
    detail::validate_line_rule(w_rule, requested_exactness, kMaximumProductPoints,
        "Tetrahedron Duffy generator: w_rule");

    // Checked interval nodes and the Duffy map establish containment up to rounding.
    std::vector<QuadPoint> points(n_u * n_v * n_w, QuadPoint::Zero());
    std::vector<double> weights(points.size());
    for (std::size_t i = 0; i < n_u; ++i) {
        const double xi = u_rule.point(i)[0];
        const double u = (1.0 + xi) / 2.0;
        const double remaining_u = detail::unit_interval_complement(xi);
        const double u_weight = u_rule.weight(i);
        const double u_scale = remaining_u * remaining_u / 8.0;
        const double outer_weight = u_weight * u_scale;
        const bool scale_outer = detail::needs_scaled_product(outer_weight, u_weight, u_scale);
        for (std::size_t j = 0; j < n_v; ++j) {
            const double eta = v_rule.point(j)[0];
            const double v = (1.0 + eta) / 2.0;
            const double remaining_v = detail::unit_interval_complement(eta);
            const double v_weight = v_rule.weight(j);
            const double weighted_v = v_weight * remaining_v;
            const double middle_weight = outer_weight * weighted_v;
            const bool rescale = scale_outer || detail::needs_scaled_product(weighted_v, v_weight, remaining_v) ||
                detail::needs_scaled_product(middle_weight, outer_weight, weighted_v);
            for (std::size_t k = 0; k < n_w; ++k) {
                const std::size_t q = (i * n_v + j) * n_w + k;
                const double w = (w_rule.point(k)[0] + 1.0) / 2.0;
                points[q] = QuadPoint{u, remaining_u * v, remaining_u * remaining_v * w};
                const double w_weight = w_rule.weight(k);
                weights[q] = rescale
                    ? detail::scaled_product({u_weight, v_weight, w_weight, remaining_u, remaining_u, remaining_v}, -3)
                    : middle_weight * w_weight;
            }
        }
    }
    // Validate transformed volume and reject nonfinite products. This does not
    // certify individual factor moments; directional errors can compensate.
    svmp::check<ConvergenceException>(detail::reproduces_measure(weights, 1.0 / 6.0, kMaximumProductPoints),
        "Tetrahedron Duffy generator: output weights do not reproduce measure 1/6");

    const int polynomial_exactness = std::min({
        u_rule.polynomial_exactness() - 2, v_rule.polynomial_exactness() - 1,
        w_rule.polynomial_exactness()});
    return QuadratureRule(
        svmp::CellFamily::Tetra, polynomial_exactness, std::move(points), std::move(weights));
}

QuadratureRule make_tetrahedron_rule(int requested_exactness, const QuadratureRule& line_rule)
{
    return make_tetrahedron_rule(requested_exactness, line_rule, line_rule, line_rule);
}

} // namespace svmp::FE::quadrature
