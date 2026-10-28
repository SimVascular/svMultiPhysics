// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file TriangleQuadrature.cpp
 * @brief Validated composition of supplied line quadrature into triangle rules.
 * @ingroup FE_Quadrature
 */

#include "FE/Quadrature/TriangleQuadrature.h"

#include "FE/Common/FEException.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace svmp::FE::quadrature {
namespace {

constexpr std::size_t kMaximumProductPoints = 128u * 128u;
constexpr double kValidationTolerance = 256.0 * std::numeric_limits<double>::epsilon();
// Scale validation sums by a power of two so even capacity*max_double stays
// within binary64 range. Stored weights and relative tolerances are unchanged.
constexpr double kValidationSumScale = 0.5 / static_cast<double>(kMaximumProductPoints);

// Return whether abs(sum(w)-measure) <= tau*max(measure,sum(abs(w))) for at most
// kMaximumProductPoints weights. Any non-finite weight fails.
bool reproduces_measure(std::span<const double> weights, double measure)
{
    double sum = 0.0;
    double correction = 0.0;
    double scale = 0.0;
    for (const double weight : weights) {
        const double term = weight * kValidationSumScale;
        // Neumaier compensation keeps long signed sums accurate.
        const double next = sum + term;
        correction += std::abs(sum) >= std::abs(term) ? (sum - next) + term : (term - next) + sum;
        sum = next;
        scale += std::abs(term);
    }
    sum += correction;
    const double scaled_measure = measure * kValidationSumScale;
    return std::isfinite(sum) && std::isfinite(scale) &&
        std::abs(sum - scaled_measure) <= kValidationTolerance * std::max(scaled_measure, scale);
}

void validate_line_rule(const QuadratureRule& rule, int required_exactness, const char* name)
{
    const std::string prefix = std::string{"Triangle Duffy generator: "} + name;
    svmp::check<InvalidArgumentException>(rule.cell_family() == svmp::CellFamily::Line,
        prefix + " must have canonical Line family");
    svmp::check<InvalidArgumentException>(rule.polynomial_exactness() >= required_exactness,
        prefix + " exactness is below the directional requirement");
    for (const QuadPoint& point : rule.points()) {
        if (std::abs(point[0]) > 1.0) {
            svmp::raise<InvalidArgumentException>(prefix + " contains a point outside [-1,1]");
        }
    }
    svmp::check<InvalidArgumentException>(reproduces_measure(rule.weights(), 2.0),
        prefix + " normalization does not reproduce measure 2");
}

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
    validate_line_rule(u_rule, requested_exactness + 1, "u_rule");
    validate_line_rule(v_rule, requested_exactness, "v_rule");

    // Factor points in [-1,1] give u,v in [0,1], so (u,(1-u)v) lies in the
    // triangle up to one rounding of 1-u.
    std::vector<QuadPoint> points(n_u * n_v, QuadPoint::Zero());
    std::vector<double> weights(points.size());
    for (std::size_t i = 0; i < n_u; ++i) {
        const double u = (u_rule.point(i)[0] + 1.0) / 2.0;
        const double remaining = 1.0 - u;
        // Two interval factors of 1/2 and the Duffy Jacobian (1-u).
        const double outer_weight = u_rule.weight(i) * (remaining / 4.0);
        for (std::size_t j = 0; j < n_v; ++j) {
            const std::size_t q = i * n_v + j;
            const double v = (v_rule.point(j)[0] + 1.0) / 2.0;
            points[q] = QuadPoint{u, remaining * v, 0.0};
            weights[q] = outer_weight * v_rule.weight(j);
        }
    }
    // Given both factor measures, this checks the u_rule first moment and
    // rejects overflowed products.
    svmp::check<ConvergenceException>(reproduces_measure(weights, 0.5),
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
