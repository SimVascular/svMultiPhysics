// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file GaussLobattoQuadrature.cpp
 * @brief Exactness-requested generation of bounded Gauss-Lobatto-Legendre line rules.
 * @ingroup FE_Quadrature
 */

#include "FE/Quadrature/GaussLobattoQuadrature.h"

#include "FE/Common/FEException.h"
#include "FE/Quadrature/GaussLobattoNodes.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace svmp::FE::quadrature {
namespace {

constexpr std::size_t kMaximumPoints = 128;
static_assert(max_gauss_lobatto_exactness() == 2 * kMaximumPoints - 3);

// Guard recurrence and Newton-update rounding; exhaustive supported-size
// sweeps qualify this scale.
constexpr double kNewtonCorrectionTolerance =
    64.0 * std::numeric_limits<double>::epsilon();
// Provide conservative O(n epsilon) accumulation headroom, qualified by those
// sweeps.
constexpr double kRuleValidationTolerance =
    32.0 * kMaximumPoints *
    std::numeric_limits<double>::epsilon();

// Preserve the failed quantity and generator context in a convergence error.
// A half-root index or iteration of -1 identifies a rule-wide validation check.
[[noreturn]] void raise_generation_failure(std::size_t num_points, int half_root_index,
    int iteration, double diagnostic_value, std::string_view detail)
{
    std::ostringstream message;
    message << "Gauss-Lobatto-Legendre generator: " << detail
            << ", num_points=" << num_points
            << ", half_root_index=" << half_root_index
            << ", diagnostic_value=" << diagnostic_value;

    svmp::raise<ConvergenceException>(
        message.str(), iteration, std::abs(diagnostic_value));
}

// Form a weight from the polynomial value already computed during refinement.
std::pair<double, double> generate_interior_root_and_weight(std::size_t num_points,
    std::size_t half_root_index, double weight_denominator_scale)
{
    const auto node = detail::refine_gauss_lobatto_node(
        num_points, half_root_index, kNewtonCorrectionTolerance);
    const double denominator =
        weight_denominator_scale * node.legendre_value * node.legendre_value;
    if (!(std::isfinite(denominator) && denominator > 0.0)) {
        raise_generation_failure(
            num_points, half_root_index, node.iterations, denominator,
            "refined root produced an invalid weight denominator");
    }

    const double weight = 2.0 / denominator;
    if (!(std::isfinite(weight) && weight > 0.0)) {
        raise_generation_failure(
            num_points, half_root_index, node.iterations, weight,
            "refined root produced an invalid quadrature weight");
    }
    return {node.coordinate, weight};
}

} // namespace

QuadratureRule make_gauss_lobatto_rule(int requested_exactness)
{
    svmp::check<InvalidArgumentException>(
        requested_exactness >= 0 &&
            requested_exactness <= max_gauss_lobatto_exactness(),
        "Gauss-Lobatto-Legendre generator: requested_exactness must be in [0, " +
            std::to_string(max_gauss_lobatto_exactness()) + ']');

    const std::size_t num_points = requested_exactness / 2 + 2;
    std::vector<QuadPoint> points(num_points, QuadPoint::Zero());
    std::vector<double> weights(points.size());

    points.front()[0] = -1.0;
    points.back()[0] = 1.0;

    const double weight_denominator_scale =
        num_points * (num_points - 1);
    const double endpoint_weight = 2.0 / weight_denominator_scale;
    weights.front() = endpoint_weight;
    weights.back() = endpoint_weight;

    const std::size_t interior_roots_to_refine = (num_points - 1) / 2;
    for (std::size_t half_root_index = 0;
         half_root_index < interior_roots_to_refine; ++half_root_index) {
        const std::size_t left_index = 1u + half_root_index;
        const std::size_t right_index =
            points.size() - 2u - half_root_index;
        const auto [root, weight] =
            generate_interior_root_and_weight(
                num_points, half_root_index, weight_denominator_scale);

        points[left_index][0] = -root;
        points[right_index][0] = root;
        weights[left_index] = weight;
        weights[right_index] = weight;
    }

    for (std::size_t point_index = 1; point_index < points.size(); ++point_index) {
        const double spacing = points[point_index][0] - points[point_index - 1u][0];
        if (!(spacing > 0.0)) {
            raise_generation_failure(
                num_points, point_index, -1, spacing,
                "generated points are not strictly increasing");
        }
    }

    // Report a failed measure instead of repairing or rescaling the weights.
    const long double weight_sum = std::accumulate(weights.begin(), weights.end(), 0.0L);
    const long double measure_error = std::abs(weight_sum - 2.0L);
    if (!(std::isfinite(weight_sum) &&
          measure_error <= kRuleValidationTolerance)) {
        raise_generation_failure(
            num_points, -1, -1, measure_error,
            "generated weights do not reproduce the reference measure");
    }

    const int polynomial_exactness = 2 * num_points - 3;
    return QuadratureRule(
        svmp::CellFamily::Line, polynomial_exactness,
        std::move(points), std::move(weights));
}

} // namespace svmp::FE::quadrature
