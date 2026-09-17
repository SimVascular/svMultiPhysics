// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file GaussLobattoQuadrature.cpp
 * @brief Exactness-requested generation of bounded Gauss-Lobatto-Legendre line rules.
 * @ingroup FE_Quadrature
 */

#include "FE/Quadrature/GaussLobattoQuadrature.h"

#include "FE/Common/FEException.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <numeric>
#include <sstream>
#include <string_view>
#include <utility>
#include <vector>

namespace svmp::FE::quadrature {
namespace {

constexpr int kMaximumPoints = 128;
static_assert(max_gauss_lobatto_exactness() == 2 * kMaximumPoints - 3);

// Defensively bound supported cosine-seeded Newton refinements for
// deterministic termination.
constexpr int kMaximumNewtonIterations = 100;
// Guard recurrence and Newton-update rounding; exhaustive supported-size
// sweeps qualify this scale.
constexpr double kNewtonCorrectionTolerance =
    64.0 * std::numeric_limits<double>::epsilon();
// Provide conservative O(n epsilon) accumulation headroom, qualified by those
// sweeps.
constexpr double kRuleValidationTolerance =
    32.0 * static_cast<double>(kMaximumPoints) *
    std::numeric_limits<double>::epsilon();

std::pair<double, double> evaluate_adjacent_legendre_values(
    int degree,
    double coordinate) noexcept
{
    double previous_value = 1.0;
    double value = coordinate;

    for (int recurrence_degree = 2; recurrence_degree <= degree; ++recurrence_degree) {
        const double next_value =
            (static_cast<double>(2 * recurrence_degree - 1) *
                 coordinate * value -
             static_cast<double>(recurrence_degree - 1) * previous_value) /
            static_cast<double>(recurrence_degree);
        previous_value = value;
        value = next_value;
    }

    return {value, previous_value};
}

[[noreturn]] void raise_generation_failure(int num_points, int half_root_index,
    int iteration, double diagnostic_value, std::string_view detail)
{
    std::ostringstream message;
    message << "Gauss-Lobatto-Legendre generator: " << detail
            << ", num_points=" << num_points
            << ", half_root_index=" << half_root_index
            << ", diagnostic_value=" << diagnostic_value;

    const double residual = std::isfinite(diagnostic_value)
                                ? std::abs(diagnostic_value)
                                : 0.0;
    svmp::raise<ConvergenceException>(message.str(), iteration, residual);
}

void require_generation(bool condition, int num_points, int half_root_index,
    int iteration, double diagnostic_value, std::string_view detail)
{
    if (!condition) {
        raise_generation_failure(
            num_points, half_root_index, iteration, diagnostic_value, detail);
    }
}

std::pair<double, double> generate_interior_root_and_weight(int num_points, int half_root_index,
    bool is_center, double weight_denominator_scale)
{
    const int polynomial_degree = num_points - 1;
    const double num_points_value = static_cast<double>(num_points);
    const double degree_value = static_cast<double>(polynomial_degree);
    const double pi = std::numbers::pi_v<double>;
    double root = std::cos(pi * static_cast<double>(half_root_index + 1) / degree_value);
    double correction = 0.0;

    // For f = x*P_m - P_(m-1), Legendre identities give
    // f' = (m+1)*P_m = n*P_m exactly.
    for (int iteration = 1; iteration <= kMaximumNewtonIterations; ++iteration) {
        const auto [polynomial_value, previous_polynomial_value] =
            evaluate_adjacent_legendre_values(polynomial_degree, root);
        require_generation(
            std::isfinite(polynomial_value) &&
                std::isfinite(previous_polynomial_value),
            num_points, half_root_index, iteration,
            polynomial_value,
            "encountered invalid adjacent Legendre values");

        const double residual =
            root * polynomial_value - previous_polynomial_value;
        const double derivative = num_points_value * polynomial_value;
        require_generation(
            std::isfinite(residual) &&
                std::isfinite(derivative) &&
                derivative != 0.0,
            num_points, half_root_index, iteration, derivative,
            "computed an invalid root-function residual or derivative");

        correction = residual / derivative;
        const double updated_root = root - correction;
        require_generation(
            std::isfinite(correction) && std::isfinite(updated_root),
            num_points, half_root_index, iteration, correction,
            "computed an invalid Newton update");
        root = updated_root;

        if (std::abs(correction) > kNewtonCorrectionTolerance) {
            continue;
        }

        if (is_center) {
            root = 0.0;
        }
        require_generation(
            root >= 0.0 && root < 1.0 && (is_center || root > 0.0),
            num_points, half_root_index, iteration, root,
            "refined root is outside the expected half interval");

        const auto [final_polynomial_value,
                    final_previous_polynomial_value] =
            evaluate_adjacent_legendre_values(polynomial_degree, root);
        require_generation(
            std::isfinite(final_polynomial_value) &&
                std::isfinite(final_previous_polynomial_value),
            num_points, half_root_index, iteration,
            final_polynomial_value,
            "refined root produced invalid adjacent Legendre values");

        const double final_residual =
            root * final_polynomial_value -
            final_previous_polynomial_value;
        const double final_derivative =
            num_points_value * final_polynomial_value;
        require_generation(
            std::isfinite(final_residual) &&
                std::isfinite(final_derivative) &&
                final_derivative != 0.0,
            num_points, half_root_index, iteration, final_derivative,
            "refined root produced an invalid residual or derivative");

        const double final_correction =
            final_residual / final_derivative;
        require_generation(
            std::isfinite(final_correction) &&
                std::abs(final_correction) <=
                    kNewtonCorrectionTolerance,
            num_points, half_root_index, iteration, final_correction,
            "refined root failed final correction validation");

        const double denominator =
            weight_denominator_scale * final_polynomial_value *
            final_polynomial_value;
        require_generation(
            std::isfinite(denominator) && denominator > 0.0,
            num_points, half_root_index, iteration, denominator,
            "refined root produced an invalid weight denominator");

        const double weight = 2.0 / denominator;
        require_generation(
            std::isfinite(weight) && weight > 0.0,
            num_points, half_root_index, iteration, weight,
            "refined root produced an invalid quadrature weight");

        return {root, weight};
    }

    raise_generation_failure(
        num_points, half_root_index, kMaximumNewtonIterations, correction,
        "Newton refinement did not converge");
}

} // namespace

QuadratureRule make_gauss_lobatto_rule(int requested_exactness)
{
    if (requested_exactness < 0 || requested_exactness > max_gauss_lobatto_exactness()) {
        std::ostringstream message;
        message << "Gauss-Lobatto-Legendre generator: "
                << "requested_exactness must be in [0, "
                << max_gauss_lobatto_exactness() << ']';
        svmp::raise<InvalidArgumentException>(message.str());
    }

    const int num_points = requested_exactness / 2 + 2;
    std::vector<QuadPoint> points(
        static_cast<std::size_t>(num_points), QuadPoint::Zero());
    std::vector<double> weights(points.size());

    points.front()[0] = -1.0;
    points.back()[0] = 1.0;

    const double weight_denominator_scale =
        static_cast<double>(num_points * (num_points - 1));
    const double endpoint_weight = 2.0 / weight_denominator_scale;
    weights.front() = endpoint_weight;
    weights.back() = endpoint_weight;

    const int interior_roots_to_refine = (num_points - 1) / 2;
    for (int half_root_index = 0;
         half_root_index < interior_roots_to_refine; ++half_root_index) {
        const std::size_t left_index =
            1u + static_cast<std::size_t>(half_root_index);
        const std::size_t right_index =
            points.size() - 2u -
            static_cast<std::size_t>(half_root_index);
        const auto [root, weight] =
            generate_interior_root_and_weight(
                num_points, half_root_index, left_index == right_index,
                weight_denominator_scale);

        points[left_index][0] = -root;
        points[right_index][0] = root;
        weights[left_index] = weight;
        weights[right_index] = weight;
    }

    for (std::size_t point_index = 1; point_index < points.size(); ++point_index) {
        const double spacing = points[point_index][0] - points[point_index - 1u][0];
        require_generation(
            spacing > 0.0,
            num_points, static_cast<int>(point_index), -1, spacing,
            "generated points are not strictly increasing");
    }

    // Report a failed measure instead of repairing or rescaling the weights.
    const long double weight_sum = std::accumulate(weights.begin(), weights.end(), 0.0L);
    const long double measure_error = std::abs(weight_sum - 2.0L);
    require_generation(
        std::isfinite(weight_sum) &&
            measure_error <=
                static_cast<long double>(kRuleValidationTolerance),
        num_points, -1, -1, static_cast<double>(measure_error),
        "generated weights do not reproduce the reference measure");

    const int polynomial_exactness = 2 * num_points - 3;
    return QuadratureRule(
        svmp::CellFamily::Line, polynomial_exactness,
        std::move(points), std::move(weights));
}

} // namespace svmp::FE::quadrature
