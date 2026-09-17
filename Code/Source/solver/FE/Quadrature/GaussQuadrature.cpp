// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file GaussQuadrature.cpp
 * @brief Exactness-requested generation of bounded Gauss-Legendre line rules.
 * @ingroup FE_Quadrature
 */

#include "FE/Quadrature/GaussQuadrature.h"

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
static_assert(max_gauss_legendre_exactness() == 2 * kMaximumPoints - 1);

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

std::pair<double, double> evaluate_legendre_with_derivative(
    int degree,
    double coordinate) noexcept
{
    double previous_value = 1.0;
    double previous_derivative = 0.0;
    double value = coordinate;
    double derivative = 1.0;
    for (int recurrence_degree = 2; recurrence_degree <= degree; ++recurrence_degree) {
        const double degree_value = static_cast<double>(recurrence_degree);
        const double recurrence_factor = static_cast<double>(2 * recurrence_degree - 1);
        const double next_value =
            (recurrence_factor * coordinate * value -
             static_cast<double>(recurrence_degree - 1) * previous_value) / degree_value;
        const double next_derivative =
            (recurrence_factor * (value + coordinate * derivative) -
             static_cast<double>(recurrence_degree - 1) * previous_derivative) / degree_value;

        previous_value = value;
        previous_derivative = derivative;
        value = next_value;
        derivative = next_derivative;
    }

    return {value, derivative};
}

[[noreturn]] void raise_generation_failure(
    int num_points,
    int root_index,
    int iteration,
    double diagnostic_value,
    std::string_view detail)
{
    std::ostringstream message;
    message << "Gauss-Legendre generator: " << detail
            << ", num_points=" << num_points
            << ", root_index=" << root_index
            << ", diagnostic_value=" << diagnostic_value;

    const double residual = std::isfinite(diagnostic_value)
                                ? std::abs(diagnostic_value)
                                : 0.0;
    svmp::raise<ConvergenceException>(message.str(), iteration, residual);
}

void require_generation(bool condition, int num_points, int root_index, int iteration,
    double diagnostic_value, std::string_view detail)
{
    if (!condition) {
        raise_generation_failure(num_points, root_index, iteration, diagnostic_value, detail);
    }
}

std::pair<double, double> generate_root_and_weight(int num_points, int root_index, bool is_center)
{
    const double pi = std::numbers::pi_v<double>;
    double root = std::cos(
        pi * (static_cast<double>(root_index) + 0.75) /
        (static_cast<double>(num_points) + 0.5));
    double correction = 0.0;

    for (int iteration = 1; iteration <= kMaximumNewtonIterations; ++iteration) {
        const auto [polynomial_value, polynomial_derivative] =
            evaluate_legendre_with_derivative(num_points, root);
        require_generation(
            std::isfinite(polynomial_value) &&
                std::isfinite(polynomial_derivative) &&
                polynomial_derivative != 0.0,
            num_points, root_index, iteration, polynomial_derivative,
            "encountered an invalid Legendre value or derivative");

        correction = polynomial_value / polynomial_derivative;
        const double updated_root = root - correction;
        require_generation(
            std::isfinite(correction) && std::isfinite(updated_root),
            num_points, root_index, iteration, correction,
            "computed an invalid Newton update");
        root = updated_root;

        if (std::abs(correction) > kNewtonCorrectionTolerance) {
            continue;
        }

        if (is_center) {
            root = 0.0;
        }
        require_generation(
            root >= 0.0 && root < 1.0 &&
                (is_center || root > 0.0),
            num_points, root_index, iteration, root,
            "refined root is outside the expected half interval");

        const auto [final_polynomial_value,
                    final_polynomial_derivative] =
            evaluate_legendre_with_derivative(num_points, root);
        require_generation(
            std::isfinite(final_polynomial_value) &&
                std::isfinite(final_polynomial_derivative) &&
                final_polynomial_derivative != 0.0,
            num_points, root_index, iteration, final_polynomial_derivative,
            "refined root produced an invalid Legendre value or derivative");

        const double final_correction =
            final_polynomial_value / final_polynomial_derivative;
        require_generation(
            std::isfinite(final_correction) &&
                std::abs(final_correction) <=
                    kNewtonCorrectionTolerance,
            num_points, root_index, iteration, final_correction,
            "refined root failed final correction validation");

        const double denominator =
            (1.0 - root) * (1.0 + root) *
            final_polynomial_derivative * final_polynomial_derivative;
        require_generation(
            std::isfinite(denominator) && denominator > 0.0,
            num_points, root_index, iteration, denominator,
            "refined root produced an invalid weight denominator");

        const double weight = 2.0 / denominator;
        require_generation(
            std::isfinite(weight) && weight > 0.0,
            num_points, root_index, iteration, weight,
            "refined root produced an invalid quadrature weight");

        return {root, weight};
    }

    raise_generation_failure(
        num_points, root_index, kMaximumNewtonIterations, correction,
        "Newton refinement did not converge");
}

} // namespace

QuadratureRule make_gauss_legendre_rule(int requested_exactness)
{
    if (requested_exactness < 0 || requested_exactness > max_gauss_legendre_exactness()) {
        std::ostringstream message;
        message << "Gauss-Legendre generator: requested_exactness must be in [0, "
                << max_gauss_legendre_exactness() << ']';
        svmp::raise<InvalidArgumentException>(message.str());
    }

    const int num_points = requested_exactness / 2 + 1;
    std::vector<QuadPoint> points(
        static_cast<std::size_t>(num_points), QuadPoint::Zero());
    std::vector<double> weights(points.size());

    const int roots_to_refine = (num_points + 1) / 2;
    for (int root_index = 0; root_index < roots_to_refine; ++root_index) {
        const std::size_t left_index = static_cast<std::size_t>(root_index);
        const std::size_t right_index = points.size() - 1u - left_index;
        const auto [root, weight] = generate_root_and_weight(
            num_points, root_index, left_index == right_index);

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

    const int polynomial_exactness = 2 * num_points - 1;
    return QuadratureRule(
        svmp::CellFamily::Line, polynomial_exactness,
        std::move(points), std::move(weights));
}

} // namespace svmp::FE::quadrature
