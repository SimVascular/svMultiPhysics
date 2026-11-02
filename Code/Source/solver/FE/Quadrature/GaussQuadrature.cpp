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

constexpr std::size_t kMaximumPoints = 128;
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
    32.0 * kMaximumPoints *
    std::numeric_limits<double>::epsilon();

// Return (P_degree(x), P'_degree(x)) for degree >= 1. Advance the Legendre
// three-term recurrence and its derivative together, starting from P_0 = 1
// and P_1 = x, without numerical differentiation.
std::pair<double, double> evaluate_legendre_with_derivative(
    std::size_t degree,
    double coordinate) noexcept
{
    double previous_value = 1.0;
    double previous_derivative = 0.0;
    double value = coordinate;
    double derivative = 1.0;
    for (std::size_t recurrence_degree = 2; recurrence_degree <= degree; ++recurrence_degree) {
        const double degree_value = recurrence_degree;
        const double recurrence_factor = 2 * recurrence_degree - 1;
        const double next_value =
            (recurrence_factor * coordinate * value -
             (recurrence_degree - 1) * previous_value) / degree_value;
        const double next_derivative =
            (recurrence_factor * (value + coordinate * derivative) -
             (recurrence_degree - 1) * previous_derivative) / degree_value;

        previous_value = value;
        previous_derivative = derivative;
        value = next_value;
        derivative = next_derivative;
    }

    return {value, derivative};
}

// Preserve the failed quantity and generator context in a convergence error.
// A root index or iteration of -1 identifies a rule-wide validation check.
[[noreturn]] void raise_generation_failure(
    std::size_t num_points,
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

    svmp::raise<ConvergenceException>(
        message.str(), iteration, std::abs(diagnostic_value));
}

// Refine a nonnegative root of P_n with cosine-seeded Newton iteration, bounded
// by the iteration and correction limits above. Recheck P_n/P'_n before forming
// w = 2 / ((1 - x*x) * P'_n(x)^2). The caller mirrors (x, w); is_center assigns
// the odd rule's center exactly to zero before the final validation.
std::pair<double, double> generate_root_and_weight(
    std::size_t num_points, std::size_t root_index, bool is_center)
{
    const double pi = std::numbers::pi_v<double>;
    double root = std::cos(
        pi * (root_index + 0.75) / (num_points + 0.5));
    double correction = 0.0;

    for (int iteration = 1; iteration <= kMaximumNewtonIterations; ++iteration) {
        const auto [polynomial_value, polynomial_derivative] =
            evaluate_legendre_with_derivative(num_points, root);
        if (!(std::isfinite(polynomial_value) &&
              std::isfinite(polynomial_derivative) &&
              polynomial_derivative != 0.0)) {
            raise_generation_failure(
                num_points, root_index, iteration, polynomial_derivative,
                "encountered an invalid Legendre value or derivative");
        }

        correction = polynomial_value / polynomial_derivative;
        const double updated_root = root - correction;
        if (!(std::isfinite(correction) && std::isfinite(updated_root))) {
            raise_generation_failure(
                num_points, root_index, iteration, correction,
                "computed an invalid Newton update");
        }
        root = updated_root;

        if (std::abs(correction) > kNewtonCorrectionTolerance) {
            continue;
        }

        if (is_center) {
            root = 0.0;
        }
        if (!(root >= 0.0 && root < 1.0 &&
              (is_center || root > 0.0))) {
            raise_generation_failure(
                num_points, root_index, iteration, root,
                "refined root is outside the expected half interval");
        }

        const auto [final_polynomial_value,
                    final_polynomial_derivative] =
            evaluate_legendre_with_derivative(num_points, root);
        if (!(std::isfinite(final_polynomial_value) &&
              std::isfinite(final_polynomial_derivative) &&
              final_polynomial_derivative != 0.0)) {
            raise_generation_failure(
                num_points, root_index, iteration, final_polynomial_derivative,
                "refined root produced an invalid Legendre value or derivative");
        }

        const double final_correction =
            final_polynomial_value / final_polynomial_derivative;
        if (!(std::isfinite(final_correction) &&
              std::abs(final_correction) <=
                  kNewtonCorrectionTolerance)) {
            raise_generation_failure(
                num_points, root_index, iteration, final_correction,
                "refined root failed final correction validation");
        }

        const double denominator =
            (1.0 - root) * (1.0 + root) *
            final_polynomial_derivative * final_polynomial_derivative;
        if (!(std::isfinite(denominator) && denominator > 0.0)) {
            raise_generation_failure(
                num_points, root_index, iteration, denominator,
                "refined root produced an invalid weight denominator");
        }

        const double weight = 2.0 / denominator;
        if (!(std::isfinite(weight) && weight > 0.0)) {
            raise_generation_failure(
                num_points, root_index, iteration, weight,
                "refined root produced an invalid quadrature weight");
        }

        return {root, weight};
    }

    raise_generation_failure(
        num_points, root_index, kMaximumNewtonIterations, correction,
        "Newton refinement did not converge");
}

} // namespace

QuadratureRule make_gauss_legendre_rule(int requested_exactness)
{
    svmp::check<InvalidArgumentException>(
        requested_exactness >= 0 &&
            requested_exactness <= max_gauss_legendre_exactness(),
        "Gauss-Legendre generator: requested_exactness must be in [0, " +
            std::to_string(max_gauss_legendre_exactness()) + ']');

    const std::size_t num_points = requested_exactness / 2 + 1;
    std::vector<QuadPoint> points(num_points, QuadPoint::Zero());
    std::vector<double> weights(points.size());

    const std::size_t roots_to_refine = (num_points + 1) / 2;
    for (std::size_t left_index = 0; left_index < roots_to_refine; ++left_index) {
        const std::size_t right_index = points.size() - 1u - left_index;
        const auto [root, weight] = generate_root_and_weight(
            num_points, left_index, left_index == right_index);

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

    const int polynomial_exactness = 2 * num_points - 1;
    return QuadratureRule(
        svmp::CellFamily::Line, polynomial_exactness,
        std::move(points), std::move(weights));
}

} // namespace svmp::FE::quadrature
