// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "FE/Quadrature/GaussLobattoNodes.h"

#include "FE/Common/FEException.h"

#include <cmath>
#include <numbers>
#include <sstream>
#include <string_view>
#include <utility>

namespace svmp::FE::quadrature::detail {
namespace {

// Bound cosine-seeded Newton refinement for deterministic termination.
constexpr int kMaximumNewtonIterations = 100;

// Return (P_degree(x), P_(degree-1)(x)) for degree >= 1 using the Legendre
// three-term recurrence, starting from P_0 = 1 and P_1 = x.
std::pair<double, double> evaluate_adjacent_legendre_values(
    std::size_t degree,
    double coordinate) noexcept
{
    double previous_value = 1.0;
    double value = coordinate;

    for (std::size_t recurrence_degree = 2; recurrence_degree <= degree; ++recurrence_degree) {
        const double next_value =
            ((2 * recurrence_degree - 1) * coordinate * value -
             (recurrence_degree - 1) * previous_value) / recurrence_degree;
        previous_value = value;
        value = next_value;
    }

    return {value, previous_value};
}

// Preserve the failed quantity and generator context in a convergence error.
[[noreturn]] void raise_generation_failure(std::size_t num_points, std::size_t half_root_index,
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

} // namespace

// Refine f(x) = x*P_(n-1)(x) - P_(n-2)(x), then recheck its correction at
// the returned node. Callers retain their existing correction tolerances.
GaussLobattoNode refine_gauss_lobatto_node(std::size_t num_points,
    std::size_t half_root_index, double correction_tolerance)
{
    if (num_points < 3 || half_root_index >= (num_points - 1) / 2) {
        svmp::raise<InvalidArgumentException>(
            "Gauss-Lobatto-Legendre nodes: invalid interior node request");
    }
    if (!(std::isfinite(correction_tolerance) && correction_tolerance > 0.0)) {
        svmp::raise<InvalidArgumentException>(
            "Gauss-Lobatto-Legendre nodes: correction tolerance must be finite and positive");
    }

    const bool is_center =
        1u + half_root_index == num_points - 2u - half_root_index;
    const std::size_t polynomial_degree = num_points - 1;
    const double num_points_value = num_points;
    const double degree_value = polynomial_degree;
    const double pi = std::numbers::pi_v<double>;
    double root = std::cos(pi * (half_root_index + 1) / degree_value);
    double correction = 0.0;

    // For f = x*P_m - P_(m-1), Legendre identities give
    // f' = (m+1)*P_m = n*P_m exactly.
    for (int iteration = 1; iteration <= kMaximumNewtonIterations; ++iteration) {
        const auto [polynomial_value, previous_polynomial_value] =
            evaluate_adjacent_legendre_values(polynomial_degree, root);
        if (!(std::isfinite(polynomial_value) &&
              std::isfinite(previous_polynomial_value))) {
            raise_generation_failure(
                num_points, half_root_index, iteration, polynomial_value,
                "encountered invalid adjacent Legendre values");
        }

        const double residual =
            root * polynomial_value - previous_polynomial_value;
        const double derivative = num_points_value * polynomial_value;
        if (!(std::isfinite(residual) &&
              std::isfinite(derivative) &&
              derivative != 0.0)) {
            raise_generation_failure(
                num_points, half_root_index, iteration, derivative,
                "computed an invalid root-function residual or derivative");
        }

        correction = residual / derivative;
        const double updated_root = root - correction;
        if (!(std::isfinite(correction) && std::isfinite(updated_root))) {
            raise_generation_failure(
                num_points, half_root_index, iteration, correction,
                "computed an invalid Newton update");
        }
        root = updated_root;

        if (std::abs(correction) > correction_tolerance) {
            continue;
        }

        if (is_center) {
            root = 0.0;
        }
        if (!(root >= 0.0 && root < 1.0 && (is_center || root > 0.0))) {
            raise_generation_failure(
                num_points, half_root_index, iteration, root,
                "refined root is outside the expected half interval");
        }

        const auto [final_polynomial_value,
                    final_previous_polynomial_value] =
            evaluate_adjacent_legendre_values(polynomial_degree, root);
        if (!(std::isfinite(final_polynomial_value) &&
              std::isfinite(final_previous_polynomial_value))) {
            raise_generation_failure(
                num_points, half_root_index, iteration, final_polynomial_value,
                "refined root produced invalid adjacent Legendre values");
        }

        const double final_residual =
            root * final_polynomial_value -
            final_previous_polynomial_value;
        const double final_derivative =
            num_points_value * final_polynomial_value;
        if (!(std::isfinite(final_residual) &&
              std::isfinite(final_derivative) &&
              final_derivative != 0.0)) {
            raise_generation_failure(
                num_points, half_root_index, iteration, final_derivative,
                "refined root produced an invalid residual or derivative");
        }

        const double final_correction =
            final_residual / final_derivative;
        if (!(std::isfinite(final_correction) &&
              std::abs(final_correction) <=
                  correction_tolerance)) {
            raise_generation_failure(
                num_points, half_root_index, iteration, final_correction,
                "refined root failed final correction validation");
        }

        return {root, final_polynomial_value, iteration};
    }

    raise_generation_failure(
        num_points, half_root_index, kMaximumNewtonIterations, correction,
        "Newton refinement did not converge");
}

} // namespace svmp::FE::quadrature::detail
