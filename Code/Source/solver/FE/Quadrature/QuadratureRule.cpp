// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file QuadratureRule.cpp
 * @brief Internal construction and structural validation for quadrature rules.
 * @ingroup FE_Quadrature
 */

#include "FE/Quadrature/QuadratureRule.h"

#include "FE/Common/FEException.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

namespace svmp::FE::quadrature {
namespace {

enum class ReferenceDomain {
    Unsupported,
    Point,
    Line,
    Triangle,
    Quad,
    Tetra,
    Hex,
    Wedge,
};

struct ReferenceCellTraits {
    ReferenceDomain domain;
    int dimension;
    double measure;
};

enum class ValidationFailure {
    None,
    UnsupportedCellFamily,
    NegativePolynomialExactness,
    InvalidTolerance,
    EmptyRule,
    PointWeightSizeMismatch,
    NonFiniteCoordinate,
    NonzeroInactiveCoordinate,
    PointOutsideReferenceCell,
    NonFiniteWeight,
    NonFiniteWeightAccumulation,
    IncorrectZerothMoment,
    UnstableStoredOrderAccumulation,
    IllConditionedWeightCancellation,
};

struct ValidationResult {
    ValidationFailure failure{ValidationFailure::None};
    std::size_t sample{std::numeric_limits<std::size_t>::max()};
};

constexpr ReferenceCellTraits unsupported_traits() noexcept
{
    return {
        ReferenceDomain::Unsupported,
        -1,
        std::numeric_limits<double>::quiet_NaN(),
    };
}

constexpr ReferenceCellTraits reference_cell_traits(svmp::CellFamily family) noexcept
{
    switch (family) {
        case svmp::CellFamily::Point:
            return {ReferenceDomain::Point, 0, 1.0};
        case svmp::CellFamily::Line:
            return {ReferenceDomain::Line, 1, 2.0};
        case svmp::CellFamily::Triangle:
            return {ReferenceDomain::Triangle, 2, 0.5};
        case svmp::CellFamily::Quad:
            return {ReferenceDomain::Quad, 2, 4.0};
        case svmp::CellFamily::Tetra:
            return {ReferenceDomain::Tetra, 3, 1.0 / 6.0};
        case svmp::CellFamily::Hex:
            return {ReferenceDomain::Hex, 3, 8.0};
        case svmp::CellFamily::Wedge:
            return {ReferenceDomain::Wedge, 3, 1.0};
        default:
            return unsupported_traits();
    }
}

ValidationResult validate_point(
    const QuadPoint& point,
    const ReferenceCellTraits& traits,
    double tolerance,
    std::size_t sample) noexcept
{
    for (std::size_t component = 0; component < 3u; ++component) {
        if (!std::isfinite(point[component])) {
            return {ValidationFailure::NonFiniteCoordinate, sample};
        }
        if (component >= static_cast<std::size_t>(traits.dimension) &&
            std::abs(point[component]) > tolerance) {
            return {ValidationFailure::NonzeroInactiveCoordinate, sample};
        }
    }

    const auto in_interval = [tolerance](double value, double lower, double upper) {
        return value >= lower - tolerance && value <= upper + tolerance;
    };

    const double x = point[0];
    const double y = point[1];
    const double z = point[2];

    bool is_contained = false;
    switch (traits.domain) {
        case ReferenceDomain::Point:
            is_contained = true;
            break;
        case ReferenceDomain::Line:
            is_contained = in_interval(x, -1.0, 1.0);
            break;
        case ReferenceDomain::Triangle:
            is_contained = x >= -tolerance && y >= -tolerance &&
                           x + y <= 1.0 + tolerance;
            break;
        case ReferenceDomain::Quad:
            is_contained = in_interval(x, -1.0, 1.0) &&
                           in_interval(y, -1.0, 1.0);
            break;
        case ReferenceDomain::Tetra:
            is_contained = x >= -tolerance && y >= -tolerance &&
                           z >= -tolerance &&
                           x + y + z <= 1.0 + tolerance;
            break;
        case ReferenceDomain::Hex:
            is_contained = in_interval(x, -1.0, 1.0) &&
                           in_interval(y, -1.0, 1.0) &&
                           in_interval(z, -1.0, 1.0);
            break;
        case ReferenceDomain::Wedge:
            is_contained = x >= -tolerance && y >= -tolerance &&
                           x + y <= 1.0 + tolerance &&
                           in_interval(z, -1.0, 1.0);
            break;
        case ReferenceDomain::Unsupported:
            break;
    }

    if (!is_contained) {
        return {ValidationFailure::PointOutsideReferenceCell, sample};
    }
    return {};
}

double add_as_binary64(double left, double right) noexcept
{
    volatile double rounded_sum = left + right;
    return rounded_sum;
}

ValidationResult validate_weights(
    const std::vector<double>& weights,
    double reference_measure,
    double tolerance) noexcept
{
    double ordered_sum = 0.0;
    long double compensated_sum = 0.0L;
    long double correction = 0.0L;
    long double absolute_sum = 0.0L;
    bool has_negative_weight = false;

    for (std::size_t sample = 0; sample < weights.size(); ++sample) {
        const double weight = weights[sample];
        if (!std::isfinite(weight)) {
            return {ValidationFailure::NonFiniteWeight, sample};
        }

        ordered_sum = add_as_binary64(ordered_sum, weight);
        if (!std::isfinite(ordered_sum)) {
            return {ValidationFailure::NonFiniteWeightAccumulation, sample};
        }
        has_negative_weight = has_negative_weight || weight < 0.0;

        const long double value = static_cast<long double>(weight);
        const long double next = compensated_sum + value;
        if (std::abs(compensated_sum) >= std::abs(value)) {
            correction += (compensated_sum - next) + value;
        } else {
            correction += (value - next) + compensated_sum;
        }
        compensated_sum = next;
        absolute_sum += std::abs(value);
        if (!std::isfinite(compensated_sum) ||
            !std::isfinite(correction) ||
            !std::isfinite(absolute_sum)) {
            return {ValidationFailure::NonFiniteWeightAccumulation, sample};
        }
    }

    const long double total = compensated_sum + correction;
    if (!std::isfinite(total)) {
        return {ValidationFailure::NonFiniteWeightAccumulation, weights.size() - 1u};
    }

    const long double expected = static_cast<long double>(reference_measure);
    const long double scale = std::max(1.0L, std::abs(expected));
    const long double requested_error_budget =
        static_cast<long double>(tolerance) * scale;
    const long double compensated_error = std::abs(total - expected);

    if (compensated_error > requested_error_budget) {
        return {ValidationFailure::IncorrectZerothMoment};
    }

    const long double ordered_error =
        std::abs(static_cast<long double>(ordered_sum) - expected);
    if (ordered_error > requested_error_budget) {
        return {ValidationFailure::UnstableStoredOrderAccumulation};
    }

    if (has_negative_weight) {
        // The actual stored-order sum is checked above. This independent guard
        // limits sensitivity to signed-weight cancellation without imposing a
        // sample-count ceiling on otherwise well-conditioned rules.
        constexpr long double cancellation_safety_factor = 8.0L;
        const long double epsilon =
            static_cast<long double>(std::numeric_limits<double>::epsilon());
        const long double stability_error_bound =
            compensated_error +
            cancellation_safety_factor * epsilon * absolute_sum;
        const long double stability_budget =
            static_cast<long double>(
                QuadratureRule::default_validation_tolerance()) * scale;
        if (!std::isfinite(stability_error_bound) ||
            stability_error_bound > stability_budget) {
            return {ValidationFailure::IllConditionedWeightCancellation};
        }
    }

    return {};
}

ValidationResult validate_rule_data(
    svmp::CellFamily family,
    int polynomial_exactness,
    const std::vector<QuadPoint>& points,
    const std::vector<double>& weights,
    double tolerance) noexcept
{
    const auto traits = reference_cell_traits(family);
    if (traits.domain == ReferenceDomain::Unsupported) {
        return {ValidationFailure::UnsupportedCellFamily};
    }
    if (polynomial_exactness < 0) {
        return {ValidationFailure::NegativePolynomialExactness};
    }
    if (!std::isfinite(tolerance) || tolerance < 0.0) {
        return {ValidationFailure::InvalidTolerance};
    }
    if (points.empty()) {
        return {ValidationFailure::EmptyRule};
    }
    if (points.size() != weights.size()) {
        return {ValidationFailure::PointWeightSizeMismatch};
    }

    for (std::size_t sample = 0; sample < points.size(); ++sample) {
        const auto result =
            validate_point(points[sample], traits, tolerance, sample);
        if (result.failure != ValidationFailure::None) {
            return result;
        }
    }

    return validate_weights(weights, traits.measure, tolerance);
}

std::string validation_failure_message(const ValidationResult& result)
{
    const auto sample_suffix = [&result]() {
        if (result.sample == std::numeric_limits<std::size_t>::max()) {
            return std::string{};
        }
        return std::string{" at sample "} + std::to_string(result.sample);
    };

    switch (result.failure) {
        case ValidationFailure::None:
            return {};
        case ValidationFailure::UnsupportedCellFamily:
            return "QuadratureRule: unsupported reference-cell family";
        case ValidationFailure::NegativePolynomialExactness:
            return "QuadratureRule: polynomial exactness must be non-negative";
        case ValidationFailure::InvalidTolerance:
            return "QuadratureRule: validation tolerance must be finite and non-negative";
        case ValidationFailure::EmptyRule:
            return "QuadratureRule: a rule must contain at least one sample";
        case ValidationFailure::PointWeightSizeMismatch:
            return "QuadratureRule: points/weights size mismatch";
        case ValidationFailure::NonFiniteCoordinate:
            return "QuadratureRule: quadrature point contains a non-finite coordinate" +
                   sample_suffix();
        case ValidationFailure::NonzeroInactiveCoordinate:
            return "QuadratureRule: quadrature point has a nonzero inactive coordinate" +
                   sample_suffix();
        case ValidationFailure::PointOutsideReferenceCell:
            return "QuadratureRule: quadrature point lies outside the canonical reference cell" +
                   sample_suffix();
        case ValidationFailure::NonFiniteWeight:
            return "QuadratureRule: quadrature weight must be finite" +
                   sample_suffix();
        case ValidationFailure::NonFiniteWeightAccumulation:
            return "QuadratureRule: weight accumulation is not finite" +
                   sample_suffix();
        case ValidationFailure::IncorrectZerothMoment:
            return "QuadratureRule: weights do not reproduce the reference-cell measure";
        case ValidationFailure::UnstableStoredOrderAccumulation:
            return "QuadratureRule: stored-order double-precision weight accumulation "
                   "does not reproduce the reference-cell measure";
        case ValidationFailure::IllConditionedWeightCancellation:
            return "QuadratureRule: weight cancellation is too ill-conditioned for "
                   "stable double-precision integration";
    }
    return "QuadratureRule: unknown validation failure";
}

} // namespace

QuadratureRule::~QuadratureRule() = default;

QuadratureRule::QuadratureRule(svmp::CellFamily family, RuleData data)
    : QuadratureRule(validate(family, std::move(data)))
{
}

QuadratureRule::QuadratureRule(ValidatedState state)
    : cell_family_(state.cell_family),
      dimension_(state.dimension),
      polynomial_exactness_(state.polynomial_exactness),
      reference_measure_(state.reference_measure),
      points_(std::move(state.points)),
      weights_(std::move(state.weights))
{
}

QuadratureRule::ValidatedState QuadratureRule::validate(
    svmp::CellFamily family,
    RuleData data)
{
    const auto validation = validate_rule_data(
        family,
        data.polynomial_exactness,
        data.points,
        data.weights,
        default_validation_tolerance());
    if (validation.failure != ValidationFailure::None) {
        svmp::raise<InvalidArgumentException>(
            validation_failure_message(validation));
    }

    const auto traits = reference_cell_traits(family);
    return {
        family,
        traits.dimension,
        data.polynomial_exactness,
        traits.measure,
        std::move(data.points),
        std::move(data.weights),
    };
}

bool QuadratureRule::is_structurally_valid(double tolerance) const noexcept
{
    return validate_rule_data(
               cell_family_,
               polynomial_exactness_,
               points_,
               weights_,
               tolerance)
               .failure == ValidationFailure::None;
}

} // namespace svmp::FE::quadrature
