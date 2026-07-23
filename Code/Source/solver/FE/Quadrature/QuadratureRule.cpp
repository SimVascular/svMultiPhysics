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
#include <optional>
#include <string>
#include <string_view>
#include <utility>

namespace svmp::FE::quadrature {
namespace {

struct ReferenceCellTraits {
    int dimension;
    double zeroth_moment;
};

constexpr double coordinate_validation_tolerance = 1.0e-12;
constexpr double moment_validation_tolerance = 1.0e-12;

struct ValidationResult {
    static constexpr std::size_t no_sample =
        std::numeric_limits<std::size_t>::max();

    std::string_view reason{};
    std::size_t sample{no_sample};

    constexpr bool valid() const noexcept
    {
        return reason.empty();
    }
};

constexpr std::optional<ReferenceCellTraits> reference_cell_traits(
    svmp::CellFamily family) noexcept
{
    switch (family) {
        case svmp::CellFamily::Point:
            return ReferenceCellTraits{0, 1.0};
        case svmp::CellFamily::Line:
            return ReferenceCellTraits{1, 2.0};
        case svmp::CellFamily::Triangle:
            return ReferenceCellTraits{2, 0.5};
        case svmp::CellFamily::Quad:
            return ReferenceCellTraits{2, 4.0};
        case svmp::CellFamily::Tetra:
            return ReferenceCellTraits{3, 1.0 / 6.0};
        case svmp::CellFamily::Hex:
            return ReferenceCellTraits{3, 8.0};
        case svmp::CellFamily::Wedge:
            return ReferenceCellTraits{3, 1.0};
        default:
            return std::nullopt;
    }
}

ValidationResult validate_point(
    const QuadPoint& point,
    svmp::CellFamily family,
    const ReferenceCellTraits& traits,
    std::size_t sample) noexcept
{
    for (std::size_t component = 0; component < 3u; ++component) {
        if (!std::isfinite(point[component])) {
            return {"quadrature point contains a non-finite coordinate", sample};
        }
        if (component >=
            static_cast<std::size_t>(traits.dimension) &&
            std::abs(point[component]) > coordinate_validation_tolerance) {
            return {"quadrature point has a nonzero inactive coordinate", sample};
        }
    }

    const auto in_interval = [](double value, double lower, double upper) {
        return value >= lower - coordinate_validation_tolerance &&
               value <= upper + coordinate_validation_tolerance;
    };

    const double x = point[0];
    const double y = point[1];
    const double z = point[2];

    bool is_contained = false;
    switch (family) {
        case svmp::CellFamily::Point:
            is_contained = true;
            break;
        case svmp::CellFamily::Line:
            is_contained = in_interval(x, -1.0, 1.0);
            break;
        case svmp::CellFamily::Triangle:
            is_contained =
                x >= -coordinate_validation_tolerance &&
                y >= -coordinate_validation_tolerance &&
                x + y <= 1.0 + coordinate_validation_tolerance;
            break;
        case svmp::CellFamily::Quad:
            is_contained = in_interval(x, -1.0, 1.0) &&
                           in_interval(y, -1.0, 1.0);
            break;
        case svmp::CellFamily::Tetra:
            is_contained =
                x >= -coordinate_validation_tolerance &&
                y >= -coordinate_validation_tolerance &&
                z >= -coordinate_validation_tolerance &&
                x + y + z <= 1.0 + coordinate_validation_tolerance;
            break;
        case svmp::CellFamily::Hex:
            is_contained = in_interval(x, -1.0, 1.0) &&
                           in_interval(y, -1.0, 1.0) &&
                           in_interval(z, -1.0, 1.0);
            break;
        case svmp::CellFamily::Wedge:
            is_contained =
                x >= -coordinate_validation_tolerance &&
                y >= -coordinate_validation_tolerance &&
                x + y <= 1.0 + coordinate_validation_tolerance &&
                in_interval(z, -1.0, 1.0);
            break;
        default:
            return {"unsupported reference-cell family"};
    }

    if (!is_contained) {
        return {
            "quadrature point lies outside the canonical reference cell",
            sample};
    }
    return {};
}

ValidationResult validate_weights(
    const std::vector<double>& weights,
    double zeroth_moment) noexcept
{
    long double compensated_sum = 0.0L;
    long double correction = 0.0L;
    long double absolute_sum = 0.0L;

    for (std::size_t sample = 0; sample < weights.size(); ++sample) {
        const double weight = weights[sample];
        if (!std::isfinite(weight)) {
            return {"quadrature weight must be finite", sample};
        }

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
            return {"weight accumulation is not finite", sample};
        }
    }

    const long double total = compensated_sum + correction;
    if (!std::isfinite(total)) {
        return {"weight accumulation is not finite", weights.size() - 1u};
    }

    const long double expected = static_cast<long double>(zeroth_moment);
    const long double scale = std::max(1.0L, std::abs(expected));
    const long double requested_error_budget =
        static_cast<long double>(moment_validation_tolerance) * scale;
    const long double compensated_error = std::abs(total - expected);

    if (compensated_error > requested_error_budget) {
        return {"weights do not reproduce the zeroth moment"};
    }

    // A single compensated accumulator is not exact for unrestricted signed
    // data. Use a conservative uncertainty estimate in the precision used by
    // validation before accepting the computed moment.
    constexpr long double validation_safety_factor = 8.0L;
    const long double validation_error_bound =
        compensated_error +
        validation_safety_factor *
            std::numeric_limits<long double>::epsilon() * absolute_sum;
    if (!std::isfinite(validation_error_bound) ||
        validation_error_bound > requested_error_budget) {
        return {
            "weight sum is too ill-conditioned to validate the zeroth moment"};
    }

    return {};
}

ValidationResult validate_rule_data(
    svmp::CellFamily family,
    const ReferenceCellTraits& traits,
    int polynomial_exactness,
    const std::vector<QuadPoint>& points,
    const std::vector<double>& weights) noexcept
{
    if (polynomial_exactness < 0) {
        return {"polynomial exactness must be non-negative"};
    }
    if (points.empty()) {
        return {"a rule must contain at least one sample"};
    }
    if (points.size() != weights.size()) {
        return {"points/weights size mismatch"};
    }

    for (std::size_t sample = 0; sample < points.size(); ++sample) {
        const auto result =
            validate_point(points[sample], family, traits, sample);
        if (!result.valid()) {
            return result;
        }
    }

    return validate_weights(weights, traits.zeroth_moment);
}

std::string validation_failure_message(const ValidationResult& result)
{
    if (result.valid()) {
        return {};
    }

    std::string message{"QuadratureRule: "};
    message.append(result.reason);
    if (result.sample != ValidationResult::no_sample) {
        message.append(" at sample ");
        message.append(std::to_string(result.sample));
    }
    return message;
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
      zeroth_moment_(state.zeroth_moment),
      points_(std::move(state.points)),
      weights_(std::move(state.weights))
{
}

QuadratureRule::ValidatedState QuadratureRule::validate(
    svmp::CellFamily family,
    RuleData data)
{
    const auto traits = reference_cell_traits(family);
    if (!traits) {
        svmp::raise<InvalidArgumentException>(
            validation_failure_message(
                {"unsupported reference-cell family"}));
    }

    const auto validation = validate_rule_data(
        family,
        *traits,
        data.polynomial_exactness,
        data.points,
        data.weights);
    if (!validation.valid()) {
        svmp::raise<InvalidArgumentException>(
            validation_failure_message(validation));
    }

    return {
        family,
        traits->dimension,
        data.polynomial_exactness,
        traits->zeroth_moment,
        std::move(data.points),
        std::move(data.weights),
    };
}

} // namespace svmp::FE::quadrature
