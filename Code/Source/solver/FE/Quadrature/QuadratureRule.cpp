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
#include <array>
#include <bit>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>

namespace svmp::FE::quadrature {
namespace {

constexpr double coordinate_validation_tolerance = 1.0e-12;
constexpr double measure_validation_tolerance = 1.0e-12;

// Exact summation decodes double bit patterns directly, so fail at compile
// time unless both the binary64 value model and object layout are supported.
static_assert(
    std::numeric_limits<double>::is_iec559 &&
        std::numeric_limits<double>::radix == 2 &&
        std::numeric_limits<double>::digits == 53 &&
        std::numeric_limits<double>::min_exponent == -1021 &&
        std::numeric_limits<double>::max_exponent == 1024 &&
        sizeof(double) == sizeof(std::uint64_t),
    "Quadrature validation requires IEEE 754 binary64 doubles");
static_assert(
    std::bit_cast<std::uint64_t>(1.0) == 0x3ff0000000000000ULL &&
        std::bit_cast<std::uint64_t>(-0.0) == 0x8000000000000000ULL &&
        std::bit_cast<std::uint64_t>(
            std::numeric_limits<double>::denorm_min()) == 1u,
    "Quadrature validation requires the standard binary64 bit layout");

/**
 * @brief Exact, order-independent sum of finite binary64 values.
 *
 * Every finite binary64 value is an integer multiple of 2^-1074. Positive and
 * negative coefficients are accumulated separately in fixed-width unsigned
 * integers, so neither cancellation nor an intermediate floating-point
 * overflow can affect the result.
 */
class ExactBinary64Sum {
public:
    /**
     * @brief Add one value exactly.
     * @return False when @p value is infinite or NaN.
     */
    bool add(double value) noexcept
    {
        const std::uint64_t bits = std::bit_cast<std::uint64_t>(value);
        const std::uint64_t exponent =
            (bits >> fraction_bit_count) & exponent_mask;
        if (exponent == exponent_mask) {
            return false;
        }

        const std::uint64_t fraction = bits & fraction_mask;
        const std::uint64_t coefficient =
            exponent == 0u ? fraction : hidden_bit | fraction;
        const std::size_t shift =
            exponent == 0u ? 0u : static_cast<std::size_t>(exponent - 1u);

        auto& magnitude =
            (bits & sign_mask) == 0u ? positive_ : negative_;
        add_shifted(magnitude, coefficient, shift);
        return true;
    }

    /**
     * @brief Test the exact sum against a binary64 target and error budget.
     */
    bool is_within(double expected, double error_budget) const noexcept
    {
        if (!(error_budget >= 0.0)) {
            return false;
        }

        ExactBinary64Sum residual = *this;
        if (!residual.add(-expected)) {
            return false;
        }

        ExactBinary64Sum budget;
        if (!budget.add(error_budget)) {
            return false;
        }

        return compare(
                   residual.absolute_magnitude(),
                   budget.absolute_magnitude()) <= 0;
    }

private:
    static constexpr std::size_t limb_bit_count = std::numeric_limits<std::uint64_t>::digits;
    static constexpr std::size_t fraction_bit_count = 52u;
    static constexpr std::uint64_t hidden_bit = std::uint64_t{1} << fraction_bit_count;
    static constexpr std::uint64_t fraction_mask = hidden_bit - 1u;
    static constexpr std::uint64_t exponent_mask = 0x7ffu;
    static constexpr std::uint64_t sign_mask = std::uint64_t{1} << (limb_bit_count - 1u);

    // A finite binary64 coefficient occupies at most this many bits when
    // scaled by 2^1074. The size_t term covers every possible vector length
    // plus the expected moment without dynamic allocation.
    static constexpr std::size_t maximum_scaled_value_bits =
        static_cast<std::size_t>(
            std::numeric_limits<double>::max_exponent -
            std::numeric_limits<double>::min_exponent +
            std::numeric_limits<double>::digits);
    static constexpr std::size_t accumulator_bit_count =
        maximum_scaled_value_bits +
        std::numeric_limits<std::size_t>::digits;
    static constexpr std::size_t accumulator_limb_count =
        (accumulator_bit_count + limb_bit_count - 1u) / limb_bit_count;

    using Magnitude =
        std::array<std::uint64_t, accumulator_limb_count>;

    static void add_word(
        Magnitude& magnitude,
        std::size_t index,
        std::uint64_t word) noexcept
    {
        while (word != 0u && index < magnitude.size()) {
            const std::uint64_t previous = magnitude[index];
            magnitude[index] += word;
            word = magnitude[index] < previous ? 1u : 0u;
            ++index;
        }
        assert(word == 0u);
    }

    static void add_shifted(
        Magnitude& magnitude,
        std::uint64_t coefficient,
        std::size_t shift) noexcept
    {
        if (coefficient == 0u) {
            return;
        }

        const std::size_t index = shift / limb_bit_count;
        const std::size_t offset = shift % limb_bit_count;
        add_word(magnitude, index, coefficient << offset);
        if (offset != 0u) {
            add_word(
                magnitude,
                index + 1u,
                coefficient >> (limb_bit_count - offset));
        }
    }

    static int compare(
        const Magnitude& left,
        const Magnitude& right) noexcept
    {
        for (std::size_t index = left.size(); index-- > 0u;) {
            if (left[index] < right[index]) {
                return -1;
            }
            if (left[index] > right[index]) {
                return 1;
            }
        }
        return 0;
    }

    static Magnitude subtract(
        const Magnitude& larger,
        const Magnitude& smaller) noexcept
    {
        Magnitude difference{};
        std::uint64_t borrow = 0u;
        for (std::size_t index = 0u; index < difference.size(); ++index) {
            const std::uint64_t after_subtraction =
                larger[index] - smaller[index];
            const bool subtraction_borrow =
                larger[index] < smaller[index];
            difference[index] = after_subtraction - borrow;
            const bool incoming_borrow = after_subtraction < borrow;
            borrow =
                subtraction_borrow || incoming_borrow ? 1u : 0u;
        }
        assert(borrow == 0u);
        return difference;
    }

    Magnitude absolute_magnitude() const noexcept
    {
        return compare(positive_, negative_) >= 0
                   ? subtract(positive_, negative_)
                   : subtract(negative_, positive_);
    }

    Magnitude positive_{};
    Magnitude negative_{};
};

constexpr int reference_dimension(svmp::CellFamily family) noexcept
{
    switch (family) {
        case svmp::CellFamily::Point:
            return 0;
        case svmp::CellFamily::Line:
            return 1;
        case svmp::CellFamily::Triangle:
        case svmp::CellFamily::Quad:
            return 2;
        case svmp::CellFamily::Tetra:
        case svmp::CellFamily::Hex:
        case svmp::CellFamily::Wedge:
            return 3;
        default:
            return -1;
    }
}

void validate_point(
    const QuadPoint& point,
    svmp::CellFamily family,
    int dimension,
    std::size_t point_index)
{
    for (std::size_t component = 0; component < 3u; ++component) {
        if (!std::isfinite(point[component])) {
            svmp::raise<InvalidArgumentException>(
                std::string{
                    "QuadratureRule: quadrature point contains a non-finite "
                    "coordinate at point index "} +
                std::to_string(point_index));
        }
        if (component >=
            static_cast<std::size_t>(dimension) &&
            std::abs(point[component]) > coordinate_validation_tolerance) {
            svmp::raise<InvalidArgumentException>(
                std::string{
                    "QuadratureRule: quadrature point has a nonzero inactive "
                    "coordinate at point index "} +
                std::to_string(point_index));
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
            svmp::raise<InvalidArgumentException>(
                "QuadratureRule: unsupported reference-cell family");
    }

    if (!is_contained) {
        svmp::raise<InvalidArgumentException>(
            std::string{
                "QuadratureRule: quadrature point lies outside the canonical "
                "reference cell at point index "} +
            std::to_string(point_index));
    }
}

void validate_weights(
    const std::vector<double>& weights,
    double reference_cell_measure)
{
    ExactBinary64Sum exact_sum;

    for (std::size_t point_index = 0;
         point_index < weights.size();
         ++point_index) {
        if (!exact_sum.add(weights[point_index])) {
            svmp::raise<InvalidArgumentException>(
                std::string{
                    "QuadratureRule: quadrature weight must be finite at point "
                    "index "} +
                std::to_string(point_index));
        }
    }

    const double scale = std::max(1.0, std::abs(reference_cell_measure));
    const double error_budget = measure_validation_tolerance * scale;
    svmp::check<InvalidArgumentException>(
        exact_sum.is_within(reference_cell_measure, error_budget),
        "QuadratureRule: weights do not reproduce the reference-cell measure");
}

} // namespace

int QuadratureRule::dimension() const noexcept
{
    const int dimension = reference_dimension(cell_family_);
    assert(dimension >= 0);
    return dimension;
}

QuadratureRule::QuadratureRule(
    svmp::CellFamily family,
    int polynomial_exactness,
    std::vector<QuadPoint> points,
    std::vector<double> weights)
    : cell_family_(family),
      polynomial_exactness_(polynomial_exactness),
      reference_cell_measure_(0.0),
      points_(std::move(points)),
      weights_(std::move(weights))
{
    switch (cell_family_) {
        case svmp::CellFamily::Point:
            reference_cell_measure_ = 1.0;
            break;
        case svmp::CellFamily::Line:
            reference_cell_measure_ = 2.0;
            break;
        case svmp::CellFamily::Triangle:
            reference_cell_measure_ = 0.5;
            break;
        case svmp::CellFamily::Quad:
            reference_cell_measure_ = 4.0;
            break;
        case svmp::CellFamily::Tetra:
            reference_cell_measure_ = 1.0 / 6.0;
            break;
        case svmp::CellFamily::Hex:
            reference_cell_measure_ = 8.0;
            break;
        case svmp::CellFamily::Wedge:
            reference_cell_measure_ = 1.0;
            break;
        default:
            svmp::raise<InvalidArgumentException>(
                "QuadratureRule: unsupported reference-cell family");
    }

    svmp::check<InvalidArgumentException>(
        polynomial_exactness_ >= 0,
        "QuadratureRule: polynomial exactness must be non-negative");
    svmp::check<InvalidArgumentException>(
        !points_.empty(),
        "QuadratureRule: a rule must contain at least one point");
    svmp::check<InvalidArgumentException>(
        points_.size() == weights_.size(),
        "QuadratureRule: points/weights size mismatch");

    const int dimension = reference_dimension(cell_family_);
    assert(dimension >= 0);

    for (std::size_t point_index = 0;
         point_index < points_.size();
         ++point_index) {
        validate_point(
            points_[point_index],
            cell_family_,
            dimension,
            point_index);
    }

    validate_weights(weights_, reference_cell_measure_);
}

} // namespace svmp::FE::quadrature
