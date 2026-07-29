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
#include <cassert>
#include <cmath>
#include <numeric>
#include <string>
#include <utility>

namespace svmp::FE::quadrature {
namespace {

constexpr double inactive_coordinate_tolerance = 1.0e-12;
constexpr double measure_validation_tolerance = 1.0e-12;

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
            std::abs(point[component]) > inactive_coordinate_tolerance) {
            svmp::raise<InvalidArgumentException>(
                std::string{
                    "QuadratureRule: quadrature point has a nonzero inactive "
                    "coordinate at point index "} +
                std::to_string(point_index));
        }
    }
}

void validate_weights(
    const std::vector<double>& weights,
    double reference_cell_measure)
{
    for (std::size_t point_index = 0;
         point_index < weights.size();
         ++point_index) {
        if (!std::isfinite(weights[point_index])) {
            svmp::raise<InvalidArgumentException>(
                std::string{
                    "QuadratureRule: quadrature weight must be finite at point "
                    "index "} +
                std::to_string(point_index));
        }
    }

    const double weight_sum =
        std::accumulate(weights.begin(), weights.end(), 0.0);
    const double scale = std::max(1.0, std::abs(reference_cell_measure));
    const double error_budget = measure_validation_tolerance * scale;
    svmp::check<InvalidArgumentException>(
        std::abs(weight_sum - reference_cell_measure) <= error_budget,
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
            dimension,
            point_index);
    }

    validate_weights(weights_, reference_cell_measure_);
}

} // namespace svmp::FE::quadrature
