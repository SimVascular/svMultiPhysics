// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file QuadratureRule.cpp
 * @brief Internal construction and structural validation for quadrature rules.
 * @ingroup FE_Quadrature
 */

#include "FE/Quadrature/QuadratureRule.h"

#include "FE/Common/FEException.h"

#include <cassert>
#include <cmath>
#include <string>
#include <utility>

namespace svmp::FE::quadrature {
namespace {

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

constexpr double reference_measure(svmp::CellFamily family) noexcept
{
    switch (family) {
        case svmp::CellFamily::Point:
            return 1.0;
        case svmp::CellFamily::Line:
            return 2.0;
        case svmp::CellFamily::Triangle:
            return 0.5;
        case svmp::CellFamily::Quad:
            return 4.0;
        case svmp::CellFamily::Tetra:
            return 1.0 / 6.0;
        case svmp::CellFamily::Hex:
            return 8.0;
        case svmp::CellFamily::Wedge:
            return 1.0;
        default:
            return -1.0;
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
            point[component] != 0.0) {
            svmp::raise<InvalidArgumentException>(
                std::string{
                    "QuadratureRule: quadrature point has a nonzero inactive "
                    "coordinate at point index "} +
                std::to_string(point_index));
        }
    }
}

void validate_weights(const std::vector<double>& weights)
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
}

} // namespace

int QuadratureRule::dimension() const noexcept
{
    const int dimension = reference_dimension(cell_family_);
    assert(dimension >= 0);
    return dimension;
}

double QuadratureRule::reference_cell_measure() const noexcept
{
    const double measure = reference_measure(cell_family_);
    assert(measure > 0.0);
    return measure;
}

QuadratureRule::QuadratureRule(
    svmp::CellFamily family,
    int polynomial_exactness,
    std::vector<QuadPoint> points,
    std::vector<double> weights)
    : cell_family_(family),
      polynomial_exactness_(polynomial_exactness),
      points_(std::move(points)),
      weights_(std::move(weights))
{
    const int dimension = reference_dimension(cell_family_);
    svmp::check<InvalidArgumentException>(
        dimension >= 0,
        "QuadratureRule: unsupported reference-cell family");

    svmp::check<InvalidArgumentException>(
        polynomial_exactness_ >= 0,
        "QuadratureRule: polynomial exactness must be non-negative");
    svmp::check<InvalidArgumentException>(
        !points_.empty(),
        "QuadratureRule: a rule must contain at least one point");
    svmp::check<InvalidArgumentException>(
        points_.size() == weights_.size(),
        "QuadratureRule: points/weights size mismatch");

    for (std::size_t point_index = 0;
         point_index < points_.size();
         ++point_index) {
        validate_point(
            points_[point_index],
            dimension,
            point_index);
    }

    validate_weights(weights_);
}

} // namespace svmp::FE::quadrature
