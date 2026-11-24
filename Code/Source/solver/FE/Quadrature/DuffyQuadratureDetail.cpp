// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "FE/Quadrature/DuffyQuadratureDetail.h"

#include "FE/Common/FEException.h"

#include <algorithm>
#include <string>

namespace svmp::FE::quadrature::detail {

double scaled_product(std::initializer_list<double> factors, int exponent)
{
    // Multiply bounded significands and restore the exponent only after the last factor.
    double significand = 1.0;
    for (const double factor : factors) {
        int shift;
        significand *= std::frexp(factor, &shift);
        exponent += shift;
    }
    return std::scalbn(significand, exponent);
}

bool reproduces_measure(std::span<const double> weights, double measure, std::size_t maximum_points)
{
    constexpr double tolerance = 256.0 * std::numeric_limits<double>::epsilon();
    // Power-of-two scaling keeps even the sum of absolute weights in binary64 range.
    const double sum_scale = 0.5 / static_cast<double>(maximum_points);
    double sum = 0.0;
    double correction = 0.0;
    double absolute_sum = 0.0;
    for (const double weight : weights) {
        const double term = weight * sum_scale;
        const double next = sum + term;
        // Neumaier compensation retains small contributions through cancellation.
        correction += std::abs(sum) >= std::abs(term) ? (sum - next) + term : (term - next) + sum;
        sum = next;
        absolute_sum += std::abs(term);
    }
    sum += correction;
    const double scaled_measure = measure * sum_scale;
    return std::isfinite(sum) && std::isfinite(absolute_sum) &&
        std::abs(sum - scaled_measure) <= tolerance * std::max(scaled_measure, absolute_sum);
}

void validate_line_rule(const QuadratureRule& rule, int required_exactness, std::size_t maximum_points,
    const char* context)
{
    const std::string prefix{context};
    svmp::check<InvalidArgumentException>(rule.cell_family() == svmp::CellFamily::Line,
        prefix + " must have canonical Line family");
    svmp::check<InvalidArgumentException>(rule.polynomial_exactness() >= required_exactness,
        prefix + " exactness is below the directional requirement");
    for (const QuadPoint& point : rule.points()) {
        if (std::abs(point[0]) > 1.0) {
            svmp::raise<InvalidArgumentException>(prefix + " contains a point outside [-1,1]");
        }
    }
    svmp::check<InvalidArgumentException>(
        reproduces_measure(rule.weights(), rule.reference_cell_measure(), maximum_points),
        prefix + " normalization does not reproduce measure 2");
}

} // namespace svmp::FE::quadrature::detail
