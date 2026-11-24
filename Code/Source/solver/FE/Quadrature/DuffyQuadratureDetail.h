// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_DUFFY_QUADRATURE_DETAIL_H
#define SVMP_FE_DUFFY_QUADRATURE_DETAIL_H

#include "FE/Quadrature/QuadratureRule.h"

#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <span>

namespace svmp::FE::quadrature::detail {

/** @cond INTERNAL */

// An interior node can round to 1 after mapping to [0,1], while this distance stays positive.
inline double unit_interval_complement(double line_node)
{
    return (1.0 - line_node) / 2.0;
}

// A later large weight can amplify rounding in a subnormal or zero intermediate.
inline bool needs_scaled_product(double product, double left, double right)
{
    return std::abs(product) < std::numeric_limits<double>::min() && left != 0.0 && right != 0.0;
}

// Return 2^exponent times the product of the finite weight and Jacobian factors.
double scaled_product(std::initializer_list<double> factors, int exponent);

// maximum_points is the provider's power-of-two cap, checked before validating any rule.
bool reproduces_measure(std::span<const double> weights, double measure, std::size_t maximum_points);
void validate_line_rule(const QuadratureRule& rule, int required_exactness, std::size_t maximum_points,
    const char* context);

/** @endcond */

} // namespace svmp::FE::quadrature::detail

#endif // SVMP_FE_DUFFY_QUADRATURE_DETAIL_H
