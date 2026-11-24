// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file TetrahedronQuadrature.cpp
 * @brief Validated composition of supplied line quadrature into tetrahedron rules.
 * @ingroup FE_Quadrature
 */

#include "FE/Quadrature/TetrahedronQuadrature.h"

#include "FE/Common/FEException.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace svmp::FE::quadrature {
namespace {

constexpr std::size_t kMaximumProductPoints = 128u * 128u * 128u;
constexpr double kValidationTolerance = 256.0 * std::numeric_limits<double>::epsilon();
// Scale validation sums by a power of two to keep bounded sums in binary64 range.
// Stored weights and relative tolerances are unchanged.
constexpr double kValidationSumScale = 0.5 / static_cast<double>(kMaximumProductPoints);

// A later large weight can amplify rounding in this subnormal or zero intermediate.
bool needs_scaled_product(double product, double left, double right)
{
    return std::abs(product) < std::numeric_limits<double>::min() && left != 0.0 && right != 0.0;
}

// Multiply bounded significands and restore the exponent only after the last factor.
double scaled_product(std::initializer_list<double> factors, int exponent)
{
    double significand = 1.0;
    for (const double factor : factors) {
        int shift;
        significand *= std::frexp(factor, &shift);
        exponent += shift;
    }
    return std::scalbn(significand, exponent);
}

// Return whether abs(sum(w)-measure) <= tau*max(measure,sum(abs(w))) for at most
// kMaximumProductPoints weights. Any nonfinite weight fails.
bool reproduces_measure(std::span<const double> weights, double measure)
{
    double sum = 0.0;
    double correction = 0.0;
    double scale = 0.0;
    for (const double weight : weights) {
        const double term = weight * kValidationSumScale;
        // Neumaier compensation keeps long signed sums accurate.
        const double next = sum + term;
        correction += std::abs(sum) >= std::abs(term) ? (sum - next) + term : (term - next) + sum;
        sum = next;
        scale += std::abs(term);
    }
    sum += correction;
    const double scaled_measure = measure * kValidationSumScale;
    return std::isfinite(sum) && std::isfinite(scale) &&
        std::abs(sum - scaled_measure) <= kValidationTolerance * std::max(scaled_measure, scale);
}

void validate_line_rule(const QuadratureRule& rule, int required_exactness, const char* name)
{
    const std::string prefix = std::string{"Tetrahedron Duffy generator: "} + name;
    svmp::check<InvalidArgumentException>(rule.cell_family() == svmp::CellFamily::Line,
        prefix + " must have canonical Line family");
    svmp::check<InvalidArgumentException>(rule.polynomial_exactness() >= required_exactness,
        prefix + " exactness is below the directional requirement");
    for (const QuadPoint& point : rule.points()) {
        if (std::abs(point[0]) > 1.0) {
            svmp::raise<InvalidArgumentException>(prefix + " contains a point outside [-1,1]");
        }
    }
    svmp::check<InvalidArgumentException>(reproduces_measure(rule.weights(), 2.0),
        prefix + " normalization does not reproduce measure 2");
}

} // namespace

QuadratureRule make_tetrahedron_centroid_rule()
{
    return QuadratureRule(svmp::CellFamily::Tetra, 1, {{0.25, 0.25, 0.25}}, {1.0 / 6.0});
}

QuadratureRule make_tetrahedron_rule(
    int requested_exactness, const QuadratureRule& u_rule,
    const QuadratureRule& v_rule, const QuadratureRule& w_rule)
{
    svmp::check<InvalidArgumentException>(
        requested_exactness >= 0 && requested_exactness <= max_tetrahedron_exactness(),
        "Tetrahedron Duffy generator: requested_exactness must be in [0, " +
            std::to_string(max_tetrahedron_exactness()) + ']');
    const std::size_t n_u = u_rule.num_points();
    const std::size_t n_v = v_rule.num_points();
    const std::size_t n_w = w_rule.num_points();
    svmp::check<InvalidArgumentException>(n_u != 0 && n_v != 0 && n_w != 0,
        "Tetrahedron Duffy generator: empty input rule");
    // Division bounds the full product and each validation sum before multiplication.
    svmp::check<InvalidArgumentException>(
        n_u <= kMaximumProductPoints / n_v / n_w,
        "Tetrahedron Duffy generator: product exceeds " +
            std::to_string(kMaximumProductPoints) + " samples");
    validate_line_rule(u_rule, requested_exactness + 2, "u_rule");
    validate_line_rule(v_rule, requested_exactness + 1, "v_rule");
    validate_line_rule(w_rule, requested_exactness, "w_rule");

    // Checked interval nodes and the Duffy map establish containment up to rounding.
    std::vector<QuadPoint> points(n_u * n_v * n_w, QuadPoint::Zero());
    std::vector<double> weights(points.size());
    for (std::size_t i = 0; i < n_u; ++i) {
        const double xi = u_rule.point(i)[0];
        const double u = (1.0 + xi) / 2.0;
        // An interior node can round to u == 1 while its distance to the face is still positive.
        const double remaining_u = (1.0 - xi) / 2.0;
        const double u_weight = u_rule.weight(i);
        const double u_scale = remaining_u * remaining_u / 8.0;
        const double outer_weight = u_weight * u_scale;
        const bool scale_outer = needs_scaled_product(outer_weight, u_weight, u_scale);
        for (std::size_t j = 0; j < n_v; ++j) {
            const double eta = v_rule.point(j)[0];
            const double v = (1.0 + eta) / 2.0;
            const double remaining_v = (1.0 - eta) / 2.0;
            const double v_weight = v_rule.weight(j);
            const double weighted_v = v_weight * remaining_v;
            const double middle_weight = outer_weight * weighted_v;
            const bool rescale = scale_outer || needs_scaled_product(weighted_v, v_weight, remaining_v) ||
                needs_scaled_product(middle_weight, outer_weight, weighted_v);
            for (std::size_t k = 0; k < n_w; ++k) {
                const std::size_t q = (i * n_v + j) * n_w + k;
                const double w = (w_rule.point(k)[0] + 1.0) / 2.0;
                points[q] = QuadPoint{u, remaining_u * v, remaining_u * remaining_v * w};
                const double w_weight = w_rule.weight(k);
                weights[q] = rescale
                    ? scaled_product({u_weight, v_weight, w_weight, remaining_u, remaining_u, remaining_v}, -3)
                    : middle_weight * w_weight;
            }
        }
    }
    // Validate transformed volume and reject nonfinite products. This does not
    // certify individual factor moments; directional errors can compensate.
    svmp::check<ConvergenceException>(reproduces_measure(weights, 1.0 / 6.0),
        "Tetrahedron Duffy generator: output weights do not reproduce measure 1/6");

    const int polynomial_exactness = std::min({
        u_rule.polynomial_exactness() - 2, v_rule.polynomial_exactness() - 1,
        w_rule.polynomial_exactness()});
    return QuadratureRule(
        svmp::CellFamily::Tetra, polynomial_exactness, std::move(points), std::move(weights));
}

QuadratureRule make_tetrahedron_rule(int requested_exactness, const QuadratureRule& line_rule)
{
    return make_tetrahedron_rule(requested_exactness, line_rule, line_rule, line_rule);
}

} // namespace svmp::FE::quadrature
