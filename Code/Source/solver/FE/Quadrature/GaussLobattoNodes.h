// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_GAUSS_LOBATTO_NODES_H
#define SVMP_FE_GAUSS_LOBATTO_NODES_H

#include <cstddef>

namespace svmp::FE::quadrature::detail {

/** @cond INTERNAL */

// Shared node refinement for interpolation and integration. The Legendre
// value is retained so Quadrature can form a weight without another recurrence.
struct GaussLobattoNode {
    double coordinate;
    double legendre_value; // P_(num_points-1) at coordinate.
    int iterations;
};

// Refine a nonnegative interior node of an n-point GLL distribution on [-1, 1].
// half_root_index runs from the positive endpoint toward the center, starting
// at zero; an odd distribution's center is returned as exactly zero.
// Requires n >= 3, index < (n-1)/2, and a finite positive correction tolerance.
// Invalid requests throw InvalidArgumentException; unsuccessful refinement
// throws ConvergenceException. Callers own endpoints, mirroring, and size caps.
[[nodiscard]] GaussLobattoNode refine_gauss_lobatto_node(
    std::size_t num_points, std::size_t half_root_index,
    double correction_tolerance);

/** @endcond */

} // namespace svmp::FE::quadrature::detail

#endif // SVMP_FE_GAUSS_LOBATTO_NODES_H
