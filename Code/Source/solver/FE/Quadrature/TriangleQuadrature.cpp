// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/** @file TriangleQuadrature.cpp
 * @brief Centroid quadrature on the unit right triangle.
 * @ingroup FE_Quadrature
 */

#include "FE/Quadrature/TriangleQuadrature.h"

namespace svmp::FE::quadrature {

QuadratureRule make_triangle_centroid_rule()
{
    return QuadratureRule(svmp::CellFamily::Triangle, 1, {{1.0 / 3.0, 1.0 / 3.0, 0.0}}, {0.5});
}

} // namespace svmp::FE::quadrature
