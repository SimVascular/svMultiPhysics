// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_TRIANGLE_QUADRATURE_H
#define SVMP_FE_TRIANGLE_QUADRATURE_H

/** @file TriangleQuadrature.h
 * @brief Centroid quadrature on the unit right triangle.
 * @ingroup FE_Quadrature
 */

#include "FE/Quadrature/QuadratureRule.h"

namespace svmp::FE::quadrature {

/** @addtogroup FE_Quadrature
 * @{
 */

/** @brief Return the degree-one centroid rule: point (1/3,1/3,0), weight 1/2.
 * @return Owning CellFamily::Triangle rule for constant or linear integration.
 */
[[nodiscard]] QuadratureRule make_triangle_centroid_rule();

/** @} */
} // namespace svmp::FE::quadrature

#endif // SVMP_FE_TRIANGLE_QUADRATURE_H
