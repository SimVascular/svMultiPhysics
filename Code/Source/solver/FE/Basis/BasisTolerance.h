/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#ifndef SVMP_FE_BASIS_BASISTOLERANCE_H
#define SVMP_FE_BASIS_BASISTOLERANCE_H

#include "Types.h"

#include <limits>

namespace svmp {
namespace FE {
namespace basis {
namespace detail {

[[nodiscard]] constexpr Real basis_abs(Real value) noexcept {
    return value < Real(0) ? -value : value;
}

[[nodiscard]] constexpr Real basis_max(Real lhs, Real rhs) noexcept {
    return lhs < rhs ? rhs : lhs;
}

[[nodiscard]] constexpr Real basis_scaled_tolerance(Real scale = Real(1),
                                                    Real multiplier = Real(64)) noexcept {
    return multiplier * std::numeric_limits<Real>::epsilon() *
           basis_max(Real(1), basis_abs(scale));
}

[[nodiscard]] constexpr bool basis_near_zero(Real value,
                                             Real scale = Real(1),
                                             Real multiplier = Real(64)) noexcept {
    return basis_abs(value) <= basis_scaled_tolerance(scale, multiplier);
}

[[nodiscard]] constexpr bool basis_nearly_equal(Real a,
                                                Real b,
                                                Real multiplier = Real(64)) noexcept {
    const Real scale = basis_max(Real(1), basis_max(basis_abs(a), basis_abs(b)));
    return basis_abs(a - b) <= basis_scaled_tolerance(scale, multiplier);
}

} // namespace detail
} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_BASISTOLERANCE_H
