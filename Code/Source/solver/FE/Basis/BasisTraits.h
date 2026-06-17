// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_BASIS_BASISTRAITS_H
#define SVMP_FE_BASIS_BASISTRAITS_H

#include "Types.h"

#include <cstddef>
#include <limits>

namespace svmp {
namespace FE {
namespace basis {

enum class BasisTopology {
    Unknown,
    Point,
    Line,
    Triangle,
    Quadrilateral,
    Tetrahedron,
    Hexahedron,
    Wedge,
};

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

[[nodiscard]] constexpr bool is_point(ElementType type) noexcept {
    return type == ElementType::Point1;
}

[[nodiscard]] constexpr bool is_line(ElementType type) noexcept {
    return type == ElementType::Line2 || type == ElementType::Line3;
}

[[nodiscard]] constexpr bool is_triangle(ElementType type) noexcept {
    return type == ElementType::Triangle3 || type == ElementType::Triangle6;
}

[[nodiscard]] constexpr bool is_quadrilateral(ElementType type) noexcept {
    return type == ElementType::Quad4 || type == ElementType::Quad8 ||
           type == ElementType::Quad9;
}

[[nodiscard]] constexpr bool is_tetrahedron(ElementType type) noexcept {
    return type == ElementType::Tetra4 || type == ElementType::Tetra10;
}

[[nodiscard]] constexpr bool is_hexahedron(ElementType type) noexcept {
    return type == ElementType::Hex8 || type == ElementType::Hex20 ||
           type == ElementType::Hex27;
}

[[nodiscard]] constexpr bool is_wedge(ElementType type) noexcept {
    return type == ElementType::Wedge6 || type == ElementType::Wedge15 ||
           type == ElementType::Wedge18;
}

[[nodiscard]] constexpr bool is_pyramid(ElementType type) noexcept {
    (void)type;
    return false;
}

[[nodiscard]] constexpr bool is_simplex(ElementType type) noexcept {
    return is_triangle(type) || is_tetrahedron(type);
}

[[nodiscard]] constexpr bool is_tensor_product(ElementType type) noexcept {
    return is_line(type) || is_quadrilateral(type) || is_hexahedron(type);
}

[[nodiscard]] constexpr int reference_dimension(ElementType type) noexcept {
    return element_dimension(type);
}

[[nodiscard]] constexpr BasisTopology topology(ElementType type) noexcept {
    if (is_point(type)) {
        return BasisTopology::Point;
    }
    if (is_line(type)) {
        return BasisTopology::Line;
    }
    if (is_triangle(type)) {
        return BasisTopology::Triangle;
    }
    if (is_quadrilateral(type)) {
        return BasisTopology::Quadrilateral;
    }
    if (is_tetrahedron(type)) {
        return BasisTopology::Tetrahedron;
    }
    if (is_hexahedron(type)) {
        return BasisTopology::Hexahedron;
    }
    if (is_wedge(type)) {
        return BasisTopology::Wedge;
    }
    return BasisTopology::Unknown;
}

[[nodiscard]] constexpr ElementType canonical_lagrange_type(ElementType type) noexcept {
    switch (type) {
        case ElementType::Line2:
        case ElementType::Line3:
            return ElementType::Line2;
        case ElementType::Triangle3:
        case ElementType::Triangle6:
            return ElementType::Triangle3;
        case ElementType::Quad4:
        case ElementType::Quad9:
            return ElementType::Quad4;
        case ElementType::Tetra4:
        case ElementType::Tetra10:
            return ElementType::Tetra4;
        case ElementType::Hex8:
        case ElementType::Hex27:
            return ElementType::Hex8;
        case ElementType::Wedge6:
        case ElementType::Wedge18:
            return ElementType::Wedge6;
        default:
            return type;
    }
}

[[nodiscard]] constexpr int complete_lagrange_alias_order(ElementType type) noexcept {
    switch (type) {
        case ElementType::Line2:
        case ElementType::Triangle3:
        case ElementType::Quad4:
        case ElementType::Tetra4:
        case ElementType::Hex8:
        case ElementType::Wedge6:
            return 1;
        case ElementType::Line3:
        case ElementType::Triangle6:
        case ElementType::Quad9:
        case ElementType::Tetra10:
        case ElementType::Hex27:
        case ElementType::Wedge18:
            return 2;
        default:
            return -1;
    }
}

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_BASISTRAITS_H
