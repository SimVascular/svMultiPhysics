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

[[nodiscard]] constexpr double basis_abs(double value) noexcept {
    return value < double(0) ? -value : value;
}

[[nodiscard]] constexpr double basis_max(double lhs, double rhs) noexcept {
    return lhs < rhs ? rhs : lhs;
}

[[nodiscard]] constexpr double basis_scaled_tolerance(double scale = double(1),
                                                    double multiplier = double(64)) noexcept {
    return multiplier * std::numeric_limits<double>::epsilon() *
           basis_max(double(1), basis_abs(scale));
}

[[nodiscard]] constexpr bool basis_near_zero(double value,
                                             double scale = double(1),
                                             double multiplier = double(64)) noexcept {
    return basis_abs(value) <= basis_scaled_tolerance(scale, multiplier);
}

[[nodiscard]] constexpr bool basis_nearly_equal(double a,
                                                double b,
                                                double multiplier = double(64)) noexcept {
    const double scale = basis_max(double(1), basis_max(basis_abs(a), basis_abs(b)));
    return basis_abs(a - b) <= basis_scaled_tolerance(scale, multiplier);
}

} // namespace detail

// Reference-cell topology is derived from the single mesh cell-family
// classification (to_mesh_family) so the basis layer never maintains a parallel
// ElementType->shape switch; adding an ElementType updates only to_mesh_family.
// ElementType::Unknown must stay Unknown here: CellFamily has no "unknown"
// member, so to_mesh_family() falls back to Point for unrecognized types.
[[nodiscard]] constexpr BasisTopology topology(ElementType type) noexcept {
    if (type == ElementType::Unknown) {
        return BasisTopology::Unknown;
    }
    switch (to_mesh_family(type)) {
        case CellFamily::Point:    return BasisTopology::Point;
        case CellFamily::Line:     return BasisTopology::Line;
        case CellFamily::Triangle: return BasisTopology::Triangle;
        case CellFamily::Quad:     return BasisTopology::Quadrilateral;
        case CellFamily::Tetra:    return BasisTopology::Tetrahedron;
        case CellFamily::Hex:      return BasisTopology::Hexahedron;
        case CellFamily::Wedge:    return BasisTopology::Wedge;
        // Pyramid/Polygon/Polyhedron are outside the current basis scope.
        default:                   return BasisTopology::Unknown;
    }
}

// The shape predicates derive from topology() so they share its single source.
[[nodiscard]] constexpr bool is_point(ElementType type) noexcept {
    return topology(type) == BasisTopology::Point;
}

[[nodiscard]] constexpr bool is_line(ElementType type) noexcept {
    return topology(type) == BasisTopology::Line;
}

[[nodiscard]] constexpr bool is_triangle(ElementType type) noexcept {
    return topology(type) == BasisTopology::Triangle;
}

[[nodiscard]] constexpr bool is_quadrilateral(ElementType type) noexcept {
    return topology(type) == BasisTopology::Quadrilateral;
}

[[nodiscard]] constexpr bool is_tetrahedron(ElementType type) noexcept {
    return topology(type) == BasisTopology::Tetrahedron;
}

[[nodiscard]] constexpr bool is_hexahedron(ElementType type) noexcept {
    return topology(type) == BasisTopology::Hexahedron;
}

[[nodiscard]] constexpr bool is_wedge(ElementType type) noexcept {
    return topology(type) == BasisTopology::Wedge;
}

[[nodiscard]] constexpr int reference_dimension(ElementType type) noexcept {
    return element_dimension(type);
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

// Reference-space dimension of a basis topology: 0 for points up to 3 for
// volume topologies; -1 for Unknown.
[[nodiscard]] constexpr int topology_dimension(BasisTopology top) noexcept {
    switch (top) {
        case BasisTopology::Point:         return 0;
        case BasisTopology::Line:          return 1;
        case BasisTopology::Triangle:
        case BasisTopology::Quadrilateral: return 2;
        case BasisTopology::Tetrahedron:
        case BasisTopology::Hexahedron:
        case BasisTopology::Wedge:         return 3;
        default:                           return -1;
    }
}

// Lowest-order named element that represents a topology. Used internally to
// drive the reference-node generators, which key on a canonical ElementType
// (and re-canonicalize it). This is the inverse of topology() for the linear
// elements and is purely an implementation detail: the node-count name never
// leaks into the public basis identity.
[[nodiscard]] constexpr ElementType lagrange_topology_representative(BasisTopology top) noexcept {
    switch (top) {
        case BasisTopology::Point:         return ElementType::Point1;
        case BasisTopology::Line:          return ElementType::Line2;
        case BasisTopology::Triangle:      return ElementType::Triangle3;
        case BasisTopology::Quadrilateral: return ElementType::Quad4;
        case BasisTopology::Tetrahedron:   return ElementType::Tetra4;
        case BasisTopology::Hexahedron:    return ElementType::Hex8;
        case BasisTopology::Wedge:         return ElementType::Wedge6;
        default:                           return ElementType::Unknown;
    }
}

// Polynomial order baked into a named Lagrange element layout: 0 for the point,
// 1 for the linear elements, 2 for the complete-quadratic aliases; -1 for types
// with no complete-Lagrange order (serendipity, pyramid, Unknown). Unlike
// complete_lagrange_alias_order this also maps Point1 -> 0, so it is the single
// source of truth the (ElementType, order) constructor validates against.
[[nodiscard]] constexpr int named_lagrange_order(ElementType type) noexcept {
    if (type == ElementType::Point1) {
        return 0;
    }
    return complete_lagrange_alias_order(type);
}

// Inverse of (topology(), order()) for the named layouts: the ElementType that a
// (topology, order, family) triple denotes, or Unknown when no named layout
// exists (order 0 on a non-point topology, any order >= 3, or a reduced family
// at an unsupported order). topology() + order() remain the authoritative
// identity; callers that want a named ElementType for a basis pass its
// topology(), order(), and basis_type() to this free helper directly.
[[nodiscard]] constexpr ElementType named_element_for(BasisTopology top, int order,
                                                      BasisType family) noexcept {
    if (family == BasisType::Serendipity) {
        switch (top) {
            case BasisTopology::Quadrilateral:
                return order == 2 ? ElementType::Quad8 : ElementType::Unknown;
            case BasisTopology::Hexahedron:
                if (order == 1) { return ElementType::Hex8; }
                if (order == 2) { return ElementType::Hex20; }
                return ElementType::Unknown;
            case BasisTopology::Wedge:
                return order == 2 ? ElementType::Wedge15 : ElementType::Unknown;
            default:
                return ElementType::Unknown;
        }
    }

    // Lagrange (and any nodal family built on the complete layouts).
    if (top == BasisTopology::Point) {
        return order == 0 ? ElementType::Point1 : ElementType::Unknown;
    }
    switch (order) {
        case 1:
            switch (top) {
                case BasisTopology::Line:          return ElementType::Line2;
                case BasisTopology::Triangle:      return ElementType::Triangle3;
                case BasisTopology::Quadrilateral: return ElementType::Quad4;
                case BasisTopology::Tetrahedron:   return ElementType::Tetra4;
                case BasisTopology::Hexahedron:    return ElementType::Hex8;
                case BasisTopology::Wedge:         return ElementType::Wedge6;
                default:                           return ElementType::Unknown;
            }
        case 2:
            switch (top) {
                case BasisTopology::Line:          return ElementType::Line3;
                case BasisTopology::Triangle:      return ElementType::Triangle6;
                case BasisTopology::Quadrilateral: return ElementType::Quad9;
                case BasisTopology::Tetrahedron:   return ElementType::Tetra10;
                case BasisTopology::Hexahedron:    return ElementType::Hex27;
                case BasisTopology::Wedge:         return ElementType::Wedge18;
                default:                           return ElementType::Unknown;
            }
        default:
            return ElementType::Unknown;
    }
}

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_BASISTRAITS_H
