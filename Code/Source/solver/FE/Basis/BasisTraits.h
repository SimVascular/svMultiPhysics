// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_BASIS_BASISTRAITS_H
#define SVMP_FE_BASIS_BASISTRAITS_H

/**
 * @file BasisTraits.h
 * @brief Reference-topology vocabulary (BasisTopology) and the internal
 *        ElementType/topology/order maps.
 */

#include "Types.h"

#include <cstddef>

namespace svmp::FE::basis {

/**
 * @brief Reference-cell topology of a basis (the shape, independent of order).
 * @ingroup FE_Basis
 *
 * @details Together with a polynomial order this is the order-agnostic identity a
 * basis is built from: the arbitrary-order constructors take a BasisTopology and an
 * order, and BasisRequest::topology selects that path. A named ElementType maps to
 * one of these through topology().
 */
enum class BasisTopology {
    Unknown,        ///< Unrecognized or uninitialized topology.
    Point,          ///< 0D point.
    Line,           ///< 1D line segment.
    Triangle,       ///< 2D triangle (simplex).
    Quadrilateral,  ///< 2D quadrilateral (tensor product).
    Tetrahedron,    ///< 3D tetrahedron (simplex).
    Hexahedron,     ///< 3D hexahedron (tensor product).
    Wedge,          ///< 3D triangular prism.
};

// The maps below are internal to the Basis module (used to build the basis classes
// and the factory); they are excluded from the public Doxygen output.
/** @cond INTERNAL */

// ---------------------------------------------------------------------------
// ElementType / BasisTopology / order mapping helpers.
//
// A basis identity is expressed three ways -- a named ElementType, a
// (BasisTopology, order) pair, and a reference dimension -- and the constexpr
// maps below convert between them. They are grouped here so the relationships
// stay in one place:
//
//   ElementType   -> BasisTopology    topology()
//   ElementType   -> ElementType      canonical_lagrange_type()   (alias -> linear representative)
//   ElementType   -> order            complete_lagrange_alias_order(), named_lagrange_order()
//   BasisTopology -> int (dimension)  topology_dimension()
//   BasisTopology -> ElementType      lagrange_topology_representative() (lowest-order representative)
//   (BasisTopology, order, family) -> ElementType   named_element_for() (inverse of topology() + order())
//
// The two ElementType -> order maps differ only at Point1:
// complete_lagrange_alias_order() returns -1 (a point is not a complete-Lagrange
// alias) while named_lagrange_order() returns 0 (the point layout's order).
// named_lagrange_order() is defined in terms of complete_lagrange_alias_order(),
// so the order-1 / order-2 alias values have a single source of truth.
// ---------------------------------------------------------------------------

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
        // BasisTopology::Unknown is a sentinel the basis constructors validate
        // and convert into a BasisElementCompatibilityException at the call site,
        // not an error raised from this constexpr noexcept classifier.
        default:                   return BasisTopology::Unknown;
    }
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
            // -1 is a sentinel for "not a complete-Lagrange alias" (serendipity
            // layouts, pyramids, Unknown), not an error: the LagrangeBasis
            // (ElementType, order) constructor compares the requested order
            // against named_lagrange_order() and raises BasisConfigurationException
            // on mismatch. These classifiers are constexpr noexcept and so cannot
            // throw themselves.
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

/** @endcond */

/**
 * @brief Named ElementType denoted by a (topology, order, family) triple.
 * @ingroup FE_Basis
 *
 * @details Inverse of topology() + order() for the named layouts: returns the
 * ElementType a basis identity denotes, or ElementType::Unknown when no named
 * layout exists (order 0 on a non-point topology, any order >= 3, or a reduced
 * family at an unsupported order). topology() + order() remain the authoritative
 * identity; callers that want a named ElementType for a basis pass its topology(),
 * order(), and basis_type() here.
 *
 * @param top Reference topology.
 * @param order Polynomial order.
 * @param family Basis family; only Serendipity is distinguished from nodal/Lagrange naming.
 * @return Named ElementType, or ElementType::Unknown when none applies.
 */
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

} // namespace svmp::FE::basis

#endif // SVMP_FE_BASIS_BASISTRAITS_H
