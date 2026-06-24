// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_BASIS_NODEORDERINGCONVENTIONS_H
#define SVMP_FE_BASIS_NODEORDERINGCONVENTIONS_H

#include "Math/Vector.h"
#include "Types.h"

#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace svmp {
namespace FE {
namespace basis {

/**
 * @brief Equispaced 1D reference coordinate on [-1, 1]: -1 + 2 i / order.
 *
 * Shared by the reference-node layout generators and the Lagrange tensor-axis
 * node initialization so the lattice formula lives in a single place.
 */
[[nodiscard]] inline constexpr double line_coord_pm_one(int i, int order) noexcept {
    if (order <= 0) {
        return double(0);
    }
    return double(-1) + double(2) * static_cast<double>(i) / static_cast<double>(order);
}

/**
 * @brief Reference Lagrange node coordinates paired with their integer lattice
 * index.
 *
 * @details `lattice[n]` is the exact integer index of `coords[n]` in the
 * element's natural index space, with every component in `[0, order]`:
 * - tensor topologies (line/quad/hex): axis indices `(i, j, k)`, unused axes `0`;
 * - simplex topologies (triangle/tetra): off-origin barycentric indices
 *   `(i, j, k)` (with `k = 0` for triangles) satisfying `i + j + k <= order`;
 * - wedge: triangle lattice `(i, j)` in the first two components and the
 *   through-axis index `r` in the third.
 *
 * Emitting the lattice alongside the coordinate lets callers consume the integer
 * index directly instead of reconstructing it from the floating-point coordinate.
 */
struct LagrangeNodeLayout {
    std::vector<math::Vector<double, 3>> coords;
    std::vector<std::array<int, 3>>      lattice;
};

class ReferenceNodeLayout {
public:
    static math::Vector<double, 3> get_node_coords(ElementType elem_type,
                                                 std::size_t local_node);
    static std::size_t num_nodes(ElementType elem_type);

    /**
     * @brief All reference node coordinates for an element type, in public layout order.
     *
     * @details Returns the complete public reference layout for @p elem_type
     * (the same coordinates get_node_coords() returns one at a time), including
     * the serendipity layouts. Prefer this single call when the whole layout is
     * needed: get_node_coords() regenerates the full list on every call.
     *
     * @param elem_type Element type to look up.
     * @return Reference node coordinates, one per node.
     */
    static std::vector<math::Vector<double, 3>> node_coords(ElementType elem_type);

    static std::vector<math::Vector<double, 3>>
    get_lagrange_node_coords(ElementType canonical_type, int order);

    /**
     * @brief Reference Lagrange nodes with their integer lattice indices.
     *
     * @details Returns the same coordinates as get_lagrange_node_coords(), paired
     * with the integer lattice index of each node (see LagrangeNodeLayout). The
     * structural invariants in the contract (size match, components in
     * `[0, order]`, simplex/wedge sum bounds) are validated before returning.
     *
     * @param canonical_type Canonical Lagrange element type (or Point1).
     * @param order Polynomial order.
     * @return Coordinates and matching lattice indices, one entry per node.
     * @throws BasisConstructionException If a structural invariant is violated.
     */
    static LagrangeNodeLayout
    get_lagrange_lattice(ElementType canonical_type, int order);

    static std::span<const std::size_t> mesh_to_basis_ordering(ElementType elem_type);
};

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_NODEORDERINGCONVENTIONS_H
