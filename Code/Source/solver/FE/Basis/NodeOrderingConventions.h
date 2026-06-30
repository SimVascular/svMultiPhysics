// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_BASIS_NODEORDERINGCONVENTIONS_H
#define SVMP_FE_BASIS_NODEORDERINGCONVENTIONS_H

#include "BasisTraits.h"
#include "Math/Vector.h"
#include "Types.h"

#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace svmp::FE::basis {

/**
 * @defgroup FE_BasisNodeOrdering Reference-node generation (internal)
 * @ingroup FE_Basis
 * @brief Reference-node generators that the basis families build on.
 *
 * @warning Internal implementation detail. Do not use these directly: obtain a
 * basis through basis_factory and read its nodes via BasisFunction::nodes().
 * These declarations are part of the internal node-ordering machinery and their
 * interface may change without notice.
 *
 * @details This is the reference-node generator the basis families build on, not
 * a consumer entry point. It is documented for FE core developers; model-level
 * code never calls it directly.
 * @{
 */

/**
 * @brief The i-th 1D tensor-axis reference node on [-1, 1] at the given order.
 *
 * @details Returns the Gauss-Lobatto-Legendre (GLL) node of index @p i for a
 * degree-@p order distribution: the endpoints are -1 and +1 and the interior
 * nodes are the roots of @f$P'_{order}@f$, so high-order tensor interpolation
 * stays well-conditioned (a logarithmic Lebesgue constant instead of the
 * exponential growth of equispaced nodes). At order 1 the nodes are
 * @f$\{-1, +1\}@f$ and at order 2 @f$\{-1, 0, +1\}@f$, so they coincide with the
 * equispaced layout for the production orders and differ only for order >= 3.
 * Returns 0 for order <= 0 when @p i is 0. Invalid indices throw.
 *
 * This is the single definition of the tensor-axis node distribution: the
 * reference-node layout generators, the Lagrange tensor-axis initialization, and
 * the serendipity edge/face/interior strata all source their 1D nodes here. The
 * LagrangeBasis and SerendipityBasis docs point back to this description of the
 * GLL distribution and its conditioning rather than restating it.
 *
 * @param i Node index in [0, order] for positive orders, or 0 for order <= 0.
 * @param order Polynomial order of the 1D distribution.
 * @return GLL node coordinate on [-1, 1].
 * @throws BasisNodeOrderingException If @p i is outside the valid range.
 */
[[nodiscard]] double line_coord_pm_one(int i, int order);

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
    std::vector<math::Vector<double, 3>> coords;   ///< Reference node coordinates, one per node.
    std::vector<std::array<int, 3>>      lattice;  ///< Integer lattice index of each node (see details above).
};

/**
 * @brief Reference-node coordinate and count lookups for an element type.
 */
class ReferenceNodeLayout {
public:
    /**
     * @brief One reference node coordinate by local index. Regenerates the full
     * layout per call; prefer node_coords() when more than one node is needed.
     *
     * @param elem_type Element type to look up.
     * @param local_node Local node index in [0, num_nodes(elem_type)).
     * @return Reference coordinate of the requested node.
     */
    static math::Vector<double, 3> node_coord_at(ElementType elem_type,
                                                 std::size_t local_node);

    /**
     * @brief Number of reference nodes in an element type's public layout.
     * @param elem_type Element type to look up.
     * @return Node count.
     */
    static std::size_t num_nodes(ElementType elem_type);

    /**
     * @brief All reference node coordinates for an element type, in public layout order.
     *
     * @details Returns the complete public reference layout for @p elem_type
     * (the same coordinates node_coord_at() returns one at a time), including
     * the serendipity layouts. Prefer this single call when the whole layout is
     * needed: node_coord_at() regenerates the full list on every call.
     *
     * @param elem_type Element type to look up.
     * @return Reference node coordinates, one per node.
     */
    static std::vector<math::Vector<double, 3>> node_coords(ElementType elem_type);

    /**
     * @brief Reference Lagrange node coordinates for a canonical type and order.
     * @param canonical_type Canonical Lagrange element type (or Point1).
     * @param order Polynomial order.
     * @return Reference node coordinates, one per node, in basis order.
     */
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

    /**
     * @brief Reference nodes for an arbitrary-order serendipity layout.
     *
     * @details Generates the stratified serendipity node set for the
     * quadrilateral or hexahedral family at the requested order: the
     * corner+edge skeleton (the leading prefix of the complete Lagrange layout
     * of the same order, in VTK boundary order) followed by the reduced face
     * and volume interior. This is the single source of serendipity node
     * geometry -- SerendipityBasis builds its mode space and coefficient table
     * on top of these coordinates for both the arbitrary-order path and the
     * named Quad8/Hex20 layouts (the order-2 instances). Wedge serendipity
     * (Wedge15) is a fixed named layout and is not generated here.
     *
     * @param topology BasisTopology::Quadrilateral or BasisTopology::Hexahedron.
     * @param order Polynomial order; must be >= 1.
     * @return Reference node coordinates in stratified (skeleton-then-interior) order.
     * @throws BasisConstructionException If @p order is below 1.
     * @throws BasisElementCompatibilityException If @p topology is not Quadrilateral or Hexahedron.
     */
    static std::vector<math::Vector<double, 3>>
    serendipity_node_coords(BasisTopology topology, int order);
};

/** @} */

} // namespace svmp::FE::basis

#endif // SVMP_FE_BASIS_NODEORDERINGCONVENTIONS_H
