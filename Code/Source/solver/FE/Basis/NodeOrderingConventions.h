/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#ifndef SVMP_FE_BASIS_NODEORDERINGCONVENTIONS_H
#define SVMP_FE_BASIS_NODEORDERINGCONVENTIONS_H

#include "Types.h"
#include "Math/Vector.h"
#include <cstddef>
#include <vector>

/**
 * @file NodeOrderingConventions.h
 * @brief Documentation of node ordering conventions for all element types
 *
 * This file provides comprehensive documentation of the node ordering
 * conventions used throughout the FE library. These orderings are consistent
 * with VTK conventions and must be matched exactly when interfacing with
 * the Mesh library.
 *
 * IMPORTANT: The FE library (Basis, Quadrature, Geometry) uses "node" to refer
 * to degrees of freedom locations on reference elements. The Mesh library uses
 * "vertex" for geometry vertices and "cell" for mesh elements. When interfacing
 * between the two, ensure consistent ordering.
 *
 * Reference Element Conventions:
 * - Line:       xi in [-1, 1]
 * - Quad:       (xi, eta) in [-1, 1] x [-1, 1]
 * - Hex:        (xi, eta, zeta) in [-1, 1]^3
 * - Triangle:   (xi, eta) in simplex with vertices (0,0), (1,0), (0,1)
 * - Tetrahedron: (xi, eta, zeta) in simplex with vertices
 *                (0,0,0), (1,0,0), (0,1,0), (0,0,1)
 * - Wedge:      Triangle base x line height, zeta in [-1, 1]
 * - Pyramid:    Quad base at z=0, apex at (0, 0, 1)
 *
 *
 * =============================================================================
 * 1D ELEMENTS
 * =============================================================================
 *
 * Line2 (Linear Line)
 * -------------------
 *   0---------1
 *   |         |
 *  xi=-1     xi=+1
 *
 * Node 0: xi = -1
 * Node 1: xi = +1
 *
 *
 * Line3 (Quadratic Line)
 * ----------------------
 *   0----2----1
 *   |    |    |
 *  xi=-1 0   xi=+1
 *
 * Node 0: xi = -1
 * Node 1: xi = +1
 * Node 2: xi =  0 (mid-edge)
 *
 *
 * =============================================================================
 * 2D QUADRILATERAL ELEMENTS
 * =============================================================================
 *
 * Quad4 (Bilinear Quadrilateral)
 * ------------------------------
 *
 *   3-----------2
 *   |           |
 *   |           |
 *   |           |
 *   0-----------1
 *
 * Node 0: (xi, eta) = (-1, -1)
 * Node 1: (xi, eta) = (+1, -1)
 * Node 2: (xi, eta) = (+1, +1)
 * Node 3: (xi, eta) = (-1, +1)
 *
 *
 * Quad8 (Serendipity Quadrilateral)
 * ---------------------------------
 *
 *   3-----6-----2
 *   |           |
 *   7           5
 *   |           |
 *   0-----4-----1
 *
 * Corners (same as Quad4):
 *   Node 0: (-1, -1)
 *   Node 1: (+1, -1)
 *   Node 2: (+1, +1)
 *   Node 3: (-1, +1)
 *
 * Mid-edge nodes:
 *   Node 4: ( 0, -1)  (edge 0-1)
 *   Node 5: (+1,  0)  (edge 1-2)
 *   Node 6: ( 0, +1)  (edge 2-3)
 *   Node 7: (-1,  0)  (edge 3-0)
 *
 *
 * Quad9 (Biquadratic Quadrilateral)
 * ---------------------------------
 *
 *   3-----6-----2
 *   |           |
 *   7     8     5
 *   |           |
 *   0-----4-----1
 *
 * Same as Quad8 plus:
 *   Node 8: (0, 0)  (center)
 *
 *
 * =============================================================================
 * 3D HEXAHEDRAL ELEMENTS
 * =============================================================================
 *
 * Hex8 (Trilinear Hexahedron)
 * ---------------------------
 *
 *        7-----------6
 *       /|          /|
 *      / |         / |
 *     4-----------5  |
 *     |  |        |  |
 *     |  3--------|--2
 *     | /         | /
 *     |/          |/
 *     0-----------1
 *
 * Bottom face (zeta = -1):
 *   Node 0: (xi, eta, zeta) = (-1, -1, -1)
 *   Node 1: (xi, eta, zeta) = (+1, -1, -1)
 *   Node 2: (xi, eta, zeta) = (+1, +1, -1)
 *   Node 3: (xi, eta, zeta) = (-1, +1, -1)
 *
 * Top face (zeta = +1):
 *   Node 4: (xi, eta, zeta) = (-1, -1, +1)
 *   Node 5: (xi, eta, zeta) = (+1, -1, +1)
 *   Node 6: (xi, eta, zeta) = (+1, +1, +1)
 *   Node 7: (xi, eta, zeta) = (-1, +1, +1)
 *
 *
 * Hex20 (Serendipity Hexahedron)
 * ------------------------------
 *
 *        7-----14-----6
 *       /|           /|
 *     15 |         13 |
 *     /  19        /  18
 *    4-----12-----5   |
 *    |   |        |   |
 *    |   3-----10-|---2
 *   16  /        17  /
 *    | 11         | 9
 *    |/           |/
 *    0------8-----1
 *
 * Corners (same as Hex8): Nodes 0-7
 *
 * Mid-edge nodes on bottom face (zeta = -1):
 *   Node 8:  ( 0, -1, -1)  (edge 0-1)
 *   Node 9:  (+1,  0, -1)  (edge 1-2)
 *   Node 10: ( 0, +1, -1)  (edge 2-3)
 *   Node 11: (-1,  0, -1)  (edge 3-0)
 *
 * Mid-edge nodes on top face (zeta = +1):
 *   Node 12: ( 0, -1, +1)  (edge 4-5)
 *   Node 13: (+1,  0, +1)  (edge 5-6)
 *   Node 14: ( 0, +1, +1)  (edge 6-7)
 *   Node 15: (-1,  0, +1)  (edge 7-4)
 *
 * Mid-edge nodes on vertical edges:
 *   Node 16: (-1, -1,  0)  (edge 0-4)
 *   Node 17: (+1, -1,  0)  (edge 1-5)
 *   Node 18: (+1, +1,  0)  (edge 2-6)
 *   Node 19: (-1, +1,  0)  (edge 3-7)
 *
 *
 * Hex27 (Triquadratic Hexahedron)
 * -------------------------------
 * Same as Hex20 plus face-center and body-center nodes:
 *
 * Face centers:
 *   Node 20: ( 0,  0, -1)  (bottom face)
 *   Node 21: ( 0,  0, +1)  (top face)
 *   Node 22: ( 0, -1,  0)  (front face)
 *   Node 23: (+1,  0,  0)  (right face)
 *   Node 24: ( 0, +1,  0)  (back face)
 *   Node 25: (-1,  0,  0)  (left face)
 *
 * Body center:
 *   Node 26: (0, 0, 0)
 *
 *
 * =============================================================================
 * 2D TRIANGULAR ELEMENTS
 * =============================================================================
 *
 * Triangle3 (Linear Triangle)
 * ---------------------------
 *
 *   2
 *   |\
 *   | \
 *   |  \
 *   |   \
 *   0----1
 *
 * Reference: (xi, eta) simplex with vertices at:
 *   Node 0: (xi, eta) = (0, 0)
 *   Node 1: (xi, eta) = (1, 0)
 *   Node 2: (xi, eta) = (0, 1)
 *
 *
 * Triangle6 (Quadratic Triangle)
 * ------------------------------
 *
 *   2
 *   |\
 *   | \
 *   5  4
 *   |   \
 *   0--3--1
 *
 * Corners: Nodes 0-2 (same as Triangle3)
 *
 * Mid-edge nodes:
 *   Node 3: (0.5,   0)  (edge 0-1)
 *   Node 4: (0.5, 0.5)  (edge 1-2)
 *   Node 5: (  0, 0.5)  (edge 2-0)
 *
 *
 * =============================================================================
 * 3D TETRAHEDRAL ELEMENTS
 * =============================================================================
 *
 * Tetrahedron4 (Linear Tetrahedron)
 * ---------------------------------
 *
 *             3
 *            /|\
 *           / | \
 *          /  |  \
 *         /   |   \
 *        /    |    \
 *       0-----|-----2
 *        \    |    /
 *         \   |   /
 *          \  |  /
 *           \ | /
 *            \|/
 *             1
 *
 * Reference: (xi, eta, zeta) simplex with vertices at:
 *   Node 0: (0, 0, 0)
 *   Node 1: (1, 0, 0)
 *   Node 2: (0, 1, 0)
 *   Node 3: (0, 0, 1)
 *
 *
 * Tetrahedron10 (Quadratic Tetrahedron)
 * -------------------------------------
 * Corners: Nodes 0-3 (same as Tet4)
 *
 * Mid-edge nodes:
 *   Node 4: (0.5,   0,   0)  (edge 0-1)
 *   Node 5: (0.5, 0.5,   0)  (edge 1-2)
 *   Node 6: (  0, 0.5,   0)  (edge 2-0)
 *   Node 7: (  0,   0, 0.5)  (edge 0-3)
 *   Node 8: (0.5,   0, 0.5)  (edge 1-3)
 *   Node 9: (  0, 0.5, 0.5)  (edge 2-3)
 *
 *
 * =============================================================================
 * 3D WEDGE (PRISM) ELEMENTS
 * =============================================================================
 *
 * Wedge6 (Linear Wedge)
 * ---------------------
 *
 *         5
 *        /|\
 *       / | \
 *      /  |  \
 *     3---|---4
 *     |   2   |
 *     |  / \  |
 *     | /   \ |
 *     |/     \|
 *     0-------1
 *
 * Reference: Triangle base at zeta = -1, top at zeta = +1
 *
 * Bottom face (zeta = -1):
 *   Node 0: (0, 0, -1)
 *   Node 1: (1, 0, -1)
 *   Node 2: (0, 1, -1)
 *
 * Top face (zeta = +1):
 *   Node 3: (0, 0, +1)
 *   Node 4: (1, 0, +1)
 *   Node 5: (0, 1, +1)
 *
 *
 * Wedge15 (Quadratic Wedge)
 * -------------------------
 * Corners: Nodes 0-5 (same as Wedge6)
 *
 * Mid-edge nodes on bottom face:
 *   Node 6:  (0.5,   0, -1)  (edge 0-1)
 *   Node 7:  (0.5, 0.5, -1)  (edge 1-2)
 *   Node 8:  (  0, 0.5, -1)  (edge 2-0)
 *
 * Mid-edge nodes on top face:
 *   Node 9:  (0.5,   0, +1)  (edge 3-4)
 *   Node 10: (0.5, 0.5, +1)  (edge 4-5)
 *   Node 11: (  0, 0.5, +1)  (edge 5-3)
 *
 * Mid-edge nodes on vertical edges:
 *   Node 12: (0, 0, 0)  (edge 0-3)
 *   Node 13: (1, 0, 0)  (edge 1-4)
 *   Node 14: (0, 1, 0)  (edge 2-5)
 *
 *
 * Wedge18 (Complete Quadratic Wedge)
 * ----------------------------------
 * Corners and mid-edges: Nodes 0-14 (same as Wedge15)
 *
 * Face-center nodes on quadrilateral faces:
 *   Node 15: (0.5, 0.0, 0.0)  (face with vertices 0-1-4-3, y = 0)
 *   Node 16: (0.5, 0.5, 0.0)  (face with vertices 1-2-5-4, x + y = 1)
 *   Node 17: (0.0, 0.5, 0.0)  (face with vertices 2-0-3-5, x = 0)
 *
 *
 * =============================================================================
 * 3D PYRAMID ELEMENTS
 * =============================================================================
 *
 * Pyramid5 (Linear Pyramid)
 * -------------------------
 *
 *           4
 *          /|\
 *         / | \
 *        /  |  \
 *       /   |   \
 *      3----|----2
 *      |    |    |
 *      |    +    |   (apex projects to center of base)
 *      |         |
 *      0---------1
 *
 * Reference: Quad base in xi-eta plane at zeta = 0, apex at zeta = 1
 *
 * Base (zeta = 0):
 *   Node 0: (-1, -1, 0)
 *   Node 1: (+1, -1, 0)
 *   Node 2: (+1, +1, 0)
 *   Node 3: (-1, +1, 0)
 *
 * Apex:
 *   Node 4: (0, 0, 1)
 *
 *
 * Pyramid13 (Quadratic Pyramid)
 * -----------------------------
 * Corners: Nodes 0-4 (same as Pyramid5)
 *
 * Mid-edge nodes on base:
 *   Node 5: ( 0, -1, 0)  (edge 0-1)
 *   Node 6: (+1,  0, 0)  (edge 1-2)
 *   Node 7: ( 0, +1, 0)  (edge 2-3)
 *   Node 8: (-1,  0, 0)  (edge 3-0)
 *
 * Mid-edge nodes to apex:
 *   Node 9:  (-0.5, -0.5, 0.5)  (edge 0-4)
 *   Node 10: (+0.5, -0.5, 0.5)  (edge 1-4)
 *   Node 11: (+0.5, +0.5, 0.5)  (edge 2-4)
 *   Node 12: (-0.5, +0.5, 0.5)  (edge 3-4)
 *
 *
 * Pyramid14 (Quadratic Rational Pyramid)
 * --------------------------------------
 *
 * This retained low-order compatibility layout matches the generated
 * complete-family quadratic Lagrange ordering for the reference pyramid with
 * base (-1,-1,0)..(1,1,0) and apex at (0,0,1). Nodes 0-12 coincide with the
 * Pyramid13 layout; node 13 is the base center.
 *
 *   Base corners (same as Pyramid5):
 *     Node 0: (-1, -1, 0)
 *     Node 1: (+1, -1, 0)
 *     Node 2: (+1, +1, 0)
 *     Node 3: (-1, +1, 0)
 *
 *   Apex:
 *     Node 4: (0, 0, 1)
 *
 *   Base mid-edges (same as Pyramid13):
 *     Node 5:  ( 0, -1, 0)   (edge 0-1)
 *     Node 6:  (+1,  0, 0)   (edge 1-2)
 *     Node 7:  ( 0, +1, 0)   (edge 2-3)
 *     Node 8:  (-1,  0, 0)   (edge 3-0)
 *
 *   Mid-edges to apex (same as Pyramid13):
 *     Node 9:  (-0.5, -0.5, 0.5)  (edge 0-4)
 *     Node 10: (+0.5, -0.5, 0.5)  (edge 1-4)
 *     Node 11: (+0.5, +0.5, 0.5)  (edge 2-4)
 *     Node 12: (-0.5, +0.5, 0.5)  (edge 3-4)
 *
 *   Base center:
 *     Node 13: (0, 0, 0)
 *
 *
 * =============================================================================
 * NOTES ON VTK COMPATIBILITY
 * =============================================================================
 *
 * The node orderings above are consistent with VTK cell types:
 *
 *   VTK_LINE           (3)  -> Line2
 *   VTK_QUADRATIC_EDGE (21) -> Line3
 *   VTK_TRIANGLE       (5)  -> Triangle3
 *   VTK_QUADRATIC_TRIANGLE (22) -> Triangle6
 *   VTK_QUAD           (9)  -> Quad4
 *   VTK_QUADRATIC_QUAD (23) -> Quad8
 *   VTK_BIQUADRATIC_QUAD (28) -> Quad9
 *   VTK_TETRA          (10) -> Tetrahedron4
 *   VTK_QUADRATIC_TETRA (24) -> Tetrahedron10
 *   VTK_HEXAHEDRON     (12) -> Hex8
 *   VTK_QUADRATIC_HEXAHEDRON (25) -> Hex20
 *   VTK_TRIQUADRATIC_HEXAHEDRON (29) -> Hex27
 *   VTK_WEDGE          (13) -> Wedge6
 *   VTK_QUADRATIC_WEDGE (26) -> Wedge15
 *   VTK_BIQUADRATIC_QUADRATIC_WEDGE (32) -> Wedge18
 *   VTK_PYRAMID        (14) -> Pyramid5
 *   VTK_QUADRATIC_PYRAMID (27) -> Pyramid13
 *
 *
 * =============================================================================
 * BARYCENTRIC COORDINATES
 * =============================================================================
 *
 * For simplex elements, barycentric coordinates (lambda_0, ..., lambda_n)
 * satisfy sum(lambda_i) = 1.
 *
 * Triangle:
 *   lambda_0 = 1 - xi - eta
 *   lambda_1 = xi
 *   lambda_2 = eta
 *
 * Tetrahedron:
 *   lambda_0 = 1 - xi - eta - zeta
 *   lambda_1 = xi
 *   lambda_2 = eta
 *   lambda_3 = zeta
 *
 */

#include <span>

namespace svmp {
namespace FE {
namespace basis {

/**
 * @brief Basis-side reference node coordinate queries
 *
 * This is intentionally named differently from `svmp::NodeOrdering` in Mesh,
 * which handles mesh-format permutations rather than reference basis layouts.
 */
class ReferenceNodeLayout {
public:
    /**
     * @brief Get reference coordinates for a node
     * @param elem_type Element type
     * @param local_node Local node index (0-based)
     *
     * Complete-family low-order Lagrange aliases (`Line2/3`, `Triangle3/6`,
     * `Quad4/9`, `Tetra4/10`, `Hex8/27`, `Wedge6/18`, `Pyramid5/14`) are
     * served by the generated arbitrary-order Lagrange ordering path. Explicit
     * hard-coded tables remain only for serendipity-only enums such as
     * `Quad8`, `Hex20`, `Wedge15`, and `Pyramid13`.
     *
     * @return Reference coordinates (xi, eta, zeta)
     */
    static math::Vector<Real, 3> get_node_coords(ElementType elem_type, std::size_t local_node);

    /**
     * @brief Get number of nodes for an element type
     *
     * The low-order complete-family Lagrange aliases share the same generated
     * ordering path used by `get_node_coords`.
     */
    static std::size_t num_nodes(ElementType elem_type);

    /**
     * @brief Generate complete-family Lagrange node coordinates for a canonical topology and order
     *
     * This covers arbitrary-order complete nodal Lagrange spaces on the
     * canonical topologies `Line2`, `Triangle3`, `Quad4`, `Tetra4`, `Hex8`,
     * `Wedge6`, and `Pyramid5`. Serendipity variants are intentionally
     * excluded.
     */
    static std::vector<math::Vector<Real, 3>>
    get_lagrange_node_coords(ElementType canonical_type, int order);

    /**
     * @brief Optional mapping from mesh/reference node order to internal basis order
     *
     * Returns an empty span when the public node order is already the basis
     * table order or no special mapping is registered.
     */
    static std::span<const std::size_t> mesh_to_basis_ordering(ElementType elem_type);

    /**
     * @brief Check if element is a simplex (triangle, tetrahedron)
     */
    static bool is_simplex(ElementType elem_type);

    /**
     * @brief Check if element uses tensor-product topology
     */
    static bool is_tensor_product(ElementType elem_type);
};

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_NODEORDERINGCONVENTIONS_H
