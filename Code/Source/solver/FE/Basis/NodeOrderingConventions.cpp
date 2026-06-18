// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "NodeOrderingConventions.h"
#include "BasisExceptions.h"
#include "BasisTraits.h"

#include <array>

namespace svmp {
namespace FE {
namespace basis {

namespace {

using Point = math::Vector<Real, 3>;

constexpr std::array<std::size_t, 20> kHex20MeshToBasisOrder = {
    0, 1, 2, 3, 4, 5, 6, 7,
    8, 13, 10, 12,
    9, 15, 11, 14,
    16, 17, 19, 18
};

Real line_coord_zero_one(int i, int order) {
    if (order <= 0) {
        return Real(0);
    }
    return static_cast<Real>(i) / static_cast<Real>(order);
}

void append_triangle_face_interior(std::vector<Point>& nodes,
                                   const Point& v0,
                                   const Point& v1,
                                   const Point& v2,
                                   int order) {
    for (int c = 1; c <= order - 2; ++c) {
        for (int b = 1; b <= order - c - 1; ++b) {
            const int a = order - b - c;
            const Real inv = Real(1) / Real(order);
            nodes.push_back(v0 * (Real(a) * inv) +
                            v1 * (Real(b) * inv) +
                            v2 * (Real(c) * inv));
        }
    }
}

std::vector<Point> generate_line_nodes(int order) {
    if (order == 0) {
        return {Point{Real(0), Real(0), Real(0)}};
    }

    std::vector<Point> nodes;
    nodes.reserve(static_cast<std::size_t>(order + 1));
    nodes.push_back(Point{Real(-1), Real(0), Real(0)});
    nodes.push_back(Point{Real(1), Real(0), Real(0)});
    for (int i = 1; i < order; ++i) {
        nodes.push_back(Point{line_coord_pm_one(i, order), Real(0), Real(0)});
    }
    return nodes;
}

std::vector<Point> generate_triangle_nodes(int order) {
    if (order == 0) {
        return {Point{Real(1) / Real(3), Real(1) / Real(3), Real(0)}};
    }

    std::vector<Point> nodes;
    nodes.reserve(static_cast<std::size_t>((order + 1) * (order + 2) / 2));
    nodes.push_back(Point{Real(0), Real(0), Real(0)});
    nodes.push_back(Point{Real(1), Real(0), Real(0)});
    nodes.push_back(Point{Real(0), Real(1), Real(0)});

    for (int m = 1; m < order; ++m) {
        nodes.push_back(Point{line_coord_zero_one(m, order), Real(0), Real(0)});
    }
    for (int m = 1; m < order; ++m) {
        nodes.push_back(Point{line_coord_zero_one(order - m, order),
                              line_coord_zero_one(m, order), Real(0)});
    }
    for (int m = 1; m < order; ++m) {
        nodes.push_back(Point{Real(0), line_coord_zero_one(order - m, order), Real(0)});
    }

    append_triangle_face_interior(nodes,
                                  Point{Real(0), Real(0), Real(0)},
                                  Point{Real(1), Real(0), Real(0)},
                                  Point{Real(0), Real(1), Real(0)},
                                  order);
    return nodes;
}

std::vector<Point> generate_quad_nodes(int order) {
    if (order == 0) {
        return {Point{Real(0), Real(0), Real(0)}};
    }

    std::vector<Point> nodes;
    nodes.reserve(static_cast<std::size_t>((order + 1) * (order + 1)));
    nodes.push_back(Point{Real(-1), Real(-1), Real(0)});
    nodes.push_back(Point{Real(1), Real(-1), Real(0)});
    nodes.push_back(Point{Real(1), Real(1), Real(0)});
    nodes.push_back(Point{Real(-1), Real(1), Real(0)});

    for (int i = 1; i < order; ++i) {
        nodes.push_back(Point{line_coord_pm_one(i, order), Real(-1), Real(0)});
    }
    for (int j = 1; j < order; ++j) {
        nodes.push_back(Point{Real(1), line_coord_pm_one(j, order), Real(0)});
    }
    for (int i = order - 1; i >= 1; --i) {
        nodes.push_back(Point{line_coord_pm_one(i, order), Real(1), Real(0)});
    }
    for (int j = order - 1; j >= 1; --j) {
        nodes.push_back(Point{Real(-1), line_coord_pm_one(j, order), Real(0)});
    }
    for (int j = 1; j < order; ++j) {
        for (int i = 1; i < order; ++i) {
            nodes.push_back(Point{line_coord_pm_one(i, order),
                                  line_coord_pm_one(j, order), Real(0)});
        }
    }
    return nodes;
}

std::vector<Point> generate_tetra_nodes(int order) {
    if (order == 0) {
        return {Point{Real(0.25), Real(0.25), Real(0.25)}};
    }

    const Point verts[] = {
        Point{Real(0), Real(0), Real(0)},
        Point{Real(1), Real(0), Real(0)},
        Point{Real(0), Real(1), Real(0)},
        Point{Real(0), Real(0), Real(1)},
    };

    std::vector<Point> nodes;
    nodes.reserve(static_cast<std::size_t>((order + 1) * (order + 2) * (order + 3) / 6));
    for (const auto& v : verts) {
        nodes.push_back(v);
    }

    const int edges[6][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};
    for (const auto& edge : edges) {
        for (int m = 1; m < order; ++m) {
            const Real t = static_cast<Real>(m) / static_cast<Real>(order);
            nodes.push_back(verts[edge[0]] * (Real(1) - t) + verts[edge[1]] * t);
        }
    }

    const int faces[4][3] = {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}};
    for (const auto& face : faces) {
        append_triangle_face_interior(nodes,
                                      verts[face[0]],
                                      verts[face[1]],
                                      verts[face[2]],
                                      order);
    }

    for (int l = 1; l <= order - 3; ++l) {
        for (int k = 1; k <= order - l - 2; ++k) {
            for (int j = 1; j <= order - l - k - 1; ++j) {
                nodes.push_back(Point{Real(j) / Real(order),
                                      Real(k) / Real(order),
                                      Real(l) / Real(order)});
            }
        }
    }
    return nodes;
}

std::vector<Point> generate_hex_nodes(int order) {
    if (order == 0) {
        return {Point{Real(0), Real(0), Real(0)}};
    }

    const Point verts[] = {
        Point{Real(-1), Real(-1), Real(-1)},
        Point{Real(1), Real(-1), Real(-1)},
        Point{Real(1), Real(1), Real(-1)},
        Point{Real(-1), Real(1), Real(-1)},
        Point{Real(-1), Real(-1), Real(1)},
        Point{Real(1), Real(-1), Real(1)},
        Point{Real(1), Real(1), Real(1)},
        Point{Real(-1), Real(1), Real(1)},
    };

    std::vector<Point> nodes;
    nodes.reserve(static_cast<std::size_t>((order + 1) * (order + 1) * (order + 1)));
    for (const auto& v : verts) {
        nodes.push_back(v);
    }

    const int edges[12][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0},
        {4, 5}, {5, 6}, {6, 7}, {7, 4},
        {0, 4}, {1, 5}, {2, 6}, {3, 7},
    };
    for (const auto& edge : edges) {
        for (int m = 1; m < order; ++m) {
            const Real t = static_cast<Real>(m) / static_cast<Real>(order);
            nodes.push_back(verts[edge[0]] * (Real(1) - t) + verts[edge[1]] * t);
        }
    }

    // Face-interior nodes, emitted in VTK face order so the layout matches the
    // VTK cell node numbering the solver ingests from .vtu meshes:
    //   -X, +X, -Y, +Y, -Z, +Z  (e.g. Hex27 face centers become nodes 20..25).
    // For order >= 3 the within-face node sequence follows the loops below; only
    // the face order is normalized to VTK, which is all the supported Hex8/20/27
    // elements require.
    // -X face (x = -1)
    for (int k = 1; k < order; ++k) {
        for (int j = order - 1; j >= 1; --j) {
            nodes.push_back(Point{Real(-1), line_coord_pm_one(j, order), line_coord_pm_one(k, order)});
        }
    }
    // +X face (x = +1)
    for (int k = 1; k < order; ++k) {
        for (int j = 1; j < order; ++j) {
            nodes.push_back(Point{Real(1), line_coord_pm_one(j, order), line_coord_pm_one(k, order)});
        }
    }
    // -Y face (y = -1)
    for (int k = 1; k < order; ++k) {
        for (int i = 1; i < order; ++i) {
            nodes.push_back(Point{line_coord_pm_one(i, order), Real(-1), line_coord_pm_one(k, order)});
        }
    }
    // +Y face (y = +1)
    for (int k = 1; k < order; ++k) {
        for (int i = order - 1; i >= 1; --i) {
            nodes.push_back(Point{line_coord_pm_one(i, order), Real(1), line_coord_pm_one(k, order)});
        }
    }
    // -Z face (z = -1)
    for (int j = 1; j < order; ++j) {
        for (int i = 1; i < order; ++i) {
            nodes.push_back(Point{line_coord_pm_one(i, order), line_coord_pm_one(j, order), Real(-1)});
        }
    }
    // +Z face (z = +1)
    for (int j = 1; j < order; ++j) {
        for (int i = 1; i < order; ++i) {
            nodes.push_back(Point{line_coord_pm_one(i, order), line_coord_pm_one(j, order), Real(1)});
        }
    }
    for (int k = 1; k < order; ++k) {
        for (int j = 1; j < order; ++j) {
            for (int i = 1; i < order; ++i) {
                nodes.push_back(Point{line_coord_pm_one(i, order),
                                      line_coord_pm_one(j, order),
                                      line_coord_pm_one(k, order)});
            }
        }
    }
    return nodes;
}

std::vector<Point> generate_wedge_nodes(int order) {
    if (order == 0) {
        return {Point{Real(1) / Real(3), Real(1) / Real(3), Real(0)}};
    }

    const Point verts[] = {
        Point{Real(0), Real(0), Real(-1)},
        Point{Real(1), Real(0), Real(-1)},
        Point{Real(0), Real(1), Real(-1)},
        Point{Real(0), Real(0), Real(1)},
        Point{Real(1), Real(0), Real(1)},
        Point{Real(0), Real(1), Real(1)},
    };

    std::vector<Point> nodes;
    nodes.reserve(static_cast<std::size_t>((order + 1) * (order + 1) * (order + 2) / 2));
    for (const auto& v : verts) {
        nodes.push_back(v);
    }

    const int edges[9][2] = {
        {0, 1}, {1, 2}, {2, 0},
        {3, 4}, {4, 5}, {5, 3},
        {0, 3}, {1, 4}, {2, 5},
    };
    for (const auto& edge : edges) {
        for (int m = 1; m < order; ++m) {
            const Real t = static_cast<Real>(m) / static_cast<Real>(order);
            nodes.push_back(verts[edge[0]] * (Real(1) - t) + verts[edge[1]] * t);
        }
    }

    append_triangle_face_interior(nodes, verts[0], verts[1], verts[2], order);
    append_triangle_face_interior(nodes, verts[3], verts[4], verts[5], order);

    for (int r = 1; r < order; ++r) {
        const Real z = line_coord_pm_one(r, order);
        for (int m = 1; m < order; ++m) {
            const Real t = static_cast<Real>(m) / static_cast<Real>(order);
            nodes.push_back(Point{t, Real(0), z});
        }
        for (int m = 1; m < order; ++m) {
            const Real t = static_cast<Real>(m) / static_cast<Real>(order);
            nodes.push_back(Point{Real(1) - t, t, z});
        }
        for (int m = 1; m < order; ++m) {
            const Real t = static_cast<Real>(m) / static_cast<Real>(order);
            nodes.push_back(Point{Real(0), Real(1) - t, z});
        }
    }

    for (int r = 1; r < order; ++r) {
        const Real z = line_coord_pm_one(r, order);
        for (int c = 1; c <= order - 2; ++c) {
            for (int b = 1; b <= order - c - 1; ++b) {
                nodes.push_back(Point{Real(b) / Real(order),
                                      Real(c) / Real(order),
                                      z});
            }
        }
    }
    return nodes;
}

std::vector<Point> complete_lagrange_nodes(ElementType canonical_type, int order) {
    svmp::throw_if<BasisNodeOrderingException>(order < 0, SVMP_HERE,
                                             "ReferenceNodeLayout requires non-negative Lagrange order");
    const ElementType type = canonical_lagrange_type(canonical_type);
    switch (type) {
        case ElementType::Point1:
            return {Point{Real(0), Real(0), Real(0)}};
        case ElementType::Line2:
            return generate_line_nodes(order);
        case ElementType::Triangle3:
            return generate_triangle_nodes(order);
        case ElementType::Quad4:
            return generate_quad_nodes(order);
        case ElementType::Tetra4:
            return generate_tetra_nodes(order);
        case ElementType::Hex8:
            return generate_hex_nodes(order);
        case ElementType::Wedge6:
            return generate_wedge_nodes(order);
        case ElementType::Pyramid5:
            svmp::raise<BasisNodeOrderingException>(SVMP_HERE,
                "ReferenceNodeLayout: pyramid node ordering is disabled");
        default:
            svmp::raise<BasisNodeOrderingException>(SVMP_HERE,
                "ReferenceNodeLayout: unsupported Lagrange topology");
    }
}

std::vector<Point> element_nodes(ElementType elem_type) {
    const int order = complete_lagrange_alias_order(elem_type);
    if (order >= 0) {
        return complete_lagrange_nodes(elem_type, order);
    }

    switch (elem_type) {
        case ElementType::Quad8: {
            auto nodes = generate_quad_nodes(2);
            nodes.resize(8u);
            return nodes;
        }
        case ElementType::Hex20: {
            auto nodes = generate_hex_nodes(2);
            nodes.resize(20u);
            return nodes;
        }
        case ElementType::Wedge15: {
            auto nodes = generate_wedge_nodes(2);
            nodes.resize(15u);
            return nodes;
        }
        case ElementType::Pyramid13:
            svmp::raise<BasisNodeOrderingException>(SVMP_HERE,
                "ReferenceNodeLayout: pyramid node ordering is disabled");
        default:
            svmp::raise<BasisNodeOrderingException>(SVMP_HERE,
                "ReferenceNodeLayout: unknown element type");
    }
}

} // namespace

math::Vector<Real, 3> ReferenceNodeLayout::get_node_coords(ElementType elem_type,
                                                           std::size_t local_node) {
    const auto nodes = element_nodes(elem_type);
    svmp::throw_if<BasisNodeOrderingException>(local_node >= nodes.size(), SVMP_HERE,
                                             "ReferenceNodeLayout::get_node_coords: node index out of range");
    return nodes[local_node];
}

std::size_t ReferenceNodeLayout::num_nodes(ElementType elem_type) {
    return element_nodes(elem_type).size();
}

std::vector<math::Vector<Real, 3>>
ReferenceNodeLayout::get_lagrange_node_coords(ElementType canonical_type, int order) {
    return complete_lagrange_nodes(canonical_type, order);
}

std::span<const std::size_t> ReferenceNodeLayout::mesh_to_basis_ordering(ElementType elem_type) {
    if (elem_type == ElementType::Hex20) {
        return std::span<const std::size_t>(kHex20MeshToBasisOrder.data(),
                                            kHex20MeshToBasisOrder.size());
    }
    return {};
}

bool ReferenceNodeLayout::is_simplex(ElementType elem_type) {
    return svmp::FE::basis::is_simplex(elem_type);
}

bool ReferenceNodeLayout::is_tensor_product(ElementType elem_type) {
    return svmp::FE::basis::is_tensor_product(elem_type);
}

} // namespace basis
} // namespace FE
} // namespace svmp
