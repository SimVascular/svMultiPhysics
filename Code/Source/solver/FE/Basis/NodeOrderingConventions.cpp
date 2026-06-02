/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#include "NodeOrderingConventions.h"
#include "Basis/BasisExceptions.h"
#include "Basis/BasisTraits.h"

#include <array>

namespace svmp {
namespace FE {
namespace basis {

namespace {

using Point = math::Vector<Real, 3>;
using RawPoint = std::array<Real, 3>;

template<std::size_t N>
using NodeTable = std::array<RawPoint, N>;

struct NodeTableView {
    const RawPoint* data{nullptr};
    std::size_t size{0};
};

inline constexpr NodeTable<2> kLine2Nodes = {{
    {Real(-1), Real(0), Real(0)},
    {Real(1), Real(0), Real(0)},
}};

inline constexpr NodeTable<3> kLine3Nodes = {{
    {Real(-1), Real(0), Real(0)},
    {Real(1), Real(0), Real(0)},
    {Real(0), Real(0), Real(0)},
}};

inline constexpr NodeTable<3> kTriangle3Nodes = {{
    {Real(0), Real(0), Real(0)},
    {Real(1), Real(0), Real(0)},
    {Real(0), Real(1), Real(0)},
}};

inline constexpr NodeTable<6> kTriangle6Nodes = {{
    {Real(0), Real(0), Real(0)},
    {Real(1), Real(0), Real(0)},
    {Real(0), Real(1), Real(0)},
    {Real(0.5), Real(0), Real(0)},
    {Real(0.5), Real(0.5), Real(0)},
    {Real(0), Real(0.5), Real(0)},
}};

inline constexpr NodeTable<4> kQuad4Nodes = {{
    {Real(-1), Real(-1), Real(0)},
    {Real(1), Real(-1), Real(0)},
    {Real(1), Real(1), Real(0)},
    {Real(-1), Real(1), Real(0)},
}};

inline constexpr NodeTable<9> kQuad9Nodes = {{
    {Real(-1), Real(-1), Real(0)},
    {Real(1), Real(-1), Real(0)},
    {Real(1), Real(1), Real(0)},
    {Real(-1), Real(1), Real(0)},
    {Real(0), Real(-1), Real(0)},
    {Real(1), Real(0), Real(0)},
    {Real(0), Real(1), Real(0)},
    {Real(-1), Real(0), Real(0)},
    {Real(0), Real(0), Real(0)},
}};

inline constexpr NodeTable<8> kQuad8Nodes = {{
    {Real(-1), Real(-1), Real(0)},
    {Real(1), Real(-1), Real(0)},
    {Real(1), Real(1), Real(0)},
    {Real(-1), Real(1), Real(0)},
    {Real(0), Real(-1), Real(0)},
    {Real(1), Real(0), Real(0)},
    {Real(0), Real(1), Real(0)},
    {Real(-1), Real(0), Real(0)},
}};

inline constexpr NodeTable<4> kTetra4Nodes = {{
    {Real(0), Real(0), Real(0)},
    {Real(1), Real(0), Real(0)},
    {Real(0), Real(1), Real(0)},
    {Real(0), Real(0), Real(1)},
}};

inline constexpr NodeTable<10> kTetra10Nodes = {{
    {Real(0), Real(0), Real(0)},
    {Real(1), Real(0), Real(0)},
    {Real(0), Real(1), Real(0)},
    {Real(0), Real(0), Real(1)},
    {Real(0.5), Real(0), Real(0)},
    {Real(0.5), Real(0.5), Real(0)},
    {Real(0), Real(0.5), Real(0)},
    {Real(0), Real(0), Real(0.5)},
    {Real(0.5), Real(0), Real(0.5)},
    {Real(0), Real(0.5), Real(0.5)},
}};

inline constexpr NodeTable<8> kHex8Nodes = {{
    {Real(-1), Real(-1), Real(-1)},
    {Real(1), Real(-1), Real(-1)},
    {Real(1), Real(1), Real(-1)},
    {Real(-1), Real(1), Real(-1)},
    {Real(-1), Real(-1), Real(1)},
    {Real(1), Real(-1), Real(1)},
    {Real(1), Real(1), Real(1)},
    {Real(-1), Real(1), Real(1)},
}};

inline constexpr NodeTable<27> kHex27Nodes = {{
    {Real(-1), Real(-1), Real(-1)},
    {Real(1), Real(-1), Real(-1)},
    {Real(1), Real(1), Real(-1)},
    {Real(-1), Real(1), Real(-1)},
    {Real(-1), Real(-1), Real(1)},
    {Real(1), Real(-1), Real(1)},
    {Real(1), Real(1), Real(1)},
    {Real(-1), Real(1), Real(1)},
    {Real(0), Real(-1), Real(-1)},
    {Real(1), Real(0), Real(-1)},
    {Real(0), Real(1), Real(-1)},
    {Real(-1), Real(0), Real(-1)},
    {Real(0), Real(-1), Real(1)},
    {Real(1), Real(0), Real(1)},
    {Real(0), Real(1), Real(1)},
    {Real(-1), Real(0), Real(1)},
    {Real(-1), Real(-1), Real(0)},
    {Real(1), Real(-1), Real(0)},
    {Real(1), Real(1), Real(0)},
    {Real(-1), Real(1), Real(0)},
    {Real(0), Real(0), Real(-1)},
    {Real(0), Real(0), Real(1)},
    {Real(0), Real(-1), Real(0)},
    {Real(1), Real(0), Real(0)},
    {Real(0), Real(1), Real(0)},
    {Real(-1), Real(0), Real(0)},
    {Real(0), Real(0), Real(0)},
}};

inline constexpr NodeTable<20> kHex20Nodes = {{
    {Real(-1), Real(-1), Real(-1)},
    {Real(1), Real(-1), Real(-1)},
    {Real(1), Real(1), Real(-1)},
    {Real(-1), Real(1), Real(-1)},
    {Real(-1), Real(-1), Real(1)},
    {Real(1), Real(-1), Real(1)},
    {Real(1), Real(1), Real(1)},
    {Real(-1), Real(1), Real(1)},
    {Real(0), Real(-1), Real(-1)},
    {Real(1), Real(0), Real(-1)},
    {Real(0), Real(1), Real(-1)},
    {Real(-1), Real(0), Real(-1)},
    {Real(0), Real(-1), Real(1)},
    {Real(1), Real(0), Real(1)},
    {Real(0), Real(1), Real(1)},
    {Real(-1), Real(0), Real(1)},
    {Real(-1), Real(-1), Real(0)},
    {Real(1), Real(-1), Real(0)},
    {Real(1), Real(1), Real(0)},
    {Real(-1), Real(1), Real(0)},
}};

// Mesh uses conventional Hex20 ordering: corners first, then edge midpoints in
// {bottom, top, vertical} groups. The quadratic Hex20 serendipity polynomial
// table uses an axis-grouped edge order. This maps public mesh/reference index
// to the internal polynomial-table index.
constexpr std::array<std::size_t, 20> kHex20MeshToBasisOrder = {
    0, 1, 2, 3, 4, 5, 6, 7,
    8, 13, 10, 12,
    9, 15, 11, 14,
    16, 17, 19, 18
};

inline constexpr NodeTable<6> kWedge6Nodes = {{
    {Real(0), Real(0), Real(-1)},
    {Real(1), Real(0), Real(-1)},
    {Real(0), Real(1), Real(-1)},
    {Real(0), Real(0), Real(1)},
    {Real(1), Real(0), Real(1)},
    {Real(0), Real(1), Real(1)},
}};

inline constexpr NodeTable<18> kWedge18Nodes = {{
    {Real(0), Real(0), Real(-1)},
    {Real(1), Real(0), Real(-1)},
    {Real(0), Real(1), Real(-1)},
    {Real(0), Real(0), Real(1)},
    {Real(1), Real(0), Real(1)},
    {Real(0), Real(1), Real(1)},
    {Real(0.5), Real(0), Real(-1)},
    {Real(0.5), Real(0.5), Real(-1)},
    {Real(0), Real(0.5), Real(-1)},
    {Real(0.5), Real(0), Real(1)},
    {Real(0.5), Real(0.5), Real(1)},
    {Real(0), Real(0.5), Real(1)},
    {Real(0), Real(0), Real(0)},
    {Real(1), Real(0), Real(0)},
    {Real(0), Real(1), Real(0)},
    {Real(0.5), Real(0), Real(0)},
    {Real(0.5), Real(0.5), Real(0)},
    {Real(0), Real(0.5), Real(0)},
}};

inline constexpr NodeTable<15> kWedge15Nodes = {{
    {Real(0), Real(0), Real(-1)},
    {Real(1), Real(0), Real(-1)},
    {Real(0), Real(1), Real(-1)},
    {Real(0), Real(0), Real(1)},
    {Real(1), Real(0), Real(1)},
    {Real(0), Real(1), Real(1)},
    {Real(0.5), Real(0), Real(-1)},
    {Real(0.5), Real(0.5), Real(-1)},
    {Real(0), Real(0.5), Real(-1)},
    {Real(0.5), Real(0), Real(1)},
    {Real(0.5), Real(0.5), Real(1)},
    {Real(0), Real(0.5), Real(1)},
    {Real(0), Real(0), Real(0)},
    {Real(1), Real(0), Real(0)},
    {Real(0), Real(1), Real(0)},
}};

inline constexpr NodeTable<5> kPyramid5Nodes = {{
    {Real(-1), Real(-1), Real(0)},
    {Real(1), Real(-1), Real(0)},
    {Real(1), Real(1), Real(0)},
    {Real(-1), Real(1), Real(0)},
    {Real(0), Real(0), Real(1)},
}};

inline constexpr NodeTable<14> kPyramid14Nodes = {{
    {Real(-1), Real(-1), Real(0)},
    {Real(1), Real(-1), Real(0)},
    {Real(1), Real(1), Real(0)},
    {Real(-1), Real(1), Real(0)},
    {Real(0), Real(0), Real(1)},
    {Real(0), Real(-1), Real(0)},
    {Real(1), Real(0), Real(0)},
    {Real(0), Real(1), Real(0)},
    {Real(-1), Real(0), Real(0)},
    {Real(-0.5), Real(-0.5), Real(0.5)},
    {Real(0.5), Real(-0.5), Real(0.5)},
    {Real(0.5), Real(0.5), Real(0.5)},
    {Real(-0.5), Real(0.5), Real(0.5)},
    {Real(0), Real(0), Real(0)},
}};

inline constexpr NodeTable<13> kPyramid13Nodes = {{
    {Real(-1), Real(-1), Real(0)},
    {Real(1), Real(-1), Real(0)},
    {Real(1), Real(1), Real(0)},
    {Real(-1), Real(1), Real(0)},
    {Real(0), Real(0), Real(1)},
    {Real(0), Real(-1), Real(0)},
    {Real(1), Real(0), Real(0)},
    {Real(0), Real(1), Real(0)},
    {Real(-1), Real(0), Real(0)},
    {Real(-0.5), Real(-0.5), Real(0.5)},
    {Real(0.5), Real(-0.5), Real(0.5)},
    {Real(0.5), Real(0.5), Real(0.5)},
    {Real(-0.5), Real(0.5), Real(0.5)},
}};

template<std::size_t N>
constexpr NodeTableView view(const NodeTable<N>& table) noexcept {
    return NodeTableView{table.data(), table.size()};
}

Point to_point(const RawPoint& raw) {
    return Point{raw[0], raw[1], raw[2]};
}

constexpr NodeTableView fixed_node_table(ElementType elem_type) noexcept {
    switch (elem_type) {
        case ElementType::Line2:     return view(kLine2Nodes);
        case ElementType::Line3:     return view(kLine3Nodes);
        case ElementType::Triangle3: return view(kTriangle3Nodes);
        case ElementType::Triangle6: return view(kTriangle6Nodes);
        case ElementType::Quad4:     return view(kQuad4Nodes);
        case ElementType::Quad8:     return view(kQuad8Nodes);
        case ElementType::Quad9:     return view(kQuad9Nodes);
        case ElementType::Tetra4:    return view(kTetra4Nodes);
        case ElementType::Tetra10:   return view(kTetra10Nodes);
        case ElementType::Hex8:      return view(kHex8Nodes);
        case ElementType::Hex20:     return view(kHex20Nodes);
        case ElementType::Hex27:     return view(kHex27Nodes);
        case ElementType::Wedge6:    return view(kWedge6Nodes);
        case ElementType::Wedge15:   return view(kWedge15Nodes);
        case ElementType::Wedge18:   return view(kWedge18Nodes);
        case ElementType::Pyramid5:  return view(kPyramid5Nodes);
        case ElementType::Pyramid13: return view(kPyramid13Nodes);
        case ElementType::Pyramid14: return view(kPyramid14Nodes);
        default:                     return {};
    }
}

constexpr NodeTableView fixed_complete_lagrange_table(ElementType canonical_type,
                                                      int order) noexcept {
    switch (canonical_type) {
        case ElementType::Line2:
            return order == 1 ? view(kLine2Nodes) :
                   order == 2 ? view(kLine3Nodes) : NodeTableView{};
        case ElementType::Triangle3:
            return order == 1 ? view(kTriangle3Nodes) :
                   order == 2 ? view(kTriangle6Nodes) : NodeTableView{};
        case ElementType::Quad4:
            return order == 1 ? view(kQuad4Nodes) :
                   order == 2 ? view(kQuad9Nodes) : NodeTableView{};
        case ElementType::Tetra4:
            return order == 1 ? view(kTetra4Nodes) :
                   order == 2 ? view(kTetra10Nodes) : NodeTableView{};
        case ElementType::Hex8:
            return order == 1 ? view(kHex8Nodes) :
                   order == 2 ? view(kHex27Nodes) : NodeTableView{};
        case ElementType::Wedge6:
            return order == 1 ? view(kWedge6Nodes) :
                   order == 2 ? view(kWedge18Nodes) : NodeTableView{};
        case ElementType::Pyramid5:
            return order == 1 ? view(kPyramid5Nodes) :
                   order == 2 ? view(kPyramid14Nodes) : NodeTableView{};
        default:
            return {};
    }
}

Real line_coord_pm_one(int i, int order) {
    if (order <= 0) {
        return Real(0);
    }
    return Real(-1) + Real(2) * static_cast<Real>(i) / static_cast<Real>(order);
}

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
            const Real la = static_cast<Real>(a) / static_cast<Real>(order);
            const Real lb = static_cast<Real>(b) / static_cast<Real>(order);
            const Real lc = static_cast<Real>(c) / static_cast<Real>(order);
            nodes.push_back(v0 * la + v1 * lb + v2 * lc);
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

    append_triangle_face_interior(
        nodes,
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
            nodes.push_back(Point{line_coord_pm_one(i, order), line_coord_pm_one(j, order), Real(0)});
        }
    }

    return nodes;
}

std::vector<Point> generate_tetra_nodes(int order) {
    if (order == 0) {
        return {Point{Real(0.25), Real(0.25), Real(0.25)}};
    }

    std::vector<Point> nodes;
    nodes.reserve(static_cast<std::size_t>((order + 1) * (order + 2) * (order + 3) / 6));

    const Point verts[] = {
        Point{Real(0), Real(0), Real(0)},
        Point{Real(1), Real(0), Real(0)},
        Point{Real(0), Real(1), Real(0)},
        Point{Real(0), Real(0), Real(1)},
    };
    for (const auto& v : verts) {
        nodes.push_back(v);
    }

    const int edges[6][2] = {
        {0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}
    };
    for (const auto& edge : edges) {
        for (int m = 1; m < order; ++m) {
            const Real t = static_cast<Real>(m) / static_cast<Real>(order);
            nodes.push_back(verts[edge[0]] * (Real(1) - t) + verts[edge[1]] * t);
        }
    }

    const int faces[4][3] = {
        {0, 1, 2},
        {0, 1, 3},
        {1, 2, 3},
        {0, 2, 3},
    };
    for (const auto& face : faces) {
        append_triangle_face_interior(
            nodes,
            verts[face[0]],
            verts[face[1]],
            verts[face[2]],
            order);
    }

    for (int l = 1; l <= order - 3; ++l) {
        for (int k = 1; k <= order - l - 2; ++k) {
            for (int j = 1; j <= order - l - k - 1; ++j) {
                const Real x = static_cast<Real>(j) / static_cast<Real>(order);
                const Real y = static_cast<Real>(k) / static_cast<Real>(order);
                const Real z = static_cast<Real>(l) / static_cast<Real>(order);
                nodes.push_back(Point{x, y, z});
            }
        }
    }

    return nodes;
}

std::vector<Point> generate_hex_nodes(int order) {
    if (order == 0) {
        return {Point{Real(0), Real(0), Real(0)}};
    }

    std::vector<Point> nodes;
    nodes.reserve(static_cast<std::size_t>((order + 1) * (order + 1) * (order + 1)));

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

    for (int j = 1; j < order; ++j) {
        for (int i = 1; i < order; ++i) {
            nodes.push_back(Point{line_coord_pm_one(i, order), line_coord_pm_one(j, order), Real(-1)});
        }
    }
    for (int j = 1; j < order; ++j) {
        for (int i = 1; i < order; ++i) {
            nodes.push_back(Point{line_coord_pm_one(i, order), line_coord_pm_one(j, order), Real(1)});
        }
    }
    for (int k = 1; k < order; ++k) {
        for (int i = 1; i < order; ++i) {
            nodes.push_back(Point{line_coord_pm_one(i, order), Real(-1), line_coord_pm_one(k, order)});
        }
    }
    for (int k = 1; k < order; ++k) {
        for (int j = 1; j < order; ++j) {
            nodes.push_back(Point{Real(1), line_coord_pm_one(j, order), line_coord_pm_one(k, order)});
        }
    }
    for (int k = 1; k < order; ++k) {
        for (int i = order - 1; i >= 1; --i) {
            nodes.push_back(Point{line_coord_pm_one(i, order), Real(1), line_coord_pm_one(k, order)});
        }
    }
    for (int k = 1; k < order; ++k) {
        for (int j = order - 1; j >= 1; --j) {
            nodes.push_back(Point{Real(-1), line_coord_pm_one(j, order), line_coord_pm_one(k, order)});
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

    std::vector<Point> nodes;
    nodes.reserve(static_cast<std::size_t>((order + 1) * (order + 1) * (order + 2) / 2));

    const Point verts[] = {
        Point{Real(0), Real(0), Real(-1)},
        Point{Real(1), Real(0), Real(-1)},
        Point{Real(0), Real(1), Real(-1)},
        Point{Real(0), Real(0), Real(1)},
        Point{Real(1), Real(0), Real(1)},
        Point{Real(0), Real(1), Real(1)},
    };
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

    append_triangle_face_interior(
        nodes, verts[0], verts[1], verts[2], order);
    append_triangle_face_interior(
        nodes, verts[3], verts[4], verts[5], order);

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
                const Real x = static_cast<Real>(b) / static_cast<Real>(order);
                const Real y = static_cast<Real>(c) / static_cast<Real>(order);
                nodes.push_back(Point{x, y, z});
            }
        }
    }

    return nodes;
}

std::vector<Point> generate_pyramid_nodes(int order) {
    if (order == 0) {
        return {Point{Real(0), Real(0), Real(0.25)}};
    }

    std::vector<Point> nodes;
    nodes.reserve(static_cast<std::size_t>((order + 1) * (order + 2) * (2 * order + 3) / 6));

    nodes.push_back(Point{Real(-1), Real(-1), Real(0)});
    nodes.push_back(Point{Real(1), Real(-1), Real(0)});
    nodes.push_back(Point{Real(1), Real(1), Real(0)});
    nodes.push_back(Point{Real(-1), Real(1), Real(0)});
    nodes.push_back(Point{Real(0), Real(0), Real(1)});

    for (int m = 1; m < order; ++m) {
        nodes.push_back(Point{line_coord_pm_one(m, order), Real(-1), Real(0)});
    }
    for (int m = 1; m < order; ++m) {
        nodes.push_back(Point{Real(1), line_coord_pm_one(m, order), Real(0)});
    }
    for (int m = order - 1; m >= 1; --m) {
        nodes.push_back(Point{line_coord_pm_one(m, order), Real(1), Real(0)});
    }
    for (int m = order - 1; m >= 1; --m) {
        nodes.push_back(Point{Real(-1), line_coord_pm_one(m, order), Real(0)});
    }

    for (int level = 1; level < order; ++level) {
        const Real z = static_cast<Real>(level) / static_cast<Real>(order);
        const Real scale = Real(1) - z;
        nodes.push_back(Point{-scale, -scale, z});
        nodes.push_back(Point{scale, -scale, z});
        nodes.push_back(Point{scale, scale, z});
        nodes.push_back(Point{-scale, scale, z});
    }

    for (int j = 1; j < order; ++j) {
        for (int i = 1; i < order; ++i) {
            nodes.push_back(Point{line_coord_pm_one(i, order), line_coord_pm_one(j, order), Real(0)});
        }
    }

    for (int level = 1; level < order - 1; ++level) {
        const int n = order - level;
        const Real z = static_cast<Real>(level) / static_cast<Real>(order);
        const Real scale = Real(1) - z;

        for (int m = 1; m < n; ++m) {
            const Real s = line_coord_pm_one(m, n) * scale;
            nodes.push_back(Point{s, -scale, z});
        }
        for (int m = 1; m < n; ++m) {
            const Real s = line_coord_pm_one(m, n) * scale;
            nodes.push_back(Point{scale, s, z});
        }
        for (int m = n - 1; m >= 1; --m) {
            const Real s = line_coord_pm_one(m, n) * scale;
            nodes.push_back(Point{s, scale, z});
        }
        for (int m = n - 1; m >= 1; --m) {
            const Real s = line_coord_pm_one(m, n) * scale;
            nodes.push_back(Point{-scale, s, z});
        }
    }

    for (int level = 1; level < order - 1; ++level) {
        const int n = order - level;
        const Real z = static_cast<Real>(level) / static_cast<Real>(order);
        const Real scale = Real(1) - z;
        for (int j = 1; j < n; ++j) {
            for (int i = 1; i < n; ++i) {
                nodes.push_back(Point{line_coord_pm_one(i, n) * scale,
                                      line_coord_pm_one(j, n) * scale,
                                      z});
            }
        }
    }

    return nodes;
}

} // namespace

math::Vector<Real, 3> ReferenceNodeLayout::get_node_coords(ElementType elem_type,
                                                     std::size_t local_node) {
    const auto table = fixed_node_table(elem_type);
    if (table.data != nullptr && local_node < table.size) {
        return to_point(table.data[local_node]);
    }

    throw BasisNodeOrderingException("Invalid element type or node index in ReferenceNodeLayout::get_node_coords",
                                     __FILE__, __LINE__, __func__);
}

std::size_t ReferenceNodeLayout::num_nodes(ElementType elem_type) {
    const auto table = fixed_node_table(elem_type);
    if (table.data != nullptr) {
        return table.size;
    }

    throw BasisNodeOrderingException("Unknown element type in ReferenceNodeLayout::num_nodes",
                                     __FILE__, __LINE__, __func__);
}

std::vector<math::Vector<Real, 3>>
ReferenceNodeLayout::get_lagrange_node_coords(ElementType canonical_type, int order) {
    if (order < 0) {
        throw BasisNodeOrderingException("ReferenceNodeLayout::get_lagrange_node_coords requires non-negative order",
                                         __FILE__, __LINE__, __func__);
    }

    const ElementType type = canonical_lagrange_type(canonical_type);
    const auto fixed_table = fixed_complete_lagrange_table(type, order);
    if (fixed_table.data != nullptr) {
        std::vector<Point> nodes;
        nodes.reserve(fixed_table.size);
        for (std::size_t i = 0; i < fixed_table.size; ++i) {
            nodes.push_back(to_point(fixed_table.data[i]));
        }
        return nodes;
    }

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
            return generate_pyramid_nodes(order);
        case ElementType::Quad8:
        case ElementType::Hex20:
        case ElementType::Wedge15:
        case ElementType::Pyramid13:
            throw BasisNodeOrderingException("ReferenceNodeLayout::get_lagrange_node_coords does not support serendipity topologies",
                                             __FILE__, __LINE__, __func__);
        default:
            throw BasisNodeOrderingException("ReferenceNodeLayout::get_lagrange_node_coords: unsupported topology",
                                             __FILE__, __LINE__, __func__);
    }
}

std::span<const std::size_t> ReferenceNodeLayout::mesh_to_basis_ordering(ElementType elem_type) {
    if (elem_type == ElementType::Hex20) {
        return std::span<const std::size_t>(
            kHex20MeshToBasisOrder.data(),
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
