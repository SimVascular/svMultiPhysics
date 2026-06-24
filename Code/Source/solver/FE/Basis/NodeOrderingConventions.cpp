// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "NodeOrderingConventions.h"
#include "BasisExceptions.h"
#include "BasisTraits.h"

#include <array>
#include <cmath>
#include <map>
#include <utility>

namespace svmp {
namespace FE {
namespace basis {

namespace {

using Point = math::Vector<double, 3>;
using Lattice = std::array<int, 3>;

// Gauss-Lobatto-Legendre nodes on [-1, 1] for a degree-`order` distribution
// (order + 1 nodes). The endpoints are -1 and +1; the interior nodes are the
// roots of P'_order, found by Newton iteration on f(x) = x P_order(x) -
// P_{order-1}(x) -- whose roots are exactly the GLL nodes -- from the
// Chebyshev-Gauss-Lobatto seed.
const std::vector<double>& gll_points(int order) {
    thread_local std::map<int, std::vector<double>> cache;
    const auto found = cache.find(order);
    if (found != cache.end()) {
        return found->second;
    }

    std::vector<double> pts(static_cast<std::size_t>(order + 1), double(0));
    if (order >= 1) {
        pts.front() = double(-1);
        pts.back() = double(1);
    }
    const double pi = std::acos(double(-1));
    const int half = order / 2;
    for (int j = 1; j <= half; ++j) {
        if (2 * j == order) {
            pts[static_cast<std::size_t>(j)] = double(0);  // exact center, even order
            continue;
        }
        double x = -std::cos(pi * static_cast<double>(j) / static_cast<double>(order));
        for (int iter = 0; iter < 100; ++iter) {
            // Legendre P_k and P'_k up to k = order at x, by the three-term
            // recurrences (regular at x = +/-1).
            double p_km1 = double(1);   // P_0
            double p_k = x;             // P_1
            double d_km1 = double(0);   // P'_0
            double d_k = double(1);     // P'_1
            for (int k = 1; k < order; ++k) {
                const double kk = static_cast<double>(k);
                const double inv = double(1) / (kk + double(1));
                const double p_kp1 =
                    ((double(2) * kk + double(1)) * x * p_k - kk * p_km1) * inv;
                const double d_kp1 =
                    ((double(2) * kk + double(1)) * (p_k + x * d_k) - kk * d_km1) * inv;
                p_km1 = p_k;
                p_k = p_kp1;
                d_km1 = d_k;
                d_k = d_kp1;
            }
            // p_k = P_order, p_km1 = P_{order-1}, d_k = P'_order, d_km1 = P'_{order-1}.
            const double f = x * p_k - p_km1;
            const double f_prime = p_k + x * d_k - d_km1;
            const double dx = f / f_prime;
            x -= dx;
            if (std::abs(dx) <= double(1e-15)) {
                break;
            }
        }
        pts[static_cast<std::size_t>(j)] = x;
    }
    for (int j = half + 1; j < order; ++j) {
        pts[static_cast<std::size_t>(j)] = -pts[static_cast<std::size_t>(order - j)];
    }

    auto inserted = cache.emplace(order, std::move(pts));
    return inserted.first->second;
}

double line_coord_zero_one(int i, int order) {
    if (order <= 0) {
        return double(0);
    }
    return static_cast<double>(i) / static_cast<double>(order);
}

// Interpolate an integer lattice index along an edge between two corner
// vertices: index = (LA * (order - m) + LB * m) / order. The division is exact
// because edge endpoints are element corners (each component is 0 or order), so
// the result is the integer lattice point at parameter m / order.
Lattice lerp_lattice(const Lattice& a, const Lattice& b, int m, int order) {
    Lattice result{0, 0, 0};
    for (std::size_t d = 0; d < 3u; ++d) {
        const int numerator = a[d] * (order - m) + b[d] * m;
        svmp::throw_if<BasisConstructionException>(
            numerator % order != 0, SVMP_HERE,
            "ReferenceNodeLayout: non-integral edge lattice index");
        result[d] = numerator / order;
    }
    return result;
}

// Barycentric combination of three corner lattice indices for a triangular
// face-interior node: index = (a * L0 + b * L1 + c * L2) / order, with
// a + b + c == order. Exact for corner inputs (components 0 or order).
Lattice combine_lattice(const Lattice& l0, const Lattice& l1, const Lattice& l2,
                        int a, int b, int c, int order) {
    Lattice result{0, 0, 0};
    for (std::size_t d = 0; d < 3u; ++d) {
        const int numerator = a * l0[d] + b * l1[d] + c * l2[d];
        svmp::throw_if<BasisConstructionException>(
            numerator % order != 0, SVMP_HERE,
            "ReferenceNodeLayout: non-integral face-interior lattice index");
        result[d] = numerator / order;
    }
    return result;
}

// Append the interior nodes of a triangular face spanned by v0, v1, v2 (with
// matching corner lattice indices l0, l1, l2), emitting both the coordinate and
// its integer lattice index. Shared by triangle interiors, tetra faces, and the
// two wedge caps.
void append_triangle_face_interior(LagrangeNodeLayout& out,
                                   const Point& v0,
                                   const Point& v1,
                                   const Point& v2,
                                   const Lattice& l0,
                                   const Lattice& l1,
                                   const Lattice& l2,
                                   int order) {
    for (int c = 1; c <= order - 2; ++c) {
        for (int b = 1; b <= order - c - 1; ++b) {
            const int a = order - b - c;
            const double inv = double(1) / double(order);
            out.coords.push_back(v0 * (double(a) * inv) +
                                 v1 * (double(b) * inv) +
                                 v2 * (double(c) * inv));
            out.lattice.push_back(combine_lattice(l0, l1, l2, a, b, c, order));
        }
    }
}

LagrangeNodeLayout generate_line_nodes(int order) {
    LagrangeNodeLayout out;
    if (order == 0) {
        out.coords.push_back(Point{double(0), double(0), double(0)});
        out.lattice.push_back(Lattice{0, 0, 0});
        return out;
    }

    out.coords.reserve(static_cast<std::size_t>(order + 1));
    out.lattice.reserve(static_cast<std::size_t>(order + 1));
    out.coords.push_back(Point{double(-1), double(0), double(0)});
    out.lattice.push_back(Lattice{0, 0, 0});
    out.coords.push_back(Point{double(1), double(0), double(0)});
    out.lattice.push_back(Lattice{order, 0, 0});
    for (int i = 1; i < order; ++i) {
        out.coords.push_back(Point{line_coord_pm_one(i, order), double(0), double(0)});
        out.lattice.push_back(Lattice{i, 0, 0});
    }
    return out;
}

LagrangeNodeLayout generate_triangle_nodes(int order) {
    LagrangeNodeLayout out;
    if (order == 0) {
        out.coords.push_back(Point{double(1) / double(3), double(1) / double(3), double(0)});
        out.lattice.push_back(Lattice{0, 0, 0});
        return out;
    }

    out.coords.reserve(static_cast<std::size_t>((order + 1) * (order + 2) / 2));
    out.lattice.reserve(static_cast<std::size_t>((order + 1) * (order + 2) / 2));
    out.coords.push_back(Point{double(0), double(0), double(0)});
    out.lattice.push_back(Lattice{0, 0, 0});
    out.coords.push_back(Point{double(1), double(0), double(0)});
    out.lattice.push_back(Lattice{order, 0, 0});
    out.coords.push_back(Point{double(0), double(1), double(0)});
    out.lattice.push_back(Lattice{0, order, 0});

    for (int m = 1; m < order; ++m) {
        out.coords.push_back(Point{line_coord_zero_one(m, order), double(0), double(0)});
        out.lattice.push_back(Lattice{m, 0, 0});
    }
    for (int m = 1; m < order; ++m) {
        out.coords.push_back(Point{line_coord_zero_one(order - m, order),
                                   line_coord_zero_one(m, order), double(0)});
        out.lattice.push_back(Lattice{order - m, m, 0});
    }
    for (int m = 1; m < order; ++m) {
        out.coords.push_back(Point{double(0), line_coord_zero_one(order - m, order), double(0)});
        out.lattice.push_back(Lattice{0, order - m, 0});
    }

    append_triangle_face_interior(out,
                                  Point{double(0), double(0), double(0)},
                                  Point{double(1), double(0), double(0)},
                                  Point{double(0), double(1), double(0)},
                                  Lattice{0, 0, 0},
                                  Lattice{order, 0, 0},
                                  Lattice{0, order, 0},
                                  order);
    return out;
}

LagrangeNodeLayout generate_quad_nodes(int order) {
    LagrangeNodeLayout out;
    if (order == 0) {
        out.coords.push_back(Point{double(0), double(0), double(0)});
        out.lattice.push_back(Lattice{0, 0, 0});
        return out;
    }

    out.coords.reserve(static_cast<std::size_t>((order + 1) * (order + 1)));
    out.lattice.reserve(static_cast<std::size_t>((order + 1) * (order + 1)));
    out.coords.push_back(Point{double(-1), double(-1), double(0)});
    out.lattice.push_back(Lattice{0, 0, 0});
    out.coords.push_back(Point{double(1), double(-1), double(0)});
    out.lattice.push_back(Lattice{order, 0, 0});
    out.coords.push_back(Point{double(1), double(1), double(0)});
    out.lattice.push_back(Lattice{order, order, 0});
    out.coords.push_back(Point{double(-1), double(1), double(0)});
    out.lattice.push_back(Lattice{0, order, 0});

    for (int i = 1; i < order; ++i) {
        out.coords.push_back(Point{line_coord_pm_one(i, order), double(-1), double(0)});
        out.lattice.push_back(Lattice{i, 0, 0});
    }
    for (int j = 1; j < order; ++j) {
        out.coords.push_back(Point{double(1), line_coord_pm_one(j, order), double(0)});
        out.lattice.push_back(Lattice{order, j, 0});
    }
    for (int i = order - 1; i >= 1; --i) {
        out.coords.push_back(Point{line_coord_pm_one(i, order), double(1), double(0)});
        out.lattice.push_back(Lattice{i, order, 0});
    }
    for (int j = order - 1; j >= 1; --j) {
        out.coords.push_back(Point{double(-1), line_coord_pm_one(j, order), double(0)});
        out.lattice.push_back(Lattice{0, j, 0});
    }
    for (int j = 1; j < order; ++j) {
        for (int i = 1; i < order; ++i) {
            out.coords.push_back(Point{line_coord_pm_one(i, order),
                                       line_coord_pm_one(j, order), double(0)});
            out.lattice.push_back(Lattice{i, j, 0});
        }
    }
    return out;
}

LagrangeNodeLayout generate_tetra_nodes(int order) {
    LagrangeNodeLayout out;
    if (order == 0) {
        out.coords.push_back(Point{double(0.25), double(0.25), double(0.25)});
        out.lattice.push_back(Lattice{0, 0, 0});
        return out;
    }

    const Point verts[] = {
        Point{double(0), double(0), double(0)},
        Point{double(1), double(0), double(0)},
        Point{double(0), double(1), double(0)},
        Point{double(0), double(0), double(1)},
    };
    const Lattice vert_lattice[] = {
        Lattice{0, 0, 0},
        Lattice{order, 0, 0},
        Lattice{0, order, 0},
        Lattice{0, 0, order},
    };

    out.coords.reserve(static_cast<std::size_t>((order + 1) * (order + 2) * (order + 3) / 6));
    out.lattice.reserve(static_cast<std::size_t>((order + 1) * (order + 2) * (order + 3) / 6));
    for (std::size_t v = 0; v < 4u; ++v) {
        out.coords.push_back(verts[v]);
        out.lattice.push_back(vert_lattice[v]);
    }

    const int edges[6][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};
    for (const auto& edge : edges) {
        for (int m = 1; m < order; ++m) {
            const double t = static_cast<double>(m) / static_cast<double>(order);
            out.coords.push_back(verts[edge[0]] * (double(1) - t) + verts[edge[1]] * t);
            out.lattice.push_back(lerp_lattice(vert_lattice[edge[0]], vert_lattice[edge[1]], m, order));
        }
    }

    const int faces[4][3] = {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}};
    for (const auto& face : faces) {
        append_triangle_face_interior(out,
                                      verts[face[0]], verts[face[1]], verts[face[2]],
                                      vert_lattice[face[0]], vert_lattice[face[1]], vert_lattice[face[2]],
                                      order);
    }

    for (int l = 1; l <= order - 3; ++l) {
        for (int k = 1; k <= order - l - 2; ++k) {
            for (int j = 1; j <= order - l - k - 1; ++j) {
                out.coords.push_back(Point{double(j) / double(order),
                                           double(k) / double(order),
                                           double(l) / double(order)});
                out.lattice.push_back(Lattice{j, k, l});
            }
        }
    }
    return out;
}

LagrangeNodeLayout generate_hex_nodes(int order) {
    LagrangeNodeLayout out;
    if (order == 0) {
        out.coords.push_back(Point{double(0), double(0), double(0)});
        out.lattice.push_back(Lattice{0, 0, 0});
        return out;
    }

    const Point verts[] = {
        Point{double(-1), double(-1), double(-1)},
        Point{double(1), double(-1), double(-1)},
        Point{double(1), double(1), double(-1)},
        Point{double(-1), double(1), double(-1)},
        Point{double(-1), double(-1), double(1)},
        Point{double(1), double(-1), double(1)},
        Point{double(1), double(1), double(1)},
        Point{double(-1), double(1), double(1)},
    };
    const Lattice vert_lattice[] = {
        Lattice{0, 0, 0},
        Lattice{order, 0, 0},
        Lattice{order, order, 0},
        Lattice{0, order, 0},
        Lattice{0, 0, order},
        Lattice{order, 0, order},
        Lattice{order, order, order},
        Lattice{0, order, order},
    };

    out.coords.reserve(static_cast<std::size_t>((order + 1) * (order + 1) * (order + 1)));
    out.lattice.reserve(static_cast<std::size_t>((order + 1) * (order + 1) * (order + 1)));
    for (std::size_t v = 0; v < 8u; ++v) {
        out.coords.push_back(verts[v]);
        out.lattice.push_back(vert_lattice[v]);
    }

    const int edges[12][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0},
        {4, 5}, {5, 6}, {6, 7}, {7, 4},
        {0, 4}, {1, 5}, {2, 6}, {3, 7},
    };
    // Edge-interior nodes at the Gauss-Lobatto-Legendre position of their lattice
    // index on each axis (line_coord_pm_one), consistent with the corner, face, and
    // interior strata and with the 1D tensor-axis nodes the evaluator uses. (A plain
    // equispaced interpolation along the edge would disagree with the GLL faces and
    // interior at order >= 3.)
    for (const auto& edge : edges) {
        for (int m = 1; m < order; ++m) {
            const Lattice idx =
                lerp_lattice(vert_lattice[edge[0]], vert_lattice[edge[1]], m, order);
            out.coords.push_back(Point{line_coord_pm_one(idx[0], order),
                                       line_coord_pm_one(idx[1], order),
                                       line_coord_pm_one(idx[2], order)});
            out.lattice.push_back(idx);
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
            out.coords.push_back(Point{double(-1), line_coord_pm_one(j, order), line_coord_pm_one(k, order)});
            out.lattice.push_back(Lattice{0, j, k});
        }
    }
    // +X face (x = +1)
    for (int k = 1; k < order; ++k) {
        for (int j = 1; j < order; ++j) {
            out.coords.push_back(Point{double(1), line_coord_pm_one(j, order), line_coord_pm_one(k, order)});
            out.lattice.push_back(Lattice{order, j, k});
        }
    }
    // -Y face (y = -1)
    for (int k = 1; k < order; ++k) {
        for (int i = 1; i < order; ++i) {
            out.coords.push_back(Point{line_coord_pm_one(i, order), double(-1), line_coord_pm_one(k, order)});
            out.lattice.push_back(Lattice{i, 0, k});
        }
    }
    // +Y face (y = +1)
    for (int k = 1; k < order; ++k) {
        for (int i = order - 1; i >= 1; --i) {
            out.coords.push_back(Point{line_coord_pm_one(i, order), double(1), line_coord_pm_one(k, order)});
            out.lattice.push_back(Lattice{i, order, k});
        }
    }
    // -Z face (z = -1)
    for (int j = 1; j < order; ++j) {
        for (int i = 1; i < order; ++i) {
            out.coords.push_back(Point{line_coord_pm_one(i, order), line_coord_pm_one(j, order), double(-1)});
            out.lattice.push_back(Lattice{i, j, 0});
        }
    }
    // +Z face (z = +1)
    for (int j = 1; j < order; ++j) {
        for (int i = 1; i < order; ++i) {
            out.coords.push_back(Point{line_coord_pm_one(i, order), line_coord_pm_one(j, order), double(1)});
            out.lattice.push_back(Lattice{i, j, order});
        }
    }
    for (int k = 1; k < order; ++k) {
        for (int j = 1; j < order; ++j) {
            for (int i = 1; i < order; ++i) {
                out.coords.push_back(Point{line_coord_pm_one(i, order),
                                           line_coord_pm_one(j, order),
                                           line_coord_pm_one(k, order)});
                out.lattice.push_back(Lattice{i, j, k});
            }
        }
    }
    return out;
}

LagrangeNodeLayout generate_wedge_nodes(int order) {
    LagrangeNodeLayout out;
    if (order == 0) {
        out.coords.push_back(Point{double(1) / double(3), double(1) / double(3), double(0)});
        out.lattice.push_back(Lattice{0, 0, 0});
        return out;
    }

    const Point verts[] = {
        Point{double(0), double(0), double(-1)},
        Point{double(1), double(0), double(-1)},
        Point{double(0), double(1), double(-1)},
        Point{double(0), double(0), double(1)},
        Point{double(1), double(0), double(1)},
        Point{double(0), double(1), double(1)},
    };
    const Lattice vert_lattice[] = {
        Lattice{0, 0, 0},
        Lattice{order, 0, 0},
        Lattice{0, order, 0},
        Lattice{0, 0, order},
        Lattice{order, 0, order},
        Lattice{0, order, order},
    };

    out.coords.reserve(static_cast<std::size_t>((order + 1) * (order + 1) * (order + 2) / 2));
    out.lattice.reserve(static_cast<std::size_t>((order + 1) * (order + 1) * (order + 2) / 2));
    for (std::size_t v = 0; v < 6u; ++v) {
        out.coords.push_back(verts[v]);
        out.lattice.push_back(vert_lattice[v]);
    }

    const int edges[9][2] = {
        {0, 1}, {1, 2}, {2, 0},
        {3, 4}, {4, 5}, {5, 3},
        {0, 3}, {1, 4}, {2, 5},
    };
    // The triangle cross-section (x, y) keeps its equispaced simplex placement; the
    // through-axis (z) uses the Gauss-Lobatto-Legendre node of the lattice index, so
    // the prism's tensor axis matches the 1D nodes the evaluator uses. (Triangle
    // edges have z lattice 0 or `order`, for which line_coord_pm_one is -1 / +1, so
    // their z is unchanged.)
    for (const auto& edge : edges) {
        for (int m = 1; m < order; ++m) {
            const double t = static_cast<double>(m) / static_cast<double>(order);
            const Lattice idx =
                lerp_lattice(vert_lattice[edge[0]], vert_lattice[edge[1]], m, order);
            Point coord = verts[edge[0]] * (double(1) - t) + verts[edge[1]] * t;
            coord[2] = line_coord_pm_one(idx[2], order);
            out.coords.push_back(coord);
            out.lattice.push_back(idx);
        }
    }

    append_triangle_face_interior(out, verts[0], verts[1], verts[2],
                                  vert_lattice[0], vert_lattice[1], vert_lattice[2], order);
    append_triangle_face_interior(out, verts[3], verts[4], verts[5],
                                  vert_lattice[3], vert_lattice[4], vert_lattice[5], order);

    for (int r = 1; r < order; ++r) {
        const double z = line_coord_pm_one(r, order);
        for (int m = 1; m < order; ++m) {
            const double t = static_cast<double>(m) / static_cast<double>(order);
            out.coords.push_back(Point{t, double(0), z});
            out.lattice.push_back(Lattice{m, 0, r});
        }
        for (int m = 1; m < order; ++m) {
            const double t = static_cast<double>(m) / static_cast<double>(order);
            out.coords.push_back(Point{double(1) - t, t, z});
            out.lattice.push_back(Lattice{order - m, m, r});
        }
        for (int m = 1; m < order; ++m) {
            const double t = static_cast<double>(m) / static_cast<double>(order);
            out.coords.push_back(Point{double(0), double(1) - t, z});
            out.lattice.push_back(Lattice{0, order - m, r});
        }
    }

    for (int r = 1; r < order; ++r) {
        const double z = line_coord_pm_one(r, order);
        for (int c = 1; c <= order - 2; ++c) {
            for (int b = 1; b <= order - c - 1; ++b) {
                out.coords.push_back(Point{double(b) / double(order),
                                           double(c) / double(order),
                                           z});
                out.lattice.push_back(Lattice{b, c, r});
            }
        }
    }
    return out;
}

LagrangeNodeLayout complete_lagrange_nodes(ElementType canonical_type, int order) {
    svmp::throw_if<BasisConfigurationException>(order < 0, SVMP_HERE,
                                              "ReferenceNodeLayout requires non-negative Lagrange order");
    const ElementType type = canonical_lagrange_type(canonical_type);
    switch (type) {
        case ElementType::Point1: {
            LagrangeNodeLayout out;
            out.coords.push_back(Point{double(0), double(0), double(0)});
            out.lattice.push_back(Lattice{0, 0, 0});
            return out;
        }
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

// Topological interior dimension of an integer lattice node: the number of
// independent directions in which the point sits in the relative interior of
// the reference cell. A vertex gives 0, an edge-interior node 1, a
// face-interior node 2, and a volume-interior node 3.
int serendipity_interior_dim(BasisTopology top, const Lattice& idx, int order) {
    const auto tensor_interior = [order](int v) { return (v > 0 && v < order) ? 1 : 0; };
    switch (top) {
        case BasisTopology::Quadrilateral:
            return tensor_interior(idx[0]) + tensor_interior(idx[1]);
        case BasisTopology::Hexahedron:
            return tensor_interior(idx[0]) + tensor_interior(idx[1]) +
                   tensor_interior(idx[2]);
        case BasisTopology::Wedge: {
            // (idx[0], idx[1]) is the triangle cross-section with implied third
            // barycentric index k; idx[2] is the tensor through-axis. A triangle
            // vertex contributes 0, a triangle edge 1, and the triangle interior 2.
            const int i = idx[0];
            const int j = idx[1];
            const int k = order - i - j;
            const bool tri_vertex = (i == order) || (j == order) || (i + j == 0);
            const bool tri_interior = (i > 0) && (j > 0) && (k > 0);
            const int tri_dim = tri_vertex ? 0 : (tri_interior ? 2 : 1);
            return tri_dim + tensor_interior(idx[2]);
        }
        default:
            return 0;
    }
}

// Build a serendipity reference layout (Quad8, Hex20, Wedge15) from the complete
// quadratic layout of the same topology. Serendipity layouts keep only the
// element's vertices and edge midpoints and drop the face- and volume-interior
// nodes; the complete-quadratic generators emit the vertex/edge nodes first, so
// the serendipity set is exactly the leading keep_count nodes.
std::vector<Point> serendipity_subset_nodes(BasisTopology top,
                                            LagrangeNodeLayout complete,
                                            std::size_t keep_count,
                                            std::size_t complete_count) {
    constexpr int kQuadraticOrder = 2;
    svmp::throw_if<BasisConstructionException>(
        complete.coords.size() != complete_count ||
            complete.lattice.size() != complete_count,
        SVMP_HERE,
        "ReferenceNodeLayout: unexpected complete-quadratic node count for serendipity layout");
    svmp::throw_if<BasisConstructionException>(
        keep_count >= complete_count, SVMP_HERE,
        "ReferenceNodeLayout: serendipity node count must be smaller than the complete layout");

    for (std::size_t n = 0; n < complete.lattice.size(); ++n) {
        const bool on_skeleton =
            serendipity_interior_dim(top, complete.lattice[n], kQuadraticOrder) <= 1;
        const bool kept = n < keep_count;
        svmp::throw_if<BasisConstructionException>(
            kept != on_skeleton, SVMP_HERE,
            "ReferenceNodeLayout: serendipity truncation does not separate skeleton nodes from interior nodes");
    }

    std::vector<Point> nodes = std::move(complete.coords);
    nodes.resize(keep_count);
    return nodes;
}

std::vector<Point> element_nodes(ElementType elem_type) {
    const int order = complete_lagrange_alias_order(elem_type);
    if (order >= 0) {
        return complete_lagrange_nodes(elem_type, order).coords;
    }

    switch (elem_type) {
        case ElementType::Quad8:
            return serendipity_subset_nodes(BasisTopology::Quadrilateral,
                                            generate_quad_nodes(2), 8u, 9u);
        case ElementType::Hex20:
            return serendipity_subset_nodes(BasisTopology::Hexahedron,
                                            generate_hex_nodes(2), 20u, 27u);
        case ElementType::Wedge15:
            return serendipity_subset_nodes(BasisTopology::Wedge,
                                            generate_wedge_nodes(2), 15u, 18u);
        case ElementType::Pyramid13:
            svmp::raise<BasisNodeOrderingException>(SVMP_HERE,
                "ReferenceNodeLayout: pyramid node ordering is disabled");
        default:
            svmp::raise<BasisNodeOrderingException>(SVMP_HERE,
                "ReferenceNodeLayout: unknown element type");
    }
}

// Structural invariants the lattice must satisfy, checked before the accessor
// hands it out. These replace the floating-point round-trip's near-equality
// guards with exact integer checks.
void validate_lattice(const LagrangeNodeLayout& layout, ElementType type, int order) {
    svmp::throw_if<BasisConstructionException>(
        layout.coords.size() != layout.lattice.size(), SVMP_HERE,
        "ReferenceNodeLayout: lattice/coordinate count mismatch");

    const BasisTopology top = topology(type);
    for (const auto& idx : layout.lattice) {
        for (std::size_t d = 0; d < 3u; ++d) {
            svmp::throw_if<BasisConstructionException>(
                idx[d] < 0 || idx[d] > order, SVMP_HERE,
                "ReferenceNodeLayout: lattice index outside [0, order]");
        }
        if (top == BasisTopology::Triangle || top == BasisTopology::Tetrahedron) {
            svmp::throw_if<BasisConstructionException>(
                idx[0] + idx[1] + idx[2] > order, SVMP_HERE,
                "ReferenceNodeLayout: simplex lattice index sum exceeds order");
        } else if (top == BasisTopology::Wedge) {
            svmp::throw_if<BasisConstructionException>(
                idx[0] + idx[1] > order, SVMP_HERE,
                "ReferenceNodeLayout: wedge triangle lattice index sum exceeds order");
        }
    }
}

} // namespace

double line_coord_pm_one(int i, int order) {
    if (order <= 0) {
        return double(0);
    }
    return gll_points(order)[static_cast<std::size_t>(i)];
}

math::Vector<double, 3> ReferenceNodeLayout::get_node_coords(ElementType elem_type,
                                                           std::size_t local_node) {
    const auto nodes = element_nodes(elem_type);
    svmp::throw_if<BasisNodeOrderingException>(local_node >= nodes.size(), SVMP_HERE,
                                             "ReferenceNodeLayout::get_node_coords: node index out of range");
    return nodes[local_node];
}

std::size_t ReferenceNodeLayout::num_nodes(ElementType elem_type) {
    return element_nodes(elem_type).size();
}

std::vector<math::Vector<double, 3>>
ReferenceNodeLayout::node_coords(ElementType elem_type) {
    return element_nodes(elem_type);
}

std::vector<math::Vector<double, 3>>
ReferenceNodeLayout::get_lagrange_node_coords(ElementType canonical_type, int order) {
    return complete_lagrange_nodes(canonical_type, order).coords;
}

LagrangeNodeLayout
ReferenceNodeLayout::get_lagrange_lattice(ElementType canonical_type, int order) {
    LagrangeNodeLayout layout = complete_lagrange_nodes(canonical_type, order);
    validate_lattice(layout, canonical_type, order);
    return layout;
}

} // namespace basis
} // namespace FE
} // namespace svmp
