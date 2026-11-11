// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/// @file test_nodal_xi.cpp
/// @brief Tests for nn::get_nodal_xi, the reference coordinates of the nodes
/// of a parent element.
///
/// These coordinates are used by post::tpost to interpolate quantities
/// computed at the corner nodes of a higher-order element to its edge, face
/// and center nodes. They must be read from the parent element and not
/// recovered from the physical node coordinates: on a curved element the edge
/// nodes lie outside the straight-sided element spanned by the corner nodes.

#include "gtest/gtest.h"

#include "consts.h"
#include "fs.h"
#include "nn.h"

#include <vector>

using namespace consts;

namespace {

struct ElementInfo {
  ElementType eType;
  int insd;
  int eNoN;
  const char* name;
};

/// @brief Element types for which post::tpost builds a lower-order function
/// space over the corner nodes (see fs::set_thood_fs).
const std::vector<ElementInfo> higher_order_elements = {
  {ElementType::TRI6,  2,  6, "TRI6" },
  {ElementType::QUD8,  2,  8, "QUD8" },
  {ElementType::QUD9,  2,  9, "QUD9" },
  {ElementType::TET10, 3, 10, "TET10"},
  {ElementType::HEX20, 3, 20, "HEX20"},
  {ElementType::HEX27, 3, 27, "HEX27"},
};

const std::vector<ElementInfo> lower_order_elements = {
  {ElementType::TRI3,  2, 3, "TRI3"},
  {ElementType::QUD4,  2, 4, "QUD4"},
  {ElementType::TET4,  3, 4, "TET4"},
  {ElementType::HEX8,  3, 8, "HEX8"},
};

/// @brief The shape functions of an element are a nodal (Lagrange) basis: N_b
/// evaluated at the reference coordinates of node a is delta_ab.
///
/// This ties get_nodal_xi to get_gnn, so a wrong coordinate or a mismatch
/// between the solver and FE Basis node orderings fails here.
void expect_nodal_basis(const ElementInfo& elem)
{
  SCOPED_TRACE(elem.name);

  Array<double> xiN(elem.insd, elem.eNoN);
  nn::get_nodal_xi(elem.insd, elem.eType, elem.eNoN, xiN);

  for (int a = 0; a < elem.eNoN; a++) {
    auto xi = xiN.col(a);
    Vector<double> N(elem.eNoN);
    Array<double> Nx(elem.insd, elem.eNoN);
    nn::get_gnn(elem.insd, elem.eType, elem.eNoN, xi, N, Nx);

    for (int b = 0; b < elem.eNoN; b++) {
      EXPECT_NEAR(N(b), (a == b) ? 1.0 : 0.0, 1.0e-12) << "node " << a << ", shape function " << b;
    }
  }
}

TEST(NodalXi, IsNodalBasisOfHigherOrderElement)
{
  for (const auto& elem : higher_order_elements) {
    expect_nodal_basis(elem);
  }
}

TEST(NodalXi, IsNodalBasisOfLowerOrderElement)
{
  for (const auto& elem : lower_order_elements) {
    expect_nodal_basis(elem);
  }
}

/// @brief The lower-order (Taylor-Hood) basis built on the corner nodes, when
/// evaluated at the nodes of the higher-order element, is a partition of unity
/// with no extrapolation that reproduces the corner values at the corner
/// nodes. This is the property post::tpost relies on to interpolate to the
/// edge, face and center nodes, and the one the old inverse map broke on
/// curved elements, where it produced weights outside [0,1].
TEST(NodalXi, LowerOrderBasisInterpolatesCornerValues)
{
  for (const auto& elem : higher_order_elements) {
    SCOPED_TRACE(elem.name);

    fsType fs;
    fs::set_thood_fs(fs, elem.eType);

    Array<double> xiN(elem.insd, elem.eNoN);
    nn::get_nodal_xi(elem.insd, elem.eType, elem.eNoN, xiN);

    for (int a = 0; a < elem.eNoN; a++) {
      auto xi = xiN.col(a);
      Vector<double> N(fs.eNoN);
      Array<double> Nx(elem.insd, fs.eNoN);
      nn::get_gnn(elem.insd, fs.eType, fs.eNoN, xi, N, Nx);

      double sum = 0.0;
      for (int b = 0; b < fs.eNoN; b++) {
        sum += N(b);
        EXPECT_GE(N(b), -1.0e-12) << "node " << a << ", corner " << b;
        EXPECT_LE(N(b), 1.0 + 1.0e-12) << "node " << a << ", corner " << b;
      }
      EXPECT_NEAR(sum, 1.0, 1.0e-12) << "node " << a;

      // At a corner node the lower-order basis is nodal as well.
      if (a < fs.eNoN) {
        for (int b = 0; b < fs.eNoN; b++) {
          EXPECT_NEAR(N(b), (a == b) ? 1.0 : 0.0, 1.0e-12) << "corner node " << a << ", shape function " << b;
        }
      }
    }
  }
}

/// @brief A mid-edge node of a quadratic element is the average of the two
/// corner nodes of its edge in the lower-order basis, so post::tpost assigns
/// it the average of the two corner values. Spot-check on TET10, whose six
/// edges are (0,1), (1,2), (0,2), (0,3), (1,3), (2,3) in solver node order.
TEST(NodalXi, Tet10EdgeNodesAverageTheirCorners)
{
  const int edges[6][2] = {{0,1}, {1,2}, {0,2}, {0,3}, {1,3}, {2,3}};

  fsType fs;
  fs::set_thood_fs(fs, ElementType::TET10);

  Array<double> xiN(3, 10);
  nn::get_nodal_xi(3, ElementType::TET10, 10, xiN);

  for (int k = 0; k < 6; k++) {
    auto xi = xiN.col(4 + k);
    Vector<double> N(fs.eNoN);
    Array<double> Nx(3, fs.eNoN);
    nn::get_gnn(3, fs.eType, fs.eNoN, xi, N, Nx);

    for (int b = 0; b < fs.eNoN; b++) {
      const bool on_edge = (b == edges[k][0]) || (b == edges[k][1]);
      EXPECT_NEAR(N(b), on_edge ? 0.5 : 0.0, 1.0e-12) << "edge node " << 4 + k << ", corner " << b;
    }
  }
}

} // namespace
