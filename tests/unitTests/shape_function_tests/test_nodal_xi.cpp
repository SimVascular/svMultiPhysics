// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/// @file test_nodal_xi.cpp
/// @brief Tests for nn::get_nodal_xi, the reference coordinates of the nodes
/// of a parent element.
///
/// get_nodal_xi reads the coordinates from the FE Basis reference node layout
/// and permutes them into solver node order. The tests check them against the
/// shape functions get_gnn evaluates at them, and against the lower-order
/// function space post::tensor_post builds over the corner nodes.

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

/// @brief Element types for which post::tensor_post builds a lower-order function
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

/// @brief Check that the shape functions of an element form a nodal basis at
/// the coordinates get_nodal_xi returns: shape function b evaluated at the
/// coordinates of node a is one for a == b and zero otherwise.
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

/// @brief The lower-order function space built on the corner nodes, evaluated
/// at the nodes of the higher-order element, weights the corner values: the
/// weights sum to one, each lies in [0,1], and at a corner node the weight of
/// that corner is one. post::tensor_post uses these weights to interpolate the
/// corner values to the edge, face and center nodes.
TEST(NodalXi, LowerOrderBasisWeightsCornerValues)
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

      if (a < fs.eNoN) {
        for (int b = 0; b < fs.eNoN; b++) {
          EXPECT_NEAR(N(b), (a == b) ? 1.0 : 0.0, 1.0e-12) << "corner node " << a << ", shape function " << b;
        }
      }
    }
  }
}

/// @brief Each mid-edge node of a TET10 weights the two corner nodes of its
/// own edge at 0.5 and the other two at zero, so post::tensor_post gives it the
/// average of those two corner values. The six edges are (0,1), (1,2), (0,2),
/// (0,3), (1,3) and (2,3) in solver node order.
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
