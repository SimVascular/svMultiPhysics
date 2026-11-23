// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file test_derivative_of_outward_normal.cpp
 * @brief Tests of nn::derivative_of_outward_normal_wrt_position.
 *
 * @par Setup
 * A mesh of one linear element (a triangle in two dimensions, a tetrahedron in
 * three) with one face, opposite to the last node, in a displaced
 * configuration. The cases combine the number of space dimensions with the
 * order of the face nodes: reversing it reverses the cross product of the
 * tangent vectors, so one of the two orders needs the normal to be flipped to
 * point outward.
 *
 * @test Elements/DerivativeOfOutwardNormal.MatchesFiniteDifferencesOfNormal
 * compares the derivative with central finite differences of nn::gnnb.
 *
 * @test Elements/DerivativeOfOutwardNormal.MatchesClosedForm compares the
 * derivative with its closed form on the segment and on the triangle.
 *
 * @par Running the tests
 * @code{.sh}
 * ./run_all_unit_tests --gtest_filter='*DerivativeOfOutwardNormal*'
 * @endcode
 */

#include "ComMod.h"
#include "SolutionStates.h"
#include "nn.h"

#include <gtest/gtest.h>

#include <array>
#include <tuple>
#include <utility>

namespace {

// Number of space dimensions, and whether the face nodes are in reverse order.
using ElementParameters = std::tuple<int, bool>;

class DerivativeOfOutwardNormal : public testing::TestWithParam<ElementParameters> {
protected:
  void SetUp() override
  {
    const auto [num_dimensions, reversed] = GetParam();
    nsd = num_dimensions;
    const int num_nodes = nsd + 1;
    const int num_face_nodes = nsd;

    com_mod.nsd = nsd;
    com_mod.mvMsh = false;

    // Node positions in the reference configuration. The last node is the
    // one that is not on the face.
    const double positions[4][3] = {{0.1, 0.2, -0.1},
                                    {1.2, 0.1, 0.3},
                                    {-0.2, 1.1, 0.2},
                                    {0.3, 0.2, 1.4}};

    com_mod.x.resize(nsd, num_nodes);
    auto& displacement = solutions.current.get_displacement();
    displacement.resize(nsd, num_nodes);
    solutions.old.get_displacement().resize(nsd, num_nodes);

    for (int a = 0; a < num_nodes; a++) {
      for (int i = 0; i < nsd; i++) {
        com_mod.x(i, a) = positions[a][i];
        displacement(i, a) = 0.03 * (i + 1) * (a + 1) * ((a % 2 == 0) ? 1.0 : -1.0);
      }
    }

    com_mod.msh.resize(1);
    auto& mesh = com_mod.msh[0];
    mesh.name = "mesh";
    mesh.lShl = false;
    mesh.eNoN = num_nodes;
    mesh.IEN.resize(num_nodes, 1);
    for (int a = 0; a < num_nodes; a++) {
      mesh.IEN(a, 0) = a;
    }

    face.name = "face";
    face.iM = 0;
    face.nEl = 1;
    face.eNoN = num_face_nodes;
    face.gE.resize(1);
    face.gE(0) = 0;
    face.IEN.resize(num_face_nodes, 1);
    for (int a = 0; a < num_face_nodes; a++) {
      face.IEN(a, 0) = a;
    }
    if (reversed) {
      std::swap(face.IEN(num_face_nodes - 2, 0), face.IEN(num_face_nodes - 1, 0));
    }

    // Shape function derivatives of the linear line and triangle.
    Nx.resize(nsd - 1, num_face_nodes);
    if (nsd == 2) {
      Nx(0, 0) = -1.0;
      Nx(0, 1) = 1.0;
    } else {
      Nx(0, 0) = -1.0;
      Nx(0, 1) = 1.0;
      Nx(1, 0) = -1.0;
      Nx(1, 2) = 1.0;
    }
  }

  /// Position of a mesh node in the displaced configuration.
  std::array<double, 3> displaced_position(const int node) const
  {
    std::array<double, 3> position{0.0, 0.0, 0.0};
    for (int i = 0; i < nsd; i++) {
      position[i] = com_mod.x(i, node) + solutions.current.get_displacement()(i, node);
    }
    return position;
  }

  const consts::MechanicalConfigurationType cfg = consts::MechanicalConfigurationType::new_timestep;
  const unsigned int displacement_index = 0;

  int nsd = 0;
  ComMod com_mod;
  faceType face;
  SolutionStates solutions;
  Array<double> Nx;
};

TEST_P(DerivativeOfOutwardNormal, MatchesFiniteDifferencesOfNormal)
{
  const auto derivative = nn::derivative_of_outward_normal_wrt_position(
      com_mod, face, 0, Nx, solutions, cfg, displacement_index);

  ASSERT_EQ(derivative.nrows(), nsd);
  ASSERT_EQ(derivative.ncols(), nsd);
  ASSERT_EQ(derivative.nslices(), face.eNoN);

  // The normal is at most quadratic in each coordinate, so central
  // differences have no truncation error, only round-off.
  const double step = 1.0e-3;
  auto& displacement = solutions.current.get_displacement();

  for (int a = 0; a < face.eNoN; a++) {
    const int node = face.IEN(a, 0);

    for (int i = 0; i < nsd; i++) {
      const double original = displacement(i, node);

      displacement(i, node) = original + step;
      const auto forward = nn::gnnb(com_mod, face, 0, 0, Nx, solutions, cfg, displacement_index);

      displacement(i, node) = original - step;
      const auto backward = nn::gnnb(com_mod, face, 0, 0, Nx, solutions, cfg, displacement_index);

      displacement(i, node) = original;

      for (int j = 0; j < nsd; j++) {
        const double finite_difference = (forward(j) - backward(j)) / (2.0 * step);

        EXPECT_NEAR(derivative(j, i, a), finite_difference, 1.0e-10)
            << "component " << j << " of the normal, coordinate " << i << " of face node " << a;
      }
    }
  }
}

TEST_P(DerivativeOfOutwardNormal, MatchesClosedForm)
{
  const auto derivative = nn::derivative_of_outward_normal_wrt_position(
      com_mod, face, 0, Nx, solutions, cfg, displacement_index);

  // Positions of the face nodes, in the order of the face, and of the node
  // that is not on the face. In two dimensions the third coordinate is zero.
  std::array<std::array<double, 3>, 3> x{};
  for (int a = 0; a < face.eNoN; a++) {
    x[a] = displaced_position(face.IEN(a, 0));
  }
  const auto interior = displaced_position(nsd);

  const auto cross = [](const std::array<double, 3>& u, const std::array<double, 3>& v) {
    return std::array<double, 3>{u[1] * v[2] - u[2] * v[1],
                                 u[2] * v[0] - u[0] * v[2],
                                 u[0] * v[1] - u[1] * v[0]};
  };
  const auto difference = [](const std::array<double, 3>& u, const std::array<double, 3>& v) {
    return std::array<double, 3>{u[0] - v[0], u[1] - v[1], u[2] - v[2]};
  };

  // The normal is 'sign' times the vector below, where the sign makes it
  // point away from the node that is not on the face:
  // - on a triangle, (x_1 - x_0) x (x_2 - x_0);
  // - on a segment, the tangent t = x_1 - x_0 rotated, (t_y, -t_x).
  std::array<double, 3> raw_normal;
  if (nsd == 3) {
    raw_normal = cross(difference(x[1], x[0]), difference(x[2], x[0]));
  } else {
    const auto tangent = difference(x[1], x[0]);
    raw_normal = {tangent[1], -tangent[0], 0.0};
  }

  const auto outward = difference(x[0], interior);
  const double sign =
      (raw_normal[0] * outward[0] + raw_normal[1] * outward[1] + raw_normal[2] * outward[2] < 0.0) ? -1.0 : 1.0;

  for (int a = 0; a < face.eNoN; a++) {
    for (int i = 0; i < nsd; i++) {
      std::array<double, 3> expected;

      if (nsd == 3) {
        // The derivative with respect to coordinate 'i' of x_a is
        // e_i x (x_b - x_c), with (a, b, c) a cyclic permutation of (0, 1, 2).
        const int b = (a + 1) % 3;
        const int c = (a + 2) % 3;
        std::array<double, 3> unit_vector{0.0, 0.0, 0.0};
        unit_vector[i] = 1.0;
        expected = cross(unit_vector, difference(x[b], x[c]));
      } else {
        // The derivative with respect to coordinate 'i' of x_1 is the i-th
        // unit vector rotated, and the one with respect to x_0 is its
        // opposite.
        const double node_sign = (a == 1) ? 1.0 : -1.0;
        if (i == 0) {
          expected = {0.0, -node_sign, 0.0};
        } else {
          expected = {node_sign, 0.0, 0.0};
        }
      }

      for (int j = 0; j < nsd; j++) {
        EXPECT_NEAR(derivative(j, i, a), sign * expected[j], 1.0e-13)
            << "component " << j << " of the normal, coordinate " << i << " of face node " << a;
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P(Elements, DerivativeOfOutwardNormal,
                         testing::Values(ElementParameters{2, false}, ElementParameters{2, true},
                                         ElementParameters{3, false}, ElementParameters{3, true}));

} // namespace
