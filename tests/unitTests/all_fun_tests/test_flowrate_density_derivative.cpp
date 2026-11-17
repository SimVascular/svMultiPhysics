// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file test_flowrate_density_derivative.cpp
 * @brief Tests of all_fun::flowrate_density_derivative_wrt_position.
 *
 * @test Elements/FlowrateDensityDerivative.MatchesFiniteDifferences compares
 * the derivative with central finite differences of the flow rate density,
 * for elements in two and three dimensions.
 *
 * @test FlowrateDensityDerivativeLinearTriangle.MatchesClosedForm compares
 * the derivative on a linear triangle with its closed form.
 *
 * @par Running the tests
 * @code{.sh}
 * ./run_all_unit_tests --gtest_filter='*FlowrateDensityDerivative*'
 * @endcode
 */

#include "all_fun.h"
#include "utils.h"

#include <gtest/gtest.h>

#include <random>
#include <tuple>

namespace {

/// Fill an array with values between -1 and 1.
void fill_random(Array<double>& array, std::mt19937& generator)
{
  std::uniform_real_distribution<double> distribution(-1.0, 1.0);

  for (int j = 0; j < array.ncols(); j++) {
    for (int i = 0; i < array.nrows(); i++) {
      array(i, j) = distribution(generator);
    }
  }
}

/// Tangent vectors at a Gauss point (nsd x insd): tangent 'k' is the sum over
/// the element nodes of Nx(k,a) times the position of node 'a'.
Array<double> tangent_vectors(const Array<double>& positions, const Array<double>& Nx)
{
  Array<double> tangents(positions.nrows(), Nx.nrows());

  for (int a = 0; a < Nx.ncols(); a++) {
    for (int k = 0; k < Nx.nrows(); k++) {
      for (int j = 0; j < positions.nrows(); j++) {
        tangents(j, k) = tangents(j, k) + Nx(k, a) * positions(j, a);
      }
    }
  }

  return tangents;
}

/// Flow rate density at a Gauss point of an element with the given node
/// positions (nsd x number of element nodes).
double flowrate_density(const Array<double>& positions, const Array<double>& Nx, const Vector<double>& velocity)
{
  return velocity * utils::cross(tangent_vectors(positions, Nx));
}

// Number of space dimensions and number of element nodes.
using ElementParameters = std::tuple<int, int>;
class FlowrateDensityDerivative : public testing::TestWithParam<ElementParameters> {};

TEST_P(FlowrateDensityDerivative, MatchesFiniteDifferences)
{
  const auto [nsd, num_nodes] = GetParam();
  const int insd = nsd - 1;

  std::mt19937 generator(0);

  Array<double> positions(nsd, num_nodes);
  Array<double> Nx(insd, num_nodes);
  Array<double> velocity_values(nsd, 1);
  fill_random(positions, generator);
  fill_random(Nx, generator);
  fill_random(velocity_values, generator);
  const Vector<double> velocity = velocity_values.col(0);

  const auto derivative = all_fun::flowrate_density_derivative_wrt_position(
      tangent_vectors(positions, Nx), Nx, velocity);

  ASSERT_EQ(derivative.nrows(), nsd);
  ASSERT_EQ(derivative.ncols(), num_nodes);

  // The flow rate density is at most quadratic in each coordinate, so central
  // differences have no truncation error, only round-off.
  const double step = 1.0e-3;

  for (int a = 0; a < num_nodes; a++) {
    for (int i = 0; i < nsd; i++) {
      Array<double> forward(positions);
      Array<double> backward(positions);
      forward(i, a) = forward(i, a) + step;
      backward(i, a) = backward(i, a) - step;

      const double finite_difference =
          (flowrate_density(forward, Nx, velocity) - flowrate_density(backward, Nx, velocity)) / (2.0 * step);

      EXPECT_NEAR(derivative(i, a), finite_difference, 1.0e-9)
          << "coordinate " << i << " of node " << a;
    }
  }
}

INSTANTIATE_TEST_SUITE_P(Elements, FlowrateDensityDerivative,
                         testing::Values(ElementParameters{2, 2}, ElementParameters{2, 3},
                                         ElementParameters{3, 3}, ElementParameters{3, 4},
                                         ElementParameters{3, 6}));

TEST(FlowrateDensityDerivativeLinearTriangle, MatchesClosedForm)
{
  // Node positions, one per column.
  const Array<double> positions{{0.1, 1.2, -0.3},
                                {0.4, 0.2, 1.5},
                                {-0.2, 0.3, 0.6}};

  // Shape function derivatives of the linear triangle.
  const Array<double> Nx{{-1.0, 1.0, 0.0},
                         {-1.0, 0.0, 1.0}};

  const Vector<double> velocity{0.7, -0.4, 1.1};

  const auto derivative = all_fun::flowrate_density_derivative_wrt_position(
      tangent_vectors(positions, Nx), Nx, velocity);

  // With nodes x_0, x_1, x_2 the flow rate density is
  // velocity . ((x_1 - x_0) x (x_2 - x_0)), and its derivative with respect
  // to x_a is (x_b - x_c) x velocity, with (a, b, c) a cyclic permutation of
  // (0, 1, 2).
  for (int a = 0; a < 3; a++) {
    const int b = (a + 1) % 3;
    const int c = (a + 2) % 3;

    Array<double> edge_and_velocity(3, 2);
    for (int i = 0; i < 3; i++) {
      edge_and_velocity(i, 0) = positions(i, b) - positions(i, c);
      edge_and_velocity(i, 1) = velocity(i);
    }
    const auto expected = utils::cross(edge_and_velocity);

    for (int i = 0; i < 3; i++) {
      EXPECT_NEAR(derivative(i, a), expected(i), 1.0e-13)
          << "coordinate " << i << " of node " << a;
    }
  }
}

} // namespace
