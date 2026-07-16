/* Copyright (c) Stanford University, The Regents of the University of
 * California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "../test_common.h"
#include "FourierInterpolation.h"

#include <cmath>
#include <numbers>

class FourierInterpolationTest : public ::testing::Test {
protected:
  void SetUp() override {}

  void TearDown() override {}

  // Generate temporal values for function f(t) = sin(t) + cos(t) + 0.1*t
  void CreateTemporalValues(int N, double x_start, double x_end,
                            Vector<double> &times, Array<double> &values) {
    const unsigned int n_components = 2;

    times.resize(N);
    values.resize(n_components, N);

    const double step = (x_end - x_start) / (N - 1);
    for (int i = 0; i < N; ++i) {
      times[i] = x_start + i * step;

      for (unsigned int j = 0; j < n_components; ++j) {
        values(j, i) = interpolated(times[i], j);
      }
    }
  }

  // Interpolated function.
  static double interpolated(const double t, const unsigned int component) {
    if (component == 0) {
      return std::sin(t) + std::cos(t) + 0.1 * t;
    } else {
      return 2.0 * std::sin(t) + 3.0 * std::cos(t) + 0.2 * t;
    }
  }

  // Derivative of the interpolated function.
  static double interpolated_derivative(const double t,
                                        const unsigned int component) {
    if (component == 0) {
      return std::cos(t) - std::sin(t) + 0.1;
    } else {
      return 2.0 * std::cos(t) - 3.0 * std::sin(t) + 0.2;
    }
  }
};

TEST_F(FourierInterpolationTest, FromTimeSeries) {
  // Construct a FourierInterpolation instance from a time series and check that
  // the dimensions of the resulting object are correct.

  int N = 100;                   // 100 timesteps
  double x_start = 0.0;          // start time
  double x_end = 2 * std::numbers::pi; // end time

  Vector<double> times;
  Array<double> values;
  CreateTemporalValues(N, x_start, x_end, times, values);

  // Compute the Fourier interpolation.
  const FourierInterpolation gt = FourierInterpolation::from_time_series(
      /* n_fourier_coefficients = */ 3, times, values,
      /* use_ramp = */ false);

  // Check the number of components and Fourier coefficients.
  ASSERT_EQ(gt.get_n_components(), 2);
  ASSERT_EQ(gt.get_n_fourier_coefficients(), 3);
}

TEST_F(FourierInterpolationTest, FromFourierCoefficients) {
  // Construct a FourierInterpolation instance from Fourier coefficients and
  // check that the dimensions of the resulting object are correct, that the
  // coefficients are set correctly, and that the object evaluates to the
  // correct function.

  const Vector<double> q_i = {0.0, 0.0};
  const Vector<double> q_s = {0.1, 0.2};
  const Array<double> c_real = {{0.0, 1.0}, {0.0, 3.0}};
  const Array<double> c_imag = {{0.0, -1.0}, {0.0, -2.0}};
  const double initial_time = 0.0;
  const double period = 2 * std::numbers::pi;

  const FourierInterpolation gt =
      FourierInterpolation::from_fourier_coefficients(q_i, q_s, c_real, c_imag,
                                                      initial_time, period);

  // Check the number of components and Fourier coefficients.
  ASSERT_EQ(gt.get_n_components(), 2);
  ASSERT_EQ(gt.get_n_fourier_coefficients(), 2);

  // Check the linear trend coefficients.
  ASSERT_EQ(gt.get_linear_trend_initial_value(0), 0.0);
  ASSERT_EQ(gt.get_linear_trend_initial_value(1), 0.0);
  ASSERT_EQ(gt.get_linear_trend_slope(0), 0.1);
  ASSERT_EQ(gt.get_linear_trend_slope(1), 0.2);

  // Check the Fourier coefficients.
  ASSERT_EQ(gt.get_coefficient_real(0, 0), 0.0);
  ASSERT_EQ(gt.get_coefficient_real(0, 1), 1.0);
  ASSERT_EQ(gt.get_coefficient_real(1, 0), 0.0);
  ASSERT_EQ(gt.get_coefficient_real(1, 1), 3.0);
  ASSERT_EQ(gt.get_coefficient_imaginary(0, 0), 0.0);
  ASSERT_EQ(gt.get_coefficient_imaginary(0, 1), -1.0);
  ASSERT_EQ(gt.get_coefficient_imaginary(1, 0), 0.0);
  ASSERT_EQ(gt.get_coefficient_imaginary(1, 1), -2.0);

  // Check that the interpolation evaluates to the correct function.
  Vector<double> value, derivative;
  for (double t = initial_time; t <= initial_time + period; t += 0.5) {
    std::tie(value, derivative) = gt.value_and_derivative(t);

    const double expected_y = interpolated(t, 0);
    const double expected_dy = interpolated_derivative(t, 0);
    const double expected_z = interpolated(t, 1);
    const double expected_dz = interpolated_derivative(t, 1);

    ASSERT_NEAR(value[0], expected_y, 1e-6);
    ASSERT_NEAR(value[1], expected_z, 1e-6);
    ASSERT_NEAR(derivative[0], expected_dy, 1e-6);
    ASSERT_NEAR(derivative[1], expected_dz, 1e-6);
  }
}

TEST_F(FourierInterpolationTest, ValueAndDerivative) {
  // Evaluate the interpolation and its derivative at a few points and compare
  // with the original function. Since the interpolated function is a sin-cos
  // combination plus a linear term, the interpolation is expected to match it
  // exactly (as long as we use at least two Fourier coefficients).

  int N = 100;                   // 100 timesteps
  double x_start = 0.0;          // start time
  double x_end = 2 * std::numbers::pi; // end time

  Vector<double> times;
  Array<double> values;
  CreateTemporalValues(N, x_start, x_end, times, values);

  // Compute the Fourier interpolation.
  const FourierInterpolation gt = FourierInterpolation::from_time_series(
      /* n_fourier_coefficients = */ 3, times, values,
      /* use_ramp = */ false);

  Vector<double> value, derivative;
  for (double t = x_start; t <= x_end; t += 0.5) {
    const double expected_y = interpolated(t, 0);
    const double expected_dy = interpolated_derivative(t, 0);

    const double expected_z = interpolated(t, 1);
    const double expected_dz = interpolated_derivative(t, 1);

    std::tie(value, derivative) = gt.value_and_derivative(t);

    ASSERT_NEAR(value[0], expected_y, 1e-2) << "Value mismatch at t = " << t;
    ASSERT_NEAR(derivative[0], expected_dy, 1e-2)
        << "Derivative mismatch at t = " << t;
    ASSERT_NEAR(value[1], expected_z, 1e-2) << "Value mismatch at t = " << t;
    ASSERT_NEAR(derivative[1], expected_dz, 1e-2)
        << "Derivative mismatch at t = " << t;
  }

  // Check periodicity.
  Vector<double> value_2, derivative_2;
  std::tie(value, derivative) = gt.value_and_derivative(x_start + 0.1);
  std::tie(value_2, derivative_2) = gt.value_and_derivative(x_end + 0.1);
  ASSERT_NEAR(value[0], value_2[0], 1e-2) << "Periodic value mismatch";
  ASSERT_NEAR(value[1], value_2[1], 1e-2) << "Periodic value mismatch";
  ASSERT_NEAR(derivative[0], derivative_2[0], 1e-2)
      << "Periodic derivative mismatch";
  ASSERT_NEAR(derivative[1], derivative_2[1], 1e-2)
      << "Periodic derivative mismatch";
}

TEST_F(FourierInterpolationTest, Coefficients) {
  // Construct the interpolation and compare the computed Fourier and linear
  // trend coefficients with the expected values. Notice that we're using a
  // period here that is different from that of the interpolated function (i.e.
  // not 2pi). Therefore, the Fourier coefficients will not be the same as those
  // of the original function.

  int N = 100;          // 100 timesteps
  double x_start = 0.0; // start time
  double x_end = 10.0;  // end time

  Vector<double> times;
  Array<double> values;
  CreateTemporalValues(N, x_start, x_end, times, values);

  // Compute the Fourier coefficients
  const FourierInterpolation gt = FourierInterpolation::from_time_series(
      /* n_fourier_coefficients = */ 16, times, values,
      /* use_ramp = */ false);

  // Check the linear trend slope
  ASSERT_NEAR(gt.get_linear_trend_slope(0), -0.13830, 1e-2)
      << "Expected slope ~-0.13830";

  // Check the real and imaginary components of the first three Fourier
  // coefficients
  ASSERT_NEAR(gt.get_coefficient_real(0, 0), 0.32094, 1e-2)
      << "Expected first real coefficient to be close to 0.32094";
  ASSERT_NEAR(gt.get_coefficient_imaginary(0, 0), 0.0, 1e-2)
      << "Expected first imaginary coefficient to be close to 0.0";
  ASSERT_NEAR(gt.get_coefficient_real(0, 1), 0.42759, 1e-2)
      << "Expected second real coefficient to be close to 0.42759";
  ASSERT_NEAR(gt.get_coefficient_imaginary(0, 1), 1.25295, 1e-2)
      << "Expected second imaginary coefficient to be close to 1.25295";
  ASSERT_NEAR(gt.get_coefficient_real(0, 2), -0.44685, 1e-2)
      << "Expected third real coefficient to be close to -0.44685";
  ASSERT_NEAR(gt.get_coefficient_imaginary(0, 2), -0.65403, 1e-2)
      << "Expected third imaginary coefficient to be close to -0.65403";
}

TEST_F(FourierInterpolationTest, Ramp) {
  // Construct a FourierInterpolation instance from a time series and check that
  // the ramp option behaves correctly.

  int N = 100;          // 100 timesteps
  double x_start = 0.0; // start time
  double x_end = 10.0;  // end time

  Vector<double> times;
  Array<double> values;
  CreateTemporalValues(N, x_start, x_end, times, values);

  // Compute the Fourier coefficients
  const FourierInterpolation gt = FourierInterpolation::from_time_series(
      /* n_fourier_coefficients = */ 16, times, values,
      /* use_ramp = */ true);

  // Check that the number of Fourier coefficients is 1 (since ramp is used)
  ASSERT_EQ(gt.get_n_fourier_coefficients(), 1)
      << "Expected number of Fourier coefficients to be 1 when using ramp";

  Vector<double> value, derivative;

  // Check that before the initial time the value is equal to the initial value
  // of the linear trend.
  std::tie(value, derivative) = gt.value_and_derivative(x_start - 1.0);
  ASSERT_NEAR(value[0], values(0, 0), 1e-6)
      << "Expected value before initial time to be equal to initial value of "
         "linear trend";
  ASSERT_NEAR(derivative[0], 0.0, 1e-6)
      << "Expected derivative before initial time to be zero";

  // Check that after the final time the value is equal to the final value of
  // the linear trend.
  std::tie(value, derivative) = gt.value_and_derivative(x_end + 1.0);
  ASSERT_NEAR(value[0], values(0, values.ncols() - 1), 1e-6)
      << "Expected value after final time to be equal to final value of linear "
         "trend";
  ASSERT_NEAR(derivative[0], 0.0, 1e-6)
      << "Expected derivative after final time to be zero";

  // Check that the value in the interpolation interval is a linear
  // interpolation of the initial and final value.
  const double expected_derivative =
      (values(0, values.ncols() - 1) - values(0, 0)) / (x_end - x_start);
  for (double t = x_start; t <= x_end; t += 0.5) {
    std::tie(value, derivative) = gt.value_and_derivative(t);

    const double expected_y =
        values(0, 0) + (values(0, values.ncols() - 1) - values(0, 0)) *
                           (t - x_start) / (x_end - x_start);
    const double expected_z =
        values(1, 0) + (values(1, values.ncols() - 1) - values(1, 0)) *
                           (t - x_start) / (x_end - x_start);

    ASSERT_NEAR(value[0], expected_y, 1e-6)
        << "Expected value in interpolation interval to be linear "
           "interpolation of initial and final value";
    ASSERT_NEAR(value[1], expected_z, 1e-6)
        << "Expected value in interpolation interval to be linear "
           "interpolation of initial and final value";
    ASSERT_NEAR(derivative[0], expected_derivative, 1e-6)
        << "Expected derivative in interpolation interval to be linear "
           "interpolation of initial and final derivative";
  }
}