// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "FourierInterpolation.h"
#include "consts.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <numbers>
#include <regex>
#include <sstream>

namespace {
/// Clean a line by removing carriage returns and leading/trailing spaces, and
/// replacing multiple spaces with a single space.
std::string clean_line(std::string line) {
  line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
  return std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1");
}

/// Read a double from a stream, raising an exception if the read fails.
double read_double(std::istream &stream, const std::string &file_name,
                   const std::string &context = "") {
  double value;
  stream >> value;

  svmp::throw_if<svmp::FileFormatException>(
      stream.fail(),
      "Cannot parse integer value" + (context != "" ? " for " + context : ""),
      file_name);

  return value;
}

/// Read an integer from a stream, raising an exception if the read fails.
int read_int(std::istream &stream, const std::string &file_name,
             const std::string &context = "") {
  int value;
  stream >> value;

  svmp::throw_if<svmp::FileFormatException>(
      stream.fail(),
      "Cannot parse integer value" + (context != "" ? " for " + context : ""),
      file_name);

  return value;
}

/// Convert a string into a vector of doubles.
std::vector<double> string_to_vector_double(const std::string &str,
                                            const std::string &file_name,
                                            unsigned int line_number) {
  std::istringstream line_string_stream(str);
  std::vector<double> values;

  while (!line_string_stream.eof()) {
    values.push_back(read_double(line_string_stream, file_name,
                                 " on line " + std::to_string(line_number)));
  }

  return values;
}

} // namespace

FourierInterpolation FourierInterpolation::from_time_series(
    const unsigned int n_fourier_coefficients, const Vector<double> &times,
    const Array<double> &values, bool use_ramp) {
  const unsigned int n_time_points = times.size();
  const unsigned int n_components = values.nrows();

  svmp::throw_if<svmp::FE::InvalidArgumentException>(
      n_time_points < 2, "At least two time points are needed to construct a "
                         "FourierInterpolation object.");

  // Check that times and values have compatible sizes.
  svmp::throw_if<svmp::FE::InvalidArgumentException>(
      values.ncols() != times.size(),
      "The number of columns of values must match the size of times.");

  // Check that times are in strictly increasing order. Notice that this implies
  // that the period is strictly positive, since it is computed as the
  // difference between the last and first time value.
  for (unsigned int i = 1; i < n_time_points; ++i) {
    svmp::throw_if<svmp::FE::InvalidArgumentException>(
        times[i] <= times[i - 1],
        "The time values must be in strictly increasing order. Got times[" +
            std::to_string(i - 1) + "] = " + std::to_string(times[i - 1]) +
            " and times[" + std::to_string(i) +
            "] = " + std::to_string(times[i]) + ".");
  }

  FourierInterpolation result;

  result.use_ramp = use_ramp;
  result.n_components = n_components;
  result.n_fourier_coefficients = use_ramp ? 1 : n_fourier_coefficients;
  result.initial_time = times[0];
  result.period = times[n_time_points - 1] - times[0];
  result.linear_trend_initial_values.resize(n_components);
  result.linear_trend_slopes.resize(n_components);
  result.fourier_coefficients_real.resize(n_components,
                                          result.n_fourier_coefficients);
  result.fourier_coefficients_imaginary.resize(n_components,
                                               result.n_fourier_coefficients);

  // Compute the linear trend part.
  for (unsigned int j = 0; j < n_components; ++j) {
    result.linear_trend_initial_values[j] = values(j, 0);
    result.linear_trend_slopes[j] =
        (values(j, n_time_points - 1) - values(j, 0)) / result.period;
  }

  // Subtract the linear trend part from the values.
  Vector<double> times_shifted = times;
  Array<double> values_shifted = values;
  for (unsigned int i = 0; i < n_time_points; ++i) {
    times_shifted[i] -= result.initial_time;
    for (unsigned int j = 0; j < n_components; ++j) {
      values_shifted(j, i) = values_shifted(j, i) -
                             result.linear_trend_initial_values[j] -
                             result.linear_trend_slopes[j] * times_shifted[i];
    }
  }

  // Compute the Fourier coefficients.
  for (unsigned int n = 0; n < result.n_fourier_coefficients; ++n) {
    const double tmp = static_cast<double>(n);
    result.fourier_coefficients_real.set_col(n, 0.0);
    result.fourier_coefficients_imaginary.set_col(n, 0.0);

    for (unsigned int i = 0; i < n_time_points - 1; ++i) {
      const double ko =
          2.0 * std::numbers::pi * tmp * times_shifted[i] / result.period;
      const double kn =
          2.0 * std::numbers::pi * tmp * times_shifted[i + 1] / result.period;

      for (unsigned int j = 0; j < result.n_components; j++) {
        if (n == 0) {
          result.fourier_coefficients_real(j, n) +=
              0.5 * (times_shifted[i + 1] - times_shifted[i]) *
              (values_shifted(j, i + 1) + values_shifted(j, i));
        } else {
          const double s = (values_shifted(j, i + 1) - values_shifted(j, i)) /
                           (times_shifted[i + 1] - times_shifted[i]);

          result.fourier_coefficients_real(j, n) +=
              s * (std::cos(kn) - std::cos(ko));
          result.fourier_coefficients_imaginary(j, n) -=
              s * (std::sin(kn) - std::sin(ko));
        }
      }
    }

    if (n == 0) {
      for (unsigned int k = 0; k < result.n_components; k++) {
        result.fourier_coefficients_real(k, n) /= result.period;
      }
    } else {
      const double tmp_2 = (std::numbers::pi * std::numbers::pi * tmp * tmp);

      for (unsigned int k = 0; k < result.n_components; k++) {
        result.fourier_coefficients_real(k, n) =
            0.5 * result.fourier_coefficients_real(k, n) * result.period /
            tmp_2;
        result.fourier_coefficients_imaginary(k, n) =
            0.5 * result.fourier_coefficients_imaginary(k, n) * result.period /
            tmp_2;
      }
    }
  }

  return result;
}

FourierInterpolation FourierInterpolation::from_time_series_file(
    const std::string &file_name, unsigned int n_components, bool use_ramp) {
  std::ifstream file(file_name);
  svmp::throw_if<svmp::FileNotFoundException>(!file.is_open(), file_name);

  // Read the header of the file.
  const unsigned int n_time_points =
      read_int(file, file_name, "number of time points");
  const unsigned int n_fourier_coefficients =
      read_int(file, file_name, "number of Fourier coefficients");

  svmp::throw_if<svmp::FileFormatException>(
      n_time_points < 2, file_name,
      "At least 2 time points are required to set up Fourier "
      "interpolation. Only " +
          std::to_string(n_time_points) + " time points were provided.");

  svmp::throw_if<svmp::FileFormatException>(
      n_fourier_coefficients == 0, file_name,
      "At least 1 Fourier coefficient is required to set up Fourier "
      "interpolation. 0 Fourier coefficients were provided.");

  // Read the time-value pairs.
  Vector<double> times(n_time_points);
  Array<double> values(n_components, n_time_points);

  std::string line;
  int line_number = 1;
  unsigned int i = 0;

  while (std::getline(file, line)) {
    line = clean_line(line);

    if (!line.empty()) {
      const std::vector<double> line_values =
          string_to_vector_double(line, file_name, line_number);

      svmp::throw_if<svmp::FileFormatException>(
          line_values.size() != 1 + n_components, file_name,
          "Error reading values for the temporal values file for line " +
              std::to_string(line_number) + ": '" + line + "'; expected " +
              std::to_string(1 + n_components) + " values, but got " +
              std::to_string(line_values.size()) + ".");

      svmp::throw_if<svmp::FileFormatException>(
          i >= n_time_points, file_name,
          "Found more than the expected " + std::to_string(n_time_points) +
              " time points in the temporal values file.");

      times[i] = line_values[0];
      for (unsigned int j = 0; j < n_components; ++j) {
        values(j, i) = line_values[j + 1];
      }
      ++i;
    }

    ++line_number;
  }

  svmp::throw_if<svmp::FileFormatException>(
      i != n_time_points, file_name,
      "Expected " + std::to_string(n_time_points) +
          " time points in the temporal values file, but found " +
          std::to_string(i) + ".");

  return FourierInterpolation::from_time_series(n_fourier_coefficients, times,
                                                values, use_ramp);
}

FourierInterpolation FourierInterpolation::from_fourier_coefficients(
    const Vector<double> &linear_trend_initial_values,
    const Vector<double> &linear_trend_slopes,
    const Array<double> &fourier_coefficients_real,
    const Array<double> &fourier_coefficients_imaginary, double initial_time,
    double period) {
  // Linear trend initial values and slopes must have the same size (which will
  // determine the number of components).
  svmp::throw_if<svmp::FE::InvalidArgumentException>(
      linear_trend_initial_values.size() != linear_trend_slopes.size(),
      "linear_trend_initial_values and linear_trend_slopes must have the same "
      "size, but their sizes are " +
          std::to_string(linear_trend_initial_values.size()) + " and " +
          std::to_string(linear_trend_slopes.size()) + ".");

  // The number of rows of the Fourier coefficients must match the number of
  // components (i.e. the size of the linear trend vectors).
  svmp::throw_if<svmp::FE::InvalidArgumentException>(
      fourier_coefficients_real.nrows() != linear_trend_initial_values.size(),
      "The number of rows of fourier_coefficients_real must match the size of "
      "linear_trend_initial_values.");

  // The number of rows of the Fourier coefficients must match the number of
  // components (i.e. the size of the linear trend vectors).
  svmp::throw_if<svmp::FE::InvalidArgumentException>(
      fourier_coefficients_imaginary.nrows() !=
          linear_trend_initial_values.size(),
      "The number of rows of fourier_coefficients_imaginary must match the "
      "size of linear_trend_initial_values.");

  // The number of columns of the real and imaginary Fourier coefficients must
  // match.
  svmp::throw_if<svmp::FE::InvalidArgumentException>(
      fourier_coefficients_real.ncols() !=
          fourier_coefficients_imaginary.ncols(),
      "The number of columns of fourier_coefficients_real and "
      "fourier_coefficients_imaginary must match.");

  // The period must be a strictly positive number.
  svmp::throw_if<svmp::FE::InvalidArgumentException>(
      period <= 0.0, "The period must be a strictly positive number.");

  FourierInterpolation result;

  result.use_ramp = false;
  result.n_fourier_coefficients = fourier_coefficients_real.ncols();
  result.n_components = linear_trend_initial_values.size();
  result.linear_trend_initial_values = linear_trend_initial_values;
  result.linear_trend_slopes = linear_trend_slopes;
  result.fourier_coefficients_real = fourier_coefficients_real;
  result.fourier_coefficients_imaginary = fourier_coefficients_imaginary;
  result.initial_time = initial_time;
  result.period = period;

  return result;
}

FourierInterpolation FourierInterpolation::from_fourier_coefficients_file(
    const std::string &file_name, unsigned int n_components) {
  std::ifstream file(file_name);
  svmp::throw_if<svmp::FileNotFoundException>(!file.is_open(), file_name);

  const double initial_time = read_double(file, file_name, "initial time");
  const double period = read_double(file, file_name, "period");

  // Read the linear trend part. Each line is expected to contain the initial
  // value and slope of the linear trend for one component, so it must hold at
  // least 2 values.
  Vector<double> linear_trend_initial_values(n_components);
  Vector<double> linear_trend_slopes(n_components);

  std::string line;
  unsigned int current_component = 0;
  unsigned int line_number = 2;

  while (current_component < n_components && std::getline(file, line)) {
    line = clean_line(line);

    if (!line.empty()) {
      const std::vector<double> values =
          string_to_vector_double(line, file_name, current_component + 2);

      svmp::throw_if<svmp::FileFormatException>(
          values.size() != 2, file_name,
          "Error reading the linear trend for component " +
              std::to_string(current_component) + ": '" + line +
              "'; expected exactly 2 values (initial value and slope), but "
              "got " +
              std::to_string(values.size()) + ".");

      linear_trend_initial_values[current_component] = values[0];
      linear_trend_slopes[current_component] = values[1];
      ++current_component;
    }

    ++line_number;
  }

  svmp::throw_if<svmp::FileFormatException>(
      current_component != n_components, file_name,
      "Expected " + std::to_string(n_components) +
          " lines of linear trend coefficients (one per component), but only "
          "found " +
          std::to_string(current_component) + ".");

  // Read the Fourier coefficients.
  const unsigned int n_fourier_coefficients =
      read_int(file, file_name, "number of Fourier coefficients");

  Array<double> fourier_coefficients_real(n_components, n_fourier_coefficients);
  Array<double> fourier_coefficients_imaginary(n_components,
                                               n_fourier_coefficients);

  // Each line is expected to hold the real parts of the coefficients for all
  // components followed by the imaginary parts, so it must hold exactly
  // 2 * n_components values.
  const unsigned int n_values_per_line = 2 * n_components;

  unsigned int current_coefficient = 0;
  while (current_coefficient < n_fourier_coefficients &&
         std::getline(file, line)) {
    line = clean_line(line);

    if (!line.empty()) {
      const std::vector<double> values =
          string_to_vector_double(line, file_name, line_number);

      svmp::throw_if<svmp::FileFormatException>(
          values.size() != n_values_per_line, file_name,
          "Error reading Fourier coefficient " +
              std::to_string(current_coefficient) + ": '" + line +
              "'; expected " + std::to_string(n_values_per_line) +
              " values (real and imaginary parts for " +
              std::to_string(n_components) + " components), but got " +
              std::to_string(values.size()) + ".");

      for (unsigned int j = 0; j < n_components; ++j) {
        fourier_coefficients_real(j, current_coefficient) = values[j];
        fourier_coefficients_imaginary(j, current_coefficient) =
            values[j + n_components];
      }

      ++current_coefficient;
    }

    ++line_number;
  }

  svmp::throw_if<svmp::FileFormatException>(
      current_coefficient != n_fourier_coefficients, file_name,
      "Expected " + std::to_string(n_fourier_coefficients) +
          " lines of Fourier coefficients, but only found " +
          std::to_string(current_coefficient) + ".");

  return FourierInterpolation::from_fourier_coefficients(
      linear_trend_initial_values, linear_trend_slopes,
      fourier_coefficients_real, fourier_coefficients_imaginary, initial_time,
      period);
}

void FourierInterpolation::distribute(const CmMod &cm_mod, const cmType &cm) {
  // Only the master knows whether the object has been initialized. Therefore,
  // we broadcast the initialization flag to all ranks, so that the following if
  // statement can be run correctly by all.
  bool initialized = defined();
  cm.bcast(cm_mod, &initialized);

  if (initialized) {
    cm.bcast(cm_mod, &use_ramp);
    cm.bcast(cm_mod, &n_fourier_coefficients);
    cm.bcast(cm_mod, &n_components);

    // All ranks but the master need to allocate the arrays before receiving
    // data.
    if (cm.slv(cm_mod)) {
      linear_trend_initial_values.resize(n_components);
      linear_trend_slopes.resize(n_components);
      fourier_coefficients_real.resize(n_components, n_fourier_coefficients);
      fourier_coefficients_imaginary.resize(n_components,
                                            n_fourier_coefficients);
    }

    cm.bcast(cm_mod, &initial_time);
    cm.bcast(cm_mod, &period);
    cm.bcast(cm_mod, linear_trend_initial_values);
    cm.bcast(cm_mod, linear_trend_slopes);
    cm.bcast(cm_mod, fourier_coefficients_real);
    cm.bcast(cm_mod, fourier_coefficients_imaginary);
  }
}

Vector<double> FourierInterpolation::value(double time) const {
  Vector<double> result(n_components);
  static Vector<double> dummy;

  evaluate_internal(time, /* evaluate_derivative = */ false, result, dummy);

  return result;
}

std::pair<Vector<double>, Vector<double>>
FourierInterpolation::value_and_derivative(double time) const {
  Vector<double> value(n_components);
  Vector<double> derivative(n_components);

  evaluate_internal(time, /* evaluate_derivative = */ true, value, derivative);

  return std::make_pair(value, derivative);
}

bool FourierInterpolation::defined() const {
  return n_fourier_coefficients != 0;
}

unsigned int FourierInterpolation::get_n_components() const {
  return n_components;
}

unsigned int FourierInterpolation::get_n_fourier_coefficients() const {
  return n_fourier_coefficients;
}

double FourierInterpolation::get_linear_trend_initial_value(
    unsigned int component) const {

  svmp::throw_if<svmp::FE::NotInitializedException>(
      !defined(),
      "Cannot get linear trend initial value of FourierInterpolation "
      "instance that has not been defined.");

  svmp::throw_if<svmp::FE::InvalidArgumentException>(
      component >= n_components,
      "Component index " + std::to_string(component) +
          " is out of bounds for FourierInterpolation instance with " +
          std::to_string(n_components) + " components.");

  return linear_trend_initial_values[component];
}

double
FourierInterpolation::get_linear_trend_slope(unsigned int component) const {

  svmp::throw_if<svmp::FE::NotInitializedException>(
      !defined(), "Cannot get linear trend initial value of "
                  "FourierInterpolation instance that has not been defined.");

  svmp::throw_if<svmp::FE::InvalidArgumentException>(
      component >= n_components,
      "Component index " + std::to_string(component) +
          " is out of bounds for FourierInterpolation instance with " +
          std::to_string(n_components) + " components.");

  return linear_trend_slopes[component];
}

double
FourierInterpolation::get_coefficient_real(unsigned int component,
                                           unsigned int frequency) const {
  svmp::throw_if<svmp::FE::NotInitializedException>(
      !defined(), "Cannot get linear trend initial value of "
                  "FourierInterpolation instance that has not been defined.");

  svmp::throw_if<svmp::FE::InvalidArgumentException>(
      component >= n_components,
      "Component index " + std::to_string(component) +
          " is out of bounds for FourierInterpolation instance with " +
          std::to_string(n_components) + " components.");

  svmp::throw_if<svmp::FE::InvalidArgumentException>(
      frequency >= n_fourier_coefficients,
      "Frequency index " + std::to_string(frequency) +
          " is out of bounds for FourierInterpolation instance with " +
          std::to_string(n_fourier_coefficients) + " Fourier coefficients.");

  return fourier_coefficients_real(component, frequency);
}

double
FourierInterpolation::get_coefficient_imaginary(unsigned int component,
                                                unsigned int frequency) const {
  svmp::throw_if<svmp::FE::NotInitializedException>(
      !defined(), "Cannot get linear trend initial value of "
                  "FourierInterpolation instance that has not been defined.");

  svmp::throw_if<svmp::FE::InvalidArgumentException>(
      component >= n_components,
      "Component index " + std::to_string(component) +
          " is out of bounds for FourierInterpolation instance with " +
          std::to_string(n_components) + " components.");

  svmp::throw_if<svmp::FE::InvalidArgumentException>(
      frequency >= n_fourier_coefficients,
      "Frequency index " + std::to_string(frequency) +
          " is out of bounds for FourierInterpolation instance with " +
          std::to_string(n_fourier_coefficients) + " Fourier coefficients.");

  return fourier_coefficients_imaginary(component, frequency);
}

void FourierInterpolation::evaluate_internal(double time,
                                             bool evaluate_derivative,
                                             Vector<double> &value,
                                             Vector<double> &derivative) const {
  svmp::throw_if<svmp::FE::NotInitializedException>(
      !defined(), "Cannot evaluate FourierInterpolation instance that has not "
                  "been defined.");

  // Shifted and rescaled time.
  // The input time is shifted by ti. Then, if using the ramp function, it is
  // clamped to the interval [0, T]. Otherwise, the time is wrapped to the
  // interval [0, T], to enable periodicity.
  const double t = use_ramp
                       ? std::max(std::min(time - initial_time, period), 0.0)
                       : std::fmod(time - initial_time, period);

  // Linear trend.
  for (unsigned int i = 0; i < n_components; ++i) {
    value[i] = linear_trend_initial_values[i] + t * linear_trend_slopes[i];

    if (evaluate_derivative) {
      if (use_ramp && (time < initial_time || time > initial_time + period)) {
        derivative[i] = 0.0;
      } else {
        derivative[i] = linear_trend_slopes[i];
      }
    }
  }

  // Fourier series.
  if (!use_ramp) {
    const double tmp = 2.0 * std::numbers::pi / period;

    // Fourier series.
    for (int i = 0; i < n_fourier_coefficients; ++i) {
      const double dk = tmp * i;
      const double K = t * dk;

      for (int j = 0; j < n_components; ++j) {
        // Using value[j] = value[j] + ... instead of value[j] += ..., because
        // the latter changes the order of operations enough to break some of
        // the tests.
        // @todo[michelebucelli] This seems pretty fragile! It happens in a few
        //   other places in this file as well, where using an increment
        //   operator (*= or +=) changes the order of operations resulting in
        //   changes in the test values. This should be investigated and the
        //   best (i.e. most accurate) operation order should be chosen.
        value[j] = value[j] + fourier_coefficients_real(j, i) * std::cos(K) -
                   fourier_coefficients_imaginary(j, i) * std::sin(K);

        if (evaluate_derivative) {
          derivative[j] -=
              (fourier_coefficients_real(j, i) * std::sin(K) +
               fourier_coefficients_imaginary(j, i) * std::cos(K)) *
              dk;
        }
      }
    }
  }
}