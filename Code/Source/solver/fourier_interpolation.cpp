// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "fourier_interpolation.h"
#include "fft.h"

#include <fstream>
#include <regex>
#include <sstream>

namespace {
/// Clean a line by removing carriage returns and leading/trailing spaces, and
/// replacing multiple spaces with a single space.
std::string clean_line(std::string line) {
  line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
  return std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1");
}
} // namespace

FourierInterpolation FourierInterpolation::from_time_series(
    const unsigned int n_fourier_coefficients,
    const std::vector<std::vector<double>> &temporal_values,
    const bool use_ramp) {
  const unsigned int n_time_points = temporal_values.size();

  // @todo[michelebucelli] This should be an appropriate exception.
  if (n_time_points < 2) {
    throw std::runtime_error("At least two time points are needed to construct "
                             "a FourierInterpolation object.");
  }

  const unsigned int n_components = temporal_values[0].size() - 1;

  // Check that all entries of temporal_values have the same number of elements.
  // @todo[michelebucelli] This should probably be an Array, in which this is
  //   enforced a priori, rather than a vector of vectors. This change would
  //   also save us the extraction for loop below.
  for (const auto &row : temporal_values) {
    if (row.size() != n_components + 1) {
      throw std::runtime_error(
          "All rows of temporal_values must have the same number of elements.");
    }
  }

  Vector<double> times(n_time_points);
  Array<double> values(n_components, n_time_points);

  for (unsigned int i = 0; i < n_time_points; ++i) {
    times[i] = temporal_values[i][0];
    for (unsigned int j = 0; j < n_components; ++j) {
      values(j, i) = temporal_values[i][j + 1];
    }
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
  for (unsigned int i = 0; i < n_time_points; ++i) {
    times[i] -= result.initial_time;
    for (unsigned int j = 0; j < n_components; ++j) {
      values(j, i) = values(j, i) - result.linear_trend_initial_values[j] -
                     result.linear_trend_slopes[j] * times[i];
    }
  }

  // Compute the Fourier coefficients.
  for (int n = 0; n < result.n_fourier_coefficients; ++n) {
    const double tmp = static_cast<double>(n);
    result.fourier_coefficients_real.set_col(n, 0.0);
    result.fourier_coefficients_imaginary.set_col(n, 0.0);

    for (int i = 0; i < n_time_points - 1; ++i) {
      const double ko = 2.0 * consts::pi * tmp * times[i] / result.period;
      const double kn = 2.0 * consts::pi * tmp * times[i + 1] / result.period;

      for (int j = 0; j < result.n_components; j++) {
        const double s =
            (values(j, i + 1) - values(j, i)) / (times[i + 1] - times[i]);

        if (n == 0) {
          result.fourier_coefficients_real(j, n) +=
              0.5 * (times[i + 1] - times[i]) *
              (values(j, i + 1) + values(j, i));
        } else {
          result.fourier_coefficients_real(j, n) +=
              s * (std::cos(kn) - std::cos(ko));
          result.fourier_coefficients_imaginary(j, n) -=
              s * (std::sin(kn) - std::sin(ko));
        }
      }
    }

    if (n == 0) {
      for (int k = 0; k < result.n_components; k++) {
        result.fourier_coefficients_real(k, n) /= result.period;
      }
    } else {
      const double tmp_2 = (consts::pi * consts::pi * tmp * tmp);

      for (int k = 0; k < result.n_components; k++) {
        // @todo[michelebucelli] Check if operator*= breaks something here.
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

FourierInterpolation
FourierInterpolation::from_time_series_file(const std::string &file_name,
                                            const unsigned int n_components,
                                            const bool use_ramp) {
  std::ifstream file(file_name);

  // @todo[michelebucelli] This should actually thrown an exception, ideally of
  //   a dedicated type such as FileNotFoundException (to be defined).
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + file_name);
  }

  // Read the header of the file.
  int n_time_points, n_fourier_coefficients;
  file >> n_time_points >> n_fourier_coefficients;

  // @todo[michelebucelli] This should also be an appropriate exception.
  if (n_time_points < 2 || n_fourier_coefficients == 0) {
    throw std::runtime_error(
        "Error reading the first line of the temporal values file '" +
        file_name + "'.");
  }

  // Read the time-value pairs.
  std::vector<std::vector<double>> values;
  double tmp;

  std::string line;
  int line_number = 1;

  while (std::getline(file, line)) {
    line = clean_line(line);
    if (line.empty())
      continue;

    std::istringstream line_string_stream(line);
    std::vector<double> line_values;

    while (!line_string_stream.eof()) {
      line_string_stream >> tmp;

      // @todo[michelebucelli] This should also be an appropriate exception.
      if (line_string_stream.fail()) {
        throw std::runtime_error(
            "Error reading values for the temporal values file '" + file_name +
            "' for line " + std::to_string(line_number) + ": '" + line +
            "'; value number " + std::to_string(line_values.size() + 1) +
            " is not a double.");
      }

      line_values.push_back(tmp);
    }

    // @todo[michelebucelli] This should also be an appropriate exception.
    if (line_values.size() != 1 + n_components) {
      throw std::runtime_error(
          "Error reading values for the temporal values file '" + file_name +
          "' for line " + std::to_string(line_number) + ": '" + line +
          "'; expected " + std::to_string(1 + n_components) +
          " values, but got " + std::to_string(line_values.size()) + ".");
    }

    values.push_back(line_values);
    ++line_number;
  }

  return FourierInterpolation::from_time_series(n_fourier_coefficients, values,
                                                use_ramp);
}

FourierInterpolation FourierInterpolation::from_fourier_coefficients(
    const Vector<double> &linear_trend_initial_values,
    const Vector<double> &linear_trend_slopes,
    const Array<double> &fourier_coefficients_real,
    const Array<double> &fourier_coefficients_imaginary,
    const double initial_time, const double period) {
  // @todo[michelebucelli] Add some correctness checks on the input arguments.

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
    const std::string &file_name, const unsigned int n_components) {
  std::ifstream file(file_name);

  // @todo[michelebucelli] This should actually thrown an exception, ideally of
  //   a dedicated type such as FileNotFoundException (to be defined).
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + file_name);
  }

  double initial_time, period;
  file >> initial_time >> period;

  // Read the linear trend part.
  Vector<double> linear_trend_initial_values(n_components);
  Vector<double> linear_trend_slopes(n_components);

  std::string line;
  double tmp;
  unsigned int current_component = 0;

  while (current_component < n_components && std::getline(file, line)) {
    line = clean_line(line);
    if (line.empty())
      continue;

    std::istringstream line_string_stream(line);
    std::vector<double> values;

    while (line_string_stream >> tmp) {
      values.push_back(tmp);
    }

    linear_trend_initial_values[current_component] = values[0];
    linear_trend_slopes[current_component] = values[1];
    ++current_component;
  }

  // Read the Fourier coefficients.
  unsigned int n_fourier_coefficients;
  file >> n_fourier_coefficients;

  Array<double> fourier_coefficients_real(n_components, n_fourier_coefficients);
  Array<double> fourier_coefficients_imaginary(n_components,
                                               n_fourier_coefficients);

  unsigned int current_coefficient = 0;
  while (current_coefficient < n_fourier_coefficients &&
         std::getline(file, line)) {
    line = clean_line(line);
    if (line.empty())
      continue;

    std::istringstream line_string_stream(line);
    std::vector<double> values;

    while (line_string_stream >> tmp) {
      values.push_back(tmp);
    }

    for (unsigned int j = 0; j < n_components; ++j) {
      fourier_coefficients_real(j, current_coefficient) = values[j];
      fourier_coefficients_imaginary(j, current_coefficient) =
          values[j + n_components];
    }

    ++current_coefficient;
  }

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

Vector<double> FourierInterpolation::value(const double time) const {
  Vector<double> result(n_components);
  static Vector<double> dummy;

  evaluate_internal(time, /* evaluate_derivative = */ false, result, dummy);

  return result;
}

std::pair<Vector<double>, Vector<double>>
FourierInterpolation::value_and_derivative(const double time) const {
  Vector<double> value(n_components);
  Vector<double> derivative(n_components);

  evaluate_internal(time, /* evaluate_derivative = */ true, value, derivative);

  return std::make_pair(value, derivative);
}

void FourierInterpolation::evaluate_internal(const double time,
                                             const bool evaluate_derivative,
                                             Vector<double> &value,
                                             Vector<double> &derivative) const {
  // @todo[michelebucelli] This should be an appropriate exception.
  if (!defined()) {
    throw std::runtime_error("Cannot evaluate FourierInterpolation instance "
                             "that has not been defined.");
  }

  // Shifted and rescaled time.
  // The input time is shifted by ti. Then, if using the ramp function, it is
  // clamped to the interval [0, T]. Otherwise, the time is wrapped to the
  // interval [0, T], to enable periodicity.
  const double t = use_ramp
                       ? std::max(std::min(time - initial_time, period), 0.0)
                       : std::fmod(time - initial_time, period);

  // Linear trend.
  for (int i = 0; i < n_components; ++i) {
    value[i] = linear_trend_initial_values[i] + t * linear_trend_slopes[i];

    if (evaluate_derivative)
      derivative[i] = linear_trend_slopes[i];
  }

  // Fourier series.
  if (!use_ramp) {
    const double tmp = 2.0 * consts::pi / period;

    // Fourier series.
    for (int i = 0; i < n_fourier_coefficients; ++i) {
      const double dk = tmp * i;
      const double K = t * dk;

      for (int j = 0; j < n_components; ++j) {
        // Using value[j] = value[j] + ... instead of value[j] += ..., because
        // the latter changes the order of operations enough to break some of
        // the tests.
        // @todo[michelebucelli] This seems pretty fragile!
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