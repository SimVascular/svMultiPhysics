// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "fcType.h"
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

fcType fcType::from_time_series_file(const std::string &file_name,
                                     const unsigned int n_dimensions,
                                     const bool is_ramp) {
  fcType result;
  result.d = n_dimensions;
  result.lrmp = is_ramp;

  std::ifstream file(file_name);

  // @todo[michelebucelli] This should actually thrown an exception, ideally of
  //   a dedicated type such as FileNotFoundException (to be defined).
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + file_name);
  }

  // Read the header of the file.
  int n_time_points;
  file >> n_time_points >> result.n;

  // @todo[michelebucelli] This should also be an appropriate exception.
  if (n_time_points < 2 || result.n == 0) {
    throw std::runtime_error(
        "Error reading the first line of the temporal values file '" +
        file_name + "'.");
  }

  if (is_ramp)
    result.n = 1;

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
    if (line_values.size() != 1 + n_dimensions) {
      throw std::runtime_error(
          "Error reading values for the temporal values file '" + file_name +
          "' for line " + std::to_string(line_number) + ": '" + line +
          "'; expected " + std::to_string(1 + n_dimensions) +
          " values, but got " + std::to_string(line_values.size()) + ".");
    }

    values.push_back(line_values);
    ++line_number;
  }

  result.qi.resize(n_dimensions);
  result.qs.resize(n_dimensions);
  result.r.resize(n_dimensions, result.n);
  result.i.resize(n_dimensions, result.n);

  fft(n_time_points, values, result);

  return result;
}

fcType fcType::from_fourier_coefficients_file(const std::string &file_name,
                                              const unsigned int n_dimensions) {
  fcType result;
  result.d = n_dimensions;
  result.lrmp = false;

  std::ifstream file(file_name);

  // @todo[michelebucelli] This should actually thrown an exception, ideally of
  //   a dedicated type such as FileNotFoundException (to be defined).
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + file_name);
  }

  // Read the initial time and period from the header.
  file >> result.ti >> result.T;

  // Read the linear trend part.
  result.qi.resize(n_dimensions);
  result.qs.resize(n_dimensions);

  std::string line;
  double tmp;
  unsigned int current_dimension = 0;

  while (current_dimension < n_dimensions && std::getline(file, line)) {
    line = clean_line(line);
    if (line.empty())
      continue;

    std::istringstream line_string_stream(line);
    std::vector<double> values;

    while (line_string_stream >> tmp) {
      values.push_back(tmp);
    }

    result.qi[current_dimension] = values[0];
    result.qs[current_dimension] = values[1];
    ++current_dimension;
  }

  // Read the Fourier coefficients.
  file >> result.n;
  result.r.resize(n_dimensions, result.n);
  result.i.resize(n_dimensions, result.n);

  unsigned int current_coefficient = 0;
  while (current_coefficient < result.n && std::getline(file, line)) {
    line = clean_line(line);
    if (line.empty())
      continue;

    std::istringstream line_string_stream(line);
    std::vector<double> values;

    while (line_string_stream >> tmp) {
      values.push_back(tmp);
    }

    for (unsigned int i = 0; i < n_dimensions; ++i) {
      result.r(i, current_coefficient) = values[i];
      result.i(i, current_coefficient) = values[i + n_dimensions];
    }

    ++current_coefficient;
  }

  return result;
}

Vector<double> fcType::value(const double time) const {
  Vector<double> result(d);
  static Vector<double> dummy;

  evaluate_internal(time, /* evaluate_derivative = */ false, result, dummy);

  return result;
}

std::pair<Vector<double>, Vector<double>>
fcType::value_and_derivative(const double time) const {
  Vector<double> value(d), derivative(d);

  evaluate_internal(time, /* evaluate_derivative = */ true, value, derivative);

  return std::make_pair(value, derivative);
}

void fcType::evaluate_internal(const double time,
                               const bool evaluate_derivative,
                               Vector<double> &value,
                               Vector<double> &derivative) const {
  // @todo[michelebucelli] This should be an appropriate exception.
  if (!defined()) {
    throw std::runtime_error(
        "Cannot evaluate fcType instance that has not been defined.");
  }

  // Shifted and rescaled time.
  // The input time is shifted by ti. Then, if using the ramp function, it is
  // clamped to the interval [0, T]. Otherwise, the time is wrapped to the
  // interval [0, T], to enable periodicity.
  const double t =
      lrmp ? std::max(std::min(time - ti, T), 0.0) : std::fmod(time - ti, T);

  // Linear trend.
  for (int i = 0; i < d; ++i) {
    value[i] = qi[i] + t * qs[i];

    if (evaluate_derivative)
      derivative[i] = qs[i];
  }

  // Fourier series.
  if (!lrmp) {
    const double tmp = 2.0 * consts::pi / T;

    // Fourier series.
    for (int i = 0; i < n; ++i) {
      const double dk = tmp * i;
      const double K = t * dk;

      for (int j = 0; j < d; ++j) {
        // Using value[j] = value[j] + ... instead of value[j] += ..., because
        // the latter changes the order of operations enough to break some of
        // the tests.
        // @todo[michelebucelli] This seems pretty fragile!
        value[j] =
            value[j] + r(j, i) * std::cos(K) - this->i(j, i) * std::sin(K);

        if (evaluate_derivative) {
          derivative[j] -=
              (r(j, i) * std::sin(K) + this->i(j, i) * std::cos(K)) * dk;
        }
      }
    }
  }
}