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

fcType fcType::from_time_series(
    const unsigned int n_fourier_coefficients,
    const std::vector<std::vector<double>> &temporal_values,
    const bool is_ramp) {
  const unsigned int n_time_points = temporal_values.size();

  // @todo[michelebucelli] This should be an appropriate exception.
  if (n_time_points < 2) {
    throw std::runtime_error(
        "At least two time points are needed to construct an fcType object.");
  }

  const unsigned int n_dimensions = temporal_values[0].size() - 1;

  // Check that all entries of temporal_values have the same number of elements.
  // @todo[michelebucelli] This should probably be an Array, in which this is
  //   enforced a priori, rather than a vector of vectors. This change would
  //   also save us the extraction for loop below.
  for (const auto &row : temporal_values) {
    if (row.size() != n_dimensions + 1) {
      throw std::runtime_error(
          "All rows of temporal_values must have the same number of elements.");
    }
  }

  Vector<double> times(n_time_points);
  Array<double> values(n_dimensions, n_time_points);

  for (unsigned int i = 0; i < n_time_points; ++i) {
    times[i] = temporal_values[i][0];
    for (unsigned int j = 0; j < n_dimensions; ++j) {
      values(j, i) = temporal_values[i][j + 1];
    }
  }

  fcType result;

  result.lrmp = is_ramp;
  result.d = n_dimensions;
  result.n = is_ramp ? 1 : n_fourier_coefficients;
  result.ti = times[0];
  result.T = times[n_time_points - 1] - times[0];
  result.qi.resize(n_dimensions);
  result.qs.resize(n_dimensions);
  result.r.resize(n_dimensions, result.n);
  result.i.resize(n_dimensions, result.n);

  // Compute the linear trend part.
  for (unsigned int j = 0; j < n_dimensions; ++j) {
    result.qi[j] = values(j, 0);
    result.qs[j] = (values(j, n_time_points - 1) - values(j, 0)) / result.T;
  }

  // Subtract the linear trend part from the values.
  for (unsigned int i = 0; i < n_time_points; ++i) {
    times[i] -= result.ti;
    for (unsigned int j = 0; j < n_dimensions; ++j) {
      values(j, i) = values(j, i) - result.qi[j] - result.qs[j] * times[i];
    }
  }

  // Compute the Fourier coefficients.
  for (int n = 0; n < result.n; ++n) {
    const double tmp = static_cast<double>(n);
    result.r.set_col(n, 0.0);
    result.i.set_col(n, 0.0);

    for (int i = 0; i < n_time_points - 1; ++i) {
      const double ko = 2.0 * consts::pi * tmp * times[i] / result.T;
      const double kn = 2.0 * consts::pi * tmp * times[i + 1] / result.T;

      for (int j = 0; j < result.d; j++) {
        const double s =
            (values(j, i + 1) - values(j, i)) / (times[i + 1] - times[i]);

        if (n == 0) {
          result.r(j, n) += 0.5 * (times[i + 1] - times[i]) *
                            (values(j, i + 1) + values(j, i));
        } else {
          result.r(j, n) += s * (std::cos(kn) - std::cos(ko));
          result.i(j, n) -= s * (std::sin(kn) - std::sin(ko));
        }
      }
    }

    if (n == 0) {
      for (int k = 0; k < result.d; k++) {
        result.r(k, n) /= result.T;
      }
    } else {
      for (int k = 0; k < result.d; k++) {
        result.r(k, n) = 0.5 * result.r(k, n) * result.T /
                         (consts::pi * consts::pi * tmp * tmp);
        result.i(k, n) = 0.5 * result.i(k, n) * result.T /
                         (consts::pi * consts::pi * tmp * tmp);
      }
    }
  }

  return result;
}

fcType fcType::from_time_series_file(const std::string &file_name,
                                     const unsigned int n_dimensions,
                                     const bool is_ramp) {
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

  return fcType::from_time_series(n_fourier_coefficients, values, is_ramp);
}

fcType fcType::from_fourier_coefficients(const Vector<double> &qi,
                                         const Vector<double> &qs,
                                         const Array<double> &r,
                                         const Array<double> &i,
                                         const double ti, const double T) {
  // @todo[michelebucelli] Add some correctness checks on the input arguments.

  fcType result;

  result.lrmp = false;
  result.n = r.ncols();
  result.d = qi.size();
  result.qi = qi;
  result.qs = qs;
  result.r = r;
  result.i = i;
  result.ti = ti;
  result.T = T;

  return result;
}

fcType fcType::from_fourier_coefficients_file(const std::string &file_name,
                                              const unsigned int n_dimensions) {
  std::ifstream file(file_name);

  // @todo[michelebucelli] This should actually thrown an exception, ideally of
  //   a dedicated type such as FileNotFoundException (to be defined).
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + file_name);
  }

  // Read the initial time and period from the header.
  double ti, T;
  file >> ti >> T;

  // Read the linear trend part.
  Vector<double> qi(n_dimensions), qs(n_dimensions);

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

    qi[current_dimension] = values[0];
    qs[current_dimension] = values[1];
    ++current_dimension;
  }

  // Read the Fourier coefficients.
  unsigned int n_fourier_coefficients;
  file >> n_fourier_coefficients;

  Array<double> r(n_dimensions, n_fourier_coefficients);
  Array<double> i(n_dimensions, n_fourier_coefficients);

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

    for (unsigned int j = 0; j < n_dimensions; ++j) {
      r(j, current_coefficient) = values[j];
      i(j, current_coefficient) = values[j + n_dimensions];
    }

    ++current_coefficient;
  }

  return fcType::from_fourier_coefficients(qi, qs, r, i, ti, T);
}

void fcType::distribute(const CmMod &cm_mod, const cmType &cm) {
  // Only the master knows whether the object has been initialized. Therefore,
  // we broadcast the initialization flag to all ranks, so that the following if
  // statement can be run correctly by all.
  bool initialized = defined();
  cm.bcast(cm_mod, &initialized);

  if (initialized) {
    cm.bcast(cm_mod, &lrmp);
    cm.bcast(cm_mod, &n);
    cm.bcast(cm_mod, &d);

    // All ranks but the master need to allocate the arrays before receiving
    // data.
    if (cm.slv(cm_mod)) {
      qi.resize(d);
      qs.resize(d);
      r.resize(d, n);
      i.resize(d, n);
    }

    cm.bcast(cm_mod, &ti);
    cm.bcast(cm_mod, &T);
    cm.bcast(cm_mod, qi);
    cm.bcast(cm_mod, qs);
    cm.bcast(cm_mod, r);
    cm.bcast(cm_mod, i);
  }
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