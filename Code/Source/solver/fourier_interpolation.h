// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef FOURIER_INTERPOLATION_H
#define FOURIER_INTERPOLATION_H

#include "Array.h"
#include "CmMod.h"
#include "Vector.h"

#include <string>
#include <utility>
#include <vector>

/// @brief Fourier coefficients that are used to specify time-dependent
/// functions.
///
/// These are used for boundary conditions and for time-dependent prescribed
/// active stress.
///
/// @todo[michelebucelli] More detailed documentation.
class FourierInterpolation {
public:
  /// @brief Construct an FourierInterpolation from a time series.
  ///
  /// @todo[michelebucelli] More detailed documentation.
  static FourierInterpolation
  from_time_series(const unsigned int n_fourier_coefficients,
                   const std::vector<std::vector<double>> &temporal_values,
                   const bool use_ramp);

  /// @brief Read a time series from file and return the corresponding
  /// FourierInterpolation instance.
  ///
  /// @todo[michelebucelli] More detailed documentation, especially on the file
  ///   format.
  static FourierInterpolation
  from_time_series_file(const std::string &file_name,
                        const unsigned int n_components, const bool use_ramp);

  /// @brief Construct an FourierInterpolation instance from Fourier
  /// coefficients.
  static FourierInterpolation
  from_fourier_coefficients(const Vector<double> &linear_trend_initial_values,
                            const Vector<double> &linear_trend_slopes,
                            const Array<double> &fourier_coefficients_real,
                            const Array<double> &fourier_coefficients_imaginary,
                            const double initial_time, const double period);

  /// @brief Read Fourier coefficients from file and return the corresponding
  /// FourierInterpolation instance.
  ///
  /// @todo[michelebucelli] More detailed documentation, especially on the file
  ///   format.
  static FourierInterpolation
  from_fourier_coefficients_file(const std::string &file_name,
                                 const unsigned int n_components);

  /// @brief Distribute the data to all parallel processes.
  ///
  /// Does nothing if the object has not been initialized (i.e. if @ref defined
  /// returns false) on the master rank.
  void distribute(const CmMod &cm_mod, const cmType &cm);

  /// @brief Return the interpolated value.
  Vector<double> value(const double time) const;

  /// @brief Return the interpolated value and time derivative.
  std::pair<Vector<double>, Vector<double>>
  value_and_derivative(const double time) const;

  /// @name Data member access.
  /// @{

  /// @brief Return whether this object has been initialized.
  bool defined() const { return n_fourier_coefficients != 0; };

  /// @brief Get the dimension of the data interpolated by this object.
  unsigned int get_n_components() const { return n_components; }

  /// @brief Get the initial value of the linear trend part for one component.
  const double
  get_linear_trend_initial_value(const unsigned int component) const {
    // @todo[michelebucelli] Add check on the index.
    return linear_trend_initial_values[component];
  }

  /// @brief Get the slope of the linear trend part for one component.
  const double get_linear_trend_slope(const unsigned int component) const {
    // @todo[michelebucelli] Add check on the index.
    return linear_trend_slopes[component];
  }

  /// @brief Get the real part of the Fourier coefficients for one component.
  const double get_coefficient_real(const unsigned int component,
                                    const unsigned int frequency) const {
    // @todo[michelebucelli] Add check on the indices.
    return fourier_coefficients_real(component, frequency);
  }

  /// @brief Get the imaginary part of the Fourier coefficients for one
  /// component.
  const double get_coefficient_imaginary(const unsigned int component,
                                         const unsigned int frequency) const {
    // @todo[michelebucelli] Add check on the indices.
    return fourier_coefficients_imaginary(component, frequency);
  }

  /// @}

private:
  /** @brief Internal evaluation function.
   *
   * Uses the inverse Fourier transform to evaluate the value, and optionally
   * the derivative, of the interpolated data. This function is not meant to be
   * used directly, but only as a backend to @ref value and @ref
   * value_and_derivative.
   *
   * The vectors value and derivative are assumed to be of size @ref d. This is
   * not checked by this function.
   *
   * @throws std::runtime_error if this FourierInterpolation instance has not
   * been initialized (i.e. if @ref defined returns false).
   *
   * @param[in] time The time at which to evaluate the interpolation.
   * @param[in] evaluate_derivative Whether to also evaluate the time
   *   derivative of the interpolation.
   * @param[out] value The interpolated value at the given time.
   * @param[out] derivative The time derivative of the interpolated value at
   *   the given time. If evaluated_derivative is false, this will not be
   *   modified or accessed.
   */
  void evaluate_internal(const double time, const bool evaluate_derivative,
                         Vector<double> &value,
                         Vector<double> &derivative) const;

  /**
   * @brief Toggle whether this is a ramp function or not.
   *
   * See the general class documentation for details on the difference between
   * ramp and periodic interpolation.
   */
  bool use_ramp = false;

  /**
   * @brief Number of Fourier coefficients.
   */
  unsigned int n_fourier_coefficients = 0;

  /**
   * @brief Number of components of the interpolated data.
   */
  unsigned int n_components = 0;

  /**
   * @brief Initial value for the linear trend.
   *
   * This is a vector with n_components entries, where each entry is the initial
   * value (i.e. the value for time = ti) of the linear trend part of the
   * interpolated data for that component.
   */
  Vector<double> linear_trend_initial_values;

  /**
   * @brief Time derivative for the linear trend.
   *
   * This is a vector with n_components entries, where each entry is the slope
   * of the linear trend part of the interpolated data for that component.
   */
  Vector<double> linear_trend_slopes;

  /**
   * @brief Period of the interpolated data.
   *
   * This is disregarded if use_ramp is true. See the general class
   * documentation for details on the difference between ramp and periodic
   * interpolation.
   */
  double period = 0.0;

  /**
   * @brief Initial time.
   */
  double initial_time = 0.0;

  /**
   * @brief Real part of the Fourier series coefficients.
   *
   * This is a 2D array with n_components rows and n_fourier_coefficients
   * columns.
   */
  Array<double> fourier_coefficients_real;

  /**
   * @brief Imaginary part of the Fourier series coefficients.
   *
   * This is a 2D array with n_components rows and n_fourier_coefficients
   * columns.
   */
  Array<double> fourier_coefficients_imaginary;
};

#endif