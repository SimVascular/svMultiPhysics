// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef FC_TYPE_H
#define FC_TYPE_H

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
class fcType {
public:
  bool defined() const { return n != 0; };

  /// @brief Read Fourier coefficients from file and return the corresponding
  /// fcType instance.
  ///
  /// @todo[michelebucelli] More detailed documentation, especially on the file
  ///   format.
  static fcType from_fourier_coefficients_file(const std::string &file_name,
                                               const unsigned int n_dimensions);

  /// @brief Read a time series from file and return the corresponding fcType
  /// instance.
  ///
  /// @todo[michelebucelli] More detailed documentation, especially on the file
  ///   format.
  static fcType from_time_series_file(const std::string &file_name,
                                      const unsigned int n_dimensions,
                                      const bool is_ramp);

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

  /// Toggle whether this is a ramp function or not.
  bool lrmp = false;

  /// Number of Fourier coefficients.
  int n = 0;

  /// Numbero of dimensions (scalar or vector).
  int d = 0;

  /// Initial value.
  Vector<double> qi;

  /// Time derivative of linear part.
  Vector<double> qs;

  /// Period.
  double T = 0.0;

  /// Initial time.
  double ti = 0.0;

  /// Imaginary part of coefficients.
  Array<double> i;

  /// Real part of coefficients.
  Array<double> r;

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
   * @throws std::runtime_error if this fcType instance has not been
   * initialized (i.e. if @ref defined returns false).
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
};

#endif