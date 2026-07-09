// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef FOURIER_INTERPOLATION_H
#define FOURIER_INTERPOLATION_H

#include "Array.h"
#include "CmMod.h"
#include "Vector.h"

#include "Core/Exception.h"
#include "FE/Common/FEException.h"

#include <string>
#include <utility>
#include <vector>

/**
 * @brief Fourier interpolation of time dependent data.
 *
 * This class implements interpolation of time dependent data through either a
 * Fourier series or a linear clamped ramp. Which of the two is used is
 * determined by the boolean member variable @ref use_ramp.
 *
 * In what follows, let @f$\{t_i, \mathbf{v}_i\}_{i=0}^{N-1}@f$ be the
 * interpolated time series, and @f$T = t_{N-1} - t_0@f$ be the period of the
 * time series. The time series is assumed to be ordered, i.e. @f$t_i <
 * t_{i+1}@f$ for all @f$i@f$.
 *
 * ## Using this class
 *
 * The FourierInterpolation class is not meant to be constructed directly, but
 * rather through one of the static methods @ref from_time_series, @ref
 * from_time_series_file, @ref from_fourier_coefficients, or @ref
 * from_fourier_coefficients_file.
 *
 * If data is only read from the master rank in a parallel setting, the method
 * @ref distribute can be used to broadcast the data to all ranks.
 *
 * Once the object has been constructed, the interpolated value at a given time
 * can be obtained through the method @ref value, and the interpolated value and
 * its time derivative can be obtained through the method @ref
 * value_and_derivative.
 *
 * ## Fourier series interpolation
 *
 * This gives rise to a periodic interpolation, with period @f$T@f$. Let @f$M@f$
 * be the user-defined number of Fourier modes. The interpolated value at time
 * @f$t@f$ is given by
 * @f[ \begin{aligned}
 * \tilde{\mathbf{v}}(t) &=
 *   \underbrace{\mathbf{q}_i + \mathbf{q}_s \tau(t)}_{\text{linear trend}} +
 *   \underbrace{
 *     \sum_{k=0}^{M-1} \left(
 *       \mathbf{c}_k^{\text{real}} \cos\left(\frac{2 \pi k \tau(t)}{T}\right) -
 *       \mathbf{c}_k^{\text{imag}} \sin\left(\frac{2 \pi k \tau(t)}{T}\right)
 *     \right)
 *   }_{\text{Fourier series}} \\
 *   \tau(t) &= (t - t_0)\;\text{mod}\;T
 * \end{aligned} @f]
 * The coefficients of the linear trend and Fourier series are determined from
 * the input time series as follows:
 * @f[ \begin{aligned}
 * \mathbf{q}_i &= \mathbf{v}_0 \\
 * \mathbf{q}_s &= \frac{\mathbf{v}_{N-1} - \mathbf{v}_0}{T} \\
 * \mathbf{c}_0^{\text{real}} &=
 *   \frac{1}{2T} \sum_{i=0}^{N-2} (\hat{t}_{i+1} - \hat{t}_i)
 *     (\hat{\mathbf{v}}_{i+1} + \hat{\mathbf{v}}_i) \\
 * \mathbf{c}_0^{\text{imag}} &= 0 \\
 * \mathbf{c}_k^{\text{real}} &= \frac{T}{2 \pi^2 k^2}
 *   \sum_{i=0}^{N-2} \frac{\hat{\mathbf{v}}_{i+1} - \hat{\mathbf{v}}_i}
 *   {\hat{t}_{i+1} - \hat{t}_i} \left(
 *     \cos\left(\frac{2 \pi k \hat{t}_{i+1}}{T}\right) -
 *     \cos\left(\frac{2 \pi k \hat{t}_i}{T}\right)
 *   \right) \\
 * \mathbf{c}_k^{\text{imag}} &= \frac{T}{2 \pi^2 k^2}
 *   -\sum_{i=0}^{N-2} \frac{\hat{\mathbf{v}}_{i+1} - \hat{\mathbf{v}}_i}
 *   {\hat{t}_{i+1} - \hat{t}_i} \left(
 *     \sin\left(\frac{2 \pi k \hat{t}_{i+1}}{T}\right) -
 *     \sin\left(\frac{2 \pi k \hat{t}_i}{T}\right)
 *   \right)
 * \end{aligned} @f]
 * The quantities @f$\hat{t}_i@f$ and @f$\hat{\mathbf{v}}_i@f$ are the time and
 * value series after subtracting the linear trend:
 * @f[ \begin{aligned}
 * \hat{t}_i &= t_i - t_0 \\
 * \hat{\mathbf{v}}_i &= \mathbf{v}_i - \mathbf{q}_i - \mathbf{q}_s \hat{t}_i
 * \end{aligned} @f]
 *
 * ## Ramp interpolation
 *
 * The interpolated value is equal to @f$\mathbf{v}_0@f$ until time @f$t_0@f$,
 * then it follows a linear ramp from @f$\mathbf{v}_0@f$ to
 * @f$\mathbf{v}_{N-1}@f$ at @f$t_{N-1}@f$, and then it remains constant at
 * @f$\mathbf{v}_{N-1}@f$ for all times greater than @f$t_{N-1}@f$. Notice in
 * particular that this interpolation is not periodic.
 *
 * The interpolated value is given by
 * @f[
 * \tilde{\mathbf{v}}(t) = \begin{cases}
 *   \mathbf{v}_0 & t < t_0 \\
 *   \mathbf{v}_0 + \frac{\mathbf{v}_{N-1}
 *     - \mathbf{v}_0}{t_{N-1} - t_0} (t - t_0)
 *     & t_0 \leq t < t_{N-1} \\
 *   \mathbf{v}_{N-1} & t \geq t_{N-1}
 * \end{cases}
 * @f]
 *
 * This is equivalent to only taking the linear trend part of the Fourier series
 * interpolation, and clamping the time to the interval @f$[t_0, t_{N-1}]@f$.
 *
 */
class FourierInterpolation {
public:
  /**
   * @brief Default constructor.
   *
   * This constructor is not meant to be used directly, and is only here to
   * facilitate storing objects of this type in STL containers. The object
   * constructed this way is not initialized and will not be usable until it is
   * assigned to a valid FourierInterpolation instance. Use the static methods
   * @ref from_time_series, @ref from_time_series_file, @ref
   * from_fourier_coefficients, or @ref from_fourier_coefficients_file to
   * construct a valid instance.
   */
  FourierInterpolation() = default;

  /**
   * @brief Construct a FourierInterpolation from a time series.
   *
   * @param[in] n_fourier_coefficients The number @f$M@f$ of Fourier modes to
   *   use in the interpolation.
   * @param[in] times The time points @f$t_i@f$ of the time series. It must be a
   *   vector in strictly ascending order, and an exception will be thrown
   *   otherwise.
   * @param[in] values The values @f$\mathbf{v}_i@f$ of the time series. It must
   *   be a 2D array with one row for each component and one column for each
   *   time point. The number of columns must match the size of @p times, and an
   *   exception will be thrown otherwise.
   * @param[in] use_ramp Whether to use a ramp function for the interpolation.
   *   See the general class documentation for the precise meaning of this
   *   choice.
   */
  static FourierInterpolation
  from_time_series(const unsigned int n_fourier_coefficients,
                   const Vector<double> &times, const Array<double> &values,
                   bool use_ramp);

  /**
   * @brief Read a time series from file and return the corresponding instance
   * of FourierInterpolation.
   *
   * The input file is expected to have the following format:
   * ```
   * <number of time points> <number of Fourier coefficients>
   * <time 0> <value 0, component 0> ... <value 0, component d-1>
   * <time 1> <value 1, component 0> ... <value 1, component d-1>
   * ...
   * <time N-1> <value N-1, component 0> ... <value N-1, component d-1>
   * ```
   *
   * @param[in] file_name The name of the file to read. An exception will be
   *   thrown if the file cannot be opened.
   * @param[in] n_components The number of components of the data to be
   *   interpolated. If any row in the file does not have this number of entries
   *   (plus one entry for the time), an exception will be thrown.
   * @param[in] use_ramp Whether to use a ramp function for the interpolation.
   *   See the general class documentation for the precise meaning of this
   *   choice. If this parameter is set to true, then the number of Fourier
   *   coefficients in the file will be ignored, and the time points except the
   *   first and last will have no effect.
   */
  static FourierInterpolation
  from_time_series_file(const std::string &file_name, unsigned int n_components,
                        bool use_ramp);

  /**
   * @brief Construct a FourierInterpolation from Fourier coefficients.
   *
   * This method bypasses the computation of the Fourier coefficients from a
   * time series, and construct the instance directly from precomputed
   * coefficients.
   *
   * This construction method always sets @ref use_ramp to false, and therefore
   * always constructs a periodic interpolation.
   *
   * @param[in] linear_trend_initial_values The initial values
   *   @f$\mathbf{q}_i@f$ of the linear trend part of the interpolation, with
   *   one entry for each component.
   * @param[in] linear_trend_slopes The slopes @f$\mathbf{q}_s@f$ of the linear
   *   trend part of the interpolation, with one entry for each component.
   * @param[in] fourier_coefficients_real The real part
   *   @f$\mathbf{c}_k^\text{real}@f$ of the Fourier interpolation. This is a 2D
   *   array, for which the first index selects the component and the second
   *   index selects the Fourier mode.
   * @param[in] fourier_coefficients_imaginary The imaginary part
   *   @f$\mathbf{c}_k^\text{imag}@f$ of the Fourier interpolation. This is a 2D
   *   array, for which the first index selects the component and the second
   *   index selects the Fourier mode.
   * @param[in] initial_time The initial time @f$t_0@f$ of the interpolation.
   * @param[in] period The period @f$T@f$ of the interpolation.
   */
  static FourierInterpolation
  from_fourier_coefficients(const Vector<double> &linear_trend_initial_values,
                            const Vector<double> &linear_trend_slopes,
                            const Array<double> &fourier_coefficients_real,
                            const Array<double> &fourier_coefficients_imaginary,
                            double initial_time, double period);

  /**
   * @brief Read Fourier coefficients from file and return the corresponding
   * FourierInterpolation instance.
   *
   * The input file is expected to have the following format:
   * ```
   * <initial time> <period>
   * <q_i, component 0> <q_s, component 0>
   * ...
   * <q_i, component d-1> <q_s, component d-1>
   * <number of Fourier coefficients>
   * <c_r[0][0]> <c_r[1][0]> ... <c_r[d-1][0]> <c_i[0][0]> <c_i[1][0]> ... <c_i[d-1][0]>
   * <c_r[0][1]> <c_r[1][1]> ... <c_r[d-1][1]> <c_i[0][1]> <c_i[1][1]> ... <c_i[d-1][1]>
   * ...
   * <c_r[0][M-1]> <c_r[1][M-1]> ... <c_r[d-1][M-1]> <c_i[0][M-1]> <c_i[1][M-1]> ... <c_i[d-1][M-1]>
   * ```
   * where <kbd>c_r</kbd> and <kbd>c_i</kbd> are the real and imaginary parts of
   * the Fourier coefficients.
   *
   * @param[in] file_name The name of the file to read. An exception will be
   *   thrown if the file cannot be opened.
   * @param[in] n_components The number of components of the data to be
   *   interpolated. This is used to check the correctness of the file format,
   *   and an exception will be thrown if the check fails.
   */
  static FourierInterpolation
  from_fourier_coefficients_file(const std::string &file_name,
                                 unsigned int n_components);

  /**
   * @brief Distribute the data to all parallel processes.
   *
   * Broadcasts the data contained in this object from the master rank to all
   * other ranks. This is necessary if only the master rank initializes the
   * object (e.g. by reading its data from a file), but all ranks need to use
   * it.
   *
   * @param[in] cm_mod The communication module to use for the broadcast.
   * @param[in] cm The communicator to use for the broadcast.
   */
  void distribute(const CmMod &cm_mod, const cmType &cm);

  /**
   * @brief Return the interpolated value at a given time.
   *
   * Refer to the general class documentation for details on how the returned
   * value is computed.
   *
   * @param[in] time The time at which to evaluate the interpolation.
   */
  Vector<double> value(double time) const;

  /**
   * @brief Return the interpolated value and its time derivative at a given
   * time.
   *
   * Refer to the general class documentation for details on how the returned
   * value is computed.
   *
   * @param[in] time The time at which to evaluate the interpolation.
   *
   * @return A pair in which the first element is the interpolated value and the
   *   second element is its time derivative.
   */
  std::pair<Vector<double>, Vector<double>>
  value_and_derivative(double time) const;

  /// @name Data member access.
  /// @{

  /// @brief Return whether this object has been initialized.
  bool defined() const;

  /// @brief Get the dimension of the data interpolated by this object.
  unsigned int get_n_components() const;

  /// @brief Get the number of Fourier coefficients used by this object.
  unsigned int get_n_fourier_coefficients() const;

  /// @brief Get the initial value of the linear trend part for one component.
  double get_linear_trend_initial_value(unsigned int component) const;

  /// @brief Get the slope of the linear trend part for one component.
  double get_linear_trend_slope(unsigned int component) const;

  /// @brief Get the real part of the Fourier coefficients for one component.
  double get_coefficient_real(unsigned int component,
                              unsigned int frequency) const;

  /// @brief Get the imaginary part of the Fourier coefficients for one
  /// component.
  double get_coefficient_imaginary(unsigned int component,
                                   unsigned int frequency) const;

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
   * @throws svmp::FE::NotInitializedException if this FourierInterpolation
   * instance has not been initialized (i.e. if @ref defined returns false).
   *
   * @param[in] time The time at which to evaluate the interpolation.
   * @param[in] evaluate_derivative Whether to also evaluate the time
   *   derivative of the interpolation.
   * @param[out] value The interpolated value at the given time.
   * @param[out] derivative The time derivative of the interpolated value at
   *   the given time. If evaluated_derivative is false, this will not be
   *   modified or accessed.
   */
  void evaluate_internal(double time, bool evaluate_derivative,
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