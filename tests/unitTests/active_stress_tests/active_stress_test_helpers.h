// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef ACTIVE_STRESS_TEST_HELPERS_H
#define ACTIVE_STRESS_TEST_HELPERS_H

#include "active_stress.h"
#include "Core/Exception.h"
#include "FE/Common/FEException.h"
#include "Vector.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

/** @brief Configure one standalone trusted-reference ActiveStress trajectory. */
struct ActiveStressTrajectoryConfiguration {
  /// End time of the simulated trajectory.
  double final_time = 0.0;

  /// Time-step size used for every trajectory update.
  double time_step = 0.0;

  /// Reference filename under UNIT_TEST_DATA_DIR.
  std::string reference_csv_filename;

  /**
   * @brief CSV column name for each model state, in model-state order.
   *
   * An empty vector selects @c s0, @c s1, and so on.
   */
  std::vector<std::string> state_reference_columns;

  /// Base tolerance in the pointwise comparison formula.
  double tolerance = 1.0e-10;
};

/**
 * @brief Run a twitch experiment on the ActiveStress model passed as template
 * argument and compare the simulated trajectory with a trusted external
 * reference solution loaded from a CSV file.
 *
 * The twitch experiment uses a prescribed calcium transient to trigger
 * contraction in the active stress model (see @ref calcium_at). Fiber stretch
 * is also prescribed (see @ref fiber_stretch_at), and fiber stretch rate is
 * obtained from the prescribed stretch trajectory.
 *
 * The active stress model is solved standalone, i.e. only one system of ODEs
 * is solved rather than a coupled three-dimensional mechanics problem.
 * Therefore, the test verifies the active stress ODE definition and time
 * stepping independently of the mechanics solver.
 *
 * ### Format of the reference solution
 *
 * The reference CSV contains a header followed by reference values at selected
 * simulation steps. The first two columns are @c step and @c Ta, containing the
 * simulation step and active tension, respectively. If the model has additional
 * internal state variables, their values follow in subsequent columns.
 *
 * The state variables are assumed by default to have the same order as in the
 * concrete model being tested. If necessary, @c state_reference_columns
 * supplies the CSV column name corresponding to each state variable, in the
 * concrete model's state order.
 *
 * Comparisons are performed only at the checkpoints listed in the reference
 * CSV; the file does not need to contain every simulated time step. At each
 * checkpoint, active tension and each state variable are compared pointwise
 * with the corresponding reference value. If a pointwise difference exceeds
 * the prescribed tolerance, GoogleTest records a failed expectation and the
 * test fails.
 *
 * Checkpoints are integer simulation-step indices and must coincide with
 * simulated time steps; intermediate checkpoints are not interpolated.
 * Checkpoint @c N is compared after the update starting at
 * @f$ t = N \Delta t @f$, which advances the interval
 * @f$[N\Delta t,(N+1)\Delta t]@f$.
 */
template <class ConcreteModel>
class ActiveStressTrajectoryTest {
  static_assert(std::is_base_of_v<ActiveStress, ConcreteModel>,
                "ConcreteModel must derive from ActiveStress");

  /// Test-only adapter that exposes the concrete model's protected node-local
  /// operations to the trajectory test.
  class TestAdapter : public ConcreteModel {
  public:
    using ConcreteModel::advance_time_step_local;
    using ConcreteModel::compute_active_tension_local;
    using ConcreteModel::init_local;
    using ConcreteModel::read_model_specific_parameters;
  };

public:
  /// Parameter type defined by the concrete active-stress model.
  using Parameters = typename ConcreteModel::Parameters;

  static_assert(std::is_base_of_v<ActiveStressModelParameters, Parameters>,
                "ConcreteModel::Parameters must derive from "
                "ActiveStressModelParameters");

  /**
   * @brief Configure a standalone trajectory test and load its reference CSV.
   *
   * @param parameters Model-specific parameters.
   * @param configuration Trajectory inputs and reference-data mapping.
   */
  ActiveStressTrajectoryTest(
      const Parameters &parameters,
      ActiveStressTrajectoryConfiguration configuration)
      : model_(make_configured_model(parameters)),
        configuration_(std::move(configuration)),
        n_steps_(number_of_steps(configuration_.final_time,
                                 configuration_.time_step))
  {
    validate_configuration();
    reference_path_ = std::string(UNIT_TEST_DATA_DIR) + "/" +
                      configuration_.reference_csv_filename;
    checkpoints_ = load_reference_data(
        reference_path_, model_->n_states, n_steps_,
        configuration_.state_reference_columns);
  }

  /** @brief Execute the common twitch protocol and compare every checkpoint. */
  void run()
  {
    Vector<double> state(model_->n_states);
    model_->init_local(state);

    size_t checkpoint_index = 0;
    for (int step = 0; step < n_steps_; ++step) {
      const double time = step * configuration_.time_step;
      const double calcium = calcium_at(time);
      const double fiber_stretch = fiber_stretch_at(time);
      const double fiber_stretch_rate =
          (fiber_stretch_at(time + configuration_.time_step) - fiber_stretch) /
          configuration_.time_step;

      model_->advance_time_step_local(time, configuration_.time_step, calcium,
                                      fiber_stretch, fiber_stretch_rate, state);
      const double active_tension =
          model_->compute_active_tension_local(state, fiber_stretch);

      compare_checkpoint_if_present(step, state, active_tension,
                                    checkpoint_index);
    }

    EXPECT_EQ(checkpoint_index, checkpoints_.size())
        << "Not all reference checkpoints were visited in " << reference_path_;
  }

private:
  /// Reference values for one simulation-step checkpoint.
  struct Checkpoint {
    /// Trajectory update index.
    int step;

    /// Expected active tension.
    double active_tension;

    /// Expected state values in concrete-model state order.
    std::vector<double> state;
  };

  /// Construct the test adapter and load its model-specific parameters.
  static std::unique_ptr<TestAdapter>
  make_configured_model(const Parameters &parameters)
  {
    auto model = std::make_unique<TestAdapter>();
    model->read_model_specific_parameters(parameters);
    return model;
  }

  /// Validate the trajectory inputs and resolve default state-column names.
  void validate_configuration()
  {
    svmp::throw_if<svmp::FE::InvalidArgumentException>(
        !(configuration_.tolerance > 0.0) ||
            !std::isfinite(configuration_.tolerance),
        "ActiveStress test tolerance must be finite and positive");
    svmp::throw_if<svmp::FE::InvalidArgumentException>(
        configuration_.reference_csv_filename.empty(),
        "ActiveStress trajectory reference CSV filename must not be empty");

    if (configuration_.state_reference_columns.empty()) {
      configuration_.state_reference_columns.reserve(model_->n_states);
      for (unsigned int state = 0; state < model_->n_states; ++state)
        configuration_.state_reference_columns.push_back(
            "s" + std::to_string(state));
    }

    svmp::throw_if<svmp::FE::InvalidArgumentException>(
        configuration_.state_reference_columns.size() != model_->n_states,
        "ActiveStress reference state-column mapping has the wrong size");

    for (size_t i = 0; i < configuration_.state_reference_columns.size(); ++i) {
      const auto &name = configuration_.state_reference_columns[i];
      svmp::throw_if<svmp::FE::InvalidArgumentException>(
          name.empty(),
          "ActiveStress reference state-column names must not be empty");
      svmp::throw_if<svmp::FE::InvalidArgumentException>(
          name == "step",
          "ActiveStress reference state-column mapping cannot use step");
      svmp::throw_if<svmp::FE::InvalidArgumentException>(
          std::find(configuration_.state_reference_columns.begin(),
                    configuration_.state_reference_columns.begin() + i,
                    name) !=
              configuration_.state_reference_columns.begin() + i,
          "ActiveStress reference state-column names must be unique");
    }
  }

  /// Validate the time configuration and return the number of simulation steps.
  static int number_of_steps(double final_time, double time_step)
  {
    svmp::throw_if<svmp::FE::InvalidArgumentException>(
        !(final_time > 0.0) || !std::isfinite(final_time),
        "ActiveStress test final time must be finite and positive");
    svmp::throw_if<svmp::FE::InvalidArgumentException>(
        !(time_step > 0.0) || !std::isfinite(time_step),
        "ActiveStress test time step must be finite and positive");

    const double step_count = final_time / time_step;
    const double rounded_step_count = std::round(step_count);
    const double scale = std::max(1.0, std::fabs(step_count));
    svmp::throw_if<svmp::FE::InvalidArgumentException>(
        std::fabs(step_count - rounded_step_count) > 1.0e-12 * scale,
        "ActiveStress test final time must be an integer multiple of the time step");
    svmp::throw_if<svmp::FE::InvalidArgumentException>(
        rounded_step_count > std::numeric_limits<int>::max(),
        "ActiveStress test has too many time steps");

    return static_cast<int>(rounded_step_count);
  }

  /// Remove leading and trailing whitespace.
  static std::string trim(const std::string &value)
  {
    size_t begin = 0;
    while (begin < value.size() &&
           std::isspace(static_cast<unsigned char>(value[begin])))
      ++begin;

    size_t end = value.size();
    while (end > begin &&
           std::isspace(static_cast<unsigned char>(value[end - 1])))
      --end;

    return value.substr(begin, end - begin);
  }

  /// Split a CSV row, trim its fields, and preserve a trailing empty field.
  static std::vector<std::string> split_csv_row(const std::string &line)
  {
    std::vector<std::string> fields;
    std::istringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ','))
      fields.push_back(trim(field));
    if (!line.empty() && line.back() == ',')
      fields.emplace_back();
    return fields;
  }

  /// Parse one finite floating-point value from a CSV field.
  static double parse_double(const std::string &field,
                             const std::string &path,
                             size_t line_number)
  {
    size_t parsed_characters = 0;
    double value = 0.0;
    try {
      value = std::stod(field, &parsed_characters);
    } catch (const std::exception &) {
      svmp::raise<svmp::ParseException>(
          "Invalid number in " + path + " at line " +
          std::to_string(line_number));
    }

    svmp::throw_if<svmp::ParseException>(
        parsed_characters != field.size() || !std::isfinite(value),
        "Invalid number in " + path + " at line " +
            std::to_string(line_number));
    return value;
  }

  /// Load and validate the reference checkpoints and state-column mapping.
  static std::vector<Checkpoint>
  load_reference_data(const std::string &path,
                      unsigned int n_states,
                      int n_steps,
                      const std::vector<std::string> &requested_state_columns)
  {
    std::ifstream file(path);
    svmp::throw_if<svmp::FileNotFoundException>(!file.is_open(), path);

    std::vector<std::string> header;
    std::vector<size_t> state_columns;
    std::vector<Checkpoint> checkpoints;
    std::string line;
    size_t line_number = 0;

    while (std::getline(file, line)) {
      ++line_number;
      const std::string stripped_line = trim(line);
      if (stripped_line.empty() || stripped_line[0] == '#')
        continue;

      if (header.empty()) {
        header = split_csv_row(stripped_line);
        svmp::throw_if<svmp::ParseException>(
            header.size() < 2 || header[0] != "step" || header[1] != "Ta",
            "ActiveStress reference CSV must begin with columns step,Ta: " +
                path);

        for (size_t i = 0; i < header.size(); ++i) {
          svmp::throw_if<svmp::ParseException>(
              header[i].empty(),
              "ActiveStress reference CSV has an empty column name: " + path);
          svmp::throw_if<svmp::ParseException>(
              std::find(header.begin(), header.begin() + i, header[i]) !=
                  header.begin() + i,
              "ActiveStress reference CSV has duplicate column " + header[i] +
                  ": " + path);
        }

        state_columns.reserve(n_states);
        for (const std::string &column_name : requested_state_columns) {
          const auto column = std::find(header.begin(), header.end(), column_name);
          svmp::throw_if<svmp::ParseException>(
              column == header.end(),
              "ActiveStress reference CSV is missing state column " +
                  column_name + ": " + path);
          state_columns.push_back(
              static_cast<size_t>(std::distance(header.begin(), column)));
        }
        continue;
      }

      const auto fields = split_csv_row(stripped_line);
      svmp::throw_if<svmp::ParseException>(
          fields.size() != header.size(),
          "Wrong column count in " + path + " at line " +
              std::to_string(line_number));

      const double step_value =
          parse_double(fields[0], path, line_number);
      const double rounded_step = std::round(step_value);
      svmp::throw_if<svmp::ParseException>(
          std::fabs(step_value - rounded_step) > 1.0e-12 ||
              rounded_step < 0.0 || rounded_step >= n_steps,
          "Invalid checkpoint step in " + path + " at line " +
              std::to_string(line_number));

      Checkpoint checkpoint;
      checkpoint.step = static_cast<int>(rounded_step);
      checkpoint.active_tension =
          parse_double(fields[1], path, line_number);
      checkpoint.state.reserve(n_states);
      for (size_t column : state_columns)
        checkpoint.state.push_back(
            parse_double(fields[column], path, line_number));

      svmp::throw_if<svmp::ParseException>(
          !checkpoints.empty() &&
              checkpoint.step <= checkpoints.back().step,
          "ActiveStress reference checkpoints must be strictly increasing: " +
              path);
      checkpoints.push_back(std::move(checkpoint));
    }

    svmp::throw_if<svmp::ParseException>(
        header.empty(), "ActiveStress reference CSV has no header: " + path);
    svmp::throw_if<svmp::ParseException>(
        checkpoints.empty(),
        "ActiveStress reference CSV has no checkpoints: " + path);
    return checkpoints;
  }

  /**
   * @brief Evaluate the prescribed calcium transient at @p time.
   *
   * The calcium concentration is
   * @f[
   * C(t) =
   * \begin{cases}
   * C_0, & t < t_0, \\[4pt]
   * C_0 + (C_{\max}-C_0)
   * \dfrac{
   * e^{-(t-t_0)/\tau_d} - e^{-(t-t_0)/\tau_r}
   * }{
   * e^{-t_p/\tau_d} - e^{-t_p/\tau_r}
   * }, & t \ge t_0,
   * \end{cases}
   * @f]
   * where
   * @f[
   * C_0 = 10^{-4}\,\mathrm{mM}, \quad
   * C_{\max} = 9\times10^{-4}\,\mathrm{mM}, \quad
   * \tau_r = 20\,\mathrm{ms}, \quad
   * \tau_d = 50\,\mathrm{ms}, \quad
   * t_0 = 10\,\mathrm{ms},
   * @f]
   * and
   * @f[
   * t_p =
   * \frac{\log(\tau_r/\tau_d)\tau_r\tau_d}
   * {\tau_r-\tau_d}.
   * @f]
   */
  static double calcium_at(double time)
  {
    constexpr double baseline_calcium = 1.0e-4;
    constexpr double peak_calcium = 9.0e-4;
    constexpr double rise_time = 20.0;
    constexpr double decay_time = 50.0;
    constexpr double onset_time = 10.0;

    if (time < onset_time)
      return baseline_calcium;

    const double peak_time =
        std::log(rise_time / decay_time) * rise_time * decay_time /
        (rise_time - decay_time);
    const double peak_raw_factor =
        std::exp(-peak_time / decay_time) -
        std::exp(-peak_time / rise_time);
    const double elapsed = time - onset_time;
    const double raw_calcium =
        std::exp(-elapsed / decay_time) - std::exp(-elapsed / rise_time);
    return baseline_calcium +
           (peak_calcium - baseline_calcium) * raw_calcium / peak_raw_factor;
  }

  /**
   * @brief Evaluate the prescribed fiber stretch at @p time.
   *
   * Fiber stretch is defined as
   * @f[
   * \lambda(t) = \frac{L(t)}{L_0},
   * @f]
   * with
   * @f[
   * L(t) =
   * \begin{cases}
   * L_0, & t < t_s, \\[4pt]
   * L_0 - \dfrac{L_0-L_{\min}}{2}
   * \left[
   * 1-\cos\left(\pi\frac{t-t_s}{t_m-t_s}\right)
   * \right],
   * & t_s \le t \le t_m, \\[8pt]
   * L_{\min} + \dfrac{L_0-L_{\min}}{2}
   * \left[
   * 1-\cos\left(\pi\frac{t-t_m}{t_r-t_m}\right)
   * \right],
   * & t_m < t \le t_r, \\[8pt]
   * L_0, & t > t_r,
   * \end{cases}
   * @f]
   * where
   * @f[
   * L_0 = 2.2\,\mu\mathrm{m}, \quad
   * L_{\min} = 2.134\,\mu\mathrm{m}, \quad
   * t_s = 30\,\mathrm{ms}, \quad
   * t_m = 150\,\mathrm{ms}, \quad
   * t_r = 350\,\mathrm{ms}.
   * @f]
   */
  double fiber_stretch_at(double time) const
  {
    constexpr double resting_sarcomere_length = 2.2;
    constexpr double minimum_sarcomere_length = 2.134;
    constexpr double shortening_onset = 30.0;
    constexpr double minimum_time = 150.0;
    constexpr double recovery_time = 350.0;

    double sarcomere_length = resting_sarcomere_length;
    if (time >= shortening_onset && time <= minimum_time) {
      const double fraction =
          (time - shortening_onset) / (minimum_time - shortening_onset);
      sarcomere_length =
          resting_sarcomere_length -
          (resting_sarcomere_length - minimum_sarcomere_length) *
              (1.0 - std::cos(std::numbers::pi * fraction)) / 2.0;
    } else if (time > minimum_time && time <= recovery_time) {
      const double fraction =
          (time - minimum_time) / (recovery_time - minimum_time);
      sarcomere_length =
          minimum_sarcomere_length +
          (resting_sarcomere_length - minimum_sarcomere_length) *
              (1.0 - std::cos(std::numbers::pi * fraction)) / 2.0;
    }

    return sarcomere_length / resting_sarcomere_length;
  }

  /// Return the pointwise tolerance @f$\mathrm{tol}(1 + |expected|)@f$.
  double scaled_tolerance(double expected) const
  {
    return configuration_.tolerance * (1.0 + std::fabs(expected));
  }

  /// Compare the checkpoint matching @p step, if one is listed.
  void compare_checkpoint_if_present(int step,
                                     const Vector<double> &state,
                                     double active_tension,
                                     size_t &checkpoint_index) const
  {
    if (checkpoint_index >= checkpoints_.size() ||
        checkpoints_[checkpoint_index].step != step)
      return;

    compare_checkpoint(checkpoints_[checkpoint_index], state, active_tension);
    ++checkpoint_index;
  }

  /// Compare active tension and every state value at one checkpoint.
  void compare_checkpoint(const Checkpoint &checkpoint,
                          const Vector<double> &state,
                          double active_tension) const
  {
    EXPECT_NEAR(active_tension, checkpoint.active_tension,
                scaled_tolerance(checkpoint.active_tension))
        << "checkpoint at trajectory update index " << checkpoint.step
        << ", active tension column Ta";

    ASSERT_EQ(state.size(), checkpoint.state.size())
        << "checkpoint at trajectory update index " << checkpoint.step
        << ", state-vector size";
    for (unsigned int index = 0; index < state.size(); ++index) {
      EXPECT_NEAR(state[index], checkpoint.state[index],
                  scaled_tolerance(checkpoint.state[index]))
          << "checkpoint at trajectory update index " << checkpoint.step
          << ", state[" << index << "] column "
          << configuration_.state_reference_columns[index];
    }
  }

  /// Configured test-only model adapter.
  std::unique_ptr<TestAdapter> model_;

  /// Validated trajectory inputs and reference-column mapping.
  ActiveStressTrajectoryConfiguration configuration_;

  /// Number of trajectory updates.
  int n_steps_;

  /// Full path to the reference CSV.
  std::string reference_path_;

  /// Parsed reference checkpoints in increasing step order.
  std::vector<Checkpoint> checkpoints_;
};

#endif // ACTIVE_STRESS_TEST_HELPERS_H
