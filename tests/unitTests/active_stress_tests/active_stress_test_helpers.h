// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef ACTIVE_STRESS_TEST_HELPERS_H
#define ACTIVE_STRESS_TEST_HELPERS_H

#include "active_stress.h"
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

/**
 * @brief Run a trusted-reference trajectory test for an ActiveStress model.
 *
 * A test-only adapter derived from the concrete model exposes the protected
 * node-local interface needed to run the trajectory experiment.
 * By default, state values are read from CSV columns @c s0 through @c sN.
 * The constructor's state-column mapping can explicitly select other columns
 * when an existing reference file uses a different, documented layout.
 *
 * Checkpoint indices are zero based: checkpoint @c N is compared with the
 * state after advance @c N, which advances the interval
 * [N * time_step, (N + 1) * time_step].
 */
template <class ConcreteModel>
class ActiveStressTrajectoryTest {
  static_assert(std::is_base_of_v<ActiveStress, ConcreteModel>,
                "ConcreteModel must derive from ActiveStress");

  class TestAdapter : public ConcreteModel {
  public:
    using ConcreteModel::advance_time_step_local;
    using ConcreteModel::compute_active_tension_local;
    using ConcreteModel::init_local;
    using ConcreteModel::read_model_specific_parameters;
  };

public:
  using Parameters = typename ConcreteModel::Parameters;

  ActiveStressTrajectoryTest(const Parameters &parameters,
                             double final_time,
                             double time_step,
                             const std::string &reference_csv_filename,
                             std::vector<std::string> state_reference_columns = {},
                             double tolerance = 1.0e-10)
      : model_(make_configured_model(parameters)),
        time_step_(time_step),
        tolerance_(tolerance),
        n_steps_(number_of_steps(final_time, time_step)),
        reference_path_(std::string(UNIT_TEST_DATA_DIR) + "/" +
                        reference_csv_filename),
        checkpoints_(load_reference_data(reference_path_, model_->n_states,
                                         n_steps_, state_reference_columns))
  {
    if (!(tolerance_ > 0.0) || !std::isfinite(tolerance_))
      throw std::invalid_argument("ActiveStress test tolerance must be finite and positive");
  }

  /** @brief Execute the common twitch protocol and compare every checkpoint. */
  void run()
  {
    Vector<double> state(model_->n_states);
    model_->init_local(state);

    size_t checkpoint_index = 0;
    for (int step = 0; step < n_steps_; ++step) {
      const double time = step * time_step_;
      const double calcium = calcium_at(time);
      const double fiber_stretch = fiber_stretch_at(time);
      const double fiber_stretch_rate =
          (fiber_stretch_at(time + time_step_) - fiber_stretch) / time_step_;

      model_->advance_time_step_local(time, time_step_, calcium, fiber_stretch,
                                      fiber_stretch_rate, state);
      const double active_tension =
          model_->compute_active_tension_local(state, fiber_stretch);

      if (checkpoint_index < checkpoints_.size() &&
          checkpoints_[checkpoint_index].step == step) {
        compare_checkpoint(checkpoints_[checkpoint_index], state,
                           active_tension);
        ++checkpoint_index;
      }
    }

    EXPECT_EQ(checkpoint_index, checkpoints_.size())
        << "Not all reference checkpoints were visited in " << reference_path_;
  }

private:
  struct Checkpoint {
    int step;
    double active_tension;
    std::vector<double> state;
  };

  static std::unique_ptr<TestAdapter>
  make_configured_model(const Parameters &parameters)
  {
    auto model = std::make_unique<TestAdapter>();
    model->read_model_specific_parameters(parameters);
    return model;
  }

  static int number_of_steps(double final_time, double time_step)
  {
    if (!(final_time > 0.0) || !std::isfinite(final_time))
      throw std::invalid_argument("ActiveStress test final time must be finite and positive");
    if (!(time_step > 0.0) || !std::isfinite(time_step))
      throw std::invalid_argument("ActiveStress test time step must be finite and positive");

    const double step_count = final_time / time_step;
    const double rounded_step_count = std::round(step_count);
    const double scale = std::max(1.0, std::fabs(step_count));
    if (std::fabs(step_count - rounded_step_count) >
        1.0e-12 * scale)
      throw std::invalid_argument(
          "ActiveStress test final time must be an integer multiple of the time step");
    if (rounded_step_count > std::numeric_limits<int>::max())
      throw std::invalid_argument("ActiveStress test has too many time steps");

    return static_cast<int>(rounded_step_count);
  }

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

  static std::vector<std::string> split_csv_row(const std::string &line)
  {
    std::vector<std::string> fields;
    std::istringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ','))
      fields.push_back(trim(field));
    return fields;
  }

  static double parse_double(const std::string &field,
                             const std::string &path,
                             size_t line_number)
  {
    size_t parsed_characters = 0;
    double value = 0.0;
    try {
      value = std::stod(field, &parsed_characters);
    } catch (const std::exception &) {
      throw std::runtime_error("Invalid number in " + path + " at line " +
                               std::to_string(line_number));
    }

    if (parsed_characters != field.size() || !std::isfinite(value))
      throw std::runtime_error("Invalid number in " + path + " at line " +
                               std::to_string(line_number));
    return value;
  }

  static std::vector<Checkpoint>
  load_reference_data(const std::string &path,
                      unsigned int n_states,
                      int n_steps,
                      const std::vector<std::string> &requested_state_columns)
  {
    std::ifstream file(path);
    if (!file.is_open())
      throw std::runtime_error("Cannot open reference CSV: " + path);

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
        if (header.size() < 2 || header[0] != "step" || header[1] != "Ta")
          throw std::runtime_error(
              "ActiveStress reference CSV must begin with columns step,Ta: " +
              path);

        std::vector<std::string> state_column_names = requested_state_columns;
        if (state_column_names.empty()) {
          state_column_names.reserve(n_states);
          for (unsigned int state = 0; state < n_states; ++state)
            state_column_names.push_back("s" + std::to_string(state));
        }
        if (state_column_names.size() != n_states)
          throw std::runtime_error(
              "ActiveStress reference state-column mapping has the wrong size: " +
              path);

        state_columns.reserve(n_states);
        for (const std::string &column_name : state_column_names) {
          const auto column = std::find(header.begin(), header.end(), column_name);
          if (column == header.end())
            throw std::runtime_error(
                "ActiveStress reference CSV is missing state column " +
                column_name + ": " + path);
          if (std::find(column + 1, header.end(), column_name) != header.end())
            throw std::runtime_error(
                "ActiveStress reference CSV has duplicate column " +
                column_name + ": " + path);
          state_columns.push_back(
              static_cast<size_t>(std::distance(header.begin(), column)));
        }
        continue;
      }

      const auto fields = split_csv_row(stripped_line);
      if (fields.size() != header.size())
        throw std::runtime_error("Wrong column count in " + path + " at line " +
                                 std::to_string(line_number));

      const double step_value =
          parse_double(fields[0], path, line_number);
      const double rounded_step = std::round(step_value);
      if (std::fabs(step_value - rounded_step) > 1.0e-12 ||
          rounded_step < 0.0 || rounded_step >= n_steps)
        throw std::runtime_error("Invalid checkpoint step in " + path +
                                 " at line " +
                                 std::to_string(line_number));

      Checkpoint checkpoint;
      checkpoint.step = static_cast<int>(rounded_step);
      checkpoint.active_tension =
          parse_double(fields[1], path, line_number);
      checkpoint.state.reserve(n_states);
      for (size_t column : state_columns)
        checkpoint.state.push_back(
            parse_double(fields[column], path, line_number));

      if (!checkpoints.empty() &&
          checkpoint.step <= checkpoints.back().step)
        throw std::runtime_error(
            "ActiveStress reference checkpoints must be strictly increasing: " +
            path);
      checkpoints.push_back(std::move(checkpoint));
    }

    if (header.empty())
      throw std::runtime_error("ActiveStress reference CSV has no header: " + path);
    if (checkpoints.empty())
      throw std::runtime_error("ActiveStress reference CSV has no checkpoints: " +
                               path);
    return checkpoints;
  }

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

  double scaled_tolerance(double expected) const
  {
    return tolerance_ * (1.0 + std::fabs(expected));
  }

  void compare_checkpoint(const Checkpoint &checkpoint,
                          const Vector<double> &state,
                          double active_tension) const
  {
    EXPECT_NEAR(active_tension, checkpoint.active_tension,
                scaled_tolerance(checkpoint.active_tension))
        << "step " << checkpoint.step << ", active tension";

    ASSERT_EQ(state.size(), checkpoint.state.size())
        << "step " << checkpoint.step << ", state-vector size";
    for (unsigned int index = 0; index < state.size(); ++index) {
      EXPECT_NEAR(state[index], checkpoint.state[index],
                  scaled_tolerance(checkpoint.state[index]))
          << "step " << checkpoint.step << ", state[" << index << "]";
    }
  }

  std::unique_ptr<TestAdapter> model_;
  double time_step_;
  double tolerance_;
  int n_steps_;
  std::string reference_path_;
  std::vector<Checkpoint> checkpoints_;
};

#endif // ACTIVE_STRESS_TEST_HELPERS_H
