// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef IONIC_MODEL_TEST_HELPERS_H
#define IONIC_MODEL_TEST_HELPERS_H

#include "ionic_model.h"
#include "Vector.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iterator>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

/** @brief Configuration for one trusted-reference IonicModel trajectory. */
struct IonicModelTrajectoryConfiguration {
  TimeIntegrationType integration_type = TimeIntegrationType::NA;
  int zone_id = 1;
  double time_step = 0.0;
  std::size_t update_count = 0;

  std::optional<std::vector<double>> initial_X_override;
  std::optional<std::vector<double>> initial_Xg_override;

  std::string reference_csv_filename;

  /**
   * @brief CSV columns corresponding to all X entries, followed by all Xg
   * entries, in their model-vector order.
   */
  std::vector<std::string> state_reference_columns;

  double tolerance = 1.0e-10;
  double stimulus = 0.0;
  double sac_coefficient = 0.0;
};

/**
 * @brief Run a trusted-reference trajectory test for an IonicModel.
 *
 * The concrete type is used only for construction and parameter typing. After
 * configuration, initialization and trajectory execution use the public
 * IonicModel interface. Reference values must use the representation exposed
 * by that interface; this helper intentionally performs no voltage, time, or
 * state transformations.
 *
 * Checkpoint N is the state after exactly N completed updates. In particular,
 * checkpoint 0 is compared immediately after init() and optional state
 * overrides, before the first call to integ(). Each update starting at N uses
 * t = N * time_step, so the resulting state is checkpoint N + 1.
 */
template <class ConcreteModel>
class IonicModelTrajectoryTest {
  static_assert(std::is_base_of_v<IonicModel, ConcreteModel>,
                "ConcreteModel must derive from IonicModel");

public:
  using Parameters = typename ConcreteModel::Parameters;

  static_assert(std::is_base_of_v<IonicModelParameters, Parameters>,
                "ConcreteModel::Parameters must derive from "
                "IonicModelParameters");

  IonicModelTrajectoryTest(
      const Parameters &parameters,
      IonicModelTrajectoryConfiguration configuration)
      : model_(make_configured_model(parameters)),
        configuration_(std::move(configuration))
  {
    validate_configuration();
    reference_path_ = std::string(UNIT_TEST_DATA_DIR) + "/" +
                      configuration_.reference_csv_filename;
    checkpoints_ = load_reference_data();
  }

  /** @brief Initialize, integrate, and compare every reference checkpoint. */
  void run()
  {
    Vector<double> X(model_->nX());
    Vector<double> Xg(model_->nG());
    model_->init(X, Xg);

    apply_override(configuration_.initial_X_override, X);
    apply_override(configuration_.initial_Xg_override, Xg);

    odeType solver;
    solver.tIntType = configuration_.integration_type;

    std::size_t checkpoint_index = 0;
    compare_checkpoint_if_present(/* completed_updates = */ 0, X, Xg,
                                  checkpoint_index);

    for (std::size_t completed_updates = 1;
         completed_updates <= configuration_.update_count;
         ++completed_updates) {
      const double time =
          static_cast<double>(completed_updates - 1) *
          configuration_.time_step;
      model_->integ(solver, configuration_.zone_id, time,
                    configuration_.time_step, configuration_.stimulus,
                    configuration_.sac_coefficient, X, Xg);
      compare_checkpoint_if_present(completed_updates, X, Xg,
                                    checkpoint_index);
    }

    EXPECT_EQ(checkpoint_index, checkpoints_.size())
        << "Not all reference checkpoints were visited in " << reference_path_;
  }

private:
  struct Checkpoint {
    std::size_t completed_updates;
    std::vector<double> state;
  };

  static std::unique_ptr<IonicModel>
  make_configured_model(const Parameters &parameters)
  {
    std::unique_ptr<IonicModel> model = std::make_unique<ConcreteModel>();
    model->read_parameters(parameters);
    return model;
  }

  void validate_configuration() const
  {
    switch (configuration_.integration_type) {
    case TimeIntegrationType::FE:
    case TimeIntegrationType::RK4:
    case TimeIntegrationType::CN2:
      break;
    default:
      throw std::invalid_argument(
          "IonicModel trajectory integration type must be FE, RK4, or CN2");
    }

    if (configuration_.zone_id <= 0)
      throw std::invalid_argument(
          "IonicModel trajectory zone_id must be positive");
    if (!(configuration_.time_step > 0.0) ||
        !std::isfinite(configuration_.time_step))
      throw std::invalid_argument(
          "IonicModel trajectory time step must be finite and positive");
    if (configuration_.update_count == 0)
      throw std::invalid_argument(
          "IonicModel trajectory update count must be positive");
    if (!std::isfinite(configuration_.time_step *
                       static_cast<double>(configuration_.update_count)))
      throw std::invalid_argument(
          "IonicModel trajectory final time must be finite");
    if (!(configuration_.tolerance > 0.0) ||
        !std::isfinite(configuration_.tolerance))
      throw std::invalid_argument(
          "IonicModel trajectory tolerance must be finite and positive");
    if (!std::isfinite(configuration_.stimulus) ||
        !std::isfinite(configuration_.sac_coefficient))
      throw std::invalid_argument(
          "IonicModel trajectory stimulus and SAC coefficient must be finite");
    if (configuration_.reference_csv_filename.empty())
      throw std::invalid_argument(
          "IonicModel trajectory reference CSV filename must not be empty");

    validate_override(configuration_.initial_X_override, model_->nX(), "X");
    validate_override(configuration_.initial_Xg_override, model_->nG(), "Xg");

    const std::size_t state_count = model_->nX() + model_->nG();
    if (configuration_.state_reference_columns.size() != state_count)
      throw std::invalid_argument(
          "IonicModel reference state-column mapping must contain all X "
          "columns followed by all Xg columns");

    for (std::size_t i = 0;
         i < configuration_.state_reference_columns.size(); ++i) {
      const auto &name = configuration_.state_reference_columns[i];
      if (name.empty())
        throw std::invalid_argument(
            "IonicModel reference state-column names must not be empty");
      if (name == "step")
        throw std::invalid_argument(
            "IonicModel reference state-column mapping cannot use step");
      if (std::find(configuration_.state_reference_columns.begin(),
                    configuration_.state_reference_columns.begin() + i,
                    name) !=
          configuration_.state_reference_columns.begin() + i)
        throw std::invalid_argument(
            "IonicModel reference state-column names must be unique");
    }
  }

  static void
  validate_override(const std::optional<std::vector<double>> &values,
                    std::size_t expected_size,
                    const std::string &state_name)
  {
    if (!values)
      return;
    if (values->size() != expected_size)
      throw std::invalid_argument("IonicModel initial " + state_name +
                                  " override has the wrong size");
    if (!std::all_of(values->begin(), values->end(),
                     [](double value) { return std::isfinite(value); }))
      throw std::invalid_argument("IonicModel initial " + state_name +
                                  " override must contain finite values");
  }

  static void
  apply_override(const std::optional<std::vector<double>> &values,
                 Vector<double> &state)
  {
    if (!values)
      return;
    for (std::size_t index = 0; index < values->size(); ++index)
      state[index] = (*values)[index];
  }

  static std::string trim(const std::string &value)
  {
    std::size_t begin = 0;
    while (begin < value.size() &&
           std::isspace(static_cast<unsigned char>(value[begin])))
      ++begin;

    std::size_t end = value.size();
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
    if (!line.empty() && line.back() == ',')
      fields.emplace_back();
    return fields;
  }

  static double parse_double(const std::string &field,
                             const std::string &path,
                             std::size_t line_number)
  {
    std::size_t parsed_characters = 0;
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

  std::vector<Checkpoint> load_reference_data() const
  {
    std::ifstream file(reference_path_);
    if (!file.is_open())
      throw std::runtime_error("Cannot open IonicModel reference CSV: " +
                               reference_path_);

    std::vector<std::string> header;
    std::vector<std::size_t> state_columns;
    std::vector<Checkpoint> checkpoints;
    std::string line;
    std::size_t line_number = 0;

    while (std::getline(file, line)) {
      ++line_number;
      const std::string stripped_line = trim(line);
      if (stripped_line.empty() || stripped_line[0] == '#')
        continue;

      if (header.empty()) {
        header = split_csv_row(stripped_line);
        validate_header(header, state_columns);
        continue;
      }

      const auto fields = split_csv_row(stripped_line);
      if (fields.size() != header.size())
        throw std::runtime_error(
            "Wrong column count in " + reference_path_ + " at line " +
            std::to_string(line_number));

      const double completed_updates_value =
          parse_double(fields[0], reference_path_, line_number);
      const double rounded_completed_updates =
          std::round(completed_updates_value);
      if (std::fabs(completed_updates_value - rounded_completed_updates) >
              1.0e-12 ||
          rounded_completed_updates < 0.0)
        throw std::runtime_error(
            "Invalid checkpoint step in " + reference_path_ + " at line " +
            std::to_string(line_number));

      if (rounded_completed_updates >
          static_cast<double>(configuration_.update_count))
        throw std::runtime_error(
            "Checkpoint step is outside the configured update range in " +
            reference_path_ + " at line " + std::to_string(line_number));

      Checkpoint checkpoint;
      checkpoint.completed_updates =
          static_cast<std::size_t>(rounded_completed_updates);
      checkpoint.state.reserve(state_columns.size());
      for (std::size_t column : state_columns)
        checkpoint.state.push_back(
            parse_double(fields[column], reference_path_, line_number));

      if (!checkpoints.empty() &&
          checkpoint.completed_updates <=
              checkpoints.back().completed_updates)
        throw std::runtime_error(
            "IonicModel reference checkpoints must be strictly increasing: " +
            reference_path_);
      checkpoints.push_back(std::move(checkpoint));
    }

    if (header.empty())
      throw std::runtime_error("IonicModel reference CSV has no header: " +
                               reference_path_);
    if (checkpoints.empty())
      throw std::runtime_error(
          "IonicModel reference CSV has no checkpoints: " + reference_path_);
    return checkpoints;
  }

  void validate_header(const std::vector<std::string> &header,
                       std::vector<std::size_t> &state_columns) const
  {
    if (header.size() < 2 || header.front() != "step")
      throw std::runtime_error(
          "IonicModel reference CSV must begin with a step column: " +
          reference_path_);

    for (std::size_t i = 0; i < header.size(); ++i) {
      if (header[i].empty())
        throw std::runtime_error(
            "IonicModel reference CSV has an empty column name: " +
            reference_path_);
      if (std::find(header.begin(), header.begin() + i, header[i]) !=
          header.begin() + i)
        throw std::runtime_error(
            "IonicModel reference CSV has duplicate column " + header[i] +
            ": " + reference_path_);
    }

    state_columns.reserve(configuration_.state_reference_columns.size());
    for (const std::string &column_name :
         configuration_.state_reference_columns) {
      const auto column =
          std::find(header.begin(), header.end(), column_name);
      if (column == header.end())
        throw std::runtime_error(
            "IonicModel reference CSV is missing state column " + column_name +
            ": " + reference_path_);
      state_columns.push_back(
          static_cast<std::size_t>(std::distance(header.begin(), column)));
    }
  }

  double scaled_tolerance(double expected) const
  {
    return configuration_.tolerance * (1.0 + std::fabs(expected));
  }

  void compare_checkpoint_if_present(std::size_t completed_updates,
                                     const Vector<double> &X,
                                     const Vector<double> &Xg,
                                     std::size_t &checkpoint_index) const
  {
    if (checkpoint_index >= checkpoints_.size() ||
        checkpoints_[checkpoint_index].completed_updates != completed_updates)
      return;

    const Checkpoint &checkpoint = checkpoints_[checkpoint_index];
    std::size_t expected_index = 0;
    for (std::size_t state = 0; state < X.size(); ++state, ++expected_index) {
      const double expected = checkpoint.state[expected_index];
      EXPECT_NEAR(X[state], expected, scaled_tolerance(expected))
          << "checkpoint after " << checkpoint.completed_updates
          << " completed updates, X[" << state << "] column "
          << configuration_.state_reference_columns[expected_index];
    }
    for (std::size_t state = 0; state < Xg.size(); ++state, ++expected_index) {
      const double expected = checkpoint.state[expected_index];
      EXPECT_NEAR(Xg[state], expected, scaled_tolerance(expected))
          << "checkpoint after " << checkpoint.completed_updates
          << " completed updates, Xg[" << state << "] column "
          << configuration_.state_reference_columns[expected_index];
    }
    ++checkpoint_index;
  }

  std::unique_ptr<IonicModel> model_;
  IonicModelTrajectoryConfiguration configuration_;
  std::string reference_path_;
  std::vector<Checkpoint> checkpoints_;
};

#endif // IONIC_MODEL_TEST_HELPERS_H
