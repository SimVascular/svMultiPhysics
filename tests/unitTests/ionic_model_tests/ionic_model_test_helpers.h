// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef IONIC_MODEL_TEST_HELPERS_H
#define IONIC_MODEL_TEST_HELPERS_H

#include "ionic_model.h"
#include "Core/Exception.h"
#include "FE/Common/FEException.h"
#include "Vector.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <iterator>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

/**
 * @brief Configure one standalone trusted-reference IonicModel trajectory.
 *
 * The configuration selects the public integration inputs and maps every
 * model state to a column in the external reference CSV.
 */
struct IonicModelTrajectoryConfiguration {
  /// Integration method passed to IonicModel::integ; must be FE, RK4, or CN2.
  TimeIntegrationType integration_type = TimeIntegrationType::NA;

  /// Positive zone identifier passed to every trajectory update.
  int zone_id = 1;

  /// Time step passed to IonicModel::integ and used to compute update times.
  double time_step = 0.0;

  /// Number of calls to IonicModel::integ.
  std::size_t update_count = 0;

  /// Optional complete X replacement applied after IonicModel::init.
  std::optional<std::vector<double>> initial_X_override;

  /// Optional complete Xg replacement applied after IonicModel::init.
  std::optional<std::vector<double>> initial_Xg_override;

  /// Reference filename under UNIT_TEST_DATA_DIR.
  std::string reference_csv_filename;

  /**
   * @brief CSV column name for every model state.
   *
   * Entries correspond in order to @c X[0] through @c X[nX-1], followed by
   * @c Xg[0] through @c Xg[nG-1]. Each entry names the CSV column containing
   * the reference value for that state; the column name need not match the
   * model state name. Every model state requires an entry, but the physical
   * CSV column order may differ from this mapping order.
   */
  std::vector<std::string> state_reference_columns;

  /// Base tolerance in the pointwise comparison formula.
  double tolerance = 1.0e-10;

  /// Constant stimulus used when stimulus_at_time is not set.
  double stimulus = 0.0;

  /**
   * @brief Optional stimulus evaluated once per update at the old-state time.
   *
   * When set, this function supplies the stimulus instead of the constant
   * value above.
   */
  std::function<double(double)> stimulus_at_time;

  /// Constant SAC coefficient passed to every trajectory update.
  double sac_coefficient = 0.0;
};

/**
 * @brief Run a standalone IonicModel trajectory and compare it with a trusted
 * external reference solution loaded from a CSV file.
 *
 * The concrete model type is given by @c ConcreteModel and is used only for
 * construction and parameter typing. The test allocates @c X and @c Xg, calls
 * IonicModel::init, applies any complete state overrides, and advances the
 * model through the public IonicModel::integ interface for the configured
 * number of updates. Only the ionic model is advanced; no coupled tissue or
 * mechanics problem is solved. Reference values must use the representation
 * exposed by the public interface; this helper intentionally performs no
 * voltage, time, or state transformations.
 *
 * If @c stimulus_at_time is set, it is evaluated once at the left endpoint of
 * each update and its result is used for the complete call to
 * IonicModel::integ. Otherwise, the constant @c stimulus value is used.
 *
 * ### Format of the reference solution
 *
 * The reference CSV contains a header followed by values at selected
 * trajectory checkpoints. The first column is @c step. Entries in
 * @c state_reference_columns map in order to @c X[0] through @c X[nX-1], then
 * @c Xg[0] through @c Xg[nG-1]; each entry names the CSV column containing the
 * corresponding state value. Those columns may occur in any physical order
 * after @c step; additional columns are ignored.
 *
 * Comparisons are performed only at the checkpoints listed in the CSV; the
 * file does not need to contain every trajectory update. Checkpoints must be
 * nonnegative integer update counts, strictly increasing, and less than or
 * equal to @c update_count. Intermediate checkpoints are not interpolated.
 *
 * Checkpoint N is the state after exactly N completed updates. In particular,
 * checkpoint 0 is compared immediately after init() and optional state
 * overrides, before the first call to integ(). Each update starting at N uses
 * t = N * time_step, so the resulting state is checkpoint N + 1.
 *
 * At each checkpoint, every @c X and @c Xg value is compared pointwise using
 * the tolerance @f$\mathrm{tol}(1+|expected|)@f$. A difference exceeding that
 * tolerance records a failed GoogleTest expectation.
 */
template <class ConcreteModel>
class IonicModelTrajectoryTest {
  static_assert(std::is_base_of_v<IonicModel, ConcreteModel>,
                "ConcreteModel must derive from IonicModel");

public:
  /// Model-specific parameter type used to construct the concrete model.
  using Parameters = typename ConcreteModel::Parameters;

  static_assert(std::is_base_of_v<IonicModelParameters, Parameters>,
                "ConcreteModel::Parameters must derive from "
                "IonicModelParameters");

  /**
   * @brief Configure the model, validate the trajectory inputs, and load the
   * reference CSV.
   *
   * @param parameters Parameters for the concrete ionic model.
   * @param configuration Trajectory inputs and reference-data mapping.
   */
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
      const double stimulus = evaluate_stimulus(time);
      model_->integ(solver, configuration_.zone_id, time,
                    configuration_.time_step, stimulus,
                    configuration_.sac_coefficient, X, Xg);
      compare_checkpoint_if_present(completed_updates, X, Xg,
                                    checkpoint_index);
    }

    EXPECT_EQ(checkpoint_index, checkpoints_.size())
        << "Not all reference checkpoints were visited in " << reference_path_;
  }

private:
  /// Reference state after a specified number of completed updates.
  struct Checkpoint {
    /// Number of completed trajectory updates.
    std::size_t completed_updates;

    /// Expected values ordered as all X entries followed by all Xg entries.
    std::vector<double> state;
  };

  /// Construct the concrete model, read its parameters, and store it as IonicModel.
  static std::unique_ptr<IonicModel>
  make_configured_model(const Parameters &parameters)
  {
    std::unique_ptr<IonicModel> model = std::make_unique<ConcreteModel>();
    model->read_parameters(parameters);
    return model;
  }

  /// Validate all trajectory inputs and state-column mappings.
  void validate_configuration() const
  {
    switch (configuration_.integration_type) {
    case TimeIntegrationType::FE:
    case TimeIntegrationType::RK4:
    case TimeIntegrationType::CN2:
      break;
    default:
      svmp::raise<svmp::FE::InvalidArgumentException>(
          "IonicModel trajectory integration type must be FE, RK4, or CN2");
    }

    svmp::throw_if<svmp::FE::InvalidArgumentException>(
        configuration_.zone_id <= 0,
        "IonicModel trajectory zone_id must be positive");
    svmp::throw_if<svmp::FE::InvalidArgumentException>(
        !(configuration_.time_step > 0.0) ||
            !std::isfinite(configuration_.time_step),
        "IonicModel trajectory time step must be finite and positive");
    svmp::throw_if<svmp::FE::InvalidArgumentException>(
        configuration_.update_count == 0,
        "IonicModel trajectory update count must be positive");
    svmp::throw_if<svmp::FE::InvalidArgumentException>(
        !std::isfinite(configuration_.time_step *
                       static_cast<double>(configuration_.update_count)),
        "IonicModel trajectory final time must be finite");
    svmp::throw_if<svmp::FE::InvalidArgumentException>(
        !(configuration_.tolerance > 0.0) ||
            !std::isfinite(configuration_.tolerance),
        "IonicModel trajectory tolerance must be finite and positive");
    svmp::throw_if<svmp::FE::InvalidArgumentException>(
        !std::isfinite(configuration_.stimulus) ||
            !std::isfinite(configuration_.sac_coefficient),
        "IonicModel trajectory stimulus and SAC coefficient must be finite");
    svmp::throw_if<svmp::FE::InvalidArgumentException>(
        configuration_.reference_csv_filename.empty(),
        "IonicModel trajectory reference CSV filename must not be empty");

    validate_override(configuration_.initial_X_override, model_->nX(), "X");
    validate_override(configuration_.initial_Xg_override, model_->nG(), "Xg");

    const std::size_t state_count = model_->nX() + model_->nG();
    svmp::throw_if<svmp::FE::InvalidArgumentException>(
        configuration_.state_reference_columns.size() != state_count,
        "IonicModel reference state-column mapping must contain all X "
        "columns followed by all Xg columns");

    for (std::size_t i = 0;
         i < configuration_.state_reference_columns.size(); ++i) {
      const auto &name = configuration_.state_reference_columns[i];
      svmp::throw_if<svmp::FE::InvalidArgumentException>(
          name.empty(),
          "IonicModel reference state-column names must not be empty");
      svmp::throw_if<svmp::FE::InvalidArgumentException>(
          name == "step",
          "IonicModel reference state-column mapping cannot use step");
      svmp::throw_if<svmp::FE::InvalidArgumentException>(
          std::find(configuration_.state_reference_columns.begin(),
                    configuration_.state_reference_columns.begin() + i,
                    name) !=
              configuration_.state_reference_columns.begin() + i,
          "IonicModel reference state-column names must be unique");
    }
  }

  /// Validate the size and finiteness of an optional complete state override.
  static void
  validate_override(const std::optional<std::vector<double>> &values,
                    std::size_t expected_size,
                    const std::string &state_name)
  {
    if (!values)
      return;
    svmp::throw_if<svmp::FE::InvalidArgumentException>(
        values->size() != expected_size,
        "IonicModel initial " + state_name + " override has the wrong size");
    svmp::throw_if<svmp::FE::InvalidArgumentException>(
        !std::all_of(values->begin(), values->end(),
                     [](double value) { return std::isfinite(value); }),
        "IonicModel initial " + state_name +
            " override must contain finite values");
  }

  /// Evaluate the callback or constant stimulus and require a finite result.
  double evaluate_stimulus(double time) const
  {
    const double stimulus = configuration_.stimulus_at_time
                                ? configuration_.stimulus_at_time(time)
                                : configuration_.stimulus;
    svmp::throw_if<svmp::FE::InvalidArgumentException>(
        !std::isfinite(stimulus),
        "IonicModel trajectory stimulus must be finite at time " +
            std::to_string(time));
    return stimulus;
  }

  /// Apply an optional complete replacement of a model state vector.
  static void
  apply_override(const std::optional<std::vector<double>> &values,
                 Vector<double> &state)
  {
    if (!values)
      return;
    for (std::size_t index = 0; index < values->size(); ++index)
      state[index] = (*values)[index];
  }

  /// Remove leading and trailing whitespace.
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

  /// Parse one complete finite floating-point value from a CSV field.
  static double parse_double(const std::string &field,
                             const std::string &path,
                             std::size_t line_number)
  {
    std::size_t parsed_characters = 0;
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

  /// Load and validate the reference checkpoints and selected state values.
  std::vector<Checkpoint> load_reference_data() const
  {
    std::ifstream file(reference_path_);
    svmp::throw_if<svmp::FileNotFoundException>(!file.is_open(),
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
      svmp::throw_if<svmp::ParseException>(
          fields.size() != header.size(),
          "Wrong column count in " + reference_path_ + " at line " +
              std::to_string(line_number));

      const double completed_updates_value =
          parse_double(fields[0], reference_path_, line_number);
      const double rounded_completed_updates =
          std::round(completed_updates_value);
      svmp::throw_if<svmp::ParseException>(
          std::fabs(completed_updates_value - rounded_completed_updates) >
                  1.0e-12 ||
              rounded_completed_updates < 0.0,
          "Invalid checkpoint step in " + reference_path_ + " at line " +
              std::to_string(line_number));

      svmp::throw_if<svmp::ParseException>(
          rounded_completed_updates >
              static_cast<double>(configuration_.update_count),
          "Checkpoint step is outside the configured update range in " +
              reference_path_ + " at line " + std::to_string(line_number));

      Checkpoint checkpoint;
      checkpoint.completed_updates =
          static_cast<std::size_t>(rounded_completed_updates);
      checkpoint.state.reserve(state_columns.size());
      for (std::size_t column : state_columns)
        checkpoint.state.push_back(
            parse_double(fields[column], reference_path_, line_number));

      svmp::throw_if<svmp::ParseException>(
          !checkpoints.empty() &&
              checkpoint.completed_updates <=
                  checkpoints.back().completed_updates,
          "IonicModel reference checkpoints must be strictly increasing: " +
              reference_path_);
      checkpoints.push_back(std::move(checkpoint));
    }

    svmp::throw_if<svmp::ParseException>(
        header.empty(),
        "IonicModel reference CSV has no header: " + reference_path_);
    svmp::throw_if<svmp::ParseException>(
        checkpoints.empty(),
        "IonicModel reference CSV has no checkpoints: " + reference_path_);
    return checkpoints;
  }

  /// Validate the CSV header and resolve the configured state-column indices.
  void validate_header(const std::vector<std::string> &header,
                       std::vector<std::size_t> &state_columns) const
  {
    svmp::throw_if<svmp::ParseException>(
        header.size() < 2 || header.front() != "step",
        "IonicModel reference CSV must begin with a step column: " +
            reference_path_);

    for (std::size_t i = 0; i < header.size(); ++i) {
      svmp::throw_if<svmp::ParseException>(
          header[i].empty(),
          "IonicModel reference CSV has an empty column name: " +
              reference_path_);
      svmp::throw_if<svmp::ParseException>(
          std::find(header.begin(), header.begin() + i, header[i]) !=
              header.begin() + i,
          "IonicModel reference CSV has duplicate column " + header[i] +
              ": " + reference_path_);
    }

    state_columns.reserve(configuration_.state_reference_columns.size());
    for (const std::string &column_name :
         configuration_.state_reference_columns) {
      const auto column =
          std::find(header.begin(), header.end(), column_name);
      svmp::throw_if<svmp::ParseException>(
          column == header.end(),
          "IonicModel reference CSV is missing state column " + column_name +
              ": " + reference_path_);
      state_columns.push_back(
          static_cast<std::size_t>(std::distance(header.begin(), column)));
    }
  }

  /// Return the pointwise tolerance @f$\mathrm{tol}(1+|expected|)@f$.
  double scaled_tolerance(double expected) const
  {
    return configuration_.tolerance * (1.0 + std::fabs(expected));
  }

  /// Compare the checkpoint matching @p completed_updates, if one is listed.
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

  /// Configured concrete model stored behind the public IonicModel interface.
  std::unique_ptr<IonicModel> model_;

  /// Validated trajectory inputs and reference-column mapping.
  IonicModelTrajectoryConfiguration configuration_;

  /// Reference CSV path resolved beneath UNIT_TEST_DATA_DIR.
  std::string reference_path_;

  /// Strictly increasing reference checkpoints loaded from the CSV.
  std::vector<Checkpoint> checkpoints_;
};

#endif // IONIC_MODEL_TEST_HELPERS_H
