// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ionic_model_test_helpers.h"
#include "ionic_ttp.h"
#include "gtest/gtest.h"

#include <string>
#include <vector>

namespace {

class TTPTestParameters : public TTP::Parameters {
public:
  void set_parameter_value(const std::string &name, double value)
  {
    parameters.at(name).set(name, true, value);
  }
};

// The independent oracles follow the TP06 equations and public svMP
// FE/Rush-Larsen split. A half-open pulse initiates one beat from the
// published phenotype-specific initial state without stretch current.
void run_ttp_trajectory(
    const TTP::Parameters &parameters,
    int zone_id,
    const std::string &reference_csv_filename,
    const std::vector<double> &initial_X_override = {},
    const std::vector<double> &initial_Xg_override = {})
{
  IonicModelTrajectoryConfiguration configuration;
  configuration.integration_type = TimeIntegrationType::FE;
  configuration.zone_id = zone_id;
  configuration.time_step = 0.005;
  configuration.update_count = 120000;
  if (!initial_X_override.empty())
    configuration.initial_X_override = initial_X_override;
  if (!initial_Xg_override.empty())
    configuration.initial_Xg_override = initial_Xg_override;
  configuration.reference_csv_filename = reference_csv_filename;
  configuration.state_reference_columns = {
      "V_mV", "Ki",  "Nai", "Cai", "Ca_ss", "Ca_SR", "R_prime",
      "Xr1",  "Xr2", "Xs",  "m",   "h",     "j",     "d",
      "f",    "f2",  "fCass", "s",  "r"};
  configuration.tolerance = 1.0e-10;
  configuration.stimulus_at_time = [](double time) {
    return 10.0 <= time && time < 11.0 ? -52.0 : 0.0;
  };
  configuration.sac_coefficient = 0.0;

  IonicModelTrajectoryTest<TTP> trajectory(parameters, configuration);
  trajectory.run();
}

} // namespace

TEST(IonicModelTrajectory, TTPEpi)
{
  TTP::Parameters parameters;
  run_ttp_trajectory(parameters, 1, "ionic_ttp_epi_trajectory.csv");
}

// The ENDO oracle uses the published TP06 initial state and G_to conductance;
// zone 2 selects the endocardial s-gate kinetics through the public interface.
TEST(IonicModelTrajectory, TTPEndo)
{
  TTPTestParameters parameters;
  parameters.set_parameter_value("G_to", 0.073);
  run_ttp_trajectory(
      parameters, 2, "ionic_ttp_endo_trajectory.csv",
      {-86.709, 138.4, 10.355, 1.3e-4, 3.6e-4, 3.715, 0.9068},
      {0.00448, 0.476, 0.0087, 0.00155, 0.7573, 0.7225,
       3.164e-5, 0.8009, 0.9778, 0.9953, 0.3212, 2.235e-8});
}

// The M oracle uses the published TP06 initial state and G_Ks conductance;
// zone 3 selects the zone-dependent gate kinetics through the public interface.
TEST(IonicModelTrajectory, TTPM)
{
  TTPTestParameters parameters;
  parameters.set_parameter_value("G_Ks", 0.098);
  run_ttp_trajectory(
      parameters, 3, "ionic_ttp_m_trajectory.csv",
      {-85.423, 138.52, 10.132, 1.53e-4, 4.2e-4, 4.272, 0.8978},
      {0.0165, 0.473, 0.0174, 0.00165, 0.749, 0.6788,
       3.288e-5, 0.7026, 0.9526, 0.9942, 0.999998, 2.347e-8});
}
