// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ionic_bueno_orovio.h"
#include "ionic_model_test_helpers.h"
#include "gtest/gtest.h"

#include <cstddef>
#include <string>

namespace {

// The independent Forward-Euler oracles follow the equations of Bueno-Orovio,
// Cherry, and Fenton (2008). The M reference uses the current svMP tau_s2
// default documented in its CSV. A half-open pulse initiates one action
// potential from the published/default resting state without stretch current.
void run_bueno_orovio_trajectory(
    int zone_id,
    std::size_t update_count,
    const std::string &reference_csv_filename)
{
  BuenoOrovio::Parameters parameters;

  IonicModelTrajectoryConfiguration configuration;
  configuration.integration_type = TimeIntegrationType::FE;
  configuration.zone_id = zone_id;
  configuration.time_step = 0.01;
  configuration.update_count = update_count;
  configuration.reference_csv_filename = reference_csv_filename;
  configuration.state_reference_columns = {"V_mV", "v", "w", "s"};
  configuration.tolerance = 1.0e-10;
  configuration.stimulus_at_time = [](double time) {
    return 10.0 <= time && time < 12.0 ? -35.714 : 0.0;
  };
  configuration.sac_coefficient = 0.0;

  IonicModelTrajectoryTest<BuenoOrovio> trajectory(parameters, configuration);
  trajectory.run();
}

} // namespace

TEST(IonicModelTrajectory, BuenoOrovioEpi)
{
  run_bueno_orovio_trajectory(
      1, 60000, "ionic_bueno_orovio_epi_trajectory.csv");
}

TEST(IonicModelTrajectory, BuenoOrovioEndo)
{
  run_bueno_orovio_trajectory(
      2, 120000, "ionic_bueno_orovio_endo_trajectory.csv");
}

TEST(IonicModelTrajectory, BuenoOrovioM)
{
  run_bueno_orovio_trajectory(
      3, 120000, "ionic_bueno_orovio_m_trajectory.csv");
}
