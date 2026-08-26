// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ionic_fitzhugh_nagumo.h"
#include "ionic_model_test_helpers.h"
#include "gtest/gtest.h"

#include <vector>

// The independent Forward-Euler oracle starts at the exact unstable
// equilibrium, applies a short positive stimulus, and follows the first
// excitation/repolarization cycle. It stops before the next autonomous
// upstroke; checkpoint N is the state after exactly N completed updates.
TEST(IonicModelTrajectory, FitzHughNagumo)
{
  FitzHughNagumo::Parameters parameters;

  IonicModelTrajectoryConfiguration configuration;
  configuration.integration_type = TimeIntegrationType::FE;
  configuration.zone_id = 1;
  configuration.time_step = 5.0e-4;
  configuration.update_count = 3000;
  configuration.initial_X_override = std::vector<double>{0.0, 0.0};
  configuration.reference_csv_filename =
      "ionic_fitzhugh_nagumo_fe_trajectory.csv";
  configuration.state_reference_columns = {"u", "w"};
  configuration.tolerance = 1.0e-12;
  configuration.stimulus_at_time = [](double time) {
    return 0.10 <= time && time < 0.12 ? 0.5 : 0.0;
  };
  configuration.sac_coefficient = 0.0;

  IonicModelTrajectoryTest<FitzHughNagumo> trajectory(parameters,
                                                       configuration);
  trajectory.run();
}
