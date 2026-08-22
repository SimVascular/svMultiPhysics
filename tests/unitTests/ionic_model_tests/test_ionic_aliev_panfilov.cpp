// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ionic_aliev_panfilov.h"
#include "ionic_model_test_helpers.h"
#include "gtest/gtest.h"

// The reference trajectory is an independent Forward-Euler evaluation of the
// split Aliev-Panfilov kinetics described by Goktepe and Kuhl (2009),
// doi:10.1002/nme.2571, using the svMultiPhysics default parameters and public
// voltage representation. Checkpoint N is the state after N completed updates.
TEST(IonicModelTrajectory, AlievPanfilov)
{
  AlievPanfilov::Parameters parameters;

  IonicModelTrajectoryConfiguration configuration;
  configuration.integration_type = TimeIntegrationType::FE;
  configuration.zone_id = 1;
  configuration.time_step = 0.129;
  configuration.update_count = 3200;
  configuration.initial_X_override = {{-30.0, 0.0}};
  configuration.reference_csv_filename =
      "ionic_aliev_panfilov_free_run.csv";
  configuration.state_reference_columns = {"V_mV", "w"};
  configuration.tolerance = 1.0e-10;
  configuration.stimulus = 0.0;
  configuration.sac_coefficient = 0.0;

  IonicModelTrajectoryTest<AlievPanfilov> trajectory(parameters,
                                                      configuration);
  trajectory.run();
}
