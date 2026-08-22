// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ionic_fitzhugh_nagumo.h"
#include "ionic_model_test_helpers.h"
#include "gtest/gtest.h"

// The external oracle is an independently validated literal Forward-Euler
// trajectory for the default svMultiPhysics FitzHugh-Nagumo parameters. Its
// checkpoints include the discrete extrema of u and w; checkpoint N is the
// state after exactly N completed updates.
TEST(IonicModelTrajectory, FitzHughNagumo)
{
  FitzHughNagumo::Parameters parameters;

  IonicModelTrajectoryConfiguration configuration;
  configuration.integration_type = TimeIntegrationType::FE;
  configuration.zone_id = 1;
  configuration.time_step = 1.0e-3;
  configuration.update_count = 2000;
  configuration.reference_csv_filename =
      "ionic_fitzhugh_nagumo_fe_trajectory.csv";
  configuration.state_reference_columns = {"u", "w"};
  configuration.tolerance = 1.0e-12;
  configuration.stimulus = 0.0;
  configuration.sac_coefficient = 0.0;

  IonicModelTrajectoryTest<FitzHughNagumo> trajectory(parameters,
                                                       configuration);
  trajectory.run();
}
