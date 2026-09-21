// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

/// @file

#include "ActiveStressNashPanfilov.h"
#include "active_stress_test_helpers.h"
#include "gtest/gtest.h"

/**
 * @test Run a standalone ActiveStressNashPanfilov twitch for 200 updates with
 * @f$\Delta t=1\,\mathrm{ms}@f$, the prescribed calcium transient, and the
 * slab-calibration parameters configured below.
 *
 * The trusted reference uses libCellML-generated code adapted from the pinned
 * Physiome Nash--Panfilov 2004 CellML model. See the
 * [generator README](../reference_generators/active_stress/nash_panfilov/README.md)
 * for provenance and reproduction. Because active tension is the model's sole
 * state, the reference @c Ta column is used for both comparisons. See
 * @ref ActiveStressNashPanfilov for the model equations and calcium adaptation.
 */
TEST(ActiveStressTrajectory, ActiveStressNashPanfilov) {
  ActiveStressNashPanfilov::Parameters params;
  params.set_scalar("epsilon_0",    0.1);
  params.set_scalar("epsilon_i",    1.0);
  params.set_scalar("xi_T",         4.0e3);
  params.set_scalar("eta_T",        1.0e2);
  params.set_scalar("calcium_rest", 1.25e-4);
  params.set_scalar("calcium_crit", 8.0e-4);

  ActiveStressTrajectoryConfiguration configuration;
  configuration.final_time = 200.0;
  configuration.time_step = 1.0;
  configuration.reference_csv_filename =
      "active_stress_nash_panfilov_twitch.csv";
  configuration.state_reference_columns = {"Ta"};

  ActiveStressTrajectoryTest<ActiveStressNashPanfilov> trajectory(params,
                                                                   configuration);
  trajectory.run();
}
