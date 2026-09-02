// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ionic_aliev_panfilov.h"
#include "ionic_model_test_helpers.h"
#include "gtest/gtest.h"

/**
 * @test Advance AlievPanfilov from its default state
 * @f$(-80\,\mathrm{mV},0.001)@f$ for 6000 Forward-Euler updates with
 * @f$\Delta t=0.1\,\mathrm{ms}@f$.
 *
 * The test applies @f$I_\mathrm{stim}=-35.714@f$ for
 * @f$10\leq t<12\,\mathrm{ms}@f$ and sets @f$K_\mathrm{sac}=0@f$. The trusted
 * reference is generated independently by
 * @c reference_generators/ionic_model/aliev_panfilov/generate_aliev_panfilov.py,
 * not by svMultiPhysics. See @ref AlievPanfilov for the equations and scaling.
 */
TEST(IonicModelTrajectory, AlievPanfilov)
{
  AlievPanfilov::Parameters parameters;

  IonicModelTrajectoryConfiguration configuration;
  configuration.integration_type = TimeIntegrationType::FE;
  configuration.zone_id = 1;
  configuration.time_step = 0.1;
  configuration.update_count = 6000;
  configuration.reference_csv_filename =
      "ionic_aliev_panfilov_stimulated_trajectory.csv";
  configuration.state_reference_columns = {"V_mV", "w"};
  configuration.tolerance = 1.0e-10;
  configuration.stimulus_at_time = [](double time) {
    return 10.0 <= time && time < 12.0 ? -35.714 : 0.0;
  };
  configuration.sac_coefficient = 0.0;

  IonicModelTrajectoryTest<AlievPanfilov> trajectory(parameters,
                                                      configuration);
  trajectory.run();
}
