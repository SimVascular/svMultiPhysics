// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ionic_fitzhugh_nagumo.h"
#include "ionic_model_test_helpers.h"
#include "gtest/gtest.h"

#include <vector>

/**
 * @test Advance FitzHughNagumo from the exact equilibrium
 * @f$(u,w)=(0,0)@f$ for 3000 Forward-Euler updates with
 * @f$\Delta t=0.0005@f$ in model-time units.
 *
 * The test applies @f$I_\mathrm{stim}=0.5@f$ for @f$0.10\leq t<0.12@f$ and
 * sets @f$K_\mathrm{sac}=0@f$. It follows the first triggered
 * excitation/repolarization cycle until @f$t=1.5@f$, before the next
 * autonomous upstroke. The trusted reference is generated independently by
 * @c reference_generators/ionic_model/fitzhugh_nagumo/generate_fitzhugh_nagumo.py,
 * not by svMultiPhysics. See @ref FitzHughNagumo for the model formulation.
 */
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
