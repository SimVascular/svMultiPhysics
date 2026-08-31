// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "active_stress_nash_panfilov.h"
#include "active_stress_test_helpers.h"
#include "gtest/gtest.h"

/**
 * @test Run a standalone NashPanfilov twitch for 200 updates with
 * @f$\Delta t=1\,\mathrm{ms}@f$, the prescribed calcium transient, and the
 * slab-calibration parameters configured below.
 *
 * The trusted reference is generated independently with Forward Euler by
 * @c reference_generators/active_stress/nash_panfilov/generate_nash_panfilov.py,
 * not by svMultiPhysics. Because active tension is the model's sole state, the
 * reference @c Ta column is used for both comparisons. See @ref NashPanfilov
 * for the model equations and calcium adaptation.
 */
TEST(ActiveStressTrajectory, NashPanfilov) {
  NashPanfilov::Parameters params;
  params.set_scalar("epsilon_0",    0.1);
  params.set_scalar("epsilon_i",    1.0);
  params.set_scalar("xi_T",         4.0e3);
  params.set_scalar("eta_T",        1.0e2);
  params.set_scalar("calcium_rest", 1.25e-4);
  params.set_scalar("calcium_crit", 8.0e-4);

  ActiveStressTrajectoryTest<NashPanfilov> trajectory(
      params, 200.0, 1.0, "active_stress_nash_panfilov_twitch.csv",
      {"Ta"});
  trajectory.run();
}
