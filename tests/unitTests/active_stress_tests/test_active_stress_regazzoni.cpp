// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "active_stress_regazzoni.h"
#include "active_stress_test_helpers.h"
#include "gtest/gtest.h"

// ---------------------------------------------------------------------------
// TwitchTrajectory
//
// Expected checkpoint values loaded from the external CSV oracle file:
//   tests/unitTests/reference_data/active_stress_regazzoni_twitch.csv
//
// Formulation source: Regazzoni, Dede', and Quarteroni (2020),
//   doi:10.1371/journal.pcbi.1008294. Its S3 numerical appendix uses Forward
//   Euler for RU dynamics and an exponential integrator for XB dynamics.
// Oracle discretization: patched authors' C++ reference implementation
//   (cardiac-activation, commit 26f05df), which instead uses Forward-Euler RU
//   substeps and implicit Euler for XB dynamics. svMultiPhysics follows that
//   C++ discretization. The oracle uses the common raised-cosine SL protocol
//   and forward-difference dSL/dt implemented by the generic runner; patch and
//   unit-conversion details are recorded in the CSV header.
// Values were NOT generated from svMP code.
//
// Active tension convention: compute_active_tension_local returns T_act [MPa].
// CSV column layout: step, Ta, s0, s1, ..., s19. States s0..s15 are the RU
// probabilities P(TL,TC,TR,CC), ordered by 8*TL + 4*TC + 2*TR + CC; s16..s19
// are [mu_P^0, mu_P^1, mu_N^0, mu_N^1], matching commit 26f05df.
// ---------------------------------------------------------------------------

TEST(ActiveStressTrajectory, Regazzoni) {
  RegazzoniActiveStress::Parameters params;
  ActiveStressTrajectoryTest<RegazzoniActiveStress> trajectory(
      params, 600.0, 1.0, "active_stress_regazzoni_twitch.csv");
  trajectory.run();
}
