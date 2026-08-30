// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "active_stress_nash_panfilov.h"
#include "active_stress_test_helpers.h"
#include "gtest/gtest.h"

// ---------------------------------------------------------------------------
// References
//   [1] Nash & Panfilov (2004)
//       doi:10.1016/j.pbiomolbio.2004.01.016
//       Introduces the active-tension ODE in Eq. 22c, driven by a
//       non-dimensional excitation variable V, with the piecewise rate switch
//       in Eq. 23.
//   [2] Goktepe & Kuhl (2009)
//       doi:10.1007/s00466-009-0434-z
//       Eq. 46 reformulates the Nash-Panfilov ODE using dimensional
//       transmembrane potential Phi and resting potential Phi_r:
//         dT/dt = epsilon(Phi) * (k_sigma * (Phi - Phi_r) - T)
//       Eq. 47 replaces the piecewise switch from [1] with a smooth Gompertz
//       rate function:
//         epsilon(Phi) = eps0 + (eps_inf - eps0)*exp(-exp(-xi*(Phi - Phi_bar)))
//       Both equations are expressed in terms of transmembrane potential Phi,
//       not calcium; see Figure 3.
//
// svMultiPhysics follows the Goktepe-Kuhl form but substitutes intracellular
// calcium concentration Ca for Phi. The parameter-role mapping is
//   eta_T        <-> k_sigma       calcium_rest <-> Phi_r
//   xi_T         <-> xi            calcium_crit <-> Phi_bar.
// This calcium substitution is an svMultiPhysics-specific adaptation not
// documented in [1] or [2].
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// TwitchTrajectory
//
// Expected checkpoint values loaded from the external CSV oracle file:
//   tests/unitTests/reference_data/active_stress_nash_panfilov_twitch.csv
//
// Oracle: independent Python Forward-Euler evaluation of the svMultiPhysics
// calcium adaptation of Eqs. 46--47 in [2], with slab-calibration parameters.
// Values were NOT generated from svMP code.
//
// Active tension convention: compute_active_tension_local returns state[0].
// The existing CSV therefore uses its Ta column for both comparisons.
// ---------------------------------------------------------------------------

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
