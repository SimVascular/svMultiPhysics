// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "active_stress_regazzoni.h"
#include "active_stress_test_helpers.h"
#include "Vector.h"
#include "gtest/gtest.h"

#include <cmath>
#include <numbers>
#include <string>

// Each test starts from a fresh RegazzoniActiveStress model with the
// reference-calibration default parameters.
// All calls to node-local virtual methods go through ActiveStress& — required
// because the derived-class overrides remain protected.
class RegazzoniActiveStressTest : public ::testing::Test {
protected:
  RegazzoniActiveStress model_concrete;
  ActiveStress &model = model_concrete;
  Vector<double> state;

  void SetUp() override {
    RegazzoniActiveStress::Parameters params;
    model.read_model_parameters(params);  // through ActiveStress&
    state = Vector<double>(RegazzoniActiveStress::n_state_variables);
    model.init_local(state);                       // through ActiveStress&
  }
};

inline double scaled_tol(double expected, double rtol = 1e-10) {
  return rtol * (1.0 + std::fabs(expected));
}

// ---------------------------------------------------------------------------
// Twitch protocol (600 steps, dt = 1 ms)
// ---------------------------------------------------------------------------
// Ca: double-exponential transient. c0 = 1e-4 mM, cmax = 9e-4 mM,
//     tau_rise = 20 ms, tau_decay = 50 ms, onset at 10 ms.
// SL: raised-cosine ramp from SL0=2.2 µm to SL_min=2.134 µm and back.
//     Shortening onset at 30 ms, minimum at 150 ms, recovery complete at 350 ms.
// fiber_stretch  = SL(t_start) / SL0          (at the start of each step)
// fiber_stretch_rate = (SL(t_start+dt) - SL(t_start)) / (dt * SL0)  [1/ms]

namespace {

constexpr double dt      = 1.0;   // ms
constexpr int    N_steps = 600;

constexpr double c0      = 1.0e-4;  // mM
constexpr double cmax    = 9.0e-4;  // mM
constexpr double tau1    = 20.0;    // ms (rise)
constexpr double tau2    = 50.0;    // ms (decay)
constexpr double t0_ca   = 10.0;    // ms onset

constexpr double SL0     = 2.2;    // µm
constexpr double SL_min  = 2.134;  // µm
constexpr double t0_sl   = 30.0;   // ms onset
constexpr double t_min   = 150.0;  // ms at minimum SL
constexpr double t_rec   = 350.0;  // ms recovery complete

double peak_raw_ca_factor() {
  const double t_peak = std::log(tau1 / tau2) * tau1 * tau2 / (tau1 - tau2);
  return std::exp(-t_peak / tau2) - std::exp(-t_peak / tau1);
}

double calcium_at(double t) {
  if (t < t0_ca) return c0;
  const double s = t - t0_ca;
  const double raw = std::exp(-s / tau2) - std::exp(-s / tau1);
  static const double prf = peak_raw_ca_factor();
  return c0 + (cmax - c0) * raw / prf;
}

double SL_at(double t) {
  if (t < t0_sl) return SL0;
  if (t <= t_min) {
    const double frac = (t - t0_sl) / (t_min - t0_sl);
    return SL0 - (SL0 - SL_min) * (1.0 - std::cos(std::numbers::pi * frac)) / 2.0;
  }
  if (t <= t_rec) {
    const double frac = (t - t_min) / (t_rec - t_min);
    return SL_min + (SL0 - SL_min) * (1.0 - std::cos(std::numbers::pi * frac)) / 2.0;
  }
  return SL0;
}

double lam_at(double t)  { return SL_at(t) / SL0; }

} // namespace

// ---------------------------------------------------------------------------

TEST_F(RegazzoniActiveStressTest, Initialization) {
  EXPECT_NEAR(state[0], 1.0, 1e-15);
  for (unsigned int i = 1; i < RegazzoniActiveStress::n_state_variables; ++i)
    EXPECT_NEAR(state[i], 0.0, 1e-15) << "state[" << i << "]";
}

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
//   C++ discretization. The oracle uses the custom raised-cosine SL protocol
//   and forward-difference dSL/dt specified above; patch and unit-conversion
//   details are recorded in the CSV header.
// Values were NOT generated from svMP code.
//
// Active tension convention: compute_active_tension_local returns T_act [MPa].
// CSV column layout: step, Ta, s0, s1, ..., s19. States s0..s15 are the RU
// probabilities P(TL,TC,TR,CC), ordered by 8*TL + 4*TC + 2*TR + CC; s16..s19
// are [mu_P^0, mu_P^1, mu_N^0, mu_N^1], matching commit 26f05df.
// ---------------------------------------------------------------------------

TEST_F(RegazzoniActiveStressTest, TwitchTrajectory) {
  const std::string csv_path =
      std::string(UNIT_TEST_DATA_DIR) + "/active_stress_regazzoni_twitch.csv";
  const auto rows = load_csv(csv_path);
  ASSERT_FALSE(rows.empty()) << "No checkpoint rows in " << csv_path;

  // Each CSV row: [step, Ta, s0..s19] — 22 columns.
  constexpr size_t n_state = RegazzoniActiveStress::n_state_variables;
  constexpr size_t expected_cols = 2 + n_state;
  for (size_t r = 0; r < rows.size(); ++r)
    ASSERT_EQ(rows[r].size(), expected_cols)
        << "CSV row " << r << " has wrong column count";

  size_t ck_idx = 0;

  auto on_step = [&](int step, const Vector<double> &s, double Ta) {
    if (ck_idx < rows.size() && static_cast<int>(rows[ck_idx][0]) == step) {
      const double expected_Ta = rows[ck_idx][1];

      EXPECT_NEAR(Ta, expected_Ta, scaled_tol(expected_Ta))
          << "step " << step << ", active tension";

      for (unsigned int sv = 0; sv < n_state; ++sv)
        EXPECT_NEAR(s[sv], rows[ck_idx][2 + sv], scaled_tol(rows[ck_idx][2 + sv]))
            << "step " << step << ", state[" << sv << "]";

      ++ck_idx;
    }
  };

  run_active_stress_trajectory(
      model, state, N_steps, dt,
      [](int step) { return calcium_at(step * dt); },
      [](int step) { return lam_at(step * dt); },
      [](int step) {
        const double lam0 = lam_at(step * dt);
        const double lam1 = lam_at((step + 1) * dt);
        return (lam1 - lam0) / dt;
      },
      on_step);

  EXPECT_EQ(ck_idx, rows.size()) << "not all oracle checkpoints were reached";
}

// ---------------------------------------------------------------------------

TEST_F(RegazzoniActiveStressTest, RUProbabilityConservation) {
  run_active_stress_trajectory(
      model, state, N_steps, dt,
      [](int step) { return calcium_at(step * dt); },
      [](int step) { return lam_at(step * dt); },
      [](int step) {
        const double lam0 = lam_at(step * dt);
        const double lam1 = lam_at((step + 1) * dt);
        return (lam1 - lam0) / dt;
      },
      [](int step, const Vector<double> &s, double Ta) {
        for (unsigned int i = 0; i < RegazzoniActiveStress::n_state_variables; ++i)
          EXPECT_TRUE(std::isfinite(s[i]))
              << "step " << step << ", state[" << i << "]";
        EXPECT_TRUE(std::isfinite(Ta)) << "step " << step << ", active tension";

        double ru_sum = 0.0;
        for (unsigned int i = 0; i < RegazzoniActiveStress::n_ru_states; ++i)
          ru_sum += s[i];
        EXPECT_NEAR(ru_sum, 1.0, 1e-10) << "step " << step;
      });
}
