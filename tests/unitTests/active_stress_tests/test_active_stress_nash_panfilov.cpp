// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "active_stress_nash_panfilov.h"
#include "Vector.h"
#include "gtest/gtest.h"

#include <cmath>

// ---------------------------------------------------------------------------
// References
//   [1] Nash & Panfilov (2004)
//       doi:10.1016/j.pbiomolbio.2004.01.016
//       Original Nash-Panfilov active stress ODE.
//   [2] Goktepe & Kuhl (2009)
//       doi:10.1007/s00466-009-0434-z
//       Gompertz sigmoid for epsilon(Ca); see Figure 3.
//
// The implementation uses the Nash-Panfilov active-stress ODE [1] with the
// Gompertz form for epsilon(Ca) from [2], using calcium as the activation
// variable in place of membrane potential u.
// ---------------------------------------------------------------------------

// Harness that exposes the protected node-local methods and parameter members for direct testing.
// No production code is modified.
struct NashPanfilovHarness : public NashPanfilov {
  using NashPanfilov::read_model_specific_parameters;
  using NashPanfilov::init_local;
  using NashPanfilov::advance_time_step_local;
  using NashPanfilov::compute_active_tension_local;
  // Expose protected parameter members so the fixture can assign slab values
  using NashPanfilov::epsilon_0;
  using NashPanfilov::epsilon_i;
  using NashPanfilov::xi_T;
  using NashPanfilov::eta_T;
  using NashPanfilov::calcium_rest;
  using NashPanfilov::calcium_crit;
};

// Each test starts from a fresh NashPanfilov model with slab-calibration parameters
// and a 1-state vector initialized by init_local().
class NashPanfilovActiveStressTest : public ::testing::Test {
protected:
  NashPanfilovHarness model;
  Vector<double> state;

  void SetUp() override {
    NashPanfilov::Parameters params;
    model.read_model_specific_parameters(params);

    // All six NashPanfilov::Parameters defaults are 1.0 (confirmed in the header).
    // Override with calibration values from
    // tests/cases/electromechanics/slab/solver_NashPanfilov.xml
    model.epsilon_0    = 0.1;
    model.epsilon_i    = 1.0;
    model.xi_T         = 4.0e3;
    model.eta_T        = 1.0e2;
    model.calcium_rest = 1.25e-4;
    model.calcium_crit = 8.0e-4;

    state = Vector<double>(1);
    model.init_local(state);
  }
};

inline double scaled_tol(double expected, double rtol = 1e-10) {
  return rtol * (1.0 + std::fabs(expected));
}

// ---------------------------------------------------------------------------
// Ca transient (200 steps, dt = 1 ms)
// ---------------------------------------------------------------------------
// Double-exponential: c0 = 1e-4 mM, cmax = 9e-4 mM,
//     tau_rise = 20 ms, tau_decay = 50 ms, onset at 10 ms.
// fiber_stretch = 1.0, fiber_stretch_rate = 0.0 (NashPanfilov ignores both).

namespace {

constexpr double dt      = 1.0;   // ms
constexpr int    N_steps = 200;

constexpr double c0      = 1.0e-4;   // mM
constexpr double cmax    = 9.0e-4;   // mM
constexpr double tau1    = 20.0;     // ms (rise)
constexpr double tau2    = 50.0;     // ms (decay)
constexpr double t0_ca   = 10.0;     // ms onset

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

} // namespace

// ---------------------------------------------------------------------------

TEST_F(NashPanfilovActiveStressTest, Initialization) {
  EXPECT_NEAR(state[0], 0.0, 1e-15);
}

// ---------------------------------------------------------------------------
// ConstantCalciumEquilibrium
//
// Closed-form discrete oracle (Oracle A): under constant Ca = Ca_const, the
// Forward Euler discretization of the Nash-Panfilov ODE [1], using the
// Gompertz epsilon(Ca) from [2], has the exact solution
//
//   T^n = beta * (1 - (1 - alpha*dt)^n)
//
// where
//   alpha = epsilon(Ca_const)     [Gompertz sigmoid [2] Fig. 3, at Ca_const]
//   beta  = eta_T * (Ca_const - calcium_rest)    [equilibrium tension T_eq]
//
// Derived from the model equations independently of svMP code; any bug in
// getf produces disagreement with the closed form.
// ---------------------------------------------------------------------------

TEST_F(NashPanfilovActiveStressTest, ConstantCalciumEquilibrium) {
  constexpr double Ca_const = cmax;  // 9e-4 mM: T_eq > 0, alpha*dt ≈ 0.56 (stable)

  const double alpha = model.epsilon_0
      + (model.epsilon_i - model.epsilon_0)
        * std::exp(-std::exp(-model.xi_T * (Ca_const - model.calcium_crit)));
  const double beta = model.eta_T * (Ca_const - model.calcium_rest);

  double prev = state[0];

  for (int step = 0; step < 30; ++step) {
    const double t_start = step * dt;
    model.advance_time_step_local(t_start, dt, Ca_const, 1.0, 0.0, state);

    const int n = step + 1;
    const double T_exact = beta * (1.0 - std::pow(1.0 - alpha * dt, n));

    EXPECT_NEAR(state[0], T_exact, scaled_tol(T_exact))
        << "step " << n << ", closed-form oracle";
    EXPECT_TRUE(std::isfinite(state[0])) << "step " << n;
    EXPECT_GE(state[0], 0.0) << "step " << n << " (T_eq > 0 => T_act non-negative)";
    EXPECT_LE(state[0], beta + scaled_tol(beta)) << "step " << n << " (no overshoot)";
    EXPECT_GE(state[0], prev) << "step " << n << " (non-decreasing toward T_eq)";
    prev = state[0];
  }
}

// ---------------------------------------------------------------------------
// TwitchTrajectory
//
// Expected checkpoint values were generated from an independent Python
// Forward-Euler implementation of the Nash-Panfilov active-stress ODE [1],
// using the Gompertz epsilon(Ca) from [2] and the slab-calibration parameters.
// The Ca transient is the double-exponential defined above. The Python script
// is not committed; checkpoint values are hardcoded here as the oracle.
//
// Active tension convention: compute_active_tension_local returns state[0]
// directly; the two should agree to machine precision.
// ---------------------------------------------------------------------------

TEST_F(NashPanfilovActiveStressTest, TwitchTrajectory) {
  struct Checkpoint { int step; double Ta; };

  const Checkpoint checkpoints[] = {
    // step 0 (t 0→1 ms; Ca < Ca_rest → T_act slightly negative)
    { 0, -2.5000016231668646e-04 },
    // step 10 (t 10→11 ms; Ca onset, T_act still rising from negative)
    { 10, -1.7154741323343946e-03 },
    // step 30 (t 30→31 ms; Ca rising, T_act strongly activated)
    { 30,  6.9827793067955765e-02 },
    // step 60 (t 60→61 ms; Ca and T_act both declining from peak)
    { 60,  6.8649748835108812e-02 },
    // step 99 (t 99→100 ms; Ca recovering, T_act declining)
    { 99,  4.0335891448303234e-02 },
    // step 149 (t 149→150 ms; late recovery)
    { 149, 1.5642696949229866e-02 },
    // step 199 (t 199→200 ms; near-baseline end)
    { 199, 4.3146770362862460e-03 },
  };

  int ck_idx = 0;
  const int n_checkpoints = static_cast<int>(sizeof(checkpoints) / sizeof(checkpoints[0]));

  for (int step = 0; step < N_steps; ++step) {
    const double t_start = step * dt;
    const double ca      = calcium_at(t_start);

    model.advance_time_step_local(t_start, dt, ca, 1.0, 0.0, state);

    if (ck_idx < n_checkpoints && checkpoints[ck_idx].step == step) {
      const Checkpoint &ck = checkpoints[ck_idx];
      const double Ta = model.compute_active_tension_local(state, 1.0);

      EXPECT_NEAR(state[0], ck.Ta, scaled_tol(ck.Ta))
          << "step " << step << ", state[0]";
      EXPECT_NEAR(Ta, ck.Ta, scaled_tol(ck.Ta))
          << "step " << step << ", compute_active_tension_local";
      EXPECT_NEAR(Ta, state[0], 1e-15)
          << "step " << step << ", compute_active_tension_local == state[0]";

      ++ck_idx;
    }
  }
  EXPECT_EQ(ck_idx, n_checkpoints) << "not all oracle checkpoints were reached";
}

// ---------------------------------------------------------------------------

TEST_F(NashPanfilovActiveStressTest, FiniteValues) {
  for (int step = 0; step < N_steps; ++step) {
    const double t_start = step * dt;
    const double ca      = calcium_at(t_start);

    model.advance_time_step_local(t_start, dt, ca, 1.0, 0.0, state);

    EXPECT_TRUE(std::isfinite(state[0])) << "step " << step << ", state[0]";
    EXPECT_TRUE(std::isfinite(model.compute_active_tension_local(state, 1.0)))
        << "step " << step << ", active tension";
  }
}
