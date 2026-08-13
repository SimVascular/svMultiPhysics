// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "active_stress_nash_panfilov.h"
#include "active_stress_test_helpers.h"
#include "Vector.h"
#include "gtest/gtest.h"

#include <cmath>
#include <string>

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

// Each test starts from a fresh NashPanfilov model with slab-calibration
// parameters configured through NashPanfilov::Parameters.
// All calls to node-local virtual methods go through ActiveStress& — required
// because the derived-class overrides remain protected.
class NashPanfilovActiveStressTest : public ::testing::Test {
protected:
  NashPanfilov model_concrete;
  ActiveStress &model = model_concrete;
  Vector<double> state;

  void SetUp() override {
    NashPanfilov::Parameters params;
    // All six NashPanfilov::Parameters defaults are 1.0.
    // Override with calibration values from
    // tests/cases/electromechanics/slab/solver_NashPanfilov.xml.
    // NashPanfilov::Parameters::set() sets value_ directly for full double
    // precision and marks value_set_ = true so check_required() passes.
    params.set("epsilon_0",    0.1);
    params.set("epsilon_i",    1.0);
    params.set("xi_T",         4.0e3);
    params.set("eta_T",        1.0e2);
    params.set("calcium_rest", 1.25e-4);
    params.set("calcium_crit", 8.0e-4);
    model.read_model_parameters(params);  // through ActiveStress&

    state = Vector<double>(1);
    model.init_local(state);                       // through ActiveStress&
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
// Forward Euler discretization of the svMultiPhysics calcium adaptation of
// Eqs. 46--47 in [2] has the exact solution
//
//   T^n = beta * (1 - (1 - alpha*dt)^n)
//
// where
//   alpha = epsilon(Ca_const)     [calcium-adapted Gompertz rate from [2]]
//   beta  = eta_T * (Ca_const - calcium_rest)    [equilibrium tension T_eq]
//
// Derived from the model equations independently of svMP code; any bug in
// getf produces disagreement with the closed form.
//
// The Gompertz epsilon and equilibrium tension are recomputed here from the
// fixture parameter values (read via the public getter get_scalar), so the
// oracle does not depend on accessing protected model members.
// ---------------------------------------------------------------------------

TEST_F(NashPanfilovActiveStressTest, ConstantCalciumEquilibrium) {
  // Read parameters back through the public Parameters getter to build the
  // oracle — no protected-member access needed.
  NashPanfilov::Parameters params;
  params.set("epsilon_0",    0.1);
  params.set("epsilon_i",    1.0);
  params.set("xi_T",         4.0e3);
  params.set("eta_T",        1.0e2);
  params.set("calcium_rest", 1.25e-4);
  params.set("calcium_crit", 8.0e-4);

  constexpr double Ca_const = cmax;  // 9e-4 mM: T_eq > 0, alpha*dt ≈ 0.56 (stable)

  const double eps0     = params.get_scalar("epsilon_0");
  const double epsi     = params.get_scalar("epsilon_i");
  const double xi_T_val = params.get_scalar("xi_T");
  const double eta_T_val= params.get_scalar("eta_T");
  const double ca_rest  = params.get_scalar("calcium_rest");
  const double ca_crit  = params.get_scalar("calcium_crit");

  const double alpha_ep = eps0 + (epsi - eps0)
      * std::exp(-std::exp(-xi_T_val * (Ca_const - ca_crit)));
  const double beta = eta_T_val * (Ca_const - ca_rest);

  double prev = state[0];

  for (int step = 0; step < 30; ++step) {
    const double t_start = step * dt;
    // Through ActiveStress& — required for public access.
    model.advance_time_step_local(t_start, dt, Ca_const, 1.0, 0.0, state);

    const int n = step + 1;
    const double T_exact = beta * (1.0 - std::pow(1.0 - alpha_ep * dt, n));

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
// Expected checkpoint values loaded from the external CSV oracle file:
//   tests/unitTests/reference_data/active_stress_nash_panfilov_twitch.csv
//
// Oracle: independent Python Forward-Euler evaluation of the svMultiPhysics
// calcium adaptation of Eqs. 46--47 in [2], with slab-calibration parameters.
// Values were NOT generated from svMP code.
//
// Active tension convention: compute_active_tension_local returns state[0]
// directly; both should agree to machine precision.
// ---------------------------------------------------------------------------

TEST_F(NashPanfilovActiveStressTest, TwitchTrajectory) {
  const std::string csv_path =
      std::string(UNIT_TEST_DATA_DIR) + "/active_stress_nash_panfilov_twitch.csv";
  const auto rows = load_csv(csv_path);
  ASSERT_FALSE(rows.empty()) << "No checkpoint rows in " << csv_path;

  size_t ck_idx = 0;

  auto on_step = [&](int step, const Vector<double> &s, double Ta) {
    if (ck_idx < rows.size() && static_cast<int>(rows[ck_idx][0]) == step) {
      const double expected_Ta = rows[ck_idx][1];

      EXPECT_NEAR(s[0], expected_Ta, scaled_tol(expected_Ta))
          << "step " << step << ", state[0]";
      EXPECT_NEAR(Ta, expected_Ta, scaled_tol(expected_Ta))
          << "step " << step << ", compute_active_tension_local";
      EXPECT_NEAR(Ta, s[0], 1e-15)
          << "step " << step << ", compute_active_tension_local == state[0]";

      ++ck_idx;
    }
  };

  run_active_stress_trajectory(
      model, state, N_steps, dt,
      [](int step) { return calcium_at(step * dt); },
      [](int /*step*/) { return 1.0; },
      [](int /*step*/) { return 0.0; },
      on_step);

  EXPECT_EQ(ck_idx, rows.size()) << "not all oracle checkpoints were reached";
}

// ---------------------------------------------------------------------------

TEST_F(NashPanfilovActiveStressTest, FiniteValues) {
  run_active_stress_trajectory(
      model, state, N_steps, dt,
      [](int step) { return calcium_at(step * dt); },
      [](int /*step*/) { return 1.0; },
      [](int /*step*/) { return 0.0; },
      [](int step, const Vector<double> &s, double Ta) {
        EXPECT_TRUE(std::isfinite(s[0])) << "step " << step << ", state[0]";
        EXPECT_TRUE(std::isfinite(Ta))   << "step " << step << ", active tension";
      });
}
