// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "active_stress_regazzoni.h"
#include "Vector.h"
#include "gtest/gtest.h"

#include <cmath>
#include <numbers>

// Harness that exposes the protected node-local methods for direct testing.
// No production code is modified.
struct RegazzoniHarness : public RegazzoniActiveStress {
  using RegazzoniActiveStress::read_model_specific_parameters;
  using RegazzoniActiveStress::init_local;
  using RegazzoniActiveStress::advance_time_step_local;
  using RegazzoniActiveStress::compute_active_tension_local;
};

// Each test starts from a fresh Regazzoni model with the reference-calibration
// parameters and a 20-state vector initialized by init_local().
class RegazzoniActiveStressTest : public ::testing::Test {
protected:
  RegazzoniHarness model;
  Vector<double> state;

  void SetUp() override {
    RegazzoniActiveStress::Parameters params;
    model.read_model_specific_parameters(params);
    state = Vector<double>(RegazzoniActiveStress::n_state_variables);
    model.init_local(state);
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
// Expected checkpoint values were generated from the patched Regazzoni C++
// reference implementation (cardiac-activation, commit 26f05df) using the
// twitch protocol above.
//
// The test below compares the svMultiPhysics RegazzoniActiveStress trajectory
// against these independent reference values using a scaled tolerance of 1e-10.
//
// Active tension convention: compute_active_tension_local returns T_act [MPa]
// directly. Oracle conversion:
//   Ta_svMP [MPa] = Ta_oracle [kPa] / 1000
// ---------------------------------------------------------------------------

TEST_F(RegazzoniActiveStressTest, TwitchTrajectory) {
  struct Checkpoint {
    int    step;
    double Ta;
    double state[RegazzoniActiveStress::n_state_variables];
  };

  const Checkpoint checkpoints[] = {
    // step 0 (after first advance, t 0→1 ms)
    { 0,
      /* Ta [MPa] */ 2.0769853711984563e-05,
      { 9.6940717188151537e-01, 2.6381424371446922e-02, 1.3503189084942508e-03,
        3.4343427180013618e-05, 1.1236588180500891e-03, 2.1053600122563298e-04,
        4.5470171062514682e-05, 6.7945663780045756e-06, 1.3503189084942508e-03,
        3.4343427180013618e-05, 1.7672866596113026e-06, 2.9934382365957969e-08,
        4.5470171062514682e-05, 6.7945663780045756e-06, 1.3909691878929928e-06,
        1.6659130249734328e-07, 1.6359975630299202e-05, 3.8979760023191681e-07,
        2.5100933115804985e-05, 5.9806222901712789e-07 } },
    // step 30 (t 30→31 ms, Ca rising)
    { 30,
      /* Ta [MPa] */ 6.3170372259402953e-03,
      { 2.6970273699972086e-01, 5.8583076955075841e-01, 9.0822857894029864e-03,
        1.7508306386496711e-02, 3.3942021959741897e-04, 7.2968975122089169e-03,
        9.2499470591179186e-04, 1.8566812618242029e-02, 9.0822857894029864e-03,
        1.7508306386496711e-02, 2.5639379411841935e-04, 2.8113908594203660e-04,
        9.2499470591179197e-04, 1.8566812618242029e-02, 2.1651243065998749e-03,
        4.1962719530947941e-02, 6.4443780653849255e-03, 1.5352064893374559e-04,
        6.1681286868829335e-03, 1.4693150647919565e-04 } },
    // step 99 (t 99→100 ms, near peak tension, shortening)
    { 99,
      /* Ta [MPa] */ 5.2391374181279100e-02,
      { 1.5732019324453009e-01, 2.0979293676465199e-01, 1.5974870723815159e-02,
        1.9784652352790824e-02, 2.0467839599633408e-04, 2.7398599272003450e-03,
        2.4988783785102964e-03, 3.3424626830984350e-02, 1.5974870723815156e-02,
        1.9784652352790824e-02, 1.4617914198968632e-03, 1.6467290361884674e-03,
        2.4988783785102964e-03, 3.3424626830984350e-02, 3.2887796467786577e-02,
        4.5057995817154656e-01, 1.0031077325241487e-01, 2.1689572017283237e-03,
        2.0243620609812121e-02, 4.0035869804356313e-04 } },
    // step 157 (t 157→158 ms, recovery beginning)
    { 157,
      /* Ta [MPa] */ 6.7187773320350211e-02,
      { 2.1301527707148393e-01, 1.3625459545317387e-01, 2.0337729296022237e-02,
        1.3890738926677160e-02, 3.0659349196029512e-04, 1.9921891163431895e-03,
        4.5754093405092769e-03, 3.0894774488222448e-02, 2.0337729296022233e-02,
        1.3890738926677160e-02, 2.0494138051662127e-03, 1.4910844091695922e-03,
        4.5754093405092769e-03, 3.0894774488222448e-02, 6.5139706850121340e-02,
        4.4035383569971798e-01, 1.1765751937916291e-01, 2.7943518198232085e-03,
        2.4046359000459934e-02, 5.6498445537614308e-04 } },
    // step 299 (t 299→300 ms, late recovery)
    { 299,
      /* Ta [MPa] */ 2.5267281291368970e-02,
      { 5.7786285533782400e-01, 1.8384587957630338e-01, 1.8534599275919040e-02,
        6.9531290368771486e-03, 8.4191991802026170e-04, 2.7210981345073609e-03,
        5.1543170439206937e-03, 1.7849846441656804e-02, 1.8534599275919040e-02,
        6.9531290368771486e-03, 7.4140640985424593e-04, 3.3804681545467166e-04,
        5.1543170439206937e-03, 1.7849846441656804e-02, 3.1002750836783124e-02,
        1.0566225937450435e-01, 3.0689769447400983e-02, 7.6141392250467059e-04,
        1.7209805063539828e-02, 4.4935336072742248e-04 } },
    // step 599 (t 599→600 ms, end of twitch)
    { 599,
      /* Ta [MPa] */ 2.9208300264620501e-03,
      { 7.4515028798360583e-01, 2.1148442974808412e-01, 7.5213053864882388e-03,
        2.1353736195492204e-03, 1.0350098623700802e-03, 2.9375949455035725e-03,
        1.5060379230515584e-03, 4.2755884037494587e-03, 7.5213053864882388e-03,
        2.1353736195492204e-03, 7.5955625924313899e-05, 2.1596502712894121e-05,
        1.5060379230515584e-03, 4.2755884037494587e-03, 2.1927274111406695e-03,
        6.2257872549735345e-03, 1.5693229395868540e-03, 3.7391150797739044e-05,
        4.2612557861625277e-03, 1.0152993604368515e-04 } },
  };

  int ck_idx = 0;
  const int n_checkpoints = sizeof(checkpoints) / sizeof(checkpoints[0]);

  for (int step = 0; step < N_steps; ++step) {
    const double t_start = step * dt;
    const double lam     = lam_at(t_start);
    const double dlam_dt = (lam_at(t_start + dt) - lam) / dt;
    const double ca      = calcium_at(t_start);

    model.advance_time_step_local(t_start, dt, ca, lam, dlam_dt, state);

    if (ck_idx < n_checkpoints && checkpoints[ck_idx].step == step) {
      const Checkpoint &ck = checkpoints[ck_idx];
      const double Ta = model.compute_active_tension_local(state, lam);

      for (unsigned int s = 0; s < RegazzoniActiveStress::n_state_variables; ++s)
        EXPECT_NEAR(state[s], ck.state[s], scaled_tol(ck.state[s]))
          << "step " << step << ", state[" << s << "]";

      EXPECT_NEAR(Ta, ck.Ta, scaled_tol(ck.Ta))
        << "step " << step << ", active tension";

      ++ck_idx;
    }
  }
}

// ---------------------------------------------------------------------------

TEST_F(RegazzoniActiveStressTest, RUProbabilityConservation) {
  for (int step = 0; step < N_steps; ++step) {
    const double t_start = step * dt;
    const double lam     = lam_at(t_start);
    const double dlam_dt = (lam_at(t_start + dt) - lam) / dt;
    const double ca      = calcium_at(t_start);

    model.advance_time_step_local(t_start, dt, ca, lam, dlam_dt, state);

    for (unsigned int i = 0; i < RegazzoniActiveStress::n_state_variables; ++i)
      EXPECT_TRUE(std::isfinite(state[i])) << "step " << step << ", state[" << i << "]";

    const double Ta = model.compute_active_tension_local(state, lam);
    EXPECT_TRUE(std::isfinite(Ta)) << "step " << step << ", active tension";

    double ru_sum = 0.0;
    for (unsigned int i = 0; i < RegazzoniActiveStress::n_ru_states; ++i)
      ru_sum += state[i];
    EXPECT_NEAR(ru_sum, 1.0, 1e-10) << "step " << step;
  }
}
