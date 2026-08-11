// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "active_stress_regazzoni.h"

#include "eigen3/Eigen/Dense"

#include <algorithm>
#include <cmath>

void RegazzoniActiveStress::read_model_specific_parameters(
    const ActiveStressModelParameters &params) {
  Kbasic = params.get_scalar("Kbasic");
  Koff = params.get_scalar("Koff");
  Q = params.get_scalar("Q");
  mu = params.get_scalar("mu");
  gamma = params.get_scalar("gamma");
  Kd0 = params.get_scalar("Kd0");
  alphaKd = params.get_scalar("alphaKd");
  if (alphaKd > 0.0)
    svmp::raise<svmp::ParseException>(
        "RegazzoniActiveStress: alphaKd must be <= 0 (positive values reduce calcium "
        "sensitivity with stretch, reversing length-dependent activation, "
        "and can produce a zero dissociation constant at physiological "
        "sarcomere lengths).");
  SL0 = params.get_scalar("SL0");
  ru_substep = params.get_scalar("ru_substep");
  kd_reference_sarcomere_length = params.get_scalar("kd_reference_sarcomere_length");

  r0 = params.get_scalar("r0");
  alpha = params.get_scalar("alpha");
  mu0_fP = params.get_scalar("mu0_fP");
  mu1_fP = params.get_scalar("mu1_fP");

  LA = params.get_scalar("LA");
  LM = params.get_scalar("LM");
  LB = params.get_scalar("LB");
  a_XB = params.get_scalar("a_XB");

  disable_force_strain_rate_feedback_ =
      params.get_bool("Disable_force_strain_rate_feedback");
}

void RegazzoniActiveStress::distribute_model_specific_parameters(
    const CmMod &cm_mod, const cmType &cm) {
  cm.bcast(cm_mod, &Kbasic);
  cm.bcast(cm_mod, &Koff);
  cm.bcast(cm_mod, &Q);
  cm.bcast(cm_mod, &mu);
  cm.bcast(cm_mod, &gamma);
  cm.bcast(cm_mod, &Kd0);
  cm.bcast(cm_mod, &alphaKd);
  cm.bcast(cm_mod, &SL0);
  cm.bcast(cm_mod, &ru_substep);
  cm.bcast(cm_mod, &kd_reference_sarcomere_length);

  cm.bcast(cm_mod, &r0);
  cm.bcast(cm_mod, &alpha);
  cm.bcast(cm_mod, &mu0_fP);
  cm.bcast(cm_mod, &mu1_fP);

  cm.bcast(cm_mod, &LA);
  cm.bcast(cm_mod, &LM);
  cm.bcast(cm_mod, &LB);
  cm.bcast(cm_mod, &a_XB);
  cm.bcast(cm_mod, &disable_force_strain_rate_feedback_);
}

void RegazzoniActiveStress::init_local(Vector<double> &state) const {
  for (unsigned int i = 0; i < n_state_variables; ++i)
    state[i] = 0.0;

  state[ru_index(0, 0, 0, 0)] = 1.0; // == state[0]
}

void RegazzoniActiveStress::advance_time_step_local(
    const double t, const double dt, const double calcium,
    const double fiber_stretch, const double fiber_stretch_rate,
    Vector<double> &state) const {
  const double sarcomere_length = SL0 * fiber_stretch;

  // Calcium/stretch-independent central-tropomyosin transition rates.
  const RUArray rates_T = ru_transition_rates_tropomyosin();

  // Troponin transition rates rates_C[CC][TC]: the calcium-binding row (CC = 0)
  // depends on calcium and sarcomere length; the unbinding row (CC = 1) does
  // not.
  const double calcium_on_rate =
      Koff /
      (Kd0 - alphaKd * (kd_reference_sarcomere_length - sarcomere_length)) *
      calcium;
  BinaryPairArray rates_C;
  rates_C[0][0] = calcium_on_rate;
  rates_C[0][1] = calcium_on_rate;
  rates_C[1][0] = Koff;
  rates_C[1][1] = Koff / mu;

  // Deserialize the 16 RU probabilities (entries 0-15). The crossbridge moments
  // (entries 16-19) are left untouched by this increment.
  RUArray state_RU;
  for (int TL = 0; TL < 2; ++TL)
    for (int TC = 0; TC < 2; ++TC)
      for (int TR = 0; TR < 2; ++TR)
        for (int CC = 0; CC < 2; ++CC)
          state_RU[TL][TC][TR][CC] = state[ru_index(TL, TC, TR, CC)];

  // Forward-Euler substepping over the outer time step. The final substep is
  // shortened so that the outer step is covered exactly.
  double time_advanced = 0.0;
  while (time_advanced <= dt - 1.0e-10) {
    const double substep = std::min(ru_substep, dt - time_advanced);
    ru_forward_euler_substep(substep, rates_T, rates_C, state_RU);
    time_advanced += substep;
  }

  // Advance the crossbridge moments (entries 16-19) from the updated RU state.
  // The velocity v = -dSL/dt / SL0 reduces to -d(lambda)/dt because SL = SL0 * lambda.
  const double velocity =
      disable_force_strain_rate_feedback_ ? 0.0 : -fiber_stretch_rate;
  XBArray state_XB;
  for (int i = 0; i < 4; ++i)
    state_XB[i] = state[xb_index(i)];
  state_XB = xb_implicit_update(dt, velocity, rates_T, state_RU, state_XB);

  // Serialize the updated RU probabilities back into the state vector.
  for (int TL = 0; TL < 2; ++TL)
    for (int TC = 0; TC < 2; ++TC)
      for (int TR = 0; TR < 2; ++TR)
        for (int CC = 0; CC < 2; ++CC)
          state[ru_index(TL, TC, TR, CC)] = state_RU[TL][TC][TR][CC];

  // Serialize the updated crossbridge moments back into the state vector.
  for (int i = 0; i < 4; ++i)
    state[xb_index(i)] = state_XB[i];
}

double RegazzoniActiveStress::compute_active_tension_local(
    const Vector<double> &state, const double fiber_stretch) const {
  const double sarcomere_length = SL0 * fiber_stretch;

  // Active tension T_act = a_XB * (μ_P^1 + μ_N^1) * φ(SL) from the
  // permissive and non-permissive XB first moments (state entries 17 and 19),
  // scaled by the single-overlap fraction and the upscaling factor a_XB.
  return a_XB * (state[xb_index(1)] + state[xb_index(3)]) *
         fraction_single_overlap(sarcomere_length);
}

RegazzoniActiveStress::RUArray
RegazzoniActiveStress::ru_transition_rates_tropomyosin() const {
  RUArray rates_T;
  for (int TL = 0; TL < 2; ++TL)
    for (int TR = 0; TR < 2; ++TR) {
      const int permissive_neighbors = TL + TR;

      // Rate of leaving the permissive central state (TC = 1).
      const double closing_rate =
          Kbasic * std::pow(gamma, 2 - permissive_neighbors);
      // Rate of leaving the non-permissive central state (TC = 0).
      const double opening_rate =
          Q * Kbasic * std::pow(gamma, permissive_neighbors);

      rates_T[TL][1][TR][0] = closing_rate;
      rates_T[TL][1][TR][1] = closing_rate;
      rates_T[TL][0][TR][0] = opening_rate / mu;
      rates_T[TL][0][TR][1] = opening_rate;
    }
  return rates_T;
}

void RegazzoniActiveStress::ru_forward_euler_substep(
    double dt, const RUArray &rates_T,
    const BinaryPairArray &rates_C, RUArray &state_RU) const {
  // Probability fluxes from central-unit transitions.
  RUArray flux_TC; // central tropomyosin
  RUArray flux_CC; // central troponin
  for (int TL = 0; TL < 2; ++TL)
    for (int TC = 0; TC < 2; ++TC)
      for (int TR = 0; TR < 2; ++TR)
        for (int CC = 0; CC < 2; ++CC) {
          flux_TC[TL][TC][TR][CC] =
              state_RU[TL][TC][TR][CC] * rates_T[TL][TC][TR][CC];
          flux_CC[TL][TC][TR][CC] =
              state_RU[TL][TC][TR][CC] * rates_C[CC][TC];
        }

  // Effective transition rates of the boundary neighbours, obtained from the
  // mean-field closure by conditioning the central-unit flux on the neighbour
  // pair state.
  BinaryPairArray rate_left;
  BinaryPairArray rate_right;
  for (int TL = 0; TL < 2; ++TL)
    for (int TC = 0; TC < 2; ++TC) {
      double flux_sum = 0.0;
      double prob_sum = 0.0;
      for (int TR = 0; TR < 2; ++TR)
        for (int CC = 0; CC < 2; ++CC) {
          flux_sum += flux_TC[TL][TC][TR][CC];
          prob_sum += state_RU[TL][TC][TR][CC];
        }
      rate_left[TL][TC] = (prob_sum > 1.0e-12) ? flux_sum / prob_sum : 0.0;
    }
  for (int TR = 0; TR < 2; ++TR)
    for (int TC = 0; TC < 2; ++TC) {
      double flux_sum = 0.0;
      double prob_sum = 0.0;
      for (int TL = 0; TL < 2; ++TL)
        for (int CC = 0; CC < 2; ++CC) {
          flux_sum += flux_TC[TL][TC][TR][CC];
          prob_sum += state_RU[TL][TC][TR][CC];
        }
      rate_right[TR][TC] = (prob_sum > 1.0e-12) ? flux_sum / prob_sum : 0.0;
    }

  // Probability fluxes from the boundary-neighbour transitions.
  // TL's only neighbour is TC on its right  → rate_right[TC][TL].
  // TR's only neighbour is TC on its left   → rate_left[TC][TR].
  // (rate_left == rate_right numerically due to mean-field LR symmetry, so the
  // result is unchanged, but the names now match the physical convention.)
  RUArray flux_TL; // left tropomyosin
  RUArray flux_TR; // right tropomyosin
  for (int TL = 0; TL < 2; ++TL)
    for (int TC = 0; TC < 2; ++TC)
      for (int TR = 0; TR < 2; ++TR)
        for (int CC = 0; CC < 2; ++CC) {
          flux_TL[TL][TC][TR][CC] =
              state_RU[TL][TC][TR][CC] * rate_right[TC][TL];
          flux_TR[TL][TC][TR][CC] =
              state_RU[TL][TC][TR][CC] * rate_left[TC][TR];
        }

  // Forward-Euler update of the 16 RU probabilities.
  for (int TL = 0; TL < 2; ++TL)
    for (int TC = 0; TC < 2; ++TC)
      for (int TR = 0; TR < 2; ++TR)
        for (int CC = 0; CC < 2; ++CC)
          state_RU[TL][TC][TR][CC] +=
              dt * (-flux_TL[TL][TC][TR][CC] + flux_TL[1 - TL][TC][TR][CC] -
                    flux_TC[TL][TC][TR][CC] + flux_TC[TL][1 - TC][TR][CC] -
                    flux_TR[TL][TC][TR][CC] + flux_TR[TL][TC][1 - TR][CC] -
                    flux_CC[TL][TC][TR][CC] + flux_CC[TL][TC][TR][1 - CC]);
}

RegazzoniActiveStress::XBArray RegazzoniActiveStress::xb_implicit_update(
    double dt, double velocity,
    const RUArray &rates_T,
    const RUArray &state_RU,
    const XBArray &state_XB) const {
  // Permissivity and the permissive/non-permissive probability fluxes from the
  // updated RU state.
  double permissivity = 0.0;
  double flux_PN = 0.0;
  double flux_NP = 0.0;
  for (int TL = 0; TL < 2; ++TL)
    for (int TR = 0; TR < 2; ++TR)
      for (int CC = 0; CC < 2; ++CC) {
        permissivity += state_RU[TL][1][TR][CC];
        flux_PN += state_RU[TL][1][TR][CC] * rates_T[TL][1][TR][CC];
        flux_NP += state_RU[TL][0][TR][CC] * rates_T[TL][0][TR][CC];
      }

  // Effective permissive->non-permissive and non-permissive->permissive rates.
  const double k_PN = (permissivity >= 1.0e-12) ? flux_PN / permissivity : 0.0;
  const double k_NP =
      ((1.0 - permissivity) >= 1.0e-12) ? flux_NP / (1.0 - permissivity) : 0.0;

  // Use the calibrated RDQ20-MF specialization of the general XB system:
  // new XBs attach only in the permissive state (f_N = 0), and both XB
  // populations share r(v) = r0 + alpha * |v|. Non-permissive moments
  // are populated by P-to-N transitions of already-attached XBs.
  const double r = r0 + alpha * std::abs(velocity);
  const double diag_P = r + k_PN;
  const double diag_N = r + k_NP;

  // Implicit-Euler system (I - dt * A) x = rhs for the four moments. The matrix
  // is zero-initialized so the structurally-zero entries are correct.
  Eigen::Matrix<double, 4, 4> system = Eigen::Matrix<double, 4, 4>::Zero();
  system(0, 0) = -diag_P;
  system(1, 1) = -diag_P;
  system(2, 2) = -diag_N;
  system(3, 3) = -diag_N;
  system(0, 2) = k_NP;
  system(1, 3) = k_NP;
  system(2, 0) = k_PN;
  system(3, 1) = k_PN;
  system(1, 0) = -velocity;
  system(3, 2) = -velocity;
  system *= -dt;
  for (int i = 0; i < 4; ++i)
    system(i, i) += 1.0;

  Eigen::Matrix<double, 4, 1> rhs;
  rhs(0) = state_XB[0] + dt * permissivity * mu0_fP;
  rhs(1) = state_XB[1] + dt * permissivity * mu1_fP;
  rhs(2) = state_XB[2];
  rhs(3) = state_XB[3];

  const Eigen::Matrix<double, 4, 1> solution =
      system.colPivHouseholderQr().solve(rhs);
  XBArray result;
  for (int i = 0; i < 4; ++i)
    result[i] = solution(i);
  return result;
}

double RegazzoniActiveStress::fraction_single_overlap(double sarcomere_length) const {
  const double SL = sarcomere_length;
  const double half_single_overlap = (LM - LB) * 0.5;

  if (SL > LA && SL <= LM)
    return (SL - LA) / half_single_overlap;
  if (SL > LM && SL <= 2.0 * LA - LB)
    return (SL + LM - 2.0 * LA) * 0.5 / half_single_overlap;
  if (SL > 2.0 * LA - LB && SL <= 2.0 * LA + LB)
    return 1.0;
  if (SL > 2.0 * LA + LB && SL <= 2.0 * LA + LM)
    return (LM + 2.0 * LA - SL) * 0.5 / half_single_overlap;
  return 0.0;
}

REGISTER_ACTIVE_STRESS_MODEL("Regazzoni", RegazzoniActiveStress);
