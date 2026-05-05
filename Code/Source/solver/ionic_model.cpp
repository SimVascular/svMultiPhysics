// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ionic_model.h"

#include "Parameters.h"
#include "mat_fun.h"

#include <iostream>

void IonicModel::set_initial_conditions(
    const IonicInitialConditionsParameters &params) {
  const auto &params_X = params.get_initial_X();
  for (auto &[label, value] : initial_X)
    value = params_X[label];
}

void IonicModel::distribute_initial_conditions(const CmMod &cm_mod,
                                               const cmType &cm) {
  for (size_t i = 0; i < initial_X.size(); ++i)
    cm.bcast(cm_mod, &initial_X[i].second);

  for (size_t i = 0; i < initial_Xg.size(); ++i)
    cm.bcast(cm_mod, &initial_Xg[i].second);
}

void IonicModel::init(const int nX, const int nG, Vector<double> &X,
                      Vector<double> &Xg) const {
  // @todo Replace with appropriate exception.
  if (initial_X.size() != X.size())
    throw std::runtime_error(
        "Initial conditions size for X does not match vector size.");

  for (size_t i = 0; i < initial_X.size(); ++i)
    X[i] = initial_X[i].second;

  if (initial_Xg.size() != Xg.size())
    throw std::runtime_error(
        "Initial conditions size for Xg does not match vector size.");

  for (size_t i = 0; i < initial_Xg.size(); ++i)
    Xg[i] = initial_Xg[i].second;
}

void IonicModel::integ_cn2(const unsigned int zone_id, const int nX,
                           const int nG, Vector<double> &X, Vector<double> &Xg,
                           const double Ts, const double Ti, const double Istim,
                           const double Ksac, Vector<int> &IPAR,
                           Vector<double> &RPAR) const {
  // @todo The nX and nG arguments can probably be removed and replaced by the
  // length of the state vectors X and Xg.

  const int itMax = IPAR(0);   // Maximum iterations for nonlinear solver.
  const double atol = RPAR(0); // Absolute tolerance for nonlinear solver.
  const double rtol = RPAR(1); // Relative tolerance for nonlinear solver.

  // Stretch-activated current.
  const double Isac = Ksac * (Vrest - X(0));

  // Rescale current time, timestep and transmembrane potential by the
  // model-specific scaling factors.
  const double dt = Ti / Tscale;
  const double t = Ts / Tscale + dt;
  X(0) = (X(0) - Voffset) / Vscale;

  const double I_stim_scaled = Istim * Tscale / Vscale;
  const double I_sac_scaled = Isac * Tscale / Vscale;

  const auto Im = mat_fun::mat_id(nX);

  // Evaluate the right-hand side function for the system at the old time.
  Vector<double> fn(nX);
  getf(zone_id, nX, nG, X, Xg, fn, I_stim_scaled, I_sac_scaled);

  int k = 0;    // Current nonlinear iteration index.
  auto Xk = X;  // Current solution. This copy is probably unnecessary.

  // @todo The following flags should be given meaningful names.
  // Flag indicating whether the maximum number of iterations was reached.
  bool l1 = false;

  // Flag indicating whether absolute tolerance is satisfied.
  bool l2 = false;

  // Flag indicating whether relative tolerance is satisfied.
  bool l3 = false;

  constexpr double eps = std::numeric_limits<double>::epsilon();

  while (true) {
    ++k;

    // Evaluate the right-hand side function for the system at the new time and
    // current nonlinear iteration.
    Vector<double> fk(nX);
    getf(zone_id, nX, nG, Xk, Xg, fk, I_stim_scaled, I_sac_scaled);

    auto rK = Xk - X - 0.5 * dt * (fk + fn);

    double rmsA = 0.0;
    double rmsR = 0.0;

    for (int i = 0; i < nX; ++i) {
      rmsA += rK(i) * rK(i);

      const double r_i = rK(i) / (Xk(i) + eps);
      rmsR += r_i * r_i;
    }

    rmsA = sqrt(rmsA / nX);
    rmsR = sqrt(rmsR / nX);

    l1 = (k > itMax);
    l2 = (rmsA <= atol);
    l3 = (rmsR <= rtol);

    if (l1 || l2 || l3)
      break;

    Array<double> JAC(nX, nX);
    getj(zone_id, nX, nG, Xk, Xg, JAC, Ksac * Tscale);

    JAC = Im - 0.5 * dt * JAC;
    JAC = mat_fun::mat_inv(JAC, nX);
    rK = mat_fun::mat_mul(JAC, rK);
    Xk = Xk - rK;
  }

  X = Xk;

  update_g(zone_id, dt, nX, nG, X, Xg);
  // @todo This call seems unnecessary, since fn is not used after this point.
  getf(zone_id, nX, nG, X, Xg, fn, I_stim_scaled, I_sac_scaled);

  // Bring the potential variable back to dimensional units.
  X(0) = X(0) * Vscale + Voffset;

  if (!l2 && !l3) {
    IPAR(1) = IPAR(1) + 1;
  }
}

void IonicModel::integ_fe(const unsigned int zone_id, const int nX,
                          const int nG, Vector<double> &X, Vector<double> &Xg,
                          const double Ts, const double Ti, const double Istim,
                          const double Ksac) const {

  // Rescale current time, timestep and transmembrane potential by the
  // model-specific scaling factors.
  const double t = Ts / Tscale;
  const double dt = Ti / Tscale;

  // Stretch-activated current.
  const double Isac = Ksac * (Vrest - X(0));

  X(0) = (X(0) - Voffset) / Vscale;

  const double I_stim_scaled = Istim * Tscale / Vscale;
  const double I_sac_scaled = Isac * Tscale / Vscale;

  Vector<double> f(nX);
  getf(zone_id, nX, nG, X, Xg, f, I_stim_scaled, I_sac_scaled);

  update_g(zone_id, dt, nX, nG, X, Xg);

  X = X + dt * f;

  // Bring the potential variable back to dimensional units.
  X(0) = X(0) * Vscale + Voffset;
}

void IonicModel::integ_rk(const unsigned int zone_id, const int nX,
                          const int nG, Vector<double> &X, Vector<double> &Xg,
                          const double Ts, const double Ti, const double Istim,
                          const double Ksac) const {
  // Stretch-activated current.
  const double Isac = Ksac * (Vrest - X(0));

  // Rescale current time, timestep and transmembrane potential by the
  // model-specific scaling factors.
  const double t = Ts / Tscale;
  const double dt = Ti / Tscale;
  X(0) = (X(0) - Voffset) / Vscale;

  const double I_stim_scaled = Istim * Tscale / Vscale;
  const double I_sac_scaled = Isac * Tscale / Vscale;

  Vector<double> Xrk(nX);
  Vector<double> frk1(nX), frk2(nX), frk3(nX), frk4(nX);

  // First RK stage.
  Xrk = X;
  getf(zone_id, nX, nG, Xrk, Xg, frk1, I_stim_scaled, I_sac_scaled);

  // Update gating variables by half a step.
  auto Xgr = Xg;
  update_g(zone_id, 0.5 * dt, nX, nG, X, Xgr);

  // Second RK stage.
  Xrk = X + 0.5 * dt * frk1;
  getf(zone_id, nX, nG, Xrk, Xgr, frk2, I_stim_scaled, I_sac_scaled);

  // Third RK stage.
  Xrk = X + 0.5 * dt * frk2;
  getf(zone_id, nX, nG, Xrk, Xgr, frk3, I_stim_scaled, I_sac_scaled);

  // Update gating variables by the whole step.
  Xgr = Xg;
  update_g(zone_id, dt, nX, nG, X, Xgr);

  // Fourth RK stage.
  Xrk = X + dt * frk3;
  getf(zone_id, nX, nG, Xrk, Xgr, frk4, I_stim_scaled, I_sac_scaled);

  X = X + dt / 6.0 * (frk1 + 2.0 * (frk2 + frk3) + frk4);
  Xg = Xgr;

  // Bring the potential variable back to dimensional units.
  X(0) = X(0) * Vscale + Voffset;
}