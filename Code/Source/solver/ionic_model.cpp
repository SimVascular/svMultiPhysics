// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ionic_model.h"

#include "ComMod.h"
#include "Parameters.h"
#include "mat_fun.h"

#include <iostream>

const std::map<std::string, TimeIntegrationType> cep_time_int_to_type = {
    {"cn", TimeIntegrationType::CN2},
    {"cn2", TimeIntegrationType::CN2},
    {"implicit", TimeIntegrationType::CN2},
    {"fe", TimeIntegrationType::FE},
    {"euler", TimeIntegrationType::FE},
    {"explicit", TimeIntegrationType::FE},
    {"rk", TimeIntegrationType::RK4},
    {"rk4", TimeIntegrationType::RK4},
    {"runge", TimeIntegrationType::RK4}};

void IonicModel::read_parameters(const IonicModelParameters &params) {
  const auto &params_X = params.get_initial_X();
  for (auto &[label, value] : initial_X)
    value = params_X[label];

  const auto &params_Xg = params.get_initial_Xg();
  for (auto &[label, value] : initial_Xg)
    value = params_Xg[label];
}

void IonicModel::distribute_parameters(const CmMod &cm_mod, const cmType &cm) {
  for (size_t i = 0; i < initial_X.size(); ++i)
    cm.bcast(cm_mod, &initial_X[i].second);

  for (size_t i = 0; i < initial_Xg.size(); ++i)
    cm.bcast(cm_mod, &initial_Xg[i].second);
}

void IonicModel::init(Vector<double> &X, Vector<double> &Xg) const {
  if (initial_X.size() != X.size()) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
        "Initial conditions size for X does not match vector size.");
  }

  for (size_t i = 0; i < initial_X.size(); ++i)
    X[i] = initial_X[i].second;

  if (initial_Xg.size() != Xg.size()) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
        "Initial conditions size for Xg does not match vector size.");
  }

  for (size_t i = 0; i < initial_Xg.size(); ++i)
    Xg[i] = initial_Xg[i].second;
}

void IonicModel::integ(const odeType &ode_solver_params, const int zone_id,
                       const double t, const double dt, const double Istim,
                       const double Ksac, Vector<double> &X,
                       Vector<double> &Xg) const {
  switch (ode_solver_params.tIntType) {
  case TimeIntegrationType::FE:
    integ_fe(zone_id, X, Xg, t, dt, Istim, Ksac);
    break;

  case TimeIntegrationType::RK4:
    integ_rk(zone_id, X, Xg, t, dt, Istim, Ksac);
    break;

  case TimeIntegrationType::CN2:
    integ_cn2(zone_id, X, Xg, t, dt, Istim, Ksac, ode_solver_params.maxItr,
              ode_solver_params.relTol, ode_solver_params.absTol);
    break;

  default:
    svmp::raise<svmp::FE::InvalidArgumentException>(
        "Unknown time integration type: " +
            std::to_string(static_cast<int>(ode_solver_params.tIntType)));
  }
}

void IonicModel::integ_cn2(const unsigned int zone_id, Vector<double> &X,
                           Vector<double> &Xg, const double Ts, const double Ti,
                           const double Istim, const double Ksac,
                           const unsigned int max_iter, const double rtol,
                           const double atol) const {
  const unsigned int nX = X.size();
  const unsigned int nG = Xg.size();

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
  const Vector<double> fn = getf(zone_id, X, Xg, I_stim_scaled, I_sac_scaled);

  int k = 0;   // Current nonlinear iteration index.
  auto Xk = X; // Current solution.

  constexpr double eps = std::numeric_limits<double>::epsilon();

  while (true) {
    ++k;

    // Evaluate the right-hand side function for the system at the new time and
    // current nonlinear iteration.
    const Vector<double> fk =
        getf(zone_id, Xk, Xg, I_stim_scaled, I_sac_scaled);

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

    if (k > max_iter || rmsA <= atol || rmsR <= rtol)
      break;

    Array<double> JAC = getj(zone_id, Xk, Xg, Ksac * Tscale);

    JAC = Im - 0.5 * dt * JAC;
    JAC = mat_fun::mat_inv(JAC, nX);
    rK = mat_fun::mat_mul(JAC, rK);
    Xk = Xk - rK;
  }

  X = Xk;

  update_g(zone_id, dt, X, Xg);

  // Bring the potential variable back to dimensional units.
  X(0) = X(0) * Vscale + Voffset;
}

void IonicModel::integ_fe(const unsigned int zone_id, Vector<double> &X,
                          Vector<double> &Xg, const double Ts, const double Ti,
                          const double Istim, const double Ksac) const {
  const unsigned int nX = X.size();
  const unsigned int nG = Xg.size();

  // Rescale current time, timestep and transmembrane potential by the
  // model-specific scaling factors.
  const double t = Ts / Tscale;
  const double dt = Ti / Tscale;

  // Stretch-activated current.
  const double Isac = Ksac * (Vrest - X(0));

  X(0) = (X(0) - Voffset) / Vscale;

  const double I_stim_scaled = Istim * Tscale / Vscale;
  const double I_sac_scaled = Isac * Tscale / Vscale;

  const Vector<double> f = getf(zone_id, X, Xg, I_stim_scaled, I_sac_scaled);

  update_g(zone_id, dt, X, Xg);

  X = X + dt * f;

  // Bring the potential variable back to dimensional units.
  X(0) = X(0) * Vscale + Voffset;
}

void IonicModel::integ_rk(const unsigned int zone_id, Vector<double> &X,
                          Vector<double> &Xg, const double Ts, const double Ti,
                          const double Istim, const double Ksac) const {
  const unsigned int nX = X.size();
  const unsigned int nG = Xg.size();

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

  // First RK stage.
  Xrk = X;
  const Vector<double> frk1 =
      getf(zone_id, Xrk, Xg, I_stim_scaled, I_sac_scaled);

  // Update gating variables by half a step.
  auto Xgr = Xg;
  update_g(zone_id, 0.5 * dt, X, Xgr);

  // Second RK stage.
  Xrk = X + 0.5 * dt * frk1;
  const Vector<double> frk2 =
      getf(zone_id, Xrk, Xgr, I_stim_scaled, I_sac_scaled);

  // Third RK stage.
  Xrk = X + 0.5 * dt * frk2;
  const Vector<double> frk3 =
      getf(zone_id, Xrk, Xgr, I_stim_scaled, I_sac_scaled);

  // Update gating variables by the whole step.
  Xgr = Xg;
  update_g(zone_id, dt, X, Xgr);

  // Fourth RK stage.
  Xrk = X + dt * frk3;
  const Vector<double> frk4 =
      getf(zone_id, Xrk, Xgr, I_stim_scaled, I_sac_scaled);

  X = X + dt / 6.0 * (frk1 + 2.0 * (frk2 + frk3) + frk4);
  Xg = Xgr;

  // Bring the potential variable back to dimensional units.
  X(0) = X(0) * Vscale + Voffset;
}

std::vector<outputType> IonicModel::get_registered_outputs() const {
  std::vector<outputType> result;

  for (const auto &[name, index] : get_output_variables()) {
    outputType out;

    out.grp = consts::OutputNameType::outGrp_ionicState;
    out.name = name;
    out.o = index;
    out.l = 1;
    out.options.spatial = true;

    result.push_back(out);
  }

  return result;
}

std::unique_ptr<IonicModel>
IonicModelFactory::create_model(const std::string &name) {
  const auto &factory_instance = get_instance();

  auto iter = factory_instance.children.find(name);
  if (iter == factory_instance.children.end()) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
        "No model with name '" + name +
                       "' was registered in the ionic model factory.");
  }

  return iter->second();
}

void IonicModelFactory::visit(
    const std::function<void(const std::string &, const IonicModel &)> &f) {
  const auto &factory_instance = get_instance();

  for (auto &[name, builder] : factory_instance.children) {
    std::unique_ptr<IonicModel> dummy = builder();
    f(name, *dummy);
  }
}
