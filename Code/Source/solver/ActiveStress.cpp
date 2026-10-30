// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ActiveStress.h"

#include "mat_fun.h"

bool supports_active_stress(const consts::EquationType eq_type) {
  return eq_type == consts::EquationType::phys_struct ||
         eq_type == consts::EquationType::phys_ustruct ||
         eq_type == consts::EquationType::phys_FSI;
}

void ActiveStress::Evaluator::update(const ActiveStress &active_stress,
                                     const Vector<int> &nodes) {
  active_stress_ = &active_stress;

  const unsigned int n_states = active_stress.n_states;

  state_.resize(n_states, nodes.size());

  // Vector::resize reallocates even when the size is unchanged.
  if (interpolated_state.size() != n_states)
    interpolated_state.resize(n_states);

  // Friend access to active_stress.states, so that gathering the state of an
  // element does not need to go through an accessor.
  for (int a = 0; a < nodes.size(); ++a)
    for (unsigned int j = 0; j < n_states; ++j)
      state_(j, a) = active_stress.states(j, nodes(a));
}

template <int nsd>
ActiveStress::ActiveTension ActiveStress::Evaluator::evaluate(
    const Eigen::Ref<const Eigen::VectorXd> &N, const mat_fun::Matrix<nsd> &F,
    const Eigen::Ref<const Eigen::Matrix<double, nsd, Eigen::Dynamic>> &fibers)
    const {
  if (active_stress_ == nullptr)
    return {};

  // The fiber stretch is only computed for the models that use it. Besides
  // saving the work, this keeps the other models usable on a mesh with no fiber
  // directions, where fibers is zero and the stretch would come out zero too.
  double fiber_stretch = 1.0;

  if (active_stress_->needs_fiber_stretch())
    fiber_stretch = (F * fibers.col(0)).norm();

  // Interpolate the nodal state to the quadrature point.
  Eigen::Map<Eigen::VectorXd>(interpolated_state.data(),
                              interpolated_state.size()) = state_ * N;

  return active_stress_->compute_tension(interpolated_state, fiber_stretch);
}

template ActiveStress::ActiveTension ActiveStress::Evaluator::evaluate<2>(
    const Eigen::Ref<const Eigen::VectorXd> &, const mat_fun::Matrix<2> &,
    const Eigen::Ref<const Eigen::Matrix<double, 2, Eigen::Dynamic>> &) const;

template ActiveStress::ActiveTension ActiveStress::Evaluator::evaluate<3>(
    const Eigen::Ref<const Eigen::VectorXd> &, const mat_fun::Matrix<3> &,
    const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>> &) const;

void ActiveStress::read_parameters(const ActiveStressParameters &params) {
  eta_f = params.get_eta_f();
  eta_s = params.get_eta_s();
  eta_n = params.get_eta_n();

  implicit_state_coupling_ = params.get_implicit_state_coupling();

  calcium_scaling_factor_ = params.get_calcium_scaling_factor();
  disable_force_strain_rate_feedback_ =
      params.get_disable_force_strain_rate_feedback();

  // With the feedback disabled the stretch rate computation can be skipped.
  // Models that never use the stretch rate must not be switched on.
  needs_fiber_stretch_rate_ &= !disable_force_strain_rate_feedback_;

  read_model_specific_parameters(
      params.get_parameters(params.get_model_name()));
}

void ActiveStress::distribute_parameters(const CmMod &cm_mod,
                                         const cmType &cm) {
  cm.bcast(cm_mod, &eta_f);
  cm.bcast(cm_mod, &eta_s);
  cm.bcast(cm_mod, &eta_n);

  cm.bcast(cm_mod, &implicit_state_coupling_);

  cm.bcast(cm_mod, &calcium_scaling_factor_);
  cm.bcast(cm_mod, &disable_force_strain_rate_feedback_);
  cm.bcast(cm_mod, &needs_fiber_stretch_rate_);

  distribute_model_specific_parameters(cm_mod, cm);
}

void ActiveStress::init(const unsigned int tnNo) {
  states.resize(n_states, tnNo);

  if (n_states > 0) {
    Vector<double> state_loc(n_states);
    init_local(state_loc);

    for (unsigned int i = 0; i < tnNo; ++i)
      for (unsigned int j = 0; j < n_states; ++j)
        states(j, i) = state_loc(j);
  }

  states_at_time_step_start.resize(n_states, tnNo);
  states_at_time_step_start = states;
}

void ActiveStress::time_advance() { states_at_time_step_start = states; }

void ActiveStress::update(const double t, const double dt,
                          const Vector<double> &calcium,
                          const Vector<double> &fiber_stretch,
                          const Vector<double> &fiber_stretch_rate) {
  time = t;

  // Advance the state from the beginning of the time step.
  for (int i = 0; i < states.ncols(); ++i) {
    Vector<double> state_loc = states_at_time_step_start.col(i);
    advance_time_step_local(
        t, dt, calcium_scaling_factor_ * calcium[i], fiber_stretch[i],
        disable_force_strain_rate_feedback_ ? 0.0 : fiber_stretch_rate[i],
        state_loc);
    states.set_col(i, state_loc);
  }
}
