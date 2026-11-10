// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PSEUDO_TRANSIENT_CONTINUATION_H
#define PSEUDO_TRANSIENT_CONTINUATION_H

#include "CmMod.h"
#include "ComMod.h"

/// @brief Pseudo-transient continuation: solve a steady nonlinear problem by
/// time stepping with a pseudo time step that grows as the residual falls, so
/// that the iteration turns into Newton's method (Kelley & Keyes, SIAM J.
/// Numer. Anal. 35 (1998) 508-523).
///
/// Each pseudo-step is one or more Newton iterations, each composed into the
/// state it solves for (see prestrain.h); with Max_iterations 1, the default,
/// it is a single one, as in the reference. The pseudo time step of the next step is
/// set from the ratio of successive first residuals, capped by max_dt; a time
/// step so large that the inertia term vanishes makes each iteration a plain
/// Newton iteration. Time is fictitious, so the loads must be steady.
///
/// Owned and driven by the Integrator. Enabled by Pseudo_transient on the
/// run's one equation, whose residual it monitors. Prestrain is the only
/// feature that uses it so far: each of its iterations starts from rest and the
/// displacement it reaches is composed into the prestrain (see prestrain.h).
class PseudoTransientContinuation {

  public:
    PseudoTransientContinuation() = default;
    explicit PseudoTransientContinuation(const PseudoTransientSettings& settings);

    bool enabled() const { return settings_.enabled; }

    /// @brief Finish a pseudo-step: report the update it reached and set the
    /// next pseudo time step from its first residual. Called once the Newton
    /// loop has exited.
    /// @param update this rank's largest nodal update composed by the step's
    /// last iteration, which tends to zero as the step converges.
    void finish_step(ComMod& com_mod, const CmMod& cm_mod, double update);

  private:
    /// @brief The pseudo time step for the next step, from the first residuals
    /// of the last two.
    double next_time_step(const double dt, const double residual_prev, const double residual) const;

    PseudoTransientSettings settings_;

    /// @brief First Newton residual of the previous pseudo-step
    double residual_prev_ = 0.0;
};

#endif
