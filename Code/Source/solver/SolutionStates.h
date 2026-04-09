// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SOLUTION_STATES_H
#define SOLUTION_STATES_H

#include "Array.h"

/**
 * @brief Represents solution variables at a single time step
 *
 * Contains the three primary solution arrays used in time integration:
 * A (time derivative), D (integrated variable), and Y (variable)
 */
struct Solution {
  Array<double>& get_acceleration() { return acceleration_; }
  const Array<double>& get_acceleration() const { return acceleration_; }

  Array<double>& get_velocity() { return velocity_; }
  const Array<double>& get_velocity() const { return velocity_; }

  Array<double>& get_displacement() { return displacement_; }
  const Array<double>& get_displacement() const { return displacement_; }

private:
  Array<double> acceleration_;  ///< Time derivative (acceleration in structural mechanics)
  Array<double> displacement_;  ///< Integrated variable (displacement in structural mechanics)
  Array<double> velocity_;      ///< Variable (velocity in structural mechanics)
};

/**
 * @brief Holds solution state at old, current, and intermediate time levels
 *
 * Contains solution arrays at three time levels for time integration:
 * - old: Previous converged solution at time n
 * - current: Current solution being computed at time n+1
 * - intermediate: Generalized-alpha intermediate level (Ag, Yg, Dg)
 */
struct SolutionStates {
  Solution old;           ///< Previous converged solution at time n (Ao, Do, Yo)
  Solution current;       ///< Current solution being computed at time n+1 (An, Dn, Yn)
  Solution intermediate;  ///< Generalized-alpha intermediate level (Ag, Yg, Dg)
};

#endif
