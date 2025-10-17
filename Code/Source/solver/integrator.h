// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "Array.h"
#include "Vector.h"
#include "Simulation.h"

/// @brief Integrator class encapsulates the Newton iteration loop for time integration.
///
/// This class handles:
/// - Solution variables (Ag, Yg, Dg)
/// - Newton iteration loop
/// - Linear system assembly and solve
/// - Boundary condition application
///
/// Related to GitHub issue #442: Encapsulate the Newton iteration in main.cpp
class Integrator {

public:
  /// @brief Constructor
  /// @param simulation Pointer to the Simulation object
  Integrator(Simulation* simulation);

  /// @brief Destructor
  ~Integrator();

  /// @brief Execute one Newton iteration loop for the current time step
  /// @return True if all equations converged, false otherwise
  bool step();

  /// @brief Get reference to solution variable Ag (time derivative of variables)
  Array<double>& get_Ag() { return Ag_; }

  /// @brief Get reference to solution variable Yg (variables)
  Array<double>& get_Yg() { return Yg_; }

  /// @brief Get reference to solution variable Dg (integrated variables)
  Array<double>& get_Dg() { return Dg_; }

private:
  /// @brief Pointer to the simulation object
  Simulation* simulation_;

  /// @brief Solution variables at generalized-alpha time levels
  /// Time derivative of variables (acceleration)
  Array<double> Ag_;

  /// Variables (velocity)
  Array<double> Yg_;

  /// Integrated variables (displacement)
  Array<double> Dg_;

  /// @brief Residual for face-based quantities
  Vector<double> res_;

  /// @brief Increment flag for faces
  Vector<int> incL_;

  /// @brief Inner iteration counter
  int inner_count_;

  /// @brief Initialize solution arrays
  void initialize_arrays();

  /// @brief Perform initiator step for Generalized alpha-Method
  void initiator_step();

  /// @brief Allocate RHS and LHS arrays
  void allocate_linear_system(eqType& eq);

  /// @brief Set body forces
  void set_body_forces();

  /// @brief Assemble global equations
  void assemble_equations();

  /// @brief Apply boundary conditions
  void apply_boundary_conditions();

  /// @brief Solve linear system
  void solve_linear_system();

  /// @brief Perform corrector step and check convergence
  /// @return True if all equations converged
  bool corrector_and_check_convergence();

  /// @brief Update residual and increment arrays
  void update_residual_arrays(eqType& eq);
};

#endif // INTEGRATOR_H
