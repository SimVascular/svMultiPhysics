/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUDAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

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
