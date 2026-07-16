// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef FSI_COUPLING_H
#define FSI_COUPLING_H

#include "ComMod.h"
#include "SolutionStates.h"

/// @brief FSI interface data exchange functions for partitioned coupling.
///
/// These functions apply residual-recovered fluid nodal forces and transfer
/// solid displacement/velocity data at the FSI interface, enabling partitioned
/// (Dirichlet-Neumann) coupling between separately solved fluid and solid
/// equations.
///
/// Related to GitHub issue #431: Implement partitioned FSI in svMultiPhysics

namespace fsi_coupling {

/// @brief Extract solid displacement at interface face nodes.
///
/// @param com_mod Common module
/// @param solid_eq The solid equation (for DOF offset)
/// @param solid_face The solid-side FSI interface face
/// @param solutions Solution states
/// @return Array(nsd, solid_face.nNo) of displacement values
Array<double> extract_solid_displacement(
    const ComMod& com_mod, const eqType& solid_eq,
    const faceType& solid_face, const SolutionStates& solutions);

/// @brief Copy generalized-alpha/Newmark coefficients from one equation to another.
///
/// Partitioned FSI uses a standalone structural subproblem, but to match the
/// monolithic FSI structural domain it must use the same time-integration
/// parameters as the FSI/fluid equation.
void copy_time_integration_parameters(const eqType& source_eq,
                                      eqType& target_eq);

/// @brief Apply velocity as strong Dirichlet BC on fluid interface nodes.
/// Directly sets Yn at the fluid equation DOF range for the face nodes.
void apply_velocity_on_fluid(
    ComMod& com_mod, const eqType& fluid_eq,
    const faceType& fluid_face,
    const Array<double>& velocity,
    SolutionStates& solutions);

/// @brief Extract force ON the solid from the unconstrained fluid residual.
///
/// The returned array contains the velocity-row residual at the fluid interface
/// converted to the external-force-on-solid convention expected by
/// apply_traction_on_solid().
///
/// @param com_mod Common module for the fluid subproblem
/// @param fluid_face The fluid-side FSI interface face
/// @param residual Unconstrained assembled fluid residual
/// @return Array(nsd, fluid_face.nNo) of consistent nodal forces
Array<double> extract_fluid_residual_force(
    const ComMod& com_mod,
    const faceType& fluid_face,
    const Array<double>& residual);

/// @brief Apply pre-computed consistent nodal forces to the solid residual.
///
/// The force array is interpreted as the external force ON the solid. It is
/// subtracted from com_mod.R at the global node locations corresponding to the
/// solid face. This should be called during the post-assembly callback of
/// step_equation() for the solid equation.
///
/// @param com_mod Common module (R is modified)
/// @param solid_eq The solid equation (for DOF offset)
/// @param solid_face The solid-side FSI interface face
/// @param force Array(nsd, solid_face.nNo) of consistent nodal forces
void apply_traction_on_solid(
    ComMod& com_mod, const eqType& solid_eq,
    const faceType& solid_face,
    const Array<double>& force);

/// @brief Apply displacement as strong Dirichlet BC on mesh interface nodes.
///
/// Directly sets the displacement in the solution arrays (An, Yn, Dn) for
/// the mesh equation DOF range at the interface face nodes.
///
/// @param com_mod Common module
/// @param mesh_eq The mesh equation (for DOF offset)
/// @param mesh_face The mesh-side FSI interface face
/// @param displacement Array(nsd, mesh_face.nNo) of displacement values
/// @param solutions Solution states (modified)
void apply_displacement_on_mesh(
    ComMod& com_mod, const eqType& mesh_eq,
    const faceType& mesh_face,
    const Array<double>& displacement,
    SolutionStates& solutions);

/// @brief Compute staged ALE fluid mesh coordinates from mesh displacement.
///
/// The returned coordinates are
/// x_stage = x_ref + theta * (mesh_Dn - mesh_Do), using the mesh equation
/// displacement block starting at mesh_s. For a generalized-alpha fluid solve,
/// theta should be the fluid equation alpha_f; before output/old-state update,
/// theta should be 1.0.
Array<double> staged_fluid_mesh_coordinates(
    const Array<double>& x_ref,
    const Array<double>& mesh_Dn,
    const Array<double>& mesh_Do,
    int mesh_s,
    int nsd,
    double theta);

} // namespace fsi_coupling

#endif // FSI_COUPLING_H
