// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "fsi_coupling.h"
#include "Integrator.h"

#include <cmath>
#include <stdexcept>

namespace fsi_coupling {

//----------------------------------------------------------------------
// extract_solid_displacement
//----------------------------------------------------------------------
Array<double> extract_solid_displacement(
    const ComMod& com_mod, const eqType& solid_eq,
    const faceType& lFa, const SolutionStates& solutions)
{
  const int nsd = com_mod.nsd;
  const int s = solid_eq.s;  // DOF offset for the solid equation
  const auto& Dn = solutions.current.get_displacement();

  Array<double> result(nsd, lFa.nNo);
  for (int a = 0; a < lFa.nNo; a++) {
    int Ac = lFa.gN(a);
    for (int i = 0; i < nsd; i++) {
      result(i, a) = Dn(i + s, Ac);
    }
  }
  return result;
}

//----------------------------------------------------------------------
// copy_time_integration_parameters
//----------------------------------------------------------------------
void copy_time_integration_parameters(const eqType& source_eq,
                                      eqType& target_eq)
{
  target_eq.roInf = source_eq.roInf;
  target_eq.am = source_eq.am;
  target_eq.af = source_eq.af;
  target_eq.beta = source_eq.beta;
  target_eq.gam = source_eq.gam;
}

//----------------------------------------------------------------------
// apply_velocity_on_fluid
//----------------------------------------------------------------------
void apply_velocity_on_fluid(
    ComMod& com_mod, const eqType& fluid_eq,
    const faceType& lFa,
    const Array<double>& velocity,
    SolutionStates& solutions)
{
  const int nsd = com_mod.nsd;
  const int s = fluid_eq.s;
  const double dt = com_mod.dt;
  const double gam = fluid_eq.gam;

  auto& An = solutions.current.get_acceleration();
  auto& Yn = solutions.current.get_velocity();
  const auto& Yo = solutions.old.get_velocity();
  const auto& Ao = solutions.old.get_acceleration();

  for (int a = 0; a < lFa.nNo; a++) {
    int Ac = lFa.gN(a);
    for (int i = 0; i < nsd; i++) {
      Yn(i + s, Ac) = velocity(i, a);
      double a_new;
      newmark::state_from_velocity(
          velocity(i, a), Yo(i + s, Ac), Ao(i + s, Ac), dt, gam, a_new);
      An(i + s, Ac) = a_new;
    }
  }
}

//----------------------------------------------------------------------
// extract_fluid_residual_force
//----------------------------------------------------------------------
Array<double> extract_fluid_residual_force(
    const ComMod& com_mod,
    const faceType& fluid_face,
    const Array<double>& residual)
{
  const int nsd = com_mod.nsd;
  if (residual.nrows() < nsd || residual.ncols() != com_mod.tnNo) {
    throw std::runtime_error(
        "[fsi_coupling] Fluid residual shape does not match fluid simulation.");
  }

  Array<double> force(nsd, fluid_face.nNo);
  for (int a = 0; a < fluid_face.nNo; a++) {
    int Ac = fluid_face.gN(a);
    for (int i = 0; i < nsd; i++) {
      force(i, a) = -residual(i, Ac);
    }
  }
  return force;
}

//----------------------------------------------------------------------
// apply_traction_on_solid
//----------------------------------------------------------------------
void apply_traction_on_solid(
    ComMod& com_mod, const eqType& solid_eq,
    const faceType& lFa,
    const Array<double>& force)
{
  // The array contains consistent nodal forces (external force on solid).
  // In svMultiPhysics, external forces are SUBTRACTED from R (for example:
  // lR -= w*N*h). So R -= force.
  for (int a = 0; a < lFa.nNo; a++) {
    int Ac = lFa.gN(a);
    for (int i = 0; i < force.nrows(); i++) {
      com_mod.R(i, Ac) -= force(i, a);
    }
  }
}

//----------------------------------------------------------------------
// apply_displacement_on_mesh
//----------------------------------------------------------------------
void apply_displacement_on_mesh(
    ComMod& com_mod, const eqType& mesh_eq,
    const faceType& lFa,
    const Array<double>& displacement,
    SolutionStates& solutions)
{
  const int nsd = com_mod.nsd;
  const int s = mesh_eq.s;
  const double dt = com_mod.dt;
  const double gam = mesh_eq.gam;
  const double beta = mesh_eq.beta;

  auto& An = solutions.current.get_acceleration();
  auto& Yn = solutions.current.get_velocity();
  auto& Dn = solutions.current.get_displacement();
  const auto& Do = solutions.old.get_displacement();
  const auto& Yo = solutions.old.get_velocity();
  const auto& Ao = solutions.old.get_acceleration();

  for (int a = 0; a < lFa.nNo; a++) {
    int Ac = lFa.gN(a);
    for (int i = 0; i < nsd; i++) {
      Dn(i + s, Ac) = displacement(i, a);
      double a_new, v_new;
      newmark::state_from_displacement(
          displacement(i, a), Do(i + s, Ac), Yo(i + s, Ac), Ao(i + s, Ac),
          dt, beta, gam, a_new, v_new);
      An(i + s, Ac) = a_new;
      Yn(i + s, Ac) = v_new;
    }
  }
}

//----------------------------------------------------------------------
// staged_fluid_mesh_coordinates
//----------------------------------------------------------------------
Array<double> staged_fluid_mesh_coordinates(
    const Array<double>& x_ref,
    const Array<double>& mesh_Dn,
    const Array<double>& mesh_Do,
    int mesh_s,
    int nsd,
    double theta)
{
  if (nsd <= 0 || mesh_s < 0 || !std::isfinite(theta)) {
    throw std::runtime_error(
        "[fsi_coupling] Invalid inputs for staged fluid mesh coordinates.");
  }
  if (x_ref.nrows() != nsd ||
      mesh_Dn.ncols() != x_ref.ncols() ||
      mesh_Do.ncols() != x_ref.ncols() ||
      mesh_Dn.nrows() < mesh_s + nsd ||
      mesh_Do.nrows() < mesh_s + nsd) {
    throw std::runtime_error(
        "[fsi_coupling] Mesh displacement shapes do not match fluid coordinates.");
  }

  Array<double> x_stage(nsd, x_ref.ncols());
  for (int a = 0; a < x_ref.ncols(); a++) {
    for (int i = 0; i < nsd; i++) {
      x_stage(i, a) = x_ref(i, a)
                    + theta * (mesh_Dn(i + mesh_s, a)
                             - mesh_Do(i + mesh_s, a));
    }
  }
  return x_stage;
}

} // namespace fsi_coupling
