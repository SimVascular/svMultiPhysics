// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "fsi_coupling.h"
#include "Integrator.h"

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
// apply_traction_on_solid
//----------------------------------------------------------------------
void apply_traction_on_solid(
    ComMod& com_mod, const eqType& solid_eq,
    const faceType& lFa,
    const Array<double>& traction)
{
  // The traction array contains consistent nodal forces (external force on solid).
  // In svMultiPhysics, external forces are SUBTRACTED from R (see b_l_elas:
  // lR -= w*N*h). So R -= traction.
  for (int a = 0; a < lFa.nNo; a++) {
    int Ac = lFa.gN(a);
    for (int i = 0; i < traction.nrows(); i++) {
      com_mod.R(i, Ac) -= traction(i, a);
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

} // namespace fsi_coupling
