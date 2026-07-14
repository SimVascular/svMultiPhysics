// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @brief Unit tests for FSI interface data exchange helpers.
 *
 * These tests avoid full mesh I/O so they can check the coupling data
 * conventions directly: solid displacement extraction, residual-based fluid
 * force extraction semantics, and solid residual force application.
 */

#include "gtest/gtest.h"

#include "fsi_coupling.h"

#include <cmath>
#include <vector>

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
static faceType make_face(const std::vector<int>& nodes)
{
  faceType face;
  face.nNo = static_cast<int>(nodes.size());
  face.gN.resize(face.nNo);
  for (int a = 0; a < face.nNo; a++) {
    face.gN(a) = nodes[a];
  }
  return face;
}

// ===========================================================================
// Tests
// ===========================================================================

/// @brief Extract solid displacement at face nodes using the solid equation offset.
TEST(FSICoupling, ExtractSolidDisplacement)
{
  ComMod com_mod;
  com_mod.nsd = 3;

  eqType solid_eq;
  solid_eq.s = 2;

  auto solid_face = make_face({0, 2, 4});

  SolutionStates solutions;
  auto& Dn = solutions.current.get_displacement();
  Dn.resize(5, 5);
  Dn = 0.0;
  for (int a = 0; a < solid_face.nNo; a++) {
    int Ac = solid_face.gN(a);
    for (int i = 0; i < com_mod.nsd; i++) {
      Dn(i + solid_eq.s, Ac) = 10.0 * (i + 1) + Ac;
    }
  }

  auto disp = fsi_coupling::extract_solid_displacement(
      com_mod, solid_eq, solid_face, solutions);

  EXPECT_EQ(disp.nrows(), com_mod.nsd);
  EXPECT_EQ(disp.ncols(), solid_face.nNo);
  for (int a = 0; a < solid_face.nNo; a++) {
    int Ac = solid_face.gN(a);
    for (int i = 0; i < com_mod.nsd; i++) {
      EXPECT_DOUBLE_EQ(disp(i, a), Dn(i + solid_eq.s, Ac));
    }
  }
}

/// @brief Partitioned solid time integration should match the FSI/fluid equation.
TEST(FSICoupling, CopyTimeIntegrationParameters)
{
  eqType source_eq;
  source_eq.roInf = 0.5;
  source_eq.am = 5.0 / 6.0;
  source_eq.af = 2.0 / 3.0;
  source_eq.beta = 49.0 / 144.0;
  source_eq.gam = 2.0 / 3.0;
  source_eq.dof = 4;
  source_eq.sym = "FS";

  eqType target_eq;
  target_eq.roInf = 0.5;
  target_eq.am = 1.0;
  target_eq.af = 2.0 / 3.0;
  target_eq.beta = 4.0 / 9.0;
  target_eq.gam = 5.0 / 6.0;
  target_eq.dof = 3;
  target_eq.sym = "ST";

  fsi_coupling::copy_time_integration_parameters(source_eq, target_eq);

  EXPECT_DOUBLE_EQ(target_eq.roInf, source_eq.roInf);
  EXPECT_DOUBLE_EQ(target_eq.am, source_eq.am);
  EXPECT_DOUBLE_EQ(target_eq.af, source_eq.af);
  EXPECT_DOUBLE_EQ(target_eq.beta, source_eq.beta);
  EXPECT_DOUBLE_EQ(target_eq.gam, source_eq.gam);
  EXPECT_EQ(target_eq.dof, 3);
  EXPECT_EQ(target_eq.sym, "ST");
}

/// @brief The residual-recovered force is the opposite of the fluid residual.
TEST(FSICoupling, ExtractFluidResidualForceUsesUnconstrainedRows)
{
  ComMod com_mod;
  com_mod.nsd = 3;
  com_mod.tnNo = 5;

  auto fluid_face = make_face({1, 3});

  Array<double> residual(4, com_mod.tnNo);
  for (int a = 0; a < residual.ncols(); a++) {
    for (int i = 0; i < residual.nrows(); i++) {
      residual(i, a) = 100.0 * (i + 1) + a;
    }
  }

  auto force = fsi_coupling::extract_fluid_residual_force(
      com_mod, fluid_face, residual);

  EXPECT_EQ(force.nrows(), com_mod.nsd);
  EXPECT_EQ(force.ncols(), fluid_face.nNo);
  for (int a = 0; a < fluid_face.nNo; a++) {
    int Ac = fluid_face.gN(a);
    for (int i = 0; i < com_mod.nsd; i++) {
      EXPECT_DOUBLE_EQ(force(i, a), -residual(i, Ac));
    }
  }

  // Simulate strong Dirichlet row clearing after extraction. The recovered force
  // must retain the unconstrained residual values.
  for (int a = 0; a < fluid_face.nNo; a++) {
    int Ac = fluid_face.gN(a);
    for (int i = 0; i < com_mod.nsd; i++) {
      residual(i, Ac) = 0.0;
    }
  }
  for (int a = 0; a < fluid_face.nNo; a++) {
    int Ac = fluid_face.gN(a);
    for (int i = 0; i < com_mod.nsd; i++) {
      EXPECT_NE(force(i, a), residual(i, Ac));
    }
  }
}

/// @brief Residual recovery gives the expected total force for a pressure load.
TEST(FSICoupling, RecoveredFluidForceMatchesKnownPressureLoad)
{
  ComMod com_mod;
  com_mod.nsd = 3;
  com_mod.tnNo = 4;

  auto fluid_face = make_face({0, 1, 2, 3});

  const double pressure = 2.5;
  const double area = 4.0;
  const double nodal_force = pressure * area / fluid_face.nNo;

  Array<double> residual(4, com_mod.tnNo);
  residual = 0.0;
  for (int a = 0; a < fluid_face.nNo; a++) {
    int Ac = fluid_face.gN(a);
    residual(0, Ac) = -nodal_force;
    residual(3, Ac) = 999.0;  // pressure row is not part of the recovered force
  }

  auto force = fsi_coupling::extract_fluid_residual_force(
      com_mod, fluid_face, residual);

  double total_force[3] = {0.0, 0.0, 0.0};
  for (int a = 0; a < fluid_face.nNo; a++) {
    for (int i = 0; i < com_mod.nsd; i++) {
      total_force[i] += force(i, a);
    }
  }

  EXPECT_DOUBLE_EQ(total_force[0], pressure * area);
  EXPECT_DOUBLE_EQ(total_force[1], 0.0);
  EXPECT_DOUBLE_EQ(total_force[2], 0.0);
}

/// @brief Applying an external force to the solid subtracts it from R.
TEST(FSICoupling, ApplyForceOnSolidSubtractsExternalForce)
{
  ComMod com_mod;
  com_mod.nsd = 3;
  com_mod.dof = 4;
  com_mod.tnNo = 5;
  com_mod.R.resize(com_mod.dof, com_mod.tnNo);
  com_mod.R = 0.0;

  eqType solid_eq;
  auto solid_face = make_face({0, 2, 4});

  Array<double> force(com_mod.nsd, solid_face.nNo);
  for (int a = 0; a < solid_face.nNo; a++) {
    for (int i = 0; i < com_mod.nsd; i++) {
      force(i, a) = static_cast<double>((i + 1) * (a + 1));
    }
  }

  fsi_coupling::apply_traction_on_solid(
      com_mod, solid_eq, solid_face, force);

  for (int a = 0; a < solid_face.nNo; a++) {
    int Ac = solid_face.gN(a);
    for (int i = 0; i < com_mod.nsd; i++) {
      EXPECT_DOUBLE_EQ(com_mod.R(i, Ac), -force(i, a));
    }
  }
}

/// @brief Fluid ALE coordinates are staged from the mesh displacement increment.
TEST(FSICoupling, StagedFluidMeshCoordinatesUseThetaIncrement)
{
  const int nsd = 3;
  const int mesh_s = 2;
  const double theta = 2.0 / 3.0;

  Array<double> x_ref(nsd, 2);
  Array<double> mesh_Do(mesh_s + nsd, 2);
  Array<double> mesh_Dn(mesh_s + nsd, 2);

  for (int a = 0; a < x_ref.ncols(); a++) {
    for (int i = 0; i < nsd; i++) {
      x_ref(i, a) = 10.0 * (i + 1) + a;
    }
    for (int i = 0; i < mesh_Do.nrows(); i++) {
      mesh_Do(i, a) = 0.25 * (i + 1) + a;
      mesh_Dn(i, a) = mesh_Do(i, a) + 0.1 * (i + 1) * (a + 1);
    }
  }

  auto x_stage = fsi_coupling::staged_fluid_mesh_coordinates(
      x_ref, mesh_Dn, mesh_Do, mesh_s, nsd, theta);

  EXPECT_EQ(x_stage.nrows(), nsd);
  EXPECT_EQ(x_stage.ncols(), x_ref.ncols());
  for (int a = 0; a < x_ref.ncols(); a++) {
    for (int i = 0; i < nsd; i++) {
      const double increment = mesh_Dn(i + mesh_s, a)
                             - mesh_Do(i + mesh_s, a);
      EXPECT_DOUBLE_EQ(x_stage(i, a), x_ref(i, a) + theta * increment);
    }
  }
}
