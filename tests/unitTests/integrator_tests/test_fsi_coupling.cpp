// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @brief Integration tests for fsi_coupling namespace functions.
 *
 * Tests the FSI interface data exchange functions: extract_fluid_traction,
 * extract_solid_displacement, apply_traction_on_solid, apply_displacement_on_mesh.
 *
 * Requires MPI and access to the FSI pipe_3d test case data files.
 */

#include "gtest/gtest.h"

#include "fsi_coupling.h"
#include "post.h"
#include "Integrator.h"
#include "Simulation.h"
#include "distribute.h"
#include "initialize.h"
#include "read_files.h"
#include "LinearAlgebra.h"
#include "set_bc.h"
#include "post.h"

#include <cmath>
#include <sys/stat.h>
#include <unistd.h>
#include <mpi.h>

#ifndef TEST_DATA_DIR
#define TEST_DATA_DIR ""
#endif

// ---------------------------------------------------------------------------
// MPI environment (same as in test_step_equation.cpp -- only one takes effect)
// ---------------------------------------------------------------------------
class MPIEnvironment_FSICoupling : public ::testing::Environment {
public:
  void SetUp() override {
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized) {
      int argc = 0;
      char** argv = nullptr;
      MPI_Init(&argc, &argv);
    }
  }
  void TearDown() override {
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (!finalized) {
      MPI_Finalize();
    }
  }
};

static testing::Environment* const mpi_env_fsi =
    testing::AddGlobalTestEnvironment(new MPIEnvironment_FSICoupling);

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
static void add_eq_la(ComMod& com_mod, eqType& lEq)
{
  lEq.linear_algebra = LinearAlgebraFactory::create_interface(lEq.linear_algebra_type);
  lEq.linear_algebra->set_preconditioner(lEq.linear_algebra_preconditioner);
  lEq.linear_algebra->initialize(com_mod, lEq);
  if (lEq.linear_algebra_assembly_type != consts::LinearAlgebraType::none) {
    lEq.linear_algebra->set_assembly(lEq.linear_algebra_assembly_type);
  }
}

static Simulation* setup_fsi_simulation()
{
  std::string fsi_dir = std::string(TEST_DATA_DIR) + "/fsi/pipe_3d";
  char orig_dir[4096];
  getcwd(orig_dir, sizeof(orig_dir));
  chdir(fsi_dir.c_str());

  auto sim = new Simulation();
  read_files_ns::read_files(sim, "solver.xml");
  distribute(sim);
  Vector<double> init_time(3);
  initialize(sim, init_time);
  for (int iEq = 0; iEq < sim->com_mod.nEq; iEq++) {
    add_eq_la(sim->com_mod, sim->com_mod.eq[iEq]);
  }

  chdir(orig_dir);
  return sim;
}

static void run_one_fsi_timestep(Simulation* sim)
{
  auto& com_mod = sim->com_mod;
  auto& integrator = sim->get_integrator();
  auto& solutions = integrator.get_solutions();

  com_mod.cTS += 1;
  com_mod.time += com_mod.dt;
  com_mod.cEq = 0;
  for (auto& eq : com_mod.eq) {
    eq.itr = 0;
    eq.ok = false;
  }

  integrator.predictor();
  set_bc::set_bc_dir(com_mod, solutions);
  integrator.step();

  solutions.old.get_acceleration() = solutions.current.get_acceleration();
  solutions.old.get_velocity() = solutions.current.get_velocity();
  if (com_mod.dFlag) {
    solutions.old.get_displacement() = solutions.current.get_displacement();
  }
}

static void teardown_sim(Simulation* sim)
{
  for (int iEq = 0; iEq < sim->com_mod.nEq; iEq++) {
    sim->com_mod.eq[iEq].linear_algebra->finalize();
  }
  delete sim;
}

static bool test_data_available()
{
  std::string path = std::string(TEST_DATA_DIR);
  if (path.empty()) return false;
  struct stat st;
  return (stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode));
}

// Find a face by name in the simulation meshes
static const faceType* find_face(const ComMod& com_mod, const std::string& name)
{
  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    for (int iFa = 0; iFa < com_mod.msh[iM].nFa; iFa++) {
      if (com_mod.msh[iM].fa[iFa].name == name) {
        return &com_mod.msh[iM].fa[iFa];
      }
    }
  }
  return nullptr;
}

static const mshType* find_mesh_for_face(const ComMod& com_mod, const std::string& face_name)
{
  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    for (int iFa = 0; iFa < com_mod.msh[iM].nFa; iFa++) {
      if (com_mod.msh[iM].fa[iFa].name == face_name) {
        return &com_mod.msh[iM];
      }
    }
  }
  return nullptr;
}

// ===========================================================================
// Tests
// ===========================================================================


/// @brief Extract solid displacement from a converged FSI solution.
TEST(FSICoupling, ExtractSolidDisplacement)
{
  if (!test_data_available()) GTEST_SKIP() << "Test data not available";

  auto sim = setup_fsi_simulation();
  auto& com_mod = sim->com_mod;
  const int nsd = com_mod.nsd;

  // Run one time step to get a non-trivial solution
  run_one_fsi_timestep(sim);

  auto& solutions = sim->get_integrator().get_solutions();
  auto& eq = com_mod.eq[0];  // FSI equation

  auto* solid_face = find_face(com_mod, "wall_inner");
  ASSERT_NE(solid_face, nullptr);

  // Extract displacement
  auto disp = fsi_coupling::extract_solid_displacement(com_mod, eq, *solid_face, solutions);

  // Verify against direct solution array access
  const auto& Dn = solutions.current.get_displacement();
  int s = eq.s;
  double max_diff = 0.0;
  for (int a = 0; a < solid_face->nNo; a++) {
    int Ac = solid_face->gN(a);
    for (int i = 0; i < nsd; i++) {
      double diff = std::abs(disp(i, a) - Dn(i + s, Ac));
      if (diff > max_diff) max_diff = diff;
    }
  }
  EXPECT_LT(max_diff, 1e-14)
      << "Extracted displacement should match solution array";

  // Verify displacement is non-zero (problem has deformation)
  double max_val = 0.0;
  for (int a = 0; a < solid_face->nNo; a++) {
    for (int i = 0; i < nsd; i++) {
      double v = std::abs(disp(i, a));
      if (v > max_val) max_val = v;
    }
  }
  EXPECT_GT(max_val, 1e-10) << "Displacement should be non-zero after FSI solve";

  teardown_sim(sim);
}

/// @brief Extract fluid traction and verify total force is reasonable.
TEST(FSICoupling, ExtractFluidTraction)
{
  if (!test_data_available()) GTEST_SKIP() << "Test data not available";

  auto sim = setup_fsi_simulation();
  auto& com_mod = sim->com_mod;
  const int nsd = com_mod.nsd;

  // Run one time step
  run_one_fsi_timestep(sim);

  auto& integrator = sim->get_integrator();
  auto& solutions = integrator.get_solutions();
  auto& eq = com_mod.eq[0];  // FSI equation

  auto* fluid_face = find_face(com_mod, "lumen_wall");
  auto* fluid_mesh = find_mesh_for_face(com_mod, "lumen_wall");
  ASSERT_NE(fluid_face, nullptr);
  ASSERT_NE(fluid_mesh, nullptr);

  // Extract consistent nodal traction forces
  com_mod.cEq = 0;  // ensure correct equation is active
  auto traction = post::compute_face_traction(
      com_mod, sim->cm_mod, *fluid_mesh, *fluid_face, eq, solutions);

  // Check dimensions
  EXPECT_EQ(traction.nrows(), nsd);
  EXPECT_EQ(traction.ncols(), fluid_face->nNo);

  // Compute total force (sum of consistent nodal forces)
  Vector<double> total_force(nsd);
  for (int a = 0; a < fluid_face->nNo; a++) {
    for (int i = 0; i < nsd; i++) {
      total_force(i) += traction(i, a);
    }
  }

  // The total force should be non-zero (pressure-driven flow in a pipe)
  double force_mag = 0.0;
  for (int i = 0; i < nsd; i++) {
    force_mag += total_force(i) * total_force(i);
  }
  force_mag = sqrt(force_mag);
  EXPECT_GT(force_mag, 1e-10)
      << "Total traction force should be non-zero for pressure-driven FSI flow";

  teardown_sim(sim);
}

