// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @brief Integration tests for Integrator::step_equation()
 *
 * Tests that step_equation() produces correct results by comparing against
 * the monolithic step() method on real solver problems.
 *
 * These tests require MPI and access to test case data files.
 * They are skipped if the test data directory is not available.
 */

#include "gtest/gtest.h"

#include "Integrator.h"
#include "Simulation.h"
#include "distribute.h"
#include "initialize.h"
#include "read_files.h"
#include "LinearAlgebra.h"
#include "set_bc.h"

#include <cmath>
#include <sys/stat.h>
#include <unistd.h>
#include <mpi.h>

// Path to test data, defined via CMake
#ifndef TEST_DATA_DIR
#define TEST_DATA_DIR ""
#endif

// ---------------------------------------------------------------------------
// MPI environment: initializes MPI once before all tests, finalizes after
// ---------------------------------------------------------------------------
class MPIEnvironment : public ::testing::Environment {
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

// Register the MPI environment (runs before any test)
static testing::Environment* const mpi_env =
    testing::AddGlobalTestEnvironment(new MPIEnvironment);

// ---------------------------------------------------------------------------
// Helper: set up a full simulation from an XML file
// ---------------------------------------------------------------------------
static Simulation* setup_simulation(const std::string& xml_path)
{
  // The solver reads mesh files relative to the XML directory,
  // so we must chdir there before calling read_files.
  std::string dir = xml_path.substr(0, xml_path.find_last_of('/'));
  std::string file = xml_path.substr(xml_path.find_last_of('/') + 1);
  char orig_dir[4096];
  getcwd(orig_dir, sizeof(orig_dir));
  chdir(dir.c_str());

  auto simulation = new Simulation();
  read_files_ns::read_files(simulation, file);
  distribute(simulation);
  Vector<double> init_time(3);
  initialize(simulation, init_time);
  for (int iEq = 0; iEq < simulation->com_mod.nEq; iEq++) {
    add_eq_linear_algebra(simulation->com_mod, simulation->com_mod.eq[iEq]);
  }

  chdir(orig_dir);
  return simulation;
}

static void teardown_simulation(Simulation* simulation)
{
  for (int iEq = 0; iEq < simulation->com_mod.nEq; iEq++) {
    simulation->com_mod.eq[iEq].linear_algebra->finalize();
  }
  delete simulation;
}

// ---------------------------------------------------------------------------
// Helper: run one time step using step()
// ---------------------------------------------------------------------------
static void run_one_timestep_step(Simulation* simulation)
{
  auto& com_mod = simulation->com_mod;
  auto& integrator = simulation->get_integrator();
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

  // Copy current -> old for next step
  solutions.old.get_acceleration() = solutions.current.get_acceleration();
  solutions.old.get_velocity() = solutions.current.get_velocity();
  if (com_mod.dFlag) {
    solutions.old.get_displacement() = solutions.current.get_displacement();
  }
  com_mod.cplBC.xo = com_mod.cplBC.xn;
}

// ---------------------------------------------------------------------------
// Helper: run one time step using step_equation() per equation
// ---------------------------------------------------------------------------
static void run_one_timestep_step_equation(Simulation* simulation,
                                           int outer_iters = 1)
{
  auto& com_mod = simulation->com_mod;
  auto& integrator = simulation->get_integrator();
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

  for (int outer = 0; outer < outer_iters; outer++) {
    for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
      integrator.step_equation(iEq);
    }
  }

  // Copy current -> old for next step
  solutions.old.get_acceleration() = solutions.current.get_acceleration();
  solutions.old.get_velocity() = solutions.current.get_velocity();
  if (com_mod.dFlag) {
    solutions.old.get_displacement() = solutions.current.get_displacement();
  }
  com_mod.cplBC.xo = com_mod.cplBC.xn;
}

// ---------------------------------------------------------------------------
// Helper: compute relative difference between two solution arrays
// ---------------------------------------------------------------------------
static double rel_diff(const Array<double>& a, const Array<double>& b)
{
  double max_diff = 0.0;
  double max_val = 0.0;
  for (int j = 0; j < a.ncols(); j++) {
    for (int i = 0; i < a.nrows(); i++) {
      double d = std::abs(a(i,j) - b(i,j));
      if (d > max_diff) max_diff = d;
      double v = std::abs(a(i,j));
      if (v > max_val) max_val = v;
    }
  }
  return (max_val > 0) ? max_diff / max_val : max_diff;
}

// ---------------------------------------------------------------------------
// Helper: check if test data directory exists
// ---------------------------------------------------------------------------
static bool test_data_available()
{
  std::string path = std::string(TEST_DATA_DIR);
  if (path.empty()) return false;
  struct stat st;
  return (stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode));
}

// ===========================================================================
// Tests
// ===========================================================================

/// @brief For a single-equation problem, step_equation(0) must produce
/// bit-identical results to step().
TEST(StepEquation, SingleEquationMatchesStep)
{
  if (!test_data_available()) GTEST_SKIP() << "Test data not available";

  std::string xml = std::string(TEST_DATA_DIR) + "/fluid/newtonian/solver.xml";

  // Run with step()
  auto sim_a = setup_simulation(xml);
  run_one_timestep_step(sim_a);
  Array<double> Yn_step = sim_a->get_integrator().get_solutions().current.get_velocity();
  teardown_simulation(sim_a);

  // Run with step_equation(0)
  auto sim_b = setup_simulation(xml);
  run_one_timestep_step_equation(sim_b, 1);
  Array<double> Yn_step_eq = sim_b->get_integrator().get_solutions().current.get_velocity();
  teardown_simulation(sim_b);

  double diff = rel_diff(Yn_step, Yn_step_eq);
  EXPECT_LT(diff, 1e-12) << "step_equation(0) should match step() for single-equation problems";
}

/// @brief For coupled FSI, sequential step_equation converges to step()
/// result when outer coupling iterations are added.
TEST(StepEquation, CoupledFSIConvergesWithOuterIterations)
{
  if (!test_data_available()) GTEST_SKIP() << "Test data not available";

  std::string xml = std::string(TEST_DATA_DIR) + "/fsi/pipe_3d/solver.xml";

  // Run with step() as reference
  auto sim_ref = setup_simulation(xml);
  run_one_timestep_step(sim_ref);
  Array<double> Yn_ref = sim_ref->get_integrator().get_solutions().current.get_velocity();
  Array<double> Dn_ref = sim_ref->get_integrator().get_solutions().current.get_displacement();
  teardown_simulation(sim_ref);

  // Run with 1 outer iteration (single pass) - should have coupling error
  auto sim_1 = setup_simulation(xml);
  run_one_timestep_step_equation(sim_1, 1);
  Array<double> Yn_1 = sim_1->get_integrator().get_solutions().current.get_velocity();
  Array<double> Dn_1 = sim_1->get_integrator().get_solutions().current.get_displacement();
  teardown_simulation(sim_1);

  double diff_vel_1 = rel_diff(Yn_ref, Yn_1);
  double diff_disp_1 = rel_diff(Dn_ref, Dn_1);

  // Run with 4 outer iterations - should converge to step() result
  auto sim_4 = setup_simulation(xml);
  run_one_timestep_step_equation(sim_4, 4);
  Array<double> Yn_4 = sim_4->get_integrator().get_solutions().current.get_velocity();
  Array<double> Dn_4 = sim_4->get_integrator().get_solutions().current.get_displacement();
  teardown_simulation(sim_4);

  double diff_vel_4 = rel_diff(Yn_ref, Yn_4);
  double diff_disp_4 = rel_diff(Dn_ref, Dn_4);

  // With 4 outer iterations, the coupling error should be orders of magnitude
  // smaller than with 1 iteration
  EXPECT_LT(diff_vel_4, diff_vel_1 * 1e-4)
      << "4 outer iterations should reduce velocity coupling error by >4 orders";
  EXPECT_LT(diff_disp_4, diff_disp_1 * 1e-4)
      << "4 outer iterations should reduce displacement coupling error by >4 orders";

  // With 4 outer iterations, the result should match step() to near machine precision
  EXPECT_LT(diff_vel_4, 1e-10)
      << "Velocity should match step() after 4 outer iterations";
  EXPECT_LT(diff_disp_4, 1e-10)
      << "Displacement should match step() after 4 outer iterations";
}

/// @brief The post_assembly callback fires on every Newton iteration.
TEST(StepEquation, PostAssemblyCallbackFires)
{
  if (!test_data_available()) GTEST_SKIP() << "Test data not available";

  std::string xml = std::string(TEST_DATA_DIR) + "/fluid/newtonian/solver.xml";

  auto sim = setup_simulation(xml);
  auto& com_mod = sim->com_mod;
  auto& integrator = sim->get_integrator();
  auto& solutions = integrator.get_solutions();

  // Set up time step
  com_mod.cTS += 1;
  com_mod.time += com_mod.dt;
  com_mod.cEq = 0;
  for (auto& eq : com_mod.eq) {
    eq.itr = 0;
    eq.ok = false;
  }

  integrator.predictor();
  set_bc::set_bc_dir(com_mod, solutions);

  // Count callback invocations
  int callback_count = 0;
  integrator.step_equation(0, [&callback_count]() {
    callback_count++;
  });

  // Callback should fire once per Newton iteration (at least minItr times)
  EXPECT_GE(callback_count, com_mod.eq[0].minItr)
      << "Callback should fire at least minItr times";
  EXPECT_LE(callback_count, com_mod.eq[0].maxItr)
      << "Callback should fire at most maxItr times";

  teardown_simulation(sim);
}
