// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "Simulation.h"
#include "Integrator.h"
#include "PartitionedFSI.h"

#include "all_fun.h"
#include "load_msh.h"

#include "mpi.h"

#include <iostream>
#include <stdexcept>

void add_eq_linear_algebra(ComMod& com_mod, eqType& lEq)
{
  lEq.linear_algebra = LinearAlgebraFactory::create_interface(lEq.linear_algebra_type);
  lEq.linear_algebra->set_preconditioner(lEq.linear_algebra_preconditioner);
  lEq.linear_algebra->initialize(com_mod, lEq);

  if (lEq.linear_algebra_assembly_type != consts::LinearAlgebraType::none) {
    lEq.linear_algebra->set_assembly(lEq.linear_algebra_assembly_type);
  }
}

Simulation::Simulation() 
{
  roInf = 0.2;
  com_mod.cm.new_cm(MPI_COMM_WORLD);

  history_file_name = "histor.dat";
}

Simulation::~Simulation() 
{
}

const mshType& Simulation::get_msh(const std::string& name)
{
  for (auto& mesh : com_mod.msh) { 
    if (mesh.name == name) {
      return mesh;
    }
  }
}

/// @brief Read solver parameters.
//
void Simulation::read_parameters(const std::string& file_name)
{
  parameters.read_xml(file_name);
}

/// @brief Set the simulation and module member data.
///
/// Replicates the README subroutine lines to set COMMOD module varliables
///
///   lPtr => list%get(nTs,"Number of time steps",1,ll=1)
//
void Simulation::set_module_parameters()
{
  // Set ComMod module varliables.
  //
  auto& general = parameters.general_simulation_parameters;

  com_mod.iniFilePath = general.simulation_initialization_file_path.value();
  com_mod.nsd = general.number_of_spatial_dimensions.value();
  com_mod.nsymd = 3*(com_mod.nsd-1);

  com_mod.nTS = general.number_of_time_steps.value();
  com_mod.nITs = general.number_of_initialization_time_steps.value();
  com_mod.startTS = general.starting_time_step.value();
  com_mod.dt = general.time_step_size.value();

  com_mod.stopTrigName = general.searched_file_name_to_trigger_stop.value();
  com_mod.ichckIEN = general.check_ien_order.value();
  com_mod.saveVTK = general.save_results_to_vtk_format.value();
  com_mod.saveName = general.name_prefix_of_saved_vtk_files.value();
  com_mod.saveName = chnl_mod.appPath + com_mod.saveName;
  com_mod.saveIncr = general.increment_in_saving_vtk_files.value();
  com_mod.saveATS = general.start_saving_after_time_step.value();
  com_mod.saveAve = general.save_averaged_results.value();
  com_mod.zeroAve = general.start_averaging_from_zero.value();
  com_mod.stFileRepl = general.overwrite_restart_file.value();
  com_mod.stFileName = chnl_mod.appPath + general.restart_file_name.value();
  com_mod.stFileIncr = general.increment_in_saving_restart_files.value();
  com_mod.rmsh.isReqd = general.simulation_requires_remeshing.value();

  auto& precomp_sol = parameters.precomputed_solution_parameters;
  com_mod.usePrecomp = precomp_sol.use_precomputed_solution.value();
  com_mod.precompFileName = precomp_sol.file_path.value();
  com_mod.precompFieldName = precomp_sol.field_name.value();
  com_mod.precompDt = precomp_sol.time_step.value();

  if ((com_mod.precompDt == 0.0) && (com_mod.usePrecomp)) {
    std::cout << "Precomputed time step size is zero. Setting to simulation time step size." << std::endl;
    com_mod.precompDt = com_mod.dt;
  }
  // Set simulation parameters.
  nTs = general.number_of_time_steps.value();
  fTmp = general.simulation_initialization_file_path.value();
  roInf = general.spectral_radius_of_infinite_time_step.value();
}

/// @brief Initialize the Integrator object after simulation setup is complete
///
/// This should be called at the end of initialize() after solution states have been
/// fully initialized. The Integrator takes ownership of these solution states.
///
/// @param solutions Solution states containing old acceleration, displacement, and velocity
void Simulation::initialize_integrator(SolutionStates&& solutions)
{
  integrator_ = std::make_unique<Integrator>(this, std::move(solutions));
}

/// @brief Get reference to the Integrator object
///
/// @return Reference to the Integrator
Integrator& Simulation::get_integrator()
{
  if (!integrator_) {
    throw std::runtime_error("Integrator not initialized. Call initialize_integrator() first.");
  }
  return *integrator_;
}

/// @brief Get pointer to PartitionedFSI object (null if not configured)
PartitionedFSI* Simulation::get_partitioned_fsi()
{
  return partitioned_fsi_.get();
}

/// @brief Initialize partitioned FSI if configured in parameters.
///
/// Parameters are only parsed on rank 0 (slaves skip read_files), so we
/// broadcast the active flag and config to all ranks before branching.
void Simulation::initialize_partitioned_fsi(const std::string& xml_file_path)
{
  auto& cm = com_mod.cm;
  auto& cm_mod_ref = cm_mod;

  // Rank 0 determines whether partitioned FSI is active and builds the config.
  // Broadcast the decision so all ranks take the same path.
  int active = 0;
  PartitionedFSIConfig config;

  if (cm.mas(cm_mod_ref)) {
    auto& pcp = parameters.partitioned_coupling_parameters;
    // Active when <Partitioned_coupling> is present and at least one equation
    // carries a partitioned role attribute.
    if (pcp.defined()) {
      for (auto* ep : parameters.equation_parameters) {
        if (ep->role.defined() && !ep->role.value().empty()) { active = 1; break; }
      }
    }
    if (active) {
      config.max_coupling_iterations = pcp.max_coupling_iterations.value();
      config.coupling_tolerance       = pcp.coupling_tolerance.value();
      config.initial_relaxation       = pcp.initial_relaxation.value();
      config.omega_max                = pcp.omega_max.value();

      std::string method = pcp.coupling_method.value();
      if (method == "constant")    config.coupling_method = CouplingMethod::constant;
      else if (method == "aitken") config.coupling_method = CouplingMethod::aitken;
      else throw std::runtime_error("[PartitionedFSI] Unknown Coupling_method: " + method);

      config.fluid_interface_face = pcp.fluid_interface_face.value();
      config.solid_interface_face = pcp.solid_interface_face.value();
    }
  }

  // Broadcast the active flag and config fields to all ranks.
  MPI_Bcast(&active, 1, MPI_INT, 0, cm.com());
  if (!active) return;

  int max_iter = config.max_coupling_iterations;
  double tol   = config.coupling_tolerance;
  double relax = config.initial_relaxation;
  double omax  = config.omega_max;
  int method_i = static_cast<int>(config.coupling_method);
  MPI_Bcast(&max_iter, 1, MPI_INT,    0, cm.com());
  MPI_Bcast(&tol,      1, MPI_DOUBLE, 0, cm.com());
  MPI_Bcast(&relax,    1, MPI_DOUBLE, 0, cm.com());
  MPI_Bcast(&omax,     1, MPI_DOUBLE, 0, cm.com());
  MPI_Bcast(&method_i, 1, MPI_INT,    0, cm.com());

  auto bcast_str = [&](std::string& s) {
    int len = static_cast<int>(s.size());
    MPI_Bcast(&len, 1, MPI_INT, 0, cm.com());
    s.resize(len);
    MPI_Bcast(s.data(), len, MPI_CHAR, 0, cm.com());
  };
  bcast_str(config.fluid_interface_face);
  bcast_str(config.solid_interface_face);

  config.max_coupling_iterations = max_iter;
  config.coupling_tolerance       = tol;
  config.initial_relaxation       = relax;
  config.omega_max                = omax;
  config.coupling_method          = static_cast<CouplingMethod>(method_i);

  partitioned_fsi_ = std::make_unique<PartitionedFSI>(this, config, xml_file_path);
}
