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
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "Simulation.h"

#include "all_fun.h"
#include "load_msh.h"

#include "mpi.h"

#include <iostream>

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
