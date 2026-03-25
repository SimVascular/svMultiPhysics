// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// The functions here reproduce the subroutines defined in svFSILS/COMMU.f.

#include "fsils.hpp"
#include "CmMod.h"

#include "mpi.h"

namespace fsi_linear_solver {

/// @brief Modifies:
///   commu.task
///   commu.nTasks
//   
void fsils_commu_create(FSILS_commuType& commu, cm_mod::MpiCommWorldType commi)
{
  // Some of these parameters are set for sequential version
  //
  commu.foC = true;
  commu.comm = commi;
  commu.nTasks = 1;
  commu.task   = 0;
  commu.master = 0;

  MPI_Comm_rank(commi, &commu.task);
  MPI_Comm_size(commi, &commu.nTasks);

  MPI_Allreduce(&commu.task, &commu.master, 1, cm_mod::mpint, MPI_MIN, commi);

  commu.masF = false;
  commu.tF = commu.task;
  if (commu.task == commu.master) {
    commu.masF = true; 
  }
}


void fsils_commu_free(FSILS_commuType& commu)
{
  //IF (.NOT.commu%foC) STOP 'COMMU is not created yet to be freed'

  commu.foC = false;
}

};


