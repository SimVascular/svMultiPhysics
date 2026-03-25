// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "bcast.h"

#include "mpi.h"

namespace bcast {

/// @brief To broadcast a variable to all processors
///
/// Reproduces code in BCAST.f.
//
void fsils_bcast(double& u, FSILS_commuType& commu)
{
  if (commu.nTasks > 1) { 
    double uG;
    MPI_Allreduce(&u, &uG, 1, cm_mod::mpreal, MPI_SUM, commu.comm);
    u = uG;
  } 
}


void fsils_bcast_v(const int n, Vector<double>& u, FSILS_commuType& commu)
{
  if (commu.nTasks > 1) { 
    Vector<double> uG(n);
    MPI_Allreduce(u.data(), uG.data(), n, cm_mod::mpreal, MPI_SUM, commu.comm);
    u = uG;
  } 
}

};


