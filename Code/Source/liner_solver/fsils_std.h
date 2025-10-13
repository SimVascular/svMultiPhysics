// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// This file reproduces FSILS_STD.h.

#ifndef FILS_STD_H
#define FILS_STD_H

#include "mpi.h"

namespace fsi_linear_solver {

// Set MPI data type names.
const decltype(MPI_LOGICAL) mplog  = MPI_LOGICAL;
const decltype(MPI_INTEGER) mpint = MPI_INTEGER;
const decltype(MPI_DOUBLE_PRECISION) mpreal = MPI_DOUBLE_PRECISION;
const decltype(MPI_CHARACTER) mpchar = MPI_CHARACTER;
//const decltype(MPI_STATUS_SIZE) mpsts = MPI_STATUS_SIZE;

};

#endif


