// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "fils_struct.hpp"

namespace bcast {

using namespace fsi_linear_solver;

void fsils_bcast(double& u, FSILS_commuType& commu);

void fsils_bcast_v(const int n, Vector<double>& u, FSILS_commuType& commu);

};
