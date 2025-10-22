// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "fils_struct.hpp"

namespace norm {

using namespace fsi_linear_solver;

double fsi_ls_norms(const int nNo, FSILS_commuType& commu, const Vector<double>& U);

double fsi_ls_normv(const int dof, const int nNo, FSILS_commuType& commu, const Array<double>& U);

};
