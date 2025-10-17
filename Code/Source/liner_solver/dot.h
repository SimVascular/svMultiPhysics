// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "fils_struct.hpp"

namespace dot {

using namespace fsi_linear_solver;

double fsils_dot_s(const int nNo, FSILS_commuType& commu, const Vector<double>& U, const Vector<double>& V);

double fsils_dot_v(const int dof, const int nNo, FSILS_commuType& commu, const Array<double>& U, const Array<double>& V);

double fsils_nc_dot_s(const int nNo, const Vector<double>& U, const Vector<double>& V);

double fsils_nc_dot_v(const int dof, const int nNo, const Array<double>& U, const Array<double>& V);

};
