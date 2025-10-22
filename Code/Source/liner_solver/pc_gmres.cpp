// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "pc_gmres.h"

#include "fsils_api.hpp"

#include "add_bc_mul.h"
#include "bcast.h"
#include "dot.h"
#include "norm.h"
#include "omp_la.h"
#include "spar_mul.h"

#include "Array3.h"

#include <math.h>

namespace pc_gmres {

/// \todo [NOTE] Not implemented.
//
void pc_gmres(fsi_linear_solver::FSILS_lhsType& lhs, fsi_linear_solver::FSILS_subLsType& ls, const int dof, 
    const Array<double>& Val, const Array<double>& R)
{
  using namespace fsi_linear_solver;

}

};
