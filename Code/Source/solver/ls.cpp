// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// Subroutines related to initializing linear solver arrays and
// function calls to svFSILS and Trilinos solver library

#include "ls.h"

#include "fsils_api.hpp"
#include "consts.h"

#include <math.h>

namespace ls_ns {

/// @brief Allocate com_mod.R and com_mod.Val arrays.
///
/// Modifies:
///    com_mod.R - Residual vector
///    com_mod.Val - LHS matrix 
///
/// Reproduces 'SUBROUTINE LSALLOC(lEq)'.
//
void ls_alloc(ComMod& com_mod, eqType& lEq)
{
  int dof = com_mod.dof;
  int tnNo = com_mod.tnNo;
  int gtnNo = com_mod.gtnNo;
  auto& lhs = com_mod.lhs;

  com_mod.R.resize(dof,tnNo);

  lEq.linear_algebra->alloc(com_mod, lEq);
}

/// @brief Modifies:    
///  com_mod.R      // Residual vector
///  com_mod.Val    // LHS matrix
///
/// Reproduces ' SUBROUTINE LSSOLVE(lEq, incL, res)'.
//
void ls_solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res) 
{
  #define n_debug_ls_solve
  #ifdef debug_ls_solve 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "lEq.sym: " << lEq.sym;
  dmsg << "lEq.useTLS: " << lEq.useTLS;
  dmsg << "lEq.assmTLS: " << lEq.assmTLS;
  #endif

  lEq.linear_algebra->solve(com_mod, lEq, incL, res);
}

};


