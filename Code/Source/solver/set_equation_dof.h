// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include <map>
#include <tuple>

/// @brief The 'equation_dof_map' map defined here sets equation dof and sym data members. 
//
using EquationDofType = std::tuple<int, std::string>; 

std::map<consts::EquationType, EquationDofType> equation_dof_map =
{
  {EquationType::phys_fluid,    std::make_tuple(nsd+1, "NS") },
  {EquationType::phys_heatF,    std::make_tuple(1,     "HF") },
  {EquationType::phys_heatS,    std::make_tuple(1,     "HS") },
  {EquationType::phys_lElas,    std::make_tuple(nsd,   "LE") },
  {EquationType::phys_struct,   std::make_tuple(nsd,   "ST") },
  {EquationType::phys_ustruct,  std::make_tuple(nsd+1, "ST") },
  {EquationType::phys_CMM,      std::make_tuple(nsd+1, "CM") },
  {EquationType::phys_shell,    std::make_tuple(nsd,   "SH") },
  {EquationType::phys_FSI,      std::make_tuple(nsd+1, "FS") },
  {EquationType::phys_mesh,     std::make_tuple(nsd,   "MS") },
  {EquationType::phys_CEP,      std::make_tuple(1,     "EP") },
  {EquationType::phys_stokes,   std::make_tuple(nsd+1, "SS") }
};

