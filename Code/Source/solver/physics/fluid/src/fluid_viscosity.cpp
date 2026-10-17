// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "read_files.h"

#include <algorithm>
#include <cctype>
#include <functional>
#include <map>
#include <stdexcept>
#include <string>

namespace read_files_ns {

#include "viscosity_props.h"

// Set the fluid viscosity material model parameters for the given domain.
void read_fluid_visc_model(Simulation* simulation, EquationParameters* eq_params, DomainParameters* domain_params, dmnType& lDmn)
{
  using namespace consts;

  FluidViscosityModelType vmodel_type;
  std::string vmodel_str;

  if (domain_params->fluid_viscosity.model.defined()) {
    vmodel_str = domain_params->fluid_viscosity.model.value();
    std::transform(vmodel_str.begin(), vmodel_str.end(), vmodel_str.begin(), ::tolower);

    try {
      vmodel_type = fluid_viscosity_model_name_to_type.at(vmodel_str);
    } catch (const std::out_of_range& exception) {
      throw std::runtime_error("Unknown fluid viscosity model '" + vmodel_str + "'.");
    }
  } else {
    vmodel_type = FluidViscosityModelType::viscType_Const;
  }

  auto& viscosity_params = domain_params->fluid_viscosity;

  try {
    set_fluid_viscosity_props[vmodel_type](simulation, viscosity_params, lDmn);
  } catch (const std::bad_function_call& exception) {
    throw std::runtime_error("[read_fluid_visc_model] Viscosity model '" + vmodel_str + "' is not supported.");
  }

  if ((lDmn.phys == EquationType::phys_stokes) && (lDmn.fluid_visc.viscType != FluidViscosityModelType::viscType_Const)) {
    throw std::runtime_error("Only constant viscosity is allowed for Stokes flow.");
  }
}

} // namespace read_files_ns
