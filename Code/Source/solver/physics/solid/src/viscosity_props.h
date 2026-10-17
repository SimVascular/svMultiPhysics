// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

using SetSolidViscosityPropertiesMapType = std::map<consts::SolidViscosityModelType, std::function<void(Simulation*, SolidViscosityParameters&,
    dmnType& lDmn)>>;

SetSolidViscosityPropertiesMapType set_solid_viscosity_props = {

//---------------------------//
//      viscType_Newtonian   //
//---------------------------//
//
{consts::SolidViscosityModelType::viscType_Newtonian, [](Simulation* simulation, SolidViscosityParameters& params, dmnType& lDmn) -> void
{
  using namespace consts;
  auto& com_mod = simulation->get_com_mod();

  lDmn.solid_visc.viscType = SolidViscosityModelType::viscType_Newtonian;
  lDmn.solid_visc.mu = params.newtonian_model.constant_value.value();
} },

//---------------------------//
//      viscType_Potential   //
//---------------------------//
//
{consts::SolidViscosityModelType::viscType_Potential, [](Simulation* simulation, SolidViscosityParameters& params, dmnType& lDmn) -> void
{
  using namespace consts;
  auto& com_mod = simulation->get_com_mod();

  lDmn.solid_visc.viscType = SolidViscosityModelType::viscType_Potential;
  lDmn.solid_visc.mu = params.potential_model.constant_value.value();
} },

};
