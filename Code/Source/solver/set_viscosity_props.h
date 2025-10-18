// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

///////////////////////////////////////////////////////////
//    S e t   V i s c o s i t y   P r o p e r t i e s    //
///////////////////////////////////////////////////////////
//
// The 'set_fluid_viscosity_props ' map defined here sets fluid 
// viscosity properties from values read in from a file.
// The 'set_solid_viscosity_props ' map defined here sets solid 
// viscosity properties from values read in from a file.

using SetFluidViscosityPropertiesMapType = std::map<consts::FluidViscosityModelType, std::function<void(Simulation*, FluidViscosityParameters&,
    dmnType& lDmn)>>;

SetFluidViscosityPropertiesMapType set_fluid_viscosity_props = {

//---------------------------//
//      viscType_Const       //
//---------------------------//
//
{consts::FluidViscosityModelType::viscType_Const, [](Simulation* simulation, FluidViscosityParameters& params, dmnType& lDmn) -> void
{
  using namespace consts;
  auto& com_mod = simulation->get_com_mod();

  lDmn.fluid_visc.viscType = FluidViscosityModelType::viscType_Const;
  lDmn.fluid_visc.mu_i = params.newtonian_model.constant_value.value();
} },

//---------------------------//
//      viscType_CY          //
//---------------------------//
//
{consts::FluidViscosityModelType::viscType_CY, [](Simulation* simulation, FluidViscosityParameters& params, dmnType& lDmn) -> void
{
  using namespace consts;
  auto& com_mod = simulation->get_com_mod();
  auto& model_params = params.carreau_yasuda_model;

  lDmn.fluid_visc.viscType = FluidViscosityModelType::viscType_CY;

  lDmn.fluid_visc.mu_i = model_params.limiting_high_shear_rate_viscosity.value();
  lDmn.fluid_visc.mu_o = model_params.limiting_low_shear_rate_viscosity.value();
  lDmn.fluid_visc.lam = model_params.shear_rate_tensor_multipler.value();
  lDmn.fluid_visc.a = model_params.shear_rate_tensor_exponent.value();
  lDmn.fluid_visc.n = model_params.power_law_index.value();

  if (lDmn.fluid_visc.mu_i > lDmn.fluid_visc.mu_o) {
    throw std::runtime_error("Unexpected inputs for Carreau-Yasuda model. " 
      "High shear-rate viscosity value (" + std::to_string(lDmn.fluid_visc.mu_i) + 
      " should be larger than low shear-rate value " + std::to_string(lDmn.fluid_visc.mu_i) + ".");
  }

} },

//---------------------------//
//      viscType_Cass        //
//---------------------------//
//
{consts::FluidViscosityModelType::viscType_Cass, [](Simulation* simulation, FluidViscosityParameters& params, dmnType& lDmn) -> void
{
  using namespace consts;
  auto& com_mod = simulation->get_com_mod();
  auto& model_params = params.cassons_model;

  lDmn.fluid_visc.viscType = FluidViscosityModelType::viscType_Cass;
  lDmn.fluid_visc.mu_i = model_params.asymptotic_viscosity();
  lDmn.fluid_visc.mu_o = model_params.yield_stress();

  if (model_params.low_shear_rate_threshold.defined()) { 
    lDmn.fluid_visc.lam = model_params.low_shear_rate_threshold();
  } else {
    lDmn.fluid_visc.lam = 0.5;
  }

} },

};





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


