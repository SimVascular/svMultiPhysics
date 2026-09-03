// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/// @brief The 'set_equation_props' map defined here sets equation 
/// properties from values read in from a file.
///
/// This replaces the 'SELECT CASE (eqName)' statement in the Fortran 'READEQ()' subroutine.
//
using SetEquationPropertiesMapType = std::map<consts::EquationType, std::function<void(Simulation*, EquationParameters*, 
  eqType&, EquationProps&, EquationOutputs&, EquationNdop&)>>;

//--------------------
// set_equation_props
//--------------------
//
SetEquationPropertiesMapType set_equation_props = {

//---------------------------//
//         phys_CEP          //
//---------------------------//
//
{consts::EquationType::phys_CEP, [](Simulation* simulation, EquationParameters* eq_params, eqType& lEq, EquationProps& propL,
      EquationOutputs& outPuts, EquationNdop& nDOP) -> void
{
  using namespace consts;
  auto& cep_mod = simulation->get_cep_mod();
  lEq.phys = consts::EquationType::phys_CEP;

  propL[0][0] = PhysicalPropertyType::fluid_density;
  propL[1][0] = PhysicalPropertyType::backflow_stab;
  propL[2][0] = PhysicalPropertyType::f_x;
  propL[3][0] = PhysicalPropertyType::f_y;

  if (simulation->com_mod.nsd == 3) {
    propL[4][0] = PhysicalPropertyType::f_z;
  }

  cep_mod.cepEq = true;

  read_domain(simulation, eq_params, lEq, propL);

  nDOP = {1, 1, 0, 0};
  outPuts[0] = OutputNameType::out_voltage;

  // Set solver parameters.
  read_ls(simulation, eq_params, SolverType::lSolver_CG, lEq);

} },

//---------------------------//
//         phys_CMM          //
//---------------------------//
//
{consts::EquationType::phys_CMM, [](Simulation* simulation, EquationParameters* eq_params, eqType& lEq, EquationProps& propL,
      EquationOutputs& outPuts, EquationNdop& nDOP) -> void
{ 
  using namespace consts;
  auto& com_mod = simulation->get_com_mod();
  lEq.phys = consts::EquationType::phys_CMM;

  bool pstEq = eq_params->prestress.defined() && eq_params->prestress.value(); 

  com_mod.cmmBdry.resize(com_mod.gtnNo);
  if (eq_params->initialize.defined()) {
    com_mod.cmmInit = true; 

    if (com_mod.nEq > 1) {
      throw std::runtime_error("More than one eqn. is not allowed while initializing CMM.");
    }

    // Determine is there is a pre-stress.
    //
    auto init_str = eq_params->initialize();
    std::transform(init_str.begin(), init_str.end(), init_str.begin(), ::tolower);

    if (std::set<std::string>{"inflate", "inf"}.count(init_str) != 0) {
      com_mod.pstEq = false;
    } else if (std::set<std::string>{"prestress", "prest"}.count(init_str) != 0) {
      com_mod.pstEq = true;
    } else {
      throw std::runtime_error("Unknown CMM initialize type '" + init_str + "'.");
    }

    // Set cmmBdry vector to be edge nodes of the wall
    for (int iM = 0; iM < com_mod.nMsh; iM++) {
      set_cmm_bdry(com_mod.msh[iM], com_mod.cmmBdry);
    }
  }

  // Set variable wall properties.
  //
  if (eq_params->variable_wall_properties.defined()) {
    com_mod.cmmVarWall = true;

    if (com_mod.varWallProps.size() == 0) {
      // varWallProps = array of size 2 x total number of nodes across all meshes and all processors; first column is thickness and second column is elastic modulus
      com_mod.varWallProps.resize(2, com_mod.gtnNo);
    }

    auto mesh_name = eq_params->variable_wall_properties.mesh_name.value();
    int iM = 0;
    int iFa = 0;
    if (com_mod.cmmInit) {
      all_fun::find_msh(com_mod.msh, mesh_name, iM);
    } else { 
      all_fun::find_face(com_mod.msh, mesh_name, iM, iFa);
    }
    auto file_path = eq_params->variable_wall_properties.wall_properties_file_path.value();
    read_wall_props_ff(com_mod, file_path, iM, iFa);
  }

  if (!com_mod.cmmInit) {
    propL[0][0] = PhysicalPropertyType::fluid_density;
    propL[1][0] = PhysicalPropertyType::backflow_stab;
    propL[2][0] = PhysicalPropertyType::solid_density;
    propL[3][0] = PhysicalPropertyType::poisson_ratio;
    propL[4][0] = PhysicalPropertyType::damping;

    if (!com_mod.cmmVarWall) {
      propL[5][0] = PhysicalPropertyType::shell_thickness;
      propL[6][0] = PhysicalPropertyType::elasticity_modulus;
    }

    propL[7][0] = PhysicalPropertyType::f_x;
    propL[8][0] = PhysicalPropertyType::f_y;
    if (simulation->com_mod.nsd == 3) {
      propL[9][0] = PhysicalPropertyType::f_z;
    }

    nDOP = {12, 4, 3, 0};
    outPuts = {
       OutputNameType::out_velocity,
       OutputNameType::out_pressure,
       OutputNameType::out_WSS,
       OutputNameType::out_displacement,
       OutputNameType::out_energyFlux,
       OutputNameType::out_traction,
       OutputNameType::out_vorticity,
       OutputNameType::out_vortex,
       OutputNameType::out_strainInv,
       OutputNameType::out_viscosity,
       OutputNameType::out_divergence,
       OutputNameType::out_acceleration
     };

  } else {
    propL[0][0] = PhysicalPropertyType::poisson_ratio;
    if (!com_mod.cmmVarWall) {
      propL[1][0] = PhysicalPropertyType::shell_thickness;
      propL[2][0] = PhysicalPropertyType::elasticity_modulus;
    }

    propL[7][0] = PhysicalPropertyType::f_x;
    propL[8][0] = PhysicalPropertyType::f_y;
    if (simulation->com_mod.nsd == 3) {
      propL[9][0] = PhysicalPropertyType::f_z;
    }

    if (pstEq) {
      nDOP = {2, 2, 0, 0};
      outPuts = { OutputNameType::out_displacement, OutputNameType::out_stress };
    } else { 
      nDOP = {1, 1, 0, 0};
      outPuts[0] = OutputNameType::out_displacement;
    }
  }

  read_domain(simulation, eq_params, lEq, propL);

  if (com_mod.cmmInit) {
    for (auto& domain : lEq.dmn) {
      domain.prop[PhysicalPropertyType::solid_density] = 0.0;
    }
  }

  // Set solver parameters.
  read_ls(simulation, eq_params, SolverType::lSolver_GMRES, lEq);

} },

//---------------------------//
//        phys_fluid         //
//---------------------------//
//
{consts::EquationType::phys_fluid, [](Simulation* simulation, EquationParameters* eq_params, eqType& lEq, EquationProps& propL, 
      EquationOutputs& outPuts, EquationNdop& nDOP) -> void 
{
  using namespace consts;
  auto& com_mod = simulation->get_com_mod();
  lEq.phys = consts::EquationType::phys_fluid;

  propL[0][0] = PhysicalPropertyType::fluid_density;
  propL[1][0] = PhysicalPropertyType::backflow_stab;
  propL[2][0] = PhysicalPropertyType::brinkman_inverse_permeability;
  propL[3][0] = PhysicalPropertyType::f_x;
  propL[4][0] = PhysicalPropertyType::f_y;

  if (simulation->com_mod.nsd == 3) {
    propL[5][0] = PhysicalPropertyType::f_z;
  }

  // Set fluid domain properties.
  read_domain(simulation, eq_params, lEq, propL);

  nDOP = {11, 2, 3, 0};

  outPuts = { 
    OutputNameType::out_velocity,
    OutputNameType::out_pressure,
    OutputNameType::out_WSS,
    OutputNameType::out_traction,
    OutputNameType::out_vorticity,
    OutputNameType::out_vortex,
    OutputNameType::out_strainInv,
    OutputNameType::out_energyFlux,
    OutputNameType::out_viscosity,
    OutputNameType::out_divergence,
    OutputNameType::out_acceleration
  };

  // Set solver parameters.
  read_ls(simulation, eq_params, SolverType::lSolver_NS, lEq);

} },

//---------------------------//
//        phys_heatF         //
//---------------------------//
//
{consts::EquationType::phys_heatF, [](Simulation* simulation, EquationParameters* eq_params, eqType& lEq, EquationProps& propL,
      EquationOutputs& outPuts, EquationNdop& nDOP) -> void
{
  using namespace consts;
  auto& com_mod = simulation->get_com_mod();
  lEq.phys = consts::EquationType::phys_heatF;

  propL[0][0] = PhysicalPropertyType::conductivity;
  propL[1][0] = PhysicalPropertyType::source_term;

  read_domain(simulation, eq_params, lEq, propL);

  nDOP = {3,1,1,0};
  outPuts = {OutputNameType::out_temperature,
             OutputNameType::out_heatFlux,
             OutputNameType::out_velocity};

  // Set solver parameters.
  read_ls(simulation, eq_params, SolverType::lSolver_GMRES, lEq);

} },

//---------------------------//
//        phys_heatS         //
//---------------------------//
//
{consts::EquationType::phys_heatS, [](Simulation* simulation, EquationParameters* eq_params, eqType& lEq, EquationProps& propL,
      EquationOutputs& outPuts, EquationNdop& nDOP) -> void
{ 
  using namespace consts;
  auto& com_mod = simulation->get_com_mod();
  lEq.phys = consts::EquationType::phys_heatS;
  
  propL[0][0] = PhysicalPropertyType::conductivity;
  propL[1][0] = PhysicalPropertyType::source_term;
  propL[2][0] = PhysicalPropertyType::solid_density;
  
  read_domain(simulation, eq_params, lEq, propL);

  nDOP = {2,1,1,0};
  outPuts = {OutputNameType::out_temperature, OutputNameType::out_heatFlux};

  // Set solver parameters.
  read_ls(simulation, eq_params, SolverType::lSolver_CG, lEq);

} },

//---------------------------//
//        phys_darcy         //
//---------------------------//
{consts::EquationType::phys_darcy, [](Simulation* simulation, EquationParameters* eq_params, eqType& lEq, EquationProps& propL,
                                      EquationOutputs& outPuts, EquationNdop& nDOP) -> void
{
    using namespace consts;
    auto& com_mod = simulation->get_com_mod();
    lEq.phys = consts::EquationType::phys_darcy;

    propL[0][0] = PhysicalPropertyType::darcy_permeability;
    propL[1][0] = PhysicalPropertyType::source_term;
    propL[2][0] = PhysicalPropertyType::fluid_density;
    propL[3][0] = PhysicalPropertyType::darcy_media_compressibility;
    propL[4][0] = PhysicalPropertyType::darcy_fluid_viscosity;

    read_domain(simulation, eq_params, lEq, propL);

    for (const auto& domain : lEq.dmn) {
      darcy::validate_material_properties(domain);
    }

    nDOP = {2,1,1,0};
    outPuts = {OutputNameType::out_darcyPressure, OutputNameType::out_darcyFlux};

    // Set solver parameters.
    read_ls(simulation, eq_params, SolverType::lSolver_CG, lEq);
} },


//---------------------------//
//         phys_FSI          //
//---------------------------//
//
{consts::EquationType::phys_FSI, [](Simulation* simulation, EquationParameters* eq_params, eqType& lEq, EquationProps& propL,
      EquationOutputs& outPuts, EquationNdop& nDOP) -> void
{ 
  using namespace consts;
  auto& com_mod = simulation->get_com_mod();
  lEq.phys = consts::EquationType::phys_FSI;
  com_mod.mvMsh = true;

  // Set the possible equations for fsi: fluid (required), struct/ustruct/lElas
  EquationPhys phys { EquationType::phys_fluid, EquationType::phys_struct, EquationType::phys_ustruct, EquationType::phys_lElas };
  
  // Set fluid properties.
  int n = 0;
  propL[0][n] = PhysicalPropertyType::fluid_density;
  propL[1][n] = PhysicalPropertyType::backflow_stab;
  propL[2][n] = PhysicalPropertyType::f_x;
  propL[3][n] = PhysicalPropertyType::f_y;
  if (simulation->com_mod.nsd == 3) {
    propL[4][n] = PhysicalPropertyType::f_z;
  }

  // Set struct properties.
  n += 1;
  propL[0][n] = PhysicalPropertyType::solid_density;
  propL[1][n] = PhysicalPropertyType::elasticity_modulus;
  propL[2][n] = PhysicalPropertyType::poisson_ratio;
  propL[3][n] = PhysicalPropertyType::damping;
  propL[4][n] = PhysicalPropertyType::f_x;
  propL[5][n] = PhysicalPropertyType::f_y;
  if (simulation->com_mod.nsd == 3) {
    propL[6][n] = PhysicalPropertyType::f_z;
  }

  // Set ustruct properties.
  n += 1;
  propL[0][n] = PhysicalPropertyType::solid_density;
  propL[1][n] = PhysicalPropertyType::elasticity_modulus;
  propL[2][n] = PhysicalPropertyType::poisson_ratio;
  propL[3][n] = PhysicalPropertyType::ctau_M;
  propL[4][n] = PhysicalPropertyType::ctau_C;
  propL[5][n] = PhysicalPropertyType::f_x;
  propL[6][n] = PhysicalPropertyType::f_y;
  if (simulation->com_mod.nsd == 3) {
    propL[7][n] = PhysicalPropertyType::f_z;
  }

  // Set lElas properties.
  n += 1;
  propL[0][n] = PhysicalPropertyType::solid_density;
  propL[1][n] = PhysicalPropertyType::elasticity_modulus;
  propL[2][n] = PhysicalPropertyType::poisson_ratio;
  propL[3][n] = PhysicalPropertyType::f_x;
  propL[4][n] = PhysicalPropertyType::f_y;
  if (simulation->com_mod.nsd == 3) {
    propL[5][n] = PhysicalPropertyType::f_z;
  }

  // Set lEq properties.
  read_domain(simulation, eq_params, lEq, propL, phys);

  nDOP = {22, 4, 2, 0};
  outPuts = {
    OutputNameType::out_velocity,
    OutputNameType::out_pressure,
    OutputNameType::out_displacement,
    OutputNameType::out_mises,

    OutputNameType::out_WSS,
    OutputNameType::out_traction,
    OutputNameType::out_vorticity,
    OutputNameType::out_vortex,
    OutputNameType::out_strainInv,
    OutputNameType::out_energyFlux,
    OutputNameType::out_viscosity,
    OutputNameType::out_absVelocity,
    OutputNameType::out_stress,
    OutputNameType::out_cauchy,
    OutputNameType::out_strain,
    OutputNameType::out_jacobian,
    OutputNameType::out_defGrad,
    OutputNameType::out_integ,
    OutputNameType::out_fibDir,
    OutputNameType::out_fibAlign,

    OutputNameType::out_divergence,
    OutputNameType::out_acceleration
  };

  // Set solver parameters.
  read_ls(simulation, eq_params, SolverType::lSolver_GMRES, lEq);

  if (com_mod.rmsh.isReqd && !com_mod.resetSim) {
    read_rmsh(simulation, eq_params);
  }

} },

//---------------------------//
//        phys_lElas         //
//---------------------------//
//
{consts::EquationType::phys_lElas, [](Simulation* simulation, EquationParameters* eq_params, eqType& lEq, EquationProps& propL,
      EquationOutputs& outPuts, EquationNdop& nDOP) -> void
{ 
  using namespace consts;
  auto& com_mod = simulation->get_com_mod();
  lEq.phys = consts::EquationType::phys_lElas;
  
  propL[0][0] = PhysicalPropertyType::solid_density;
  propL[1][0] = PhysicalPropertyType::elasticity_modulus;
  propL[2][0] = PhysicalPropertyType::poisson_ratio;
  propL[3][0] = PhysicalPropertyType::f_x;
  propL[4][0] = PhysicalPropertyType::f_y;
  if (simulation->com_mod.nsd == 3) {
    propL[5][0] = PhysicalPropertyType::f_z;
  }

  read_domain(simulation, eq_params, lEq, propL);

  if (eq_params->prestress.defined() && eq_params->prestress.value()) { 
    nDOP = {3,2,0,0};
    outPuts = {OutputNameType::out_displacement, OutputNameType::out_stress, OutputNameType::out_strain};
  } else {
    nDOP = {8,2,0,0};
    outPuts = {
      OutputNameType::out_displacement, OutputNameType::out_mises, OutputNameType::out_stress,
      OutputNameType::out_strain, OutputNameType::out_velocity, OutputNameType::out_acceleration,
      OutputNameType::out_integ, OutputNameType::out_jacobian 
    };
  }

  // Set solver parameters.
  read_ls(simulation, eq_params, SolverType::lSolver_CG, lEq);

} },

//---------------------------//
//          phys_mesh        //
//---------------------------//
//
{consts::EquationType::phys_mesh, [](Simulation* simulation, EquationParameters* eq_params, eqType& lEq, EquationProps& propL,
      EquationOutputs& outPuts, EquationNdop& nDOP) -> void
{
  using namespace consts;
  auto& com_mod = simulation->get_com_mod();
  lEq.phys = consts::EquationType::phys_mesh;

  propL[0][0] = PhysicalPropertyType::solid_density;
  propL[1][0] = PhysicalPropertyType::elasticity_modulus;
  propL[2][0] = PhysicalPropertyType::poisson_ratio;
  propL[3][0] = PhysicalPropertyType::f_x;
  propL[4][0] = PhysicalPropertyType::f_y;
  if (simulation->com_mod.nsd == 3) {
    propL[5][0] = PhysicalPropertyType::f_z;
  }

  read_domain(simulation, eq_params, lEq, propL);

  for (auto& domain : lEq.dmn) {
      domain.prop[PhysicalPropertyType::solid_density] = 0.0;
      domain.prop[PhysicalPropertyType::elasticity_modulus] = 1.0;
  }

  nDOP = {3, 1, 0, 0};
  outPuts = {OutputNameType::out_displacement, OutputNameType::out_velocity, OutputNameType::out_acceleration };

  lEq.ls.relTol = 0.2;

  // Set solver parameters.
  read_ls(simulation, eq_params, SolverType::lSolver_CG, lEq);

} },

//---------------------------//
//          phys_shell       //
//---------------------------//
//
{consts::EquationType::phys_shell, [](Simulation* simulation, EquationParameters* eq_params, eqType& lEq, EquationProps& propL,
      EquationOutputs& outPuts, EquationNdop& nDOP) -> void
{ 
  using namespace consts;
  auto& com_mod = simulation->get_com_mod();
  lEq.phys = consts::EquationType::phys_shell;
  com_mod.shlEq = true;
  
  propL[0][0] = PhysicalPropertyType::solid_density;
  propL[1][0] = PhysicalPropertyType::damping;
  propL[2][0] = PhysicalPropertyType::elasticity_modulus;
  propL[3][0] = PhysicalPropertyType::poisson_ratio;
  propL[4][0] = PhysicalPropertyType::shell_thickness;
  propL[5][0] = PhysicalPropertyType::f_x;
  propL[6][0] = PhysicalPropertyType::f_y;
  propL[7][0] = PhysicalPropertyType::f_z;
  
  read_domain(simulation, eq_params, lEq, propL);
  
  nDOP = {9,1,0,0};
  outPuts = {
    OutputNameType::out_displacement, 
    OutputNameType::out_stress, 
    OutputNameType::out_strain, 
    OutputNameType::out_jacobian, 
    OutputNameType::out_defGrad, 
    OutputNameType::out_velocity, 
    OutputNameType::out_integ,
    OutputNameType::out_CGstrain,
    OutputNameType::out_CGInv1
  };

  // Set solver parameters.
  read_ls(simulation, eq_params, SolverType::lSolver_CG, lEq);

} },

//---------------------------//
//        phys_stokes        //
//---------------------------//
//
{consts::EquationType::phys_stokes, [](Simulation* simulation, EquationParameters* eq_params, eqType& lEq, EquationProps& propL,
      EquationOutputs& outPuts, EquationNdop& nDOP) -> void
{
  using namespace consts;
  auto& com_mod = simulation->get_com_mod();
  lEq.phys = consts::EquationType::phys_stokes;

  propL[0][0] = PhysicalPropertyType::ctau_M;
  propL[1][0] = PhysicalPropertyType::f_x;
  propL[2][0] = PhysicalPropertyType::f_y;
  if (simulation->com_mod.nsd == 3) {
    propL[3][0] = PhysicalPropertyType::f_z;
  }
  read_domain(simulation, eq_params, lEq, propL);

  nDOP = {8, 2, 3, 0};
  outPuts = {
    OutputNameType::out_velocity,
    OutputNameType::out_pressure,
    OutputNameType::out_WSS,
    OutputNameType::out_vorticity,
    OutputNameType::out_traction,
    OutputNameType::out_strainInv,
    OutputNameType::out_viscosity,
    OutputNameType::out_divergence
  };

  // Set solver parameters.
  read_ls(simulation, eq_params, SolverType::lSolver_GMRES, lEq);

} },

//---------------------------//
//          phys_struct      //
//---------------------------//
//
{consts::EquationType::phys_struct, [](Simulation* simulation, EquationParameters* eq_params, eqType& lEq, EquationProps& propL,
      EquationOutputs& outPuts, EquationNdop& nDOP) -> void
{
  using namespace consts;
  auto& com_mod = simulation->get_com_mod();
  lEq.phys = consts::EquationType::phys_struct;

  propL[0][0] = PhysicalPropertyType::solid_density;
  propL[1][0] = PhysicalPropertyType::damping;
  propL[2][0] = PhysicalPropertyType::elasticity_modulus;
  propL[3][0] = PhysicalPropertyType::poisson_ratio;
  propL[4][0] = PhysicalPropertyType::f_x;
  propL[5][0] = PhysicalPropertyType::f_y;
  if (simulation->com_mod.nsd == 3) {
    propL[6][0] = PhysicalPropertyType::f_z;
  }

  read_domain(simulation, eq_params, lEq, propL);

  if (eq_params->prestress.defined() && eq_params->prestress.value()) { 
    nDOP = {4,2,0,0};
    outPuts = {OutputNameType::out_displacement, OutputNameType::out_stress, OutputNameType::out_cauchy, OutputNameType::out_strain};
    //simulation->com_mod.pstEq = true;
  } else {
    nDOP = {17, 2, 0, 0};
    outPuts = {OutputNameType::out_displacement,
               OutputNameType::out_mises,
               OutputNameType::out_stress,
               OutputNameType::out_cauchy,
               OutputNameType::out_strain,
               OutputNameType::out_jacobian,
               OutputNameType::out_defGrad,
               OutputNameType::out_integ,
               OutputNameType::out_fibDir,
               OutputNameType::out_fibAlign,
               OutputNameType::out_velocity,
               OutputNameType::out_acceleration,
               OutputNameType::out_fibStretch,
               OutputNameType::out_fibStretchRate,
               OutputNameType::out_activeTensionFibers,
               OutputNameType::out_activeTensionSheets,
               OutputNameType::out_activeTensionNormal};
  }

  // Set solver parameters.
  read_ls(simulation, eq_params, SolverType::lSolver_CG, lEq);

} },

//---------------------------//
//        phys_ustruct       //
//---------------------------//
//
{consts::EquationType::phys_ustruct, [](Simulation* simulation, EquationParameters* eq_params, eqType& lEq, EquationProps& propL,
      EquationOutputs& outPuts, EquationNdop& nDOP) -> void
{ 
  using namespace consts;
  auto& com_mod = simulation->get_com_mod();

  lEq.phys = consts::EquationType::phys_ustruct;
  com_mod.sstEq = true;
  
  propL[0][0] = PhysicalPropertyType::solid_density;
  propL[1][0] = PhysicalPropertyType::elasticity_modulus;
  propL[2][0] = PhysicalPropertyType::poisson_ratio;
  propL[3][0] = PhysicalPropertyType::ctau_M;
  propL[4][0] = PhysicalPropertyType::ctau_C;
  propL[5][0] = PhysicalPropertyType::f_x;
  propL[6][0] = PhysicalPropertyType::f_y;
  if (simulation->com_mod.nsd == 3) {
    propL[7][0] = PhysicalPropertyType::f_z;
  }

  read_domain(simulation, eq_params, lEq, propL);

  nDOP = {19, 2, 0, 0};
  outPuts = {OutputNameType::out_displacement,
             OutputNameType::out_mises,
             OutputNameType::out_stress,
             OutputNameType::out_cauchy,
             OutputNameType::out_strain,
             OutputNameType::out_jacobian,
             OutputNameType::out_defGrad,
             OutputNameType::out_integ,
             OutputNameType::out_fibDir,
             OutputNameType::out_fibAlign,
             OutputNameType::out_velocity,
             OutputNameType::out_pressure,
             OutputNameType::out_acceleration,
             OutputNameType::out_divergence,
             OutputNameType::out_fibStretch,
             OutputNameType::out_fibStretchRate,
             OutputNameType::out_activeTensionFibers,
             OutputNameType::out_activeTensionSheets,
             OutputNameType::out_activeTensionNormal};

  // Set solver parameters.
  read_ls(simulation, eq_params, SolverType::lSolver_GMRES, lEq);

} },
};

