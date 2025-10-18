// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef READ_FILES_H 
#define READ_FILES_H 

#include "Simulation.h"

#include <string>

/// @brief Define some types used to pass data to functions.
///
/// \todo [TODO:DaveP] maxOutput=5 is is defined in consts but in the Fortran READEQ
/// subroutine is defined as maxOutput=22.
//
namespace read_files_ns {

  const int maxOutput = 22;
  using EquationNdop = std::array<int, 4>;
  using EquationOutputs = std::array<consts::OutputNameType, maxOutput>;
  using EquationPhys = std::vector<consts::EquationType>;
  using EquationProps = std::array<std::array<consts::PhysicalProperyType, consts::maxNProp>, 20>;

  void face_match(ComMod& com_mod, faceType& lFa, faceType& gFa, Vector<int>& ptr);

  void read_bc(Simulation* simulation, EquationParameters* eq_params, eqType& lEq, BoundaryConditionParameters* bc_params, bcType& lBc);

  void read_bct(ComMod& com_mod, MBType& lMB, faceType& lFa, const std::string& fName);

  void read_bf(ComMod& com_mod, BodyForceParameters* bf_params, bfType& lBf);

  void read_cplbc_initialization_file(const std::string& file_name, cplBCType& cplBC);

  void read_domain(Simulation* simulation, EquationParameters* eq_params, eqType& lEq, EquationProps& propL,  EquationPhys = {});

  void read_eq(Simulation* simulation, EquationParameters* params, eqType& eq);

  void read_files(Simulation* simulation, const std::string& file_name);

  void read_fourier_coeff_values_file(const std::string& file_name, bcType& lBc);
  void read_fourier_coeff_values_file(const std::string& file_name, bfType& lBf);

  void read_ls(Simulation* simulation, EquationParameters* eq_params, consts::SolverType solver_type, eqType& lEq);

  void read_mat_model(Simulation* simulation, EquationParameters* eq_params, DomainParameters* domain_params, dmnType& lDmn);

  void read_outputs(Simulation* simulation, EquationParameters* eq_params, eqType& lEq, EquationNdop& nDOP,  EquationOutputs& outPuts);

  void read_rmsh(Simulation* simulation, EquationParameters* eq_param);

  void read_spatial_values(const ComMod& com_mod, const mshType& msh, const faceType& lFa, const std::string& file_name, bcType& lBc);

  void read_temporal_values(const std::string& file_name, bcType& lBc);
  void read_temporal_values(const std::string& file_name, bfType& lBf);

  void read_temp_spat_values(const ComMod& com_mod, const mshType& msh, const faceType& lFa, 
      const std::string& file_name, bcType& lBc);
  void read_temp_spat_values(const ComMod& com_mod, const mshType& msh, const std::string& file_name, bfType& lBf);

  void read_trac_bcff(ComMod& com_mod, MBType& lMB, faceType& lFa, const std::string& file_name);

  void read_fluid_visc_model(Simulation* simulation, EquationParameters* eq_params, DomainParameters* domain_params, dmnType& lDmn);

  void read_solid_visc_model(Simulation* simulation, EquationParameters* eq_params, DomainParameters* domain_params, dmnType& lDmn);

  void read_wall_props_ff(ComMod& com_mod, const std::string& file_path, const int iM, const int iFa);

  void set_cmm_bdry(mshType& lM, Vector<int>& bNds);

  void set_equation_properties(Simulation* simulation, EquationParameters* eq_params, eqType& lEq, EquationProps& propL, 
    EquationOutputs& outPuts, EquationNdop& nDOP);


};

#endif

