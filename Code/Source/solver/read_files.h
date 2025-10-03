/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

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

  //void copyBC(ComMod& com_mod, const bcType& oBc, bcType& nBc);

  void addCapBC(ComMod& com_mod, eqType& lEq, int iBc);


};

#endif

