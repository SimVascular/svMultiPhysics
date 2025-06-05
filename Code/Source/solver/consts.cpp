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

#include "consts.h" 

namespace consts {

/// @brief Reproduces the 'SELECT CASE (lM%eType)' statements in the
/// Fortran 'PARTMSH' subroutine. 
//
const std::map<ElementType,int> element_type_to_elem_nonb = 
{
  {ElementType::LIN1,  1 },
  {ElementType::LIN2,  1 },
  {ElementType::TRI3,  2 },
  {ElementType::TRI6,  3 },
  {ElementType::QUD4,  2 },
  {ElementType::QUD8,  3 },
  {ElementType::QUD9,  3 },
  {ElementType::TET4,  3 },
  {ElementType::TET10, 6 },
  {ElementType::HEX8,  4 },
  {ElementType::HEX20, 8 },
  {ElementType::HEX27, 9 },
  {ElementType::WDG,   3 }
};

/// @brief Map containing the element type as a string. 
///
const std::map<ElementType,std::string> element_type_to_string = 
{
  {ElementType::NA,    "unknown" },
  {ElementType::PNT,   "point" },
  {ElementType::LIN1,  "line1" },
  {ElementType::LIN2,  "line2" },
  {ElementType::TRI3,  "tri3" },
  {ElementType::TRI6,  "tri6" },
  {ElementType::QUD4,  "quad4" },
  {ElementType::QUD8,  "quad8" },
  {ElementType::QUD9,  "qud9" },
  {ElementType::TET4,  "tet4" },
  {ElementType::TET10, "tet10" },
  {ElementType::HEX8,  "hex8" },
  {ElementType::HEX20, "hex20" },
  {ElementType::HEX27, "hex27" },
  {ElementType::WDG,   "wedge" }
};

/// @brief Map containing the dimension for each element type. 
///
const std::map<ElementType,int> element_dimension = 
{
  {ElementType::PNT,   0 },
  {ElementType::LIN1,  1 },
  {ElementType::LIN2,  1 },
  {ElementType::TRI3,  2 },
  {ElementType::TRI6,  2 },
  {ElementType::QUD4,  2 },
  {ElementType::QUD8,  2 },
  {ElementType::QUD9,  2 },
  {ElementType::TET4,  3 },
  {ElementType::TET10, 3 },
  {ElementType::HEX8,  3 },
  {ElementType::HEX20, 3 },
  {ElementType::HEX27, 3 },
  {ElementType::WDG,   3 }
};

/// @brief Map for constitutive_model string to ConstitutiveModelType. 
///
/// Type of constitutive model (isochoric) for structure equation:
/// St.Venant-Kirchhoff, modified St.Venant-Kirchhoff, NeoHookean,
/// Mooney-Rivlin, modified Holzapfel-Gasser-Ogden with dispersion,
/// Linear model (S = mu*I), Guccione (1995), Holzapfel & Ogden model
//
const std::map<std::string,ConstitutiveModelType> constitutive_model_name_to_type = 
{
  {"lin", ConstitutiveModelType::stIso_lin}, 
  {"linear", ConstitutiveModelType::stIso_lin},

  {"stVK", ConstitutiveModelType::stIso_StVK},
  {"stVenantKirchhoff", ConstitutiveModelType::stIso_StVK},

  {"m-stVK", ConstitutiveModelType::stIso_mStVK},
  {"modified-stVK", ConstitutiveModelType::stIso_mStVK},
  {"modified-stVenantKirchhoff", ConstitutiveModelType::stIso_mStVK},

  {"nHK", ConstitutiveModelType::stIso_nHook}, 
  {"nHK91", ConstitutiveModelType::stIso_nHook},
  {"neoHookean", ConstitutiveModelType::stIso_nHook},
  {"neoHookeanSimo91", ConstitutiveModelType::stIso_nHook},

  {"MR", ConstitutiveModelType::stIso_MR}, 
  {"Mooney-Rivlin", ConstitutiveModelType::stIso_MR},

  {"HGO", ConstitutiveModelType::stIso_HGO}, 

  {"Guccione", ConstitutiveModelType::stIso_Gucci},
  {"Gucci", ConstitutiveModelType::stIso_Gucci},

  {"HO", ConstitutiveModelType::stIso_HO}, 
  {"HolzapfelOgden", ConstitutiveModelType::stIso_HO},

  {"HO_ma", ConstitutiveModelType::stIso_HO_ma}, 
  {"HolzapfelOgden-ModifiedAnisotropy", ConstitutiveModelType::stIso_HO_ma},

  {"quad", ConstitutiveModelType::stVol_Quad},
  {"Quad", ConstitutiveModelType::stVol_Quad},
  {"quadratic", ConstitutiveModelType::stVol_Quad},
  {"Quadratic",ConstitutiveModelType::stVol_Quad},

  {"ST91", ConstitutiveModelType::stVol_ST91},
  {"Simo-Taylor91", ConstitutiveModelType::stVol_ST91},

  {"M94", ConstitutiveModelType::stVol_M94},
  {"Miehe94", ConstitutiveModelType::stVol_M94},

};

/// @brief Map for contact model string name to ContacteModelType
//
const std::map<std::string,ContactModelType> contact_model_name_to_type =
{
  {"penalty", ContactModelType::cntctM_penalty},
  {"potential", ContactModelType::cntctM_potential},
};


/// @brief Map for fluid viscosity model string to FluidViscosityModelType. 
//
const std::map<std::string,FluidViscosityModelType> fluid_viscosity_model_name_to_type
{
  {"const", FluidViscosityModelType::viscType_Const},
  {"constant", FluidViscosityModelType::viscType_Const},
  {"newtonian", FluidViscosityModelType::viscType_Const},

  {"cy", FluidViscosityModelType::viscType_CY},
  {"carreau-yasuda", FluidViscosityModelType::viscType_CY},

  {"cass", FluidViscosityModelType::viscType_Cass},
  {"cassons", FluidViscosityModelType::viscType_Cass},

};

/// @brief Map for solid viscosity model string to SolidViscosityModelType.
//
const std::map<std::string,SolidViscosityModelType> solid_viscosity_model_name_to_type
{
  {"Newtonian", SolidViscosityModelType::viscType_Newtonian},
  {"newtonian", SolidViscosityModelType::viscType_Newtonian},
  {"Newt", SolidViscosityModelType::viscType_Newtonian},
  {"newt", SolidViscosityModelType::viscType_Newtonian},

  {"Potential", SolidViscosityModelType::viscType_Potential},
  {"potential", SolidViscosityModelType::viscType_Potential},
  {"Pot", SolidViscosityModelType::viscType_Potential},
  {"pot", SolidViscosityModelType::viscType_Potential},
};


/// @brief Map number of element nodes to element type.
///
/// Note that this is not an ambiguous mapping between, 
/// need to handle some overlap in number of element nodes.
//
static const std::map<int,ElementType> num_nodes_to_type = {
    {1, ElementType::PNT},
    {2, ElementType::LIN1}, 
    {2, ElementType::LIN2}, 
    {3, ElementType::TRI3},
    {6, ElementType::TRI6},
    {4, ElementType::QUD4},
    {8, ElementType::QUD8},
    {9, ElementType::QUD9},
    {4, ElementType::TET4},
    {10, ElementType::TET10}, 
    {8, ElementType::HEX8},
    {20, ElementType::HEX20},
    {27, ElementType::HEX27},
    {6, ElementType::WDG}, 
    {0, ElementType::NRB}
  };

const std::map<std::string,CplBCType> cplbc_name_to_type = {
    {"E", CplBCType::cplBC_E},
    {"explicit", CplBCType::cplBC_E},
    {"I", CplBCType::cplBC_I},
    {"implicit", CplBCType::cplBC_I},
    {"SI", CplBCType::cplBC_SI},
    {"semi-implicit", CplBCType::cplBC_SI}
  };

/// @brief Map equation name to a type
//
const std::map<std::string,EquationType> equation_name_to_type = {
    {"advection_diffusion", EquationType::phys_heatF},
    {"AD", EquationType::phys_heatF},
    {"heatF", EquationType::phys_heatF},
    {"dyeTransport", EquationType::phys_heatF},
    {"scalarTransport", EquationType::phys_heatF},

    {"cardiac_electro_physiology", EquationType::phys_CEP},
    {"CEP", EquationType::phys_CEP},

    {"coupled_momentum", EquationType::phys_CMM},
    {"CMM", EquationType::phys_CMM},

    {"fluid", EquationType::phys_fluid},

    {"fluid-solid-interaction", EquationType::phys_FSI},
    {"FSI", EquationType::phys_FSI},

    {"linear_elasticity", EquationType::phys_lElas},
    {"lElas", EquationType::phys_lElas},

    {"mesh", EquationType::phys_mesh},

    {"shell", EquationType::phys_shell},

    {"solid_heat", EquationType::phys_heatS},
    {"heatS", EquationType::phys_heatS},
    {"laplace", EquationType::phys_heatS},
    {"poisson", EquationType::phys_heatS},

    {"stokes", EquationType::phys_stokes},

    {"structural", EquationType::phys_struct},
    {"struct", EquationType::phys_struct},

    {"structural_velocity_pressure", EquationType::phys_ustruct},
    {"ustruct", EquationType::phys_ustruct},

  };

// Map from output type string to OutputType.
//
const std::map<std::string,OutputType> output_type_name_to_type = {
  {"B_INT", OutputType::boundary_integral}, 
  {"Boundary_integral", OutputType::boundary_integral},
  {"Spatial", OutputType::spatial}, 
  {"V_INT", OutputType::volume_integral}, 
  {"Volume_integral", OutputType::volume_integral}
};

const std::map<std::string,MeshGeneratorType> mesh_generator_name_to_type = {
    {"Tetgen", MeshGeneratorType::RMSH_TETGEN},
    {"Meshsim", MeshGeneratorType::RMSH_MESHSIM}
};

/// @brief The list of Trilinos preconditioners. 
const std::set<PreconditionerType> trilinos_preconditioners = {
  PreconditionerType::PREC_TRILINOS_DIAGONAL,
  PreconditionerType::PREC_TRILINOS_BLOCK_JACOBI,
  PreconditionerType::PREC_TRILINOS_ILU,
  PreconditionerType::PREC_TRILINOS_ILUT,
  PreconditionerType::PREC_TRILINOS_IC,
  PreconditionerType::PREC_TRILINOS_ICT,
  PreconditionerType::PREC_TRILINOS_ML
};

/// @brief The list of FSILS preconditioners. 
const std::set<PreconditionerType> fsils_preconditioners = {
  PreconditionerType::PREC_FSILS,
  PreconditionerType::PREC_RCS
};

/// @brief The list of PETSc preconditioners. 
const std::set<PreconditionerType> petsc_preconditioners = {
  PreconditionerType::PREC_PETSC_JACOBI,
  PreconditionerType::PREC_PETSC_RCS
};

/// @brief Map for preconditioner type string to PreconditionerType enum
//
const std::map<std::string,PreconditionerType> preconditioner_name_to_type =
{
  {"none", PreconditionerType::PREC_NONE},

  {"fsils", PreconditionerType::PREC_FSILS},
  {"rcs", PreconditionerType::PREC_RCS},
  {"row-column-scaling", PreconditionerType::PREC_RCS},

  {"trilinos-diagonal", PreconditionerType::PREC_TRILINOS_DIAGONAL},
  {"trilinos-blockjacobi", PreconditionerType::PREC_TRILINOS_BLOCK_JACOBI},
  {"trilinos-ilu", PreconditionerType::PREC_TRILINOS_ILU},
  {"trilinos-ilut", PreconditionerType::PREC_TRILINOS_ILUT},
  {"trilinos-ic", PreconditionerType::PREC_TRILINOS_IC},
  {"trilinos-ict", PreconditionerType::PREC_TRILINOS_ICT},
  {"trilinos-ml", PreconditionerType::PREC_TRILINOS_ML},

  {"petsc-jacobi", PreconditionerType::PREC_PETSC_JACOBI},
  {"petsc-rcs", PreconditionerType::PREC_PETSC_RCS}
};

/// @brief Map for PreconditionerType enum to a string name.
//
const std::map<PreconditionerType, std::string> preconditioner_type_to_name {
  {PreconditionerType::PREC_FSILS, "fsils"}, 
  {PreconditionerType::PREC_NONE, "none"}, 
  {PreconditionerType::PREC_RCS, "row-column-scaling"}, 
  {PreconditionerType::PREC_TRILINOS_DIAGONAL, "trilinos-diagonal"}, 
  {PreconditionerType::PREC_TRILINOS_BLOCK_JACOBI, "trilinos-blockjacobi"}, 
  {PreconditionerType::PREC_TRILINOS_ILU, "trilinos-ilu"}, 
  {PreconditionerType::PREC_TRILINOS_ILUT, "trilinos-ilut"}, 
  {PreconditionerType::PREC_TRILINOS_IC, "trilinos-ic"}, 
  {PreconditionerType::PREC_TRILINOS_IC, "trilinos-ict"}, 
  {PreconditionerType::PREC_TRILINOS_ML, "trilinos-ml"},
  {PreconditionerType::PREC_PETSC_JACOBI, "petsc-jacobi"},
  {PreconditionerType::PREC_PETSC_RCS, "petsc-rcs"}
};

/// @brief Map solver type string to SolverType enum. 
//
const std::map<std::string,SolverType> solver_name_to_type 
{
  {"bi-partitioned", SolverType::lSolver_NS},
  {"ns", SolverType::lSolver_NS},
  {"bpn", SolverType::lSolver_NS},
  {"bipn", SolverType::lSolver_NS},

  {"gmres", SolverType::lSolver_GMRES},

  {"conjugate-gradient", SolverType::lSolver_CG},
  {"cg", SolverType::lSolver_CG},

  {"bi-conjugate-gradient", SolverType::lSolver_BICGS},
  {"bicg", SolverType::lSolver_BICGS},
  {"bicgs", SolverType::lSolver_BICGS}
};


};

