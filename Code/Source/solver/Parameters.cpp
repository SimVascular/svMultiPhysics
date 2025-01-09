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

//
// The class methods defined here are used to process svFSIplus simulation parameters 
// read in from an XML-format file. 
//
// XML files are parsed using tinyxml2 (https://github.com/leethomason/tinyxml2).
//
//-----------------
// Section objects 
//-----------------
// Groups of related parameters making up the sections of an XML simulation 
// file are stored in section objects. For example
//
//   - GeneralSimulationParameters - General parameter section
//   - MeshParameters - Mesh parameter section
//   - EquationParameters - Equation parameter section
//   - ProjectionParameters  - Projection parameter section
//
// These section objects may also contain objects representing the sub-sections 
// defined for each section. 
//
// The name and default value for each parameter is defined in a section object's 
// constructor using the 'ParameterLists::set_parameter()' method. 
//
// Parameter values are set using the 'set_values()' method which contains calls to tinyxml2  
// to parse parameter values from an XML file. The XML elements within a section are 
// extracted in a while loop. Sub-sections or data will need checked and processed. 
// The 'ParameterLists::set_parameter_value()' method is used to set the value of a parameter 
// from a string. See MeshParameters::set_values() for an example.
//
// If a section does not contain any sub-sections then all parameters can be parsed automatically. 
// See LinearSolverParameters::set_values() for an example. 
//
#include "Parameters.h"
#include "consts.h"
#include "LinearAlgebra.h"
#include "ustruct.h"

#include <iostream>
#include <regex>
#include <set>
#include <sstream>
#include <limits>
#include <math.h>

/// @brief Set paramaters using a function pointing to the 'ParameterLists::set_parameter_value' method.
//
// Subsection names given in 'sub_sections' are ignored and processed elsewhere.
//
void xml_util_set_parameters( std::function<void(const std::string&, const std::string&)> fn, tinyxml2::XMLElement* xml_elem,
    const std::string& error_msg, std::set<std::string> sub_sections = std::set<std::string>())
{
  auto item = xml_elem->FirstChildElement();

  while (item != nullptr) {
    auto name = std::string(item->Value());

    if (sub_sections.count(name) == 0) {
      if (item->GetText() != nullptr) {
        auto value = item->GetText();
        try {
          fn(name, value);
        } catch (const std::bad_function_call& exception) {
          throw std::runtime_error(error_msg + name + "'.");
        }
      } else {
        throw std::runtime_error(error_msg + name + "'.");
      }
    }

    item = item->NextSiblingElement();
  }
}

//////////////////////////////////////////////////////////
//                   Parameters                         //
//////////////////////////////////////////////////////////

const std::string Parameters::FSI_FILE = "svMultiPhysicsFile";

const std::set<std::string> Parameters::constitutive_model_names = {
  "none",
  "nHK",
};

const std::set<std::string> Parameters::equation_names = {
  "none",
  "fluid",
  "struct",
};

Parameters::Parameters() 
{
}

void Parameters::get_logging_levels(int& verbose, int& warning, int& debug)
{ 
  /*
  verbose = general_simulation_parameters_.verbose;
  warning = general_simulation_parameters_.warning;
  debug = general_simulation_parameters_.debug;
  */
}

void Parameters::print_parameters()
{ 
  general_simulation_parameters.print_parameters();

  for (auto& mesh : mesh_parameters) { 
      mesh->print_parameters();
  }

  for (auto& equation : equation_parameters) { 
      equation->print_parameters();
  }
}

/// @brief Set the simulation parameter values given in an XML format file.
void Parameters::read_xml(std::string file_name)
{ 
  tinyxml2::XMLDocument doc;
  
  auto error = doc.LoadFile(file_name.c_str());
  
  auto root_element = doc.FirstChildElement(FSI_FILE.c_str());
  if (root_element == nullptr) {
    throw std::runtime_error("The following error occured while reading the XML file '" + file_name + "'.\n" + 
        "[svFSI] ERROR " + std::string(doc.ErrorStr()));
  }

  // Get general parameters.
  general_simulation_parameters.set_values(root_element);

  // Set Contact values.
  set_contact_values(root_element);

  // Set Add_mesh values.
  set_mesh_values(root_element);

  // Set Precomputed_solution values.
  set_precomputed_solution_values(root_element);

  // Set mesh projection parameters.
  set_projection_values(root_element);

  // Set Add_equation values.
  set_equation_values(root_element);
}

void Parameters::set_contact_values(tinyxml2::XMLElement* root_element)
{
  auto item = root_element->FirstChildElement(ContactParameters::xml_element_name_.c_str());
  
  if (item == nullptr) {
    return;
  }

  contact_parameters.set_values(item);
}

void Parameters::set_equation_values(tinyxml2::XMLElement* root_element)
{
  auto add_eq_item = root_element->FirstChildElement(EquationParameters::xml_element_name_.c_str());

  while (add_eq_item) {
    const char* eq_type;
    auto result = add_eq_item->QueryStringAttribute("type", &eq_type);

    auto eq_params = new EquationParameters();
    eq_params->type.set(std::string(eq_type));
    eq_params->set_values(add_eq_item);
    equation_parameters.push_back(eq_params);

    add_eq_item = add_eq_item->NextSiblingElement(EquationParameters::xml_element_name_.c_str());
  }
}

void Parameters::set_mesh_values(tinyxml2::XMLElement* root_element)
{
  auto add_mesh_item = root_element->FirstChildElement(MeshParameters::xml_element_name_.c_str());

  while (add_mesh_item) {
    const char* mesh_name;
    auto result = add_mesh_item->QueryStringAttribute("name", &mesh_name);

    MeshParameters* mesh_params = new MeshParameters();
    mesh_params->name.set(std::string(mesh_name));
    mesh_params->set_values(add_mesh_item);
    mesh_parameters.push_back(mesh_params);

    add_mesh_item = add_mesh_item->NextSiblingElement(MeshParameters::xml_element_name_.c_str());
  }
}

void Parameters::set_precomputed_solution_values(tinyxml2::XMLElement* root_element)
{
  auto add_pre_sol_item = root_element->FirstChildElement(PrecomputedSolutionParameters::xml_element_name_.c_str());
  if (add_pre_sol_item == nullptr) { 
    return;
  }

  precomputed_solution_parameters.set_values(add_pre_sol_item);
}

void Parameters::set_projection_values(tinyxml2::XMLElement* root_element)
{
  auto add_proj_item = root_element->FirstChildElement(ProjectionParameters::xml_element_name_.c_str());

  while (add_proj_item) {
    const char* proj_name;
    auto result = add_proj_item->QueryStringAttribute("name", &proj_name);

    ProjectionParameters* proj_params = new ProjectionParameters();
    proj_params->name.set(std::string(proj_name));
    proj_params->set_values(add_proj_item);
    projection_parameters.push_back(proj_params);

    add_proj_item = add_proj_item->NextSiblingElement(ProjectionParameters::xml_element_name_.c_str());
  }
}

//////////////////////////////////////////////////////////
//                BodyForceParameters                   //
//////////////////////////////////////////////////////////

// Body force over a mesh using the "Add_BF" command.

/// @brief Define the XML element name for boundary condition parameters.
const std::string BodyForceParameters::xml_element_name_ = "Add_BF";

BodyForceParameters::BodyForceParameters()
{
  // A parameter that must be defined.
  bool required = true;

  mesh_name = Parameter<std::string>("mesh", "", required);

  set_parameter("Fourier_coefficients_file_path", "", !required, fourier_coefficients_file_path);
  // [TODO:DaveP] I'm not sure if this is required in the Add_BF element.
  //set_parameter("Ramp_function", false, !required, ramp_function);
  set_parameter("Spatial_values_file_path", "", !required, spatial_values_file_path);

  set_parameter("Temporal_and_spatial_values_file_path", "", !required, temporal_and_spatial_values_file_path);
  set_parameter("Temporal_values_file_path", "", !required, temporal_values_file_path);
  set_parameter("Time_dependence", "Steady", !required, time_dependence);
  set_parameter("Type", "", required, type);

  set_parameter("Value", 0.0, !required, value);
}

void BodyForceParameters::print_parameters()
{
  std::cout << std::endl;
  std::cout << "---------------------" << std::endl;
  std::cout << "Body Force Parameters" << std::endl;
  std::cout << "---------------------" << std::endl;
  std::cout << mesh_name.name() << ": " << mesh_name.value() << std::endl;

  auto params_name_value = get_parameter_list();
  for (auto& [ key, value ] : params_name_value) { 
    std::cout << key << ": " << value << std::endl;
  }
}

void BodyForceParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  using namespace tinyxml2;
  std::string error_msg = "Unknown " + xml_element_name_ + " XML element '";

  // Get the 'type' from the <LS type=TYPE> element.
  const char* smesh;
  auto result = xml_elem->QueryStringAttribute("mesh", &smesh);
  if (smesh == nullptr) {
    throw std::runtime_error("No MESH given in the XML <Add_BF mesh=MESH> element.");
  }
  mesh_name.set(std::string(smesh));
  //auto item = xml_elem->FirstChildElement();

  using std::placeholders::_1;
  using std::placeholders::_2;

  std::function<void(const std::string&, const std::string&)> ftpr =
      std::bind( &BodyForceParameters::set_parameter_value, *this, _1, _2);

  xml_util_set_parameters(ftpr, xml_elem, error_msg);
}

//////////////////////////////////////////////////////////
//            BoundaryConditionParameters               //
//////////////////////////////////////////////////////////

// The BoundaryConditionParameters stores paramaters for various
// type of boundary conditions under the Add_BC XML element.

/// @brief Define the XML element name for equation boundary condition parameters.
const std::string BoundaryConditionParameters::xml_element_name_ = "Add_BC";
const std::string BoundaryConditionRCRParameters::xml_element_name_ = "RCR_values";

/// @brief RCR values for Neumann BC type.
BoundaryConditionRCRParameters::BoundaryConditionRCRParameters()
{
  // A parameter that must be defined.
  bool required = true;

  set_parameter("Capacitance", 0.0, required, capacitance);

  set_parameter("Distal_resistance", 0.0, required, distal_resistance);
  set_parameter("Distal_pressure", 0.0, !required, distal_pressure);

  set_parameter("Initial_pressure", 0.0, !required, initial_pressure);

  set_parameter("Proximal_resistance", 0.0, required, proximal_resistance);
}

void BoundaryConditionRCRParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  std::string error_msg = "Unknown " + xml_element_name_ + " XML element '"; 
  using std::placeholders::_1;
  using std::placeholders::_2;
  std::function<void(const std::string&, const std::string&)> ftpr =
      std::bind( &BoundaryConditionRCRParameters::set_parameter_value, *this, _1, _2);
  xml_util_set_parameters(ftpr, xml_elem, error_msg);

  value_set = true;
}

void BoundaryConditionRCRParameters::print_parameters()
{
  std::cout << std::endl; 
  std::cout << "---------------------------------" << std::endl;
  std::cout << "Boundary Condition RCR Parameters" << std::endl;
  std::cout << "---------------------------------" << std::endl;

  auto params_name_value = get_parameter_list();
  for (auto& [ key, value ] : params_name_value) {
    std::cout << key << ": " << value << std::endl;
  }
}

BoundaryConditionParameters::BoundaryConditionParameters()
{
  // A parameter that must be defined.
  bool required = true;

  // Set name from  Add_BC name="" XML element.
  name = Parameter<std::string>("name", "", required);

  set_parameter("Apply_along_normal_direction", false, !required, apply_along_normal_direction);
  set_parameter("Bct_file_path", "", !required, bct_file_path);

  set_parameter("Damping", 1.0, !required, damping);
  set_parameter("Distal_pressure", {}, !required, distal_pressure);

  set_parameter("Effective_direction", {}, !required, effective_direction);
  //set_parameter("Effective_direction", {0,0,0}, !required, effective_direction);
  set_parameter("Follower_pressure_load", false, !required, follower_pressure_load);
  set_parameter("Fourier_coefficients_file_path", "", !required, fourier_coefficients_file_path);

  set_parameter("Impose_flux", false, !required, impose_flux);
  set_parameter("Impose_on_state_variable_integral", false, !required, impose_on_state_variable_integral);
  set_parameter("Initial_displacements_file_path", "", !required, initial_displacements_file_path);

  set_parameter("Penalty_parameter", 0.0, !required, penalty_parameter);
  set_parameter("Penalty_parameter_normal", 0.0, !required, penalty_parameter_normal);
  set_parameter("Penalty_parameter_tangential", 0.0, !required, penalty_parameter_tangential);
  set_parameter("Prestress_file_path", "", !required, prestress_file_path);
  set_parameter("Profile", "Flat", !required, profile);

  set_parameter("Ramp_function", false, !required, ramp_function);

  set_parameter("CST_shell_bc_type", "", !required, cst_shell_bc_type);
  set_parameter("Spatial_profile_file_path", "", !required, spatial_profile_file_path);
  set_parameter("Spatial_values_file_path", "", !required, spatial_values_file_path);
  set_parameter("Stiffness", 1.0, !required, stiffness);

  set_parameter("Temporal_and_spatial_values_file_path", "", !required, temporal_and_spatial_values_file_path);
  set_parameter("Temporal_values_file_path", "", !required, temporal_values_file_path);
  set_parameter("Time_dependence", "Steady", !required, time_dependence);
  set_parameter("Traction_values_file_path", "", !required, traction_values_file_path);
  set_parameter("Traction_multiplier", 1.0, !required, traction_multiplier);
  set_parameter("Type", "", required, type);

  set_parameter("Undeforming_neu_face", false, !required, undeforming_neu_face);
  set_parameter("Value", 0.0, !required, value);

  set_parameter("Weakly_applied", false, !required, weakly_applied);
  set_parameter("Zero_out_perimeter", false, !required, zero_out_perimeter);
}

void BoundaryConditionParameters::print_parameters()
{ 
  std::cout << std::endl; 
  std::cout << "-----------------------------" << std::endl;
  std::cout << "Boundary Condition Parameters" << std::endl;
  std::cout << "-----------------------------" << std::endl;
  std::cout << name.name() << ": " << name.value() << std::endl;

  auto params_name_value = get_parameter_list();
  for (auto& [ key, value ] : params_name_value) {
    std::cout << key << ": " << value << std::endl;
  }

  rcr.print_parameters();
}

void BoundaryConditionParameters::set_values(tinyxml2::XMLElement* xml_elem)
{ 
  using namespace tinyxml2;
  std::string error_msg = "Unknown " + xml_element_name_ + " XML element '"; 

  // Get the 'name' from the <Add_BC name=NAME> element.
  const char* sname;
  auto result = xml_elem->QueryStringAttribute("name", &sname);
  if (sname == nullptr) {
    throw std::runtime_error("No NAME given in the XML <Add_BC name=NAME> element.");
  }
  name.set(std::string(sname));

  auto item = xml_elem->FirstChildElement();

  while (item != nullptr) {
    auto name = std::string(item->Value());

    if (name == BoundaryConditionRCRParameters::xml_element_name_) {
      rcr.set_values(item);
    }
   
    else if (item->GetText() != nullptr) {
      auto value = item->GetText();
      try {
        set_parameter_value(name, value);
      } catch (const std::bad_function_call& exception) {
        throw std::runtime_error(error_msg + name + "'.");
      }
    } else {
      throw std::runtime_error(error_msg + name + "'.");
    }

    item = item->NextSiblingElement();
  }
}

//////////////////////////////////////////////////////////
//            ConstitutiveModelParameters               //
//////////////////////////////////////////////////////////

/// @brief Process parameters for various constitutive models.
///
/// Define the XML element name for constitutive parameters.
const std::string ConstitutiveModelParameters::xml_element_name_ = "Constitutive_model";

// [TODO] Should use the types defined in consts.h.
const std::string ConstitutiveModelParameters::GUCCIONE_MODEL = "Guccione";
const std::string ConstitutiveModelParameters::HGO_MODEL = "HGO";
const std::string ConstitutiveModelParameters::HOLZAPFEL_OGDEN_MODEL = "HolzapfelOgden";
const std::string ConstitutiveModelParameters::HOLZAPFEL_OGDEN_MA_MODEL = "HolzapfelOgden-ModifiedAnisotropy";
const std::string ConstitutiveModelParameters::LEE_SACKS = "Lee-Sacks";
const std::string ConstitutiveModelParameters::NEOHOOKEAN_MODEL = "neoHookean";
const std::string ConstitutiveModelParameters::STVENANT_KIRCHHOFF_MODEL = "stVenantKirchhoff";

/// @brief Supported constitutive model types and their aliases.
const std::map<std::string, std::string> ConstitutiveModelParameters::constitutive_model_types = {
  { ConstitutiveModelParameters::GUCCIONE_MODEL, ConstitutiveModelParameters::GUCCIONE_MODEL},
  { "Gucci",                                     ConstitutiveModelParameters::GUCCIONE_MODEL},

  {ConstitutiveModelParameters::HGO_MODEL, ConstitutiveModelParameters::HGO_MODEL},

  {ConstitutiveModelParameters::HOLZAPFEL_OGDEN_MODEL, ConstitutiveModelParameters::HOLZAPFEL_OGDEN_MODEL},

  {ConstitutiveModelParameters::HOLZAPFEL_OGDEN_MA_MODEL, ConstitutiveModelParameters::HOLZAPFEL_OGDEN_MA_MODEL},

  {ConstitutiveModelParameters::LEE_SACKS, ConstitutiveModelParameters::LEE_SACKS},

  {ConstitutiveModelParameters::NEOHOOKEAN_MODEL, ConstitutiveModelParameters::NEOHOOKEAN_MODEL},
  {"nHK", ConstitutiveModelParameters::NEOHOOKEAN_MODEL},

  {ConstitutiveModelParameters::STVENANT_KIRCHHOFF_MODEL, ConstitutiveModelParameters::STVENANT_KIRCHHOFF_MODEL},
  {"stVK",                                                ConstitutiveModelParameters::STVENANT_KIRCHHOFF_MODEL},
}; 

/// @brief Define a map to set the parameters for each constitutive model.
using CmpType = ConstitutiveModelParameters*;
using CmpXmlType = tinyxml2::XMLElement*;
using SetConstitutiveModelParamMapType = std::map<std::string, std::function<void(CmpType, CmpXmlType)>>;

SetConstitutiveModelParamMapType SetConstitutiveModelParamMap = {
  {ConstitutiveModelParameters::GUCCIONE_MODEL, [](CmpType cp, CmpXmlType params) -> void {cp->guccione.set_values(params);}},
  {ConstitutiveModelParameters::HGO_MODEL, [](CmpType cp, CmpXmlType params) -> void {cp->holzapfel_gasser_ogden.set_values(params);}},
  {ConstitutiveModelParameters::HOLZAPFEL_OGDEN_MODEL, [](CmpType cp, CmpXmlType params) -> void {cp->holzapfel.set_values(params);}},
  {ConstitutiveModelParameters::HOLZAPFEL_OGDEN_MA_MODEL, [](CmpType cp, CmpXmlType params) -> void {cp->holzapfel.set_values(params);}},
  {ConstitutiveModelParameters::LEE_SACKS, [](CmpType cp, CmpXmlType params) -> void {cp->lee_sacks.set_values(params);}},
  {ConstitutiveModelParameters::NEOHOOKEAN_MODEL, [](CmpType cp, CmpXmlType params) -> void {cp->neo_hookean.set_values(params);}},
  {ConstitutiveModelParameters::STVENANT_KIRCHHOFF_MODEL, [](CmpType cp, CmpXmlType params) -> void {cp->stvenant_kirchhoff.set_values(params);}},
};

/// @brief Define a map to print parameters for each constitutive model.
using PrintConstitutiveModelParamMapType = std::map<std::string, std::function<void(CmpType)>>;

PrintConstitutiveModelParamMapType PrintConstitutiveModelParamMap = {
  {ConstitutiveModelParameters::GUCCIONE_MODEL, [](CmpType cp) -> void {cp->guccione.print_parameters();}},
  {ConstitutiveModelParameters::HGO_MODEL, [](CmpType cp) -> void {cp->holzapfel_gasser_ogden.print_parameters();}},
  {ConstitutiveModelParameters::HOLZAPFEL_OGDEN_MODEL, [](CmpType cp) -> void {cp->holzapfel.print_parameters();}},
  {ConstitutiveModelParameters::HOLZAPFEL_OGDEN_MA_MODEL, [](CmpType cp) -> void {cp->holzapfel.print_parameters();}},
  {ConstitutiveModelParameters::LEE_SACKS, [](CmpType cp) -> void {cp->lee_sacks.print_parameters();}},
  {ConstitutiveModelParameters::NEOHOOKEAN_MODEL, [](CmpType cp) -> void {cp->neo_hookean.print_parameters();}},
  {ConstitutiveModelParameters::STVENANT_KIRCHHOFF_MODEL, [](CmpType cp) -> void {cp->stvenant_kirchhoff.print_parameters();}},
};


GuccioneParameters::GuccioneParameters()
{
  // A parameter that must be defined.
  bool required = true;

  set_parameter("bf", 0.0, required, bf);
  set_parameter("bfs", 0.0, required, bfs);
  set_parameter("bt", 0.0, required, bt);
  set_parameter("c", 0.0, required, c);

  set_xml_element_name("Constitutive_model type=Guccione");
}

void GuccioneParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  std::string error_msg = "Unknown Constitutive_model type=Guccione XML element '"; 

  using std::placeholders::_1;
  using std::placeholders::_2;
  std::function<void(const std::string&, const std::string&)> ftpr =
      std::bind( &GuccioneParameters::set_parameter_value, *this, _1, _2);

  xml_util_set_parameters(ftpr, xml_elem, error_msg);

  value_set = true;
}

void GuccioneParameters::print_parameters()
{
  std::cout << "Guccione: " << std::endl;
  auto params_name_value = get_parameter_list(); 
  for (auto& [ key, value ] : params_name_value) {
    std::cout << key << ": " << value << std::endl;
  }
}

HolzapfelParameters::HolzapfelParameters()
{ 
  // A parameter that must be defined.
  bool required = true;

  set_parameter("a", 0.0, required, a);
  set_parameter("b", 0.0, required, b);

  set_parameter("a4f", 0.0, required, a4f);
  set_parameter("b4f", 0.0, required, b4f);

  set_parameter("a4s", 0.0, required, a4s);
  set_parameter("b4s", 0.0, required, b4s);

  set_parameter("afs", 0.0, required, afs);
  set_parameter("bfs", 0.0, required, bfs);

  set_parameter("k", 0.0, required, k);
  
  set_xml_element_name("Constitutive_model type=Holzapfel");
}

void HolzapfelParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
 std::string error_msg = "Unknown Constitutive_model type=Holzapfel XML element '";
  
  using std::placeholders::_1;
  using std::placeholders::_2;
  std::function<void(const std::string&, const std::string&)> ftpr =
      std::bind( &HolzapfelParameters::set_parameter_value, *this, _1, _2);
  
  xml_util_set_parameters(ftpr, xml_elem, error_msg);
  
  value_set = true;
}

void HolzapfelParameters::print_parameters()
{
  auto params_name_value = get_parameter_list(); 
  for (auto& [ key, value ] : params_name_value) {
    std::cout << key << ": " << value << std::endl;
  }
}

HolzapfelGasserOgdenParameters::HolzapfelGasserOgdenParameters()
{ 
  // A parameter that must be defined.
  bool required = true;

  set_parameter("a4", 0.0, required, a4);
  set_parameter("b4", 0.0, required, b4);
  set_parameter("a6", 0.0, required, a6);
  set_parameter("b6", 0.0, required, b6);
  set_parameter("kappa", 0.0, required, kappa);

  set_xml_element_name("Constitutive_model type=HGO");
}

void HolzapfelGasserOgdenParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  std::string error_msg = "Unknown Constitutive_model type=HGO XML element '"; 

  using std::placeholders::_1;
  using std::placeholders::_2;
  std::function<void(const std::string&, const std::string&)> ftpr = 
      std::bind( &HolzapfelGasserOgdenParameters::set_parameter_value, *this, _1, _2);

  xml_util_set_parameters(ftpr, xml_elem, error_msg);
  value_set = true;
}

void HolzapfelGasserOgdenParameters::print_parameters()
{
  auto params_name_value = get_parameter_list();
  for (auto& [ key, value ] : params_name_value) {
    std::cout << key << ": " << value << std::endl;
  }
}

LeeSacksParameters::LeeSacksParameters()
{
  // A parameter that must be defined.
  bool required = true;

  set_parameter("a", 0.0, required, a);
  set_parameter("a0", 0.0, required, a);
  set_parameter("b1", 0.0, required, b1);
  set_parameter("b2", 0.0, required, b2);
  set_parameter("mu0", 0.0, required, mu0);

  set_xml_element_name("Constitutive_model type=Lee-Sacks");
}

void LeeSacksParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  std::string error_msg = "Unknown Constitutive_model type=Lee-Sacks XML element '";

  using std::placeholders::_1;
  using std::placeholders::_2;
  std::function<void(const std::string&, const std::string&)> ftpr =
      std::bind( &LeeSacksParameters::set_parameter_value, *this, _1, _2);

  xml_util_set_parameters(ftpr, xml_elem, error_msg);

  value_set = true;
}

void LeeSacksParameters::print_parameters()
{
  std::cout << "Lee-Sacks: " << std::endl;
  auto params_name_value = get_parameter_list();
  for (auto& [ key, value ] : params_name_value) {
    std::cout << key << ": " << value << std::endl;
  }
}

MooneyRivlinParameters::MooneyRivlinParameters()
{
  // A parameter that must be defined.
  bool required = true;

  set_parameter("c1", 0.0, required, c1);
  set_parameter("c2", 0.0, required, c2);

  set_xml_element_name("Constitutive_model type=Mooney-Rivlin");
}

void MooneyRivlinParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  std::string error_msg = "Unknown Constitutive_model type=Mooney-Rivlin  XML element '";

  using std::placeholders::_1;
  using std::placeholders::_2;
  std::function<void(const std::string&, const std::string&)> ftpr =
      std::bind( &MooneyRivlinParameters::set_parameter_value, *this, _1, _2);

  xml_util_set_parameters(ftpr, xml_elem, error_msg);

  value_set = true;
}

void MooneyRivlinParameters::print_parameters()
{
  std::cout << "MooneyRivlin: " << std::endl;
  std::cout << c1.name() << ": " << c1.value() << std::endl;;
  std::cout << c2.name() << ": " << c2.value() << std::endl;;
}

/// @brief There are no parameters associated with a Neohookean model.
NeoHookeanParameters::NeoHookeanParameters()
{
}

void NeoHookeanParameters::set_values(tinyxml2::XMLElement* con_params)
{
  value_set = true;
}

void NeoHookeanParameters::print_parameters()
{
}

/// @brief There are no parameters associated with a StVenantKirchhoff model.
StVenantKirchhoffParameters::StVenantKirchhoffParameters()
{
  value_set = true;
}

void StVenantKirchhoffParameters::set_values(tinyxml2::XMLElement* con_params)
{
  value_set = true;
}

void StVenantKirchhoffParameters::print_parameters()
{
}

ConstitutiveModelParameters::ConstitutiveModelParameters()
{
  // A parameter that must be defined.
  bool required = true;

  // Type from <Constitutive_model type=TYPE > XML element.
  type = Parameter<std::string>("type", "", required);
}

void ConstitutiveModelParameters::print_parameters()
{
  if (!value_set) {
    return;
  }
  std::cout << std::endl;
  std::cout << "-----------------------------" << std::endl;
  std::cout << "Constitutive Model Parameters" << std::endl;
  std::cout << "-----------------------------" << std::endl;
  std::cout << type.name() << ": '" << type.value() << "'" << std::endl;

  PrintConstitutiveModelParamMap[type.value()](this);
}

void ConstitutiveModelParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  using namespace tinyxml2;
  std::string error_msg = "Unknown " + xml_element_name_ + " XML element '"; 

  // Get the 'type' from the <Constitutive_model type= > element.
  const char* stype;
  auto result = xml_elem->QueryStringAttribute("type", &stype);
  if (stype == nullptr) {
    throw std::runtime_error("No TYPE given in the XML <Constitutive_model type=TYPE> element.");
  }
  type.set(std::string(stype));

  // Check constitutive model type.
  if (constitutive_model_types.count(type.value()) == 0) {
    auto msg_1 = "Unknown constitutive model '" + type.value() + " in the '" + xml_elem->Name() + "' XML element.\n";
    std::string msg_2 = "Valid types are:";
    for (auto& name : constitutive_model_types) {
      msg_2 += " " + name.first;
    }
    msg_2 += "\n";
    auto msg = msg_1 + msg_2;
    throw std::runtime_error(msg);
  }
  auto model_type = constitutive_model_types.at(type.value());
  type.set(model_type);

  // Set parameters for the given constitutive model.
  SetConstitutiveModelParamMap[model_type](this, xml_elem);

  value_set = true;
}

/// @brief Check if a constitutive model is valid for the given equation.
//
void ConstitutiveModelParameters::check_constitutive_model(const Parameter<std::string>& eq_type_str)
{
  auto eq_type = consts::equation_name_to_type.at(eq_type_str.value());
  auto model = consts::constitutive_model_name_to_type.at(type.value());

  if (eq_type == consts::EquationType::phys_ustruct) {
    if (! ustruct::constitutive_model_is_valid(model)) {
      throw std::runtime_error("The " + type.value() + " constitutive model is not valid for ustruct equations.");
    }
  }
}

//////////////////////////////////////////////////////////
//                  CoupleCplBCParameters               //
//////////////////////////////////////////////////////////

/// @brief Couple to reduced-order models.
///
/// Define the XML element name for equation Couple_to_genBC parameters.
const std::string CoupleCplBCParameters::xml_element_name_ = "Couple_to_cplBC";

CoupleCplBCParameters::CoupleCplBCParameters()
{
  // A parameter that must be defined.
  bool required = true;

  // Define attributes.
  type = Parameter<std::string>("type", "", required);

  set_parameter("File_name_for_0D_3D_communication", "", required, file_name_for_0D_3D_communication);
  set_parameter("File_name_for_saving_unknowns", "", required, file_name_for_saving_unknowns);
  set_parameter("Number_of_unknowns", 0, required, number_of_unknowns);
  set_parameter("Number_of_user_defined_outputs", 0, required, number_of_user_defined_outputs);
  set_parameter("Unknowns_initialization_file_path", "", !required, unknowns_initialization_file_path);
  set_parameter("ZeroD_code_file_path", "", required, zerod_code_file_path);
}

void CoupleCplBCParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  std::string error_msg = "Unknown Couple_to_cplBC type=TYPE XML element '";

  // Get the 'type' from the <Couple_to_cplBC type=TYPE> element.
  const char* stype;
  auto result = xml_elem->QueryStringAttribute("type", &stype);
  if (stype == nullptr) {
    throw std::runtime_error("No TYPE given in the XML <Stimulus=TYPE> element.");
  }
  type.set(std::string(stype));
  auto item = xml_elem->FirstChildElement();

  using std::placeholders::_1;
  using std::placeholders::_2;
  std::function<void(const std::string&, const std::string&)> ftpr =
      std::bind( &CoupleCplBCParameters::set_parameter_value, *this, _1, _2);

  xml_util_set_parameters(ftpr, xml_elem, error_msg);

  value_set = true;
}

void CoupleCplBCParameters::print_parameters()
{
  if (!value_set) { 
    return;
  }
  std::cout << std::endl;
  std::cout << "----------------------" << std::endl;
  std::cout << "CoupleCplBC Parameters" << std::endl;
  std::cout << "----------------------" << std::endl;
  std::cout << type.name() << ": " << type.value() << std::endl;

  auto params_name_value = get_parameter_list();
  for (auto& [ key, value ] : params_name_value) {
    std::cout << key << ": " << value << std::endl;
  }
}

//////////////////////////////////////////////////////////
//                  CoupleGenBCParameters               //
//////////////////////////////////////////////////////////

/// @brief Coupling to GenBC.
///
/// Define the XML element name for equation Couple_to_genBC parameters.
const std::string CoupleGenBCParameters::xml_element_name_ = "Couple_to_genBC";

CoupleGenBCParameters::CoupleGenBCParameters()
{
  // A parameter that must be defined.
  bool required = true;

  type = Parameter<std::string>("type", "", required);

  set_parameter("ZeroD_code_file_path", "", required, zerod_code_file_path);
};

void CoupleGenBCParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  std::string error_msg = "Unknown Couple_to_genBC type=TYPE XML element '";
  
  // Get the 'type' from the <Couple_to_genBC type=TYPE> element.
  const char* stype;
  auto result = xml_elem->QueryStringAttribute("type", &stype);
  if (stype == nullptr) {
    throw std::runtime_error("No TYPE given in the XML <Couple_to_genBC type=TYPE> element.");
  }
  type.set(std::string(stype));
  auto item = xml_elem->FirstChildElement();
  
  using std::placeholders::_1;
  using std::placeholders::_2;
  std::function<void(const std::string&, const std::string&)> ftpr = 
      std::bind( &CoupleGenBCParameters::set_parameter_value, *this, _1, _2);
  
  xml_util_set_parameters(ftpr, xml_elem, error_msg);
  
  value_set = true;
}


//////////////////////////////////////////////////////////
//                  CoupleSvZeroDParameters               //
//////////////////////////////////////////////////////////

// Coupling to svZeroD.

// Define the XML element name for equation Couple_to_svZeroD parameters.
const std::string CoupleSvZeroDParameters::xml_element_name_ = "Couple_to_svZeroD";

CoupleSvZeroDParameters::CoupleSvZeroDParameters()
{
  // A parameter that must be defined.
  bool required = true;

  type = Parameter<std::string>("type", "", required);
};

void CoupleSvZeroDParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  std::string error_msg = "Unknown Couple_to_svZeroD type=TYPE XML element '";
  
  // Get the 'type' from the <Couple_to_genBC type=TYPE> element.
  const char* stype;
  auto result = xml_elem->QueryStringAttribute("type", &stype);
  if (stype == nullptr) {
    throw std::runtime_error("No TYPE given in the XML <Couple_to_svZeroD type=TYPE> element.");
  }
  type.set(std::string(stype));
  auto item = xml_elem->FirstChildElement();
  
  using std::placeholders::_1;
  using std::placeholders::_2;
  std::function<void(const std::string&, const std::string&)> ftpr = 
      std::bind( &CoupleSvZeroDParameters::set_parameter_value, *this, _1, _2);
  
  xml_util_set_parameters(ftpr, xml_elem, error_msg);
  
  value_set = true;
}


//////////////////////////////////////////////////////////
//                  OutputParameters                    //
//////////////////////////////////////////////////////////

// The OutputParameters class stores parameters for the
// Output XML sub-element under Add_equation element.

/// @brief Define the XML element name for equation output parameters.
const std::string OutputParameters::xml_element_name_ = "Output";

OutputParameters::OutputParameters()
{
  // A parameter that must be defined.
  bool required = true;

  // Type for the <Output type="Spatial" > xml element. 
  type = Parameter<std::string>("type", "", required);
}

void OutputParameters::print_parameters()
{ 
  std::cout << std::endl;
  std::cout << "-----------------" << std::endl;
  std::cout << "Output Parameters" << std::endl;
  std::cout << "-----------------" << std::endl;
  std::cout << type.name() << ": " << type.value() << std::endl;
}

void OutputParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  std::string msg("[OutputParameters::set_values] ");
  std::string error_msg = "Unknown " + xml_element_name_ + " XML element '"; 

  const char* stype;
  auto result = xml_elem->QueryStringAttribute("type", &stype);
  type.set(std::string(stype));

  // Get values from XML file.
  //
  // The 'Alias' type needs special processing, storing
  // alias output names in a map.
  //
  if (type.value() == "Alias") {
    auto item = xml_elem->FirstChildElement();
    while (item != nullptr) {
      auto name = std::string(item->Name());
      auto value = std::string(item->GetText());
      Parameter<std::string> param(name, "", false);
      param.set(value);
      alias_list.emplace_back(param);
      item = item->NextSiblingElement();
    }
  } else {
    auto item = xml_elem->FirstChildElement();
    while (item != nullptr) {
      auto name = std::string(item->Name());
      auto value = std::string(item->GetText());
      Parameter<bool> param(name, false, false);
      param.set(value);
      output_list.emplace_back(param);
      item = item->NextSiblingElement();
    }
  }
}

/// @brief Get the value of an alias by name.
std::string OutputParameters::get_alias_value(const std::string& name)
{
  for (auto& param : alias_list) {
    if (param.name_ == name) {
      return param.value();
    }
  }

  return "";
}

/// @brief Get the value of an output by name.
bool OutputParameters::get_output_value(const std::string& name)
{
  for (auto& param : output_list) {
    if (param.name_ == name) {
      return param.value();
    }
  }

  return false;
}

//////////////////////////////////////////////////////////
//               VariableWallPropsParameters            //
//////////////////////////////////////////////////////////

/// @brief The VariableWallPropsParameters class stores parameters for
/// variable wall properties for the CMM equation.
///
/// Define the XML element name for variable wall parameters.
const std::string VariableWallPropsParameters::xml_element_name_ = "Variable_wall_properties";

VariableWallPropsParameters::VariableWallPropsParameters()
{
  // A parameter that must be defined.
  bool required = true;

  mesh_name = Parameter<std::string>("mesh_name", "", required);

  set_parameter("Wall_properties_file_path", "", required, wall_properties_file_path);
}

void VariableWallPropsParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  using namespace tinyxml2;
  std::string error_msg = "Unknown " + xml_element_name_ + " XML element '";

  // Get the 'type' from the <Variable_wall_properties mesh_name=NAME> element.
  const char* sname;
  auto result = xml_elem->QueryStringAttribute("mesh_name", &sname);
  if (sname == nullptr) {
    throw std::runtime_error("No TYPE given in the XML <Variable_wall_properties mesh_name=NAME> element.");
  }
  mesh_name.set(std::string(sname));
  auto item = xml_elem->FirstChildElement();

  using std::placeholders::_1;
  using std::placeholders::_2;

  std::function<void(const std::string&, const std::string&)> ftpr =
      std::bind( &VariableWallPropsParameters::set_parameter_value, *this, _1, _2);

  xml_util_set_parameters(ftpr, xml_elem, error_msg);

  value_set = true;
}

//////////////////////////////////////////////////////////
//                  FluidViscosityParameters                 //
//////////////////////////////////////////////////////////

/// @brief Process parameters for various fluid viscosity models.
///
/// Define the XML element name for viscosiity parameters.
const std::string FluidViscosityParameters::xml_element_name_ = "Viscosity";

const std::string FluidViscosityParameters::CONSTANT_MODEL = "Constant";
const std::string FluidViscosityParameters::CARREAU_YASUDA_MODEL = "Carreau-Yasuda";
const std::string FluidViscosityParameters::CASSONS_MODEL = "Cassons";

const std::set<std::string> FluidViscosityParameters::model_names = {
  FluidViscosityParameters::CONSTANT_MODEL,
  FluidViscosityParameters::CARREAU_YASUDA_MODEL,
  FluidViscosityParameters::CASSONS_MODEL
};

/// @brief Define a map to set parameters for each fluid viscosity model.
using FVpType = FluidViscosityParameters*;
using XmlType = tinyxml2::XMLElement*;
using SetFluidViscosityParamMapType = std::map<std::string, std::function<void(FVpType, XmlType)>>;
SetFluidViscosityParamMapType SetFluidViscosityModelParamsMap = {
  {FluidViscosityParameters::CARREAU_YASUDA_MODEL, [](FVpType vp, XmlType params) -> void { vp->carreau_yasuda_model.set_values(params); }},
  {FluidViscosityParameters::CASSONS_MODEL, [](FVpType vp, XmlType params) -> void { vp->cassons_model.set_values(params); }},
  {FluidViscosityParameters::CONSTANT_MODEL, [](FVpType vp, XmlType params) -> void { vp->newtonian_model.set_values(params); }},
};

/// @brief Define a map to print parameters for each fluid viscosity model.
using PrintFluidViscosityParamaMapType = std::map<std::string, std::function<void(FVpType)>>;
PrintFluidViscosityParamaMapType PrintFluidViscosityModelParamsMap = {
  {FluidViscosityParameters::CARREAU_YASUDA_MODEL, [](FVpType vp) -> void { vp->carreau_yasuda_model.print_parameters(); }},
  {FluidViscosityParameters::CASSONS_MODEL, [](FVpType vp) -> void { vp->cassons_model.print_parameters(); }},
  {FluidViscosityParameters::CONSTANT_MODEL, [](FVpType vp) -> void { vp->newtonian_model.print_parameters(); }},
};

FluidViscosityNewtonianParameters::FluidViscosityNewtonianParameters()
{
  // A parameter that must be defined.
  bool required = true;

  set_parameter("Value", 0.0, !required, constant_value);
}

void FluidViscosityNewtonianParameters::print_parameters()
{
  std::cout << constant_value.name_ << ": " << constant_value.value_ << std::endl;
}

void FluidViscosityNewtonianParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  std::string error_msg = "Unknown Constitutive_model type=Newtonian XML element '";

  using std::placeholders::_1;
  using std::placeholders::_2;
  std::function<void(const std::string&, const std::string&)> ftpr =
      std::bind( &FluidViscosityNewtonianParameters::set_parameter_value, *this, _1, _2);

  xml_util_set_parameters(ftpr, xml_elem, error_msg);
}

FluidViscosityCarreauYasudaParameters::FluidViscosityCarreauYasudaParameters()
{
  // A parameter that must be defined.
  bool required = true;

  set_parameter("Limiting_high_shear_rate_viscosity", 0.0, !required, limiting_high_shear_rate_viscosity);
  set_parameter("Limiting_low_shear_rate_viscosity", 0.0, !required, limiting_low_shear_rate_viscosity);

  set_parameter("Power_law_index", 0.0, !required, power_law_index);

  set_parameter("Shear_rate_tensor_multiplier", 0.0, !required, shear_rate_tensor_multipler);
  set_parameter("Shear_rate_tensor_exponent", 0.0, !required, shear_rate_tensor_exponent); 
}

void FluidViscosityCarreauYasudaParameters::print_parameters()
{
  std::cout << limiting_high_shear_rate_viscosity.name_ << ": " << limiting_high_shear_rate_viscosity.value_ << std::endl;
  std::cout << limiting_low_shear_rate_viscosity.name_ << ": " << limiting_low_shear_rate_viscosity.value_ << std::endl;
  std::cout << power_law_index.name_ << ": " << power_law_index.value_ << std::endl;
  std::cout << shear_rate_tensor_exponent.name_ << ": " << shear_rate_tensor_exponent.value_ << std::endl;
  std::cout << shear_rate_tensor_multipler.name_ << ": " << shear_rate_tensor_multipler.value_ << std::endl;
}

void FluidViscosityCarreauYasudaParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  std::string error_msg = "Unknown Constitutive_model type=CarreauYasuda XML element '";
  using std::placeholders::_1;
  using std::placeholders::_2;
  std::function<void(const std::string&, const std::string&)> ftpr =
      std::bind( &FluidViscosityCarreauYasudaParameters::set_parameter_value, *this, _1, _2);

  xml_util_set_parameters(ftpr, xml_elem, error_msg);
}

FluidViscosityCassonsParameters::FluidViscosityCassonsParameters()
{
  // A parameter that must be defined.
  bool required = true;

  set_parameter("Asymptotic_viscosity_parameter", 0.0, !required, asymptotic_viscosity);
  set_parameter("Low_shear_rate_threshold", 0.0, !required, low_shear_rate_threshold);
  set_parameter("Yield_stress_parameter", 0.0, !required, yield_stress);
}

void FluidViscosityCassonsParameters::print_parameters()
{
  std::cout << asymptotic_viscosity.name_ << ": " << asymptotic_viscosity.value_ << std::endl;
  std::cout << low_shear_rate_threshold.name_ << ": " << low_shear_rate_threshold.value_ << std::endl;
  std::cout << yield_stress.name_ << ": " << yield_stress.value_ << std::endl;
}

void FluidViscosityCassonsParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  std::string error_msg = "Unknown Constitutive_model type=Cassons XML element '";

  using std::placeholders::_1;
  using std::placeholders::_2;
  std::function<void(const std::string&, const std::string&)> ftpr =
      std::bind( &FluidViscosityCassonsParameters::set_parameter_value, *this, _1, _2);

  xml_util_set_parameters(ftpr, xml_elem, error_msg);
}

FluidViscosityParameters::FluidViscosityParameters()
{
  // A parameter that must be defined.
  bool required = true;

  // Stores model from the <Viscosity model= > XML element 
  model = Parameter<std::string>("model", "", required);
}

void FluidViscosityParameters::print_parameters()
{ 
  std::cout << std::endl;
  std::cout << "--------------------" << std::endl;
  std::cout << "Viscosity Parameters" << std::endl;
  std::cout << "--------------------" << std::endl;
  std::cout << model.name() << ": '" << model.value() << "'" << std::endl;

  // Print parameters for the given fluid_viscosity model.
  PrintFluidViscosityModelParamsMap[model.value_](this);
}

void FluidViscosityParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  using namespace tinyxml2;

  const char* smodel;
  auto result = xml_elem->QueryStringAttribute("model", &smodel);

  if (smodel == nullptr) {
    throw std::runtime_error("No MODEL given in the <Viscosity model=MODEL > XML element."); 
  }
  model.set(std::string(smodel));

  // Check fluid_viscosity model name.
  if (model_names.count(model.value()) == 0) { 
      throw std::runtime_error("Unknown fluid viscosity model '" + model.value() + 
        "' in '" + xml_elem->Name() + "'.");
  }

  // Set parameters for the given fluid_viscosity model.
  SetFluidViscosityModelParamsMap[model.value()](this, xml_elem);
}

//////////////////////////////////////////////////////////
//                  SolidViscosityParameters                 //
//////////////////////////////////////////////////////////

/// @brief Process parameters for various solid viscosity models.
///
/// Define the XML element name for viscosiity parameters.
const std::string SolidViscosityParameters::xml_element_name_ = "Viscosity";

const std::string SolidViscosityParameters::NEWTONIAN_MODEL = "Newtonian";
const std::string SolidViscosityParameters::POTENTIAL_MODEL = "Potential";

const std::set<std::string> SolidViscosityParameters::model_names = {
  SolidViscosityParameters::NEWTONIAN_MODEL,
  SolidViscosityParameters::POTENTIAL_MODEL
};

/// @brief Define a map to set parameters for each solid viscosity model.
using SVpType = SolidViscosityParameters*;
using XmlType = tinyxml2::XMLElement*;
using SetSolidViscosityParamMapType = std::map<std::string, std::function<void(SVpType, XmlType)>>;
SetSolidViscosityParamMapType SetSolidViscosityModelParamsMap = {
  {SolidViscosityParameters::NEWTONIAN_MODEL, [](SVpType vp, XmlType params) -> void { vp->newtonian_model.set_values(params); }},
  {SolidViscosityParameters::POTENTIAL_MODEL, [](SVpType vp, XmlType params) -> void { vp->potential_model.set_values(params); }},
};

/// @brief Define a map to print parameters for each solid viscosity model.
using PrintSolidViscosityParamaMapType = std::map<std::string, std::function<void(SVpType)>>;
PrintSolidViscosityParamaMapType PrintSolidViscosityModelParamsMap = {
  {SolidViscosityParameters::NEWTONIAN_MODEL, [](SVpType vp) -> void { vp->newtonian_model.print_parameters(); }},
  {SolidViscosityParameters::POTENTIAL_MODEL, [](SVpType vp) -> void { vp->potential_model.print_parameters(); }},
};

SolidViscosityNewtonianParameters::SolidViscosityNewtonianParameters()
{
  // A parameter that must be defined.
  bool required = true;

  set_parameter("Value", 0.0, !required, constant_value);
}

void SolidViscosityNewtonianParameters::print_parameters()
{
  std::cout << constant_value.name_ << ": " << constant_value.value_ << std::endl;
}

void SolidViscosityNewtonianParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  std::string error_msg = "Unknown Constitutive_model type=Newtonian XML element '";

  using std::placeholders::_1;
  using std::placeholders::_2;
  std::function<void(const std::string&, const std::string&)> ftpr =
      std::bind( &SolidViscosityNewtonianParameters::set_parameter_value, *this, _1, _2);

  xml_util_set_parameters(ftpr, xml_elem, error_msg);
}

SolidViscosityPotentialParameters::SolidViscosityPotentialParameters()
{
  // A parameter that must be defined.
  bool required = true;

  set_parameter("Value", 0.0, !required, constant_value);
}

void SolidViscosityPotentialParameters::print_parameters()
{
  std::cout << constant_value.name_ << ": " << constant_value.value_ << std::endl;
}

void SolidViscosityPotentialParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  std::string error_msg = "Unknown Constitutive_model type=Potential XML element '";
  using std::placeholders::_1;
  using std::placeholders::_2;
  std::function<void(const std::string&, const std::string&)> ftpr =
      std::bind( &SolidViscosityPotentialParameters::set_parameter_value, *this, _1, _2);

  xml_util_set_parameters(ftpr, xml_elem, error_msg);
}

SolidViscosityParameters::SolidViscosityParameters()
{
  // A parameter that must be defined.
  bool required = true;

  // Stores model from the <Viscosity model= > XML element 
  model = Parameter<std::string>("model", "", required);
}

void SolidViscosityParameters::print_parameters()
{ 
  std::cout << std::endl;
  std::cout << "--------------------" << std::endl;
  std::cout << "Viscosity Parameters" << std::endl;
  std::cout << "--------------------" << std::endl;
  std::cout << model.name() << ": '" << model.value() << "'" << std::endl;

  // Print parameters for the given solid_viscosity model.
  PrintSolidViscosityModelParamsMap[model.value_](this);
}

void SolidViscosityParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  using namespace tinyxml2;

  const char* smodel;
  auto result = xml_elem->QueryStringAttribute("model", &smodel);

  if (smodel == nullptr) {
    throw std::runtime_error("No MODEL given in the <Viscosity model=MODEL > XML element."); 
  }
  model.set(std::string(smodel));

  // Check solid viscosity model name.
  if (model_names.count(model.value()) == 0) { 
      throw std::runtime_error("Unknown solid viscosity model '" + model.value() + 
        "' in '" + xml_elem->Name() + "'.");
  }

  // Set parameters for the given solid viscosity model.
  SetSolidViscosityModelParamsMap[model.value()](this, xml_elem);
}

//////////////////////////////////////////////////////////
//                  DomainParameters                    //
//////////////////////////////////////////////////////////

// Process parameters for the XML // 'Domain' element to 
// specify properties for solving equations.

/// @brief Define the XML element name for domain parameters.
const std::string DomainParameters::xml_element_name_ = "Domain";

DomainParameters::DomainParameters()
{
  // A parameter that must be defined.
  bool required = true;

  // Set value from <Domain id=> element. 
  id = Parameter<std::string>("id", "", required);

  set_parameter("Absolute_tolerance", 1e-6, !required, absolute_tolerance);
  set_parameter("Anisotropic_conductivity", {}, !required, anisotropic_conductivity);
  set_parameter("Backflow_stabilization_coefficient", 0.2, !required, backflow_stabilization_coefficient);

  set_parameter("Conductivity", 0.0, !required, conductivity);
  //set_parameter("Constitutive_model", "", !required, constitutive_model);
  set_parameter("Continuity_stabilization_coefficient", 0.0, !required, continuity_stabilization_coefficient);

  set_parameter("Density", 0.5, !required, density);
  set_parameter("Dilational_penalty_model", "", !required, dilational_penalty_model);

  set_parameter("Elasticity_modulus", 1.0e7, !required, elasticity_modulus);
  set_parameter("Electrophysiology_model", "", !required, electrophysiology_model);
  set_parameter("Equation", "", !required, equation);

  set_parameter("Feedback_parameter_for_stretch_activated_currents", 0.5, !required, feedback_parameter_for_stretch_activated_currents);
  set_parameter("Fluid_density", 0.5, !required, fluid_density);
  set_parameter("Force_x", 0.0, !required, force_x);
  set_parameter("Force_y", 0.0, !required, force_y);
  set_parameter("Force_z", 0.0, !required, force_z);

  set_parameter("Isotropic_conductivity", 0.0, !required, isotropic_conductivity);

  set_parameter("Mass_damping", 0.0, !required, mass_damping);
  set_parameter("Maximum_iterations", 5, !required, maximum_iterations);
  set_parameter("Momentum_stabilization_coefficient", 0.0, !required, momentum_stabilization_coefficient);
  set_parameter("Myocardial_zone", "epicardium", !required, myocardial_zone);

  set_parameter("G_Na", 14.838, !required, G_Na);
  set_parameter("G_CaL", 3.98E-5, !required, G_CaL);
  set_parameter("G_Kr", 0.153, !required, G_Kr);
  set_parameter("G_Ks", 0.392, !required, G_Ks);
  set_parameter("G_to", 0.294, !required, G_to);

  set_parameter("tau_fi", 0.110, !required, tau_fi);
  set_parameter("tau_si", 1.88750, !required, tau_si);

  set_parameter("ODE_solver", "euler", !required, ode_solver);

  set_parameter("Penalty_parameter", 0.0, !required, penalty_parameter);
  set_parameter("Poisson_ratio", 0.3, !required, poisson_ratio);

  set_parameter("Relative_tolerance", 1e-4, !required, relative_tolerance);
  set_parameter("Shell_thickness", 0.0, !required, shell_thickness);
  set_parameter("Solid_density", 0.5, !required, solid_density);
  set_parameter("Source_term", 0.0, !required, source_term);
  set_parameter("Time_step_for_integration", 0.0, !required, time_step_for_integration);

  set_parameter("Inverse_darcy_permeability", 0.0, !required, inverse_darcy_permeability);
}

void DomainParameters::print_parameters()
{
  std::cout << std::endl;
  std::cout << "-----------------" << std::endl;
  std::cout << "Domain Parameters" << std::endl;
  std::cout << "-----------------" << std::endl;
  std::cout << id.name() << ": " << id.value() << std::endl;

  auto params_name_value = get_parameter_list();
  for (auto& [ key, value ] : params_name_value) { 
    std::cout << key << ": " << value << std::endl;
  }

  constitutive_model.print_parameters();

  fiber_reinforcement_stress.print_parameters();

  stimulus.print_parameters();

  fluid_viscosity.print_parameters();

  solid_viscosity.print_parameters();
}

void DomainParameters::set_values(tinyxml2::XMLElement* domain_elem)
{
  using namespace tinyxml2;
  std::string error_msg = "Unknown " + xml_element_name_ + " XML element '"; 

  const char* sid;
  auto result = domain_elem->QueryStringAttribute("id", &sid);
  if (sid == nullptr) {
    throw std::runtime_error("No ID found in the  <Domain id=ID> XML element.");
  }
  id.set(std::string(sid));

  auto item = domain_elem->FirstChildElement();
  
  // Parse XML elements for varius sub-elements (e.g. Viscosity).
  //
  while (item != nullptr) {
    auto name = std::string(item->Value());
  
    if (name == ConstitutiveModelParameters::xml_element_name_) {
      constitutive_model.set_values(item);

    } else if (name == FiberReinforcementStressParameters::xml_element_name_) {
      fiber_reinforcement_stress.set_values(item);

    } else if (name == StimulusParameters::xml_element_name_) {
      stimulus.set_values(item);

    } else if (name == FluidViscosityParameters::xml_element_name_ || name == SolidViscosityParameters::xml_element_name_) {
      auto eq_type = consts::equation_name_to_type.at(equation.value());
      if (eq_type == consts::EquationType::phys_fluid || eq_type == consts::EquationType::phys_CMM || eq_type == consts::EquationType::phys_stokes) {
        fluid_viscosity.set_values(item);
      } else if (eq_type == consts::EquationType::phys_struct || eq_type == consts::EquationType::phys_ustruct) {
        solid_viscosity.set_values(item);
      }
      else {
        throw std::runtime_error("Viscosity model not supported for equation '" + equation.value() + "'.");
      }
  
    } else if (item->GetText() != nullptr) {
      auto value = item->GetText();
      try {
        set_parameter_value(name, value);
      } catch (const std::bad_function_call& exception) {
        throw std::runtime_error(error_msg + name + "'.");
      }

    } else {
      throw std::runtime_error(error_msg + name + "'.");
    }
  
    item = item->NextSiblingElement();
  }

/*

  // Check values for some parameters..
  //
  if (Parameters::constitutive_model_names.count(constitutive_model.value()) == 0) {
    throw std::runtime_error("Unknown constitutive model '" + constitutive_model.value_ + "' for '" + constitutive_model.name_ + 
      "' in '" + domain_params->Name() + "'.");
  }

  if (Parameters::equation_names.count(equation.value()) == 0) {
    throw std::runtime_error("Unknown equation name '" + equation.value() + "' for '" + equation.name_ + 
      "' in '" + domain_params->Name() + "'.");
  }
*/
}

//////////////////////////////////////////////////////////
//            FiberReinforcementStressParameters        //
//////////////////////////////////////////////////////////

// Process parameters for the fiber reinforcement stress 'Fiber_reinforcement_stress` 
// XML element.

/// @brief Define the XML element name for fiber reinforcement stress parameters.
const std::string FiberReinforcementStressParameters::xml_element_name_ = "Fiber_reinforcement_stress";

FiberReinforcementStressParameters::FiberReinforcementStressParameters()
{
  // A parameter that must be defined.
  bool required = true;

  // Define attributes.
  type = Parameter<std::string>("type", "", required);

  set_parameter("Ramp_function", false, !required, ramp_function);
  set_parameter("Temporal_values_file_path", "", !required, temporal_values_file_path);
  set_parameter("Value", 0.0, !required, value);
}

void FiberReinforcementStressParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  using namespace tinyxml2;
  std::string error_msg = "Unknown " + xml_element_name_ + " XML element '";

  // Get the 'type' from the <LS type=TYPE> element.
  const char* stype;
  auto result = xml_elem->QueryStringAttribute("type", &stype);
  if (stype == nullptr) {
    throw std::runtime_error("No TYPE given in the XML <Stimulus=TYPE> element.");
  }
  type.set(std::string(stype));
  auto item = xml_elem->FirstChildElement();

  using std::placeholders::_1;
  using std::placeholders::_2;

  std::function<void(const std::string&, const std::string&)> ftpr =
      std::bind( &FiberReinforcementStressParameters::set_parameter_value, *this, _1, _2);

  xml_util_set_parameters(ftpr, xml_elem, error_msg);

  value_set = true;
}

void FiberReinforcementStressParameters::print_parameters()
{
  if (!value_set) {
    return;
  }
  std::cout << std::endl;
  std::cout << "-----------------------------------" << std::endl;
  std::cout << "FiberReinforcementStress Parameters" << std::endl;
  std::cout << "-----------------------------------" << std::endl;
  std::cout << type.name() << ": " << type.value() << std::endl;

  auto params_name_value = get_parameter_list();
  for (auto& [ key, value ] : params_name_value) {
    std::cout << key << ": " << value << std::endl;
  }
}

//////////////////////////////////////////////////////////
//                 StimulusParameters                   //
//////////////////////////////////////////////////////////

// The StimulusParameters class stores parameters for 
// 'Stimulus' XML element used to parameters for 
// pacemaker cells.

/// @brief Define the XML element name for equation output parameters.
const std::string StimulusParameters::xml_element_name_ = "Stimulus";

StimulusParameters::StimulusParameters()
{
  // A parameter that must be defined.
  bool required = true;

  // Define attributes.
  type = Parameter<std::string>("type", "", required);

  set_parameter("Amplitude", 0.0, !required, amplitude);
  set_parameter("Cycle_length", 0.0, !required, cycle_length);
  set_parameter("Duration", 0.0, !required, duration);
  set_parameter("Start_time", 0.0, !required, start_time);
}

void StimulusParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  std::string error_msg = "Unknown " + xml_element_name_ + " XML element '";

  // Get the 'type' from the <LS type=TYPE> element.
  const char* stype;
  auto result = xml_elem->QueryStringAttribute("type", &stype);
  if (stype == nullptr) {
    throw std::runtime_error("No TYPE given in the XML <Stimulus=TYPE> element.");
  }
  type.set(std::string(stype));
  auto item = xml_elem->FirstChildElement();
  
  using std::placeholders::_1;
  using std::placeholders::_2;

  std::function<void(const std::string&, const std::string&)> ftpr =
      std::bind( &StimulusParameters::set_parameter_value, *this, _1, _2);

  xml_util_set_parameters(ftpr, xml_elem, error_msg);

  value_set = true;
}

void StimulusParameters::print_parameters()
{
  if (!value_set) {
    return;
  }
  std::cout << std::endl;
  std::cout << "-------------------" << std::endl;
  std::cout << "Stimulus Parameters" << std::endl;
  std::cout << "-------------------" << std::endl;
  std::cout << type.name() << ": " << type.value() << std::endl;

  auto params_name_value = get_parameter_list();
  for (auto& [ key, value ] : params_name_value) { 
    std::cout << key << ": " << value << std::endl;
  }
}

//////////////////////////////////////////////////////////
//                 ECGLeadsParameters                   //
//////////////////////////////////////////////////////////

/// @brief Define the XML element name for ECG leads parameters.
const std::string ECGLeadsParameters::xml_element_name_ = "ECGLeads";

ECGLeadsParameters::ECGLeadsParameters()
{
  // A parameter that must be defined.
  bool required = true;

  // Define attributes.
  set_parameter("X_coords_file_path", "", !required, x_coords_file_path);
  set_parameter("Y_coords_file_path", "", !required, y_coords_file_path);
  set_parameter("Z_coords_file_path", "", !required, z_coords_file_path);
}

void ECGLeadsParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  std::string error_msg = "Unknown " + xml_element_name_ + " XML element '";

  using std::placeholders::_1;
  using std::placeholders::_2;

  std::function<void(const std::string&, const std::string&)> ftpr =
      std::bind( &ECGLeadsParameters::set_parameter_value, *this, _1, _2);

  xml_util_set_parameters(ftpr, xml_elem, error_msg);

  value_set = true;
}

void ECGLeadsParameters::print_parameters()
{
  if (!value_set) {
    return;
  }
  std::cout << std::endl;
  std::cout << "--------------------"  << std::endl;
  std::cout << "ECG Leads Parameters" << std::endl;
  std::cout << "--------------------"  << std::endl;

  auto params_name_value = get_parameter_list();
  for (auto& [ key, value ] : params_name_value) { 
    std::cout << key << ": " << value << std::endl;
  }
}

//////////////////////////////////////////////////////////
//                  ContactParameters                   //
//////////////////////////////////////////////////////////

/// @brief Process parameters for the 'Contact' XML element
/// used to specify parameters for contact computation.
///
/// Define the XML element name for contact parameters.
const std::string ContactParameters::xml_element_name_ = "Contact";

ContactParameters::ContactParameters()
{
  set_xml_element_name(xml_element_name_);

  // A parameter that must be defined.
  bool required = true;

  // Contact model.
  model = Parameter<std::string>("model", "", required);

  // Define contact parameters.
  //
  set_parameter("Closest_gap_to_activate_penalty", 1.0, !required, closest_gap_to_activate_penalty);
  set_parameter("Desired_separation", 0.05, !required, desired_separation);
  set_parameter("Min_norm_of_face_normals", 0.7, !required, min_norm_of_face_normals);
  set_parameter("Penalty_constant", 1e5, !required, penalty_constant);
}

void ContactParameters::print_parameters()
{
  std::cout << std::endl;
  std::cout << "-------------------" << std::endl;
  std::cout << "Contact Parameters" << std::endl;
  std::cout << "-------------------" << std::endl;
  std::cout << model.name() << ": " << model.value() << std::endl;

  auto params_name_value = get_parameter_list();
  for (auto& [ key, value ] : params_name_value) {
    std::cout << key << ": " << value << std::endl;
  }
}

void ContactParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  using namespace tinyxml2;
  std::string error_msg = "Unknown " + xml_element_name_ + " XML element '";

  // Get the 'type' from the <Add_projection name=NAME> element.
  const char* mname;
  auto result = xml_elem->QueryStringAttribute("model", &mname);
  if (mname == nullptr) {
    throw std::runtime_error("No MODEL given in the XML <Contact model=MODEL> element.");
  }
  model.set(std::string(mname));

  using std::placeholders::_1;
  using std::placeholders::_2;

  std::function<void(const std::string&, const std::string&)> ftpr =
      std::bind( &ProjectionParameters::set_parameter_value, *this, _1, _2);

  xml_util_set_parameters(ftpr, xml_elem, error_msg);
}

//////////////////////////////////////////////////////////
//                  EquationParameters                  //
//////////////////////////////////////////////////////////

// Process parameters for the 'Add_equation' XML element 
// used to specify an equation to be solved (e.g. fluid).

/// @brief Define the XML element name for equation parameters.
const std::string EquationParameters::xml_element_name_ = "Add_equation";

EquationParameters::EquationParameters()
{
  set_xml_element_name(xml_element_name_);

  // A parameter that must be defined.
  bool required = true;

  // Equation type.
  type = Parameter<std::string>("type", "", required);

  // Define equation parameters.
  //
  set_parameter("Coupled", false, !required, coupled);

  set_parameter("Initialize", "", !required, initialize);
  set_parameter("Initialize_RCR_from_flow", false, !required, initialize_rcr_from_flow);

  set_parameter("Max_iterations", 1, !required, max_iterations);
  set_parameter("Min_iterations", 1, !required, min_iterations);

  set_parameter("Prestress", false, !required, prestress);

  set_parameter("Tolerance", 0.5, !required, tolerance);
  set_parameter("Use_taylor_hood_type_basis", false, !required, use_taylor_hood_type_basis);
}

void EquationParameters::print_parameters()
{ 
  std::cout << std::endl;
  std::cout << "-------------------" << std::endl;
  std::cout << "Equation Parameters" << std::endl;
  std::cout << "-------------------" << std::endl;
  std::cout << type.name() << ": " << type.value() << std::endl;

  auto params_name_value = get_parameter_list();
  for (auto& [ key, value ] : params_name_value) { 
    std::cout << key << ": " << value << std::endl;
  }

  if (domains.size() == 0) {
    default_domain->print_parameters();
  } else {
    for (auto& domain : domains) {
      domain->print_parameters();
    }
  } 

  //stimulus.print_parameters();

  for (auto& output : outputs) {
    output->print_parameters();
  }

  linear_solver.print_parameters();

  couple_to_cplBC.print_parameters();

  for (auto& bc : boundary_conditions) {
    bc->print_parameters();
  }

  for (auto& bf : body_forces) {
    bf->print_parameters();
  }

  ecg_leads.print_parameters();
}

void EquationParameters::set_values(tinyxml2::XMLElement* eq_elem)
{
  using namespace tinyxml2;
  default_domain = new DomainParameters();
  auto item = eq_elem->FirstChildElement();

  // Parse XML sub-elements.
  //
  // [TODO] Replace all of these ifs with a string/lambda map.
  //
  while (item != nullptr) {
    auto name = std::string(item->Value());

    if (name == BodyForceParameters::xml_element_name_) {
      auto bf_params = new BodyForceParameters();
      bf_params->set_values(item);
      body_forces.push_back(bf_params);

    } else if (name == BoundaryConditionParameters::xml_element_name_) {
      auto bc_params = new BoundaryConditionParameters();
      bc_params->set_values(item);
      boundary_conditions.push_back(bc_params);

    } else if (name == ConstitutiveModelParameters::xml_element_name_) {
      default_domain->constitutive_model.set_values(item);
      default_domain->constitutive_model.check_constitutive_model(type);

    } else if (name == CoupleCplBCParameters::xml_element_name_) {
      couple_to_cplBC.set_values(item);

    } else if (name == CoupleGenBCParameters::xml_element_name_) {
      couple_to_genBC.set_values(item);

    } else if (name == CoupleSvZeroDParameters::xml_element_name_) {
      couple_to_svZeroD.set_values(item);

    } else if (name == DomainParameters::xml_element_name_) {
      auto domain_params = new DomainParameters();
      domain_params->set_values(item);
      domains.push_back(domain_params);

    } else if (name == FiberReinforcementStressParameters::xml_element_name_) {
      default_domain->fiber_reinforcement_stress.set_values(item);

    } else if (name == LinearSolverParameters::xml_element_name_) {
      linear_solver.set_values(item);

    } else if (name == OutputParameters::xml_element_name_) {
      auto output_params = new OutputParameters();
      output_params->set_values(item);
      outputs.push_back(output_params);

    } else if (name == RemesherParameters::xml_element_name_) {
      remesher.set_values(item);

    } else if (name == StimulusParameters::xml_element_name_) {
      default_domain->stimulus.set_values(item);

    } else if (name == FluidViscosityParameters::xml_element_name_ || name == SolidViscosityParameters::xml_element_name_) {
      auto eq_type = consts::equation_name_to_type.at(type.value());
      if (eq_type == consts::EquationType::phys_fluid || eq_type == consts::EquationType::phys_CMM || eq_type == consts::EquationType::phys_stokes) {
        default_domain->fluid_viscosity.set_values(item);
      } else if (eq_type == consts::EquationType::phys_struct || eq_type == consts::EquationType::phys_ustruct) {
        default_domain->solid_viscosity.set_values(item);
      } else {
        throw std::runtime_error("Viscosity model not supported for equation '" + type.value() + "'.");
      }

    } else if (name == ECGLeadsParameters::xml_element_name_) {
      ecg_leads.set_values(item);

    } else if (name == VariableWallPropsParameters::xml_element_name_) {
      variable_wall_properties.set_values(item);

    } else if (item->GetText() != nullptr) {
      auto value = item->GetText();

      // A parameter can be an EqutionParameter or DomainParameter.
      //
      try {
        set_parameter_value(name, value);
      } catch (const std::exception &exception) {

        try {
          default_domain->set_parameter_value(name, value);
        } catch (const std::bad_function_call& exception) {
          throw std::runtime_error("Unknown " + xml_element_name_ + " XML element '" + name + ".");
        }
      }

    } else {
      throw std::runtime_error("Unknown " + xml_element_name_ + " XML element '" + name + ".");
    }

    item = item->NextSiblingElement();
  }

  /*
  if (domains.size() == 0) { 
    auto domain_params = new DomainParameters();
    domain_params->set_values(item);
    domains.push_back(domain_params);
  }
  */
}

//////////////////////////////////////////////////////////
//               GeneralSimulationParameters            //
//////////////////////////////////////////////////////////

/// @brief Process paramaters for the 'GeneralSimulationParameters' XML element.
GeneralSimulationParameters::GeneralSimulationParameters()
{
  int int_inf = std::numeric_limits<int>::infinity();

  // Define the XML element name for general simulation parameters.
  xml_element_name = "GeneralSimulationParameters";
  set_xml_element_name(xml_element_name);

  // A parameter that must be defined.
  bool required = true;

  set_parameter("Check_IEN_order", true, !required, check_ien_order);
  set_parameter("Continue_previous_simulation", false, required, continue_previous_simulation);
  set_parameter("Convert_BIN_to_VTK_format", false, !required, convert_bin_to_vtk_format);

  set_parameter("Debug", false, !required, debug);

  set_parameter("Increment_in_saving_restart_files", 0, !required, increment_in_saving_restart_files);
  set_parameter("Increment_in_saving_VTK_files", 0, !required, increment_in_saving_vtk_files);

  set_parameter("Name_prefix_of_saved_VTK_files", "", !required, name_prefix_of_saved_vtk_files);
  set_parameter("Number_of_initialization_time_steps", 0, !required, number_of_initialization_time_steps, {0,int_inf});
  set_parameter("Number_of_spatial_dimensions", 3, !required, number_of_spatial_dimensions);
  set_parameter("Number_of_time_steps", 0, required, number_of_time_steps, {0,int_inf});

  set_parameter("Overwrite_restart_file", false, !required, overwrite_restart_file);

  set_parameter("Restart_file_name", "stFile", !required, restart_file_name);

  set_parameter("Save_averaged_results", false, !required, save_averaged_results);
  set_parameter("Save_results_in_folder", "", !required, save_results_in_folder);
  set_parameter("Save_results_to_VTK_format", false, required, save_results_to_vtk_format);
  set_parameter("Searched_file_name_to_trigger_stop", "", !required, searched_file_name_to_trigger_stop);
  set_parameter("Simulation_initialization_file_path", "", !required, simulation_initialization_file_path);
  set_parameter("Simulation_requires_remeshing", false, !required, simulation_requires_remeshing);
  set_parameter("Spectral_radius_of_infinite_time_step", 0.5, required, spectral_radius_of_infinite_time_step);
  set_parameter("Start_averaging_from_zero", false, !required, start_averaging_from_zero);
  set_parameter("Start_saving_after_time_step", 0, required, start_saving_after_time_step);
  set_parameter("Starting time step", 0, !required, starting_time_step);

  set_parameter("Time_step_size", 0.0, required, time_step_size);
  set_parameter("Verbose", false, !required, verbose);
  set_parameter("Warning", false, !required, warning);
}

void GeneralSimulationParameters::print_parameters()
{
  std::cout << std::endl;
  std::cout << "-----------------------------" << std::endl;
  std::cout << "General Simulation Parameters" << std::endl;
  std::cout << "-----------------------------" << std::endl;

  auto params_name_value = get_parameter_list();
  for (auto& [ key, value ] : params_name_value) { 
    std::cout << key << ": " << value << std::endl;
  }
}

/// @brief Set general parameters values from XML.
void GeneralSimulationParameters::set_values(tinyxml2::XMLElement* xml_element)
{
  using namespace tinyxml2;

  // Set parameter values from the XML elements.
  //
  auto general_params = xml_element->FirstChildElement(xml_element_name.c_str());
  auto item = general_params->FirstChildElement();

  while (item != nullptr) {
    std::string name = std::string(item->Value());
    auto value = item->GetText();
    try {
      set_parameter_value(name, value);
    } catch (const std::bad_function_call& exception) {
      throw std::runtime_error("Unknown XML GeneralSimulationParameters element '" + name + ".");
    }
    item = item->NextSiblingElement();
  }

  // Check that required parameters have been set.
  check_required();
}

//////////////////////////////////////////////////////////
//             F a c e P a r a m e t e r s              //
//////////////////////////////////////////////////////////

/// @brief Process parameters for the 'Add_face' XML element.
///
/// Define the XML element name for face parameters.
const std::string FaceParameters::xml_element_name_ = "Add_face";

FaceParameters::FaceParameters()
{
  set_xml_element_name(xml_element_name_);

  // A parameter that must be defined.
  bool required = true;

  name = Parameter<std::string>("name", "", required);

  set_parameter("End_nodes_face_file_path", "", !required, end_nodes_face_file_path);
  set_parameter("Face_file_path", "", !required, face_file_path);

  set_parameter("Quadrature_modifier_TRI3", (2.0/3.0), !required, quadrature_modifier_TRI3);
}

void FaceParameters::print_parameters()
{
  std::cout << std::endl;
  std::cout << "---------------" << std::endl;
  std::cout << "Face Parameters" << std::endl;
  std::cout << "---------------" << std::endl;
  std::cout << name.name() << ": " << name.value() << std::endl;
  std::cout << face_file_path.name() << ": " << face_file_path.value() << std::endl;
}

void FaceParameters::set_values(tinyxml2::XMLElement* face_elem)
{
  using namespace tinyxml2;

  std::string error_msg = "Unknown " + xml_element_name_ + " XML element '"; 
  const char* face_name;
  auto result = face_elem->QueryStringAttribute("name", &face_name);
  name.set(std::string(face_name));
  auto item = face_elem->FirstChildElement();

  while (item != nullptr) {
    auto name = std::string(item->Value());
    auto value = item->GetText();

    if (value == nullptr) { 
      throw std::runtime_error(error_msg + name + "'.");
    }

    try {
      set_parameter_value(name, value);
    } catch (const std::bad_function_call& exception) {
      throw std::runtime_error(error_msg + name + "'.");
    }

    item = item->NextSiblingElement();
  }
}

//////////////////////////////////////////////////////////
//           R e m e s h e r P a r a m e t e r s        //
//////////////////////////////////////////////////////////

// Process parameters for the 'Remesher' XML element used for remeshing.

/// @brief Define the XML element name for mesh parameters.
const std::string RemesherParameters::xml_element_name_ = "Remesher";

RemesherParameters::RemesherParameters()
{
  bool required = true;

  type = Parameter<std::string>("type", "", required);

  set_parameter("Min_dihedral_angle", 10.0, !required, min_dihedral_angle);
  set_parameter("Max_radius_ratio", 1.15, !required, max_radius_ratio);
  set_parameter("Remesh_frequency", 100, !required, remesh_frequency);
  set_parameter("Frequency_for_copying_data", 10, !required, frequency_for_copying_data);
}

void RemesherParameters::print_parameters()
{
  std::cout << std::endl;
  std::cout << "-------------------" << std::endl;
  std::cout << "Remesher Parameters" << std::endl;
  std::cout << "-------------------" << std::endl;

  std::cout << type.name() << ": " << type.value() << std::endl;
}

void RemesherParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  std::string error_msg = "Unknown " + xml_element_name + " XML element '";

  // Get the 'type' from the <Remesher type=TYPE> element.
  const char* stype;
  auto result = xml_elem->QueryStringAttribute("type", &stype);
  if (stype == nullptr) {
    throw std::runtime_error("No TYPE given in the XML <Remesher type=TYPE> element.");
  }
  type.set(std::string(stype));
  values_set_ = true;

  // Iterate over sub-elements.
  //
  auto item = xml_elem->FirstChildElement();

  while (item != nullptr) {
    auto name = std::string(item->Value());

    if (name == "Max_edge_size") {
      const char* name;
      const char* value;
      auto result = item->QueryStringAttribute("name", &name);
      if (name == nullptr) {
        throw std::runtime_error("No NAME given in the XML Remesher <Max_edge_size name=NAME  value=VALUE> element.");
      }

      result = item->QueryStringAttribute("value", &value);
      if (value == nullptr) {
        throw std::runtime_error("No VALUE given in the XML Remesher <Max_edge_size name=NAME  value=VALUE> element.");
      }
      auto svalue = std::string(value);

      try {
        double dvalue = std::stod(svalue);
        max_edge_sizes_[std::string(name)] = dvalue;
      } catch (...) {
        throw std::runtime_error("VALUE=" + svalue + 
            " is not a valid float in the XML Remesher <Max_edge_size name=NAME  value=VALUE> element.");
      }

    } else if (item->GetText() != nullptr) {
      auto value = item->GetText();
      try {
        set_parameter_value(name, value);
      } catch (const std::bad_function_call& exception) {
        throw std::runtime_error(error_msg + name + "'.");
      }

    } else {
      throw std::runtime_error(error_msg + name + "'.");
    }

    item = item->NextSiblingElement();
  }
}

//////////////////////////////////////////////////////////
//             M e s h P a r a m e t e r s              //
//////////////////////////////////////////////////////////

// Process parameters for the 'Add_mesh' XML element used for defining mesh elements.

/// @brief Define the XML element name for mesh parameters.
const std::string MeshParameters::xml_element_name_ = "Add_mesh";

MeshParameters::MeshParameters()
{
  bool required = true;

  // Mesh name from Add_mesh element.
  name = Parameter<std::string>("name", "", required);

  // Parameters under Add_mesh element.
  //
  set_parameter("Domain", 0,  !required, domain_id);
  set_parameter("Domain_file_path", "", !required, domain_file_path);

  //set_parameter("Fiber_direction", {}, !required, fiber_direction);
  set_parameter("Fiber_direction_file_path", {}, !required, fiber_direction_file_paths);

  set_parameter("Mesh_file_path", "", !required, mesh_file_path);
  set_parameter("Mesh_scale_factor", 1.0, !required, mesh_scale_factor);
  set_parameter("Prestress_file_path", "", !required, prestress_file_path);

  set_parameter("Initial_displacements_file_path", "", !required, initial_displacements_file_path);
  set_parameter("Initial_pressures_file_path", "", !required, initial_pressures_file_path);
  set_parameter("Initial_velocities_file_path", "", !required, initial_velocities_file_path);

  set_parameter("Set_mesh_as_fibers", false, !required, set_mesh_as_fibers);
  set_parameter("Set_mesh_as_shell", false, !required, set_mesh_as_shell);

  set_parameter("Quadrature_modifier_TET4", (5.0+3.0*sqrt(5.0))/20.0, !required, quadrature_modifier_TET4);
}

void MeshParameters::print_parameters()
{
  std::cout << std::endl;
  std::cout << "---------------" << std::endl;
  std::cout << "Mesh Parameters" << std::endl;
  std::cout << "---------------" << std::endl;
  std::cout << name.name() << ": " << name.value() << std::endl;

  auto params_name_value = get_parameter_list();
  for (auto& [ key, value ] : params_name_value) { 
    std::cout << key << ": " << value << std::endl;
  }

  for (auto& dir : fiber_directions) { 
    std::cout << dir.name() << ": " << dir.svalue() << std::endl;
  }

  for (auto& face : face_parameters) {
    face->print_parameters();
  }
}

void MeshParameters::set_values(tinyxml2::XMLElement* mesh_elem)
{
  using namespace tinyxml2;
  std::string error_msg = "Unknown " + xml_element_name_ + " XML element '"; 
  auto item = mesh_elem->FirstChildElement();

  while (item != nullptr) {
    auto name = std::string(item->Value());

    // Add_face sub-element.
    if (name == FaceParameters::xml_element_name_) {
      auto face_params = new FaceParameters();
      face_params->set_values(item);
      face_parameters.push_back(face_params);

    // There may be multiple 'Fiber_direction' elements so store
    // them as a list of VectorParameter<double>. 
    //
    } else if (name == "Fiber_direction") {
      auto value = item->GetText();
      VectorParameter<double> dir("Fiber_direction", {}, false, {});
      dir.set(value);
      fiber_directions.push_back(dir);

    // Just a simple element. 
    } else if (item->GetText() != nullptr) {
      auto value = item->GetText();
      try {
        set_parameter_value(name, value);
      } catch (const std::bad_function_call& exception) {
        throw std::runtime_error(error_msg + name + "'.");
      }
    } else {
      throw std::runtime_error(error_msg + name + "'.");
    }

    item = item->NextSiblingElement();
  }
}

/////////////////////////////////////////////////////////////////////////////
//        P r e c o m p u t e d S o l u t i o n P a r a m e t e r s       //
/////////////////////////////////////////////////////////////////////////////

// The PrecomputedSolutionParameters class stores parameters for the
// 'Precomputed_solution' XML element used to read in the data from a 
// precomputed solution for the simulation state.

const std::string PrecomputedSolutionParameters::xml_element_name_ = "Precomputed_solution";

PrecomputedSolutionParameters::PrecomputedSolutionParameters()
{
  // A parameter that must be defined.
  bool required = true;

  set_parameter("Field_name", "", required, field_name);
  set_parameter("File_path", "", required, file_path);
  set_parameter("Time_step", 0.0, !required, time_step);
  set_parameter("Use_precomputed_solution", false, !required, use_precomputed_solution);
}

void PrecomputedSolutionParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  using namespace tinyxml2;
  std::string error_msg = "Unknown " + xml_element_name_ + " XML element '";
  using std::placeholders::_1;
  using std::placeholders::_2;

  std::function<void(const std::string&, const std::string&)> ftpr =
      std::bind( &PrecomputedSolutionParameters::set_parameter_value, *this, _1, _2);

  xml_util_set_parameters(ftpr, xml_elem, error_msg);
}

//////////////////////////////////////////////////////////
//        P r o j e c t i o n P a r a m e t e r s       //
//////////////////////////////////////////////////////////

// The ProjectionParameters class stores parameters for the
// 'Add_projection' XML element used for fluid-structure interaction simulations.

/// @brief Define the XML element name for mesh parameters.
const std::string ProjectionParameters::xml_element_name_ = "Add_projection";

ProjectionParameters::ProjectionParameters()
{
  // A parameter that must be defined.
  bool required = true;

  name = Parameter<std::string>("name", "", required);

  set_parameter("Project_from_face", "", required, project_from_face);
  set_parameter("Projection_tolerance", 0.0, !required, projection_tolerance);
}

void ProjectionParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  using namespace tinyxml2;
  std::string error_msg = "Unknown " + xml_element_name_ + " XML element '";

  // Get the 'type' from the <Add_projection name=NAME> element.
  const char* sname;
  auto result = xml_elem->QueryStringAttribute("name", &sname);
  if (sname == nullptr) {
    throw std::runtime_error("No TYPE given in the XML <Add_projection name=NAME> element.");
  }
  name.set(std::string(sname));

  using std::placeholders::_1;
  using std::placeholders::_2;

  std::function<void(const std::string&, const std::string&)> ftpr =
      std::bind( &ProjectionParameters::set_parameter_value, *this, _1, _2);

  xml_util_set_parameters(ftpr, xml_elem, error_msg);
}

//////////////////////////////////////////////////////////
//                 LinearAlgebraParameters              //
//////////////////////////////////////////////////////////

// The LinearAlgebraParameters class stores parameters for
// the 'Linear_algebra' XML element.

/// @brief Define the XML element name for equation output parameters.
const std::string LinearAlgebraParameters::xml_element_name_ = "Linear_algebra";

LinearAlgebraParameters::LinearAlgebraParameters()
{
  // A parameter that must be defined.
  bool required = true;

  auto alg_type = LinearAlgebra::type_to_name.at(consts::LinearAlgebraType::fsils);
  type = Parameter<std::string>("type", alg_type, required);

  set_parameter("Configuration_file", "", !required, configuration_file);

  auto prec_type = consts::preconditioner_type_to_name.at(consts::PreconditionerType::PREC_NONE);
  set_parameter("Preconditioner", prec_type, !required, preconditioner);

  auto assemble_type = LinearAlgebra::type_to_name.at(consts::LinearAlgebraType::none);
  set_parameter("Assembly", assemble_type, !required, assembly);
}

void LinearAlgebraParameters::print_parameters()
{ 
  std::cout << std::endl;
  std::cout << "-------------------------" << std::endl;
  std::cout << "Linear Algebra Parameters" << std::endl;
  std::cout << "-------------------------" << std::endl;
	  
  std::cout << type.name() << ": " << type.value() << std::endl;
  
  auto params_name_value = get_parameter_list(); 
  for (auto& [ key, value ] : params_name_value) {
    std::cout << key << ": " << value << std::endl;
  }
}

void LinearAlgebraParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  std::string error_msg = "Unknown " + xml_element_name + " XML element '";

  // Get the 'type' from the <Linear_algebra type=TYPE> element.
  const char* stype;
  auto result = xml_elem->QueryStringAttribute("type", &stype);
  if (stype == nullptr) {
    throw std::runtime_error("No TYPE given in the XML <Linear_algebra type=TYPE> element.");
  }
  type.set(std::string(stype));

  // Check Linear_algebra type=TYPE> element.
  //
  // TYPE in (fsils petsc trilinos)
  //
  if (LinearAlgebra::name_to_type.count(type.value()) == 0) {
    std::string valid_types = "";
    std::for_each(LinearAlgebra::name_to_type.begin(), LinearAlgebra::name_to_type.end(), 
        [&valid_types](std::pair<const std::string, const consts::LinearAlgebraType> p) {valid_types += p.first+" ";}); 
    throw std::runtime_error("Unknown TYPE '" + type.value() + 
        "' given in the XML <Linear_algebra type=TYPE> element.\nValid types are: " + valid_types);
  }

  // Create a function pointer 'fptr' to 'LinearAlgebraParameters::set_parameter_value'.
  //
  using std::placeholders::_1;
  using std::placeholders::_2;
  std::function<void(const std::string&, const std::string&)> ftpr =
      std::bind( &LinearAlgebraParameters::set_parameter_value, *this, _1, _2);
      
  // Parse XML and set parameter values.
  xml_util_set_parameters(ftpr, xml_elem, error_msg);

  // Check preconditioner type.
  if (consts::preconditioner_name_to_type.count(preconditioner.value()) == 0) {
    std::string valid_types = "";
    std::for_each(consts::preconditioner_name_to_type.begin(), consts::preconditioner_name_to_type.end(),
        [&valid_types](std::pair<const std::string, const consts::PreconditionerType> p) {valid_types += p.first+" ";});
    throw std::runtime_error("Unknown TYPE '" + preconditioner() + 
        "' given in the XML <Linear_algebra> <Preconditioner> element.\nValid types are: " + valid_types);
  }     

  check_input_parameters();

  values_set_ = true;
}

/// @brief Check the validity of the input parameters.
void LinearAlgebraParameters::check_input_parameters()
{
  auto linear_algebra_type = LinearAlgebra::name_to_type.at(type());
  auto prec_cond_type = consts::preconditioner_name_to_type.at(preconditioner.value()); 
  auto assembly_type = LinearAlgebra::name_to_type.at(assembly.value()); 

  auto linear_algebra = LinearAlgebraFactory::create_interface(linear_algebra_type);
  linear_algebra->check_options(prec_cond_type, assembly_type);
}

//////////////////////////////////////////////////////////
//                 LinearSolverParameters               //
//////////////////////////////////////////////////////////

// The LinearSolverParameters class stores parameters for
// the 'LS' XML element.

/// @brief Define the XML element name for equation output parameters.
const std::string LinearSolverParameters::xml_element_name_ = "LS";

LinearSolverParameters::LinearSolverParameters()
{
  // A parameter that must be defined.
  bool required = true;

  // Define type attribute.
  type = Parameter<std::string>("type", "", required);

  set_parameter("Absolute_tolerance", 1.0e-10, !required, absolute_tolerance);

  set_parameter("Krylov_space_dimension", 50, !required, krylov_space_dimension);

  set_parameter("Max_iterations", 1000, !required, max_iterations);

  set_parameter("NS_CG_max_iterations", 1000, !required, ns_cg_max_iterations);
  set_parameter("NS_CG_tolerance", 1.0e-2, !required, ns_cg_tolerance);
  set_parameter("NS_GM_max_iterations", 1000, !required, ns_gm_max_iterations);
  set_parameter("NS_GM_tolerance", 1.0e-2, !required, ns_gm_tolerance);

  //set_parameter("Preconditioner", "", !required, preconditioner);

  set_parameter("Tolerance", 0.5, !required, tolerance);
}

void LinearSolverParameters::print_parameters()
{
  std::cout << std::endl;
  std::cout << "------------------------" << std::endl;
  std::cout << "Linear Solver Parameters" << std::endl;
  std::cout << "------------------------" << std::endl;

  std::cout << type.name() << ": " << type.value() << std::endl;

  auto params_name_value = get_parameter_list();
  for (auto& [ key, value ] : params_name_value) {
    std::cout << key << ": " << value << std::endl;
  }
}

void LinearSolverParameters::set_values(tinyxml2::XMLElement* xml_elem)
{
  std::string error_msg = "Unknown " + xml_element_name + " XML element '";

  // Get the 'type' from the <LS type=TYPE> element.
  const char* stype;
  auto result = xml_elem->QueryStringAttribute("type", &stype);
  if (stype == nullptr) {
    throw std::runtime_error("No TYPE given in the XML <LStype=TYPE> element.");
  }

  type.set(std::string(stype));

  using std::placeholders::_1;
  using std::placeholders::_2;

  // Create a function pointer 'fptr' to 'LinearSolverParameters::set_parameter_value'.
  std::function<void(const std::string&, const std::string&)> ftpr =
      std::bind( &LinearSolverParameters::set_parameter_value, *this, _1, _2);

  // Parse XML and set parameter values.
  std::set<std::string> sub_sections = {LinearAlgebraParameters::xml_element_name_};
  xml_util_set_parameters(ftpr, xml_elem, error_msg, sub_sections);

  // Set subsection values.
  //
  auto item = xml_elem->FirstChildElement();

  while (item != nullptr) {
    auto name = std::string(item->Value());
    if (name == LinearAlgebraParameters::xml_element_name_) {
      linear_algebra.set_values(item);
    }
    item = item->NextSiblingElement();
  }

}

