// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PARAMETERS_H 
#define PARAMETERS_H 

#include <any>
#include <functional>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

#include "tinyxml2.h"

template<typename T>

/// @brief The Parameter class template is used to store a named 
/// paramater and its scalar value as a basic type: bool, double, 
/// int and string.
///
/// The classes defined here are used to process svFSIplus simulation parameters read in 
/// from an Extensible Markup Language (XML) format file. XML is a simple text-based format 
/// for representing structured information. 
///
/// An XML document is formed as an element tree. The XML tree starts at a root element and 
/// branches from the root to sub-elements. All elements can have sub-elements:
///
/// \code{.cpp}
/// <svMultiPhysicsFile>
///   <element>
///     <subelement>.....</subelement>
///   </element>
/// </svMultiPhysicsFile>
/// \endcode
///
/// The elements in the svFSIplus simulation file are represented by sections of
/// related parameters. Sub-elements are refered to as sub-sections.
///
///-----------------
/// Parameters class 
///-----------------
/// The Parameters class is the top level class. It contains objects used to store 
/// parameters for the sections making up an XML simulation parameters file
///
///   1) General     (GeneralSimulationParameters)
///   2) Mesh        (MeshParameters)
///   3) Equation    (EquationParameters)
///   4) Projection  (ProjectionParameters)
///   5) RIS Projection  (RIS ProjectionParameters)
///
/// Each object contains methods to parse the XML file for the parameters defined for it.
/// These section objects may also contain objects representing the sub-sections defined 
/// for each section.
///
///-----------------
/// Section objects 
///-----------------
/// Each section object contains objects representing parameters. A parameter's name and value 
/// is stored using either a Parameter and VectorParamater template objects. A parameter 
/// value is stored as a basic type: bool, double, int and string. 
///
/// Parameter objects in a section class are named using the same XML element name with the 1st 
/// character lower case.
///
/// Example: GeneralSimulationParameters class 
///
/// \code{.cpp}
///   Parameter<bool> verbose;                                  // <Verbose> 
///   Parameter<double> spectral_radius_of_infinite_time_step;  // <Spectral_radius_of_infinite_time_step>
///   Parameter<double> time_step_size;                         // <Time_step_size> 
/// \endcode
///
/// The name and default value for each parameter is defined in a section object's constructor.
///
/// Parameter values are set using the set_values() method which contains calls to tinyxml2  
/// to parse parameter values from an XML file.
///
/// Each section object inherits from the ParameterLists class. This class provides methods to
/// store parameters in a map and process iterated over them for setting values and other operations.
class Parameter
{
  public:
    Parameter() {};

    Parameter(const std::string& name, T value, bool required, std::vector<T> range = {}) :
      value_(value), name_(name), required_(required)
    { 
      value_ = value;
      range_ = range;
    };

    std::string name() const { return name_; };
    T value() const { return value_; };
    T operator()() const { return value_; };
    bool defined() const { return value_set_; };

    /// @brief Get the value of a parameter as a string.
    std::string svalue()
    {
      std::ostringstream str_stream;
      str_stream << value_;
      return str_stream.str();
    }

    friend std::ostream& operator << (std::ostream& out, const Parameter<T>& param)
    {
      out << param.value();
      return out;
    }

    /// @brief Set the parameter name and value, and if it is required.
    void set(const std::string& name, bool required, T value) { 
      name_ = name;
      required_ = required;
      value_ = value;
    }

    /// @brief Set the parameter value from a string.
    void set(const std::string& str_value)
    {
      if (str_value == "") {
        value_ = T{0};
      }

      auto str = str_value;
      std::string::iterator end_pos = std::remove(str.begin(), str.end(), ' ');
      str.erase(end_pos, str.end());

      std::istringstream str_stream(str);
      if (!(str_stream >> value_)) {
        std::istringstream str_stream(str);
        if (!(str_stream >> std::boolalpha >> value_)) {
          throw std::runtime_error("Incorrect value '" + str + "' for '" + name_ + "'.");
        }
      }

      value_set_ = true;
    }

    bool check_required_set()
    {
      if (!required_) {
        return true;
      }
      return value_set_;
    }

    T value_ = T{0};
    std::string name_ = "";
    bool required_ = false;
    bool value_set_ = false;
    std::vector<T> range_;
};

/// @brief The VectorParameter class template is used to store a named 
/// paramater and its vector of values as a basic type: bool, double, 
/// int and string.
template<typename T>
class VectorParameter
{
  public:
    VectorParameter() {};

    VectorParameter(const std::string& name, const std::vector<T>& value, bool required, std::vector<T> range = {}) :
      value_(value), name_(name), required_(required)
    { 
      value_ = value;
      range_ = range;
    };

    std::string name() const { return name_; };
    std::vector<T> value() const { return value_; };
    bool defined() const { return value_set_; };
    int size() const { return value_.size(); };

    std::vector<T> operator()() const { return value_; };
    const double& operator[](const int i) const { return value_[i]; };

    /// @brief Get the string representation of the parameter value.
    std::string svalue() 
    {
      std::string str;

      if constexpr (std::is_same<T, std::string>::value) {
        for (auto v : value_) {
          str += " " + v + " ";
        }
      } else {
        for (auto v : value_) {
          str += " " + std::to_string(v);
        }
      }

      return str;
    }

    friend std::ostream& operator << (std::ostream& out, const VectorParameter<T>& param)
    {
      for (int i = 0; i < param.size(); i++) {
        out << param.value_[i];
       }
       return out;
    }

    /// @brief Set the parameter name and value, and if it is required.
    void set(const std::string& name, bool required, const std::vector<T>& value) 
    { 
      name_ = name;
      required_ = required;
      value_ = value;
    }

    /// @brief Set the parameter value from a string.
    void set(const std::string& str_value)
    {
      if (str_value == "") {
        return;
      }

      std::string error_msg = "Improper vector format '" + str_value + "' found in '" + name_ + "'." + " Vector format is: (x,y,z)";
      std::regex sep("\\(|\\)|\\,");
      auto str = std::regex_replace(str_value, sep, " ");

      if constexpr (std::is_same<T, std::string>::value) {
        std::stringstream ssin(str);
        std::string value;
        while (ssin >> value) {
          value_.push_back(value);
        }
      } else {
        T value;
        std::istringstream ssin(str);
        while (ssin >> value) {
          value_.push_back(value);
        }
      }
    }

    bool check_required_set()
    {
      if (!required_) {
        return true;
      }
      return value_set_;
    }

    std::vector<T> value_;
    std::string name_;
    bool required_ = false;
    bool value_set_ = false;
    std::vector<T> range_;
};

/// @brief struct to define a row of CANN model parameter table
struct CANNRow {
    Parameter<int> invariant_index;
    VectorParameter<int> activation_functions;  // Fixed size (3 values)
    VectorParameter<double> weights;           // Fixed size (3 values)
};

/// @brief Defines parameter name and value, and stores them in
/// maps for settng values from XML.
class ParameterLists
{
  public:

    ParameterLists() { }

    void set_xml_element_name(const std::string& name) 
    {
      xml_element_name = name;
    }

    /// @brief Set the name, default value and the parameter required flag.
    void set_parameter(const std::string& name, const bool value, bool required, Parameter<bool>& param) 
    {
      param.set(name, required, value); 
      params_map[name] = &param;
    }

    void set_parameter(const std::string& name, const double value, bool required, Parameter<double>& param)
    {
      param.set(name, required, value); 
      params_map[name] = &param;
    }

    void set_parameter(const std::string& name, std::initializer_list<double> value, bool required, VectorParameter<double>& param)
    {
      param.set(name, required, value); 
      params_map[name] = &param;
    }

    void set_parameter(const std::string& name, std::initializer_list<int> value, bool required, VectorParameter<int>& param)
    {
      param.set(name, required, value); 
      params_map[name] = &param;
    }

    void set_parameter(const std::string& name, std::initializer_list<std::string> value, bool required, 
        VectorParameter<std::string>& param)
    {
      param.set(name, required, value); 
      params_map[name] = &param;
    }

    void set_parameter(const std::string& name, const int value, bool required, Parameter<int>& param, std::vector<int> range = {}) 
    {
      param.set(name, required, value); 
      params_map[name] = &param;
    }

    void set_parameter(const std::string& name, const std::string& value, bool required, Parameter<std::string>& param) 
    {
      param.set(name, required, value); 
      params_map[name] = &param;
    }

    /// @brief set_parameter function to handle CANNRow
    void set_parameter_value_CANN(const std::string& name, const std::string& value) 
    {
      if (params_map.count(name) == 0) {
        throw std::runtime_error("Unknown " + xml_element_name + " XML element '" + name + "'.");
      }

      auto& param_variant = params_map[name];

      // Check for Activation_functions
      if (name == "Activation_functions") {
        if (auto* vec_param = std::get_if<VectorParameter<int>*>(&param_variant)) {
          (*vec_param)->value_.clear();  // Clear the vector before setting
          (*vec_param)->set(value);  // Set the new value
        } else {
          throw std::runtime_error("Activation_functions is not a VectorParameter<int>.");
        }
      }
      // Check for Weights
      else if (name == "Weights") {
        if (auto* vec_param = std::get_if<VectorParameter<double>*>(&param_variant)) {
          (*vec_param)->value_.clear();  // Clear the vector before setting
          (*vec_param)->set(value);  // Set the new value
        } else {
          throw std::runtime_error("Weights is not a VectorParameter<double>.");
        }
      }
      // Default: everything else
      else {
        std::visit([&](auto&& p) -> void {
          p->set(value);
        }, param_variant);
      }
    }



    /// @brief Set the value of a paramter from a string.
    void set_parameter_value(const std::string& name, const std::string& value) 
    {
      if (params_map.count(name) == 0) {
        throw std::runtime_error("Unknown " + xml_element_name + " XML element '" + name + "'.");
      }

      std::visit([value](auto&& p) { p->set(value); }, params_map[name]);
    }

    /// @brief Check if any required parameters have not been set.
    void check_required()
    {
      bool unset_found = false;

      for (auto& [ key, param ] : params_map) {
        if (std::visit([](auto&& p) {
          return !p->check_required_set();
        }, param)) { 
          throw std::runtime_error(xml_element_name + " XML element '" + key + "' has not been set.");
        }
      }
    }

    /// @brief Get the defined parameters as a map of strings.
    std::map<std::string,std::string> get_parameter_list()
    {
      std::map<std::string,std::string> params;

      for (auto& [ key, param ] : params_map) {
        std::visit([&params](auto&& p) {
          params[p->name()] = p->svalue();
        }, param); 
      }

      return params;
    }

    /// @brief Print the parameters.
    void print_parameter_list()
    {
      for (auto& [ key, param ] : params_map) {
        std::cout << key << ": ";
        std::visit([](auto& p) {
          std::cout << p->name_ << std::endl;
          std::cout << p->svalue() << std::endl;
        }, param);
      }
    }

    /// @brief Map used for storing parameters by name / Parameter template union.
    std::map<std::string, std::variant<Parameter<bool>*, Parameter<double>*, Parameter<int>*, 
        Parameter<std::string>*, VectorParameter<double>*, VectorParameter<int>*,
        VectorParameter<std::string>* >> params_map;

    std::string xml_element_name = "";
};

//----------------------
// IncludeParameterFile
//----------------------
// The IncludeParameterFile class is used to read and set the 
// root element of external XML file.
//
class IncludeParametersFile 
{
  public:
    IncludeParametersFile(const char* file_name);
    tinyxml2::XMLDocument document;
    tinyxml2::XMLElement* root_element = nullptr;
    static std::string NAME;
};

//////////////////////////////////////////////////////////
//            ConstitutiveModelParameters               //
//////////////////////////////////////////////////////////

// The following classes are used to store parameters for
// various constitutive models.

class LeeSacksParameters : public ParameterLists
{
  public:
    LeeSacksParameters();
    bool defined() const { return value_set; };
    void set_values(tinyxml2::XMLElement* con_model_params);
    void print_parameters();
    Parameter<double> a;
    Parameter<double> a0;
    Parameter<double> b1;
    Parameter<double> b2;
    Parameter<double> mu0;
    bool value_set = false;
};

class GuccioneParameters : public ParameterLists
{
  public:
    GuccioneParameters();
    bool defined() const { return value_set; };
    void set_values(tinyxml2::XMLElement* con_model_params);
    void print_parameters();
    Parameter<double> bf;
    Parameter<double> bfs;
    Parameter<double> bt;
    Parameter<double> c;
    bool value_set = false;
};

//---------------------
// HolzapfelParameters
//---------------------
class HolzapfelParameters : public ParameterLists
{
  public:
    HolzapfelParameters();
    bool defined() const { return value_set; };
    void set_values(tinyxml2::XMLElement* con_model_params);
    void print_parameters();

    Parameter<double> a;
    Parameter<double> b;
    Parameter<double> a4f;
    Parameter<double> b4f;
    Parameter<double> a4s;
    Parameter<double> b4s;
    Parameter<double> afs;
    Parameter<double> bfs;
    Parameter<double> k;

    bool value_set = false;
};

class HolzapfelGasserOgdenParameters : public ParameterLists
{ 
  public:
    HolzapfelGasserOgdenParameters();
    bool defined() const { return value_set; };
    void set_values(tinyxml2::XMLElement* con_model_params);
    void print_parameters();
    Parameter<double> a4;
    Parameter<double> b4;
    Parameter<double> a6;
    Parameter<double> b6;
    Parameter<double> kappa;
    bool value_set = false;
};

class MooneyRivlinParameters : public ParameterLists
{ 
  public:
    MooneyRivlinParameters();
    bool defined() const { return value_set; };
    void set_values(tinyxml2::XMLElement* con_model_params);
    void print_parameters();
    Parameter<double> c1;
    Parameter<double> c2;
    bool value_set = false;
};

class NeoHookeanParameters : public ParameterLists
{
  public:
    NeoHookeanParameters();
    void set_values(tinyxml2::XMLElement* modl_params);
    void print_parameters();
    bool value_set = false;
};

class StVenantKirchhoffParameters : public ParameterLists
{ 
  public:
    StVenantKirchhoffParameters();
    void set_values(tinyxml2::XMLElement* modl_params);
    void print_parameters();
    bool value_set = false;
};

/// @brief The CANNRowParameters class is used to store the parameters for
/// each row of the CANN table for the xml element "Add_row"
class CANNRowParameters : public ParameterLists
{
  public:
    CANNRowParameters();

    void print_parameters();
    void set_values(tinyxml2::XMLElement* xml_elem);

    static const std::string xml_element_name_;

    Parameter<std::string> row_name; // to identify each row
    CANNRow row;

};

/// @brief The CANNParameters class stores the parameters table rows
/// for xml element "Constitutive_model type=CANN". Each row is handled
/// with xml element "Add_row"
///
/// \code {.xml}
/// <Constitutive_model type="CANN">
/// <Add_row row_name="1">
///   <Invariant_num> 1 </Invariant_num>
///   <Activation_functions> (1,1,1) </Activation_functions>
///   <Weights> (1.0,1.0,1.0) </Weights>
///  </Add_row>
///  </Constitutive_model>
/// \endcode
class CANNParameters : public ParameterLists
{ 
  public:
    CANNParameters();
    ~CANNParameters();
    bool defined() const { return value_set; };
    void set_values(tinyxml2::XMLElement* con_model_params);
    void print_parameters();
    
    std::vector<CANNRowParameters*> rows;  // Store multiple rows

    bool value_set = false;
};

/// @brief The ConstitutiveModelParameters class store parameters
/// for various constitutive models.
class ConstitutiveModelParameters : public ParameterLists
{
  public:
    ConstitutiveModelParameters();
    void print_parameters();
    void check_constitutive_model(const Parameter<std::string>& eq_type);
    bool defined() const { return value_set; };
    void set_values(tinyxml2::XMLElement* modl_params);
    static const std::string xml_element_name_;

    // Model types supported.
    static const std::string GUCCIONE_MODEL;
    static const std::string HGO_MODEL;
    static const std::string HOLZAPFEL_OGDEN_MODEL;
    static const std::string HOLZAPFEL_OGDEN_MA_MODEL;
    static const std::string LEE_SACKS;
    static const std::string NEOHOOKEAN_MODEL;
    static const std::string STVENANT_KIRCHHOFF_MODEL;
    static const std::string CANN_MODEL;
    static const std::map<std::string, std::string> constitutive_model_types;

    // Constitutive model type.
    Parameter<std::string> type;

    GuccioneParameters guccione;
    HolzapfelParameters holzapfel;
    HolzapfelGasserOgdenParameters holzapfel_gasser_ogden;
    LeeSacksParameters lee_sacks;
    MooneyRivlinParameters mooney_rivlin;
    NeoHookeanParameters neo_hookean;
    StVenantKirchhoffParameters stvenant_kirchhoff;
    CANNParameters cann;

    bool value_set = false;
};

/// @brief Couple to reduced-order models.
class CoupleCplBCParameters : public ParameterLists
{
  public:
    CoupleCplBCParameters();

    static const std::string xml_element_name_;

    bool defined() const { return value_set; };
    void set_values(tinyxml2::XMLElement* xml_elem);
    void print_parameters();

    // attribute.
    Parameter<std::string> type;

    Parameter<std::string> file_name_for_0D_3D_communication;
    Parameter<std::string> file_name_for_saving_unknowns;
    Parameter<int> number_of_unknowns;
    Parameter<int> number_of_user_defined_outputs;
    Parameter<std::string> unknowns_initialization_file_path;

    Parameter<std::string> zerod_code_file_path;

    bool value_set = false;
};

/// @brief Coupling to GenBC.
class CoupleGenBCParameters : public ParameterLists
{
  public:
    CoupleGenBCParameters();

    static const std::string xml_element_name_;

    bool defined() const { return value_set; };
    void set_values(tinyxml2::XMLElement* xml_elem);

    // attributes.
    Parameter<std::string> type;

    // String parameters.
    Parameter<std::string> zerod_code_file_path;

    bool value_set = false;
};

//----------------------------------
// svZeroDSolverInterfaceParameters
//----------------------------------
//
class svZeroDSolverInterfaceParameters : public ParameterLists
{
  public:
    svZeroDSolverInterfaceParameters();

    static const std::string xml_element_name_;

    bool defined() const { return value_set; };
    void set_values(tinyxml2::XMLElement* xml_elem);

    Parameter<std::string> configuration_file;
    Parameter<std::string> coupling_type;

    Parameter<double> initial_flows;
    Parameter<double> initial_pressures;

    Parameter<std::string> shared_library;

    bool value_set = false;
};

/// @brief Body force over a mesh using the "Add_BF" command.
///
/// \code {.xml}
/// <Add_BF mesh="msh" >
///   <Type> volumetric </Type>
///   <Time_dependence> general </Time_dependence>
///   <Temporal_and_spatial_values_file_path> bforce.dat </Temporal_and_spatial_values_file_path>
/// </Add_BF>
/// \endcode
class BodyForceParameters : public ParameterLists
{
  public:
    BodyForceParameters();
    void print_parameters();
    void set_values(tinyxml2::XMLElement* xml_elem);
    static const std::string xml_element_name_;

    // Attributes.
    Parameter<std::string> mesh_name;

    // Boolean parameters.
    Parameter<bool> ramp_function;

    // Double parameters.
    Parameter<double> value;

    // String parameters.
    Parameter<std::string> fourier_coefficients_file_path;
    Parameter<std::string> spatial_values_file_path;
    Parameter<std::string> temporal_and_spatial_values_file_path;
    Parameter<std::string> temporal_values_file_path;
    Parameter<std::string> time_dependence;
    Parameter<std::string> type;
};

/// @brief RCR values for Neumann BC type.
///
/// \code {.xml}
/// <RCR_values>
///   <Proximal_resistance> 121.0 </Proximal_resistance>
///   <Capacitance> 1.5e-5 </Capacitance>
///   <Distal_resistance> 1212.0 </Distal_resistance>
/// </RCR_values>
/// \endcode
class BoundaryConditionRCRParameters : public ParameterLists
{
  public:
    BoundaryConditionRCRParameters();

    static const std::string xml_element_name_;

    void set_values(tinyxml2::XMLElement* xml_elem);
    void print_parameters();

    Parameter<double> capacitance;
    Parameter<double> distal_pressure;
    Parameter<double> distal_resistance;
    Parameter<double> initial_pressure;
    Parameter<double> proximal_resistance;

    bool value_set = false;
};

/// @brief The BoundaryConditionParameters stores paramaters for various
/// type of boundary conditions under the Add_BC XML element.
class BoundaryConditionParameters : public ParameterLists
{
  public:
    BoundaryConditionParameters();
    void print_parameters();
    void set_values(tinyxml2::XMLElement* bc_params);
    static const std::string xml_element_name_;

    // RCR parameters sub-element.
    BoundaryConditionRCRParameters rcr;

    // Add_BC name= attribute.
    Parameter<std::string> name;

    // Add_BC XML elements.
    //
    Parameter<bool> apply_along_normal_direction;
    Parameter<std::string> bct_file_path;

    Parameter<double> damping;
    Parameter<double> distal_pressure;
    VectorParameter<int> effective_direction;
    Parameter<bool> follower_pressure_load;
    Parameter<std::string> fourier_coefficients_file_path;

    Parameter<bool> impose_flux;
    Parameter<bool> impose_on_state_variable_integral;
    Parameter<std::string> initial_displacements_file_path;

    Parameter<double> penalty_parameter;
    Parameter<double> penalty_parameter_normal;
    Parameter<double> penalty_parameter_tangential;
    Parameter<std::string> prestress_file_path;
    Parameter<std::string> profile;
    Parameter<bool> ramp_function;

    Parameter<std::string> cst_shell_bc_type;
    Parameter<std::string> spatial_profile_file_path;
    Parameter<std::string> spatial_values_file_path;
    Parameter<double> stiffness;
    Parameter<std::string> svzerod_solver_block;

    Parameter<std::string> temporal_and_spatial_values_file_path;
    Parameter<std::string> temporal_values_file_path;
    Parameter<std::string> time_dependence;
    Parameter<std::string> traction_values_file_path;
    Parameter<double> traction_multiplier;
    Parameter<std::string> type;

    Parameter<bool> undeforming_neu_face;
    Parameter<double> value;
    Parameter<bool> weakly_applied;
    Parameter<bool> zero_out_perimeter;

    Parameter<double> resistance;
};

/// @brief The OutputParameters class stores parameters for the
/// Output XML element under Add_equation.
///
/// \code {.xml}
/// <Output type="Volume_integral" >
///   <Temperature> true </Temperature>
/// </Output>
/// \endcode
class OutputParameters : public ParameterLists
{
  public:
    OutputParameters();

    static const std::string xml_element_name_;

    void print_parameters();
    void set_values(tinyxml2::XMLElement* xml_elem);
    bool get_output_value(const std::string& name);
    std::string get_alias_value(const std::string& name);

    Parameter<std::string> type;

    // List of output names.
    std::vector<Parameter<bool>> output_list;

    // List of alias output names.
    std::vector<Parameter<std::string>> alias_list;
};

/// @brief The PrecomputedSolutionParameters class stores parameters for the
/// 'Precomputed_solution' XML element used to read in the data from a precomputed solution 
/// for the simulation state
/// \code {.xml}
/// <Precomputed_solution>
///   <Project_from_face> lumen_wall </Project_from_face>
/// </Precomputed_solution>
/// \endcode

class PrecomputedSolutionParameters: public ParameterLists
{
  public:
    PrecomputedSolutionParameters();

    void set_values(tinyxml2::XMLElement* xml_elem);

    static const std::string xml_element_name_;

    Parameter<std::string> field_name;
    Parameter<std::string> file_path;
    Parameter<double> time_step;
    Parameter<bool> use_precomputed_solution;
};

/// @brief The ProjectionParameters class stores parameters for the
/// 'Add_projection' XML element used for fluid-structure interaction 
/// simulations.
/// \code {.xml}
/// <Add_projection name="wall_inner" >
///   <Project_from_face> lumen_wall </Project_from_face>
/// </Add_projection>
/// \endcode
class ProjectionParameters : public ParameterLists
{
  public:
    ProjectionParameters();

    void set_values(tinyxml2::XMLElement* xml_elem);

    static const std::string xml_element_name_;

     Parameter<std::string> name;

     Parameter<std::string> project_from_face;
     Parameter<double> projection_tolerance;
};

/// @brief The VariableWallPropsParameters class stores parameters for
/// variable wall properties for the CMM equation.
class VariableWallPropsParameters : public ParameterLists
{
  public:
    VariableWallPropsParameters();
    static const std::string xml_element_name_;
    bool defined() const { return value_set; };
    void set_values(tinyxml2::XMLElement* xml_elemnt);

    Parameter<std::string> mesh_name;
    Parameter<std::string> wall_properties_file_path;
    bool value_set = false;
};


//////////////////////////////////////////////////////////
//                 FluidViscosity                            //
//////////////////////////////////////////////////////////

// The following classes are used to store parameters for
// various fluid viscosity models.

class FluidViscosityNewtonianParameters : public ParameterLists
{
  public:
    FluidViscosityNewtonianParameters();
    void print_parameters();
    void set_values(tinyxml2::XMLElement* equation_params);
    Parameter<double> constant_value;
};

class FluidViscosityCarreauYasudaParameters : public ParameterLists
{
  public:
    FluidViscosityCarreauYasudaParameters();
    void print_parameters();
    void set_values(tinyxml2::XMLElement* xml_elem);

    Parameter<double> limiting_high_shear_rate_viscosity;
    Parameter<double> limiting_low_shear_rate_viscosity;
    Parameter<double> power_law_index;
    Parameter<double> shear_rate_tensor_multipler;
    Parameter<double> shear_rate_tensor_exponent;
};

class FluidViscosityCassonsParameters : public ParameterLists
{
  public:
    FluidViscosityCassonsParameters();
    void print_parameters();
    void set_values(tinyxml2::XMLElement* xml_elem);
    Parameter<double> asymptotic_viscosity; 
    Parameter<double> yield_stress;
    Parameter<double> low_shear_rate_threshold;
};

class FluidViscosityParameters : public ParameterLists
{
  public:
    FluidViscosityParameters();

    static const std::string xml_element_name_;

    static const std::string CONSTANT_MODEL;
    static const std::string CARREAU_YASUDA_MODEL;
    static const std::string CASSONS_MODEL;
    static const std::set<std::string> model_names;

    void print_parameters();
    void set_values(tinyxml2::XMLElement* xml_elem);

    Parameter<std::string> model;

    FluidViscosityNewtonianParameters newtonian_model;
    FluidViscosityCarreauYasudaParameters carreau_yasuda_model;
    FluidViscosityCassonsParameters cassons_model;
};

//////////////////////////////////////////////////////////
//                 SolidViscosity                            //
//////////////////////////////////////////////////////////

// The following classes are used to store parameters for
// various solid viscosity models.

class SolidViscosityNewtonianParameters : public ParameterLists
{
  public:
    SolidViscosityNewtonianParameters();
    void print_parameters();
    void set_values(tinyxml2::XMLElement* equation_params);
    Parameter<double> constant_value;
};

class SolidViscosityPotentialParameters : public ParameterLists
{
  public:
    SolidViscosityPotentialParameters();
    void print_parameters();
    void set_values(tinyxml2::XMLElement* equation_params);
    Parameter<double> constant_value;
};

class SolidViscosityParameters : public ParameterLists
{
  public:
    SolidViscosityParameters();

    static const std::string xml_element_name_;

    static const std::string NEWTONIAN_MODEL;
    static const std::string POTENTIAL_MODEL;
    static const std::set<std::string> model_names;

    void print_parameters();
    void set_values(tinyxml2::XMLElement* xml_elem);

    Parameter<std::string> model;

    SolidViscosityNewtonianParameters newtonian_model;
    SolidViscosityPotentialParameters potential_model;
};


/// @brief The LinearAlgebraParameters class stores parameters for
/// the 'Linear_algebra' XML element.
class LinearAlgebraParameters : public ParameterLists
{
  public:
    static const std::string xml_element_name_;
    LinearAlgebraParameters();
    void check_input_parameters();
    void print_parameters();
    void set_values(tinyxml2::XMLElement* fsi_file);
    bool defined() const { return values_set_; };

    bool values_set_ = false;
    Parameter<std::string> type;

    Parameter<std::string> assembly;
    Parameter<std::string> configuration_file;
    Parameter<std::string> preconditioner;
};

/// @brief The LinearSolverParameters class stores parameters for
/// the 'LS' XML element.
class LinearSolverParameters : public ParameterLists
{
  public:
    LinearSolverParameters();

    void print_parameters();
    void set_values(tinyxml2::XMLElement* fsi_file);

    static const std::string xml_element_name_;

    Parameter<std::string> type;

    Parameter<double> absolute_tolerance;
    Parameter<int> krylov_space_dimension;

    Parameter<int> max_iterations;
    Parameter<int> ns_cg_max_iterations;
    Parameter<double> ns_cg_tolerance;
    Parameter<int> ns_gm_max_iterations; 
    Parameter<double> ns_gm_tolerance;

    //Parameter<std::string> preconditioner;

    Parameter<double> tolerance;

    LinearAlgebraParameters linear_algebra;
};

/// @brief The StimulusParameters class stores parameters for 
/// 'Stimulus' XML element used to parameters for 
/// pacemaker cells.
///
/// \code {.xml}
/// <Stimulus type="Istim" >
///   <Amplitude> -52.0 </Amplitude>
///   <Start_time> 0.0 </Start_time>
///   <Duration> 1.0 </Duration>
///   <Cycle_length> 10000.0 </Cycle_length>
/// </Stimulus>
/// \endcode
class StimulusParameters : public ParameterLists
{ 
  public:
    StimulusParameters();

    static const std::string xml_element_name_;
    
    bool defined() const { return value_set; };
    void print_parameters();
    void set_values(tinyxml2::XMLElement* xml_elem);
    
    Parameter<std::string> type;
    
    Parameter<double> amplitude;
    Parameter<double> cycle_length;
    Parameter<double> duration;
    Parameter<double> start_time;
    
    bool value_set = false;
};

class ECGLeadsParameters : public ParameterLists
{ 
  public:
    ECGLeadsParameters();

    static const std::string xml_element_name_;
    
    bool defined() const { return value_set; };
    void print_parameters();
    void set_values(tinyxml2::XMLElement* xml_elem);
    
    Parameter<std::string> x_coords_file_path;
    Parameter<std::string> y_coords_file_path;
    Parameter<std::string> z_coords_file_path;
    
    bool value_set = false;
};

/// @brief The DirectionalDistributionParameters class stores directional
/// distribution parameters for active stress.
///
/// \code {.xml}
/// <Directional_distribution>
///   <Fiber_direction> 1.0 </Fiber_direction>
///   <Sheet_direction> 0.0 </Sheet_direction>
///   <Sheet_normal_direction> 0.0 </Sheet_normal_direction>
/// </Directional_distribution>
/// \endcode
class DirectionalDistributionParameters : public ParameterLists
{
  public:
    DirectionalDistributionParameters();

    static const std::string xml_element_name_;

    bool defined() const { return value_set; };
    void print_parameters();
    void set_values(tinyxml2::XMLElement* xml_elem);
    void validate() const;  // Validate directional fractions

    Parameter<double> fiber_direction;
    Parameter<double> sheet_direction;
    Parameter<double> sheet_normal_direction;

    bool value_set = false;
};

/// @brief The FiberReinforcementStressParameters class stores fiber
/// reinforcement stress parameters for the 'Fiber_reinforcement_stress` 
/// XML element.
///
/// \code {.xml}
/// <Fiber_reinforcement_stress type="Unsteady" >
///   <Temporal_values_file_path> fib_stress.dat </Temporal_values_file_path>
///   <Ramp_function> true </Ramp_function>
///   <Directional_distribution>
///     <Fiber_direction> 0.7 </Fiber_direction>
///     <Sheet_direction> 0.2 </Sheet_direction>
///     <Sheet_normal_direction> 0.1 </Sheet_normal_direction>
///   </Directional_distribution>
/// </Fiber_reinforcement_stress>
/// \endcode
class FiberReinforcementStressParameters : public ParameterLists
{
  public:
    FiberReinforcementStressParameters();

    static const std::string xml_element_name_;

    bool defined() const { return value_set; };
    void print_parameters();
    void set_values(tinyxml2::XMLElement* xml_elem);

    Parameter<std::string> type;

    Parameter<bool> ramp_function;
    Parameter<std::string> temporal_values_file_path;
    Parameter<double> value;

    // Directional stress distribution parameters
    DirectionalDistributionParameters directional_distribution;

    bool value_set = false;
};

/// @brief The DomainParameters class stores parameters for the XML
/// 'Domain' element to specify properties for solving equations.
///
/// \code {.xml}
/// <Domain id="1" >
///   <Equation> fluid </Equation>
///   <Density> 1.06 </Density>
///   <Viscosity model="Constant" >
///     <Value> 0.04 </Value>
///   </Viscosity>
///   <Backflow_stabilization_coefficient> 0.2 </Backflow_stabilization_coefficient>
/// </Domain>
/// \endcode
class DomainParameters : public ParameterLists
{
  public:
    DomainParameters();

    static const std::string xml_element_name_;

    void print_parameters();
    void set_values(tinyxml2::XMLElement* xml_elem, bool from_external_xml = false);

    // Parameters for sub-elements under the Domain element.
    ConstitutiveModelParameters constitutive_model;
    FiberReinforcementStressParameters fiber_reinforcement_stress;
    StimulusParameters stimulus;
    FluidViscosityParameters fluid_viscosity;
    SolidViscosityParameters solid_viscosity;

    // Attributes.
    Parameter<std::string> id;

    Parameter<double> absolute_tolerance;
    VectorParameter<double> anisotropic_conductivity;
    Parameter<double> backflow_stabilization_coefficient;

    Parameter<double> conductivity;
    //Parameter<std::string> constitutive_model_name;
    Parameter<double> continuity_stabilization_coefficient;

    Parameter<double> density;
    Parameter<std::string> dilational_penalty_model;

    Parameter<std::string> equation;
    Parameter<double> elasticity_modulus;
    Parameter<std::string> electrophysiology_model;

    Parameter<double> feedback_parameter_for_stretch_activated_currents;
    Parameter<double> fluid_density;
    Parameter<double> force_x;
    Parameter<double> force_y;
    Parameter<double> force_z;

    Parameter<std::string> include_xml;
    Parameter<double> isotropic_conductivity;

    Parameter<double> mass_damping;
    Parameter<int> maximum_iterations;
    Parameter<double> momentum_stabilization_coefficient;
    Parameter<std::string> myocardial_zone;

    Parameter<double> G_Na;
    Parameter<double> G_CaL;
    Parameter<double> G_Kr;
    Parameter<double> G_Ks;
    Parameter<double> G_to;

    Parameter<double> tau_fi;
    Parameter<double> tau_si;

    Parameter<std::string> ode_solver;
    Parameter<double> penalty_parameter;
    Parameter<double> poisson_ratio;
    Parameter<double> relative_tolerance;

    Parameter<double> shell_thickness;
    Parameter<double> solid_density;
    Parameter<double> source_term;
    Parameter<double> time_step_for_integration;
    
    // Inverse of Darcy permeability. Default value of 0.0 for Navier-Stokes and non-zero for Navier-Stokes-Brinkman
    Parameter<double> inverse_darcy_permeability;
};

/// @brief The RemesherParameters class stores parameters for the 
/// 'Remesher' XML element used for remeshing.
///
/// \code {.xml}
/// <Remesher type="Tetgen" >
///   <Max_edge_size name="lumen" value="0.7"> </Max_edge_size>
///   <Max_edge_size name="wall"  value="0.5"> </Max_edge_size>
///   <Min_dihedral_angle> 10.0 </Min_dihedral_angle>
///   <Max_radius_ratio> 1.1 </Max_radius_ratio>
///   <Remesh_frequency> 1000 </Remesh_frequency>
///   <Frequency_for_copying_data> 1 </Frequency_for_copying_data>
/// </Remesher>
/// \endcode
class RemesherParameters : public ParameterLists
{
  public:
    RemesherParameters();

    static const std::string xml_element_name_;
    bool values_set_ = false;

    bool defined() const { return values_set_; };
    void print_parameters();
    double get_edge_size(const std::string& name) const { return max_edge_sizes_.at(name); } 
    bool has_edge_size(const std::string& name) const { return max_edge_sizes_.count(name) == 1; }
    void set_values(tinyxml2::XMLElement* mesh_elem);

    // Values given in the 'Max_edge_size' element.
    std::map<std::string, double> max_edge_sizes_;

    Parameter<std::string> type;
    Parameter<double> min_dihedral_angle;
    Parameter<double> max_radius_ratio; 
    Parameter<int> remesh_frequency;
    Parameter<int> frequency_for_copying_data;
};

/// @brief The ContactParameters class stores parameters for the 'Contact''
/// XML element used to specify parameter values for contact
/// computations. 
class ContactParameters : public ParameterLists
{
  public:
    ContactParameters();

    static const std::string xml_element_name_;

    void print_parameters();
    void set_values(tinyxml2::XMLElement* xml_elem);

    Parameter<double> closest_gap_to_activate_penalty;

    Parameter<double> desired_separation;

    Parameter<double> min_norm_of_face_normals;

    Parameter<std::string> model;

    Parameter<double> penalty_constant;
};

/// @brief The EquationParameters class stores parameters for the 'Add_equation'
/// XML element used to specify an equation to be solved (e.g. fluid).
///
/// \code {.xml}
/// <Add_equation type="FSI" >
///   <Coupled> true </Coupled>
///   <Min_iterations> 1 </Min_iterations>
///   <Max_iterations> 1 </Max_iterations>
///   .
///   .
///   .
/// </Add_equation>
/// \endcode
class EquationParameters : public ParameterLists
{
  public:
    EquationParameters();

    static const std::string xml_element_name_;

    void print_parameters();
    void set_values(tinyxml2::XMLElement* xml_elem, DomainParameters* default_domain=nullptr);

    Parameter<double> backflow_stabilization_coefficient;

    Parameter<double> conductivity;
    Parameter<double> continuity_stabilization_coefficient;
    Parameter<bool> coupled;

    Parameter<double> density;
    Parameter<std::string> dilational_penalty_model;

    Parameter<double> elasticity_modulus;

    Parameter<std::string> include_xml;
    Parameter<std::string> initialize;
    Parameter<bool> initialize_rcr_from_flow;

    Parameter<int> max_iterations;
    Parameter<int> min_iterations;
    Parameter<double> momentum_stabilization_coefficient;

    Parameter<double> penalty_parameter;
    Parameter<double> poisson_ratio;
    Parameter<bool> prestress;

    Parameter<double> source_term;
    Parameter<double> tolerance;

    Parameter<std::string> type;
    Parameter<bool> use_taylor_hood_type_basis;

    // Inverse of Darcy permeability. Default value of 0.0 for Navier-Stokes and non-zero for Navier-Stokes-Brinkman
    Parameter<double> inverse_darcy_permeability;

    // Sub-element parameters.
    //
    std::vector<BodyForceParameters*> body_forces;

    std::vector<BoundaryConditionParameters*> boundary_conditions;

    CoupleCplBCParameters couple_to_cplBC;
    CoupleGenBCParameters couple_to_genBC;

    svZeroDSolverInterfaceParameters svzerodsolver_interface_parameters;

    DomainParameters* default_domain = nullptr;

    std::vector<DomainParameters*> domains;

    LinearSolverParameters linear_solver;

    std::vector<OutputParameters*> outputs;

    RemesherParameters remesher;

    VariableWallPropsParameters variable_wall_properties;

    FluidViscosityParameters fluid_viscosity;

    SolidViscosityParameters solid_viscosity;

    ECGLeadsParameters ecg_leads;

};

/// @brief The GeneralSimulationParameters class stores paramaters for the
/// 'GeneralSimulationParameters' XML element.
///
/// \code {.xml}
/// <GeneralSimulationParameters>
///   <Continue_previous_simulation> 0 </Continue_previous_simulation>
///   <Number_of_spatial_dimensions> 3 </Number_of_spatial_dimensions>
///   <Number_of_time_steps> 1 </Number_of_time_steps>
///   <Time_step_size> 1e-4 </Time_step_size>
///   <Spectral_radius_of_infinite_time_step> 0.50 </Spectral_radius_of_infinite_time_step>
///   <Searched_file_name_to_trigger_stop> STOP_SIM </Searched_file_name_to_trigger_stop>
///   <Save_results_to_VTK_format> true </Save_results_to_VTK_format>
///   <Name_prefix_of_saved_VTK_files> result </Name_prefix_of_saved_VTK_files>
///   <Increment_in_saving_VTK_files> 1 </Increment_in_saving_VTK_files>
///   <Start_saving_after_time_step> 1 </Start_saving_after_time_step>
///   <Increment_in_saving_restart_files> 1 </Increment_in_saving_restart_files>
///   <Convert_BIN_to_VTK_format> 0 </Convert_BIN_to_VTK_format>
///   <Verbose> 1 </Verbose>
///   <Warning> 0 </Warning>
///   <Debug> 0 </Debug>
///   <Simulation_requires_remeshing> true </Simulation_requires_remeshing>
/// </GeneralSimulationParameters>
/// \endcode
class GeneralSimulationParameters : public ParameterLists 
{
  public:
    GeneralSimulationParameters();

    void print_parameters();
    void set_values(tinyxml2::XMLElement* xml_element, bool from_external_xml = false);

    std::string xml_element_name;

    Parameter<bool> check_ien_order;
    Parameter<bool> continue_previous_simulation;
    Parameter<bool> convert_bin_to_vtk_format;
    Parameter<bool> debug;
    Parameter<bool> overwrite_restart_file;
    Parameter<bool> save_averaged_results;
    Parameter<bool> save_results_to_vtk_format;
    Parameter<bool> simulation_requires_remeshing;
    Parameter<bool> start_averaging_from_zero;
    Parameter<bool> verbose;
    Parameter<bool> warning;

    Parameter<double> spectral_radius_of_infinite_time_step;
    Parameter<double> time_step_size;

    Parameter<std::string> include_xml;
    Parameter<int> increment_in_saving_restart_files;
    Parameter<int> increment_in_saving_vtk_files;
    Parameter<int> number_of_spatial_dimensions;
    Parameter<int> number_of_initialization_time_steps;
    Parameter<int> start_saving_after_time_step;
    Parameter<int> starting_time_step;
    Parameter<int> number_of_time_steps;

    Parameter<std::string> name_prefix_of_saved_vtk_files;
    Parameter<std::string> restart_file_name; 
    Parameter<std::string> searched_file_name_to_trigger_stop; 
    Parameter<std::string> save_results_in_folder; 
    Parameter<std::string> simulation_initialization_file_path;
};

/// @brief The FaceParameters class is used to store parameters for the
/// 'Add_face' XML element.
class FaceParameters : public ParameterLists
{
  public:
    FaceParameters();

    void print_parameters();
    void set_values(tinyxml2::XMLElement* xml_elem);

    static const std::string xml_element_name_;

    Parameter<std::string> end_nodes_face_file_path;
    Parameter<std::string> face_file_path;
    Parameter<std::string> name;

    Parameter<double> quadrature_modifier_TRI3;
};

/// @brief The MeshParameters class is used to store paramaters for the
/// 'Add_mesh' XML element.
///
/// \code {.xml}
/// <Add_mesh name="lumen" >
///   <Mesh_file_path> mesh/lumen/mesh-complete.mesh.vtu  </Mesh_file_path>
///
///   <Add_face name="lumen_inlet">
///       <Face_file_path> mesh/lumen/mesh-surfaces/lumen_inlet.vtp </Face_file_path>
///   </Add_face>
///
///   <Add_face name="lumen_outlet">
///       <Face_file_path> mesh/lumen/mesh-surfaces/lumen_outlet.vtp </Face_file_path>
///   </Add_face>
///
///   <Add_face name="lumen_wall">
///       <Face_file_path> mesh/lumen/mesh-surfaces/lumen_wall.vtp </Face_file_path>
///   </Add_face>
///
///   <Domain> 0 </Domain>
/// 
/// </Add_mesh>
/// \endcode
class MeshParameters : public ParameterLists
{
  public:
    MeshParameters();

    static const std::string xml_element_name_;

    void print_parameters();
    void set_values(tinyxml2::XMLElement* mesh_elem, bool from_external_xml = false);
    std::string get_name() const { return name.value(); };
    std::string get_path() const { return mesh_file_path.value(); };

    std::vector<FaceParameters*> face_parameters;

    // Add_mesh name= 
    Parameter<std::string> name;

    // Parameters under Add_mesh 
    //
    Parameter<int> domain_id;
    Parameter<std::string> domain_file_path;

    VectorParameter<std::string> fiber_direction_file_paths;
    //Parameter<std::string> fiber_direction_file_path;
    std::vector<VectorParameter<double>> fiber_directions;
    //VectorParameter<double> fiber_direction;

    Parameter<std::string> include_xml;
    Parameter<std::string> initial_displacements_file_path;
    Parameter<std::string> initial_pressures_file_path;
    Parameter<bool> initialize_rcr_from_flow;
    Parameter<std::string> initial_velocities_file_path;

    Parameter<std::string> mesh_file_path;
    Parameter<double> mesh_scale_factor;
    Parameter<std::string> prestress_file_path;

    Parameter<bool> set_mesh_as_fibers;
    Parameter<bool> set_mesh_as_shell;

    Parameter<double> quadrature_modifier_TET4;
};

//////////////////////////////////////////////////////////
//         Resistive immersed surfaces method           //
//////////////////////////////////////////////////////////

/// @brief The RISProjectionParameters class stores parameters for the
/// 'Add_RIS_projection' XML element used for RIS valve simulations.
/// \code {.xml}
/// <Add_RIS_projection name="left_ris" >
///   <Project_from_face> right_ris </Project_from_face>
///   <Resistance> 1.e6 </Resistance>
/// </Add_RIS_projection>
/// \endcode
class RISProjectionParameters : public ParameterLists
{
  public:
    RISProjectionParameters();

    void set_values(tinyxml2::XMLElement* xml_elem);

    static const std::string xml_element_name_;

    Parameter<std::string> name;

    Parameter<std::string> project_from_face;
    Parameter<double> resistance;
    Parameter<double> projection_tolerance;
};

/// @brief The URISFaceParameters class is used to store parameters for the
/// 'Add_URIS_face' XML element.
class URISFaceParameters : public ParameterLists
{
  public:
    URISFaceParameters();

    void print_parameters();
    void set_values(tinyxml2::XMLElement* xml_elem);

    static const std::string xml_element_name_;

    Parameter<std::string> name; // Name of the valve surface

    Parameter<std::string> face_file_path; // File path for the valve surface
    Parameter<std::string> open_motion_file_path; // File path for the open motion of the valve
    Parameter<std::string> close_motion_file_path; // File path for the close motion of the valve

};

/// @brief The URISMeshParameters class is used to store paramaters for the
/// 'Add_URIS_mesh' XML element.
///
/// \code {.xml}
/// <Add_uris_mesh name="MV" > 
///   <Add_uris_face name="LCC" > 
///     <Face_file_path> meshes/uris_face.vtu </Face_file_path>
///     <Open_motion_file_path> meshes/uris_facemotion_open.dat </Open_motion_file_path>
///     <Close_motion_file_path> meshes/uris_facemotion_close.dat </Close_motion_file_path>
///   </Add_uris_face>
///   <Mesh_scale_factor> 1.0 </Mesh_scale_factor>
///   <Thickness> 0.25 </Thickness>
///   <Resistance> 1.0e5 </Resistance>
///   <Positive_flow_normal_file> meshes/normal.dat </Positive_flow_normal_file>
/// </Add_uris_mesh>
/// \endcode
class URISMeshParameters : public ParameterLists
{
  public:
    URISMeshParameters();

    static const std::string xml_element_name_;

    void print_parameters();
    void set_values(tinyxml2::XMLElement* mesh_elem);
    std::string get_name() const { return name.value(); };
    // std::string get_path() const { return mesh_file_path.value(); };

    std::vector<URISFaceParameters*> URIS_face_parameters;

    // Add_mesh name
    Parameter<std::string> name; // Name of the valve mesh

    // Parameters under Add_URIS_mesh 
    Parameter<double> mesh_scale_factor; // Scale factor for the mesh
    Parameter<double> thickness; // Thickness of the valve
    Parameter<double> close_thickness; // Thickness of the valve when it is closed
    Parameter<double> resistance; // Resistance of the valve
    Parameter<double> resistance_close; // Resistance of the valve when it is closed
    Parameter<bool> valve_starts_as_closed; // Whether the valve starts as closed
    Parameter<std::string> positive_flow_normal_file_path; // File path for the positive flow normal

};



/// @brief The Parameters class stores parameter values read in from a solver input file.
class Parameters {

  public:
    Parameters();

    static const std::set<std::string> constitutive_model_names;
    static const std::set<std::string> equation_names;
    static const std::string FSI_FILE;

    void get_logging_levels(int& verbose, int& warning, int& debug);
    void print_parameters();
    void read_xml(std::string file_name);

    void set_contact_values(tinyxml2::XMLElement* root_element);
    void set_equation_values(tinyxml2::XMLElement* root_element);
    void set_mesh_values(tinyxml2::XMLElement* root_element);
    void set_precomputed_solution_values(tinyxml2::XMLElement* root_element);
    void set_projection_values(tinyxml2::XMLElement* root_element);
    void set_svzerodsolver_interface_values(tinyxml2::XMLElement* root_element);

    void set_RIS_projection_values(tinyxml2::XMLElement* root_element);
    void set_URIS_mesh_values(tinyxml2::XMLElement* root_element);

    // Objects representing each parameter section of XML file.
    ContactParameters contact_parameters;
    GeneralSimulationParameters general_simulation_parameters;
    std::vector<MeshParameters*> mesh_parameters;
    std::vector<EquationParameters*> equation_parameters;
    std::vector<ProjectionParameters*> projection_parameters;
    PrecomputedSolutionParameters precomputed_solution_parameters;

    std::vector<RISProjectionParameters*> RIS_projection_parameters;
    std::vector<URISMeshParameters*> URIS_mesh_parameters;

};

#endif

