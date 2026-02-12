// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// The classes defined here duplicate the data structures in the Fortran CEPMOD module
// defined in CEPMOD.f. 

// This module defines data structures for cardiac electrophysiology
// model equation. It also interfaces with individual modules for
// the cellular activation model.


#ifndef CEP_MOD_H 
#define CEP_MOD_H 

#include "CepModAp.h"
#include "CepModBo.h"
#include "CepModFn.h"
#include "CepModTtp.h"
#include "consts.h"

#include "Array.h"
#include "Vector.h"
#include <map>

/// @brief Type of cardiac electrophysiology models.
enum class ElectrophysiologyModelType {
  NA = 100, 
  AP = 101,
  BO = 102, 
  FN = 103, 
  TTP = 104
};

extern const std::map<std::string,ElectrophysiologyModelType> cep_model_name_to_type;

/// @brief Print ElectrophysiologyModelType as a string.
static std::ostream &operator << ( std::ostream& strm, ElectrophysiologyModelType type)
{
  const std::map<ElectrophysiologyModelType, std::string> names = { 
    {ElectrophysiologyModelType::NA, "NA"}, 
    {ElectrophysiologyModelType::AP,"AP"}, 
    {ElectrophysiologyModelType::BO, "BO"}, 
    {ElectrophysiologyModelType::FN, "FN"}, 
    {ElectrophysiologyModelType::TTP, "TTP"}, 
  };
  return strm << names.at(type);
}

/// @brief Time integration scheme.
enum class TimeIntegratioType {
  NA = 200, 
  FE = 201,
  RK4 = 202, 
  CN2 = 203
};

extern const std::map<std::string,TimeIntegratioType> cep_time_int_to_type;

static std::ostream &operator << ( std::ostream& strm, TimeIntegratioType type)
{
  const std::map<TimeIntegratioType, std::string> names = { 
    {TimeIntegratioType::NA, "NA"}, 
    {TimeIntegratioType::FE, "FE"}, 
    {TimeIntegratioType::RK4, "RK4"}, 
    {TimeIntegratioType::CN2, "CN2"}, 
  };
  return strm << names.at(type);
}

/// @brief Time integration scheme and related parameters
class odeType {
  public:
    odeType() {};

    /// @brief Time integration method type
    TimeIntegratioType tIntType = TimeIntegratioType::NA;
    //int tIntType = tIntType_NA;

    /// @brief Max. iterations for Newton-Raphson method
    int maxItr = 5;

    /// @brief Absolute tolerance
    double absTol = 1.E-8;

    /// @brief Relative tolerance
    double relTol = 1.E-4;
};

/// @brief External stimulus type
class stimType
{
  public:
    /// @brief start time
    double Ts = 0.0;

    /// @brief duration of stimulus
    double Td = 0.0;

    /// @brief cycle length
    double CL = 0.0;

    /// @brief stimulus amplitude
    double A = 0.0;
};

/// @brief ECG leads type
class ecgLeadsType
{
  public:
    /// @brief Number of leads
    int num_leads = 0;

    /// @brief x coordinates
    Vector<double> x_coords;

    /// @brief y coordinates
    Vector<double> y_coords;

    /// @brief z coordinates
    Vector<double> z_coords;

    /// @brief Pseudo ECG over each lead
    Vector<double> pseudo_ECG;

    /// @brief Output files
    std::vector<std::string> out_files;
};

/// @brief Cardiac electrophysiology model type
class cepModelType
{
  public:
    cepModelType();
    ~cepModelType();

    /// @brief Type of cardiac electrophysiology model
    ElectrophysiologyModelType cepType = ElectrophysiologyModelType::NA;

    /// @brief Number of state variables
    int nX = 0;

    /// @brief Number of gating variables
    int nG = 0;

    /// @brief  Number of fiber directions
    int nFn = 0;

    /// @brief  Myocardium zone id, default to epicardium.
    int imyo = 1;

    /// @brief  Time step for integration
    double dt = 0.0;

    /// @brief  Constant for stretch-activated-currents
    double Ksac = 0.0;

    /// @brief  Isotropic conductivity
    double Diso = 0.0;

    /// @brief  Anisotropic conductivity
    Vector<double> Dani;

    /// @brief  External stimulus
    stimType Istim;

    /// @brief  Time integration options
    odeType odes;

    /// @brief Interface for Aliev-Panfilov cellular activation model
    CepModAp ap;

    /// @brief Interface for Bueno-Orovio cellular activation model
    CepModBo bo;

    /// @brief Interface for Fitzhugh-Nagumo cellular activation model
    CepModFn fn;

    /// @brief Interface for Tusscher-Panfilov cellular activation model
    CepModTtp ttp;
};

/// @brief Cardiac electromechanics model type
class cemModelType
{
  public:
    /// @brief  Whether electrophysiology and mechanics are coupled
    bool cpld = false;
    //bool cpld = .FALSE.

    /// @brief  Whether active stress formulation is employed
    bool aStress = false;
    //bool aStress = .FALSE.

    /// @brief  Whether active strain formulation is employed
    bool aStrain = false;
    //bool aStrain = .FALSE.

    /// @brief  Local variable integrated in time
    ///    := activation force for active stress model
    ///    := fiber stretch for active strain model
    Vector<double> Ya;
};

class CepMod 
{
  public:

    /// @brief Whether cardiac electrophysiology is solved
    bool cepEq;

    /// @brief Max. dof in cellular activation model
    int nXion = 0;

    /// @brief Unknowns stored at all nodes
    Array<double> Xion;

    /// @brief Cardiac electromechanics type
    cemModelType cem;

    /// @brief Interface for Aliev-Panfilov cellular activation model.
    CepModAp ap;

    /// @brief Interface for ABueno-Orovio cellular activation model.
    CepModBo bo;

    /// @brief Interface for Fitzhugh-Nagumo cellular activation model.
    CepModFn fn;

    /// @brief Interface for Tusscher-Panfilov cellular activation model.
    CepModTtp ttp;

    /// @brief ECG leads
    ecgLeadsType ecgleads;
};

#endif

