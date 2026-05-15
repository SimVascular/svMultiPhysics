// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// The classes defined here duplicate the data structures in the Fortran CEPMOD module
// defined in CEPMOD.f. 

// This module defines data structures for cardiac electrophysiology
// model equation. It also interfaces with individual modules for
// the cellular activation model.

#ifndef CEP_MOD_H
#define CEP_MOD_H

#include "consts.h"
#include "ionic_model.h"

#include "Array.h"
#include "Vector.h"
#include <map>
#include <memory>

/// @brief Type of cardiac electrophysiology models.
enum class ElectrophysiologyModelType {
  NA = 100, 
  AP = 101,
  BO = 102, 
  FN = 103, 
  TTP = 104
};

extern const std::map<ElectrophysiologyModelType, std::string> cep_model_type_to_name;
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

    /// @brief Ionic model instance.
    std::shared_ptr<IonicModel> ionic_model;
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

    /// @brief ECG leads
    ecgLeadsType ecgleads;
};

#endif

