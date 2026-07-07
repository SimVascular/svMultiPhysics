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

class ComMod;
class CmMod;
class cmType;
class StimulusParameters;

/// @brief External stimulus type
class stimType
{
  public:
    /// @brief Spatial bounds for a CEP stimulus region.
    class SpatialBounds
    {
      public:
        /// @brief Set box bounds.
        void set_box(const Vector<double>& min, const Vector<double>& max);

        /// @brief Set sphere bounds.
        void set_sphere(const Vector<double>& center, const double radius);

        /// @brief Return true if the point lies inside all active spatial bounds.
        bool contains(const Vector<double>& x) const;

        /// @brief Broadcast spatial bounds to all MPI ranks.
        void distribute(const CmMod& cm_mod, const cmType& cm);

      private:
        /// @brief True if a box region has been set.
        bool has_box = false;
        /// @brief True if a sphere region has been set.
        bool has_sphere = false;

        /// @brief Minimum corner of the box region.
        Vector<double> box_min;
        /// @brief Maximum corner of the box region.
        Vector<double> box_max;
        /// @brief Center of the sphere region.
        Vector<double> sphere_center;
        /// @brief Radius of the sphere region.
        double sphere_radius = 0.0;

        /// @brief Return true if x lies inside the box. Assumes has_box is true.
        bool inside_box(const Vector<double>& x) const;
        /// @brief Return true if x lies inside the sphere. Assumes has_sphere is true.
        bool inside_sphere(const Vector<double>& x) const;
    };

    /// @brief Return the applied stimulus value at a point and time.
    double operator()(const double time, const Vector<double>& x) const;

    /// @brief Set stimulus parameters from parsed XML parameters.
    void read_parameters(const StimulusParameters& params, const int nsd, const double default_cycle_length);

    /// @brief Broadcast stimulus parameters to all ranks.
    void distribute(const CmMod& cm_mod, const cmType& cm);

  private:
    /// @brief Time at which the stimulus begins within each cycle.
    double start_time = 0.0;
    /// @brief Duration of the stimulus pulse within each cycle.
    double duration = 0.0;
    /// @brief Length of one stimulus cycle.
    double cycle_length = 0.0;
    /// @brief Amplitude of the applied stimulus.
    double amplitude = 0.0;

    /// @brief Spatial region to which the stimulus is applied.
    SpatialBounds spatial_bounds;

    /// @brief Return true if the stimulus is active at the given time.
    bool is_active(const double time) const;
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

