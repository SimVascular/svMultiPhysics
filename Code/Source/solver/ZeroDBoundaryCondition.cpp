// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "ComMod.h"
#include "SimulationLogger.h"
#include "ZeroDBoundaryCondition.h"
#include "all_fun.h"
#include "consts.h"
#include "utils.h"

#define n_debug_zerod_bc

static std::string bc_type_to_string(consts::BoundaryConditionType type)
{
    switch (type) {
        case consts::BoundaryConditionType::bType_Dir: return "Dirichlet";
        case consts::BoundaryConditionType::bType_Neu: return "Neumann";
        default: return "Unsupported";
    }
}

ZeroDBoundaryCondition::ZeroDBoundaryCondition(consts::BoundaryConditionType bc_type, const faceType& face, SimulationLogger& logger)
    : face_(&face)
    , logger_(&logger)
    , bc_type_(bc_type)
{
    if (bc_type_ != consts::BoundaryConditionType::bType_Dir &&
        bc_type_ != consts::BoundaryConditionType::bType_Neu) {
        throw std::runtime_error("[ZeroDBoundaryCondition] bc_type must be bType_Dir or bType_Neu.");
    }
    #ifdef debug_zerod_bc
    std::cout << "ZeroDBoundaryCondition: BC type set to '" << bc_type_to_string(bc_type_) << "'" << std::endl;
    #endif
}

ZeroDBoundaryCondition::ZeroDBoundaryCondition(consts::BoundaryConditionType bc_type, const std::string& cap_face_vtp_file, const faceType& face, SimulationLogger& logger)
    : cap_face_vtp_file_(cap_face_vtp_file)
    , face_(&face)
    , logger_(&logger)
    , bc_type_(bc_type)
{
    if (bc_type_ != consts::BoundaryConditionType::bType_Dir &&
        bc_type_ != consts::BoundaryConditionType::bType_Neu) {
        throw std::runtime_error("[ZeroDBoundaryCondition] bc_type must be bType_Dir or bType_Neu.");
    }
    #ifdef debug_zerod_bc
    std::cout << "ZeroDBoundaryCondition: Cap face VTP file set to '" << cap_face_vtp_file_ << "'" << std::endl;
    std::cout << "ZeroDBoundaryCondition: BC type set to '" << bc_type_to_string(bc_type_) << "'" << std::endl;
    #endif
}

void ZeroDBoundaryCondition::load_cap_face_vtp(const std::string& vtp_file_path)
{
    cap_face_vtp_file_ = vtp_file_path;
    #ifdef debug_zerod_bc
    std::cout << "ZeroDBoundaryCondition: Loading cap VTP file '" << vtp_file_path << "'" << std::endl;
    #endif
    // TODO: Load VTP file and extract geometric data
    // This would typically use VtkVtpData to read the file
}

// =========================================================================
// svZeroD block configuration
// =========================================================================

void ZeroDBoundaryCondition::set_block_name(const std::string& block_name)
{
    block_name_ = block_name;
    #ifdef debug_zerod_bc
    std::cout << "ZeroDBoundaryCondition: Block name set to '" << block_name_ << "'" << std::endl;
    #endif
}

const std::string& ZeroDBoundaryCondition::get_block_name() const
{
    return block_name_;
}

void ZeroDBoundaryCondition::set_face_name(const std::string& face_name)
{
    face_name_ = face_name;
}

const std::string& ZeroDBoundaryCondition::get_face_name() const
{
    return face_name_;
}

void ZeroDBoundaryCondition::set_solution_ids(int flow_id, int pressure_id, double in_out_sign)
{
    flow_sol_id_ = flow_id;
    pressure_sol_id_ = pressure_id;
    in_out_sign_ = in_out_sign;
    
    #ifdef debug_zerod_bc
    std::cout << "ZeroDBoundaryCondition: Solution IDs set - flow: " << flow_sol_id_ 
              << ", pressure: " << pressure_sol_id_ << ", sign: " << in_out_sign_ << std::endl;
    #endif
}

int ZeroDBoundaryCondition::get_flow_sol_id() const
{
    return flow_sol_id_;
}

int ZeroDBoundaryCondition::get_pressure_sol_id() const
{
    return pressure_sol_id_;
}

double ZeroDBoundaryCondition::get_in_out_sign() const
{
    return in_out_sign_;
}

// =========================================================================
// Flowrate computation and access
// =========================================================================

void ZeroDBoundaryCondition::set_follower_pressure_load(bool flwP)
{
    follower_pressure_load_ = flwP;
}

bool ZeroDBoundaryCondition::get_follower_pressure_load() const
{
    return follower_pressure_load_;
}

/// @brief Compute flowrates at the boundary face at old and new timesteps
///
/// This replicates the flowrate computation done in set_bc::calc_der_cpl_bc and
/// set_bc::set_bc_cpl for coupled Neumann boundary conditions.
///
/// The flowrate is computed as the integral of velocity dotted with the face normal.
/// For struct/ustruct physics, the integral is computed on the deformed configuration.
/// For fluid/FSI/CMM physics, the integral is computed on the reference configuration.
void ZeroDBoundaryCondition::compute_flowrates(ComMod& com_mod, const CmMod& cm_mod, consts::EquationType phys)
{
    using namespace consts;
    
    #ifdef debug_zerod_bc
    std::cout << "ZeroDBoundaryCondition::compute_flowrates called" << std::endl;
    #endif
    
    if (face_ == nullptr) {
        throw std::runtime_error("[ZeroDBoundaryCondition::compute_flowrates] Face is not set.");
    }
    
    int nsd = com_mod.nsd;
    
    // Determine mechanical configuration based on physics type
    MechanicalConfigurationType cfg_o = MechanicalConfigurationType::reference;
    MechanicalConfigurationType cfg_n = MechanicalConfigurationType::reference;
    
    // For struct or ustruct, use old and new configurations to compute flowrate integral
    if ((phys == EquationType::phys_struct) || (phys == EquationType::phys_ustruct)) {
        // Must use follower pressure load for 0D coupling with struct/ustruct
        if (!follower_pressure_load_) {
            throw std::runtime_error("[ZeroDBoundaryCondition::compute_flowrates] Follower pressure load must be used for 0D coupling with struct/ustruct");
        }
        cfg_o = MechanicalConfigurationType::old_timestep;
        cfg_n = MechanicalConfigurationType::new_timestep;
    }
    // For fluid, FSI, or CMM, use reference configuration to compute flowrate integral
    else if ((phys == EquationType::phys_fluid) || (phys == EquationType::phys_FSI) || (phys == EquationType::phys_CMM)) {
        cfg_o = MechanicalConfigurationType::reference;
        cfg_n = MechanicalConfigurationType::reference;
    }
    else {
        throw std::runtime_error("[ZeroDBoundaryCondition::compute_flowrates] Invalid physics type for 0D coupling");
    }
    
    // Compute flowrates by integrating velocity over face
    // The all_fun::integ function with indices 0 to nsd-1 integrates the velocity vector
    // dotted with the face normal, giving the volumetric flowrate
    Qo_ = all_fun::integ(com_mod, cm_mod, *face_, com_mod.Yo, 0, nsd-1, false, cfg_o);
    Qn_ = all_fun::integ(com_mod, cm_mod, *face_, com_mod.Yn, 0, nsd-1, false, cfg_n);
    
    #ifdef debug_zerod_bc
    std::cout << "ZeroDBoundaryCondition::compute_flowrates - Qo: " << Qo_ << ", Qn: " << Qn_ << std::endl;
    #endif
}

/// @brief Compute average pressures at the boundary face at old and new timesteps
///
/// This replicates the pressure computation done in set_bc::calc_der_cpl_bc and
/// set_bc::set_bc_cpl for coupled Dirichlet boundary conditions.
///
/// The pressure is computed as the average pressure over the face by integrating
/// pressure (at index nsd in the solution vector) and dividing by the face area.
void ZeroDBoundaryCondition::compute_pressures(ComMod& com_mod, const CmMod& cm_mod)
{
    using namespace consts;
    
    #ifdef debug_zerod_bc
    std::cout << "ZeroDBoundaryCondition::compute_pressures called" << std::endl;
    #endif
    
    if (face_ == nullptr) {
        throw std::runtime_error("[ZeroDBoundaryCondition::compute_pressures] Face is not set.");
    }
    
    int nsd = com_mod.nsd;
    double area = face_->area;
    
    if (utils::is_zero(area)) {
        throw std::runtime_error("[ZeroDBoundaryCondition::compute_pressures] Face area is zero.");
    }
    
    // Compute average pressures by integrating pressure over face and dividing by area
    // The all_fun::integ function with index nsd integrates the pressure scalar
    Po_ = all_fun::integ(com_mod, cm_mod, *face_, com_mod.Yo, nsd) / area;
    Pn_ = all_fun::integ(com_mod, cm_mod, *face_, com_mod.Yn, nsd) / area;
    
    #ifdef debug_zerod_bc
    std::cout << "ZeroDBoundaryCondition::compute_pressures - Po: " << Po_ << ", Pn: " << Pn_ << std::endl;
    #endif
}

double ZeroDBoundaryCondition::get_Qo() const
{
    return Qo_;
}

double ZeroDBoundaryCondition::get_Qn() const
{
    return Qn_;
}

void ZeroDBoundaryCondition::set_flowrates(double Qo, double Qn)
{
    Qo_ = Qo;
    Qn_ = Qn;
}

void ZeroDBoundaryCondition::perturb_flowrate(double diff)
{
    Qn_ += diff;
}

// =========================================================================
// Pressure access
// =========================================================================

void ZeroDBoundaryCondition::set_pressure(double pressure)
{
    pressure_ = pressure;
    #ifdef debug_zerod_bc
    std::cout << "ZeroDBoundaryCondition::set_pressure - pressure: " << pressure_ << std::endl;
    #endif
}

double ZeroDBoundaryCondition::get_pressure() const
{
    return pressure_;
}

double ZeroDBoundaryCondition::get_Po() const
{
    return Po_;
}

double ZeroDBoundaryCondition::get_Pn() const
{
    return Pn_;
}

// =========================================================================
// State management
// =========================================================================

ZeroDBoundaryCondition::State ZeroDBoundaryCondition::save_state() const
{
    return State{Qn_, pressure_};
}

void ZeroDBoundaryCondition::restore_state(const State& state)
{
    Qn_ = state.Qn;
    pressure_ = state.pressure;
}

// =========================================================================
// Utility methods
// =========================================================================

const faceType* ZeroDBoundaryCondition::get_face() const
{
    return face_;
}

bool ZeroDBoundaryCondition::is_initialized() const
{
    return (face_ != nullptr);
}
