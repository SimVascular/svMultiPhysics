// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "ComMod.h"
#include "SimulationLogger.h"
#include "CoupledBoundaryCondition.h"
#include "all_fun.h"
#include "consts.h"
#include "utils.h"
#include "VtkData.h"
#include "nn.h"
#include "fils_struct.hpp"
#include <optional>
#include <unordered_map>
#include <vtkCellType.h>
#include <fstream>


CoupledBoundaryCondition::CoupledBoundaryCondition(const CoupledBoundaryCondition& other)
    : face_(other.face_)
    , cap_face_vtp_file_(other.cap_face_vtp_file_)
    , logger_(other.logger_)
    , bc_type_(other.bc_type_)
    , block_name_(other.block_name_)
    , face_name_(other.face_name_)
    , Qo_(other.Qo_)
    , Qn_(other.Qn_)
    , Po_(other.Po_)
    , Pn_(other.Pn_)
    , pressure_(other.pressure_)
    , flow_sol_id_(other.flow_sol_id_)
    , pressure_sol_id_(other.pressure_sol_id_)
    , in_out_sign_(other.in_out_sign_)
    , follower_pressure_load_(other.follower_pressure_load_)
    , has_cap_(other.has_cap_)
    , cap_global_node_ids_(other.cap_global_node_ids_)
    , cap_area_computed_(other.cap_area_computed_)
    , cap_gnNo_to_tnNo_(other.cap_gnNo_to_tnNo_)
    , cap_valM_(other.cap_valM_)
    , cap_initial_normals_(other.cap_initial_normals_)
{
    // Deep copy the cap_face_ unique_ptr if it exists
    // Note: faceType uses default copy constructor, which should properly copy all arrays
    if (other.cap_face_ != nullptr) {
        try {
            cap_face_ = std::make_unique<faceType>(*other.cap_face_);
            // Validate the copied faceType
            if (cap_face_ != nullptr) {
                // Basic validation - ensure sizes are consistent
                if (cap_face_->nNo > 0 && cap_face_->gN.size() != cap_face_->nNo) {
                    throw std::runtime_error("[CoupledBoundaryCondition::copy constructor] Invalid cap_face_ after copy: gN.size()=" +
                                            std::to_string(cap_face_->gN.size()) + " != nNo=" + std::to_string(cap_face_->nNo));
                }
                if (cap_face_->nEl > 0 && cap_face_->IEN.ncols() != cap_face_->nEl) {
                    throw std::runtime_error("[CoupledBoundaryCondition::copy constructor] Invalid cap_face_ after copy: IEN.ncols()=" +
                                            std::to_string(cap_face_->IEN.ncols()) + " != nEl=" + std::to_string(cap_face_->nEl));
                }
            }
        } catch (const std::exception& e) {
            throw std::runtime_error("[CoupledBoundaryCondition::copy constructor] Failed to copy cap_face_: " + std::string(e.what()));
        }
    } else {
        cap_face_.reset();
    }
}

CoupledBoundaryCondition& CoupledBoundaryCondition::operator=(const CoupledBoundaryCondition& other)
{
    if (this != &other) {
        face_ = other.face_;
        cap_face_vtp_file_ = other.cap_face_vtp_file_;
        logger_ = other.logger_;
        bc_type_ = other.bc_type_;
        block_name_ = other.block_name_;
        face_name_ = other.face_name_;
        Qo_ = other.Qo_;
        Qn_ = other.Qn_;
        Po_ = other.Po_;
        Pn_ = other.Pn_;
        pressure_ = other.pressure_;
        flow_sol_id_ = other.flow_sol_id_;
        pressure_sol_id_ = other.pressure_sol_id_;
        in_out_sign_ = other.in_out_sign_;
        follower_pressure_load_ = other.follower_pressure_load_;
        has_cap_ = other.has_cap_;
        cap_global_node_ids_ = other.cap_global_node_ids_;
    cap_area_computed_ = other.cap_area_computed_;
    cap_gnNo_to_tnNo_ = other.cap_gnNo_to_tnNo_;
    cap_valM_ = other.cap_valM_;
    cap_initial_normals_ = other.cap_initial_normals_;
        
        // Deep copy the cap_face_ unique_ptr if it exists
        // Note: faceType uses default copy constructor, which should properly copy all arrays
        if (other.cap_face_ != nullptr) {
            try {
                cap_face_ = std::make_unique<faceType>(*other.cap_face_);
            } catch (const std::exception& e) {
                throw std::runtime_error("[CoupledBoundaryCondition::operator=] Failed to copy cap_face_: " + std::string(e.what()));
            }
        } else {
            cap_face_.reset();
        }
    }
    return *this;
}

CoupledBoundaryCondition::CoupledBoundaryCondition(CoupledBoundaryCondition&& other) noexcept
    : face_(other.face_)
    , cap_face_vtp_file_(std::move(other.cap_face_vtp_file_))
    , logger_(other.logger_)
    , bc_type_(other.bc_type_)
    , block_name_(std::move(other.block_name_))
    , face_name_(std::move(other.face_name_))
    , Qo_(other.Qo_)
    , Qn_(other.Qn_)
    , Po_(other.Po_)
    , Pn_(other.Pn_)
    , pressure_(other.pressure_)
    , flow_sol_id_(other.flow_sol_id_)
    , pressure_sol_id_(other.pressure_sol_id_)
    , in_out_sign_(other.in_out_sign_)
    , follower_pressure_load_(other.follower_pressure_load_)
    , has_cap_(other.has_cap_)
    , cap_face_(std::move(other.cap_face_))
    , cap_global_node_ids_(std::move(other.cap_global_node_ids_))
    , cap_area_computed_(other.cap_area_computed_)
    , cap_gnNo_to_tnNo_(std::move(other.cap_gnNo_to_tnNo_))
    , cap_valM_(std::move(other.cap_valM_))
    , cap_initial_normals_(std::move(other.cap_initial_normals_))
{
    // Reset moved-from object to valid state
    other.face_ = nullptr;
    other.logger_ = nullptr;
    other.has_cap_ = false;
    other.cap_area_computed_ = false;
    other.flow_sol_id_ = -1;
    other.pressure_sol_id_ = -1;
    other.Qo_ = 0.0;
    other.Qn_ = 0.0;
    other.Po_ = 0.0;
    other.Pn_ = 0.0;
    other.pressure_ = 0.0;
}

CoupledBoundaryCondition& CoupledBoundaryCondition::operator=(CoupledBoundaryCondition&& other) noexcept
{
    if (this != &other) {
        face_ = other.face_;
        cap_face_vtp_file_ = std::move(other.cap_face_vtp_file_);
        logger_ = other.logger_;
        bc_type_ = other.bc_type_;
        block_name_ = std::move(other.block_name_);
        face_name_ = std::move(other.face_name_);
        Qo_ = other.Qo_;
        Qn_ = other.Qn_;
        Po_ = other.Po_;
        Pn_ = other.Pn_;
        pressure_ = other.pressure_;
        flow_sol_id_ = other.flow_sol_id_;
        pressure_sol_id_ = other.pressure_sol_id_;
        in_out_sign_ = other.in_out_sign_;
        follower_pressure_load_ = other.follower_pressure_load_;
        has_cap_ = other.has_cap_;
        cap_face_ = std::move(other.cap_face_);
        cap_global_node_ids_ = std::move(other.cap_global_node_ids_);
    cap_area_computed_ = other.cap_area_computed_;
    cap_gnNo_to_tnNo_ = std::move(other.cap_gnNo_to_tnNo_);
    cap_valM_ = std::move(other.cap_valM_);
    cap_initial_normals_ = std::move(other.cap_initial_normals_);
        
        // Reset moved-from object to valid state
        other.face_ = nullptr;
        other.logger_ = nullptr;
        other.has_cap_ = false;
        other.cap_area_computed_ = false;
        other.flow_sol_id_ = -1;
        other.pressure_sol_id_ = -1;
        other.Qo_ = 0.0;
        other.Qn_ = 0.0;
        other.Po_ = 0.0;
        other.Pn_ = 0.0;
        other.pressure_ = 0.0;
    }
    return *this;
}

CoupledBoundaryCondition::CoupledBoundaryCondition(consts::BoundaryConditionType bc_type, const faceType& face, const std::string& face_name,
                                               const std::string& block_name, SimulationLogger& logger)
    : face_(&face)
    , logger_(&logger)
    , bc_type_(bc_type)
    , block_name_(block_name)
    , face_name_(face_name)
{
    if (bc_type_ != consts::BoundaryConditionType::bType_Dir &&
        bc_type_ != consts::BoundaryConditionType::bType_Neu) {
        throw std::runtime_error("[CoupledBoundaryCondition] bc_type must be bType_Dir or bType_Neu.");
    }
}

CoupledBoundaryCondition::CoupledBoundaryCondition(consts::BoundaryConditionType bc_type, const faceType& face, const std::string& face_name,
                                               const std::string& block_name, const std::string& cap_face_vtp_file, SimulationLogger& logger)
    : cap_face_vtp_file_(cap_face_vtp_file)
    , face_(&face)
    , logger_(&logger)
    , bc_type_(bc_type)
    , block_name_(block_name)
    , face_name_(face_name)
{
    if (bc_type_ != consts::BoundaryConditionType::bType_Dir &&
        bc_type_ != consts::BoundaryConditionType::bType_Neu) {
        throw std::runtime_error("[CoupledBoundaryCondition] bc_type must be bType_Dir or bType_Neu.");
    }
    // Load the cap VTP file if provided
    if (!cap_face_vtp_file_.empty()) {
        load_cap_face_vtp(cap_face_vtp_file_);
    }
}

// =========================================================================
// svZeroD block configuration
// =========================================================================

void CoupledBoundaryCondition::set_block_name(const std::string& block_name)
{
    block_name_ = block_name;
}

const std::string& CoupledBoundaryCondition::get_block_name() const
{
    return block_name_;
}

void CoupledBoundaryCondition::set_face_name(const std::string& face_name)
{
    face_name_ = face_name;
}

void CoupledBoundaryCondition::set_solution_ids(int flow_id, int pressure_id, double in_out_sign)
{
    flow_sol_id_ = flow_id;
    pressure_sol_id_ = pressure_id;
    in_out_sign_ = in_out_sign;
}

int CoupledBoundaryCondition::get_flow_sol_id() const
{
    return flow_sol_id_;
}

int CoupledBoundaryCondition::get_pressure_sol_id() const
{
    return pressure_sol_id_;
}

double CoupledBoundaryCondition::get_in_out_sign() const
{
    return in_out_sign_;
}

// =========================================================================
// Flowrate computation and access
// =========================================================================

void CoupledBoundaryCondition::set_follower_pressure_load(bool flwP)
{
    follower_pressure_load_ = flwP;
}

bool CoupledBoundaryCondition::get_follower_pressure_load() const
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
void CoupledBoundaryCondition::compute_flowrates(ComMod& com_mod, const CmMod& cm_mod, consts::EquationType phys)
{
    using namespace consts;
    
    
    if (face_ == nullptr) {
        throw std::runtime_error("[CoupledBoundaryCondition::compute_flowrates] Face is not set.");
    }
    
    int nsd = com_mod.nsd;
    
    // Determine mechanical configuration based on physics type
    MechanicalConfigurationType cfg_o = MechanicalConfigurationType::reference;
    MechanicalConfigurationType cfg_n = MechanicalConfigurationType::reference;
    
    // For struct or ustruct, use old and new configurations to compute flowrate integral
    if ((phys == EquationType::phys_struct) || (phys == EquationType::phys_ustruct)) {
        // Must use follower pressure load for 0D coupling with struct/ustruct
        if (!follower_pressure_load_) {
            throw std::runtime_error("[CoupledBoundaryCondition::compute_flowrates] Follower pressure load must be used for 0D coupling with struct/ustruct");
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
        throw std::runtime_error("[CoupledBoundaryCondition::compute_flowrates] Invalid physics type for 0D coupling");
    }
    
    // Compute flowrates by integrating velocity over face
    // The all_fun::integ function with indices 0 to nsd-1 integrates the velocity vector
    // dotted with the face normal, giving the volumetric flowrate
    Qo_ = all_fun::integ(com_mod, cm_mod, *face_, com_mod.Yo, 0, nsd-1, false, cfg_o);
    Qn_ = all_fun::integ(com_mod, cm_mod, *face_, com_mod.Yn, 0, nsd-1, false, cfg_n);
    
    // Add cap contribution. prepare_cap_gathered_data (no-op in serial) gathers Yo and Yn in one go for parallel.
    if (has_cap()) {
        prepare_cap_gathered_data(com_mod, cm_mod, com_mod.Yo, com_mod.Yn, 0, nsd, cfg_o, cfg_n);
        double Qo_cap = integrate_over_cap(com_mod, cm_mod, com_mod.Yo, 0, nsd-1, cfg_o);
        double Qn_cap = integrate_over_cap(com_mod, cm_mod, com_mod.Yn, 0, nsd-1, cfg_n);
        Qo_ += Qo_cap;
        Qn_ += Qn_cap;
    }
    
}

/// @brief Compute average pressures at the boundary face at old and new timesteps
///
/// This replicates the pressure computation done in set_bc::calc_der_cpl_bc and
/// set_bc::set_bc_cpl for coupled Dirichlet boundary conditions.
///
/// The pressure is computed as the average pressure over the face by integrating
/// pressure (at index nsd in the solution vector) and dividing by the face area.
void CoupledBoundaryCondition::compute_pressures(ComMod& com_mod, const CmMod& cm_mod)
{
    using namespace consts;
    
    
    if (face_ == nullptr) {
        throw std::runtime_error("[CoupledBoundaryCondition::compute_pressures] Face is not set.");
    }
    
    int nsd = com_mod.nsd;
    double area = face_->area;
    
    Po_ = all_fun::integ(com_mod, cm_mod, *face_, com_mod.Yo, nsd) / area;
    Pn_ = all_fun::integ(com_mod, cm_mod, *face_, com_mod.Yn, nsd) / area;
    
}

double CoupledBoundaryCondition::get_Qo() const
{
    return Qo_;
}

double CoupledBoundaryCondition::get_Qn() const
{
    return Qn_;
}

void CoupledBoundaryCondition::set_flowrates(double Qo, double Qn)
{
    Qo_ = Qo;
    Qn_ = Qn;
}

void CoupledBoundaryCondition::perturb_flowrate(double diff)
{
    Qn_ += diff;
}

// =========================================================================
// Pressure access
// =========================================================================

void CoupledBoundaryCondition::set_pressure(double pressure)
{
    pressure_ = pressure;
}

double CoupledBoundaryCondition::get_pressure() const
{
    return pressure_;
}

double CoupledBoundaryCondition::get_Po() const
{
    return Po_;
}

double CoupledBoundaryCondition::get_Pn() const
{
    return Pn_;
}

// =========================================================================
// State management
// =========================================================================

CoupledBoundaryCondition::State CoupledBoundaryCondition::save_state() const
{
    return State{Qn_, pressure_};
}

void CoupledBoundaryCondition::restore_state(const State& state)
{
    Qn_ = state.Qn;
    pressure_ = state.pressure;
}

// =========================================================================
// Utility methods
// =========================================================================

void CoupledBoundaryCondition::distribute(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, const faceType& face)
{
    #define n_debug_coupled_distribute

    // Update face pointer to local face
    face_ = &face;

    const bool is_slave = cm.slv(cm_mod);
    
    // Distribute BC type (Dirichlet or Neumann)
    int bc_type_int = static_cast<int>(bc_type_);
    cm.bcast(cm_mod, &bc_type_int);
    if (is_slave) {
        bc_type_ = static_cast<consts::BoundaryConditionType>(bc_type_int);
    }
    
    // Distribute block name
    cm.bcast(cm_mod, block_name_);
    
    // Distribute face name
    cm.bcast(cm_mod, face_name_);
    
    // Distribute follower pressure load flag
    cm.bcast(cm_mod, &follower_pressure_load_);
    
    // Distribute solution IDs
    cm.bcast(cm_mod, &flow_sol_id_);
    cm.bcast(cm_mod, &pressure_sol_id_);
    cm.bcast(cm_mod, &in_out_sign_);

    // Distribute cap flag so all ranks agree (master loaded cap_face_; slaves have has_cap_ set here)
    cm.bcast(cm_mod, &has_cap_);
}

void CoupledBoundaryCondition::set_face(const faceType& face)
{
    face_ = &face;
}

const faceType* CoupledBoundaryCondition::get_face() const
{
    return face_;
}

bool CoupledBoundaryCondition::is_initialized() const
{
    return (face_ != nullptr);
}

// =========================================================================
// Cap: helpers and implementation 
// =========================================================================

// Map VTK cell types to ElementType (same as vtk_xml_parser.cpp)
static consts::ElementType vtk_cell_type_to_element_type(int vtk_cell_type)
{
    using namespace consts;
    switch (vtk_cell_type) {
        case VTK_TRIANGLE: return ElementType::TRI3;
        case VTK_QUADRATIC_TRIANGLE: return ElementType::TRI6;
        case VTK_QUAD: return ElementType::QUD4;
        case VTK_QUADRATIC_QUAD: return ElementType::QUD8;
        case VTK_BIQUADRATIC_QUAD: return ElementType::QUD9;
        case VTK_LINE: return ElementType::LIN1;
        case VTK_TETRA: return ElementType::TET4;
        case VTK_QUADRATIC_TETRA: return ElementType::TET10;
        case VTK_HEXAHEDRON: return ElementType::HEX8;
        case VTK_QUADRATIC_HEXAHEDRON: return ElementType::HEX20;
        case VTK_TRIQUADRATIC_HEXAHEDRON: return ElementType::HEX27;
        case VTK_WEDGE: return ElementType::WDG;
        default:
            throw std::runtime_error("[CoupledBoundaryCondition] Unsupported VTK cell type " + 
                                    std::to_string(vtk_cell_type) + " in cap VTP file.");
    }
}

template <typename GetXlFn, typename ComputeJacNFn, typename OnGaussFn>
void for_each_cap_gauss_point(const faceType* cap_face, int nsd, int insd,
                              GetXlFn&& get_xl,
                              ComputeJacNFn&& compute_jac_n,
                              OnGaussFn&& on_gauss)
{
    for (int e = 0; e < cap_face->nEl; e++) {
        Array<double> xl = get_xl(e);
        for (int g = 0; g < cap_face->nG; g++) {
            auto [Jac, n] = compute_jac_n(xl, e, g, nsd, insd);
            on_gauss(e, g, Jac, n);
        }
    }
}

void CoupledBoundaryCondition::load_cap_face_vtp(const std::string& vtp_file_path)
{
    // Safety check: ensure face_ is set before loading cap
    if (face_ == nullptr) {
        throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Cannot load cap: face_ is null. Call set_face() first.");
    }
    
    cap_face_vtp_file_ = vtp_file_path;
    has_cap_ = false;
    cap_face_.reset();
    cap_gnNo_to_tnNo_.clear();  // Clear the mapping since we're reloading
    cap_area_computed_ = false;  // Reset area computation flag
    cap_initial_normals_.resize(0, 0);  // Clear initial normals

    if (vtp_file_path.empty()) {
        return;
    }
    
    // Check if file exists and is readable
    std::ifstream file_check(vtp_file_path);
    if (!file_check.good()) {
        throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Cannot open cap VTP file '" + vtp_file_path + "' for reading.");
    }
    file_check.close();
    
    
    // Load the VTP file - wrap in try-catch to handle any exceptions
    // Note: VtkVtpData constructor with file_name calls read_file() which may throw or crash
    VtkVtpData vtp_data;
    try {
        // Construct VtkVtpData and read file in one step
        // The constructor with file_name will call read_file() internally
        vtp_data = VtkVtpData(vtp_file_path, true);
    } catch (const std::exception& e) {
        throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Failed to construct VtkVtpData from file '" + 
                                vtp_file_path + "': " + e.what());
    } catch (...) {
        throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Unknown error constructing VtkVtpData from file '" + 
                                vtp_file_path + "'. This may indicate a crash in the VTK library.");
    }
    
    int nNo = 0;
    try {
        nNo = vtp_data.num_points();
    } catch (const std::exception& e) {
        throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Failed to get number of points from VTP file '" + 
                                vtp_file_path + "': " + e.what());
    } catch (...) {
        throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Unknown error getting number of points from file '" + 
                                vtp_file_path + "'");
    }
    if (nNo == 0) {
        throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Cap VTP file '" + vtp_file_path + "' does not contain any points.");
    }
    
    int num_elems = 0;
    try {
        num_elems = vtp_data.num_elems();
    } catch (const std::exception& e) {
        throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Failed to get number of elements from VTP file '" + 
                                vtp_file_path + "': " + e.what());
    } catch (...) {
        throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Unknown error getting number of elements from file '" + 
                                vtp_file_path + "'");
    }
    if (num_elems == 0) {
        throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Cap VTP file '" + vtp_file_path + "' does not contain any elements.");
    }
    
    // Get GlobalNodeID from VTP file
    bool has_global_node_id = false;
    try {
        has_global_node_id = vtp_data.has_point_data("GlobalNodeID");
    } catch (const std::exception& e) {
        throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Failed to check for GlobalNodeID in VTP file '" + 
                                vtp_file_path + "': " + e.what());
    } catch (...) {
        throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Unknown error checking for GlobalNodeID in file '" + 
                                vtp_file_path + "'");
    }
    if (!has_global_node_id) {
        throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Cap VTP file '" + vtp_file_path + "' does not contain 'GlobalNodeID' point data.");
    }
    
    try {
        cap_global_node_ids_.resize(nNo);
        vtp_data.copy_point_data("GlobalNodeID", cap_global_node_ids_);
    } catch (const std::exception& e) {
        throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Failed to copy GlobalNodeID from VTP file '" + 
                                vtp_file_path + "': " + e.what());
    } catch (...) {
        throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Unknown error copying GlobalNodeID from file '" + 
                                vtp_file_path + "'");
    }
    
    // Get connectivity from VTP file
    Array<int> conn;
    int eNoN = 0;
    int vtk_cell_type = 0;
    try {
        conn = vtp_data.get_connectivity();
        eNoN = vtp_data.np_elem();
        vtk_cell_type = vtp_data.elem_type();  // This returns VTK cell type (e.g., 5 for VTK_TRIANGLE)
    } catch (const std::exception& e) {
        throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Failed to get connectivity from VTP file '" + 
                                vtp_file_path + "': " + e.what());
    } catch (...) {
        throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Unknown error getting connectivity from file '" + 
                                vtp_file_path + "'");
    }
    if (eNoN <= 0) {
        throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Invalid number of nodes per element: " + std::to_string(eNoN));
    }
    
    // Convert VTK cell type to ElementType enum
    consts::ElementType eType = vtk_cell_type_to_element_type(vtk_cell_type);
    
    // Create a temporary faceType structure for the cap
    cap_face_ = std::make_unique<faceType>();
    cap_face_->name = face_name_ + "_cap";
    cap_face_->iM = face_->iM;  // Use the same mesh index as the main face
    cap_face_->nNo = nNo;
    cap_face_->nEl = num_elems;
    cap_face_->gnEl = num_elems;
    cap_face_->eNoN = eNoN;
    cap_face_->eType = eType;
    
    // Set global node IDs (map GlobalNodeID from VTP to mesh global node indices)
    // Verify cap_global_node_ids_ is properly sized
    if (cap_global_node_ids_.size() != nNo) {
        throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] cap_global_node_ids_ size mismatch: " +
                                std::to_string(cap_global_node_ids_.size()) + " != " + std::to_string(nNo));
    }
    cap_face_->gN.resize(nNo);
    for (int a = 0; a < nNo; a++) {
        cap_face_->gN(a) = cap_global_node_ids_(a) - 1;  // Convert from 1-based to 0-based
    }
    
    // Map connectivity from local cap node indices to global mesh node indices
    // The VTP connectivity contains local node indices (0 to nNo-1)
    cap_face_->IEN.resize(eNoN, num_elems);
    
    for (int e = 0; e < num_elems; e++) {
        for (int a = 0; a < eNoN; a++) {
            int local_node_idx = conn(a, e);  // Local node index in cap (0 to nNo-1)
            if (local_node_idx < 0 || local_node_idx >= nNo) {
                throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Invalid local node index " + 
                                        std::to_string(local_node_idx) + " in cap connectivity (element " + 
                                        std::to_string(e) + ", node " + std::to_string(a) + ", nNo=" + std::to_string(nNo) + ").");
            }
            cap_face_->IEN(a, e) = cap_face_->gN(local_node_idx);  // Map to global mesh node index
        }
    }
    
    // Load initial element normals from VTP file if available
    // The normals are stored as cell data array called "Normals"
    try {
        bool has_normals = vtp_data.has_cell_data("Normals");
        if (has_normals) {
            auto [num_comp, num_tuples] = vtp_data.get_cell_data_dimensions("Normals");

            if (num_comp == 0 || num_tuples == 0) {
                // Provide more detailed error message
                std::string error_msg = "[CoupledBoundaryCondition::load_cap_face_vtp] Normals array exists but has zero size. ";
                error_msg += "num_components=" + std::to_string(num_comp) + ", num_tuples=" + std::to_string(num_tuples);
                error_msg += ", expected num_elems=" + std::to_string(num_elems);
                error_msg += ". The array may not be a numeric type or may be empty.";
                throw std::runtime_error(error_msg);
            }
            
            if (num_tuples != num_elems) {
                throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Normals array size mismatch: " +
                                        std::to_string(num_tuples) + " != " + std::to_string(num_elems));
            }
            
            if (num_comp != 2 && num_comp != 3) {
                throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Invalid number of components in Normals array: " +
                                        std::to_string(num_comp) + " (expected 2 or 3)");
            }
            
            // Resize array to correct dimensions
            cap_initial_normals_.resize(num_comp, num_elems);
            vtp_data.copy_cell_data("Normals", cap_initial_normals_);
            
            // Verify that data was actually copied
            if (cap_initial_normals_.nrows() != num_comp || cap_initial_normals_.ncols() != num_elems) {
                throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Failed to copy Normals data. "
                                        "Expected size: " + std::to_string(num_comp) + "x" + std::to_string(num_elems) +
                                        ", Actual size: " + std::to_string(cap_initial_normals_.nrows()) + "x" + 
                                        std::to_string(cap_initial_normals_.ncols()));
            }
        } else {
            // Normals not found - initialize to empty array
            cap_initial_normals_.resize(0, 0);
        }
    } catch (const std::exception& e) {
        throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Failed to load Normals from VTP file '" + 
                                vtp_file_path + "': " + e.what());
    } catch (...) {
        throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Unknown error loading Normals from file '" + 
                                vtp_file_path + "'");
    }

    has_cap_ = true;
}

// =========================================================================
// Cap surface integration
// =========================================================================

void CoupledBoundaryCondition::initialize_cap_integration(ComMod& com_mod, const CmMod& cm_mod)
{
    using namespace consts;
    
    if (!cap_face_ready()) {
        return;
    }
    if (cap_face_->nEl == 0 || cap_face_->nNo == 0) {
        throw std::runtime_error("[CoupledBoundaryCondition::initialize_cap_integration] Cap face is not properly initialized.");
    }
    
    // Skip if already initialized
    if (cap_face_->nG > 0 && cap_face_->w.size() > 0 && cap_face_->N.nrows() > 0 && 
        cap_face_->Nx.nslices() == cap_face_->nG && !cap_gnNo_to_tnNo_.empty()) {
        return;
    }
    
    int nsd = com_mod.nsd;
    int insd = nsd - 1;
    
    
    // Build reverse map from gnNo to tnNo using com_mod.ltg
    // com_mod.ltg maps from tnNo to gnNo, so we reverse it
    // IEN contains gnNo indices from load_cap_face_vtp, we'll map them on-the-fly during integration
    // Note: com_mod.ltg might not be set up yet when called from read_files.cpp.
    // If it's not ready, we'll skip initialization and it will be done later when needed.
    auto& msh = com_mod.msh[face_->iM];
    if (com_mod.ltg.size() == 0 || com_mod.tnNo == 0) {
        // ltg not ready yet - will be initialized later when ltg is set up
        return;
    }
    if (com_mod.ltg.size() != com_mod.tnNo) {
        throw std::runtime_error("[CoupledBoundaryCondition::initialize_cap_integration] com_mod.ltg size (" + 
                                std::to_string(com_mod.ltg.size()) + ") does not match com_mod.tnNo (" + 
                                std::to_string(com_mod.tnNo) + ")");
    }
    
    cap_gnNo_to_tnNo_.clear();
    cap_gnNo_to_tnNo_.reserve(com_mod.tnNo);  // Reserve space for better performance
    for (int a = 0; a < com_mod.tnNo; a++) {
        int gnNo_idx = com_mod.ltg(a);
        // Check for valid gnNo index
        if (gnNo_idx < 0 || gnNo_idx >= msh.gnNo) {
            throw std::runtime_error("[CoupledBoundaryCondition::initialize_cap_integration] Invalid gnNo index " + 
                                    std::to_string(gnNo_idx) + " in com_mod.ltg at position " + 
                                    std::to_string(a) + " (msh.gnNo=" + std::to_string(msh.gnNo) + ")");
        }
        cap_gnNo_to_tnNo_[gnNo_idx] = a;
    }
    
    try {
        // Determine nG based on element type (similar to select_eleb)
        // For common surface elements:
        if (cap_face_->eType == ElementType::TRI3) {
            // Use 1-point centroid rule for TRI3 (more efficient than 3-point rule)
            cap_face_->nG = 1;
        } else if (cap_face_->eType == ElementType::QUD4) {
            cap_face_->nG = 4;  // 4 Gauss points for QUD4
        } else if (cap_face_->eType == ElementType::TRI6) {
            cap_face_->nG = 3;  // 3 Gauss points for TRI6
        } else {
            // Default: use 1 Gauss point per element node
            cap_face_->nG = cap_face_->eNoN;
        }
        
        // Allocate and initialize Gauss points and weights
        cap_face_->w.resize(cap_face_->nG);
        cap_face_->xi.resize(insd, cap_face_->nG);
        
        // For TRI3 with 1 point, manually set up centroid rule (get_gip is hardcoded for 3 points)
        if (cap_face_->eType == ElementType::TRI3 && cap_face_->nG == 1) {
            cap_face_->w(0) = 0.5;  // Weight for centroid rule on reference triangle
            cap_face_->xi(0, 0) = 1.0/3.0;  // Centroid coordinate
            cap_face_->xi(1, 0) = 1.0/3.0;  // Centroid coordinate
        } else {
            // For other element types, use get_gip
            nn::get_gip(insd, cap_face_->eType, cap_face_->nG, cap_face_->w, cap_face_->xi);
        }
        
        // Allocate and compute shape functions
        cap_face_->N.resize(cap_face_->eNoN, cap_face_->nG);
        cap_face_->Nx.resize(insd, cap_face_->eNoN, cap_face_->nG);
        for (int g = 0; g < cap_face_->nG; g++) {
            nn::get_gnn(insd, cap_face_->eType, cap_face_->eNoN, g, cap_face_->xi, cap_face_->N, cap_face_->Nx);
        }
    } catch (const std::exception& e) {
        throw std::runtime_error("[CoupledBoundaryCondition::initialize_cap_integration] Failed to initialize cap face shape functions: " + 
                                std::string(e.what()));
    }
}

// =========================================================================
// Cap integration helper functions
// =========================================================================

Array<double> CoupledBoundaryCondition::update_cap_element_position(ComMod& com_mod, int e, 
                                                                  const std::unordered_map<int, int>& gnNo_to_tnNo,
                                                                  consts::MechanicalConfigurationType cfg)
{
    using namespace consts;
    
    // Caller must ensure cap_face_ready(); e and IEN are validated here
    if (cap_face_->IEN.nrows() == 0 || cap_face_->IEN.ncols() == 0) {
        throw std::runtime_error("[CoupledBoundaryCondition::update_cap_element_position] cap_face_->IEN is not allocated.");
    }
    if (e < 0 || e >= cap_face_->IEN.ncols()) {
        throw std::runtime_error("[CoupledBoundaryCondition::update_cap_element_position] Element index e=" + 
                                std::to_string(e) + " is out of bounds (IEN.ncols()=" + std::to_string(cap_face_->IEN.ncols()) + ").");
    }
    if (cap_face_->eNoN <= 0) {
        throw std::runtime_error("[CoupledBoundaryCondition::update_cap_element_position] cap_face_->eNoN is invalid: " + 
                                std::to_string(cap_face_->eNoN));
    }
    
    int nsd = com_mod.nsd;
    Array<double> xl(nsd, cap_face_->eNoN);
    
    for (int a = 0; a < cap_face_->eNoN; a++) {
        // IEN contains gnNo index, map it to tnNo index
        int gnNo_idx = cap_face_->IEN(a, e);
        auto it = gnNo_to_tnNo.find(gnNo_idx);
        if (it == gnNo_to_tnNo.end()) {
            throw std::runtime_error("[CoupledBoundaryCondition::update_cap_element_position] IEN entry (element " + 
                                    std::to_string(e) + ", node " + std::to_string(a) + 
                                    ") contains invalid gnNo index " + std::to_string(gnNo_idx) + 
                                    " not found in com_mod.ltg mapping.");
        }
        int Ac = it->second;  // tnNo index
        
        // Bounds checking
        if (Ac < 0 || Ac >= com_mod.tnNo) {
            std::string msg = "[CoupledBoundaryCondition::update_cap_element_position] Invalid node index Ac=" + 
                             std::to_string(Ac) + " (tnNo=" + std::to_string(com_mod.tnNo) + 
                             ") at element " + std::to_string(e) + ", local node " + std::to_string(a);
            throw std::runtime_error(msg);
        }
        
        // Bounds check for com_mod.x
        if (Ac >= com_mod.x.ncols() || nsd > com_mod.x.nrows()) {
            throw std::runtime_error("[CoupledBoundaryCondition::update_cap_element_position] Invalid bounds for com_mod.x: Ac=" + 
                                    std::to_string(Ac) + ", x.ncols()=" + std::to_string(com_mod.x.ncols()) + 
                                    ", nsd=" + std::to_string(nsd) + ", x.nrows()=" + std::to_string(com_mod.x.nrows()));
        }
        
        xl.set_col(a, com_mod.x.col(Ac));
        
        // Apply displacement if needed based on configuration
        if (cfg == MechanicalConfigurationType::old_timestep) {
            // Bounds check for com_mod.Do
            if (Ac >= com_mod.Do.ncols() || nsd > com_mod.Do.nrows()) {
                throw std::runtime_error("[CoupledBoundaryCondition::update_cap_element_position] Invalid bounds for com_mod.Do: Ac=" + 
                                        std::to_string(Ac) + ", Do.ncols()=" + std::to_string(com_mod.Do.ncols()) + 
                                        ", nsd=" + std::to_string(nsd) + ", Do.nrows()=" + std::to_string(com_mod.Do.nrows()));
            }
            for (int i = 0; i < nsd; i++) {
                xl(i, a) += com_mod.Do(i, Ac);
            }
        } else if (cfg == MechanicalConfigurationType::new_timestep) {
            // Bounds check for com_mod.Dn
            if (Ac >= com_mod.Dn.ncols() || nsd > com_mod.Dn.nrows()) {
                throw std::runtime_error("[CoupledBoundaryCondition::update_cap_element_position] Invalid bounds for com_mod.Dn: Ac=" + 
                                        std::to_string(Ac) + ", Dn.ncols()=" + std::to_string(com_mod.Dn.ncols()) + 
                                        ", nsd=" + std::to_string(nsd) + ", Dn.nrows()=" + std::to_string(com_mod.Dn.nrows()));
            }
            for (int i = 0; i < nsd; i++) {
                xl(i, a) += com_mod.Dn(i, Ac);
            }
        }
    }
    
    return xl;
}

Array<double> CoupledBoundaryCondition::update_cap_element_position(int e, consts::MechanicalConfigurationType cfg,
                                                                  const Array<double>& cap_x, const Array<double>& cap_Do, const Array<double>& cap_Dn,
                                                                  const std::unordered_map<int, int>& gnNo_to_capIdx)
{
    if (cap_face_->IEN.nrows() == 0 || cap_face_->IEN.ncols() == 0) {
        throw std::runtime_error("[CoupledBoundaryCondition::update_cap_element_position(gathered)] cap_face_->IEN not ready.");
    }
    int nsd = cap_x.nrows();
    Array<double> xl(nsd, cap_face_->eNoN);
    for (int a = 0; a < cap_face_->eNoN; a++) {
        int gnNo_idx = cap_face_->IEN(a, e);
        auto it = gnNo_to_capIdx.find(gnNo_idx);
        if (it == gnNo_to_capIdx.end()) {
            throw std::runtime_error("[CoupledBoundaryCondition::update_cap_element_position(gathered)] IEN gnNo " +
                                    std::to_string(gnNo_idx) + " not in gnNo_to_capIdx.");
        }
        int cap_idx = it->second;
        for (int i = 0; i < nsd; i++) {
            xl(i, a) = cap_x(i, cap_idx);
        }
        if (cfg == consts::MechanicalConfigurationType::old_timestep) {
            for (int i = 0; i < nsd; i++) xl(i, a) += cap_Do(i, cap_idx);
        } else if (cfg == consts::MechanicalConfigurationType::new_timestep) {
            for (int i = 0; i < nsd; i++) xl(i, a) += cap_Dn(i, cap_idx);
        }
    }
    return xl;
}

void CoupledBoundaryCondition::gather_cap_node_data_to_master(ComMod& com_mod, const CmMod& cm_mod,
                                                            const Array<double>& s_old, const Array<double>& s_new,
                                                            int l, int s_comps, consts::MechanicalConfigurationType cfg_o, consts::MechanicalConfigurationType cfg_n,
                                                            int cap_nNo, int nsd,
                                                            Array<double>& cap_x, Array<double>& cap_Do, Array<double>& cap_Dn,
                                                            Array<double>& cap_Yo, Array<double>& cap_Yn)
{
    (void)cfg_o;
    (void)cfg_n;
    auto& cm = com_mod.cm;
    const int root = cm_mod.master;
    const int nProcs = cm.np();
    if (nProcs == 1) {
        if (!cap_face_ready()) return;
        for (int a = 0; a < cap_nNo; a++) {
            auto it = cap_gnNo_to_tnNo_.find(cap_face_->gN(a));
            if (it == cap_gnNo_to_tnNo_.end()) continue;
            int Ac = it->second;
            for (int i = 0; i < nsd; i++) {
                cap_x(i, a) = com_mod.x(i, Ac);
                cap_Do(i, a) = com_mod.Do(i, Ac);
                cap_Dn(i, a) = com_mod.Dn(i, Ac);
            }
            for (int i = 0; i < s_comps; i++) {
                cap_Yo(i, a) = s_old(l + i, Ac);
                cap_Yn(i, a) = s_new(l + i, Ac);
            }
        }
        return;
    }
    // Build local map gnNo -> Ac from this rank's ltg so every rank can send its owned cap nodes
    std::unordered_map<int, int> gnNo_to_Ac;
    for (int Ac = 0; Ac < com_mod.tnNo; Ac++)
        gnNo_to_Ac[com_mod.ltg(Ac)] = Ac;

    const int per_node = 1 + 3 * nsd + 2 * s_comps;
    Vector<double> send_buf;
    int n_owned = 0;
    for (int a = 0; a < cap_nNo; a++) {
        if (gnNo_to_Ac.find(cap_gN_broadcast_(a)) != gnNo_to_Ac.end()) n_owned++;
    }
    send_buf.resize(n_owned * per_node);
    int idx = 0;
    for (int a = 0; a < cap_nNo; a++) {
        auto it = gnNo_to_Ac.find(cap_gN_broadcast_(a));
        if (it == gnNo_to_Ac.end()) continue;
        int Ac = it->second;
        send_buf(idx++) = static_cast<double>(a);
        for (int i = 0; i < nsd; i++) send_buf(idx++) = com_mod.x(i, Ac);
        for (int i = 0; i < nsd; i++) send_buf(idx++) = com_mod.Do(i, Ac);
        for (int i = 0; i < nsd; i++) send_buf(idx++) = com_mod.Dn(i, Ac);
        for (int i = 0; i < s_comps; i++) send_buf(idx++) = s_old(l + i, Ac);
        for (int i = 0; i < s_comps; i++) send_buf(idx++) = s_new(l + i, Ac);
    }
    const int my_send_count = static_cast<int>(send_buf.size());
    Vector<int> send_count_vec(1);
    send_count_vec(0) = my_send_count;
    Vector<int> recv_counts(nProcs);
    cm.gather(cm_mod, send_count_vec, recv_counts, root);
    Vector<double> recv_buf;
    Vector<int> displs(nProcs);
    if (cm.idcm() == root) {
        int total = 0;
        for (int i = 0; i < nProcs; i++) {
            displs(i) = total;
            total += recv_counts(i);
        }
        recv_buf.resize(total);
    }
    cm.gatherv(cm_mod, send_buf, recv_buf, recv_counts, displs, root);
    if (cm.idcm() == root) {
        int pos = 0;
        while (pos < static_cast<int>(recv_buf.size())) {
            int cap_idx = static_cast<int>(recv_buf(pos++));
            for (int i = 0; i < nsd; i++) cap_x(i, cap_idx) = recv_buf(pos++);
            for (int i = 0; i < nsd; i++) cap_Do(i, cap_idx) = recv_buf(pos++);
            for (int i = 0; i < nsd; i++) cap_Dn(i, cap_idx) = recv_buf(pos++);
            for (int i = 0; i < s_comps; i++) cap_Yo(i, cap_idx) = recv_buf(pos++);
            for (int i = 0; i < s_comps; i++) cap_Yn(i, cap_idx) = recv_buf(pos++);
        }
    }
}

void CoupledBoundaryCondition::prepare_cap_gathered_data(ComMod& com_mod, const CmMod& cm_mod,
                                                       const Array<double>& Yo, const Array<double>& Yn,
                                                       int l, int s_comps, consts::MechanicalConfigurationType cfg_o, consts::MechanicalConfigurationType cfg_n)
{
    auto& cm = com_mod.cm;
    if (cm.seq()) {
        if (!cap_face_ready()) return;
        const int cap_nNo = cap_face_->nNo;
        if (cap_nNo == 0) return;
        cap_nNo_gathered_ = cap_nNo;
        const int nsd = com_mod.nsd;
        cap_x_gathered_.resize(nsd, cap_nNo);
        cap_Do_gathered_.resize(nsd, cap_nNo);
        cap_Dn_gathered_.resize(nsd, cap_nNo);
        cap_Yo_gathered_.resize(s_comps, cap_nNo);
        cap_Yn_gathered_.resize(s_comps, cap_nNo);
        for (int a = 0; a < cap_nNo; a++) {
            auto it = cap_gnNo_to_tnNo_.find(cap_face_->gN(a));
            if (it == cap_gnNo_to_tnNo_.end()) continue;
            int Ac = it->second;
            for (int i = 0; i < nsd; i++) {
                cap_x_gathered_(i, a) = com_mod.x(i, Ac);
                cap_Do_gathered_(i, a) = com_mod.Do(i, Ac);
                cap_Dn_gathered_(i, a) = com_mod.Dn(i, Ac);
            }
            for (int i = 0; i < s_comps; i++) {
                cap_Yo_gathered_(i, a) = Yo(l + i, Ac);
                cap_Yn_gathered_(i, a) = Yn(l + i, Ac);
            }
        }
        return;
    }
    int cap_nNo = 0;
    if (cm.idcm() == cm_mod.master && cap_face_ready()) {
        cap_nNo = cap_face_->nNo;
    }
    cm.bcast(cm_mod, &cap_nNo);
    cap_nNo_gathered_ = cap_nNo;
    if (cap_nNo == 0) {
        return;
    }
    // Broadcast cap global node IDs so all ranks can determine which cap nodes they own and send
    cap_gN_broadcast_.resize(cap_nNo);
    if (cm.idcm() == cm_mod.master && cap_face_ready()) {
        for (int a = 0; a < cap_nNo; a++)
            cap_gN_broadcast_(a) = cap_face_->gN(a);
    }
    cm.bcast(cm_mod, cap_gN_broadcast_);

    const int nsd = com_mod.nsd;
    if (cm.idcm() == cm_mod.master) {
        cap_x_gathered_.resize(nsd, cap_nNo);
        cap_Do_gathered_.resize(nsd, cap_nNo);
        cap_Dn_gathered_.resize(nsd, cap_nNo);
        cap_Yo_gathered_.resize(s_comps, cap_nNo);
        cap_Yn_gathered_.resize(s_comps, cap_nNo);
    }
    gather_cap_node_data_to_master(com_mod, cm_mod, Yo, Yn, l, s_comps, cfg_o, cfg_n, cap_nNo, nsd,
                                   cap_x_gathered_, cap_Do_gathered_, cap_Dn_gathered_, cap_Yo_gathered_, cap_Yn_gathered_);
}

std::pair<double, Vector<double>> CoupledBoundaryCondition::compute_cap_jacobian_and_normal(const Array<double>& xl,
                                                                                           int e, int g, int nsd, int insd)
{
    if (xl.nrows() != nsd || xl.ncols() != cap_face_->eNoN) {
        throw std::runtime_error("[CoupledBoundaryCondition::compute_cap_jacobian_and_normal] xl has wrong dimensions: " +
                                std::to_string(xl.nrows()) + "x" + std::to_string(xl.ncols()) + 
                                " (expected " + std::to_string(nsd) + "x" + std::to_string(cap_face_->eNoN) + ").");
    }
    
    // Get shape function derivatives at this Gauss point
    Array<double> Nx_g = cap_face_->Nx.slice(g);  // (insd x eNoN)
    
    // Compute tangent vectors using basis functions
    // For 3D surface: compute ∂x/∂ξ₁ and ∂x/∂ξ₂
    // For 2D curve: compute ∂x/∂ξ
    Array<double> xXi(nsd, insd);  // Tangent vectors
    xXi = 0.0;
    
    for (int a = 0; a < cap_face_->eNoN; a++) {
        for (int i = 0; i < insd; i++) {
            for (int j = 0; j < nsd; j++) {
                xXi(j, i) += xl(j, a) * Nx_g(i, a);
            }
        }
    }
    
    // Compute Jacobian (area element) from tangent vectors
    double Jac = 0.0;
    Vector<double> n(nsd);
    if (nsd == 3 && insd == 2) {
        // 3D surface: Jacobian = ||∂x/∂ξ₁ × ∂x/∂ξ₂||
        n = utils::cross(xXi);
        Jac = sqrt(utils::norm(n));  // utils::norm returns squared norm, need sqrt
    } else if (nsd == 2 && insd == 1) {
        // 2D curve: Jacobian = ||∂x/∂ξ||
        Jac = sqrt(utils::norm(xXi.col(0)));  // utils::norm returns squared norm, need sqrt
        // For 2D, compute normal from tangent
        n(0) = -xXi(1, 0);
        n(1) = xXi(0, 0);
    } else {
        throw std::runtime_error("[CoupledBoundaryCondition::compute_cap_jacobian_and_normal] Unsupported nsd/insd combination: " + 
                                std::to_string(nsd) + "/" + std::to_string(insd));
    }
    
    if (utils::is_zero(Jac)) {
        throw std::runtime_error("[CoupledBoundaryCondition::compute_cap_jacobian_and_normal] Zero Jacobian at Gauss point " + 
                                std::to_string(g));
    }
    
    // Normalize normal vector
    n = n / Jac;
    
    // Check if initial normals are available and ensure consistency
    if (cap_initial_normals_.ncols() > 0 && cap_initial_normals_.nrows() == nsd) {
        // Check bounds for element index
        if (e < 0 || e >= cap_initial_normals_.ncols()) {
            throw std::runtime_error("[CoupledBoundaryCondition::compute_cap_jacobian_and_normal] Element index e=" + 
                                    std::to_string(e) + " is out of bounds for cap_initial_normals_ (ncols=" + 
                                    std::to_string(cap_initial_normals_.ncols()) + ").");
        }
        
        // Get initial normal for this element
        Vector<double> n0(nsd);
        for (int i = 0; i < nsd; i++) {
            n0(i) = cap_initial_normals_(i, e);
        }
        
        // Normalize initial normal (in case it's not normalized)
        double n0_norm = sqrt(utils::norm(n0));
        if (!utils::is_zero(n0_norm)) {
            n0 = n0 / n0_norm;
            
            // Check if calculated normal and initial normal point in opposite directions
            // Dot product < 0 means they point in opposite directions
            double dot_product = 0.0;
            for (int i = 0; i < nsd; i++) {
                dot_product += n(i) * n0(i);
            }
            
            // If dot product is negative, flip the calculated normal
            if (dot_product < 0.0) {
                n = -n;
            }
        }
    }
    
    return std::make_pair(Jac, n);
}

double CoupledBoundaryCondition::integrate_over_cap(ComMod& com_mod, const CmMod& cm_mod, const Array<double>& s, 
                                                   int l, std::optional<int> u, consts::MechanicalConfigurationType cfg)
{
    using namespace consts;
    auto& cm = com_mod.cm;
    int nsd = com_mod.nsd;
    int insd = nsd - 1;
    int u_val = u.has_value() ? u.value() : l;
    bool is_scalar = (u_val == l);
    const int s_comps = is_scalar ? 1 : nsd;

    auto integrate_kernel = [&](const faceType* cap_face, auto&& get_xl, auto&& get_value) -> double
    {
        double result = 0.0;
        auto compute_jac_n = [&](const Array<double>& xl, int e, int g, int l_nsd, int l_insd) {
            return compute_cap_jacobian_and_normal(xl, e, g, l_nsd, l_insd);
        };
        auto on_gauss = [&](int e, int g, double Jac, const Vector<double>& n) {
                double sHat = 0.0;
                for (int a = 0; a < cap_face->eNoN; a++) {
                    const double Na = cap_face->N(a, g);
                    if (is_scalar) {
                        sHat += Na * get_value(e, a, 0);
                    } else {
                        for (int i = 0; i < s_comps; i++) {
                            sHat += Na * get_value(e, a, i) * n(i);
                        }
                    }
                }
                result += cap_face->w(g) * Jac * sHat;
        };
        for_each_cap_gauss_point(cap_face, nsd, insd, get_xl, compute_jac_n, on_gauss);
        return result;
    };

    // Serial: only master has cap; early return if this rank has no cap
    if (cm.seq()) {
        if (!cap_face_ready()) {
            return 0.0;
        }
        faceType* cap_face = cap_face_.get();
        // initialize_cap_integration() must have been called first (e.g. from baf_ini).
        auto get_xl = [&](int e) {
            return update_cap_element_position(com_mod, e, cap_gnNo_to_tnNo_, cfg);
        };
        auto get_value = [&](int e, int a, int comp) -> double {
            int gnNo_idx = cap_face->IEN(a, e);
            auto it = cap_gnNo_to_tnNo_.find(gnNo_idx);
            if (it == cap_gnNo_to_tnNo_.end()) {
                throw std::runtime_error("[CoupledBoundaryCondition::integrate_over_cap] IEN entry references gnNo not found in cap_gnNo_to_tnNo_.");
            }
            int Ac = it->second;
            if (Ac < 0 || Ac >= com_mod.tnNo) {
                throw std::runtime_error("[CoupledBoundaryCondition::integrate_over_cap] Invalid tnNo index while integrating cap.");
            }
            if (l + comp >= s.nrows() || Ac >= s.ncols()) {
                throw std::runtime_error("[CoupledBoundaryCondition::integrate_over_cap] Array bounds exceeded while integrating cap.");
            }
            return s(l + comp, Ac);
        };
        double result = integrate_kernel(cap_face, get_xl, get_value);
        return result;
    }

    // Parallel: use pre-gathered data (caller must call prepare_cap_gathered_data first). Only master has cap and does the integration.
    if (cap_nNo_gathered_ == 0) {
        return 0.0;
    }
    const int cap_nNo = cap_nNo_gathered_;
    const Array<double>& cap_s_use = (cfg == MechanicalConfigurationType::new_timestep) ? cap_Yn_gathered_ : cap_Yo_gathered_;
    double result = 0.0;
    if (cm.idcm() == cm_mod.master) {
        faceType* cap_face = cap_face_.get();
        std::unordered_map<int, int> gnNo_to_capIdx;
        gnNo_to_capIdx.reserve(cap_nNo);
        for (int a = 0; a < cap_nNo; a++) {
            gnNo_to_capIdx[cap_face->gN(a)] = a;
        }
        auto get_xl = [&](int e) {
            return update_cap_element_position(e, cfg, cap_x_gathered_, cap_Do_gathered_, cap_Dn_gathered_, gnNo_to_capIdx);
        };
        auto get_value = [&](int e, int a, int comp) -> double {
            int gnNo_idx = cap_face->IEN(a, e);
            auto it = gnNo_to_capIdx.find(gnNo_idx);
            if (it == gnNo_to_capIdx.end()) {
                throw std::runtime_error("[CoupledBoundaryCondition::integrate_over_cap] IEN entry references gnNo not found in gathered cap map.");
            }
            int cap_idx = it->second;
            if (comp >= cap_s_use.nrows() || cap_idx >= cap_s_use.ncols()) {
                throw std::runtime_error("[CoupledBoundaryCondition::integrate_over_cap] Gathered array bounds exceeded while integrating cap.");
            }
            return cap_s_use(comp, cap_idx);
        };
        result = integrate_kernel(cap_face, get_xl, get_value);
    }
    cm.bcast(cm_mod, &result);
    return result;
}

// =========================================================================
// Cap valM computation (precomputed normal integrals)
// =========================================================================

void CoupledBoundaryCondition::compute_cap_valM(ComMod& com_mod, const CmMod& cm_mod, consts::MechanicalConfigurationType cfg)
{
    using namespace consts;
    
    if (!cap_face_ready()) {
        cap_valM_.resize(0, 0);
        return;
    }
    faceType* cap_face = cap_face_.get();
    // initialize_cap_integration() must have been called first (e.g. from baf_ini).
    int nsd = com_mod.nsd;
    int cap_nNo = cap_face->nNo;
    
    // Initialize cap_valM_ to zero: nsd x cap_nNo (stored with cap face-local indices, like face.valM)
    cap_valM_.resize(nsd, cap_nNo);
    cap_valM_ = 0.0;
    
    // Create mapping from gnNo to cap face-local index for efficient lookup
    std::unordered_map<int, int> gnNo_to_cap_local;
    for (int a = 0; a < cap_nNo; a++) {
        int gnNo = cap_face->gN(a);
        gnNo_to_cap_local[gnNo] = a;
    }

    auto& cm = com_mod.cm;
    const bool use_gathered = !cm.seq() && cap_nNo_gathered_ > 0;

    auto get_xl = [&](int e) {
        if (use_gathered) {
            return update_cap_element_position(e, cfg, cap_x_gathered_, cap_Do_gathered_, cap_Dn_gathered_, gnNo_to_cap_local);
        }
        return update_cap_element_position(com_mod, e, cap_gnNo_to_tnNo_, cfg);
    };
    auto compute_jac_n = [&](const Array<double>& xl, int e, int g, int l_nsd, int l_insd) {
        return compute_cap_jacobian_and_normal(xl, e, g, l_nsd, l_insd);
    };
    auto on_gauss = [&](int e, int g, double Jac, const Vector<double>& n) {
        // Accumulate ∫ N_A n_i dΓ at each node
        for (int a = 0; a < cap_face->eNoN; a++) {
            int gnNo_idx = cap_face->IEN(a, e);
            auto it = gnNo_to_cap_local.find(gnNo_idx);
            if (it == gnNo_to_cap_local.end()) {
                throw std::runtime_error("[CoupledBoundaryCondition::compute_cap_valM] IEN entry (element " +
                                        std::to_string(e) + ", node " + std::to_string(a) +
                                        ") contains invalid gnNo index " + std::to_string(gnNo_idx) +
                                        " not found in cap face nodes.");
            }
            int cap_a = it->second;
            if (cap_a < 0 || cap_a >= cap_nNo) {
                throw std::runtime_error("[CoupledBoundaryCondition::compute_cap_valM] Invalid cap face-local index cap_a=" +
                                        std::to_string(cap_a) + " (cap_nNo=" + std::to_string(cap_nNo) + ")");
            }
            for (int i = 0; i < nsd; i++) {
                cap_valM_(i, cap_a) += cap_face->N(a, g) * cap_face->w(g) * Jac * n(i);
            }
        }
    };
    for_each_cap_gauss_point(cap_face, nsd, nsd - 1, get_xl, compute_jac_n, on_gauss);

    // Reduce across processors only when using local data (each rank had partial cap_gnNo_to_tnNo_); when using gathered data only master computed
    if (!cm.seq() && !use_gathered) {
        for (int i = 0; i < nsd; i++) {
            Vector<double> row = cap_valM_.row(i);
            row = cm.reduce(cm_mod, row);
            cap_valM_.set_row(i, row);
        }
    }
}

// =========================================================================
// Copy cap data to linear solver face structure
// =========================================================================

void CoupledBoundaryCondition::copy_cap_data_to_linear_solver_face(ComMod& com_mod, const CmMod& cm_mod, 
                                                                   fsi_linear_solver::FSILS_faceType& lhs_face,
                                                                   consts::MechanicalConfigurationType cfg)
{
    lhs_face.has_cap = has_cap();
    if (!has_cap()) {
        lhs_face.cap_val.resize(0, 0);
        lhs_face.cap_valM.resize(0, 0);
        lhs_face.cap_glob.resize(0);
        return;
    }

    const int nsd = com_mod.nsd;
    const bool serial = com_mod.cm.seq();

    if (serial) {
        // Serial: only this rank has cap_face_; fill lhs_face directly.
        if (!cap_face_ready()) {
            lhs_face.cap_val.resize(0, 0);
            lhs_face.cap_valM.resize(0, 0);
            lhs_face.cap_glob.resize(0);
            return;
        }
        compute_cap_valM(com_mod, cm_mod, cfg);
        const faceType* cap_face = cap_face_.get();
        int cap_nNo = cap_face->nNo;
        lhs_face.cap_val = cap_valM_;
        lhs_face.cap_valM.resize(nsd, cap_nNo);
        lhs_face.cap_valM = 0.0;
        lhs_face.cap_glob.resize(cap_nNo);
        for (int a = 0; a < cap_nNo; a++) {
            int gnNo = cap_face->gN(a);
            int localIdx = -1;
            for (int i = 0; i < com_mod.tnNo; i++) {
                if (com_mod.ltg(i) == gnNo) {
                    localIdx = i;
                    break;
                }
            }
            lhs_face.cap_glob(a) = (localIdx >= 0) ? com_mod.lhs.map(localIdx) : -1;
        }
        return;
    }

    // Parallel: broadcast full cap data from master, then each rank builds local subset.
    int cap_nNo = 0;
    Vector<int> cap_gN_all;
    Array<double> cap_val_all;

    if (cap_face_ready()) {
        compute_cap_valM(com_mod, cm_mod, cfg);
        const faceType* cap_face = cap_face_.get();
        if (cap_face->nNo > 0) {
            cap_nNo = cap_face->nNo;
            cap_gN_all.resize(cap_nNo);
            cap_val_all.resize(nsd, cap_nNo);
            for (int a = 0; a < cap_nNo; a++) {
                cap_gN_all(a) = cap_face->gN(a);
                for (int i = 0; i < nsd; i++)
                    cap_val_all(i, a) = cap_valM_(i, a);
            }
        }
    }

    com_mod.cm.bcast(cm_mod, &cap_nNo);
    if (cap_nNo == 0) {
        lhs_face.cap_val.resize(0, 0);
        lhs_face.cap_valM.resize(0, 0);
        lhs_face.cap_glob.resize(0);
        return;
    }
    // Only resize on ranks that don't have the data; Vector::resize() deallocates and
    // zero-initializes, which would wipe the sender's cap_gN_all/cap_val_all.
    const bool i_am_sender = cap_face_ready();
    if (!i_am_sender) {
        cap_gN_all.resize(cap_nNo);
        cap_val_all.resize(nsd, cap_nNo);
    }
    com_mod.cm.bcast(cm_mod, cap_gN_all);
    com_mod.cm.bcast(cm_mod, cap_val_all);

    int n_owned = 0;
    for (int a = 0; a < cap_nNo; a++) {
        int gnNo = cap_gN_all(a);
        for (int i = 0; i < com_mod.tnNo; i++) {
            if (com_mod.ltg(i) == gnNo) {
                n_owned++;
                break;
            }
        }
    }
    lhs_face.cap_glob.resize(n_owned);
    lhs_face.cap_val.resize(nsd, n_owned);
    lhs_face.cap_valM.resize(nsd, n_owned);
    lhs_face.cap_valM = 0.0;

    int idx = 0;
    for (int a = 0; a < cap_nNo; a++) {
        int gnNo = cap_gN_all(a);
        int localIdx = -1;
        for (int i = 0; i < com_mod.tnNo; i++) {
            if (com_mod.ltg(i) == gnNo) {
                localIdx = i;
                break;
            }
        }
        if (localIdx >= 0) {
            lhs_face.cap_glob(idx) = com_mod.lhs.map(localIdx);
            for (int i = 0; i < nsd; i++)
                lhs_face.cap_val(i, idx) = cap_val_all(i, a);
            idx++;
        }
    }
}
