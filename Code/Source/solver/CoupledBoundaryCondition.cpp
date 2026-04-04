// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "ComMod.h"
#include "SimulationLogger.h"
#include "CoupledBoundaryCondition.h"
#include "all_fun.h"
#include "consts.h"
#include "utils.h"
#include "fils_struct.hpp"
#include "VtkData.h"
#include "nn.h"
#include <fstream>
#include <unordered_set>
#include <vtkCellType.h>


// =========================================================================
// CoupledBoundaryCondition
// =========================================================================

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
    , phys_(other.phys_)
    , flowrate_cfg_o_(other.flowrate_cfg_o_)
    , flowrate_cfg_n_(other.flowrate_cfg_n_)
    , has_cap_(other.has_cap_)
    , cap_(other.cap_)
{
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
        phys_ = other.phys_;
        flowrate_cfg_o_ = other.flowrate_cfg_o_;
        flowrate_cfg_n_ = other.flowrate_cfg_n_;
        has_cap_ = other.has_cap_;
        cap_ = other.cap_;
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
    , phys_(other.phys_)
    , flowrate_cfg_o_(other.flowrate_cfg_o_)
    , flowrate_cfg_n_(other.flowrate_cfg_n_)
    , has_cap_(other.has_cap_)
    , cap_(std::move(other.cap_))
{
    other.face_ = nullptr;
    other.logger_ = nullptr;
    other.has_cap_ = false;
    other.flow_sol_id_ = -1;
    other.pressure_sol_id_ = -1;
    other.Qo_ = 0.0;
    other.Qn_ = 0.0;
    other.Po_ = 0.0;
    other.Pn_ = 0.0;
    other.pressure_ = 0.0;
    other.phys_ = consts::EquationType::phys_NA;
    other.flowrate_cfg_o_ = consts::MechanicalConfigurationType::reference;
    other.flowrate_cfg_n_ = consts::MechanicalConfigurationType::reference;
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
        phys_ = other.phys_;
        flowrate_cfg_o_ = other.flowrate_cfg_o_;
        flowrate_cfg_n_ = other.flowrate_cfg_n_;
        has_cap_ = other.has_cap_;
        cap_ = std::move(other.cap_);

        other.face_ = nullptr;
        other.logger_ = nullptr;
        other.has_cap_ = false;
        other.flow_sol_id_ = -1;
        other.pressure_sol_id_ = -1;
        other.Qo_ = 0.0;
        other.Qn_ = 0.0;
        other.Po_ = 0.0;
        other.Pn_ = 0.0;
        other.pressure_ = 0.0;
        other.phys_ = consts::EquationType::phys_NA;
        other.flowrate_cfg_o_ = consts::MechanicalConfigurationType::reference;
        other.flowrate_cfg_n_ = consts::MechanicalConfigurationType::reference;
    }
    return *this;
}

CoupledBoundaryCondition::CoupledBoundaryCondition(consts::BoundaryConditionType bc_type, const faceType& face, const std::string& face_name,
                                               const std::string& block_name, consts::EquationType phys, bool follower_pressure_load,
                                               SimulationLogger& logger)
    : face_(&face)
    , logger_(&logger)
    , bc_type_(bc_type)
    , block_name_(block_name)
    , face_name_(face_name)
    , phys_(phys)
{
    set_flowrate_mechanical_configurations(phys, follower_pressure_load);
}

CoupledBoundaryCondition::CoupledBoundaryCondition(consts::BoundaryConditionType bc_type, const faceType& face, const std::string& face_name,
                                               const std::string& block_name, const std::string& cap_face_vtp_file,
                                               consts::EquationType phys, bool follower_pressure_load, SimulationLogger& logger)
    : cap_face_vtp_file_(cap_face_vtp_file)
    , face_(&face)
    , logger_(&logger)
    , bc_type_(bc_type)
    , block_name_(block_name)
    , face_name_(face_name)
    , phys_(phys)
{
    // Load the cap VTP file if provided
    if (!cap_face_vtp_file_.empty()) {
        load_cap_face_vtp(cap_face_vtp_file_);
    }
    set_flowrate_mechanical_configurations(phys, follower_pressure_load);
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

void CoupledBoundaryCondition::set_flowrate_mechanical_configurations(consts::EquationType phys, bool follower_pressure_load)
{
    using namespace consts;
    phys_ = phys;
    follower_pressure_load_ = follower_pressure_load;
    if ((phys == EquationType::phys_struct) || (phys == EquationType::phys_ustruct)) {
        flowrate_cfg_o_ = MechanicalConfigurationType::old_timestep;
        flowrate_cfg_n_ = MechanicalConfigurationType::new_timestep;
    } else {
        flowrate_cfg_o_ = MechanicalConfigurationType::reference;
        flowrate_cfg_n_ = MechanicalConfigurationType::reference;
    }
}

/// @brief Compute flowrates at the boundary face at old and new timesteps
///
/// This replicates the flowrate computation done in set_bc::calc_der_cpl_bc and
/// set_bc::set_bc_cpl for coupled Neumann boundary conditions.
///
/// The flowrate is computed as the integral of velocity dotted with the face normal.
/// For struct/ustruct physics, the integral is computed on the deformed configuration.
/// For fluid/FSI/CMM physics, the integral is computed on the reference configuration.
void CoupledBoundaryCondition::compute_flowrates(ComMod& com_mod, const CmMod& cm_mod)
{
    if (face_ == nullptr) {
        throw std::runtime_error("[CoupledBoundaryCondition::compute_flowrates] Face is not set.");
    }
    
    int nsd = com_mod.nsd;
    
    Qo_ = all_fun::integ(com_mod, cm_mod, *face_, com_mod.Yo, 0, nsd-1, false, flowrate_cfg_o_);
    Qn_ = all_fun::integ(com_mod, cm_mod, *face_, com_mod.Yn, 0, nsd-1, false, flowrate_cfg_n_);
    
    if (has_cap_) {
        const auto [Qo_cap, Qn_cap] =
            calculate_cap_contribution(com_mod, cm_mod, nsd, flowrate_cfg_o_, flowrate_cfg_n_);
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

    // In the constructor, the face pointer is set to the global face, which is then partitioned and distributed among all processes.
    // Here, we update the face pointer to the local face.
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

    int phys_int = static_cast<int>(phys_);
    cm.bcast(cm_mod, &phys_int);
    int cfg_o_int = static_cast<int>(flowrate_cfg_o_);
    int cfg_n_int = static_cast<int>(flowrate_cfg_n_);
    cm.bcast(cm_mod, &cfg_o_int);
    cm.bcast(cm_mod, &cfg_n_int);
    if (is_slave) {
        phys_ = static_cast<consts::EquationType>(phys_int);
        flowrate_cfg_o_ = static_cast<consts::MechanicalConfigurationType>(cfg_o_int);
        flowrate_cfg_n_ = static_cast<consts::MechanicalConfigurationType>(cfg_n_int);
    }
    
    // Distribute solution IDs
    cm.bcast(cm_mod, &flow_sol_id_);
    cm.bcast(cm_mod, &pressure_sol_id_);
    cm.bcast(cm_mod, &in_out_sign_);

    // Distribute cap flag so all ranks agree (master loaded cap_face_; slaves have has_cap_ set here)
    cm.bcast(cm_mod, &has_cap_);

    // Keep a cap object on every rank when has_cap_ is true so collective
    // cap gather/integration preparation can be called everywhere. Only
    // ranks with loaded geometry have cap_->is_ready() == true.
    if (has_cap_) {
        if (!cap_.has_value()) {
            cap_.emplace();
        }
    } else {
        cap_.reset();
    }
}

const faceType* CoupledBoundaryCondition::get_face() const
{
    return face_;
}

bool CoupledBoundaryCondition::is_initialized() const
{
    return (face_ != nullptr);
}

bool CoupledBoundaryCondition::cap_face_ready() const
{
    return has_cap_ && cap_.has_value() && cap_->is_ready();
}


// =========================================================================
// Cap surface loading and integration
// =========================================================================

void CoupledBoundaryCondition::load_cap_face_vtp(const std::string& vtp_file_path)
{
    if (face_ == nullptr) {
        throw std::runtime_error("[CoupledBoundaryCondition::load_cap_face_vtp] Cannot load cap: face_ is null (BC not constructed with a face or distribute() not run).");
    }

    cap_face_vtp_file_ = vtp_file_path;
    has_cap_ = false;
    cap_.reset();

    if (vtp_file_path.empty()) {
        return;
    }

    cap_.emplace();
    try {
        cap_->load_from_vtp(vtp_file_path, *face_, face_name_);
    } catch (...) {
        cap_.reset();
        has_cap_ = false;
        throw;
    }
    has_cap_ = true;
}

std::pair<double, double> CoupledBoundaryCondition::calculate_cap_contribution(ComMod& com_mod, const CmMod& cm_mod,
    int nsd, consts::MechanicalConfigurationType cfg_o,
    consts::MechanicalConfigurationType cfg_n)
{
    auto& cm = com_mod.cm;
    double Qo_cap = 0.0;
    double Qn_cap = 0.0;
    const bool i_am_master = cm.mas(cm_mod);
    const bool serial_run = (cm.np() == 1);

    // Gather Yo and Yn on the cap surface
    if (!cap_) {
        throw std::runtime_error("[CoupledBoundaryCondition::calculate_cap_contribution] Cap is enabled (has_cap) but no CappingSurface on this rank.");
    }
    cap_->prepare_gathered_data(com_mod, cm_mod, com_mod.Yo, com_mod.Yn, 0, nsd, cfg_o, cfg_n);

    // Integrate over the cap surface
    if (serial_run || i_am_master) {
        if (cap_ && cap_->is_ready()) {
            Qo_cap = cap_->integrate_over(com_mod, cm_mod, com_mod.Yo, 0, nsd - 1, cfg_o);
            Qn_cap = cap_->integrate_over(com_mod, cm_mod, com_mod.Yn, 0, nsd - 1, cfg_n);
        }
    }

    // Broadcast Qo_cap and Qn_cap to all ranks
    if (!serial_run) {
        cm.bcast(cm_mod, &Qo_cap);
        cm.bcast(cm_mod, &Qn_cap);
    }
    return {Qo_cap, Qn_cap};
}

// =========================================================================
// CappingSurface (definitions; class at bottom of CoupledBoundaryCondition.h)
// =========================================================================

namespace {

consts::ElementType vtk_cell_type_to_element_type(int vtk_cell_type)
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
            throw std::runtime_error("[CappingSurface] Unsupported VTK cell type " +
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

} // namespace

CappingSurface::CappingSurface(const CappingSurface& other)
    : global_node_ids_(other.global_node_ids_)
    , area_computed_(other.area_computed_)
    , gnNo_to_tnNo_(other.gnNo_to_tnNo_)
    , valM_(other.valM_)
    , initial_normals_(other.initial_normals_)
    , x_gathered_(other.x_gathered_)
    , Do_gathered_(other.Do_gathered_)
    , Dn_gathered_(other.Dn_gathered_)
    , Yo_gathered_(other.Yo_gathered_)
    , Yn_gathered_(other.Yn_gathered_)
    , nNo_gathered_(other.nNo_gathered_)
    , gN_broadcast_(other.gN_broadcast_)
{
    if (other.face_ != nullptr) {
        try {
            face_ = std::make_unique<faceType>(*other.face_);
            if (face_ != nullptr) {
                if (face_->nNo > 0 && face_->gN.size() != face_->nNo) {
                    throw std::runtime_error("[CappingSurface::copy constructor] Invalid face_: gN.size()=" +
                                            std::to_string(face_->gN.size()) + " != nNo=" + std::to_string(face_->nNo));
                }
                if (face_->nEl > 0 && face_->IEN.ncols() != face_->nEl) {
                    throw std::runtime_error("[CappingSurface::copy constructor] Invalid face_: IEN.ncols()=" +
                                            std::to_string(face_->IEN.ncols()) + " != nEl=" + std::to_string(face_->nEl));
                }
            }
        } catch (const std::exception& e) {
            throw std::runtime_error("[CappingSurface::copy constructor] Failed to copy face_: " + std::string(e.what()));
        }
    } else {
        face_.reset();
    }
}

CappingSurface& CappingSurface::operator=(const CappingSurface& other)
{
    if (this != &other) {
        global_node_ids_ = other.global_node_ids_;
        area_computed_ = other.area_computed_;
        gnNo_to_tnNo_ = other.gnNo_to_tnNo_;
        valM_ = other.valM_;
        initial_normals_ = other.initial_normals_;
        x_gathered_ = other.x_gathered_;
        Do_gathered_ = other.Do_gathered_;
        Dn_gathered_ = other.Dn_gathered_;
        Yo_gathered_ = other.Yo_gathered_;
        Yn_gathered_ = other.Yn_gathered_;
        nNo_gathered_ = other.nNo_gathered_;
        gN_broadcast_ = other.gN_broadcast_;
        if (other.face_ != nullptr) {
            try {
                face_ = std::make_unique<faceType>(*other.face_);
            } catch (const std::exception& e) {
                throw std::runtime_error("[CappingSurface::operator=] Failed to copy face_: " + std::string(e.what()));
            }
        } else {
            face_.reset();
        }
    }
    return *this;
}

void CappingSurface::load_from_vtp(const std::string& vtp_file_path, const faceType& coupled_face,
                                   const std::string& coupled_face_name)
{
    face_.reset();
    gnNo_to_tnNo_.clear();
    area_computed_ = false;
    initial_normals_.resize(0, 0);

    if (vtp_file_path.empty()) {
        return;
    }

    std::ifstream file_check(vtp_file_path);
    if (!file_check.good()) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Cannot open cap VTP file '" + vtp_file_path +
                                "' for reading.");
    }
    file_check.close();

    VtkVtpData vtp_data;
    try {
        vtp_data = VtkVtpData(vtp_file_path, true);
    } catch (const std::exception& e) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Failed to construct VtkVtpData from file '" +
                                vtp_file_path + "': " + e.what());
    } catch (...) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Unknown error constructing VtkVtpData from file '" +
                                vtp_file_path + "'. This may indicate a crash in the VTK library.");
    }

    int nNo = 0;
    try {
        nNo = vtp_data.num_points();
    } catch (const std::exception& e) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Failed to get number of points from VTP file '" +
                                vtp_file_path + "': " + e.what());
    } catch (...) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Unknown error getting number of points from file '" +
                                vtp_file_path + "'");
    }
    if (nNo == 0) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Cap VTP file '" + vtp_file_path +
                                "' does not contain any points.");
    }

    int num_elems = 0;
    try {
        num_elems = vtp_data.num_elems();
    } catch (const std::exception& e) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Failed to get number of elements from VTP file '" +
                                vtp_file_path + "': " + e.what());
    } catch (...) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Unknown error getting number of elements from file '" +
                                vtp_file_path + "'");
    }
    if (num_elems == 0) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Cap VTP file '" + vtp_file_path +
                                "' does not contain any elements.");
    }

    bool has_global_node_id = false;
    try {
        has_global_node_id = vtp_data.has_point_data("GlobalNodeID");
    } catch (const std::exception& e) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Failed to check for GlobalNodeID in VTP file '" +
                                vtp_file_path + "': " + e.what());
    } catch (...) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Unknown error checking for GlobalNodeID in file '" +
                                vtp_file_path + "'");
    }
    if (!has_global_node_id) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Cap VTP file '" + vtp_file_path +
                                "' does not contain 'GlobalNodeID' point data.");
    }

    try {
        global_node_ids_.resize(nNo);
        vtp_data.copy_point_data("GlobalNodeID", global_node_ids_);
    } catch (const std::exception& e) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Failed to copy GlobalNodeID from VTP file '" +
                                vtp_file_path + "': " + e.what());
    } catch (...) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Unknown error copying GlobalNodeID from file '" +
                                vtp_file_path + "'");
    }

    Array<int> conn;
    int eNoN = 0;
    int vtk_cell_type = 0;
    try {
        conn = vtp_data.get_connectivity();
        eNoN = vtp_data.np_elem();
        vtk_cell_type = vtp_data.elem_type();
    } catch (const std::exception& e) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Failed to get connectivity from VTP file '" +
                                vtp_file_path + "': " + e.what());
    } catch (...) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Unknown error getting connectivity from file '" +
                                vtp_file_path + "'");
    }
    if (eNoN <= 0) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Invalid number of nodes per element: " +
                                std::to_string(eNoN));
    }

    consts::ElementType eType = vtk_cell_type_to_element_type(vtk_cell_type);

    face_ = std::make_unique<faceType>();
    face_->name = coupled_face_name + "_cap";
    face_->iM = coupled_face.iM;
    face_->nNo = nNo;
    face_->nEl = num_elems;
    face_->gnEl = num_elems;
    face_->eNoN = eNoN;
    face_->eType = eType;

    if (global_node_ids_.size() != nNo) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] global_node_ids_ size mismatch: " +
                                std::to_string(global_node_ids_.size()) + " != " + std::to_string(nNo));
    }
    face_->gN.resize(nNo);
    for (int a = 0; a < nNo; a++) {
        face_->gN(a) = global_node_ids_(a) - 1;
    }

    face_->IEN.resize(eNoN, num_elems);

    for (int e = 0; e < num_elems; e++) {
        for (int a = 0; a < eNoN; a++) {
            int local_node_idx = conn(a, e);
            if (local_node_idx < 0 || local_node_idx >= nNo) {
                throw std::runtime_error("[CappingSurface::load_from_vtp] Invalid local node index " +
                                        std::to_string(local_node_idx) + " in cap connectivity (element " +
                                        std::to_string(e) + ", node " + std::to_string(a) + ", nNo=" + std::to_string(nNo) + ").");
            }
            face_->IEN(a, e) = face_->gN(local_node_idx);
        }
    }

    try {
        bool has_normals = vtp_data.has_cell_data("Normals");
        if (has_normals) {
            auto [num_comp, num_tuples] = vtp_data.get_cell_data_dimensions("Normals");

            if (num_comp == 0 || num_tuples == 0) {
                std::string error_msg = "[CappingSurface::load_from_vtp] Normals array exists but has zero size. ";
                error_msg += "num_components=" + std::to_string(num_comp) + ", num_tuples=" + std::to_string(num_tuples);
                error_msg += ", expected num_elems=" + std::to_string(num_elems);
                error_msg += ". The array may not be a numeric type or may be empty.";
                throw std::runtime_error(error_msg);
            }

            if (num_tuples != num_elems) {
                throw std::runtime_error("[CappingSurface::load_from_vtp] Normals array size mismatch: " +
                                        std::to_string(num_tuples) + " != " + std::to_string(num_elems));
            }

            if (num_comp != 2 && num_comp != 3) {
                throw std::runtime_error("[CappingSurface::load_from_vtp] Invalid number of components in Normals array: " +
                                        std::to_string(num_comp) + " (expected 2 or 3)");
            }

            initial_normals_.resize(num_comp, num_elems);
            vtp_data.copy_cell_data("Normals", initial_normals_);

            if (initial_normals_.nrows() != num_comp || initial_normals_.ncols() != num_elems) {
                throw std::runtime_error("[CappingSurface::load_from_vtp] Failed to copy Normals data. "
                                        "Expected size: " + std::to_string(num_comp) + "x" + std::to_string(num_elems) +
                                        ", Actual size: " + std::to_string(initial_normals_.nrows()) + "x" +
                                        std::to_string(initial_normals_.ncols()));
            }
        } else {
            initial_normals_.resize(0, 0);
        }
    } catch (const std::exception& e) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Failed to load Normals from VTP file '" +
                                vtp_file_path + "': " + e.what());
    } catch (...) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Unknown error loading Normals from file '" +
                                vtp_file_path + "'");
    }

    std::unordered_set<int> face_gn;
    face_gn.reserve(static_cast<size_t>(coupled_face.nNo));
    for (int a = 0; a < coupled_face.nNo; a++) {
        face_gn.insert(coupled_face.gN(a));
    }
    bool any_shared = false;
    for (int a = 0; a < face_->nNo; a++) {
        if (face_gn.find(face_->gN(a)) != face_gn.end()) {
            any_shared = true;
            break;
        }
    }
    if (!any_shared) {
        throw std::runtime_error(
            "[CappingSurface::load_from_vtp] Cap VTP file '" + vtp_file_path +
            "' has no GlobalNodeID entries in common with coupled face '" + coupled_face_name +
            "'. The cap must share at least one mesh node with that face.");
    }
}

void CappingSurface::initialize_integration(ComMod& com_mod, const CmMod& cm_mod)
{
    (void)cm_mod;
    using namespace consts;

    if (!is_ready()) {
        return;
    }
    if (face_->nEl == 0 || face_->nNo == 0) {
        throw std::runtime_error("[CappingSurface::initialize_integration] Cap face is not properly initialized.");
    }

    if (face_->nG > 0 && face_->w.size() > 0 && face_->N.nrows() > 0 &&
        face_->Nx.nslices() == face_->nG && !gnNo_to_tnNo_.empty()) {
        return;
    }

    int nsd = com_mod.nsd;
    int insd = nsd - 1;

    auto& msh = com_mod.msh[face_->iM];
    if (com_mod.ltg.size() == 0 || com_mod.tnNo == 0) {
        return;
    }
    if (com_mod.ltg.size() != com_mod.tnNo) {
        throw std::runtime_error("[CappingSurface::initialize_integration] com_mod.ltg size (" +
                                std::to_string(com_mod.ltg.size()) + ") does not match com_mod.tnNo (" +
                                std::to_string(com_mod.tnNo) + ")");
    }

    gnNo_to_tnNo_.clear();
    gnNo_to_tnNo_.reserve(com_mod.tnNo);
    for (int a = 0; a < com_mod.tnNo; a++) {
        int gnNo_idx = com_mod.ltg(a);
        if (gnNo_idx < 0 || gnNo_idx >= msh.gnNo) {
            throw std::runtime_error("[CappingSurface::initialize_integration] Invalid gnNo index " +
                                    std::to_string(gnNo_idx) + " in com_mod.ltg at position " +
                                    std::to_string(a) + " (msh.gnNo=" + std::to_string(msh.gnNo) + ")");
        }
        gnNo_to_tnNo_[gnNo_idx] = a;
    }

    try {
        if (face_->eType == ElementType::TRI3) {
            face_->nG = 1;
        } else if (face_->eType == ElementType::QUD4) {
            face_->nG = 4;
        } else if (face_->eType == ElementType::TRI6) {
            face_->nG = 3;
        } else {
            face_->nG = face_->eNoN;
        }

        face_->w.resize(face_->nG);
        face_->xi.resize(insd, face_->nG);

        if (face_->eType == ElementType::TRI3 && face_->nG == 1) {
            face_->w(0) = 0.5;
            face_->xi(0, 0) = 1.0 / 3.0;
            face_->xi(1, 0) = 1.0 / 3.0;
        } else {
            nn::get_gip(insd, face_->eType, face_->nG, face_->w, face_->xi);
        }

        face_->N.resize(face_->eNoN, face_->nG);
        face_->Nx.resize(insd, face_->eNoN, face_->nG);
        for (int g = 0; g < face_->nG; g++) {
            nn::get_gnn(insd, face_->eType, face_->eNoN, g, face_->xi, face_->N, face_->Nx);
        }
    } catch (const std::exception& e) {
        throw std::runtime_error("[CappingSurface::initialize_integration] Failed to initialize cap face shape functions: " +
                                std::string(e.what()));
    }
}

Array<double> CappingSurface::update_element_position(ComMod& com_mod, int e,
                                                      const std::unordered_map<int, int>& gnNo_to_tnNo,
                                                      consts::MechanicalConfigurationType cfg)
{
    using namespace consts;

    if (face_->IEN.nrows() == 0 || face_->IEN.ncols() == 0) {
        throw std::runtime_error("[CappingSurface::update_element_position] face_->IEN is not allocated.");
    }
    if (e < 0 || e >= face_->IEN.ncols()) {
        throw std::runtime_error("[CappingSurface::update_element_position] Element index e=" +
                                std::to_string(e) + " is out of bounds (IEN.ncols()=" + std::to_string(face_->IEN.ncols()) + ").");
    }
    if (face_->eNoN <= 0) {
        throw std::runtime_error("[CappingSurface::update_element_position] face_->eNoN is invalid: " +
                                std::to_string(face_->eNoN));
    }

    int nsd = com_mod.nsd;
    Array<double> xl(nsd, face_->eNoN);

    for (int a = 0; a < face_->eNoN; a++) {
        int gnNo_idx = face_->IEN(a, e);
        auto it = gnNo_to_tnNo.find(gnNo_idx);
        if (it == gnNo_to_tnNo.end()) {
            throw std::runtime_error("[CappingSurface::update_element_position] IEN entry (element " +
                                    std::to_string(e) + ", node " + std::to_string(a) +
                                    ") contains invalid gnNo index " + std::to_string(gnNo_idx) +
                                    " not found in com_mod.ltg mapping.");
        }
        int Ac = it->second;

        if (Ac < 0 || Ac >= com_mod.tnNo) {
            std::string msg = "[CappingSurface::update_element_position] Invalid node index Ac=" +
                             std::to_string(Ac) + " (tnNo=" + std::to_string(com_mod.tnNo) +
                             ") at element " + std::to_string(e) + ", local node " + std::to_string(a);
            throw std::runtime_error(msg);
        }

        if (Ac >= com_mod.x.ncols() || nsd > com_mod.x.nrows()) {
            throw std::runtime_error("[CappingSurface::update_element_position] Invalid bounds for com_mod.x: Ac=" +
                                    std::to_string(Ac) + ", x.ncols()=" + std::to_string(com_mod.x.ncols()) +
                                    ", nsd=" + std::to_string(nsd) + ", x.nrows()=" + std::to_string(com_mod.x.nrows()));
        }

        xl.set_col(a, com_mod.x.col(Ac));

        if (cfg == MechanicalConfigurationType::old_timestep) {
            if (Ac >= com_mod.Do.ncols() || nsd > com_mod.Do.nrows()) {
                throw std::runtime_error("[CappingSurface::update_element_position] Invalid bounds for com_mod.Do: Ac=" +
                                        std::to_string(Ac) + ", Do.ncols()=" + std::to_string(com_mod.Do.ncols()) +
                                        ", nsd=" + std::to_string(nsd) + ", Do.nrows()=" + std::to_string(com_mod.Do.nrows()));
            }
            for (int i = 0; i < nsd; i++) {
                xl(i, a) += com_mod.Do(i, Ac);
            }
        } else if (cfg == MechanicalConfigurationType::new_timestep) {
            if (Ac >= com_mod.Dn.ncols() || nsd > com_mod.Dn.nrows()) {
                throw std::runtime_error("[CappingSurface::update_element_position] Invalid bounds for com_mod.Dn: Ac=" +
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

Array<double> CappingSurface::update_element_position(int e, consts::MechanicalConfigurationType cfg,
                                                      const Array<double>& cap_x, const Array<double>& cap_Do,
                                                      const Array<double>& cap_Dn,
                                                      const std::unordered_map<int, int>& gnNo_to_capIdx)
{
    if (face_->IEN.nrows() == 0 || face_->IEN.ncols() == 0) {
        throw std::runtime_error("[CappingSurface::update_element_position(gathered)] face_->IEN not ready.");
    }
    int nsd = cap_x.nrows();
    Array<double> xl(nsd, face_->eNoN);
    for (int a = 0; a < face_->eNoN; a++) {
        int gnNo_idx = face_->IEN(a, e);
        auto it = gnNo_to_capIdx.find(gnNo_idx);
        if (it == gnNo_to_capIdx.end()) {
            throw std::runtime_error("[CappingSurface::update_element_position(gathered)] IEN gnNo " +
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

void CappingSurface::gather_node_data_to_master(ComMod& com_mod, const CmMod& cm_mod,
                                                const Array<double>& s_old, const Array<double>& s_new,
                                                int l, int s_comps, consts::MechanicalConfigurationType cfg_o,
                                                consts::MechanicalConfigurationType cfg_n,
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
        if (!is_ready()) return;
        for (int a = 0; a < cap_nNo; a++) {
            auto it = gnNo_to_tnNo_.find(face_->gN(a));
            if (it == gnNo_to_tnNo_.end()) continue;
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
    std::unordered_map<int, int> gnNo_to_Ac;
    for (int Ac = 0; Ac < com_mod.tnNo; Ac++)
        gnNo_to_Ac[com_mod.ltg(Ac)] = Ac;

    const int per_node = 1 + 3 * nsd + 2 * s_comps;
    Vector<double> send_buf;
    int n_owned = 0;
    for (int a = 0; a < cap_nNo; a++) {
        if (gnNo_to_Ac.find(gN_broadcast_(a)) != gnNo_to_Ac.end()) n_owned++;
    }
    send_buf.resize(n_owned * per_node);
    int idx = 0;
    for (int a = 0; a < cap_nNo; a++) {
        auto it = gnNo_to_Ac.find(gN_broadcast_(a));
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

void CappingSurface::prepare_gathered_data(ComMod& com_mod, const CmMod& cm_mod,
                                           const Array<double>& Yo, const Array<double>& Yn,
                                           int l, int s_comps, consts::MechanicalConfigurationType cfg_o,
                                           consts::MechanicalConfigurationType cfg_n)
{
    auto& cm = com_mod.cm;
    if (cm.seq()) {
        if (!is_ready()) return;
        const int cap_nNo = face_->nNo;
        if (cap_nNo == 0) return;
        nNo_gathered_ = cap_nNo;
        const int nsd = com_mod.nsd;
        x_gathered_.resize(nsd, cap_nNo);
        Do_gathered_.resize(nsd, cap_nNo);
        Dn_gathered_.resize(nsd, cap_nNo);
        Yo_gathered_.resize(s_comps, cap_nNo);
        Yn_gathered_.resize(s_comps, cap_nNo);
        for (int a = 0; a < cap_nNo; a++) {
            auto it = gnNo_to_tnNo_.find(face_->gN(a));
            if (it == gnNo_to_tnNo_.end()) continue;
            int Ac = it->second;
            for (int i = 0; i < nsd; i++) {
                x_gathered_(i, a) = com_mod.x(i, Ac);
                Do_gathered_(i, a) = com_mod.Do(i, Ac);
                Dn_gathered_(i, a) = com_mod.Dn(i, Ac);
            }
            for (int i = 0; i < s_comps; i++) {
                Yo_gathered_(i, a) = Yo(l + i, Ac);
                Yn_gathered_(i, a) = Yn(l + i, Ac);
            }
        }
        return;
    }
    int cap_nNo = 0;
    if (cm.idcm() == cm_mod.master && is_ready()) {
        cap_nNo = face_->nNo;
    }
    cm.bcast(cm_mod, &cap_nNo);
    nNo_gathered_ = cap_nNo;
    if (cap_nNo == 0) {
        return;
    }
    gN_broadcast_.resize(cap_nNo);
    if (cm.idcm() == cm_mod.master && is_ready()) {
        for (int a = 0; a < cap_nNo; a++)
            gN_broadcast_(a) = face_->gN(a);
    }
    cm.bcast(cm_mod, gN_broadcast_);

    const int nsd = com_mod.nsd;
    if (cm.idcm() == cm_mod.master) {
        x_gathered_.resize(nsd, cap_nNo);
        Do_gathered_.resize(nsd, cap_nNo);
        Dn_gathered_.resize(nsd, cap_nNo);
        Yo_gathered_.resize(s_comps, cap_nNo);
        Yn_gathered_.resize(s_comps, cap_nNo);
    }
    gather_node_data_to_master(com_mod, cm_mod, Yo, Yn, l, s_comps, cfg_o, cfg_n, cap_nNo, nsd,
                              x_gathered_, Do_gathered_, Dn_gathered_, Yo_gathered_, Yn_gathered_);
}

std::pair<double, Vector<double>> CappingSurface::compute_jacobian_and_normal(const Array<double>& xl,
                                                                              int e, int g, int nsd, int insd)
{
    if (xl.nrows() != nsd || xl.ncols() != face_->eNoN) {
        throw std::runtime_error("[CappingSurface::compute_jacobian_and_normal] xl has wrong dimensions: " +
                                std::to_string(xl.nrows()) + "x" + std::to_string(xl.ncols()) +
                                " (expected " + std::to_string(nsd) + "x" + std::to_string(face_->eNoN) + ").");
    }

    Array<double> Nx_g = face_->Nx.slice(g);

    Array<double> xXi(nsd, insd);
    xXi = 0.0;

    for (int a = 0; a < face_->eNoN; a++) {
        for (int i = 0; i < insd; i++) {
            for (int j = 0; j < nsd; j++) {
                xXi(j, i) += xl(j, a) * Nx_g(i, a);
            }
        }
    }

    double Jac = 0.0;
    Vector<double> n(nsd);
    if (nsd == 3 && insd == 2) {
        n = utils::cross(xXi);
        Jac = sqrt(utils::norm(n));
    } else if (nsd == 2 && insd == 1) {
        Jac = sqrt(utils::norm(xXi.col(0)));
        n(0) = -xXi(1, 0);
        n(1) = xXi(0, 0);
    } else {
        throw std::runtime_error("[CappingSurface::compute_jacobian_and_normal] Unsupported nsd/insd combination: " +
                                std::to_string(nsd) + "/" + std::to_string(insd));
    }

    if (utils::is_zero(Jac)) {
        throw std::runtime_error("[CappingSurface::compute_jacobian_and_normal] Zero Jacobian at Gauss point " +
                                std::to_string(g));
    }

    n = n / Jac;

    if (initial_normals_.ncols() > 0 && initial_normals_.nrows() == nsd) {
        if (e < 0 || e >= initial_normals_.ncols()) {
            throw std::runtime_error("[CappingSurface::compute_jacobian_and_normal] Element index e=" +
                                    std::to_string(e) + " is out of bounds for initial_normals_ (ncols=" +
                                    std::to_string(initial_normals_.ncols()) + ").");
        }

        Vector<double> n0(nsd);
        for (int i = 0; i < nsd; i++) {
            n0(i) = initial_normals_(i, e);
        }

        double n0_norm = sqrt(utils::norm(n0));
        if (!utils::is_zero(n0_norm)) {
            n0 = n0 / n0_norm;

            double dot_product = 0.0;
            for (int i = 0; i < nsd; i++) {
                dot_product += n(i) * n0(i);
            }

            if (dot_product < 0.0) {
                n = -n;
            }
        }
    }

    return std::make_pair(Jac, n);
}

double CappingSurface::integrate_over(ComMod& com_mod, const CmMod& cm_mod, const Array<double>& s,
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
            return compute_jacobian_and_normal(xl, e, g, l_nsd, l_insd);
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

    if (cm.seq()) {
        if (!is_ready()) {
            return 0.0;
        }
        faceType* cap_face = face_.get();
        auto get_xl = [&](int e) {
            return update_element_position(com_mod, e, gnNo_to_tnNo_, cfg);
        };
        auto get_value = [&](int e, int a, int comp) -> double {
            int gnNo_idx = cap_face->IEN(a, e);
            auto it = gnNo_to_tnNo_.find(gnNo_idx);
            if (it == gnNo_to_tnNo_.end()) {
                throw std::runtime_error("[CappingSurface::integrate_over] IEN entry references gnNo not found in gnNo_to_tnNo_.");
            }
            int Ac = it->second;
            if (Ac < 0 || Ac >= com_mod.tnNo) {
                throw std::runtime_error("[CappingSurface::integrate_over] Invalid tnNo index while integrating cap.");
            }
            if (l + comp >= s.nrows() || Ac >= s.ncols()) {
                throw std::runtime_error("[CappingSurface::integrate_over] Array bounds exceeded while integrating cap.");
            }
            return s(l + comp, Ac);
        };
        double result = integrate_kernel(cap_face, get_xl, get_value);
        return result;
    }

    if (nNo_gathered_ == 0) {
        return 0.0;
    }
    const int cap_nNo = nNo_gathered_;
    const Array<double>& cap_s_use = (cfg == MechanicalConfigurationType::new_timestep) ? Yn_gathered_ : Yo_gathered_;
    double result = 0.0;
    if (cm.idcm() == cm_mod.master) {
        faceType* cap_face = face_.get();
        std::unordered_map<int, int> gnNo_to_capIdx;
        gnNo_to_capIdx.reserve(cap_nNo);
        for (int a = 0; a < cap_nNo; a++) {
            gnNo_to_capIdx[cap_face->gN(a)] = a;
        }
        auto get_xl = [&](int e) {
            return update_element_position(e, cfg, x_gathered_, Do_gathered_, Dn_gathered_, gnNo_to_capIdx);
        };
        auto get_value = [&](int e, int a, int comp) -> double {
            int gnNo_idx = cap_face->IEN(a, e);
            auto it = gnNo_to_capIdx.find(gnNo_idx);
            if (it == gnNo_to_capIdx.end()) {
                throw std::runtime_error("[CappingSurface::integrate_over] IEN entry references gnNo not found in gathered cap map.");
            }
            int cap_idx = it->second;
            if (comp >= cap_s_use.nrows() || cap_idx >= cap_s_use.ncols()) {
                throw std::runtime_error("[CappingSurface::integrate_over] Gathered array bounds exceeded while integrating cap.");
            }
            return cap_s_use(comp, cap_idx);
        };
        result = integrate_kernel(cap_face, get_xl, get_value);
    }
    return result;
}

void CappingSurface::compute_valM(ComMod& com_mod, const CmMod& cm_mod, consts::MechanicalConfigurationType cfg)
{
    using namespace consts;

    if (!is_ready()) {
        valM_.resize(0, 0);
        return;
    }
    faceType* cap_face = face_.get();
    int nsd = com_mod.nsd;
    int cap_nNo = cap_face->nNo;

    valM_.resize(nsd, cap_nNo);
    valM_ = 0.0;

    std::unordered_map<int, int> gnNo_to_cap_local;
    for (int a = 0; a < cap_nNo; a++) {
        int gnNo = cap_face->gN(a);
        gnNo_to_cap_local[gnNo] = a;
    }

    auto& cm = com_mod.cm;
    const bool use_gathered = !cm.seq() && nNo_gathered_ > 0;

    auto get_xl = [&](int e) {
        if (use_gathered) {
            return update_element_position(e, cfg, x_gathered_, Do_gathered_, Dn_gathered_, gnNo_to_cap_local);
        }
        return update_element_position(com_mod, e, gnNo_to_tnNo_, cfg);
    };
    auto compute_jac_n = [&](const Array<double>& xl, int e, int g, int l_nsd, int l_insd) {
        return compute_jacobian_and_normal(xl, e, g, l_nsd, l_insd);
    };
    auto on_gauss = [&](int e, int g, double Jac, const Vector<double>& n) {
        for (int a = 0; a < cap_face->eNoN; a++) {
            int gnNo_idx = cap_face->IEN(a, e);
            auto it = gnNo_to_cap_local.find(gnNo_idx);
            if (it == gnNo_to_cap_local.end()) {
                throw std::runtime_error("[CappingSurface::compute_valM] IEN entry (element " +
                                        std::to_string(e) + ", node " + std::to_string(a) +
                                        ") contains invalid gnNo index " + std::to_string(gnNo_idx) +
                                        " not found in cap face nodes.");
            }
            int cap_a = it->second;
            if (cap_a < 0 || cap_a >= cap_nNo) {
                throw std::runtime_error("[CappingSurface::compute_valM] Invalid cap face-local index cap_a=" +
                                        std::to_string(cap_a) + " (cap_nNo=" + std::to_string(cap_nNo) + ")");
            }
            for (int i = 0; i < nsd; i++) {
                valM_(i, cap_a) += cap_face->N(a, g) * cap_face->w(g) * Jac * n(i);
            }
        }
    };
    for_each_cap_gauss_point(cap_face, nsd, nsd - 1, get_xl, compute_jac_n, on_gauss);

    if (!cm.seq() && !use_gathered) {
        for (int i = 0; i < nsd; i++) {
            Vector<double> row = valM_.row(i);
            row = cm.reduce(cm_mod, row);
            valM_.set_row(i, row);
        }
    }
}

void CappingSurface::copy_to_linear_solver_face(ComMod& com_mod, const CmMod& cm_mod,
                                                fsi_linear_solver::FSILS_faceType& lhs_face,
                                                consts::MechanicalConfigurationType cfg)
{
    const int nsd = com_mod.nsd;
    const bool serial = com_mod.cm.seq();

    if (serial) {
        if (!is_ready()) {
            lhs_face.cap_val.resize(0, 0);
            lhs_face.cap_valM.resize(0, 0);
            lhs_face.cap_glob.resize(0);
            return;
        }
        compute_valM(com_mod, cm_mod, cfg);
        const faceType* cap_face = face_.get();
        int cap_nNo = cap_face->nNo;
        lhs_face.cap_val = valM_;
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

    int cap_nNo = 0;
    Vector<int> cap_gN_all;
    Array<double> cap_val_all;

    if (is_ready()) {
        compute_valM(com_mod, cm_mod, cfg);
        const faceType* cap_face = face_.get();
        if (cap_face->nNo > 0) {
            cap_nNo = cap_face->nNo;
            cap_gN_all.resize(cap_nNo);
            cap_val_all.resize(nsd, cap_nNo);
            for (int a = 0; a < cap_nNo; a++) {
                cap_gN_all(a) = cap_face->gN(a);
                for (int i = 0; i < nsd; i++)
                    cap_val_all(i, a) = valM_(i, a);
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
    const bool i_am_sender = is_ready();
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
