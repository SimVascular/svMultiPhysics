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
#include <unordered_map>
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
    , owns_cap_(other.owns_cap_)
    , cap_mesh_global_node_ids_(other.cap_mesh_global_node_ids_)
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
        owns_cap_ = other.owns_cap_;
        cap_mesh_global_node_ids_ = other.cap_mesh_global_node_ids_;
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
    , owns_cap_(other.owns_cap_)
    , cap_mesh_global_node_ids_(std::move(other.cap_mesh_global_node_ids_))
    , cap_(std::move(other.cap_))
{
    other.face_ = nullptr;
    other.logger_ = nullptr;
    other.has_cap_ = false;
    other.owns_cap_ = false;
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
        owns_cap_ = other.owns_cap_;
        cap_mesh_global_node_ids_ = std::move(other.cap_mesh_global_node_ids_);
        cap_ = std::move(other.cap_);

        other.face_ = nullptr;
        other.logger_ = nullptr;
        other.has_cap_ = false;
        other.owns_cap_ = false;
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

/// @brief Constructor for a coupled boundary condition
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

/// @brief Constructor for a coupled boundary condition with a cap
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

    // Distribute cap mesh global node IDs
    int cap_nn = static_cast<int>(cap_mesh_global_node_ids_.size());
    cm.bcast(cm_mod, &cap_nn);
    if (is_slave) {
        cap_mesh_global_node_ids_.resize(cap_nn);
    }
    if (cap_nn > 0) {
        cm.bcast(cm_mod, cap_mesh_global_node_ids_);
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


// =========================================================================
// Cap surface loading and integration
// =========================================================================

/// @brief Load the cap face VTP file and associate it with this boundary condition
/// @param vtp_file_path Path to the cap face VTP file
void CoupledBoundaryCondition::load_cap_face_vtp(const std::string& vtp_file_path)
{
    cap_face_vtp_file_ = vtp_file_path;
    has_cap_ = false;
    owns_cap_ = false;
    cap_.reset();
    cap_mesh_global_node_ids_.resize(0);

    if (vtp_file_path.empty()) {
        return;
    }

    cap_.emplace();
    try {
        cap_->load_from_vtp(vtp_file_path, *face_, face_name_);
    } catch (...) {
        cap_.reset();
        has_cap_ = false;
        owns_cap_ = false;
        cap_mesh_global_node_ids_.resize(0);
        throw;
    }
    has_cap_ = true;
    owns_cap_ = true;
    const faceType* cf = cap_->face();
    if (cf != nullptr && cf->nNo > 0) {
        cap_mesh_global_node_ids_.resize(cf->nNo);
        for (int a = 0; a < cf->nNo; a++) {
            cap_mesh_global_node_ids_(a) = cf->gN(a);
        }
    }
}

namespace {

void gather_global_mesh_state_serial(ComMod& com_mod, const Array<double>& Yo, const Array<double>& Yn, int s_comps,
                                     int nsd, int gtnNo, int tnNo, CapGlobalMeshState& out)
{
    if (tnNo <= 0) {
        return;
    }
    out.gtnNo = gtnNo;
    out.x.resize(nsd, gtnNo);
    out.Do.resize(nsd, gtnNo);
    out.Dn.resize(nsd, gtnNo);
    out.x = 0.0;
    out.Do = 0.0;
    out.Dn = 0.0;
    if (s_comps > 0) {
        out.Yo.resize(s_comps, gtnNo);
        out.Yn.resize(s_comps, gtnNo);
        out.Yo = 0.0;
        out.Yn = 0.0;
    }
    for (int Ac = 0; Ac < tnNo; Ac++) {
        int g = com_mod.ltg(Ac);
        if (g < 0 || g >= gtnNo) {
            continue;
        }
        for (int i = 0; i < nsd; i++) {
            out.x(i, g) = com_mod.x(i, Ac);
            out.Do(i, g) = com_mod.Do(i, Ac);
            out.Dn(i, g) = com_mod.Dn(i, Ac);
        }
        for (int i = 0; i < s_comps; i++) {
            out.Yo(i, g) = Yo(i, Ac);
            out.Yn(i, g) = Yn(i, Ac);
        }
    }
}

void gather_global_mesh_state_parallel(ComMod& com_mod, const CmMod& cm_mod, cmType& cm, const Array<double>& Yo,
                                       const Array<double>& Yn, int s_comps, int nsd, int gtnNo, int tnNo, int root,
                                       int nProcs, CapGlobalMeshState& out)
{
    const int per_node = 1 + 3 * nsd + 2 * s_comps;
    const int nPack = tnNo > 0 ? tnNo : 0;
    Vector<double> send_buf(static_cast<size_t>(nPack * per_node));
    int idx = 0;
    for (int Ac = 0; Ac < tnNo; Ac++) {
        int g = com_mod.ltg(Ac);
        send_buf(idx++) = static_cast<double>(g);
        for (int i = 0; i < nsd; i++) {
            send_buf(idx++) = com_mod.x(i, Ac);
        }
        for (int i = 0; i < nsd; i++) {
            send_buf(idx++) = com_mod.Do(i, Ac);
        }
        for (int i = 0; i < nsd; i++) {
            send_buf(idx++) = com_mod.Dn(i, Ac);
        }
        for (int i = 0; i < s_comps; i++) {
            send_buf(idx++) = Yo(i, Ac);
        }
        for (int i = 0; i < s_comps; i++) {
            send_buf(idx++) = Yn(i, Ac);
        }
    }

    Vector<int> send_count_vec(1);
    send_count_vec(0) = static_cast<int>(send_buf.size());
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

    if (cm.idcm() != root) {
        return;
    }

    out.gtnNo = gtnNo;
    out.x.resize(nsd, gtnNo);
    out.Do.resize(nsd, gtnNo);
    out.Dn.resize(nsd, gtnNo);
    out.x = 0.0;
    out.Do = 0.0;
    out.Dn = 0.0;
    if (s_comps > 0) {
        out.Yo.resize(s_comps, gtnNo);
        out.Yn.resize(s_comps, gtnNo);
        out.Yo = 0.0;
        out.Yn = 0.0;
    }

    int pos = 0;
    while (pos + per_node <= static_cast<int>(recv_buf.size())) {
        int g = static_cast<int>(recv_buf(pos++));
        if (g < 0 || g >= gtnNo) {
            pos += 3 * nsd + 2 * s_comps;
            continue;
        }
        for (int i = 0; i < nsd; i++) {
            out.x(i, g) = recv_buf(pos++);
        }
        for (int i = 0; i < nsd; i++) {
            out.Do(i, g) = recv_buf(pos++);
        }
        for (int i = 0; i < nsd; i++) {
            out.Dn(i, g) = recv_buf(pos++);
        }
        for (int i = 0; i < s_comps; i++) {
            out.Yo(i, g) = recv_buf(pos++);
        }
        for (int i = 0; i < s_comps; i++) {
            out.Yn(i, g) = recv_buf(pos++);
        }
    }
}

} // namespace

/// Assemble \ref CapGlobalMeshState on the MPI root from local mesh arrays (\c ltg maps local node → global column).
void CoupledBoundaryCondition::gather_global_mesh_state(ComMod& com_mod, const CmMod& cm_mod, const Array<double>& Yo,
                                                        const Array<double>& Yn, int s_comps, CapGlobalMeshState& out)
{
    out.clear();
    auto& cm = com_mod.cm;
    const int nsd = com_mod.nsd;
    const int gtnNo = com_mod.gtnNo;
    const int tnNo = com_mod.tnNo;
    const int root = cm_mod.master;
    const int nProcs = cm.np();

    // Serial: no MPI, just gather the global mesh state locally.
    if (cm.seq()) {
        gather_global_mesh_state_serial(com_mod, Yo, Yn, s_comps, nsd, gtnNo, tnNo, out);
        return;
    }

    // Parallel: gather the global mesh state on the MPI root.
    gather_global_mesh_state_parallel(com_mod, cm_mod, cm, Yo, Yn, s_comps, nsd, gtnNo, tnNo, root, nProcs, out);
}

/// @brief Initialize the cap quadrature.
/// @param com_mod The com_mod object.
/// @param cm_mod The cm_mod object.
void CoupledBoundaryCondition::initialize_cap(ComMod& com_mod)
{
    if (!has_cap_ || !owns_cap_ || !cap_) {
        return;
    }
    if (com_mod.nsd != 3) {
        throw std::runtime_error("[CoupledBoundaryCondition::initialize_cap] Cap surface requires nsd=3 (TRI3 surface in 3D).");
    }
    if (cap_->face()) {
        cap_->init_cap_face_quadrature(com_mod);
        cap_->valM_.resize(3, cap_->face()->nNo);
        cap_->valM_ = 0.0;
    }
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

    if (!has_cap_) {
        return {Qo_cap, Qn_cap};
    }

    CapGlobalMeshState st;
    gather_global_mesh_state(com_mod, cm_mod, com_mod.Yo, com_mod.Yn, nsd, st);

    if ((serial_run || i_am_master) && owns_cap_ && cap_ && cap_->face()) {
        Qo_cap = cap_->integrate_velocity_flux(st, false, 0, cfg_o);
        Qn_cap = cap_->integrate_velocity_flux(st, true, 0, cfg_n);
    }

    if (!serial_run) {
        cm.bcast(cm_mod, &Qo_cap);
        cm.bcast(cm_mod, &Qn_cap);
    }
    return {Qo_cap, Qn_cap};
}

// =========================================================================
// CappingSurface (definitions; class in CoupledBoundaryCondition.h)
// =========================================================================

namespace {

/// Load cap VTP and validate point/element counts and presence of GlobalNodeID (used from \c CappingSurface::load_from_vtp).
VtkVtpData load_cap_vtp(const std::string& vtp_file_path)
{
    VtkVtpData vtp_data;
    try {
        vtp_data = VtkVtpData(vtp_file_path, true);
    } catch (const std::exception& e) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Failed to construct VtkVtpData from file '" +
                                vtp_file_path + "': " + e.what());
    }

    int nNo = 0;
    try {
        nNo = vtp_data.num_points();
    } catch (const std::exception& e) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Failed to get number of points from VTP file '" +
                                vtp_file_path + "': " + e.what());
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
    }
    if (!has_global_node_id) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Cap VTP file '" + vtp_file_path +
                                "' does not contain 'GlobalNodeID' point data.");
    }

    return vtp_data;
}

/// Optionally load cell "Normals" from cap VTP into \p initial_normals (used from \c CappingSurface::load_from_vtp).
void load_cap_vtp_normals(VtkVtpData& vtp_data, const std::string& vtp_file_path, int num_elems,
                                    Array<double>& initial_normals)
{
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

            initial_normals.resize(num_comp, num_elems);
            vtp_data.copy_cell_data("Normals", initial_normals);

            if (initial_normals.nrows() != num_comp || initial_normals.ncols() != num_elems) {
                throw std::runtime_error("[CappingSurface::load_from_vtp] Failed to copy Normals data. "
                                        "Expected size: " + std::to_string(num_comp) + "x" + std::to_string(num_elems) +
                                        ", Actual size: " + std::to_string(initial_normals.nrows()) + "x" +
                                        std::to_string(initial_normals.ncols()));
            }
        } else {
            initial_normals.resize(0, 0);
        }
    } catch (const std::exception& e) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Failed to load Normals from VTP file '" +
                                vtp_file_path + "': " + e.what());
    }
}

/// Throws if cap and coupled boundary faces share no global mesh node (used from \c CappingSurface::load_from_vtp).
void check_cap_shares_nodes_with_coupled_face(const faceType& coupled_face, const faceType& cap_face,
                                                       const std::string& vtp_file_path,
                                                       const std::string& coupled_face_name)
{
    std::unordered_set<int> coupled_gn;
    coupled_gn.reserve(static_cast<size_t>(coupled_face.nNo));
    for (int a = 0; a < coupled_face.nNo; a++) {
        coupled_gn.insert(coupled_face.gN(a));
    }
    for (int a = 0; a < cap_face.nNo; a++) {
        if (coupled_gn.find(cap_face.gN(a)) != coupled_gn.end()) {
            return;
        }
    }
    throw std::runtime_error(
        "[CappingSurface::load_from_vtp] Cap VTP file '" + vtp_file_path +
        "' has no GlobalNodeID entries in common with coupled face '" + coupled_face_name +
        "'. The cap must share at least one mesh node with that face.");
}

template <typename GetXlFn, typename ComputeJacNFn, typename OnGaussFn>
void for_each_cap_gauss_point(const faceType* cap_face,
                              GetXlFn&& get_xl,
                              ComputeJacNFn&& compute_jac_n,
                              OnGaussFn&& on_gauss)
{
    for (int e = 0; e < cap_face->nEl; e++) {
        Array<double> xl = get_xl(e);
        for (int g = 0; g < cap_face->nG; g++) {
            auto [Jac, n] = compute_jac_n(xl, e, g);
            on_gauss(e, g, Jac, n);
        }
    }
}

} // namespace

CappingSurface::CappingSurface(const CappingSurface& other)
    : global_node_ids_(other.global_node_ids_)
    , valM_(other.valM_)
    , initial_normals_(other.initial_normals_)
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
        valM_ = other.valM_;
        initial_normals_ = other.initial_normals_;
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
    
    // Check if the VTP file exists.
    std::ifstream file_check(vtp_file_path);
    if (!file_check.good()) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Cannot open cap VTP file '" + vtp_file_path +
                                "' for reading.");
    }
    file_check.close();

    // Load the VTP file and validate the header.
    VtkVtpData vtp_data = load_cap_vtp(vtp_file_path);
    const int nNo = vtp_data.num_points();
    const int num_elems = vtp_data.num_elems();

    // Check if the VTP file contains the GlobalNodeID point data and 
    // copy it to the global_node_ids_ array.
    try {
        global_node_ids_.resize(nNo);
        vtp_data.copy_point_data("GlobalNodeID", global_node_ids_);
    } catch (const std::exception& e) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Failed to copy GlobalNodeID from VTP file '" +
                                vtp_file_path + "': " + e.what());
    }

    // Get the connectivity from the VTP file and validate the number of nodes per element.
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
    }
    if (vtk_cell_type != VTK_TRIANGLE) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Unsupported cap cell type " +
                                std::to_string(vtk_cell_type) + ". Only VTK_TRIANGLE (TRI3) is supported.");
    }
    if (eNoN != 3) {
        throw std::runtime_error("[CappingSurface::load_from_vtp] Invalid nodes-per-element for triangle cap: " +
                                std::to_string(eNoN) + " (expected 3).");
    }

    // Create the face object (cap path assumes TRI3 everywhere).
    face_ = std::make_unique<faceType>();
    face_->name = coupled_face_name + "_cap";
    face_->iM = coupled_face.iM;
    face_->nNo = nNo;
    face_->nEl = num_elems;
    face_->gnEl = num_elems;
    face_->eNoN = 3;
    face_->eType = consts::ElementType::TRI3;

    // Copy the global node IDs to the face object.
    face_->gN.resize(nNo);
    for (int a = 0; a < nNo; a++) {
        face_->gN(a) = global_node_ids_(a) - 1;
    }

    // Copy the connectivity to the face object.
    face_->IEN.resize(face_->eNoN, num_elems);

    for (int e = 0; e < num_elems; e++) {
        for (int a = 0; a < face_->eNoN; a++) {
            int local_node_idx = conn(a, e);
            face_->IEN(a, e) = face_->gN(local_node_idx);
        }
    }

    // Load the normals from the VTP file.
    load_cap_vtp_normals(vtp_data, vtp_file_path, num_elems, initial_normals_);

    // Check if the cap face shares any nodes with the coupled face.
    check_cap_shares_nodes_with_coupled_face(coupled_face, *face_, vtp_file_path, coupled_face_name);
}

/// @brief Initializes the cap face quadrature (assumes triangular elements).
/// @param com_mod The com_mod object.
void CappingSurface::init_cap_face_quadrature(const ComMod& com_mod)
{
    using namespace consts;
    int nsd = com_mod.nsd;

    try {
        if (nsd != cap_nsd_) {
            throw std::runtime_error("[CappingSurface::init_cap_face_quadrature] Cap surface requires nsd=3.");
        }
        face_->nG = 1;

        face_->w.resize(face_->nG);
        face_->xi.resize(cap_insd_, face_->nG);

        face_->w(0) = 0.5;
        face_->xi(0, 0) = 1.0 / 3.0;
        face_->xi(1, 0) = 1.0 / 3.0;

        face_->N.resize(face_->eNoN, face_->nG);
        face_->Nx.resize(cap_insd_, face_->eNoN, face_->nG);
        for (int g = 0; g < face_->nG; g++) {
            nn::get_gnn(cap_insd_, face_->eType, face_->eNoN, g, face_->xi, face_->N, face_->Nx);
        }
    } catch (const std::exception& e) {
        throw std::runtime_error("[CappingSurface::init_cap_face_quadrature] Failed to initialize cap face shape functions: " +
                                std::string(e.what()));
    }
}

/// @brief Updates the element position in global coordinates.
/// @param e The element index.
/// @param cfg The mechanical configuration type.
/// @param mesh_x The mesh x coordinates.
/// @param mesh_Do The mesh Do coordinates.
/// @param mesh_Dn The mesh Dn coordinates.
/// @param gtnNo_cols The number of global nodes.
/// @return The element position in global coordinates.
Array<double> CappingSurface::update_element_position_global(int e, consts::MechanicalConfigurationType cfg,
                                                             const Array<double>& mesh_x, const Array<double>& mesh_Do,
                                                             const Array<double>& mesh_Dn, int gtnNo_cols) const
{
    using namespace consts;

    if (mesh_x.nrows() < cap_nsd_ || mesh_Do.nrows() < cap_nsd_ || mesh_Dn.nrows() < cap_nsd_) {
        throw std::runtime_error("[CappingSurface::update_element_position_global] Mesh arrays must have at least 3 rows.");
    }
    Array<double> xl(cap_nsd_, face_->eNoN);

    for (int a = 0; a < face_->eNoN; a++) {
        int g = face_->IEN(a, e);
        for (int i = 0; i < cap_nsd_; i++) {
            xl(i, a) = mesh_x(i, g);
        }
        if (cfg == MechanicalConfigurationType::old_timestep) {
            for (int i = 0; i < cap_nsd_; i++) {
                xl(i, a) += mesh_Do(i, g);
            }
        } else if (cfg == MechanicalConfigurationType::new_timestep) {
            for (int i = 0; i < cap_nsd_; i++) {
                xl(i, a) += mesh_Dn(i, g);
            }
        }
    }

    return xl;
}

std::pair<double, Vector<double>> CappingSurface::compute_jacobian_and_normal(const Array<double>& xl, int e, int g)
{
    if (xl.nrows() != cap_nsd_ || xl.ncols() != face_->eNoN) {
        throw std::runtime_error("[CappingSurface::compute_jacobian_and_normal] xl has wrong dimensions: " +
                                std::to_string(xl.nrows()) + "x" + std::to_string(xl.ncols()) +
                                " (expected " + std::to_string(cap_nsd_) + "x" + std::to_string(face_->eNoN) + ").");
    }

    // Get the shape function derivatives for the Gauss point.
    Array<double> Nx_g = face_->Nx.rslice(g);
    Array<double> xXi(cap_nsd_, cap_insd_);
    xXi = 0.0;

    // Compute the Jacobian matrix of the element.
    for (int a = 0; a < face_->eNoN; a++) {
        for (int i = 0; i < cap_insd_; i++) {
            for (int j = 0; j < cap_nsd_; j++) {
                xXi(j, i) += xl(j, a) * Nx_g(i, a);
            }
        }
    }

    // Compute the Jacobian and normal vector.
    double Jac = 0.0;
    Vector<double> n(cap_nsd_);
    n = utils::cross(xXi);
    Jac = sqrt(utils::norm(n));

    if (utils::is_zero(Jac)) {
        throw std::runtime_error("[CappingSurface::compute_jacobian_and_normal] Zero Jacobian at Gauss point " +
                                std::to_string(g));
    }

    n = n / Jac;

    // Check if the initial normals are provided and if they are valid.
    if (initial_normals_.ncols() > 0 && initial_normals_.nrows() == cap_nsd_) {
        if (e < 0 || e >= initial_normals_.ncols()) {
            throw std::runtime_error("[CappingSurface::compute_jacobian_and_normal] Element index e=" +
                                    std::to_string(e) + " is out of bounds for initial_normals_ (ncols=" +
                                    std::to_string(initial_normals_.ncols()) + ").");
        }

        Vector<double> n0(cap_nsd_);
        for (int i = 0; i < cap_nsd_; i++) {
            n0(i) = initial_normals_(i, e);
        }

        double n0_norm = sqrt(utils::norm(n0));
        if (!utils::is_zero(n0_norm)) {
            n0 = n0 / n0_norm;

            double dot_product = 0.0;
            for (int i = 0; i < cap_nsd_; i++) {
                dot_product += n(i) * n0(i);
            }

            if (dot_product < 0.0) {
                n = -n;
            }
        }
    }

    return std::make_pair(Jac, n);
}

double CappingSurface::integrate_velocity_flux(const CapGlobalMeshState& st, bool use_Yn_velocity, int l_vel,
                                                 consts::MechanicalConfigurationType cfg)
{
    using namespace consts;
    if (!face_ || st.gtnNo <= 0) {
        return 0.0;
    }
    const int s_comps = cap_nsd_;
    const Array<double>& cap_vel = use_Yn_velocity ? st.Yn : st.Yo;
    if (cap_vel.nrows() < l_vel + s_comps || cap_vel.ncols() < st.gtnNo ||
        st.x.nrows() < cap_nsd_ || st.Do.nrows() < cap_nsd_ || st.Dn.nrows() < cap_nsd_ || st.x.ncols() < st.gtnNo) {
        return 0.0;
    }

    faceType* cap_face = face_.get();

    auto integrate_kernel = [&](const faceType* cf, auto&& get_xl, auto&& get_value) -> double {
        double result = 0.0;
        auto compute_jac_n = [&](const Array<double>& xl, int e, int g) {
            return compute_jacobian_and_normal(xl, e, g);
        };
        auto on_gauss = [&](int e, int g, double Jac, const Vector<double>& n) {
            double sHat = 0.0;
            for (int a = 0; a < cf->eNoN; a++) {
                const double Na = cf->N(a, g);
                for (int i = 0; i < s_comps; i++) {
                    sHat += Na * get_value(e, a, i) * n(i);
                }
            }
            result += cf->w(g) * Jac * sHat;
        };
        for_each_cap_gauss_point(cf, get_xl, compute_jac_n, on_gauss);
        return result;
    };

    auto get_xl = [&](int e) { return update_element_position_global(e, cfg, st.x, st.Do, st.Dn, st.gtnNo); };
    auto get_value = [&](int e, int a, int comp) -> double {
        int g = cap_face->IEN(a, e);
        if (g < 0 || g >= st.gtnNo) {
            throw std::runtime_error("[CappingSurface::integrate_velocity_flux] IEN global id " + std::to_string(g) +
                                    " out of range.");
        }
        if (l_vel + comp >= cap_vel.nrows()) {
            throw std::runtime_error("[CappingSurface::integrate_velocity_flux] Velocity array row out of range.");
        }
        return cap_vel(l_vel + comp, g);
    };
    return integrate_kernel(cap_face, get_xl, get_value);
}


/// @brief Computes the cap contribution to the linear solver face.
/// @param com_mod The com_mod object.
/// @param cfg The mechanical configuration type.
/// @param st The cap global mesh state.
void CappingSurface::compute_valM(ComMod& com_mod, consts::MechanicalConfigurationType cfg,
                                  const CapGlobalMeshState& st)
{
    using namespace consts;
    (void)com_mod;

    int cap_nNo = face_->nNo;

    valM_ = 0.0;

    // Map global node IDs to cap face-local indices
    std::unordered_map<int, int> gnNo_to_cap_local;
    for (int a = 0; a < cap_nNo; a++) {
        int gnNo = face_->gN(a);
        gnNo_to_cap_local[gnNo] = a;
    }

    // Get the element position in global coordinates
    auto get_xl = [&](int e) {
        return update_element_position_global(e, cfg, st.x, st.Do, st.Dn, st.gtnNo);
    };
    auto compute_jac_n = [&](const Array<double>& xl, int e, int g) {
        return compute_jacobian_and_normal(xl, e, g);
    };
    auto on_gauss = [&](int e, int g, double Jac, const Vector<double>& n) {
        for (int a = 0; a < face_->eNoN; a++) {
            int gnNo_idx = face_->IEN(a, e);
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
            for (int i = 0; i < cap_nsd_; i++) {
                valM_(i, cap_a) += face_->N(a, g) * face_->w(g) * Jac * n(i);
            }
        }
    };
    for_each_cap_gauss_point(face_.get(), get_xl, compute_jac_n, on_gauss);
}

namespace {


/// @brief Broadcasts the cap contribution to the linear solver face to all ranks.
/// @param com_mod The com_mod object.
/// @param cm_mod The cm_mod object.
/// @param lhs_face The linear solver face.
/// @param cap_nNo The number of cap nodes.
/// @param cap_gN_all The global node IDs of the cap nodes.
/// @param cap_val_all The values of the cap nodes.
void bcast_cap_lhs_contribution(ComMod& com_mod, const CmMod& cm_mod,
                                fsi_linear_solver::FSILS_faceType& lhs_face, int& cap_nNo,
                                Vector<int>& cap_gN_all, Array<double>& cap_val_all)
{
    auto& cm = com_mod.cm;
    const bool i_am_sender = cm.mas(cm_mod);
    int nsd = com_mod.nsd;

    // Resize the cap contribution to the linear solver face (in all ranks)
    if (!i_am_sender) {
        cap_gN_all.resize(cap_nNo);
        cap_val_all.resize(nsd, cap_nNo);
    }
    cm.bcast(cm_mod, cap_gN_all);
    cm.bcast(cm_mod, cap_val_all);

    // Count the number of owned cap nodes
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

    // Resize the cap contribution to the linear solver face (in all ranks)
    lhs_face.cap_glob.resize(n_owned);
    lhs_face.cap_val.resize(nsd, n_owned);
    lhs_face.cap_valM.resize(nsd, n_owned);
    lhs_face.cap_valM = 0.0;

    // Fill the cap contribution to the linear solver face (in all ranks)
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
            for (int i = 0; i < nsd; i++) {
                lhs_face.cap_val(i, idx) = cap_val_all(i, a);
            }
            idx++;
        }
    }
}

} // namespace

/// @brief Calculates the cap contribution to the linear solver face and broadcasts it to all ranks.
void CoupledBoundaryCondition::copy_cap_surface_to_linear_solver_face(ComMod& com_mod, const CmMod& cm_mod,
                                                                      fsi_linear_solver::FSILS_faceType& lhs_face,
                                                                      consts::MechanicalConfigurationType cfg)
{
    const int nsd = com_mod.nsd;

    CapGlobalMeshState gstate;
    gather_global_mesh_state(com_mod, cm_mod, com_mod.Yo, com_mod.Yn, 0, gstate);

    int cap_nNo = 0;
    Vector<int> cap_gN_all;
    Array<double> cap_val_all;

    // Calculate the cap contribution to the linear solver face (in master rank)
    if (owns_cap_ && cap_) {
        const faceType* cap_face = cap_->face();
        cap_->compute_valM(com_mod, cfg, gstate);
        cap_nNo = cap_face->nNo;
        cap_gN_all.resize(cap_nNo);
        cap_val_all.resize(nsd, cap_nNo);
        for (int a = 0; a < cap_nNo; a++) {
            cap_gN_all(a) = cap_face->gN(a);
            for (int i = 0; i < nsd; i++) {
                cap_val_all(i, a) = cap_->valM()(i, a);
            }
        }
    }

    // Broadcast the cap contribution to all ranks
    bcast_cap_lhs_contribution(com_mod, cm_mod, lhs_face, cap_nNo, cap_gN_all, cap_val_all);
}
