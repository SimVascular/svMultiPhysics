// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "ComMod.h"
#include "CoupledBoundaryCondition.h"
#include "all_fun.h"
#include "utils.h"
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
    , cap_n_no_(other.cap_n_no_)
    , cap_mesh_global_node_ids_(other.cap_mesh_global_node_ids_)
    , cap_g_to_cap_col_(other.cap_g_to_cap_col_)
    , cap_(other.cap_)
    , cap_global_mesh_state_(other.cap_global_mesh_state_)
    , cm_mod_(other.cm_mod_)
{
}

CoupledBoundaryCondition& CoupledBoundaryCondition::operator=(const CoupledBoundaryCondition& other)
{
    if (this != &other) {
        face_ = other.face_;
        cap_face_vtp_file_ = other.cap_face_vtp_file_;
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
        cap_n_no_ = other.cap_n_no_;
        cap_mesh_global_node_ids_ = other.cap_mesh_global_node_ids_;
        cap_g_to_cap_col_ = other.cap_g_to_cap_col_;
        cap_ = other.cap_;
        cap_global_mesh_state_ = other.cap_global_mesh_state_;
        cm_mod_ = other.cm_mod_;
    }
    return *this;
}

CoupledBoundaryCondition::CoupledBoundaryCondition(CoupledBoundaryCondition&& other) noexcept
    : face_(other.face_)
    , cap_face_vtp_file_(std::move(other.cap_face_vtp_file_))
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
    , cap_n_no_(other.cap_n_no_)
    , cap_mesh_global_node_ids_(std::move(other.cap_mesh_global_node_ids_))
    , cap_g_to_cap_col_(std::move(other.cap_g_to_cap_col_))
    , cap_(std::move(other.cap_))
    , cap_global_mesh_state_(std::move(other.cap_global_mesh_state_))
    , cm_mod_(other.cm_mod_)
{
    other.face_ = nullptr;
    other.has_cap_ = false;
    other.owns_cap_ = false;
    other.cap_n_no_ = 0;
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
    other.cap_global_mesh_state_.clear();
}

CoupledBoundaryCondition& CoupledBoundaryCondition::operator=(CoupledBoundaryCondition&& other) noexcept
{
    if (this != &other) {
        face_ = other.face_;
        cap_face_vtp_file_ = std::move(other.cap_face_vtp_file_);
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
        cap_n_no_ = other.cap_n_no_;
        cap_mesh_global_node_ids_ = std::move(other.cap_mesh_global_node_ids_);
        cap_g_to_cap_col_ = std::move(other.cap_g_to_cap_col_);
        cap_ = std::move(other.cap_);
        cap_global_mesh_state_ = std::move(other.cap_global_mesh_state_);
        cm_mod_ = other.cm_mod_;

        other.face_ = nullptr;
        other.has_cap_ = false;
        other.owns_cap_ = false;
        other.cap_n_no_ = 0;
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
        other.cap_global_mesh_state_.clear();
    }
    return *this;
}

/// @brief Constructor for a coupled boundary condition
CoupledBoundaryCondition::CoupledBoundaryCondition(consts::BoundaryConditionType bc_type, const faceType& face, const std::string& face_name,
                                               const std::string& block_name, consts::EquationType phys, bool follower_pressure_load)
    : face_(&face)
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
                                               consts::EquationType phys, bool follower_pressure_load)
    : cap_face_vtp_file_(cap_face_vtp_file)
    , face_(&face)
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

const std::string& CoupledBoundaryCondition::get_block_name() const
{
    return block_name_;
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

/// @brief Set the flowrate mechanical configurations.
/// @param phys The equation type.
/// @param follower_pressure_load The follower pressure load flag.
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
void CoupledBoundaryCondition::compute_flowrates(ComMod& com_mod, const CmMod& cm_mod, const SolutionStates& solutions)
{
    int nsd = com_mod.nsd;
    const auto& Yo = solutions.old.get_velocity();
    const auto& Yn = solutions.current.get_velocity();
    
    Qo_ = all_fun::integ(com_mod, cm_mod, *face_, Yo, 0, solutions,
                         std::optional<int>(nsd - 1), false, flowrate_cfg_o_);
    Qn_ = all_fun::integ(com_mod, cm_mod, *face_, Yn, 0, solutions,
                         std::optional<int>(nsd - 1), false, flowrate_cfg_n_);
    
    if (has_cap_) {
        const auto [Qo_cap, Qn_cap] =
            calculate_cap_contribution(com_mod, cm_mod, solutions, flowrate_cfg_o_, flowrate_cfg_n_);
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
void CoupledBoundaryCondition::compute_pressures(ComMod& com_mod, const CmMod& cm_mod, const SolutionStates& solutions)
{
    using namespace consts;

    int nsd = com_mod.nsd;
    double area = face_->area;
    const auto& Yo = solutions.old.get_velocity();
    const auto& Yn = solutions.current.get_velocity();
    
    Po_ = all_fun::integ(com_mod, cm_mod, *face_, Yo, nsd, solutions,
                         std::nullopt, false, flowrate_cfg_o_) / area;
    Pn_ = all_fun::integ(com_mod, cm_mod, *face_, Yn, nsd, solutions,
                         std::nullopt, false, flowrate_cfg_n_) / area;
    
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

/// @brief Rebuild the map from global node ID to cap column index.
void CoupledBoundaryCondition::rebuild_cap_global_to_col_map()
{
    cap_g_to_cap_col_.clear();
    cap_g_to_cap_col_.reserve(static_cast<size_t>(cap_mesh_global_node_ids_.size()));
    for (int a = 0; a < static_cast<int>(cap_mesh_global_node_ids_.size()); a++) {
        cap_g_to_cap_col_[cap_mesh_global_node_ids_(a)] = a;
    }
}

/// @brief Distribute the boundary condition metadata to all processes.
/// @param com_mod The com_mod object.
/// @param cm_mod The cm_mod object.
/// @param cm The cmType object.
/// @param face The faceType object.
void CoupledBoundaryCondition::distribute(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, const faceType& face)
{
    #define n_debug_coupled_distribute

    // In the constructor, the face pointer is set to the global face, which is then partitioned and distributed among all processes.
    // Here, we update the face pointer to the local face.
    face_ = &face;
    cm_mod_ = cm_mod;

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
    
    // Distribute follower pressure load flag and physical configuration
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

    // Distribute cap surface global node IDs (same count on every rank; master filled at load).
    int cap_nn = static_cast<int>(cap_mesh_global_node_ids_.size());
    cm.bcast(cm_mod, &cap_nn);
    if (is_slave) {
        cap_mesh_global_node_ids_.resize(cap_nn);
    }
    if (cap_nn > 0) {
        cm.bcast(cm_mod, cap_mesh_global_node_ids_);
    }
    cap_n_no_ = cap_nn;
    rebuild_cap_global_to_col_map();
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
    cap_n_no_ = 0;
    cap_mesh_global_node_ids_.resize(0);
    cap_g_to_cap_col_.clear();

    if (vtp_file_path.empty()) {
        return;
    }
    
    // Load the cap face from the VTP file.
    cap_.emplace();
    try {
        cap_->load_from_vtp(vtp_file_path, *face_, face_name_);
    } catch (...) {
        cap_.reset();
        has_cap_ = false;
        owns_cap_ = false;
        cap_n_no_ = 0;
        cap_mesh_global_node_ids_.resize(0);
        cap_g_to_cap_col_.clear();
        cap_global_mesh_state_.clear();
        throw;
    }

    // Set the cap face flag.
    has_cap_ = true;
    owns_cap_ = true;

    // Set the cap mesh global node IDs.
    const faceType* cf = cap_->face();
    if (cf != nullptr && cf->nNo > 0) {
        cap_mesh_global_node_ids_.resize(cf->nNo);
        for (int a = 0; a < cf->nNo; a++) {
            cap_mesh_global_node_ids_(a) = cf->gN(a);
        }
        cap_n_no_ = cf->nNo;
    } else {
        cap_mesh_global_node_ids_.resize(0);
        cap_n_no_ = 0;
    }

    // Rebuild the map from global node ID to cap column index.
    rebuild_cap_global_to_col_map();
}


/// @brief Initialize the cap quadrature.
/// @param com_mod The com_mod object.
void CoupledBoundaryCondition::initialize_cap(ComMod& com_mod)
{
    if (!has_cap_ || !owns_cap_ || !cap_) {
        return;
    }
    if (com_mod.nsd != 3) {
        throw std::runtime_error("[CoupledBoundaryCondition::initialize_cap] Cap surface requires nsd=3 (TRI3 surface in 3D).");
    }

    // Initialize the cap face quadrature.
    cap_->init_cap_face_quadrature(com_mod);

    // Initialize cap contribution storage.
    cap_->initialize_valM();

    // Reset the cap global mesh state buffers.
    const int nsd = com_mod.nsd;
    const int n_cap = cap_n_no_;
    cap_global_mesh_state_.n_cap = n_cap;
    cap_global_mesh_state_.x.resize(nsd, n_cap);
    cap_global_mesh_state_.Do.resize(nsd, n_cap);
    cap_global_mesh_state_.Dn.resize(nsd, n_cap);
    cap_global_mesh_state_.Yo.resize(nsd, n_cap);
    cap_global_mesh_state_.Yn.resize(nsd, n_cap);
    cap_global_mesh_state_.x = 0.0;
    cap_global_mesh_state_.Do = 0.0;
    cap_global_mesh_state_.Dn = 0.0;
    cap_global_mesh_state_.Yo = 0.0;
    cap_global_mesh_state_.Yn = 0.0;
}

namespace {

/// @brief Gathers cap-node mesh state on the serial rank (columns 0..n_cap-1 match \p cap_gn order).
void gather_global_mesh_state_serial(ComMod& com_mod, const SolutionStates& solutions, bool gather_Y, int nsd, int tnNo,
                                     const Vector<int>& cap_gn, CapGlobalMeshState& out)
{
    const auto& Do = solutions.old.get_displacement();
    const auto& Dn = solutions.current.get_displacement();
    const auto& Yo = solutions.old.get_velocity();
    const auto& Yn = solutions.current.get_velocity();

    const int n_cap = static_cast<int>(cap_gn.size());

    // Reset the cap global mesh state buffers.
    out.x = 0.0;
    out.Do = 0.0;
    out.Dn = 0.0;
    if (gather_Y) {
        out.Yo = 0.0;
        out.Yn = 0.0;
    }
    
    // Gather the cap global mesh state.
    for (int a = 0; a < n_cap; a++) {
        const int gn = cap_gn(a);
        for (int Ac = 0; Ac < tnNo; Ac++) {
            if (com_mod.ltg(Ac) != gn) {
                continue;
            }
            for (int i = 0; i < nsd; i++) {
                out.x(i, a) = com_mod.x(i, Ac);
                out.Do(i, a) = Do(i, Ac);
                out.Dn(i, a) = Dn(i, Ac);
            }
            if (gather_Y) {
                for (int i = 0; i < nsd; i++) {
                    out.Yo(i, a) = Yo(i, Ac);
                    out.Yn(i, a) = Yn(i, Ac);
                }
            }
            break;
        }
    }
}


/// @brief Gathers cap-node mesh state on the MPI root (columns 0..n_cap-1 match \p cap_gn order).
void gather_global_mesh_state_parallel(ComMod& com_mod, const CmMod& cm_mod, cmType& cm, const SolutionStates& solutions,
                                        bool gather_Y, int nsd, int tnNo, int root, int nProcs, const Vector<int>& cap_gn,
                                        const std::unordered_map<int, int>& g_to_cap_col,
                                        CapGlobalMeshState& out)
{
    const auto& Do = solutions.old.get_displacement();
    const auto& Dn = solutions.current.get_displacement();
    const auto& Yo = solutions.old.get_velocity();
    const auto& Yn = solutions.current.get_velocity();

    // Count the number of cap nodes.
    const int n_cap = static_cast<int>(cap_gn.size());
    std::unordered_set<int> cap_gn_set;
    cap_gn_set.reserve(static_cast<size_t>(n_cap));
    for (int a = 0; a < n_cap; a++) {
        cap_gn_set.insert(cap_gn(a));
    }

    // Count the number of values to pack for each node.
    const int per_node = 1 + 3 * nsd + (gather_Y ? 2 * nsd : 0);
    int n_pack = 0;
    for (int Ac = 0; Ac < tnNo; Ac++) {
        const int g = com_mod.ltg(Ac);
        if (cap_gn_set.find(g) != cap_gn_set.end()) {
            n_pack++;
        }
    }

    // Pack the cap global mesh state.
    Vector<double> send_buf(static_cast<size_t>(n_pack * per_node));
    int idx = 0;
    for (int Ac = 0; Ac < tnNo; Ac++) {
        const int g = com_mod.ltg(Ac);
        if (cap_gn_set.find(g) == cap_gn_set.end()) {
            continue;
        }
        send_buf(idx++) = static_cast<double>(g);
        for (int i = 0; i < nsd; i++) {
            send_buf(idx++) = com_mod.x(i, Ac);
        }
        for (int i = 0; i < nsd; i++) {
            send_buf(idx++) = Do(i, Ac);
        }
        for (int i = 0; i < nsd; i++) {
            send_buf(idx++) = Dn(i, Ac);
        }
        if (gather_Y) {
            for (int i = 0; i < nsd; i++) {
                send_buf(idx++) = Yo(i, Ac);
            }
            for (int i = 0; i < nsd; i++) {
                send_buf(idx++) = Yn(i, Ac);
            }
        }
    }

    // Gather the cap global mesh state.
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

    // If not the root rank, return.
    if (cm.idcm() != root) {
        return;
    }

    // Initialize the cap global mesh state.
    out.x = 0.0;
    out.Do = 0.0;
    out.Dn = 0.0;
    if (gather_Y) {
        out.Yo = 0.0;
        out.Yn = 0.0;
    }

    // Unpack the cap global mesh state.
    int pos = 0;
    while (pos + per_node <= static_cast<int>(recv_buf.size())) {
        const int g = static_cast<int>(recv_buf(pos++));
        const auto it = g_to_cap_col.find(g);
        if (it == g_to_cap_col.end()) {
            pos += 3 * nsd + (gather_Y ? 2 * nsd : 0);
            continue;
        }
        const int col = it->second;
        for (int i = 0; i < nsd; i++) {
            out.x(i, col) = recv_buf(pos++);
        }
        for (int i = 0; i < nsd; i++) {
            out.Do(i, col) = recv_buf(pos++);
        }
        for (int i = 0; i < nsd; i++) {
            out.Dn(i, col) = recv_buf(pos++);
        }
        if (gather_Y) {
            for (int i = 0; i < nsd; i++) {
                out.Yo(i, col) = recv_buf(pos++);
            }
            for (int i = 0; i < nsd; i++) {
                out.Yn(i, col) = recv_buf(pos++);
            }
        }
    }
}

} // namespace

/// @brief Gathers cap-node mesh state (see \ref cap_mesh_global_node_ids_).
/// @param com_mod The com_mod object.
/// @param cm_mod The cm_mod object.
/// @param solutions Old/current displacement and velocity for gather.
/// @param gather_Y If true, gather Yo/Yn with \c nsd rows per node; if false, geometry only (x, Do, Dn).
void CoupledBoundaryCondition::gather_global_mesh_state(ComMod& com_mod, const CmMod& cm_mod,
                                                        const SolutionStates& solutions, bool gather_Y) const
{
    auto& cm = com_mod.cm;
    const int nsd = com_mod.nsd;
    const int tnNo = com_mod.tnNo;
    const int root = cm_mod.master;
    const int nProcs = cm.np();

    if (cm.seq()) {
        gather_global_mesh_state_serial(com_mod, solutions, gather_Y, nsd, tnNo, cap_mesh_global_node_ids_,
                                        cap_global_mesh_state_);
        return;
    }

    gather_global_mesh_state_parallel(com_mod, cm_mod, cm, solutions, gather_Y, nsd, tnNo, root, nProcs,
                                        cap_mesh_global_node_ids_, cap_g_to_cap_col_, cap_global_mesh_state_);
}

/// @brief Calculates the cap contribution to the linear solver face.
/// @param com_mod The com_mod object.
/// @param cm_mod The cm_mod object.
/// @param cfg_o The old mechanical configuration type.
/// @param cfg_n The new mechanical configuration type.
/// @return The cap contribution.
std::pair<double, double> CoupledBoundaryCondition::calculate_cap_contribution(ComMod& com_mod, const CmMod& cm_mod,
    const SolutionStates& solutions,
    consts::MechanicalConfigurationType cfg_o, consts::MechanicalConfigurationType cfg_n)
{
    auto& cm = com_mod.cm;
    double Qo_cap = 0.0;
    double Qn_cap = 0.0;
    const bool i_am_master = cm.mas(cm_mod);
    const bool serial_run = (cm.np() == 1);

    gather_global_mesh_state(com_mod, cm_mod, solutions, true);

    // Integrate the velocity flux over the cap surface (in master rank)
    if ((serial_run || i_am_master) && owns_cap_ && cap_ && cap_->face()) {
        Qo_cap = cap_->integrate_velocity_flux(cap_global_mesh_state_, false, cfg_o);
        Qn_cap = cap_->integrate_velocity_flux(cap_global_mesh_state_, true, cfg_n);
    }

    // Broadcast the cap contribution to all ranks.
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

/// @brief Loads the cap VTP file and validates the number of points and elements.
/// @param vtp_file_path The path to the VTP file.
/// @return The VTP data.
VtkVtpData load_cap_vtp(const std::string& vtp_file_path)
{
    VtkVtpData vtp_data;
    try {
        vtp_data = VtkVtpData(vtp_file_path, true);
    } catch (const std::exception& e) {
        throw CappingSurfaceVtpException("Failed to construct VtkVtpData from file '" + vtp_file_path +
                                         "': " + e.what());
    }

    int nNo = 0;
    try {
        nNo = vtp_data.num_points();
    } catch (const std::exception& e) {
        throw CappingSurfaceVtpException("Failed to get number of points from VTP file '" + vtp_file_path +
                                         "': " + e.what());
    }
    if (nNo == 0) {
        throw CappingSurfaceVtpException("Cap VTP file '" + vtp_file_path + "' does not contain any points.");
    }

    int num_elems = 0;
    try {
        num_elems = vtp_data.num_elems();
    } catch (const std::exception& e) {
        throw CappingSurfaceVtpException("Failed to get number of elements from VTP file '" + vtp_file_path +
                                         "': " + e.what());
    }
    if (num_elems == 0) {
        throw CappingSurfaceVtpException("Cap VTP file '" + vtp_file_path + "' does not contain any elements.");
    }

    bool has_global_node_id = false;
    try {
        has_global_node_id = vtp_data.has_point_data("GlobalNodeID");
    } catch (const std::exception& e) {
        throw CappingSurfaceVtpException("Failed to check for GlobalNodeID in VTP file '" + vtp_file_path +
                                         "': " + e.what());
    }
    if (!has_global_node_id) {
        throw CappingSurfaceVtpException("Cap VTP file '" + vtp_file_path +
                                         "' does not contain 'GlobalNodeID' point data.");
    }

    return vtp_data;
}

/// @brief Loads the cell "Normals" from the cap VTP file.
/// @param vtp_data The VTP data.
/// @param vtp_file_path The path to the VTP file.
/// @param num_elems The number of elements.
/// @param normals Output array for cell normals (num_comp x num_elems).
void load_cap_vtp_normals(VtkVtpData& vtp_data, const std::string& vtp_file_path, int num_elems,
                                    Array<double>& normals)
{
    try {
        bool has_normals = vtp_data.has_cell_data("Normals");
        if (has_normals) {
            auto [num_comp, num_tuples] = vtp_data.get_cell_data_dimensions("Normals");

            if (num_comp == 0 || num_tuples == 0) {
                std::string error_msg = "Normals array exists but has zero size. ";
                error_msg += "num_components=" + std::to_string(num_comp) + ", num_tuples=" + std::to_string(num_tuples);
                error_msg += ", expected num_elems=" + std::to_string(num_elems);
                error_msg += ". The array may not be a numeric type or may be empty.";
                throw CappingSurfaceVtpException(error_msg);
            }

            if (num_tuples != num_elems) {
                throw CappingSurfaceVtpException("Normals array size mismatch: " + std::to_string(num_tuples) +
                                                 " != " + std::to_string(num_elems));
            }

            if (num_comp != 2 && num_comp != 3) {
                throw CappingSurfaceVtpException("Invalid number of components in Normals array: " +
                                                 std::to_string(num_comp) + " (expected 2 or 3)");
            }

            normals.resize(num_comp, num_elems);
            vtp_data.copy_cell_data("Normals", normals);

            if (normals.nrows() != num_comp || normals.ncols() != num_elems) {
                throw CappingSurfaceVtpException("Failed to copy Normals data. Expected size: " +
                                                 std::to_string(num_comp) + "x" + std::to_string(num_elems) +
                                                 ", Actual size: " + std::to_string(normals.nrows()) + "x" +
                                                 std::to_string(normals.ncols()));
            }
        } else {
            normals.resize(0, 0);
        }
    } catch (const CappingSurfaceBaseException&) {
        throw;
    } catch (const std::exception& e) {
        throw CappingSurfaceVtpException("Failed to load Normals from VTP file '" + vtp_file_path + "': " + e.what());
    }
}

/// @brief Checks if the cap face shares any nodes with the coupled face.
/// @param coupled_face The coupled face.
/// @param cap_face The cap face.
/// @param vtp_file_path The path to the VTP file.
/// @param coupled_face_name The name of the coupled face.
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
    throw CappingSurfaceCouplingTopologyException(vtp_file_path, coupled_face_name);
}

} // namespace

CappingSurface::CappingSurface(const CappingSurface& other)
    : global_node_ids_(other.global_node_ids_)
    , valM_(other.valM_)
    , normals_(other.normals_)
{
    if (other.face_) {
        face_ = std::make_unique<faceType>(*other.face_);
    }
}

CappingSurface& CappingSurface::operator=(const CappingSurface& other)
{
    if (this != &other) {
        global_node_ids_ = other.global_node_ids_;
        valM_ = other.valM_;
        normals_ = other.normals_;
        if (other.face_) {
            face_ = std::make_unique<faceType>(*other.face_);
        } else {
            face_.reset();
        }
    }
    return *this;
}

/// @brief Loads the cap VTP file and creates the cap face.
/// @param vtp_file_path The path to the VTP file.
/// @param coupled_face The coupled face.
/// @param coupled_face_name The name of the coupled face.
void CappingSurface::load_from_vtp(const std::string& vtp_file_path, const faceType& coupled_face,
                                   const std::string& coupled_face_name)
{
    
    // Check if the VTP file exists.
    std::ifstream file_check(vtp_file_path);
    if (!file_check.good()) {
        throw CappingSurfaceFileException(vtp_file_path);
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
        throw CappingSurfaceVtpException("Failed to copy GlobalNodeID from VTP file '" + vtp_file_path + "': " +
                                         e.what());
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
        throw CappingSurfaceVtpException("Failed to get connectivity from VTP file '" + vtp_file_path + "': " +
                                         e.what());
    }
    if (vtk_cell_type != VTK_TRIANGLE) {
        throw CappingSurfaceUnsupportedCellException(vtk_cell_type);
    }
    if (eNoN != 3) {
        throw CappingSurfaceInvalidElementNodesException(eNoN, 3);
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
    load_cap_vtp_normals(vtp_data, vtp_file_path, num_elems, normals_);

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
            throw CappingSurfaceBaseException("[CappingSurface::init_cap_face_quadrature] Cap surface requires nsd=3.");
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
    } catch (const CappingSurfaceBaseException&) {
        throw;
    } catch (const std::exception& e) {
        throw CappingSurfaceQuadratureException(std::string(e.what()));
    }
}

/// @brief Initialize the cap contribution storage.
void CappingSurface::initialize_valM()
{
    if (!face_) {
        valM_.resize(0, 0);
        return;
    }
    valM_.resize(cap_nsd_, face_->nNo);
    valM_ = 0.0;
}

/// @brief Updates the element position in global coordinates.
/// @param e The element index.
/// @param cfg The mechanical configuration type.
/// @param mesh_x The mesh x coordinates.
/// @param mesh_Do The mesh Do coordinates.
/// @param mesh_Dn The mesh Dn coordinates.
/// @return The element position in global coordinates.
Array<double> CappingSurface::update_element_position_global(int e, consts::MechanicalConfigurationType cfg,
                                                             const Array<double>& mesh_x, const Array<double>& mesh_Do,
                                                             const Array<double>& mesh_Dn,
                                                             const std::unordered_map<int, int>& gn_to_cap_local) const
{
    using namespace consts;

    Array<double> xl(cap_nsd_, face_->eNoN);

    for (int a = 0; a < face_->eNoN; a++) {
        const int gn = face_->IEN(a, e);
        const int col = gn_to_cap_local.at(gn);
        for (int i = 0; i < cap_nsd_; i++) {
            xl(i, a) = mesh_x(i, col);
        }
        if (cfg == MechanicalConfigurationType::old_timestep) {
            for (int i = 0; i < cap_nsd_; i++) {
                xl(i, a) += mesh_Do(i, col);
            }
        } else if (cfg == MechanicalConfigurationType::new_timestep) {
            for (int i = 0; i < cap_nsd_; i++) {
                xl(i, a) += mesh_Dn(i, col);
            }
        }
    }

    return xl;
}

/// @brief Computes the Jacobian and normal vector for a given element and Gauss point.
/// @param xl The element position in global coordinates.
/// @param e The element index.
/// @param g The Gauss point index.
/// @return The Jacobian and normal vector.
std::pair<double, Vector<double>> CappingSurface::compute_jacobian_and_normal(const Array<double>& xl, int e, int g) const
{
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
        throw CappingSurfaceBaseException("[CappingSurface::compute_jacobian_and_normal] Zero Jacobian at Gauss point " +
                                              std::to_string(g));
    }

    n = n / Jac;

    // Check if the initial normals are provided and if they are valid.
    if (normals_.ncols() > 0 && normals_.nrows() == cap_nsd_) {

        Vector<double> n0(cap_nsd_);
        for (int i = 0; i < cap_nsd_; i++) {
            n0(i) = normals_(i, e);
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

/// @brief Integrates the velocity flux over the cap surface.
/// @param st The cap global mesh state.
/// @param use_Yn_velocity Whether to use the Yn velocity.
/// @param cfg The mechanical configuration type.
/// @return The velocity flux.
double CappingSurface::integrate_velocity_flux(const CapGlobalMeshState& st, bool use_Yn_velocity,
                                                 consts::MechanicalConfigurationType cfg)
{
    const Array<double>& cap_vel = use_Yn_velocity ? st.Yn : st.Yo;

    std::unordered_map<int, int> gn_to_cap_local;
    gn_to_cap_local.reserve(static_cast<size_t>(face_->nNo));
    for (int a = 0; a < face_->nNo; a++) {
        gn_to_cap_local[face_->gN(a)] = a;
    }

    double result = 0.0;
    for (int e = 0; e < face_->nEl; e++) {
        Array<double> xl = update_element_position_global(e, cfg, st.x, st.Do, st.Dn, gn_to_cap_local);
        for (int g = 0; g < face_->nG; g++) {
            auto [Jac, n] = compute_jacobian_and_normal(xl, e, g);
            double sHat = 0.0;
            for (int a = 0; a < face_->eNoN; a++) {
                const int gn = face_->IEN(a, e);
                const int col = gn_to_cap_local.at(gn);
                const double Na = face_->N(a, g);
                for (int i = 0; i < cap_nsd_; i++) {
                    sHat += Na * cap_vel(i, col) * n(i);
                }
            }
            result += face_->w(g) * Jac * sHat;
        }
    }
    return result;
}


/// @brief Computes the cap contribution to the linear solver face.
/// @param cfg The mechanical configuration type.
/// @param st The cap global mesh state.
void CappingSurface::compute_valM(consts::MechanicalConfigurationType cfg, const CapGlobalMeshState& st) const
{
    int cap_nNo = face_->nNo;

    valM_ = 0.0;

    // Map global node IDs to cap face-local indices
    std::unordered_map<int, int> gnNo_to_cap_local;
    for (int a = 0; a < cap_nNo; a++) {
        int gnNo = face_->gN(a);
        gnNo_to_cap_local[gnNo] = a;
    }

    for (int e = 0; e < face_->nEl; e++) {
        Array<double> xl = update_element_position_global(e, cfg, st.x, st.Do, st.Dn, gnNo_to_cap_local);
        for (int g = 0; g < face_->nG; g++) {
            auto [Jac, n] = compute_jacobian_and_normal(xl, e, g);
            for (int a = 0; a < face_->eNoN; a++) {
                const int gnNo_idx = face_->IEN(a, e);
                const int cap_a = gnNo_to_cap_local.at(gnNo_idx);
                for (int i = 0; i < cap_nsd_; i++) {
                    valM_(i, cap_a) += face_->N(a, g) * face_->w(g) * Jac * n(i);
                }
            }
        }
    }
}

namespace {


/// @brief Broadcasts the cap contribution to the linear solver face to all ranks.
/// @param com_mod The com_mod object.
/// @param cm_mod The cm_mod object.
/// @param lhs_face The linear solver face.
/// @param cap_nNo The number of cap nodes (same on all ranks: \c cap_n_no_ from \c distribute()).
/// @param cap_gN_all The global node IDs of the cap nodes.
/// @param cap_val_all The values of the cap nodes.
void bcast_cap_lhs_contribution(ComMod& com_mod, const CmMod& cm_mod,
                                fsi_linear_solver::FSILS_faceType& lhs_face, int cap_nNo,
                                Vector<int>& cap_gN_all, Array<double>& cap_val_all)
{
    auto& cm = com_mod.cm;
    const bool i_am_sender = cm.mas(cm_mod);
    int nsd = com_mod.nsd;

    if (!i_am_sender && cap_nNo > 0) {
        cap_gN_all.resize(cap_nNo);
        cap_val_all.resize(nsd, cap_nNo);
    }

    if (cap_nNo > 0) {
        cm.bcast(cm_mod, cap_gN_all);
        cm.bcast(cm_mod, cap_val_all);
    }

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
/// @param com_mod The com_mod object.
/// @param lhs_face The linear solver face.
/// @param cfg The mechanical configuration type.
void CoupledBoundaryCondition::copy_cap_surface_to_linear_solver_face(ComMod& com_mod,
                                                                      fsi_linear_solver::FSILS_faceType& lhs_face,
                                                                      consts::MechanicalConfigurationType cfg,
                                                                      const SolutionStates& solutions) const
{
    const int nsd = com_mod.nsd;

    const int cap_nNo = cap_n_no_;
    Vector<int> cap_gN_all;
    Array<double> cap_val_all;

    // Calculate the cap contribution to the linear solver face (in master rank)
    if (owns_cap_ && cap_) {
        const faceType* cap_face = cap_->face();
        cap_->compute_valM(cfg, cap_global_mesh_state_);
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
    bcast_cap_lhs_contribution(com_mod, cm_mod_, lhs_face, cap_nNo, cap_gN_all, cap_val_all);
}

/// @brief Broadcasts the coupled Neumann pressure to all ranks.
/// @param cm_mod The cm_mod object.
/// @param cm The cm object.
void CoupledBoundaryCondition::bcast_coupled_neumann_pressure(const CmMod& cm_mod, cmType& cm)
{
    if (cm.seq()) {
        return;
    }
    using namespace consts;
    if (get_bc_type() != BoundaryConditionType::bType_Neu) {
        return;
    }
    double pr = 0.0;
    if (cm.mas(cm_mod)) {
        pr = get_pressure();
    }
    cm.bcast(cm_mod, &pr);
    set_pressure(pr);
}
