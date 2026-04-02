// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef COUPLED_BOUNDARY_CONDITION_H
#define COUPLED_BOUNDARY_CONDITION_H

#include <string>
#include <vector>
#include <memory>
#include <optional>
#include <unordered_map>
#include <utility>
#include "consts.h"
#include "Array.h"
#include "Vector.h"

// Forward declarations to avoid heavy includes
class LPNSolverInterface;
class faceType;
class ComMod;
class CmMod;
class SimulationLogger;

namespace fsi_linear_solver {
    class FSILS_faceType;
}

class CappingSurface;

/// @brief Object-oriented Coupled boundary condition on a cap face
///
/// This class provides an interface for:
///  - loading a cap face VTP,
///  - computing flowrates on the face for coupling, and
///  - getting/setting pressure values from/to a 0D solver.
///
/// The class manages its own coupling data. svZeroD_subroutines accesses
/// coupled boundary conditions by iterating through com_mod.eq[].bc[].
class CoupledBoundaryCondition {
protected:
    /// @brief Data members for BC
    const faceType* face_ = nullptr;         ///< Face associated with the BC (not owned by CoupledBoundaryCondition)
    std::string cap_face_vtp_file_;          ///< Path to VTP file (empty if no cap)
    const SimulationLogger* logger_ = nullptr;  ///< Logger for warnings/info (not owned by CoupledBoundaryCondition)

    /// @brief 3D boundary condition type (Dirichlet or Neumann) for this Coupled BC.
    consts::BoundaryConditionType bc_type_ = consts::BoundaryConditionType::bType_Neu;

    /// @brief svZeroD coupling data
    std::string block_name_;                 ///< Block name in svZeroDSolver configuration
    std::string face_name_;                  ///< Face name from the mesh
    
    /// @brief Flowrate data
    double Qo_ = 0.0;                        ///< Flowrate at old timestep (t_n)
    double Qn_ = 0.0;                        ///< Flowrate at new timestep (t_{n+1})
    
    /// @brief Pressure data  
    double Po_ = 0.0;                        ///< Pressure at old timestep (for completeness)
    double Pn_ = 0.0;                        ///< Pressure at new timestep (for completeness)
    double pressure_ = 0.0;                  ///< Current pressure value from 0D solver (result)
    
    /// @brief svZeroD solution IDs
    int flow_sol_id_ = -1;                   ///< ID in svZeroD solution vector for flow
    int pressure_sol_id_ = -1;               ///< ID in svZeroD solution vector for pressure
    double in_out_sign_ = 1.0;               ///< Sign for inlet/outlet (+1 inlet to 0D model, -1 outlet)
    
    /// @brief Configuration for flowrate computation
    bool follower_pressure_load_ = false;   ///< Whether to use follower pressure load (for struct/ustruct)
    
    /// @brief True if this BC uses a chamber cap (broadcast in distribute so all ranks agree).
    bool has_cap_ = false;
    /// @brief Cap geometry and integration; set when a cap VTP path is given and load succeeds.
    std::unique_ptr<CappingSurface> cap_;

public:
    /// @brief Default constructor - creates an uninitialized object
    CoupledBoundaryCondition() = default;

    /// @brief Destructor (out-of-line so std::unique_ptr<CappingSurface> is valid with forward declare)
    ~CoupledBoundaryCondition();
    
    /// @brief Copy constructor
    CoupledBoundaryCondition(const CoupledBoundaryCondition& other);
    
    /// @brief Copy assignment operator
    CoupledBoundaryCondition& operator=(const CoupledBoundaryCondition& other);
    
    /// @brief Move constructor
    CoupledBoundaryCondition(CoupledBoundaryCondition&& other) noexcept;
    
    /// @brief Move assignment operator
    CoupledBoundaryCondition& operator=(CoupledBoundaryCondition&& other) noexcept;

    /// @brief Construct with a face association (no VTP data loaded)
    /// @param bc_type The 3D boundary condition type (must be bType_Dir or bType_Neu)
    /// @param face Face associated with this BC
    /// @param face_name Face name from the mesh
    /// @param block_name Block name in svZeroDSolver configuration
    /// @param logger Simulation logger used to write warnings
    CoupledBoundaryCondition(consts::BoundaryConditionType bc_type, const faceType& face, const std::string& face_name,
                          const std::string& block_name, SimulationLogger& logger);

    /// @brief Construct and optionally point to a cap face VTP file
    /// @param bc_type The 3D boundary condition type (must be bType_Dir or bType_Neu)
    /// @param face Face associated with this BC
    /// @param face_name Face name from the mesh
    /// @param block_name Block name in svZeroDSolver configuration
    /// @param cap_face_vtp_file Path to the cap face VTP file
    /// @param logger Simulation logger used to write warnings
    CoupledBoundaryCondition(consts::BoundaryConditionType bc_type, const faceType& face, const std::string& face_name,
                          const std::string& block_name, const std::string& cap_face_vtp_file, SimulationLogger& logger);

    /// @brief Get the 3D BC type for this Coupled boundary condition.
    consts::BoundaryConditionType get_bc_type() const { return bc_type_; }

    /// @brief Load the cap face VTP file and associate it with this boundary condition
    /// @param vtp_file_path Path to the cap face VTP file
    void load_cap_face_vtp(const std::string& vtp_file_path);

    // =========================================================================
    // svZeroD block configuration
    // =========================================================================
    
    /// @brief Set the svZeroD block name
    /// @param block_name Block name in svZeroDSolver configuration
    void set_block_name(const std::string& block_name);
    
    /// @brief Get the svZeroD block name
    /// @return Block name
    const std::string& get_block_name() const;
    
    /// @brief Set the face name
    /// @param face_name Face name from the mesh
    void set_face_name(const std::string& face_name);
    
    /// @brief Set the svZeroD solution IDs for flow and pressure
    /// @param flow_id Flow solution ID
    /// @param pressure_id Pressure solution ID
    /// @param in_out_sign Sign for inlet/outlet
    void set_solution_ids(int flow_id, int pressure_id, double in_out_sign);
    
    /// @brief Get the flow solution ID
    int get_flow_sol_id() const;
    
    /// @brief Get the pressure solution ID
    int get_pressure_sol_id() const;
    
    /// @brief Get the inlet/outlet sign
    double get_in_out_sign() const;

    // =========================================================================
    // Flowrate computation and access
    // =========================================================================

    /// @brief Set the follower pressure load flag
    /// @param flwP Whether to use follower pressure load
    void set_follower_pressure_load(bool flwP);
    
    /// @brief Get the follower pressure load flag
    /// @return Whether to use follower pressure load
    bool get_follower_pressure_load() const;

    /// @brief Compute flowrates at the boundary face at old and new timesteps
    /// @param com_mod ComMod reference containing simulation data
    /// @param cm_mod CmMod reference for communication
    /// @param phys Current physics type (struct, ustruct, fluid, etc.)
    void compute_flowrates(ComMod& com_mod, const CmMod& cm_mod, consts::EquationType phys);

    /// @brief Gather cap surface nodal data for integration (serial path or MPI bcast/gatherv).
    /// Call on every rank when has_cap() so cap MPI collectives match; no-op when cap disabled.
    void prepare_cap_gathered_data(ComMod& com_mod, const CmMod& cm_mod,
                                   const Array<double>& Yo, const Array<double>& Yn,
                                   int l, int s_comps, consts::MechanicalConfigurationType cfg_o,
                                   consts::MechanicalConfigurationType cfg_n);

    /// @brief Extra volumetric flux through the cap (old/new timestep); {0,0} if no cap; MPI-safe on all ranks.
    std::pair<double, double> calculate_cap_contribution(ComMod& com_mod, const CmMod& cm_mod, int nsd,
                                                         consts::MechanicalConfigurationType cfg_o,
                                                         consts::MechanicalConfigurationType cfg_n);
    
    /// @brief Compute average pressures at the boundary face at old and new timesteps (for Dirichlet BCs)
    /// @param com_mod ComMod reference containing simulation data
    /// @param cm_mod CmMod reference for communication
    void compute_pressures(ComMod& com_mod, const CmMod& cm_mod);

    /// @brief Get the flowrate at old timestep
    /// @return Flowrate at t_n
    double get_Qo() const;
    
    /// @brief Get the flowrate at new timestep
    /// @return Flowrate at t_{n+1}
    double get_Qn() const;
    
    /// @brief Set the flowrates directly
    /// @param Qo Flowrate at old timestep
    /// @param Qn Flowrate at new timestep
    void set_flowrates(double Qo, double Qn);
    
    /// @brief Perturb the new timestep flowrate by a given amount
    /// @param diff Perturbation to add to Qn
    void perturb_flowrate(double diff);

    // =========================================================================
    // Pressure access (result from 0D solver)
    // =========================================================================

    /// @brief Set the pressure value from 0D solver
    /// @param pressure Pressure value to be applied as Neumann BC
    void set_pressure(double pressure);
    
    /// @brief Get the current pressure value
    /// @return Current pressure value from 0D solver
    double get_pressure() const;
    
    /// @brief Get the pressure at old timestep
    /// @return Pressure at t_n
    double get_Po() const;
    
    /// @brief Get the pressure at new timestep
    /// @return Pressure at t_{n+1}
    double get_Pn() const;
    
    // =========================================================================
    // State management for derivative computation
    // =========================================================================
    
    /// @brief State struct for saving/restoring Qn and pressure
    struct State {
        double Qn = 0.0;
        double pressure = 0.0;
    };
    
    /// @brief Save current state (Qn and pressure)
    /// @return Current state
    State save_state() const;
    
    /// @brief Restore state from a saved state
    /// @param state State to restore
    void restore_state(const State& state);

    // =========================================================================
    // Utility methods
    // =========================================================================

    /// @brief Distribute BC metadata from master to slave processes
    /// @param com_mod Reference to ComMod object
    /// @param cm_mod Reference to CmMod object for MPI communication
    /// @param cm Reference to cmType object for MPI communication
    /// @param face Face associated with the BC (after distribution)
    void distribute(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, const faceType& face);
    
    /// @brief Set the associated face (used during distribution)
    /// @param face Reference to the face
    void set_face(const faceType& face);
    
    /// @brief Get the associated face
    /// @return Pointer to the associated face (may be nullptr if not set)
    const faceType* get_face() const;
    
    /// @brief Check if the BC is properly initialized
    /// @return true if face is set, false otherwise
    bool is_initialized() const;
    
    /// @brief Check if this BC has a cap (broadcast in distribute so all ranks agree).
    bool has_cap() const { return has_cap_; }

    /// @brief Cap object on this rank after load; may be nullptr while has_cap() is true (e.g. rank without mesh).
    CappingSurface* capping_surface() noexcept { return cap_.get(); }
    const CappingSurface* capping_surface() const noexcept { return cap_.get(); }

    /// @brief True if this rank holds loaded cap surface data (mesh / quadrature).
    bool cap_face_ready() const;
};


/// @brief Capping surface geometry and integration for a coupled boundary.
///
/// Loaded from a cap VTP (GlobalNodeID connectivity), used for extra surface
/// flux contribution and linear-operator cap blocks. Definitions in CoupledBoundaryCondition.cpp.
class CappingSurface {
public:
    CappingSurface() = default;
    CappingSurface(const CappingSurface& other);
    CappingSurface& operator=(const CappingSurface& other);
    CappingSurface(CappingSurface&& other) noexcept = default;
    CappingSurface& operator=(CappingSurface&& other) noexcept = default;

    void load_from_vtp(const std::string& vtp_file_path, const faceType& coupled_face,
                       const std::string& coupled_face_name);

    void initialize_integration(ComMod& com_mod, const CmMod& cm_mod);

    bool is_ready() const { return face_ != nullptr; }

    const faceType* face() const { return face_.get(); }
    faceType* face() { return face_.get(); }

    const Array<double>& valM() const { return valM_; }

    void prepare_gathered_data(ComMod& com_mod, const CmMod& cm_mod,
                               const Array<double>& Yo, const Array<double>& Yn,
                               int l, int s_comps, consts::MechanicalConfigurationType cfg_o,
                               consts::MechanicalConfigurationType cfg_n);

    double integrate_over(ComMod& com_mod, const CmMod& cm_mod, const Array<double>& s,
                          int l, std::optional<int> u, consts::MechanicalConfigurationType cfg);

    void compute_valM(ComMod& com_mod, const CmMod& cm_mod, consts::MechanicalConfigurationType cfg);

    void copy_to_linear_solver_face(ComMod& com_mod, const CmMod& cm_mod,
                                    fsi_linear_solver::FSILS_faceType& lhs_face,
                                    consts::MechanicalConfigurationType cfg);

private:
    std::unique_ptr<faceType> face_;
    Vector<int> global_node_ids_;
    bool area_computed_ = false;
    std::unordered_map<int, int> gnNo_to_tnNo_;
    Array<double> valM_;
    Array<double> initial_normals_;

    Array<double> x_gathered_;
    Array<double> Do_gathered_;
    Array<double> Dn_gathered_;
    Array<double> Yo_gathered_;
    Array<double> Yn_gathered_;
    int nNo_gathered_ = 0;
    Vector<int> gN_broadcast_;

    void gather_node_data_to_master(ComMod& com_mod, const CmMod& cm_mod,
                                    const Array<double>& s_old, const Array<double>& s_new,
                                    int l, int s_comps, consts::MechanicalConfigurationType cfg_o,
                                    consts::MechanicalConfigurationType cfg_n,
                                    int cap_nNo, int nsd,
                                    Array<double>& cap_x, Array<double>& cap_Do, Array<double>& cap_Dn,
                                    Array<double>& cap_Yo, Array<double>& cap_Yn);

    Array<double> update_element_position(ComMod& com_mod, int e,
                                          const std::unordered_map<int, int>& gnNo_to_tnNo,
                                          consts::MechanicalConfigurationType cfg);

    Array<double> update_element_position(int e, consts::MechanicalConfigurationType cfg,
                                          const Array<double>& cap_x, const Array<double>& cap_Do,
                                          const Array<double>& cap_Dn,
                                          const std::unordered_map<int, int>& gnNo_to_capIdx);

    std::pair<double, Vector<double>> compute_jacobian_and_normal(const Array<double>& xl,
                                                                  int e, int g, int nsd, int insd);
};

#endif // COUPLED_BOUNDARY_CONDITION_H
