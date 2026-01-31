// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef ZEROD_BOUNDARY_CONDITION_H
#define ZEROD_BOUNDARY_CONDITION_H

#include <string>
#include <vector>
#include "consts.h"
#include "Array.h"
#include "Vector.h"

// Forward declarations to avoid heavy includes
class LPNSolverInterface;
class faceType;
class ComMod;
class CmMod;
class SimulationLogger;

/// @brief Object-oriented 0D-coupled boundary condition on a cap face
///
/// This class provides an interface for:
///  - loading a cap face VTP,
///  - computing flowrates on the face for coupling, and
///  - getting/setting pressure values from/to a 0D solver.
///
/// The class manages its own coupling data. svZeroD_subroutines accesses
/// ZeroD boundary conditions by iterating through com_mod.eq[].bc[].
class ZeroDBoundaryCondition {
protected:
    /// @brief Data members for BC
    const faceType* face_ = nullptr;         ///< Face associated with the BC (not owned by ZeroDBoundaryCondition)
    std::string cap_face_vtp_file_;          ///< Path to VTP file (empty if no cap)
    const SimulationLogger* logger_ = nullptr;  ///< Logger for warnings/info (not owned by ZeroDBoundaryCondition)

    /// @brief 3D boundary condition type (Dirichlet or Neumann) for this 0D-coupled BC.
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
    double in_out_sign_ = 1.0;               ///< Sign for inlet/outlet (+1 inlet to LPN, -1 outlet)
    
    /// @brief Configuration for flowrate computation
    bool follower_pressure_load_ = false;   ///< Whether to use follower pressure load (for struct/ustruct)

public:
    /// @brief Default constructor - creates an uninitialized object
    ZeroDBoundaryCondition() = default;

    /// @brief Construct with a face association (no VTP data loaded)
    /// @param bc_type The 3D boundary condition type (must be bType_Dir or bType_Neu)
    /// @param face Face associated with this BC
    /// @param logger Simulation logger used to write warnings
    ZeroDBoundaryCondition(consts::BoundaryConditionType bc_type, const faceType& face, SimulationLogger& logger);

    /// @brief Construct and optionally point to a cap face VTP file
    /// @param bc_type The 3D boundary condition type (must be bType_Dir or bType_Neu)
    /// @param cap_face_vtp_file Path to the cap face VTP file
    /// @param face Face associated with this BC
    /// @param logger Simulation logger used to write warnings
    ZeroDBoundaryCondition(consts::BoundaryConditionType bc_type, const std::string& cap_face_vtp_file, const faceType& face, SimulationLogger& logger);

    /// @brief Get the 3D BC type for this 0D-coupled boundary condition.
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
    
    /// @brief Get the face name
    /// @return Face name
    const std::string& get_face_name() const;
    
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

    /// @brief Get the associated face
    /// @return Pointer to the associated face (may be nullptr if not set)
    const faceType* get_face() const;
    
    /// @brief Check if the BC is properly initialized
    /// @return true if face is set, false otherwise
    bool is_initialized() const;
};

#endif // ZEROD_BOUNDARY_CONDITION_H
