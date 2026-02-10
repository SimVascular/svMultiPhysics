// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef ZEROD_BOUNDARY_CONDITION_H
#define ZEROD_BOUNDARY_CONDITION_H

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
    
    /// @brief Cap surface data (for caps that don't share elements with mesh)
    bool cap_loaded_ = false;                ///< Whether cap VTP has been loaded
    std::unique_ptr<faceType> cap_face_;     ///< Temporary face structure for cap integration
    Vector<int> cap_global_node_ids_;        ///< GlobalNodeID from cap VTP file
    bool cap_area_computed_ = false;         ///< Whether cap area has been computed and printed
    std::unordered_map<int, int> cap_gnNo_to_tnNo_;  ///< Mapping from global node numbers to total node numbers for cap
    Array<double> cap_valM_;                 ///< Precomputed cap normal integrals (nsd x cap_face->nNo), stored with cap face-local indices
                                            ///< Similar to face.valM which is stored with face-local indices
                                            ///< If size is 0, no cap contribution (either no cap or not computed)
    Array<double> cap_initial_normals_;     ///< Initial element normals from VTP file (nsd x num_elems)
                                            ///< Used to ensure calculated normals point in the same direction

public:
    /// @brief Default constructor - creates an uninitialized object
    ZeroDBoundaryCondition() = default;
    
    /// @brief Copy constructor
    ZeroDBoundaryCondition(const ZeroDBoundaryCondition& other);
    
    /// @brief Copy assignment operator
    ZeroDBoundaryCondition& operator=(const ZeroDBoundaryCondition& other);
    
    /// @brief Move constructor
    ZeroDBoundaryCondition(ZeroDBoundaryCondition&& other) noexcept;
    
    /// @brief Move assignment operator
    ZeroDBoundaryCondition& operator=(ZeroDBoundaryCondition&& other) noexcept;

    /// @brief Construct with a face association (no VTP data loaded)
    /// @param bc_type The 3D boundary condition type (must be bType_Dir or bType_Neu)
    /// @param face Face associated with this BC
    /// @param face_name Face name from the mesh
    /// @param block_name Block name in svZeroDSolver configuration
    /// @param logger Simulation logger used to write warnings
    ZeroDBoundaryCondition(consts::BoundaryConditionType bc_type, const faceType& face, const std::string& face_name,
                          const std::string& block_name, SimulationLogger& logger);

    /// @brief Construct and optionally point to a cap face VTP file
    /// @param bc_type The 3D boundary condition type (must be bType_Dir or bType_Neu)
    /// @param face Face associated with this BC
    /// @param face_name Face name from the mesh
    /// @param block_name Block name in svZeroDSolver configuration
    /// @param cap_face_vtp_file Path to the cap face VTP file
    /// @param logger Simulation logger used to write warnings
    ZeroDBoundaryCondition(consts::BoundaryConditionType bc_type, const faceType& face, const std::string& face_name,
                          const std::string& block_name, const std::string& cap_face_vtp_file, SimulationLogger& logger);

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
    
    /// @brief Initialize cap integration requirements (Gauss points, weights, shape functions)
    /// @param com_mod ComMod reference containing simulation data
    /// @param cm_mod CmMod reference for communication
    /// @note This should be called after the cap is loaded and before any integration operations
    void initialize_cap_integration(ComMod& com_mod, const CmMod& cm_mod);
    
    /// @brief Check if cap is loaded and valid
    /// @return true if cap is loaded, false otherwise
    bool has_cap() const { return cap_loaded_ && cap_face_ != nullptr; }
    
    /// @brief Compute the area of the cap surface
    /// @param com_mod ComMod reference containing simulation data
    /// @param cm_mod CmMod reference for communication
    /// @return Area of the cap surface (0.0 if no cap is loaded)
    double compute_cap_area(ComMod& com_mod, const CmMod& cm_mod);
    
    /// @brief Compute and update cap_valM (precomputed cap normal integrals)
    /// Similar to fsi_ls_upd, this computes ∫ N_A n_i dΓ over the cap surface
    /// and stores the result in cap_valM_ for efficient use in matrix-vector products
    /// @param com_mod ComMod reference containing simulation data
    /// @param cm_mod CmMod reference for communication
    /// @param cfg Mechanical configuration type (reference, old_timestep, or new_timestep)
    void compute_cap_valM(ComMod& com_mod, const CmMod& cm_mod, consts::MechanicalConfigurationType cfg);
    
    /// @brief Get the precomputed cap_valM array
    /// @return Reference to cap_valM_ array (nsd x cap_face->nNo)
    /// @note Check cap_valM_.size() > 0 to see if it has been computed and contains data
    const Array<double>& get_cap_valM() const { return cap_valM_; }
    
    /// @brief Get cap face structure (for accessing cap_face->nNo and cap_face->gN)
    /// @return Pointer to cap_face_ (nullptr if no cap)
    const faceType* get_cap_face() const { return cap_face_.get(); }
    
    /// @brief Copy cap data to linear solver face structure for use in preconditioner
    /// This method computes cap_valM if not already done, then copies it to the linear solver
    /// face structure's cap_val field, initializes cap_valM to zero, and sets up cap_glob mapping
    /// @param com_mod ComMod reference containing simulation data
    /// @param cm_mod CmMod reference for communication
    /// @param lhs_face Reference to the linear solver face structure (FSILS_faceType)
    /// @param cfg Mechanical configuration type (reference, old_timestep, or new_timestep)
    void copy_cap_data_to_linear_solver_face(ComMod& com_mod, const CmMod& cm_mod, 
                                              fsi_linear_solver::FSILS_faceType& lhs_face,
                                              consts::MechanicalConfigurationType cfg);
    
    // =========================================================================
    // Cap surface integration (private helper methods)
    // =========================================================================
    
private:
    /// @brief Integrate a variable over the cap surface
    /// @param com_mod ComMod reference containing simulation data
    /// @param cm_mod CmMod reference for communication
    /// @param s Array containing variable values (nsd x tnNo for vectors, or scalar)
    /// @param l Lower index of s
    /// @param u Upper index of s (optional, defaults to l)
    /// @param cfg Mechanical configuration type
    /// @return Integral value
    double integrate_over_cap(ComMod& com_mod, const CmMod& cm_mod, const Array<double>& s, 
                               int l, std::optional<int> u, consts::MechanicalConfigurationType cfg);
    
    /// @brief Update cap element nodal coordinates based on configuration
    /// @param com_mod ComMod reference containing simulation data
    /// @param e Element index
    /// @param gnNo_to_tnNo Mapping from global node numbers to total node numbers
    /// @param cfg Mechanical configuration type
    /// @return Element nodal coordinates in the specified configuration
    Array<double> update_cap_element_position(ComMod& com_mod, int e, 
                                              const std::unordered_map<int, int>& gnNo_to_tnNo,
                                              consts::MechanicalConfigurationType cfg);
    
    /// @brief Compute Jacobian and normal vector at a Gauss point
    /// @param xl Element nodal coordinates
    /// @param e Element index (for accessing initial normal)
    /// @param g Gauss point index
    /// @param nsd Number of spatial dimensions
    /// @param insd Integration space dimensions
    /// @return Pair of (Jacobian, normalized normal vector)
    std::pair<double, Vector<double>> compute_cap_jacobian_and_normal(const Array<double>& xl, 
                                                                       int e, int g, int nsd, int insd);
    
    /// @brief Integrate scalar field over cap element at a Gauss point
    /// @param com_mod ComMod reference containing simulation data
    /// @param s Array containing scalar values
    /// @param l Index of scalar in s
    /// @param e Element index
    /// @param g Gauss point index
    /// @param gnNo_to_tnNo Mapping from global node numbers to total node numbers
    /// @return Scalar integrand value
    double integrate_scalar_at_gauss_point(ComMod& com_mod, const Array<double>& s, int l,
                                           int e, int g, const std::unordered_map<int, int>& gnNo_to_tnNo);
    
    /// @brief Integrate vector field over cap element at a Gauss point
    /// @param com_mod ComMod reference containing simulation data
    /// @param s Array containing vector values
    /// @param l Lower index of vector in s
    /// @param nsd Number of spatial dimensions
    /// @param e Element index
    /// @param g Gauss point index
    /// @param n Normalized normal vector
    /// @param gnNo_to_tnNo Mapping from global node numbers to total node numbers
    /// @return Vector integrand value (s dot n)
    double integrate_vector_at_gauss_point(ComMod& com_mod, const Array<double>& s, int l, int nsd,
                                           int e, int g, const Vector<double>& n,
                                           const std::unordered_map<int, int>& gnNo_to_tnNo);
};

#endif // ZEROD_BOUNDARY_CONDITION_H
