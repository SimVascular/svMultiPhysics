// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef COUPLED_BOUNDARY_CONDITION_H
#define COUPLED_BOUNDARY_CONDITION_H

#include <string>
#include <memory>
#include <exception>
#include <optional>
#include <utility>
#include "consts.h"
#include "CmMod.h"
#include "SolutionStates.h"

// Forward declarations to avoid heavy includes
class LPNSolverInterface;
class faceType;
class ComMod;

namespace fsi_linear_solver {
    class FSILS_faceType;
}

/// @brief Mesh-wide nodal state (columns indexed by com_mod.ltg / gtnNo).
struct CapGlobalMeshState {
    int gtnNo = 0;
    Array<double> x;
    Array<double> Do;
    Array<double> Dn;
    Array<double> Yo;
    Array<double> Yn;

    void clear()
    {
        gtnNo = 0;
        x.resize(0, 0);
        Do.resize(0, 0);
        Dn.resize(0, 0);
        Yo.resize(0, 0);
        Yn.resize(0, 0);
    }
};

/// @brief Base exception for capping surface (cap VTP) errors.
///
/// These indicate fatal errors while loading or using a cap surface. They are not
/// expected to be recovered; callers may catch \ref CappingSurfaceBaseException to
/// handle all cap-related failures.
class CappingSurfaceBaseException : public std::exception {
public:
    explicit CappingSurfaceBaseException(std::string msg) : message_(std::move(msg)) {}

    const char* what() const noexcept override { return message_.c_str(); }

private:
    std::string message_;
};

/// @brief Cap VTP file cannot be opened.
class CappingSurfaceFileException : public CappingSurfaceBaseException {
public:
    explicit CappingSurfaceFileException(const std::string& path)
        : CappingSurfaceBaseException("[CappingSurface::load_from_vtp] Cannot open cap VTP file '" + path +
                                      "' for reading.") {}
};

/// @brief VTP read, parse, or validation error during cap load.
class CappingSurfaceVtpException : public CappingSurfaceBaseException {
public:
    explicit CappingSurfaceVtpException(const std::string& detail)
        : CappingSurfaceBaseException("[CappingSurface::load_from_vtp] " + detail) {}
};

/// @brief Cap mesh shares no nodes with the coupled boundary face.
class CappingSurfaceCouplingTopologyException : public CappingSurfaceBaseException {
public:
    explicit CappingSurfaceCouplingTopologyException(const std::string& vtp_path, const std::string& coupled_face_name)
        : CappingSurfaceBaseException("[CappingSurface::load_from_vtp] Cap VTP file '" + vtp_path +
                                      "' has no GlobalNodeID entries in common with coupled face '" +
                                      coupled_face_name + "'. The cap must share at least one mesh node with that face.") {}
};

/// @brief Cap VTP uses an unsupported cell type (only TRI3 is supported).
class CappingSurfaceUnsupportedCellException : public CappingSurfaceBaseException {
public:
    explicit CappingSurfaceUnsupportedCellException(int vtk_cell_type)
        : CappingSurfaceBaseException("[CappingSurface::load_from_vtp] Unsupported cap cell type " +
                                      std::to_string(vtk_cell_type) + ". Only VTK_TRIANGLE (TRI3) is supported.") {}
};

/// @brief Cap VTP connectivity does not match expected TRI3 topology.
class CappingSurfaceInvalidElementNodesException : public CappingSurfaceBaseException {
public:
    explicit CappingSurfaceInvalidElementNodesException(int eNoN, int expected)
        : CappingSurfaceBaseException("[CappingSurface::load_from_vtp] Invalid nodes-per-element for triangle cap: " +
                                      std::to_string(eNoN) + " (expected " + std::to_string(expected) + ").") {}
};

/// @brief Geometry or Jacobian error during cap surface integration.
class CappingSurfaceGeometryException : public CappingSurfaceBaseException {
public:
    explicit CappingSurfaceGeometryException(const std::string& detail) : CappingSurfaceBaseException(detail) {}
};

/// @brief Cap face quadrature (shape functions on TRI3) setup failed.
class CappingSurfaceQuadratureException : public CappingSurfaceBaseException {
public:
    explicit CappingSurfaceQuadratureException(const std::string& nested)
        : CappingSurfaceBaseException("[CappingSurface::init_cap_face_quadrature] Failed to initialize cap face shape "
                                      "functions: " + nested) {}
};

/// @brief Inconsistent cap connectivity or assembly indexing.
class CappingSurfaceAssemblyException : public CappingSurfaceBaseException {
public:
    explicit CappingSurfaceAssemblyException(const std::string& detail) : CappingSurfaceBaseException(detail) {}
};

/// @brief Capping surface geometry and integration for a coupled boundary.
class CappingSurface {
    public:
        /// @brief Default constructor.
        CappingSurface() = default;
        /// @brief Copy constructor.
        CappingSurface(const CappingSurface& other);
        /// @brief Copy assignment operator.
        CappingSurface& operator=(const CappingSurface& other);
        /// @brief Move constructor.
        CappingSurface(CappingSurface&& other) noexcept = default;
        /// @brief Move assignment operator.
        CappingSurface& operator=(CappingSurface&& other) noexcept = default;
    
        /// @brief Load the cap face from a VTP file.
        void load_from_vtp(const std::string& vtp_file_path, const faceType& coupled_face,
                           const std::string& coupled_face_name);
    
        /// @brief Initialize the cap face quadrature.
        void init_cap_face_quadrature(const ComMod& com_mod);

        /// @brief Initialize the cap contribution storage.
        void initialize_valM();
    
        /// Surface velocity flux through the cap using \a st columns indexed by cap IEN / GlobalNodeID (master / serial).
        double integrate_velocity_flux(const CapGlobalMeshState& st, bool use_Yn_velocity,
            consts::MechanicalConfigurationType cfg);

        /// @brief Compute the cap contribution to the linear solver face.
        void compute_valM(consts::MechanicalConfigurationType cfg, const CapGlobalMeshState& st);

        /// @brief Get the cap face.
        faceType* face() { return face_.get(); }

        /// @brief Get the cap contribution.
        const Array<double>& valM() const { return valM_; }

    private:
        /// @brief The cap face.
        std::unique_ptr<faceType> face_;
        /// @brief The global node IDs.
        Vector<int> global_node_ids_;
        /// @brief The normals.
        Array<double> normals_;

        /// @brief The number of spatial dimensions (3D).
        static constexpr int cap_nsd_ = 3;  
        /// @brief The number of independent spatial dimensions (2D).
        static constexpr int cap_insd_ = 2; 

        /// @brief The cap contribution.
        Array<double> valM_;

        /// @brief Update the element position in global coordinates.
        Array<double> update_element_position_global(int e, consts::MechanicalConfigurationType cfg,
                                                     const Array<double>& mesh_x, const Array<double>& mesh_Do,
                                                     const Array<double>& mesh_Dn, int gtnNo_cols) const;
    
        /// @brief Compute the Jacobian and normal vector for a given element and Gauss point.
        std::pair<double, Vector<double>> compute_jacobian_and_normal(const Array<double>& xl, int e, int g);
    };

/// @brief Object-oriented Coupled boundary condition on a cap face
///
/// This class provides an interface for:
///  - loading a cap face VTP,
///  - computing flowrates on the face for coupling, and
///  - getting/setting pressure values from/to a 0D solver.
///
/// The class manages its own coupling data. svZeroD interface code accesses
/// coupled boundary conditions by iterating through com_mod.eq[].bc[].
class CoupledBoundaryCondition {
private:
    /// @brief Data members for BC
    const faceType* face_ = nullptr;         ///< Face associated with the BC (not owned by CoupledBoundaryCondition)
    std::string cap_face_vtp_file_;          ///< Path to VTP file (empty if no cap)

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
    consts::EquationType phys_ = consts::EquationType::phys_NA;  ///< Equation physics for this coupled BC (set at construction)
    consts::MechanicalConfigurationType flowrate_cfg_o_ = consts::MechanicalConfigurationType::reference;
    consts::MechanicalConfigurationType flowrate_cfg_n_ = consts::MechanicalConfigurationType::reference;
    
    /// @brief True if this BC uses a chamber cap (broadcast in distribute so all ranks agree).
    bool has_cap_ = false;
    /// @brief True on ranks that hold \ref cap_ mesh/quadrature (MPI master when \ref has_cap_; true in serial when cap loaded).
    bool owns_cap_ = false;
    /// @brief Global mesh node IDs for cap surface nodes (same on all ranks after distribute; master fills at load).
    Vector<int> cap_mesh_global_node_ids_;
    /// @brief Cap geometry on ranks with \ref owns_cap_; empty on non-owning MPI ranks.
    std::optional<CappingSurface> cap_;
    /// @brief Reused global-column mesh state for cap gather (allocated on ranks that unpack; includes Yo/Yn rows).
    CapGlobalMeshState cap_global_mesh_state_;

    /// @brief Default CmMod used for cap gather/bcast (master rank id and MPI datatypes; stateless for this use).
    static const CmMod cm_mod_;

    /// Fill \ref cap_global_mesh_state_. Serial: all ranks. Parallel: root only (slaves skip buffer allocation).
    void gather_global_mesh_state(ComMod& com_mod, const CmMod& cm_mod, const SolutionStates& solutions, bool gather_Y);


public:
    /// @brief Default constructor - creates an uninitialized object
    CoupledBoundaryCondition() = default;

    /// @brief Destructor
    ~CoupledBoundaryCondition() = default;
    
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
    /// @param phys Equation physics for this boundary (struct, fluid, FSI, etc.)
    /// @param follower_pressure_load Follower pressure load flag (struct/ustruct); false for fluid-like physics
    CoupledBoundaryCondition(consts::BoundaryConditionType bc_type, const faceType& face, const std::string& face_name,
                          const std::string& block_name, consts::EquationType phys, bool follower_pressure_load);

    /// @brief Construct and optionally point to a cap face VTP file
    /// @param bc_type The 3D boundary condition type (must be bType_Dir or bType_Neu)
    /// @param face Face associated with this BC
    /// @param face_name Face name from the mesh
    /// @param block_name Block name in svZeroDSolver configuration
    /// @param cap_face_vtp_file Path to the cap face VTP file
    /// @param phys Equation physics for this boundary (struct, fluid, FSI, etc.)
    /// @param follower_pressure_load Follower pressure load flag (struct/ustruct); false for fluid-like physics
    CoupledBoundaryCondition(consts::BoundaryConditionType bc_type, const faceType& face, const std::string& face_name,
                          const std::string& block_name, const std::string& cap_face_vtp_file,
                          consts::EquationType phys, bool follower_pressure_load);

    /// @brief Get the 3D BC type for this Coupled boundary condition.
    consts::BoundaryConditionType get_bc_type() const { return bc_type_; }

    // =========================================================================
    // svZeroD block configuration
    // =========================================================================
    
    /// @brief Get the svZeroD block name
    /// @return Block name
    const std::string& get_block_name() const;
    
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

    /// @brief Set follower load flag and mechanical configs used for flowrate integration (also run from the face constructors).
    void set_flowrate_mechanical_configurations(consts::EquationType phys, bool follower_pressure_load);

    /// @brief Compute flowrates at the boundary face at old and new timesteps
    /// @param com_mod ComMod reference containing simulation data
    /// @param cm_mod CmMod reference for communication
    void compute_flowrates(ComMod& com_mod, const CmMod& cm_mod, const SolutionStates& solutions);

    /// @brief Initialize cap quadrature on the master (call from \c baf_ini after partition).
    void initialize_cap(ComMod& com_mod);

    /// @brief Gather mesh geometry, compute cap \a valM on master, and copy to FSILS face (all ranks enter MPI gather).
    void copy_cap_surface_to_linear_solver_face(ComMod& com_mod, fsi_linear_solver::FSILS_faceType& lhs_face,
                                                consts::MechanicalConfigurationType cfg,
                                                const SolutionStates& solutions);

    /// @brief Extra volumetric flux through the cap (old/new timestep); {0,0} if no cap; MPI-safe on all ranks.
    std::pair<double, double> calculate_cap_contribution(ComMod& com_mod, const CmMod& cm_mod,
                                                         const SolutionStates& solutions,
                                                         consts::MechanicalConfigurationType cfg_o,
                                                         consts::MechanicalConfigurationType cfg_n);
    
    /// @brief Compute average pressures at the boundary face at old and new timesteps (for Dirichlet BCs)
    /// @param com_mod ComMod reference containing simulation data
    /// @param cm_mod CmMod reference for communication
    void compute_pressures(ComMod& com_mod, const CmMod& cm_mod, const SolutionStates& solutions);

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

    /// @brief Load the cap face VTP file and associate it with this boundary condition
    /// @param vtp_file_path Path to the cap face VTP file
    void load_cap_face_vtp(const std::string& vtp_file_path);

    /// @brief Check if this BC has a cap (broadcast in distribute so all ranks agree).
    bool has_cap() const { return has_cap_; }

    /// @brief True if this rank stores the cap mesh / quadrature in \ref cap_.
    bool owns_cap() const { return owns_cap_; }

    /// @brief Master reads Neumann pressure, one scalar \c MPI_Bcast, all ranks set pressure (svZeroD sync).
    void bcast_coupled_neumann_pressure(const CmMod& cm_mod, cmType& cm);

};

#endif // COUPLED_BOUNDARY_CONDITION_H
