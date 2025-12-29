// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef ZEROD_BOUNDARY_CONDITION_H
#define ZEROD_BOUNDARY_CONDITION_H

#include <string>
#include <vector>

// Forward declaration to avoid heavy includes
class LPNSolverInterface;
class faceType;
class ComMod;

/// @brief Skeleton for a 0D-coupled boundary condition on a cap face
///
/// This class provides an interface for:
///  - loading a cap face VTP,
///  - computing per-face quantities required for coupling, and
///  - communicating with a 0D (lumped parameter) solver interface.
///
/// Implementation is intentionally omitted here; this header only declares the API.
class ZeroDBoundaryCondition {
protected:

    /// @brief Data members for BC
    const faceType* face_ = nullptr;         ///< Face associated with the BC (not owned by ZeroDBoundaryCondition)
    std::string cap_face_vtp_file_;              ///< Path to VTP file (empty if no cap)
    const SimulationLogger* logger_ = nullptr;  ///< Logger for warnings/info (not owned by ZeroDBoundaryCondition)

public:
    /// @brief Default constructor - creates an uninitialized object
    ZeroDBoundaryCondition() = default;

    /// @brief Construct with a face association (no VTP data loaded)
    /// @param face Face associated with this BC
    /// @param logger Simulation logger used to write warnings
    ZeroDBoundaryCondition(const faceType& face, SimulationLogger& logger);

    /// @brief Construct and optionally point to a cap face VTP file
    /// @param cap_face_vtp_file Path to the cap face VTP file
    /// @param face Face associated with this BC
    /// @param logger Simulation logger used to write warnings
    ZeroDBoundaryCondition(const std::string& cap_face_vtp_file, const faceType& face, SimulationLogger& logger);

    /// @brief Load the cap face VTP file and associate it with this boundary condition
    /// @param vtp_file_path Path to the cap face VTP file
    /// @note This only declares the API; implementation is provided elsewhere.
    void load_cap_face_vtp(const std::string& vtp_file_path);

//     /// @brief Compute geometric/physical quantities on the cap face needed for coupling
//     /// Typical quantities may include:
//     ///  - face area or lumped areas per node,
//     ///  - outward normal integration factors,
//     ///  - centroid(s), and
//     ///  - mappings between global/local node IDs and VTP arrays.
//     /// @note This only declares the API; implementation is provided elsewhere.
//     void compute_face_quantities();

//     /// @brief Set the interface to the external 0D solver
//     /// @param interface Pointer to an already-initialized 0D solver interface
//     /// @note The class does not take ownership of the pointer.
//     void set_solver_interface(LPNSolverInterface* interface) noexcept;

//     /// @brief Send 3D-side boundary data (e.g., flow) to the 0D solver
//     /// @param time Current simulation time
//     /// @note This only declares the API; implementation is provided elsewhere.
//     void send_to_zerod(double time);

//     /// @brief Receive 0D-side boundary data (e.g., pressure) from the 0D solver
//     /// @param time Current simulation time
//     /// @note This only declares the API; implementation is provided elsewhere.
//     void receive_from_zerod(double time);

//     /// @brief Convenience method to perform a full exchange with the 0D solver
//     /// @param time Current simulation time
//     /// @note This only declares the API; implementation is provided elsewhere.
//     void exchange_with_zerod(double time);

// protected:
//     /// @brief No specific array validation required for 0D coupling template
//     void validate_array_value(const std::string&, double) const override {}

// private:
//     /// @brief Optional path to the cap face VTP used for geometric data
//     std::string cap_vtp_file_path_;

//     /// @brief External 0D solver interface (not owned)
//     LPNSolverInterface* solver_interface_ = nullptr;
};

#endif // ZEROD_BOUNDARY_CONDITION_H


