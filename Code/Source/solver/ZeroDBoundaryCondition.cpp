// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "ComMod.h"
#include "SimulationLogger.h"
#include "ZeroDBoundaryCondition.h"

#define n_debug_zerod_bc

ZeroDBoundaryCondition::ZeroDBoundaryCondition(const faceType& face, SimulationLogger& logger)
    : face_(&face)
    , logger_(&logger)
{
}

ZeroDBoundaryCondition::ZeroDBoundaryCondition(const std::string& cap_face_vtp_file, const faceType& face, SimulationLogger& logger)
    : cap_face_vtp_file_(cap_face_vtp_file)
    , face_(&face)
    , logger_(&logger)
{
    std::cout << "ZeroDBoundaryCondition: Cap face VTP file set to '" << cap_face_vtp_file_ << "'" << std::endl;
}

void ZeroDBoundaryCondition::load_cap_face_vtp(const std::string& vtp_file_path)
{
    cap_face_vtp_file_ = vtp_file_path;
    std::cout << "ZeroDBoundaryCondition: No cap '" << "'" << std::endl;
    // TODO: Load VTP file and extract geometric data
    // This would typically use VtkVtpData to read the file
}

// void ZeroDBoundaryCondition::compute_face_quantities()
// {
//     // TODO: Compute geometric/physical quantities on the cap face
//     // Typical quantities include:
//     //  - face area or lumped areas per node
//     //  - outward normal integration factors
//     //  - centroid(s)
//     //  - mappings between global/local node IDs and VTP arrays
// }

// void ZeroDBoundaryCondition::set_solver_interface(LPNSolverInterface* interface) noexcept
// {
//     solver_interface_ = interface;
// }

// void ZeroDBoundaryCondition::send_to_zerod(double time)
// {
//     // TODO: Send 3D-side boundary data (e.g., flow) to the 0D solver
//     // This would communicate with solver_interface_
//     (void)time; // Suppress unused parameter warning
// }

// void ZeroDBoundaryCondition::receive_from_zerod(double time)
// {
//     // TODO: Receive 0D-side boundary data (e.g., pressure) from the 0D solver
//     // This would communicate with solver_interface_
//     (void)time; // Suppress unused parameter warning
// }

// void ZeroDBoundaryCondition::exchange_with_zerod(double time)
// {
//     send_to_zerod(time);
//     receive_from_zerod(time);
// }

