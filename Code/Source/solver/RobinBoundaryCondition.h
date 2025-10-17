// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef ROBIN_BOUNDARY_CONDITION_H
#define ROBIN_BOUNDARY_CONDITION_H

#include "BoundaryCondition.h"
#include <string>
#include <map>
#include <vector>

/// @brief Class to handle Robin boundary condition with potentially spatially variable arrays
/// 
/// This class extends the generic BoundaryCondition class to handle Robin boundary conditions, which require
/// stiffness and damping arrays. While it supports any number of named arrays through its base class,
/// it provides specific validation and convenience methods for stiffness and damping values.
///
/// Example usage:
/// ```cpp
/// // Read multiple arrays from VTP file
/// std::vector<std::string> array_names = {"Stiffness", "Damping"};
/// RobinBoundaryCondition bc(vtp_file_path, array_names, face);
///
/// // Access values
/// double stiffness = bc.get_stiffness(node_id);  // Convenience method
/// double damping = bc.get_damping(node_id);      // Convenience method
///
/// // Create with uniform values
/// std::map<std::string, double> uniform_values = {
///     {"Stiffness", 1.0},
///     {"Damping", 0.5},
/// };
/// RobinBoundaryCondition bc(uniform_values, face);
/// ```

#define debug_robin_bc
class RobinBoundaryCondition : public BoundaryCondition {
public:
    /// @brief Default constructor - creates an empty RobinBoundaryCondition
    RobinBoundaryCondition() : BoundaryCondition() {}

    /// @brief Constructor - reads stiffness and damping from VTP file
    /// @param vtp_file_path Path to VTP file containing Stiffness and Damping point arrays
    /// @param normal_only Flag to apply only along normal direction
    /// @param face Face associated with the Robin BC
    /// @param logger Simulation logger used to write warnings
    /// @throws BoundaryConditionFileException if file cannot be read
    /// @throws BoundaryConditionVtpArrayException if arrays are missing
    /// @throws BoundaryConditionValidationException if values are invalid
    RobinBoundaryCondition(const std::string& vtp_file_path, bool normal_only, const faceType& face, SimulationLogger& logger)
        : BoundaryCondition(vtp_file_path, std::vector<std::string>{"Stiffness", "Damping"}, StringBoolMap{{"normal_direction_only", normal_only}}, face, logger) {}


    /// @brief Constructor for uniform values
    /// @param uniform_stiffness Uniform stiffness value for all nodes
    /// @param uniform_damping Uniform damping value for all nodes
    /// @param normal_only Flag to apply only along normal direction
    /// @param face Face associated with the Robin BC
    /// @param logger Simulation logger used to write warnings
    /// @throws BoundaryConditionValidationException if values are invalid
    RobinBoundaryCondition(double uniform_stiffness, double uniform_damping, bool normal_only, const faceType& face, SimulationLogger& logger);
 
    /// @brief Apply only along normal direction (getter)
    /// @return true if BC should be applied only along normal direction
    /// @throws BoundaryConditionFlagException if "normal_direction_only" flag not found
    bool normal_direction_only() const { return this->get_flag("normal_direction_only"); }
 
    /// @brief Get stiffness value for a specific node (convenience method)
    /// @param node_id Node index on the face
    /// @return Stiffness value for the node
    /// @throws BoundaryConditionArrayException if "Stiffness" array not found
    /// @throws BoundaryConditionNodeIdException if node_id is out of range
    double get_stiffness(int node_id) const {
        return get_value("Stiffness", node_id);
    }
    
    /// @brief Get damping value for a specific node (convenience method)
    /// @param node_id Node index on the face
    /// @return Damping value for the node
    /// @throws BoundaryConditionArrayException if "Damping" array not found
    /// @throws BoundaryConditionNodeIdException if node_id is out of range
    double get_damping(int node_id) const {
        return get_value("Damping", node_id);
    }

    /// @brief Assemble the Robin BC into the global residual vector and stiffness matrix
    /// Currently not implemented
    /// @return 0
    double assemble() const { return 0; }

protected:
    /// @brief Validate array values for Robin BC
    /// @param array_name Name of the array being validated
    /// @param value Value to validate
    /// @throws BoundaryConditionValidationException if validation fails
    void validate_array_value(const std::string& array_name, double value) const override {
        if (value < 0.0) {
            throw BoundaryConditionValidationException(array_name, value);
        }
    }
};

#endif // ROBIN_BOUNDARY_CONDITION_H
