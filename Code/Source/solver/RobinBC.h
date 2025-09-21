/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef ROBIN_BC_H
#define ROBIN_BC_H

#include "BC.h"
#include <string>
#include <map>
#include <vector>

/// @brief Class to handle Robin boundary condition with potentially spatially variable arrays
/// 
/// This class extends the generic BC class to handle Robin boundary conditions, which require
/// stiffness and damping arrays. While it supports any number of named arrays through its base class,
/// it provides specific validation and convenience methods for stiffness and damping values.
///
/// Example usage:
/// ```cpp
/// // Read multiple arrays from VTP file
/// std::vector<std::string> array_names = {"Stiffness", "Damping"};
/// RobinBC bc(vtp_file_path, array_names, face);
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
/// RobinBC bc(uniform_values, face);
/// ```

#define debug_robin_bc
class RobinBC : public BC {
public:
    /// @brief Default constructor - creates an empty RobinBC
    RobinBC() : BC() {}

    /// @brief Constructor - reads data from VTP file
    /// @param vtp_file_path Path to VTP file containing arrays
    /// @param array_names Names of arrays to read from VTP file
    /// @param face Face associated with the Robin BC
    /// @throws std::runtime_error if file cannot be read or arrays are missing
    RobinBC(const std::string& vtp_file_path, const std::vector<std::string>& array_names, const faceType& face)
        : BC(vtp_file_path, array_names, face) {}

    /// @brief Legacy constructor - reads stiffness and damping from VTP file
    /// @param vtp_file_path Path to VTP file containing Stiffness and Damping point arrays
    /// @param face Face associated with the Robin BC
    /// @throws std::runtime_error if file cannot be read or arrays are missing
    RobinBC(const std::string& vtp_file_path, const faceType& face) 
        : RobinBC(vtp_file_path, std::vector<std::string>{"Stiffness", "Damping"}, face) {}

    /// @brief Constructor for uniform values
    /// @param uniform_values Map of array names to uniform values
    /// @param face Face associated with the Robin BC
    RobinBC(const std::map<std::string, double>& uniform_values, const faceType& face)
        : BC(uniform_values, face) {}

    /// @brief Legacy constructor for uniform values
    /// @param uniform_stiffness Uniform stiffness value for all nodes
    /// @param uniform_damping Uniform damping value for all nodes
    /// @param face Face associated with the Robin BC
    RobinBC(double uniform_stiffness, double uniform_damping, const faceType& face) {
        std::map<std::string, double> uniform_values = {
            {"Stiffness", uniform_stiffness},
            {"Damping", uniform_damping}
        };
        BC::init_uniform(uniform_values, face);
    }

    /// @brief Initialize from VTP file
    /// @param vtp_file_path Path to VTP file containing Stiffness and Damping point arrays
    /// @param face Face associated with the BC
    /// @throws std::runtime_error if file cannot be read or arrays are missing
    void init_from_vtp(const std::string& vtp_file_path, const faceType& face);

    /// @brief Initialize with uniform values
    /// @param uniform_stiffness Uniform stiffness value for all nodes
    /// @param uniform_damping Uniform damping value for all nodes
    /// @param face Face associated with the BC
    void init_uniform(double uniform_stiffness, double uniform_damping, const faceType& face);

    /// @brief Get stiffness value for a specific node (convenience method)
    /// @param node_id Node index on the face
    /// @return Stiffness value for the node
    double get_stiffness(int node_id) const { return get_value("Stiffness", node_id); }
    
    /// @brief Get damping value for a specific node (convenience method)
    /// @param node_id Node index on the face
    /// @return Damping value for the node
    double get_damping(int node_id) const { return get_value("Damping", node_id); }

    /// @brief Assemble the Robin BC into the global residual vector and stiffness matrix
    /// Currently not implemented
    /// @return 0
    double assemble() const { return 0; }

protected:
    /// @brief Validate array values for Robin BC
    /// @param array_name Name of the array being validated
    /// @param value Value to validate
    /// @throws std::runtime_error if validation fails
    void validate_array_value(const std::string& array_name, double value) const override {
        if (value < 0.0) {
            throw std::runtime_error("Negative value for array '" + array_name + "'");
        }
    }
};

#endif // ROBIN_BC_H