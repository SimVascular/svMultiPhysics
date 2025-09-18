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

#ifndef ROBIN_BC_DATA_H
#define ROBIN_BC_DATA_H

#include "Array.h"
#include "Vector.h"
#include "VtkData.h"
#include <string>
#include <memory>
#include <map>

// Forward declarations
class faceType;


/// @brief Class to handle Robin boundary condition data with per-node stiffness and damping
/// 
/// This class reads stiffness and damping arrays from a VTP file and provides
/// efficient access to per-node values during boundary condition assembly.
class VariableRobinBCData {
public:
    /// @brief Default constructor - initializes to empty state
    VariableRobinBCData() : num_nodes_(0), from_vtp_(false) {}

    /// @brief Constructor - reads data from VTP file
    /// @param vtp_file_path Path to VTP file containing Stiffness and Damping point arrays
    /// @param face Reference to the face for validation
    /// @throws std::runtime_error if file cannot be read or arrays are missing
    VariableRobinBCData(const std::string& vtp_file_path, const faceType& face);
    
    /// @brief Constructor for uniform values
    /// @param uniform_stiffness Uniform stiffness value for all nodes
    /// @param uniform_damping Uniform damping value for all nodes
    /// @param num_nodes Number of nodes on the face
    VariableRobinBCData(double uniform_stiffness, double uniform_damping, int num_nodes);
    
    /// @brief Destructor
    ~VariableRobinBCData() = default;
    
    /// @brief Get stiffness value for a specific node
    /// @param node_id Node index on the face
    /// @return Stiffness value for the node
    double get_stiffness(int node_id) const;
    
    /// @brief Get damping value for a specific node
    /// @param node_id Node index on the face
    /// @return Damping value for the node
    double get_damping(int node_id) const;
    
    /// @brief Get number of nodes
    /// @return Number of nodes on the face
    int get_num_nodes() const { return num_nodes_; }

    /// @brief Get local array index for a global node ID
    /// @param global_node_id The global node ID defined on the face
    /// @return Local array index for stiffness_array_ and damping_array_
    /// @throws std::runtime_error if global_node_id is not found in the map
    int get_local_index(int global_node_id) const;
    
    /// @brief Check if data is loaded from VTP file
    /// @return true if loaded from VTP, false if using uniform values
    bool is_from_vtp() const { return from_vtp_; }
    
    /// @brief Get the VTP file path (empty if using uniform values)
    /// @return VTP file path
    const std::string& get_vtp_path() const { return vtp_file_path_; }

private:
    /// @brief Load stiffness and damping arrays from VTP file
    /// @param vtp_file_path Path to VTP file
    /// @param face Reference to face for validation
    void load_from_vtp(const std::string& vtp_file_path, const faceType& face);
    
    /// @brief Initialize uniform values
    /// @param stiffness Uniform stiffness value
    /// @param damping Uniform damping value
    /// @param num_nodes Number of nodes
    void initialize_uniform(double stiffness, double damping, int num_nodes);
    
    /// @brief Validate that VTP data matches face geometry
    /// @param vtp_data VTP data object
    /// @param face Reference to face
    void validate_vtp_data(const VtkVtpData& vtp_data, const faceType& face);
    
    int num_nodes_;                    ///< Number of nodes on the face
    Vector<double> stiffness_array_;   ///< Stiffness values for each node
    Vector<double> damping_array_;     ///< Damping values for each node
    bool from_vtp_;                    ///< Flag indicating if data is from VTP file
    std::string vtp_file_path_;        ///< Path to VTP file (empty if uniform)
    std::map<int, int> global_to_local_map_; ///< Maps global node IDs to local array indices
};

#endif // ROBIN_BC_DATA_H
