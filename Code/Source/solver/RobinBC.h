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
#include "CmMod.h"

// Forward declarations
class faceType;


/// @brief Class to handle Robin boundary condition with potentially spatially variable stiffness and damping
/// 
/// This class can read stiffness and damping arrays from a VTP file, distribute them to all processes, and provides
/// efficient access to per-node values during boundary condition assembly.
class RobinBC {
public:
    /// @brief Default constructor deleted - RobinBC must be constructed with a face
    RobinBC() = delete;

    /// @brief Constructor - reads data from VTP file
    /// @param vtp_file_path Path to VTP file containing Stiffness and Damping point arrays
    /// @param face Face associated with the Robin BC
    /// @throws std::runtime_error if file cannot be read or arrays are missing
    RobinBC(const std::string& vtp_file_path, const faceType& face);
    
    /// @brief Constructor for uniform values
    /// @param uniform_stiffness Uniform stiffness value for all nodes
    /// @param uniform_damping Uniform damping value for all nodes
    /// @param face Face associated with the Robin BC
    RobinBC(double uniform_stiffness, double uniform_damping, const faceType& face);
    
    /// @brief Destructor
    ~RobinBC() = default;
    
    /// @brief Get stiffness value for a specific node
    /// @param node_id Node index on the face
    /// @return Stiffness value for the node
    double get_stiffness(int node_id) const;
    
    /// @brief Get damping value for a specific node
    /// @param node_id Node index on the face
    /// @return Damping value for the node
    double get_damping(int node_id) const;
    
    /// @brief Get global number of nodes
    /// @return Global number of nodes on the face
    int get_global_num_nodes() const { return global_num_nodes_; }

    /// @brief Get local number of nodes
    /// @return Local number of nodes on the face on this processor
    int get_local_num_nodes() const { return local_num_nodes_; }

    /// @brief Get local array index for a global node ID
    /// @param global_node_id The global node ID defined on the face
    /// @return Local array index for stiffness_array_ and damping_array_
    /// @throws std::runtime_error if global_node_id is not found in the map
    int get_local_index(int global_node_id) const;
    
    /// @brief Check if data is loaded from VTP file
    /// @return true if loaded from VTP, false if using uniform values
    bool is_from_vtp() const { return spatially_variable; }
    
    /// @brief Get the VTP file path (empty if using uniform values)
    /// @return VTP file path
    const std::string& get_vtp_path() const { return vtp_file_path_; }

    /// @brief Distribute Robin BC data from the master process to the slave processes. This is called by the distribute function. Follows distribution of face.
    /// @param cm_mod Reference to CmMod object for MPI communication
    /// @param cm Reference to cmType object for MPI communication
    void distribute(const CmMod& cm_mod, const cmType& cm);

private:
    /// @brief Initialize stiffness and damping arrays from VTP file. This is called by the constructor. Should only be called by the  master process.
    /// @param vtp_file_path Path to VTP file
    /// @param face Face associated with the Robin BC
    void initialize_from_vtp(const std::string& vtp_file_path);
    
    /// @brief Initialize uniform values. This is called by the constructor.
    /// @param stiffness Uniform stiffness value
    /// @param damping Uniform damping value
    /// @param num_nodes Number of nodes
    void initialize_uniform(double stiffness, double damping, int num_nodes);

    /// @brief Build map from mesh global node IDs to indices in the global arrays
    /// @param vtp_data VTP data object containing point coordinates
    /// @return Map from mesh global node IDs to indices in global arrays
    /// @throws std::runtime_error if a matching point cannot be found for any face node
    std::map<int, int> build_global_node_map(const VtkVtpData& vtp_data) const;
    
    const faceType& face_;             ///< Face associated with the Robin BC (reference cannot be null)
    int global_num_nodes_;             ///< Global number of nodes on the face
    int local_num_nodes_;              ///< Local number of nodes on this processor
    Vector<double> local_stiffness_;   ///< Local stiffness values for each node on this processor
    Vector<double> local_damping_;     ///< Local damping values for each node on this processor
    Vector<double> global_stiffness_;  ///< Global stiffness values (only populated on master)
    Vector<double> global_damping_;    ///< Global damping values (only populated on master)
    bool spatially_variable;                    ///< Flag indicating if data is from VTP file
    std::string vtp_file_path_;        ///< Path to VTP file (empty if uniform)
    std::map<int, int> global_node_map_; ///< Maps global node IDs to VTP point indices
    VtkVtpData vtp_data_;                 ///< VTP data object
};

#endif // ROBIN_BC_DATA_H
