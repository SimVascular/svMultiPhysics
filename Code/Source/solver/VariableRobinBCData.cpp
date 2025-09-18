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

#include "VariableRobinBCData.h"
#include "ComMod.h"
#include "DebugMsg.h"
#include <stdexcept>
#include <iostream>
#include <cstdio>

VariableRobinBCData::VariableRobinBCData(const std::string& vtp_file_path, const faceType& face)
{
    #define debug_robin_bc_data
    #ifdef debug_robin_bc_data
    DebugMsg dmsg(__func__, 0);
    dmsg << "Constructor with VTP file path" << std::endl;
    dmsg << "Loading Robin BC from VTP file: " << vtp_file_path << std::endl;
    #endif

    num_nodes_ = face.nNo;
    from_vtp_ = true;
    vtp_file_path_ = vtp_file_path;
    load_from_vtp(vtp_file_path, face);
}

VariableRobinBCData::VariableRobinBCData(double uniform_stiffness, double uniform_damping, int num_nodes)
{
    #ifdef debug_robin_bc_data
    DebugMsg dmsg(__func__, 0);
    dmsg << "Constructor with uniform values" << std::endl;
    dmsg << "Initializing uniform Robin BC values" << std::endl;
    #endif

    num_nodes_ = num_nodes;
    from_vtp_ = false;
    vtp_file_path_ = "";
    initialize_uniform(uniform_stiffness, uniform_damping, num_nodes);
}

void VariableRobinBCData::load_from_vtp(const std::string& vtp_file_path, const faceType& face)
{
    #ifdef debug_robin_bc_data
    DebugMsg dmsg(__func__, 0);
    dmsg << "Loading Robin BC from VTP file: " << vtp_file_path << std::endl;
    #endif

    // Check if file exists
    if (FILE *file = fopen(vtp_file_path.c_str(), "r")) {
        fclose(file);
        #ifdef debug_robin_bc_data
        dmsg << "File exists and is readable" << std::endl;
        #endif
    } else {
        throw std::runtime_error("Robin BC VTP file '" + vtp_file_path + "' cannot be read.");
    }
    
    // Read the VTP file
    #ifdef debug_robin_bc_data
    dmsg << "About to create VtkVtpData object" << std::endl;
    #endif
    VtkVtpData vtp_data(vtp_file_path, true);
    #ifdef debug_robin_bc_data
    dmsg << "VtkVtpData object created successfully" << std::endl;
    #endif
    
    // Validate the VTP data matches the face geometry
    validate_vtp_data(vtp_data, face);
    
    // Check for required point arrays
    if (!vtp_data.has_point_data("Stiffness")) {
        throw std::runtime_error("VTP file '" + vtp_file_path + "' does not contain 'Stiffness' point array.");
    }
    
    if (!vtp_data.has_point_data("Damping")) {
        throw std::runtime_error("VTP file '" + vtp_file_path + "' does not contain 'Damping' point array.");
    }
    
    // Load stiffness array
    auto stiffness_data = vtp_data.get_point_data("Stiffness");
    if (stiffness_data.nrows() != num_nodes_ || stiffness_data.ncols() != 1) {
        throw std::runtime_error("'Stiffness' array in VTP file '" + vtp_file_path + 
                                "' has incorrect dimensions. Expected " + std::to_string(num_nodes_) + 
                                " x 1, got " + std::to_string(stiffness_data.nrows()) + " x " + 
                                std::to_string(stiffness_data.ncols()) + ".");
    }
    
    // Load damping array
    auto damping_data = vtp_data.get_point_data("Damping");
    if (damping_data.nrows() != num_nodes_ || damping_data.ncols() != 1) {
        throw std::runtime_error("'Damping' array in VTP file '" + vtp_file_path + 
                                "' has incorrect dimensions. Expected " + std::to_string(num_nodes_) + 
                                " x 1, got " + std::to_string(damping_data.nrows()) + " x " + 
                                std::to_string(damping_data.ncols()) + ".");
    }
    
    // Copy data to internal arrays and build global-to-local map
    stiffness_array_.resize(num_nodes_);
    damping_array_.resize(num_nodes_);
    global_to_local_map_.clear();
    
    for (int i = 0; i < num_nodes_; i++) {
        stiffness_array_(i) = stiffness_data(i, 0);
        damping_array_(i) = damping_data(i, 0);
        
        // Store mapping from global node ID to local index
        global_to_local_map_[face.gN(i)] = i;
        
        // Validate that values are non-negative
        if (stiffness_array_(i) < 0.0) {
            throw std::runtime_error("Negative stiffness value at node " + std::to_string(i) + 
                                    " in VTP file '" + vtp_file_path + "'.");
        }
        if (damping_array_(i) < 0.0) {
            throw std::runtime_error("Negative damping value at node " + std::to_string(i) + 
                                    " in VTP file '" + vtp_file_path + "'.");
        }
    }
}

void VariableRobinBCData::initialize_uniform(double stiffness, double damping, int num_nodes)
{
    if (stiffness < 0.0) {
        throw std::runtime_error("Uniform stiffness value cannot be negative.");
    }
    if (damping < 0.0) {
        throw std::runtime_error("Uniform damping value cannot be negative.");
    }
    if (num_nodes <= 0) {
        throw std::runtime_error("Number of nodes must be positive.");
    }
    
    stiffness_array_.resize(num_nodes);
    damping_array_.resize(num_nodes);
    global_to_local_map_.clear();
    
    for (int i = 0; i < num_nodes; i++) {
        stiffness_array_(i) = stiffness;
        damping_array_(i) = damping;
        // For uniform values, we assume global IDs are sequential starting from 0
        global_to_local_map_[i] = i;
    }
}

void VariableRobinBCData::validate_vtp_data(const VtkVtpData& vtp_data, const faceType& face)
{
    #ifdef debug_robin_bc_data
    DebugMsg dmsg(__func__, 0);
    dmsg << "Validating VTP data" << std::endl;
    #endif
    int vtp_num_points = vtp_data.num_points();
    if (vtp_num_points != face.nNo) {
        throw std::runtime_error("Number of points in VTP file (" + std::to_string(vtp_num_points) + 
                                ") does not match number of nodes on face '" + face.name + 
                                "' (" + std::to_string(face.nNo) + ").");
    }
    
    int vtp_num_elems = vtp_data.num_elems();
    if (vtp_num_elems != face.nEl) {
        throw std::runtime_error("Number of elements in VTP file (" + std::to_string(vtp_num_elems) + 
                                ") does not match number of elements on face '" + face.name + 
                                "' (" + std::to_string(face.nEl) + ").");
    }

    // Get points from VTP file, (3,N) array
    Array<double> vtp_points = vtp_data.get_points();
    
    // Compare each point's coordinates
    const double tolerance = 1e-12;  // Small tolerance for floating point comparison
    for (int i = 0; i < face.nNo; i++) {
        double face_x = face.x(0,i);
        double face_y = face.x(1,i);
        double face_z = face.x(2,i);
        double vtp_x = vtp_points(0,i);
        double vtp_y = vtp_points(1,i);
        double vtp_z = vtp_points(2,i);
        double dx = std::abs(vtp_x - face_x);
        double dy = std::abs(vtp_y - face_y);
        double dz = std::abs(vtp_z - face_z);
        
        if (dx > tolerance || dy > tolerance || dz > tolerance) {
            throw std::runtime_error("Node " + std::to_string(i) + " in VTP file is not colocated with " +
                                   "corresponding node in face '" + face.name + "'. " +
                                   "VTP: (" + std::to_string(vtp_x) + ", " + 
                                   std::to_string(vtp_y) + ", " + 
                                   std::to_string(vtp_z) + ") " +
                                   "Face: (" + std::to_string(face_x) + ", " + 
                                   std::to_string(face_y) + ", " + 
                                   std::to_string(face_z) + ")");
        }
    }
}

double VariableRobinBCData::get_stiffness(int node_id) const
{
    if (node_id < 0 || node_id >= num_nodes_) {
        throw std::runtime_error("Node ID " + std::to_string(node_id) + 
                                " is out of range [0, " + std::to_string(num_nodes_ - 1) + "].");
    }
    return stiffness_array_(node_id);
}

double VariableRobinBCData::get_damping(int node_id) const
{
    if (node_id < 0 || node_id >= num_nodes_) {
        throw std::runtime_error("Node ID " + std::to_string(node_id) + 
                                " is out of range [0, " + std::to_string(num_nodes_ - 1) + "].");
    }
    return damping_array_(node_id);
}

int VariableRobinBCData::get_local_index(int global_node_id) const
{
    auto it = global_to_local_map_.find(global_node_id);
    if (it == global_to_local_map_.end()) {
        throw std::runtime_error("Global node ID " + std::to_string(global_node_id) + 
                               " not found in global-to-local map.");
    }
    return it->second;
}
