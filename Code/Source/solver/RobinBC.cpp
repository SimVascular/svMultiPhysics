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

#include "RobinBCData.h"
#include "ComMod.h"
#include "DebugMsg.h"
#include <stdexcept>
#include <iostream>
#include <cstdio>

RobinBC::RobinBC(const std::string& vtp_file_path, const faceType& face)
{
    #define debug_robin_bc_data
    #ifdef debug_robin_bc_data
    DebugMsg dmsg(__func__, 0);
    dmsg << "Constructor with VTP file path" << std::endl;
    dmsg << "Loading Robin BC from VTP file: " << vtp_file_path << std::endl;
    #endif

    global_num_nodes_ = 0;
    local_num_nodes_ = 0;
    from_vtp_ = true;
    vtp_file_path_ = vtp_file_path;
    global_to_local_map_.clear();
    vtp_point_map_.clear();
    initialize_from_vtp(vtp_file_path, face);
}

RobinBC::RobinBC(double uniform_stiffness, double uniform_damping, const faceType& face)
{
    #ifdef debug_robin_bc_data
    DebugMsg dmsg(__func__, 0);
    dmsg << "Constructor with uniform values" << std::endl;
    dmsg << "Initializing uniform Robin BC values" << std::endl;
    #endif

    if (uniform_stiffness < 0.0) {
        throw std::runtime_error("Uniform stiffness value cannot be negative.");
    }
    if (uniform_damping < 0.0) {
        throw std::runtime_error("Uniform damping value cannot be negative.");
    }

    global_num_nodes_ = face.nNo;
    local_num_nodes_ = 0;
    from_vtp_ = false;
    vtp_file_path_ = "";
    global_to_local_map_.clear();
    vtp_point_map_.clear();

    global_stiffness_.resize(1);
    global_damping_.resize(1);
    global_stiffness_(0) = uniform_stiffness;
    global_damping_(0) = uniform_damping;

}

void RobinBC::initialize_from_vtp(const std::string& vtp_file_path, const faceType& face)
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
    
    // Check for required point arrays
    if (!vtp_data.has_point_data("Stiffness")) {
        throw std::runtime_error("VTP file '" + vtp_file_path + "' does not contain 'Stiffness' point array.");
    }
    
    if (!vtp_data.has_point_data("Damping")) {
        throw std::runtime_error("VTP file '" + vtp_file_path + "' does not contain 'Damping' point array.");
    }
    
    // Load global stiffness and damping arrays
    auto stiffness_data = vtp_data.get_point_data("Stiffness");
    auto damping_data = vtp_data.get_point_data("Damping");
    
    // Store the total number of points in the VTP file
    global_num_nodes_ = vtp_data.num_points();
    if (global_num_nodes_ != face.nNo) {
        throw std::runtime_error("Number of nodes in VTP file '" + vtp_file_path + "' does not match number of nodes on face '" + face.name + "'.");
    }
    
    if (stiffness_data.nrows() != global_num_nodes_ || stiffness_data.ncols() != 1) {
        throw std::runtime_error("'Stiffness' array in VTP file '" + vtp_file_path + 
                                "' has incorrect dimensions. Expected " + std::to_string(global_num_nodes_) + 
                                " x 1, got " + std::to_string(stiffness_data.nrows()) + " x " + 
                                std::to_string(stiffness_data.ncols()) + ".");
    }
    
    if (damping_data.nrows() != global_num_nodes_ || damping_data.ncols() != 1) {
        throw std::runtime_error("'Damping' array in VTP file '" + vtp_file_path + 
                                "' has incorrect dimensions. Expected " + std::to_string(global_num_nodes_) + 
                                " x 1, got " + std::to_string(damping_data.nrows()) + " x " + 
                                std::to_string(damping_data.ncols()) + ".");
    }
    
    // Store global data (master process only)
    global_stiffness_.resize(global_num_nodes_);
    global_damping_.resize(global_num_nodes_);
    
    // Copy VTP data to global arrays
    for (int i = 0; i < global_num_nodes_; i++) {
        global_stiffness_(i) = stiffness_data(i, 0);
        global_damping_(i) = damping_data(i, 0);
        
        // Validate that values are non-negative
        if (global_stiffness_(i) < 0.0) {
            throw std::runtime_error("Negative stiffness value at node " + std::to_string(i) + 
                                    " in VTP file '" + vtp_file_path + "'.");
        }
        if (global_damping_(i) < 0.0) {
            throw std::runtime_error("Negative damping value at node " + std::to_string(i) + 
                                    " in VTP file '" + vtp_file_path + "'.");
        }
    }
    
    // Build point map from VTP to face nodes for this processor's portion
    vtp_point_map_ = build_point_map(vtp_data, face);
    
    
    #ifdef debug_robin_bc_data
    dmsg << "Loaded Robin BC data for " << global_num_nodes_ << " total nodes" << std::endl;
    #endif
}

std::map<int, int> RobinBC::build_point_map(const VtkVtpData& vtp_data, const faceType& face) const
{
    #ifdef debug_robin_bc_data
    DebugMsg dmsg(__func__, 0);
    dmsg << "Building point map" << std::endl;
    dmsg << "VTP file has " << vtp_data.num_points() << " points" << std::endl;
    dmsg << "Face '" << face.name << "' has " << face.nNo << " nodes on this processor" << std::endl;
    dmsg << "Face '" << face.name << "' has global node IDs: " << face.gN << std::endl;
    dmsg << "face.x '" << face.name << "' has shape (" << face.x.nrows() << ", " << face.x.ncols() << ")" << std::endl;
    #endif

    // Get points from VTP file, (3,N) array
    Array<double> vtp_points = vtp_data.get_points();
    
    // Build map from VTP point index to face node index by matching coordinates
    const double tolerance = 1e-12;  // Small tolerance for floating point comparison
    std::map<int, int> vtp_to_face_map;  // Maps VTP point index to face node index

    #ifdef debug_robin_bc_data
    dmsg << "Building point map" << std::endl;
    #endif
    
    for (int i = 0; i < face.nNo; i++) {
        #ifdef debug_robin_bc_data
        dmsg << "Processing face node " << i << std::endl;
        #endif
        double face_x = face.x(0,i);
        double face_y = face.x(1,i);
        double face_z = face.x(2,i);
        #ifdef debug_robin_bc_data
        dmsg << "Face node " << i << " has coordinates (" << face_x << ", " << face_y << ", " << face_z << ")" << std::endl;
        #endif
        bool found_match = false;
        for (int j = 0; j < vtp_data.num_points(); j++) {
            double vtp_x = vtp_points(0,j);
            double vtp_y = vtp_points(1,j);
            double vtp_z = vtp_points(2,j);
            double dx = std::abs(vtp_x - face_x);
            double dy = std::abs(vtp_y - face_y);
            double dz = std::abs(vtp_z - face_z);
            
            if (dx <= tolerance && dy <= tolerance && dz <= tolerance) {
                vtp_to_face_map[j] = i;
                found_match = true;
                #ifdef debug_robin_bc_data
                dmsg << "Found matching VTP point " << j << " for face node " << i << std::endl;
                #endif
                break;
            }
        }

        if (!found_match) {
            throw std::runtime_error("Could not find matching point in VTP file for node " + 
                                   std::to_string(i) + " in face '" + face.name + "' at position (" +
                                   std::to_string(face_x) + ", " + std::to_string(face_y) + ", " +
                                   std::to_string(face_z) + ")");
        }
    }

    #ifdef debug_robin_bc_data
    dmsg << "Found matching VTP points for all " << face.nNo << " face nodes" << std::endl;
    #endif

    return vtp_to_face_map;
}

double RobinBC::get_stiffness(int node_id) const
{
    if (node_id < 0 || node_id >= global_num_nodes_) {
        throw std::runtime_error("Node ID " + std::to_string(node_id) + 
                                " is out of range [0, " + std::to_string(global_num_nodes_ - 1) + "].");
    }
    return local_stiffness_(node_id);
}

double RobinBC::get_damping(int node_id) const
{
    if (node_id < 0 || node_id >= global_num_nodes_) {
        throw std::runtime_error("Node ID " + std::to_string(node_id) + 
                                " is out of range [0, " + std::to_string(global_num_nodes_ - 1) + "].");
    }
    return local_damping_(node_id);
}

int RobinBC::get_local_index(int global_node_id) const
{
    auto it = global_to_local_map_.find(global_node_id);
    if (it == global_to_local_map_.end()) {
        throw std::runtime_error("Global node ID " + std::to_string(global_node_id) + 
                               " not found in global-to-local map.");
    }
    return it->second;
}

void RobinBC::distribute(const CmMod& cm_mod, const cmType& cm, const faceType& face, const Vector<int>& gmtl)
{
    #ifdef debug_robin_bc_data
    DebugMsg dmsg(__func__, cm.idcm());
    dmsg << "Distributing Robin BC data" << std::endl;
    #endif

    bool is_slave = cm.slv(cm_mod);

    // First communicate whether data is from VTP
    #ifdef debug_robin_bc_data
    dmsg << "Distributing whether data is from VTP" << std::endl;
    #endif
    cm.bcast(cm_mod, &from_vtp_);

    // Communicate VTP file path if data is from VTP
    if (from_vtp_) {
        #ifdef debug_robin_bc_data
        dmsg << "Distributing VTP file path" << std::endl;
        #endif
        cm.bcast(cm_mod, vtp_file_path_);
    }

    // Communicate total number of nodes in VTP file
    #ifdef debug_robin_bc_data
    dmsg << "Distributing total number of nodes" << std::endl;
    #endif
    cm.bcast(cm_mod, &global_num_nodes_);

    if (from_vtp_) {
        // Broadcast global arrays from master to all processes
        if (is_slave) {
            global_stiffness_.resize(global_num_nodes_);
            global_damping_.resize(global_num_nodes_);
        }

        #ifdef debug_robin_bc_data
        dmsg << "Broadcasting global stiffness and damping arrays" << std::endl;
        #endif
        cm.bcast(cm_mod, global_stiffness_);
        cm.bcast(cm_mod, global_damping_);

        // Each process needs to build its own point map based on its portion of the face
        #ifdef debug_robin_bc_data
        dmsg << "Reading VTP file on each process to build point map" << std::endl;
        #endif
        VtkVtpData vtp_data(vtp_file_path_, true);
        vtp_point_map_ = build_point_map(vtp_data, face);

        // Initialize local arrays based on this processor's portion of the face
        #ifdef debug_robin_bc_data
        dmsg << "Initializing local arrays for this processor's portion" << std::endl;
        #endif
        local_stiffness_.resize(face.nNo);
        local_damping_.resize(face.nNo);
        global_to_local_map_.clear();

        // Map global data to local arrays for this processor's portion
        for (const auto& pair : vtp_point_map_) {
            int vtp_idx = pair.first;   // Index in VTP file
            int face_idx = pair.second;  // Local index in face

            local_stiffness_(face_idx) = global_stiffness_(vtp_idx);
            local_damping_(face_idx) = global_damping_(vtp_idx);

            // Store mapping from global node ID to local index
            global_to_local_map_[face.gN(face_idx)] = face_idx;
        }

        // Clear global arrays on slave processes to save memory
        if (is_slave) {
            global_stiffness_.clear();
            global_damping_.clear();
        }
    } else {
        // For uniform values, just broadcast the first value from master's arrays
        double uniform_stiffness = 0.0;
        double uniform_damping = 0.0;
        
        if (!is_slave) {
            uniform_stiffness = local_stiffness_(0);
            uniform_damping = local_damping_(0);
        }
        
        cm.bcast(cm_mod, &uniform_stiffness);
        cm.bcast(cm_mod, &uniform_damping);

        // Initialize uniform values on each process
        initialize_uniform(uniform_stiffness, uniform_damping, face.nNo);
    }

    #ifdef debug_robin_bc_data
    dmsg << "Finished distributing Robin BC data" << std::endl;
    dmsg << "Number of nodes on this processor: " << face.nNo << std::endl;
    #endif
}
