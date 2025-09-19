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

#include "RobinBC.h"
#include "ComMod.h"
#include "DebugMsg.h"
#include <stdexcept>
#include <iostream>
#include <cstdio>

RobinBC::RobinBC(const std::string& vtp_file_path, const faceType& face)
    : face_(face), global_num_nodes_(0), local_num_nodes_(0), spatially_variable(true), vtp_file_path_(vtp_file_path)
{
    #define debug_robin_bc_data
    #ifdef debug_robin_bc_data
    DebugMsg dmsg(__func__, 0);
    dmsg << "Constructor with VTP file path" << std::endl;
    dmsg << "Loading Robin BC from VTP file: " << vtp_file_path << std::endl;
    #endif

    // Clear all containers before initialization
    global_node_map_.clear();
    local_stiffness_.clear();
    local_damping_.clear();
    global_stiffness_.clear();
    global_damping_.clear();

    // initialize_from_vtp will properly size and populate the vectors
    initialize_from_vtp(vtp_file_path, face);
}

RobinBC::RobinBC(double uniform_stiffness, double uniform_damping, const faceType& face)
    : face_(face), global_num_nodes_(face.nNo), local_num_nodes_(0), spatially_variable(false), vtp_file_path_("")
{
    // Validate input parameters first
    if (uniform_stiffness < 0.0) {
        throw std::runtime_error("Uniform stiffness value cannot be negative.");
    }
    if (uniform_damping < 0.0) {
        throw std::runtime_error("Uniform damping value cannot be negative.");
    }

    #ifdef debug_robin_bc_data
    DebugMsg dmsg(__func__, 0);
    dmsg << "Constructor with uniform values" << std::endl;
    dmsg << "Initializing uniform Robin BC values" << std::endl;
    #endif

    // Clear all maps
    global_node_map_.clear();

    // Clear global arrays (not used in uniform case)
    global_stiffness_.clear();
    global_damping_.clear();

    // Initialize localarrays with uniform values
    local_stiffness_.resize(1);
    local_damping_.resize(1);
    local_stiffness_(0) = uniform_stiffness;
    local_damping_(0) = uniform_damping;
}

void RobinBC::initialize_from_vtp(const std::string& vtp_file_path)
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
    vtp_data_ = VtkVtpData(vtp_file_path, true);
    #ifdef debug_robin_bc_data
    dmsg << "VtkVtpData object created successfully" << std::endl;
    #endif
    
    // Check for required point arrays
    if (!vtp_data_.has_point_data("Stiffness")) {
        throw std::runtime_error("VTP file '" + vtp_file_path + "' does not contain 'Stiffness' point array.");
    }
    
    if (!vtp_data_.has_point_data("Damping")) {
        throw std::runtime_error("VTP file '" + vtp_file_path + "' does not contain 'Damping' point array.");
    }
    
    // Load global stiffness and damping arrays
    auto stiffness_data = vtp_data_.get_point_data("Stiffness");
    auto damping_data = vtp_data_.get_point_data("Damping");
    
    // Store the total number of points in the VTP file
    global_num_nodes_ = vtp_data_.num_points();
    if (global_num_nodes_ != face_.nNo) {
        throw std::runtime_error("Number of nodes in VTP file '" + vtp_file_path + "' does not match number of nodes on face '" + face_.name + "'.");
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

    // Master process: build mapping from mesh global node IDs from face to VTP indices
    global_node_map_ = build_global_node_map(vtp_data_);
    
    
    #ifdef debug_robin_bc_data
    dmsg << "Loaded Robin BC data for " << global_num_nodes_ << " total nodes" << std::endl;
    #endif
}

double RobinBC::get_stiffness(int node_id) const
{
    if (spatially_variable) {
        if (node_id < 0 || node_id >= global_num_nodes_) {
            throw std::runtime_error("Node ID " + std::to_string(node_id) + 
                                        " is out of range [0, " + std::to_string(global_num_nodes_ - 1) + "].");
            }
            return local_stiffness_(node_id);
    } else {
        return local_stiffness_(0);
    }
}

double RobinBC::get_damping(int node_id) const
{
    if (spatially_variable) {
    if (node_id < 0 || node_id >= global_num_nodes_) {
        throw std::runtime_error("Node ID " + std::to_string(node_id) + 
                                " is out of range [0, " + std::to_string(global_num_nodes_ - 1) + "].");
    }
        return local_damping_(node_id);
    } else {
        return local_damping_(0);
    }
}

int RobinBC::get_local_index(int global_node_id) const
{
    auto it = global_node_map_.find(global_node_id);
    if (it == global_node_map_.end()) {
        throw std::runtime_error("Global node ID " + std::to_string(global_node_id) + 
                               " not found in global-to-local map.");
    }
    return it->second;
}

std::map<int, int> RobinBC::build_global_node_map(const VtkVtpData& vtp_data) const
{
    #ifdef debug_robin_bc_data
    DebugMsg dmsg(__func__, 0);
    dmsg << "Building global node map" << std::endl;
    #endif

    // Get points from VTP file
    Array<double> vtp_points = vtp_data.get_points();
    
    // Build map from mesh global node IDs to VTP point indices
    const double tolerance = 1e-12;
    std::map<int, int> global_node_to_vtp_idx;

    for (int i = 0; i < face_.nNo; i++) {
        double face_x = face_.x(0,i);
        double face_y = face_.x(1,i);
        double face_z = face_.x(2,i);
        
        bool found_match = false;
        for (int j = 0; j < vtp_data.num_points(); j++) {
            double vtp_x = vtp_points(0,j);
            double vtp_y = vtp_points(1,j);
            double vtp_z = vtp_points(2,j);
            
            if (std::abs(vtp_x - face_x) <= tolerance && 
                std::abs(vtp_y - face_y) <= tolerance && 
                std::abs(vtp_z - face_z) <= tolerance) {
                global_node_to_vtp_idx[face_.gN(i)] = j;
                found_match = true;
                break;
            }
        }

        if (!found_match) {
            throw std::runtime_error("Could not find matching point in VTP file for node " + 
                                   std::to_string(i) + " (global ID " + std::to_string(face_.gN(i)) + 
                                   ") in face '" + face_.name + "' at position (" +
                                   std::to_string(face_x) + ", " + std::to_string(face_y) + ", " +
                                   std::to_string(face_z) + ")");
        }
    }

    return global_node_to_vtp_idx;
}

void RobinBC::distribute(const CmMod& cm_mod, const cmType& cm)
{
    #ifdef debug_robin_bc_data
    DebugMsg dmsg(__func__, cm.idcm());
    dmsg << "Distributing Robin BC data" << std::endl;
    #endif

    bool is_slave = cm.slv(cm_mod);

    // First communicate whether data is spatially variable
    #ifdef debug_robin_bc_data
    dmsg << "Distributing whether data is spatially variable" << std::endl;
    #endif
    cm.bcast(cm_mod, &spatially_variable);

    // Communicate VTP file path if data is spatially variable
    if (spatially_variable) {
        #ifdef debug_robin_bc_data
        dmsg << "Distributing VTP file path" << std::endl;
        #endif
        cm.bcast(cm_mod, vtp_file_path_);
    }

    // Communicate global number of nodes
    #ifdef debug_robin_bc_data
    dmsg << "Distributing total number of nodes" << std::endl;
    #endif
    cm.bcast(cm_mod, &global_num_nodes_);

    if (spatially_variable) {
        // Initialize local storage
        local_num_nodes_ = face.nNo;
        local_stiffness_.resize(local_num_nodes_);
        local_damping_.resize(local_num_nodes_);

        // Each process collects its global node IDs
        Vector<int> local_global_ids(local_num_nodes_);
        for (int i = 0; i < local_num_nodes_; i++) {
            local_global_ids(i) = face_.gN(i);
        }

        // Gather number of nodes from each process to master
        Vector<int> proc_num_nodes(cm.np());
        cm.gather(cm_mod, &local_num_nodes_, 1, proc_num_nodes.data(), 1, 0);

        // Calculate displacements for gatherv/scatterv
        Vector<int> displs(cm.np());
        int total_num_nodes = 0;
        for (int i = 0; i < cm.np(); i++) {
            displs(i) = total_num_nodes;
            total_num_nodes += proc_num_nodes(i);
        }

        // Master process: gather all global IDs and prepare data to scatter back to all processes
        Vector<int> all_global_ids;
        Vector<double> all_stiffness, all_damping;
        if (!is_slave) {
            // Resize receive buffers based on total number of nodes
            all_global_ids.resize(total_num_nodes);
            all_stiffness.resize(total_num_nodes);
            all_damping.resize(total_num_nodes);

            // Gather all global IDs to master using gatherv
            cm.gatherv(cm_mod, local_global_ids, all_global_ids, proc_num_nodes, displs, 0);

            // Look up data for all nodes
            for (int i = 0; i < total_num_nodes; i++) {
                int vtp_idx = global_node_map_[all_global_ids(i)];
                all_stiffness(i) = global_stiffness_(vtp_idx);
                all_damping(i) = global_damping_(vtp_idx);
            }

            // Clear global arrays to save memory
            global_stiffness_.clear();
            global_damping_.clear();
        } else {
            // Slave processes: send global IDs to master
            cm.gatherv(cm_mod, local_global_ids, all_global_ids, proc_num_nodes, displs, 0);
        }

        // Scatter data back to all processes using scatterv
        cm.scatterv(cm_mod, all_stiffness, proc_num_nodes, displs, local_stiffness_, 0);
        cm.scatterv(cm_mod, all_damping, proc_num_nodes, displs, local_damping_, 0);

    } else {
        // For uniform values, just broadcast the single values
        double uniform_stiffness = 0.0;
        double uniform_damping = 0.0;
        
        if (!is_slave) {
            uniform_stiffness = local_stiffness_(0);
            uniform_damping = local_damping_(0);
        }
        
        cm.bcast(cm_mod, &uniform_stiffness);
        cm.bcast(cm_mod, &uniform_damping);

        // Initialize uniform values on each process
        local_num_nodes_ = face_.nNo;
        local_stiffness_.resize(1);
        local_damping_.resize(1);
        local_stiffness_(0) = uniform_stiffness;
        local_damping_(0) = uniform_damping;
    }

    #ifdef debug_robin_bc_data
    dmsg << "Finished distributing Robin BC data" << std::endl;
    dmsg << "Number of nodes on this processor: " << local_num_nodes_ << std::endl;
    #endif
}
