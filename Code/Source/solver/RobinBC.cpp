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
#include "Vector.h"
#include <stdexcept>
#include <iostream>
#include <cstdio>
#include <vector>

#define debug_robin_bc_data

RobinBC::RobinBC(const std::string& vtp_file_path, const faceType& face)
{
    init_from_vtp(vtp_file_path, face);
}

RobinBC::RobinBC(double uniform_stiffness, double uniform_damping, const faceType& face)
{
    init_uniform(uniform_stiffness, uniform_damping, face);
}

RobinBC::RobinBC(const RobinBC& other)
    : face_(other.face_)
    , global_num_nodes_(other.global_num_nodes_)
    , local_num_nodes_(other.local_num_nodes_)
    , local_stiffness_(other.local_stiffness_)
    , local_damping_(other.local_damping_)
    , global_stiffness_(other.global_stiffness_)
    , global_damping_(other.global_damping_)
    , spatially_variable(other.spatially_variable)
    , vtp_file_path_(other.vtp_file_path_)
    , global_node_map_(other.global_node_map_)
    , defined_(other.defined_)
{
    if (other.vtp_data_) {
        vtp_data_ = std::make_unique<VtkVtpData>(*other.vtp_data_);
    }
}

RobinBC& RobinBC::operator=(const RobinBC& other)
{
    if (this != &other) {
        face_ = other.face_;
        global_num_nodes_ = other.global_num_nodes_;
        local_num_nodes_ = other.local_num_nodes_;
        local_stiffness_ = other.local_stiffness_;
        local_damping_ = other.local_damping_;
        global_stiffness_ = other.global_stiffness_;
        global_damping_ = other.global_damping_;
        spatially_variable = other.spatially_variable;
        vtp_file_path_ = other.vtp_file_path_;
        global_node_map_ = other.global_node_map_;
        defined_ = other.defined_;

        if (other.vtp_data_) {
            vtp_data_ = std::make_unique<VtkVtpData>(*other.vtp_data_);
        } else {
            vtp_data_.reset();
        }
    }
    return *this;
}

RobinBC::RobinBC(RobinBC&& other) noexcept
    : face_(other.face_)
    , global_num_nodes_(other.global_num_nodes_)
    , local_num_nodes_(other.local_num_nodes_)
    , local_stiffness_(std::move(other.local_stiffness_))
    , local_damping_(std::move(other.local_damping_))
    , global_stiffness_(std::move(other.global_stiffness_))
    , global_damping_(std::move(other.global_damping_))
    , spatially_variable(other.spatially_variable)
    , vtp_file_path_(std::move(other.vtp_file_path_))
    , global_node_map_(std::move(other.global_node_map_))
    , vtp_data_(std::move(other.vtp_data_))
    , defined_(other.defined_)
{
    other.face_ = nullptr;
    other.global_num_nodes_ = 0;
    other.local_num_nodes_ = 0;
    other.spatially_variable = false;
    other.defined_ = false;
}

RobinBC& RobinBC::operator=(RobinBC&& other) noexcept
{
    if (this != &other) {
        face_ = other.face_;
        global_num_nodes_ = other.global_num_nodes_;
        local_num_nodes_ = other.local_num_nodes_;
        local_stiffness_ = std::move(other.local_stiffness_);
        local_damping_ = std::move(other.local_damping_);
        global_stiffness_ = std::move(other.global_stiffness_);
        global_damping_ = std::move(other.global_damping_);
        spatially_variable = other.spatially_variable;
        vtp_file_path_ = std::move(other.vtp_file_path_);
        global_node_map_ = std::move(other.global_node_map_);
        vtp_data_ = std::move(other.vtp_data_);
        defined_ = other.defined_;

        other.face_ = nullptr;
        other.global_num_nodes_ = 0;
        other.local_num_nodes_ = 0;
        other.spatially_variable = false;
        other.defined_ = false;
    }
    return *this;
}

void RobinBC::init_uniform(double uniform_stiffness, double uniform_damping, const faceType& face)
{
    face_ = &face;
    global_num_nodes_ = face_->nNo;
    local_num_nodes_ = 0;
    spatially_variable = false;
    vtp_file_path_ = "";
    defined_ = true;

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

void RobinBC::init_from_vtp(const std::string& vtp_file_path, const faceType& face)
{
    // Note that this function is only called by the master process

    face_ = &face;
    global_num_nodes_ = face_->nNo;
    local_num_nodes_ = 0;
    spatially_variable = true;
    vtp_file_path_ = vtp_file_path;
    defined_ = true;

    // Clear all containers before initialization
    global_node_map_.clear();
    local_stiffness_.clear();
    local_damping_.clear();
    global_stiffness_.clear();
    global_damping_.clear();

    #ifdef debug_robin_bc_data
    DebugMsg dmsg(__func__, 0);
    dmsg << "Constructor with VTP file path" << std::endl;
    #endif

    // read_data_from_vtp_file will properly size and populate the vectors
    std::vector<std::string> array_names = {"Stiffness", "Damping"};
    auto vtp_data = read_data_from_vtp_file(vtp_file_path, array_names);

    // Store global data
    global_stiffness_.resize(global_num_nodes_);
    global_damping_.resize(global_num_nodes_);
    
    // Copy VTP data to global arrays
    for (int i = 0; i < global_num_nodes_; i++) {
        global_stiffness_(i) = vtp_data["Stiffness"](i, 0);
        global_damping_(i) = vtp_data["Damping"](i, 0);
        
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

    // In case we are running sequentially, we need to fill the local arrays 
    // and the global node map as well, because distribute is not called in sequential mode.
    local_stiffness_.resize(global_num_nodes_);
    local_damping_.resize(global_num_nodes_);
    for (int i = 0; i < global_num_nodes_; i++) {
        local_stiffness_(i) = global_stiffness_(i);
        local_damping_(i) = global_damping_(i);
    }

    // Initialize global node map, in case we are running sequentially
    global_node_map_.clear();
    for (int i = 0; i < global_num_nodes_; i++) {
        global_node_map_[face_->gN(i)] = i;
    }
}

RobinBC::StringArrayMap RobinBC::read_data_from_vtp_file(const std::string& vtp_file_path, const std::vector<std::string>& array_names)
{
    #ifdef debug_robin_bc_data
    DebugMsg dmsg(__func__, 0);
    dmsg << "Loading data from VTP file: " << vtp_file_path << std::endl;
    #endif

    // Check if file exists
    if (FILE *file = fopen(vtp_file_path.c_str(), "r")) {
        fclose(file);
        #ifdef debug_robin_bc_data
        dmsg << "File exists and is readable" << std::endl;
        #endif
    } else {
        throw std::runtime_error("VTP file '" + vtp_file_path + "' cannot be read.");
    }
    
    // Read the VTP file
    vtp_data_ = std::make_unique<VtkVtpData>(vtp_file_path, true);

    #ifdef debug_robin_bc_data
    dmsg << "VtkVtpData object created successfully" << std::endl;
    #endif
    
    // Check if the number of nodes in the VTP file matches the number of nodes on the face
    if (global_num_nodes_ != face_->nNo) {
        throw std::runtime_error("Number of nodes in VTP file '" + vtp_file_path + "' does not match number of nodes on face '" + face_->name + "'.");
    }

    // Create map to store results
    StringArrayMap result;

    // Check for and load each requested array
    for (const auto& array_name : array_names) {
        // Check if array exists
        if (!vtp_data_->has_point_data(array_name)) {
            throw std::runtime_error("VTP file '" + vtp_file_path + "' does not contain '" + array_name + "' point array.");
        }

        // Load array data
        auto array_data = vtp_data_->get_point_data(array_name);

        // Validate dimensions
        if (array_data.nrows() != global_num_nodes_ || array_data.ncols() != 1) {
            throw std::runtime_error("'" + array_name + "' array in VTP file '" + vtp_file_path + 
                                   "' has incorrect dimensions. Expected " + std::to_string(global_num_nodes_) + 
                                   " x 1, got " + std::to_string(array_data.nrows()) + " x " + 
                                   std::to_string(array_data.ncols()) + ".");
        }

        // Store array in result map
        result[array_name] = array_data;

        #ifdef debug_robin_bc_data
        dmsg << "Successfully loaded " << array_name << " data" << std::endl;
        dmsg << array_name << " data size: " << array_data.nrows() << " x " << array_data.ncols() << std::endl;
        #endif
    }

    #ifdef debug_robin_bc_data
    dmsg << "Finished loading VTP data" << std::endl;
    #endif

    return result;
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
    if (spatially_variable) {
        auto it = global_node_map_.find(global_node_id);
        if (it == global_node_map_.end()) {
            throw std::runtime_error("Global node ID " + std::to_string(global_node_id) + 
                                   " not found in global-to-local map.");
        }
        return it->second;
    } else {
        return 0;
    }
}

int RobinBC::find_vtp_point_index(double x, double y, double z,
                                const Array<double>& vtp_points) const
{
    const double tolerance = 1e-12;
    const int num_points = vtp_points.ncols();
    
    // Simple linear search through all points
    for (int i = 0; i < num_points; i++) {
        // Compare coordinates with tolerance
        if (std::abs(vtp_points(0,i) - x) <= tolerance &&
            std::abs(vtp_points(1,i) - y) <= tolerance &&
            std::abs(vtp_points(2,i) - z) <= tolerance) {

            #define n_debug_robin_bc_data_find_vtp_point_index
            #ifdef debug_robin_bc_data_find_vtp_point_index
            DebugMsg dmsg(__func__, 0);
            dmsg << "Found VTP point index for node at position (" << x << ", " << y << ", " << z << ")" << std::endl;
            dmsg << "VTP point index: " << i << std::endl;
            #endif

            return i;  // Found match
        }
    }
    
    throw std::runtime_error("Could not find matching point in VTP file for node at position (" +
                          std::to_string(x) + ", " + std::to_string(y) + ", " +
                          std::to_string(z) + ")");
}

void RobinBC::distribute(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, const faceType& face)
{
    #define n_debug_robin_bc_data_distribute
    #ifdef debug_robin_bc_data_distribute
    DebugMsg dmsg(__func__, cm.idcm());
    dmsg << "Distributing Robin BC data" << std::endl;
    #endif

    // Before this point, only the master process had a face_ pointer, which contained
    // the entire (global) face. Within the distribute::distribute() function,
    // the global face was partitioned and distributed among all processes. Each
    // local face contains only a portion of the global face, corresponding to 
    // the portion of the volume mesh owned by this process. Here, we update the
    // face_ pointer to the local face.
    face_ = &face;

    // Number of nodes on the face on this processor
    local_num_nodes_ = face_->nNo;

    bool is_slave = cm.slv(cm_mod);

    // First communicate whether data is spatially variable
    cm.bcast(cm_mod, &spatially_variable);

    // Communicate VTP file path if data is spatially variable 
    // (not necessary, but we do it for consistency)
    if (spatially_variable) {
        cm.bcast(cm_mod, vtp_file_path_);
    }

    // Communicate global number of nodes
    // (not necessary, but we do it for consistency)
    cm.bcast(cm_mod, &global_num_nodes_);

    // Communicate stiffness and damping values needed by each process
    if (spatially_variable) {

        // Setup
        if (face_ == nullptr) {
            throw std::runtime_error("face_ is nullptr during distribute");
        }
        local_stiffness_.resize(local_num_nodes_);
        local_damping_.resize(local_num_nodes_);

        // Each processor collects the global node IDs and nodal positions of its 
        // associated face portion
        Vector<int> local_global_ids = face_->gN;
        Array<double> local_positions(3, local_num_nodes_);
        for (int i = 0; i < local_num_nodes_; i++) {
            local_positions.set_col(i, com_mod.x.col(face_->gN(i)));
        }

        #ifdef debug_robin_bc_data_distribute
        dmsg << "Number of face nodes on this processor: " << local_num_nodes_ << std::endl;
        dmsg << "Local global IDs: " << local_global_ids << std::endl;
        dmsg << "Local positions: " << local_positions << std::endl;
        #endif

        // Gather number of face nodes from each processor to master
        Vector<int> proc_num_nodes(cm.np());
        cm.gather(cm_mod, &local_num_nodes_, 1, proc_num_nodes.data(), 1, 0);

        // Calculate displacements for gatherv/scatterv and compute total number of nodes
        // total_num_nodes is the total number of face nodes across all processors.
        Vector<int> displs(cm.np());
        int total_num_nodes = 0;
        for (int i = 0; i < cm.np(); i++) {
            displs(i) = total_num_nodes;
            total_num_nodes += proc_num_nodes(i);
        }

        // Master process: gather the nodal positions of face nodes from all processors,
        // get the corresponding stiffness and damping values by matching the positions to the VTP points,
        // and scatter the data back to all processors.
        Array<double> all_positions;
        Vector<double> all_stiffness, all_damping;
        if (!is_slave) {

            // Resize receive buffers based on total number of nodes
            all_positions.resize(3, total_num_nodes);
            all_stiffness.resize(total_num_nodes);
            all_damping.resize(total_num_nodes);

            // Gather all positions to master using gatherv
            for (int d = 0; d < 3; d++) {
                Vector<double> local_pos_d(local_num_nodes_);
                Vector<double> all_pos_d(total_num_nodes);
                for (int i = 0; i < local_num_nodes_; i++) {
                    local_pos_d(i) = local_positions(d,i);
                }
                cm.gatherv(cm_mod, local_pos_d, all_pos_d, proc_num_nodes, displs, 0);
                for (int i = 0; i < total_num_nodes; i++) {
                    all_positions(d,i) = all_pos_d(i);
                }
            }

            // Get VTP points for position matching
            Array<double> vtp_points = vtp_data_->get_points();
            
            // Look up data for all nodes using point matching
            for (int i = 0; i < total_num_nodes; i++) {
                int vtp_idx = find_vtp_point_index(all_positions(0,i), all_positions(1,i), all_positions(2,i), vtp_points);
                all_stiffness(i) = global_stiffness_(vtp_idx);
                all_damping(i) = global_damping_(vtp_idx);

                #ifdef debug_robin_bc_data_distribute
                dmsg << "Global stiffness value at node " << i << ": " << all_stiffness(i) << std::endl;
                dmsg << "Global damping value at node " << i << ": " << all_damping(i) << std::endl;
                #endif
            }

            // Clear global arrays to save memory
            global_stiffness_.clear();
            global_damping_.clear();

        } else {
            // Slave processes: send positions to master
            for (int d = 0; d < 3; d++) {
                Vector<double> local_pos_d(local_num_nodes_);
                for (int i = 0; i < local_num_nodes_; i++) {
                    local_pos_d(i) = local_positions(d,i);
                }
                Vector<double> dummy_recv(total_num_nodes);
                cm.gatherv(cm_mod, local_pos_d, dummy_recv, proc_num_nodes, displs, 0);
            }
        }

        // Scatter data back to all processes using scatterv
        cm.scatterv(cm_mod, all_stiffness, proc_num_nodes, displs, local_stiffness_, 0);
        cm.scatterv(cm_mod, all_damping, proc_num_nodes, displs, local_damping_, 0);

        // Build mapping from face global node IDs to local array indices so we can
        // get data from a global node ID
        global_node_map_.clear();
        for (int i = 0; i < local_num_nodes_; i++) {
            global_node_map_[local_global_ids(i)] = i;
        }

        // Check if local arrays and node positions are consistent
        #ifdef debug_robin_bc_data_distribute
        dmsg << "Checking if local arrays and node positions are consistent" << std::endl;
        for (int i = 0; i < local_num_nodes_; i++) {
            dmsg << "Local global ID: " << local_global_ids(i) << std::endl;
            dmsg << "Local index: " << get_local_index(local_global_ids(i)) << std::endl;
            dmsg << "Local position: " << com_mod.x.col(local_global_ids(i)) << std::endl;
            dmsg << "Local stiffness: " << local_stiffness_(i) << std::endl;
            dmsg << "Local damping: " << local_damping_(i) << std::endl;
        }
        #endif

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
        local_stiffness_.resize(1);
        local_damping_.resize(1);
        local_stiffness_(0) = uniform_stiffness;
        local_damping_(0) = uniform_damping;
    }

    // Mark as defined
    defined_ = true;

    #ifdef debug_robin_bc_data_distribute
    dmsg << "Finished distributing Robin BC data" << std::endl;
    dmsg << "Number of face nodes on this processor: " << local_num_nodes_ << std::endl;
    #endif
}