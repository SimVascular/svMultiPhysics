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

#include "BC.h"
#include "ComMod.h"
#include "DebugMsg.h"
#include "Vector.h"
#include <stdexcept>
#include <iostream>
#include <cstdio>
#include <vector>

#define debug_bc

BC::BC(const std::string& vtp_file_path, const std::vector<std::string>& array_names, const faceType& face)
{
    init_from_vtp(vtp_file_path, array_names, face);
}

BC::BC(const std::map<std::string, double>& uniform_values, const faceType& face)
{
    init_uniform(uniform_values, face);
}

BC::BC(const BC& other)
    : face_(other.face_)
    , global_num_nodes_(other.global_num_nodes_)
    , local_num_nodes_(other.local_num_nodes_)
    , array_names_(other.array_names_)
    , local_data_(other.local_data_)
    , global_data_(other.global_data_)
    , spatially_variable(other.spatially_variable)
    , vtp_file_path_(other.vtp_file_path_)
    , global_node_map_(other.global_node_map_)
    , defined_(other.defined_)
{
    if (other.vtp_data_) {
        vtp_data_ = std::make_unique<VtkVtpData>(*other.vtp_data_);
    }
}

BC& BC::operator=(const BC& other)
{
    if (this != &other) {
        face_ = other.face_;
        global_num_nodes_ = other.global_num_nodes_;
        local_num_nodes_ = other.local_num_nodes_;
        array_names_ = other.array_names_;
        local_data_ = other.local_data_;
        global_data_ = other.global_data_;
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

BC::BC(BC&& other) noexcept
    : face_(other.face_)
    , global_num_nodes_(other.global_num_nodes_)
    , local_num_nodes_(other.local_num_nodes_)
    , array_names_(std::move(other.array_names_))
    , local_data_(std::move(other.local_data_))
    , global_data_(std::move(other.global_data_))
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

BC& BC::operator=(BC&& other) noexcept
{
    if (this != &other) {
        face_ = other.face_;
        global_num_nodes_ = other.global_num_nodes_;
        local_num_nodes_ = other.local_num_nodes_;
        array_names_ = std::move(other.array_names_);
        local_data_ = std::move(other.local_data_);
        global_data_ = std::move(other.global_data_);
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

void BC::init_uniform(const std::map<std::string, double>& uniform_values, const faceType& face)
{
    face_ = &face;
    global_num_nodes_ = face_->nNo;
    local_num_nodes_ = 0;
    spatially_variable = false;
    vtp_file_path_ = "";
    defined_ = true;

    // Store array names
    array_names_.clear();
    for (const auto& [name, value] : uniform_values) {
        array_names_.push_back(name);
        
        // Validate value
        validate_array_value(name, value);
    }

    #ifdef debug_bc
    DebugMsg dmsg(__func__, 0);
    dmsg << "Constructor with uniform values" << std::endl;
    dmsg << "Initializing uniform BC values" << std::endl;
    #endif

    // Initialize local data with uniform values
    for (const auto& [name, value] : uniform_values) {
        Vector<double> uniform_array(1);
        uniform_array(0) = value;
        local_data_[name] = Array<double>(1, 1);
        local_data_[name].set_col(0, uniform_array);
    }
}

void BC::init_from_vtp(const std::string& vtp_file_path, const std::vector<std::string>& array_names, const faceType& face)
{
    // Note that this function is only called by the master process

    #ifdef debug_bc
    DebugMsg dmsg(__func__, 0);
    dmsg << "Initializing BC from VTP file" << std::endl;
    dmsg << "VTP file path: " << vtp_file_path << std::endl;
    dmsg << "Face name: " << face.name << std::endl;
    dmsg << "Face nNo: " << face.nNo << std::endl;
    #endif

    #ifdef debug_bc
    dmsg << "Setting face pointer" << std::endl;
    #endif
    face_ = &face;
    if (!face_) {
        throw std::runtime_error("Face pointer is null after assignment");
    }
    #ifdef debug_bc
    dmsg << "Getting number of nodes from face" << std::endl;
    #endif
    global_num_nodes_ = face_->nNo;
    local_num_nodes_ = 0;
    spatially_variable = true;
    vtp_file_path_ = vtp_file_path;
    array_names_ = array_names;
    defined_ = true;

    #ifdef debug_bc
    dmsg << "Constructor with VTP file path" << std::endl;
    dmsg << "VTP file path: " << vtp_file_path << std::endl;
    dmsg << "Face name: " << face.name << std::endl;
    dmsg << "Array names: " << array_names[0] << " and " << array_names[1] << std::endl;
    #endif

    // Read data from VTP file
    global_data_ = read_data_from_vtp_file(vtp_file_path, array_names);

    // Validate values
    for (const auto& [name, data] : global_data_) {
        for (int i = 0; i < global_num_nodes_; i++) {
            validate_array_value(name, data(i, 0));
        }
    }

    // In case we are running sequentially, we need to fill the local arrays 
    // and the global node map as well, because distribute is not called in sequential mode.
    local_data_ = global_data_;

    // Initialize global node map, in case we are running sequentially
    global_node_map_.clear();
    for (int i = 0; i < global_num_nodes_; i++) {
        global_node_map_[face_->gN(i)] = i;
    }
}

BC::StringArrayMap BC::read_data_from_vtp_file(const std::string& vtp_file_path, const std::vector<std::string>& array_names)
{
    #ifdef debug_bc
    DebugMsg dmsg(__func__, 0);
    dmsg << "Loading data from VTP file: " << vtp_file_path << std::endl;
    dmsg << "Array names: " << array_names[0] << " and " << array_names[1] << std::endl;
    #endif

    // Check if file exists
    #ifdef debug_bc
    dmsg << "Checking if file exists: " << vtp_file_path << std::endl;
    #endif
    if (FILE *file = fopen(vtp_file_path.c_str(), "r")) {
        fclose(file);
        #ifdef debug_bc
        dmsg << "File exists and is readable" << std::endl;
        #endif
    } else {
        throw BCFileException(vtp_file_path);
    }
    
    // Read the VTP file
    #ifdef debug_bc
    dmsg << "Creating VtkVtpData object" << std::endl;
    #endif
    try {
        vtp_data_ = std::make_unique<VtkVtpData>(vtp_file_path, true);
        #ifdef debug_bc
        dmsg << "VtkVtpData object created successfully" << std::endl;
        #endif
    } catch (const std::exception& e) {
        throw BCFileException(vtp_file_path);
    }
    
    // Check if the number of nodes in the VTP file matches the number of nodes on the face
    if (global_num_nodes_ != face_->nNo) {
        throw BCNodeCountException(vtp_file_path, face_->name);
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

        #ifdef debug_bc
        dmsg << "Successfully loaded " << array_name << " data" << std::endl;
        dmsg << array_name << " data size: " << array_data.nrows() << " x " << array_data.ncols() << std::endl;
        #endif
    }

    #ifdef debug_bc
    dmsg << "Finished loading VTP data" << std::endl;
    #endif

    return result;
}

double BC::get_value(const std::string& array_name, int node_id) const
{
    // Check if array exists
    auto it = local_data_.find(array_name);
    if (it == local_data_.end()) {
        throw std::runtime_error("Array '" + array_name + "' not found.");
    }

    // Check node ID range
    if (node_id < 0 || node_id >= global_num_nodes_) {
        throw std::runtime_error("Node ID " + std::to_string(node_id) + 
                                " is out of range [0, " + std::to_string(global_num_nodes_ - 1) + "].");
    }

    // Return value
    if (spatially_variable) {
        return it->second(node_id, 0);
    } else {
        return it->second(0, 0);
    }
}

int BC::get_local_index(int global_node_id) const
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

int BC::find_vtp_point_index(double x, double y, double z,
                                const Array<double>& vtp_points) const
{
    const int num_points = vtp_points.ncols();
    Vector<double> target_point{x, y, z};
    
    // Simple linear search through all points in the VTP file
    for (int i = 0; i < num_points; i++) {
        auto vtp_point = vtp_points.col(i);
        auto diff = vtp_point - target_point;
        
        // Compute distance as sqrt of dot product of diff with itself
        double dist = sqrt(diff.dot(diff));
        
        // Compare squared distance with squared tolerance to avoid sqrt
        if (dist <= POINT_MATCH_TOLERANCE) {
            #define n_debug_bc_find_vtp_point_index
            #ifdef debug_bc_find_vtp_point_index
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

void BC::distribute(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, const faceType& face)
{
    #define n_debug_bc_distribute
    #ifdef debug_bc_distribute
    DebugMsg dmsg(__func__, cm.idcm());
    dmsg << "Distributing BC data" << std::endl;
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

    // Communicate array names
    int num_arrays = array_names_.size();
    cm.bcast(cm_mod, &num_arrays);
    if (is_slave) {
        array_names_.resize(num_arrays);
    }
    for (int i = 0; i < num_arrays; i++) {
        if (!is_slave) {
            std::string& array_name = array_names_[i];
            cm.bcast(cm_mod, array_name);
        } else {
            std::string array_name;
            cm.bcast(cm_mod, array_name);
            array_names_[i] = array_name;
        }
    }

    // Communicate array values needed by each process
    if (spatially_variable) {
        // Setup
        if (face_ == nullptr) {
            throw std::runtime_error("face_ is nullptr during distribute");
        }

        // Each processor collects the global node IDs and nodal positions of its 
        // associated face portion
        Vector<int> local_global_ids = face_->gN;
        Array<double> local_positions(3, local_num_nodes_);
        for (int i = 0; i < local_num_nodes_; i++) {
            local_positions.set_col(i, com_mod.x.col(face_->gN(i)));
        }

        #ifdef debug_bc_distribute
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
        // get the corresponding array values by matching the positions to the VTP points,
        // and scatter the data back to all processors.
        Array<double> all_positions;
        std::map<std::string, Vector<double>> all_values;
        if (!is_slave) {
            // Resize receive buffers based on total number of nodes
            all_positions.resize(3, total_num_nodes);

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
            for (const auto& array_name : array_names_) {
                all_values[array_name].resize(total_num_nodes);
                for (int i = 0; i < total_num_nodes; i++) {
                    int vtp_idx = find_vtp_point_index(all_positions(0,i), all_positions(1,i), all_positions(2,i), vtp_points);
                    all_values[array_name](i) = global_data_[array_name](vtp_idx, 0);
                }
            }

            // Clear global data to save memory
            global_data_.clear();

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
        local_data_.clear();
        for (const auto& array_name : array_names_) {
            Vector<double> local_values(local_num_nodes_);
            cm.scatterv(cm_mod, all_values[array_name], proc_num_nodes, displs, local_values, 0);
            local_data_[array_name] = Array<double>(local_num_nodes_, 1);
            local_data_[array_name].set_col(0, local_values);
        }

        // Build mapping from face global node IDs to local array indices so we can
        // get data from a global node ID
        global_node_map_.clear();
        for (int i = 0; i < local_num_nodes_; i++) {
            global_node_map_[local_global_ids(i)] = i;
        }

        // Check if local arrays and node positions are consistent
        #ifdef debug_bc_distribute
        dmsg << "Checking if local arrays and node positions are consistent" << std::endl;
        for (int i = 0; i < local_num_nodes_; i++) {
            dmsg << "Local global ID: " << local_global_ids(i) << std::endl;
            dmsg << "Local index: " << get_local_index(local_global_ids(i)) << std::endl;
            dmsg << "Local position: " << com_mod.x.col(local_global_ids(i)) << std::endl;
            for (const auto& array_name : array_names_) {
                dmsg << "Local " << array_name << ": " << local_data_[array_name](i, 0) << std::endl;
            }
        }
        #endif

    } else {
        // For uniform values, just broadcast the single values
        if (!is_slave) {
            for (const auto& array_name : array_names_) {
                double uniform_value = local_data_[array_name](0, 0);
                cm.bcast(cm_mod, &uniform_value);
            }
        } else {
            local_data_.clear();
            for (const auto& array_name : array_names_) {
                double uniform_value;
                cm.bcast(cm_mod, &uniform_value);
                local_data_[array_name] = Array<double>(1, 1);
                local_data_[array_name](0, 0) = uniform_value;
            }
        }
    }

    // Mark as defined
    defined_ = true;

    #ifdef debug_bc_distribute
    dmsg << "Finished distributing BC data" << std::endl;
    dmsg << "Number of face nodes on this processor: " << local_num_nodes_ << std::endl;
    #endif
}
