// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "BoundaryCondition.h"
#include "ComMod.h"
#include "DebugMsg.h"
#include "Vector.h"
#include <stdexcept>
#include <iostream>
#include <cstdio>
#include <vector>
#include <utility>

#define n_debug_bc

BoundaryCondition::BoundaryCondition(const std::string& vtp_file_path, const std::vector<std::string>& array_names, const StringBoolMap& flags, const faceType& face, SimulationLogger& logger)
    : face_(&face)
    , global_num_nodes_(face.nNo)
    , local_num_nodes_(0)
    , array_names_(array_names)
    , spatially_variable(true)
    , vtp_file_path_(vtp_file_path)
    , flags_(flags)
    , logger_(&logger)
{
    try {
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

        global_node_map_.clear();
        for (int i = 0; i < global_num_nodes_; i++) {
            global_node_map_[face_->gN(i)] = i;
        }
    } catch (const std::exception& e) {
        // Constructor failed - resources are automatically cleaned up
        throw;  // Re-throw the exception
    }
}

BoundaryCondition::BoundaryCondition(const StringDoubleMap& uniform_values, const StringBoolMap& flags, const faceType& face, SimulationLogger& logger)
    : face_(&face)
    , global_num_nodes_(face.nNo)
    , local_num_nodes_(0)
    , spatially_variable(false)
    , vtp_file_path_("")
    , flags_(flags)
    , logger_(&logger)
{
    try {
        // Store array names, validate and store values
        array_names_.clear();
        for (const auto& [name, value] : uniform_values) {
            array_names_.push_back(name);
            validate_array_value(name, value);
            local_data_[name] = Array<double>(1, 1);
            local_data_[name](0, 0) = value;
        }
    } catch (const std::exception& e) {
        // Constructor failed - resources are automatically cleaned up
        throw;  // Re-throw the exception
    }
}

BoundaryCondition::BoundaryCondition(const BoundaryCondition& other)
    : face_(other.face_)
    , global_num_nodes_(other.global_num_nodes_)
    , local_num_nodes_(other.local_num_nodes_)
    , array_names_(other.array_names_)
    , local_data_(other.local_data_)
    , global_data_(other.global_data_)
    , spatially_variable(other.spatially_variable)
    , vtp_file_path_(other.vtp_file_path_)
    , flags_(other.flags_)
    , global_node_map_(other.global_node_map_)
    , logger_(other.logger_)
{
    if (other.vtp_data_) {
        vtp_data_ = std::make_unique<VtkVtpData>(*other.vtp_data_);
    }
}

void swap(BoundaryCondition& lhs, BoundaryCondition& rhs) noexcept {
    using std::swap;
    swap(lhs.face_, rhs.face_);
    swap(lhs.global_num_nodes_, rhs.global_num_nodes_);
    swap(lhs.local_num_nodes_, rhs.local_num_nodes_);
    swap(lhs.array_names_, rhs.array_names_);
    swap(lhs.local_data_, rhs.local_data_);
    swap(lhs.global_data_, rhs.global_data_);
    swap(lhs.spatially_variable, rhs.spatially_variable);
    swap(lhs.vtp_file_path_, rhs.vtp_file_path_);
    swap(lhs.flags_, rhs.flags_);
    swap(lhs.global_node_map_, rhs.global_node_map_);
    swap(lhs.vtp_data_, rhs.vtp_data_);
    swap(lhs.logger_, rhs.logger_);
}

BoundaryCondition& BoundaryCondition::operator=(BoundaryCondition other) {
    swap(*this, other);
    return *this;
}

BoundaryCondition::BoundaryCondition(BoundaryCondition&& other) noexcept
    : face_(other.face_)
    , global_num_nodes_(other.global_num_nodes_)
    , local_num_nodes_(other.local_num_nodes_)
    , array_names_(std::move(other.array_names_))
    , local_data_(std::move(other.local_data_))
    , global_data_(std::move(other.global_data_))
    , spatially_variable(other.spatially_variable)
    , vtp_file_path_(std::move(other.vtp_file_path_))
    , flags_(std::move(other.flags_))
    , global_node_map_(std::move(other.global_node_map_))
    , vtp_data_(std::move(other.vtp_data_))
    , logger_(other.logger_)
{
    other.face_ = nullptr;
    other.global_num_nodes_ = 0;
    other.local_num_nodes_ = 0;
    other.spatially_variable = false;
}


BoundaryCondition::StringArrayMap BoundaryCondition::read_data_from_vtp_file(const std::string& vtp_file_path, const std::vector<std::string>& array_names)
{
    #ifdef debug_bc
    DebugMsg dmsg(__func__, 0);
    dmsg << "Loading data from VTP file: " << vtp_file_path << std::endl;
    dmsg << "Array names: " << array_names[0] << " and " << array_names[1] << std::endl;
    #endif

    // Check if file exists
    if (FILE *file = fopen(vtp_file_path.c_str(), "r")) {
        fclose(file);
    } else {
        throw BoundaryConditionFileException(vtp_file_path);
    }
    
    // Read the VTP file
    try {
        vtp_data_ = std::make_unique<VtkVtpData>(vtp_file_path, true);
    } catch (const std::exception& e) {
        throw BoundaryConditionFileException(vtp_file_path);
    }
    
    if (global_num_nodes_ != face_->nNo) {
        throw BoundaryConditionNodeCountException(vtp_file_path, face_->name);
    }

    // Read in the data from the VTP file
    StringArrayMap result;
    for (const auto& array_name : array_names) {
        if (!vtp_data_->has_point_data(array_name)) {
            throw BoundaryConditionVtpArrayException(vtp_file_path, array_name);
        }

        auto array_data = vtp_data_->get_point_data(array_name);

        if (array_data.nrows() != global_num_nodes_ || array_data.ncols() != 1) {
            throw BoundaryConditionVtpArrayDimensionException(vtp_file_path, array_name, 
                                                           global_num_nodes_, 1, 
                                                           array_data.nrows(), array_data.ncols());
        }

        // Store array in result map
        result[array_name] = array_data;

        #ifdef debug_bc
        dmsg << "Successfully loaded " << array_name << " data" << std::endl;
        dmsg << array_name << " data size: " << array_data.nrows() << " x " << array_data.ncols() << std::endl;
        #endif
    }

    logger_ -> log_message("[BoundaryCondition] Loaded from VTP file");
    logger_ -> log_message("\t File path:", vtp_file_path);
    logger_ -> log_message("\t Arrays:", array_names);
    logger_ -> log_message("\t Face:", face_->name);

    return result;
}

double BoundaryCondition::get_value(const std::string& array_name, int node_id) const
{
    auto it = local_data_.find(array_name);
    if (it == local_data_.end()) {
        throw BoundaryConditionArrayException(array_name);
    }

    if (node_id < 0 || node_id >= global_num_nodes_) {
        throw BoundaryConditionNodeIdException(node_id, global_num_nodes_);
    }

    // Return value
    if (spatially_variable) {
        return it->second(node_id, 0);
    } else {
        return it->second(0, 0);
    }
}

bool BoundaryCondition::get_flag(const std::string& name) const
{
    auto it = flags_.find(name);
    if (it == flags_.end()) {
        #ifdef debug_bc
        DebugMsg dmsg(__func__, 0);
        dmsg << "Flag '" << name << "' not found. Available flags: " << flags_to_string() << std::endl;
        #endif
        throw BoundaryConditionFlagException(name);
    }
    return it->second;
}

int BoundaryCondition::get_local_index(int global_node_id) const
{
    if (spatially_variable) {
        auto it = global_node_map_.find(global_node_id);
        if (it == global_node_map_.end()) {
            throw BoundaryConditionGlobalNodeIdException(global_node_id);
        }
        return it->second;
    } else {
        return 0;
    }
}

void BoundaryCondition::distribute(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, const faceType& face)
{
    #define n_debug_distribute
    #ifdef debug_distribute
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
    const bool is_slave = cm.slv(cm_mod);
    distribute_metadata(cm_mod, cm, is_slave);
    if (spatially_variable) {
        distribute_spatially_variable(com_mod, cm_mod, cm, is_slave);
    } else {
        distribute_uniform(cm_mod, cm, is_slave);
    }
    distribute_flags(cm_mod, cm, is_slave);

    #ifdef debug_distribute
    dmsg << "Finished distributing BC data" << std::endl;
    dmsg << "Number of face nodes on this processor: " << local_num_nodes_ << std::endl;
    #endif
}

void BoundaryCondition::distribute_metadata(const CmMod& cm_mod, const cmType& cm, bool is_slave)
{
    cm.bcast(cm_mod, &spatially_variable);

    // Not necessary, but we do it for consistency
    if (spatially_variable) {
        cm.bcast(cm_mod, vtp_file_path_);
    }

    // Not necessary, but we do it for consistency
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
}

void BoundaryCondition::distribute_spatially_variable(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, bool is_slave)
{
    #define n_debug_distribute_spatially_variable
    #ifdef debug_distribute_spatially_variable
    DebugMsg dmsg(__func__, 0);
    #endif

    if (face_ == nullptr) {
        throw BoundaryConditionNullFaceException();
    }
    // Each processor collects the global node IDs and nodal positions of its 
    // associated face portion
    Vector<int> local_global_ids = face_->gN;
    Array<double> local_positions(3, local_num_nodes_);
    for (int i = 0; i < local_num_nodes_; i++) {
        local_positions.set_col(i, com_mod.x.col(face_->gN(i)));
    }

    #ifdef debug_distribute_spatially_variable
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
        
        // Get mesh scale factor from the face's mesh
        double mesh_scale_factor = 1.0; // Default scale factor
        if (face_ != nullptr) {
            mesh_scale_factor = com_mod.msh[face_->iM].scF;
            #ifdef debug_distribute_spatially_variable
            dmsg << "Mesh scale factor: " << mesh_scale_factor << std::endl;
            #endif
        }
        
        // Look up data for all nodes using point matching
        for (const auto& array_name : array_names_) {
            all_values[array_name].resize(total_num_nodes);
            for (int i = 0; i < total_num_nodes; i++) {
                int vtp_idx = find_vtp_point_index(all_positions(0,i), all_positions(1,i), all_positions(2,i), vtp_points, mesh_scale_factor);
                all_values[array_name](i) = global_data_[array_name](vtp_idx, 0);
            }
        }
        
        // Clear global data to save memory
        global_data_.clear();
    } else {
        // Slave processes: send node positions to master
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

    #ifdef debug_distribute_spatially_variable
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
}

void BoundaryCondition::distribute_uniform(const CmMod& cm_mod, const cmType& cm, bool is_slave)
{
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

void BoundaryCondition::distribute_flags(const CmMod& cm_mod, const cmType& cm, bool is_slave)
{
    if (cm.seq()) return;
    int num_flags = 0;
    if (!is_slave) {
        num_flags = static_cast<int>(flags_.size());
    }
    cm.bcast(cm_mod, &num_flags);
    if (is_slave) {
        flags_.clear();
    }
    for (int i = 0; i < num_flags; i++) {
        std::string key;
        bool val = false;
        if (!is_slave) {
            auto it = std::next(flags_.begin(), i);
            key = it->first;
            val = it->second;
            cm.bcast(cm_mod, key);
            cm.bcast(cm_mod, &val);
        } else {
            cm.bcast(cm_mod, key);
            cm.bcast(cm_mod, &val);
            flags_[key] = val;
        }
    }
}

int BoundaryCondition::find_vtp_point_index(double x, double y, double z,
                                const Array<double>& vtp_points, double mesh_scale_factor) const
{
    const int num_points = vtp_points.ncols();
    
    // Scale down the target coordinates to match the unscaled VTP coordinates
    // The simulation coordinates are scaled by mesh_scale_factor, but VTP coordinates are not
    Vector<double> target_point{x / mesh_scale_factor, y / mesh_scale_factor, z / mesh_scale_factor};

    // Simple linear search through all points in the VTP file
    for (int i = 0; i < num_points; i++) {
        auto vtp_point = vtp_points.col(i);
        auto diff = vtp_point - target_point;
        double distance = sqrt(diff.dot(diff));

        if (distance <= POINT_MATCH_TOLERANCE) {
            #define n_debug_bc_find_vtp_point_index
            #ifdef debug_bc_find_vtp_point_index
            DebugMsg dmsg(__func__, 0);
            dmsg << "Found VTP point index for node at position (" << x << ", " << y << ", " << z << ")" << std::endl;
            dmsg << "Scaled target position (" << target_point(0) << ", " << target_point(1) << ", " << target_point(2) << ")" << std::endl;
            dmsg << "VTP point index: " << i << std::endl;
            #endif

            return i;
        }
    }

    throw BoundaryConditionPointNotFoundException(x, y, z);
}

std::string BoundaryCondition::flags_to_string() const {
    std::string result;
    for (const auto& [name, value] : flags_) {
        result += name + ": " + (value ? "true" : "false") + ", ";
    }
    return result;
}