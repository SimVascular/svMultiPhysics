// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

#include "Array.h"
#include "Vector.h"
#include "CmMod.h"
#include "VtkData.h"
#include <string>
#include "SimulationLogger.h"
#include <memory>
#include <map>
#include <vector>
#include <stdexcept>
#include <utility>

// Forward declarations. These are needed because including ComMod.h causes a 
// circular header dependency.
class faceType;
class ComMod;

/// @brief Base class for boundary conditions with spatially variable arrays
/// 
/// This class provides common functionality for boundary conditions that need to
/// read and manage arrays of values from VTP files or uniform values. It handles
/// distribution of data across processes and provides efficient access to values.
///
/// This class is intended to be subclassed by specific boundary condition types.
///
/// Development note: this class is intended to eventually be an object-oriented
/// replacement of the existing bcType, although it is not yet complete.
///
/// Example usage:
/// ```cpp
/// class MyBoundaryCondition : public BoundaryCondition {
/// public:
///     MyBoundaryCondition(const std::string& vtp_file_path, const std::vector<std::string>& array_names, const faceType& face)
///         : BoundaryCondition(vtp_file_path, array_names, face) {}
///
///     MyBoundaryCondition(const std::map<std::string, double>& uniform_values, const faceType& face)
///         : BoundaryCondition(uniform_values, face) {}
///
///     // Add BC-specific functionality here
/// };
/// ```
class BoundaryCondition {
protected:
    /// @brief Type alias for map of array names to array data
    using StringArrayMap = std::map<std::string, Array<double>>;
    using StringBoolMap = std::map<std::string, bool>;
    using StringDoubleMap = std::map<std::string, double>;

    /// @brief Data members for BC
    const faceType* face_ = nullptr;         ///< Face associated with the BC (not owned by BoundaryCondition)
    int global_num_nodes_ = 0;               ///< Global number of nodes on the face
    int local_num_nodes_ = 0;                ///< Local number of nodes on this processor
    std::vector<std::string> array_names_;   ///< Names of arrays to read from VTP file
    StringArrayMap local_data_;              ///< Local array values for each node on this processor
    StringArrayMap global_data_;             ///< Global array values (only populated on master)
    StringBoolMap flags_;                    ///< Named boolean flags for BC behavior
    bool spatially_variable = false;         ///< Flag indicating if data is from VTP file
    std::string vtp_file_path_;              ///< Path to VTP file (empty if uniform)
    std::map<int, int> global_node_map_;     ///< Maps global node IDs to local array indices
    std::unique_ptr<VtkVtpData> vtp_data_;   ///< VTP data object
    const SimulationLogger* logger_ = nullptr;  ///< Logger for warnings/info (not owned by BoundaryCondition)

public:
    /// @brief Tolerance for point matching in VTP files
    static constexpr double POINT_MATCH_TOLERANCE = 1e-12;

    /// @brief Default constructor - creates an empty BC
    BoundaryCondition() = default;

    /// @brief Constructor - reads data from VTP file
    /// @param vtp_file_path Path to VTP file containing arrays
    /// @param array_names Names of arrays to read from VTP file
    /// @param face Face associated with the BC
    /// @param logger Simulation logger used to write warnings
    /// @throws std::runtime_error if file cannot be read or arrays are missing
    BoundaryCondition(const std::string& vtp_file_path, const std::vector<std::string>& array_names, const StringBoolMap& flags, const faceType& face, SimulationLogger& logger);

    /// @brief Constructor for uniform values
    /// @param uniform_values Map of array names to uniform values
    /// @param face Face associated with the BC
    /// @param logger Simulation logger used to write warnings
    BoundaryCondition(const StringDoubleMap& uniform_values, const StringBoolMap& flags, const faceType& face, SimulationLogger& logger);

    /// @brief Copy constructor
    BoundaryCondition(const BoundaryCondition& other);

    /// @brief Unified assignment operator (handles both copy and move)
    BoundaryCondition& operator=(BoundaryCondition other);

    /// @brief Move constructor
    BoundaryCondition(BoundaryCondition&& other) noexcept;

    /// @brief Virtual destructor
    virtual ~BoundaryCondition() noexcept = default;

    /// @brief Swap function for copy-and-swap idiom (friend function)
    friend void swap(BoundaryCondition& lhs, BoundaryCondition& rhs) noexcept;


    /// @brief Get value for a specific array and node
    /// @param array_name Name of the array
    /// @param node_id Node index on the face
    /// @return Value for the array at the specified node
    /// @throws std::runtime_error if array_name is not found
    double get_value(const std::string& array_name, int node_id) const;

    /// @brief Get a boolean flag by name
    /// @param name Name of the flag
    /// @return Value of the flag
    /// @throws std::runtime_error if flag is not found
    bool get_flag(const std::string& name) const;

    /// @brief Get a string representation of the flags
    /// @return String representation of the flags
    std::string flags_to_string() const;

    /// @brief Get global number of nodes
    /// @return Global number of nodes on the face
    int get_global_num_nodes() const noexcept {
        return global_num_nodes_;
    }

    /// @brief Get local number of nodes
    /// @return Local number of nodes on the face on this processor
    int get_local_num_nodes() const noexcept {
        return local_num_nodes_;
    }

    /// @brief Get local array index for a global node ID
    /// @param global_node_id The global node ID defined on the face
    /// @return Local array index for data arrays
    /// @throws std::runtime_error if global_node_id is not found in the map
    int get_local_index(int global_node_id) const;

    /// @brief Check if data is loaded from VTP file
    /// @return true if loaded from VTP, false if using uniform values
    bool is_from_vtp() const noexcept {
        return spatially_variable;
    }

    /// @brief Get the VTP file path (empty if using uniform values)
    /// @return VTP file path
    const std::string& get_vtp_path() const noexcept {
        return vtp_file_path_;
    }

    /// @brief Check if this BC has been properly initialized with data
    /// @return true if BC has data (either global or local arrays are populated)
    bool is_initialized() const noexcept {
        return !global_data_.empty() || !local_data_.empty();
    }

    /// @brief Distribute BC data from the master process to the slave processes
    /// @param com_mod Reference to ComMod object for global coordinates
    /// @param cm_mod Reference to CmMod object for MPI communication
    /// @param cm Reference to cmType object for MPI communication
    /// @param face Face associated with the BC
    void distribute(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, const faceType& face);

protected:
    /// @brief Read data from VTP file
    /// @param vtp_file_path Path to VTP file
    /// @param array_names Names of arrays to read
    /// @return Map of array names to array data
    StringArrayMap read_data_from_vtp_file(const std::string& vtp_file_path, const std::vector<std::string>& array_names);

    /// @brief Find index of a point in the VTP points array
    /// @param x X coordinate
    /// @param y Y coordinate
    /// @param z Z coordinate
    /// @param vtp_points VTP points array
    /// @param mesh_scale_factor Scale factor applied to mesh coordinates
    /// @return Index of the matching point in the VTP array
    /// @throws std::runtime_error if no matching point is found
    int find_vtp_point_index(double x, double y, double z, const Array<double>& vtp_points, double mesh_scale_factor) const;

    /// @brief Hook for derived classes to validate array values
    /// @param array_name Name of the array being validated
    /// @param value Value to validate
    /// @throws std::runtime_error if validation fails
    virtual void validate_array_value(const std::string& array_name, double value) const {}

     // ---- distribute helpers ----
     void distribute_metadata(const CmMod& cm_mod, const cmType& cm, bool is_slave);
     void distribute_spatially_variable(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, bool is_slave);
     void distribute_uniform(const CmMod& cm_mod, const cmType& cm, bool is_slave);
     void distribute_flags(const CmMod& cm_mod, const cmType& cm, bool is_slave);
 
};

/// @brief Base exception class for BC errors
/// 
/// These exceptions indicate fatal errors that should terminate the solver.
/// They should not be caught and handled gracefully - the solver should exit.
class BoundaryConditionBaseException : public std::exception {
public:
    /// @brief Constructor
    /// @param msg Error message
    explicit BoundaryConditionBaseException(const std::string& msg) : message_(msg) {}

    /// @brief Get error message
    /// @return Error message
    const char* what() const noexcept override {
        return message_.c_str();
    }

private:
    std::string message_;
};

/// @brief Exception thrown when VTP file cannot be read or is invalid
class BoundaryConditionFileException : public BoundaryConditionBaseException {
public:
    /// @brief Constructor
    /// @param file Path to VTP file
    explicit BoundaryConditionFileException(const std::string& file)
        : BoundaryConditionBaseException("Failed to open or read the VTP file '" + file + "'") {}
};

/// @brief Exception thrown when node count mismatch between VTP and face
class BoundaryConditionNodeCountException : public BoundaryConditionBaseException {
public:
    /// @brief Constructor
    /// @param vtp_file VTP file path
    /// @param face_name Face name
    explicit BoundaryConditionNodeCountException(const std::string& vtp_file, const std::string& face_name)
        : BoundaryConditionBaseException("Number of nodes in VTP file '" + vtp_file +
                                       "' does not match number of nodes on face '" + face_name + "'") {}
};

/// @brief Exception thrown when array validation fails
class BoundaryConditionValidationException : public BoundaryConditionBaseException {
public:
    /// @brief Constructor
    /// @param array_name Name of array that failed validation
    /// @param value Value that failed validation
    explicit BoundaryConditionValidationException(const std::string& array_name, double value) 
        : BoundaryConditionBaseException("Invalid value " + std::to_string(value) + 
                                       " for array '" + array_name + "'") {}
};

/// @brief Exception thrown when a requested flag is not defined
class BoundaryConditionFlagException : public BoundaryConditionBaseException {
public:
    explicit BoundaryConditionFlagException(const std::string& flag_name)
        : BoundaryConditionBaseException("BoundaryCondition flag not found: '" + flag_name + "'") {}
};

/// @brief Exception thrown when a requested array is not found
class BoundaryConditionArrayException : public BoundaryConditionBaseException {
public:
    explicit BoundaryConditionArrayException(const std::string& array_name)
        : BoundaryConditionBaseException("BoundaryCondition array not found: '" + array_name + "'") {}
};

/// @brief Exception thrown when BoundaryCondition is not properly initialized
class BoundaryConditionNotInitializedException : public BoundaryConditionBaseException {
public:
    explicit BoundaryConditionNotInitializedException()
        : BoundaryConditionBaseException("BoundaryCondition not properly initialized - no data available") {}
};

/// @brief Exception thrown when a node ID is out of range
class BoundaryConditionNodeIdException : public BoundaryConditionBaseException {
public:
    explicit BoundaryConditionNodeIdException(int node_id, int max_node_id)
        : BoundaryConditionBaseException("Node ID " + std::to_string(node_id) + 
                                       " is out of range [0, " + std::to_string(max_node_id - 1) + "]") {}
};

/// @brief Exception thrown when a global node ID is not found in the global-to-local map
class BoundaryConditionGlobalNodeIdException : public BoundaryConditionBaseException {
public:
    explicit BoundaryConditionGlobalNodeIdException(int global_node_id)
        : BoundaryConditionBaseException("Global node ID " + std::to_string(global_node_id) + 
                                       " not found in global-to-local map") {}
};

/// @brief Exception thrown when a VTP file doesn't contain a required array
class BoundaryConditionVtpArrayException : public BoundaryConditionBaseException {
public:
    explicit BoundaryConditionVtpArrayException(const std::string& vtp_file, const std::string& array_name)
        : BoundaryConditionBaseException("VTP file '" + vtp_file + "' does not contain '" + array_name + "' point array") {}
};

/// @brief Exception thrown when a VTP array has incorrect dimensions
class BoundaryConditionVtpArrayDimensionException : public BoundaryConditionBaseException {
public:
    explicit BoundaryConditionVtpArrayDimensionException(const std::string& vtp_file, const std::string& array_name, 
                                                       int expected_rows, int expected_cols, int actual_rows, int actual_cols)
        : BoundaryConditionBaseException("'" + array_name + "' array in VTP file '" + vtp_file +
                                       "' has incorrect dimensions. Expected " + std::to_string(expected_rows) +
                                       " x " + std::to_string(expected_cols) + ", got " + std::to_string(actual_rows) +
                                       " x " + std::to_string(actual_cols)) {}
};

/// @brief Exception thrown when face_ is nullptr during distribute
class BoundaryConditionNullFaceException : public BoundaryConditionBaseException {
public:
    explicit BoundaryConditionNullFaceException()
        : BoundaryConditionBaseException("face_ is nullptr during distribute") {}
};

/// @brief Exception thrown when a point cannot be found in VTP file
class BoundaryConditionPointNotFoundException : public BoundaryConditionBaseException {
public:
    explicit BoundaryConditionPointNotFoundException(double x, double y, double z)
        : BoundaryConditionBaseException("Could not find matching point in VTP file for node at position (" +
                                       std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")") {}
};

#endif // BOUNDARY_CONDITION_H
