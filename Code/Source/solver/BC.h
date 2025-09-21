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

#ifndef BC_H
#define BC_H

#include "Array.h"
#include "Vector.h"
#include "CmMod.h"
#include "VtkData.h"
#include <string>
#include <memory>
#include <map>
#include <vector>

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
/// class MyBC : public BC {
/// public:
///     MyBC(const std::string& vtp_file_path, const std::vector<std::string>& array_names, const faceType& face)
///         : BC(vtp_file_path, array_names, face) {}
///
///     MyBC(const std::map<std::string, double>& uniform_values, const faceType& face)
///         : BC(uniform_values, face) {}
///
///     // Add BC-specific functionality here
/// };
/// ```
class BC {
public:
    /// @brief Type alias for map of array names to array data
    using StringArrayMap = std::map<std::string, Array<double>>;

    /// @brief Default constructor - creates an empty BC
    BC() : face_(nullptr), global_num_nodes_(0), local_num_nodes_(0), spatially_variable(false), defined_(false) {}

    /// @brief Constructor - reads data from VTP file
    /// @param vtp_file_path Path to VTP file containing arrays
    /// @param array_names Names of arrays to read from VTP file
    /// @param face Face associated with the BC
    /// @throws std::runtime_error if file cannot be read or arrays are missing
    BC(const std::string& vtp_file_path, const std::vector<std::string>& array_names, const faceType& face);

    /// @brief Constructor for uniform values
    /// @param uniform_values Map of array names to uniform values
    /// @param face Face associated with the BC
    BC(const std::map<std::string, double>& uniform_values, const faceType& face);

    /// @brief Copy constructor
    BC(const BC& other);

    /// @brief Copy assignment operator
    BC& operator=(const BC& other);

    /// @brief Move constructor
    BC(BC&& other) noexcept;

    /// @brief Move assignment operator
    BC& operator=(BC&& other) noexcept;

    /// @brief Virtual destructor
    virtual ~BC() = default;

    /// @brief Initialize from VTP file
    /// @param vtp_file_path Path to VTP file containing arrays
    /// @param array_names Names of arrays to read from VTP file
    /// @param face Face associated with the BC
    /// @throws std::runtime_error if file cannot be read or arrays are missing
    void init_from_vtp(const std::string& vtp_file_path, const std::vector<std::string>& array_names, const faceType& face);

    /// @brief Initialize with uniform values
    /// @param uniform_values Map of array names to uniform values
    /// @param face Face associated with the BC
    void init_uniform(const std::map<std::string, double>& uniform_values, const faceType& face);

    /// @brief Get value for a specific array and node
    /// @param array_name Name of the array
    /// @param node_id Node index on the face
    /// @return Value for the array at the specified node
    double get_value(const std::string& array_name, int node_id) const;

    /// @brief Get global number of nodes
    /// @return Global number of nodes on the face
    int get_global_num_nodes() const { return global_num_nodes_; }

    /// @brief Get local number of nodes
    /// @return Local number of nodes on the face on this processor
    int get_local_num_nodes() const { return local_num_nodes_; }

    /// @brief Get local array index for a global node ID
    /// @param global_node_id The global node ID defined on the face
    /// @return Local array index for data arrays
    /// @throws std::runtime_error if global_node_id is not found in the map
    int get_local_index(int global_node_id) const;

    /// @brief Check if data is loaded from VTP file
    /// @return true if loaded from VTP, false if using uniform values
    bool is_from_vtp() const { return spatially_variable; }

    /// @brief Get the VTP file path (empty if using uniform values)
    /// @return VTP file path
    const std::string& get_vtp_path() const { return vtp_file_path_; }

    /// @brief Check if this BC is properly defined
    /// @return true if BC has been initialized with either VTP data or uniform values
    bool is_defined() const { return defined_; }

    /// @brief Distribute BC data from the master process to the slave processes
    /// @param com_mod Reference to ComMod object for global coordinates
    /// @param cm_mod Reference to CmMod object for MPI communication
    /// @param cm Reference to cmType object for MPI communication
    /// @param face Face associated with the BC
    void distribute(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, const faceType& face);

protected:
    /// @brief Data members for BC
    const faceType* face_;             ///< Face associated with the BC (can be null)
    int global_num_nodes_;             ///< Global number of nodes on the face
    int local_num_nodes_;              ///< Local number of nodes on this processor
    std::vector<std::string> array_names_; ///< Names of arrays to read from VTP file
    StringArrayMap local_data_;        ///< Local array values for each node on this processor
    StringArrayMap global_data_;       ///< Global array values (only populated on master)
    bool spatially_variable;           ///< Flag indicating if data is from VTP file
    std::string vtp_file_path_;        ///< Path to VTP file (empty if uniform)
    std::map<int, int> global_node_map_; ///< Maps global node IDs to local array indices
    std::unique_ptr<VtkVtpData> vtp_data_;  ///< VTP data object
    bool defined_;                     ///< Whether this BC has been properly initialized

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
    /// @return Index of the matching point in the VTP array
    /// @throws std::runtime_error if no matching point is found
    int find_vtp_point_index(double x, double y, double z, const Array<double>& vtp_points) const;

    /// @brief Hook for derived classes to validate array values
    /// @param array_name Name of the array being validated
    /// @param value Value to validate
    /// @throws std::runtime_error if validation fails
    virtual void validate_array_value(const std::string& array_name, double value) const {}
};

#endif // BC_H
