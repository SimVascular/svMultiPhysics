// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "ComMod.h"
#include "SimulationLogger.h"
#include "ZeroDBoundaryCondition.h"
#include "all_fun.h"
#include "consts.h"
#include "utils.h"
#include "VtkData.h"
#include "nn.h"
#include "fils_struct.hpp"
#include <optional>
#include <unordered_map>
#include <vtkCellType.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>
#include <fstream>


// Map VTK cell types to ElementType (same as vtk_xml_parser.cpp)
static consts::ElementType vtk_cell_type_to_element_type(int vtk_cell_type)
{
    using namespace consts;
    switch (vtk_cell_type) {
        case VTK_TRIANGLE: return ElementType::TRI3;
        case VTK_QUADRATIC_TRIANGLE: return ElementType::TRI6;
        case VTK_QUAD: return ElementType::QUD4;
        case VTK_QUADRATIC_QUAD: return ElementType::QUD8;
        case VTK_BIQUADRATIC_QUAD: return ElementType::QUD9;
        case VTK_LINE: return ElementType::LIN1;
        case VTK_TETRA: return ElementType::TET4;
        case VTK_QUADRATIC_TETRA: return ElementType::TET10;
        case VTK_HEXAHEDRON: return ElementType::HEX8;
        case VTK_QUADRATIC_HEXAHEDRON: return ElementType::HEX20;
        case VTK_TRIQUADRATIC_HEXAHEDRON: return ElementType::HEX27;
        case VTK_WEDGE: return ElementType::WDG;
        default:
            throw std::runtime_error("[ZeroDBoundaryCondition] Unsupported VTK cell type " + 
                                    std::to_string(vtk_cell_type) + " in cap VTP file.");
    }
}

static std::string bc_type_to_string(consts::BoundaryConditionType type)
{
    switch (type) {
        case consts::BoundaryConditionType::bType_Dir: return "Dirichlet";
        case consts::BoundaryConditionType::bType_Neu: return "Neumann";
        default: return "Unsupported";
    }
}

// Helper functions to access cell data from VTP file
namespace {
    /// @brief Check if a cell data array exists in a VTP file
    bool has_cell_data_in_vtp(const std::string& vtp_file_path, const std::string& data_name)
    {
        auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        reader->SetFileName(vtp_file_path.c_str());
        reader->Update();
        auto polydata = reader->GetOutput();
        
        if (polydata == nullptr) {
            return false;
        }
        
        int num_arrays = polydata->GetCellData()->GetNumberOfArrays();
        for (int i = 0; i < num_arrays; i++) {
            if (!strcmp(polydata->GetCellData()->GetArrayName(i), data_name.c_str())) {
                return true;
            }
        }
        
        return false;
    }
    
    /// @brief Get dimensions of a cell data array from a VTP file
    /// @return Pair of (num_components, num_tuples), or (0, 0) if not found
    std::pair<int, int> get_cell_data_dimensions_from_vtp(const std::string& vtp_file_path, const std::string& data_name)
    {
        auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        reader->SetFileName(vtp_file_path.c_str());
        reader->Update();
        auto polydata = reader->GetOutput();
        
        if (polydata == nullptr) {
            return std::make_pair(0, 0);
        }
        
        auto cell_data = polydata->GetCellData();
        if (cell_data == nullptr) {
            return std::make_pair(0, 0);
        }
        
        auto vtk_array = cell_data->GetArray(data_name.c_str());
        if (vtk_array == nullptr) {
            return std::make_pair(0, 0);
        }
        
        // Try to cast to vtkDoubleArray first
        auto vtk_data = vtkDoubleArray::SafeDownCast(vtk_array);
        if (vtk_data == nullptr) {
            // If not double, try to get dimensions from the base array
            // This handles cases where the array might be float or another numeric type
            int num_comp = vtk_array->GetNumberOfComponents();
            int num_tuples = vtk_array->GetNumberOfTuples();
            return std::make_pair(num_comp, num_tuples);
        }
        
        int num_comp = vtk_data->GetNumberOfComponents();
        int num_tuples = vtk_data->GetNumberOfTuples();
        
        return std::make_pair(num_comp, num_tuples);
    }
    
    /// @brief Copy cell data from a VTP file to an Array
    void copy_cell_data_from_vtp(const std::string& vtp_file_path, const std::string& data_name, Array<double>& mesh_data)
    {
        auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        reader->SetFileName(vtp_file_path.c_str());
        reader->Update();
        auto polydata = reader->GetOutput();
        
        if (polydata == nullptr) {
            return;
        }
        
        auto cell_data = polydata->GetCellData();
        if (cell_data == nullptr) {
            return;
        }
        
        auto vtk_array = cell_data->GetArray(data_name.c_str());
        if (vtk_array == nullptr) {
            return;
        }
        
        int num_data = vtk_array->GetNumberOfTuples();
        if (num_data == 0) {
            return;
        }
        
        int num_comp = vtk_array->GetNumberOfComponents();
        
        // Set the data - handle both double and float arrays
        auto vtk_data = vtkDoubleArray::SafeDownCast(vtk_array);
        if (vtk_data != nullptr) {
            // Double array - direct copy
            for (int i = 0; i < num_data; i++) {
                auto tuple = vtk_data->GetTuple(i);
                for (int j = 0; j < num_comp; j++) {
                    mesh_data(j, i) = tuple[j];
                }
            }
        } else {
            // Other numeric type - convert to double
            for (int i = 0; i < num_data; i++) {
                for (int j = 0; j < num_comp; j++) {
                    mesh_data(j, i) = vtk_array->GetComponent(i, j);
                }
            }
        }
    }
}

ZeroDBoundaryCondition::ZeroDBoundaryCondition(const ZeroDBoundaryCondition& other)
    : face_(other.face_)
    , cap_face_vtp_file_(other.cap_face_vtp_file_)
    , logger_(other.logger_)
    , bc_type_(other.bc_type_)
    , block_name_(other.block_name_)
    , face_name_(other.face_name_)
    , Qo_(other.Qo_)
    , Qn_(other.Qn_)
    , Po_(other.Po_)
    , Pn_(other.Pn_)
    , pressure_(other.pressure_)
    , flow_sol_id_(other.flow_sol_id_)
    , pressure_sol_id_(other.pressure_sol_id_)
    , in_out_sign_(other.in_out_sign_)
    , follower_pressure_load_(other.follower_pressure_load_)
    , cap_loaded_(other.cap_loaded_)
    , cap_global_node_ids_(other.cap_global_node_ids_)
    , cap_area_computed_(other.cap_area_computed_)
    , cap_gnNo_to_tnNo_(other.cap_gnNo_to_tnNo_)
    , cap_valM_(other.cap_valM_)
    , cap_initial_normals_(other.cap_initial_normals_)
{
    // Deep copy the cap_face_ unique_ptr if it exists
    // Note: faceType uses default copy constructor, which should properly copy all arrays
    if (other.cap_face_ != nullptr) {
        try {
            cap_face_ = std::make_unique<faceType>(*other.cap_face_);
            // Validate the copied faceType
            if (cap_face_ != nullptr) {
                // Basic validation - ensure sizes are consistent
                if (cap_face_->nNo > 0 && cap_face_->gN.size() != cap_face_->nNo) {
                    throw std::runtime_error("[ZeroDBoundaryCondition::copy constructor] Invalid cap_face_ after copy: gN.size()=" +
                                            std::to_string(cap_face_->gN.size()) + " != nNo=" + std::to_string(cap_face_->nNo));
                }
                if (cap_face_->nEl > 0 && cap_face_->IEN.ncols() != cap_face_->nEl) {
                    throw std::runtime_error("[ZeroDBoundaryCondition::copy constructor] Invalid cap_face_ after copy: IEN.ncols()=" +
                                            std::to_string(cap_face_->IEN.ncols()) + " != nEl=" + std::to_string(cap_face_->nEl));
                }
            }
        } catch (const std::exception& e) {
            throw std::runtime_error("[ZeroDBoundaryCondition::copy constructor] Failed to copy cap_face_: " + std::string(e.what()));
        }
    } else {
        cap_face_.reset();
    }
}

ZeroDBoundaryCondition& ZeroDBoundaryCondition::operator=(const ZeroDBoundaryCondition& other)
{
    if (this != &other) {
        face_ = other.face_;
        cap_face_vtp_file_ = other.cap_face_vtp_file_;
        logger_ = other.logger_;
        bc_type_ = other.bc_type_;
        block_name_ = other.block_name_;
        face_name_ = other.face_name_;
        Qo_ = other.Qo_;
        Qn_ = other.Qn_;
        Po_ = other.Po_;
        Pn_ = other.Pn_;
        pressure_ = other.pressure_;
        flow_sol_id_ = other.flow_sol_id_;
        pressure_sol_id_ = other.pressure_sol_id_;
        in_out_sign_ = other.in_out_sign_;
        follower_pressure_load_ = other.follower_pressure_load_;
        cap_loaded_ = other.cap_loaded_;
        cap_global_node_ids_ = other.cap_global_node_ids_;
    cap_area_computed_ = other.cap_area_computed_;
    cap_gnNo_to_tnNo_ = other.cap_gnNo_to_tnNo_;
    cap_valM_ = other.cap_valM_;
    cap_initial_normals_ = other.cap_initial_normals_;
        
        // Deep copy the cap_face_ unique_ptr if it exists
        // Note: faceType uses default copy constructor, which should properly copy all arrays
        if (other.cap_face_ != nullptr) {
            try {
                cap_face_ = std::make_unique<faceType>(*other.cap_face_);
            } catch (const std::exception& e) {
                throw std::runtime_error("[ZeroDBoundaryCondition::operator=] Failed to copy cap_face_: " + std::string(e.what()));
            }
        } else {
            cap_face_.reset();
        }
    }
    return *this;
}

ZeroDBoundaryCondition::ZeroDBoundaryCondition(ZeroDBoundaryCondition&& other) noexcept
    : face_(other.face_)
    , cap_face_vtp_file_(std::move(other.cap_face_vtp_file_))
    , logger_(other.logger_)
    , bc_type_(other.bc_type_)
    , block_name_(std::move(other.block_name_))
    , face_name_(std::move(other.face_name_))
    , Qo_(other.Qo_)
    , Qn_(other.Qn_)
    , Po_(other.Po_)
    , Pn_(other.Pn_)
    , pressure_(other.pressure_)
    , flow_sol_id_(other.flow_sol_id_)
    , pressure_sol_id_(other.pressure_sol_id_)
    , in_out_sign_(other.in_out_sign_)
    , follower_pressure_load_(other.follower_pressure_load_)
    , cap_loaded_(other.cap_loaded_)
    , cap_face_(std::move(other.cap_face_))
    , cap_global_node_ids_(std::move(other.cap_global_node_ids_))
    , cap_area_computed_(other.cap_area_computed_)
    , cap_gnNo_to_tnNo_(std::move(other.cap_gnNo_to_tnNo_))
    , cap_valM_(std::move(other.cap_valM_))
    , cap_initial_normals_(std::move(other.cap_initial_normals_))
{
    // Reset moved-from object to valid state
    other.face_ = nullptr;
    other.logger_ = nullptr;
    other.cap_loaded_ = false;
    other.cap_area_computed_ = false;
    other.flow_sol_id_ = -1;
    other.pressure_sol_id_ = -1;
    other.Qo_ = 0.0;
    other.Qn_ = 0.0;
    other.Po_ = 0.0;
    other.Pn_ = 0.0;
    other.pressure_ = 0.0;
}

ZeroDBoundaryCondition& ZeroDBoundaryCondition::operator=(ZeroDBoundaryCondition&& other) noexcept
{
    if (this != &other) {
        face_ = other.face_;
        cap_face_vtp_file_ = std::move(other.cap_face_vtp_file_);
        logger_ = other.logger_;
        bc_type_ = other.bc_type_;
        block_name_ = std::move(other.block_name_);
        face_name_ = std::move(other.face_name_);
        Qo_ = other.Qo_;
        Qn_ = other.Qn_;
        Po_ = other.Po_;
        Pn_ = other.Pn_;
        pressure_ = other.pressure_;
        flow_sol_id_ = other.flow_sol_id_;
        pressure_sol_id_ = other.pressure_sol_id_;
        in_out_sign_ = other.in_out_sign_;
        follower_pressure_load_ = other.follower_pressure_load_;
        cap_loaded_ = other.cap_loaded_;
        cap_face_ = std::move(other.cap_face_);
        cap_global_node_ids_ = std::move(other.cap_global_node_ids_);
    cap_area_computed_ = other.cap_area_computed_;
    cap_gnNo_to_tnNo_ = std::move(other.cap_gnNo_to_tnNo_);
    cap_valM_ = std::move(other.cap_valM_);
    cap_initial_normals_ = std::move(other.cap_initial_normals_);
        
        // Reset moved-from object to valid state
        other.face_ = nullptr;
        other.logger_ = nullptr;
        other.cap_loaded_ = false;
        other.cap_area_computed_ = false;
        other.flow_sol_id_ = -1;
        other.pressure_sol_id_ = -1;
        other.Qo_ = 0.0;
        other.Qn_ = 0.0;
        other.Po_ = 0.0;
        other.Pn_ = 0.0;
        other.pressure_ = 0.0;
    }
    return *this;
}

ZeroDBoundaryCondition::ZeroDBoundaryCondition(consts::BoundaryConditionType bc_type, const faceType& face, const std::string& face_name,
                                               const std::string& block_name, SimulationLogger& logger)
    : face_(&face)
    , logger_(&logger)
    , bc_type_(bc_type)
    , block_name_(block_name)
    , face_name_(face_name)
{
    if (bc_type_ != consts::BoundaryConditionType::bType_Dir &&
        bc_type_ != consts::BoundaryConditionType::bType_Neu) {
        throw std::runtime_error("[ZeroDBoundaryCondition] bc_type must be bType_Dir or bType_Neu.");
    }
}

ZeroDBoundaryCondition::ZeroDBoundaryCondition(consts::BoundaryConditionType bc_type, const faceType& face, const std::string& face_name,
                                               const std::string& block_name, const std::string& cap_face_vtp_file, SimulationLogger& logger)
    : cap_face_vtp_file_(cap_face_vtp_file)
    , face_(&face)
    , logger_(&logger)
    , bc_type_(bc_type)
    , block_name_(block_name)
    , face_name_(face_name)
{
    if (bc_type_ != consts::BoundaryConditionType::bType_Dir &&
        bc_type_ != consts::BoundaryConditionType::bType_Neu) {
        throw std::runtime_error("[ZeroDBoundaryCondition] bc_type must be bType_Dir or bType_Neu.");
    }
    // Load the cap VTP file if provided
    if (!cap_face_vtp_file_.empty()) {
        load_cap_face_vtp(cap_face_vtp_file_);
    }
}

void ZeroDBoundaryCondition::load_cap_face_vtp(const std::string& vtp_file_path)
{
    // Safety check: ensure face_ is set before loading cap
    if (face_ == nullptr) {
        throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Cannot load cap: face_ is null. Call set_face() first.");
    }
    
    cap_face_vtp_file_ = vtp_file_path;
    cap_loaded_ = false;
    cap_face_.reset();
    cap_gnNo_to_tnNo_.clear();  // Clear the mapping since we're reloading
    cap_area_computed_ = false;  // Reset area computation flag
    cap_initial_normals_.resize(0, 0);  // Clear initial normals
    
    if (vtp_file_path.empty()) {
        return;
    }
    
    // Check if file exists and is readable
    std::ifstream file_check(vtp_file_path);
    if (!file_check.good()) {
        throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Cannot open cap VTP file '" + vtp_file_path + "' for reading.");
    }
    file_check.close();
    
    
    // Load the VTP file - wrap in try-catch to handle any exceptions
    // Note: VtkVtpData constructor with file_name calls read_file() which may throw or crash
    VtkVtpData vtp_data;
    try {
        // Construct VtkVtpData and read file in one step
        // The constructor with file_name will call read_file() internally
        vtp_data = VtkVtpData(vtp_file_path, true);
    } catch (const std::exception& e) {
        throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Failed to construct VtkVtpData from file '" + 
                                vtp_file_path + "': " + e.what());
    } catch (...) {
        throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Unknown error constructing VtkVtpData from file '" + 
                                vtp_file_path + "'. This may indicate a crash in the VTK library.");
    }
    
    int nNo = 0;
    try {
        nNo = vtp_data.num_points();
    } catch (const std::exception& e) {
        throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Failed to get number of points from VTP file '" + 
                                vtp_file_path + "': " + e.what());
    } catch (...) {
        throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Unknown error getting number of points from file '" + 
                                vtp_file_path + "'");
    }
    if (nNo == 0) {
        throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Cap VTP file '" + vtp_file_path + "' does not contain any points.");
    }
    
    int num_elems = 0;
    try {
        num_elems = vtp_data.num_elems();
    } catch (const std::exception& e) {
        throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Failed to get number of elements from VTP file '" + 
                                vtp_file_path + "': " + e.what());
    } catch (...) {
        throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Unknown error getting number of elements from file '" + 
                                vtp_file_path + "'");
    }
    if (num_elems == 0) {
        throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Cap VTP file '" + vtp_file_path + "' does not contain any elements.");
    }
    
    // Get GlobalNodeID from VTP file
    bool has_global_node_id = false;
    try {
        has_global_node_id = vtp_data.has_point_data("GlobalNodeID");
    } catch (const std::exception& e) {
        throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Failed to check for GlobalNodeID in VTP file '" + 
                                vtp_file_path + "': " + e.what());
    } catch (...) {
        throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Unknown error checking for GlobalNodeID in file '" + 
                                vtp_file_path + "'");
    }
    if (!has_global_node_id) {
        throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Cap VTP file '" + vtp_file_path + "' does not contain 'GlobalNodeID' point data.");
    }
    
    try {
        cap_global_node_ids_.resize(nNo);
        vtp_data.copy_point_data("GlobalNodeID", cap_global_node_ids_);
    } catch (const std::exception& e) {
        throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Failed to copy GlobalNodeID from VTP file '" + 
                                vtp_file_path + "': " + e.what());
    } catch (...) {
        throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Unknown error copying GlobalNodeID from file '" + 
                                vtp_file_path + "'");
    }
    
    // Get connectivity from VTP file
    Array<int> conn;
    int eNoN = 0;
    int vtk_cell_type = 0;
    try {
        conn = vtp_data.get_connectivity();
        eNoN = vtp_data.np_elem();
        vtk_cell_type = vtp_data.elem_type();  // This returns VTK cell type (e.g., 5 for VTK_TRIANGLE)
    } catch (const std::exception& e) {
        throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Failed to get connectivity from VTP file '" + 
                                vtp_file_path + "': " + e.what());
    } catch (...) {
        throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Unknown error getting connectivity from file '" + 
                                vtp_file_path + "'");
    }
    if (eNoN <= 0) {
        throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Invalid number of nodes per element: " + std::to_string(eNoN));
    }
    
    // Convert VTK cell type to ElementType enum
    consts::ElementType eType = vtk_cell_type_to_element_type(vtk_cell_type);
    
    // Create a temporary faceType structure for the cap
    cap_face_ = std::make_unique<faceType>();
    cap_face_->name = face_name_ + "_cap";
    cap_face_->iM = face_->iM;  // Use the same mesh index as the main face
    cap_face_->nNo = nNo;
    cap_face_->nEl = num_elems;
    cap_face_->gnEl = num_elems;
    cap_face_->eNoN = eNoN;
    cap_face_->eType = eType;
    
    // Set global node IDs (map GlobalNodeID from VTP to mesh global node indices)
    // Verify cap_global_node_ids_ is properly sized
    if (cap_global_node_ids_.size() != nNo) {
        throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] cap_global_node_ids_ size mismatch: " +
                                std::to_string(cap_global_node_ids_.size()) + " != " + std::to_string(nNo));
    }
    cap_face_->gN.resize(nNo);
    for (int a = 0; a < nNo; a++) {
        // Bounds check before accessing arrays
        if (a < 0 || a >= cap_global_node_ids_.size()) {
            throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Index out of bounds for cap_global_node_ids_: " +
                                    std::to_string(a) + " (size=" + std::to_string(cap_global_node_ids_.size()) + ")");
        }
        if (a < 0 || a >= cap_face_->gN.size()) {
            throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Index out of bounds for gN: " +
                                    std::to_string(a) + " (size=" + std::to_string(cap_face_->gN.size()) + ")");
        }
        cap_face_->gN(a) = cap_global_node_ids_(a) - 1;  // Convert from 1-based to 0-based
    }
    
    // Map connectivity from local cap node indices to global mesh node indices
    // The VTP connectivity contains local node indices (0 to nNo-1)
    // We need to map these to global mesh node indices using gN
    cap_face_->IEN.resize(eNoN, num_elems);
    
    // Verify gN is properly sized before using it
    if (cap_face_->gN.size() != nNo) {
        throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] cap_face_->gN size mismatch: " +
                                std::to_string(cap_face_->gN.size()) + " != " + std::to_string(nNo));
    }
    
    for (int e = 0; e < num_elems; e++) {
        // Bounds check for element index
        if (e < 0 || e >= conn.ncols()) {
            throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Invalid element index " + 
                                    std::to_string(e) + " (conn.ncols()=" + std::to_string(conn.ncols()) + ")");
        }
        for (int a = 0; a < eNoN; a++) {
            // Bounds check for node index in connectivity
            if (a < 0 || a >= conn.nrows()) {
                throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Invalid node index in connectivity " + 
                                        std::to_string(a) + " (conn.nrows()=" + std::to_string(conn.nrows()) + ")");
            }
            int local_node_idx = conn(a, e);  // Local node index in cap (0 to nNo-1)
            if (local_node_idx < 0 || local_node_idx >= nNo) {
                throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Invalid local node index " + 
                                        std::to_string(local_node_idx) + " in cap connectivity (element " + 
                                        std::to_string(e) + ", node " + std::to_string(a) + ", nNo=" + std::to_string(nNo) + ").");
            }
            // Bounds check for gN access
            if (local_node_idx < 0 || local_node_idx >= cap_face_->gN.size()) {
                throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Index out of bounds for gN access: " + 
                                        std::to_string(local_node_idx) + " (gN.size()=" + std::to_string(cap_face_->gN.size()) + ")");
            }
            cap_face_->IEN(a, e) = cap_face_->gN(local_node_idx);  // Map to global mesh node index
        }
    }
    
    // Load initial element normals from VTP file if available
    // The normals are stored as cell data array called "Normals"
    try {
        bool has_normals = has_cell_data_in_vtp(vtp_file_path, "Normals");
        if (has_normals) {
            // Get dimensions from VTK array
            auto [num_comp, num_tuples] = get_cell_data_dimensions_from_vtp(vtp_file_path, "Normals");
            if (num_comp == 0 || num_tuples == 0) {
                // Provide more detailed error message
                std::string error_msg = "[ZeroDBoundaryCondition::load_cap_face_vtp] Normals array exists but has zero size. ";
                error_msg += "num_components=" + std::to_string(num_comp) + ", num_tuples=" + std::to_string(num_tuples);
                error_msg += ", expected num_elems=" + std::to_string(num_elems);
                error_msg += ". The array may not be a numeric type or may be empty.";
                throw std::runtime_error(error_msg);
            }
            
            if (num_tuples != num_elems) {
                throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Normals array size mismatch: " +
                                        std::to_string(num_tuples) + " != " + std::to_string(num_elems));
            }
            
            if (num_comp != 2 && num_comp != 3) {
                throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Invalid number of components in Normals array: " +
                                        std::to_string(num_comp) + " (expected 2 or 3)");
            }
            
            // Resize array to correct dimensions
            cap_initial_normals_.resize(num_comp, num_elems);
            copy_cell_data_from_vtp(vtp_file_path, "Normals", cap_initial_normals_);
            
            // Verify that data was actually copied
            if (cap_initial_normals_.nrows() != num_comp || cap_initial_normals_.ncols() != num_elems) {
                throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Failed to copy Normals data. "
                                        "Expected size: " + std::to_string(num_comp) + "x" + std::to_string(num_elems) +
                                        ", Actual size: " + std::to_string(cap_initial_normals_.nrows()) + "x" + 
                                        std::to_string(cap_initial_normals_.ncols()));
            }
        } else {
            // Normals not found - initialize to empty array
            cap_initial_normals_.resize(0, 0);
        }
    } catch (const std::exception& e) {
        throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Failed to load Normals from VTP file '" + 
                                vtp_file_path + "': " + e.what());
    } catch (...) {
        throw std::runtime_error("[ZeroDBoundaryCondition::load_cap_face_vtp] Unknown error loading Normals from file '" + 
                                vtp_file_path + "'");
    }
    
    cap_loaded_ = true;
    
}

// =========================================================================
// svZeroD block configuration
// =========================================================================

void ZeroDBoundaryCondition::set_block_name(const std::string& block_name)
{
    block_name_ = block_name;
}

const std::string& ZeroDBoundaryCondition::get_block_name() const
{
    return block_name_;
}

void ZeroDBoundaryCondition::set_face_name(const std::string& face_name)
{
    face_name_ = face_name;
}

const std::string& ZeroDBoundaryCondition::get_face_name() const
{
    return face_name_;
}

void ZeroDBoundaryCondition::set_solution_ids(int flow_id, int pressure_id, double in_out_sign)
{
    flow_sol_id_ = flow_id;
    pressure_sol_id_ = pressure_id;
    in_out_sign_ = in_out_sign;
}

int ZeroDBoundaryCondition::get_flow_sol_id() const
{
    return flow_sol_id_;
}

int ZeroDBoundaryCondition::get_pressure_sol_id() const
{
    return pressure_sol_id_;
}

double ZeroDBoundaryCondition::get_in_out_sign() const
{
    return in_out_sign_;
}

// =========================================================================
// Flowrate computation and access
// =========================================================================

void ZeroDBoundaryCondition::set_follower_pressure_load(bool flwP)
{
    follower_pressure_load_ = flwP;
}

bool ZeroDBoundaryCondition::get_follower_pressure_load() const
{
    return follower_pressure_load_;
}

/// @brief Compute flowrates at the boundary face at old and new timesteps
///
/// This replicates the flowrate computation done in set_bc::calc_der_cpl_bc and
/// set_bc::set_bc_cpl for coupled Neumann boundary conditions.
///
/// The flowrate is computed as the integral of velocity dotted with the face normal.
/// For struct/ustruct physics, the integral is computed on the deformed configuration.
/// For fluid/FSI/CMM physics, the integral is computed on the reference configuration.
void ZeroDBoundaryCondition::compute_flowrates(ComMod& com_mod, const CmMod& cm_mod, consts::EquationType phys)
{
    using namespace consts;
    
    
    if (face_ == nullptr) {
        throw std::runtime_error("[ZeroDBoundaryCondition::compute_flowrates] Face is not set.");
    }
    
    int nsd = com_mod.nsd;
    
    // Determine mechanical configuration based on physics type
    MechanicalConfigurationType cfg_o = MechanicalConfigurationType::reference;
    MechanicalConfigurationType cfg_n = MechanicalConfigurationType::reference;
    
    // For struct or ustruct, use old and new configurations to compute flowrate integral
    if ((phys == EquationType::phys_struct) || (phys == EquationType::phys_ustruct)) {
        // Must use follower pressure load for 0D coupling with struct/ustruct
        if (!follower_pressure_load_) {
            throw std::runtime_error("[ZeroDBoundaryCondition::compute_flowrates] Follower pressure load must be used for 0D coupling with struct/ustruct");
        }
        cfg_o = MechanicalConfigurationType::old_timestep;
        cfg_n = MechanicalConfigurationType::new_timestep;
    }
    // For fluid, FSI, or CMM, use reference configuration to compute flowrate integral
    else if ((phys == EquationType::phys_fluid) || (phys == EquationType::phys_FSI) || (phys == EquationType::phys_CMM)) {
        cfg_o = MechanicalConfigurationType::reference;
        cfg_n = MechanicalConfigurationType::reference;
    }
    else {
        throw std::runtime_error("[ZeroDBoundaryCondition::compute_flowrates] Invalid physics type for 0D coupling");
    }
    
    // Compute flowrates by integrating velocity over face
    // The all_fun::integ function with indices 0 to nsd-1 integrates the velocity vector
    // dotted with the face normal, giving the volumetric flowrate
    Qo_ = all_fun::integ(com_mod, cm_mod, *face_, com_mod.Yo, 0, nsd-1, false, cfg_o);
    Qn_ = all_fun::integ(com_mod, cm_mod, *face_, com_mod.Yn, 0, nsd-1, false, cfg_n);
    
    // Add cap contribution if cap exists
    if (has_cap()) {
        double Qo_cap = integrate_over_cap(com_mod, cm_mod, com_mod.Yo, 0, nsd-1, cfg_o);
        double Qn_cap = integrate_over_cap(com_mod, cm_mod, com_mod.Yn, 0, nsd-1, cfg_n);
        Qo_ += Qo_cap;
        Qn_ += Qn_cap;
    }
    
}

/// @brief Compute average pressures at the boundary face at old and new timesteps
///
/// This replicates the pressure computation done in set_bc::calc_der_cpl_bc and
/// set_bc::set_bc_cpl for coupled Dirichlet boundary conditions.
///
/// The pressure is computed as the average pressure over the face by integrating
/// pressure (at index nsd in the solution vector) and dividing by the face area.
void ZeroDBoundaryCondition::compute_pressures(ComMod& com_mod, const CmMod& cm_mod)
{
    using namespace consts;
    
    
    if (face_ == nullptr) {
        throw std::runtime_error("[ZeroDBoundaryCondition::compute_pressures] Face is not set.");
    }
    
    int nsd = com_mod.nsd;
    double area = face_->area;
    
    Po_ = all_fun::integ(com_mod, cm_mod, *face_, com_mod.Yo, nsd) / area;
    Pn_ = all_fun::integ(com_mod, cm_mod, *face_, com_mod.Yn, nsd) / area;
    
}

double ZeroDBoundaryCondition::get_Qo() const
{
    return Qo_;
}

double ZeroDBoundaryCondition::get_Qn() const
{
    return Qn_;
}

void ZeroDBoundaryCondition::set_flowrates(double Qo, double Qn)
{
    Qo_ = Qo;
    Qn_ = Qn;
}

void ZeroDBoundaryCondition::perturb_flowrate(double diff)
{
    Qn_ += diff;
}

// =========================================================================
// Pressure access
// =========================================================================

void ZeroDBoundaryCondition::set_pressure(double pressure)
{
    pressure_ = pressure;
}

double ZeroDBoundaryCondition::get_pressure() const
{
    return pressure_;
}

double ZeroDBoundaryCondition::get_Po() const
{
    return Po_;
}

double ZeroDBoundaryCondition::get_Pn() const
{
    return Pn_;
}

// =========================================================================
// State management
// =========================================================================

ZeroDBoundaryCondition::State ZeroDBoundaryCondition::save_state() const
{
    return State{Qn_, pressure_};
}

void ZeroDBoundaryCondition::restore_state(const State& state)
{
    Qn_ = state.Qn;
    pressure_ = state.pressure;
}

// =========================================================================
// Utility methods
// =========================================================================

void ZeroDBoundaryCondition::distribute(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, const faceType& face)
{
    #define n_debug_zerod_distribute

    // Update face pointer to local face
    face_ = &face;

    const bool is_slave = cm.slv(cm_mod);
    
    // Distribute BC type (Dirichlet or Neumann)
    int bc_type_int = static_cast<int>(bc_type_);
    cm.bcast(cm_mod, &bc_type_int);
    if (is_slave) {
        bc_type_ = static_cast<consts::BoundaryConditionType>(bc_type_int);
    }
    
    // Distribute block name
    cm.bcast(cm_mod, block_name_);
    
    // Distribute face name
    cm.bcast(cm_mod, face_name_);
    
    // Distribute follower pressure load flag
    cm.bcast(cm_mod, &follower_pressure_load_);
    
    // Distribute solution IDs
    cm.bcast(cm_mod, &flow_sol_id_);
    cm.bcast(cm_mod, &pressure_sol_id_);
    cm.bcast(cm_mod, &in_out_sign_);
    
    // Distribute cap VTP file path
    cm.bcast(cm_mod, cap_face_vtp_file_);
    
    // Load cap on all processes if not already loaded
    // Note: This must happen after face_ is set (line 513) because load_cap_face_vtp uses face_->iM
    // On slave processes, cap won't be loaded yet, so we need to load it
    // On master process, cap is already loaded in constructor, but we reload to ensure consistency
    if (!cap_face_vtp_file_.empty()) {
        if (face_ == nullptr) {
            throw std::runtime_error("[ZeroDBoundaryCondition::distribute] Cannot load cap: face_ is null.");
        }
        // Always reload on slave processes, reload on master to ensure consistency after face_ update
        try {
            load_cap_face_vtp(cap_face_vtp_file_);
        } catch (const std::exception& e) {
            throw std::runtime_error("[ZeroDBoundaryCondition::distribute] Failed to load cap VTP file '" + 
                                    cap_face_vtp_file_ + "' on process " + std::to_string(cm.idcm()) + ": " + e.what());
        }
    } else {
        // If cap file path is empty, clear any existing cap
        cap_loaded_ = false;
        cap_face_.reset();
        cap_gnNo_to_tnNo_.clear();
        cap_area_computed_ = false;
        cap_initial_normals_.resize(0, 0);
    }
    
    // Note: Cap integration initialization happens in initialize() after distribute()
    // sets up com_mod.ltg for each process. The gnNo_to_tnNo mapping will be rebuilt in initialize().
    
}

void ZeroDBoundaryCondition::set_face(const faceType& face)
{
    face_ = &face;
}

const faceType* ZeroDBoundaryCondition::get_face() const
{
    return face_;
}

bool ZeroDBoundaryCondition::is_initialized() const
{
    return (face_ != nullptr);
}

// =========================================================================
// Cap surface integration
// =========================================================================

void ZeroDBoundaryCondition::initialize_cap_integration(ComMod& com_mod, const CmMod& cm_mod)
{
    using namespace consts;
    
    if (!has_cap()) {
        return;
    }
    
    // Safety check: ensure cap_face_ is properly initialized
    if (cap_face_ == nullptr || cap_face_->nEl == 0 || cap_face_->nNo == 0) {
        throw std::runtime_error("[ZeroDBoundaryCondition::initialize_cap_integration] Cap face is not properly initialized.");
    }
    
    // Skip if already initialized
    if (cap_face_->nG > 0 && cap_face_->w.size() > 0 && cap_face_->N.nrows() > 0 && 
        cap_face_->Nx.nslices() == cap_face_->nG && !cap_gnNo_to_tnNo_.empty()) {
        return;
    }
    
    int nsd = com_mod.nsd;
    int insd = nsd - 1;
    
    
    // Build reverse map from gnNo to tnNo using com_mod.ltg
    // com_mod.ltg maps from tnNo to gnNo, so we reverse it
    // IEN contains gnNo indices from load_cap_face_vtp, we'll map them on-the-fly during integration
    // Note: com_mod.ltg might not be set up yet when called from read_files.cpp.
    // If it's not ready, we'll skip initialization and it will be done later when needed.
    auto& msh = com_mod.msh[face_->iM];
    if (com_mod.ltg.size() == 0 || com_mod.tnNo == 0) {
        // ltg not ready yet - will be initialized later when ltg is set up
        return;
    }
    if (com_mod.ltg.size() != com_mod.tnNo) {
        throw std::runtime_error("[ZeroDBoundaryCondition::initialize_cap_integration] com_mod.ltg size (" + 
                                std::to_string(com_mod.ltg.size()) + ") does not match com_mod.tnNo (" + 
                                std::to_string(com_mod.tnNo) + ")");
    }
    
    cap_gnNo_to_tnNo_.clear();
    cap_gnNo_to_tnNo_.reserve(com_mod.tnNo);  // Reserve space for better performance
    for (int a = 0; a < com_mod.tnNo; a++) {
        int gnNo_idx = com_mod.ltg(a);
        // Check for valid gnNo index
        if (gnNo_idx < 0 || gnNo_idx >= msh.gnNo) {
            throw std::runtime_error("[ZeroDBoundaryCondition::initialize_cap_integration] Invalid gnNo index " + 
                                    std::to_string(gnNo_idx) + " in com_mod.ltg at position " + 
                                    std::to_string(a) + " (msh.gnNo=" + std::to_string(msh.gnNo) + ")");
        }
        cap_gnNo_to_tnNo_[gnNo_idx] = a;
    }
    
    try {
        // Determine nG based on element type (similar to select_eleb)
        // For common surface elements:
        if (cap_face_->eType == ElementType::TRI3) {
            // Use 1-point centroid rule for TRI3 (more efficient than 3-point rule)
            cap_face_->nG = 1;
        } else if (cap_face_->eType == ElementType::QUD4) {
            cap_face_->nG = 4;  // 4 Gauss points for QUD4
        } else if (cap_face_->eType == ElementType::TRI6) {
            cap_face_->nG = 3;  // 3 Gauss points for TRI6
        } else {
            // Default: use 1 Gauss point per element node
            cap_face_->nG = cap_face_->eNoN;
        }
        
        // Allocate and initialize Gauss points and weights
        cap_face_->w.resize(cap_face_->nG);
        cap_face_->xi.resize(insd, cap_face_->nG);
        
        // For TRI3 with 1 point, manually set up centroid rule (get_gip is hardcoded for 3 points)
        if (cap_face_->eType == ElementType::TRI3 && cap_face_->nG == 1) {
            cap_face_->w(0) = 0.5;  // Weight for centroid rule on reference triangle
            cap_face_->xi(0, 0) = 1.0/3.0;  // Centroid coordinate
            cap_face_->xi(1, 0) = 1.0/3.0;  // Centroid coordinate
        } else {
            // For other element types, use get_gip
            nn::get_gip(insd, cap_face_->eType, cap_face_->nG, cap_face_->w, cap_face_->xi);
        }
        
        // Allocate and compute shape functions
        cap_face_->N.resize(cap_face_->eNoN, cap_face_->nG);
        cap_face_->Nx.resize(insd, cap_face_->eNoN, cap_face_->nG);
        for (int g = 0; g < cap_face_->nG; g++) {
            nn::get_gnn(insd, cap_face_->eType, cap_face_->eNoN, g, cap_face_->xi, cap_face_->N, cap_face_->Nx);
        }
    } catch (const std::exception& e) {
        throw std::runtime_error("[ZeroDBoundaryCondition::initialize_cap_integration] Failed to initialize cap face shape functions: " + 
                                std::string(e.what()));
    }
}

// =========================================================================
// Cap integration helper functions
// =========================================================================

Array<double> ZeroDBoundaryCondition::update_cap_element_position(ComMod& com_mod, int e, 
                                                                  const std::unordered_map<int, int>& gnNo_to_tnNo,
                                                                  consts::MechanicalConfigurationType cfg)
{
    using namespace consts;
    
    // Safety checks
    if (cap_face_ == nullptr) {
        throw std::runtime_error("[ZeroDBoundaryCondition::update_cap_element_position] cap_face_ is null.");
    }
    if (cap_face_->IEN.nrows() == 0 || cap_face_->IEN.ncols() == 0) {
        throw std::runtime_error("[ZeroDBoundaryCondition::update_cap_element_position] cap_face_->IEN is not allocated.");
    }
    if (e < 0 || e >= cap_face_->IEN.ncols()) {
        throw std::runtime_error("[ZeroDBoundaryCondition::update_cap_element_position] Element index e=" + 
                                std::to_string(e) + " is out of bounds (IEN.ncols()=" + std::to_string(cap_face_->IEN.ncols()) + ").");
    }
    if (cap_face_->eNoN <= 0) {
        throw std::runtime_error("[ZeroDBoundaryCondition::update_cap_element_position] cap_face_->eNoN is invalid: " + 
                                std::to_string(cap_face_->eNoN));
    }
    
    int nsd = com_mod.nsd;
    Array<double> xl(nsd, cap_face_->eNoN);
    
    for (int a = 0; a < cap_face_->eNoN; a++) {
        // Bounds check for IEN access
        if (a < 0 || a >= cap_face_->IEN.nrows()) {
            throw std::runtime_error("[ZeroDBoundaryCondition::update_cap_element_position] Node index a=" + 
                                    std::to_string(a) + " is out of bounds (IEN.nrows()=" + std::to_string(cap_face_->IEN.nrows()) + ").");
        }
        
        // IEN contains gnNo index, map it to tnNo index
        int gnNo_idx = cap_face_->IEN(a, e);
        auto it = gnNo_to_tnNo.find(gnNo_idx);
        if (it == gnNo_to_tnNo.end()) {
            throw std::runtime_error("[ZeroDBoundaryCondition::update_cap_element_position] IEN entry (element " + 
                                    std::to_string(e) + ", node " + std::to_string(a) + 
                                    ") contains invalid gnNo index " + std::to_string(gnNo_idx) + 
                                    " not found in com_mod.ltg mapping.");
        }
        int Ac = it->second;  // tnNo index
        
        // Bounds checking
        if (Ac < 0 || Ac >= com_mod.tnNo) {
            std::string msg = "[ZeroDBoundaryCondition::update_cap_element_position] Invalid node index Ac=" + 
                             std::to_string(Ac) + " (tnNo=" + std::to_string(com_mod.tnNo) + 
                             ") at element " + std::to_string(e) + ", local node " + std::to_string(a);
            throw std::runtime_error(msg);
        }
        
        // Bounds check for com_mod.x
        if (Ac >= com_mod.x.ncols() || nsd > com_mod.x.nrows()) {
            throw std::runtime_error("[ZeroDBoundaryCondition::update_cap_element_position] Invalid bounds for com_mod.x: Ac=" + 
                                    std::to_string(Ac) + ", x.ncols()=" + std::to_string(com_mod.x.ncols()) + 
                                    ", nsd=" + std::to_string(nsd) + ", x.nrows()=" + std::to_string(com_mod.x.nrows()));
        }
        
        xl.set_col(a, com_mod.x.col(Ac));
        
        // Apply displacement if needed based on configuration
        if (cfg == MechanicalConfigurationType::old_timestep) {
            // Bounds check for com_mod.Do
            if (Ac >= com_mod.Do.ncols() || nsd > com_mod.Do.nrows()) {
                throw std::runtime_error("[ZeroDBoundaryCondition::update_cap_element_position] Invalid bounds for com_mod.Do: Ac=" + 
                                        std::to_string(Ac) + ", Do.ncols()=" + std::to_string(com_mod.Do.ncols()) + 
                                        ", nsd=" + std::to_string(nsd) + ", Do.nrows()=" + std::to_string(com_mod.Do.nrows()));
            }
            for (int i = 0; i < nsd; i++) {
                xl(i, a) += com_mod.Do(i, Ac);
            }
        } else if (cfg == MechanicalConfigurationType::new_timestep) {
            // Bounds check for com_mod.Dn
            if (Ac >= com_mod.Dn.ncols() || nsd > com_mod.Dn.nrows()) {
                throw std::runtime_error("[ZeroDBoundaryCondition::update_cap_element_position] Invalid bounds for com_mod.Dn: Ac=" + 
                                        std::to_string(Ac) + ", Dn.ncols()=" + std::to_string(com_mod.Dn.ncols()) + 
                                        ", nsd=" + std::to_string(nsd) + ", Dn.nrows()=" + std::to_string(com_mod.Dn.nrows()));
            }
            for (int i = 0; i < nsd; i++) {
                xl(i, a) += com_mod.Dn(i, Ac);
            }
        }
    }
    
    return xl;
}

std::pair<double, Vector<double>> ZeroDBoundaryCondition::compute_cap_jacobian_and_normal(const Array<double>& xl, 
                                                                                           int e, int g, int nsd, int insd)
{
    // Safety checks
    if (cap_face_ == nullptr) {
        throw std::runtime_error("[ZeroDBoundaryCondition::compute_cap_jacobian_and_normal] cap_face_ is null.");
    }
    if (cap_face_->eNoN <= 0) {
        throw std::runtime_error("[ZeroDBoundaryCondition::compute_cap_jacobian_and_normal] cap_face_->eNoN is invalid: " + 
                                std::to_string(cap_face_->eNoN));
    }
    if (cap_face_->Nx.nslices() == 0) {
        throw std::runtime_error("[ZeroDBoundaryCondition::compute_cap_jacobian_and_normal] cap_face_->Nx is not allocated.");
    }
    if (e < 0 || e >= cap_face_->nEl) {
        throw std::runtime_error("[ZeroDBoundaryCondition::compute_cap_jacobian_and_normal] Element index e=" + 
                                std::to_string(e) + " is out of bounds (nEl=" + std::to_string(cap_face_->nEl) + ").");
    }
    if (g < 0 || g >= cap_face_->Nx.nslices()) {
        throw std::runtime_error("[ZeroDBoundaryCondition::compute_cap_jacobian_and_normal] Gauss point index g=" + 
                                std::to_string(g) + " is out of bounds (nG=" + std::to_string(cap_face_->Nx.nslices()) + ").");
    }
    if (cap_face_->Nx.nrows() != insd || cap_face_->Nx.ncols() != cap_face_->eNoN) {
        throw std::runtime_error("[ZeroDBoundaryCondition::compute_cap_jacobian_and_normal] cap_face_->Nx has wrong dimensions: " +
                                std::to_string(cap_face_->Nx.nrows()) + "x" + std::to_string(cap_face_->Nx.ncols()) + 
                                " (expected " + std::to_string(insd) + "x" + std::to_string(cap_face_->eNoN) + ").");
    }
    if (xl.nrows() != nsd || xl.ncols() != cap_face_->eNoN) {
        throw std::runtime_error("[ZeroDBoundaryCondition::compute_cap_jacobian_and_normal] xl has wrong dimensions: " +
                                std::to_string(xl.nrows()) + "x" + std::to_string(xl.ncols()) + 
                                " (expected " + std::to_string(nsd) + "x" + std::to_string(cap_face_->eNoN) + ").");
    }
    
    // Get shape function derivatives at this Gauss point
    Array<double> Nx_g = cap_face_->Nx.slice(g);  // (insd x eNoN)
    
    // Compute tangent vectors using basis functions
    // For 3D surface: compute ∂x/∂ξ₁ and ∂x/∂ξ₂
    // For 2D curve: compute ∂x/∂ξ
    Array<double> xXi(nsd, insd);  // Tangent vectors
    xXi = 0.0;
    
    for (int a = 0; a < cap_face_->eNoN; a++) {
        for (int i = 0; i < insd; i++) {
            for (int j = 0; j < nsd; j++) {
                xXi(j, i) += xl(j, a) * Nx_g(i, a);
            }
        }
    }
    
    // Compute Jacobian (area element) from tangent vectors
    double Jac = 0.0;
    Vector<double> n(nsd);
    if (nsd == 3 && insd == 2) {
        // 3D surface: Jacobian = ||∂x/∂ξ₁ × ∂x/∂ξ₂||
        // Verify xXi dimensions before calling cross
        if (xXi.nrows() != 3 || xXi.ncols() != 2) {
            throw std::runtime_error("[ZeroDBoundaryCondition::compute_cap_jacobian_and_normal] xXi has wrong dimensions: " +
                                    std::to_string(xXi.nrows()) + "x" + std::to_string(xXi.ncols()) + 
                                    " (expected 3x2 for 3D surface).");
        }
        n = utils::cross(xXi);
        Jac = sqrt(utils::norm(n));  // utils::norm returns squared norm, need sqrt
    } else if (nsd == 2 && insd == 1) {
        // 2D curve: Jacobian = ||∂x/∂ξ||
        Jac = sqrt(utils::norm(xXi.col(0)));  // utils::norm returns squared norm, need sqrt
        // For 2D, compute normal from tangent
        n(0) = -xXi(1, 0);
        n(1) = xXi(0, 0);
    } else {
        throw std::runtime_error("[ZeroDBoundaryCondition::compute_cap_jacobian_and_normal] Unsupported nsd/insd combination: " + 
                                std::to_string(nsd) + "/" + std::to_string(insd));
    }
    
    if (utils::is_zero(Jac)) {
        throw std::runtime_error("[ZeroDBoundaryCondition::compute_cap_jacobian_and_normal] Zero Jacobian at Gauss point " + 
                                std::to_string(g));
    }
    
    // Normalize normal vector
    n = n / Jac;
    
    // Check if initial normals are available and ensure consistency
    if (cap_initial_normals_.ncols() > 0 && cap_initial_normals_.nrows() == nsd) {
        // Check bounds for element index
        if (e < 0 || e >= cap_initial_normals_.ncols()) {
            throw std::runtime_error("[ZeroDBoundaryCondition::compute_cap_jacobian_and_normal] Element index e=" + 
                                    std::to_string(e) + " is out of bounds for cap_initial_normals_ (ncols=" + 
                                    std::to_string(cap_initial_normals_.ncols()) + ").");
        }
        
        // Get initial normal for this element
        Vector<double> n0(nsd);
        for (int i = 0; i < nsd; i++) {
            n0(i) = cap_initial_normals_(i, e);
        }
        
        // Normalize initial normal (in case it's not normalized)
        double n0_norm = sqrt(utils::norm(n0));
        if (!utils::is_zero(n0_norm)) {
            n0 = n0 / n0_norm;
            
            // Check if calculated normal and initial normal point in opposite directions
            // Dot product < 0 means they point in opposite directions
            double dot_product = 0.0;
            for (int i = 0; i < nsd; i++) {
                dot_product += n(i) * n0(i);
            }
            
            // If dot product is negative, flip the calculated normal
            if (dot_product < 0.0) {
                n = -n;
            }
        }
    }
    
    return std::make_pair(Jac, n);
}

double ZeroDBoundaryCondition::integrate_scalar_at_gauss_point(ComMod& com_mod, const Array<double>& s, int l,
                                                                int e, int g, const std::unordered_map<int, int>& gnNo_to_tnNo)
{
    // Safety checks
    if (cap_face_ == nullptr) {
        throw std::runtime_error("[ZeroDBoundaryCondition::integrate_scalar_at_gauss_point] cap_face_ is null.");
    }
    if (cap_face_->N.nrows() == 0 || cap_face_->N.ncols() == 0) {
        throw std::runtime_error("[ZeroDBoundaryCondition::integrate_scalar_at_gauss_point] cap_face_->N is not allocated.");
    }
    if (g < 0 || g >= cap_face_->N.ncols()) {
        throw std::runtime_error("[ZeroDBoundaryCondition::integrate_scalar_at_gauss_point] Gauss point index g=" + 
                                std::to_string(g) + " is out of bounds (N.ncols()=" + std::to_string(cap_face_->N.ncols()) + ").");
    }
    if (cap_face_->IEN.nrows() == 0 || cap_face_->IEN.ncols() == 0) {
        throw std::runtime_error("[ZeroDBoundaryCondition::integrate_scalar_at_gauss_point] cap_face_->IEN is not allocated.");
    }
    if (e < 0 || e >= cap_face_->IEN.ncols()) {
        throw std::runtime_error("[ZeroDBoundaryCondition::integrate_scalar_at_gauss_point] Element index e=" + 
                                std::to_string(e) + " is out of bounds (IEN.ncols()=" + std::to_string(cap_face_->IEN.ncols()) + ").");
    }
    
    double sHat = 0.0;
    
    for (int a = 0; a < cap_face_->eNoN; a++) {
        // Bounds check for IEN access
        if (a < 0 || a >= cap_face_->IEN.nrows()) {
            throw std::runtime_error("[ZeroDBoundaryCondition::integrate_scalar_at_gauss_point] Node index a=" + 
                                    std::to_string(a) + " is out of bounds (IEN.nrows()=" + std::to_string(cap_face_->IEN.nrows()) + ").");
        }
        // Bounds check for N access
        if (a < 0 || a >= cap_face_->N.nrows()) {
            throw std::runtime_error("[ZeroDBoundaryCondition::integrate_scalar_at_gauss_point] Node index a=" + 
                                    std::to_string(a) + " is out of bounds (N.nrows()=" + std::to_string(cap_face_->N.nrows()) + ").");
        }
        
        // IEN contains gnNo index, map it to tnNo index
        int gnNo_idx = cap_face_->IEN(a, e);
        auto it = gnNo_to_tnNo.find(gnNo_idx);
        if (it == gnNo_to_tnNo.end()) {
            throw std::runtime_error("[ZeroDBoundaryCondition::integrate_scalar_at_gauss_point] IEN entry (element " + 
                                    std::to_string(e) + ", node " + std::to_string(a) + 
                                    ") contains invalid gnNo index " + std::to_string(gnNo_idx) + 
                                    " not found in com_mod.ltg mapping.");
        }
        int Ac = it->second;  // tnNo index
        
        // Bounds checking
        if (Ac < 0 || Ac >= com_mod.tnNo) {
            std::string msg = "[ZeroDBoundaryCondition::integrate_scalar_at_gauss_point] Invalid node index Ac=" + 
                             std::to_string(Ac);
            throw std::runtime_error(msg);
        }
        // Check array bounds for s
        if (l >= s.nrows() || Ac >= s.ncols()) {
            std::string msg = std::string("[ZeroDBoundaryCondition::integrate_scalar_at_gauss_point] Array s bounds exceeded: ") +
                             std::string("s.nrows()=") + std::to_string(s.nrows()) + std::string(", s.ncols()=") + std::to_string(s.ncols()) +
                             std::string(", l=") + std::to_string(l) + std::string(", Ac=") + std::to_string(Ac);
            throw std::runtime_error(msg);
        }
        sHat += cap_face_->N(a, g) * s(l, Ac);
    }
    
    return sHat;
}

double ZeroDBoundaryCondition::integrate_vector_at_gauss_point(ComMod& com_mod, const Array<double>& s, int l, int nsd,
                                                                int e, int g, const Vector<double>& n,
                                                                const std::unordered_map<int, int>& gnNo_to_tnNo)
{
    // Safety checks
    if (cap_face_ == nullptr) {
        throw std::runtime_error("[ZeroDBoundaryCondition::integrate_vector_at_gauss_point] cap_face_ is null.");
    }
    if (cap_face_->N.nrows() == 0 || cap_face_->N.ncols() == 0) {
        throw std::runtime_error("[ZeroDBoundaryCondition::integrate_vector_at_gauss_point] cap_face_->N is not allocated.");
    }
    if (g < 0 || g >= cap_face_->N.ncols()) {
        throw std::runtime_error("[ZeroDBoundaryCondition::integrate_vector_at_gauss_point] Gauss point index g=" + 
                                std::to_string(g) + " is out of bounds (N.ncols()=" + std::to_string(cap_face_->N.ncols()) + ").");
    }
    if (cap_face_->IEN.nrows() == 0 || cap_face_->IEN.ncols() == 0) {
        throw std::runtime_error("[ZeroDBoundaryCondition::integrate_vector_at_gauss_point] cap_face_->IEN is not allocated.");
    }
    if (e < 0 || e >= cap_face_->IEN.ncols()) {
        throw std::runtime_error("[ZeroDBoundaryCondition::integrate_vector_at_gauss_point] Element index e=" + 
                                std::to_string(e) + " is out of bounds (IEN.ncols()=" + std::to_string(cap_face_->IEN.ncols()) + ").");
    }
    if (n.size() != nsd) {
        throw std::runtime_error("[ZeroDBoundaryCondition::integrate_vector_at_gauss_point] Normal vector size mismatch: " +
                                std::to_string(n.size()) + " != " + std::to_string(nsd));
    }
    
    double sHat = 0.0;
    
    for (int a = 0; a < cap_face_->eNoN; a++) {
        // Bounds check for IEN access
        if (a < 0 || a >= cap_face_->IEN.nrows()) {
            throw std::runtime_error("[ZeroDBoundaryCondition::integrate_vector_at_gauss_point] Node index a=" + 
                                    std::to_string(a) + " is out of bounds (IEN.nrows()=" + std::to_string(cap_face_->IEN.nrows()) + ").");
        }
        // Bounds check for N access
        if (a < 0 || a >= cap_face_->N.nrows()) {
            throw std::runtime_error("[ZeroDBoundaryCondition::integrate_vector_at_gauss_point] Node index a=" + 
                                    std::to_string(a) + " is out of bounds (N.nrows()=" + std::to_string(cap_face_->N.nrows()) + ").");
        }
        
        // IEN contains gnNo index, map it to tnNo index
        int gnNo_idx = cap_face_->IEN(a, e);
        auto it = gnNo_to_tnNo.find(gnNo_idx);
        if (it == gnNo_to_tnNo.end()) {
            throw std::runtime_error("[ZeroDBoundaryCondition::integrate_vector_at_gauss_point] IEN entry (element " + 
                                    std::to_string(e) + ", node " + std::to_string(a) + 
                                    ") contains invalid gnNo index " + std::to_string(gnNo_idx) + 
                                    " not found in com_mod.ltg mapping.");
        }
        int Ac = it->second;  // tnNo index
        
        // Bounds checking
        if (Ac < 0 || Ac >= com_mod.tnNo) {
            std::string msg = "[ZeroDBoundaryCondition::integrate_vector_at_gauss_point] Invalid node index Ac=" + 
                             std::to_string(Ac);
            throw std::runtime_error(msg);
        }
        // Check array bounds for s
        if (l + nsd - 1 >= s.nrows() || Ac >= s.ncols()) {
            std::string msg = std::string("[ZeroDBoundaryCondition::integrate_vector_at_gauss_point] Array s bounds exceeded: ") +
                             std::string("s.nrows()=") + std::to_string(s.nrows()) + std::string(", s.ncols()=") + std::to_string(s.ncols()) +
                             std::string(", l=") + std::to_string(l) + std::string(", nsd=") + std::to_string(nsd) + std::string(", Ac=") + std::to_string(Ac);
            throw std::runtime_error(msg);
        }
        for (int i = 0; i < nsd; i++) {
            sHat += cap_face_->N(a, g) * s(l + i, Ac) * n(i);
        }
    }
    
    return sHat;
}

double ZeroDBoundaryCondition::integrate_over_cap(ComMod& com_mod, const CmMod& cm_mod, const Array<double>& s, 
                                                   int l, std::optional<int> u, consts::MechanicalConfigurationType cfg)
{
    using namespace consts;
    
    if (!has_cap()) {
        return 0.0;
    }
    
    // Get a local pointer to cap_face_ to avoid it being reset between check and use
    // This prevents race conditions if cap_face_ is modified concurrently
    faceType* cap_face = cap_face_.get();
    if (cap_face == nullptr) {
        throw std::runtime_error("[ZeroDBoundaryCondition::integrate_over_cap] Cap face is null.");
    }
    
    // Safety check: ensure cap_face_ is properly initialized
    if (cap_face->nEl == 0 || cap_face->nNo == 0) {
        throw std::runtime_error("[ZeroDBoundaryCondition::integrate_over_cap] Cap face is not properly initialized.");
    }
    
    // Ensure integration requirements are initialized
    if (cap_face->nG == 0 || cap_face->w.size() == 0 || cap_face->N.nrows() == 0 || cap_gnNo_to_tnNo_.empty()) {
        throw std::runtime_error("[ZeroDBoundaryCondition::integrate_over_cap] Cap integration not initialized. Call initialize_cap_integration() first.");
    }
    if (cap_face->Nx.nslices() == 0 || cap_face->Nx.nslices() != cap_face->nG) {
        throw std::runtime_error(std::string("[ZeroDBoundaryCondition::integrate_over_cap] cap_face_->Nx is not properly allocated: ") +
                                "nG=" + std::to_string(cap_face->nG) + ", Nx.nslices()=" + std::to_string(cap_face->Nx.nslices()) + 
                                ". Call initialize_cap_integration() first.");
    }
    
    int nsd = com_mod.nsd;
    int insd = nsd - 1;
    
    
    // Determine if scalar or vector integration
    int u_val = u.has_value() ? u.value() : l;
    bool is_scalar = (u_val == l);
    
    // Custom integration for cap (since cap doesn't share elements with mesh)
    // We compute the normal directly from cap geometry without using gnnb
    double result = 0.0;
    
    // Loop over cap elements
    for (int e = 0; e < cap_face->nEl; e++) {
        // Update cap element position based on configuration
        Array<double> xl = update_cap_element_position(com_mod, e, cap_gnNo_to_tnNo_, cfg);
        
        // Loop over Gauss points
        for (int g = 0; g < cap_face->nG; g++) {
            // Bounds check for w access
            if (g < 0 || g >= cap_face->w.size()) {
                throw std::runtime_error("[ZeroDBoundaryCondition::integrate_over_cap] Gauss point index g=" + 
                                        std::to_string(g) + " is out of bounds (w.size()=" + std::to_string(cap_face->w.size()) + ").");
            }
            
            // Compute Jacobian and normal vector
            auto [Jac, n] = compute_cap_jacobian_and_normal(xl, e, g, nsd, insd);
            
            // Compute integrand at this Gauss point
            double sHat = 0.0;
            if (is_scalar) {
                // Scalar integration: ∫ s dA
                sHat = integrate_scalar_at_gauss_point(com_mod, s, l, e, g, cap_gnNo_to_tnNo_);
            } else {
                // Vector integration: ∫ (s dot n) dA
                sHat = integrate_vector_at_gauss_point(com_mod, s, l, nsd, e, g, n, cap_gnNo_to_tnNo_);
            }
            
            // Add contribution to integral
            result += cap_face->w(g) * Jac * sHat;
        }
    }
    
    // Reduce across processors if needed
    auto& cm = com_mod.cm;
    if (!cm.seq()) {
        result = cm.reduce(cm_mod, result);
    }
    
    return result;
}

// =========================================================================
// Cap area calculation
// =========================================================================

double ZeroDBoundaryCondition::compute_cap_area(ComMod& com_mod, const CmMod& cm_mod)
{
    using namespace consts;
    
    if (!has_cap()) {
        return 0.0;
    }
    
    // Get a local pointer to cap_face_ to avoid it being reset between check and use
    faceType* cap_face = cap_face_.get();
    if (cap_face == nullptr) {
        throw std::runtime_error("[ZeroDBoundaryCondition::compute_cap_area] Cap face is null.");
    }
    
    // Safety check: ensure cap_face_ is properly initialized
    if (cap_face->nEl == 0 || cap_face->nNo == 0) {
        throw std::runtime_error("[ZeroDBoundaryCondition::compute_cap_area] Cap face is not properly initialized.");
    }
    
    // Ensure integration requirements are initialized
    if (cap_face->nG == 0 || cap_face->w.size() == 0 || cap_face->N.nrows() == 0 || cap_gnNo_to_tnNo_.empty()) {
        throw std::runtime_error("[ZeroDBoundaryCondition::compute_cap_area] Cap integration not initialized. Call initialize_cap_integration() first.");
    }
    
    int nsd = com_mod.nsd;
    int insd = nsd - 1;
    
    // Compute area by integrating 1 over the cap surface (using reference configuration)
    double area = 0.0;
    
    // Loop over cap elements
    for (int e = 0; e < cap_face->nEl; e++) {
        // Update cap element position (using reference configuration for area calculation)
        Array<double> xl = update_cap_element_position(com_mod, e, cap_gnNo_to_tnNo_, MechanicalConfigurationType::reference);
        
        // Loop over Gauss points
        for (int g = 0; g < cap_face->nG; g++) {
            // Compute Jacobian and normal vector
            auto [Jac, n] = compute_cap_jacobian_and_normal(xl, e, g, nsd, insd);
            
            // Integrate 1 over the surface: ∫ 1 dA = Σ w_g * Jac_g
            area += cap_face->w(g) * Jac;
        }
    }
    
    // Reduce across processors if needed
    auto& cm = com_mod.cm;
    if (!cm.seq()) {
        area = cm.reduce(cm_mod, area);
    }
    
    
    return area;
}

// =========================================================================
// Cap valM computation (precomputed normal integrals)
// =========================================================================

void ZeroDBoundaryCondition::compute_cap_valM(ComMod& com_mod, const CmMod& cm_mod, consts::MechanicalConfigurationType cfg)
{
    using namespace consts;
    
    if (!has_cap()) {
        cap_valM_.resize(0, 0);
        return;
    }
    
    // Get a local pointer to cap_face_ to avoid it being reset between check and use
    faceType* cap_face = cap_face_.get();
    if (cap_face == nullptr) {
        throw std::runtime_error("[ZeroDBoundaryCondition::compute_cap_valM] Cap face is null.");
    }
    
    // Safety check: ensure cap_face_ is properly initialized
    if (cap_face->nEl == 0 || cap_face->nNo == 0) {
        throw std::runtime_error("[ZeroDBoundaryCondition::compute_cap_valM] Cap face is not properly initialized.");
    }
    
    // Ensure integration requirements are initialized
    if (cap_face->nG == 0 || cap_face->w.size() == 0 || cap_face->N.nrows() == 0 || cap_gnNo_to_tnNo_.empty()) {
        throw std::runtime_error("[ZeroDBoundaryCondition::compute_cap_valM] Cap integration not initialized. Call initialize_cap_integration() first.");
    }
    if (cap_face->Nx.nslices() == 0 || cap_face->Nx.nslices() != cap_face->nG) {
        throw std::runtime_error(std::string("[ZeroDBoundaryCondition::compute_cap_valM] cap_face_->Nx is not properly allocated: ") +
                                "nG=" + std::to_string(cap_face->nG) + ", Nx.nslices()=" + std::to_string(cap_face->Nx.nslices()) + 
                                ". Call initialize_cap_integration() first.");
    }
    
    int nsd = com_mod.nsd;
    int cap_nNo = cap_face->nNo;
    
    // Initialize cap_valM_ to zero: nsd x cap_nNo (stored with cap face-local indices, like face.valM)
    cap_valM_.resize(nsd, cap_nNo);
    cap_valM_ = 0.0;
    
    // Create mapping from gnNo to cap face-local index for efficient lookup
    std::unordered_map<int, int> gnNo_to_cap_local;
    for (int a = 0; a < cap_nNo; a++) {
        int gnNo = cap_face->gN(a);
        gnNo_to_cap_local[gnNo] = a;
    }
    
    // Loop over cap elements
    for (int e = 0; e < cap_face->nEl; e++) {
        // Update cap element position based on configuration
        Array<double> xl = update_cap_element_position(com_mod, e, cap_gnNo_to_tnNo_, cfg);
        
        // Loop over Gauss points
        for (int g = 0; g < cap_face->nG; g++) {
            // Compute Jacobian and normal vector
            auto [Jac, n] = compute_cap_jacobian_and_normal(xl, e, g, nsd, nsd - 1);
            
            // Accumulate ∫ N_A n_i dΓ at each node
            // Similar to fsi_ls_upd, but accumulate by cap face-local node index
            for (int a = 0; a < cap_face->eNoN; a++) {
                // Get the global node index from IEN
                int gnNo_idx = cap_face->IEN(a, e);
                
                // Map gnNo to cap face-local index
                auto it = gnNo_to_cap_local.find(gnNo_idx);
                if (it == gnNo_to_cap_local.end()) {
                    throw std::runtime_error("[ZeroDBoundaryCondition::compute_cap_valM] IEN entry (element " + 
                                            std::to_string(e) + ", node " + std::to_string(a) + 
                                            ") contains invalid gnNo index " + std::to_string(gnNo_idx) + 
                                            " not found in cap face nodes.");
                }
                int cap_a = it->second;  // cap face-local index
                
                // Bounds checking
                if (cap_a < 0 || cap_a >= cap_nNo) {
                    throw std::runtime_error("[ZeroDBoundaryCondition::compute_cap_valM] Invalid cap face-local index cap_a=" + 
                                            std::to_string(cap_a) + " (cap_nNo=" + std::to_string(cap_nNo) + ")");
                }
                
                // Accumulate: cap_valM_(i,cap_a) += N(a,g) * w(g) * Jac * n(i)
                // Note: n is normalized, so we multiply by Jac to get the surface element
                // Stored with cap face-local indices, matching GenBC implementation
                for (int i = 0; i < nsd; i++) {
                    cap_valM_(i, cap_a) += cap_face->N(a, g) * cap_face->w(g) * Jac * n(i);
                }
            }
        }
    }
    
    // Reduce across processors if needed
    auto& cm = com_mod.cm;
    if (!cm.seq()) {
        for (int i = 0; i < nsd; i++) {
            Vector<double> row = cap_valM_.row(i);
            row = cm.reduce(cm_mod, row);
            cap_valM_.set_row(i, row);
        }
    }
    
}

// =========================================================================
// Copy cap data to linear solver face structure
// =========================================================================

void ZeroDBoundaryCondition::copy_cap_data_to_linear_solver_face(ComMod& com_mod, const CmMod& cm_mod, 
                                                                   fsi_linear_solver::FSILS_faceType& lhs_face,
                                                                   consts::MechanicalConfigurationType cfg)
{
    if (!has_cap()) {
        // Clear all cap-related fields if no cap is present
        lhs_face.cap_val.resize(0, 0);
        lhs_face.cap_valM.resize(0, 0);
        lhs_face.cap_glob.resize(0);
        return;
    }
    
    // Compute cap_valM if not already computed
    compute_cap_valM(com_mod, cm_mod, cfg);
    
    const faceType* cap_face = cap_face_.get();
    if (cap_face == nullptr || cap_face->nNo == 0) {
        lhs_face.cap_val.resize(0, 0);
        lhs_face.cap_valM.resize(0, 0);
        lhs_face.cap_glob.resize(0);
        return;
    }
    
    int nsd = com_mod.nsd;
    int cap_nNo = cap_face->nNo;
    
    // Copy cap_valM_ to cap_val (unpreconditioned values)
    // cap_val will be preconditioned later in precond_diag to produce cap_valM
    lhs_face.cap_val = cap_valM_;
    
    // Initialize cap_valM to zero (will be computed in precond_diag)
    lhs_face.cap_valM.resize(nsd, cap_nNo);
    lhs_face.cap_valM = 0.0;
    
    // Set up cap_glob mapping: cap face-local index -> linear solver index
    // Similar to face.glob, but for cap nodes
    lhs_face.cap_glob.resize(cap_nNo);
    for (int a = 0; a < cap_nNo; a++) {
        int gnNo = cap_face->gN(a);
        int Ac = com_mod.lhs.map(gnNo);
        lhs_face.cap_glob(a) = Ac;
    }
}
