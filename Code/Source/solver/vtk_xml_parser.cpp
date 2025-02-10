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

// The functions defined here replicate the Fortran functions defined in vtkXMLParser.f90.
//
// Volume mesh and face data is read using VTK readers to read VTU and VTP files. Data 
// read in is stored directly into mshType and faceType objects.

#include "vtk_xml_parser.h" 
#include "Array.h" 
#include "Array3.h"

#include <vtkDoubleArray.h>
#include "vtkCellData.h"
#include <vtkGenericCell.h>
#include <vtkGeometryFilter.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <cctype>

namespace vtk_xml_parser {

// VTK file extensions.
const std::string VtkFileExtentions::VTK_VTU_EXTENSION = "vtu";
const std::string VtkFileExtentions::VTK_VTP_EXTENSION = "vtp";

/// Map used to convert VTK cell types to number of nodes per element.
std::map<unsigned char,int> vtk_cell_to_elem {
  {VTK_HEXAHEDRON, 8},
  {VTK_LINE, 2},
  {VTK_QUAD, 4},
  {VTK_TETRA, 4},
  {VTK_TRIANGLE, 3},
  {VTK_WEDGE, 6},
  {VTK_QUADRATIC_TRIANGLE, 6},
  {VTK_BIQUADRATIC_TRIANGLE, 7},
  {VTK_QUADRATIC_QUAD, 8},
  {VTK_BIQUADRATIC_QUAD, 9},
  {VTK_QUADRATIC_TETRA, 10},
  {VTK_QUADRATIC_HEXAHEDRON, 20},
  {VTK_TRIQUADRATIC_HEXAHEDRON, 27}
};

/// Map used to convert VTK cell types to consts::ElementType. 
std::map<unsigned char,consts::ElementType> vtk_cell_to_elem_type {
  {VTK_HEXAHEDRON, consts::ElementType::HEX8}, 
  {VTK_LINE, consts::ElementType::LIN1},
  {VTK_QUAD, consts::ElementType::QUD4},
  {VTK_TETRA, consts::ElementType::TET4},
  {VTK_TRIANGLE, consts::ElementType::TRI3},
  {VTK_WEDGE, consts::ElementType::WDG}, 
  {VTK_QUADRATIC_TRIANGLE, consts::ElementType::TRI6},
  {VTK_QUADRATIC_QUAD, consts::ElementType::QUD8},
  {VTK_BIQUADRATIC_QUAD, consts::ElementType::QUD9}, 
  {VTK_QUADRATIC_TETRA, consts::ElementType::TET10},
  {VTK_QUADRATIC_HEXAHEDRON, consts::ElementType::HEX20},
  {VTK_TRIQUADRATIC_HEXAHEDRON, consts::ElementType::HEX27}
};

std::map<unsigned char, std::vector<std::vector<int>>> vtk_cell_ordering {
  {VTK_HEXAHEDRON, {{0,3,2,1},
                    {4,5,6,7},
                    {0,1,5,4},
                    {1,2,6,5},
                    {2,3,7,6},
                    {3,0,4,7}}},
  {VTK_LINE, {{0},
              {1}}},
  {VTK_QUAD, {{0,1},
              {1,2},
              {2,3},
              {3,0}}},
  {VTK_TETRA, {{0,1,2},
               {0,1,3},
               {1,2,3},
               {2,0,3}}},
  {VTK_TRIANGLE, {{0,1},
                  {1,2},
                  {2,0}}},
  {VTK_WEDGE, {{0,1,2},
               {3,4,5},
               {0,1,4,3},
               {1,2,5,4},
               {2,0,3,5}}},
  {VTK_QUADRATIC_TRIANGLE, {{0,3,1},
                            {1,4,2},
                            {2,5,0}}},
  {VTK_BIQUADRATIC_TRIANGLE, {{0,3,1},
                              {1,4,2},
                              {2,5,0}}},
  {VTK_QUADRATIC_QUAD, {{0,4,1},
                        {1,5,2},
                        {2,6,3},
                        {3,7,0}}},
  {VTK_BIQUADRATIC_QUAD, {{0,4,1},
                          {1,5,2},
                          {2,6,3},
                          {3,7,0}}},
  {VTK_QUADRATIC_TETRA, {{0,1,2,4,5,6},
                       {0,3,1,7,8,4},
                       {1,3,2,8,9,5},
                       {2,3,0,9,7,6}}},
  {VTK_QUADRATIC_HEXAHEDRON, {{0,11,3,10,2,9,1,8},
                              {4,12,5,13,6,14,7,15},
                              {0,8,1,17,5,12,4,16},
                              {1,9,2,18,6,13,5,17},
                              {2,10,3,19,7,14,6,18},
                              {3,11,0,16,4,15,7,19}}},
  {VTK_TRIQUADRATIC_HEXAHEDRON, {{0,11,3,10,2,9,1,8,24},
                                {4,12,5,13,6,14,7,15,25},
                                {0,8,1,17,5,12,4,16,22},
                                {1,9,2,18,6,13,5,17,21},
                                {2,10,3,19,7,14,6,18,23},
                                {3,11,0,16,4,15,7,19,20}}}
};
/// Names of data arrays store in VTK mesh files.
const std::string NODE_IDS_NAME("GlobalNodeID");
const std::string ELEMENT_IDS_NAME("GlobalElementID");

/////////////////////////////////////////////////////////////////
//             I n t e r n a l  U t i l i t i e s              //
/////////////////////////////////////////////////////////////////

/// @brief Get the mesh nodes per element and ordering. 
///
void get_mesh_ordering(const int num_elems, vtkSmartPointer<vtkUnsignedCharArray> cell_types, int& np_elem, 
    std::vector<std::vector<int>>& ordering)
{
  int num_hex = 0;
  int num_line = 0;
  int num_quad = 0;
  int num_tet = 0;
  int num_tri = 0;
  int num_unknown = 0;
  int num_wedge = 0;
  int num_quadratic_tri = 0;
  int num_biquadratic_tri = 0;
  int num_quadratic_quad = 0;
  int num_biquadratic_quad = 0;
  int num_quadratic_tetra = 0;
  int num_quadratic_hexahedron = 0;
  int num_triquadratic_hexahedron = 0;

  for (int i = 0; i < num_elems; i++) {
    switch (cell_types->GetValue(i)) {
      case VTK_HEXAHEDRON:
        num_hex += 1;
      break;

      case VTK_LINE:
        num_line += 1;
      break;

      case VTK_QUAD:
        num_quad += 1;
      break;

      case VTK_TETRA:
        num_tet += 1;
      break;

      case VTK_TRIANGLE:
        num_tri += 1;
      break;

      case VTK_WEDGE:
        num_wedge += 1;
      break;

      case VTK_QUADRATIC_TRIANGLE:
        num_quadratic_tri += 1;
      break;

      case VTK_BIQUADRATIC_TRIANGLE:
        num_biquadratic_tri += 1;
      break;

      case VTK_QUADRATIC_QUAD:
        num_quadratic_quad += 1;
      break;

      case VTK_BIQUADRATIC_QUAD:
        num_biquadratic_quad += 1;
      break;

      case VTK_QUADRATIC_TETRA:
        num_quadratic_tetra += 1;
      break;

      case VTK_QUADRATIC_HEXAHEDRON:
        num_quadratic_hexahedron += 1;
      break;

      case VTK_TRIQUADRATIC_HEXAHEDRON:
        num_triquadratic_hexahedron += 1;
      break;

      default:
        num_unknown += 1;
      break;
    }
  }

  #ifdef debug_get_mesh_ordering
  std::cout << "[store_element_conn] num_line: " << num_line <<  std::endl;
  std::cout << "[store_element_conn] num_quad: " << num_quad <<  std::endl;
  std::cout << "[store_element_conn] num_hex: " << num_hex <<  std::endl;
  std::cout << "[store_element_conn] num_tet: " << num_tet <<  std::endl;
  std::cout << "[store_element_conn] num_tri: " << num_tri <<  std::endl;
  std::cout << "[store_element_conn] num_unknown: " << num_unknown <<  std::endl;
  #endif

  np_elem = 0;
  ordering.clear();

  if (num_line != 0) {
    np_elem = vtk_cell_to_elem[VTK_LINE];
    ordering = vtk_cell_ordering[VTK_LINE];
  } if (num_hex != 0) {
    np_elem = vtk_cell_to_elem[VTK_HEXAHEDRON];
    ordering = vtk_cell_ordering[VTK_HEXAHEDRON];
  } if (num_quad != 0) {
    np_elem = vtk_cell_to_elem[VTK_QUAD];
    ordering = vtk_cell_ordering[VTK_QUAD];
  } if (num_tet != 0) {
    np_elem = vtk_cell_to_elem[VTK_TETRA];
    ordering = vtk_cell_ordering[VTK_TETRA];
  } if (num_tri != 0) {
    np_elem = vtk_cell_to_elem[VTK_TRIANGLE];
    ordering = vtk_cell_ordering[VTK_TRIANGLE];
  } if (num_wedge != 0) {
    np_elem = vtk_cell_to_elem[VTK_WEDGE];
    ordering = vtk_cell_ordering[VTK_WEDGE];
  }

  // For higher-order elements with mid-side nodes. 
  //
  if (num_quadratic_tri != 0) {
    np_elem = vtk_cell_to_elem[VTK_QUADRATIC_TRIANGLE];
    ordering = vtk_cell_ordering[VTK_QUADRATIC_TRIANGLE];
  } if (num_biquadratic_tri != 0){
    np_elem = vtk_cell_to_elem[VTK_BIQUADRATIC_TRIANGLE];
    ordering = vtk_cell_ordering[VTK_BIQUADRATIC_TRIANGLE];
  } if (num_quadratic_quad != 0) {
    np_elem = vtk_cell_to_elem[VTK_QUADRATIC_QUAD];
    ordering = vtk_cell_ordering[VTK_QUADRATIC_QUAD];
  } if (num_biquadratic_quad != 0) {
    np_elem = vtk_cell_to_elem[VTK_BIQUADRATIC_QUAD];
    ordering = vtk_cell_ordering[VTK_BIQUADRATIC_QUAD];
  } if (num_quadratic_tetra != 0) {
    np_elem = vtk_cell_to_elem[VTK_QUADRATIC_TETRA];
    ordering = vtk_cell_ordering[VTK_QUADRATIC_TETRA];
  } if (num_quadratic_hexahedron != 0) {
    np_elem = vtk_cell_to_elem[VTK_QUADRATIC_HEXAHEDRON];
    ordering = vtk_cell_ordering[VTK_QUADRATIC_HEXAHEDRON];
  } if (num_triquadratic_hexahedron != 0) {
    np_elem = vtk_cell_to_elem[VTK_TRIQUADRATIC_HEXAHEDRON];
    ordering = vtk_cell_ordering[VTK_TRIQUADRATIC_HEXAHEDRON];
  }
}

/// @brief Store element connectivity into the Face object.
///
/// Face variables set
///   face.nEl - number of elements
///   face.eNoN - number of noders per element
///   face.IEN - element connectivity (num_nodes_per_elem, num_elems)
///
/// Note: The face.IEN array is allocated as in the Fortran code as (face.eNoN, face.nEl).
//
void store_element_conn(vtkSmartPointer<vtkPolyData> vtk_polydata, faceType& face)
{
  #define n_debug_store_element_conn
  #ifdef debug_store_element_conn
  std::cout << "[store_element_conn(polydata)] " << std::endl;
  std::cout << "[store_element_conn(polydata)] ========== store_element_conn(polydata) =========" << std::endl;
  #endif

  vtkObject::GlobalWarningDisplayOff();

  // Get the number of nodes per cell.
  auto cell = vtkGenericCell::New();
  vtk_polydata->GetCell(0, cell);
  int np_elem = cell->GetNumberOfPoints();
  #ifdef debug_store_element_conn
  std::cout << "[store_element_conn(polydata)] np_elem: " << np_elem << std::endl;
  std::cout << "[store_element_conn(polydata)] cell type: " << vtk_polydata->GetCellType(0) << std::endl;
  #endif

  // Allocate connectivity array.
  auto num_elems = vtk_polydata->GetNumberOfCells();
  face.nEl = num_elems;
  face.eNoN = np_elem; 
  face.IEN = Array<int>(np_elem, num_elems);

  for (int i = 0; i < num_elems; i++) {
    vtk_polydata->GetCell(i, cell);
    auto num_cell_pts = cell->GetNumberOfPoints();
    for (int j = 0; j < num_cell_pts; j++) {
      auto id = cell->PointIds->GetId(j);
      face.IEN(j,i) = id;
    }
  }
}

void store_element_conn(vtkSmartPointer<vtkUnstructuredGrid> vtk_ugrid, faceType& face)
{
  // Get the number of nodes per cell.
  auto cell = vtkGenericCell::New();
  vtk_ugrid->GetCell(0, cell);
  int np_elem = cell->GetNumberOfPoints();

  auto num_elems = vtk_ugrid->GetNumberOfCells();
  face.nEl = num_elems;
  face.eNoN = np_elem; 
  face.IEN = Array<int>(np_elem, num_elems);

  for (int i = 0; i < num_elems; i++) {
    vtk_ugrid->GetCell(i, cell);
    auto num_cell_pts = cell->GetNumberOfPoints();
    for (int j = 0; j < num_cell_pts; j++) {
      auto id = cell->PointIds->GetId(j);
      face.IEN(j,i) = id;
    }
  }
}

void store_element_conn(vtkSmartPointer<vtkPolyData> vtk_polydata, mshType& mesh)
{
  #ifdef debug_store_element_conn
  std::cout << "[store_element_conn(polydata,mesh)] " << std::endl;
  std::cout << "[store_element_conn(polydata,mesh)] ========== store_element_conn(polydata,mesh) =========" << std::endl;
  #endif

  // Get the number of nodes per cell.
  auto cell = vtkGenericCell::New();
  vtk_polydata->GetCell(0, cell);
  int np_elem = cell->GetNumberOfPoints();

  // Allocate connectivity array.
  auto num_elems = vtk_polydata->GetNumberOfCells();
  mesh.gnEl = num_elems;
  mesh.eNoN = np_elem;
  mesh.gIEN = Array<int>(np_elem, num_elems);
  #ifdef debug_store_element_conn
  std::cout << "[store_element_conn(polydata,mesh)] num_elems: " << num_elems << std::endl;
  std::cout << "[store_element_conn(polydata,mesh)] np_elem: " << np_elem << std::endl;
  #endif

  for (int i = 0; i < num_elems; i++) {
    //int elem_id = elem_ids->GetValue(i);
    vtk_polydata->GetCell(i, cell);
    auto num_cell_pts = cell->GetNumberOfPoints();
    for (int j = 0; j < num_cell_pts; j++) {
      auto id = cell->PointIds->GetId(j);
      mesh.gIEN(j,i) = id;
    }
  }
}

/// @brief Store VTK vtkUnstructuredGrid cell connectivity into a mesh element connectivity. 
///
/// Mesh variables set
///   mesh.gnEl - number of elements
///   mesh.eNoN - number of noders per element
///   mesh.gIEN - element connectivity (num_nodes_per_elem, num_elems)
///
/// Note: The mesh.gIEN array is allocated as in the Fortran code as (np_elem, num_elems).
//
void store_element_conn(vtkSmartPointer<vtkUnstructuredGrid> vtk_ugrid, mshType& mesh)
{
  #ifdef debug_store_element_conn
  std::cout << "[store_element_conn(ugrid)] " << std::endl;
  std::cout << "[store_element_conn(ugrid)] ========== store_element_conn(ugrid) =========" << std::endl;
  #endif
  vtkSmartPointer<vtkUnsignedCharArray> cell_types = vtk_ugrid->GetCellTypesArray();
  auto num_elems = vtk_ugrid->GetNumberOfCells();
  #ifdef debug_store_element_conn
  std::cout << "[store_element_conn(ugrid)] num_elems: " << num_elems << std::endl;
  #endif

  // Get the mesh nodes per element and ordering. 
  //
  int np_elem = 0; 
  std::vector<std::vector<int>> ordering;
  get_mesh_ordering(num_elems, cell_types, np_elem, ordering);

  // For generic higher-order elements.
  //
  if (np_elem == 0) { 
    auto cell = vtkGenericCell::New();
    vtk_ugrid->GetCell(0, cell);
    np_elem = cell->GetNumberOfPoints();
  }

  mesh.gnEl = num_elems;
  mesh.eNoN = np_elem; 
  mesh.gIEN = Array<int>(np_elem, num_elems);
  mesh.ordering = ordering;

  #ifdef debug_store_element_conn
  std::cout << "[store_element_conn] np_elem: " << np_elem <<  std::endl;
  std::cout << "[store_element_conn] Number of elements: " << num_elems <<  std::endl;
  std::cout << "[store_element_conn] Number of nodes per element: " << np_elem <<  std::endl;
  #endif

  auto cell = vtkGenericCell::New();
  for (int i = 0; i < num_elems; i++) {
    vtk_ugrid->GetCell(i, cell);
    auto dim = cell->GetCellDimension();
    auto num_cell_pts = cell->GetNumberOfPoints();
    if (num_cell_pts != np_elem) { 
      throw std::runtime_error("[store_element_conn] Error in VTK mesh data for mesh '" + mesh.name + "'.");
    }
    for (int j = 0; j < num_cell_pts; j++) {
      auto id = cell->PointIds->GetId(j);
      //mesh.gIEN(i,j) = id;
      mesh.gIEN(j,i) = id;
    }
  }
}

/// @brief Store element IDs into a mshType object. 
///
/// \todo [NOTE] Are element IDs used?
//
void store_element_ids(vtkSmartPointer<vtkUnstructuredGrid> vtk_ugrid, mshType& mesh)
{
  auto elem_ids = vtkIntArray::SafeDownCast(vtk_ugrid->GetCellData()->GetArray(ELEMENT_IDS_NAME.c_str()));
  if (elem_ids == nullptr) { 
    throw std::runtime_error("No '" + ELEMENT_IDS_NAME + "' data of type Int32 found in VTK mesh.");
  }
  int num_elem_ids = elem_ids->GetNumberOfTuples();
  for (int i = 0; i < num_elem_ids; i++) { 
  }
}

/// @brief Store element IDs into a faceType object.
///
/// Face data set
///   face.gE
//
void store_element_ids(vtkSmartPointer<vtkPolyData> vtk_polydata, faceType& face)
{
  auto elem_ids = vtkIntArray::SafeDownCast(vtk_polydata->GetCellData()->GetArray(ELEMENT_IDS_NAME.c_str()));
  if (elem_ids == nullptr) {
    throw std::runtime_error("No '" + ELEMENT_IDS_NAME + "' data of type Int32 found in VTK mesh.");
    return;
  }
  #ifdef debug_store_element_ids
  std::cout << "[store_element_ids] Allocate face.gE ... " <<  std::endl;
  #endif
  int num_elem_ids = elem_ids->GetNumberOfTuples();
  face.gE = Vector<int>(num_elem_ids);
  // [NOTE] It is not clear how these IDs are used but if they
  // index into arrays or vectors then they need to be offset by -1.
  for (int i = 0; i < num_elem_ids; i++) {
    face.gE(i) = elem_ids->GetValue(i) - 1;
  }
}

/// @brief Store element IDs into a faceType object.
///
/// Face data set
///   face.gE
//
void store_element_ids(vtkSmartPointer<vtkUnstructuredGrid> vtk_ugrid, faceType& face)
{
  auto elem_ids = vtkIntArray::SafeDownCast(vtk_ugrid->GetCellData()->GetArray(ELEMENT_IDS_NAME.c_str()));
  if (elem_ids == nullptr) {
    throw std::runtime_error("No '" + ELEMENT_IDS_NAME + "' data of type Int32 found in VTK mesh.");
    return;
  }
  #ifdef debug_store_element_ids
  std::cout << "[store_element_ids] Allocate face.gE ... " <<  std::endl;
  #endif
  int num_elem_ids = elem_ids->GetNumberOfTuples();
  face.gE = Vector<int>(num_elem_ids);
  // [NOTE] It is not clear how these IDs are used but if they
  // index into arrays or vectors then they need to be offset by -1.
  for (int i = 0; i < num_elem_ids; i++) {
    face.gE(i) = elem_ids->GetValue(i) - 1;
  }
}

/// @brief Store nodal coordinates from the VTK points into a face.
///
/// Face data set
///   face.nNo - number of nodes
///   face.x - node coordinates
//
void store_nodal_coords(vtkPoints* points, faceType& face)
{
  vtkIdType num_nodes = points->GetNumberOfPoints();
  face.nNo = num_nodes;
  face.x = Array<double>(3, num_nodes);
  double point[3];
  for (int i = 0; i < num_nodes; i++) {
    points->GetPoint(i, point);
    face.x(0,i) = point[0];
    face.x(1,i) = point[1];
    face.x(2,i) = point[2];
  }
}

/// @brief Store VTK points into a mesh nodal coordinates.
///
/// Mesh variables set
///   mesh.gnNo - number of nodes
///   mesh.x - node coordinates
//
void store_nodal_coords(vtkPoints* points, mshType& mesh)
{
  vtkIdType num_nodes = points->GetNumberOfPoints();
  mesh.gnNo = num_nodes;
  mesh.x = Array<double>(3, num_nodes);

  for (int i = 0; i < num_nodes; i++) {
    double point[3];
    points->GetPoint(i, point);
    mesh.x(0,i) = point[0];
    mesh.x(1,i) = point[1];
    mesh.x(2,i) = point[2];
  }
}

/// @brief Store VTK node IDs data array into a mesh nodal IDs.
///
/// Note: This array appeats to be optional. 
///
/// Mesh variables set
///   mesh.gN - nodal IDs
//
void store_nodal_ids(vtkSmartPointer<vtkUnstructuredGrid> vtk_ugrid, mshType& mesh)
{
  vtkIdType num_nodes = vtk_ugrid->GetNumberOfPoints();
  auto node_ids = vtkIntArray::SafeDownCast(vtk_ugrid->GetPointData()->GetArray(NODE_IDS_NAME.c_str()));
  if (node_ids == nullptr) {
    return;
  }
  mesh.gN = Vector<int>(num_nodes);
  for (int i = 0; i < num_nodes; i++) {
    mesh.gN(i) = node_ids->GetValue(i);
    //std::cout << "[store_nodal_ids] mesh.gN(" << i << "): " << mesh.gN(i) << std::endl;
  }
}

/// @brief Store VTK node IDs data array into a face nodal IDs.
///
void store_nodal_ids(vtkSmartPointer<vtkUnstructuredGrid> vtk_ugrid, faceType& face)
{
  vtkIdType num_nodes = vtk_ugrid->GetNumberOfPoints();
  auto node_ids = vtkIntArray::SafeDownCast(vtk_ugrid->GetPointData()->GetArray(NODE_IDS_NAME.c_str()));
  if (node_ids == nullptr) {
    return;
  }
  face.gN = Vector<int>(num_nodes);
  for (int i = 0; i < num_nodes; i++) {
    face.gN(i) = node_ids->GetValue(i);
    //std::cout << "[store_nodal_ids] mesh.gN(" << i << "): " << mesh.gN(i) << std::endl;
  }
}

/// @brief Store VTK node IDs data array into a face nodal IDs.
///
/// Note: This array appeats to be optional. 
///
/// Face variables set
///   face.gN - nodal IDs
//
void store_nodal_ids(vtkSmartPointer<vtkPolyData> vtk_polydata, faceType& face)
{
  vtkIdType num_nodes = vtk_polydata->GetNumberOfPoints();
  auto node_ids = vtkIntArray::SafeDownCast(vtk_polydata->GetPointData()->GetArray(NODE_IDS_NAME.c_str()));
  if (node_ids == nullptr) {
    return; 
  }
  face.gN = Vector<int>(num_nodes);
  for (int i = 0; i < num_nodes; i++) {
    // [NOTE] It seems that face node IDs are 1 based.
    face.gN(i) = node_ids->GetValue(i) - 1;
    //std::cout << "[store_nodal_ids] face.gN(" << i << "): " << face.gN(i) << std::endl;
  }
}

void store_nodal_ids(vtkSmartPointer<vtkPolyData> vtk_polydata, mshType& mesh)
{
  vtkIdType num_nodes = vtk_polydata->GetNumberOfPoints();
  auto node_ids = vtkIntArray::SafeDownCast(vtk_polydata->GetPointData()->GetArray(NODE_IDS_NAME.c_str()));
  if (node_ids == nullptr) {
    return;
  }
  mesh.gN = Vector<int>(num_nodes);
  for (int i = 0; i < num_nodes; i++) {
    // [NOTE] It seems that face node IDs are 1 based.
    mesh.gN(i) = node_ids->GetValue(i) - 1;
    //std::cout << "[store_nodal_ids] face.gN(" << i << "): " << face.gN(i) << std::endl;
  }
}


/////////////////////////////////////////////////////////////////
//             E x p o s e d    U t i l i t i e s              //
/////////////////////////////////////////////////////////////////

/// @brief Read fiber direction data from a VTK VTU file and copy it into a mesh..
///
/// Data is stored in the mesh for all fiber direction files.
///
/// Mesh variables set
///   mesh.fN - Fiber orientations stored at the element level.  
///             mesh.fN shape is Array<double>(num_fiber_files*nsd, num_elems)
///
/// Arguments:
///   file_name - The name of the VTK VTU file storing fiber direction data
///   data_name - The name of the VTK Cell Data Array storing fiber direction data
///   idx - The id used to identify a fiber file (0, 1, ..., number of fiber files - 1)
///   mesh - The mesh to store the data into.
//
void load_fiber_direction_vtu(const std::string& file_name, const std::string& data_name, const int idx, 
    const int nsd, mshType& mesh)
{
  #ifdef debug_load_fiber_direction_vtu
  std::cout << "[load_fiber_direction_vtu] " << std::endl;
  std::cout << "[load_fiber_direction_vtu] ===== vtk_xml_parser::load_fiber_direction_vtu ===== " << std::endl;
  #endif
  using namespace vtk_xml_parser;

  if (FILE *file = fopen(file_name.c_str(), "r")) {
      fclose(file);
  } else {
    throw std::runtime_error("The fiber direction VTK file '" + file_name + "' can't be read.");
  }

  auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(file_name.c_str());
  reader->Update();
  vtkSmartPointer<vtkUnstructuredGrid> vtk_ugrid = reader->GetOutput();

  vtkIdType num_nodes = vtk_ugrid->GetNumberOfPoints();
  if (num_nodes == 0) {
    throw std::runtime_error("Failed reading the VTK file '" + file_name + "'.");
  }

  vtkIdType num_elems = vtk_ugrid->GetNumberOfCells();
  if (mesh.gnEl != num_elems) {
    throw std::runtime_error("The number of elements (" + std::to_string(num_elems) + 
        ") in the fiber direction VTK file '" + file_name + "' is not equal to the number of elements (" 
        + std::to_string(mesh.gnEl) + ") for the mesh named '" + mesh.name + "'.");
  }

  // Get the 3-component fiber orientation data.
  auto fiber_data = vtkDoubleArray::SafeDownCast(vtk_ugrid->GetCellData()->GetArray(data_name.c_str()));
  if (fiber_data == nullptr) { 
    throw std::runtime_error("No '" + data_name + "' data found in the fiber direction VTK file '" + file_name + "'");
  }

  // Set the fiber orientations data.
  int offset = idx * nsd;
  #ifdef debug_load_fiber_direction_vtu
  std::cout << "[load_fiber_direction_vtu] idx: " << idx << std::endl;
  std::cout << "[load_fiber_direction_vtu] offset: " << offset << std::endl;
  #endif

  for (int e = 0; e < mesh.gnEl; e++) {
    auto fiber_dir = fiber_data->GetTuple(e);
    for (int i = 0; i < nsd; i++) {
      mesh.fN(i+offset, e) = fiber_dir[i];
    }
    //std::cout << "[load_fiber_direction_vtu] e mesh.fN: " << e+1 << " " << mesh.fN.col(e) << std::endl;
  }
}

/// @brief Store a surface mesh read from a VTK .vtp file into a Face object.
//
void load_vtp(const std::string& file_name, faceType& face)
{
  #define n_debug_load_vtp 
  #ifdef debug_load_vtp 
  std::cout << "[load_vtp] " << std::endl;
  std::cout << "[load_vtp] ===== vtk_xml_parser.cpp::load_vtp ===== " << std::endl;
  #endif
  auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(file_name.c_str());
  reader->Update();
  vtkSmartPointer<vtkPolyData> vtk_polydata = reader->GetOutput();

  vtkIdType num_nodes = vtk_polydata->GetNumberOfPoints();
  if (num_nodes == 0) {
    throw std::runtime_error("Failed reading the VTK file '" + file_name + "'.");
  }

  vtkIdType num_elems = vtk_polydata->GetNumberOfCells();
  auto cell = vtkGenericCell::New();
  vtk_polydata->GetCell(0, cell);
  int np_elem = cell->GetNumberOfPoints();
  int elem_type = cell->GetCellType();

  // Set the mesh element type from the vtk cell type.
  try {
    face.eType = vtk_cell_to_elem_type[elem_type];
  } catch (const std::bad_function_call& exception) {
    throw std::runtime_error("[load_vtp] The VTK cell type " + std::to_string(elem_type) + " is not supported.");
  }
  #ifdef debug_load_vtp 
  std::cout << "[load_vtp] Number of nodes: " << num_nodes << std::endl;
  std::cout << "[load_vtp] Number of elements: " << num_elems << std::endl;
  #endif

  // Store nodal coordinates.
  auto points = vtk_polydata->GetPoints();
  store_nodal_coords(points, face);

  // Store nodal IDs.
  store_nodal_ids(vtk_polydata, face);

  // Store element connectivity.
  store_element_conn(vtk_polydata, face);

  // Store element IDs.
  store_element_ids(vtk_polydata, face);

  #ifdef debug_load_vtp 
  std::cout << "[load_vtp] Done. " << std::endl;
  #endif
}

/// @brief Store a surface mesh read from a VTK .vtp file into a Mesh object.
//
void load_vtp(const std::string& file_name, mshType& mesh) 
{
  #ifdef n_debug_load_vtp 
  std::cout << "[load_vtp] " << std::endl;
  std::cout << "[load_vtp] ===== vtk_xml_parser.cpp::load_vtp ===== " << std::endl;
  #endif
  auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(file_name.c_str());
  reader->Update();
  vtkSmartPointer<vtkPolyData> vtk_polydata = reader->GetOutput();

  vtkIdType num_nodes = vtk_polydata->GetNumberOfPoints();
  if (num_nodes == 0) {
    throw std::runtime_error("Failed reading the VTK file '" + file_name + "'.");
  }

  vtkIdType num_elems = vtk_polydata->GetNumberOfCells();
  auto cell = vtkGenericCell::New();
  vtk_polydata->GetCell(0, cell);
  int np_elem = cell->GetNumberOfPoints();
  int elem_type = cell->GetCellType();

  // Set the mesh element type from the vtk cell type.
  try {
    mesh.eType = vtk_cell_to_elem_type[elem_type];
  } catch (const std::bad_function_call& exception) {
    throw std::runtime_error("[load_vtp] The VTK cell type " + std::to_string(elem_type) + " is not supported.");
  }

  #ifdef debug_load_vtp 
  std::cout << "[load_vtp] Number of nodes: " << num_nodes << std::endl;
  std::cout << "[load_vtp] Number of elements: " << num_elems << std::endl;
  #endif

  // Store nodal coordinates.
  auto points = vtk_polydata->GetPoints();
  store_nodal_coords(points, mesh);

  // Store nodal IDs.
  store_nodal_ids(vtk_polydata, mesh);

  // Store element connectivity.
  store_element_conn(vtk_polydata, mesh);

  // Store element IDs.
  //store_element_ids(vtk_polydata, mesh);

  #ifdef debug_load_vtp 
  std::cout << "[load_vtp] Done. " << std::endl;
  #endif
}

/// @brief Read a mesh from a .vtu file.
///
/// Mesh variables set
///   mesh.gnNo - number of nodes
///   mesh.x - node coordinates
///   mesh.gN - node IDs 
///   mesh.gnEl - number of elements
///   mesh.eNoN - number of noders per element
///   mesh.gIEN - element connectivity (num_nodes_per_elem, num_elems)
///   mesh.eType - element type
///
/// Replicates 'subroutine loadVTK(vtk,fName,istat)' defined in vtkXMLParser.f90.
//
void load_vtu(const std::string& file_name, mshType& mesh)
{
  #define n_debug_load_vtu 
  #ifdef debug_load_vtu 
  std::cout << "[load_vtu] " << std::endl;
  std::cout << "[load_vtu] ===== vtk_xml_parser::load_vtu ===== " << std::endl;
  std::cout << "[load_vtu] file_name: " << file_name << std::endl;
  #endif

  auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(file_name.c_str());
  reader->Update();
  vtkSmartPointer<vtkUnstructuredGrid> vtk_ugrid = reader->GetOutput();

  vtkIdType num_nodes = vtk_ugrid->GetNumberOfPoints();
  if (num_nodes == 0) {
    throw std::runtime_error("Failed reading the VTK file '" + file_name + "'.");
  }

  vtkIdType num_elems = vtk_ugrid->GetNumberOfCells();
  auto cell = vtkGenericCell::New();
  vtk_ugrid->GetCell(0, cell);
  int np_elem = cell->GetNumberOfPoints();
  int elem_type = cell->GetCellType();

  // Set the mesh element type from the vtk cell type.
  try {
    mesh.eType = vtk_cell_to_elem_type[elem_type];
  } catch (const std::bad_function_call& exception) {
    throw std::runtime_error("[load_vtu] The VTK cell type " + std::to_string(elem_type) + " is not supported.");
  }

  #ifdef debug_load_vtu 
  std::cout << "[load_vtu] Number of nodes: " << num_nodes << std::endl;
  std::cout << "[load_vtu] Number of elements: " << num_elems << std::endl;
  std::cout << "[load_vtu] Number of nodes per element: " << np_elem << std::endl;
  std::cout << "[load_vtu] Element type: " << elem_type << std::endl;
  #endif

  // Store nodal coordinates.
  auto points = vtk_ugrid->GetPoints();
  store_nodal_coords(points, mesh);

  // Store nodal IDs.
  store_nodal_ids(vtk_ugrid, mesh);

  // Store element connectivity.
  store_element_conn(vtk_ugrid, mesh);
}

/// @brief Store a surface mesh read from a VTK .vtu file into a Face object.
//
void load_vtu(const std::string& file_name, faceType& face)
{
  #define n_debug_load_vtu 
  #ifdef debug_load_vtu 
  std::cout << "[load_vtu] " << std::endl;
  std::cout << "[load_vtu] ===== vtk_xml_parser::load_vtu ===== " << std::endl;
  std::cout << "[load_vtu] file_name: " << file_name << std::endl;
  #endif

  auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(file_name.c_str());
  reader->Update();
  vtkSmartPointer<vtkUnstructuredGrid> vtk_ugrid = reader->GetOutput();

  vtkIdType num_nodes = vtk_ugrid->GetNumberOfPoints();
  if (num_nodes == 0) {
    throw std::runtime_error("Failed reading the VTK file '" + file_name + "'.");
  }

  vtkIdType num_elems = vtk_ugrid->GetNumberOfCells();
  auto cell = vtkGenericCell::New();
  vtk_ugrid->GetCell(0, cell);
  int np_elem = cell->GetNumberOfPoints();
  int elem_type = cell->GetCellType();

  // Set the mesh element type from the vtk cell type.
  try {
    face.eType = vtk_cell_to_elem_type[elem_type];
  } catch (const std::bad_function_call& exception) {
    throw std::runtime_error("[load_vtu] The VTK cell type " + std::to_string(elem_type) + " is not supported.");
  }

  #ifdef debug_load_vtu 
  std::cout << "[load_vtu] Number of nodes: " << num_nodes << std::endl;
  std::cout << "[load_vtu] Number of elements: " << num_elems << std::endl;
  std::cout << "[load_vtu] Number of nodes per element: " << np_elem << std::endl;
  #endif

  // Store nodal coordinates.
  auto points = vtk_ugrid->GetPoints();
  store_nodal_coords(points, face);

  // Store nodal IDs.
  store_nodal_ids(vtk_ugrid, face);

  // Store element connectivity.
  store_element_conn(vtk_ugrid, face);

  // Store element IDs.
  store_element_ids(vtk_ugrid, face);
}


/// @brief Read a time series field from a VTK .vtu file.
///
/// Mesh variables set
///   mesh.Ys - time series field data (num_components, num_nodes, num_time_steps)
///
//
void load_time_varying_field_vtu(const std::string file_name, const std::string field_name, mshType& mesh)
{
  #define n_debug_load_vtu
  #ifdef debug_load_vtu
  std::cout << "[load_vtu] " << std::endl;
  std::cout << "[load_vtu] ===== vtk_xml_parser::load_time_varying_field_vtu ===== " << std::endl;
  std::cout << "[load_vtu] file_name: " << file_name << std::endl;
  #endif

    auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(file_name.c_str());
    reader->Update();
    vtkSmartPointer<vtkUnstructuredGrid> vtk_ugrid = reader->GetOutput();
    vtkIdType num_nodes = vtk_ugrid->GetNumberOfPoints();
    int array_count = 0;
    std::vector<std::pair<std::string, int>> array_names;

    if (num_nodes == 0) {
      throw std::runtime_error("Failed reading the VTK file '" + file_name + "'.");
    }
    // Store all array names
    for (int i = 0; i < vtk_ugrid->GetPointData()->GetNumberOfArrays(); i++) {
      std::string array_name = vtk_ugrid->GetPointData()->GetArrayName(i);
      size_t pos = array_name.find(field_name.c_str());
      if (pos != std::string::npos) {
        auto not_digit = [](char c) { return !std::isdigit(c); };
        auto it = std::find_if(array_name.rbegin(), array_name.rend(), not_digit);
        std::string time_step = std::string(it.base(), array_name.end());
        array_count++;
        if (!time_step.empty()) {
          array_names.push_back({array_name, std::stoi(time_step)});
        } else {
          array_names.push_back({array_name, 0});
        }
      }
    }
    // Check if there are any fields present in the VTK file
    if (array_count == 0) {
      throw std::runtime_error("No '" + field_name + "' data found in the VTK file '" + file_name + "'.");
    }

    // Order all array names by time step
    std::sort(array_names.begin(), array_names.end(), [](const std::pair<std::string, int>& a, const std::pair<std::string, int>& b) {
      return a.second < b.second;
    });
    // Get the expected number of state-variable components
    int num_components = vtk_ugrid->GetPointData()->GetArray(array_names[0].first.c_str())->GetNumberOfComponents();
    mesh.Ys.resize(num_components, num_nodes, array_count);

    for (int i = 0; i < array_count; i++) {
      auto array = vtk_ugrid->GetPointData()->GetArray(array_names[i].first.c_str());
      if (array == nullptr) {
        throw std::runtime_error("No '" + array_names[i].first + "' data found in the VTK file '" + file_name + "'.");
      }
      if (array->GetNumberOfComponents() != num_components) {
        throw std::runtime_error("The number of components in the field '" + array_names[i].first + "' is not equal to the number of components in the first field.");
      }
      for (int j = 0; j < num_nodes; j++) {
        for (int k = 0; k < num_components; k++) {
          mesh.Ys(k, j, i) = array->GetComponent(j, k);
        }
      }
    }
}
} // namespace vtk_utils

