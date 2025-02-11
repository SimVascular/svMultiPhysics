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

#include "VtkData.h"
#include "Array.h"

#include <vtkDoubleArray.h>
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include <vtkGenericCell.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <string>
#include <map>

/////////////////////////////////////////////////////////////////
//        I n t e r n a l   I m p l e m e n t a t i o n        //
/////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////
//                    V t k V t p D a t a                      //
/////////////////////////////////////////////////////////////////

class VtkVtpData::VtkVtpDataImpl {
  public:
    VtkVtpDataImpl(); 
    void read_file(const std::string& file_name);
    void set_connectivity(const int nsd, const Array<int>& conn, const int pid);
    void set_point_data(const std::string& data_name, const Vector<int>& data);
    void set_points(const Array<double>& points);
    void write(const std::string& file_name);

    vtkSmartPointer<vtkPolyData> vtk_polydata;
    int elem_type;
    int num_elems;
    int np_elem;
    int num_points;
};

VtkVtpData::VtkVtpDataImpl::VtkVtpDataImpl()
{
  vtk_polydata = vtkSmartPointer<vtkPolyData>::New();
}

void VtkVtpData::VtkVtpDataImpl::read_file(const std::string& file_name)
{
  auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(file_name.c_str());
  reader->Update();
  vtk_polydata= reader->GetOutput();
  num_elems = vtk_polydata->GetNumberOfCells();
  num_points = vtk_polydata->GetPoints()->GetNumberOfPoints();
  if (num_points == 0) {
    throw std::runtime_error("Error reading the VTK VTP file '" + file_name + "'.");
  }

  // Get the cell type.
  auto cell = vtkGenericCell::New();
  vtk_polydata->GetCell(0, cell);
  np_elem = cell->GetNumberOfPoints();
  elem_type = cell->GetCellType();
}

void VtkVtpData::VtkVtpDataImpl::set_connectivity(const int nsd, const Array<int>& conn, const int pid)
{
  //std::cout << "[VtkVtpData.set_connectivity] " << std::endl;
  //std::cout << "[VtkVtpData.set_connectivity] vtk_polydata: " << vtk_polydata << std::endl;
  //std::cout << "[VtkVtpData.set_connectivity] nsd: " << nsd << std::endl;
  int num_elems = conn.ncols();
  int np_elem = conn.nrows();
  unsigned char vtk_cell_type;
  //std::cout << "[VtkVtpData.set_connectivity] num_elems: " << num_elems << std::endl;
  //std::cout << "[VtkVtpData.set_connectivity] np_elem: " << np_elem << std::endl;

  if (nsd == 2) {
    if (np_elem == 4) {
      vtk_cell_type = VTK_QUAD;
    }
    else if (np_elem == 3) {
      vtk_cell_type = VTK_TRIANGLE;
    }
    else if (np_elem == 6) {
        vtk_cell_type = VTK_QUADRATIC_TRIANGLE;
    }
    else if (np_elem == 8) {
        vtk_cell_type = VTK_QUADRATIC_QUAD;
    }
    else if (np_elem == 9) {
        vtk_cell_type = VTK_BIQUADRATIC_QUAD;
    }

  } else if (nsd == 3) {
    if (np_elem == 4) {
      vtk_cell_type = VTK_QUAD;
    }
    else if (np_elem == 3) {
      vtk_cell_type = VTK_TRIANGLE;
      //std::cout << "[VtkVtpData.set_connectivity] vtk_cell_type = VTK_TRIANGLE " << std::endl;
    }
    else if (np_elem == 6) {
        vtk_cell_type = VTK_QUADRATIC_TRIANGLE;
    }
    else if (np_elem == 8) {
        vtk_cell_type = VTK_HEXAHEDRON;
    }
    else if (np_elem == 10) {
        vtk_cell_type = VTK_QUADRATIC_TETRA;
    }
    else if (np_elem == 20) {
        vtk_cell_type = VTK_QUADRATIC_HEXAHEDRON;
    }
    else if (np_elem == 27) {
        vtk_cell_type = VTK_TRIQUADRATIC_HEXAHEDRON;
    }
  }
  
  if (np_elem == 2) {
    vtk_cell_type = VTK_LINE;
  }

  auto elem_nodes = vtkSmartPointer<vtkIdList>::New();
  elem_nodes->Allocate(np_elem);
  elem_nodes->Initialize();
  elem_nodes->SetNumberOfIds(np_elem);

  auto elem_ids = vtkSmartPointer<vtkIntArray>::New();
  elem_ids->SetNumberOfComponents(1);
  elem_ids->Allocate(num_elems,1000);
  elem_ids->SetNumberOfTuples(num_elems);
  elem_ids->SetName("GlobalElementID");

  vtkSmartPointer<vtkCellArray> element_cells = vtkSmartPointer<vtkCellArray>::New();

  for (int i = 0; i < num_elems; i++) {
    //std::cout << "[VtkVtpData.set_connectivity] ---------- i " << i << std::endl;
    for (int j = 0; j < np_elem; j++) {
      //std::cout << "[VtkVtpData.set_connectivity] ----- j " << j << std::endl;
      //std::cout << "[VtkVtpData.set_connectivity] conn(j,i): " << conn(j,i) << std::endl;
      elem_nodes->SetId(j, conn(j,i));
    }
    element_cells->InsertNextCell(elem_nodes);
    //vtk_polydata->InsertNextCell(vtk_cell_type, elem_nodes);
    elem_ids->SetTuple1(i,i+1);
  }

  vtk_polydata->SetPolys(element_cells);
  vtk_polydata->GetCellData()->AddArray(elem_ids);
}

void VtkVtpData::VtkVtpDataImpl::set_point_data(const std::string& data_name, const Vector<int>& data)
{
  int num_vals = data.size();
  auto data_array = vtkSmartPointer<vtkIntArray>::New();
  data_array->SetNumberOfComponents(1);
  data_array->Allocate(num_vals);
  data_array->SetName(data_name.c_str());

  for (int i = 0; i < num_vals; i++) {
    data_array->InsertNextTuple1(data(i));
  }

  vtk_polydata->GetPointData()->AddArray(data_array);
}

/// @brief Set the 3D point (coordinate) data for the polydata.
//
void VtkVtpData::VtkVtpDataImpl::set_points(const Array<double>& points)
{
  //std::cout << "[VtkVtpData.set_points] vtk_polydata: " << vtk_polydata << std::endl;
  int num_coords = points.ncols();
  if (num_coords == 0) { 
    throw std::runtime_error("Error in VTK VTP set_points: the number of points is zero.");
  }
  //std::cout << "[VtkVtpData.set_points] num_coords: " << num_coords << std::endl;

  auto node_coords = vtkSmartPointer<vtkPoints>::New();
  node_coords->Allocate(num_coords, 1000);
  node_coords->SetNumberOfPoints(num_coords);

  auto node_ids = vtkSmartPointer<vtkIntArray>::New();
  node_ids->SetNumberOfComponents(1);
  node_ids->Allocate(num_coords,1000);
  node_ids->SetNumberOfTuples(num_coords);
  node_ids->SetName("GlobalNodeID");

  for (int i = 0; i < num_coords; i++ ) {
    node_coords->SetPoint(i, points(0,i), points(1,i), points(2,i));
    node_ids->SetTuple1(i,i+1);
  }

  vtk_polydata->SetPoints(node_coords);
  vtk_polydata->GetPointData()->AddArray(node_ids);
}

void VtkVtpData::VtkVtpDataImpl::write(const std::string& file_name)
{
  //std::cout << "[VtkVtpData.write] file_name: " << file_name << std::endl;
  auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetInputDataObject(vtk_polydata);
  writer->SetFileName(file_name.c_str());
  writer->Write();
}

/////////////////////////////////////////////////////////////////
//                    V t k V t u D a t a                      //
/////////////////////////////////////////////////////////////////

class VtkVtuData::VtkVtuDataImpl {
  public:
    void create_grid();
    void read_file(const std::string& file_name);
    void set_connectivity(const int nsd, const Array<int>& conn, const int pid);

    void set_element_data(const std::string& data_name, const Array<double>& data);
    void set_element_data(const std::string& data_name, const Array<int>& data);

    void set_point_data(const std::string& data_name, const Array<double>& data);
    void set_point_data(const std::string& data_name, const Array<int>& data);
    void set_point_data(const std::string& data_name, const Vector<int>& data);

    void set_points(const Array<double>& points);
    void write(const std::string& file_name);

    template<typename T1, typename T2>
    void set_element_data(const std::string& data_name, const T1& data, T2& data_array)
    {
      int num_vals = data.ncols();
      int num_comp = data.nrows();
      data_array->SetNumberOfComponents(num_comp);
      data_array->Allocate(num_vals,1000);
      data_array->SetNumberOfTuples(num_vals);
      data_array->SetName(data_name.c_str());
      for (int i = 0; i < num_vals; i++) {
        for (int j = 0; j < num_comp; j++) {
          data_array->SetComponent(i, j, data(j,i));
        }
      }
      vtk_ugrid->GetCellData()->AddArray(data_array);
    };

    vtkSmartPointer<vtkUnstructuredGrid> vtk_ugrid;
    int elem_type;
    int num_elems;
    int np_elem;
    int num_points;
};

void VtkVtuData::VtkVtuDataImpl::create_grid()
{
  vtk_ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
}

void VtkVtuData::VtkVtuDataImpl::read_file(const std::string& file_name)
{
  auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(file_name.c_str());
  reader->Update();
  vtk_ugrid = reader->GetOutput();
  num_elems = vtk_ugrid->GetNumberOfCells();
  num_points = vtk_ugrid->GetPoints()->GetNumberOfPoints();
  if (num_points == 0) {
    throw std::runtime_error("Error reading the VTK VTU file '" + file_name + "'.");
  }

  // Get the cell type.
  auto cell = vtkGenericCell::New();
  vtk_ugrid->GetCell(0, cell);
  np_elem = cell->GetNumberOfPoints();
  elem_type = cell->GetCellType();
}

void VtkVtuData::VtkVtuDataImpl::set_connectivity(const int nsd, const Array<int>& conn, const int pid)
{
  int num_elems = conn.ncols();
  int np_elem = conn.nrows();
  int num_coords = vtk_ugrid->GetPoints()->GetNumberOfPoints();
  unsigned char vtk_cell_type;
  /*
  std::cout << "[VtkVtuData.set_connectivity] " << std::endl;
  std::cout << "[VtkVtuData.set_connectivity] nsd: " << nsd << std::endl;
  std::cout << "[VtkVtuData.set_connectivity] num_elems: " << num_elems << std::endl;
  std::cout << "[VtkVtuData.set_connectivity] np_elem: " << np_elem << std::endl;
  std::cout << "[VtkVtuData.set_connectivity] num_coords: " << num_coords << std::endl;
  */

  if (nsd == 2) {

    if (np_elem == 3) {
      vtk_cell_type = VTK_TRIANGLE;

    } else if (np_elem == 4) {
      vtk_cell_type = VTK_QUAD;

    } else if (np_elem == 6) {
      vtk_cell_type = VTK_QUADRATIC_TRIANGLE;

    } else if (np_elem == 8) {
      vtk_cell_type = VTK_QUADRATIC_QUAD;

    } else if (np_elem == 9) {
        vtk_cell_type = VTK_BIQUADRATIC_QUAD;
    }

  } else if (nsd == 3) {

    if (np_elem == 3) {
      vtk_cell_type = VTK_TRIANGLE;

    } else if (np_elem == 4) {
      vtk_cell_type = VTK_TETRA;

    } else if (np_elem == 6) {
      vtk_cell_type = VTK_WEDGE;

    } else if (np_elem == 8) {
      vtk_cell_type = VTK_HEXAHEDRON;

    } else if (np_elem == 10) {
      vtk_cell_type = VTK_QUADRATIC_TETRA;

    } else if (np_elem == 20) {
      vtk_cell_type = VTK_QUADRATIC_HEXAHEDRON;

    } else if (np_elem == 27) {
      vtk_cell_type = VTK_TRIQUADRATIC_HEXAHEDRON;
    }

  }
  
  if (np_elem == 2) {
    vtk_cell_type = VTK_LINE;
  }

  auto elem_nodes = vtkSmartPointer<vtkIdList>::New();
  elem_nodes->Allocate(np_elem);
  elem_nodes->Initialize();
  elem_nodes->SetNumberOfIds(np_elem);
  //std::cout << "[VtkVtuData.set_connectivity] Set conn ... " << std::endl;

  for (int i = 0; i < num_elems; i++) {
    //std::cout << "[VtkVtuData.set_connectivity] " << i << ": ";
    for (int j = 0; j < np_elem; j++) {
      //std::cout << conn(j,i) << " "; 
      if ((conn(j,i) < 0) || (conn(j,i) >= num_coords)) {
        throw std::runtime_error("[VtkVtuData.set_connectivity] Element " + std::to_string(i+1) +
            " has the non-valid node ID " + std::to_string(conn(j,i)) + ".");
      }
      elem_nodes->SetId(j, conn(j,i));
    }
    //std::cout << std::endl; 
    vtk_ugrid->InsertNextCell(vtk_cell_type, elem_nodes);
  }
}

//------------------
// set_element_data
//------------------
//
/*
void VtkVtuData::VtkVtuDataImpl::set_element_data(const std::string& data_name, const Array<double>& data)
{
  int num_vals = data.num_cols();
  int num_comp = data.num_rows();

  auto data_array = vtkSmartPointer<vtkDoubleArray>::New();
  data_array->SetNumberOfComponents(num_comp);
  data_array->Allocate(num_vals,1000);
  data_array->SetNumberOfTuples(num_vals);
  data_array->SetName(data_name.c_str());

  for (int i = 0; i < num_vals; i++) {
    for (int j = 0; j < num_comp; j++) {
      data_array->SetComponent(i, j, data(j,i));
    }
  }

  vtk_ugrid->GetCellData()->AddArray(data_array);
}

//------------------
// set_element_data
//------------------
//
void VtkVtuData::VtkVtuDataImpl::set_element_data(const std::string& data_name, const Array<int>& data)
{
  int num_vals = data.num_cols();
  int num_comp = data.num_rows();

  auto data_array = vtkSmartPointer<vtkIntArray>::New();
  data_array->SetNumberOfComponents(num_comp);
  data_array->Allocate(num_vals,1000);
  data_array->SetNumberOfTuples(num_vals);
  data_array->SetName(data_name.c_str());

  for (int i = 0; i < num_vals; i++) {
    for (int j = 0; j < num_comp; j++) {
      data_array->SetComponent(i, j, data(j,i));
    }
  }

  vtk_ugrid->GetCellData()->AddArray(data_array);
}
*/

//----------------
// set_point_data
//----------------
//
void VtkVtuData::VtkVtuDataImpl::set_point_data(const std::string& data_name, const Array<double>& data)
{
  int num_vals = data.ncols();
  int num_comp = data.nrows();

  auto data_array = vtkSmartPointer<vtkDoubleArray>::New();
  data_array->SetNumberOfComponents(num_comp);
  data_array->Allocate(num_vals,1000);
  data_array->SetNumberOfTuples(num_vals);
  data_array->SetName(data_name.c_str());

  for (int i = 0; i < num_vals; i++) {
    for (int j = 0; j < num_comp; j++) {
      data_array->SetComponent(i, j, data(j,i));
    }
  }

  vtk_ugrid->GetPointData()->AddArray(data_array);
}

void VtkVtuData::VtkVtuDataImpl::set_point_data(const std::string& data_name, const Array<int>& data)
{
  int num_vals = data.ncols();
  int num_comp = data.nrows();

  auto data_array = vtkSmartPointer<vtkIntArray>::New();
  data_array->SetNumberOfComponents(num_comp);
  data_array->Allocate(num_vals,1000);
  data_array->SetNumberOfTuples(num_vals);
  data_array->SetName(data_name.c_str());

  for (int i = 0; i < num_vals; i++) {
    for (int j = 0; j < num_comp; j++) {
      data_array->SetComponent(i, j, data(j,i));
    }
  }

  vtk_ugrid->GetPointData()->AddArray(data_array);
}

void VtkVtuData::VtkVtuDataImpl::set_point_data(const std::string& data_name, const Vector<int>& data)
{
  throw std::runtime_error("[VtkVtuData] set_point_data for Vector<int> not implemented.");
}

/// @brief Set the 3D points (coordinates) data for the unstructured grid.
//
void VtkVtuData::VtkVtuDataImpl::set_points(const Array<double>& points)
{
  int num_coords = points.ncols();
  auto node_coords = vtkSmartPointer<vtkPoints>::New();
  node_coords->Allocate(num_coords ,1000);
  node_coords->SetNumberOfPoints(num_coords);
  //std::cout << "[VtkVtuData.set_points] " << std::endl;
  //std::cout << "[VtkVtuData.set_points] num_coords: " << num_coords << std::endl;

  for (int i = 0; i < num_coords; i++ ) {
    node_coords->SetPoint(i, points(0,i), points(1,i), points(2,i)); 
  }

  vtk_ugrid->SetPoints(node_coords);
}

void VtkVtuData::VtkVtuDataImpl::write(const std::string& file_name)
{
  auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetInputDataObject(vtk_ugrid);
  writer->SetFileName(file_name.c_str());
  writer->Write();
}

/////////////////////////////////////////////////////////////////
//          V t k D a t a     I m p l e m e n t a t i o n      //
/////////////////////////////////////////////////////////////////

VtkData::VtkData()
{
}

VtkData::~VtkData()
{
}

VtkData* VtkData::create_reader(const std::string& file_name)
{
  auto file_ext = file_name.substr(file_name.find_last_of(".") + 1);
  if (file_ext == "vtp") {
    return new VtkVtpData(file_name);
  } else if (file_ext == "vtu") {
    return new VtkVtuData(file_name);
  }
}

VtkData* VtkData::create_writer(const std::string& file_name)
{
  auto file_ext = file_name.substr(file_name.find_last_of(".") + 1);
  bool reader = false;
  if (file_ext == "vtp") {
    return new VtkVtpData(file_name, reader);
  } else if (file_ext == "vtu") {
    return new VtkVtuData(file_name, reader);
  }
}

//void VtkData::write(const std::string& file_name)
//{
//}

/////////////////////////////////////////////////////////////////
//      V t k V t p D a t a     I m p l e m e n t a t i o n    //
/////////////////////////////////////////////////////////////////

VtkVtpData::VtkVtpData()
{
  impl = new VtkVtpDataImpl;
}

VtkVtpData::VtkVtpData(const std::string& file_name, bool reader)
{
  this->file_name = file_name;
  impl = new VtkVtpDataImpl;
  if (reader) {
    read_file(file_name); 
   }
}

VtkVtpData::~VtkVtpData()
{
  delete impl;
}

Array<int> VtkVtpData::get_connectivity()
{
  int num_elems = impl->num_elems; 
  int np_elem = impl->np_elem; 
  Array<int> conn(np_elem, num_elems);

  auto cell = vtkGenericCell::New();
  for (int i = 0; i < num_elems; i++) {
    impl->vtk_polydata->GetCell(i, cell);
    auto num_cell_pts = cell->GetNumberOfPoints();
    for (int j = 0; j < num_cell_pts; j++) {
      auto id = cell->PointIds->GetId(j);
      conn(j,i) = id;
    }
  }

  return  conn;
}

/// @brief Copy an array of point data from an polydata mesh into the given Array.
//
void VtkVtpData::copy_point_data(const std::string& data_name, Array<double>& mesh_data)
{
  auto vtk_data = vtkDoubleArray::SafeDownCast(impl->vtk_polydata->GetPointData()->GetArray(data_name.c_str()));
  if (vtk_data == nullptr) { 
    return;
  }

  int num_data = vtk_data->GetNumberOfTuples();
  if (num_data == 0) { 
    return; 
  }

  int num_comp = vtk_data->GetNumberOfComponents();

  // Set the data.
  for (int i = 0; i < num_data; i++) {
    auto tuple = vtk_data->GetTuple(i);
    for (int j = 0; j < num_comp; j++) {
      mesh_data(j, i) = tuple[j];
    }
  }
}

void VtkVtpData::copy_point_data(const std::string& data_name, Vector<double>& mesh_data)
{
  auto vtk_data = vtkDoubleArray::SafeDownCast(impl->vtk_polydata->GetPointData()->GetArray(data_name.c_str()));
  if (vtk_data == nullptr) { 
    return;
  }

  int num_data = vtk_data->GetNumberOfTuples();
  if (num_data == 0) { 
    return; 
  }

  int num_comp = vtk_data->GetNumberOfComponents();

  // Set the data.
  for (int i = 0; i < num_data; i++) {
    mesh_data(i) = vtk_data->GetValue(i);
  }
}

void VtkVtpData::copy_point_data(const std::string& data_name, Vector<int>& mesh_data)
{
  auto vtk_data = vtkIntArray::SafeDownCast(impl->vtk_polydata->GetPointData()->GetArray(data_name.c_str()));
  if (vtk_data == nullptr) {
    return;
  }

  int num_data = vtk_data->GetNumberOfTuples();
  if (num_data == 0) {
    return;
  }

  int num_comp = vtk_data->GetNumberOfComponents();

  // Set the data.
  for (int i = 0; i < num_data; i++) {
    mesh_data(i) = vtk_data->GetValue(i);
  }
}

/// @brief Copy points into the given array.
//
void VtkVtpData::copy_points(Array<double>& points)
{
  auto vtk_points = impl->vtk_polydata->GetPoints();
  auto num_points = vtk_points->GetNumberOfPoints();
  Array<double> points_array(3, num_points);

  double point[3];
  for (int i = 0; i < num_points; i++) {
    vtk_points->GetPoint(i, point);
    points(0,i) = point[0];
    points(1,i) = point[1];
    points(2,i) = point[2];
  }

  return;
}

/// @brief Get an array of point data from an unstructured grid.
//
Array<double> VtkVtpData::get_point_data(const std::string& data_name)
{
  auto vtk_data = vtkDoubleArray::SafeDownCast(impl->vtk_polydata->GetPointData()->GetArray(data_name.c_str()));
  if (vtk_data == nullptr) { 
    return Array<double>();
  }

  int num_data = vtk_data->GetNumberOfTuples();
  if (num_data == 0) { 
    return Array<double>();
  }

  int num_comp = vtk_data->GetNumberOfComponents();

  // Set the data.
  Array<double> data(num_data, num_comp);
  for (int i = 0; i < num_data; i++) {
    auto tuple = vtk_data->GetTuple(i);
    for (int j = 0; j < num_comp; j++) {
      data(i, j) = tuple[j];
    }
  }

  return data;
}

/// @brief Get a list of point data names.
std::vector<std::string> VtkVtpData::get_point_data_names()
{
  std::vector<std::string> data_names; 
  int num_arrays = impl->vtk_polydata->GetPointData()->GetNumberOfArrays();

  for (int i = 0; i < num_arrays; i++) {
    auto array_name = impl->vtk_polydata->GetPointData()->GetArrayName(i);
    data_names.push_back(array_name); 
  }

  return data_names; 
}

/// @brief Get an array of point data from an unstructured grid.
//
Array<double> VtkVtpData::get_points()
{
  auto vtk_points = impl->vtk_polydata->GetPoints();
  auto num_points = vtk_points->GetNumberOfPoints();
  Array<double> points_array(3, num_points);

  double point[3];
  for (int i = 0; i < num_points; i++) {
    vtk_points->GetPoint(i, point);
    points_array(0,i) = point[0];
    points_array(1,i) = point[1];
    points_array(2,i) = point[2];
  }

  return points_array;
}

bool VtkVtpData::has_point_data(const std::string& data_name)
{
  int num_arrays = impl->vtk_polydata->GetPointData()->GetNumberOfArrays();

  for (int i = 0; i < num_arrays; i++) {
    if (!strcmp(impl->vtk_polydata->GetPointData()->GetArrayName(i), data_name.c_str())) {
      return true;
    }
  }

  return false;
}

int VtkVtpData::elem_type() 
{ 
  return impl->elem_type; 
}

int VtkVtpData::num_elems() 
{ 
  return impl->num_elems; 
}

int VtkVtpData::np_elem() 
{ 
  return impl->np_elem; 
}

int VtkVtpData::num_points() 
{ 
  return impl->num_points; 
}

void VtkVtpData::read_file(const std::string& file_name) 
{ 
  impl->read_file(file_name); 
}

void VtkVtpData::set_connectivity(const int nsd, const Array<int>& conn, const int pid)
{
  impl->set_connectivity(nsd, conn, pid);
}

void VtkVtpData::set_element_data(const std::string& data_name, const Array<double>& data)
{
  throw std::runtime_error("[VtkVtpData] set_element_data not implemented.");
}

void VtkVtpData::set_element_data(const std::string& data_name, const Array<int>& data) 
{
  throw std::runtime_error("[VtkVtpData] set_element_data not implemented.");
}

void VtkVtpData::set_point_data(const std::string& data_name, const Array<double>& data)
{
  throw std::runtime_error("[VtkVtpData] set_point_data for Array<double> not implemented.");
}

void VtkVtpData::set_point_data(const std::string& data_name, const Array<int>& data)
{
  throw std::runtime_error("[VtkVtpData] set_point_data Array<int> not implemented.");
}

void VtkVtpData::set_point_data(const std::string& data_name, const Vector<int>& data)
{
  impl->set_point_data(data_name, data);
}

void VtkVtpData::set_points(const Array<double>& points)
{
  impl->set_points(points);
}

void VtkVtpData::write()
{
  impl->write(file_name);
}

/////////////////////////////////////////////////////////////////
//      V t k V t u D a t a     I m p l e m e n t a t i o n    //
/////////////////////////////////////////////////////////////////


VtkVtuData::VtkVtuData()
{
  impl = new VtkVtuDataImpl;
}

VtkVtuData::VtkVtuData(const std::string& file_name, bool reader)
{
  this->file_name = file_name;
  impl = new VtkVtuDataImpl;
  if (reader) {
    read_file(file_name); 
  } else {
    impl->create_grid();
  }
}

VtkVtuData::~VtkVtuData()
{
  delete impl;
}

Array<int> VtkVtuData::get_connectivity()
{
  int num_elems = impl->num_elems; 
  int np_elem = impl->np_elem; 

  Array<int> conn(np_elem, num_elems);

  auto cell = vtkGenericCell::New();
  for (int i = 0; i < num_elems; i++) {
    impl->vtk_ugrid->GetCell(i, cell);
    auto num_cell_pts = cell->GetNumberOfPoints();
    for (int j = 0; j < num_cell_pts; j++) {
      auto id = cell->PointIds->GetId(j);
      conn(j,i) = id;
    }
  }
  return conn;
}

/// @brief Get a list of point data names.
std::vector<std::string> VtkVtuData::get_point_data_names()
{
  std::vector<std::string> data_names;
  int num_arrays = impl->vtk_ugrid->GetPointData()->GetNumberOfArrays();

  for (int i = 0; i < num_arrays; i++) {
    auto array_name = impl->vtk_ugrid->GetPointData()->GetArrayName(i);
    data_names.push_back(array_name);
  }

  return data_names;
}

/// @brief Copy an array of point data from an unstructured grid into the given Array.
//
void VtkVtuData::copy_point_data(const std::string& data_name, Array<double>& mesh_data)
{

  auto vtk_data = vtkDoubleArray::SafeDownCast(impl->vtk_ugrid->GetPointData()->GetArray(data_name.c_str()));
  if (vtk_data == nullptr) { 
    return;
  }

  int num_data = vtk_data->GetNumberOfTuples();
  if (num_data == 0) { 
    return; 
  }

  int num_comp = vtk_data->GetNumberOfComponents();

  // Set the data.
  for (int i = 0; i < num_data; i++) {
    auto tuple = vtk_data->GetTuple(i);
    for (int j = 0; j < num_comp; j++) {
      mesh_data(j, i) = tuple[j];
    }
  }
}

void VtkVtuData::copy_point_data(const std::string& data_name, Vector<double>& mesh_data)
{
  auto vtk_data = vtkDoubleArray::SafeDownCast(impl->vtk_ugrid->GetPointData()->GetArray(data_name.c_str()));
  if (vtk_data == nullptr) {
    return;
  }

  int num_data = vtk_data->GetNumberOfTuples();
  if (num_data == 0) {
    return;
  }

  int num_comp = vtk_data->GetNumberOfComponents();

  // Set the data.
  for (int i = 0; i < num_data; i++) {
    mesh_data[i] = vtk_data->GetValue(i);
  }
}

void VtkVtuData::copy_point_data(const std::string& data_name, Vector<int>& mesh_data)
{
  auto vtk_data = vtkIntArray::SafeDownCast(impl->vtk_ugrid->GetPointData()->GetArray(data_name.c_str()));
  if (vtk_data == nullptr) {
    return;
  }

  int num_data = vtk_data->GetNumberOfTuples();
  if (num_data == 0) {
    return;
  }

  int num_comp = vtk_data->GetNumberOfComponents();

  // Set the data.
  for (int i = 0; i < num_data; i++) {
    mesh_data[i] = vtk_data->GetValue(i);
  }
}

/// @brief Copy points into the given array.
//
void VtkVtuData::copy_points(Array<double>& points)
{ 
  auto vtk_points = impl->vtk_ugrid->GetPoints();
  auto num_points = vtk_points->GetNumberOfPoints();
  Array<double> points_array(3, num_points);
  
  double point[3];
  for (int i = 0; i < num_points; i++) {
    vtk_points->GetPoint(i, point);
    points(0,i) = point[0];
    points(1,i) = point[1];
    points(2,i) = point[2];
  } 
  
  return;
}

bool VtkVtuData::has_point_data(const std::string& data_name)
{
  int num_arrays = impl->vtk_ugrid->GetPointData()->GetNumberOfArrays();

  for (int i = 0; i < num_arrays; i++) {
    if (!strcmp(impl->vtk_ugrid->GetPointData()->GetArrayName(i), data_name.c_str())) {
      return true;
    }
  }

  return false;
}

/// @brief Get an array of point data from an unstructured grid.
//
Array<double> VtkVtuData::get_point_data(const std::string& data_name)
{
  auto vtk_data = vtkDoubleArray::SafeDownCast(impl->vtk_ugrid->GetPointData()->GetArray(data_name.c_str()));
  if (vtk_data == nullptr) { 
    return Array<double>();
  }

  int num_data = vtk_data->GetNumberOfTuples();
  if (num_data == 0) { 
    return Array<double>();
  }

  int num_comp = vtk_data->GetNumberOfComponents();

  // Set the data.
  Array<double> data(num_data, num_comp);
  for (int i = 0; i < num_data; i++) {
    auto tuple = vtk_data->GetTuple(i);
    for (int j = 0; j < num_comp; j++) {
      data(i, j) = tuple[j];
    }
  }

  return data;
}

Array<double> VtkVtuData::get_points()
{
  auto vtk_points = impl->vtk_ugrid->GetPoints();
  auto num_points = vtk_points->GetNumberOfPoints();
  Array<double> points_array(3, num_points);

  double point[3];
  for (int i = 0; i < num_points; i++) {
    vtk_points->GetPoint(i, point);
    points_array(0,i) = point[0];
    points_array(1,i) = point[1];
    points_array(2,i) = point[2];
  }

  return points_array;
}

int VtkVtuData::elem_type() 
{ 
  return impl->elem_type; 
}

int VtkVtuData::num_elems() 
{ 
  return impl->num_elems; 
}

int VtkVtuData::np_elem() 
{ 
  return impl->np_elem; 
}

int VtkVtuData::num_points() 
{ 
  return impl->num_points; 
}


void VtkVtuData::read_file(const std::string& file_name) 
{ 
  impl->read_file(file_name); 
}

void VtkVtuData::set_connectivity(const int nsd, const Array<int>& conn, const int pid)
{
  impl->set_connectivity(nsd, conn, pid);
}

void VtkVtuData::set_element_data(const std::string& data_name, const Array<double>& data)
{
  auto data_array = vtkSmartPointer<vtkDoubleArray>::New();
  impl->set_element_data(data_name, data, data_array); 
  //impl->set_element_data(data_name, data);
}

void VtkVtuData::set_element_data(const std::string& data_name, const Array<int>& data)
{
  auto data_array = vtkSmartPointer<vtkIntArray>::New();
  impl->set_element_data(data_name, data, data_array); 
  //impl->set_element_data(data_name, data);
}

void VtkVtuData::set_point_data(const std::string& data_name, const Array<double>& data)
{
  impl->set_point_data(data_name, data);
}

void VtkVtuData::set_point_data(const std::string& data_name, const Array<int>& data)
{
  impl->set_point_data(data_name, data);
}

void VtkVtuData::set_point_data(const std::string& data_name, const Vector<int>& data)
{
  impl->set_point_data(data_name, data);
}

void VtkVtuData::set_points(const Array<double>& points)
{
  impl->set_points(points);
}

void VtkVtuData::write()
{
  impl->write(file_name);
}

