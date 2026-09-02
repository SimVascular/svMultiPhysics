// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "VtkData.h"

#include "Core/Exception.h"
#include "FE/Common/FEException.h"

#include <cstring>
#include <string>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkErrorCode.h>
#include <vtkFieldData.h>
#include <vtkGenericCell.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

void VtkData::read_file(const std::string &file_name) {
  file_name_ = file_name;

  read_file_internal(file_name);

  // Extract metadata (number of cells, points, cell type).
  {
    // A failed read leaves the data object empty, without any points.
    const auto points = vtk_data->GetPoints();
    num_points_ = (points == nullptr) ? 0 : points->GetNumberOfPoints();
    svmp::check<svmp::FileFormatException>(num_points_ != 0, file_name,
                                           "The file has no points.");

    num_elems_ = vtk_data->GetNumberOfCells();
    svmp::check<svmp::FileFormatException>(num_elems_ != 0, file_name,
                                           "The file has no elements.");

    // Get the cell type.
    auto cell = vtkSmartPointer<vtkGenericCell>::New();
    vtk_data->GetCell(0, cell);
    num_points_per_elem_ = cell->GetNumberOfPoints();
    elem_type_ = cell->GetCellType();
  }
}

Array<int> VtkData::get_connectivity() const {
  Array<int> connectivity(num_points_per_elem_, num_elems_);

  auto cell = vtkGenericCell::New();
  for (int i = 0; i < num_elems_; i++) {
    vtk_data->GetCell(i, cell);

    const int num_cell_pts = cell->GetNumberOfPoints();

    for (int j = 0; j < num_cell_pts; ++j) {
      connectivity(j, i) = cell->PointIds->GetId(j);
    }
  }

  return connectivity;
}

Array<double> VtkData::get_points() const {
  auto vtk_points = vtk_data->GetPoints();
  auto num_points = vtk_points->GetNumberOfPoints();

  Array<double> points_array(3, num_points);

  double point[3];
  for (int i = 0; i < num_points; i++) {
    vtk_points->GetPoint(i, point);
    points_array(0, i) = point[0];
    points_array(1, i) = point[1];
    points_array(2, i) = point[2];
  }

  return points_array;
}

int VtkData::num_elems() const { return num_elems_; }

int VtkData::elem_type() const { return elem_type_; }

int VtkData::num_points_per_elem() const { return num_points_per_elem_; }

int VtkData::num_points() const { return num_points_; }

void VtkData::set_element_data(const std::string &data_name,
                               const Array<double> &data) {
  const int num_vals = data.ncols();
  const int num_components = data.nrows();

  svmp::check<svmp::FE::InvalidArgumentException>(
      num_vals == vtk_data->GetNumberOfCells(),
      "The element data array named '" + data_name + "' holds " +
          std::to_string(num_vals) + " values, while the mesh has " +
          std::to_string(vtk_data->GetNumberOfCells()) + " elements.");

  auto data_array = vtkSmartPointer<vtkDoubleArray>::New();
  data_array->SetNumberOfComponents(num_components);
  data_array->Allocate(num_vals, 1000);
  data_array->SetNumberOfTuples(num_vals);
  data_array->SetName(data_name.c_str());
  for (int i = 0; i < num_vals; ++i) {
    for (int j = 0; j < num_components; ++j) {
      data_array->SetComponent(i, j, data(j, i));
    }
  }

  vtk_data->GetCellData()->AddArray(data_array);
}

void VtkData::set_element_data(const std::string &data_name,
                               const Array<int> &data) {
  const int num_vals = data.ncols();
  const int num_components = data.nrows();

  svmp::check<svmp::FE::InvalidArgumentException>(
      num_vals == vtk_data->GetNumberOfCells(),
      "The element data array named '" + data_name + "' holds " +
          std::to_string(num_vals) + " values, while the mesh has " +
          std::to_string(vtk_data->GetNumberOfCells()) + " elements.");

  auto data_array = vtkSmartPointer<vtkIntArray>::New();
  data_array->SetNumberOfComponents(num_components);
  data_array->Allocate(num_vals, 1000);
  data_array->SetNumberOfTuples(num_vals);
  data_array->SetName(data_name.c_str());
  for (int i = 0; i < num_vals; ++i) {
    for (int j = 0; j < num_components; ++j) {
      data_array->SetComponent(i, j, data(j, i));
    }
  }

  vtk_data->GetCellData()->AddArray(data_array);
}

void VtkData::set_element_data(const std::string &data_name,
                               const Vector<int> &data) {
  const int num_vals = data.size();

  svmp::check<svmp::FE::InvalidArgumentException>(
      num_vals == vtk_data->GetNumberOfCells(),
      "The element data array named '" + data_name + "' holds " +
          std::to_string(num_vals) + " values, while the mesh has " +
          std::to_string(vtk_data->GetNumberOfCells()) + " elements.");

  auto data_array = vtkSmartPointer<vtkIntArray>::New();
  data_array->SetNumberOfComponents(1);
  data_array->Allocate(num_vals);
  data_array->SetName(data_name.c_str());

  for (int i = 0; i < num_vals; ++i) {
    data_array->InsertNextTuple1(data(i));
  }

  vtk_data->GetCellData()->AddArray(data_array);
}

void VtkData::set_point_data(const std::string &data_name,
                             const Array<double> &data) {
  const int num_vals = data.ncols();
  const int num_comp = data.nrows();

  svmp::check<svmp::FE::InvalidArgumentException>(
      num_vals == vtk_data->GetNumberOfPoints(),
      "The point data array named '" + data_name + "' holds " +
          std::to_string(num_vals) + " values, while the mesh has " +
          std::to_string(vtk_data->GetNumberOfPoints()) + " points.");

  auto data_array = vtkSmartPointer<vtkDoubleArray>::New();
  data_array->SetNumberOfComponents(num_comp);
  data_array->Allocate(num_vals, 1000);
  data_array->SetNumberOfTuples(num_vals);
  data_array->SetName(data_name.c_str());

  for (int i = 0; i < num_vals; i++) {
    for (int j = 0; j < num_comp; j++) {
      data_array->SetComponent(i, j, data(j, i));
    }
  }

  vtk_data->GetPointData()->AddArray(data_array);
}

void VtkData::set_point_data(const std::string &data_name,
                             const Array<int> &data) {
  int num_vals = data.ncols();
  int num_comp = data.nrows();

  svmp::check<svmp::FE::InvalidArgumentException>(
      num_vals == vtk_data->GetNumberOfPoints(),
      "The point data array named '" + data_name + "' holds " +
          std::to_string(num_vals) + " values, while the mesh has " +
          std::to_string(vtk_data->GetNumberOfPoints()) + " points.");

  auto data_array = vtkSmartPointer<vtkIntArray>::New();
  data_array->SetNumberOfComponents(num_comp);
  data_array->Allocate(num_vals, 1000);
  data_array->SetNumberOfTuples(num_vals);
  data_array->SetName(data_name.c_str());

  for (int i = 0; i < num_vals; i++) {
    for (int j = 0; j < num_comp; j++) {
      data_array->SetComponent(i, j, data(j, i));
    }
  }

  vtk_data->GetPointData()->AddArray(data_array);
}

void VtkData::set_point_data(const std::string &data_name,
                             const Vector<int> &data) {
  const int num_vals = data.size();

  svmp::check<svmp::FE::InvalidArgumentException>(
      num_vals == vtk_data->GetNumberOfPoints(),
      "The point data array named '" + data_name + "' holds " +
          std::to_string(num_vals) + " values, while the mesh has " +
          std::to_string(vtk_data->GetNumberOfPoints()) + " points.");

  auto data_array = vtkSmartPointer<vtkIntArray>::New();
  data_array->SetNumberOfComponents(1);
  data_array->Allocate(num_vals);
  data_array->SetName(data_name.c_str());

  for (int i = 0; i < num_vals; i++) {
    data_array->InsertNextTuple1(data(i));
  }

  vtk_data->GetPointData()->AddArray(data_array);
}

void VtkData::set_points(const Array<double> &points) {
  const int num_coords = points.ncols();
  svmp::check<svmp::FE::InvalidArgumentException>(
      num_coords != 0, "The number of points is zero.");

  svmp::check<svmp::FE::InvalidArgumentException>(
      points.nrows() >= 3,
      "The point coordinates are given as an array of " +
          std::to_string(points.nrows()) +
          " rows, while three coordinates per point are needed.");

  auto node_coords = vtkSmartPointer<vtkPoints>::New();
  node_coords->Allocate(num_coords, 1000);
  node_coords->SetNumberOfPoints(num_coords);

  for (int i = 0; i < num_coords; i++) {
    node_coords->SetPoint(i, points(0, i), points(1, i), points(2, i));
  }

  vtk_data->SetPoints(node_coords);
}

void VtkData::set_connectivity(const int nsd, const Array<int> &conn) {
  int num_elems = conn.ncols();
  int np_elem = conn.nrows();

  const vtkIdType num_points = vtk_data->GetNumberOfPoints();

  auto elem_nodes = vtkSmartPointer<vtkIdList>::New();
  elem_nodes->Allocate(np_elem);
  elem_nodes->Initialize();
  elem_nodes->SetNumberOfIds(np_elem);

  for (int i = 0; i < num_elems; i++) {
    for (int j = 0; j < np_elem; j++) {
      const int node_id = conn(j, i);

      // The check is written out to keep the error message from being composed
      // for every point of every element.
      if (node_id < 0 || node_id >= num_points) {
        svmp::raise<svmp::FE::InvalidArgumentException>(
            "Element " + std::to_string(i) + " refers to point " +
            std::to_string(node_id) + ", which is not among the " +
            std::to_string(num_points) + " points of the mesh.");
      }

      elem_nodes->SetId(j, node_id);
    }

    insert_cell(cell_type(nsd, np_elem), elem_nodes);
  }
}

void VtkData::set_time_value(const double time) {
  auto time_array = vtkSmartPointer<vtkDoubleArray>::New();

  time_array->SetName("TimeValue");
  time_array->SetNumberOfComponents(1);
  time_array->SetNumberOfTuples(1);
  time_array->SetValue(0, time);

  vtk_data->GetFieldData()->AddArray(time_array);
}

bool VtkData::has_cell_data(const std::string &data_name) const {
  const int num_arrays = vtk_data->GetCellData()->GetNumberOfArrays();

  for (int i = 0; i < num_arrays; i++) {
    const char *array_name = vtk_data->GetCellData()->GetArrayName(i);
    if (array_name != nullptr && !strcmp(array_name, data_name.c_str())) {
      return true;
    }
  }

  return false;
}

bool VtkData::has_point_data(const std::string &data_name) const {
  const int num_arrays = vtk_data->GetPointData()->GetNumberOfArrays();

  for (int i = 0; i < num_arrays; i++) {
    const char *array_name = vtk_data->GetPointData()->GetArrayName(i);
    if (array_name != nullptr && !strcmp(array_name, data_name.c_str())) {
      return true;
    }
  }

  return false;
}

void VtkData::copy_points(Array<double> &points) const {
  auto vtk_points = vtk_data->GetPoints();
  auto num_points = vtk_points->GetNumberOfPoints();

  svmp::check<svmp::FE::InvalidArgumentException>(
      points.nrows() >= 3 && points.ncols() >= num_points,
      "The " + std::to_string(num_points) + " points of the VTK file '" +
          file_name_ + "' do not fit in an array of " +
          std::to_string(points.nrows()) + " rows and " +
          std::to_string(points.ncols()) + " columns.");

  double point[3];
  for (int i = 0; i < num_points; i++) {
    vtk_points->GetPoint(i, point);
    points(0, i) = point[0];
    points(1, i) = point[1];
    points(2, i) = point[2];
  }
}

void VtkData::copy_point_data(const std::string &data_name,
                              Array<double> &mesh_data) const {
  const auto vtk_array = vtk_data->GetPointData()->GetArray(data_name.c_str());
  svmp::check<svmp::FE::InvalidArgumentException>(
      vtk_array != nullptr,
      "There is no double-valued point data array named '" + data_name +
          "' in the VTK file '" + file_name_ + "'.");

  const int num_data = vtk_array->GetNumberOfTuples();
  const int num_comp = vtk_array->GetNumberOfComponents();

  svmp::check<svmp::FE::InvalidArgumentException>(
      num_data <= mesh_data.ncols() && num_comp <= mesh_data.nrows(),
      "The point data array named '" + data_name + "' of the VTK file '" +
          file_name_ + "' has " + std::to_string(num_comp) +
          " components and " + std::to_string(num_data) +
          " tuples, which do not fit in an array of " +
          std::to_string(mesh_data.nrows()) + " rows and " +
          std::to_string(mesh_data.ncols()) + " columns.");

  // Set the data.
  for (int i = 0; i < num_data; i++) {
    const auto tuple = vtk_array->GetTuple(i);
    for (int j = 0; j < num_comp; j++) {
      mesh_data(j, i) = tuple[j];
    }
  }
}

void VtkData::copy_point_data(const std::string &data_name,
                              Vector<double> &mesh_data) const {
  const auto vtk_array = vtkDoubleArray::SafeDownCast(
      vtk_data->GetPointData()->GetArray(data_name.c_str()));
  svmp::check<svmp::FE::InvalidArgumentException>(
      vtk_array != nullptr,
      "There is no double-valued point data array named '" + data_name +
          "' in the VTK file '" + file_name_ + "'.");

  int num_data = vtk_array->GetNumberOfTuples();

  svmp::check<svmp::FE::InvalidArgumentException>(
      num_data <= mesh_data.size(),
      "The point data array named '" + data_name + "' of the VTK file '" +
          file_name_ + "' has " + std::to_string(num_data) +
          " values, which do not fit in a vector of size " +
          std::to_string(mesh_data.size()) + ".");

  // Set the data.
  for (int i = 0; i < num_data; i++) {
    mesh_data(i) = vtk_array->GetValue(i);
  }
}

void VtkData::copy_point_data(const std::string &data_name,
                              Vector<int> &mesh_data) const {
  const auto vtk_array = vtkIntArray::SafeDownCast(
      vtk_data->GetPointData()->GetArray(data_name.c_str()));
  svmp::check<svmp::FE::InvalidArgumentException>(
      vtk_array != nullptr, "There is no int-valued point data array named '" +
                                data_name + "' in the VTK file '" + file_name_ +
                                "'.");

  int num_data = vtk_array->GetNumberOfTuples();

  svmp::check<svmp::FE::InvalidArgumentException>(
      num_data <= mesh_data.size(),
      "The point data array named '" + data_name + "' of the VTK file '" +
          file_name_ + "' has " + std::to_string(num_data) +
          " values, which do not fit in a vector of size " +
          std::to_string(mesh_data.size()) + ".");

  // Set the data.
  for (int i = 0; i < num_data; i++) {
    mesh_data(i) = vtk_array->GetValue(i);
  }
}

void VtkData::copy_cell_data(const std::string &data_name,
                             Array<double> &mesh_data) const {
  const auto vtk_array = vtk_data->GetCellData()->GetArray(data_name.c_str());
  svmp::check<svmp::FE::InvalidArgumentException>(
      vtk_array != nullptr,
      "There is no double-valued element data array named '" + data_name +
          "' in the VTK file '" + file_name_ + "'.");

  const int num_data = vtk_array->GetNumberOfTuples();
  const int num_comp = vtk_array->GetNumberOfComponents();

  svmp::check<svmp::FE::InvalidArgumentException>(
      num_data <= mesh_data.ncols() && num_comp <= mesh_data.nrows(),
      "The element data array named '" + data_name + "' of the VTK file '" +
          file_name_ + "' has " + std::to_string(num_comp) +
          " components and " + std::to_string(num_data) +
          " tuples, which do not fit in an array of " +
          std::to_string(mesh_data.nrows()) + " rows and " +
          std::to_string(mesh_data.ncols()) + " columns.");

  // Set the data.
  for (int i = 0; i < num_data; i++) {
    const auto tuple = vtk_array->GetTuple(i);
    for (int j = 0; j < num_comp; j++) {
      mesh_data(j, i) = tuple[j];
    }
  }
}

void VtkData::copy_cell_data(const std::string &data_name,
                             Vector<double> &mesh_data) const {
  const auto vtk_array = vtkDoubleArray::SafeDownCast(
      vtk_data->GetCellData()->GetArray(data_name.c_str()));
  svmp::check<svmp::FE::InvalidArgumentException>(
      vtk_array != nullptr,
      "There is no double-valued element data array named '" + data_name +
          "' in the VTK file '" + file_name_ + "'.");

  const int num_data = vtk_array->GetNumberOfTuples();

  svmp::check<svmp::FE::InvalidArgumentException>(
      num_data <= mesh_data.size(),
      "The element data array named '" + data_name + "' of the VTK file '" +
          file_name_ + "' has " + std::to_string(num_data) +
          " values, which do not fit in a vector of size " +
          std::to_string(mesh_data.size()) + ".");

  // Set the data.
  for (int i = 0; i < num_data; i++) {
    mesh_data(i) = vtk_array->GetValue(i);
  }
}

void VtkData::copy_cell_data(const std::string &data_name,
                             Vector<int> &mesh_data) const {
  const auto vtk_array = vtkIntArray::SafeDownCast(
      vtk_data->GetCellData()->GetArray(data_name.c_str()));
  svmp::check<svmp::FE::InvalidArgumentException>(
      vtk_array != nullptr,
      "There is no int-valued element data array named '" + data_name +
          "' in the VTK file '" + file_name_ + "'.");

  const int num_data = vtk_array->GetNumberOfTuples();

  svmp::check<svmp::FE::InvalidArgumentException>(
      num_data <= mesh_data.size(),
      "The element data array named '" + data_name + "' of the VTK file '" +
          file_name_ + "' has " + std::to_string(num_data) +
          " values, which do not fit in a vector of size " +
          std::to_string(mesh_data.size()) + ".");

  // Set the data.
  for (int i = 0; i < num_data; i++) {
    mesh_data(i) = vtk_array->GetValue(i);
  }
}

Array<double> VtkData::get_point_data(const std::string &data_name) const {
  auto vtk_array = vtkDoubleArray::SafeDownCast(
      vtk_data->GetPointData()->GetArray(data_name.c_str()));
  svmp::check<svmp::FE::InvalidArgumentException>(
      vtk_array != nullptr,
      "There is no double-valued point data array named '" + data_name +
          "' in the VTK file '" + file_name_ + "'.");

  int num_data = vtk_array->GetNumberOfTuples();
  int num_comp = vtk_array->GetNumberOfComponents();

  // Set the data.
  Array<double> data(num_data, num_comp);
  for (int i = 0; i < num_data; i++) {
    auto tuple = vtk_array->GetTuple(i);
    for (int j = 0; j < num_comp; j++) {
      data(i, j) = tuple[j];
    }
  }

  return data;
}

std::vector<std::string> VtkData::get_point_data_names() const {
  std::vector<std::string> data_names;

  const int num_arrays = vtk_data->GetPointData()->GetNumberOfArrays();
  for (int i = 0; i < num_arrays; i++) {
    data_names.push_back(vtk_data->GetPointData()->GetArrayName(i));
  }

  return data_names;
}

std::pair<int, int>
VtkData::get_cell_data_dimensions(const std::string &data_name) const {
  auto vtk_array = vtk_data->GetCellData()->GetArray(data_name.c_str());
  if (vtk_array == nullptr) {
    return std::make_pair(0, 0);
  }

  return std::make_pair(vtk_array->GetNumberOfComponents(),
                        vtk_array->GetNumberOfTuples());
}

VtkData *VtkData::create_reader(const std::string &file_name) {
  auto file_ext = file_name.substr(file_name.find_last_of(".") + 1);
  if (file_ext == "vtp") {
    return new VtkVtpData(file_name);
  } else if (file_ext == "vtu") {
    return new VtkVtuData(file_name);
  }

  svmp::raise<svmp::FE::InvalidArgumentException>(
      "Error in VtkData::create_reader: the file '" + file_name +
      "' has the extension '" + file_ext + "', which is not 'vtp' or 'vtu'.");
}

VtkData *VtkData::create_writer(const std::string &file_name) {
  auto file_ext = file_name.substr(file_name.find_last_of(".") + 1);
  bool reader = false;
  if (file_ext == "vtp") {
    return new VtkVtpData(file_name, reader);
  } else if (file_ext == "vtu") {
    return new VtkVtuData(file_name, reader);
  }

  svmp::raise<svmp::FE::InvalidArgumentException>(
      "Error in VtkData::create_writer: the file '" + file_name +
      "' has the extension '" + file_ext + "', which is not 'vtp' or 'vtu'.");
}

VtkVtpData::VtkVtpData() { create_grid(); }

VtkVtpData::VtkVtpData(const std::string &file_name, bool reader) {
  file_name_ = file_name;

  if (reader) {
    read_file(file_name);
  } else {
    create_grid();
  }
}

void VtkVtpData::create_grid() {
  vtk_polydata = vtkSmartPointer<vtkPolyData>::New();
  vtk_data = vtk_polydata;

  elem_type_ = -1;
  num_elems_ = 0;
  num_points_per_elem_ = 0;
  num_points_ = 0;
}

void VtkVtpData::write() const {
  auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetInputDataObject(vtk_polydata);
  writer->SetFileName(file_name_.c_str());

  const int status = writer->Write();
  svmp::check<svmp::CoreException>(
      status != 0,
      "Error writing the VTK file '" + file_name_ + "': " +
          vtkErrorCode::GetStringFromErrorCode(writer->GetErrorCode()) + ".",
      svmp::StatusCode::IOError);
}

void VtkVtpData::read_file_internal(const std::string &file_name) {
  auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(file_name.c_str());
  reader->Update();

  svmp::check<svmp::FileFormatException>(
      reader->GetErrorCode() == vtkErrorCode::NoError, file_name,
      std::string("Error reading VTK file: ") +
          vtkErrorCode::GetStringFromErrorCode(reader->GetErrorCode()) + ".");

  vtk_polydata = reader->GetOutput();
  vtk_data = vtk_polydata;
}

int VtkVtpData::cell_type(int nsd, int np_elem) const {
  if (np_elem == 2) {
    return VTK_LINE;
  }

  // A vtkPolyData holds surface elements only, whose type is determined by the
  // number of points alone.
  switch (np_elem) {
  case 3:
    return VTK_TRIANGLE;
  case 4:
    return VTK_QUAD;
  case 6:
    return VTK_QUADRATIC_TRIANGLE;
  case 8:
    return VTK_QUADRATIC_QUAD;
  case 9:
    return VTK_BIQUADRATIC_QUAD;
  }

  svmp::raise<svmp::FE::InvalidArgumentException>(
      "Error in VtkVtpData::cell_type: no surface cell type for an element "
      "with " +
      std::to_string(np_elem) + " points in " + std::to_string(nsd) +
      " dimensions.");
}

void VtkVtpData::insert_cell(int vtk_cell_type,
                             vtkSmartPointer<vtkIdList> elem_nodes) {
  // A vtkPolyData derives the type of a polygon from its number of points, so
  // vtk_cell_type is not needed here.
  vtk_polydata->GetPolys()->InsertNextCell(elem_nodes);
}

VtkVtuData::VtkVtuData() { create_grid(); }

VtkVtuData::VtkVtuData(const std::string &file_name, bool reader) {
  file_name_ = file_name;

  if (reader) {
    read_file(file_name);
  } else {
    create_grid();
  }
}

void VtkVtuData::create_grid() {
  vtk_ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtk_data = vtk_ugrid;

  elem_type_ = -1;
  num_elems_ = 0;
  num_points_per_elem_ = 0;
  num_points_ = 0;
}

void VtkVtuData::write() const {
  auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetInputDataObject(vtk_ugrid);
  writer->SetFileName(file_name_.c_str());

  const int status = writer->Write();
  svmp::check<svmp::CoreException>(
      status != 0,
      "Error writing the VTK file '" + file_name_ + "': " +
          vtkErrorCode::GetStringFromErrorCode(writer->GetErrorCode()) + ".",
      svmp::StatusCode::IOError);
}

void VtkVtuData::read_file_internal(const std::string &file_name) {
  auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(file_name.c_str());
  reader->Update();

  svmp::check<svmp::FileFormatException>(
      reader->GetErrorCode() == vtkErrorCode::NoError, file_name,
      std::string("Error reading VTK file: ") +
          vtkErrorCode::GetStringFromErrorCode(reader->GetErrorCode()) + ".");

  vtk_ugrid = reader->GetOutput();
  vtk_data = vtk_ugrid;
}

int VtkVtuData::cell_type(int nsd, int np_elem) const {
  if (np_elem == 2) {
    return VTK_LINE;
  }

  if (nsd == 2) {
    switch (np_elem) {
    case 3:
      return VTK_TRIANGLE;
    case 4:
      return VTK_QUAD;
    case 6:
      return VTK_QUADRATIC_TRIANGLE;
    case 8:
      return VTK_QUADRATIC_QUAD;
    case 9:
      return VTK_BIQUADRATIC_QUAD;
    }
  } else if (nsd == 3) {
    switch (np_elem) {
    case 3:
      return VTK_TRIANGLE;
    case 4:
      return VTK_TETRA;
    case 6:
      return VTK_WEDGE;
    case 8:
      return VTK_HEXAHEDRON;
    case 10:
      return VTK_QUADRATIC_TETRA;
    case 20:
      return VTK_QUADRATIC_HEXAHEDRON;
    case 27:
      return VTK_TRIQUADRATIC_HEXAHEDRON;
    }
  }

  svmp::raise<svmp::FE::InvalidArgumentException>(
      "Error in VtkVtuData::cell_type: no cell type for an element with " +
      std::to_string(np_elem) + " points in " + std::to_string(nsd) +
      " dimensions.");
}

void VtkVtuData::insert_cell(int vtk_cell_type,
                             vtkSmartPointer<vtkIdList> elem_nodes) {
  vtk_ugrid->InsertNextCell(vtk_cell_type, elem_nodes);
}
