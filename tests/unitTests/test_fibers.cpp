// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "ComMod.h"
#include "Core/Exception.h"
#include "Fibers.h"
#include "Parameters.h"

#include <gtest/gtest.h>

#include <vtkDoubleArray.h>
#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <cmath>
#include <filesystem>
#include <numbers>
#include <string>
#include <vector>

namespace {

mshType make_single_element_mesh()
{
  mshType mesh;
  mesh.gnEl = 1;
  mesh.gnNo = 1;
  mesh.nG = 1;
  mesh.eNoN = 1;
  mesh.name = "unit-test-mesh";
  mesh.N.resize(1, 1);
  mesh.N(0, 0) = 1.0;
  mesh.gIEN.resize(1, 1);
  mesh.gIEN(0, 0) = 0;
  return mesh;
}

mshType make_two_node_line_mesh(const int num_quadrature_points = 1)
{
  mshType mesh;
  mesh.name = "two-node-line";
  mesh.gnEl = 1;
  mesh.gnNo = 2;
  mesh.eNoN = 2;
  mesh.nG = num_quadrature_points;
  mesh.gIEN.resize(mesh.eNoN, mesh.gnEl);
  mesh.gIEN(0, 0) = 0;
  mesh.gIEN(1, 0) = 1;
  mesh.N.resize(mesh.eNoN, mesh.nG);
  for (int g = 0; g < mesh.nG; ++g) {
    mesh.N(0, g) = 0.5;
    mesh.N(1, g) = 0.5;
  }
  return mesh;
}

void write_nodal_fiber_file(
    const std::string& file_name,
    const std::vector<std::vector<double>>& directions)
{
  auto points = vtkSmartPointer<vtkPoints>::New();
  auto fibers = vtkSmartPointer<vtkDoubleArray>::New();
  fibers->SetName("FIB_DIR");
  fibers->SetNumberOfComponents(3);
  for (int node = 0; node < static_cast<int>(directions.size()); ++node) {
    points->InsertNextPoint(static_cast<double>(node), 0.0, 0.0);
    const auto& direction = directions[node];
    fibers->InsertNextTuple3(
      direction[0], direction[1], direction.size() > 2 ? direction[2] : 0.0);
  }

  auto line = vtkSmartPointer<vtkLine>::New();
  line->GetPointIds()->SetId(0, 0);
  line->GetPointIds()->SetId(1, 1);

  auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  grid->SetPoints(points);
  grid->InsertNextCell(line->GetCellType(), line->GetPointIds());
  grid->GetPointData()->AddArray(fibers);

  auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName(file_name.c_str());
  writer->SetInputData(grid);
  ASSERT_EQ(writer->Write(), 1);
}

} // namespace

TEST(Fibers, LoadsRawThenCompletesConstantDirections)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{1.0, 0.0, 0.0}, false);
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{0.0, 1.0, 0.0}, false);

  auto mesh = make_single_element_mesh();
  Fibers fibers;

  EXPECT_FALSE(fibers.has_fibers());
  EXPECT_EQ(fibers.get_num_directions(), 0);
  {
    const auto zeros = Fibers::zero_fibers(3);
    ASSERT_EQ(zeros.nrows(), 3);
    ASSERT_EQ(zeros.ncols(), 1);
    EXPECT_DOUBLE_EQ(zeros(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(zeros(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(zeros(2, 0), 0.0);

  }
  fibers.read(parameters, mesh, 3);

  ASSERT_TRUE(fibers.has_fibers());
  EXPECT_EQ(fibers.get_num_directions(), 2);

  auto quadrature_fibers = fibers.get_quadrature_fibers(3, 0, 0);
  ASSERT_EQ(quadrature_fibers.nrows(), 3);
  ASSERT_EQ(quadrature_fibers.ncols(), 2);
  EXPECT_DOUBLE_EQ(quadrature_fibers(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(2, 0), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(0, 1), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(2, 1), 0.0);

  const auto* raw_data = quadrature_fibers.data();
  const auto basis = fibers.get_quadrature_fibers(3, 0, 0, true);
  const auto raw_after = fibers.get_quadrature_fibers(3, 0, 0);
  ASSERT_EQ(raw_after.ncols(), 2);
  EXPECT_EQ(raw_after.data(), raw_data);
  EXPECT_DOUBLE_EQ(raw_after(0, 0), 1.0);
  ASSERT_EQ(basis.ncols(), 3);
  EXPECT_DOUBLE_EQ(basis(2, 2), 1.0);
  const auto second_basis = fibers.get_quadrature_fibers(3, 0, 0, true);
  EXPECT_NE(second_basis.data(), basis.data());
  EXPECT_DOUBLE_EQ(second_basis(2, 2), 1.0);

  fibers.clear();
  EXPECT_FALSE(fibers.has_fibers());
  EXPECT_EQ(fibers.get_num_directions(), 0);
}

TEST(Fibers, AcceptsRawCountRejectedByOrthonormalization)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{1.0, 0.0, 0.0}, false);

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  EXPECT_NO_THROW(fibers.read(parameters, mesh, 3));
  EXPECT_EQ(fibers.get_num_directions(), 1);
  EXPECT_THROW(fibers.get_quadrature_fibers(3, 0, 0, true),
    svmp::ParseException);
}

TEST(Fibers, PreservesConstantDirectionsWhenCompletingBasis)
{
  constexpr double fiber_x = 1.0 + 5.0e-7;
  constexpr double sheet_y = 1.0 - 5.0e-7;
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{fiber_x, 0.0, 0.0}, false);
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{0.0, sheet_y, 0.0}, false);

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  fibers.read(parameters, mesh, 3);

  const auto basis = fibers.get_quadrature_fibers(3, 0, 0, true);
  EXPECT_DOUBLE_EQ(basis(0, 0), fiber_x);
  EXPECT_DOUBLE_EQ(basis(1, 1), sheet_y);
  EXPECT_DOUBLE_EQ(basis(2, 2), fiber_x * sheet_y);
}

TEST(Fibers, CompletesTwoDimensionalConstantBasisByRotation)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{0.6, 0.8}, false);

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  fibers.read(parameters, mesh, 2);

  const auto basis = fibers.get_quadrature_fibers(2, 0, 0, true);
  EXPECT_DOUBLE_EQ(basis(0, 0), 0.6);
  EXPECT_DOUBLE_EQ(basis(1, 0), 0.8);
  EXPECT_DOUBLE_EQ(basis(0, 1), -0.8);
  EXPECT_DOUBLE_EQ(basis(1, 1), 0.6);
}

TEST(Fibers, AcceptsNoDirections)
{
  MeshParameters parameters;
  auto mesh = make_single_element_mesh();
  Fibers fibers;

  EXPECT_NO_THROW(fibers.read(parameters, mesh, 3));
  EXPECT_FALSE(fibers.has_fibers());
  EXPECT_EQ(fibers.get_num_directions(), 0);
}

TEST(Fibers, AcceptsArbitraryRawDirectionCount)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{1.0, 0.0, 0.0}, false);
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{0.0, 1.0, 0.0}, false);
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{0.0, 0.0, 1.0}, false);
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{-1.0, 0.0, 0.0}, false);

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  EXPECT_NO_THROW(fibers.read(parameters, mesh, 3));
  EXPECT_EQ(fibers.get_num_directions(), 4);
  EXPECT_EQ(fibers.get_quadrature_fibers(3, 0, 0).ncols(), 4);
  EXPECT_THROW(fibers.get_quadrature_fibers(3, 0, 0, true),
    svmp::ParseException);
}

TEST(Fibers, RejectsNonUnitConstantDirection)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{3.0, 4.0, 0.0}, false);
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{0.0, 1.0, 0.0}, false);

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  EXPECT_THROW(fibers.read(parameters, mesh, 3), svmp::ParseException);
}

TEST(Fibers, RejectsNonOrthogonalConstantDirections)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{1.0, 0.0, 0.0}, false);
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{0.6, 0.8, 0.0}, false);

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  EXPECT_NO_THROW(fibers.read(parameters, mesh, 3));
  const auto raw = fibers.get_quadrature_fibers(3, 0, 0);
  EXPECT_DOUBLE_EQ(raw(0, 1), 0.6);
  EXPECT_DOUBLE_EQ(raw(1, 1), 0.8);
  EXPECT_THROW(fibers.get_quadrature_fibers(3, 0, 0, true),
    svmp::ParseException);
}

TEST(Fibers, LoadsCellDirectionsFromVtu)
{
  MeshParameters parameters;
  parameters.fiber_direction_file_paths_cells.set(
      "Fiber_direction_from_file_cells", false,
      {std::string(UNIT_TEST_DATA_DIR) + "/fibers_cell_data.vtu",
       std::string(UNIT_TEST_DATA_DIR) + "/fibers_cell_data_sheet.vtu"});

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  fibers.read(parameters, mesh, 3);

  ASSERT_TRUE(fibers.has_fibers());
  ASSERT_EQ(fibers.get_num_directions(), 2);
  const auto quadrature_fibers = fibers.get_quadrature_fibers(3, 0, 0);
  ASSERT_EQ(quadrature_fibers.nrows(), 3);
  ASSERT_EQ(quadrature_fibers.ncols(), 2);
  EXPECT_DOUBLE_EQ(quadrature_fibers(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(1, 0), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(2, 0), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(0, 1), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(2, 1), 0.0);

  const auto basis = fibers.get_quadrature_fibers(3, 0, 0, true);
  ASSERT_EQ(basis.ncols(), 3);
  EXPECT_DOUBLE_EQ(basis(0, 0), quadrature_fibers(0, 0));
  EXPECT_DOUBLE_EQ(basis(1, 1), quadrature_fibers(1, 1));
  EXPECT_DOUBLE_EQ(basis(2, 2), 1.0);
}

TEST(Fibers, LoadsPointDirectionsFromVtu)
{
  MeshParameters parameters;
  parameters.fiber_direction_file_paths_points.set(
      "Fiber_direction_from_file_points", false,
      {std::string(UNIT_TEST_DATA_DIR) + "/fibers_point_data.vtu",
       std::string(UNIT_TEST_DATA_DIR) + "/fibers_point_data_sheet.vtu"});

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  fibers.read(parameters, mesh, 3);
  ASSERT_TRUE(fibers.has_fibers());
  ASSERT_EQ(fibers.get_num_directions(), 2);

  const auto raw = fibers.get_quadrature_fibers(3, 0, 0);
  ASSERT_EQ(raw.ncols(), 2);
  const auto quadrature_fibers =
    fibers.get_quadrature_fibers(3, 0, 0, true);
  ASSERT_EQ(quadrature_fibers.nrows(), 3);
  ASSERT_EQ(quadrature_fibers.ncols(), 3);
  EXPECT_DOUBLE_EQ(quadrature_fibers(0, 0), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(1, 0), 1.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(2, 0), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(0, 1), 1.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(1, 1), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(2, 1), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(0, 2), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(1, 2), 0.0);
  EXPECT_DOUBLE_EQ(quadrature_fibers(2, 2), -1.0);

  Vector<double> basis_functions(1);
  basis_functions[0] = 1.0;
  const auto reinterpolated =
    fibers.get_quadrature_fibers(3, basis_functions, 0, true);
  ASSERT_EQ(reinterpolated.ncols(), 3);
  EXPECT_DOUBLE_EQ(reinterpolated(2, 2), -1.0);
}

TEST(Fibers, InterpolatesNodalFibersAtWeightedQuadraturePoints)
{
  const std::string file_name =
      ::testing::TempDir() + "weighted_nodal_fibers.vtu";
  write_nodal_fiber_file(file_name, {{1.0, 0.0}, {0.0, 1.0}});

  MeshParameters parameters;
  parameters.fiber_direction_file_paths_points.set(
      "Fiber_direction_from_file_points", false, {file_name});

  auto mesh = make_two_node_line_mesh(2);
  mesh.N(0, 0) = 0.6;
  mesh.N(1, 0) = 0.4;
  mesh.N(0, 1) = 0.25;
  mesh.N(1, 1) = 0.75;

  Fibers fibers;
  fibers.read(parameters, mesh, 2);

  ASSERT_EQ(fibers.get_num_directions(), 1);
  EXPECT_NEAR(fibers.get_quadrature_fibers(2, 0, 0)(0, 0), 1.0, 1.0e-12);
  EXPECT_NEAR(fibers.get_quadrature_fibers(2, 0, 0)(1, 0), 0.0, 1.0e-12);
  EXPECT_NEAR(fibers.get_quadrature_fibers(2, 1, 0)(0, 0), 0.0, 1.0e-12);
  EXPECT_NEAR(fibers.get_quadrature_fibers(2, 1, 0)(1, 0), 1.0, 1.0e-12);

  std::filesystem::remove(file_name);
}

TEST(Fibers, DyadicInterpolationDoesNotCancelOppositeDirections)
{
  const std::string file_name =
      ::testing::TempDir() + "opposite_nodal_fibers.vtu";
  write_nodal_fiber_file(file_name, {{1.0, 0.0}, {-1.0, 0.0}});

  MeshParameters parameters;
  parameters.fiber_direction_file_paths_points.set(
      "Fiber_direction_from_file_points", false, {file_name});

  auto mesh = make_two_node_line_mesh();
  Fibers fibers;
  fibers.read(parameters, mesh, 2);

  const auto direction = fibers.get_quadrature_fibers(2, 0, 0);
  EXPECT_NEAR(std::abs(direction(0, 0)), 1.0, 1.0e-12);
  EXPECT_NEAR(direction(1, 0), 0.0, 1.0e-12);

  std::filesystem::remove(file_name);
}

TEST(Fibers, DyadicInterpolationReturnsPrincipalDirection)
{
  const double sqrt_one_half = std::sqrt(0.5);
  const std::string file_name =
      ::testing::TempDir() + "principal_nodal_fibers.vtu";
  write_nodal_fiber_file(file_name, {
    {1.0, 0.0},
    {sqrt_one_half, sqrt_one_half}
  });

  MeshParameters parameters;
  parameters.fiber_direction_file_paths_points.set(
      "Fiber_direction_from_file_points", false, {file_name});

  auto mesh = make_two_node_line_mesh();
  Fibers fibers;
  fibers.read(parameters, mesh, 2);

  const auto direction = fibers.get_quadrature_fibers(2, 0, 0);
  EXPECT_NEAR(direction(0, 0), std::cos(std::numbers::pi / 8.0), 1.0e-12);
  EXPECT_NEAR(direction(1, 0), std::sin(std::numbers::pi / 8.0), 1.0e-12);

  std::filesystem::remove(file_name);
}

TEST(Fibers, RejectsCellTagWithPointData)
{
  MeshParameters parameters;
  parameters.fiber_direction_file_paths_cells.set(
      "Fiber_direction_from_file_cells", false,
      {std::string(UNIT_TEST_DATA_DIR) + "/fibers_point_data.vtu",
       std::string(UNIT_TEST_DATA_DIR) + "/fibers_point_data_sheet.vtu"});

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  EXPECT_THROW(fibers.read(parameters, mesh, 3), svmp::FileFormatException);
}

TEST(Fibers, RejectsPointTagWithCellData)
{
  MeshParameters parameters;
  parameters.fiber_direction_file_paths_points.set(
      "Fiber_direction_from_file_points", false,
      {std::string(UNIT_TEST_DATA_DIR) + "/fibers_cell_data.vtu",
       std::string(UNIT_TEST_DATA_DIR) + "/fibers_cell_data_sheet.vtu"});

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  EXPECT_THROW(fibers.read(parameters, mesh, 3), svmp::FileFormatException);
}

TEST(Fibers, RejectsMixedInputModes)
{
  MeshParameters parameters;
  parameters.fiber_directions.emplace_back("Fiber_direction_constant",
      std::vector<double>{1.0, 0.0, 0.0}, false);
  parameters.fiber_direction_file_paths_cells.set(
      "Fiber_direction_from_file_cells", false,
      {std::string(UNIT_TEST_DATA_DIR) + "/fibers_cell_data.vtu"});

  auto mesh = make_single_element_mesh();
  Fibers fibers;
  EXPECT_THROW(fibers.read(parameters, mesh, 3), svmp::ParseException);
}
