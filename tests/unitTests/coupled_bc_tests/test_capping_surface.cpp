// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file test_capping_surface.cpp
 * @brief Tests of CappingSurface: orientation of the cap elements and
 * derivative of the flow rate through the cap.
 *
 * @par Setup
 * A cap of two triangles that are not coplanar, written to a VTP file and
 * loaded from it. With the node order (0, 1, 2) and (1, 3, 2), the normal of
 * both triangles has a positive z component. The cases combine:
 * - whether the file lists the second triangle in the reverse order, so that
 *   the two triangles of the file are oriented differently;
 * - whether the normals that the file gives for the triangles point along +z
 *   or -z.
 *
 * @test Cases/CappingSurfaceTest.ElementsFollowFileNormals checks that the
 * node order of each loaded element gives a normal along the one of the file.
 *
 * @test Cases/CappingSurfaceTest.FlowRateOfVelocityAlongFileNormalsIsPositive
 * checks the flow rate of a uniform velocity along the normals of the file.
 *
 * @test Cases/CappingSurfaceTest.VelocityPartMatchesFiniteDifferences compares
 * the flow rate derivative, without the displacement term, with central
 * finite differences of the flow rate with respect to the node velocities.
 *
 * @test Cases/CappingSurfaceTest.DisplacementPartMatchesFiniteDifferences
 * compares the displacement term of the flow rate derivative with central
 * finite differences of the flow rate with respect to the node displacements.
 *
 * @test Cases/CappingSurfaceTest.DisplacementPartScalesWithItsFactor checks
 * that the displacement term is proportional to displacement_per_velocity.
 *
 * @par Running the tests
 * @code{.sh}
 * ./run_all_unit_tests --gtest_filter='*CappingSurface*'
 * @endcode
 */

#include "ComMod.h"
#include "CoupledBoundaryCondition.h"
#include "consts.h"

#include <gtest/gtest.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>

#include <array>
#include <filesystem>
#include <string>
#include <tuple>

namespace {

// Whether the file lists the second triangle in the reverse order, and the z
// component of the normals of the file.
using CapParameters = std::tuple<bool, double>;

constexpr int num_nodes = 4;
constexpr int num_elements = 2;

// Positions of the cap nodes in the reference configuration.
constexpr double positions[num_nodes][3] = {{0.0, 0.0, 0.0},
                                            {1.0, 0.0, 0.1},
                                            {0.0, 1.0, -0.1},
                                            {1.1, 0.9, 0.3}};

// The GlobalNodeID of the first cap node. The following nodes are numbered
// consecutively.
constexpr int first_global_node_id = 11;

/// Write the cap to a VTP file.
void write_cap_file(const std::string& file_name, const bool reverse_second_element, const double normal_z)
{
  auto points = vtkSmartPointer<vtkPoints>::New();
  auto global_node_ids = vtkSmartPointer<vtkIntArray>::New();
  global_node_ids->SetName("GlobalNodeID");

  for (int a = 0; a < num_nodes; a++) {
    points->InsertNextPoint(positions[a][0], positions[a][1], positions[a][2]);
    global_node_ids->InsertNextValue(first_global_node_id + a);
  }

  const vtkIdType first_element[3] = {0, 1, 2};
  const vtkIdType second_element[3] = {1, 3, 2};
  const vtkIdType reversed_second_element[3] = {1, 2, 3};

  auto elements = vtkSmartPointer<vtkCellArray>::New();
  elements->InsertNextCell(3, first_element);
  elements->InsertNextCell(3, reverse_second_element ? reversed_second_element : second_element);

  auto normals = vtkSmartPointer<vtkDoubleArray>::New();
  normals->SetName("Normals");
  normals->SetNumberOfComponents(3);
  for (int e = 0; e < num_elements; e++) {
    normals->InsertNextTuple3(0.0, 0.0, normal_z);
  }

  auto polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);
  polydata->SetPolys(elements);
  polydata->GetPointData()->AddArray(global_node_ids);
  polydata->GetCellData()->AddArray(normals);

  auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(file_name.c_str());
  writer->SetInputData(polydata);
  ASSERT_EQ(writer->Write(), 1);
}

class CappingSurfaceTest : public testing::TestWithParam<CapParameters> {
protected:
  void SetUp() override
  {
    const auto [reverse_second_element, file_normal_z] = GetParam();
    normal_z = file_normal_z;

    const std::string file_name = ::testing::TempDir() + "capping_surface_" +
                                  (reverse_second_element ? "reversed_" : "") + (normal_z > 0.0 ? "up" : "down") +
                                  ".vtp";
    write_cap_file(file_name, reverse_second_element, normal_z);

    // The coupled face has to share a node with the cap.
    faceType coupled_face;
    coupled_face.iM = 0;
    coupled_face.nNo = 1;
    coupled_face.gN.resize(1);
    coupled_face.gN(0) = first_global_node_id - 1;

    ComMod com_mod;
    com_mod.nsd = 3;

    cap.load_from_vtp(file_name, coupled_face, "face");
    cap.init_cap_face_quadrature(com_mod);

    std::filesystem::remove(file_name);

    // Positions, displacements and velocities of the cap nodes.
    state.n_cap = num_nodes;
    state.x.resize(3, num_nodes);
    state.Do.resize(3, num_nodes);
    state.Dn.resize(3, num_nodes);
    state.Yo.resize(3, num_nodes);
    state.Yn.resize(3, num_nodes);

    for (int a = 0; a < num_nodes; a++) {
      for (int i = 0; i < 3; i++) {
        state.x(i, a) = positions[a][i];
        state.Dn(i, a) = 0.02 * (i + 1) * (a + 1) * ((a % 2 == 0) ? 1.0 : -1.0);
        state.Yn(i, a) = 0.3 * (i + 2) - 0.2 * (a + 1) * (i - 1);
      }
    }
  }

  /// Flow rate through the cap in the new configuration.
  double flowrate() { return cap.integrate_velocity_flux(state, true, cfg); }

  /// Central finite difference of the flow rate with respect to an entry of
  /// an array of the cap state.
  double flowrate_finite_difference(Array<double>& values, const int i, const int a)
  {
    const double step = 1.0e-3;
    const double original = values(i, a);

    values(i, a) = original + step;
    const double forward = flowrate();

    values(i, a) = original - step;
    const double backward = flowrate();

    values(i, a) = original;

    return (forward - backward) / (2.0 * step);
  }

  const consts::MechanicalConfigurationType cfg = consts::MechanicalConfigurationType::new_timestep;

  double normal_z = 0.0;
  CappingSurface cap;
  CapGlobalMeshState state;
};

TEST_P(CappingSurfaceTest, ElementsFollowFileNormals)
{
  const faceType* face = cap.face();
  ASSERT_EQ(face->nEl, num_elements);

  for (int e = 0; e < num_elements; e++) {
    // The elements refer to their nodes by global node ID, starting from zero.
    std::array<std::array<double, 3>, 3> x;
    for (int a = 0; a < 3; a++) {
      const int node = face->IEN(a, e) - (first_global_node_id - 1);
      ASSERT_GE(node, 0);
      ASSERT_LT(node, num_nodes);
      x[a] = {positions[node][0], positions[node][1], positions[node][2]};
    }

    // The z component of (x_1 - x_0) x (x_2 - x_0), the normal that follows
    // the right-hand rule on the node order.
    const double element_normal_z =
        (x[1][0] - x[0][0]) * (x[2][1] - x[0][1]) - (x[1][1] - x[0][1]) * (x[2][0] - x[0][0]);

    EXPECT_GT(element_normal_z * normal_z, 0.0) << "element " << e;
  }
}

TEST_P(CappingSurfaceTest, FlowRateOfVelocityAlongFileNormalsIsPositive)
{
  // A uniform velocity of unit magnitude along the normals of the file.
  CapGlobalMeshState uniform_state = state;
  uniform_state.Yn = 0.0;
  for (int a = 0; a < num_nodes; a++) {
    uniform_state.Yn(2, a) = normal_z;
  }

  // The flow rate is the area of the projection of the cap on the xy plane:
  // one half of the z components of the cross products of the edges of the
  // two triangles, which are both equal to 1.
  const double flowrate =
      cap.integrate_velocity_flux(uniform_state, true, consts::MechanicalConfigurationType::reference);

  EXPECT_NEAR(flowrate, 1.0, 1.0e-12);
}

TEST_P(CappingSurfaceTest, VelocityPartMatchesFiniteDifferences)
{
  const auto derivative = cap.compute_flowrate_derivative(cfg, state, 0.0);

  ASSERT_EQ(derivative.nrows(), 3);
  ASSERT_EQ(derivative.ncols(), num_nodes);

  // The flow rate is linear in the velocity, so central differences have no
  // truncation error, only round-off.
  for (int a = 0; a < num_nodes; a++) {
    for (int i = 0; i < 3; i++) {
      EXPECT_NEAR(derivative(i, a), flowrate_finite_difference(state.Yn, i, a), 1.0e-10)
          << "coordinate " << i << " of cap node " << a;
    }
  }
}

TEST_P(CappingSurfaceTest, DisplacementPartMatchesFiniteDifferences)
{
  const auto without_displacement = cap.compute_flowrate_derivative(cfg, state, 0.0);
  const auto with_displacement = cap.compute_flowrate_derivative(cfg, state, 1.0);

  // The flow rate is quadratic in the displacement, so central differences
  // have no truncation error, only round-off.
  for (int a = 0; a < num_nodes; a++) {
    for (int i = 0; i < 3; i++) {
      EXPECT_NEAR(with_displacement(i, a) - without_displacement(i, a),
                  flowrate_finite_difference(state.Dn, i, a), 1.0e-10)
          << "coordinate " << i << " of cap node " << a;
    }
  }
}

TEST_P(CappingSurfaceTest, DisplacementPartScalesWithItsFactor)
{
  const double displacement_per_velocity = 0.37;

  const auto without_displacement = cap.compute_flowrate_derivative(cfg, state, 0.0);
  const auto with_unit_factor = cap.compute_flowrate_derivative(cfg, state, 1.0);
  const auto derivative = cap.compute_flowrate_derivative(cfg, state, displacement_per_velocity);

  for (int a = 0; a < num_nodes; a++) {
    for (int i = 0; i < 3; i++) {
      const double expected = without_displacement(i, a) +
                              displacement_per_velocity * (with_unit_factor(i, a) - without_displacement(i, a));

      EXPECT_NEAR(derivative(i, a), expected, 1.0e-13) << "coordinate " << i << " of cap node " << a;
    }
  }
}

INSTANTIATE_TEST_SUITE_P(Cases, CappingSurfaceTest,
                         testing::Values(CapParameters{false, 1.0}, CapParameters{false, -1.0},
                                         CapParameters{true, 1.0}, CapParameters{true, -1.0}));

} // namespace
