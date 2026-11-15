// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "ActiveStress.h"
#include "CepMod.h"
#include "ComMod.h"
#include "Core/Exception.h"
#include "Parameters.h"
#include "consts.h"
#include "mat_fun.h"
#include "prestrain.h"
#include "solid/sv_struct.h"

#include <gtest/gtest.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkLine.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <random>
#include <string>
#include <vector>

namespace {

using Matrix3 = mat_fun::Matrix<3>;
using Vector3 = Eigen::Vector3d;

mshType make_single_element_mesh(const int num_quadrature_points = 1)
{
  mshType mesh;
  mesh.name = "unit-test-mesh";
  mesh.gnEl = 1;
  mesh.nEl = 1;
  mesh.nG = num_quadrature_points;
  return mesh;
}

/// One cell array Prestrain_F_g<g> per entry of @p gauss_point_values, giving
/// every cell the same F0, component i + j*nsd for row i and column j.
void write_prestrain_file(
    const std::string& file_name, const int num_cells,
    const std::vector<std::vector<double>>& gauss_point_values)
{
  auto points = vtkSmartPointer<vtkPoints>::New();
  auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  for (int cell = 0; cell < num_cells; ++cell) {
    points->InsertNextPoint(static_cast<double>(cell), 0.0, 0.0);
    points->InsertNextPoint(static_cast<double>(cell) + 0.5, 0.0, 0.0);
  }
  grid->SetPoints(points);
  for (int cell = 0; cell < num_cells; ++cell) {
    auto line = vtkSmartPointer<vtkLine>::New();
    line->GetPointIds()->SetId(0, 2 * cell);
    line->GetPointIds()->SetId(1, 2 * cell + 1);
    grid->InsertNextCell(line->GetCellType(), line->GetPointIds());
  }

  for (int g = 0; g < static_cast<int>(gauss_point_values.size()); ++g) {
    const auto& values = gauss_point_values[g];
    auto array = vtkSmartPointer<vtkDoubleArray>::New();
    array->SetName(("Prestrain_F_g" + std::to_string(g)).c_str());
    array->SetNumberOfComponents(static_cast<int>(values.size()));
    for (int cell = 0; cell < num_cells; ++cell) {
      array->InsertNextTuple(values.data());
    }
    grid->GetCellData()->AddArray(array);
  }

  auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName(file_name.c_str());
  writer->SetInputData(grid);
  ASSERT_EQ(writer->Write(), 1);
}

/// A single linear simplex on the mesh coordinates @p corners (nsd x nsd+1),
/// with @p num_quadrature_points Gauss points sharing its constant gradients.
mshType make_simplex_mesh(ComMod& com_mod, const Array<double>& corners,
    const int num_quadrature_points)
{
  const int nsd = corners.nrows();
  const int eNoN = corners.ncols();
  com_mod.nsd = nsd;
  com_mod.x = corners;

  mshType mesh;
  mesh.name = "unit-test-simplex";
  mesh.gnEl = 1;
  mesh.nEl = 1;
  mesh.eNoN = eNoN;
  mesh.nG = num_quadrature_points;
  mesh.lShpF = true;
  mesh.IEN.resize(eNoN, 1);
  for (int a = 0; a < eNoN; ++a) {
    mesh.IEN(a, 0) = a;
  }
  mesh.Nx.resize(nsd, eNoN, num_quadrature_points);
  mesh.Nx = 0.0;
  for (int g = 0; g < num_quadrature_points; ++g) {
    for (int i = 0; i < nsd; ++i) {
      mesh.Nx(i, 0, g) = -1.0;
      mesh.Nx(i, i + 1, g) = 1.0;
    }
  }
  return mesh;
}

/// The nodal displacement u = A x of an affine map with gradient @p A.
Array<double> affine_displacement(const Eigen::MatrixXd& A, const Array<double>& x)
{
  Array<double> u(x.nrows(), x.ncols());
  mat_fun::eigen_view_mutable(u) = A * mat_fun::eigen_view(x);
  return u;
}

/// One tetrahedron of a struct equation with no inertia, whose residual and
/// tangent the struct_3d kernel integrates on the stress-free configuration
/// of a prestrain, with the mesh fibers pulled back as the assembly loops do.
class StructElement
{
public:
  StructElement()
  {
    using namespace consts;
    mat_fun::ten_init(3);
    com_mod.nsd = 3;
    com_mod.dof = 3;
    com_mod.tDof = 3;
    com_mod.nsymd = 6;
    com_mod.cEq = 0;
    com_mod.cDmn = 0;
    com_mod.dt = 1.0e-3;
    com_mod.eq.resize(1);
    auto& eq = com_mod.eq[0];
    eq.s = 0;
    eq.af = 0.6;
    eq.am = 0.7;
    eq.beta = 0.3;
    eq.gam = 0.6;
    eq.maxItr = 1;
    eq.dmn.resize(1);
    auto& dmn = eq.dmn[0];
    dmn.phys = EquationType::phys_struct;
    dmn.prop[PhysicalPropertyType::solid_density] = 0.0;
    dmn.prop[PhysicalPropertyType::damping] = 0.0;
    dmn.prop[PhysicalPropertyType::f_x] = 0.0;
    dmn.prop[PhysicalPropertyType::f_y] = 0.0;
    dmn.prop[PhysicalPropertyType::f_z] = 0.0;
    dmn.stM.volType = ConstitutiveModelType::stVol_ST91;
    dmn.stM.Kpen = 5.0e4;
    use_holzapfel_ogden();

    const double corners[4][3] = {
      {0.02, -0.01, 0.03}, {1.1, 0.05, -0.02}, {-0.04, 0.95, 0.06}, {0.07, 0.03, 1.05}};
    for (int a = 0; a < 4; ++a) {
      for (int i = 0; i < 3; ++i) {
        x(i, a) = corners[a][i];
      }
    }
    Matrix3 edges;
    for (int k = 0; k < 3; ++k) {
      edges.col(k) = mat_fun::eigen_view<3>(x).col(k + 1) - mat_fun::eigen_view<3>(x).col(0);
    }
    Eigen::Matrix<double, 3, 4> parent_gradients;
    parent_gradients << -1, 1, 0, 0,
                        -1, 0, 1, 0,
                        -1, 0, 0, 1;
    mat_fun::eigen_view_mutable(Nx) = edges.inverse().transpose() * parent_gradients;
    w = std::abs(edges.determinant()) / 6.0;
    N = 0.25;
    use_orthonormal_fibers();
  }

  void use_holzapfel_ogden()
  {
    auto& stM = com_mod.eq[0].dmn[0].stM;
    stM.isoType = consts::ConstitutiveModelType::stIso_HO;
    stM.a = 590.0;
    stM.b = 8.0;
    stM.aff = 1840.0;
    stM.bff = 16.0;
    stM.ass = 270.0;
    stM.bss = 11.0;
    stM.afs = 220.0;
    stM.bfs = 9.0;
    stM.khs = 100.0;
  }

  void use_guccione()
  {
    auto& stM = com_mod.eq[0].dmn[0].stM;
    stM.isoType = consts::ConstitutiveModelType::stIso_Gucci;
    stM.C10 = 880.0;
    stM.bff = 8.0;
    stM.bss = 6.0;
    stM.bfs = 12.0;
  }

  void use_neo_hookean()
  {
    auto& stM = com_mod.eq[0].dmn[0].stM;
    stM.isoType = consts::ConstitutiveModelType::stIso_nHook;
    stM.C10 = 500.0;
  }

  void use_holzapfel_gasser_ogden()
  {
    auto& stM = com_mod.eq[0].dmn[0].stM;
    stM.isoType = consts::ConstitutiveModelType::stIso_HGO;
    stM.C10 = 500.0;
    stM.aff = 1000.0;
    stM.bff = 5.0;
    stM.ass = 800.0;
    stM.bss = 4.0;
    stM.kap = 0.1;
  }

  void use_orthonormal_fibers()
  {
    const Vector3 fiber = Vector3(1.0, 0.3, 0.1).normalized();
    const Vector3 sheet = fiber.cross(Vector3(0.0, 0.0, 1.0)).normalized();
    fibers.resize(3, 3);
    mat_fun::eigen_view_mutable(fibers) << fiber, sheet, fiber.cross(sheet);
  }

  void use_oblique_fibers()
  {
    const Vector3 fiber = Vector3(1.0, 0.3, 0.1).normalized();
    const Vector3 sheet =
      (fiber.cross(Vector3(0.0, 0.0, 1.0)).normalized() + 0.4 * fiber).normalized();
    fibers.resize(3, 2);
    mat_fun::eigen_view_mutable(fibers) << fiber, sheet;
  }

  void evaluate(const Matrix3& F0, const Array<double>& u, Array<double>& lR, Array3<double>& lK)
  {
    Array<double> al(3, 4), yl(3, 4), bfl(3, 4), pS0l(6, 4);
    al = 0.0;
    yl = 0.0;
    bfl = 0.0;
    pS0l = 0.0;
    Vector<double> pSl(6);
    lR.resize(3, 4);
    lR = 0.0;
    lK.resize(9, 4, 4);
    lK = 0.0;

    Array<double> Nx0 = Nx;
    mat_fun::eigen_view_mutable(Nx0) = F0.transpose() * mat_fun::eigen_view(Nx);
    Array<double> pulled_back_fibers = fibers;
    mat_fun::eigen_view_mutable(pulled_back_fibers) =
      (F0.inverse() * mat_fun::eigen_view(fibers)).colwise().normalized();
    ActiveStress::Evaluator active_stress_evaluator;
    struct_ns::struct_3d(com_mod, cep_mod, 4, w / F0.determinant(), N, Nx0, al, yl, u, bfl, F0,
        mat_fun::eigen_view(pulled_back_fibers), pS0l, pSl, active_stress_evaluator, lR, lK, true);
  }

  ComMod com_mod;
  CepMod cep_mod;
  Array<double> x = Array<double>(3, 4);
  Array<double> Nx = Array<double>(3, 4);
  Array<double> fibers;
  Vector<double> N = Vector<double>(4);
  double w = 0.0;
};

Matrix3 random_prestrain(const int seed)
{
  std::mt19937 generator(seed);
  std::normal_distribution<double> normal(0.0, 1.0);
  Matrix3 F0 = Matrix3::Identity();
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      F0(i, j) += 0.15 * normal(generator);
    }
  }
  return F0;
}

double max_abs(const Array<double>& A)
{
  return mat_fun::eigen_view(A).cwiseAbs().maxCoeff();
}

} // namespace

TEST(Prestrain, AcceptsNoFile)
{
  MeshParameters parameters;
  auto mesh = make_single_element_mesh();
  Prestrain prestrain;

  EXPECT_FALSE(prestrain.has_prestrain());
  EXPECT_NO_THROW(prestrain.read(parameters, mesh, 3));
  EXPECT_FALSE(prestrain.has_prestrain());
}

TEST(Prestrain, SetsIdentityAtEveryGaussPoint)
{
  Prestrain prestrain;
  prestrain.set_identity(2, 3, 2);

  ASSERT_TRUE(prestrain.has_prestrain());
  for (int e = 0; e < 2; ++e) {
    for (int g = 0; g < 3; ++g) {
      EXPECT_TRUE(prestrain.get(2, g, e).isIdentity(0.0));
    }
    EXPECT_TRUE(prestrain.get_mean(2, e).isIdentity(0.0));
  }
}

TEST(Prestrain, LoadsGaussPointPrestrainFromVtu)
{
  const std::string file_name = ::testing::TempDir() + "gauss_point_prestrain.vtu";
  write_prestrain_file(file_name, 1, {
    {1.1, 0.2, 0.3, 0.4, 0.9, 0.5, 0.6, 0.7, 1.2},
    {1.3, 0.0, 0.1, 0.2, 1.1, 0.0, 0.0, 0.3, 0.8}
  });

  MeshParameters parameters;
  parameters.prestrain_file_path.set(file_name);

  auto mesh = make_single_element_mesh(2);
  Prestrain prestrain;
  prestrain.read(parameters, mesh, 3);

  ASSERT_TRUE(prestrain.has_prestrain());
  const auto F0_g0 = prestrain.get(3, 0, 0);
  ASSERT_EQ(F0_g0.rows(), 3);
  ASSERT_EQ(F0_g0.cols(), 3);
  EXPECT_DOUBLE_EQ(F0_g0(0, 0), 1.1);
  EXPECT_DOUBLE_EQ(F0_g0(1, 0), 0.2);
  EXPECT_DOUBLE_EQ(F0_g0(2, 0), 0.3);
  EXPECT_DOUBLE_EQ(F0_g0(0, 1), 0.4);
  EXPECT_DOUBLE_EQ(F0_g0(1, 1), 0.9);
  EXPECT_DOUBLE_EQ(F0_g0(2, 1), 0.5);
  EXPECT_DOUBLE_EQ(F0_g0(0, 2), 0.6);
  EXPECT_DOUBLE_EQ(F0_g0(1, 2), 0.7);
  EXPECT_DOUBLE_EQ(F0_g0(2, 2), 1.2);

  const auto F0_g1 = prestrain.get(3, 1, 0);
  EXPECT_DOUBLE_EQ(F0_g1(0, 0), 1.3);
  EXPECT_DOUBLE_EQ(F0_g1(1, 2), 0.3);

  const auto mean = prestrain.get_mean(3, 0);
  EXPECT_DOUBLE_EQ(mean(0, 0), 1.2);
  EXPECT_DOUBLE_EQ(mean(1, 0), 0.1);
  EXPECT_DOUBLE_EQ(mean(1, 2), 0.5);

  std::filesystem::remove(file_name);
}

TEST(Prestrain, RejectsMissingFile)
{
  MeshParameters parameters;
  parameters.prestrain_file_path.set(::testing::TempDir() + "missing_prestrain.vtu");

  auto mesh = make_single_element_mesh();
  Prestrain prestrain;
  EXPECT_THROW(prestrain.read(parameters, mesh, 3), svmp::FileNotFoundException);
}

TEST(Prestrain, RejectsElementCountMismatch)
{
  const std::string file_name = ::testing::TempDir() + "two_cell_prestrain.vtu";
  write_prestrain_file(file_name, 2, {{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}});

  MeshParameters parameters;
  parameters.prestrain_file_path.set(file_name);

  auto mesh = make_single_element_mesh();
  Prestrain prestrain;
  EXPECT_THROW(prestrain.read(parameters, mesh, 3), svmp::FileFormatException);

  std::filesystem::remove(file_name);
}

TEST(Prestrain, RejectsMissingGaussPointArray)
{
  const std::string file_name = ::testing::TempDir() + "one_gauss_point_prestrain.vtu";
  write_prestrain_file(file_name, 1, {{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}});

  MeshParameters parameters;
  parameters.prestrain_file_path.set(file_name);

  auto mesh = make_single_element_mesh(2);
  Prestrain prestrain;
  EXPECT_THROW(prestrain.read(parameters, mesh, 3), svmp::FileFormatException);

  std::filesystem::remove(file_name);
}

TEST(Prestrain, ComposesDisplacementGradientOnTheLeft)
{
  ComMod com_mod;
  Array<double> corners(3, 4);
  mat_fun::eigen_view_mutable(corners) << 0.1, 1.2, 0.0, 0.2,
                                         0.0, 0.1, 0.9, 0.1,
                                         0.2, 0.0, 0.1, 1.1;
  const auto mesh = make_simplex_mesh(com_mod, corners, 2);

  Eigen::Matrix3d A, B;
  A << 0.10, 0.20, -0.05,
       0.00, -0.10, 0.30,
       0.15, 0.05, 0.20;
  B << -0.20, 0.00, 0.10,
       0.25, 0.10, 0.00,
       0.00, -0.15, 0.05;

  Prestrain prestrain;
  prestrain.set_identity(3, mesh.nG, mesh.nEl);
  prestrain.compose(com_mod, mesh, affine_displacement(A, corners));
  prestrain.compose(com_mod, mesh, affine_displacement(B, corners));

  const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  const Eigen::Matrix3d expected = (I + B) * (I + A);
  for (int g = 0; g < mesh.nG; ++g) {
    EXPECT_TRUE(prestrain.get(3, g, 0).isApprox(expected, 1.0e-12));
  }
  EXPECT_FALSE(expected.isApprox((I + A) * (I + B), 1.0e-3));
}

TEST(Prestrain, ComposesTwoDimensionalDisplacementGradientOnTheLeft)
{
  ComMod com_mod;
  Array<double> corners(2, 3);
  mat_fun::eigen_view_mutable(corners) << 0.1, 1.2, 0.0,
                                         0.0, 0.1, 0.9;
  const auto mesh = make_simplex_mesh(com_mod, corners, 1);

  Eigen::Matrix2d A, B;
  A << 0.10, 0.20,
       -0.05, 0.30;
  B << -0.20, 0.10,
       0.25, 0.05;

  Prestrain prestrain;
  prestrain.set_identity(2, mesh.nG, mesh.nEl);
  prestrain.compose(com_mod, mesh, affine_displacement(A, corners));
  prestrain.compose(com_mod, mesh, affine_displacement(B, corners));

  const Eigen::Matrix2d I = Eigen::Matrix2d::Identity();
  EXPECT_TRUE(prestrain.get(2, 0, 0).isApprox((I + B) * (I + A), 1.0e-12));
}

TEST(Prestrain, InternalForcesAreObjectiveUnderRigidRotation)
{
  StructElement element;
  const Matrix3 F0 = random_prestrain(11);

  const Matrix3 Q = Eigen::AngleAxisd(0.7, Vector3(0.3, -1.0, 0.5).normalized()).toRotationMatrix();
  Array<double> u(3, 4), u_rotated(3, 4);
  u = 0.0;
  mat_fun::eigen_view_mutable(u_rotated) =
    (Q - Matrix3::Identity()) * mat_fun::eigen_view(element.x);

  Array<double> lR, lR_rotated;
  Array3<double> lK;
  element.evaluate(F0, u, lR, lK);
  element.evaluate(F0, u_rotated, lR_rotated, lK);

  const double error =
    (mat_fun::eigen_view(lR_rotated) - Q * mat_fun::eigen_view(lR)).cwiseAbs().maxCoeff();
  EXPECT_LT(error / max_abs(lR), 1.0e-13);
}

TEST(Prestrain, ComposedTangentMatchesFiniteDifferences)
{
  struct Case {
    const char* name;
    void (StructElement::*material)();
    void (StructElement::*fibers)();
  };
  const Case cases[] = {
    {"Holzapfel-Ogden, orthonormal fibers",
      &StructElement::use_holzapfel_ogden, &StructElement::use_orthonormal_fibers},
    {"Guccione, orthonormal fibers",
      &StructElement::use_guccione, &StructElement::use_orthonormal_fibers},
    {"neo-Hookean, orthonormal fibers",
      &StructElement::use_neo_hookean, &StructElement::use_orthonormal_fibers},
    {"Holzapfel-Ogden, oblique fibers",
      &StructElement::use_holzapfel_ogden, &StructElement::use_oblique_fibers},
    {"Holzapfel-Gasser-Ogden, oblique fibers",
      &StructElement::use_holzapfel_gasser_ogden, &StructElement::use_oblique_fibers},
  };

  for (const auto& test_case : cases) {
    StructElement element;
    (element.*test_case.material)();
    (element.*test_case.fibers)();
    const Matrix3 F0 = random_prestrain(11);
    const auto& eq = element.com_mod.eq[0];
    const double afu = eq.af * eq.beta * element.com_mod.dt * element.com_mod.dt;
    Array<double> rest(3, 4);
    rest = 0.0;

    Array<double> lR;
    Array3<double> lK_forward, lK;
    element.com_mod.prestrainEq = false;
    element.evaluate(F0, rest, lR, lK_forward);
    element.com_mod.prestrainEq = true;
    element.evaluate(F0, rest, lR, lK);

    const auto composed_residual = [&](const Array<double>& du) {
      const Matrix3 gradu = mat_fun::eigen_view<3>(du) * mat_fun::eigen_view<3>(element.Nx).transpose();
      Array<double> R;
      Array3<double> K;
      element.evaluate((Matrix3::Identity() + gradu) * F0, rest, R, K);
      return R;
    };

    const double h = 1.0e-6;
    double error_forward = 0.0, error = 0.0, scale = 0.0;
    for (int b = 0; b < 4; ++b) {
      for (int j = 0; j < 3; ++j) {
        Array<double> du_plus(3, 4), du_minus(3, 4);
        du_plus = 0.0;
        du_minus = 0.0;
        du_plus(j, b) = h;
        du_minus(j, b) = -h;
        const Array<double> R_plus = composed_residual(du_plus);
        const Array<double> R_minus = composed_residual(du_minus);
        for (int a = 0; a < 4; ++a) {
          for (int i = 0; i < 3; ++i) {
            const double fd = afu * (R_plus(i, a) - R_minus(i, a)) / (2.0 * h);
            error_forward = std::max(error_forward, std::abs(fd - lK_forward(i * 3 + j, a, b)));
            error = std::max(error, std::abs(fd - lK(i * 3 + j, a, b)));
            scale = std::max(scale, std::abs(fd));
          }
        }
      }
    }
    EXPECT_LT(error / scale, 1.0e-7) << test_case.name;
    EXPECT_GT(error_forward / scale, 1.0e-3) << test_case.name;
  }
}
