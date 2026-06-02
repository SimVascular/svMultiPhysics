/**
 * @file test_LagrangeBasis.cpp
 * @brief Unit tests for Lagrange basis functions
 */

#include <gtest/gtest.h>
#include "FE/Basis/BasisFactory.h"
#include "FE/Basis/LagrangeBasis.h"
#include "FE/Basis/NodeOrderingConventions.h"
#include "FE/Basis/SerendipityBasis.h"
#include "fs.h"
#include "nn.h"
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <math.h>
#include <numeric>
#include <string>
#include <vector>

namespace legacy_solver_nn {
using namespace consts;
#include "nn_elem_gip.h"
#include "nn_elem_gnn.h"
#include "nn_elem_gnnxx.h"
} // namespace legacy_solver_nn

using svmp::FE::basis::LagrangeBasis;
using svmp::FE::ElementType;
using svmp::FE::Real;
using svmp::FE::basis::Gradient;
using svmp::FE::basis::Hessian;
using svmp::FE::basis::ReferenceNodeLayout;

namespace {

using Point = svmp::FE::math::Vector<Real, 3>;

struct SolverBasisAdapterCase {
    consts::ElementType type;
    consts::ElementType quadrature_type;
    int insd;
    int eNoN;
    int nG;
};

std::vector<SolverBasisAdapterCase> solver_basis_adapter_cases() {
    using consts::ElementType;
    return {
        {ElementType::LIN1, ElementType::LIN1, 1, 2, 2},
        {ElementType::LIN2, ElementType::LIN2, 1, 3, 3},
        {ElementType::TRI3, ElementType::TRI3, 2, 3, 3},
        {ElementType::TRI6, ElementType::TRI6, 2, 6, 7},
        {ElementType::QUD4, ElementType::QUD4, 2, 4, 4},
        {ElementType::QUD8, ElementType::QUD9, 2, 8, 9},
        {ElementType::QUD9, ElementType::QUD9, 2, 9, 9},
        {ElementType::TET4, ElementType::TET4, 3, 4, 4},
        {ElementType::TET10, ElementType::TET10, 3, 10, 15},
        {ElementType::HEX8, ElementType::HEX8, 3, 8, 8},
        {ElementType::HEX20, ElementType::HEX20, 3, 20, 27},
        {ElementType::HEX27, ElementType::HEX27, 3, 27, 27},
        {ElementType::WDG, ElementType::WDG, 3, 6, 6},
    };
}

std::vector<SolverBasisAdapterCase> solver_face_basis_adapter_cases() {
    using consts::ElementType;
    return {
        {ElementType::LIN1, ElementType::LIN1, 1, 2, 2},
        {ElementType::LIN2, ElementType::LIN2, 1, 3, 3},
        {ElementType::TRI3, ElementType::TRI3, 2, 3, 3},
        {ElementType::TRI6, ElementType::TRI6, 2, 6, 7},
        {ElementType::QUD4, ElementType::QUD4, 2, 4, 4},
        {ElementType::QUD8, ElementType::QUD8, 2, 8, 9},
        {ElementType::QUD9, ElementType::QUD9, 2, 9, 9},
    };
}

std::vector<SolverBasisAdapterCase> solver_hessian_adapter_cases() {
    return solver_basis_adapter_cases();
}

std::vector<SolverBasisAdapterCase> solver_legacy_hessian_parity_cases() {
    using consts::ElementType;
    return {
        {ElementType::TRI6, ElementType::TRI6, 2, 6, 7},
        {ElementType::QUD9, ElementType::QUD9, 2, 9, 9},
        {ElementType::TET10, ElementType::TET10, 3, 10, 15},
    };
}

int packed_hessian_components(int insd) {
    if (insd == 1) {
        return 1;
    }
    if (insd == 2) {
        return 3;
    }
    return 6;
}

void fill_legacy_quadrature(const SolverBasisAdapterCase& c,
                            Vector<double>& w,
                            Array<double>& xi) {
    mshType mesh;
    mesh.eType = c.quadrature_type;
    mesh.eNoN = c.eNoN;
    mesh.nG = c.nG;
    mesh.w.resize(c.nG);
    mesh.xi.resize(c.insd, c.nG);
    legacy_solver_nn::set_element_gauss_int_data.at(c.quadrature_type)(mesh);
    w = mesh.w;
    xi = mesh.xi;
}

faceType initialized_face_for_case(const SolverBasisAdapterCase& c) {
    faceType face;
    face.eType = c.type;
    face.eNoN = c.eNoN;
    face.nG = c.nG;
    face.w.resize(c.nG);
    face.xi.resize(c.insd, c.nG);
    legacy_solver_nn::set_face_gauss_int_data.at(c.quadrature_type)(face);
    face.N.resize(c.eNoN, c.nG);
    face.Nx.resize(c.insd, c.eNoN, c.nG);
    return face;
}

void expect_arrays_near(const Array<double>& actual,
                        const Array<double>& expected,
                        double tol) {
    ASSERT_EQ(actual.nrows(), expected.nrows());
    ASSERT_EQ(actual.ncols(), expected.ncols());
    for (int col = 0; col < actual.ncols(); ++col) {
        for (int row = 0; row < actual.nrows(); ++row) {
            EXPECT_NEAR(actual(row, col), expected(row, col), tol)
                << "row=" << row << ", col=" << col;
        }
    }
}

void expect_vectors_near(const Vector<double>& actual,
                         const Vector<double>& expected,
                         double tol) {
    ASSERT_EQ(actual.size(), expected.size());
    for (int i = 0; i < actual.size(); ++i) {
        EXPECT_NEAR(actual(i), expected(i), tol) << "index=" << i;
    }
}

void expect_array3_near(const Array3<double>& actual,
                        const Array3<double>& expected,
                        double tol) {
    ASSERT_EQ(actual.nrows(), expected.nrows());
    ASSERT_EQ(actual.ncols(), expected.ncols());
    ASSERT_EQ(actual.nslices(), expected.nslices());
    for (int slice = 0; slice < actual.nslices(); ++slice) {
        for (int col = 0; col < actual.ncols(); ++col) {
            for (int row = 0; row < actual.nrows(); ++row) {
                EXPECT_NEAR(actual(row, col, slice), expected(row, col, slice), tol)
                    << "row=" << row << ", col=" << col << ", slice=" << slice;
            }
        }
    }
}

void fill_array3(Array3<double>& values, double value) {
    for (int slice = 0; slice < values.nslices(); ++slice) {
        for (int col = 0; col < values.ncols(); ++col) {
            for (int row = 0; row < values.nrows(); ++row) {
                values(row, col, slice) = value;
            }
        }
    }
}

void expect_face_partition_identities(const SolverBasisAdapterCase& c,
                                      const faceType& face,
                                      int g,
                                      double tol) {
    double partition = 0.0;
    std::array<double, 3> gradient_sum{0.0, 0.0, 0.0};

    for (int a = 0; a < c.eNoN; ++a) {
        EXPECT_TRUE(std::isfinite(face.N(a, g)))
            << "element=" << static_cast<int>(c.type)
            << ", node=" << a
            << ", g=" << g;
        partition += face.N(a, g);

        for (int d = 0; d < c.insd; ++d) {
            EXPECT_TRUE(std::isfinite(face.Nx(d, a, g)))
                << "element=" << static_cast<int>(c.type)
                << ", d=" << d
                << ", node=" << a
                << ", g=" << g;
            gradient_sum[static_cast<std::size_t>(d)] += face.Nx(d, a, g);
        }
    }

    EXPECT_NEAR(partition, 1.0, tol)
        << "element=" << static_cast<int>(c.type) << ", g=" << g;
    for (int d = 0; d < c.insd; ++d) {
        EXPECT_NEAR(gradient_sum[static_cast<std::size_t>(d)], 0.0, tol)
            << "element=" << static_cast<int>(c.type) << ", d=" << d << ", g=" << g;
    }
}

bool array3_has_nonzero_component(const Array3<double>& values,
                                  int row,
                                  double tol) {
    for (int slice = 0; slice < values.nslices(); ++slice) {
        for (int col = 0; col < values.ncols(); ++col) {
            if (std::abs(values(row, col, slice)) > tol) {
                return true;
            }
        }
    }
    return false;
}

Array<double> single_point_xi(const SolverBasisAdapterCase& c,
                              const Array<double>& xi,
                              int g) {
    Array<double> point(c.insd, 1);
    for (int d = 0; d < c.insd; ++d) {
        point(d, 0) = xi(d, g);
    }
    return point;
}

std::vector<double> finite_difference_solver_second_derivative(
    const SolverBasisAdapterCase& c,
    const Array<double>& point,
    int gradient_component,
    int coordinate_component,
    double eps) {
    Array<double> xi_plus = point;
    Array<double> xi_minus = point;
    xi_plus(coordinate_component, 0) += eps;
    xi_minus(coordinate_component, 0) -= eps;

    Array<double> N_plus(c.eNoN, 1);
    Array<double> N_minus(c.eNoN, 1);
    Array3<double> Nx_plus(c.insd, c.eNoN, 1);
    Array3<double> Nx_minus(c.insd, c.eNoN, 1);

    nn::get_gnn(c.insd, c.type, c.eNoN, 0, xi_plus, N_plus, Nx_plus);
    nn::get_gnn(c.insd, c.type, c.eNoN, 0, xi_minus, N_minus, Nx_minus);

    std::vector<double> values(static_cast<std::size_t>(c.eNoN));
    for (int a = 0; a < c.eNoN; ++a) {
        values[static_cast<std::size_t>(a)] =
            (Nx_plus(gradient_component, a, 0) - Nx_minus(gradient_component, a, 0)) /
            (2.0 * eps);
    }
    return values;
}

void expect_packed_hessian_component_matches_finite_difference(
    const SolverBasisAdapterCase& c,
    const Array<double>& point,
    const Array3<double>& Nxx,
    int g,
    int packed_row,
    int first_derivative_component,
    int second_derivative_component,
    double tol) {
    const double eps = 2e-6;
    const auto numerical = finite_difference_solver_second_derivative(
        c, point, first_derivative_component, second_derivative_component, eps);
    for (int a = 0; a < c.eNoN; ++a) {
        EXPECT_NEAR(Nxx(packed_row, a, g), numerical[static_cast<std::size_t>(a)], tol)
            << "element=" << static_cast<int>(c.type)
            << ", packed_row=" << packed_row
            << ", node=" << a
            << ", g=" << g;
    }

    if (first_derivative_component != second_derivative_component) {
        const auto symmetric_numerical = finite_difference_solver_second_derivative(
            c, point, second_derivative_component, first_derivative_component, eps);
        for (int a = 0; a < c.eNoN; ++a) {
            EXPECT_NEAR(Nxx(packed_row, a, g),
                        symmetric_numerical[static_cast<std::size_t>(a)],
                        tol)
                << "element=" << static_cast<int>(c.type)
                << ", symmetry packed_row=" << packed_row
                << ", node=" << a
                << ", g=" << g;
        }
    }
}

void expect_solver_hessian_matches_gradient_finite_difference(
    const SolverBasisAdapterCase& c,
    const Array<double>& xi,
    int g,
    const Array3<double>& Nxx,
    double tol) {
    const Array<double> point = single_point_xi(c, xi, g);

    expect_packed_hessian_component_matches_finite_difference(c, point, Nxx, g, 0, 0, 0, tol);
    if (c.insd >= 2) {
        expect_packed_hessian_component_matches_finite_difference(c, point, Nxx, g, 1, 1, 1, tol);
    }
    if (c.insd == 2) {
        expect_packed_hessian_component_matches_finite_difference(c, point, Nxx, g, 2, 0, 1, tol);
    } else if (c.insd >= 3) {
        expect_packed_hessian_component_matches_finite_difference(c, point, Nxx, g, 2, 2, 2, tol);
        expect_packed_hessian_component_matches_finite_difference(c, point, Nxx, g, 3, 0, 1, tol);
        expect_packed_hessian_component_matches_finite_difference(c, point, Nxx, g, 4, 1, 2, tol);
        expect_packed_hessian_component_matches_finite_difference(c, point, Nxx, g, 5, 0, 2, tol);
    }
}

void expect_partition_hessian_identity(const SolverBasisAdapterCase& c,
                                       const Array3<double>& Nxx,
                                       int g,
                                       double tol) {
    for (int row = 0; row < Nxx.nrows(); ++row) {
        double sum = 0.0;
        for (int a = 0; a < c.eNoN; ++a) {
            sum += Nxx(row, a, g);
        }
        EXPECT_NEAR(sum, 0.0, tol)
            << "element=" << static_cast<int>(c.type)
            << ", packed_row=" << row
            << ", g=" << g;
    }
}

void expect_all_hessians_zero(const SolverBasisAdapterCase& c,
                              const Array3<double>& Nxx,
                              int g,
                              double tol) {
    for (int row = 0; row < Nxx.nrows(); ++row) {
        for (int a = 0; a < c.eNoN; ++a) {
            EXPECT_NEAR(Nxx(row, a, g), 0.0, tol)
                << "element=" << static_cast<int>(c.type)
                << ", packed_row=" << row
                << ", node=" << a
                << ", g=" << g;
        }
    }
}

mshType initialized_mesh_for_case(const SolverBasisAdapterCase& c, bool force_lShpF) {
    mshType mesh;
    mesh.nFs = 1;
    mesh.eType = c.type;
    mesh.eNoN = c.eNoN;
    mesh.nG = c.nG;
    mesh.lShpF = force_lShpF;
    mesh.w.resize(c.nG);
    mesh.xi.resize(c.insd, c.nG);
    mesh.N.resize(c.eNoN, c.nG);
    mesh.Nx.resize(c.insd, c.eNoN, c.nG);
    mesh.xib.resize(2, c.insd);
    mesh.Nb.resize(2, c.eNoN);

    nn::get_gip(c.insd, c.quadrature_type, c.nG, mesh.w, mesh.xi);
    for (int g = 0; g < c.nG; ++g) {
        nn::get_gnn(c.insd, c.type, c.eNoN, g, mesh.xi, mesh.N, mesh.Nx);
    }
    nn::get_nn_bnds(c.insd, c.type, c.eNoN, mesh.xib, mesh.Nb);
    return mesh;
}

enum class PyramidFace {
    Base,
    South,
    East,
    North,
    West
};

enum class PyramidEdge {
    BaseSouth,
    BaseEast,
    BaseNorth,
    BaseWest,
    VerticalSW,
    VerticalSE,
    VerticalNE,
    VerticalNW
};

struct LagrangeAccuracyCase {
    ElementType type;
    int order;
    std::vector<Point> points;
};

std::size_t expected_lagrange_size(ElementType type, int order) {
    switch (type) {
        case ElementType::Point1:
            return 1u;
        case ElementType::Line2:
        case ElementType::Line3:
            return static_cast<std::size_t>(order + 1);
        case ElementType::Triangle3:
        case ElementType::Triangle6:
            return static_cast<std::size_t>(order + 1) * static_cast<std::size_t>(order + 2) / 2;
        case ElementType::Quad4:
        case ElementType::Quad9:
            return static_cast<std::size_t>(order + 1) * static_cast<std::size_t>(order + 1);
        case ElementType::Tetra4:
        case ElementType::Tetra10:
            return static_cast<std::size_t>(order + 1) *
                   static_cast<std::size_t>(order + 2) *
                   static_cast<std::size_t>(order + 3) / 6;
        case ElementType::Hex8:
        case ElementType::Hex27:
            return static_cast<std::size_t>(order + 1) *
                   static_cast<std::size_t>(order + 1) *
                   static_cast<std::size_t>(order + 1);
        case ElementType::Wedge6:
        case ElementType::Wedge18:
            return static_cast<std::size_t>(order + 1) *
                   static_cast<std::size_t>(order + 1) *
                   static_cast<std::size_t>(order + 2) / 2;
        case ElementType::Pyramid5:
        case ElementType::Pyramid14:
            return static_cast<std::size_t>(order + 1) *
                   static_cast<std::size_t>(order + 2) *
                   static_cast<std::size_t>(2 * order + 3) / 6;
        default:
            return 0u;
    }
}

int expected_dimension(ElementType type) {
    switch (type) {
        case ElementType::Point1:
            return 0;
        case ElementType::Line2:
        case ElementType::Line3:
            return 1;
        case ElementType::Triangle3:
        case ElementType::Triangle6:
        case ElementType::Quad4:
        case ElementType::Quad9:
            return 2;
        default:
            return 3;
    }
}

bool points_close(const Point& a,
                  const Point& b,
                  Real tol = Real(1e-12)) {
    return std::abs(a[0] - b[0]) <= tol &&
           std::abs(a[1] - b[1]) <= tol &&
           std::abs(a[2] - b[2]) <= tol;
}

std::vector<Point> reference_node_coords(ElementType type) {
    switch (type) {
        case ElementType::Line2:
            return {
                Point{Real(-1), Real(0), Real(0)},
                Point{Real(1), Real(0), Real(0)},
            };
        case ElementType::Line3:
            return {
                Point{Real(-1), Real(0), Real(0)},
                Point{Real(1), Real(0), Real(0)},
                Point{Real(0), Real(0), Real(0)},
            };
        case ElementType::Triangle3:
            return {
                Point{Real(0), Real(0), Real(0)},
                Point{Real(1), Real(0), Real(0)},
                Point{Real(0), Real(1), Real(0)},
            };
        case ElementType::Triangle6:
            return {
                Point{Real(0), Real(0), Real(0)},
                Point{Real(1), Real(0), Real(0)},
                Point{Real(0), Real(1), Real(0)},
                Point{Real(0.5), Real(0), Real(0)},
                Point{Real(0.5), Real(0.5), Real(0)},
                Point{Real(0), Real(0.5), Real(0)},
            };
        case ElementType::Quad4:
            return {
                Point{Real(-1), Real(-1), Real(0)},
                Point{Real(1), Real(-1), Real(0)},
                Point{Real(1), Real(1), Real(0)},
                Point{Real(-1), Real(1), Real(0)},
            };
        case ElementType::Quad8:
            return {
                Point{Real(-1), Real(-1), Real(0)},
                Point{Real(1), Real(-1), Real(0)},
                Point{Real(1), Real(1), Real(0)},
                Point{Real(-1), Real(1), Real(0)},
                Point{Real(0), Real(-1), Real(0)},
                Point{Real(1), Real(0), Real(0)},
                Point{Real(0), Real(1), Real(0)},
                Point{Real(-1), Real(0), Real(0)},
            };
        case ElementType::Quad9:
            return {
                Point{Real(-1), Real(-1), Real(0)},
                Point{Real(1), Real(-1), Real(0)},
                Point{Real(1), Real(1), Real(0)},
                Point{Real(-1), Real(1), Real(0)},
                Point{Real(0), Real(-1), Real(0)},
                Point{Real(1), Real(0), Real(0)},
                Point{Real(0), Real(1), Real(0)},
                Point{Real(-1), Real(0), Real(0)},
                Point{Real(0), Real(0), Real(0)},
            };
        case ElementType::Tetra4:
            return {
                Point{Real(0), Real(0), Real(0)},
                Point{Real(1), Real(0), Real(0)},
                Point{Real(0), Real(1), Real(0)},
                Point{Real(0), Real(0), Real(1)},
            };
        case ElementType::Tetra10:
            return {
                Point{Real(0), Real(0), Real(0)},
                Point{Real(1), Real(0), Real(0)},
                Point{Real(0), Real(1), Real(0)},
                Point{Real(0), Real(0), Real(1)},
                Point{Real(0.5), Real(0), Real(0)},
                Point{Real(0.5), Real(0.5), Real(0)},
                Point{Real(0), Real(0.5), Real(0)},
                Point{Real(0), Real(0), Real(0.5)},
                Point{Real(0.5), Real(0), Real(0.5)},
                Point{Real(0), Real(0.5), Real(0.5)},
            };
        case ElementType::Hex8:
            return {
                Point{Real(-1), Real(-1), Real(-1)},
                Point{Real(1), Real(-1), Real(-1)},
                Point{Real(1), Real(1), Real(-1)},
                Point{Real(-1), Real(1), Real(-1)},
                Point{Real(-1), Real(-1), Real(1)},
                Point{Real(1), Real(-1), Real(1)},
                Point{Real(1), Real(1), Real(1)},
                Point{Real(-1), Real(1), Real(1)},
            };
        case ElementType::Hex20:
            return {
                Point{Real(-1), Real(-1), Real(-1)},
                Point{Real(1), Real(-1), Real(-1)},
                Point{Real(1), Real(1), Real(-1)},
                Point{Real(-1), Real(1), Real(-1)},
                Point{Real(-1), Real(-1), Real(1)},
                Point{Real(1), Real(-1), Real(1)},
                Point{Real(1), Real(1), Real(1)},
                Point{Real(-1), Real(1), Real(1)},
                Point{Real(0), Real(-1), Real(-1)},
                Point{Real(1), Real(0), Real(-1)},
                Point{Real(0), Real(1), Real(-1)},
                Point{Real(-1), Real(0), Real(-1)},
                Point{Real(0), Real(-1), Real(1)},
                Point{Real(1), Real(0), Real(1)},
                Point{Real(0), Real(1), Real(1)},
                Point{Real(-1), Real(0), Real(1)},
                Point{Real(-1), Real(-1), Real(0)},
                Point{Real(1), Real(-1), Real(0)},
                Point{Real(1), Real(1), Real(0)},
                Point{Real(-1), Real(1), Real(0)},
            };
        case ElementType::Hex27:
            return {
                Point{Real(-1), Real(-1), Real(-1)},
                Point{Real(1), Real(-1), Real(-1)},
                Point{Real(1), Real(1), Real(-1)},
                Point{Real(-1), Real(1), Real(-1)},
                Point{Real(-1), Real(-1), Real(1)},
                Point{Real(1), Real(-1), Real(1)},
                Point{Real(1), Real(1), Real(1)},
                Point{Real(-1), Real(1), Real(1)},
                Point{Real(0), Real(-1), Real(-1)},
                Point{Real(1), Real(0), Real(-1)},
                Point{Real(0), Real(1), Real(-1)},
                Point{Real(-1), Real(0), Real(-1)},
                Point{Real(0), Real(-1), Real(1)},
                Point{Real(1), Real(0), Real(1)},
                Point{Real(0), Real(1), Real(1)},
                Point{Real(-1), Real(0), Real(1)},
                Point{Real(-1), Real(-1), Real(0)},
                Point{Real(1), Real(-1), Real(0)},
                Point{Real(1), Real(1), Real(0)},
                Point{Real(-1), Real(1), Real(0)},
                Point{Real(0), Real(0), Real(-1)},
                Point{Real(0), Real(0), Real(1)},
                Point{Real(0), Real(-1), Real(0)},
                Point{Real(1), Real(0), Real(0)},
                Point{Real(0), Real(1), Real(0)},
                Point{Real(-1), Real(0), Real(0)},
                Point{Real(0), Real(0), Real(0)},
            };
        case ElementType::Wedge6:
            return {
                Point{Real(0), Real(0), Real(-1)},
                Point{Real(1), Real(0), Real(-1)},
                Point{Real(0), Real(1), Real(-1)},
                Point{Real(0), Real(0), Real(1)},
                Point{Real(1), Real(0), Real(1)},
                Point{Real(0), Real(1), Real(1)},
            };
        case ElementType::Wedge15:
            return {
                Point{Real(0), Real(0), Real(-1)},
                Point{Real(1), Real(0), Real(-1)},
                Point{Real(0), Real(1), Real(-1)},
                Point{Real(0), Real(0), Real(1)},
                Point{Real(1), Real(0), Real(1)},
                Point{Real(0), Real(1), Real(1)},
                Point{Real(0.5), Real(0), Real(-1)},
                Point{Real(0.5), Real(0.5), Real(-1)},
                Point{Real(0), Real(0.5), Real(-1)},
                Point{Real(0.5), Real(0), Real(1)},
                Point{Real(0.5), Real(0.5), Real(1)},
                Point{Real(0), Real(0.5), Real(1)},
                Point{Real(0), Real(0), Real(0)},
                Point{Real(1), Real(0), Real(0)},
                Point{Real(0), Real(1), Real(0)},
            };
        case ElementType::Wedge18:
            return {
                Point{Real(0), Real(0), Real(-1)},
                Point{Real(1), Real(0), Real(-1)},
                Point{Real(0), Real(1), Real(-1)},
                Point{Real(0), Real(0), Real(1)},
                Point{Real(1), Real(0), Real(1)},
                Point{Real(0), Real(1), Real(1)},
                Point{Real(0.5), Real(0), Real(-1)},
                Point{Real(0.5), Real(0.5), Real(-1)},
                Point{Real(0), Real(0.5), Real(-1)},
                Point{Real(0.5), Real(0), Real(1)},
                Point{Real(0.5), Real(0.5), Real(1)},
                Point{Real(0), Real(0.5), Real(1)},
                Point{Real(0), Real(0), Real(0)},
                Point{Real(1), Real(0), Real(0)},
                Point{Real(0), Real(1), Real(0)},
                Point{Real(0.5), Real(0), Real(0)},
                Point{Real(0.5), Real(0.5), Real(0)},
                Point{Real(0), Real(0.5), Real(0)},
            };
        case ElementType::Pyramid5:
            return {
                Point{Real(-1), Real(-1), Real(0)},
                Point{Real(1), Real(-1), Real(0)},
                Point{Real(1), Real(1), Real(0)},
                Point{Real(-1), Real(1), Real(0)},
                Point{Real(0), Real(0), Real(1)},
            };
        case ElementType::Pyramid13:
            return {
                Point{Real(-1), Real(-1), Real(0)},
                Point{Real(1), Real(-1), Real(0)},
                Point{Real(1), Real(1), Real(0)},
                Point{Real(-1), Real(1), Real(0)},
                Point{Real(0), Real(0), Real(1)},
                Point{Real(0), Real(-1), Real(0)},
                Point{Real(1), Real(0), Real(0)},
                Point{Real(0), Real(1), Real(0)},
                Point{Real(-1), Real(0), Real(0)},
                Point{Real(-0.5), Real(-0.5), Real(0.5)},
                Point{Real(0.5), Real(-0.5), Real(0.5)},
                Point{Real(0.5), Real(0.5), Real(0.5)},
                Point{Real(-0.5), Real(0.5), Real(0.5)},
            };
        case ElementType::Pyramid14:
            return {
                Point{Real(-1), Real(-1), Real(0)},
                Point{Real(1), Real(-1), Real(0)},
                Point{Real(1), Real(1), Real(0)},
                Point{Real(-1), Real(1), Real(0)},
                Point{Real(0), Real(0), Real(1)},
                Point{Real(0), Real(-1), Real(0)},
                Point{Real(1), Real(0), Real(0)},
                Point{Real(0), Real(1), Real(0)},
                Point{Real(-1), Real(0), Real(0)},
                Point{Real(-0.5), Real(-0.5), Real(0.5)},
                Point{Real(0.5), Real(-0.5), Real(0.5)},
                Point{Real(0.5), Real(0.5), Real(0.5)},
                Point{Real(-0.5), Real(0.5), Real(0.5)},
                Point{Real(0), Real(0), Real(0)},
            };
        default:
            return {};
    }
}

void expect_nodes_match_node_ordering(ElementType canonical_type,
                                      int order,
                                      ElementType node_ordering_type) {
    LagrangeBasis basis(canonical_type, order);
    const auto& nodes = basis.nodes();

    ASSERT_EQ(nodes.size(), ReferenceNodeLayout::num_nodes(node_ordering_type));
    ASSERT_EQ(nodes.size(), basis.size());

    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const auto expected = ReferenceNodeLayout::get_node_coords(node_ordering_type, i);
        EXPECT_NEAR(nodes[i][0], expected[0], 1e-14);
        EXPECT_NEAR(nodes[i][1], expected[1], 1e-14);
        EXPECT_NEAR(nodes[i][2], expected[2], 1e-14);

        std::vector<Real> vals;
        basis.evaluate_values(expected, vals);
        ASSERT_EQ(vals.size(), nodes.size());
        for (std::size_t j = 0; j < vals.size(); ++j) {
            const double expected_delta = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(vals[j], expected_delta, 1e-12);
        }
    }
}

void expect_alias_matches_canonical(ElementType alias_type,
                                    ElementType canonical_type,
                                    int canonical_order,
                                    const std::vector<Point>& points,
                                    Real tol = Real(1e-12)) {
    LagrangeBasis alias(alias_type, canonical_order);
    LagrangeBasis canonical(canonical_type, canonical_order);

    ASSERT_EQ(alias.element_type(), canonical.element_type());
    ASSERT_EQ(alias.order(), canonical.order());
    ASSERT_EQ(alias.size(), canonical.size());
    ASSERT_EQ(alias.nodes().size(), canonical.nodes().size());

    for (std::size_t i = 0; i < alias.nodes().size(); ++i) {
        EXPECT_NEAR(alias.nodes()[i][0], canonical.nodes()[i][0], tol);
        EXPECT_NEAR(alias.nodes()[i][1], canonical.nodes()[i][1], tol);
        EXPECT_NEAR(alias.nodes()[i][2], canonical.nodes()[i][2], tol);
    }

    for (const auto& xi : points) {
        std::vector<Real> alias_values;
        std::vector<Real> canonical_values;
        std::vector<Gradient> alias_gradients;
        std::vector<Gradient> canonical_gradients;
        std::vector<Hessian> alias_hessians;
        std::vector<Hessian> canonical_hessians;

        alias.evaluate_values(xi, alias_values);
        canonical.evaluate_values(xi, canonical_values);
        alias.evaluate_gradients(xi, alias_gradients);
        canonical.evaluate_gradients(xi, canonical_gradients);
        alias.evaluate_hessians(xi, alias_hessians);
        canonical.evaluate_hessians(xi, canonical_hessians);

        ASSERT_EQ(alias_values.size(), canonical_values.size());
        ASSERT_EQ(alias_gradients.size(), canonical_gradients.size());
        ASSERT_EQ(alias_hessians.size(), canonical_hessians.size());

        for (std::size_t i = 0; i < alias_values.size(); ++i) {
            EXPECT_NEAR(alias_values[i], canonical_values[i], tol);
            for (int d = 0; d < canonical.dimension(); ++d) {
                const std::size_t sd = static_cast<std::size_t>(d);
                EXPECT_NEAR(alias_gradients[i][sd], canonical_gradients[i][sd], tol);
                for (int e = 0; e < canonical.dimension(); ++e) {
                    const std::size_t se = static_cast<std::size_t>(e);
                    EXPECT_NEAR(alias_hessians[i](sd, se), canonical_hessians[i](sd, se), Real(5) * tol);
                }
            }
        }
    }
}

std::vector<Point> sample_points_for(ElementType type) {
    switch (type) {
        case ElementType::Line2:
        case ElementType::Line3:
            return {
                Point{Real(-0.7), Real(0), Real(0)},
                Point{Real(0.1), Real(0), Real(0)},
                Point{Real(0.65), Real(0), Real(0)}
            };
        case ElementType::Triangle3:
        case ElementType::Triangle6:
            return {
                Point{Real(0.15), Real(0.2), Real(0)},
                Point{Real(0.25), Real(0.1), Real(0)},
                Point{Real(0.2), Real(0.3), Real(0)}
            };
        case ElementType::Quad4:
        case ElementType::Quad9:
            return {
                Point{Real(0.2), Real(-0.35), Real(0)},
                Point{Real(-0.4), Real(0.25), Real(0)},
                Point{Real(0.55), Real(0.1), Real(0)}
            };
        case ElementType::Tetra4:
        case ElementType::Tetra10:
            return {
                Point{Real(0.1), Real(0.2), Real(0.15)},
                Point{Real(0.2), Real(0.1), Real(0.25)},
                Point{Real(0.15), Real(0.15), Real(0.2)}
            };
        case ElementType::Hex8:
        case ElementType::Hex27:
            return {
                Point{Real(0.2), Real(-0.3), Real(0.25)},
                Point{Real(-0.5), Real(0.4), Real(-0.2)},
                Point{Real(0.1), Real(0.15), Real(0.6)}
            };
        case ElementType::Wedge6:
        case ElementType::Wedge18:
            return {
                Point{Real(0.2), Real(0.25), Real(0.0)},
                Point{Real(0.1), Real(0.2), Real(-0.45)},
                Point{Real(0.3), Real(0.15), Real(0.5)}
            };
        case ElementType::Pyramid5:
        case ElementType::Pyramid14:
            return {
                Point{Real(0.0), Real(0.0), Real(0.25)},
                Point{Real(0.15), Real(-0.1), Real(0.3)},
                Point{Real(-0.1), Real(0.2), Real(0.4)}
            };
        default:
            return {Point{Real(0), Real(0), Real(0)}};
    }
}

std::vector<Point> boundary_stress_points_for(ElementType type);

std::vector<Point> dense_sample_points_for(ElementType type) {
    const auto interior = sample_points_for(type);
    const auto boundary = boundary_stress_points_for(type);

    std::vector<Point> points;
    points.reserve(interior.size() + boundary.size());
    points.insert(points.end(), interior.begin(), interior.end());
    points.insert(points.end(), boundary.begin(), boundary.end());

    if (type == ElementType::Pyramid5 || type == ElementType::Pyramid14) {
        points.push_back(Point{Real(0.0), Real(0.0), Real(0.85)});
        points.push_back(Point{Real(0.02), Real(-0.015), Real(0.95)});
    }
    return points;
}

std::vector<Point> boundary_stress_points_for(ElementType type) {
    switch (type) {
        case ElementType::Line2:
        case ElementType::Line3:
            return {
                Point{Real(-0.999), Real(0), Real(0)},
                Point{Real(-0.75), Real(0), Real(0)},
                Point{Real(0.0), Real(0), Real(0)},
                Point{Real(0.8), Real(0), Real(0)},
                Point{Real(0.999), Real(0), Real(0)}
            };
        case ElementType::Triangle3:
        case ElementType::Triangle6:
            return {
                Point{Real(1e-6), Real(1e-6), Real(0)},
                Point{Real(0.98), Real(0.01), Real(0)},
                Point{Real(0.01), Real(0.98), Real(0)},
                Point{Real(0.25), Real(1e-4), Real(0)},
                Point{Real(0.49), Real(0.49), Real(0)}
            };
        case ElementType::Quad4:
        case ElementType::Quad9:
            return {
                Point{Real(-0.99), Real(-0.99), Real(0)},
                Point{Real(0.99), Real(-0.99), Real(0)},
                Point{Real(0.99), Real(0.99), Real(0)},
                Point{Real(-0.99), Real(0.99), Real(0)},
                Point{Real(0.0), Real(0.95), Real(0)}
            };
        case ElementType::Tetra4:
        case ElementType::Tetra10:
            return {
                Point{Real(1e-6), Real(1e-6), Real(1e-6)},
                Point{Real(0.97), Real(0.01), Real(0.01)},
                Point{Real(0.01), Real(0.97), Real(0.01)},
                Point{Real(0.01), Real(0.01), Real(0.97)},
                Point{Real(0.32), Real(0.33), Real(0.01)}
            };
        case ElementType::Hex8:
        case ElementType::Hex27:
            return {
                Point{Real(-0.99), Real(-0.99), Real(-0.99)},
                Point{Real(0.99), Real(-0.99), Real(0.99)},
                Point{Real(0.99), Real(0.99), Real(-0.99)},
                Point{Real(-0.99), Real(0.99), Real(0.99)},
                Point{Real(0.0), Real(0.0), Real(0.95)}
            };
        case ElementType::Wedge6:
        case ElementType::Wedge18:
            return {
                Point{Real(1e-6), Real(1e-6), Real(-0.99)},
                Point{Real(0.98), Real(0.01), Real(-0.99)},
                Point{Real(0.01), Real(0.98), Real(0.99)},
                Point{Real(0.49), Real(0.49), Real(0.0)},
                Point{Real(0.25), Real(1e-4), Real(0.95)}
            };
        case ElementType::Pyramid5:
        case ElementType::Pyramid14:
            return {
                Point{Real(0.0), Real(0.0), Real(0.95)},
                Point{Real(0.01), Real(-0.01), Real(0.98)},
                Point{Real(0.6), Real(-0.6), Real(0.2)},
                Point{Real(0.79), Real(0.0), Real(0.2)},
                Point{Real(0.0), Real(0.79), Real(0.2)}
            };
        default:
            return {Point{Real(0), Real(0), Real(0)}};
    }
}

Real monomial_value(const Point& xi, int px, int py, int pz) {
    return std::pow(xi[0], px) * std::pow(xi[1], py) * std::pow(xi[2], pz);
}

void expect_gradients_match_finite_difference(const LagrangeAccuracyCase& c,
                                              Real eps,
                                              Real tol) {
    LagrangeBasis basis(c.type, c.order);

    for (const auto& xi : c.points) {
        std::vector<Gradient> gradients;
        basis.evaluate_gradients(xi, gradients);
        ASSERT_EQ(gradients.size(), basis.size());

        for (int d = 0; d < basis.dimension(); ++d) {
            Point xp = xi;
            Point xm = xi;
            xp[d] += eps;
            xm[d] -= eps;

            std::vector<Real> values_p;
            std::vector<Real> values_m;
            basis.evaluate_values(xp, values_p);
            basis.evaluate_values(xm, values_m);

            ASSERT_EQ(values_p.size(), basis.size());
            ASSERT_EQ(values_m.size(), basis.size());
            for (std::size_t i = 0; i < basis.size(); ++i) {
                const Real fd = (values_p[i] - values_m[i]) / (Real(2) * eps);
                EXPECT_NEAR(gradients[i][d], fd, tol)
                    << "type=" << static_cast<int>(c.type)
                    << ", order=" << c.order
                    << ", dim=" << d
                    << ", basis_i=" << i
                    << ", xi=(" << xi[0] << "," << xi[1] << "," << xi[2] << ")";
            }
        }
    }
}

void expect_polynomial_reproduction(const LagrangeAccuracyCase& c,
                                    const std::vector<std::array<int, 3>>& exponents,
                                    Real tol) {
    LagrangeBasis basis(c.type, c.order);
    const auto& nodes = basis.nodes();
    ASSERT_EQ(nodes.size(), basis.size());

    for (const auto& exp : exponents) {
        std::vector<Real> coeffs(basis.size(), Real(0));
        for (std::size_t i = 0; i < basis.size(); ++i) {
            coeffs[i] = monomial_value(nodes[i], exp[0], exp[1], exp[2]);
        }

        for (const auto& xi : c.points) {
            std::vector<Real> values;
            basis.evaluate_values(xi, values);
            ASSERT_EQ(values.size(), basis.size());

            Real interpolated = Real(0);
            for (std::size_t i = 0; i < basis.size(); ++i) {
                interpolated += coeffs[i] * values[i];
            }

            const Real exact = monomial_value(xi, exp[0], exp[1], exp[2]);
            EXPECT_NEAR(interpolated, exact, tol)
                << "type=" << static_cast<int>(c.type)
                << ", order=" << c.order
                << ", monomial=(" << exp[0] << "," << exp[1] << "," << exp[2] << ")"
                << ", xi=(" << xi[0] << "," << xi[1] << "," << xi[2] << ")";
        }
    }
}

template<typename Container>
void expect_all_finite(const Container& values) {
    for (const auto& value : values) {
        for (std::size_t d = 0; d < 3; ++d) {
            EXPECT_TRUE(std::isfinite(value[d]));
        }
    }
}

void expect_hessians_finite(const std::vector<Hessian>& hessians,
                            int dimension) {
    for (const auto& H : hessians) {
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                EXPECT_TRUE(std::isfinite(H(static_cast<std::size_t>(i),
                                            static_cast<std::size_t>(j))));
            }
        }
    }
}

void expect_partition_gradient_hessian_sums(const LagrangeBasis& basis,
                                            const std::vector<Point>& points,
                                            Real value_tol,
                                            Real derivative_tol) {
    for (const auto& xi : points) {
        std::vector<Real> values;
        std::vector<Gradient> gradients;
        std::vector<Hessian> hessians;
        basis.evaluate_values(xi, values);
        basis.evaluate_gradients(xi, gradients);
        basis.evaluate_hessians(xi, hessians);

        ASSERT_EQ(values.size(), basis.size());
        ASSERT_EQ(gradients.size(), basis.size());
        ASSERT_EQ(hessians.size(), basis.size());

        Real value_sum = Real(0);
        Gradient gradient_sum{};
        Hessian hessian_sum{};
        for (std::size_t i = 0; i < basis.size(); ++i) {
            value_sum += values[i];
            for (int d = 0; d < basis.dimension(); ++d) {
                const std::size_t sd = static_cast<std::size_t>(d);
                gradient_sum[sd] += gradients[i][sd];
                for (int e = 0; e < basis.dimension(); ++e) {
                    const std::size_t se = static_cast<std::size_t>(e);
                    hessian_sum(sd, se) += hessians[i](sd, se);
                }
            }
        }

        EXPECT_NEAR(value_sum, Real(1), value_tol)
            << "Element type " << static_cast<int>(basis.element_type())
            << ", order " << basis.order()
            << ", xi=(" << xi[0] << "," << xi[1] << "," << xi[2] << ")";

        for (int d = 0; d < basis.dimension(); ++d) {
            const std::size_t sd = static_cast<std::size_t>(d);
            EXPECT_NEAR(gradient_sum[sd], Real(0), derivative_tol)
                << "Gradient sum mismatch for element type " << static_cast<int>(basis.element_type())
                << ", order " << basis.order()
                << ", dim " << d;
            for (int e = 0; e < basis.dimension(); ++e) {
                const std::size_t se = static_cast<std::size_t>(e);
                EXPECT_NEAR(hessian_sum(sd, se), Real(0), derivative_tol)
                    << "Hessian sum mismatch for element type " << static_cast<int>(basis.element_type())
                    << ", order " << basis.order()
                    << ", component (" << d << "," << e << ")";
            }
        }
    }
}

bool is_on_pyramid_face(const Point& point,
                        PyramidFace face,
                        Real tol = Real(1e-12)) {
    const Real scale = Real(1) - point[2];
    switch (face) {
        case PyramidFace::Base:
            return std::abs(point[2]) <= tol;
        case PyramidFace::South:
            return std::abs(point[1] + scale) <= tol;
        case PyramidFace::East:
            return std::abs(point[0] - scale) <= tol;
        case PyramidFace::North:
            return std::abs(point[1] - scale) <= tol;
        case PyramidFace::West:
            return std::abs(point[0] + scale) <= tol;
    }
    return false;
}

Point map_pyramid_face_to_reference(PyramidFace face,
                                    const Point& point) {
    const Real scale = Real(1) - point[2];
    switch (face) {
        case PyramidFace::Base:
            return Point{point[0], point[1], Real(0)};
        case PyramidFace::South:
            return Point{(scale - point[0]) / Real(2), point[2], Real(0)};
        case PyramidFace::East:
            return Point{(scale + point[1]) / Real(2), point[2], Real(0)};
        case PyramidFace::North:
            return Point{(scale + point[0]) / Real(2), point[2], Real(0)};
        case PyramidFace::West:
            return Point{(scale - point[1]) / Real(2), point[2], Real(0)};
    }
    return Point{};
}

std::vector<Point> sample_points_for_pyramid_face(PyramidFace face) {
    switch (face) {
        case PyramidFace::Base:
            return {
                Point{Real(0.15), Real(-0.2), Real(0)},
                Point{Real(-0.55), Real(0.35), Real(0)}
            };
        case PyramidFace::South:
            return {
                Point{Real(-0.2), Real(-0.8), Real(0.2)},
                Point{Real(0.05), Real(-0.35), Real(0.65)}
            };
        case PyramidFace::East:
            return {
                Point{Real(0.8), Real(-0.25), Real(0.2)},
                Point{Real(0.3), Real(0.08), Real(0.7)}
            };
        case PyramidFace::North:
            return {
                Point{Real(0.25), Real(0.8), Real(0.2)},
                Point{Real(-0.08), Real(0.35), Real(0.65)}
            };
        case PyramidFace::West:
            return {
                Point{Real(-0.8), Real(0.2), Real(0.2)},
                Point{Real(-0.3), Real(-0.05), Real(0.7)}
            };
    }
    return {};
}

bool is_on_pyramid_edge(const Point& point,
                        PyramidEdge edge,
                        Real tol = Real(1e-12)) {
    const Real scale = Real(1) - point[2];
    switch (edge) {
        case PyramidEdge::BaseSouth:
            return std::abs(point[2]) <= tol && std::abs(point[1] + Real(1)) <= tol;
        case PyramidEdge::BaseEast:
            return std::abs(point[2]) <= tol && std::abs(point[0] - Real(1)) <= tol;
        case PyramidEdge::BaseNorth:
            return std::abs(point[2]) <= tol && std::abs(point[1] - Real(1)) <= tol;
        case PyramidEdge::BaseWest:
            return std::abs(point[2]) <= tol && std::abs(point[0] + Real(1)) <= tol;
        case PyramidEdge::VerticalSW:
            return std::abs(point[0] + scale) <= tol && std::abs(point[1] + scale) <= tol;
        case PyramidEdge::VerticalSE:
            return std::abs(point[0] - scale) <= tol && std::abs(point[1] + scale) <= tol;
        case PyramidEdge::VerticalNE:
            return std::abs(point[0] - scale) <= tol && std::abs(point[1] - scale) <= tol;
        case PyramidEdge::VerticalNW:
            return std::abs(point[0] + scale) <= tol && std::abs(point[1] - scale) <= tol;
    }
    return false;
}

Point map_pyramid_edge_to_reference(PyramidEdge edge,
                                    const Point& point) {
    switch (edge) {
        case PyramidEdge::BaseSouth:
        case PyramidEdge::BaseNorth:
            return Point{point[0], Real(0), Real(0)};
        case PyramidEdge::BaseEast:
        case PyramidEdge::BaseWest:
            return Point{point[1], Real(0), Real(0)};
        case PyramidEdge::VerticalSW:
        case PyramidEdge::VerticalSE:
        case PyramidEdge::VerticalNE:
        case PyramidEdge::VerticalNW:
            return Point{Real(2) * point[2] - Real(1), Real(0), Real(0)};
    }
    return Point{};
}

std::vector<Point> sample_points_for_pyramid_edge(PyramidEdge edge) {
    switch (edge) {
        case PyramidEdge::BaseSouth:
            return {Point{Real(-0.65), Real(-1), Real(0)}, Point{Real(0.35), Real(-1), Real(0)}};
        case PyramidEdge::BaseEast:
            return {Point{Real(1), Real(-0.45), Real(0)}, Point{Real(1), Real(0.55), Real(0)}};
        case PyramidEdge::BaseNorth:
            return {Point{Real(-0.55), Real(1), Real(0)}, Point{Real(0.45), Real(1), Real(0)}};
        case PyramidEdge::BaseWest:
            return {Point{Real(-1), Real(-0.55), Real(0)}, Point{Real(-1), Real(0.45), Real(0)}};
        case PyramidEdge::VerticalSW:
            return {Point{Real(-0.75), Real(-0.75), Real(0.25)}, Point{Real(-0.3), Real(-0.3), Real(0.7)}};
        case PyramidEdge::VerticalSE:
            return {Point{Real(0.75), Real(-0.75), Real(0.25)}, Point{Real(0.3), Real(-0.3), Real(0.7)}};
        case PyramidEdge::VerticalNE:
            return {Point{Real(0.75), Real(0.75), Real(0.25)}, Point{Real(0.3), Real(0.3), Real(0.7)}};
        case PyramidEdge::VerticalNW:
            return {Point{Real(-0.75), Real(0.75), Real(0.25)}, Point{Real(-0.3), Real(0.3), Real(0.7)}};
    }
    return {};
}

std::vector<int> map_pyramid_nodes_to_lower_basis_nodes(
    const std::vector<Point>& pyramid_nodes,
    const std::vector<Point>& lower_basis_nodes,
    const std::function<bool(const Point&)>& selector,
    const std::function<Point(const Point&)>& mapper) {
    std::vector<int> mapping(pyramid_nodes.size(), -1);
    std::size_t face_count = 0;
    for (std::size_t i = 0; i < pyramid_nodes.size(); ++i) {
        if (!selector(pyramid_nodes[i])) {
            continue;
        }

        ++face_count;
        const Point mapped = mapper(pyramid_nodes[i]);
        bool found = false;
        for (std::size_t j = 0; j < lower_basis_nodes.size(); ++j) {
            if (points_close(mapped, lower_basis_nodes[j])) {
                mapping[i] = static_cast<int>(j);
                found = true;
                break;
            }
        }
        EXPECT_TRUE(found)
            << "Failed to match pyramid trace node at (" << pyramid_nodes[i][0] << ","
            << pyramid_nodes[i][1] << "," << pyramid_nodes[i][2] << ")";
    }

    EXPECT_EQ(face_count, lower_basis_nodes.size());
    return mapping;
}

void expect_pyramid_face_trace_matches_lower_basis(int order,
                                                   PyramidFace face,
                                                   Real tol = Real(2e-10)) {
    LagrangeBasis pyramid(ElementType::Pyramid5, order);
    const bool base_face = face == PyramidFace::Base;
    LagrangeBasis lower(base_face ? ElementType::Quad4 : ElementType::Triangle3, order);

    const auto mapping = map_pyramid_nodes_to_lower_basis_nodes(
        pyramid.nodes(),
        lower.nodes(),
        [&](const Point& point) { return is_on_pyramid_face(point, face); },
        [&](const Point& point) { return map_pyramid_face_to_reference(face, point); });

    for (const auto& face_point : sample_points_for_pyramid_face(face)) {
        std::vector<Real> pyramid_values;
        std::vector<Real> lower_values;
        pyramid.evaluate_values(face_point, pyramid_values);
        lower.evaluate_values(map_pyramid_face_to_reference(face, face_point), lower_values);

        ASSERT_EQ(pyramid_values.size(), pyramid.size());
        ASSERT_EQ(lower_values.size(), lower.size());

        for (std::size_t i = 0; i < pyramid.size(); ++i) {
            if (mapping[i] >= 0) {
                EXPECT_NEAR(pyramid_values[i], lower_values[static_cast<std::size_t>(mapping[i])], tol)
                    << "Face trace mismatch for order " << order
                    << ", face " << static_cast<int>(face)
                    << ", basis " << i;
            } else {
                EXPECT_NEAR(pyramid_values[i], Real(0), tol)
                    << "Off-face pyramid basis should vanish on face for order " << order
                    << ", face " << static_cast<int>(face)
                    << ", basis " << i;
            }
        }
    }
}

void expect_pyramid_edge_trace_matches_line_basis(int order,
                                                  PyramidEdge edge,
                                                  Real tol = Real(2e-10)) {
    LagrangeBasis pyramid(ElementType::Pyramid5, order);
    LagrangeBasis line(ElementType::Line2, order);

    const auto mapping = map_pyramid_nodes_to_lower_basis_nodes(
        pyramid.nodes(),
        line.nodes(),
        [&](const Point& point) { return is_on_pyramid_edge(point, edge); },
        [&](const Point& point) { return map_pyramid_edge_to_reference(edge, point); });

    for (const auto& edge_point : sample_points_for_pyramid_edge(edge)) {
        std::vector<Real> pyramid_values;
        std::vector<Real> line_values;
        pyramid.evaluate_values(edge_point, pyramid_values);
        line.evaluate_values(map_pyramid_edge_to_reference(edge, edge_point), line_values);

        ASSERT_EQ(pyramid_values.size(), pyramid.size());
        ASSERT_EQ(line_values.size(), line.size());

        for (std::size_t i = 0; i < pyramid.size(); ++i) {
            if (mapping[i] >= 0) {
                EXPECT_NEAR(pyramid_values[i], line_values[static_cast<std::size_t>(mapping[i])], tol)
                    << "Edge trace mismatch for order " << order
                    << ", edge " << static_cast<int>(edge)
                    << ", basis " << i;
            } else {
                EXPECT_NEAR(pyramid_values[i], Real(0), tol)
                    << "Off-edge pyramid basis should vanish on edge for order " << order
                    << ", edge " << static_cast<int>(edge)
                    << ", basis " << i;
            }
        }
    }
}

struct StridedOutputRequest {
    bool values;
    bool gradients;
    bool hessians;
};

void expect_strided_matches_pointwise(ElementType type,
                                      int order,
                                      const StridedOutputRequest& request) {
    LagrangeBasis basis(type, order);
    const auto points = dense_sample_points_for(type);
    const std::size_t stride = points.size() + 3u;
    constexpr Real sentinel = Real(-12345.25);

    std::vector<Real> values(request.values ? basis.size() * stride : 0u, sentinel);
    std::vector<Real> gradients(request.gradients ? basis.size() * 3u * stride : 0u, sentinel);
    std::vector<Real> hessians(request.hessians ? basis.size() * 9u * stride : 0u, sentinel);

    basis.evaluate_at_quadrature_points_strided(
        points,
        stride,
        request.values ? values.data() : nullptr,
        request.gradients ? gradients.data() : nullptr,
        request.hessians ? hessians.data() : nullptr);

    const Real tol = (type == ElementType::Pyramid5 || type == ElementType::Pyramid14)
        ? Real(5e-10)
        : Real(1e-12);

    for (std::size_t q = 0; q < points.size(); ++q) {
        if (request.values) {
            std::vector<Real> expected;
            basis.evaluate_values(points[q], expected);
            ASSERT_EQ(expected.size(), basis.size());
            for (std::size_t d = 0; d < basis.size(); ++d) {
                EXPECT_NEAR(values[d * stride + q], expected[d], tol)
                    << "type=" << static_cast<int>(type)
                    << ", order=" << order
                    << ", dof=" << d
                    << ", q=" << q;
            }
        }

        if (request.gradients) {
            std::vector<Gradient> expected;
            basis.evaluate_gradients(points[q], expected);
            ASSERT_EQ(expected.size(), basis.size());
            for (std::size_t d = 0; d < basis.size(); ++d) {
                for (std::size_t c = 0; c < 3u; ++c) {
                    EXPECT_NEAR(gradients[(d * 3u + c) * stride + q], expected[d][c], tol)
                        << "type=" << static_cast<int>(type)
                        << ", order=" << order
                        << ", dof=" << d
                        << ", component=" << c
                        << ", q=" << q;
                }
            }
        }

        if (request.hessians) {
            std::vector<Hessian> expected;
            basis.evaluate_hessians(points[q], expected);
            ASSERT_EQ(expected.size(), basis.size());
            for (std::size_t d = 0; d < basis.size(); ++d) {
                for (std::size_t r = 0; r < 3u; ++r) {
                    for (std::size_t c = 0; c < 3u; ++c) {
                        EXPECT_NEAR(hessians[(d * 9u + r * 3u + c) * stride + q],
                                    expected[d](r, c),
                                    Real(4) * tol)
                            << "type=" << static_cast<int>(type)
                            << ", order=" << order
                            << ", dof=" << d
                            << ", hessian=(" << r << "," << c << ")"
                            << ", q=" << q;
                    }
                }
            }
        }
    }

    const auto expect_padding_untouched = [&](const std::vector<Real>& buffer,
                                              std::size_t rows) {
        for (std::size_t row = 0; row < rows; ++row) {
            for (std::size_t q = points.size(); q < stride; ++q) {
                EXPECT_EQ(buffer[row * stride + q], sentinel)
                    << "type=" << static_cast<int>(type)
                    << ", order=" << order
                    << ", row=" << row
                    << ", padding q=" << q;
            }
        }
    };

    if (request.values) {
        expect_padding_untouched(values, basis.size());
    }
    if (request.gradients) {
        expect_padding_untouched(gradients, basis.size() * 3u);
    }
    if (request.hessians) {
        expect_padding_untouched(hessians, basis.size() * 9u);
    }
}

void expect_raw_to_matches_vector_evaluation(ElementType type, int order) {
    LagrangeBasis basis(type, order);
    const Real tol = (type == ElementType::Pyramid5 || type == ElementType::Pyramid14)
        ? Real(5e-10)
        : Real(1e-12);

    for (const auto& point : sample_points_for(type)) {
        std::vector<Real> values;
        std::vector<Gradient> gradients;
        std::vector<Hessian> hessians;
        basis.evaluate_all(point, values, gradients, hessians);

        std::vector<Real> raw_values(basis.size());
        std::vector<Real> raw_gradients(basis.size() * 3u);
        std::vector<Real> raw_hessians(basis.size() * 9u);
        basis.evaluate_values_to(point, raw_values.data());
        basis.evaluate_gradients_to(point, raw_gradients.data());
        basis.evaluate_hessians_to(point, raw_hessians.data());

        for (std::size_t i = 0; i < basis.size(); ++i) {
            EXPECT_NEAR(raw_values[i], values[i], tol)
                << "type=" << static_cast<int>(type) << ", order=" << order << ", dof=" << i;
            for (std::size_t c = 0; c < 3u; ++c) {
                EXPECT_NEAR(raw_gradients[i * 3u + c], gradients[i][c], tol)
                    << "type=" << static_cast<int>(type)
                    << ", order=" << order
                    << ", dof=" << i
                    << ", gradient component=" << c;
            }
            for (std::size_t r = 0; r < 3u; ++r) {
                for (std::size_t c = 0; c < 3u; ++c) {
                    EXPECT_NEAR(raw_hessians[i * 9u + r * 3u + c], hessians[i](r, c), Real(4) * tol)
                        << "type=" << static_cast<int>(type)
                        << ", order=" << order
                        << ", dof=" << i
                        << ", hessian=(" << r << "," << c << ")";
                }
            }
        }
    }
}

} // namespace

TEST(SolverBasisAdapter, ShapeValuesGradientsAndMeshOverloadMatchLegacy) {
    constexpr double tol = 2e-12;

    for (const auto& c : solver_basis_adapter_cases()) {
        SCOPED_TRACE("element=" + std::to_string(static_cast<int>(c.type)));
        Vector<double> weights;
        Array<double> xi;
        fill_legacy_quadrature(c, weights, xi);

        Array<double> legacy_N(c.eNoN, c.nG);
        Array<double> adapter_N(c.eNoN, c.nG);
        Array3<double> legacy_Nx(c.insd, c.eNoN, c.nG);
        Array3<double> adapter_Nx(c.insd, c.eNoN, c.nG);
        auto legacy_shape = legacy_solver_nn::get_element_shape_data.find(c.type);

        faceType legacy_face;
        if (legacy_shape == legacy_solver_nn::get_element_shape_data.end()) {
            ASSERT_EQ(c.type, consts::ElementType::QUD8);
            legacy_face.eType = c.type;
            legacy_face.eNoN = c.eNoN;
            legacy_face.nG = c.nG;
            legacy_face.xi = xi;
            legacy_face.N.resize(c.eNoN, c.nG);
            legacy_face.Nx.resize(c.insd, c.eNoN, c.nG);
        }

        for (int g = 0; g < c.nG; ++g) {
            if (legacy_shape != legacy_solver_nn::get_element_shape_data.end()) {
                legacy_shape->second(c.insd, c.eNoN, g, xi, legacy_N, legacy_Nx);
            } else {
                legacy_solver_nn::set_face_shape_data.at(c.type)(g, legacy_face);
            }
            nn::get_gnn(c.insd, c.type, c.eNoN, g, xi, adapter_N, adapter_Nx);

            double partition = 0.0;
            std::array<double, 3> gradient_sum{0.0, 0.0, 0.0};
            for (int a = 0; a < c.eNoN; ++a) {
                partition += adapter_N(a, g);
                for (int d = 0; d < c.insd; ++d) {
                    gradient_sum[static_cast<std::size_t>(d)] += adapter_Nx(d, a, g);
                }
            }

            EXPECT_NEAR(partition, 1.0, tol)
                << "element=" << static_cast<int>(c.type) << ", g=" << g;
            for (int d = 0; d < c.insd; ++d) {
                EXPECT_NEAR(gradient_sum[static_cast<std::size_t>(d)], 0.0, tol)
                    << "element=" << static_cast<int>(c.type) << ", d=" << d << ", g=" << g;
            }
        }

        if (legacy_shape == legacy_solver_nn::get_element_shape_data.end()) {
            legacy_N = legacy_face.N;
            legacy_Nx = legacy_face.Nx;
        }

        expect_arrays_near(adapter_N, legacy_N, tol);
        expect_array3_near(adapter_Nx, legacy_Nx, tol);

        mshType mesh;
        mesh.eType = c.type;
        mesh.eNoN = c.eNoN;
        mesh.nG = c.nG;
        mesh.xi = xi;
        mesh.N.resize(c.eNoN, c.nG);
        mesh.Nx.resize(c.insd, c.eNoN, c.nG);
        for (int g = 0; g < c.nG; ++g) {
            nn::get_gnn(g, mesh);
        }

        expect_arrays_near(mesh.N, legacy_N, tol);
        expect_array3_near(mesh.Nx, legacy_Nx, tol);
    }
}

TEST(SolverFaceBasisAdapter, ShapeValuesGradientsAndDispatchMatchLegacyFaceTable) {
    constexpr double tol = 2e-12;

    int covered = 0;
    for (const auto& c : solver_face_basis_adapter_cases()) {
        SCOPED_TRACE("face element=" + std::to_string(static_cast<int>(c.type)));

        faceType legacy_face = initialized_face_for_case(c);
        faceType basis_face = initialized_face_for_case(c);

        for (int g = 0; g < c.nG; ++g) {
            legacy_solver_nn::set_face_shape_data.at(c.type)(g, legacy_face);
            nn::get_gnn(nullptr, g, basis_face);
            expect_face_partition_identities(c, basis_face, g, tol);
        }

        expect_arrays_near(basis_face.N, legacy_face.N, tol);
        expect_array3_near(basis_face.Nx, legacy_face.Nx, tol);
        ++covered;
    }

    EXPECT_EQ(covered, 7);
}

TEST(SolverFaceBasisAdapter, MappedFacesFailClosedWithoutLegacyFallback) {
    using consts::ElementType;

    SolverBasisAdapterCase c{ElementType::LIN1, ElementType::LIN1, 1, 3, 2};
    faceType face = initialized_face_for_case(c);

    try {
        nn::get_gnn(nullptr, 0, face);
        FAIL() << "Expected mapped face dispatch to reject mismatched eNoN";
    } catch (const svmp::FE::basis::BasisEvaluationException& exception) {
        const std::string message = exception.what();
        EXPECT_NE(message.find("legacy fallback was not attempted"), std::string::npos)
            << message;
    }
}

TEST(SolverFaceBasisAdapter, PointFaceRemainsLegacyValuePath) {
    faceType face;
    face.eType = consts::ElementType::PNT;
    face.eNoN = 1;
    face.nG = 1;
    face.N.resize(1, 1);
    face.Nx.resize(1, 1, 1);
    face.N(0, 0) = -7.0;
    face.Nx(0, 0, 0) = 42.0;

    nn::get_gnn(nullptr, 0, face);

    EXPECT_DOUBLE_EQ(face.N(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(face.Nx(0, 0, 0), 42.0);
}

TEST(SolverFaceBasisAdapter, UnsupportedFacesThrowClearErrors) {
    faceType nrb_face;
    nrb_face.eType = consts::ElementType::NRB;
    nrb_face.eNoN = 1;
    nrb_face.nG = 1;
    nrb_face.N.resize(1, 1);
    nrb_face.Nx.resize(1, 1, 1);
    EXPECT_THROW(nn::get_gnn(nullptr, 0, nrb_face), svmp::FE::NotImplementedException);

    faceType unknown_face;
    unknown_face.eType = consts::ElementType::NA;
    unknown_face.eNoN = 1;
    unknown_face.nG = 1;
    unknown_face.N.resize(1, 1);
    unknown_face.Nx.resize(1, 1, 1);
    EXPECT_THROW(nn::get_gnn(nullptr, 0, unknown_face), svmp::FE::InvalidElementException);
}

TEST(SolverBasisAdapter, QuadraturePathsRemainLegacyCompatible) {
    constexpr double tol = 0.0;

    for (const auto& c : solver_basis_adapter_cases()) {
        auto mesh_it = legacy_solver_nn::set_element_gauss_int_data.find(c.type);
        if (mesh_it != legacy_solver_nn::set_element_gauss_int_data.end()) {
            mshType legacy_mesh;
            legacy_mesh.eType = c.type;
            legacy_mesh.eNoN = c.eNoN;
            legacy_mesh.nG = c.nG;
            legacy_mesh.w.resize(c.nG);
            legacy_mesh.xi.resize(c.insd, c.nG);
            mesh_it->second(legacy_mesh);

            mshType adapter_mesh;
            adapter_mesh.eType = c.type;
            adapter_mesh.eNoN = c.eNoN;
            adapter_mesh.nG = c.nG;
            adapter_mesh.w.resize(c.nG);
            adapter_mesh.xi.resize(c.insd, c.nG);
            nn::get_gip(adapter_mesh);

            expect_vectors_near(adapter_mesh.w, legacy_mesh.w, tol);
            expect_arrays_near(adapter_mesh.xi, legacy_mesh.xi, tol);
        }

        auto scalar_it = legacy_solver_nn::get_element_gauss_int_data.find(c.type);
        if (scalar_it != legacy_solver_nn::get_element_gauss_int_data.end()) {
            Vector<double> legacy_w(c.nG);
            Vector<double> adapter_w(c.nG);
            Array<double> legacy_xi(c.insd, c.nG);
            Array<double> adapter_xi(c.insd, c.nG);

            scalar_it->second(c.insd, c.nG, legacy_w, legacy_xi);
            nn::get_gip(c.insd, c.type, c.nG, adapter_w, adapter_xi);

            expect_vectors_near(adapter_w, legacy_w, tol);
            expect_arrays_near(adapter_xi, legacy_xi, tol);
        }
    }

    mshType legacy_tet;
    legacy_tet.eType = consts::ElementType::TET4;
    legacy_tet.eNoN = 4;
    legacy_tet.nG = 4;
    legacy_tet.qmTET4 = 0.25;
    legacy_tet.w.resize(4);
    legacy_tet.xi.resize(3, 4);
    legacy_solver_nn::set_element_gauss_int_data.at(consts::ElementType::TET4)(legacy_tet);

    mshType adapter_tet;
    adapter_tet.eType = consts::ElementType::TET4;
    adapter_tet.eNoN = 4;
    adapter_tet.nG = 4;
    adapter_tet.qmTET4 = 0.25;
    adapter_tet.w.resize(4);
    adapter_tet.xi.resize(3, 4);
    nn::get_gip(adapter_tet);

    expect_vectors_near(adapter_tet.w, legacy_tet.w, tol);
    expect_arrays_near(adapter_tet.xi, legacy_tet.xi, tol);
}

TEST(SolverBasisAdapter, HessiansCoverEveryMappedScalarVolumeElement) {
    constexpr double partition_tol = 2e-10;
    constexpr double finite_difference_tol = 2e-5;
    constexpr double zero_tol = 2e-12;

    int covered = 0;
    for (const auto& c : solver_hessian_adapter_cases()) {
        SCOPED_TRACE("element=" + std::to_string(static_cast<int>(c.type)));
        Vector<double> weights;
        Array<double> xi;
        fill_legacy_quadrature(c, weights, xi);

        const int ind2 = packed_hessian_components(c.insd);
        Array3<double> adapter_Nxx(ind2, c.eNoN, c.nG);
        fill_array3(adapter_Nxx, std::numeric_limits<double>::quiet_NaN());

        for (int g = 0; g < c.nG; ++g) {
            nn::get_gn_nxx(c.insd, ind2, c.type, c.eNoN, g, xi, adapter_Nxx);
            expect_partition_hessian_identity(c, adapter_Nxx, g, partition_tol);
            expect_solver_hessian_matches_gradient_finite_difference(
                c, xi, g, adapter_Nxx, finite_difference_tol);

            if (c.type == consts::ElementType::LIN1 ||
                c.type == consts::ElementType::TRI3 ||
                c.type == consts::ElementType::TET4) {
                expect_all_hessians_zero(c, adapter_Nxx, g, zero_tol);
            }
        }

        if (c.type == consts::ElementType::QUD4) {
            EXPECT_TRUE(array3_has_nonzero_component(adapter_Nxx, 2, zero_tol));
        } else if (c.type == consts::ElementType::HEX8) {
            EXPECT_TRUE(array3_has_nonzero_component(adapter_Nxx, 3, zero_tol));
            EXPECT_TRUE(array3_has_nonzero_component(adapter_Nxx, 4, zero_tol));
            EXPECT_TRUE(array3_has_nonzero_component(adapter_Nxx, 5, zero_tol));
        } else if (c.type == consts::ElementType::WDG) {
            EXPECT_TRUE(array3_has_nonzero_component(adapter_Nxx, 5, zero_tol));
        }
        ++covered;
    }

    EXPECT_EQ(covered, 13);
}

TEST(SolverBasisAdapter, HessianPackingMatchesLegacyWhereLegacyIsApproved) {
    constexpr double tol = 2e-12;

    for (const auto& c : solver_legacy_hessian_parity_cases()) {
        Vector<double> weights;
        Array<double> xi;
        fill_legacy_quadrature(c, weights, xi);

        const int ind2 = packed_hessian_components(c.insd);
        Array3<double> legacy_Nxx(ind2, c.eNoN, c.nG);
        Array3<double> adapter_Nxx(ind2, c.eNoN, c.nG);

        for (int g = 0; g < c.nG; ++g) {
            legacy_solver_nn::get_element_2nd_derivs.at(c.type)(
                c.insd, ind2, c.eNoN, g, xi, legacy_Nxx);
            nn::get_gn_nxx(c.insd, ind2, c.type, c.eNoN, g, xi, adapter_Nxx);
        }

        expect_array3_near(adapter_Nxx, legacy_Nxx, tol);
    }
}

TEST(SolverBasisAdapter, Qud8HessiansDoNotUseLegacyFallback) {
    using consts::ElementType;
    SolverBasisAdapterCase c{ElementType::QUD8, ElementType::QUD9, 2, 8, 9};

    Vector<double> weights;
    Array<double> xi;
    fill_legacy_quadrature(c, weights, xi);

    const int ind2 = packed_hessian_components(c.insd);
    Array3<double> legacy_Nxx(ind2, c.eNoN, c.nG);
    Array3<double> adapter_Nxx(ind2, c.eNoN, c.nG);
    fill_array3(legacy_Nxx, 0.0);
    fill_array3(adapter_Nxx, 0.0);

    for (int g = 0; g < c.nG; ++g) {
        legacy_solver_nn::get_element_2nd_derivs.at(c.type)(
            c.insd, ind2, c.eNoN, g, xi, legacy_Nxx);
        nn::get_gn_nxx(c.insd, ind2, c.type, c.eNoN, g, xi, adapter_Nxx);
    }

    double max_abs_difference = 0.0;
    for (int g = 0; g < c.nG; ++g) {
        for (int a = 0; a < c.eNoN; ++a) {
            for (int row = 0; row < ind2; ++row) {
                max_abs_difference = std::max(
                    max_abs_difference,
                    std::abs(adapter_Nxx(row, a, g) - legacy_Nxx(row, a, g)));
            }
        }
    }

    EXPECT_GT(max_abs_difference, 1e-8);
}

TEST(SolverBasisAdapter, UnsupportedHessianFamiliesRemainNoOp) {
    Array<double> xi(1, 1);
    xi(0, 0) = 0.0;
    Array3<double> Nxx(1, 1, 1);

    for (const auto unsupported : {consts::ElementType::NRB, consts::ElementType::PNT}) {
        fill_array3(Nxx, 42.0);
        nn::get_gn_nxx(1, 1, unsupported, 1, 0, xi, Nxx);
        EXPECT_DOUBLE_EQ(Nxx(0, 0, 0), 42.0)
            << "element=" << static_cast<int>(unsupported);
    }
}

TEST(SolverBasisAdapter, InitFsMshPopulatesMappedHessiansWithoutLShpFGate) {
    using consts::ElementType;
    const SolverBasisAdapterCase cases[] = {
        {ElementType::QUD4, ElementType::QUD4, 2, 4, 4},
        {ElementType::HEX8, ElementType::HEX8, 3, 8, 8},
        {ElementType::HEX20, ElementType::HEX20, 3, 20, 27},
        {ElementType::HEX27, ElementType::HEX27, 3, 27, 27},
        {ElementType::WDG, ElementType::WDG, 3, 6, 6},
    };

    for (const auto& c : cases) {
        SCOPED_TRACE("element=" + std::to_string(static_cast<int>(c.type)));
        ComMod com_mod;
        com_mod.nsd = c.insd;
        mshType mesh = initialized_mesh_for_case(c, true);

        fs::init_fs_msh(com_mod, mesh);

        ASSERT_EQ(mesh.fs.size(), 1u);
        ASSERT_EQ(mesh.fs[0].Nxx.nrows(), packed_hessian_components(c.insd));
        if (c.type == ElementType::QUD4) {
            EXPECT_TRUE(array3_has_nonzero_component(mesh.fs[0].Nxx, 2, 2e-12));
        } else if (c.type == ElementType::HEX8) {
            EXPECT_TRUE(array3_has_nonzero_component(mesh.fs[0].Nxx, 3, 2e-12));
        } else if (c.type == ElementType::WDG) {
            EXPECT_TRUE(array3_has_nonzero_component(mesh.fs[0].Nxx, 5, 2e-12));
        } else {
            bool has_nonzero = false;
            for (int row = 0; row < mesh.fs[0].Nxx.nrows(); ++row) {
                has_nonzero = has_nonzero ||
                    array3_has_nonzero_component(mesh.fs[0].Nxx, row, 2e-12);
            }
            EXPECT_TRUE(has_nonzero);
        }
    }
}

TEST(LagrangeBasis, QuadPartitionOfUnity) {
    LagrangeBasis basis(ElementType::Quad4, 1);
    svmp::FE::math::Vector<Real, 3> xi{0.2, -0.3, 0.0};

    std::vector<Real> values;
    basis.evaluate_values(xi, values);

    double sum = std::accumulate(values.begin(), values.end(), 0.0);
    EXPECT_NEAR(sum, 1.0, 1e-12);
}

TEST(LagrangeBasis, LineGradientLinear) {
    LagrangeBasis basis(ElementType::Line2, 1);
    svmp::FE::math::Vector<Real, 3> xi{0.0, 0.0, 0.0};
    std::vector<Gradient> grad;
    basis.evaluate_gradients(xi, grad);

    ASSERT_EQ(grad.size(), 2u);
    EXPECT_NEAR(grad[0][0], -0.5, 1e-12);
    EXPECT_NEAR(grad[1][0], 0.5, 1e-12);
}

TEST(LagrangeBasis, TrianglePartitionOfUnity) {
    LagrangeBasis basis(ElementType::Triangle3, 1);
    svmp::FE::math::Vector<Real, 3> xi{0.2, 0.3, 0.0};
    std::vector<Real> values;
    basis.evaluate_values(xi, values);

    double sum = std::accumulate(values.begin(), values.end(), 0.0);
    EXPECT_NEAR(sum, 1.0, 1e-12);
}

TEST(LagrangeBasis, SizeFormulasPerElement) {
    for (int order = 0; order <= 3; ++order) {
        {
            LagrangeBasis line(ElementType::Line2, order);
            EXPECT_EQ(line.size(), static_cast<std::size_t>(order + 1));
        }
        {
            LagrangeBasis quad(ElementType::Quad4, order);
            const std::size_t n1d = static_cast<std::size_t>(order + 1);
            EXPECT_EQ(quad.size(), n1d * n1d);
        }
        {
            LagrangeBasis hex(ElementType::Hex8, order);
            const std::size_t n1d = static_cast<std::size_t>(order + 1);
            EXPECT_EQ(hex.size(), n1d * n1d * n1d);
        }
        {
            LagrangeBasis tri(ElementType::Triangle3, order);
            const std::size_t expected =
                static_cast<std::size_t>(order + 1) *
                static_cast<std::size_t>(order + 2) / 2;
            EXPECT_EQ(tri.size(), expected);
        }
        {
            LagrangeBasis tet(ElementType::Tetra4, order);
            const std::size_t expected =
                static_cast<std::size_t>(order + 1) *
                static_cast<std::size_t>(order + 2) *
                static_cast<std::size_t>(order + 3) / 6;
            EXPECT_EQ(tet.size(), expected);
        }
    }
}

TEST(LagrangeBasis, KroneckerDeltaAtNodes) {
    const std::vector<std::pair<ElementType, int>> cases = {
        {ElementType::Line2, 1},
        {ElementType::Quad4, 1},
        {ElementType::Triangle3, 1},
        {ElementType::Tetra4, 1},
        {ElementType::Hex8, 1},
        {ElementType::Triangle3, 2},
        {ElementType::Tetra4, 2},
        {ElementType::Quad4, 2},
        {ElementType::Hex8, 2},
        {ElementType::Wedge6, 2}
    };

    for (const auto& c : cases) {
        LagrangeBasis basis(c.first, c.second);
        const auto& nodes = basis.nodes();
        ASSERT_EQ(nodes.size(), basis.size());

        for (std::size_t i = 0; i < nodes.size(); ++i) {
            std::vector<Real> vals;
            basis.evaluate_values(nodes[i], vals);
            ASSERT_EQ(vals.size(), nodes.size());
            for (std::size_t j = 0; j < nodes.size(); ++j) {
                if (i == j) {
                    EXPECT_NEAR(vals[j], 1.0, 1e-12);
                } else {
                    EXPECT_NEAR(vals[j], 0.0, 1e-12);
                }
            }
        }
    }
}

TEST(LagrangeBasis, MatchesNodeOrderingConventionsForLinearAndQuadratic) {
    // Tensor-product elements
    expect_nodes_match_node_ordering(ElementType::Line2, 1, ElementType::Line2);
    expect_nodes_match_node_ordering(ElementType::Line2, 2, ElementType::Line3);
    expect_nodes_match_node_ordering(ElementType::Quad4, 1, ElementType::Quad4);
    expect_nodes_match_node_ordering(ElementType::Quad4, 2, ElementType::Quad9);
    expect_nodes_match_node_ordering(ElementType::Hex8, 1, ElementType::Hex8);
    expect_nodes_match_node_ordering(ElementType::Hex8, 2, ElementType::Hex27);

    // Simplex elements
    expect_nodes_match_node_ordering(ElementType::Triangle3, 1, ElementType::Triangle3);
    expect_nodes_match_node_ordering(ElementType::Triangle3, 2, ElementType::Triangle6);
    expect_nodes_match_node_ordering(ElementType::Tetra4, 1, ElementType::Tetra4);
    expect_nodes_match_node_ordering(ElementType::Tetra4, 2, ElementType::Tetra10);

    // Mixed topology
    expect_nodes_match_node_ordering(ElementType::Wedge6, 1, ElementType::Wedge6);
    expect_nodes_match_node_ordering(ElementType::Wedge6, 2, ElementType::Wedge18);

    // Pyramid
    expect_nodes_match_node_ordering(ElementType::Pyramid5, 1, ElementType::Pyramid5);
    expect_nodes_match_node_ordering(ElementType::Pyramid14, 2, ElementType::Pyramid14);
}

TEST(LagrangeBasis, WedgeAndPyramidPartitionOfUnity) {
    {
        LagrangeBasis wedge(ElementType::Wedge6, 1);
        svmp::FE::math::Vector<Real, 3> xi{Real(0.2), Real(0.1), Real(0.3)};
        std::vector<Real> vals;
        wedge.evaluate_values(xi, vals);
        const double sum = std::accumulate(vals.begin(), vals.end(), 0.0);
        EXPECT_NEAR(sum, 1.0, 1e-12);
    }

    {
        LagrangeBasis wedge_q(ElementType::Wedge18, 2);
        svmp::FE::math::Vector<Real, 3> xi{Real(0.2), Real(0.1), Real(-0.25)};
        std::vector<Real> vals;
        wedge_q.evaluate_values(xi, vals);
        const double sum = std::accumulate(vals.begin(), vals.end(), 0.0);
        EXPECT_NEAR(sum, 1.0, 1e-12);

        // Wedge18 should report 18 nodes in ReferenceNodeLayout
        EXPECT_EQ(ReferenceNodeLayout::num_nodes(ElementType::Wedge18), 18u);
        // Corner nodes should match Wedge6 vertices
        auto v0 = ReferenceNodeLayout::get_node_coords(ElementType::Wedge18, 0);
        auto v1 = ReferenceNodeLayout::get_node_coords(ElementType::Wedge18, 1);
        auto v2 = ReferenceNodeLayout::get_node_coords(ElementType::Wedge18, 2);
        EXPECT_NEAR(v0[0], Real(0), 1e-14);
        EXPECT_NEAR(v0[1], Real(0), 1e-14);
        EXPECT_NEAR(v0[2], Real(-1), 1e-14);
        EXPECT_NEAR(v1[0], Real(1), 1e-14);
        EXPECT_NEAR(v1[1], Real(0), 1e-14);
        EXPECT_NEAR(v1[2], Real(-1), 1e-14);
        EXPECT_NEAR(v2[0], Real(0), 1e-14);
        EXPECT_NEAR(v2[1], Real(1), 1e-14);
        EXPECT_NEAR(v2[2], Real(-1), 1e-14);
    }

    {
        LagrangeBasis pyr(ElementType::Pyramid5, 1);
        svmp::FE::math::Vector<Real, 3> xi{Real(0.1), Real(-0.2), Real(0.4)};
        std::vector<Real> vals;
        pyr.evaluate_values(xi, vals);
        const double sum = std::accumulate(vals.begin(), vals.end(), 0.0);
        EXPECT_NEAR(sum, 1.0, 1e-12);
    }
}

TEST(LagrangeBasis, NonTensorStridedEvaluationMatchesPointwise) {
    const std::vector<std::pair<ElementType, int>> cases = {
        {ElementType::Triangle3, 3},
        {ElementType::Tetra4, 3},
        {ElementType::Wedge6, 3},
        {ElementType::Pyramid5, 3},
    };
    const std::vector<StridedOutputRequest> requests = {
        {true, false, false},
        {false, true, false},
        {false, false, true},
        {true, true, false},
        {true, false, true},
        {false, true, true},
        {true, true, true},
    };

    for (const auto& [type, order] : cases) {
        for (const auto& request : requests) {
            SCOPED_TRACE(static_cast<int>(type));
            SCOPED_TRACE(order);
            SCOPED_TRACE(request.values ? "values" : "no values");
            SCOPED_TRACE(request.gradients ? "gradients" : "no gradients");
            SCOPED_TRACE(request.hessians ? "hessians" : "no hessians");
            expect_strided_matches_pointwise(type, order, request);
        }
    }
}

TEST(LagrangeBasis, RawOutputSinksMatchVectorEvaluationAcrossTopologies) {
    const std::vector<std::pair<ElementType, int>> cases = {
        {ElementType::Line2, 4},
        {ElementType::Quad4, 3},
        {ElementType::Hex8, 3},
        {ElementType::Triangle3, 4},
        {ElementType::Tetra4, 3},
        {ElementType::Wedge6, 3},
        {ElementType::Pyramid5, 3},
    };

    for (const auto& [type, order] : cases) {
        SCOPED_TRACE(static_cast<int>(type));
        SCOPED_TRACE(order);
        expect_raw_to_matches_vector_evaluation(type, order);
    }
}

TEST(LagrangeBasis, CanonicalConstructorsSupportArbitraryOrders) {
    const struct Case {
        ElementType type;
        int max_order;
    } cases[] = {
        {ElementType::Line2, 8},
        {ElementType::Triangle3, 6},
        {ElementType::Quad4, 6},
        {ElementType::Tetra4, 5},
        {ElementType::Hex8, 5},
        {ElementType::Wedge6, 5},
        {ElementType::Pyramid5, 5},
    };

    for (const auto& c : cases) {
        for (int order = 0; order <= c.max_order; ++order) {
            LagrangeBasis basis(c.type, order);
            EXPECT_EQ(basis.element_type(), c.type);
            EXPECT_EQ(basis.order(), order);
            EXPECT_EQ(basis.dimension(), expected_dimension(c.type));
            EXPECT_EQ(basis.size(), expected_lagrange_size(c.type, order));
            EXPECT_EQ(basis.nodes().size(), basis.size());
        }
    }
}

TEST(LagrangeBasis, AliasVariantsNormalizeToCanonicalPaths) {
    expect_alias_matches_canonical(
        ElementType::Line3, ElementType::Line2, 2, sample_points_for(ElementType::Line2));
    expect_alias_matches_canonical(
        ElementType::Triangle6, ElementType::Triangle3, 2, sample_points_for(ElementType::Triangle3));
    expect_alias_matches_canonical(
        ElementType::Quad9, ElementType::Quad4, 2, sample_points_for(ElementType::Quad4));
    expect_alias_matches_canonical(
        ElementType::Tetra10, ElementType::Tetra4, 2, sample_points_for(ElementType::Tetra4));
    expect_alias_matches_canonical(
        ElementType::Hex27, ElementType::Hex8, 2, sample_points_for(ElementType::Hex8));
    expect_alias_matches_canonical(
        ElementType::Wedge18, ElementType::Wedge6, 2, sample_points_for(ElementType::Wedge6));
    expect_alias_matches_canonical(
        ElementType::Pyramid14, ElementType::Pyramid5, 2, sample_points_for(ElementType::Pyramid5),
        Real(2e-10));
}

TEST(LagrangeBasis, SerendipityVariantsRemainRejected) {
    EXPECT_THROW((void)LagrangeBasis(ElementType::Quad8, 2), svmp::FE::FEException);
    EXPECT_THROW((void)LagrangeBasis(ElementType::Hex20, 2), svmp::FE::FEException);
    EXPECT_THROW((void)LagrangeBasis(ElementType::Wedge15, 2), svmp::FE::FEException);
    EXPECT_THROW((void)LagrangeBasis(ElementType::Pyramid13, 2), svmp::FE::FEException);
}

TEST(LagrangeBasis, GeneratedNodeOrderingIsDeterministicAcrossOrders) {
    const struct Case {
        ElementType type;
        int max_order;
    } cases[] = {
        {ElementType::Line2, 8},
        {ElementType::Triangle3, 6},
        {ElementType::Quad4, 6},
        {ElementType::Tetra4, 5},
        {ElementType::Hex8, 5},
        {ElementType::Wedge6, 5},
        {ElementType::Pyramid5, 5},
    };

    for (const auto& c : cases) {
        for (int order = 0; order <= c.max_order; ++order) {
            const auto generated_a = ReferenceNodeLayout::get_lagrange_node_coords(c.type, order);
            const auto generated_b = ReferenceNodeLayout::get_lagrange_node_coords(c.type, order);
            ASSERT_EQ(generated_a.size(), expected_lagrange_size(c.type, order));
            ASSERT_EQ(generated_a.size(), generated_b.size());
            for (std::size_t i = 0; i < generated_a.size(); ++i) {
                EXPECT_TRUE(points_close(generated_a[i], generated_b[i]));
            }
        }
    }
}

TEST(LagrangeBasis, NodeOrderingMatchesReferenceCoordinateOracles) {
    const std::array<ElementType, 18> cases = {
        ElementType::Line2, ElementType::Line3,
        ElementType::Triangle3, ElementType::Triangle6,
        ElementType::Quad4, ElementType::Quad8, ElementType::Quad9,
        ElementType::Tetra4, ElementType::Tetra10,
        ElementType::Hex8, ElementType::Hex20, ElementType::Hex27,
        ElementType::Wedge6, ElementType::Wedge15, ElementType::Wedge18,
        ElementType::Pyramid5, ElementType::Pyramid13, ElementType::Pyramid14,
    };

    for (ElementType type : cases) {
        const auto expected = reference_node_coords(type);
        ASSERT_FALSE(expected.empty());
        ASSERT_EQ(ReferenceNodeLayout::num_nodes(type), expected.size());
        for (std::size_t i = 0; i < expected.size(); ++i) {
            const auto actual = ReferenceNodeLayout::get_node_coords(type, i);
            EXPECT_TRUE(points_close(actual, expected[i]))
                << "Element type " << static_cast<int>(type)
                << ", node " << i;
        }
    }
}

TEST(LagrangeBasis, GeneratedLowOrderOrderingMatchesPublicAliasPaths) {
    const struct Case {
        ElementType type;
        int order;
        ElementType public_alias;
    } cases[] = {
        {ElementType::Line2, 1, ElementType::Line2},
        {ElementType::Line2, 2, ElementType::Line3},
        {ElementType::Triangle3, 1, ElementType::Triangle3},
        {ElementType::Triangle3, 2, ElementType::Triangle6},
        {ElementType::Quad4, 1, ElementType::Quad4},
        {ElementType::Quad4, 2, ElementType::Quad9},
        {ElementType::Tetra4, 1, ElementType::Tetra4},
        {ElementType::Tetra4, 2, ElementType::Tetra10},
        {ElementType::Hex8, 1, ElementType::Hex8},
        {ElementType::Hex8, 2, ElementType::Hex27},
        {ElementType::Wedge6, 1, ElementType::Wedge6},
        {ElementType::Wedge6, 2, ElementType::Wedge18},
        {ElementType::Pyramid5, 1, ElementType::Pyramid5},
        {ElementType::Pyramid5, 2, ElementType::Pyramid14},
    };

    for (const auto& c : cases) {
        const auto generated = ReferenceNodeLayout::get_lagrange_node_coords(c.type, c.order);
        ASSERT_EQ(generated.size(), ReferenceNodeLayout::num_nodes(c.public_alias));
        for (std::size_t i = 0; i < generated.size(); ++i) {
            const auto public_alias = ReferenceNodeLayout::get_node_coords(c.public_alias, i);
            EXPECT_TRUE(points_close(generated[i], public_alias));
        }
    }
}

TEST(LagrangeBasis, KroneckerDeltaAcrossCanonicalTopologiesAndOrders) {
    const struct Case {
        ElementType type;
        int max_order;
    } cases[] = {
        {ElementType::Line2, 8},
        {ElementType::Triangle3, 6},
        {ElementType::Quad4, 6},
        {ElementType::Tetra4, 5},
        {ElementType::Hex8, 5},
        {ElementType::Wedge6, 5},
        {ElementType::Pyramid5, 5},
    };

    for (const auto& c : cases) {
        for (int order = 0; order <= c.max_order; ++order) {
            LagrangeBasis basis(c.type, order);
            ASSERT_EQ(basis.size(), expected_lagrange_size(c.type, order));

            std::vector<Real> values;
            for (std::size_t node_i = 0; node_i < basis.size(); ++node_i) {
                basis.evaluate_values(basis.nodes()[node_i], values);
                ASSERT_EQ(values.size(), basis.size());
                for (std::size_t basis_i = 0; basis_i < basis.size(); ++basis_i) {
                    EXPECT_NEAR(values[basis_i], basis_i == node_i ? Real(1) : Real(0), Real(2e-10))
                        << "Element type " << static_cast<int>(c.type)
                        << ", order " << order
                        << ", node " << node_i
                        << ", basis " << basis_i;
                }
            }
        }
    }
}

TEST(LagrangeBasis, PartitionGradientAndHessianSumsAcrossCanonicalTopologiesAndOrders) {
    const struct Case {
        ElementType type;
        int max_order;
        Real tol;
    } cases[] = {
        {ElementType::Line2, 8, Real(1e-11)},
        {ElementType::Triangle3, 6, Real(1e-10)},
        {ElementType::Quad4, 6, Real(1e-10)},
        {ElementType::Tetra4, 5, Real(2e-10)},
        {ElementType::Hex8, 5, Real(2e-10)},
        {ElementType::Wedge6, 5, Real(5e-10)},
        {ElementType::Pyramid5, 5, Real(5e-7)},
    };

    for (const auto& c : cases) {
        for (int order = 0; order <= c.max_order; ++order) {
            LagrangeBasis basis(c.type, order);
            expect_partition_gradient_hessian_sums(basis, dense_sample_points_for(c.type), c.tol, c.tol);
        }
    }
}

TEST(LagrangeBasis, SimplexAxisScratchDynamicFallbackForHighOrder) {
    const struct Case {
        ElementType type;
        int order;
        Point point;
        Real tolerance;
    } cases[] = {
        {ElementType::Triangle3, 13, Point{Real(0.19), Real(0.31), Real(0)}, Real(1e-8)},
        {ElementType::Tetra4, 13, Point{Real(0.13), Real(0.17), Real(0.19)}, Real(1e-7)},
    };

    for (const auto& c : cases) {
        LagrangeBasis basis(c.type, c.order);
        std::vector<Real> values;
        std::vector<Gradient> gradients;
        std::vector<Hessian> hessians;
        basis.evaluate_all(c.point, values, gradients, hessians);

        ASSERT_EQ(values.size(), basis.size());
        ASSERT_EQ(gradients.size(), basis.size());
        ASSERT_EQ(hessians.size(), basis.size());

        Real value_sum = Real(0);
        Gradient gradient_sum{};
        Hessian hessian_sum{};
        for (std::size_t i = 0; i < basis.size(); ++i) {
            value_sum += values[i];
            for (std::size_t d = 0; d < 3u; ++d) {
                gradient_sum[d] += gradients[i][d];
                for (std::size_t e = 0; e < 3u; ++e) {
                    hessian_sum(d, e) += hessians[i](d, e);
                }
            }
        }

        EXPECT_NEAR(value_sum, Real(1), c.tolerance);
        for (std::size_t d = 0; d < 3u; ++d) {
            EXPECT_NEAR(gradient_sum[d], Real(0), c.tolerance);
            for (std::size_t e = 0; e < 3u; ++e) {
                EXPECT_NEAR(hessian_sum(d, e), Real(0), Real(10) * c.tolerance);
            }
        }
    }
}

TEST(LagrangeBasis, HighOrderAxisNearNodeMaintainsPartitionAndDerivativeSums) {
    const int order = 16;
    const LagrangeBasis basis(ElementType::Line2, order);
    const Real node = Real(-1) + Real(2 * 5) / static_cast<Real>(order);
    const Point point{node + Real(1e-7), Real(0), Real(0)};

    std::vector<Real> values;
    std::vector<Gradient> gradients;
    std::vector<Hessian> hessians;
    basis.evaluate_all(point, values, gradients, hessians);
    ASSERT_EQ(values.size(), basis.size());

    Real value_sum = Real(0);
    Real gradient_sum = Real(0);
    Real hessian_sum = Real(0);
    for (std::size_t i = 0; i < basis.size(); ++i) {
        value_sum += values[i];
        gradient_sum += gradients[i][0];
        hessian_sum += hessians[i](0, 0);
    }

    EXPECT_NEAR(value_sum, Real(1), Real(1e-12));
    EXPECT_NEAR(gradient_sum, Real(0), Real(1e-8));
    EXPECT_NEAR(hessian_sum, Real(0), Real(1e-5));
}

TEST(LagrangeBasis, PyramidFaceTracesMatchLowerDimensionalLagrangeBases) {
    const PyramidFace faces[] = {
        PyramidFace::Base,
        PyramidFace::South,
        PyramidFace::East,
        PyramidFace::North,
        PyramidFace::West,
    };

    for (int order = 1; order <= 5; ++order) {
        for (const auto face : faces) {
            expect_pyramid_face_trace_matches_lower_basis(
                order, face, face == PyramidFace::Base ? Real(2e-10) : Real(5e-10));
        }
    }
}

TEST(LagrangeBasis, PyramidEdgeTracesMatchLineLagrangeBasis) {
    const PyramidEdge edges[] = {
        PyramidEdge::BaseSouth,
        PyramidEdge::BaseEast,
        PyramidEdge::BaseNorth,
        PyramidEdge::BaseWest,
        PyramidEdge::VerticalSW,
        PyramidEdge::VerticalSE,
        PyramidEdge::VerticalNE,
        PyramidEdge::VerticalNW,
    };

    for (int order = 1; order <= 5; ++order) {
        for (const auto edge : edges) {
            expect_pyramid_edge_trace_matches_line_basis(order, edge, Real(5e-10));
        }
    }
}

TEST(LagrangeBasis, Pyramid14RationalNodalAndPartition) {
    using svmp::FE::basis::ReferenceNodeLayout;

    LagrangeBasis basis(ElementType::Pyramid14, 2);
    EXPECT_EQ(basis.dimension(), 3);
    EXPECT_EQ(basis.size(), 14u);

    // Kronecker nodal property at all Pyramid14 nodes
    for (std::size_t i = 0; i < basis.size(); ++i) {
        auto xi = ReferenceNodeLayout::get_node_coords(ElementType::Pyramid14, i);
        std::vector<Real> vals;
        basis.evaluate_values(xi, vals);
        ASSERT_EQ(vals.size(), basis.size());
        for (std::size_t j = 0; j < basis.size(); ++j) {
            const double expected = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(vals[j], expected, 1e-12);
        }
    }

    // Partition of unity at an interior point
    svmp::FE::math::Vector<Real, 3> xi{Real(0.1), Real(-0.2), Real(0.3)};
    std::vector<Real> vals;
    basis.evaluate_values(xi, vals);
    const double sum = std::accumulate(vals.begin(), vals.end(), 0.0);
    EXPECT_NEAR(sum, 1.0, 1e-12);
}

TEST(LagrangeBasis, Pyramid14GradientSumZero) {
    LagrangeBasis basis(ElementType::Pyramid14, 2);
    svmp::FE::math::Vector<Real, 3> xi{Real(0.15), Real(-0.1), Real(0.3)};

    std::vector<Gradient> grads;
    basis.evaluate_gradients(xi, grads);
    ASSERT_EQ(grads.size(), basis.size());

    Gradient sum{};
    for (const auto& g : grads) {
        sum[0] += g[0];
        sum[1] += g[1];
        sum[2] += g[2];
    }
    EXPECT_NEAR(sum[0], 0.0, 1e-8);
    EXPECT_NEAR(sum[1], 0.0, 1e-8);
    EXPECT_NEAR(sum[2], 0.0, 1e-8);
}

TEST(LagrangeBasis, HigherOrderP4KroneckerAndPartition) {
    struct Case {
        ElementType type;
        int order;
        svmp::FE::math::Vector<Real, 3> xi;
    };

    const std::vector<Case> cases = {
        {ElementType::Line2, 4, {Real(0.11), Real(0), Real(0)}},
        {ElementType::Quad4, 4, {Real(0.2), Real(-0.3), Real(0)}},
        {ElementType::Triangle3, 4, {Real(0.2), Real(0.3), Real(0)}},
        {ElementType::Hex8, 4, {Real(0.2), Real(-0.3), Real(0.4)}},
    };

    for (const auto& c : cases) {
        LagrangeBasis basis(c.type, c.order);

        // Partition of unity at an interior point
        std::vector<Real> values;
        basis.evaluate_values(c.xi, values);
        const double sum = std::accumulate(values.begin(), values.end(), 0.0);
        EXPECT_NEAR(sum, 1.0, 1e-12);

        // Kronecker delta property at all nodes
        const auto& nodes = basis.nodes();
        ASSERT_EQ(nodes.size(), basis.size());
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            basis.evaluate_values(nodes[i], values);
            ASSERT_EQ(values.size(), nodes.size());
            for (std::size_t j = 0; j < nodes.size(); ++j) {
                const double expected = (i == j) ? 1.0 : 0.0;
                EXPECT_NEAR(values[j], expected, 1e-12);
            }
        }
    }
}

TEST(LagrangeBasis, Pyramid14InterpolatesQuadraticPolynomials) {
    using svmp::FE::basis::ReferenceNodeLayout;

    LagrangeBasis basis(ElementType::Pyramid14, 2);
    const std::size_t n = basis.size();

    // Precompute nodal coordinates
    std::vector<svmp::FE::math::Vector<Real,3>> nodes;
    nodes.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        nodes.push_back(ReferenceNodeLayout::get_node_coords(ElementType::Pyramid14, i));
    }

    auto interpolate_and_check = [&](auto f, Real tol) {
        // Nodal coefficients
        std::vector<Real> coeffs(n);
        for (std::size_t i = 0; i < n; ++i) {
            const auto& x = nodes[i];
            coeffs[i] = f(x[0], x[1], x[2]);
        }

        // Test at a few interior points
        const svmp::FE::math::Vector<Real,3> test_pts[] = {
            {Real(0.1), Real(-0.2), Real(0.2)},
            {Real(-0.2), Real(0.15), Real(0.4)},
            {Real(0.05), Real(0.05), Real(0.3)}
        };

        for (const auto& xi : test_pts) {
            std::vector<Real> vals;
            basis.evaluate_values(xi, vals);
            ASSERT_EQ(vals.size(), n);

            Real u_interp = Real(0);
            for (std::size_t i = 0; i < n; ++i) {
                u_interp += coeffs[i] * vals[i];
            }

            const Real u_exact = f(xi[0], xi[1], xi[2]);
            EXPECT_NEAR(u_interp, u_exact, tol);
        }
    };

    // Constant, linear and quadratic monomials
    interpolate_and_check([](Real, Real, Real) { return Real(1); }, Real(1e-12));
    interpolate_and_check([](Real x, Real, Real) { return x; }, Real(1e-11));
    interpolate_and_check([](Real, Real y, Real) { return y; }, Real(1e-11));
    interpolate_and_check([](Real, Real, Real z) { return z; }, Real(1e-11));
    interpolate_and_check([](Real x, Real y, Real) { return x * y; }, Real(1e-10));
    interpolate_and_check([](Real x, Real, Real z) { return x * z; }, Real(1e-10));
    interpolate_and_check([](Real, Real y, Real z) { return y * z; }, Real(1e-10));
    interpolate_and_check([](Real x, Real, Real) { return x * x; }, Real(1e-10));
    interpolate_and_check([](Real, Real y, Real) { return y * y; }, Real(1e-10));
    interpolate_and_check([](Real, Real, Real z) { return z * z; }, Real(1e-10));
}

TEST(LagrangeBasis, Pyramid14GradientMatchesLinearFunctionGradient) {
    using svmp::FE::basis::ReferenceNodeLayout;

    LagrangeBasis basis(ElementType::Pyramid14, 2);
    const std::size_t n = basis.size();

    // Nodal coordinates and coefficients for f(x,y,z) = ax + by + cz
    const Real a = Real(1.2);
    const Real b = Real(-0.7);
    const Real c = Real(0.5);

    std::vector<Real> coeffs(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto x = ReferenceNodeLayout::get_node_coords(ElementType::Pyramid14, i);
        coeffs[i] = a * x[0] + b * x[1] + c * x[2];
    }

    const svmp::FE::math::Vector<Real,3> xi{Real(0.1), Real(-0.15), Real(0.35)};

    std::vector<Gradient> grads;
    basis.evaluate_gradients(xi, grads);
    ASSERT_EQ(grads.size(), n);

    Gradient g_interp{};
    for (std::size_t i = 0; i < n; ++i) {
        g_interp[0] += coeffs[i] * grads[i][0];
        g_interp[1] += coeffs[i] * grads[i][1];
        g_interp[2] += coeffs[i] * grads[i][2];
    }

    EXPECT_NEAR(g_interp[0], a, 1e-6);
    EXPECT_NEAR(g_interp[1], b, 1e-6);
    EXPECT_NEAR(g_interp[2], c, 1e-6);
}

TEST(LagrangeBasis, PyramidApexValuesRemainExactAcrossRepresentativeOrders) {
    const struct Case {
        ElementType type;
        int order;
    } cases[] = {
        {ElementType::Pyramid5, 1},
        {ElementType::Pyramid14, 2},
        {ElementType::Pyramid5, 4},
    };

    const svmp::FE::math::Vector<Real, 3> apex{Real(0), Real(0), Real(1)};
    for (const auto& c : cases) {
        LagrangeBasis basis(c.type, c.order);
        std::vector<Real> values;
        basis.evaluate_values(apex, values);
        ASSERT_EQ(values.size(), basis.size());

        const auto& nodes = basis.nodes();
        auto apex_it = std::find_if(
            nodes.begin(), nodes.end(),
            [](const auto& node) {
                return std::abs(node[0]) <= Real(1e-14) &&
                       std::abs(node[1]) <= Real(1e-14) &&
                       std::abs(node[2] - Real(1)) <= Real(1e-14);
            });
        ASSERT_NE(apex_it, nodes.end());
        const std::size_t apex_index = static_cast<std::size_t>(
            std::distance(nodes.begin(), apex_it));

        Real sum = Real(0);
        for (std::size_t i = 0; i < values.size(); ++i) {
            EXPECT_TRUE(std::isfinite(static_cast<double>(values[i])));
            sum += values[i];
            const Real expected = (i == apex_index) ? Real(1) : Real(0);
            EXPECT_NEAR(values[i], expected, 1e-12)
                << "order " << c.order << ", basis " << i;
        }
        EXPECT_NEAR(sum, Real(1), 1e-12);
    }
}

TEST(LagrangeBasis, PyramidGradientAtExactApexThrowsWhenLimitIsNotUnique) {
    const struct Case {
        ElementType type;
        int order;
    } cases[] = {
        {ElementType::Pyramid5, 1},
        {ElementType::Pyramid14, 2},
        {ElementType::Pyramid5, 4},
    };

    const svmp::FE::math::Vector<Real, 3> apex{Real(0), Real(0), Real(1)};
    for (const auto& c : cases) {
        LagrangeBasis basis(c.type, c.order);
        std::vector<Gradient> gradients;
        EXPECT_THROW(basis.evaluate_gradients(apex, gradients), svmp::FE::basis::BasisEvaluationException)
            << "order " << c.order;
    }
}

TEST(LagrangeBasis, PyramidApexValuesMatchDirectionalNearApexLimits) {
    const struct Case {
        ElementType type;
        int order;
        Real tol;
    } cases[] = {
        {ElementType::Pyramid5, 1, Real(3e-6)},
        {ElementType::Pyramid14, 2, Real(4e-6)},
        {ElementType::Pyramid5, 4, Real(1e-5)},
    };

    const std::array<std::array<Real, 2>, 4> directions = {{
        {Real(0), Real(0)},
        {Real(0.35), Real(-0.25)},
        {Real(-0.50), Real(0.45)},
        {Real(0.20), Real(0.60)},
    }};
    const Real t = Real(1e-6);
    const svmp::FE::math::Vector<Real, 3> apex{Real(0), Real(0), Real(1)};

    for (const auto& c : cases) {
        LagrangeBasis basis(c.type, c.order);
        std::vector<Real> apex_values;
        basis.evaluate_values(apex, apex_values);

        for (const auto& direction : directions) {
            const svmp::FE::math::Vector<Real, 3> xi{
                t * direction[0],
                t * direction[1],
                Real(1) - t
            };

            std::vector<Real> values;
            basis.evaluate_values(xi, values);
            ASSERT_EQ(values.size(), apex_values.size());

            for (std::size_t i = 0; i < values.size(); ++i) {
                EXPECT_NEAR(values[i], apex_values[i], c.tol)
                    << "order " << c.order
                    << ", basis " << i
                    << ", direction (" << direction[0] << ", " << direction[1] << ")";
            }
        }
    }
}

TEST(LagrangeBasis, PyramidNearApexGradientShowsDirectionalSpread) {
    const struct Case {
        ElementType type;
        int order;
        Real min_spread;
    } cases[] = {
        {ElementType::Pyramid5, 1, Real(5e-2)},
        {ElementType::Pyramid14, 2, Real(5e-2)},
    };

    const std::array<std::array<Real, 2>, 4> directions = {{
        {Real(0), Real(0)},
        {Real(0.45), Real(-0.30)},
        {Real(-0.35), Real(0.40)},
        {Real(0.25), Real(0.55)},
    }};
    const Real t = Real(1e-6);

    for (const auto& c : cases) {
        LagrangeBasis basis(c.type, c.order);
        double max_spread = 0.0;

        std::vector<std::vector<Gradient>> directional_gradients;
        directional_gradients.reserve(directions.size());
        for (const auto& direction : directions) {
            const svmp::FE::math::Vector<Real, 3> xi{
                t * direction[0],
                t * direction[1],
                Real(1) - t
            };

            std::vector<Gradient> gradients;
            basis.evaluate_gradients(xi, gradients);
            directional_gradients.push_back(std::move(gradients));
        }

        for (std::size_t i = 0; i < basis.size(); ++i) {
            for (int d = 0; d < 3; ++d) {
                double min_value = std::numeric_limits<double>::infinity();
                double max_value = -std::numeric_limits<double>::infinity();
                for (const auto& gradients : directional_gradients) {
                    const double value = static_cast<double>(gradients[i][static_cast<std::size_t>(d)]);
                    min_value = std::min(min_value, value);
                    max_value = std::max(max_value, value);
                }
                max_spread = std::max(max_spread, max_value - min_value);
            }
        }

        EXPECT_GT(max_spread, static_cast<double>(c.min_spread))
            << "order " << c.order;
    }
}

TEST(LagrangeBasis, GradientSumZeroQuadAndTet) {
    const std::vector<std::pair<ElementType, svmp::FE::math::Vector<Real, 3>>> cases = {
        {ElementType::Quad4, svmp::FE::math::Vector<Real, 3>{Real(0.2), Real(-0.1), Real(0)}},
        {ElementType::Tetra4, svmp::FE::math::Vector<Real, 3>{Real(0.1), Real(0.2), Real(0.1)}}
    };

    for (const auto& c : cases) {
        LagrangeBasis basis(c.first, 1);
        std::vector<Gradient> grads;
        basis.evaluate_gradients(c.second, grads);

        ASSERT_EQ(grads.size(), basis.size());
        Gradient sum{};
        for (const auto& g : grads) {
            sum[0] += g[0];
            sum[1] += g[1];
            sum[2] += g[2];
        }
        EXPECT_NEAR(sum[0], 0.0, 1e-12);
        EXPECT_NEAR(sum[1], 0.0, 1e-12);
        EXPECT_NEAR(sum[2], 0.0, 1e-12);
    }
}

TEST(LagrangeBasis, HexPartitionAndGradientSumZeroOrderThree) {
    LagrangeBasis basis(ElementType::Hex8, 3);
    svmp::FE::math::Vector<Real, 3> xi{Real(0.1), Real(-0.2), Real(0.25)};

    std::vector<Real> values;
    basis.evaluate_values(xi, values);
    const double sum = std::accumulate(values.begin(), values.end(), 0.0);
    EXPECT_NEAR(sum, 1.0, 1e-12);

    std::vector<Gradient> grads;
    basis.evaluate_gradients(xi, grads);
    Gradient gsum{};
    for (const auto& g : grads) {
        gsum[0] += g[0];
        gsum[1] += g[1];
        gsum[2] += g[2];
    }
    EXPECT_NEAR(gsum[0], 0.0, 1e-10);
    EXPECT_NEAR(gsum[1], 0.0, 1e-10);
    EXPECT_NEAR(gsum[2], 0.0, 1e-10);
}

TEST(LagrangeBasis, OracleLine3ValuesGradientsAndHessians) {
    LagrangeBasis basis(ElementType::Line3, 2);
    const Point xi{Real(0.2), Real(0), Real(0)};

    std::vector<Real> values;
    std::vector<Gradient> gradients;
    std::vector<Hessian> hessians;
    basis.evaluate_values(xi, values);
    basis.evaluate_gradients(xi, gradients);
    basis.evaluate_hessians(xi, hessians);

    ASSERT_EQ(values.size(), 3u);
    ASSERT_EQ(gradients.size(), 3u);
    ASSERT_EQ(hessians.size(), 3u);

    const Real expected_values[] = {Real(-2) / Real(25), Real(3) / Real(25), Real(24) / Real(25)};
    const Real expected_gradients[] = {Real(-3) / Real(10), Real(7) / Real(10), Real(-2) / Real(5)};
    const Real expected_hessians[] = {Real(1), Real(1), Real(-2)};

    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(values[i], expected_values[i], 1e-14);
        EXPECT_NEAR(gradients[i][0], expected_gradients[i], 1e-14);
        EXPECT_NEAR(hessians[i](0, 0), expected_hessians[i], 1e-14);
    }
}

TEST(LagrangeBasis, OracleTriangle3ValuesGradientsAndHessians) {
    LagrangeBasis basis(ElementType::Triangle3, 1);
    const Point xi{Real(0.2), Real(0.3), Real(0)};

    std::vector<Real> values;
    std::vector<Gradient> gradients;
    std::vector<Hessian> hessians;
    basis.evaluate_values(xi, values);
    basis.evaluate_gradients(xi, gradients);
    basis.evaluate_hessians(xi, hessians);

    ASSERT_EQ(values.size(), 3u);
    const Point expected_gradients[] = {
        Point{Real(-1), Real(-1), Real(0)},
        Point{Real(1), Real(0), Real(0)},
        Point{Real(0), Real(1), Real(0)}
    };
    const Real expected_values[] = {Real(0.5), Real(0.2), Real(0.3)};

    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(values[i], expected_values[i], 1e-14);
        EXPECT_NEAR(gradients[i][0], expected_gradients[i][0], 1e-14);
        EXPECT_NEAR(gradients[i][1], expected_gradients[i][1], 1e-14);
        for (int a = 0; a < 2; ++a) {
            for (int b = 0; b < 2; ++b) {
                EXPECT_NEAR(hessians[i](static_cast<std::size_t>(a), static_cast<std::size_t>(b)),
                            Real(0), 1e-14);
            }
        }
    }
}

TEST(LagrangeBasis, OracleQuad4ValuesGradientsAndHessians) {
    LagrangeBasis basis(ElementType::Quad4, 1);
    const Point xi{Real(0.2), Real(-0.4), Real(0)};

    std::vector<Real> values;
    std::vector<Gradient> gradients;
    std::vector<Hessian> hessians;
    basis.evaluate_values(xi, values);
    basis.evaluate_gradients(xi, gradients);
    basis.evaluate_hessians(xi, hessians);

    ASSERT_EQ(values.size(), 4u);
    const Real expected_values[] = {Real(7) / Real(25), Real(21) / Real(50),
                                    Real(9) / Real(50), Real(3) / Real(25)};
    const Point expected_gradients[] = {
        Point{Real(-7) / Real(20), Real(-1) / Real(5), Real(0)},
        Point{Real(7) / Real(20), Real(-3) / Real(10), Real(0)},
        Point{Real(3) / Real(20), Real(3) / Real(10), Real(0)},
        Point{Real(-3) / Real(20), Real(1) / Real(5), Real(0)}
    };
    const Real expected_hxy[] = {Real(1) / Real(4), Real(-1) / Real(4),
                                 Real(1) / Real(4), Real(-1) / Real(4)};

    for (std::size_t i = 0; i < 4; ++i) {
        EXPECT_NEAR(values[i], expected_values[i], 1e-14);
        EXPECT_NEAR(gradients[i][0], expected_gradients[i][0], 1e-14);
        EXPECT_NEAR(gradients[i][1], expected_gradients[i][1], 1e-14);
        EXPECT_NEAR(hessians[i](0, 0), Real(0), 1e-14);
        EXPECT_NEAR(hessians[i](1, 1), Real(0), 1e-14);
        EXPECT_NEAR(hessians[i](0, 1), expected_hxy[i], 1e-14);
        EXPECT_NEAR(hessians[i](1, 0), expected_hxy[i], 1e-14);
    }
}

TEST(LagrangeBasis, OracleWedge6ValuesGradientsAndHessians) {
    LagrangeBasis basis(ElementType::Wedge6, 1);
    const Point xi{Real(0.2), Real(0.25), Real(-0.3)};

    std::vector<Real> values;
    std::vector<Gradient> gradients;
    std::vector<Hessian> hessians;
    basis.evaluate_values(xi, values);
    basis.evaluate_gradients(xi, gradients);
    basis.evaluate_hessians(xi, hessians);

    ASSERT_EQ(values.size(), 6u);
    const Real expected_values[] = {
        Real(143) / Real(400), Real(13) / Real(100), Real(13) / Real(80),
        Real(77) / Real(400), Real(7) / Real(100), Real(7) / Real(80)
    };
    const Point expected_gradients[] = {
        Point{Real(-13) / Real(20), Real(-13) / Real(20), Real(-11) / Real(40)},
        Point{Real(13) / Real(20), Real(0), Real(-1) / Real(10)},
        Point{Real(0), Real(13) / Real(20), Real(-1) / Real(8)},
        Point{Real(-7) / Real(20), Real(-7) / Real(20), Real(11) / Real(40)},
        Point{Real(7) / Real(20), Real(0), Real(1) / Real(10)},
        Point{Real(0), Real(7) / Real(20), Real(1) / Real(8)}
    };
    const Point expected_hxz[] = {
        Point{Real(1) / Real(2), Real(1) / Real(2), Real(0)},
        Point{Real(-1) / Real(2), Real(0), Real(0)},
        Point{Real(0), Real(-1) / Real(2), Real(0)},
        Point{Real(-1) / Real(2), Real(-1) / Real(2), Real(0)},
        Point{Real(1) / Real(2), Real(0), Real(0)},
        Point{Real(0), Real(1) / Real(2), Real(0)}
    };

    for (std::size_t i = 0; i < 6; ++i) {
        EXPECT_NEAR(values[i], expected_values[i], 1e-14);
        EXPECT_NEAR(gradients[i][0], expected_gradients[i][0], 1e-14);
        EXPECT_NEAR(gradients[i][1], expected_gradients[i][1], 1e-14);
        EXPECT_NEAR(gradients[i][2], expected_gradients[i][2], 1e-14);
        EXPECT_NEAR(hessians[i](0, 0), Real(0), 1e-14);
        EXPECT_NEAR(hessians[i](1, 1), Real(0), 1e-14);
        EXPECT_NEAR(hessians[i](2, 2), Real(0), 1e-14);
        EXPECT_NEAR(hessians[i](0, 1), Real(0), 1e-14);
        EXPECT_NEAR(hessians[i](1, 0), Real(0), 1e-14);
        EXPECT_NEAR(hessians[i](0, 2), expected_hxz[i][0], 1e-14);
        EXPECT_NEAR(hessians[i](2, 0), expected_hxz[i][0], 1e-14);
        EXPECT_NEAR(hessians[i](1, 2), expected_hxz[i][1], 1e-14);
        EXPECT_NEAR(hessians[i](2, 1), expected_hxz[i][1], 1e-14);
    }
}

TEST(LagrangeBasis, DeterministicBoundarySweepMaintainsPartitionAndFiniteDerivatives) {
    const std::vector<std::pair<ElementType, int>> cases = {
        {ElementType::Line2, 1},
        {ElementType::Line3, 2},
        {ElementType::Triangle3, 1},
        {ElementType::Triangle6, 2},
        {ElementType::Quad4, 1},
        {ElementType::Quad9, 2},
        {ElementType::Tetra4, 1},
        {ElementType::Tetra10, 2},
        {ElementType::Hex8, 1},
        {ElementType::Hex27, 2},
        {ElementType::Wedge6, 1},
        {ElementType::Wedge18, 2},
        {ElementType::Pyramid5, 1},
        {ElementType::Pyramid14, 2},
    };

    for (const auto& [type, order] : cases) {
        LagrangeBasis basis(type, order);
        for (const auto& xi : boundary_stress_points_for(type)) {
            std::vector<Real> values;
            std::vector<Gradient> gradients;
            std::vector<Hessian> hessians;
            basis.evaluate_values(xi, values);
            basis.evaluate_gradients(xi, gradients);
            basis.evaluate_hessians(xi, hessians);

            ASSERT_EQ(values.size(), basis.size());
            ASSERT_EQ(gradients.size(), basis.size());
            ASSERT_EQ(hessians.size(), basis.size());

            Real sum = Real(0);
            for (Real value : values) {
                EXPECT_TRUE(std::isfinite(value));
                sum += value;
            }
            expect_all_finite(gradients);
            expect_hessians_finite(hessians, basis.dimension());
            EXPECT_NEAR(sum, Real(1), type == ElementType::Pyramid5 || type == ElementType::Pyramid14
                                       ? Real(1e-8)
                                       : Real(1e-12))
                << "type=" << static_cast<int>(type)
                << ", order=" << order
                << ", xi=(" << xi[0] << "," << xi[1] << "," << xi[2] << ")";
        }
    }
}

TEST(LagrangeBasis, FiniteDifferenceGradientsAcrossSupportedLinearShapes) {
    const std::vector<LagrangeAccuracyCase> cases = {
        {ElementType::Line2, 1, sample_points_for(ElementType::Line2)},
        {ElementType::Triangle3, 1, sample_points_for(ElementType::Triangle3)},
        {ElementType::Quad4, 1, sample_points_for(ElementType::Quad4)},
        {ElementType::Tetra4, 1, sample_points_for(ElementType::Tetra4)},
        {ElementType::Hex8, 1, sample_points_for(ElementType::Hex8)},
        {ElementType::Wedge6, 1, sample_points_for(ElementType::Wedge6)},
        {ElementType::Pyramid5, 1, sample_points_for(ElementType::Pyramid5)},
    };

    for (const auto& c : cases) {
        expect_gradients_match_finite_difference(c, Real(1e-6), Real(1e-6));
    }
}

TEST(LagrangeBasis, FiniteDifferenceGradientsAcrossSupportedQuadraticShapes) {
    const std::vector<LagrangeAccuracyCase> cases = {
        {ElementType::Line3, 2, sample_points_for(ElementType::Line3)},
        {ElementType::Triangle6, 2, sample_points_for(ElementType::Triangle6)},
        {ElementType::Quad9, 2, sample_points_for(ElementType::Quad9)},
        {ElementType::Tetra10, 2, sample_points_for(ElementType::Tetra10)},
        {ElementType::Hex27, 2, sample_points_for(ElementType::Hex27)},
        {ElementType::Wedge18, 2, sample_points_for(ElementType::Wedge18)},
        {ElementType::Pyramid14, 2, sample_points_for(ElementType::Pyramid14)},
    };

    for (const auto& c : cases) {
        expect_gradients_match_finite_difference(c, Real(1e-6), Real(2e-6));
    }
}

TEST(LagrangeBasis, LinearPolynomialReproductionAcrossSupportedLinearShapes) {
    const std::vector<LagrangeAccuracyCase> cases = {
        {ElementType::Line2, 1, sample_points_for(ElementType::Line2)},
        {ElementType::Triangle3, 1, sample_points_for(ElementType::Triangle3)},
        {ElementType::Quad4, 1, sample_points_for(ElementType::Quad4)},
        {ElementType::Tetra4, 1, sample_points_for(ElementType::Tetra4)},
        {ElementType::Hex8, 1, sample_points_for(ElementType::Hex8)},
        {ElementType::Wedge6, 1, sample_points_for(ElementType::Wedge6)},
        {ElementType::Pyramid5, 1, sample_points_for(ElementType::Pyramid5)},
    };

    const std::vector<std::array<int, 3>> exponents = {
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1},
    };

    for (const auto& c : cases) {
        const std::vector<std::array<int, 3>> relevant(
            exponents.begin(),
            exponents.begin() + static_cast<std::ptrdiff_t>(c.type == ElementType::Line2 ? 2 :
                                                            (c.type == ElementType::Triangle3 ||
                                                             c.type == ElementType::Quad4) ? 3 : 4));
        expect_polynomial_reproduction(c, relevant, Real(1e-12));
    }
}

TEST(LagrangeBasis, QuadraticPolynomialReproductionAcrossSupportedQuadraticShapes) {
    const std::vector<LagrangeAccuracyCase> cases = {
        {ElementType::Line3, 2, sample_points_for(ElementType::Line3)},
        {ElementType::Triangle6, 2, sample_points_for(ElementType::Triangle6)},
        {ElementType::Quad9, 2, sample_points_for(ElementType::Quad9)},
        {ElementType::Tetra10, 2, sample_points_for(ElementType::Tetra10)},
        {ElementType::Hex27, 2, sample_points_for(ElementType::Hex27)},
        {ElementType::Wedge18, 2, sample_points_for(ElementType::Wedge18)},
        {ElementType::Pyramid14, 2, sample_points_for(ElementType::Pyramid14)},
    };

    const std::vector<std::array<int, 3>> line_exponents = {
        {0, 0, 0}, {1, 0, 0}, {2, 0, 0}
    };
    const std::vector<std::array<int, 3>> surface_exponents = {
        {0, 0, 0}, {1, 0, 0}, {0, 1, 0},
        {2, 0, 0}, {1, 1, 0}, {0, 2, 0}
    };
    const std::vector<std::array<int, 3>> volume_exponents = {
        {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
        {2, 0, 0}, {1, 1, 0}, {0, 2, 0},
        {1, 0, 1}, {0, 1, 1}, {0, 0, 2}
    };

    for (const auto& c : cases) {
        if (c.type == ElementType::Line3) {
            expect_polynomial_reproduction(c, line_exponents, Real(1e-12));
        } else if (c.type == ElementType::Triangle6 || c.type == ElementType::Quad9) {
            expect_polynomial_reproduction(c, surface_exponents, Real(1e-11));
        } else {
            expect_polynomial_reproduction(c, volume_exponents, Real(2e-10));
        }
    }
}

TEST(LagrangeBasis, HighOrderTensorLagrangeMaintainsPartitionAndDerivativeSums) {
    const std::vector<LagrangeAccuracyCase> cases = {
        {ElementType::Line2, 8, {Point{-0.875, 0, 0}, Point{0.125, 0, 0}, Point{1, 0, 0}}},
        {ElementType::Quad4, 7, {Point{0.2, -0.35, 0}, Point{-1, 0.5, 0}, Point{0.5, 1, 0}}},
        {ElementType::Hex8, 6, {Point{0.1, -0.2, 0.3}, Point{-1, 0.5, 1}, Point{0.75, -1, -0.5}}},
    };

    for (const auto& c : cases) {
        LagrangeBasis basis(c.type, c.order);
        expect_partition_gradient_hessian_sums(basis, c.points, Real(2e-12), Real(2e-8));
    }
}

TEST(LagrangeBasis, HighOrderTensorLagrangeReproducesTensorPolynomials) {
    const LagrangeAccuracyCase line{ElementType::Line2,
                                    8,
                                    {Point{-0.73, 0, 0}, Point{-0.1, 0, 0}, Point{0.64, 0, 0}}};
    expect_polynomial_reproduction(line,
                                   {{0, 0, 0}, {1, 0, 0}, {4, 0, 0}, {8, 0, 0}},
                                   Real(1e-11));

    const LagrangeAccuracyCase quad{ElementType::Quad4,
                                    7,
                                    {Point{-0.6, -0.2, 0}, Point{0.15, 0.45, 0}, Point{0.8, -0.55, 0}}};
    expect_polynomial_reproduction(quad,
                                   {{0, 0, 0}, {7, 0, 0}, {0, 7, 0}, {4, 3, 0}},
                                   Real(5e-10));

    const LagrangeAccuracyCase hex{ElementType::Hex8,
                                   6,
                                   {Point{-0.4, 0.2, -0.3}, Point{0.35, -0.55, 0.25}, Point{0.75, 0.4, -0.65}}};
    expect_polynomial_reproduction(hex,
                                   {{0, 0, 0}, {6, 0, 0}, {0, 6, 0}, {0, 0, 6}, {3, 2, 4}},
                                   Real(2e-9));
}
