/**
 * @file test_BasisHessians.cpp
 * @brief Analytical Hessian coverage for the migrated Lagrange basis.
 */

#include <gtest/gtest.h>

#include "FE/Basis/BasisFactory.h"
#include "FE/Basis/LagrangeBasis.h"
#include "FE/Basis/SerendipityBasis.h"

#include <array>
#include <limits>
#include <vector>

using namespace svmp::FE;
using namespace svmp::FE::basis;

namespace {

// The exact Hessian identities below -- the partition sum (sum_i Hess N_i = 0) and
// symmetry (Hess N_i = Hess N_i^T) -- have a floating-point round-off residual at
// every order and family, so they share one tolerance. The finite-difference-vs-
// analytic comparisons keep their own, larger, order-dependent tolerances because
// finite-difference error grows with order.
constexpr double kHessianInvariantTol = double(1e-12);

void numerical_gradient_helper(const BasisFunction& basis,
                               const math::Vector<double, 3>& xi,
                               std::vector<Gradient>& gradients,
                               double eps = double(1e-6))
{
    std::vector<double> base;
    basis.evaluate_values(xi, base);
    gradients.assign(base.size(), Gradient::Zero());

    for (int d = 0; d < basis.dimension(); ++d) {
        const std::size_t sd = static_cast<std::size_t>(d);
        math::Vector<double, 3> xi_p = xi;
        math::Vector<double, 3> xi_m = xi;
        xi_p[sd] += eps;
        xi_m[sd] -= eps;

        std::vector<double> v_p;
        std::vector<double> v_m;
        basis.evaluate_values(xi_p, v_p);
        basis.evaluate_values(xi_m, v_m);

        for (std::size_t n = 0; n < base.size(); ++n) {
            gradients[n][sd] = (v_p[n] - v_m[n]) / (double(2) * eps);
        }
    }
}

void numerical_hessian_helper(const BasisFunction& basis,
                              const math::Vector<double, 3>& xi,
                              std::vector<Hessian>& hessians,
                              double eps = double(1e-5))
{
    hessians.assign(basis.size(), Hessian::Zero());
    const int dim = basis.dimension();

    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            math::Vector<double, 3> xi_p = xi;
            math::Vector<double, 3> xi_m = xi;
            const std::size_t sj = static_cast<std::size_t>(j);
            xi_p[sj] += eps;
            xi_m[sj] -= eps;

            std::vector<Gradient> g_p;
            std::vector<Gradient> g_m;
            basis.evaluate_gradients(xi_p, g_p);
            basis.evaluate_gradients(xi_m, g_m);

            for (std::size_t n = 0; n < basis.size(); ++n) {
                const std::size_t si = static_cast<std::size_t>(i);
                hessians[n](si, sj) = (g_p[n][si] - g_m[n][si]) / (double(2) * eps);
            }
        }
    }
}

std::vector<math::Vector<double, 3>> sample_points_for(BasisTopology topology) {
    switch (topology) {
        case BasisTopology::Line:
            return {{double(-0.35), double(0), double(0)}, {double(0.2), double(0), double(0)}};
        case BasisTopology::Triangle:
            return {{double(0.15), double(0.2), double(0)}, {double(0.25), double(0.1), double(0)}};
        case BasisTopology::Quadrilateral:
            return {{double(0.2), double(-0.3), double(0)}, {double(-0.45), double(0.25), double(0)}};
        case BasisTopology::Tetrahedron:
            return {{double(0.12), double(0.18), double(0.16)}, {double(0.2), double(0.1), double(0.18)}};
        case BasisTopology::Hexahedron:
            return {{double(0.1), double(-0.2), double(0.3)}, {double(-0.35), double(0.25), double(-0.15)}};
        case BasisTopology::Wedge:
            return {{double(0.18), double(0.22), double(-0.2)}, {double(0.12), double(0.16), double(0.1)}};
        default:
            return {{double(0), double(0), double(0)}};
    }
}

void expect_gradients_match_numerical(const BasisFunction& basis,
                                      const std::vector<math::Vector<double, 3>>& points,
                                      double tol,
                                      double eps = double(1e-6))
{
    for (const auto& xi : points) {
        std::vector<Gradient> analytical;
        std::vector<Gradient> numerical;
        basis.evaluate_gradients(xi, analytical);
        numerical_gradient_helper(basis, xi, numerical, eps);

        ASSERT_EQ(analytical.size(), numerical.size());
        for (std::size_t n = 0; n < analytical.size(); ++n) {
            for (int d = 0; d < basis.dimension(); ++d) {
                const std::size_t sd = static_cast<std::size_t>(d);
                EXPECT_NEAR(analytical[n][sd], numerical[n][sd], tol)
                    << "basis " << n << ", component " << d
                    << ", element " << static_cast<int>(named_element_for(basis.topology(), basis.order(), basis.basis_type()))
                    << ", order " << basis.order();
            }
        }
    }
}

void expect_hessians_match_numerical(const BasisFunction& basis,
                                     const std::vector<math::Vector<double, 3>>& points,
                                     double tol,
                                     double eps = double(1e-5))
{
    for (const auto& xi : points) {
        std::vector<Hessian> analytical;
        std::vector<Hessian> numerical;
        basis.evaluate_hessians(xi, analytical);
        numerical_hessian_helper(basis, xi, numerical, eps);

        ASSERT_EQ(analytical.size(), numerical.size());
        for (std::size_t n = 0; n < analytical.size(); ++n) {
            for (int i = 0; i < basis.dimension(); ++i) {
                for (int j = 0; j < basis.dimension(); ++j) {
                    const std::size_t si = static_cast<std::size_t>(i);
                    const std::size_t sj = static_cast<std::size_t>(j);
                    EXPECT_NEAR(analytical[n](si, sj), numerical[n](si, sj), tol)
                        << "basis " << n << ", component (" << i << "," << j
                        << "), element " << static_cast<int>(named_element_for(basis.topology(), basis.order(), basis.basis_type()))
                        << ", order " << basis.order();
                }
            }
        }
    }
}

void expect_partition_hessian_sum_zero(const BasisFunction& basis,
                                       const math::Vector<double, 3>& xi,
                                       double tol)
{
    std::vector<Hessian> hessians;
    basis.evaluate_hessians(xi, hessians);

    Hessian sum = Hessian::Zero();
    for (const auto& hessian : hessians) {
        for (std::size_t r = 0; r < 3u; ++r) {
            for (std::size_t c = 0; c < 3u; ++c) {
                sum(r, c) += hessian(r, c);
            }
        }
    }

    for (int r = 0; r < basis.dimension(); ++r) {
        for (int c = 0; c < basis.dimension(); ++c) {
            EXPECT_NEAR(sum(static_cast<std::size_t>(r), static_cast<std::size_t>(c)),
                        double(0),
                        tol)
                << "element " << static_cast<int>(named_element_for(basis.topology(), basis.order(), basis.basis_type()))
                << ", order " << basis.order();
        }
    }
}

void expect_hessians_symmetric(const BasisFunction& basis,
                               const math::Vector<double, 3>& xi,
                               double tol)
{
    std::vector<Hessian> hessians;
    basis.evaluate_hessians(xi, hessians);

    for (const auto& hessian : hessians) {
        for (int r = 0; r < basis.dimension(); ++r) {
            for (int c = r + 1; c < basis.dimension(); ++c) {
                const std::size_t sr = static_cast<std::size_t>(r);
                const std::size_t sc = static_cast<std::size_t>(c);
                EXPECT_NEAR(hessian(sr, sc), hessian(sc, sr), tol);
            }
        }
    }
}

void expect_inactive_z_derivatives_zero(const BasisFunction& basis,
                                        const std::vector<math::Vector<double, 3>>& points,
                                        double tol)
{
    ASSERT_EQ(basis.dimension(), 2);
    for (const auto& xi : points) {
        std::vector<Gradient> gradients;
        std::vector<Hessian> hessians;
        basis.evaluate_gradients(xi, gradients);
        basis.evaluate_hessians(xi, hessians);

        ASSERT_EQ(gradients.size(), basis.size());
        ASSERT_EQ(hessians.size(), basis.size());
        for (std::size_t n = 0; n < basis.size(); ++n) {
            EXPECT_NEAR(gradients[n][2], double(0), tol)
                << "basis " << n << ", element "
                << static_cast<int>(named_element_for(basis.topology(), basis.order(), basis.basis_type()))
                << ", order " << basis.order();
            for (std::size_t d = 0; d < 3u; ++d) {
                EXPECT_NEAR(hessians[n](2, d), double(0), tol)
                    << "basis " << n << ", component (2," << d
                    << "), element " << static_cast<int>(named_element_for(basis.topology(), basis.order(), basis.basis_type()))
                    << ", order " << basis.order();
                EXPECT_NEAR(hessians[n](d, 2), double(0), tol)
                    << "basis " << n << ", component (" << d
                    << ",2), element " << static_cast<int>(named_element_for(basis.topology(), basis.order(), basis.basis_type()))
                    << ", order " << basis.order();
            }
        }
    }
}

std::vector<math::Vector<double, 3>> serendipity_sample_points(BasisTopology topology) {
    if (topology == BasisTopology::Quadrilateral) {
        return {{double(0.17), double(-0.31), double(0)}, {double(-0.45), double(0.25), double(0)}};
    }
    if (topology == BasisTopology::Hexahedron) {
        return {{double(0.2), double(-0.1), double(0.3)}, {double(-0.35), double(0.25), double(-0.15)}};
    }
    return {{double(0.2), double(0.3), double(0.1)}, {double(0.12), double(0.16), double(-0.2)}};  // wedge
}

} // namespace

TEST(BasisHessians, LagrangeCanonicalTopologiesMatchNumericalHessians) {
    const struct Case {
        BasisTopology topology;
        int order;
        double tol;
        double eps;
    } cases[] = {
        {BasisTopology::Line, 3, double(1e-7), double(1e-5)},
        {BasisTopology::Triangle, 3, double(2e-6), double(1e-5)},
        {BasisTopology::Quadrilateral, 3, double(1e-6), double(1e-5)},
        {BasisTopology::Tetrahedron, 2, double(1e-6), double(1e-5)},
        {BasisTopology::Hexahedron, 2, double(1e-6), double(1e-5)},
        {BasisTopology::Wedge, 2, double(1e-5), double(1e-5)},
    };

    for (const auto& c : cases) {
        LagrangeBasis basis(c.topology, c.order);
        expect_hessians_match_numerical(basis, sample_points_for(c.topology), c.tol, c.eps);
    }
}

TEST(BasisHessians, LagrangeHessiansSumToZeroAndAreSymmetric) {
    const struct Case {
        BasisTopology topology;
        int order;
        math::Vector<double, 3> xi;
    } cases[] = {
        {BasisTopology::Line, 3, {double(0.15), double(0), double(0)}},
        {BasisTopology::Triangle, 3, {double(0.2), double(0.25), double(0)}},
        {BasisTopology::Quadrilateral, 3, {double(0.3), double(-0.2), double(0)}},
        {BasisTopology::Tetrahedron, 2, {double(0.15), double(0.2), double(0.1)}},
        {BasisTopology::Hexahedron, 2, {double(0.1), double(-0.2), double(0.3)}},
        {BasisTopology::Wedge, 2, {double(0.2), double(0.15), double(-0.3)}},
    };

    for (const auto& c : cases) {
        LagrangeBasis basis(c.topology, c.order);
        expect_partition_hessian_sum_zero(basis, c.xi, kHessianInvariantTol);
        expect_hessians_symmetric(basis, c.xi, kHessianInvariantTol);
    }
}

TEST(BasisHessians, SerendipityHessiansSumToZeroAndAreSymmetric) {
    const struct Case {
        ElementType type;
        int order;
        math::Vector<double, 3> xi;
    } cases[] = {
        {ElementType::Quad8, 2, {double(0.17), double(-0.31), double(0)}},
        {ElementType::Hex20, 2, {double(0.2), double(-0.1), double(0.3)}},
        {ElementType::Wedge15, 2, {double(0.2), double(0.3), double(0.1)}},
    };

    for (const auto& c : cases) {
        SerendipityBasis basis(c.type, c.order);
        expect_partition_hessian_sum_zero(basis, c.xi, kHessianInvariantTol);
        expect_hessians_symmetric(basis, c.xi, kHessianInvariantTol);
    }
}

TEST(BasisHessians, SolverMappedVolumeSelectionsSatisfyInvariants) {
    // Mirrors the full default element set pinned in
    // BasisFactoryDefaults.SelectionsArePinnedForAllSupportedElements (including
    // both wedge defaults: Wedge15 serendipity and Wedge18 Lagrange), so every
    // family the solver adapter can map is exercised for the Hessian invariants.
    const struct Case {
        ElementType type;
        BasisType basis_type;
        int order;
        math::Vector<double, 3> xi;
    } cases[] = {
        {ElementType::Line2, BasisType::Lagrange, 1, {double(0.15), double(0), double(0)}},
        {ElementType::Line3, BasisType::Lagrange, 2, {double(-0.25), double(0), double(0)}},
        {ElementType::Triangle3, BasisType::Lagrange, 1, {double(0.2), double(0.25), double(0)}},
        {ElementType::Triangle6, BasisType::Lagrange, 2, {double(0.2), double(0.25), double(0)}},
        {ElementType::Quad4, BasisType::Lagrange, 1, {double(0.3), double(-0.2), double(0)}},
        {ElementType::Quad8, BasisType::Serendipity, 2, {double(0.17), double(-0.31), double(0)}},
        {ElementType::Quad9, BasisType::Lagrange, 2, {double(0.3), double(-0.2), double(0)}},
        {ElementType::Tetra4, BasisType::Lagrange, 1, {double(0.15), double(0.2), double(0.1)}},
        {ElementType::Tetra10, BasisType::Lagrange, 2, {double(0.15), double(0.2), double(0.1)}},
        {ElementType::Hex8, BasisType::Lagrange, 1, {double(0.1), double(-0.2), double(0.3)}},
        {ElementType::Hex20, BasisType::Serendipity, 2, {double(0.2), double(-0.1), double(0.3)}},
        {ElementType::Hex27, BasisType::Lagrange, 2, {double(0.1), double(-0.2), double(0.3)}},
        {ElementType::Wedge6, BasisType::Lagrange, 1, {double(0.2), double(0.15), double(-0.3)}},
        {ElementType::Wedge15, BasisType::Serendipity, 2, {double(0.2), double(0.3), double(0.1)}},
        {ElementType::Wedge18, BasisType::Lagrange, 2, {double(0.2), double(0.15), double(-0.3)}},
    };

    for (const auto& c : cases) {
        auto basis = basis_factory::create(BasisRequest{c.type, c.basis_type, c.order});
        ASSERT_NE(basis, nullptr) << "element=" << static_cast<int>(c.type);
        expect_partition_hessian_sum_zero(*basis, c.xi, kHessianInvariantTol);
        expect_hessians_symmetric(*basis, c.xi, kHessianInvariantTol);
    }
}

// Gradients must match centered finite differences of values. This is the only
// check that ties the gradient code path back to the value code path; partition
// sums and Hessian-vs-FD(gradient) comparisons cannot catch a systematic error
// shared by the first- and second-derivative recurrences.
TEST(BasisGradients, LagrangeCanonicalTopologiesMatchNumericalGradients) {
    const struct Case {
        BasisTopology topology;
        int order;
        double tol;
    } cases[] = {
        {BasisTopology::Line, 3, double(1e-8)},
        {BasisTopology::Triangle, 3, double(1e-7)},
        {BasisTopology::Quadrilateral, 3, double(1e-7)},
        {BasisTopology::Tetrahedron, 2, double(1e-7)},
        {BasisTopology::Hexahedron, 2, double(1e-7)},
        {BasisTopology::Wedge, 2, double(1e-7)},
    };

    for (const auto& c : cases) {
        LagrangeBasis basis(c.topology, c.order);
        expect_gradients_match_numerical(basis, sample_points_for(c.topology), c.tol);
    }
}

// The serendipity coefficient tables (Hex20 20x20, Wedge15 15x15) and the quad
// inverse-Vandermonde path each differentiate values through hand-written code
// that is independent of the value evaluation. Partition sums only verify that
// the constant function differentiates to zero, and symmetry is assigned
// structurally, so neither can detect a wrong derivative formula. Finite
// differences of values are the authoritative check.
TEST(BasisGradients, SerendipityFamiliesMatchNumericalGradients) {
    // Arbitrary-order quadrilateral serendipity (topology path).
    const struct QuadCase { int order; double tol; } quad_cases[] = {
        {1, double(1e-8)}, {3, double(1e-7)}, {4, double(5e-7)}, {6, double(2e-6)},
    };
    for (const auto& c : quad_cases) {
        SerendipityBasis basis(BasisTopology::Quadrilateral, c.order);
        expect_gradients_match_numerical(
            basis, serendipity_sample_points(BasisTopology::Quadrilateral), c.tol);
    }

    // Arbitrary-order hexahedral serendipity (topology path).
    const struct HexCase { int order; double tol; } hex_cases[] = {
        {1, double(1e-8)}, {2, double(1e-7)}, {3, double(5e-7)},
        {4, double(1e-6)}, {5, double(5e-6)},
    };
    for (const auto& c : hex_cases) {
        SerendipityBasis basis(BasisTopology::Hexahedron, c.order);
        expect_gradients_match_numerical(
            basis, serendipity_sample_points(BasisTopology::Hexahedron), c.tol);
    }

    // Named fixed serendipity layouts.
    const struct NamedCase { ElementType type; int order; double tol; } named_cases[] = {
        {ElementType::Quad8, 2, double(1e-7)},
        {ElementType::Hex8, 1, double(1e-8)},
        {ElementType::Hex20, 2, double(1e-7)},
        {ElementType::Wedge15, 2, double(1e-7)},
    };
    for (const auto& c : named_cases) {
        SerendipityBasis basis(c.type, c.order);
        expect_gradients_match_numerical(
            basis, serendipity_sample_points(basis.topology()), c.tol);
    }
}

TEST(BasisGradients, QuadrilateralSerendipityInactiveZDerivativesRemainZero) {
    // All quadrilateral serendipity, including the production order 2, exercised
    // through the arbitrary-order topology path.
    for (const int order : {1, 2, 4, 6, 10}) {
        SerendipityBasis basis(BasisTopology::Quadrilateral, order);
        expect_inactive_z_derivatives_zero(
            basis,
            serendipity_sample_points(BasisTopology::Quadrilateral),
            double(1e-12));
    }
}

TEST(BasisHessians, SerendipityFamiliesMatchNumericalHessians) {
    // Arbitrary-order quadrilateral serendipity (topology path).
    const struct QuadCase { int order; double tol; } quad_cases[] = {
        {1, double(1e-6)}, {3, double(1e-6)}, {4, double(2e-6)}, {6, double(5e-6)},
    };
    for (const auto& c : quad_cases) {
        SerendipityBasis basis(BasisTopology::Quadrilateral, c.order);
        expect_hessians_match_numerical(
            basis, serendipity_sample_points(BasisTopology::Quadrilateral), c.tol);
    }

    // Arbitrary-order hexahedral serendipity (topology path).
    const struct HexCase { int order; double tol; } hex_cases[] = {
        {1, double(1e-6)}, {2, double(1e-6)}, {3, double(2e-6)},
        {4, double(5e-6)}, {5, double(1e-5)},
    };
    for (const auto& c : hex_cases) {
        SerendipityBasis basis(BasisTopology::Hexahedron, c.order);
        expect_hessians_match_numerical(
            basis, serendipity_sample_points(BasisTopology::Hexahedron), c.tol);
    }

    // Named fixed serendipity layouts.
    const struct NamedCase { ElementType type; int order; double tol; } named_cases[] = {
        {ElementType::Quad8, 2, double(1e-6)},
        {ElementType::Hex8, 1, double(1e-6)},
        {ElementType::Hex20, 2, double(1e-6)},
        {ElementType::Wedge15, 2, double(1e-6)},
    };
    for (const auto& c : named_cases) {
        SerendipityBasis basis(c.type, c.order);
        expect_hessians_match_numerical(
            basis, serendipity_sample_points(basis.topology()), c.tol);
    }
}
