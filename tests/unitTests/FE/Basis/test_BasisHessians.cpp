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

void numerical_gradient_helper(const BasisFunction& basis,
                               const math::Vector<Real, 3>& xi,
                               std::vector<Gradient>& gradients,
                               Real eps = Real(1e-6))
{
    std::vector<Real> base;
    basis.evaluate_values(xi, base);
    gradients.assign(base.size(), Gradient::Zero());

    for (int d = 0; d < basis.dimension(); ++d) {
        const std::size_t sd = static_cast<std::size_t>(d);
        math::Vector<Real, 3> xi_p = xi;
        math::Vector<Real, 3> xi_m = xi;
        xi_p[sd] += eps;
        xi_m[sd] -= eps;

        std::vector<Real> v_p;
        std::vector<Real> v_m;
        basis.evaluate_values(xi_p, v_p);
        basis.evaluate_values(xi_m, v_m);

        for (std::size_t n = 0; n < base.size(); ++n) {
            gradients[n][sd] = (v_p[n] - v_m[n]) / (Real(2) * eps);
        }
    }
}

void numerical_hessian_helper(const BasisFunction& basis,
                              const math::Vector<Real, 3>& xi,
                              std::vector<Hessian>& hessians,
                              Real eps = Real(1e-5))
{
    hessians.assign(basis.size(), Hessian::Zero());
    const int dim = basis.dimension();

    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            math::Vector<Real, 3> xi_p = xi;
            math::Vector<Real, 3> xi_m = xi;
            const std::size_t sj = static_cast<std::size_t>(j);
            xi_p[sj] += eps;
            xi_m[sj] -= eps;

            std::vector<Gradient> g_p;
            std::vector<Gradient> g_m;
            basis.evaluate_gradients(xi_p, g_p);
            basis.evaluate_gradients(xi_m, g_m);

            for (std::size_t n = 0; n < basis.size(); ++n) {
                const std::size_t si = static_cast<std::size_t>(i);
                hessians[n](si, sj) = (g_p[n][si] - g_m[n][si]) / (Real(2) * eps);
            }
        }
    }
}

std::vector<math::Vector<Real, 3>> sample_points_for(ElementType type) {
    switch (type) {
        case ElementType::Line2:
            return {{Real(-0.35), Real(0), Real(0)}, {Real(0.2), Real(0), Real(0)}};
        case ElementType::Triangle3:
            return {{Real(0.15), Real(0.2), Real(0)}, {Real(0.25), Real(0.1), Real(0)}};
        case ElementType::Quad4:
            return {{Real(0.2), Real(-0.3), Real(0)}, {Real(-0.45), Real(0.25), Real(0)}};
        case ElementType::Tetra4:
            return {{Real(0.12), Real(0.18), Real(0.16)}, {Real(0.2), Real(0.1), Real(0.18)}};
        case ElementType::Hex8:
            return {{Real(0.1), Real(-0.2), Real(0.3)}, {Real(-0.35), Real(0.25), Real(-0.15)}};
        case ElementType::Wedge6:
            return {{Real(0.18), Real(0.22), Real(-0.2)}, {Real(0.12), Real(0.16), Real(0.1)}};
        default:
            return {{Real(0), Real(0), Real(0)}};
    }
}

void expect_gradients_match_numerical(const BasisFunction& basis,
                                      const std::vector<math::Vector<Real, 3>>& points,
                                      Real tol,
                                      Real eps = Real(1e-6))
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
                    << ", element " << static_cast<int>(basis.element_type())
                    << ", order " << basis.order();
            }
        }
    }
}

void expect_hessians_match_numerical(const BasisFunction& basis,
                                     const std::vector<math::Vector<Real, 3>>& points,
                                     Real tol,
                                     Real eps = Real(1e-5))
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
                        << "), element " << static_cast<int>(basis.element_type())
                        << ", order " << basis.order();
                }
            }
        }
    }
}

void expect_partition_hessian_sum_zero(const BasisFunction& basis,
                                       const math::Vector<Real, 3>& xi,
                                       Real tol)
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
                        Real(0),
                        tol)
                << "element " << static_cast<int>(basis.element_type())
                << ", order " << basis.order();
        }
    }
}

void expect_hessians_symmetric(const BasisFunction& basis,
                               const math::Vector<Real, 3>& xi,
                               Real tol)
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

std::vector<math::Vector<Real, 3>> serendipity_sample_points(ElementType type) {
    if (type == ElementType::Quad4 || type == ElementType::Quad8) {
        return {{Real(0.17), Real(-0.31), Real(0)}, {Real(-0.45), Real(0.25), Real(0)}};
    }
    if (type == ElementType::Hex8 || type == ElementType::Hex20) {
        return {{Real(0.2), Real(-0.1), Real(0.3)}, {Real(-0.35), Real(0.25), Real(-0.15)}};
    }
    return {{Real(0.2), Real(0.3), Real(0.1)}, {Real(0.12), Real(0.16), Real(-0.2)}};
}

} // namespace

TEST(BasisHessians, LagrangeCanonicalTopologiesMatchNumericalHessians) {
    const struct Case {
        ElementType type;
        int order;
        Real tol;
        Real eps;
    } cases[] = {
        {ElementType::Line2, 3, Real(1e-7), Real(1e-5)},
        {ElementType::Triangle3, 3, Real(2e-6), Real(1e-5)},
        {ElementType::Quad4, 3, Real(1e-6), Real(1e-5)},
        {ElementType::Tetra4, 2, Real(1e-6), Real(1e-5)},
        {ElementType::Hex8, 2, Real(1e-6), Real(1e-5)},
        {ElementType::Wedge6, 2, Real(1e-5), Real(1e-5)},
    };

    for (const auto& c : cases) {
        LagrangeBasis basis(c.type, c.order);
        expect_hessians_match_numerical(basis, sample_points_for(c.type), c.tol, c.eps);
    }
}

TEST(BasisHessians, LagrangeHessiansSumToZeroAndAreSymmetric) {
    const struct Case {
        ElementType type;
        int order;
        math::Vector<Real, 3> xi;
        Real tol;
    } cases[] = {
        {ElementType::Line2, 3, {Real(0.15), Real(0), Real(0)}, Real(1e-12)},
        {ElementType::Triangle3, 3, {Real(0.2), Real(0.25), Real(0)}, Real(1e-10)},
        {ElementType::Quad4, 3, {Real(0.3), Real(-0.2), Real(0)}, Real(1e-12)},
        {ElementType::Tetra4, 2, {Real(0.15), Real(0.2), Real(0.1)}, Real(1e-10)},
        {ElementType::Hex8, 2, {Real(0.1), Real(-0.2), Real(0.3)}, Real(1e-12)},
        {ElementType::Wedge6, 2, {Real(0.2), Real(0.15), Real(-0.3)}, Real(1e-10)},
    };

    for (const auto& c : cases) {
        LagrangeBasis basis(c.type, c.order);
        expect_partition_hessian_sum_zero(basis, c.xi, Real(10) * c.tol);
        expect_hessians_symmetric(basis, c.xi, c.tol);
    }
}

TEST(BasisHessians, SerendipityHessiansSumToZeroAndAreSymmetric) {
    const struct Case {
        ElementType type;
        int order;
        math::Vector<Real, 3> xi;
        Real tol;
    } cases[] = {
        {ElementType::Quad8, 2, {Real(0.17), Real(-0.31), Real(0)}, Real(1e-10)},
        {ElementType::Hex20, 2, {Real(0.2), Real(-0.1), Real(0.3)}, Real(1e-10)},
        {ElementType::Wedge15, 2, {Real(0.2), Real(0.3), Real(0.1)}, Real(1e-10)},
    };

    for (const auto& c : cases) {
        SerendipityBasis basis(c.type, c.order);
        expect_partition_hessian_sum_zero(basis, c.xi, c.tol);
        expect_hessians_symmetric(basis, c.xi, c.tol);
    }
}

TEST(BasisHessians, SolverMappedVolumeSelectionsSatisfyInvariants) {
    const struct Case {
        ElementType type;
        BasisType basis_type;
        int order;
        math::Vector<Real, 3> xi;
        Real tol;
    } cases[] = {
        {ElementType::Line2, BasisType::Lagrange, 1, {Real(0.15), Real(0), Real(0)}, Real(1e-12)},
        {ElementType::Line3, BasisType::Lagrange, 2, {Real(-0.25), Real(0), Real(0)}, Real(1e-12)},
        {ElementType::Triangle3, BasisType::Lagrange, 1, {Real(0.2), Real(0.25), Real(0)}, Real(1e-12)},
        {ElementType::Triangle6, BasisType::Lagrange, 2, {Real(0.2), Real(0.25), Real(0)}, Real(1e-12)},
        {ElementType::Quad4, BasisType::Lagrange, 1, {Real(0.3), Real(-0.2), Real(0)}, Real(1e-12)},
        {ElementType::Quad8, BasisType::Serendipity, 2, {Real(0.17), Real(-0.31), Real(0)}, Real(1e-10)},
        {ElementType::Quad9, BasisType::Lagrange, 2, {Real(0.3), Real(-0.2), Real(0)}, Real(1e-12)},
        {ElementType::Tetra4, BasisType::Lagrange, 1, {Real(0.15), Real(0.2), Real(0.1)}, Real(1e-12)},
        {ElementType::Tetra10, BasisType::Lagrange, 2, {Real(0.15), Real(0.2), Real(0.1)}, Real(1e-10)},
        {ElementType::Hex8, BasisType::Lagrange, 1, {Real(0.1), Real(-0.2), Real(0.3)}, Real(1e-12)},
        {ElementType::Hex20, BasisType::Serendipity, 2, {Real(0.2), Real(-0.1), Real(0.3)}, Real(1e-10)},
        {ElementType::Hex27, BasisType::Lagrange, 2, {Real(0.1), Real(-0.2), Real(0.3)}, Real(1e-12)},
        {ElementType::Wedge6, BasisType::Lagrange, 1, {Real(0.2), Real(0.15), Real(-0.3)}, Real(1e-12)},
    };

    for (const auto& c : cases) {
        auto basis = basis_factory::create(BasisRequest{c.type, c.basis_type, c.order});
        expect_partition_hessian_sum_zero(*basis, c.xi, c.tol);
        expect_hessians_symmetric(*basis, c.xi, c.tol);
    }
}

// Gradients must match centered finite differences of values. This is the only
// check that ties the gradient code path back to the value code path; partition
// sums and Hessian-vs-FD(gradient) comparisons cannot catch a systematic error
// shared by the first- and second-derivative recurrences.
TEST(BasisGradients, LagrangeCanonicalTopologiesMatchNumericalGradients) {
    const struct Case {
        ElementType type;
        int order;
        Real tol;
    } cases[] = {
        {ElementType::Line2, 3, Real(1e-8)},
        {ElementType::Triangle3, 3, Real(1e-7)},
        {ElementType::Quad4, 3, Real(1e-7)},
        {ElementType::Tetra4, 2, Real(1e-7)},
        {ElementType::Hex8, 2, Real(1e-7)},
        {ElementType::Wedge6, 2, Real(1e-7)},
    };

    for (const auto& c : cases) {
        LagrangeBasis basis(c.type, c.order);
        expect_gradients_match_numerical(basis, sample_points_for(c.type), c.tol);
    }
}

// The serendipity coefficient tables (Hex20 20x20, Wedge15 15x15) and the quad
// inverse-Vandermonde path each differentiate values through hand-written code
// that is independent of the value evaluation. Partition sums only verify that
// the constant function differentiates to zero, and symmetry is assigned
// structurally, so neither can detect a wrong derivative formula. Finite
// differences of values are the authoritative check.
TEST(BasisGradients, SerendipityFamiliesMatchNumericalGradients) {
    const struct Case {
        ElementType type;
        int order;
        Real tol;
    } cases[] = {
        {ElementType::Quad4, 1, Real(1e-8)},
        {ElementType::Quad8, 2, Real(1e-7)},
        {ElementType::Quad4, 3, Real(1e-7)},
        {ElementType::Quad4, 4, Real(5e-7)},
        {ElementType::Hex8, 1, Real(1e-8)},
        {ElementType::Hex20, 2, Real(1e-7)},
        {ElementType::Wedge15, 2, Real(1e-7)},
    };

    for (const auto& c : cases) {
        SerendipityBasis basis(c.type, c.order);
        expect_gradients_match_numerical(basis, serendipity_sample_points(c.type), c.tol);
    }
}

TEST(BasisHessians, SerendipityFamiliesMatchNumericalHessians) {
    const struct Case {
        ElementType type;
        int order;
        Real tol;
    } cases[] = {
        {ElementType::Quad4, 1, Real(1e-6)},
        {ElementType::Quad8, 2, Real(1e-6)},
        {ElementType::Quad4, 3, Real(1e-6)},
        {ElementType::Quad4, 4, Real(5e-6)},
        {ElementType::Hex8, 1, Real(1e-6)},
        {ElementType::Hex20, 2, Real(1e-6)},
        {ElementType::Wedge15, 2, Real(1e-6)},
    };

    for (const auto& c : cases) {
        SerendipityBasis basis(c.type, c.order);
        expect_hessians_match_numerical(basis, serendipity_sample_points(c.type), c.tol);
    }
}
