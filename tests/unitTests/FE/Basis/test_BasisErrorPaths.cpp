/**
 * @file test_BasisErrorPaths.cpp
 * @brief Error-path coverage for the Lagrange-focused Basis subset.
 */

#include <gtest/gtest.h>

#include "FE/Basis/BasisExceptions.h"
#include "FE/Basis/BasisFactory.h"
#include "FE/Basis/BasisFunction.h"
#include "FE/Basis/LagrangeBasis.h"
#include "FE/Basis/NodeOrderingConventions.h"
#include "FE/Basis/SerendipityBasis.h"

#include <vector>

using namespace svmp::FE;
using namespace svmp::FE::basis;

namespace {

class MinimalScalarBasis : public BasisFunction {
public:
    BasisType basis_type() const noexcept override { return BasisType::Lagrange; }
    ElementType element_type() const noexcept override { return ElementType::Line2; }
    int dimension() const noexcept override { return 1; }
    int order() const noexcept override { return 1; }
    std::size_t size() const noexcept override { return 2u; }

    void evaluate_values(const math::Vector<Real, 3>&,
                         std::vector<Real>& values) const override
    {
        values.assign(size(), Real(0));
    }
};

// Quadratic scalar basis with exact analytic derivatives, used to verify the
// protected numerical_gradient/numerical_hessian development helpers. Centered
// differences are exact (up to roundoff) on quadratics, so any mismatch is a
// bug in the helpers themselves.
class ExactQuadraticBasis : public BasisFunction {
public:
    using BasisFunction::numerical_gradient;
    using BasisFunction::numerical_hessian;

    BasisType basis_type() const noexcept override { return BasisType::Custom; }
    ElementType element_type() const noexcept override { return ElementType::Hex8; }
    int dimension() const noexcept override { return 3; }
    int order() const noexcept override { return 2; }
    std::size_t size() const noexcept override { return 2u; }

    void evaluate_values(const math::Vector<Real, 3>& xi,
                         std::vector<Real>& values) const override
    {
        const Real x = xi[0];
        const Real y = xi[1];
        const Real z = xi[2];
        values.resize(size());
        values[0] = Real(1) + Real(2) * x - y + Real(0.5) * z +
                    x * x + Real(0.75) * y * y - Real(0.25) * z * z +
                    Real(0.2) * x * y - Real(0.3) * x * z + Real(0.4) * y * z;
        values[1] = Real(3) - x + Real(2) * y + z +
                    Real(0.5) * x * x - y * y + z * z +
                    x * y + x * z - y * z;
    }

    void evaluate_gradients(const math::Vector<Real, 3>& xi,
                            std::vector<Gradient>& gradients) const override
    {
        const Real x = xi[0];
        const Real y = xi[1];
        const Real z = xi[2];
        gradients.assign(size(), Gradient::Zero());
        gradients[0][0] = Real(2) + Real(2) * x + Real(0.2) * y - Real(0.3) * z;
        gradients[0][1] = Real(-1) + Real(1.5) * y + Real(0.2) * x + Real(0.4) * z;
        gradients[0][2] = Real(0.5) - Real(0.5) * z - Real(0.3) * x + Real(0.4) * y;
        gradients[1][0] = Real(-1) + x + y + z;
        gradients[1][1] = Real(2) - Real(2) * y + x - z;
        gradients[1][2] = Real(1) + Real(2) * z + x - y;
    }

    void exact_hessians(std::vector<Hessian>& hessians) const
    {
        hessians.assign(size(), Hessian::Zero());
        hessians[0] = make_symmetric_hessian(Real(2), Real(1.5), Real(-0.5),
                                             Real(0.2), Real(-0.3), Real(0.4));
        hessians[1] = make_symmetric_hessian(Real(1), Real(-2), Real(2),
                                             Real(1), Real(1), Real(-1));
    }
};

class CompleteFallbackBasis : public BasisFunction {
public:
    BasisType basis_type() const noexcept override { return BasisType::Lagrange; }
    ElementType element_type() const noexcept override { return ElementType::Triangle3; }
    int dimension() const noexcept override { return 2; }
    int order() const noexcept override { return 1; }
    std::size_t size() const noexcept override { return 2u; }

    void evaluate_values(const math::Vector<Real, 3>& xi,
                         std::vector<Real>& values) const override
    {
        values.resize(size());
        values[0] = Real(1) + xi[0];
        values[1] = Real(2) + xi[1];
    }

    void evaluate_gradients(const math::Vector<Real, 3>&,
                            std::vector<Gradient>& gradients) const override
    {
        gradients.assign(size(), Gradient::Zero());
        gradients[0][0] = Real(1);
        gradients[1][1] = Real(1);
    }

    void evaluate_hessians(const math::Vector<Real, 3>& xi,
                           std::vector<Hessian>& hessians) const override
    {
        hessians.assign(size(), Hessian::Zero());
        for (std::size_t d = 0; d < hessians.size(); ++d) {
            for (std::size_t r = 0; r < 3u; ++r) {
                for (std::size_t c = 0; c < 3u; ++c) {
                    hessians[d](r, c) = Real(100) * static_cast<Real>(d + 1u) +
                                        Real(10) * static_cast<Real>(r) +
                                        static_cast<Real>(c) + xi[2];
                }
            }
        }
    }
};

} // namespace

TEST(BasisErrorPaths, LagrangeInvalidRequestsThrowBasisExceptions) {
    EXPECT_THROW(LagrangeBasis(ElementType::Unknown, 1),
                 BasisElementCompatibilityException);
    EXPECT_THROW(LagrangeBasis(ElementType::Line2, -1),
                 BasisConfigurationException);
    EXPECT_THROW(LagrangeBasis(ElementType::Quad8, 2),
                 BasisElementCompatibilityException);
}

TEST(BasisErrorPaths, SerendipityInvalidRequestsThrowBasisExceptions) {
    EXPECT_THROW(SerendipityBasis(ElementType::Unknown, 2),
                 BasisElementCompatibilityException);
    EXPECT_THROW(SerendipityBasis(ElementType::Quad8, 3),
                 BasisConfigurationException);
    EXPECT_THROW(SerendipityBasis(ElementType::Pyramid13, 2),
                 BasisElementCompatibilityException);
    EXPECT_THROW(SerendipityBasis(ElementType::Pyramid14, 2),
                 BasisElementCompatibilityException);
}

TEST(BasisErrorPaths, BasisFactoryRejectsNonC0Continuity) {
    BasisRequest c1_request{ElementType::Line2, BasisType::Lagrange, 1};
    c1_request.continuity = Continuity::C1;
    EXPECT_THROW((void)basis_factory::create(c1_request), BasisConfigurationException);

    BasisRequest l2_request{ElementType::Quad8, BasisType::Serendipity, 2};
    l2_request.continuity = Continuity::L2;
    EXPECT_THROW((void)basis_factory::create(l2_request), BasisConfigurationException);
}

TEST(BasisErrorPaths, BasisFactoryInvalidRequestsThrowBasisExceptions) {
    EXPECT_THROW((void)basis_factory::create(
                     BasisRequest{ElementType::Line2, BasisType::Lagrange}),
                 BasisConfigurationException);
    EXPECT_THROW((void)basis_factory::create(
                     BasisRequest{ElementType::Line2, BasisType::Lagrange, -1}),
                 BasisConfigurationException);
    EXPECT_THROW((void)basis_factory::create(
                     BasisRequest{ElementType::Line2, BasisType::Bernstein, 1}),
                 BasisConfigurationException);
    EXPECT_THROW((void)basis_factory::create(
                     BasisRequest{ElementType::Pyramid5, BasisType::Lagrange, 1}),
                 BasisElementCompatibilityException);

    BasisRequest vector_req{ElementType::Line2, BasisType::Lagrange, 1};
    vector_req.field_type = FieldType::Vector;
    EXPECT_THROW((void)basis_factory::create(vector_req), BasisConfigurationException);

    auto serendipity = basis_factory::create(
        BasisRequest{ElementType::Quad8, BasisType::Serendipity, 2});
    ASSERT_NE(serendipity, nullptr);
    EXPECT_EQ(serendipity->basis_type(), BasisType::Serendipity);
}

TEST(BasisErrorPaths, BasisExceptionsUseCommonStatusCodes) {
    try {
        throw BasisConfigurationException("invalid config", __FILE__, __LINE__, __func__);
    } catch (const FEException& e) {
        EXPECT_EQ(e.status(), svmp::StatusCode::InvalidArgument);
    }

    try {
        throw BasisConstructionException("construction failure", __FILE__, __LINE__, __func__);
    } catch (const FEException& e) {
        EXPECT_EQ(e.status(), svmp::StatusCode::InternalError);
    }
}

TEST(BasisErrorPaths, NodeOrderingInvalidNodeThrows) {
    EXPECT_THROW((void)ReferenceNodeLayout::get_node_coords(ElementType::Quad8, 99u),
                 BasisNodeOrderingException);
    EXPECT_THROW((void)ReferenceNodeLayout::get_lagrange_node_coords(ElementType::Quad8, 2),
                 BasisNodeOrderingException);
    EXPECT_THROW((void)ReferenceNodeLayout::num_nodes(ElementType::Pyramid5),
                 BasisNodeOrderingException);
}

TEST(BasisErrorPaths, BasisFunctionDefaultsThrowForMissingDerivatives) {
    MinimalScalarBasis basis;
    const math::Vector<Real, 3> xi{Real(0), Real(0), Real(0)};
    std::vector<Gradient> gradients;
    std::vector<Hessian> hessians;

    EXPECT_THROW(basis.evaluate_gradients(xi, gradients), BasisEvaluationException);
    EXPECT_THROW(basis.evaluate_hessians(xi, hessians), BasisEvaluationException);
}

TEST(BasisErrorPaths, NumericalDerivativeHelpersMatchAnalyticDerivatives) {
    ExactQuadraticBasis basis;
    const math::Vector<Real, 3> xi{Real(0.2), Real(-0.35), Real(0.4)};

    std::vector<Gradient> exact_gradients;
    basis.evaluate_gradients(xi, exact_gradients);

    std::vector<Gradient> approx_gradients;
    basis.numerical_gradient(xi, approx_gradients);
    ASSERT_EQ(approx_gradients.size(), basis.size());
    for (std::size_t n = 0; n < basis.size(); ++n) {
        for (int d = 0; d < basis.dimension(); ++d) {
            const std::size_t sd = static_cast<std::size_t>(d);
            EXPECT_NEAR(approx_gradients[n][sd], exact_gradients[n][sd], Real(1e-8))
                << "basis=" << n << " component=" << d;
        }
    }

    std::vector<Hessian> exact_hessians;
    basis.exact_hessians(exact_hessians);

    std::vector<Hessian> approx_hessians;
    basis.numerical_hessian(xi, approx_hessians);
    ASSERT_EQ(approx_hessians.size(), basis.size());
    for (std::size_t n = 0; n < basis.size(); ++n) {
        for (int r = 0; r < basis.dimension(); ++r) {
            for (int c = 0; c < basis.dimension(); ++c) {
                const std::size_t sr = static_cast<std::size_t>(r);
                const std::size_t sc = static_cast<std::size_t>(c);
                EXPECT_NEAR(approx_hessians[n](sr, sc), exact_hessians[n](sr, sc),
                            Real(1e-8))
                    << "basis=" << n << " component=(" << r << "," << c << ")";
            }
        }
    }
}

TEST(BasisErrorPaths, BasisFunctionFallbackWritesRawLayouts) {
    CompleteFallbackBasis basis;
    const math::Vector<Real, 3> point{Real(0.25), Real(0.5), Real(-0.25)};

    std::vector<Real> flat_values(basis.size());
    std::vector<Real> flat_gradients(basis.size() * 3u);
    std::vector<Real> flat_hessians(basis.size() * 9u);
    basis.evaluate_values_to(point, flat_values.data());
    basis.evaluate_gradients_to(point, flat_gradients.data());
    basis.evaluate_hessians_to(point, flat_hessians.data());

    std::vector<Real> expected_values;
    std::vector<Gradient> expected_gradients;
    std::vector<Hessian> expected_hessians;
    basis.evaluate_all(point, expected_values, expected_gradients, expected_hessians);
    for (std::size_t d = 0; d < basis.size(); ++d) {
        EXPECT_EQ(flat_values[d], expected_values[d]);
        for (std::size_t c = 0; c < 3u; ++c) {
            EXPECT_EQ(flat_gradients[d * 3u + c], expected_gradients[d][c]);
        }
        for (std::size_t r = 0; r < 3u; ++r) {
            for (std::size_t c = 0; c < 3u; ++c) {
                EXPECT_EQ(flat_hessians[d * 9u + r * 3u + c], expected_hessians[d](r, c));
            }
        }
    }
}
