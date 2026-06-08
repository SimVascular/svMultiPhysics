/**
 * @file test_BasisErrorPaths.cpp
 * @brief Error-path coverage for the migrated Lagrange-focused Basis subset.
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
        gradients.assign(size(), Gradient{});
        gradients[0][0] = Real(1);
        gradients[1][1] = Real(1);
    }

    void evaluate_hessians(const math::Vector<Real, 3>& xi,
                           std::vector<Hessian>& hessians) const override
    {
        hessians.assign(size(), Hessian{});
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

TEST(BasisErrorPaths, BasisFunctionFallbackWritesRawLayouts) {
    CompleteFallbackBasis basis;
    const math::Vector<Real, 3> point{Real(0.25), Real(0.5), Real(-0.25)};
    prewarm_basis_function_scratch(basis.size());

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
