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

#include <algorithm>
#include <span>
#include <string>
#include <vector>

using namespace svmp::FE;
using namespace svmp::FE::basis;

namespace {

// Build a symmetric 3x3 Hessian from its six independent components. Local to
// this test; the production basis evaluators fill Hessians directly.
[[nodiscard]] Hessian make_symmetric_hessian(double xx,
                                             double yy,
                                             double zz,
                                             double xy,
                                             double xz,
                                             double yz) {
    Hessian hessian = Hessian::Zero();
    hessian(0, 0) = xx;
    hessian(1, 1) = yy;
    hessian(2, 2) = zz;
    hessian(0, 1) = hessian(1, 0) = xy;
    hessian(0, 2) = hessian(2, 0) = xz;
    hessian(1, 2) = hessian(2, 1) = yz;
    return hessian;
}

class MinimalScalarBasis : public BasisFunction {
public:
    BasisType basis_type() const noexcept override { return BasisType::Lagrange; }
    BasisTopology topology() const noexcept override { return BasisTopology::Line; }
    int dimension() const noexcept override { return 1; }
    int order() const noexcept override { return 1; }
    std::size_t size() const noexcept override { return 2u; }

    void evaluate_values_to(const math::Vector<double, 3>&,
                            std::span<double> values_out) const override
    {
        std::fill(values_out.begin(), values_out.end(), double(0));
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
    BasisTopology topology() const noexcept override { return BasisTopology::Hexahedron; }
    int dimension() const noexcept override { return 3; }
    int order() const noexcept override { return 2; }
    std::size_t size() const noexcept override { return 2u; }

    void evaluate_values_to(const math::Vector<double, 3>& xi,
                            std::span<double> values_out) const override
    {
        const double x = xi[0];
        const double y = xi[1];
        const double z = xi[2];
        values_out[0] = double(1) + double(2) * x - y + double(0.5) * z +
                        x * x + double(0.75) * y * y - double(0.25) * z * z +
                        double(0.2) * x * y - double(0.3) * x * z + double(0.4) * y * z;
        values_out[1] = double(3) - x + double(2) * y + z +
                        double(0.5) * x * x - y * y + z * z +
                        x * y + x * z - y * z;
    }

    void evaluate_gradients_to(const math::Vector<double, 3>& xi,
                               std::span<Gradient> gradients_out) const override
    {
        const double x = xi[0];
        const double y = xi[1];
        const double z = xi[2];
        gradients_out[0] = Gradient::Zero();
        gradients_out[1] = Gradient::Zero();
        gradients_out[0][0] = double(2) + double(2) * x + double(0.2) * y - double(0.3) * z;
        gradients_out[0][1] = double(-1) + double(1.5) * y + double(0.2) * x + double(0.4) * z;
        gradients_out[0][2] = double(0.5) - double(0.5) * z - double(0.3) * x + double(0.4) * y;
        gradients_out[1][0] = double(-1) + x + y + z;
        gradients_out[1][1] = double(2) - double(2) * y + x - z;
        gradients_out[1][2] = double(1) + double(2) * z + x - y;
    }

    void exact_hessians(std::vector<Hessian>& hessians) const
    {
        hessians.assign(size(), Hessian::Zero());
        hessians[0] = make_symmetric_hessian(double(2), double(1.5), double(-0.5),
                                             double(0.2), double(-0.3), double(0.4));
        hessians[1] = make_symmetric_hessian(double(1), double(-2), double(2),
                                             double(1), double(1), double(-1));
    }
};

// Basis that implements only the span primitives and deliberately does not
// override the combined evaluate_all_to. It therefore exercises the base class's
// vector overloads and the default combined evaluator, both of which must forward
// to these primitives.
class SpanPrimitiveBasis : public BasisFunction {
public:
    BasisType basis_type() const noexcept override { return BasisType::Lagrange; }
    BasisTopology topology() const noexcept override { return BasisTopology::Triangle; }
    int dimension() const noexcept override { return 2; }
    int order() const noexcept override { return 1; }
    std::size_t size() const noexcept override { return 2u; }

    void evaluate_values_to(const math::Vector<double, 3>& xi,
                            std::span<double> values_out) const override
    {
        values_out[0] = double(1) + xi[0];
        values_out[1] = double(2) + xi[1];
    }

    void evaluate_gradients_to(const math::Vector<double, 3>&,
                               std::span<Gradient> gradients_out) const override
    {
        gradients_out[0] = Gradient::Zero();
        gradients_out[1] = Gradient::Zero();
        gradients_out[0][0] = double(1);
        gradients_out[1][1] = double(1);
    }

    void evaluate_hessians_to(const math::Vector<double, 3>& xi,
                              std::span<Hessian> hessians_out) const override
    {
        for (std::size_t d = 0; d < size(); ++d) {
            hessians_out[d] = Hessian::Zero();
            for (std::size_t r = 0; r < 3u; ++r) {
                for (std::size_t c = 0; c < 3u; ++c) {
                    hessians_out[d](r, c) = double(100) * static_cast<double>(d + 1u) +
                                            double(10) * static_cast<double>(r) +
                                            static_cast<double>(c) + xi[2];
                }
            }
        }
    }
};

void expect_source_location(const svmp::ExceptionBase& e)
{
    EXPECT_NE(e.context().file().find("test_BasisErrorPaths.cpp"), std::string::npos);
    EXPECT_GT(e.context().line(), 0);
    EXPECT_FALSE(e.context().function().empty());
}

// The core helpers raise both FE-subsystem exceptions and Core exceptions (for
// example not_implemented() defaults to svmp::NotImplementedException), so this
// catches their common base, svmp::ExceptionBase.
template <class Thrower>
void expect_core_helper_preserves_source_location(Thrower&& thrower)
{
    try {
        thrower();
        FAIL() << "Expected an svmp::ExceptionBase";
    } catch (const svmp::ExceptionBase& e) {
        expect_source_location(e);
    }
}

} // namespace

TEST(BasisErrorPaths, LagrangeInvalidRequestsThrowBasisExceptions) {
    EXPECT_THROW(LagrangeBasis(ElementType::Unknown, 1),
                 BasisElementCompatibilityException);
    EXPECT_THROW(LagrangeBasis(ElementType::Line2, -1),
                 BasisConfigurationException);
    EXPECT_THROW(LagrangeBasis(ElementType::Quad8, 2),
                 BasisElementCompatibilityException);
    EXPECT_NO_THROW((void)LagrangeBasis(BasisTopology::Point, 0));
    EXPECT_THROW((void)LagrangeBasis(BasisTopology::Point, 1),
                 BasisConfigurationException);
}

// A named Lagrange element layout fixes its polynomial order: the matching order
// is accepted and any other order is rejected. Arbitrary orders must be
// requested through the BasisTopology overload, never by over-/under-specifying
// a node-count-named element.
TEST(BasisErrorPaths, NamedLagrangeElementsRejectNonBakedOrders) {
    const std::vector<std::pair<ElementType, int>> named = {
        {ElementType::Point1, 0},
        {ElementType::Line2, 1},     {ElementType::Line3, 2},
        {ElementType::Triangle3, 1}, {ElementType::Triangle6, 2},
        {ElementType::Quad4, 1},     {ElementType::Quad9, 2},
        {ElementType::Tetra4, 1},    {ElementType::Tetra10, 2},
        {ElementType::Hex8, 1},      {ElementType::Hex27, 2},
        {ElementType::Wedge6, 1},    {ElementType::Wedge18, 2},
    };

    for (const auto& [type, baked] : named) {
        EXPECT_NO_THROW((void)LagrangeBasis(type, baked))
            << "element=" << static_cast<int>(type);
        EXPECT_THROW((void)LagrangeBasis(type, baked + 1), BasisConfigurationException)
            << "element=" << static_cast<int>(type);
        if (baked > 0) {
            EXPECT_THROW((void)LagrangeBasis(type, baked - 1), BasisConfigurationException)
                << "element=" << static_cast<int>(type);
        }
    }
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
    EXPECT_THROW(SerendipityBasis(ElementType::Hex8, 2),
                 BasisConfigurationException);
    EXPECT_THROW(SerendipityBasis(ElementType::Hex20, 1),
                 BasisConfigurationException);
    EXPECT_THROW(SerendipityBasis(ElementType::Hex20, 3),
                 BasisConfigurationException);
    EXPECT_THROW(SerendipityBasis(ElementType::Wedge15, 1),
                 BasisConfigurationException);
    EXPECT_THROW(SerendipityBasis(ElementType::Wedge15, 3),
                 BasisConfigurationException);

    // Order 0 and negative orders are rejected for every serendipity layout; a
    // named element is pinned to its inferred order and is never floored up to it.
    EXPECT_THROW(SerendipityBasis(ElementType::Quad8, 0),
                 BasisConfigurationException);
    EXPECT_THROW(SerendipityBasis(ElementType::Hex8, 0),
                 BasisConfigurationException);
    EXPECT_THROW(SerendipityBasis(ElementType::Hex20, 0),
                 BasisConfigurationException);
    EXPECT_THROW(SerendipityBasis(ElementType::Wedge15, 0),
                 BasisConfigurationException);
    EXPECT_THROW(SerendipityBasis(ElementType::Hex8, -1),
                 BasisConfigurationException);
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
    // NURBS is a declared but unimplemented family, so the scalar factory rejects
    // it as outside the Lagrange/Serendipity scope.
    EXPECT_THROW((void)basis_factory::create(
                     BasisRequest{ElementType::Line2, BasisType::NURBS, 1}),
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
        svmp::raise<BasisConfigurationException>("invalid config");
    } catch (const FEException& e) {
        EXPECT_EQ(e.status(), svmp::StatusCode::InvalidArgument);
    }

    try {
        svmp::raise<BasisConstructionException>("construction failure");
    } catch (const FEException& e) {
        EXPECT_EQ(e.status(), svmp::StatusCode::InternalError);
    }
}

TEST(BasisErrorPaths, CoreHelpersPreserveSourceLocation) {
    expect_core_helper_preserves_source_location([] {
        svmp::raise<BasisEvaluationException>("raise location");
    });

    expect_core_helper_preserves_source_location([] {
        svmp::throw_if<BasisEvaluationException>(
            true, "throw_if location");
    });

    expect_core_helper_preserves_source_location([] {
        svmp::check<BasisEvaluationException>(
            false, "check location");
    });

    expect_core_helper_preserves_source_location([] {
        const int* ptr = nullptr;
        svmp::check_not_null<BasisEvaluationException>(
            ptr, "check_not_null location");
    });

    expect_core_helper_preserves_source_location([] {
        svmp::check_index<BasisEvaluationException>(1, 1);
    });

    expect_core_helper_preserves_source_location([] {
        svmp::not_implemented<NotImplementedException>(
            "test feature");
    });
}

TEST(BasisErrorPaths, NodeOrderingInvalidNodeThrows) {
    EXPECT_THROW((void)line_coord_pm_one(-1, 1),
                 BasisNodeOrderingException);
    EXPECT_THROW((void)line_coord_pm_one(2, 1),
                 BasisNodeOrderingException);
    EXPECT_THROW((void)line_coord_pm_one(-1, 0),
                 BasisNodeOrderingException);
    EXPECT_THROW((void)line_coord_pm_one(1, 0),
                 BasisNodeOrderingException);
    EXPECT_THROW((void)ReferenceNodeLayout::node_coord_at(ElementType::Quad8, 99u),
                 BasisNodeOrderingException);
    EXPECT_THROW((void)ReferenceNodeLayout::get_lagrange_node_coords(ElementType::Quad8, 2),
                 BasisNodeOrderingException);
    EXPECT_THROW((void)ReferenceNodeLayout::num_nodes(ElementType::Pyramid5),
                 BasisNodeOrderingException);
}

TEST(BasisErrorPaths, BasisFunctionDefaultsThrowForMissingDerivatives) {
    MinimalScalarBasis basis;
    const math::Vector<double, 3> xi{double(0), double(0), double(0)};
    std::vector<Gradient> gradients;
    std::vector<Hessian> hessians;

    EXPECT_THROW(basis.evaluate_gradients(xi, gradients), BasisEvaluationException);
    EXPECT_THROW(basis.evaluate_hessians(xi, hessians), BasisEvaluationException);
}

TEST(BasisErrorPaths, NumericalDerivativeHelpersMatchAnalyticDerivatives) {
    ExactQuadraticBasis basis;
    const math::Vector<double, 3> xi{double(0.2), double(-0.35), double(0.4)};

    // On a quadratic, centered differences are exact except for the round-off
    // floor ~ eps_machine/step. The tolerances below are a few times
    // those floors -- tight enough that a wrong difference or analytic formula
    // (which would give an O(step) or O(1) error) cannot slip through.
    std::vector<Gradient> exact_gradients;
    basis.evaluate_gradients(xi, exact_gradients);

    std::vector<Gradient> approx_gradients;
    basis.numerical_gradient(xi, approx_gradients);
    ASSERT_EQ(approx_gradients.size(), basis.size());
    for (std::size_t n = 0; n < basis.size(); ++n) {
        for (int d = 0; d < basis.dimension(); ++d) {
            const std::size_t sd = static_cast<std::size_t>(d);
            EXPECT_NEAR(approx_gradients[n][sd], exact_gradients[n][sd], double(3e-9))
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
                            double(2e-10))
                    << "basis=" << n << " component=(" << r << "," << c << ")";
            }
        }
    }
}

TEST(BasisErrorPaths, BasisFunctionVectorOverloadsForwardToSpanPrimitives) {
    SpanPrimitiveBasis basis;
    const math::Vector<double, 3> point{double(0.25), double(0.5), double(-0.25)};

    // Reference results taken directly from the span primitives the basis defines.
    std::vector<double> span_values(basis.size());
    std::vector<Gradient> span_gradients(basis.size());
    std::vector<Hessian> span_hessians(basis.size());
    basis.evaluate_values_to(point, span_values);
    basis.evaluate_gradients_to(point, span_gradients);
    basis.evaluate_hessians_to(point, span_hessians);

    // The base-class vector overloads must size their outputs and forward to the
    // span primitives; evaluate_all() goes through the default combined evaluator.
    std::vector<double> values;
    basis.evaluate_values(point, values);
    std::vector<double> all_values;
    std::vector<Gradient> all_gradients;
    std::vector<Hessian> all_hessians;
    basis.evaluate_all(point, all_values, all_gradients, all_hessians);

    ASSERT_EQ(values.size(), basis.size());
    ASSERT_EQ(all_values.size(), basis.size());
    ASSERT_EQ(all_gradients.size(), basis.size());
    ASSERT_EQ(all_hessians.size(), basis.size());
    for (std::size_t d = 0; d < basis.size(); ++d) {
        EXPECT_EQ(values[d], span_values[d]);
        EXPECT_EQ(all_values[d], span_values[d]);
        for (std::size_t c = 0; c < 3u; ++c) {
            EXPECT_EQ(all_gradients[d][c], span_gradients[d][c]);
        }
        for (std::size_t r = 0; r < 3u; ++r) {
            for (std::size_t c = 0; c < 3u; ++c) {
                EXPECT_EQ(all_hessians[d](r, c), span_hessians[d](r, c));
            }
        }
    }
}
