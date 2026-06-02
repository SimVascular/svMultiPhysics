/**
 * @file test_BasisCacheFactory.cpp
 * @brief Tests for the migrated Basis cache and factory subset.
 */

#include <gtest/gtest.h>

#include "FE/Basis/BasisCache.h"
#include "FE/Basis/BasisFactory.h"
#include "FE/Basis/LagrangeBasis.h"
#include "FE/Basis/SerendipityBasis.h"
#include "FE/Quadrature/QuadratureRule.h"

#include <memory>
#include <vector>

using namespace svmp::FE;
using namespace svmp::FE::basis;
using namespace svmp::FE::quadrature;

namespace {

class CustomQuadratureRule final : public QuadratureRule {
public:
    CustomQuadratureRule(svmp::CellFamily family,
                         int dimension,
                         int order,
                         std::vector<QuadPoint> points,
                         std::vector<Real> weights)
        : QuadratureRule(family, dimension, order)
    {
        set_data(std::move(points), std::move(weights));
    }
};

CustomQuadratureRule line_rule() {
    return CustomQuadratureRule(
        svmp::CellFamily::Line, 1, 3,
        {
            QuadPoint{Real(-0.5), Real(0), Real(0)},
            QuadPoint{Real(0.5), Real(0), Real(0)}
        },
        {Real(1), Real(1)});
}

CustomQuadratureRule quad_rule(Real first_weight = Real(1)) {
    return CustomQuadratureRule(
        svmp::CellFamily::Quad, 2, 3,
        {
            QuadPoint{Real(-0.5), Real(-0.5), Real(0)},
            QuadPoint{Real(0.5), Real(-0.25), Real(0)},
            QuadPoint{Real(0.0), Real(0.5), Real(0)}
        },
        {first_weight, Real(1), Real(2)});
}

class TestCustomScalarBasis final : public BasisFunction {
public:
    explicit TestCustomScalarBasis(int tag)
        : tag_(tag)
    {
    }

    BasisType basis_type() const noexcept override { return BasisType::Custom; }
    ElementType element_type() const noexcept override { return ElementType::Line2; }
    int dimension() const noexcept override { return 1; }
    int order() const noexcept override { return 1; }
    std::size_t size() const noexcept override { return 2u; }

    std::string cache_identity() const override {
        return BasisFunction::cache_identity() + "|tag=" + std::to_string(tag_);
    }

    void evaluate_values(const math::Vector<Real, 3>& xi,
                         std::vector<Real>& values) const override
    {
        values.resize(2u);
        const Real shift = Real(tag_) * Real(0.125);
        values[0] = Real(0.5) * (Real(1) - xi[0]) + shift;
        values[1] = Real(0.5) * (Real(1) + xi[0]) - shift;
    }

    void evaluate_gradients(const math::Vector<Real, 3>&,
                            std::vector<Gradient>& gradients) const override
    {
        gradients.assign(2u, Gradient{});
        gradients[0][0] = Real(-0.5);
        gradients[1][0] = Real(0.5);
    }

private:
    int tag_{0};
};

class StructuredIdentityScalarBasis final : public BasisFunction {
public:
    explicit StructuredIdentityScalarBasis(int tag)
        : tag_(tag)
    {
    }

    BasisType basis_type() const noexcept override { return BasisType::Custom; }
    ElementType element_type() const noexcept override { return ElementType::Line2; }
    int dimension() const noexcept override { return 1; }
    int order() const noexcept override { return 1; }
    std::size_t size() const noexcept override { return 2u; }

    bool cache_identity_words(std::vector<std::uint64_t>& words) const override {
        words.push_back(0x7374727563746964ULL);
        words.push_back(static_cast<std::uint64_t>(tag_));
        return true;
    }

    std::string cache_identity() const override {
        ++string_identity_calls;
        return BasisFunction::cache_identity() + "|structured-tag=" + std::to_string(tag_);
    }

    void evaluate_values(const math::Vector<Real, 3>& xi,
                         std::vector<Real>& values) const override
    {
        values.resize(2u);
        values[0] = Real(1) - xi[0] + Real(tag_);
        values[1] = xi[0] - Real(tag_);
    }

    mutable std::size_t string_identity_calls{0};

private:
    int tag_{0};
};

} // namespace

TEST(BasisFactory, CreatesLagrangeAndSerendipityBases) {
    auto lagrange = basis_factory::create(
        BasisRequest{ElementType::Line2, BasisType::Lagrange, 2});
    ASSERT_NE(lagrange, nullptr);
    EXPECT_EQ(lagrange->basis_type(), BasisType::Lagrange);
    EXPECT_EQ(lagrange->element_type(), ElementType::Line2);
    EXPECT_EQ(lagrange->order(), 2);

    auto serendipity = basis_factory::create(
        BasisRequest{ElementType::Quad8, BasisType::Serendipity, 2});
    ASSERT_NE(serendipity, nullptr);
    EXPECT_EQ(serendipity->basis_type(), BasisType::Serendipity);
    EXPECT_EQ(serendipity->element_type(), ElementType::Quad8);
    EXPECT_EQ(serendipity->size(), 8u);
}

TEST(BasisFactory, RejectsOutOfScopeAndInvalidRequests) {
    EXPECT_THROW(
        (void)basis_factory::create(BasisRequest{ElementType::Line2, BasisType::Lagrange}),
        BasisConfigurationException);
    EXPECT_THROW(
        (void)basis_factory::create(
            BasisRequest{ElementType::Line2, BasisType::Lagrange, -1}),
        BasisConfigurationException);
    EXPECT_THROW(
        (void)basis_factory::create(
            BasisRequest{ElementType::Line2, BasisType::Bernstein, 1}),
        BasisConfigurationException);
    EXPECT_THROW(
        (void)basis_factory::create(
            BasisRequest{ElementType::Line2,
                         BasisType::Lagrange,
                         1,
                         Continuity::H_div,
                         FieldType::Vector}),
        BasisConfigurationException);
}

TEST(BasisFactory, SupportsCustomFactoryRegistration) {
    basis_factory::clear_custom_registry_for_tests();
    basis_factory::register_custom(
        "test-custom",
        [](const BasisRequest& req) {
            const int tag = req.order.value_or(0);
            return std::make_shared<TestCustomScalarBasis>(tag);
        });

    BasisRequest req{ElementType::Line2, BasisType::Custom, 7};
    req.custom_id = "test-custom";
    auto custom = basis_factory::create(req);
    ASSERT_NE(custom, nullptr);
    EXPECT_EQ(custom->basis_type(), BasisType::Custom);
    EXPECT_EQ(custom->size(), 2u);

    basis_factory::unregister_custom("test-custom");
    EXPECT_THROW((void)basis_factory::create(req), BasisConfigurationException);
    basis_factory::clear_custom_registry_for_tests();
}

TEST(BasisCache, ReusesEntriesForSameBasisAndQuadratureCoordinates) {
    LagrangeBasis basis(ElementType::Line2, 2);
    const auto quad = line_rule();

    auto& cache = BasisCache::instance();
    cache.clear();
    const auto& entry1 = cache.get_or_compute(basis, quad, true, true);
    const auto& entry2 = cache.get_or_compute(basis, quad, true, true);

    EXPECT_EQ(&entry1, &entry2);
    EXPECT_EQ(entry1.num_qpts, quad.num_points());
    EXPECT_EQ(entry1.num_dofs, basis.size());
    ASSERT_EQ(entry1.scalar_values.size(), basis.size() * quad.num_points());
    ASSERT_EQ(entry1.gradients.size(), basis.size() * 3u * quad.num_points());
    ASSERT_EQ(entry1.hessians.size(), basis.size() * 9u * quad.num_points());
    EXPECT_EQ(cache.size(), 1u);
}

TEST(BasisCache, ReusesCoordinateIdenticalQuadratureRulesIgnoringWeights) {
    SerendipityBasis basis(ElementType::Quad8, 2);
    const auto quad_a = quad_rule(Real(1));
    const auto quad_b = quad_rule(Real(0.25));

    auto& cache = BasisCache::instance();
    cache.clear();
    const auto& entry_a = cache.get_or_compute(basis, quad_a, true, false);
    const auto& entry_b = cache.get_or_compute(basis, quad_b, true, false);

    EXPECT_EQ(&entry_a, &entry_b);
    EXPECT_EQ(cache.size(), 1u);
}

TEST(BasisCache, SeparatesStringIdentityCustomBases) {
    TestCustomScalarBasis custom_a(1);
    TestCustomScalarBasis custom_b(2);
    const auto quad = line_rule();

    auto& cache = BasisCache::instance();
    cache.clear();
    const auto& entry_a = cache.get_or_compute(custom_a, quad, false, false);
    const auto& entry_b = cache.get_or_compute(custom_b, quad, false, false);

    EXPECT_NE(&entry_a, &entry_b);
    EXPECT_NE(entry_a.scalar_values, entry_b.scalar_values);
    EXPECT_EQ(cache.size(), 2u);
}

TEST(BasisCache, StructuredIdentityAvoidsStringFallbackAndSeparatesBases) {
    StructuredIdentityScalarBasis custom_a(1);
    StructuredIdentityScalarBasis custom_b(2);
    const auto quad = line_rule();

    auto& cache = BasisCache::instance();
    cache.clear();
    const auto& entry_a = cache.get_or_compute(custom_a, quad, false, false);
    const auto& entry_b = cache.get_or_compute(custom_b, quad, false, false);

    EXPECT_NE(&entry_a, &entry_b);
    EXPECT_EQ(custom_a.string_identity_calls, 0u);
    EXPECT_EQ(custom_b.string_identity_calls, 0u);
    EXPECT_EQ(cache.size(), 2u);
}

