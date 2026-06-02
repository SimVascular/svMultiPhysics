/**
 * @file test_ConstexprBasis.cpp
 * @brief Compile-time and lightweight runtime checks for migrated Basis helpers.
 */

#include "FE/Basis/BasisTolerance.h"
#include "FE/Basis/BasisTraits.h"
#include "FE/Basis/LagrangeBasis.h"
#include "FE/Basis/LagrangeBasisFast.h"
#include "FE/Basis/NodeOrderingConventions.h"

#include <gtest/gtest.h>

#include <array>
#include <limits>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace svmp {
namespace FE {
namespace basis {
namespace {

static_assert(is_line(ElementType::Line2));
static_assert(is_line(ElementType::Line3));
static_assert(is_triangle(ElementType::Triangle6));
static_assert(is_quadrilateral(ElementType::Quad8));
static_assert(is_tetrahedron(ElementType::Tetra10));
static_assert(is_hexahedron(ElementType::Hex20));
static_assert(is_wedge(ElementType::Wedge18));
static_assert(is_pyramid(ElementType::Pyramid14));
static_assert(is_simplex(ElementType::Triangle3));
static_assert(is_simplex(ElementType::Tetra4));
static_assert(!is_simplex(ElementType::Wedge6));
static_assert(is_tensor_product(ElementType::Line2));
static_assert(is_tensor_product(ElementType::Quad9));
static_assert(is_tensor_product(ElementType::Hex27));
static_assert(!is_tensor_product(ElementType::Pyramid5));
static_assert(reference_dimension(ElementType::Pyramid14) == 3);
static_assert(canonical_lagrange_type(ElementType::Hex27) == ElementType::Hex8);
static_assert(canonical_lagrange_type(ElementType::Pyramid13) == ElementType::Pyramid13);
static_assert(complete_lagrange_alias_order(ElementType::Wedge18) == 2);
static_assert(complete_lagrange_alias_order(ElementType::Hex20) == -1);
static_assert(line_lagrange_size(2) == 3u);
static_assert(triangle_lagrange_size(2) == 6u);
static_assert(quad_lagrange_size(2) == 9u);
static_assert(tetra_lagrange_size(2) == 10u);
static_assert(hex_lagrange_size(2) == 27u);
static_assert(wedge_lagrange_size(2) == 18u);
static_assert(pyramid_lagrange_size(2) == 14u);
static_assert(detail::basis_abs(Real(-2)) == Real(2));
static_assert(detail::basis_max(Real(2), Real(3)) == Real(3));
static_assert(detail::basis_near_zero(std::numeric_limits<Real>::epsilon() * Real(32)));
static_assert(detail::basis_nearly_equal(
    Real(1),
    Real(1) + std::numeric_limits<Real>::epsilon() * Real(32)));

constexpr auto kLineFastValues = [] {
    math::Vector<Real, 3> xi{Real(0), Real(0), Real(0)};
    std::array<Real, LagrangeLineFast<1>::n_dofs> values{};
    LagrangeLineFast<1>::evaluate(xi, values);
    return values;
}();
static_assert(kLineFastValues[0] == Real(0.5));
static_assert(kLineFastValues[1] == Real(0.5));

constexpr auto kLineP2FastHessians = [] {
    math::Vector<Real, 3> xi{Real(0), Real(0), Real(0)};
    std::array<Hessian, LagrangeLineFast<2>::n_dofs> hessians{};
    LagrangeLineFast<2>::evaluate_hessians(xi, hessians);
    return hessians;
}();
static_assert(kLineP2FastHessians[0](0, 0) == Real(1));
static_assert(kLineP2FastHessians[1](0, 0) == Real(1));
static_assert(kLineP2FastHessians[2](0, 0) == Real(-2));

constexpr auto kTriP2FastValues = [] {
    math::Vector<Real, 3> xi{Real(0.25), Real(0.25), Real(0)};
    std::array<Real, LagrangeTriFast<2>::n_dofs> values{};
    LagrangeTriFast<2>::evaluate(xi, values);
    return values;
}();
static_assert(kTriP2FastValues[0] == Real(0));
static_assert(kTriP2FastValues[3] == Real(0.5));
static_assert(kTriP2FastValues[4] == Real(0.25));

template<typename Basis>
constexpr bool overrides_scalar_strided_v =
    !std::is_same_v<decltype(&Basis::evaluate_at_quadrature_points_strided),
                    decltype(&BasisFunction::evaluate_at_quadrature_points_strided)>;

template<typename FastBasis>
void expect_fast_matches_lagrange(ElementType type,
                                  int order,
                                  const std::vector<math::Vector<Real, 3>>& points)
{
    LagrangeBasis basis(type, order);
    for (const auto& xi : points) {
        std::vector<Real> expected_values;
        std::vector<Gradient> expected_gradients;
        std::vector<Hessian> expected_hessians;
        basis.evaluate_all(xi, expected_values, expected_gradients, expected_hessians);

        std::array<Real, FastBasis::n_dofs> values{};
        std::array<Gradient, FastBasis::n_dofs> gradients{};
        std::array<Hessian, FastBasis::n_dofs> hessians{};
        FastBasis::evaluate(xi, values);
        FastBasis::evaluate_gradients(xi, gradients);
        FastBasis::evaluate_hessians(xi, hessians);

        ASSERT_EQ(expected_values.size(), values.size());
        for (std::size_t i = 0; i < values.size(); ++i) {
            EXPECT_NEAR(values[i], expected_values[i], Real(1e-14));
            for (std::size_t d = 0; d < 3u; ++d) {
                EXPECT_NEAR(gradients[i][d], expected_gradients[i][d], Real(1e-14));
                for (std::size_t e = 0; e < 3u; ++e) {
                    EXPECT_NEAR(hessians[i](d, e), expected_hessians[i](d, e), Real(1e-14));
                }
            }
        }
    }
}

TEST(ConstexprBasis, FixedNodeTableSizes) {
    const std::vector<std::pair<ElementType, std::size_t>> expected = {
        {ElementType::Line2, 2u},
        {ElementType::Line3, 3u},
        {ElementType::Triangle3, 3u},
        {ElementType::Triangle6, 6u},
        {ElementType::Quad4, 4u},
        {ElementType::Quad8, 8u},
        {ElementType::Quad9, 9u},
        {ElementType::Tetra4, 4u},
        {ElementType::Tetra10, 10u},
        {ElementType::Hex8, 8u},
        {ElementType::Hex20, 20u},
        {ElementType::Hex27, 27u},
        {ElementType::Wedge6, 6u},
        {ElementType::Wedge15, 15u},
        {ElementType::Wedge18, 18u},
        {ElementType::Pyramid5, 5u},
        {ElementType::Pyramid13, 13u},
        {ElementType::Pyramid14, 14u},
    };

    for (const auto& [type, size] : expected) {
        EXPECT_EQ(ReferenceNodeLayout::num_nodes(type), size);
    }
}

TEST(ConstexprBasis, BasisToleranceScalesWithRealPrecision) {
    const Real eps = std::numeric_limits<Real>::epsilon();
    EXPECT_GT(detail::basis_scaled_tolerance(), eps);
    EXPECT_TRUE(detail::basis_near_zero(eps * Real(32)));
    EXPECT_FALSE(detail::basis_near_zero(eps * Real(128)));
    EXPECT_TRUE(detail::basis_nearly_equal(Real(1), Real(1) + eps * Real(32)));
    EXPECT_FALSE(detail::basis_nearly_equal(Real(1), Real(1) + eps * Real(128)));
}

TEST(ConstexprBasis, LagrangeOverridesStridedEvaluation) {
    EXPECT_TRUE(overrides_scalar_strided_v<LagrangeBasis>);
}

TEST(ConstexprBasis, FastSidecarsMatchRuntimeLagrangeBasis) {
    expect_fast_matches_lagrange<LagrangeLineFast<1>>(
        ElementType::Line2, 1,
        {{Real(-0.2), Real(0), Real(0)}, {Real(0.35), Real(0), Real(0)}});
    expect_fast_matches_lagrange<LagrangeLineFast<2>>(
        ElementType::Line2, 2,
        {{Real(-0.2), Real(0), Real(0)}, {Real(0.35), Real(0), Real(0)}});
    expect_fast_matches_lagrange<LagrangeQuadFast<1>>(
        ElementType::Quad4, 1,
        {{Real(-0.2), Real(0.3), Real(0)}, {Real(0.35), Real(-0.45), Real(0)}});
    expect_fast_matches_lagrange<LagrangeHexFast<1>>(
        ElementType::Hex8, 1,
        {{Real(-0.2), Real(0.3), Real(0.1)}, {Real(0.35), Real(-0.45), Real(0.25)}});
    expect_fast_matches_lagrange<LagrangeTriFast<1>>(
        ElementType::Triangle3, 1,
        {{Real(0.2), Real(0.3), Real(0)}, {Real(0.1), Real(0.6), Real(0)}});
    expect_fast_matches_lagrange<LagrangeTriFast<2>>(
        ElementType::Triangle3, 2,
        {{Real(0.2), Real(0.3), Real(0)}, {Real(0.1), Real(0.6), Real(0)}});
    expect_fast_matches_lagrange<LagrangeTetFast<1>>(
        ElementType::Tetra4, 1,
        {{Real(0.2), Real(0.3), Real(0.1)}, {Real(0.1), Real(0.2), Real(0.4)}});
    expect_fast_matches_lagrange<LagrangeTetFast<2>>(
        ElementType::Tetra4, 2,
        {{Real(0.2), Real(0.3), Real(0.1)}, {Real(0.1), Real(0.2), Real(0.4)}});
}

TEST(ConstexprBasis, CompleteAliasTablesMatchGeneratedLagrangeNodes) {
    const std::vector<std::tuple<ElementType, ElementType, int>> aliases = {
        {ElementType::Line2, ElementType::Line2, 1},
        {ElementType::Line3, ElementType::Line2, 2},
        {ElementType::Triangle3, ElementType::Triangle3, 1},
        {ElementType::Triangle6, ElementType::Triangle3, 2},
        {ElementType::Quad4, ElementType::Quad4, 1},
        {ElementType::Quad9, ElementType::Quad4, 2},
        {ElementType::Tetra4, ElementType::Tetra4, 1},
        {ElementType::Tetra10, ElementType::Tetra4, 2},
        {ElementType::Hex8, ElementType::Hex8, 1},
        {ElementType::Hex27, ElementType::Hex8, 2},
        {ElementType::Wedge6, ElementType::Wedge6, 1},
        {ElementType::Wedge18, ElementType::Wedge6, 2},
        {ElementType::Pyramid5, ElementType::Pyramid5, 1},
        {ElementType::Pyramid14, ElementType::Pyramid5, 2},
    };

    for (const auto& [alias, canonical_type, order] : aliases) {
        const auto nodes = ReferenceNodeLayout::get_lagrange_node_coords(canonical_type, order);
        ASSERT_EQ(nodes.size(), ReferenceNodeLayout::num_nodes(alias));
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            const auto direct = ReferenceNodeLayout::get_node_coords(alias, i);
            EXPECT_EQ(nodes[i][0], direct[0]);
            EXPECT_EQ(nodes[i][1], direct[1]);
            EXPECT_EQ(nodes[i][2], direct[2]);
        }
    }
}

} // namespace
} // namespace basis
} // namespace FE
} // namespace svmp
