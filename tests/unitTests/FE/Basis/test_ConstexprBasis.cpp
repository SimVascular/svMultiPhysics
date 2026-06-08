/**
 * @file test_ConstexprBasis.cpp
 * @brief Compile-time and lightweight runtime checks for reduced Basis helpers.
 */

#include "FE/Basis/BasisExceptions.h"
#include "FE/Basis/BasisTraits.h"
#include "FE/Basis/NodeOrderingConventions.h"

#include <gtest/gtest.h>

#include <limits>
#include <tuple>
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
static_assert(!is_pyramid(ElementType::Pyramid5));
static_assert(!is_pyramid(ElementType::Pyramid14));
static_assert(is_simplex(ElementType::Triangle3));
static_assert(is_simplex(ElementType::Tetra4));
static_assert(!is_simplex(ElementType::Wedge6));
static_assert(is_tensor_product(ElementType::Line2));
static_assert(is_tensor_product(ElementType::Quad9));
static_assert(is_tensor_product(ElementType::Hex27));
static_assert(!is_tensor_product(ElementType::Wedge6));
static_assert(topology(ElementType::Pyramid5) == BasisTopology::Unknown);
static_assert(canonical_lagrange_type(ElementType::Hex27) == ElementType::Hex8);
static_assert(canonical_lagrange_type(ElementType::Pyramid13) == ElementType::Pyramid13);
static_assert(complete_lagrange_alias_order(ElementType::Wedge18) == 2);
static_assert(complete_lagrange_alias_order(ElementType::Pyramid14) == -1);
static_assert(line_lagrange_size(2) == 3u);
static_assert(triangle_lagrange_size(2) == 6u);
static_assert(quad_lagrange_size(2) == 9u);
static_assert(tetra_lagrange_size(2) == 10u);
static_assert(hex_lagrange_size(2) == 27u);
static_assert(wedge_lagrange_size(2) == 18u);
static_assert(complete_lagrange_alias_size(ElementType::Pyramid14) == 0u);
static_assert(detail::basis_abs(Real(-2)) == Real(2));
static_assert(detail::basis_max(Real(2), Real(3)) == Real(3));
static_assert(detail::basis_near_zero(std::numeric_limits<Real>::epsilon() * Real(32)));
static_assert(detail::basis_nearly_equal(
    Real(1),
    Real(1) + std::numeric_limits<Real>::epsilon() * Real(32)));

TEST(ConstexprBasis, FixedNodeTableSizesForSupportedLayouts) {
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
    };

    for (const auto& [type, size] : expected) {
        EXPECT_EQ(ReferenceNodeLayout::num_nodes(type), size);
    }
}

TEST(ConstexprBasis, TraitToleranceScalesWithRealPrecision) {
    const Real eps = std::numeric_limits<Real>::epsilon();
    EXPECT_GT(detail::basis_scaled_tolerance(), eps);
    EXPECT_TRUE(detail::basis_near_zero(eps * Real(32)));
    EXPECT_FALSE(detail::basis_near_zero(eps * Real(128)));
    EXPECT_TRUE(detail::basis_nearly_equal(Real(1), Real(1) + eps * Real(32)));
    EXPECT_FALSE(detail::basis_nearly_equal(Real(1), Real(1) + eps * Real(128)));
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

TEST(ConstexprBasis, PyramidNodeOrderingIsOutsideCurrentScope) {
    EXPECT_THROW((void)ReferenceNodeLayout::num_nodes(ElementType::Pyramid5),
                 BasisNodeOrderingException);
    EXPECT_THROW((void)ReferenceNodeLayout::num_nodes(ElementType::Pyramid13),
                 BasisNodeOrderingException);
    EXPECT_THROW((void)ReferenceNodeLayout::get_lagrange_node_coords(ElementType::Pyramid5, 1),
                 BasisNodeOrderingException);
}

} // namespace
} // namespace basis
} // namespace FE
} // namespace svmp
