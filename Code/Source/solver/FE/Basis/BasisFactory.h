// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_BASIS_BASISFACTORY_H
#define SVMP_FE_BASIS_BASISFACTORY_H

/**
 * @file BasisFactory.h
 * @brief Runtime creation of basis families
 */

#include "BasisFunction.h"
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace svmp {
namespace FE {
namespace basis {

struct BasisRequest {
    // Named mesh element layout for default/mesh-compatible bases. Leave Unknown
    // when requesting an arbitrary-order basis by reference topology.
    ElementType element_type{ElementType::Unknown};
    BasisType basis_type{BasisType::Lagrange};
    std::optional<int> order{};
    Continuity continuity{Continuity::C0};
    FieldType field_type{FieldType::Scalar};
    std::vector<double> knot_vector{};
    std::vector<double> weights{};
    std::vector<int> axis_orders{};
    std::vector<std::vector<double>> axis_knot_vectors{};
    std::vector<std::vector<double>> axis_weights{};
    std::vector<int> tensor_extents{};
    std::string custom_id{};
    // Reference topology for arbitrary-order bases. This field is intentionally
    // last so existing aggregate initializers for named elements keep their
    // positional meaning.
    BasisTopology topology{BasisTopology::Unknown};
};

namespace basis_factory {

/**
 * @brief Create a basis from a runtime request.
 *
 * @details A request must identify exactly one construction target: set
 * BasisRequest::element_type for a named mesh-node layout, or set
 * BasisRequest::topology for an arbitrary-order reference-topology basis.
 * Setting neither target, or setting both, is rejected. Named element requests
 * keep the element's fixed polynomial order contract; topology requests are the
 * arbitrary-order path.
 *
 * @param req Basis family, target, and order request.
 * @return Shared basis instance.
 */
[[nodiscard]] std::shared_ptr<BasisFunction> create(const BasisRequest& req);

/**
 * @brief Return the default basis request (family and order) for an element type.
 *
 * @details This is the single source of truth for which basis family and
 * polynomial order a given element type uses by default: serendipity node
 * layouts (Quad8, Hex20, Wedge15) select the quadratic serendipity family,
 * and every complete Lagrange element selects the Lagrange family at the
 * order given by its node layout. Solver-facing adapters should translate
 * their element names to ElementType and delegate the basis choice here
 * rather than tabulating family/order themselves.
 *
 * @param element_type Element type to select a default basis for.
 * @return Basis request suitable for create().
 * @throws BasisElementCompatibilityException If no default basis is defined
 *         for the element type.
 */
[[nodiscard]] BasisRequest default_basis_request(ElementType element_type);

/**
 * @brief Create the default basis for an element type.
 *
 * @details Equivalent to create(default_basis_request(element_type)).
 *
 * @param element_type Element type to create a default basis for.
 * @return Shared basis instance.
 */
[[nodiscard]] std::shared_ptr<BasisFunction> create_default_for(ElementType element_type);

} // namespace basis_factory

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_BASISFACTORY_H
