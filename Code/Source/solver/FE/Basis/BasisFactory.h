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

namespace svmp::FE::basis {

/**
 * @brief Runtime description of a basis to construct.
 * @ingroup FE_Basis
 *
 * @details A request identifies exactly one construction target -- a named
 * ElementType layout, or a reference BasisTopology with an explicit order -- plus
 * the family and field policy; basis_factory::create() validates and builds from
 * it. The spline/NURBS fields are reserved for future families and are unused by
 * the scalar Lagrange/Serendipity factory.
 */
struct BasisRequest {
    ElementType element_type{ElementType::Unknown};  ///< Named element layout, or Unknown to request by topology.
    BasisType basis_type{BasisType::Lagrange};       ///< Basis family to construct.
    std::optional<int> order{};                      ///< Polynomial order; required by the factory.
    Continuity continuity{Continuity::C0};           ///< Inter-element continuity (Lagrange/Serendipity are C0).
    FieldType field_type{FieldType::Scalar};         ///< Field type (Lagrange/Serendipity support Scalar).
    std::vector<double> knot_vector{};               ///< Reserved for spline/NURBS families; unused here.
    std::vector<double> weights{};                   ///< Reserved for rational (NURBS) families; unused here.
    std::vector<int> axis_orders{};                  ///< Reserved for per-axis tensor spline orders; unused here.
    std::vector<std::vector<double>> axis_knot_vectors{};  ///< Reserved for per-axis spline knots; unused here.
    std::vector<std::vector<double>> axis_weights{};       ///< Reserved for per-axis rational weights; unused here.
    std::vector<int> tensor_extents{};               ///< Reserved for tensor-product extents; unused here.
    std::string custom_id{};                         ///< Optional identifier for Custom families.
    /// Reference topology for arbitrary-order requests, or Unknown to request by
    /// element_type.
    BasisTopology topology{BasisTopology::Unknown};
    // Implementation note (kept out of the rendered docs): topology is declared
    // last so existing aggregate initializers for named elements keep their
    // positional meaning.
};

namespace basis_factory {

/**
 * @brief Create a basis from a runtime request.
 * @ingroup FE_Basis
 *
 * @details A request must identify exactly one construction target: set
 * BasisRequest::element_type for a named mesh-node layout, or set
 * BasisRequest::topology for an arbitrary-order reference-topology basis.
 * Setting neither target, or setting both, is rejected. Named element requests
 * keep the element's fixed polynomial order contract; topology requests are the
 * arbitrary-order path.
 *
 * @param req Basis family, target, and order request.
 * @return Unique basis instance. Move it into a std::shared_ptr at the call site
 *         if shared ownership is needed.
 */
[[nodiscard]] std::unique_ptr<BasisFunction> create(const BasisRequest& req);

/**
 * @brief Return the default basis request (family and order) for an element type.
 * @ingroup FE_Basis
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
 * @ingroup FE_Basis
 *
 * @details Equivalent to create(default_basis_request(element_type)).
 *
 * @param element_type Element type to create a default basis for.
 * @return Unique basis instance. Move it into a std::shared_ptr at the call site
 *         if shared ownership is needed.
 */
[[nodiscard]] std::unique_ptr<BasisFunction> create_default_for(ElementType element_type);

} // namespace basis_factory

} // namespace svmp::FE::basis

#endif // SVMP_FE_BASIS_BASISFACTORY_H
