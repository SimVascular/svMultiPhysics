// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "BasisFactory.h"

#include "BasisTraits.h"
#include "LagrangeBasis.h"
#include "SerendipityBasis.h"

namespace svmp {
namespace FE {
namespace basis {

namespace {

enum class RequestTarget {
    NamedElement,
    Topology,
};

int require_basis_order(const BasisRequest& req,
                        const char* missing_message,
                        const char* negative_message) {
    svmp::throw_if<BasisConfigurationException>(!req.order.has_value(), SVMP_HERE,
                                              missing_message);
    svmp::throw_if<BasisConfigurationException>(*req.order < 0, SVMP_HERE,
                                              negative_message);
    return *req.order;
}

RequestTarget require_single_request_target(const BasisRequest& req) {
    const bool has_named_element = req.element_type != ElementType::Unknown;
    const bool has_topology = req.topology != BasisTopology::Unknown;
    svmp::throw_if<BasisConfigurationException>(
        !has_named_element && !has_topology, SVMP_HERE,
        "BasisFactory: request must specify either a named element_type or a reference topology");
    svmp::throw_if<BasisConfigurationException>(
        has_named_element && has_topology, SVMP_HERE,
        "BasisFactory: request must specify element_type or topology, not both");
    return has_topology ? RequestTarget::Topology : RequestTarget::NamedElement;
}

void require_scalar_c0_request(const BasisRequest& req) {
    svmp::throw_if<BasisConfigurationException>(
        req.field_type != FieldType::Scalar, SVMP_HERE,
        "BasisFactory: Lagrange/Serendipity bases support scalar fields only");
    svmp::throw_if<BasisConfigurationException>(
        req.continuity != Continuity::C0, SVMP_HERE,
        "BasisFactory: Lagrange/Serendipity bases support C0 continuity only");
}

std::shared_ptr<BasisFunction> create_lagrange(const BasisRequest& req) {
    require_scalar_c0_request(req);
    const int order = require_basis_order(
        req,
        "BasisFactory: Lagrange creation requires an explicit order",
        "BasisFactory: Lagrange requires non-negative order");
    if (require_single_request_target(req) == RequestTarget::Topology) {
        return std::make_shared<LagrangeBasis>(req.topology, order);
    }
    return std::make_shared<LagrangeBasis>(req.element_type, order);
}

std::shared_ptr<BasisFunction> create_serendipity(const BasisRequest& req) {
    require_scalar_c0_request(req);
    const int order = require_basis_order(
        req,
        "BasisFactory: Serendipity creation requires an explicit order",
        "BasisFactory: Serendipity requires non-negative order");
    if (require_single_request_target(req) == RequestTarget::Topology) {
        return std::make_shared<SerendipityBasis>(req.topology, order);
    }
    return std::make_shared<SerendipityBasis>(req.element_type, order);
}

} // namespace

namespace basis_factory {

std::shared_ptr<BasisFunction> create(const BasisRequest& req) {
    switch (req.basis_type) {
        case BasisType::Lagrange:
            return create_lagrange(req);
        case BasisType::Serendipity:
            return create_serendipity(req);
        default:
            svmp::raise<BasisConfigurationException>(SVMP_HERE,
                "BasisFactory: requested basis family is outside the scalar Lagrange/Serendipity scope");
    }
}

BasisRequest default_basis_request(ElementType element_type) {
    switch (element_type) {
        // Reduced serendipity node layouts have no complete Lagrange basis at
        // their node count; they always use the quadratic serendipity space.
        case ElementType::Quad8:
        case ElementType::Hex20:
        case ElementType::Wedge15:
            return BasisRequest{element_type, BasisType::Serendipity, 2};
        case ElementType::Point1:
            return BasisRequest{element_type, BasisType::Lagrange, 0};
        default: {
            const int order = complete_lagrange_alias_order(element_type);
            if (order >= 0) {
                return BasisRequest{element_type, BasisType::Lagrange, order};
            }
            svmp::raise<BasisElementCompatibilityException>(SVMP_HERE,
                "BasisFactory: no default basis is defined for the requested element type");
        }
    }
}

std::shared_ptr<BasisFunction> create_default_for(ElementType element_type) {
    return create(default_basis_request(element_type));
}

} // namespace basis_factory

} // namespace basis
} // namespace FE
} // namespace svmp
