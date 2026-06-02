/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#include "BasisFactory.h"

#include "LagrangeBasis.h"
#include "SerendipityBasis.h"

#include <mutex>
#include <unordered_map>
#include <utility>

namespace svmp {
namespace FE {
namespace basis {

namespace {

using CustomRegistryMap =
    std::unordered_map<std::string, basis_factory::CustomFactory>;

CustomRegistryMap& custom_registry() {
    static CustomRegistryMap registry;
    return registry;
}

std::mutex& custom_registry_mutex() {
    static std::mutex mutex;
    return mutex;
}

int require_basis_order(const BasisRequest& req,
                        const char* missing_message,
                        const char* negative_message) {
    if (!req.order.has_value()) {
        throw BasisConfigurationException(missing_message,
                                          __FILE__, __LINE__, __func__);
    }
    if (*req.order < 0) {
        throw BasisConfigurationException(negative_message,
                                          __FILE__, __LINE__, __func__);
    }
    return *req.order;
}

void require_scalar_c0_request(const BasisRequest& req) {
    if (req.field_type != FieldType::Scalar) {
        throw BasisConfigurationException(
            "BasisFactory: Lagrange/Serendipity bases currently support scalar fields only",
            __FILE__, __LINE__, __func__);
    }
    if (req.continuity != Continuity::C0) {
        throw BasisConfigurationException(
            "BasisFactory: migrated Lagrange/Serendipity scope supports C0 continuity only",
            __FILE__, __LINE__, __func__);
    }
}

std::shared_ptr<BasisFunction> create_lagrange(const BasisRequest& req) {
    require_scalar_c0_request(req);
    const int order = require_basis_order(
        req,
        "BasisFactory: Lagrange creation requires an explicit order",
        "BasisFactory: Lagrange requires non-negative order");
    return std::make_shared<LagrangeBasis>(req.element_type, order);
}

std::shared_ptr<BasisFunction> create_serendipity(const BasisRequest& req) {
    require_scalar_c0_request(req);
    const int order = require_basis_order(
        req,
        "BasisFactory: Serendipity creation requires an explicit order",
        "BasisFactory: Serendipity requires non-negative order");
    return std::make_shared<SerendipityBasis>(req.element_type, order);
}

std::shared_ptr<BasisFunction> create_custom(const BasisRequest& req) {
    if (req.custom_id.empty()) {
        throw BasisConfigurationException(
            "BasisFactory: custom basis requests require custom_id",
            __FILE__, __LINE__, __func__);
    }

    basis_factory::CustomFactory factory;
    {
        std::lock_guard<std::mutex> lock(custom_registry_mutex());
        const auto it = custom_registry().find(req.custom_id);
        if (it == custom_registry().end()) {
            throw BasisConfigurationException(
                "BasisFactory: no custom basis factory registered for id '" +
                    req.custom_id + "'",
                __FILE__, __LINE__, __func__);
        }
        factory = it->second;
    }

    auto basis = factory(req);
    if (!basis) {
        throw BasisConstructionException(
            "BasisFactory: custom factory returned null basis for id '" +
                req.custom_id + "'",
            __FILE__, __LINE__, __func__);
    }
    return basis;
}

} // namespace

namespace basis_factory {

std::shared_ptr<BasisFunction> create(const BasisRequest& req) {
    switch (req.basis_type) {
        case BasisType::Lagrange:
            return create_lagrange(req);
        case BasisType::Serendipity:
            return create_serendipity(req);
        case BasisType::Custom:
            return create_custom(req);
        default:
            throw BasisConfigurationException(
                "BasisFactory: requested basis family is outside the migrated Lagrange/Serendipity scope",
                __FILE__, __LINE__, __func__);
    }
}

void register_custom(std::string custom_id, CustomFactory factory) {
    if (custom_id.empty()) {
        throw BasisConfigurationException(
            "BasisFactory: custom factory id must not be empty",
            __FILE__, __LINE__, __func__);
    }
    if (!factory) {
        throw BasisConfigurationException(
            "BasisFactory: custom factory must be callable",
            __FILE__, __LINE__, __func__);
    }

    std::lock_guard<std::mutex> lock(custom_registry_mutex());
    custom_registry()[std::move(custom_id)] = std::move(factory);
}

void unregister_custom(const std::string& custom_id) {
    std::lock_guard<std::mutex> lock(custom_registry_mutex());
    custom_registry().erase(custom_id);
}

void clear_custom_registry_for_tests() {
    std::lock_guard<std::mutex> lock(custom_registry_mutex());
    custom_registry().clear();
}

} // namespace basis_factory

} // namespace basis
} // namespace FE
} // namespace svmp
