/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#ifndef SVMP_FE_BASIS_BASISFACTORY_H
#define SVMP_FE_BASIS_BASISFACTORY_H

/**
 * @file BasisFactory.h
 * @brief Runtime creation of basis families
 */

#include "BasisFunction.h"
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace svmp {
namespace FE {
namespace basis {

struct BasisRequest {
    ElementType element_type;
    BasisType basis_type;
    std::optional<int> order{};
    Continuity continuity{Continuity::C0};
    FieldType field_type{FieldType::Scalar};
    std::vector<Real> knot_vector{};
    std::vector<Real> weights{};
    std::vector<int> axis_orders{};
    std::vector<std::vector<Real>> axis_knot_vectors{};
    std::vector<std::vector<Real>> axis_weights{};
    std::vector<int> tensor_extents{};
    std::string custom_id{};
};

namespace basis_factory {

using CustomFactory = std::function<std::shared_ptr<BasisFunction>(const BasisRequest&)>;

[[nodiscard]] std::shared_ptr<BasisFunction> create(const BasisRequest& req);
void register_custom(std::string custom_id, CustomFactory factory);
void unregister_custom(const std::string& custom_id);
void clear_custom_registry_for_tests();

} // namespace basis_factory

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_BASISFACTORY_H
