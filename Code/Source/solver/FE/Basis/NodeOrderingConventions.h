/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#ifndef SVMP_FE_BASIS_NODEORDERINGCONVENTIONS_H
#define SVMP_FE_BASIS_NODEORDERINGCONVENTIONS_H

#include "Math/Vector.h"
#include "Types.h"

#include <cstddef>
#include <span>
#include <vector>

namespace svmp {
namespace FE {
namespace basis {

class ReferenceNodeLayout {
public:
    static math::Vector<Real, 3> get_node_coords(ElementType elem_type,
                                                 std::size_t local_node);
    static std::size_t num_nodes(ElementType elem_type);

    static std::vector<math::Vector<Real, 3>>
    get_lagrange_node_coords(ElementType canonical_type, int order);

    static std::span<const std::size_t> mesh_to_basis_ordering(ElementType elem_type);
    static bool is_simplex(ElementType elem_type);
    static bool is_tensor_product(ElementType elem_type);
};

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_NODEORDERINGCONVENTIONS_H
