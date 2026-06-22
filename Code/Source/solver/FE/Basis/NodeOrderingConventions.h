// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

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

/// \brief Equispaced 1D reference coordinate on [-1, 1]: -1 + 2 i / order.
///
/// Shared by the reference-node layout generators and the Lagrange tensor-axis
/// node initialization so the lattice formula lives in a single place.
[[nodiscard]] inline constexpr double line_coord_pm_one(int i, int order) noexcept {
    if (order <= 0) {
        return double(0);
    }
    return double(-1) + double(2) * static_cast<double>(i) / static_cast<double>(order);
}

class ReferenceNodeLayout {
public:
    static math::Vector<double, 3> get_node_coords(ElementType elem_type,
                                                 std::size_t local_node);
    static std::size_t num_nodes(ElementType elem_type);

    static std::vector<math::Vector<double, 3>>
    get_lagrange_node_coords(ElementType canonical_type, int order);

    static std::span<const std::size_t> mesh_to_basis_ordering(ElementType elem_type);
};

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_NODEORDERINGCONVENTIONS_H
