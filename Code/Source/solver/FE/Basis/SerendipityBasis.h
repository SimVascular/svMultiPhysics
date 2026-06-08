// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_BASIS_SERENDIPITYBASIS_H
#define SVMP_FE_BASIS_SERENDIPITYBASIS_H

/**
 * @file SerendipityBasis.h
 * @brief Reduced-degree-of-freedom serendipity bases
 */

#include "BasisFunction.h"

#include <array>

namespace svmp {
namespace FE {
namespace basis {

class SerendipityBasis : public BasisFunction {
public:
    SerendipityBasis(ElementType type, int order, bool geometry_mode = false);

    BasisType basis_type() const noexcept override { return BasisType::Serendipity; }
    ElementType element_type() const noexcept override { return element_type_; }
    int dimension() const noexcept override { return dimension_; }
    int order() const noexcept override { return order_; }
    std::size_t size() const noexcept override { return size_; }
    const std::vector<math::Vector<Real, 3>>& nodes() const noexcept { return nodes_; }

    void evaluate_values(const math::Vector<Real, 3>& xi,
                         std::vector<Real>& values) const override;

    void evaluate_gradients(const math::Vector<Real, 3>& xi,
                            std::vector<Gradient>& gradients) const override;

    void evaluate_hessians(const math::Vector<Real, 3>& xi,
                           std::vector<Hessian>& hessians) const override;

private:
    ElementType element_type_;
    int dimension_;
    int order_;
    std::size_t size_;
    std::vector<math::Vector<Real, 3>> nodes_;
    std::vector<std::array<int, 2>> quad_monomial_exponents_;
    // Row-major inverse Vandermonde, indexed as [monomial, basis].
    std::vector<Real> quad_inv_vandermonde_;

    // When true, this basis is used purely for geometry mapping and may use
    // reduced polynomial order (e.g., Hex20 geometry as Hex8).
    bool geometry_mode_;
};

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_SERENDIPITYBASIS_H
