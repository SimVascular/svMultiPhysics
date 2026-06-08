/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#ifndef SVMP_FE_BASIS_LAGRANGEBASIS_H
#define SVMP_FE_BASIS_LAGRANGEBASIS_H

#include "BasisFunction.h"
#include "BasisTraits.h"

#include <array>
#include <cstddef>

namespace svmp {
namespace FE {
namespace basis {

class LagrangeBasis : public BasisFunction {
public:
    using TensorNodeIndex = std::array<std::size_t, 3>;
    using SimplexExponent = std::array<int, 4>;
    using WedgeNodeIndex = std::array<std::size_t, 2>;

    LagrangeBasis(ElementType type, int order);

    BasisType basis_type() const noexcept override { return BasisType::Lagrange; }
    ElementType element_type() const noexcept override { return element_type_; }
    int dimension() const noexcept override { return dimension_; }
    int order() const noexcept override { return order_; }
    std::size_t size() const noexcept override { return nodes_.size(); }

    const std::vector<math::Vector<Real, 3>>& nodes() const noexcept { return nodes_; }

    void evaluate_values(const math::Vector<Real, 3>& xi,
                         std::vector<Real>& values) const final;
    void evaluate_gradients(const math::Vector<Real, 3>& xi,
                            std::vector<Gradient>& gradients) const final;
    void evaluate_hessians(const math::Vector<Real, 3>& xi,
                           std::vector<Hessian>& hessians) const final;
    void evaluate_all(const math::Vector<Real, 3>& xi,
                      std::vector<Real>& values,
                      std::vector<Gradient>& gradients,
                      std::vector<Hessian>& hessians) const final;

    void evaluate_values_to(const math::Vector<Real, 3>& xi,
                            Real* SVMP_RESTRICT values_out) const final;
    void evaluate_gradients_to(const math::Vector<Real, 3>& xi,
                               Real* SVMP_RESTRICT gradients_out) const final;
    void evaluate_hessians_to(const math::Vector<Real, 3>& xi,
                              Real* SVMP_RESTRICT hessians_out) const final;

private:
    ElementType element_type_;
    BasisTopology topology_{BasisTopology::Unknown};
    int dimension_{0};
    int order_{0};

    std::vector<Real> nodes_1d_;
    std::vector<math::Vector<Real, 3>> nodes_;
    std::vector<TensorNodeIndex> tensor_indices_;
    std::vector<SimplexExponent> simplex_exponents_;
    std::vector<WedgeNodeIndex> wedge_indices_;

    void init_nodes();
    void build_point_nodes();
    void build_tensor_product_nodes(int dimensions);
    void build_simplex_nodes();
    void build_wedge_nodes();
    void init_equispaced_1d_nodes();

    void evaluate_all_to(const math::Vector<Real, 3>& xi,
                         Real* SVMP_RESTRICT values_out,
                         Real* SVMP_RESTRICT gradients_out,
                         Real* SVMP_RESTRICT hessians_out) const;
};

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_LAGRANGEBASIS_H
