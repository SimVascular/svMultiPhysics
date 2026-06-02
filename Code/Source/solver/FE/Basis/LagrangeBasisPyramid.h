#ifndef SVMP_FE_BASIS_LAGRANGEBASISPYRAMID_H
#define SVMP_FE_BASIS_LAGRANGEBASISPYRAMID_H

// Private declarations for the rational pyramid Lagrange helper implemented in
// LagrangeBasisPyramid.cpp. This header is intentionally small so the large
// construction and apex-classification code stays out of LagrangeBasis.cpp.

#include "BasisFunction.h"

#include <cstddef>
#include <vector>

namespace svmp {
namespace FE {
namespace basis {
namespace detail {
namespace lagrange_pyramid {

const std::vector<math::Vector<Real, 3>>& nodes(int order);

void prewarm_scratch(int order, std::size_t max_qpts = 0);

void evaluate_values(int order,
                     const math::Vector<Real, 3>& xi,
                     std::vector<Real>& values);
void evaluate_gradients(int order,
                        const math::Vector<Real, 3>& xi,
                        std::vector<Gradient>& gradients);
void evaluate_hessians(int order,
                       const math::Vector<Real, 3>& xi,
                       std::vector<Hessian>& hessians);
void evaluate_all(int order,
                  const math::Vector<Real, 3>& xi,
                  std::vector<Real>& values,
                  std::vector<Gradient>& gradients,
                  std::vector<Hessian>& hessians);

void evaluate_values_to(int order,
                        const math::Vector<Real, 3>& xi,
                        Real* SVMP_RESTRICT values_out);
void evaluate_gradients_to(int order,
                           const math::Vector<Real, 3>& xi,
                           Real* SVMP_RESTRICT gradients_out);
void evaluate_hessians_to(int order,
                          const math::Vector<Real, 3>& xi,
                          Real* SVMP_RESTRICT hessians_out);
void evaluate_all_to(int order,
                     const math::Vector<Real, 3>& xi,
                     Real* SVMP_RESTRICT values_out,
                     Real* SVMP_RESTRICT gradients_out,
                     Real* SVMP_RESTRICT hessians_out);

void evaluate_at_quadrature_points_strided(
    int order,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out);

} // namespace lagrange_pyramid
} // namespace detail
} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_LAGRANGEBASISPYRAMID_H
