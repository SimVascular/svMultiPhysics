#ifndef SVMP_FE_BASIS_LAGRANGEBASISSIMPLEX_H
#define SVMP_FE_BASIS_LAGRANGEBASISSIMPLEX_H

// Private declarations for simplex Lagrange evaluation helpers implemented in
// LagrangeBasisSimplex.cpp.

#include "BasisFunction.h"

#include <array>
#include <cstddef>
#include <vector>

namespace svmp {
namespace FE {
namespace basis {
namespace detail {

void evaluate_triangle_simplex_basis(const std::vector<std::array<int, 4>>& simplex_exponents,
                                     int order,
                                     const math::Vector<Real, 3>& xi,
                                     std::vector<Real>* values,
                                     std::vector<Gradient>* gradients,
                                     std::vector<Hessian>* hessians);

void evaluate_triangle_simplex_basis_to(const std::vector<std::array<int, 4>>& simplex_exponents,
                                        int order,
                                        const math::Vector<Real, 3>& xi,
                                        Real* SVMP_RESTRICT values_out,
                                        Real* SVMP_RESTRICT gradients_out,
                                        Real* SVMP_RESTRICT hessians_out);

void evaluate_triangle_simplex_basis_strided(
    const std::vector<std::array<int, 4>>& simplex_exponents,
    int order,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out);

void evaluate_triangle_simplex_basis_wedge_components_strided(
    const std::vector<std::array<int, 4>>& simplex_exponents,
    int order,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_xy_out,
    Real* SVMP_RESTRICT hessians_xx_xy_yy_out);

void evaluate_tetrahedron_simplex_basis(const std::vector<std::array<int, 4>>& simplex_exponents,
                                        int order,
                                        const math::Vector<Real, 3>& xi,
                                        std::vector<Real>* values,
                                        std::vector<Gradient>* gradients,
                                        std::vector<Hessian>* hessians);

void evaluate_tetrahedron_simplex_basis_to(const std::vector<std::array<int, 4>>& simplex_exponents,
                                           int order,
                                           const math::Vector<Real, 3>& xi,
                                           Real* SVMP_RESTRICT values_out,
                                           Real* SVMP_RESTRICT gradients_out,
                                           Real* SVMP_RESTRICT hessians_out);

void evaluate_tetrahedron_simplex_basis_strided(
    const std::vector<std::array<int, 4>>& simplex_exponents,
    int order,
    const std::vector<math::Vector<Real, 3>>& points,
    std::size_t output_stride,
    Real* SVMP_RESTRICT values_out,
    Real* SVMP_RESTRICT gradients_out,
    Real* SVMP_RESTRICT hessians_out);

} // namespace detail
} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_LAGRANGEBASISSIMPLEX_H
