/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#ifndef SVMP_FE_BASIS_BASISFUNCTION_H
#define SVMP_FE_BASIS_BASISFUNCTION_H

#include "BasisExceptions.h"
#include "Math/Matrix.h"
#include "Math/Vector.h"
#include "Types.h"

#include <cstddef>
#include <vector>

namespace svmp {
namespace FE {
namespace basis {

using Gradient = math::Vector<Real, 3>;
using Hessian  = math::Matrix<Real, 3, 3>;

void prewarm_basis_function_scratch(std::size_t max_size);

[[nodiscard]] inline Hessian make_symmetric_hessian(Real xx,
                                                    Real yy,
                                                    Real zz,
                                                    Real xy,
                                                    Real xz,
                                                    Real yz) {
    Hessian hessian{};
    hessian(0, 0) = xx;
    hessian(1, 1) = yy;
    hessian(2, 2) = zz;
    hessian(0, 1) = hessian(1, 0) = xy;
    hessian(0, 2) = hessian(2, 0) = xz;
    hessian(1, 2) = hessian(2, 1) = yz;
    return hessian;
}

inline void store_hessian(const Hessian& hessian, Real* dst) noexcept {
    dst[0] = hessian(0, 0);
    dst[1] = hessian(0, 1);
    dst[2] = hessian(0, 2);
    dst[3] = hessian(1, 0);
    dst[4] = hessian(1, 1);
    dst[5] = hessian(1, 2);
    dst[6] = hessian(2, 0);
    dst[7] = hessian(2, 1);
    dst[8] = hessian(2, 2);
}

[[nodiscard]] inline Hessian load_hessian(const Real* src) noexcept {
    Hessian hessian{};
    hessian(0, 0) = src[0];
    hessian(0, 1) = src[1];
    hessian(0, 2) = src[2];
    hessian(1, 0) = src[3];
    hessian(1, 1) = src[4];
    hessian(1, 2) = src[5];
    hessian(2, 0) = src[6];
    hessian(2, 1) = src[7];
    hessian(2, 2) = src[8];
    return hessian;
}

inline void add_scaled_hessian(Hessian& target,
                               const Hessian& source,
                               Real scale) noexcept {
    for (std::size_t r = 0; r < 3u; ++r) {
        for (std::size_t c = 0; c < 3u; ++c) {
            target(r, c) += scale * source(r, c);
        }
    }
}

class BasisFunction {
public:
    virtual ~BasisFunction() = default;

    virtual BasisType basis_type() const noexcept = 0;
    virtual ElementType element_type() const noexcept = 0;
    virtual int dimension() const noexcept = 0;
    virtual int order() const noexcept = 0;
    virtual std::size_t size() const noexcept = 0;

    virtual void evaluate_values(const math::Vector<Real, 3>& xi,
                                 std::vector<Real>& values) const = 0;
    virtual void evaluate_gradients(const math::Vector<Real, 3>& xi,
                                    std::vector<Gradient>& gradients) const;
    virtual void evaluate_hessians(const math::Vector<Real, 3>& xi,
                                   std::vector<Hessian>& hessians) const;
    virtual void evaluate_all(const math::Vector<Real, 3>& xi,
                              std::vector<Real>& values,
                              std::vector<Gradient>& gradients,
                              std::vector<Hessian>& hessians) const;

    virtual void evaluate_values_to(const math::Vector<Real, 3>& xi,
                                    Real* SVMP_RESTRICT values_out) const;
    virtual void evaluate_gradients_to(const math::Vector<Real, 3>& xi,
                                       Real* SVMP_RESTRICT gradients_out) const;
    virtual void evaluate_hessians_to(const math::Vector<Real, 3>& xi,
                                      Real* SVMP_RESTRICT hessians_out) const;

protected:
    void numerical_gradient(const math::Vector<Real, 3>& xi,
                            std::vector<Gradient>& gradients,
                            Real eps = Real(1e-6)) const;
    void numerical_hessian(const math::Vector<Real, 3>& xi,
                           std::vector<Hessian>& hessians,
                           Real eps = Real(1e-5)) const;
};

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_BASISFUNCTION_H
