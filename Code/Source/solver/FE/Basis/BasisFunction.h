// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_BASIS_BASISFUNCTION_H
#define SVMP_FE_BASIS_BASISFUNCTION_H

#include "BasisExceptions.h"
#include "Math/Matrix.h"
#include "Math/Vector.h"
#include "Types.h"

#include <cstddef>
#include <vector>

/// \defgroup FE FE Library
/// \brief Finite-element interfaces and utilities used by the solver.
///
/// The FE library groups basis functions, math utilities, assembly interfaces,
/// and related support code that can be built and consumed as a coherent
/// finite-element component.

/// \defgroup FE_Basis Basis
/// \ingroup FE
/// \brief Basis-function interfaces, concrete basis families, and reference-node conventions.

namespace svmp {
namespace FE {
namespace basis {

/// \brief Gradient vector type used by basis evaluators.
using Gradient = math::Vector<Real, 3>;

/// \brief Hessian matrix type used by basis evaluators.
using Hessian  = math::Matrix<Real, 3, 3>;

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

/// \brief Abstract interface for finite-element basis-function families.
/// \ingroup FE_Basis
///
/// BasisFunction defines the common query and evaluation API used by solver
/// code that does not need to know the concrete basis implementation. Derived
/// classes provide values at minimum and can override analytical gradients,
/// Hessians, combined evaluation, and flat-buffer output paths.
class BasisFunction {
public:
    /// \brief Destroy a basis function through the abstract interface.
    virtual ~BasisFunction() = default;

    /// \brief Return the concrete basis family.
    /// \return Basis family identifier.
    virtual BasisType basis_type() const noexcept = 0;

    /// \brief Return the canonical element type represented by this basis.
    /// \return Element type used for node layout and evaluation.
    virtual ElementType element_type() const noexcept = 0;

    /// \brief Return the reference-space dimension of the basis.
    /// \return Reference dimension, from zero for points through three for volume elements.
    virtual int dimension() const noexcept = 0;

    /// \brief Return the polynomial order represented by this basis.
    /// \return Effective polynomial order after any element-family normalization.
    virtual int order() const noexcept = 0;

    /// \brief Return the number of basis functions and reference nodes.
    /// \return Basis function count.
    virtual std::size_t size() const noexcept = 0;

    /// \brief Evaluate basis function values at a reference coordinate.
    /// \param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
    /// \param values Receives one value per basis function.
    virtual void evaluate_values(const math::Vector<Real, 3>& xi,
                                 std::vector<Real>& values) const = 0;

    /// \brief Evaluate basis gradients at a reference coordinate.
    /// \param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
    /// \param gradients Receives one three-component gradient per basis function.
    /// \throws BasisEvaluationException If gradients are not available for the basis.
    virtual void evaluate_gradients(const math::Vector<Real, 3>& xi,
                                    std::vector<Gradient>& gradients) const;

    /// \brief Evaluate basis Hessians at a reference coordinate.
    /// \param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
    /// \param hessians Receives one 3-by-3 Hessian per basis function.
    /// \throws BasisEvaluationException If Hessians are not available for the basis.
    virtual void evaluate_hessians(const math::Vector<Real, 3>& xi,
                                   std::vector<Hessian>& hessians) const;

    /// \brief Evaluate values, gradients, and Hessians together.
    /// \param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
    /// \param values Receives one value per basis function.
    /// \param gradients Receives one three-component gradient per basis function.
    /// \param hessians Receives one 3-by-3 Hessian per basis function.
    virtual void evaluate_all(const math::Vector<Real, 3>& xi,
                              std::vector<Real>& values,
                              std::vector<Gradient>& gradients,
                              std::vector<Hessian>& hessians) const;

    /// \brief Evaluate basis values into a flat caller-provided buffer.
    /// \param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
    /// \param values_out Output buffer with at least size() entries.
    virtual void evaluate_values_to(const math::Vector<Real, 3>& xi,
                                    Real* SVMP_RESTRICT values_out) const;

    /// \brief Evaluate basis gradients into a flat caller-provided buffer.
    /// \param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
    /// \param gradients_out Output buffer with node-major layout: node * 3 + component.
    virtual void evaluate_gradients_to(const math::Vector<Real, 3>& xi,
                                       Real* SVMP_RESTRICT gradients_out) const;

    /// \brief Evaluate basis Hessians into a flat caller-provided buffer.
    /// \param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
    /// \param hessians_out Output buffer with node-major row-major layout: node * 9 + row * 3 + col.
    virtual void evaluate_hessians_to(const math::Vector<Real, 3>& xi,
                                      Real* SVMP_RESTRICT hessians_out) const;

protected:
    /// \brief Approximate gradients by centered finite differences of values.
    ///
    /// \details This helper exists as a development and fallback utility for
    /// basis implementations that do not yet provide analytical gradients. It
    /// is useful for prototyping new basis families and for checking analytical
    /// derivative formulas in tests. Production element assembly should prefer
    /// analytical gradients when available because finite differences introduce
    /// truncation/roundoff sensitivity and require multiple value evaluations
    /// per reference coordinate.
    void numerical_gradient(const math::Vector<Real, 3>& xi,
                            std::vector<Gradient>& gradients,
                            Real eps = Real(1e-6)) const;

    /// \brief Approximate Hessians by centered finite differences of gradients.
    ///
    /// \details This helper exists for the same reason as numerical_gradient:
    /// it provides a simple reference implementation for prototyping and
    /// derivative verification when analytical second derivatives are not yet
    /// implemented. It depends on evaluate_gradients(), so it is only available
    /// for basis implementations that can already provide gradients. Analytical
    /// Hessians should be used in performance-sensitive solver paths because
    /// finite-difference Hessians amplify numerical error and require repeated
    /// gradient evaluations.
    void numerical_hessian(const math::Vector<Real, 3>& xi,
                           std::vector<Hessian>& hessians,
                           Real eps = Real(1e-5)) const;
};

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_BASISFUNCTION_H
