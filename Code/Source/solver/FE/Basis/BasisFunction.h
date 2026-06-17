// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_BASIS_BASISFUNCTION_H
#define SVMP_FE_BASIS_BASISFUNCTION_H

#include "BasisExceptions.h"
#include "Math/Matrix.h"
#include "Math/Vector.h"
#include "Types.h"

#include <cstddef>
#include <span>
#include <vector>

/// \defgroup FE_Basis Basis
/// \ingroup FE
/// \brief Basis-function interfaces, concrete basis families, and reference-node conventions.
///
/// \details
/// ## Scope
///
/// The Basis module owns reference-element shape functions. It provides the
/// number of basis functions and the values and derivatives,
/// \f$N_i\f$, \f$\partial N_i / \partial \xi_j\f$, and
/// \f$\partial^2 N_i / \partial \xi_j \partial \xi_k\f$ at reference
/// points. It does not own mesh storage, quadrature selection, field
/// formulation policy, or transformation of derivatives to physical
/// coordinates. Those decisions stay with the solver layer that has the mesh,
/// material model, and equation context.
///
/// The main pieces are:
/// - BasisFunction (BasisFunction.h): the abstract query and evaluation
///   contract for code that does not need to know the concrete family.
/// - \ref FE_LagrangeBasis "LagrangeBasis" and
///   \ref FE_SerendipityBasis "SerendipityBasis": the implemented nodal
///   families, including analytical first and second derivatives in reference
///   coordinates.
/// - basis_factory (BasisFactory.h): runtime construction from a BasisRequest.
///   basis_factory::default_basis_request() centralizes the family/order that
///   matches each supported element's public node layout.
/// - ReferenceNodeLayout (NodeOrderingConventions.h): canonical reference-node
///   coordinates and the output ordering used by every basis evaluator.
/// - BasisTraits.h and BasisExceptions.h: topology classification,
///   compile-time helpers, and module-specific exception types.
///
/// ## Object and evaluation contract
///
/// A basis object is immutable after construction. It represents one reference
/// topology, basis family, and effective polynomial order, and can be shared
/// safely across evaluations. Construction may build node lattices or invert
/// interpolation matrices, so callers should construct through basis_factory
/// and cache one instance for each distinct basis request instead of rebuilding
/// inside element loops.
///
/// Every evaluator takes a three-component reference coordinate. For
/// lower-dimensional elements, only the first dimension() components are
/// active. Returned gradients always have three components and Hessians are
/// always 3-by-3 matrices; inactive reference directions are expected to be
/// zero for conforming lower-dimensional bases. The std::vector overloads are
/// convenient for setup, tests, and adapter code. The *_to overloads write to
/// caller-owned spans and are the allocation-free path for assembly.
///
/// Outputs are in ReferenceNodeLayout basis order, not necessarily the mesh or
/// solver's native node order. A caller that stores elements in another local
/// ordering must apply the appropriate permutation at the boundary between the
/// basis module and that storage format.
///
/// ## Inputs and ownership
///
/// Constructing and evaluating a basis combines several independent choices:
///
/// - **Element topology comes from the mesh.** The mesh cell type is translated
///   to ElementType, which defines the reference topology and public node
///   layout. This is structural information, not a complete discretization
///   policy.
/// - **Geometry interpolation follows the mesh nodes.** The basis used for the
///   reference-to-physical map must be compatible with the element's node
///   count and ordering. For that case, callers normally use
///   basis_factory::create_default_for(element_type), which selects the
///   Lagrange or serendipity space associated with that element layout. A
///   Tetra10 mesh therefore implies a quadratic geometry map; a Hex20 mesh
///   implies the supported Hex20 serendipity geometry basis.
/// - **Field approximation is chosen by the formulation.** Field bases do not
///   have to match the geometry map. Mixed formulations, stabilized methods,
///   enrichment, and convergence studies may use different families or orders
///   for different fields on the same mesh topology. Those bases should be
///   requested explicitly with basis_factory::create() and a BasisRequest
///   naming the desired family and order.
/// - **Evaluation points come from the caller.** Quadrature rules, probe
///   points, interpolation targets, and error-sampling locations are outside
///   this module. The basis only evaluates at the reference coordinates it is
///   given.
///
/// \dot "Basis inputs and responsibilities"
/// digraph fe_basis_information_flow {
///   rankdir=LR;
///   node [shape=box, fontname=Helvetica, fontsize=10];
///   mesh     [label="Mesh element type"];
///   request  [label="BasisRequest\nfamily + order"];
///   topology [label="Reference topology\nand node layout"];
///   basis    [label="Basis object", style=filled, fillcolor=lightgray];
///   points   [label="Reference points"];
///   outputs  [label="Reference values\nand derivatives"];
///   mesh -> topology;
///   request -> basis;
///   topology -> basis;
///   basis -> outputs;
///   points -> outputs;
/// }
/// \enddot
///
/// ## Reference scope and the solver adapter
///
/// The solver-facing adapter in nn.cpp is the boundary between this reference
/// basis contract and legacy solver storage. It translates solver element
/// enums to ElementType, obtains cached default bases for mesh/face shape
/// tables, permutes from ReferenceNodeLayout order into solver node order, and
/// stores N, Nx, and, where needed, packed Nxx at Gauss points. At that stage
/// Nx and Nxx are still derivatives with respect to reference coordinates.
/// Physical-coordinate derivatives are formed later, for a particular
/// configuration and element geometry, by composing the cached reference data
/// with the mapping Jacobian (nn::gnn for first derivatives and nn::gn_nxx for
/// second derivatives).

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
    Hessian hessian = Hessian::Zero();
    hessian(0, 0) = xx;
    hessian(1, 1) = yy;
    hessian(2, 2) = zz;
    hessian(0, 1) = hessian(1, 0) = xy;
    hessian(0, 2) = hessian(2, 0) = xz;
    hessian(1, 2) = hessian(2, 1) = yz;
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
/// Hessians, combined evaluation, and span output paths. The interface
/// is deliberately limited to reference-space quantities; callers own node
/// ordering translation, physical mapping, and any field-level discretization
/// policy.
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

    /// \brief Return the reference interpolation nodes in basis ordering.
    ///
    /// \details Nodal families return one reference-element coordinate per basis
    /// function, in the same order as the evaluator outputs. Bases that do not
    /// define interpolation nodes (non-nodal families, or abstract base usage)
    /// return an empty vector. The returned reference is valid for the lifetime
    /// of the basis object.
    ///
    /// \return Reference node coordinates: size() entries for nodal families,
    ///         empty otherwise.
    virtual const std::vector<math::Vector<Real, 3>>& nodes() const noexcept;

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

    /// \brief Evaluate basis values into caller-provided storage.
    /// \param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
    /// \param values_out Output span with at least size() entries.
    virtual void evaluate_values_to(const math::Vector<Real, 3>& xi,
                                    std::span<Real> values_out) const;

    /// \brief Evaluate basis gradients into caller-provided storage.
    /// \param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
    /// \param gradients_out Output span with at least size() entries.
    virtual void evaluate_gradients_to(const math::Vector<Real, 3>& xi,
                                       std::span<Gradient> gradients_out) const;

    /// \brief Evaluate basis Hessians into caller-provided storage.
    /// \param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
    /// \param hessians_out Output span with at least size() entries.
    virtual void evaluate_hessians_to(const math::Vector<Real, 3>& xi,
                                      std::span<Hessian> hessians_out) const;

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
