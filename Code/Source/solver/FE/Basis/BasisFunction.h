// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_BASIS_BASISFUNCTION_H
#define SVMP_FE_BASIS_BASISFUNCTION_H

#include "BasisExceptions.h"
#include "BasisTraits.h"
#include "Math/Matrix.h"
#include "Math/Vector.h"
#include "Types.h"

#include <cstddef>
#include <span>
#include <vector>

/**
 * @defgroup FE_Basis Basis
 * @ingroup FE
 * @brief Basis-function interfaces, concrete basis families, and reference-node conventions.
 *
 * @details
 * ## Scope
 *
 * The Basis module owns reference-element shape functions. It provides the
 * number of basis functions and the values and derivatives,
 * @f$N_i@f$, @f$\partial N_i / \partial \xi_j@f$, and
 * @f$\partial^2 N_i / \partial \xi_j \partial \xi_k@f$ at reference
 * points. It does not own mesh storage, quadrature selection, field
 * formulation policy, or transformation of derivatives to physical
 * coordinates. Those decisions stay with the solver layer that has the mesh,
 * material model, and equation context.
 *
 * The main pieces are:
 * - BasisFunction (BasisFunction.h): the abstract query and evaluation
 *   contract for code that does not need to know the concrete family.
 * - @ref FE_LagrangeBasis "LagrangeBasis" and
 *   @ref FE_SerendipityBasis "SerendipityBasis": the implemented nodal
 *   families, including analytical first and second derivatives in reference
 *   coordinates.
 * - basis_factory (BasisFactory.h): runtime construction from a BasisRequest.
 *   basis_factory::default_basis_request() centralizes the family/order that
 *   matches each supported element's public node layout.
 * - ReferenceNodeLayout (NodeOrderingConventions.h): canonical reference-node
 *   coordinates and the output ordering used by every basis evaluator.
 * - BasisTraits.h and BasisExceptions.h: topology classification,
 *   compile-time helpers, and module-specific exception types.
 *
 * ## Object and evaluation contract
 *
 * A basis object is immutable after construction. It represents one reference
 * topology, basis family, and effective polynomial order, and can be shared
 * safely across evaluations. Construction may build node lattices or invert
 * interpolation matrices, so callers should construct through basis_factory
 * and cache one instance for each distinct basis request instead of rebuilding
 * inside element loops.
 *
 * Every evaluator takes a three-component reference coordinate. For
 * lower-dimensional elements, only the first dimension() components are
 * active. Returned gradients always have three components and Hessians are
 * always 3-by-3 matrices; inactive reference directions are expected to be
 * zero for conforming lower-dimensional bases. The std::vector overloads are
 * convenient for setup, tests, and adapter code. The *_to overloads write to
 * caller-owned spans; the concrete nodal families (LagrangeBasis,
 * SerendipityBasis) compute directly into the span and so provide the
 * allocation-free path for assembly. The base-class defaults instead evaluate
 * into a temporary and copy into the span, so a basis that implements only the
 * vector form still works through the span API, just without the allocation
 * savings.
 *
 * Outputs are in ReferenceNodeLayout basis order, not necessarily the mesh or
 * solver's native node order. A caller that stores elements in another local
 * ordering must apply the appropriate permutation at the boundary between the
 * basis module and that storage format.
 *
 * ## Inputs and ownership
 *
 * Constructing and evaluating a basis combines several independent choices:
 *
 * - **Element topology comes from the mesh.** The mesh cell type is translated
 *   to ElementType, which defines the reference topology and public node
 *   layout. This is structural information, not a complete discretization
 *   policy.
 * - **Geometry interpolation follows the mesh nodes.** The basis used for the
 *   reference-to-physical map must be compatible with the element's node
 *   count and ordering. For that case, callers normally use
 *   basis_factory::create_default_for(element_type), which selects the
 *   Lagrange or serendipity space associated with that element layout. A
 *   Tetra10 mesh therefore implies a quadratic geometry map; a Hex20 mesh
 *   implies the supported Hex20 serendipity geometry basis.
 * - **Field approximation is chosen by the formulation.** Field bases do not
 *   have to match the geometry map. Mixed formulations, stabilized methods,
 *   enrichment, and convergence studies may use different families or orders
 *   for different fields on the same mesh topology. Those bases should be
 *   requested explicitly with basis_factory::create() and a BasisRequest
 *   naming the desired family, topology, and order.
 * - **Evaluation points come from the caller.** Quadrature rules, probe
 *   points, interpolation targets, and error-sampling locations are outside
 *   this module. The basis only evaluates at the reference coordinates it is
 *   given.
 *
 * @dot "Basis inputs and responsibilities"
 * digraph fe_basis_information_flow {
 *   rankdir=LR;
 *   node [shape=box, fontname=Helvetica, fontsize=10];
 *   mesh     [label="Mesh element type"];
 *   request  [label="BasisRequest\nfamily + order"];
 *   topology [label="Reference topology\nand node layout"];
 *   basis    [label="Basis object", style=filled, fillcolor=lightgray];
 *   points   [label="Reference points"];
 *   outputs  [label="Reference values\nand derivatives"];
 *   mesh -> topology;
 *   request -> basis;
 *   topology -> basis;
 *   basis -> outputs;
 *   points -> outputs;
 * }
 * @enddot
 *
 * ## Reference scope and the solver adapter
 *
 * The solver-facing adapter in nn.cpp is the boundary between this reference
 * basis contract and legacy solver storage. It translates solver element
 * enums to ElementType, obtains cached default bases for mesh/face shape
 * tables, permutes from ReferenceNodeLayout order into solver node order, and
 * stores N, Nx, and, where needed, packed Nxx at Gauss points. At that stage
 * Nx and Nxx are still derivatives with respect to reference coordinates.
 * Physical-coordinate derivatives are formed later, for a particular
 * configuration and element geometry, by composing the cached reference data
 * with the mapping Jacobian (nn::gnn for first derivatives and nn::gn_nxx for
 * second derivatives).
 */

namespace svmp {
namespace FE {
namespace basis {

/** @brief Gradient vector type used by basis evaluators. */
using Gradient = math::Vector<double, 3>;

/** @brief Hessian matrix type used by basis evaluators. */
using Hessian  = math::Matrix<double, 3, 3>;

/**
 * @brief Throw BasisEvaluationException when an output span is smaller than the
 * basis size. \p label is the full "Class::method" context used in the message,
 * so each basis family passes its own qualified name.
 */
void require_span_size(std::size_t actual, std::size_t expected, const char* label);

/**
 * @brief Check a requested output span unless it is empty, following the
 * "skip this output" convention used by the combined evaluators.
 */
template <typename T>
void require_requested_span_size(std::span<T> output,
                                 std::size_t expected,
                                 const char* label) {
    if (!output.empty()) {
        require_span_size(output.size(), expected, label);
    }
}

/**
 * @brief Abstract interface for finite-element basis-function families.
 * @ingroup FE_Basis
 *
 * BasisFunction defines the common query and evaluation API used by solver
 * code that does not need to know the concrete basis implementation. Derived
 * classes provide values at minimum and can override analytical gradients,
 * Hessians, combined evaluation, and span output paths. The interface
 * is deliberately limited to reference-space quantities; callers own node
 * ordering translation, physical mapping, and any field-level discretization
 * policy.
 */
class BasisFunction {
public:
    /** @brief Destroy a basis function through the abstract interface. */
    virtual ~BasisFunction() = default;

    /**
     * @brief Return the concrete basis family.
     * @return Basis family identifier.
     */
    virtual BasisType basis_type() const noexcept = 0;

    /**
     * @brief Return the reference topology of this basis.
     *
     * @details Together with order(), this is the authoritative identity of a
     * basis: a topology plus a polynomial order, with no node-count assumption.
     * Arbitrary-order bases are constructed from a BasisTopology and an order;
     * named ElementType layouts (Hex8, Hex27, ...) are a fixed-order shorthand
     * that maps to the same (topology, order) pair.
     *
     * @return Reference topology.
     */
    virtual BasisTopology topology() const noexcept = 0;

    /**
     * @brief Return the reference-space dimension of the basis.
     * @return Reference dimension, from zero for points through three for volume elements.
     */
    virtual int dimension() const noexcept = 0;

    /**
     * @brief Return the polynomial order represented by this basis.
     * @return Effective polynomial order after any element-family normalization.
     */
    virtual int order() const noexcept = 0;

    /**
     * @brief Return the number of basis functions and reference nodes.
     * @return Basis function count.
     */
    virtual std::size_t size() const noexcept = 0;

    /**
     * @brief Return the reference interpolation nodes in basis ordering.
     *
     * @details Nodal families return one reference-element coordinate per basis
     * function, in the same order as the evaluator outputs. Bases that do not
     * define interpolation nodes (non-nodal families, or abstract base usage)
     * return an empty vector. The returned reference is valid for the lifetime
     * of the basis object.
     *
     * @return Reference node coordinates: size() entries for nodal families,
     *         empty otherwise.
     */
    virtual const std::vector<math::Vector<double, 3>>& nodes() const noexcept;

    /**
     * @brief Evaluate basis function values at a reference coordinate.
     * @param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
     * @param values Receives one value per basis function.
     */
    virtual void evaluate_values(const math::Vector<double, 3>& xi,
                                 std::vector<double>& values) const = 0;

    /**
     * @brief Evaluate basis gradients at a reference coordinate.
     * @param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
     * @param gradients Receives one three-component gradient per basis function.
     * @throws BasisEvaluationException If gradients are not available for the basis.
     */
    virtual void evaluate_gradients(const math::Vector<double, 3>& xi,
                                    std::vector<Gradient>& gradients) const;

    /**
     * @brief Evaluate basis Hessians at a reference coordinate.
     * @param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
     * @param hessians Receives one 3-by-3 Hessian per basis function.
     * @throws BasisEvaluationException If Hessians are not available for the basis.
     */
    virtual void evaluate_hessians(const math::Vector<double, 3>& xi,
                                   std::vector<Hessian>& hessians) const;

    /**
     * @brief Evaluate values, gradients, and Hessians together.
     * @param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
     * @param values Receives one value per basis function.
     * @param gradients Receives one three-component gradient per basis function.
     * @param hessians Receives one 3-by-3 Hessian per basis function.
     */
    virtual void evaluate_all(const math::Vector<double, 3>& xi,
                              std::vector<double>& values,
                              std::vector<Gradient>& gradients,
                              std::vector<Hessian>& hessians) const;

    /**
     * @brief Evaluate basis values into caller-provided storage.
     * @param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
     * @param values_out Output span with at least size() entries.
     * @note The base-class default evaluates into a temporary and copies; nodal
     *       families override this to write directly into the span.
     */
    virtual void evaluate_values_to(const math::Vector<double, 3>& xi,
                                    std::span<double> values_out) const;

    /**
     * @brief Evaluate basis gradients into caller-provided storage.
     * @param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
     * @param gradients_out Output span with at least size() entries.
     * @note The base-class default evaluates into a temporary and copies; nodal
     *       families override this to write directly into the span.
     */
    virtual void evaluate_gradients_to(const math::Vector<double, 3>& xi,
                                       std::span<Gradient> gradients_out) const;

    /**
     * @brief Evaluate basis Hessians into caller-provided storage.
     * @param xi Reference coordinate. Lower-dimensional elements use the active prefix components.
     * @param hessians_out Output span with at least size() entries.
     * @note The base-class default evaluates into a temporary and copies; nodal
     *       families override this to write directly into the span.
     */
    virtual void evaluate_hessians_to(const math::Vector<double, 3>& xi,
                                      std::span<Hessian> hessians_out) const;

protected:
    /**
     * @brief Approximate gradients by centered finite differences of values.
     *
     * @details This helper is primarily a verification utility for tests: it
     * provides a basis-independent reference that checks a concrete basis's
     * analytical evaluate_gradients() against centered finite differences of
     * evaluate_values(). It lives on the base class so any BasisFunction can be
     * checked uniformly, and having no production caller is by design — every
     * shipped basis supplies analytical gradients. Centered differences add
     * truncation/roundoff sensitivity and require multiple value evaluations
     * per reference coordinate, so analytical gradients are always preferred
     * outside this testing context.
     */
    void numerical_gradient(const math::Vector<double, 3>& xi,
                            std::vector<Gradient>& gradients,
                            double eps = double(1e-6)) const;

    /**
     * @brief Approximate Hessians by centered finite differences of gradients.
     *
     * @details Companion verification utility to numerical_gradient: it checks
     * a basis's analytical evaluate_hessians() against centered finite
     * differences of evaluate_gradients(). Because it differentiates gradients,
     * it is only meaningful for bases that already provide them. Like
     * numerical_gradient it is test-support rather than a production fallback —
     * finite-difference Hessians amplify numerical error and require repeated
     * gradient evaluations, so analytical Hessians are used everywhere outside
     * tests.
     */
    void numerical_hessian(const math::Vector<double, 3>& xi,
                           std::vector<Hessian>& hessians,
                           double eps = double(1e-5)) const;
};

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_BASISFUNCTION_H
