// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_BASIS_BASISEXCEPTIONS_H
#define SVMP_FE_BASIS_BASISEXCEPTIONS_H

#include "FEException.h"

namespace svmp::FE::basis {

/**
 * @defgroup FE_BasisExceptions Exceptions
 * @ingroup FE_Basis
 * @brief Basis-module exception hierarchy.
 *
 * @details Every Basis exception derives from BasisException (and thus FEException),
 * so a caller can catch a specific basis failure or the FE base type. See
 * BasisException for why the module raises these basis-specific types rather than
 * the generic FE exceptions.
 * @{
 */

/**
 * @brief Base exception type for errors originating in the Basis module
 *
 * @details The Basis module raises these basis-specific types -- rather than the
 * generic FE exceptions in FEException.h -- so a caller can catch a basis failure
 * precisely: an unsupported element/order pairing, a non-unisolvent node set, an
 * out-of-range reference-node index. They all derive from FEException, so code
 * that only wants "some FE error" can still catch the base type, and they carry
 * the same StatusCode and source location as the rest of the hierarchy.
 */
class BasisException : public FEException {
public:
    explicit BasisException(const std::string& message,
                            StatusCode status = StatusCode::Unknown)
        : FEException(message, status) {}
};

/**
 * @brief Invalid Basis request or configuration
 *
 * @details Raised when a request is malformed before any geometry is built: a
 * missing or negative polynomial order, a named element layout paired with an
 * explicit order that does not match its fixed order, or a field type or
 * continuity the scalar Lagrange/Serendipity factory does not support. Example:
 * constructing a LagrangeBasis for Tetra10 at order 1, when that layout is fixed
 * at order 2.
 */
SVMP_DEFINE_EXCEPTION(BasisConfigurationException, BasisException,
                      StatusCode::InvalidArgument);

/**
 * @brief Requested element topology is incompatible with the basis family
 *
 * @details Raised when the family cannot represent the requested topology or
 * named layout. Example: requesting wedge serendipity through the arbitrary-order
 * topology path (only the named Wedge15 layout is supported), or requesting a
 * basis on ElementType::Unknown.
 */
SVMP_DEFINE_EXCEPTION(BasisElementCompatibilityException, BasisException,
                      StatusCode::InvalidArgument);

/**
 * @brief Basis evaluation request cannot be satisfied
 *
 * @details Raised at evaluation time rather than construction time. Example: an
 * output span smaller than size(), or requesting analytical gradients or Hessians
 * from a basis that does not provide them.
 */
SVMP_DEFINE_EXCEPTION(BasisEvaluationException, BasisException,
                      StatusCode::InvalidArgument);

/**
 * @brief Public-to-canonical node ordering or coordinate lookup failure
 *
 * @details Raised when a node index or coordinate lookup falls outside the
 * reference layout. Example: requesting a tensor-axis node index outside
 * [0, order] from line_coord_pm_one.
 */
SVMP_DEFINE_EXCEPTION(BasisNodeOrderingException, BasisException,
                      StatusCode::InvalidArgument);

/**
 * @brief Internal basis construction or transform setup failure
 *
 * @details Signals a violated internal invariant during setup (StatusCode::
 * InternalError) rather than bad user input. Example: a generated Lagrange node
 * lattice whose index components fall outside [0, order] in get_lagrange_lattice.
 */
SVMP_DEFINE_EXCEPTION(BasisConstructionException, BasisException,
                      StatusCode::InternalError);

/** @} */

} // namespace svmp::FE::basis

#endif // SVMP_FE_BASIS_BASISEXCEPTIONS_H
