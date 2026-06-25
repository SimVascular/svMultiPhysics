// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_BASIS_BASISEXCEPTIONS_H
#define SVMP_FE_BASIS_BASISEXCEPTIONS_H

#include "FEException.h"

namespace svmp {
namespace FE {
namespace basis {

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
    BasisException(const std::string& message,
                   const char* file,
                   int line,
                   const char* function,
                   StatusCode status = StatusCode::Unknown)
        : FEException(message, status, file, line, function) {}
};

/**
 * @brief Invalid Basis request or configuration
 */
class BasisConfigurationException : public BasisException {
public:
    BasisConfigurationException(const std::string& message,
                                const char* file,
                                int line,
                                const char* function)
        : BasisException(message, file, line, function, StatusCode::InvalidArgument) {}
};

/**
 * @brief Requested element topology is incompatible with the basis family
 */
class BasisElementCompatibilityException : public BasisException {
public:
    BasisElementCompatibilityException(const std::string& message,
                                       const char* file,
                                       int line,
                                       const char* function)
        : BasisException(message, file, line, function, StatusCode::InvalidArgument) {}
};

/**
 * @brief Basis evaluation request cannot be satisfied
 */
class BasisEvaluationException : public BasisException {
public:
    BasisEvaluationException(const std::string& message,
                             const char* file,
                             int line,
                             const char* function)
        : BasisException(message, file, line, function, StatusCode::InvalidArgument) {}
};

/**
 * @brief Public-to-canonical node ordering or coordinate lookup failure
 */
class BasisNodeOrderingException : public BasisException {
public:
    BasisNodeOrderingException(const std::string& message,
                               const char* file,
                               int line,
                               const char* function)
        : BasisException(message, file, line, function, StatusCode::InvalidArgument) {}
};

/**
 * @brief Internal basis construction or transform setup failure
 */
class BasisConstructionException : public BasisException {
public:
    BasisConstructionException(const std::string& message,
                               const char* file,
                               int line,
                               const char* function)
        : BasisException(message, file, line, function, StatusCode::InternalError) {}
};

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_BASISEXCEPTIONS_H
