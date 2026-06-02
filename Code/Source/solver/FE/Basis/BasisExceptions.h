/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#ifndef SVMP_FE_BASIS_BASISEXCEPTIONS_H
#define SVMP_FE_BASIS_BASISEXCEPTIONS_H

#include "FEException.h"

namespace svmp {
namespace FE {
namespace basis {

/**
 * @brief Base exception type for errors originating in the Basis module
 */
class BasisException : public FEException {
public:
    BasisException(const std::string& message,
                   const char* file = "",
                   int line = 0,
                   const char* function = "",
                   StatusCode status = StatusCode::Unknown)
        : FEException(message, status, file, line, function) {}
};

/**
 * @brief Invalid Basis request or configuration
 */
class BasisConfigurationException : public BasisException {
public:
    BasisConfigurationException(const std::string& message,
                                const char* file = "",
                                int line = 0,
                                const char* function = "")
        : BasisException(message, file, line, function, StatusCode::InvalidArgument) {}
};

/**
 * @brief Requested element topology is incompatible with the basis family
 */
class BasisElementCompatibilityException : public BasisException {
public:
    BasisElementCompatibilityException(const std::string& message,
                                       const char* file = "",
                                       int line = 0,
                                       const char* function = "")
        : BasisException(message, file, line, function, StatusCode::InvalidArgument) {}
};

/**
 * @brief Basis evaluation request cannot be satisfied
 */
class BasisEvaluationException : public BasisException {
public:
    BasisEvaluationException(const std::string& message,
                             const char* file = "",
                             int line = 0,
                             const char* function = "")
        : BasisException(message, file, line, function, StatusCode::InvalidArgument) {}
};

/**
 * @brief Public-to-canonical node ordering or coordinate lookup failure
 */
class BasisNodeOrderingException : public BasisException {
public:
    BasisNodeOrderingException(const std::string& message,
                               const char* file = "",
                               int line = 0,
                               const char* function = "")
        : BasisException(message, file, line, function, StatusCode::InvalidArgument) {}
};

/**
 * @brief Internal basis construction or transform setup failure
 */
class BasisConstructionException : public BasisException {
public:
    BasisConstructionException(const std::string& message,
                               const char* file = "",
                               int line = 0,
                               const char* function = "")
        : BasisException(message, file, line, function, StatusCode::InternalError) {}
};

#define BASIS_CHECK_CONFIG(condition, message)                                                 \
    do {                                                                                       \
        if (!(condition)) {                                                                    \
            throw ::svmp::FE::basis::BasisConfigurationException((message),                    \
                                                                  __FILE__, __LINE__, __func__); \
        }                                                                                      \
    } while (false)

#define BASIS_CHECK_COMPAT(condition, message)                                                 \
    do {                                                                                       \
        if (!(condition)) {                                                                    \
            throw ::svmp::FE::basis::BasisElementCompatibilityException((message),             \
                                                                         __FILE__, __LINE__, __func__); \
        }                                                                                      \
    } while (false)

#define BASIS_CHECK_EVAL(condition, message)                                                   \
    do {                                                                                       \
        if (!(condition)) {                                                                    \
            throw ::svmp::FE::basis::BasisEvaluationException((message),                       \
                                                               __FILE__, __LINE__, __func__);  \
        }                                                                                      \
    } while (false)

#define BASIS_CHECK_NODE_ORDER(condition, message)                                             \
    do {                                                                                       \
        if (!(condition)) {                                                                    \
            throw ::svmp::FE::basis::BasisNodeOrderingException((message),                     \
                                                                 __FILE__, __LINE__, __func__); \
        }                                                                                      \
    } while (false)

#define BASIS_CHECK_CONSTRUCTION(condition, message)                                           \
    do {                                                                                       \
        if (!(condition)) {                                                                    \
            throw ::svmp::FE::basis::BasisConstructionException((message),                     \
                                                                 __FILE__, __LINE__, __func__); \
        }                                                                                      \
    } while (false)

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_BASISEXCEPTIONS_H
