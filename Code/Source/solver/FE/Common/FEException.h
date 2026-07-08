// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_EXCEPTION_H
#define SVMP_FE_EXCEPTION_H

/**
 * @file FEException.h
 * @brief Exception hierarchy for error handling in the FE library
 *
 * This header defines FE-specific exception types that derive from the shared
 * solver exception infrastructure in Exception.h. FEException marks
 * failures from finite-element assembly, backend, DOF, and element operations.
 */

#include "Exception.h"

#include <sstream>
#include <string>
#include <utility>

namespace svmp::FE {

/**
 * @defgroup FE_CommonExceptions Exceptions
 * @ingroup FE_Common
 * @brief FE exception hierarchy.
 *
 * @details All FE-specific exceptions derive from FEException, which itself
 * derives from the shared solver ExceptionBase. Specialized subclasses carry
 * structured context (element type, DOF index, backend name and error code,
 * iteration counts, Jacobian determinants) so call sites can report
 * actionable diagnostics.
 *
 * Throw FE exceptions through the canonical core helpers in Core/Exception.h:
 *
 * @code
 * svmp::raise<ExceptionT>(message);
 * svmp::throw_if<ExceptionT>(failure_condition, message);
 * svmp::check<ExceptionT>(valid_condition, message);
 * svmp::check_not_null<ExceptionT>(ptr, message);
 * svmp::check_index<ExceptionT>(index, size);
 * svmp::not_implemented(message);
 * @endcode
 *
 * check() raises when its (success) condition is false; throw_if() raises when
 * its (failure) condition is true. FE owns exception types; helper spelling is
 * owned by the core layer.
 * @{
 */

/**
 * @brief Base exception type for errors originating in the FE library
 *
 * Carries a status code and source location alongside the message. Derived
 * classes select an appropriate StatusCode and may attach additional
 * structured context.
 */
class FEException : public ExceptionBase {
public:
    /**
     * @brief Construct with a message and optional status code.
     * @param message Human-readable error description.
     * @param status Status code classifying the failure.
     *
     * @details The source location is stamped by svmp::raise(); construct FE
     * exceptions through the core helpers rather than passing file/line/function.
     */
    explicit FEException(const std::string& message,
                         StatusCode status = StatusCode::Unknown)
        : ExceptionBase(message, status, "FE Exception")
    {
    }

    /**
     * @brief Status code classifying the failure.
     * @return The status code recorded at construction.
     */
    StatusCode status() const noexcept { return status_code(); }
};

/**
 * @brief An argument failed validation
 */
SVMP_DEFINE_EXCEPTION(InvalidArgumentException, FEException,
                      StatusCode::InvalidArgument);

/**
 * @brief Unsupported or malformed element request
 *
 * Records the offending element type so error reports can name it.
 */
class InvalidElementException : public FEException {
public:
    /**
     * @brief Construct with a message and optional element-type context.
     * @param message Human-readable error description.
     * @param element_type Name of the offending element type; appended to the message when non-empty.
     */
    InvalidElementException(const std::string& message,
                            std::string element_type = "")
        : FEException(build_message(message, element_type),
                      StatusCode::InvalidArgument),
          element_type_(std::move(element_type))
    {
    }

    /**
     * @brief Name of the offending element type.
     * @return Element-type name; empty when not provided.
     */
    const std::string& element_type() const noexcept { return element_type_; }

private:
    static std::string build_message(const std::string& message,
                                     const std::string& element_type)
    {
        if (element_type.empty()) {
            return message;
        }

        return message + " (Element type: " + element_type + ")";
    }

    std::string element_type_;
};

/**
 * @brief Degree-of-freedom numbering or lookup failure
 *
 * Records the offending DOF index so error reports can name it.
 */
class DofException : public FEException {
public:
    /**
     * @brief Construct with a message and optional DOF-index context.
     * @param message Human-readable error description.
     * @param dof_index Offending DOF index; appended to the message unless it equals invalid_dof_index().
     */
    DofException(const std::string& message,
                 long long dof_index = invalid_dof_index())
        : FEException(build_message(message, dof_index),
                      StatusCode::InvalidArgument),
          dof_index_(dof_index)
    {
    }

    /**
     * @brief Offending DOF index.
     * @return DOF index; invalid_dof_index() when not provided.
     */
    long long dof_index() const noexcept { return dof_index_; }
    /**
     * @brief Sentinel meaning "no DOF index attached".
     * @return The sentinel value -1.
     */
    static constexpr long long invalid_dof_index() noexcept { return -1; }

private:
    static std::string build_message(const std::string& message,
                                     long long dof_index)
    {
        if (dof_index == invalid_dof_index()) {
            return message;
        }

        return message + " (DOF index: " + std::to_string(dof_index) + ")";
    }

    long long dof_index_ = invalid_dof_index();
};

/**
 * @brief Global assembly failure
 */
SVMP_DEFINE_EXCEPTION(AssemblyException, FEException, StatusCode::InvalidState);

/**
 * @brief Failure reported by a linear-algebra or solver backend
 *
 * Records the backend name and its native error code so error reports can
 * identify the failing dependency.
 */
class BackendException : public FEException {
public:
    /**
     * @brief Construct with a message and optional backend context.
     * @param message Human-readable error description.
     * @param backend_name Name of the failing backend; appended to the message when non-empty.
     * @param error_code Backend-native error code; appended to the message when nonzero.
     */
    BackendException(const std::string& message,
                     std::string backend_name = "",
                     int error_code = 0)
        : FEException(build_message(message, backend_name, error_code),
                      StatusCode::DependencyError),
          backend_name_(std::move(backend_name)),
          error_code_(error_code)
    {
    }

    /**
     * @brief Name of the failing backend.
     * @return Backend name; empty when not provided.
     */
    const std::string& backend_name() const noexcept { return backend_name_; }
    /**
     * @brief Backend-native error code.
     * @return Error code; zero when not provided.
     */
    int error_code() const noexcept { return error_code_; }

private:
    static std::string build_message(const std::string& message,
                                     const std::string& backend_name,
                                     int error_code)
    {
        std::ostringstream oss;
        oss << message;
        if (!backend_name.empty() || error_code != 0) {
            oss << " (";
            if (!backend_name.empty()) {
                oss << "Backend: " << backend_name;
            }
            if (error_code != 0) {
                if (!backend_name.empty()) {
                    oss << ", ";
                }
                oss << "Error code: " << error_code;
            }
            oss << ")";
        }
        return oss.str();
    }

    std::string backend_name_;
    int error_code_ = 0;
};

/**
 * @brief Requested feature is not implemented.
 *
 * @details Alias for svmp::NotImplementedException (Core/Exception.h), the single
 * not-implemented type used across the solver and the default raised by
 * svmp::not_implemented(). Kept in the FE namespace for source compatibility; it
 * derives from CoreException, not FEException.
 */
using NotImplementedException = svmp::NotImplementedException;

/**
 * @brief Required initialization step has not been performed
 */
class NotInitializedException : public FEException {
public:
  /**
   * @brief Construct from the name of the uninitialized feature.
   * @param feature Description of the missing initialization.
   */
  explicit NotInitializedException(const std::string& feature)
      : FEException("Missing initialization: " + feature,
                    StatusCode::InvalidState)
  {
  }
};

/**
 * @brief Iterative process failed to converge
 *
 * Records the iteration count and final residual so error reports can show
 * how far the iteration progressed.
 */
class ConvergenceException : public FEException {
public:
    /**
     * @brief Construct with a message and optional iteration context.
     * @param message Human-readable error description.
     * @param iteration Iteration at which the failure was detected; appended to the message when non-negative.
     * @param residual Final residual; appended to the message when positive.
     */
    ConvergenceException(const std::string& message,
                         int iteration = -1,
                         double residual = 0.0)
        : FEException(build_message(message, iteration, residual),
                      StatusCode::InvalidState),
          iteration_(iteration),
          residual_(residual)
    {
    }

    /**
     * @brief Iteration at which the failure was detected.
     * @return Iteration count; -1 when not provided.
     */
    int iteration() const noexcept { return iteration_; }
    /**
     * @brief Final residual value.
     * @return Residual; 0.0 when not provided.
     */
    double residual() const noexcept { return residual_; }

private:
    static std::string build_message(const std::string& message,
                                     int iteration,
                                     double residual)
    {
        std::ostringstream oss;
        oss << message;
        if (iteration >= 0) {
            oss << " (Iteration: " << iteration;
            if (residual > 0.0) {
                oss << ", Residual: " << residual;
            }
            oss << ")";
        }
        return oss.str();
    }

    int iteration_ = -1;
    double residual_ = 0.0;
};

/**
 * @brief Element geometric mapping is singular or inverted
 *
 * Records the offending Jacobian determinant so error reports can show the
 * degeneracy.
 */
class SingularMappingException : public FEException {
public:
    /**
     * @brief Construct with a message and the offending Jacobian determinant.
     * @param message Human-readable error description.
     * @param jacobian_det Jacobian determinant at the failure point; appended to the message.
     */
    SingularMappingException(const std::string& message,
                             double jacobian_det = 0.0)
        : FEException(build_message(message, jacobian_det),
                      StatusCode::InvalidState),
          jacobian_det_(jacobian_det)
    {
    }

    /**
     * @brief Jacobian determinant at the failure point.
     * @return The determinant recorded at construction.
     */
    double jacobian_det() const noexcept { return jacobian_det_; }

private:
    static std::string build_message(const std::string& message, double det)
    {
        return message + " (Jacobian determinant: " + std::to_string(det) +
               ")";
    }

    double jacobian_det_ = 0.0;
};

/** @} */

} // namespace svmp::FE

#endif // SVMP_FE_EXCEPTION_H
