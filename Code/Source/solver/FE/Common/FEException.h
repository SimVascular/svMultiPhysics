// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_EXCEPTION_H
#define SVMP_FE_EXCEPTION_H

/**
 * @file FEException.h
 * @brief Exception hierarchy for error handling in the FE library
 *
 * This header defines FE-specific exception types that derive from the shared
 * solver exception infrastructure in Exception.h.
 */

#include "Exception.h"

#include <sstream>
#include <string>
#include <utility>

namespace svmp {
namespace FE {

class FEException : public ExceptionBase {
public:
    FEException(const std::string& message,
                StatusCode status = StatusCode::Unknown,
                const char* file = "",
                int line = 0,
                const char* function = "")
        : ExceptionBase(prefix_status_message(message, status),
                        status,
                        file,
                        line,
                        function),
          status_(status)
    {
        rebuild_what(subsystem_label());
    }

    FEException(const std::string& message,
                const char* file,
                int line,
                const char* function = "")
        : FEException(message, StatusCode::Unknown, file, line, function)
    {
    }

    StatusCode status() const noexcept { return status_; }

protected:
    const char* subsystem_label() const noexcept override
    {
        return "FE Exception";
    }

private:
    static std::string prefix_status_message(const std::string& message,
                                             StatusCode status)
    {
        if (status == StatusCode::Success || status == StatusCode::Unknown) {
            return message;
        }

        return std::string("[") + status_code_to_string(status) + "] " +
               message;
    }

    StatusCode status_ = StatusCode::Unknown;
};

class InvalidArgumentException : public FEException {
public:
    InvalidArgumentException(const std::string& message,
                             const char* file = "",
                             int line = 0,
                             const char* function = "")
        : FEException(message, StatusCode::InvalidArgument, file, line,
                      function)
    {
    }
};

class InvalidElementException : public FEException {
public:
    InvalidElementException(const std::string& message,
                            std::string element_type = "",
                            const char* file = "",
                            int line = 0,
                            const char* function = "")
        : FEException(build_message(message, element_type),
                      StatusCode::InvalidArgument,
                      file,
                      line,
                      function),
          element_type_(std::move(element_type))
    {
    }

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

class DofException : public FEException {
public:
    DofException(const std::string& message,
                 long long dof_index = invalid_dof_index(),
                 const char* file = "",
                 int line = 0,
                 const char* function = "")
        : FEException(build_message(message, dof_index),
                      StatusCode::InvalidArgument,
                      file,
                      line,
                      function),
          dof_index_(dof_index)
    {
    }

    long long dof_index() const noexcept { return dof_index_; }
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

class AssemblyException : public FEException {
public:
    AssemblyException(const std::string& message,
                      const char* file = "",
                      int line = 0,
                      const char* function = "")
        : FEException(message, StatusCode::InvalidState, file, line, function)
    {
    }
};

class BackendException : public FEException {
public:
    BackendException(const std::string& message,
                     std::string backend_name = "",
                     int error_code = 0,
                     const char* file = "",
                     int line = 0,
                     const char* function = "")
        : FEException(build_message(message, backend_name, error_code),
                      StatusCode::DependencyError,
                      file,
                      line,
                      function),
          backend_name_(std::move(backend_name)),
          error_code_(error_code)
    {
    }

    const std::string& backend_name() const noexcept { return backend_name_; }
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

class NotImplementedException : public FEException {
public:
    NotImplementedException(const std::string& feature,
                            const char* file = "",
                            int line = 0,
                            const char* function = "")
        : FEException("Feature not implemented: " + feature,
                      StatusCode::NotImplemented,
                      file,
                      line,
                      function)
    {
    }
};

class ConvergenceException : public FEException {
public:
    ConvergenceException(const std::string& message,
                         int iteration = -1,
                         double residual = 0.0,
                         const char* file = "",
                         int line = 0,
                         const char* function = "")
        : FEException(build_message(message, iteration, residual),
                      StatusCode::InvalidState,
                      file,
                      line,
                      function),
          iteration_(iteration),
          residual_(residual)
    {
    }

    int iteration() const noexcept { return iteration_; }
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

class SingularMappingException : public FEException {
public:
    SingularMappingException(const std::string& message,
                             double jacobian_det = 0.0,
                             const char* file = "",
                             int line = 0,
                             const char* function = "")
        : FEException(build_message(message, jacobian_det),
                      StatusCode::InvalidState,
                      file,
                      line,
                      function),
          jacobian_det_(jacobian_det)
    {
    }

    double jacobian_det() const noexcept { return jacobian_det_; }

private:
    static std::string build_message(const std::string& message, double det)
    {
        return message + " (Jacobian determinant: " + std::to_string(det) +
               ")";
    }

    double jacobian_det_ = 0.0;
};

template <class ExceptionT, class... Args>
[[noreturn]] inline void raise(SourceLocation location, Args&&... args)
{
    ::svmp::raise<ExceptionT>(location, std::forward<Args>(args)...);
}

template <class ExceptionT = FEException, class... Args>
inline void throw_if(bool condition, SourceLocation location, Args&&... args)
{
    if (condition) {
        ::svmp::FE::raise<ExceptionT>(location, std::forward<Args>(args)...);
    }
}

template <class ExceptionT = InvalidArgumentException, class... Args>
inline void check_arg(bool condition, SourceLocation location, Args&&... args)
{
    ::svmp::check_arg<ExceptionT>(condition, location,
                                  std::forward<Args>(args)...);
}

template <class ExceptionT = InvalidArgumentException, class PointerT,
          class... Args>
inline PointerT* check_not_null(PointerT* ptr, SourceLocation location,
                                Args&&... args)
{
    return ::svmp::check_not_null<ExceptionT>(ptr, location,
                                              std::forward<Args>(args)...);
}

template <class ExceptionT = InvalidArgumentException, class IndexT,
          class SizeT>
inline void check_index(IndexT index, SizeT size, SourceLocation location)
{
    const long long fe_check_index_value = static_cast<long long>(index);
    const long long fe_check_size_value = static_cast<long long>(size);

    ::svmp::FE::check_arg<ExceptionT>(
        fe_check_index_value >= 0 &&
            fe_check_index_value < fe_check_size_value,
        location,
        "Index " + std::to_string(fe_check_index_value) +
            " out of bounds [0, " + std::to_string(fe_check_size_value) + ")");
}

[[noreturn]] inline void not_implemented(const std::string& feature,
                                         SourceLocation location)
{
    ::svmp::FE::raise<NotImplementedException>(location, feature);
}

} // namespace FE
} // namespace svmp

#endif // SVMP_FE_EXCEPTION_H
