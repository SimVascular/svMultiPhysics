// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_EXCEPTION_H
#define SVMP_FE_EXCEPTION_H

/**
 * @file FEException.h
 * @brief Exception hierarchy for error handling in the FE library
 *
 * This header defines a comprehensive exception system for the FE library
 * with support for detailed error messages, stack traces, MPI-aware error
 * handling, and structured error codes.
 */

#include "Types.h"
#include "FEConfig.h"
#include <exception>
#include <string>
#include <sstream>
#include <vector>
#include <memory>

#if FE_HAS_MPI
#include <mpi.h>
#include <iostream>
#endif

// Platform-specific includes for stack traces
#if defined(__GNUC__) && !defined(_WIN32)
#include <execinfo.h>
#include <cxxabi.h>
#endif

namespace svmp {
namespace FE {

// ============================================================================
// Base Exception Class
// ============================================================================

/**
 * @brief Base exception class for all FE library exceptions
 *
 * Provides rich error information including:
 * - Detailed error messages
 * - Source file and line information
 * - Optional stack traces
 * - MPI rank information in parallel runs
 * - Nested exception support
 */
class FEException : public std::exception {
public:
    /**
     * @brief Construct exception with message
     */
    FEException(const std::string& message,
                FEStatus status = FEStatus::Unknown)
        : message_(message),
          status_(status),
          file_(""),
          line_(0),
          function_(""),
          mpi_rank_(-1) {
        capture_context();
        build_what();
    }

    /**
     * @brief Construct exception with source location
     */
    FEException(const std::string& message,
                const char* file,
                int line,
                const char* function = "",
                FEStatus status = FEStatus::Unknown)
        : message_(message),
          status_(status),
          file_(file),
          line_(line),
          function_(function),
          mpi_rank_(-1) {
        capture_context();
        build_what();
    }

    /**
     * @brief Copy constructor
     */
    FEException(const FEException&) = default;

    /**
     * @brief Virtual destructor
     */
    virtual ~FEException() noexcept = default;

    /**
     * @brief Get exception message
     */
    virtual const char* what() const noexcept override {
        return what_.c_str();
    }

    /**
     * @brief Get error status code
     */
    FEStatus status() const noexcept {
        return status_;
    }

    /**
     * @brief Get source file where exception was thrown
     */
    const std::string& file() const noexcept {
        return file_;
    }

    /**
     * @brief Get line number where exception was thrown
     */
    int line() const noexcept {
        return line_;
    }

    /**
     * @brief Get function name where exception was thrown
     */
    const std::string& function() const noexcept {
        return function_;
    }

    /**
     * @brief Get MPI rank (if applicable)
     */
    int mpi_rank() const noexcept {
        return mpi_rank_;
    }

    /**
     * @brief Get stack trace (if available)
     */
    const std::vector<std::string>& stack_trace() const noexcept {
        return stack_trace_;
    }

    /**
     * @brief Add nested exception context
     */
    void add_context(const std::string& context) {
        message_ = context + "\n  -> " + message_;
        build_what();
    }

protected:
    std::string message_;
    FEStatus status_;
    std::string file_;
    int line_;
    std::string function_;
    int mpi_rank_;
    std::vector<std::string> stack_trace_;
    std::string what_;

    /**
     * @brief Capture MPI rank and stack trace
     */
    void capture_context() {
        // Get MPI rank if available
        #if FE_HAS_MPI
        int initialized = 0;
        MPI_Initialized(&initialized);
        if (initialized) {
            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
        }
        #endif

        // Capture stack trace if in debug mode
        #if FE_DEBUG_MODE
        capture_stack_trace();
        #endif
    }

    /**
     * @brief Capture stack trace (platform-specific)
     */
    void capture_stack_trace() {
        #if defined(__GNUC__) && !defined(_WIN32)
        constexpr int MAX_FRAMES = 32;
        void* frames[MAX_FRAMES];
        int n_frames = backtrace(frames, MAX_FRAMES);

        char** symbols = backtrace_symbols(frames, n_frames);
        if (symbols) {
            if (n_frames > 0) {
                stack_trace_.reserve(static_cast<std::size_t>(n_frames));
            }
            for (int i = 1; i < n_frames; ++i) {  // Skip this function
                std::string symbol(symbols[i]);

                // Try to demangle C++ names
                size_t start = symbol.find('(');
                size_t end = symbol.find('+', start);
                if (start != std::string::npos && end != std::string::npos) {
                    std::string mangled = symbol.substr(start + 1, end - start - 1);
                    int status;
                    char* demangled = abi::__cxa_demangle(mangled.c_str(), nullptr, nullptr, &status);
                    if (status == 0 && demangled) {
                        symbol.replace(start + 1, end - start - 1, demangled);
                        free(demangled);
                    }
                }

                stack_trace_.push_back(symbol);
            }
            free(symbols);
        }
        #endif
    }

    /**
     * @brief Build the complete error message
     */
    void build_what() {
        std::ostringstream oss;

        // Error type and status
        oss << "[FE Exception] " << status_to_string(status_);

        // MPI rank if applicable
        if (mpi_rank_ >= 0) {
            oss << " (Rank " << mpi_rank_ << ")";
        }

        oss << "\n";

        // Location information
        if (!file_.empty()) {
            oss << "  Location: " << file_ << ":" << line_;
            if (!function_.empty()) {
                oss << " in " << function_ << "()";
            }
            oss << "\n";
        }

        // Error message
        oss << "  Message: " << message_ << "\n";

        // Stack trace in debug mode
        #if FE_DEBUG_MODE
        if (!stack_trace_.empty()) {
            oss << "  Stack trace:\n";
            for (size_t i = 0; i < stack_trace_.size() && i < 10; ++i) {
                oss << "    #" << i << " " << stack_trace_[i] << "\n";
            }
        }
        #endif

        what_ = oss.str();
    }
};

// ============================================================================
// Specific Exception Types
// ============================================================================

/**
 * @brief Exception for invalid arguments
 */
class InvalidArgumentException : public FEException {
public:
    InvalidArgumentException(const std::string& message,
                            const char* file = "",
                            int line = 0,
                            const char* function = "")
        : FEException(message, file, line, function, FEStatus::InvalidArgument) {}
};

/**
 * @brief Exception for invalid elements
 */
class InvalidElementException : public FEException {
public:
    InvalidElementException(const std::string& message,
                           ElementType elem_type,
                           const char* file = "",
                           int line = 0)
        : FEException(build_message(message, elem_type), file, line, "", FEStatus::InvalidElement),
          element_type_(elem_type) {}

    ElementType element_type() const { return element_type_; }

private:
    ElementType element_type_;

    static std::string build_message(const std::string& msg, ElementType elem) {
        return msg + " (Element type: " + std::to_string(static_cast<int>(elem)) + ")";
    }
};

/**
 * @brief Exception for DOF-related errors
 */
class DofException : public FEException {
public:
    DofException(const std::string& message,
                 GlobalIndex dof_index = INVALID_GLOBAL_INDEX,
                 const char* file = "",
                 int line = 0)
        : FEException(build_message(message, dof_index), file, line),
          dof_index_(dof_index) {}

    GlobalIndex dof_index() const { return dof_index_; }

private:
    GlobalIndex dof_index_;

    static std::string build_message(const std::string& msg, GlobalIndex dof) {
        if (dof != INVALID_GLOBAL_INDEX) {
            return msg + " (DOF index: " + std::to_string(dof) + ")";
        }
        return msg;
    }
};

/**
 * @brief Exception for assembly errors
 */
class AssemblyException : public FEException {
public:
    AssemblyException(const std::string& message,
                      const char* file = "",
                      int line = 0)
        : FEException(message, file, line, "", FEStatus::AssemblyError) {}
};

/**
 * @brief Exception for backend-specific errors
 */
class BackendException : public FEException {
public:
    BackendException(const std::string& message,
                     config::Backend backend,
                     int error_code = 0,
                     const char* file = "",
                     int line = 0)
        : FEException(build_message(message, backend, error_code), file, line, "", FEStatus::BackendError),
          backend_(backend),
          error_code_(error_code) {}

    config::Backend backend() const { return backend_; }
    int error_code() const { return error_code_; }

private:
    config::Backend backend_;
    int error_code_;

    static std::string build_message(const std::string& msg, config::Backend backend, int code) {
        std::ostringstream oss;
        oss << msg << " (Backend: ";
        switch(backend) {
            case config::Backend::Trilinos: oss << "Trilinos"; break;
            case config::Backend::PETSc: oss << "PETSc"; break;
            case config::Backend::Eigen: oss << "Eigen"; break;
            default: oss << "Unknown";
        }
        if (code != 0) {
            oss << ", Error code: " << code;
        }
        oss << ")";
        return oss.str();
    }
};

/**
 * @brief Exception for not implemented features
 */
class NotImplementedException : public FEException {
public:
    NotImplementedException(const std::string& feature,
                           const char* file = "",
                           int line = 0,
                           const char* function = "")
        : FEException("Feature not implemented: " + feature, file, line, function, FEStatus::NotImplemented) {}
};

/**
 * @brief Exception for convergence failures
 */
class ConvergenceException : public FEException {
public:
    ConvergenceException(const std::string& message,
                        int iteration = -1,
                        Real residual = 0.0,
                        const char* file = "",
                        int line = 0)
        : FEException(build_message(message, iteration, residual), file, line, "", FEStatus::ConvergenceError),
          iteration_(iteration),
          residual_(residual) {}

    int iteration() const { return iteration_; }
    Real residual() const { return residual_; }

private:
    int iteration_;
    Real residual_;

    static std::string build_message(const std::string& msg, int iter, Real res) {
        std::ostringstream oss;
        oss << msg;
        if (iter >= 0) {
            oss << " (Iteration: " << iter;
            if (res > 0) {
                oss << ", Residual: " << res;
            }
            oss << ")";
        }
        return oss.str();
    }
};

/**
 * @brief Exception for singular mappings
 */
class SingularMappingException : public FEException {
public:
    SingularMappingException(const std::string& message,
                            Real jacobian_det = 0.0,
                            const char* file = "",
                            int line = 0)
        : FEException(build_message(message, jacobian_det), file, line, "", FEStatus::SingularMapping),
          jacobian_det_(jacobian_det) {}

    Real jacobian_det() const { return jacobian_det_; }

private:
    Real jacobian_det_;

    static std::string build_message(const std::string& msg, Real det) {
        return msg + " (Jacobian determinant: " + std::to_string(det) + ")";
    }
};

// ============================================================================
// Exception Throwing Macros
// ============================================================================

/**
 * @brief Throw exception with automatic source location
 */
#define FE_THROW(ExceptionType, message) \
    throw ExceptionType(message, __FILE__, __LINE__, __FUNCTION__)

/**
 * @brief Conditional throw with source location (3 args: condition, ExceptionType, message)
 */
#define FE_THROW_IF_3(condition, ExceptionType, message) \
    do { \
        if (FE_UNLIKELY(condition)) { \
            FE_THROW(ExceptionType, message); \
        } \
    } while(0)

/**
 * @brief Conditional throw with FEException (2 args: condition, message)
 */
#define FE_THROW_IF_2(condition, message) \
    do { \
        if (FE_UNLIKELY(condition)) { \
            FE_THROW(FEException, message); \
        } \
    } while(0)

/**
 * @brief Helper macro to select between 2 and 3 argument versions
 */
#define FE_THROW_IF_SELECT(_1, _2, _3, NAME, ...) NAME

/**
 * @brief Conditional throw with source location
 *
 * Can be called with 2 arguments (condition, message) using FEException,
 * or 3 arguments (condition, ExceptionType, message) for specific exception types.
 */
#define FE_THROW_IF(...) \
    FE_THROW_IF_SELECT(__VA_ARGS__, FE_THROW_IF_3, FE_THROW_IF_2)(__VA_ARGS__)

/**
 * @brief Check and throw InvalidArgumentException
 */
#define FE_CHECK_ARG(condition, message) \
    FE_THROW_IF(!(condition), InvalidArgumentException, message)

/**
 * @brief Check for null pointers
 */
#define FE_CHECK_NOT_NULL(ptr, name) \
    FE_THROW_IF((ptr) == nullptr, InvalidArgumentException, \
                std::string(name) + " is null")

/**
 * @brief Check index bounds
 */
#define FE_CHECK_INDEX(index, size) \
    FE_THROW_IF((index) < 0 || (index) >= (size), \
                InvalidArgumentException, \
                "Index " + std::to_string(index) + " out of bounds [0, " + \
                std::to_string(size) + ")")

/**
 * @brief Throw NotImplementedException
 */
#define FE_NOT_IMPLEMENTED(feature) \
    FE_THROW(NotImplementedException, feature)

// ============================================================================
// MPI Error Handling
// ============================================================================

#if FE_HAS_MPI

/**
 * @brief MPI-aware error handler
 *
 * Ensures all ranks abort cleanly when an exception occurs
 */
class MPIErrorHandler {
public:
    static void abort_all_ranks(const FEException& e) {
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        // Print error on rank that threw
        if (rank == e.mpi_rank() || e.mpi_rank() < 0) {
            std::cerr << e.what() << std::endl;
        }

        // Ensure all ranks abort
        MPI_Abort(MPI_COMM_WORLD, static_cast<int>(e.status()));
    }

    /**
     * @brief Install MPI error handler
     */
    static void install() {
        std::set_terminate([]() {
            try {
                std::rethrow_exception(std::current_exception());
            } catch (const FEException& e) {
                abort_all_ranks(e);
            } catch (const std::exception& e) {
                std::cerr << "Unhandled exception: " << e.what() << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            } catch (...) {
                std::cerr << "Unknown exception" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        });
    }
};

#endif // FE_HAS_MPI

} // namespace FE
} // namespace svmp

#endif // SVMP_FE_EXCEPTION_H
