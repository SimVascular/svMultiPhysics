// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_CORE_EXCEPTION_H
#define SVMP_CORE_EXCEPTION_H

#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#if !defined(NDEBUG) || defined(DEBUG) || defined(_DEBUG)
#  define SVMP_EXCEPTION_DEBUG_MODE 1
#else
#  define SVMP_EXCEPTION_DEBUG_MODE 0
#endif

namespace svmp {

enum class StatusCode : std::uint8_t {
    Success = 0,
    InvalidArgument,
    InvalidState,
    ParseError,
    IOError,
    ResourceExhausted,
    DependencyError,
    MPIError,
    NotImplemented,
    UnsupportedOperation,
    InternalError,
    Unknown = 255
};

inline const char* status_code_to_string(StatusCode status) noexcept
{
    switch (status) {
        case StatusCode::Success:
            return "Success";
        case StatusCode::InvalidArgument:
            return "Invalid argument";
        case StatusCode::InvalidState:
            return "Invalid state";
        case StatusCode::ParseError:
            return "Parse error";
        case StatusCode::IOError:
            return "I/O error";
        case StatusCode::ResourceExhausted:
            return "Resource exhausted";
        case StatusCode::DependencyError:
            return "Dependency error";
        case StatusCode::MPIError:
            return "MPI error";
        case StatusCode::NotImplemented:
            return "Not implemented";
        case StatusCode::UnsupportedOperation:
            return "Unsupported operation";
        case StatusCode::InternalError:
            return "Internal error";
        default:
            return "Unknown error";
    }
}

struct SourceLocation {
    const char* file;
    int line;
    const char* function;
};

class StackTraceFrame final {
public:
    StackTraceFrame() = default;

    StackTraceFrame(std::string symbol, std::string module, std::string file,
                    int line, std::uintptr_t address) noexcept
        : symbol_(std::move(symbol)),
          module_(std::move(module)),
          file_(std::move(file)),
          line_(line),
          address_(address)
    {
    }

    const std::string& symbol() const noexcept { return symbol_; }
    const std::string& module() const noexcept { return module_; }
    const std::string& file() const noexcept { return file_; }
    int line() const noexcept { return line_; }
    std::uintptr_t address() const noexcept { return address_; }

private:
    std::string symbol_;
    std::string module_;
    std::string file_;
    int line_ = 0;
    std::uintptr_t address_ = 0;
};

class StackTrace final {
public:
    bool empty() const noexcept { return frames_.empty(); }
    std::size_t size() const noexcept { return frames_.size(); }
    const std::vector<StackTraceFrame>& frames() const noexcept { return frames_; }
    void add_frame(StackTraceFrame frame) { frames_.push_back(std::move(frame)); }

private:
    std::vector<StackTraceFrame> frames_;
};

struct PlatformCapabilities final {
    bool has_stack_trace;
    bool has_symbol_resolution;
    bool has_demangling;
};

class PlatformSupport final {
public:
    static PlatformCapabilities capabilities() noexcept;
    static int query_mpi_rank() noexcept;
    static StackTrace capture_stack_trace();
    static std::string demangle_symbol(const char* symbol);
    static void finalize_mpi_if_needed() noexcept;
    static void abort_mpi_if_needed(int exit_code) noexcept;
};
} // namespace svmp

#define SVMP_CORE_EXCEPTION_INCLUDE_PLATFORM_SUPPORT
#include "PlatformSupport.inl"
#undef SVMP_CORE_EXCEPTION_INCLUDE_PLATFORM_SUPPORT

namespace svmp {

class ExceptionContext final {
public:
    StatusCode status_code() const noexcept { return status_code_; }
    const std::string& file() const noexcept { return file_; }
    int line() const noexcept { return line_; }
    const std::string& function() const noexcept { return function_; }
    int mpi_rank() const noexcept { return mpi_rank_; }
    const StackTrace& stack_trace() const noexcept { return stack_trace_; }

    void set_status_code(StatusCode status_code) noexcept
    {
        status_code_ = status_code;
    }

    void set_source_location(const char* file, int line, const char* function)
    {
        file_ = (file == nullptr) ? std::string() : std::string(file);
        line_ = line;
        function_ =
            (function == nullptr) ? std::string() : std::string(function);
    }

    void set_mpi_rank(int mpi_rank) noexcept { mpi_rank_ = mpi_rank; }
    void set_stack_trace(StackTrace stack_trace)
    {
        stack_trace_ = std::move(stack_trace);
    }

private:
    StatusCode status_code_ = StatusCode::Unknown;
    std::string file_;
    int line_ = 0;
    std::string function_;
    int mpi_rank_ = -1;
    StackTrace stack_trace_;
};

namespace ExceptionFormatter {

inline std::string format(const ExceptionContext& context,
                          const std::string& message,
                          std::string_view subsystem_label = "Exception")
{
    if (subsystem_label.empty()) {
        subsystem_label = "Exception";
    }

    std::ostringstream oss;

    oss << "[" << subsystem_label << "] "
        << status_code_to_string(context.status_code());
    if (context.mpi_rank() >= 0) {
        oss << " (Rank " << context.mpi_rank() << ")";
    }
    oss << "\n";

    if (!context.file().empty()) {
        oss << "  Location: " << context.file() << ":" << context.line();
        if (!context.function().empty()) {
            oss << " in " << context.function() << "()";
        }
        oss << "\n";
    }

    oss << "  Message: " << message << "\n";

    if (!context.stack_trace().empty()) {
        oss << "  Stack trace:\n";
        std::size_t frame_index = 0;
        for (const auto& frame : context.stack_trace().frames()) {
            oss << "    #" << frame_index++ << " ";
            if (!frame.symbol().empty()) {
                oss << frame.symbol();
            } else {
                std::ostringstream address;
                address << "0x" << std::hex << frame.address();
                oss << address.str();
            }

            if (!frame.module().empty()) {
                oss << " [" << frame.module() << "]";
            }

            if (!frame.file().empty()) {
                oss << " (" << frame.file();
                if (frame.line() > 0) {
                    oss << ":" << frame.line();
                }
                oss << ")";
            }

            oss << "\n";
        }
    }

    return oss.str();
}

} // namespace ExceptionFormatter

class ExceptionBase;

namespace ExceptionRuntime {

    inline int query_mpi_rank() noexcept
    {
        return PlatformSupport::query_mpi_rank();
    }

    inline StackTrace capture_stack_trace()
    {
        return PlatformSupport::capture_stack_trace();
    }

    inline void finalize_mpi_if_needed() noexcept
    {
        PlatformSupport::finalize_mpi_if_needed();
    }

    void install_terminate_handler();

    inline void report_unhandled_exception(const std::exception& exception)
    {
        std::cerr << exception.what() << std::endl;
    }

    inline void abort_mpi_if_needed(int exit_code) noexcept
    {
        PlatformSupport::abort_mpi_if_needed(exit_code);
    }

} // namespace ExceptionRuntime

class ExceptionBase : public std::exception {
public:
    const char* what() const noexcept override { return what_.c_str(); }
    StatusCode status_code() const noexcept { return context_.status_code(); }
    const std::string& message() const noexcept { return message_; }
    const ExceptionContext& context() const noexcept { return context_; }

    void add_context(const std::string& context)
    {
        message_ = context + "\n  -> " + message_;
        rebuild_what();
    }

    /// @brief Record the originating source location and refresh what().
    ///
    /// @details Called by raise() after construction, so exception constructors
    /// do not need to accept (and forward) the file/line/function themselves.
    void set_source_location(const SourceLocation& location)
    {
        context_.set_source_location(location.file, location.line,
                                     location.function);
        rebuild_what();
    }

    virtual ~ExceptionBase() noexcept = default;

protected:
    ExceptionBase(std::string message, StatusCode status,
                  std::string_view subsystem_label, const char* file = "",
                  int line = 0, const char* function = "")
        : message_(std::move(message)),
          subsystem_label_(subsystem_label.empty() ? std::string_view("Exception")
                                                   : subsystem_label)
    {
        context_.set_status_code(status);
        context_.set_source_location(file, line, function);
        context_.set_mpi_rank(ExceptionRuntime::query_mpi_rank());
#if SVMP_EXCEPTION_DEBUG_MODE
        context_.set_stack_trace(ExceptionRuntime::capture_stack_trace());
#endif
        rebuild_what();
    }

    void rebuild_what()
    {
        what_ = ExceptionFormatter::format(context_, message_, subsystem_label_);
    }

    std::string message_;
    ExceptionContext context_;
    std::string_view subsystem_label_;
    std::string what_;
};

class CoreException : public ExceptionBase {
public:
    CoreException(const std::string& message,
                  StatusCode status = StatusCode::Unknown,
                  const char* file = "",
                  int line = 0,
                  const char* function = "")
        : ExceptionBase(message, status, "Core Exception", file, line, function)
    {
    }
};

/**
 * @brief Define a simple, message-only exception type in one line.
 *
 * @details Expands to a class @p Name deriving from @p Base with a single
 * `explicit Name(const std::string& message)` constructor that records @p Status.
 * Use it for exceptions that carry only a message; write the class by hand when it
 * needs extra structured context (members and accessors). The source location is
 * stamped by raise(), so no file/line/function constructor is needed.
 *
 * @p Base must be an ExceptionBase-derived type whose constructor accepts
 * `(const std::string&, StatusCode)` (CoreException, FEException, and the
 * subsystem bases do).
 *
 * @code
 * SVMP_DEFINE_EXCEPTION(ParseException, CoreException, StatusCode::ParseError);
 * @endcode
 */
#define SVMP_DEFINE_EXCEPTION(Name, Base, Status)                              \
    class Name : public Base {                                                 \
    public:                                                                    \
        explicit Name(const std::string& message)                             \
            : Base(message, (Status))                                          \
        {                                                                      \
        }                                                                      \
    }

/// @brief A parsing or input-format error.
SVMP_DEFINE_EXCEPTION(ParseException, CoreException, StatusCode::ParseError);

/// @brief A required dependency is missing or failed to load.
SVMP_DEFINE_EXCEPTION(DependencyException, CoreException,
                      StatusCode::DependencyError);

/// @brief A requested operation or feature is not implemented.
///
/// @details The default exception raised by not_implemented().
SVMP_DEFINE_EXCEPTION(NotImplementedException, CoreException,
                      StatusCode::NotImplemented);

/// @brief An index is outside its valid range.
///
/// @details The default exception raised by check_index(); the status code is
/// InvalidArgument because an out-of-range index is a caller error.
SVMP_DEFINE_EXCEPTION(IndexOutOfRangeException, CoreException,
                      StatusCode::InvalidArgument);

inline void ExceptionRuntime::install_terminate_handler()
{
    std::set_terminate([]() {
        try {
            const std::exception_ptr current = std::current_exception();
            if (current != nullptr) {
                std::rethrow_exception(current);
            }
        } catch (const std::exception& exception) {
            ExceptionRuntime::report_unhandled_exception(exception);
        } catch (...) {
            std::cerr << "[Unhandled Exception] Unknown non-std exception"
                      << std::endl;
        }

        ExceptionRuntime::abort_mpi_if_needed(EXIT_FAILURE);
        std::abort();
    });
}

/**
 * @brief A diagnostic message bundled with the source location where it was written.
 *
 * @details The core helpers (raise(), check(), throw_if(), check_not_null()) take
 * a Diagnostic in place of an explicit source location. Its file/line/function
 * arguments default to the compiler builtins __builtin_FILE()/__builtin_LINE()/
 * __builtin_FUNCTION(), which capture the caller's location, so a string literal
 * or std::string passed at the call site is implicitly wrapped into a Diagnostic
 * that records exactly where the call appears -- callers pass no explicit location:
 * @code
 * svmp::check<MyException>(ptr != nullptr, "pointer must not be null");
 * @endcode
 */
class Diagnostic {
public:
    /**
     * @brief Wrap a message, capturing the caller's source location by default.
     * @param message The diagnostic message.
     * @param file Source file; defaults to the caller's via __builtin_FILE().
     * @param line Source line; defaults to the caller's via __builtin_LINE().
     * @param function Function; defaults to the caller's via __builtin_FUNCTION().
     */
    Diagnostic(const char* message,
               const char* file = __builtin_FILE(),
               int line = __builtin_LINE(),
               const char* function = __builtin_FUNCTION())
        : message_(message), location_{file, line, function}
    {
    }
    /**
     * @brief Wrap a message, capturing the caller's source location by default.
     * @param message The diagnostic message.
     * @param file Source file; defaults to the caller's via __builtin_FILE().
     * @param line Source line; defaults to the caller's via __builtin_LINE().
     * @param function Function; defaults to the caller's via __builtin_FUNCTION().
     */
    Diagnostic(std::string message,
               const char* file = __builtin_FILE(),
               int line = __builtin_LINE(),
               const char* function = __builtin_FUNCTION())
        : message_(std::move(message)), location_{file, line, function}
    {
    }

    /**
     * @brief The diagnostic message.
     * @return The stored message.
     */
    const std::string& message() const noexcept { return message_; }
    /**
     * @brief The source location captured when the Diagnostic was constructed.
     * @return The stored source location.
     */
    const SourceLocation& location() const noexcept { return location_; }

private:
    std::string message_;
    SourceLocation location_;
};

/**
 * @brief Construct @p ExceptionT from the diagnostic message and @p args, stamp
 * the source location, and throw it.
 *
 * @details @p diagnostic carries the message and the source location captured at
 * the call site; @p args are forwarded to the exception constructor after the
 * message. The location is recorded via ExceptionBase::set_source_location(), so
 * exception types never need a file/line/function constructor -- a `(message)`
 * (plus any structured-context) constructor is enough.
 */
template <class ExceptionT, class... Args>
[[noreturn]] void raise(Diagnostic diagnostic, Args&&... args)
{
    static_assert(std::is_base_of_v<ExceptionBase, ExceptionT>,
                  "raise<>() requires an svmp::ExceptionBase-derived exception type");
    ExceptionT exception(diagnostic.message(), std::forward<Args>(args)...);
    exception.set_source_location(diagnostic.location());
    throw exception;
}

/**
 * @brief Raise @p ExceptionT when @p condition is false (a required condition).
 *
 * @details The general success-condition check, used for argument, state, and
 * invariant validation: `check<E>(ptr != nullptr, "...")`. throw_if() is the
 * logical inverse -- it raises when its condition is true.
 */
template <class ExceptionT, class... Args>
void check(bool condition, Diagnostic diagnostic, Args&&... args)
{
    if (!condition) {
        raise<ExceptionT>(std::move(diagnostic), std::forward<Args>(args)...);
    }
}

/**
 * @brief Raise @p ExceptionT when @p ptr is null.
 */
template <class ExceptionT, class PointerT, class... Args>
void check_not_null(PointerT ptr, Diagnostic diagnostic, Args&&... args)
{
    if (ptr == nullptr) {
        raise<ExceptionT>(std::move(diagnostic), std::forward<Args>(args)...);
    }
}

/**
 * @brief Raise @p ExceptionT when @p condition is true.
 *
 * @details The logical inverse of check(): check() raises when its condition is
 * false (a required condition); throw_if() raises when its condition is true (a
 * failure condition). The two are not interchangeable.
 */
template <class ExceptionT, class... Args>
void throw_if(bool condition, Diagnostic diagnostic, Args&&... args)
{
    if (condition) {
        raise<ExceptionT>(std::move(diagnostic), std::forward<Args>(args)...);
    }
}

/**
 * @brief Raise an exception when @p index is outside [0, @p size).
 *
 * @details @p ExceptionT defaults to IndexOutOfRangeException; supply a different
 * type only when a subsystem needs its own exception. The bounds message and the
 * source location are generated automatically.
 */
template <class ExceptionT = IndexOutOfRangeException, class IndexT, class SizeT>
void check_index(IndexT index, SizeT size,
                 const char* file = __builtin_FILE(),
                 int line = __builtin_LINE(),
                 const char* function = __builtin_FUNCTION())
{
    const long long index_value = static_cast<long long>(index);
    const long long size_value = static_cast<long long>(size);
    check<ExceptionT>(
        index_value >= 0 && index_value < size_value,
        Diagnostic("Index " + std::to_string(index_value) + " out of bounds [0, " +
                       std::to_string(size_value) + ")",
                   file, line, function));
}

/**
 * @brief Raise an exception reporting an unimplemented feature, with a message.
 *
 * @details @p ExceptionT defaults to NotImplementedException, so most call sites
 * pass only a message:
 * @code
 * svmp::not_implemented("GPU assembly is not supported");
 * @endcode
 * Pass a different exception type explicitly only when a subsystem needs one.
 */
template <class ExceptionT = NotImplementedException>
[[noreturn]] void not_implemented(std::string message,
                                  const char* file = __builtin_FILE(),
                                  int line = __builtin_LINE(),
                                  const char* function = __builtin_FUNCTION())
{
    raise<ExceptionT>(Diagnostic(std::move(message), file, line, function));
}

} // namespace svmp

#if SVMP_EXCEPTION_DEBUG_MODE
#define SVMP_DEBUG_CHECK(ExceptionT, condition, ...)                         \
    do {                                                                     \
        if (!(condition)) {                                                  \
            ::svmp::raise<ExceptionT>(__VA_ARGS__);                          \
        }                                                                    \
    } while (false)
#else
#define SVMP_DEBUG_CHECK(ExceptionT, condition, ...)                         \
    do {                                                                     \
    } while (false)
#endif

#endif // SVMP_CORE_EXCEPTION_H
