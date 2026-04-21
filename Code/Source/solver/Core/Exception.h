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

#include "Core/PlatformSupport.h"

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

class ExceptionFormatter final {
public:
    static std::string format(const ExceptionContext& context,
                              const std::string& message,
                              const char* subsystem_label)
    {
        std::ostringstream oss;

        oss << "[" << ((subsystem_label == nullptr) ? "Exception" : subsystem_label)
            << "] " << status_code_to_string(context.status_code());
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
};

class ExceptionBase;

class ExceptionRuntime final {
public:
    static int query_mpi_rank() noexcept
    {
        return PlatformSupport::query_mpi_rank();
    }

    static StackTrace capture_stack_trace()
    {
        return PlatformSupport::capture_stack_trace();
    }

    static void finalize_mpi_if_needed() noexcept
    {
        PlatformSupport::finalize_mpi_if_needed();
    }

    static void install_terminate_handler();

    static void report_unhandled_exception(const std::exception& exception)
    {
        std::cerr << exception.what() << std::endl;
    }

    static void abort_mpi_if_needed(int exit_code) noexcept
    {
        PlatformSupport::abort_mpi_if_needed(exit_code);
    }
};

class ExceptionBase : public std::exception {
public:
    const char* what() const noexcept override { return what_.c_str(); }
    StatusCode status_code() const noexcept { return context_.status_code(); }
    const std::string& message() const noexcept { return message_; }
    const ExceptionContext& context() const noexcept { return context_; }

    void add_context(const std::string& context)
    {
        message_ = context + "\n  -> " + message_;
        rebuild_what(subsystem_label());
    }

    virtual ~ExceptionBase() noexcept = default;

protected:
    ExceptionBase(std::string message, StatusCode status, const char* file,
                  int line, const char* function)
        : message_(std::move(message))
    {
        context_.set_status_code(status);
        context_.set_source_location(file, line, function);
        context_.set_mpi_rank(ExceptionRuntime::query_mpi_rank());
#if SVMP_EXCEPTION_DEBUG_MODE
        context_.set_stack_trace(ExceptionRuntime::capture_stack_trace());
#endif
    }

    virtual const char* subsystem_label() const noexcept = 0;

    void rebuild_what(const char* subsystem_label)
    {
        what_ = ExceptionFormatter::format(context_, message_, subsystem_label);
    }

    std::string message_;
    ExceptionContext context_;
    std::string what_;
};

class CoreException : public ExceptionBase {
public:
    CoreException(const std::string& message,
                  StatusCode status = StatusCode::Unknown,
                  const char* file = "",
                  int line = 0,
                  const char* function = "")
        : ExceptionBase(message, status, file, line, function)
    {
        rebuild_what(subsystem_label());
    }

protected:
    const char* subsystem_label() const noexcept override
    {
        return "Core Exception";
    }
};

class ParseException : public CoreException {
public:
    ParseException(const std::string& message,
                   const char* file = "",
                   int line = 0,
                   const char* function = "")
        : CoreException(message, StatusCode::ParseError, file, line, function)
    {
    }
};

class DependencyException : public CoreException {
public:
    DependencyException(const std::string& message,
                        const char* file = "",
                        int line = 0,
                        const char* function = "")
        : CoreException(message, StatusCode::DependencyError, file, line,
                        function)
    {
    }
};

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

template <class ExceptionT, class... Args>
[[noreturn]] void raise(SourceLocation location, Args&&... args)
{
    throw ExceptionT(std::forward<Args>(args)..., location.file, location.line,
                     location.function);
}

template <class ExceptionT, class... Args>
void check(bool condition, SourceLocation location, Args&&... args)
{
    if (!condition) {
        raise<ExceptionT>(location, std::forward<Args>(args)...);
    }
}

template <class ExceptionT, class... Args>
void check_arg(bool condition, SourceLocation location, Args&&... args)
{
    if (!condition) {
        raise<ExceptionT>(location, std::forward<Args>(args)...);
    }
}

template <class ExceptionT, class PointerT, class... Args>
PointerT* check_not_null(PointerT* ptr, SourceLocation location, Args&&... args)
{
    if (ptr == nullptr) {
        raise<ExceptionT>(location, std::forward<Args>(args)...);
    }
    return ptr;
}

} // namespace svmp

#define SVMP_HERE ::svmp::SourceLocation{__FILE__, __LINE__, __func__}

#endif // SVMP_CORE_EXCEPTION_H
