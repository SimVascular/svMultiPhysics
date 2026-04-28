// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_CORE_PLATFORM_SUPPORT_INL
#define SVMP_CORE_PLATFORM_SUPPORT_INL

#ifndef SVMP_CORE_EXCEPTION_INCLUDE_PLATFORM_SUPPORT
#error "PlatformSupport.inl is private; include Core/Exception.h instead."
#endif

#include <mpi.h>

#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <string>

#if defined(_WIN32)
#  include <windows.h>
#elif defined(__GNUC__) || defined(__clang__)
#  include <cxxabi.h>
#  include <execinfo.h>
#endif

namespace svmp {

inline PlatformCapabilities PlatformSupport::capabilities() noexcept
{
#if defined(_WIN32)
    return PlatformCapabilities{true, false, false};
#elif defined(__GNUC__) || defined(__clang__)
    return PlatformCapabilities{true, true, true};
#else
    return PlatformCapabilities{false, false, false};
#endif
}

inline int PlatformSupport::query_mpi_rank() noexcept
{
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized) {
        return -1;
    }

    int finalized = 0;
    MPI_Finalized(&finalized);
    if (finalized) {
        return -1;
    }

    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

inline std::string PlatformSupport::demangle_symbol(const char* symbol)
{
    if (symbol == nullptr) {
        return std::string();
    }

#if defined(__GNUC__) || defined(__clang__)
    int status = 0;
    char* demangled = abi::__cxa_demangle(symbol, nullptr, nullptr, &status);
    if (status == 0 && demangled != nullptr) {
        std::string result(demangled);
        std::free(demangled);
        return result;
    }
#endif

    return std::string(symbol);
}

inline StackTrace PlatformSupport::capture_stack_trace()
{
    StackTrace trace;

#if defined(_WIN32)
    constexpr unsigned long max_frames = 32;
    void* frames[max_frames] = {};
    const USHORT count = CaptureStackBackTrace(0, max_frames, frames, nullptr);

    for (USHORT i = 1; i < count; ++i) {
        const std::uintptr_t address =
            reinterpret_cast<std::uintptr_t>(frames[i]);
        std::ostringstream symbol;
        symbol << "0x" << std::hex << address;
        trace.add_frame(StackTraceFrame(symbol.str(), std::string(),
                                        std::string(), 0, address));
    }
#elif defined(__GNUC__) || defined(__clang__)
    constexpr int max_frames = 32;
    void* frames[max_frames] = {};
    const int count = backtrace(frames, max_frames);
    char** symbols = backtrace_symbols(frames, count);

    if (symbols != nullptr) {
        for (int i = 2; i < count; ++i) {
            const std::uintptr_t address =
                reinterpret_cast<std::uintptr_t>(frames[i]);
            std::string symbol(symbols[i]);

            const std::size_t start = symbol.find('(');
            const std::size_t end = symbol.find('+', start);
            if (start != std::string::npos && end != std::string::npos &&
                end > start + 1) {
                const std::string mangled =
                    symbol.substr(start + 1, end - start - 1);
                const std::string demangled =
                    PlatformSupport::demangle_symbol(mangled.c_str());
                if (!demangled.empty()) {
                    symbol.replace(start + 1, end - start - 1, demangled);
                }
            }

            trace.add_frame(
                StackTraceFrame(symbol, std::string(), std::string(), 0, address));
        }
        std::free(symbols);
    }
#endif

    return trace;
}

inline void PlatformSupport::finalize_mpi_if_needed() noexcept
{
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized) {
        return;
    }

    int finalized = 0;
    MPI_Finalized(&finalized);
    if (!finalized) {
        MPI_Finalize();
    }
}

inline void PlatformSupport::abort_mpi_if_needed(int exit_code) noexcept
{
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized) {
        return;
    }

    int finalized = 0;
    MPI_Finalized(&finalized);
    if (finalized) {
        return;
    }

    int size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size > 1) {
        MPI_Abort(MPI_COMM_WORLD, exit_code);
    }
}

} // namespace svmp

#endif // SVMP_CORE_PLATFORM_SUPPORT_INL
