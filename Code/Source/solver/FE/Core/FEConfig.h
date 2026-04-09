// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_CONFIG_H
#define SVMP_FE_CONFIG_H

/**
 * @file FEConfig.h
 * @brief Compile-time configuration and feature settings for the FE library
 *
 * This header centralizes all compile-time configuration options,
 * feature flags, and platform-specific settings for the finite element
 * library. Settings can be overridden via CMake or compiler flags.
 */

#include "Types.h"
#include <cstddef>

// ============================================================================
// Build Configuration Detection
// ============================================================================

// Debug/Release mode detection
#if !defined(NDEBUG) || defined(DEBUG) || defined(_DEBUG)
    #define FE_DEBUG_MODE 1
#else
    #define FE_DEBUG_MODE 0
#endif

#if FE_DEBUG_MODE
#include <cassert>
#endif

// MPI support detection.
//
// FE can be built in two modes:
//  - Integrated with Mesh: Mesh defines `MESH_HAS_MPI` when compiled with MPI.
//  - Standalone FE with MPI: FE's CMake sets `MESH_ENABLE_MPI` for compatibility.
//
// Normalize to a numeric macro suitable for `#if FE_HAS_MPI` use.
#ifdef FE_HAS_MPI
#  undef FE_HAS_MPI
#endif
#if defined(MESH_HAS_MPI) || defined(MESH_ENABLE_MPI)
#  define FE_HAS_MPI 1
#else
#  define FE_HAS_MPI 0
#endif

// OpenMP support detection
#ifdef _OPENMP
    #define FE_HAS_OPENMP 1
    #include <omp.h>
#else
    #define FE_HAS_OPENMP 0
#endif

namespace svmp {
namespace FE {
namespace config {

// ============================================================================
// Dimension Configuration
// ============================================================================

/**
 * @brief Maximum spatial dimension supported
 *
 * Controls template instantiations and static array sizes.
 * Can be overridden at compile time with -DFE_MAX_DIM=N
 */
#ifndef FE_MAX_DIM
    constexpr int MAX_SPATIAL_DIM = 3;
#else
    constexpr int MAX_SPATIAL_DIM = FE_MAX_DIM;
#endif

/**
 * @brief Enable specific dimensions
 *
 * Controls which template instantiations are generated.
 * Disabling unused dimensions reduces binary size.
 */
#ifndef FE_ENABLE_1D
    #define FE_ENABLE_1D 1
#endif

#ifndef FE_ENABLE_2D
    #define FE_ENABLE_2D 1
#endif

#ifndef FE_ENABLE_3D
    #define FE_ENABLE_3D 1
#endif

// ============================================================================
// Precision Configuration
// ============================================================================

/**
 * @brief Floating point precision settings
 */
enum class Precision {
    Single,    // 32-bit float
    Double,    // 64-bit double (default)
    Quad       // 128-bit quad precision (if available)
};

/**
 * @brief Default precision for FE computations
 *
 * Can be overridden with -DFE_USE_SINGLE_PRECISION
 */
#ifdef FE_USE_SINGLE_PRECISION
    constexpr Precision DEFAULT_PRECISION = Precision::Single;
    using DefaultReal = float;
#elif defined(FE_USE_QUAD_PRECISION)
    constexpr Precision DEFAULT_PRECISION = Precision::Quad;
    using DefaultReal = long double;
#else
    constexpr Precision DEFAULT_PRECISION = Precision::Double;
    using DefaultReal = double;
#endif

// ============================================================================
// Backend Configuration
// ============================================================================

/**
 * @brief Linear algebra backend selection
 */
enum class Backend {
    None,      // No backend (testing only)
    Trilinos,  // Trilinos (Epetra/Tpetra)
    PETSc,     // PETSc
    Eigen,     // Eigen (for small/local problems)
    Custom     // User-provided backend
};

/**
 * @brief Default backend selection
 *
 * Set via CMake based on available libraries
 */
#if defined(FE_USE_TRILINOS)
    constexpr Backend DEFAULT_BACKEND = Backend::Trilinos;
    #define FE_HAS_TRILINOS 1
#elif defined(FE_USE_PETSC)
    constexpr Backend DEFAULT_BACKEND = Backend::PETSc;
    #define FE_HAS_PETSC 1
#elif defined(FE_USE_EIGEN)
    constexpr Backend DEFAULT_BACKEND = Backend::Eigen;
    #define FE_HAS_EIGEN 1
#else
    constexpr Backend DEFAULT_BACKEND = Backend::None;
#endif

// ============================================================================
// Performance Configuration
// ============================================================================

/**
 * @brief Cache line size in bytes
 *
 * Used for memory alignment and padding to avoid false sharing
 */
#ifndef FE_CACHE_LINE_SIZE
    constexpr std::size_t CACHE_LINE_SIZE = 64;
#else
    constexpr std::size_t CACHE_LINE_SIZE = FE_CACHE_LINE_SIZE;
#endif

/**
 * @brief Enable SIMD vectorization
 *
 * Controls whether vectorized implementations are used
 */
#ifndef FE_ENABLE_VECTORIZATION
    #if defined(__AVX512F__)
        #define FE_VECTORIZATION_LEVEL 512
        #define FE_ENABLE_VECTORIZATION 1
    #elif defined(__AVX2__)
        #define FE_VECTORIZATION_LEVEL 256
        #define FE_ENABLE_VECTORIZATION 1
    #elif defined(__SSE4_2__)
        #define FE_VECTORIZATION_LEVEL 128
        #define FE_ENABLE_VECTORIZATION 1
    #else
        #define FE_VECTORIZATION_LEVEL 0
        #define FE_ENABLE_VECTORIZATION 0
    #endif
#endif

/**
 * @brief Enable aggressive loop unrolling
 */
#ifndef FE_ENABLE_LOOP_UNROLL
    #define FE_ENABLE_LOOP_UNROLL 1
#endif

/**
 * @brief Assembly kernel caching
 *
 * Cache element matrices/vectors to reduce recomputation
 */
#ifndef FE_ENABLE_ASSEMBLY_CACHE
    #define FE_ENABLE_ASSEMBLY_CACHE 1
#endif

/**
 * @brief Maximum elements to process in vectorized batches
 */
constexpr std::size_t VECTOR_BATCH_SIZE = 32;

// ============================================================================
// Memory Configuration
// ============================================================================

/**
 * @brief Memory allocation strategies
 */
enum class AllocStrategy {
    Standard,     // Standard new/delete
    Pooled,       // Memory pool allocator
    Aligned,      // Aligned allocation for SIMD
    Custom        // User-provided allocator
};

/**
 * @brief Default memory allocation strategy
 */
#ifdef FE_USE_MEMORY_POOL
    constexpr AllocStrategy DEFAULT_ALLOC_STRATEGY = AllocStrategy::Pooled;
#elif FE_ENABLE_VECTORIZATION
    constexpr AllocStrategy DEFAULT_ALLOC_STRATEGY = AllocStrategy::Aligned;
#else
    constexpr AllocStrategy DEFAULT_ALLOC_STRATEGY = AllocStrategy::Standard;
#endif

/**
 * @brief Stack array size threshold
 *
 * Arrays smaller than this use stack allocation
 */
constexpr std::size_t STACK_ARRAY_THRESHOLD = 256;

// ============================================================================
// Debug/Diagnostic Configuration
// ============================================================================

/**
 * @brief Bounds checking level
 */
enum class BoundsCheckLevel {
    None,         // No bounds checking (release mode)
    Basic,        // Basic index validation
    Full          // Full bounds and consistency checking
};

#if FE_DEBUG_MODE
    constexpr BoundsCheckLevel BOUNDS_CHECK_LEVEL = BoundsCheckLevel::Full;
    constexpr bool ENABLE_ASSERTIONS = true;
    constexpr bool ENABLE_NAN_CHECK = true;
#else
    constexpr BoundsCheckLevel BOUNDS_CHECK_LEVEL = BoundsCheckLevel::None;
    constexpr bool ENABLE_ASSERTIONS = false;
    constexpr bool ENABLE_NAN_CHECK = false;
#endif

/**
 * @brief Performance profiling
 */
#ifndef FE_ENABLE_PROFILING
    #define FE_ENABLE_PROFILING FE_DEBUG_MODE
#endif

/**
 * @brief Detailed timing of assembly operations
 */
#ifndef FE_ENABLE_TIMING
    #define FE_ENABLE_TIMING 0
#endif

// ============================================================================
// Feature Toggles
// ============================================================================


/**
 * @brief Maximum polynomial order supported
 */
#ifndef FE_MAX_POLYNOMIAL_ORDER
    constexpr int MAX_POLYNOMIAL_ORDER = 10;
#else
    constexpr int MAX_POLYNOMIAL_ORDER = FE_MAX_POLYNOMIAL_ORDER;
#endif

// ============================================================================
// Quadrature Configuration
// ============================================================================


/**
 * @brief Cache quadrature rules
 */
#ifndef FE_CACHE_QUADRATURE
    #define FE_CACHE_QUADRATURE 1
#endif

// ============================================================================
// Compiler Hints and Attributes
// ============================================================================

// Likely/unlikely branch hints
#if defined(__GNUC__) || defined(__clang__)
    #define FE_LIKELY(x)   __builtin_expect(!!(x), 1)
    #define FE_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
    #define FE_LIKELY(x)   (x)
    #define FE_UNLIKELY(x) (x)
#endif

// Force inline
#if defined(_MSC_VER)
    #define FE_ALWAYS_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
    #define FE_ALWAYS_INLINE __attribute__((always_inline)) inline
#else
    #define FE_ALWAYS_INLINE inline
#endif

// Restrict pointer aliasing
#if defined(__GNUC__) || defined(__clang__)
    #define FE_RESTRICT __restrict__
#elif defined(_MSC_VER)
    #define FE_RESTRICT __restrict
#else
    #define FE_RESTRICT
#endif

// Memory alignment
#if defined(__GNUC__) || defined(__clang__)
    #define FE_ALIGN(x) __attribute__((aligned(x)))
#elif defined(_MSC_VER)
    #define FE_ALIGN(x) __declspec(align(x))
#else
    #define FE_ALIGN(x)
#endif

// ============================================================================
// Assertion Macros
// ============================================================================

#if FE_DEBUG_MODE
    #define FE_ASSERT(cond) assert(cond)
    #define FE_ASSERT_MSG(cond, msg) assert((cond) && (msg))
#else
    #define FE_ASSERT(cond) ((void)0)
    #define FE_ASSERT_MSG(cond, msg) ((void)0)
#endif

// Bounds checking macro
#if BOUNDS_CHECK_LEVEL == 2  // Full
    #define FE_CHECK_BOUNDS(index, size) \
        FE_ASSERT_MSG((index) >= 0 && (index) < (size), "Index out of bounds")
#elif BOUNDS_CHECK_LEVEL == 1  // Basic
    #define FE_CHECK_BOUNDS(index, size) \
        FE_ASSERT((index) < (size))
#else
    #define FE_CHECK_BOUNDS(index, size) ((void)0)
#endif

// ============================================================================
// Configuration Summary
// ============================================================================

/**
 * @brief Print configuration summary
 *
 * Useful for debugging and verifying build configuration
 */
inline void print_config() {
    #define PRINT_BOOL(x) ((x) ? "ON" : "OFF")

    printf("FE Library Configuration:\n");
    printf("  Debug Mode: %s\n", PRINT_BOOL(FE_DEBUG_MODE));
    printf("  MPI Support: %s\n", PRINT_BOOL(FE_HAS_MPI));
    printf("  OpenMP Support: %s\n", PRINT_BOOL(FE_HAS_OPENMP));
    printf("  Max Dimension: %d\n", MAX_SPATIAL_DIM);
    printf("  Vectorization: %s (Level: %d)\n",
           PRINT_BOOL(FE_ENABLE_VECTORIZATION), FE_VECTORIZATION_LEVEL);
    printf("  Cache Line Size: %zu bytes\n", CACHE_LINE_SIZE);
    printf("  Profiling: %s\n", PRINT_BOOL(FE_ENABLE_PROFILING));
    printf("  Matrix-Free: %s\n", PRINT_BOOL(FE_ENABLE_MATRIX_FREE));
    printf("  Max Polynomial Order: %d\n", MAX_POLYNOMIAL_ORDER);

    #undef PRINT_BOOL
}

} // namespace config
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_CONFIG_H
