// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_TYPES_H
#define SVMP_FE_TYPES_H

/**
 * @file Types.h
 * @brief Fundamental type definitions for the finite element library
 *
 * This header provides core type aliases, enumerations, and strong type
 * definitions used throughout the FE library. It establishes a consistent
 * type system that integrates with the Mesh library while maintaining
 * independence from backend-specific types.
 */

// The Mesh library is an optional, external module. When the build enables it
// (SVMP_FE_WITH_MESH), FE imports the Mesh scalar/index types so the two libraries
// share a vocabulary; otherwise FE compiles standalone using the fallback
// definitions below (e.g. svmp::CellFamily and the Mesh* aliases). The Mesh
// headers are not part of this repository.
#if defined(SVMP_FE_WITH_MESH) && SVMP_FE_WITH_MESH
#  include "Mesh/Core/MeshTypes.h"
/** Nonzero when FE shares scalar/index types with the Mesh library. */
#  define SVMP_FE_HAS_MESH_TYPES 1
#else
// Build FE without Mesh types unless explicitly enabled.
/** Nonzero when FE shares scalar/index types with the Mesh library. */
#  define SVMP_FE_HAS_MESH_TYPES 0
#endif

#if !SVMP_FE_HAS_MESH_TYPES
namespace svmp {
#ifndef SVMP_CELL_FAMILY_DEFINED
/** Guard marking that svmp::CellFamily has been defined. */
#define SVMP_CELL_FAMILY_DEFINED 1
/**
 * @brief Minimal fallback for svmp::CellFamily when the Mesh library is unavailable
 * @ingroup FE_CommonTypes
 *
 * Keeps FE compilation self-contained while preserving the same namespace
 * and enumerator set as the Mesh library's cell-family classification.
 */
enum class CellFamily {
    Point,
    Line,
    Triangle,
    Quad,
    Tetra,
    Hex,
    Wedge,
    Pyramid,
    Polygon,
    Polyhedron
};
#endif
} // namespace svmp
#endif
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <type_traits>
#include <limits>

/**
 * @defgroup FE_Common Common
 * @ingroup FE
 * @brief Shared vocabulary types, constants, and exception infrastructure used by every FE module.
 *
 * @details The Common module collects the foundational definitions that the
 * rest of the FE library builds on: index and scalar type aliases; element,
 * basis, quadrature, and field enumerations; sentinel constants and strong
 * type wrappers; and the FE exception hierarchy together with its
 * argument-checking helpers.
 */

namespace svmp::FE {

/**
 * @defgroup FE_CommonTypes Types
 * @ingroup FE_Common
 * @brief Core type aliases, enumerations, constants, geometric types, and compile-time traits.
 *
 * @details This group documents the index and identifier types used for
 * element-local and global numbering, the element/basis/quadrature/field
 * enumerations shared across modules, sentinel constants, reference- and
 * physical-space geometric aliases, and the strong-type utilities that
 * prevent accidental mixing of conceptually distinct values.
 * @{
 */

// ============================================================================
// Index Types
// ============================================================================

/**
 * @brief Local index type for element-level operations
 *
 * Used for local node numbering within elements, local DOF indices,
 * and other element-local indexing. Unsigned for safety.
 */
using LocalIndex = std::uint32_t;

/**
 * @brief Global index type for distributed DOF numbering
 *
 * Signed 64-bit for compatibility with PETSc and Trilinos.
 * Negative values can indicate special conditions or invalid indices.
 *
 * @note Kept as a plain integer alias rather than a StrongType wrapper: this is
 * the raw interop type handed directly to PETSc/Trilinos, where a wrapper would
 * force an unwrap at every call. Type safety for DOF indices is provided by
 * DofIndex (below), the strong wrapper around a GlobalIndex.
 */
using GlobalIndex = std::int64_t;

/**
 * @brief Field identifier type
 *
 * Used to distinguish between different physical fields in multi-field problems.
 */
using FieldId = std::uint16_t;

/**
 * @brief Block identifier for block-structured systems
 */
using BlockId = std::uint16_t;

// Import mesh library scalar/index types when available (optional dependency).
#if SVMP_FE_HAS_MESH_TYPES
using MeshIndex = svmp::index_t;        ///< Local mesh entity index, shared with the Mesh library.
using MeshOffset = svmp::offset_t;      ///< Offset type for mesh connectivity arrays.
using MeshGlobalId = svmp::gid_t;       ///< Global mesh entity identifier.
#else
using MeshIndex = std::int32_t;         ///< Local mesh entity index, shared with the Mesh library.
using MeshOffset = std::int64_t;        ///< Offset type for mesh connectivity arrays.
using MeshGlobalId = std::int64_t;      ///< Global mesh entity identifier.
#endif

// ============================================================================
// Constants
// ============================================================================

/** Sentinel for an unset or out-of-range local index. */
constexpr LocalIndex INVALID_LOCAL_INDEX = std::numeric_limits<LocalIndex>::max();
/** Sentinel for an unset or out-of-range global index. */
constexpr GlobalIndex INVALID_GLOBAL_INDEX = -1;
/** Sentinel FieldId meaning "uninitialized / no field". */
constexpr FieldId INVALID_FIELD_ID = std::numeric_limits<FieldId>::max();
/**
 * Sentinel FieldId for geometry-only quantities (no DOF dependence).
 * Uses first registered field's space for quadrature, but logically decoupled
 * from any specific field's DOFs.
 */
constexpr FieldId GEOMETRY_FIELD_ID = std::numeric_limits<FieldId>::max() - 1;
/** Sentinel for an unset or out-of-range block identifier. */
constexpr BlockId INVALID_BLOCK_ID = std::numeric_limits<BlockId>::max();

/**
 * @brief Sentinel FieldId representing "the current solution state" in tangent forms.
 *
 * When differentiating a residual form to obtain the tangent (Jacobian), undifferentiated
 * TrialFunction occurrences are rewritten to StateField nodes. Those that represent the
 * block's own primary unknown (rather than a named external field) use this sentinel
 * FieldId. The assembler maps it to the current solution coefficients at each quadrature
 * point, regardless of which physics or field variables are involved.
 *
 * This is distinct from INVALID_FIELD_ID, which means "uninitialized / no field."
 * CURRENT_SOLUTION_FIELD_ID uses the same numeric value for backward compatibility
 * with existing KernelIR encodings, but carries explicit semantic intent.
 */
constexpr FieldId CURRENT_SOLUTION_FIELD_ID = std::numeric_limits<FieldId>::max();

/** Preferred cache-line/SIMD alignment for performance-critical arrays. */
inline constexpr std::size_t kFEPreferredAlignmentBytes = 64u;

/** Alignment for small fixed-size math objects that are commonly passed by value. */
inline constexpr std::size_t kFEFixedObjectAlignmentBytes = 32u;

// ============================================================================
// Field Value Entry (for point evaluation of field-dependent expressions)
// ============================================================================

/** Maximum number of components in a FieldValueEntry (3x3 tensor). */
constexpr int MAX_FIELD_VALUE_COMPONENTS = 9;

/**
 * @brief Field value at an evaluation point — scalar, vector, or tensor.
 *
 * Used by PointEvaluator and the auxiliary assembly path to supply FE
 * field values at entity locations (e.g., nodal DOF values for
 * Node-scoped auxiliary models with Lagrange Kronecker delta).
 */
struct FieldValueEntry {
    FieldId field{INVALID_FIELD_ID};                  ///< Field this value belongs to.
    int n_components{0};                              ///< Number of valid entries in components.
    double components[MAX_FIELD_VALUE_COMPONENTS]{};    ///< Component values, row-major for tensors.
};

// ============================================================================
// Element Type Enumerations
// ============================================================================

/**
 * @brief Reference element types supported by the FE library
 *
 * Maps to svmp::CellFamily from the Mesh library but provides
 * FE-specific categorization including higher-order variants.
 *
 * @note The enum is consumed by name (the switches in to_mesh_family() and
 * element_dimension() and the basis classifiers); nothing depends on the
 * underlying numeric values, so they are left implicit. Entries are grouped by
 * polynomial order (linear, quadratic) plus a special section.
 */
enum class ElementType : std::uint8_t {
    // Linear elements
    Line2,       ///< 2-node line
    Triangle3,   ///< 3-node triangle
    Quad4,       ///< 4-node quadrilateral
    Tetra4,      ///< 4-node tetrahedron
    Hex8,        ///< 8-node hexahedron
    Wedge6,      ///< 6-node wedge/prism
    Pyramid5,    ///< 5-node pyramid

    // Quadratic elements
    Line3,       ///< 3-node line
    Triangle6,   ///< 6-node triangle
    Quad9,       ///< 9-node quadrilateral (bi-quadratic)
    Quad8,       ///< 8-node quadrilateral (serendipity)
    Tetra10,     ///< 10-node tetrahedron
    Hex27,       ///< 27-node hexahedron (tri-quadratic)
    Hex20,       ///< 20-node hexahedron (serendipity)
    Wedge15,     ///< 15-node wedge
    Wedge18,     ///< 18-node wedge (complete quadratic)
    Pyramid13,   ///< 13-node pyramid
    Pyramid14,   ///< 14-node pyramid

    // Special elements
    Point1,      ///< 1-node point element

    Unknown      ///< Unrecognized or uninitialized element type
};

/**
 * @brief Quadrature rule types
 */
enum class QuadratureType : std::uint8_t {
    GaussLegendre,     ///< Standard Gaussian quadrature
    GaussLobatto,      ///< Includes endpoints (for spectral elements)
    Newton,            ///< Newton-Cotes rules
    Reduced,           ///< Order-based reduced integration for locking
    PositionBased,     ///< Position-based reduced integration (legacy compatible)
    Composite,         ///< Composite rules for adaptivity
    Custom             ///< User-defined quadrature points
};

/**
 * @brief Basis function families
 */
enum class BasisType : std::uint8_t {
    Lagrange,          ///< Standard nodal Lagrange basis
    NURBS,             ///< Non-uniform rational B-splines (reserved; not yet implemented)
    Serendipity,       ///< Serendipity elements
    Custom             ///< User-defined basis
};

/**
 * @brief Field types for function spaces
 */
enum class FieldType : std::uint8_t {
    Scalar,            ///< Scalar field (temperature, pressure)
    Vector,            ///< Vector field (velocity, displacement)
    Tensor,            ///< Tensor field (stress, strain)
    SymmetricTensor,   ///< Symmetric tensor field
    Mixed              ///< Mixed/composite field
};

/**
 * @brief Continuity requirements for function spaces
 */
enum class Continuity : std::uint8_t {
    C0,                ///< Continuous (standard FEM)
    C1,                ///< C1 continuous (for plates/shells)
    L2,                ///< L2 (discontinuous)
    H_div,             ///< H(div) conforming
    H_curl,            ///< H(curl) conforming
    Custom             ///< User-defined continuity requirement
};

/**
 * @brief Assembly strategies
 */
enum class AssemblyStrategy : std::uint8_t {
    ElementByElement,  ///< Traditional element loop
    Vectorized,        ///< SIMD vectorized assembly
    MatrixFree,        ///< Matrix-free operators
    Hybrid             ///< Mixed strategy
};

// ============================================================================
// Geometric Types
// ============================================================================

/**
 * @brief Point in reference element coordinates
 * @tparam Dim Reference-space dimension
 */
template<int Dim>
using ReferencePoint = std::array<double, static_cast<std::size_t>(Dim)>;

/**
 * @brief Point in physical coordinates
 */
using PhysicalPoint = std::array<double, 3>;

/**
 * @brief Jacobian matrix type
 * @tparam SpatialDim Physical-space dimension (rows)
 * @tparam ReferenceDim Reference-space dimension (columns)
 */
template<int SpatialDim, int ReferenceDim = SpatialDim>
using Jacobian = std::array<std::array<double, static_cast<std::size_t>(ReferenceDim)>, static_cast<std::size_t>(SpatialDim)>;

// ============================================================================
// Strong Type Wrappers (C++17 idiom for type safety)
// ============================================================================

/**
 * @brief Strong type wrapper template for type-safe programming
 *
 * Prevents accidental mixing of conceptually different types that have
 * the same underlying representation.
 *
 * @tparam T Underlying value type
 * @tparam Tag Empty tag type that distinguishes otherwise identical wrappers
 */
template<typename T, typename Tag>
class StrongType {
public:
    /** @brief Underlying value type. */
    using ValueType = T;

    /** @brief Value-initialize the wrapped value. */
    constexpr StrongType() noexcept(std::is_nothrow_default_constructible_v<T>)
        : value_{} {}

    /**
     * @brief Wrap an explicit value.
     * @param value Value to store.
     */
    constexpr explicit StrongType(T value) noexcept(std::is_nothrow_move_constructible_v<T>)
        : value_(std::move(value)) {}

    /**
     * @brief Access the wrapped value.
     * @return Reference to the wrapped value.
     */
    constexpr T& get() noexcept { return value_; }
    /**
     * @brief Access the wrapped value.
     * @return Reference to the wrapped value.
     */
    constexpr const T& get() const noexcept { return value_; }

    /**
     * @brief Explicitly convert back to the underlying type.
     * @return Copy of the wrapped value.
     */
    constexpr explicit operator T() const noexcept { return value_; }

    /**
     * @brief Compare wrapped values for equality.
     * @param other Wrapper to compare against.
     * @return True when the wrapped values are equal.
     */
    constexpr bool operator==(const StrongType& other) const noexcept {
        return value_ == other.value_;
    }
    /**
     * @brief Compare wrapped values for inequality.
     * @param other Wrapper to compare against.
     * @return True when the wrapped values differ.
     */
    constexpr bool operator!=(const StrongType& other) const noexcept {
        return value_ != other.value_;
    }
    /**
     * @brief Order by wrapped value.
     * @param other Wrapper to compare against.
     * @return True when this wrapped value orders before the other.
     */
    constexpr bool operator<(const StrongType& other) const noexcept {
        return value_ < other.value_;
    }

private:
    T value_;
};

// Specific strong types for common use cases
struct QuadraturePointTag {};   ///< Tag type for quadrature-point indices.
struct QuadratureWeightTag {};  ///< Tag type for quadrature weights.
struct BasisValueTag {};        ///< Tag type for basis-function values.
struct BasisGradientTag {};     ///< Tag type for basis-function gradients.
struct DofTag {};               ///< Tag type for global DOF indices.

/** Type-safe index of a quadrature point within a rule. */
using QuadraturePointIndex = StrongType<LocalIndex, QuadraturePointTag>;
/** Type-safe quadrature weight value. */
using QuadratureWeight = StrongType<double, QuadratureWeightTag>;

/**
 * @brief DOF-specific index type
 *
 * @details A StrongType over GlobalIndex that prevents mixing DOF indices with
 * other indices: conversion back to GlobalIndex is explicit (via get()), so a
 * DofIndex cannot silently decay to a raw integer. Over the base StrongType it
 * adds two DOF-specific conveniences -- it default-constructs to the invalid
 * sentinel (-1) and exposes is_valid() -- for distributed DOF numbering where a
 * negative value marks an unset or non-local DOF.
 */
class DofIndex : public StrongType<GlobalIndex, DofTag> {
public:
    using StrongType::StrongType;

    /** @brief Construct an invalid DOF index (the negative sentinel). */
    constexpr DofIndex() noexcept : StrongType(GlobalIndex{-1}) {}

    /**
     * @brief Check whether this index refers to a valid DOF.
     * @return True when the stored value is non-negative.
     */
    constexpr bool is_valid() const noexcept { return get() >= 0; }
};

// ============================================================================
// Type Traits
// ============================================================================

/**
 * @brief Check if a type is a valid index type
 */
template<typename T>
struct is_index_type : std::false_type {};

template<>
struct is_index_type<LocalIndex> : std::true_type {};

template<>
struct is_index_type<GlobalIndex> : std::true_type {};

template<>
struct is_index_type<DofIndex> : std::true_type {};

/** Convenience variable template for is_index_type. */
template<typename T>
inline constexpr bool is_index_type_v = is_index_type<T>::value;

/**
 * @brief Check if a type represents a field type
 */
template<typename T>
struct is_field_type : std::false_type {};

template<>
struct is_field_type<FieldType> : std::true_type {};

/** Convenience variable template for is_field_type. */
template<typename T>
inline constexpr bool is_field_type_v = is_field_type<T>::value;

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * @brief Convert FE ElementType to Mesh CellFamily
 * @param elem Element type to classify.
 * @return Cell family of the element's linear topology; Point for unknown types.
 */
constexpr svmp::CellFamily to_mesh_family(ElementType elem) noexcept {
    switch(elem) {
        case ElementType::Line2:
        case ElementType::Line3:
            return svmp::CellFamily::Line;

        case ElementType::Triangle3:
        case ElementType::Triangle6:
            return svmp::CellFamily::Triangle;

        case ElementType::Quad4:
        case ElementType::Quad8:
        case ElementType::Quad9:
            return svmp::CellFamily::Quad;

        case ElementType::Tetra4:
        case ElementType::Tetra10:
            return svmp::CellFamily::Tetra;

        case ElementType::Hex8:
        case ElementType::Hex20:
        case ElementType::Hex27:
            return svmp::CellFamily::Hex;

        case ElementType::Wedge6:
        case ElementType::Wedge15:
        case ElementType::Wedge18:
            return svmp::CellFamily::Wedge;

        case ElementType::Pyramid5:
        case ElementType::Pyramid13:
        case ElementType::Pyramid14:
            return svmp::CellFamily::Pyramid;

        case ElementType::Point1:
            return svmp::CellFamily::Point;

        default:
            return svmp::CellFamily::Point;  // Fallback
    }
}

/**
 * @brief Get spatial dimension of element type
 * @param elem Element type to query.
 * @return Reference dimension from 0 (point) to 3 (volume); -1 for unknown types.
 */
constexpr int element_dimension(ElementType elem) noexcept {
    switch(elem) {
        case ElementType::Point1:
            return 0;
        case ElementType::Line2:
        case ElementType::Line3:
            return 1;
        case ElementType::Triangle3:
        case ElementType::Triangle6:
        case ElementType::Quad4:
        case ElementType::Quad8:
        case ElementType::Quad9:
            return 2;
        case ElementType::Tetra4:
        case ElementType::Tetra10:
        case ElementType::Hex8:
        case ElementType::Hex20:
        case ElementType::Hex27:
        case ElementType::Wedge6:
        case ElementType::Wedge15:
        case ElementType::Wedge18:
        case ElementType::Pyramid5:
        case ElementType::Pyramid13:
        case ElementType::Pyramid14:
            return 3;
        default:
            return -1;
    }
}

/** @} */

} // namespace svmp::FE

#endif // SVMP_FE_TYPES_H
