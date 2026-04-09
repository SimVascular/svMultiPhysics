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

namespace svmp {
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
} // namespace svmp

#include <cstdint>
#include <array>
#include <string>
#include <type_traits>
#include <limits>

namespace svmp {
namespace FE {

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
 */
using GlobalIndex = std::int64_t;

/**
 * @brief DOF-specific index type
 *
 * Strong type alias to prevent mixing DOF indices with other indices.
 * Provides type safety at compile time.
 */
struct DofIndex {
    GlobalIndex value;

    constexpr explicit DofIndex(GlobalIndex v = -1) noexcept : value(v) {}
    constexpr operator GlobalIndex() const noexcept { return value; }
    constexpr bool is_valid() const noexcept { return value >= 0; }
};

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

// Import mesh library scalar/index types when available.
using MeshIndex = std::int32_t;
using MeshOffset = std::int64_t;
using MeshGlobalId = std::int64_t;
using Real = double;

// ============================================================================
// Constants
// ============================================================================

constexpr LocalIndex INVALID_LOCAL_INDEX = std::numeric_limits<LocalIndex>::max();
constexpr GlobalIndex INVALID_GLOBAL_INDEX = -1;
constexpr FieldId INVALID_FIELD_ID = std::numeric_limits<FieldId>::max();
/// Sentinel FieldId for geometry-only quantities (no DOF dependence).
/// Uses first registered field's space for quadrature, but logically decoupled
/// from any specific field's DOFs.
constexpr FieldId GEOMETRY_FIELD_ID = std::numeric_limits<FieldId>::max() - 1;
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

// ============================================================================
// Field Value Entry (for point evaluation of field-dependent expressions)
// ============================================================================

/// Maximum number of components in a FieldValueEntry (3x3 tensor).
constexpr int MAX_FIELD_VALUE_COMPONENTS = 9;

/**
 * @brief Field value at an evaluation point — scalar, vector, or tensor.
 *
 * Used by PointEvaluator and the auxiliary assembly path to supply FE
 * field values at entity locations (e.g., nodal DOF values for
 * Node-scoped auxiliary models with Lagrange Kronecker delta).
 */
struct FieldValueEntry {
    FieldId field{INVALID_FIELD_ID};
    int n_components{0};
    Real components[MAX_FIELD_VALUE_COMPONENTS]{};
};

// ============================================================================
// Element Type Enumerations
// ============================================================================

/**
 * @brief Reference element types supported by the FE library
 *
 * Maps to svmp::CellFamily from the Mesh library but provides
 * FE-specific categorization including higher-order variants.
 */
enum class ElementType : std::uint8_t {
    // Linear elements
    Line2      = 0,   // 2-node line
    Triangle3  = 1,   // 3-node triangle
    Quad4      = 2,   // 4-node quadrilateral
    Tetra4     = 3,   // 4-node tetrahedron
    Hex8       = 4,   // 8-node hexahedron
    Wedge6     = 5,   // 6-node wedge/prism
    Pyramid5   = 6,   // 5-node pyramid

    // Quadratic elements
    Line3      = 10,  // 3-node line
    Triangle6  = 11,  // 6-node triangle
    Quad9      = 12,  // 9-node quadrilateral (bi-quadratic)
    Quad8      = 13,  // 8-node quadrilateral (serendipity)
    Tetra10    = 14,  // 10-node tetrahedron
    Hex27      = 15,  // 27-node hexahedron (tri-quadratic)
    Hex20      = 16,  // 20-node hexahedron (serendipity)
    Wedge15    = 17,  // 15-node wedge
    Wedge18    = 18,  // 18-node wedge (complete quadratic)
    Pyramid13  = 19,  // 13-node pyramid
    Pyramid14  = 20,  // 14-node pyramid

    // Special elements
    Point1     = 30,  // 1-node point element

    Unknown    = 255
};

/**
 * @brief Quadrature rule types
 */
enum class QuadratureType : std::uint8_t {
    GaussLegendre,     // Standard Gaussian quadrature
    GaussLobatto,      // Includes endpoints (for spectral elements)
    Newton,            // Newton-Cotes rules
    Reduced,           // Order-based reduced integration for locking
    PositionBased,     // Position-based reduced integration (legacy compatible)
    Composite,         // Composite rules for adaptivity
    Custom             // User-defined quadrature points
};

/**
 * @brief Basis function families
 */
enum class BasisType : std::uint8_t {
    Lagrange,          // Standard nodal Lagrange basis
    Hierarchical,      // Hierarchical/modal basis
    Bernstein,         // Bernstein polynomials
    NURBS,             // Non-uniform rational B-splines
    BSpline,           // Non-rational B-spline basis
    Spectral,          // Spectral element basis
    Serendipity,       // Serendipity elements
    Hermite,           // Hermite C1 continuity basis
    RaviartThomas,     // H(div) Raviart-Thomas family
    Nedelec,           // H(curl) Nedelec edge elements
    BDM,               // H(div) Brezzi-Douglas-Marini family
    Bubble,            // Interior bubble functions for enrichment
    Custom             // User-defined basis
};

/**
 * @brief Field types for function spaces
 */
enum class FieldType : std::uint8_t {
    Scalar,            // Scalar field
    Vector,            // Vector field
    Tensor,            // Tensor field
    SymmetricTensor,   // Symmetric tensor field
    Mixed              // Mixed/composite field
};

/**
 * @brief Continuity requirements for function spaces
 */
enum class Continuity : std::uint8_t {
    C0,                // Continuous (standard FEM)
    C1,                // C1 continuous (for plates/shells)
    L2,                // L2 (discontinuous)
    H_div,             // H(div) conforming
    H_curl,            // H(curl) conforming
    Custom
};

/**
 * @brief Assembly strategies
 */
enum class AssemblyStrategy : std::uint8_t {
    ElementByElement,  // Traditional element loop
    Vectorized,        // SIMD vectorized assembly
    MatrixFree,        // Matrix-free operators
    Hybrid             // Mixed strategy
};

/**
 * @brief Status codes for FE operations
 */
enum class FEStatus : std::uint8_t {
    Success           = 0,
    InvalidArgument   = 1,
    InvalidElement    = 2,
    SingularMapping   = 3,
    QuadratureError   = 4,
    AssemblyError     = 5,
    BackendError      = 6,
    NotImplemented    = 7,
    ConvergenceError  = 8,
    AllocationError   = 9,
    MPIError          = 10,
    IOError           = 11,
    Unknown           = 255
};

// ============================================================================
// Geometric Types
// ============================================================================

/**
 * @brief Point in reference element coordinates
 */
template<int Dim>
using ReferencePoint = std::array<Real, static_cast<std::size_t>(Dim)>;

/**
 * @brief Point in physical coordinates
 */
using PhysicalPoint = std::array<Real, 3>;

/**
 * @brief Jacobian matrix type
 */
template<int SpatialDim, int ReferenceDim = SpatialDim>
using Jacobian = std::array<std::array<Real, static_cast<std::size_t>(ReferenceDim)>, static_cast<std::size_t>(SpatialDim)>;

// ============================================================================
// Strong Type Wrappers
// ============================================================================

/**
 * @brief Strong type wrapper template for type-safe programming
 *
 * Prevents accidental mixing of conceptually different types that have
 * the same underlying representation.
 */
template<typename T, typename Tag>
class StrongType {
public:
    using ValueType = T;

    constexpr StrongType() noexcept(std::is_nothrow_default_constructible_v<T>)
        : value_{} {}

    constexpr explicit StrongType(T value) noexcept(std::is_nothrow_move_constructible_v<T>)
        : value_(std::move(value)) {}

    constexpr T& get() noexcept { return value_; }
    constexpr const T& get() const noexcept { return value_; }

    // Explicit conversion
    constexpr explicit operator T() const noexcept { return value_; }

    // Comparison operators
    constexpr bool operator==(const StrongType& other) const noexcept {
        return value_ == other.value_;
    }
    constexpr bool operator!=(const StrongType& other) const noexcept {
        return value_ != other.value_;
    }
    constexpr bool operator<(const StrongType& other) const noexcept {
        return value_ < other.value_;
    }

private:
    T value_;
};

// Specific strong types for common use cases
struct QuadraturePointTag {};
struct QuadratureWeightTag {};
struct BasisValueTag {};
struct BasisGradientTag {};

using QuadraturePointIndex = StrongType<LocalIndex, QuadraturePointTag>;
using QuadratureWeight = StrongType<Real, QuadratureWeightTag>;

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

template<typename T>
inline constexpr bool is_index_type_v = is_index_type<T>::value;

/**
 * @brief Check if a type represents a field type
 */
template<typename T>
struct is_field_type : std::false_type {};

template<>
struct is_field_type<FieldType> : std::true_type {};

template<typename T>
inline constexpr bool is_field_type_v = is_field_type<T>::value;

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * @brief Convert FE ElementType to Mesh CellFamily
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

/**
 * @brief Convert status code to string for error reporting
 */
inline const char* status_to_string(FEStatus status) noexcept {
    switch(status) {
        case FEStatus::Success:          return "Success";
        case FEStatus::InvalidArgument:  return "Invalid argument";
        case FEStatus::InvalidElement:   return "Invalid element";
        case FEStatus::SingularMapping:  return "Singular mapping";
        case FEStatus::QuadratureError:  return "Quadrature error";
        case FEStatus::AssemblyError:    return "Assembly error";
        case FEStatus::BackendError:     return "Backend error";
        case FEStatus::NotImplemented:   return "Not implemented";
        case FEStatus::ConvergenceError: return "Convergence error";
        case FEStatus::AllocationError:  return "Allocation error";
        case FEStatus::MPIError:         return "MPI error";
        case FEStatus::IOError:          return "I/O error";
        default:                         return "Unknown error";
    }
}

} // namespace FE
} // namespace svmp

#endif // SVMP_FE_TYPES_H
