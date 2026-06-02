/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#ifndef SVMP_FE_QUADRATURE_RULE_H
#define SVMP_FE_QUADRATURE_RULE_H

/**
 * @file QuadratureRule.h
 * @brief Abstracted quadrature rule representation for FE integration
 *
 * This header defines the base class for all quadrature rules used by the
 * finite element infrastructure. Rules are expressed in reference element
 * space only; mapping to physical space is handled by the Geometry module.
 *
 * The interface is intentionally lightweight and header-only to avoid coupling
 * Quadrature to other modules while remaining compatible with the Mesh library
 * through shared type aliases provided by FE/Common/Types.h.
 */

#include "Types.h"
#include "FEException.h"
#include "Math/Vector.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace svmp {
namespace FE {
namespace quadrature {

/// Convenience alias for quadrature point representation in reference space
using QuadPoint = math::Vector<Real, 3>;

struct QuadraturePointFingerprint {
    int dimension{0};
    std::size_t num_points{0};
    std::uint64_t points_hash_a{0};
    std::uint64_t points_hash_b{0};
};

/**
 * @brief Base class for quadrature rules over reference elements
 *
 * Derived classes populate the point/weight data via the protected setters.
 * The class performs lightweight consistency checks (size agreement, basic
 * reference-measure validation) but leaves element-specific checks to callers.
 */
class QuadratureRule {
public:
    virtual ~QuadratureRule() = default;

    /// Number of quadrature points
    std::size_t num_points() const noexcept { return points_.size(); }

    /// Polynomial exactness degree reported by the rule
    int order() const noexcept { return order_; }

    /// Spatial dimension of the reference domain
    int dimension() const noexcept { return dimension_; }

    /// Cell family the rule integrates over (line, tri, quad, ...)
    svmp::CellFamily cell_family() const noexcept { return cell_family_; }

    /// Access a single quadrature point (no bounds checking)
    QuadPoint point(std::size_t i) const noexcept { return points_[i]; }

    /// Access a single quadrature weight (no bounds checking)
    Real weight(std::size_t i) const noexcept { return weights_[i]; }

    /// Bulk accessors
    const std::vector<QuadPoint>& points() const noexcept { return points_; }
    const std::vector<Real>& weights() const noexcept { return weights_; }

    /// Cached coordinate-only fingerprint for consumers whose values depend on
    /// reference points but not quadrature weights.
    QuadraturePointFingerprint point_fingerprint() const noexcept { return point_fingerprint_; }

    /// Stable semantic identity used by BasisCache
    virtual std::string cache_identity() const;

    /**
     * @brief Validate rule data for basic consistency
     * @param tol Relative tolerance for weight sum check
     * @return True if rule passes size and weight checks
     */
    virtual bool is_valid(Real tol = 1e-12) const;

    /**
     * @brief Reference-domain measure for the element family
     *
     * Length/area/volume of the canonical reference element:
     * - Line   [-1,1]            -> 2
     * - Quad   [-1,1]^2          -> 4
     * - Hex    [-1,1]^3          -> 8
     * - Tri    (0,0)-(1,0)-(0,1) -> 0.5
     * - Tet    simplex at origin -> 1/6
     * - Wedge  (triangle x line) -> 1
     * - Pyramid (x,y in [-1,1], z in [0,1]) -> 4/3
     */
    Real reference_measure() const noexcept;

protected:
    QuadratureRule(svmp::CellFamily family, int dimension, int order = 0)
        : cell_family_(family), dimension_(dimension), order_(order) {}

    /// Assign point and weight storage (sizes must match)
    void set_data(std::vector<QuadPoint> pts, std::vector<Real> wts);

    /// Override computed order in derived classes
    void set_order(int ord) noexcept { order_ = ord; }

private:
    std::string build_cache_identity() const;
    QuadraturePointFingerprint build_point_fingerprint() const noexcept;

    svmp::CellFamily cell_family_;
    int dimension_;
    int order_;
    std::vector<QuadPoint> points_;
    std::vector<Real> weights_;
    std::string cache_identity_;
    QuadraturePointFingerprint point_fingerprint_;
};

// --------------------------------------------------------------------------------
// Inline implementations
// --------------------------------------------------------------------------------

inline void QuadratureRule::set_data(std::vector<QuadPoint> pts, std::vector<Real> wts) {
    if (pts.size() != wts.size()) {
        throw FEException("QuadratureRule: points/weights size mismatch",
                          StatusCode::InvalidArgument,
                          __FILE__, __LINE__, __func__);
    }
    points_ = std::move(pts);
    weights_ = std::move(wts);
    point_fingerprint_ = build_point_fingerprint();
    cache_identity_ = build_cache_identity();
}

inline bool QuadratureRule::is_valid(Real tol) const {
    if (points_.empty() || points_.size() != weights_.size()) {
        return false;
    }
    Real sum_w = Real(0);
    for (Real w : weights_) {
        if (!std::isfinite(w)) {
            return false;
        }
        sum_w += w;
    }
    const Real ref = reference_measure();
    const Real denom = std::max(Real(1), std::abs(ref));
    return std::abs(sum_w - ref) <= tol * denom;
}

inline std::string QuadratureRule::cache_identity() const {
    if (!cache_identity_.empty()) {
        return cache_identity_;
    }
    return build_cache_identity();
}

inline std::string QuadratureRule::build_cache_identity() const {
    std::ostringstream oss;
    oss << "dim=" << dimension_
        << "|npts=" << points_.size();

    oss << std::setprecision(std::numeric_limits<Real>::max_digits10);
    for (const auto& pt : points_) {
        oss << "|pt=" << pt[0] << ',' << pt[1] << ',' << pt[2];
    }
    return oss.str();
}

inline QuadraturePointFingerprint QuadratureRule::build_point_fingerprint() const noexcept {
    auto real_bits = [](Real value) noexcept {
        static_assert(sizeof(Real) <= sizeof(std::uint64_t),
                      "Quadrature point fingerprints assume Real fits in 64 bits");
        std::uint64_t bits = 0;
        std::memcpy(&bits, &value, sizeof(Real));
        return bits;
    };
    auto mix_hash = [](std::uint64_t& seed, std::uint64_t value) noexcept {
        seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6u) + (seed >> 2u);
    };

    QuadraturePointFingerprint fingerprint;
    fingerprint.dimension = dimension_;
    fingerprint.num_points = points_.size();
    fingerprint.points_hash_a = 1469598103934665603ULL;
    fingerprint.points_hash_b = 1099511628211ULL;

    mix_hash(fingerprint.points_hash_a, static_cast<std::uint64_t>(fingerprint.dimension));
    mix_hash(fingerprint.points_hash_a, static_cast<std::uint64_t>(fingerprint.num_points));
    mix_hash(fingerprint.points_hash_b, static_cast<std::uint64_t>(fingerprint.num_points));
    mix_hash(fingerprint.points_hash_b, static_cast<std::uint64_t>(fingerprint.dimension));
    for (const auto& point : points_) {
        for (std::size_t component = 0; component < 3u; ++component) {
            const std::uint64_t bits = real_bits(point[component]);
            mix_hash(fingerprint.points_hash_a, bits);
            mix_hash(fingerprint.points_hash_b, bits ^ (0xbf58476d1ce4e5b9ULL + component));
        }
    }
    return fingerprint;
}

inline Real QuadratureRule::reference_measure() const noexcept {
    switch (cell_family_) {
        case svmp::CellFamily::Line:      return Real(2);
        case svmp::CellFamily::Quad:      return Real(4);
        case svmp::CellFamily::Hex:       return Real(8);
        case svmp::CellFamily::Triangle:  return Real(0.5);
        case svmp::CellFamily::Tetra:     return Real(1.0 / 6.0);
        case svmp::CellFamily::Wedge:     return Real(1.0);     // 0.5 area * length 2
        case svmp::CellFamily::Pyramid:   return Real(4.0 / 3.0);
        case svmp::CellFamily::Point:     return Real(1.0);
        default:                          return Real(1.0);
    }
}

} // namespace quadrature
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_QUADRATURE_RULE_H
