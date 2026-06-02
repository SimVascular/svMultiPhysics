/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#ifndef SVMP_FE_BASIS_BASISCACHE_H
#define SVMP_FE_BASIS_BASISCACHE_H

/**
 * @file BasisCache.h
 * @brief Cache for basis evaluations at quadrature points
 */

#include "BasisFunction.h"
#include "Quadrature/QuadratureRule.h"
#include <cstddef>
#include <condition_variable>
#include <exception>
#include <functional>
#include <memory>
#include <mutex>
#include <cstdint>
#include <shared_mutex>
#include <span>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

namespace svmp {
namespace FE {
namespace basis {

struct QuadratureCacheKey {
    int dimension{0};
    std::size_t num_points{0};
    // Quadrature coordinates are intentionally fingerprinted from their exact
    // Real bit patterns. Values such as -0.0 and +0.0 therefore produce
    // distinct cache keys unless a future API explicitly normalizes them. The
    // key intentionally ignores weights and rule class because basis values only
    // depend on reference coordinates; bit-identical point sets share entries.
    std::uint64_t points_hash_a{0};
    std::uint64_t points_hash_b{0};

    bool operator==(const QuadratureCacheKey& other) const noexcept {
        return dimension == other.dimension &&
               num_points == other.num_points &&
               points_hash_a == other.points_hash_a &&
               points_hash_b == other.points_hash_b;
    }
};

struct StructuralBasisKey {
    BasisType basis_type{BasisType::Custom};
    ElementType element_type{ElementType::Unknown};
    int dimension{0};
    int order{0};
    std::size_t num_dofs{0};
    bool vector_valued{false};
    QuadratureCacheKey quadrature;
    bool with_gradients{false};
    bool with_hessians{false};

    bool operator==(const StructuralBasisKey& other) const noexcept {
        return basis_type == other.basis_type &&
               element_type == other.element_type &&
               dimension == other.dimension &&
               order == other.order &&
               num_dofs == other.num_dofs &&
               vector_valued == other.vector_valued &&
               quadrature == other.quadrature &&
               with_gradients == other.with_gradients &&
               with_hessians == other.with_hessians;
    }
};

struct ParameterizedBasisKey {
    StructuralBasisKey structural;
    bool uses_structured_identity{false};
    std::uint64_t identity_hash_a{0};
    std::uint64_t identity_hash_b{0};
    std::vector<std::uint64_t> basis_identity_words;
    std::string basis_identity;

    bool operator==(const ParameterizedBasisKey& other) const noexcept {
        return structural == other.structural &&
               uses_structured_identity == other.uses_structured_identity &&
               identity_hash_a == other.identity_hash_a &&
               identity_hash_b == other.identity_hash_b &&
               basis_identity_words == other.basis_identity_words &&
               basis_identity == other.basis_identity;
    }
};

struct BasisCacheKey {
    std::variant<StructuralBasisKey, ParameterizedBasisKey> value;

    bool operator==(const BasisCacheKey& other) const noexcept {
        return value == other.value;
    }
};

struct BasisCacheKeyHash {
    std::size_t operator()(const BasisCacheKey& key) const noexcept {
        std::size_t seed = 0;
        auto combine = [&seed](std::size_t value) noexcept {
            seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6u) + (seed >> 2u);
        };

        auto hash_structural = [&](const StructuralBasisKey& structural) noexcept {
            combine(std::hash<int>()(structural.quadrature.dimension));
            combine(std::hash<std::size_t>()(structural.quadrature.num_points));
            combine(std::hash<std::uint64_t>()(structural.quadrature.points_hash_a));
            combine(std::hash<std::uint64_t>()(structural.quadrature.points_hash_b));
            combine(std::hash<int>()(static_cast<int>(structural.basis_type)));
            combine(std::hash<int>()(static_cast<int>(structural.element_type)));
            combine(std::hash<int>()(structural.dimension));
            combine(std::hash<int>()(structural.order));
            combine(std::hash<std::size_t>()(structural.num_dofs));
            unsigned flags = 0u;
            flags |= structural.vector_valued ? 1u : 0u;
            flags |= structural.with_gradients ? 2u : 0u;
            flags |= structural.with_hessians ? 4u : 0u;
            combine(std::hash<unsigned>()(flags));
        };

        std::visit([&](const auto& active_key) {
            using ActiveKey = std::decay_t<decltype(active_key)>;
            if constexpr (std::is_same_v<ActiveKey, StructuralBasisKey>) {
                combine(0x5354525543544b45ULL);
                hash_structural(active_key);
            } else {
                combine(0x504152414d4b4559ULL);
                hash_structural(active_key.structural);
                combine(active_key.uses_structured_identity ? 1u : 0u);
                combine(std::hash<std::uint64_t>()(active_key.identity_hash_a));
                combine(std::hash<std::uint64_t>()(active_key.identity_hash_b));
            }
        }, key.value);
        return seed;
    }
};

struct BasisCacheEntry {
    std::size_t num_qpts{0};
    std::size_t num_dofs{0};
    // Scalar basis values in dof-major SoA layout: [dof * num_qpts + qp].
    std::vector<Real> scalar_values;
    // Scalar reference gradients in dof/component/qpt SoA layout:
    // [(dof * 3 + component) * num_qpts + qp].
    std::vector<Real> gradients;
    // Scalar reference Hessians in dof/component/qpt SoA layout:
    // [(dof * 9 + row * 3 + col) * num_qpts + qp].
    std::vector<Real> hessians;

    // Vector basis values in dof/component/qpt SoA layout:
    // [(dof * 3 + component) * num_qpts + qp].
    std::vector<Real> vector_values_xyz;
    // Vector basis reference Jacobians in dof/component/derivative/qpt layout:
    // [(dof * 9 + component * 3 + derivative) * num_qpts + qp].
    std::vector<Real> vector_jacobians;
    // Vector basis curls in dof/component/qpt SoA layout.
    std::vector<Real> vector_curls_xyz;
    // Vector basis divergences in dof/qpt SoA layout.
    std::vector<Real> vector_divergence;

    // The object-returning accessors below are convenience helpers for tests,
    // diagnostics, and occasional scalar use. Hot loops should prefer the SoA
    // span accessors so they do not reconstruct Gradient, Hessian, or matrix
    // objects per DOF and quadrature point.

    [[nodiscard]] Real scalarValue(std::size_t dof, std::size_t qp) const noexcept {
        return scalar_values[dof * num_qpts + qp];
    }

    [[nodiscard]] std::span<const Real> scalarValuesForDof(std::size_t dof) const noexcept {
        if (num_qpts == 0) return {};
        return std::span<const Real>(scalar_values.data() + dof * num_qpts, num_qpts);
    }

    [[nodiscard]] Real gradientValue(std::size_t dof,
                                     std::size_t component,
                                     std::size_t qp) const noexcept {
        return gradients[(dof * 3u + component) * num_qpts + qp];
    }

    [[nodiscard]] Gradient gradientVector(std::size_t dof, std::size_t qp) const noexcept {
        Gradient out{};
        for (std::size_t component = 0; component < 3u; ++component) {
            out[component] = gradientValue(dof, component, qp);
        }
        return out;
    }

    [[nodiscard]] std::span<const Real> gradientsForDofComponent(std::size_t dof,
                                                                  std::size_t component) const noexcept {
        if (num_qpts == 0) return {};
        return std::span<const Real>(gradients.data() + (dof * 3u + component) * num_qpts, num_qpts);
    }

    [[nodiscard]] std::span<const Real> gradientsForDof(std::size_t dof) const noexcept {
        if (num_qpts == 0) return {};
        return std::span<const Real>(gradients.data() + dof * 3u * num_qpts, 3u * num_qpts);
    }

    [[nodiscard]] Real hessianValue(std::size_t dof,
                                    std::size_t row,
                                    std::size_t col,
                                    std::size_t qp) const noexcept {
        return hessians[(dof * 9u + row * 3u + col) * num_qpts + qp];
    }

    [[nodiscard]] Hessian hessianMatrix(std::size_t dof, std::size_t qp) const noexcept {
        Hessian out{};
        for (std::size_t row = 0; row < 3u; ++row) {
            for (std::size_t col = 0; col < 3u; ++col) {
                out(row, col) = hessianValue(dof, row, col, qp);
            }
        }
        return out;
    }

    [[nodiscard]] std::span<const Real> hessiansForDofComponent(std::size_t dof,
                                                                 std::size_t row,
                                                                 std::size_t col) const noexcept {
        if (num_qpts == 0) return {};
        return std::span<const Real>(hessians.data() + (dof * 9u + row * 3u + col) * num_qpts, num_qpts);
    }

    [[nodiscard]] std::span<const Real> hessiansForDof(std::size_t dof) const noexcept {
        if (num_qpts == 0) return {};
        return std::span<const Real>(hessians.data() + dof * 9u * num_qpts, 9u * num_qpts);
    }

    [[nodiscard]] Real vectorValue(std::size_t dof,
                                   std::size_t component,
                                   std::size_t qp) const noexcept {
        return vector_values_xyz[(dof * 3u + component) * num_qpts + qp];
    }

    [[nodiscard]] math::Vector<Real, 3> vectorValue(std::size_t dof,
                                                     std::size_t qp) const noexcept {
        math::Vector<Real, 3> out{};
        for (std::size_t component = 0; component < 3u; ++component) {
            out[component] = vectorValue(dof, component, qp);
        }
        return out;
    }

    [[nodiscard]] std::span<const Real> vectorValuesForDofComponent(std::size_t dof,
                                                                     std::size_t component) const noexcept {
        if (num_qpts == 0) return {};
        return std::span<const Real>(vector_values_xyz.data() + (dof * 3u + component) * num_qpts, num_qpts);
    }

    [[nodiscard]] std::span<const Real> vectorValuesForDof(std::size_t dof) const noexcept {
        if (num_qpts == 0 || vector_values_xyz.empty()) return {};
        return std::span<const Real>(vector_values_xyz.data() + dof * 3u * num_qpts, 3u * num_qpts);
    }

    [[nodiscard]] Real vectorJacobianValue(std::size_t dof,
                                           std::size_t component,
                                           std::size_t derivative,
                                           std::size_t qp) const noexcept {
        return vector_jacobians[(dof * 9u + component * 3u + derivative) * num_qpts + qp];
    }

    [[nodiscard]] VectorJacobian vectorJacobianMatrix(std::size_t dof,
                                                       std::size_t qp) const noexcept {
        VectorJacobian out{};
        for (std::size_t component = 0; component < 3u; ++component) {
            for (std::size_t derivative = 0; derivative < 3u; ++derivative) {
                out(component, derivative) =
                    vectorJacobianValue(dof, component, derivative, qp);
            }
        }
        return out;
    }

    [[nodiscard]] std::span<const Real> vectorJacobiansForDofComponentDerivative(
        std::size_t dof,
        std::size_t component,
        std::size_t derivative) const noexcept {
        if (num_qpts == 0 || vector_jacobians.empty()) return {};
        return std::span<const Real>(
            vector_jacobians.data() + (dof * 9u + component * 3u + derivative) * num_qpts,
            num_qpts);
    }

    [[nodiscard]] std::span<const Real> vectorJacobiansForDof(std::size_t dof) const noexcept {
        if (num_qpts == 0 || vector_jacobians.empty()) return {};
        return std::span<const Real>(vector_jacobians.data() + dof * 9u * num_qpts, 9u * num_qpts);
    }

    [[nodiscard]] Real vectorCurlValue(std::size_t dof,
                                       std::size_t component,
                                       std::size_t qp) const noexcept {
        return vector_curls_xyz[(dof * 3u + component) * num_qpts + qp];
    }

    [[nodiscard]] math::Vector<Real, 3> vectorCurl(std::size_t dof,
                                                    std::size_t qp) const noexcept {
        math::Vector<Real, 3> out{};
        for (std::size_t component = 0; component < 3u; ++component) {
            out[component] = vectorCurlValue(dof, component, qp);
        }
        return out;
    }

    [[nodiscard]] std::span<const Real> vectorCurlsForDofComponent(std::size_t dof,
                                                                    std::size_t component) const noexcept {
        if (num_qpts == 0 || vector_curls_xyz.empty()) return {};
        return std::span<const Real>(vector_curls_xyz.data() + (dof * 3u + component) * num_qpts, num_qpts);
    }

    [[nodiscard]] std::span<const Real> vectorCurlsForDof(std::size_t dof) const noexcept {
        if (num_qpts == 0 || vector_curls_xyz.empty()) return {};
        return std::span<const Real>(vector_curls_xyz.data() + dof * 3u * num_qpts, 3u * num_qpts);
    }

    [[nodiscard]] Real vectorDivergenceValue(std::size_t dof,
                                             std::size_t qp) const noexcept {
        return vector_divergence[dof * num_qpts + qp];
    }

    [[nodiscard]] std::span<const Real> vectorDivergenceForDof(std::size_t dof) const noexcept {
        if (num_qpts == 0 || vector_divergence.empty()) return {};
        return std::span<const Real>(vector_divergence.data() + dof * num_qpts, num_qpts);
    }
};

class BasisCacheHandle {
public:
    BasisCacheHandle() = default;

    [[nodiscard]] const BasisCacheEntry& entry() const {
        BASIS_CHECK_CONFIG(entry_ != nullptr,
                           "BasisCacheHandle: attempted to access an empty handle");
        return *entry_;
    }

    [[nodiscard]] bool valid() const noexcept { return entry_ != nullptr; }
    explicit operator bool() const noexcept { return valid(); }

private:
    friend class BasisCache;

    explicit BasisCacheHandle(std::shared_ptr<const BasisCacheEntry> entry)
        : entry_(std::move(entry)) {}

    std::shared_ptr<const BasisCacheEntry> entry_;
};

class BasisCache {
public:
    static BasisCache& instance();

    const BasisCacheEntry& get_or_compute(
        const BasisFunction& basis,
        const quadrature::QuadratureRule& quad,
        bool gradients = true,
        bool hessians = false);

    /**
     * @brief Compute an entry without consulting, publishing to, or waiting on
     * the shared cache.
     */
    BasisCacheEntry compute_uncached(
        const BasisFunction& basis,
        const quadrature::QuadratureRule& quad,
        bool gradients = true,
        bool hessians = false) const;

    /**
     * @brief Eagerly populate the cache for the given (basis, quadrature) key
     *
     * Pays the compute cost up front so that subsequent get_or_compute calls
     * for the same key hit the warm-cache path immediately. Equivalent to
     * calling get_or_compute and discarding the return value.
     *
     * Returns the inserted (or pre-existing) entry for convenience.
     */
    const BasisCacheEntry& prewarm(
        const BasisFunction& basis,
        const quadrature::QuadratureRule& quad,
        bool gradients = true,
        bool hessians = false);

    /**
     * @brief Eagerly populate the cache and return a hot-loop handle.
     *
     * The returned handle owns a shared reference to the completed entry. Access
     * through BasisCacheHandle::entry() performs no key construction, hashing,
     * map lookup, or cache mutex acquisition. Calling clear() removes the entry
     * from the global lookup map but does not invalidate existing handles.
     */
    BasisCacheHandle prewarm_handle(
        const BasisFunction& basis,
        const quadrature::QuadratureRule& quad,
        bool gradients = true,
        bool hessians = false);

    /**
     * @brief Remove completed cache entries.
     *
     * This is a soft clear: computations that were already in flight before
     * clear() was called are allowed to publish their completed entry afterward.
     * This preserves the returned-reference lifetime contract for concurrent
     * get_or_compute() callers while still dropping all entries that had already
     * completed at the time of the call.
     */
    void clear();
    std::size_t size() const;

private:
    struct InFlightComputation {
        std::mutex mutex;
        std::condition_variable ready_cv;
        bool ready{false};
        std::shared_ptr<BasisCacheEntry> entry;
        std::exception_ptr exception;
    };

    struct CacheSlot {
        std::shared_ptr<BasisCacheEntry> entry;
        std::shared_ptr<InFlightComputation> pending;
    };

    BasisCache() = default;
    BasisCache(const BasisCache&) = delete;
    BasisCache& operator=(const BasisCache&) = delete;

    BasisCacheEntry compute(const BasisFunction& basis,
                            const quadrature::QuadratureRule& quad,
                            bool gradients,
                            bool hessians) const;

    std::shared_ptr<const BasisCacheEntry> get_or_compute_shared(
        const BasisFunction& basis,
        const quadrature::QuadratureRule& quad,
        bool gradients,
        bool hessians);

    mutable std::shared_mutex mutex_;
    std::unordered_map<BasisCacheKey, CacheSlot, BasisCacheKeyHash> slots_;
};

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_BASISCACHE_H
