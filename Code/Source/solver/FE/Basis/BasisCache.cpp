/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#include "BasisCache.h"
#include <utility>

namespace svmp {
namespace FE {
namespace basis {

namespace {

QuadratureCacheKey make_quadrature_cache_key(const quadrature::QuadratureRule& quad) noexcept {
    const auto fingerprint = quad.point_fingerprint();
    return QuadratureCacheKey{fingerprint.dimension,
                              fingerprint.num_points,
                              fingerprint.points_hash_a,
                              fingerprint.points_hash_b};
}

void mix_hash_word(std::uint64_t word,
                   std::uint64_t& hash_a,
                   std::uint64_t& hash_b) noexcept {
    hash_a ^= word + 0x9e3779b97f4a7c15ULL + (hash_a << 6u) + (hash_a >> 2u);
    hash_b ^= (word + 0xbf58476d1ce4e5b9ULL) + (hash_b << 7u) + (hash_b >> 3u);
}

std::pair<std::uint64_t, std::uint64_t>
identity_fingerprint(const std::string& identity) noexcept {
    std::uint64_t hash_a = 0xa4093822299f31d0ULL;
    std::uint64_t hash_b = 0x082efa98ec4e6c89ULL;
    mix_hash_word(static_cast<std::uint64_t>(identity.size()), hash_a, hash_b);
    for (const char c : identity) {
        mix_hash_word(static_cast<std::uint64_t>(static_cast<unsigned char>(c)), hash_a, hash_b);
    }
    return {hash_a, hash_b};
}

BasisCacheKey make_basis_cache_key(const BasisFunction& basis,
                                   const quadrature::QuadratureRule& quad,
                                   bool gradients,
                                   bool hessians) {
    StructuralBasisKey structural_key{
        basis.basis_type(),
        basis.element_type(),
        basis.dimension(),
        basis.order(),
        basis.size(),
        basis.is_vector_valued(),
        make_quadrature_cache_key(quad),
        gradients,
        hessians
    };

    BasisCacheKey key;
    const bool uses_basis_identity = !basis.cache_identity_is_structural();
    if (!uses_basis_identity) {
        key.value = structural_key;
        return key;
    }

    std::vector<std::uint64_t> basis_identity_words;
    const bool uses_structured_identity = basis.cache_identity_words(basis_identity_words);
    if (!uses_structured_identity) {
        basis_identity_words.clear();
    }
    const std::string basis_identity =
        uses_structured_identity ? std::string{} : basis.cache_identity();
    BasisIdentityFingerprint cached_identity_hash{};
    const bool has_cached_identity_hash =
        uses_structured_identity &&
        basis.cache_identity_fingerprint(cached_identity_hash.hash_a,
                                         cached_identity_hash.hash_b);
    const auto identity_hash = uses_structured_identity
        ? has_cached_identity_hash
              ? std::pair<std::uint64_t, std::uint64_t>{
                    cached_identity_hash.hash_a,
                    cached_identity_hash.hash_b}
              : [&basis_identity_words] {
                    const auto fingerprint =
                        compute_basis_identity_fingerprint(basis_identity_words);
                    return std::pair<std::uint64_t, std::uint64_t>{
                        fingerprint.hash_a,
                        fingerprint.hash_b};
                }()
        : identity_fingerprint(basis_identity);
    key.value = ParameterizedBasisKey{
        structural_key,
        uses_structured_identity,
        identity_hash.first,
        identity_hash.second,
        std::move(basis_identity_words),
        basis_identity
    };
    return key;
}

} // namespace

BasisCache& BasisCache::instance() {
    static BasisCache cache;
    return cache;
}

const BasisCacheEntry& BasisCache::get_or_compute(
    const BasisFunction& basis,
    const quadrature::QuadratureRule& quad,
    bool gradients,
    bool hessians) {
    return *get_or_compute_shared(basis, quad, gradients, hessians);
}

std::shared_ptr<const BasisCacheEntry> BasisCache::get_or_compute_shared(
    const BasisFunction& basis,
    const quadrature::QuadratureRule& quad,
    bool gradients,
    bool hessians) {
    const BasisCacheKey key = make_basis_cache_key(basis, quad, gradients, hessians);

    // Warm path: shared (reader) lock allows concurrent cache hits.
    {
        std::shared_lock<std::shared_mutex> read_lock(mutex_);
        auto it = slots_.find(key);
        if (it != slots_.end() && it->second.entry) {
            return it->second.entry;
        }
    }

    std::shared_ptr<InFlightComputation> in_flight;
    bool owner = false;
    {
        std::unique_lock<std::shared_mutex> write_lock(mutex_);
        auto& slot = slots_[key];
        if (slot.entry) {
            return slot.entry;
        }

        if (!slot.pending) {
            in_flight = std::make_shared<InFlightComputation>();
            slot.pending = in_flight;
            owner = true;
        } else {
            in_flight = slot.pending;
        }
    }

    if (!owner) {
        std::unique_lock<std::mutex> wait_lock(in_flight->mutex);
        in_flight->ready_cv.wait(wait_lock, [&in_flight] { return in_flight->ready; });
        if (in_flight->exception) {
            std::rethrow_exception(in_flight->exception);
        }
        return in_flight->entry;
    }

    try {
        auto entry = std::make_shared<BasisCacheEntry>(compute(basis, quad, gradients, hessians));
        {
            std::unique_lock<std::shared_mutex> write_lock(mutex_);
            auto slot_it = slots_.find(key);
            if (slot_it == slots_.end()) {
                slot_it = slots_.emplace(key, CacheSlot{}).first;
            }
            auto& slot = slot_it->second;
            if (slot.entry) {
                entry = slot.entry;
            } else {
                slot.entry = entry;
            }
            if (slot.pending == in_flight) {
                slot.pending.reset();
            }
        }
        {
            std::lock_guard<std::mutex> ready_lock(in_flight->mutex);
            in_flight->entry = entry;
            in_flight->ready = true;
        }
        in_flight->ready_cv.notify_all();
        return entry;
    } catch (...) {
        {
            std::lock_guard<std::mutex> ready_lock(in_flight->mutex);
            in_flight->exception = std::current_exception();
            in_flight->ready = true;
        }
        {
            std::unique_lock<std::shared_mutex> write_lock(mutex_);
            auto slot_it = slots_.find(key);
            if (slot_it != slots_.end() && slot_it->second.pending == in_flight) {
                slot_it->second.pending.reset();
                if (!slot_it->second.entry) {
                    slots_.erase(slot_it);
                }
            }
        }
        in_flight->ready_cv.notify_all();
        throw;
    }
}

const BasisCacheEntry& BasisCache::prewarm(
    const BasisFunction& basis,
    const quadrature::QuadratureRule& quad,
    bool gradients,
    bool hessians) {
    return get_or_compute(basis, quad, gradients, hessians);
}

BasisCacheHandle BasisCache::prewarm_handle(
    const BasisFunction& basis,
    const quadrature::QuadratureRule& quad,
    bool gradients,
    bool hessians) {
    return BasisCacheHandle(get_or_compute_shared(basis, quad, gradients, hessians));
}

BasisCacheEntry BasisCache::compute_uncached(
    const BasisFunction& basis,
    const quadrature::QuadratureRule& quad,
    bool gradients,
    bool hessians) const {
    return compute(basis, quad, gradients, hessians);
}

void BasisCache::clear() {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    for (auto it = slots_.begin(); it != slots_.end();) {
        if (it->second.pending) {
            it->second.entry.reset();
            ++it;
        } else {
            it = slots_.erase(it);
        }
    }
}

std::size_t BasisCache::size() const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    std::size_t completed = 0;
    for (const auto& [key, slot] : slots_) {
        (void)key;
        if (slot.entry) {
            ++completed;
        }
    }
    return completed;
}

BasisCacheEntry BasisCache::compute(const BasisFunction& basis,
                                    const quadrature::QuadratureRule& quad,
                                    bool gradients,
                                    bool hessians) const {
    BasisCacheEntry entry;
    const auto& points = quad.points();
    entry.num_qpts = points.size();
    entry.num_dofs = basis.size();

    const bool vector_basis = basis.is_vector_valued();
    if (!vector_basis) {
        entry.scalar_values.assign(entry.num_dofs * entry.num_qpts, Real(0));
        if (gradients) {
            entry.gradients.assign(entry.num_dofs * 3u * entry.num_qpts, Real(0));
        }
        if (hessians) {
            entry.hessians.assign(entry.num_dofs * 9u * entry.num_qpts, Real(0));
        }
    } else {
        entry.vector_values_xyz.assign(entry.num_dofs * 3u * entry.num_qpts, Real(0));
        if (gradients && basis.supports_vector_jacobians()) {
            entry.vector_jacobians.assign(entry.num_dofs * 9u * entry.num_qpts, Real(0));
        }
        if (gradients && basis.supports_curl()) {
            entry.vector_curls_xyz.assign(entry.num_dofs * 3u * entry.num_qpts, Real(0));
        }
        if (gradients && basis.supports_divergence()) {
            entry.vector_divergence.assign(entry.num_dofs * entry.num_qpts, Real(0));
        }
    }

    if (vector_basis) {
        if (entry.num_dofs > 0 && entry.num_qpts > 0) {
            basis.evaluate_vector_at_quadrature_points(
                points,
                entry.vector_values_xyz.data(),
                entry.vector_jacobians.empty() ? nullptr : entry.vector_jacobians.data(),
                entry.vector_curls_xyz.empty() ? nullptr : entry.vector_curls_xyz.data(),
                entry.vector_divergence.empty() ? nullptr : entry.vector_divergence.data());
        }
        return entry;
    }

    if (entry.num_dofs > 0 && entry.num_qpts > 0) {
        basis.fill_scalar_cache_entry(points,
                                      entry.num_qpts,
                                      entry.scalar_values.data(),
                                      gradients ? entry.gradients.data() : nullptr,
                                      hessians ? entry.hessians.data() : nullptr);
    }

    return entry;
}
} // namespace basis
} // namespace FE
} // namespace svmp
