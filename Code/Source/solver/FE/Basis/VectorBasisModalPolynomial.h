/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#ifndef SVMP_FE_BASIS_VECTORBASISMODALPOLYNOMIAL_H
#define SVMP_FE_BASIS_VECTORBASISMODALPOLYNOMIAL_H

#include "Types.h"

#include <algorithm>
#include <array>
#include <vector>

namespace svmp {
namespace FE {
namespace basis {

struct VectorBasisModalTerm {
    int component{0}; // 0=x, 1=y, 2=z
    int px{0};
    int py{0};
    int pz{0};
    Real coefficient{Real(1)};
};

struct VectorBasisModalPolynomial {
    std::array<VectorBasisModalTerm, 4> terms{};
    int num_terms{0};
};

inline bool modal_terms_equal(const VectorBasisModalTerm& lhs,
                              const VectorBasisModalTerm& rhs) noexcept {
    return lhs.component == rhs.component &&
           lhs.px == rhs.px &&
           lhs.py == rhs.py &&
           lhs.pz == rhs.pz &&
           lhs.coefficient == rhs.coefficient;
}

inline bool modal_polynomials_equal(const VectorBasisModalPolynomial& lhs,
                                    const VectorBasisModalPolynomial& rhs) noexcept {
    if (lhs.num_terms != rhs.num_terms) {
        return false;
    }
    for (int term = 0; term < lhs.num_terms; ++term) {
        const auto index = static_cast<std::size_t>(term);
        if (!modal_terms_equal(lhs.terms[index], rhs.terms[index])) {
            return false;
        }
    }
    return true;
}

inline bool append_unique_modal_polynomial(
    std::vector<VectorBasisModalPolynomial>& polynomials,
    const VectorBasisModalPolynomial& polynomial) {
    const auto found = std::find_if(
        polynomials.begin(),
        polynomials.end(),
        [&](const VectorBasisModalPolynomial& existing) {
            return modal_polynomials_equal(existing, polynomial);
        });
    if (found != polynomials.end()) {
        return false;
    }
    polynomials.push_back(polynomial);
    return true;
}

} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_VECTORBASISMODALPOLYNOMIAL_H
