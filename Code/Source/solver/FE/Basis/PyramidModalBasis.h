#ifndef SVMP_FE_BASIS_PYRAMIDMODALBASIS_H
#define SVMP_FE_BASIS_PYRAMIDMODALBASIS_H

// Shared rational/modal pyramid helpers for scalar complete-family and spectral
// pyramid bases. The degenerate z=1 top plane is evaluated by its apex limit;
// callers that reject non-apex top-plane queries must validate before calling.

#include "BasisFunction.h"
#include "BasisTolerance.h"
#include <algorithm>
#include <cmath>
#include <vector>

namespace svmp {
namespace FE {
namespace basis {
namespace pyramid_modal {

struct Term {
    int px{0};
    int py{0};
    int pz{0};
    int denom_power{0};
};

struct EvaluationPoint {
    Real x{Real(0)};
    Real y{Real(0)};
    Real z{Real(0)};
    Real t{Real(1)};
    bool top_plane{false};
    std::vector<Real> x_powers;
    std::vector<Real> y_powers;
    std::vector<Real> z_powers;
    std::vector<Real> t_powers;
};

inline std::vector<Term> build_terms(int order) {
    std::vector<Term> terms;
    terms.reserve(static_cast<std::size_t>((order + 1) * (order + 2) *
                                           (2 * order + 3) / 6));
    for (int pz = 0; pz <= order; ++pz) {
        const int n = order - pz;
        for (int py = 0; py <= n; ++py) {
            for (int px = 0; px <= n; ++px) {
                terms.push_back(Term{px, py, pz, std::min(px, py)});
            }
        }
    }
    return terms;
}

inline bool on_degenerate_top_plane(const math::Vector<Real, 3>& xi,
                                    Real tolerance = detail::basis_scaled_tolerance()) {
    return std::abs(Real(1) - xi[2]) <= tolerance;
}

inline void fill_powers(Real base, int max_power, std::vector<Real>& powers) {
    powers.assign(static_cast<std::size_t>(max_power + 1), Real(1));
    for (int p = 1; p <= max_power; ++p) {
        powers[static_cast<std::size_t>(p)] =
            powers[static_cast<std::size_t>(p - 1)] * base;
    }
}

inline void prepare_evaluation_point(const math::Vector<Real, 3>& xi,
                                     int max_px,
                                     int max_py,
                                     int max_pz,
                                     int max_denom_power,
                                     EvaluationPoint& point) {
    point.x = xi[0];
    point.y = xi[1];
    point.z = xi[2];
    point.t = Real(1) - point.z;
    point.top_plane = on_degenerate_top_plane(xi);

    fill_powers(point.x, std::max(max_px, 0), point.x_powers);
    fill_powers(point.y, std::max(max_py, 0), point.y_powers);
    fill_powers(point.z, std::max(max_pz, 0), point.z_powers);
    if (point.top_plane) [[unlikely]] {
        point.t_powers.assign(1u, Real(1));
    } else {
        fill_powers(point.t, std::max(max_denom_power + 2, 0), point.t_powers);
    }
}

inline void prepare_evaluation_point(const std::vector<Term>& terms,
                                     const math::Vector<Real, 3>& xi,
                                     EvaluationPoint& point) {
    int max_px = 0;
    int max_py = 0;
    int max_pz = 0;
    int max_denom_power = 0;
    for (const Term& term : terms) {
        max_px = std::max(max_px, term.px);
        max_py = std::max(max_py, term.py);
        max_pz = std::max(max_pz, term.pz);
        max_denom_power = std::max(max_denom_power, term.denom_power);
    }
    prepare_evaluation_point(xi, max_px, max_py, max_pz, max_denom_power, point);
}

inline void evaluate_term(const Term& term,
                          const EvaluationPoint& point,
                          Real& value,
                          Gradient* gradient = nullptr,
                          Hessian* hessian = nullptr) {
    const auto pow_x = [&](int p) -> Real {
        return point.x_powers[static_cast<std::size_t>(p)];
    };
    const auto pow_y = [&](int p) -> Real {
        return point.y_powers[static_cast<std::size_t>(p)];
    };
    const auto pow_z = [&](int p) -> Real {
        return point.z_powers[static_cast<std::size_t>(p)];
    };
    const auto pow_t = [&](int p) -> Real {
        return point.t_powers[static_cast<std::size_t>(p)];
    };

    if (point.top_plane) [[unlikely]] {
        if (term.px == 0 && term.py == 0) {
            value = pow_z(term.pz);
        } else {
            value = Real(0);
        }
        if (gradient != nullptr) {
            *gradient = Gradient{};
            if (term.px == 0 && term.py == 0 && term.pz > 0) {
                (*gradient)[2] = static_cast<Real>(term.pz) * pow_z(term.pz - 1);
            }
        }
        if (hessian != nullptr) {
            *hessian = Hessian{};
            if (term.px == 0 && term.py == 0 && term.pz > 1) {
                (*hessian)(2, 2) =
                    static_cast<Real>(term.pz * (term.pz - 1)) *
                    pow_z(term.pz - 2);
            }
        }
        return;
    }

    const Real base = pow_x(term.px) * pow_y(term.py) * pow_z(term.pz);
    const Real denom = pow_t(term.denom_power);
    value = base / denom;

    if (gradient != nullptr) {
        *gradient = Gradient{};
        if (term.px > 0) {
            (*gradient)[0] =
                static_cast<Real>(term.px) * pow_x(term.px - 1) *
                pow_y(term.py) * pow_z(term.pz) / denom;
        }
        if (term.py > 0) {
            (*gradient)[1] =
                static_cast<Real>(term.py) * pow_x(term.px) *
                pow_y(term.py - 1) * pow_z(term.pz) / denom;
        }

        Real gz = Real(0);
        if (term.pz > 0) {
            gz += static_cast<Real>(term.pz) * pow_x(term.px) *
                  pow_y(term.py) * pow_z(term.pz - 1) / denom;
        }
        if (term.denom_power > 0) {
            gz += static_cast<Real>(term.denom_power) * base / pow_t(term.denom_power + 1);
        }
        (*gradient)[2] = gz;
    }

    if (hessian == nullptr) {
        return;
    }

    *hessian = Hessian{};
    if (term.px > 1) {
        (*hessian)(0, 0) =
            static_cast<Real>(term.px * (term.px - 1)) *
            pow_x(term.px - 2) * pow_y(term.py) * pow_z(term.pz) / denom;
    }
    if (term.py > 1) {
        (*hessian)(1, 1) =
            static_cast<Real>(term.py * (term.py - 1)) *
            pow_x(term.px) * pow_y(term.py - 2) * pow_z(term.pz) / denom;
    }
    if (term.px > 0 && term.py > 0) {
        const Real hxy =
            static_cast<Real>(term.px * term.py) *
            pow_x(term.px - 1) * pow_y(term.py - 1) * pow_z(term.pz) / denom;
        (*hessian)(0, 1) = hxy;
        (*hessian)(1, 0) = hxy;
    }

    if (term.px > 0) {
        Real hxz =
            static_cast<Real>(term.px) * pow_x(term.px - 1) *
            pow_y(term.py) / denom;
        if (term.pz > 0) {
            hxz *= static_cast<Real>(term.pz) * pow_z(term.pz - 1);
        } else {
            hxz = Real(0);
        }
        if (term.denom_power > 0) {
            hxz += static_cast<Real>(term.px * term.denom_power) *
                   pow_x(term.px - 1) * pow_y(term.py) *
                   pow_z(term.pz) / pow_t(term.denom_power + 1);
        }
        (*hessian)(0, 2) = hxz;
        (*hessian)(2, 0) = hxz;
    }

    if (term.py > 0) {
        Real hyz =
            static_cast<Real>(term.py) * pow_x(term.px) *
            pow_y(term.py - 1) / denom;
        if (term.pz > 0) {
            hyz *= static_cast<Real>(term.pz) * pow_z(term.pz - 1);
        } else {
            hyz = Real(0);
        }
        if (term.denom_power > 0) {
            hyz += static_cast<Real>(term.py * term.denom_power) *
                   pow_x(term.px) * pow_y(term.py - 1) *
                   pow_z(term.pz) / pow_t(term.denom_power + 1);
        }
        (*hessian)(1, 2) = hyz;
        (*hessian)(2, 1) = hyz;
    }

    Real hzz = Real(0);
    if (term.pz > 1) {
        hzz += static_cast<Real>(term.pz * (term.pz - 1)) *
               pow_x(term.px) * pow_y(term.py) * pow_z(term.pz - 2) / denom;
    }
    if (term.pz > 0 && term.denom_power > 0) {
        hzz += static_cast<Real>(2 * term.pz * term.denom_power) *
               pow_x(term.px) * pow_y(term.py) *
               pow_z(term.pz - 1) / pow_t(term.denom_power + 1);
    }
    if (term.denom_power > 0) {
        hzz += static_cast<Real>(term.denom_power * (term.denom_power + 1)) *
               base / pow_t(term.denom_power + 2);
    }
    (*hessian)(2, 2) = hzz;
}

inline void evaluate_term(const Term& term,
                          const math::Vector<Real, 3>& xi,
                          Real& value,
                          Gradient* gradient = nullptr,
                          Hessian* hessian = nullptr) {
    EvaluationPoint point;
    prepare_evaluation_point(
        xi, term.px, term.py, term.pz, term.denom_power, point);
    evaluate_term(term, point, value, gradient, hessian);
}

} // namespace pyramid_modal
} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_PYRAMIDMODALBASIS_H
