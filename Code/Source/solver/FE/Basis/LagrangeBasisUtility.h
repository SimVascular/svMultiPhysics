#ifndef SVMP_FE_BASIS_LAGRANGEBASISUTILITY_H
#define SVMP_FE_BASIS_LAGRANGEBASISUTILITY_H

// Private helper for LagrangeBasis internals.
// This header is only intended to be included after the FE basis scalar types
// are already available.

namespace svmp {
namespace FE {
namespace basis {
namespace detail {

inline constexpr Real equispaced_pm_one_coord(int i, int order) {
    if (order <= 0) {
        return Real(0);
    }
    return Real(-1) + Real(2) * static_cast<Real>(i) / static_cast<Real>(order);
}

} // namespace detail
} // namespace basis
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_BASIS_LAGRANGEBASISUTILITY_H
