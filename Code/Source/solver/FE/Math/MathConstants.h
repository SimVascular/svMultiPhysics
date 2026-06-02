#ifndef SVMP_FE_MATH_CONSTANTS_H
#define SVMP_FE_MATH_CONSTANTS_H

/**
 * @file MathConstants.h
 * @brief Mathematical constants and numerical tolerances for FE computations
 *
 * This header provides mathematical constants (π, e, √2, etc.) and numerical
 * tolerances used throughout the FE library. All constants are templated
 * to support different precision types.
 */

#include <cmath>
#include <limits>
#include <type_traits>
#include <algorithm>

namespace svmp {
namespace FE {
namespace math {

/**
 * @brief Mathematical constants templated by type
 * @tparam T The numeric type (float, double, long double)
 */
template<typename T>
struct Constants {
    static_assert(std::is_floating_point_v<T>,
                  "Constants only defined for floating-point types");

    // Mathematical constants
    static constexpr T pi           = T(3.14159265358979323846264338327950288419716939937510L);
    static constexpr T two_pi       = T(6.28318530717958647692528676655900576839433879875021L);
    static constexpr T half_pi      = T(1.57079632679489661923132169163975144209858469968755L);
    static constexpr T quarter_pi   = T(0.78539816339744830961566084581987572104929234984378L);
    static constexpr T inv_pi       = T(0.31830988618379067153776752674502872406891929148091L);
    static constexpr T inv_two_pi   = T(0.15915494309189533576888376337251436203445964574046L);

    static constexpr T e            = T(2.71828182845904523536028747135266249775724709369995L);
    static constexpr T log2e        = T(1.44269504088896340735992468100189213742664595415299L);
    static constexpr T log10e       = T(0.43429448190325182765112891891660508229439700580367L);
    static constexpr T ln2          = T(0.69314718055994530941723212145817656807550013436026L);
    static constexpr T ln10         = T(2.30258509299404568401799145468436420760110148862877L);

    static constexpr T sqrt2        = T(1.41421356237309504880168872420969807856967187537694L);
    static constexpr T sqrt3        = T(1.73205080756887729352744634150587236694280525381038L);
    static constexpr T inv_sqrt2    = T(0.70710678118654752440084436210484903928483593768847L);
    static constexpr T inv_sqrt3    = T(0.57735026918962576450914878050195745564760175127013L);

    // Golden ratio
    static constexpr T phi          = T(1.61803398874989484820458683436563811772030917980576L);

    // Degrees to radians conversion
    static constexpr T deg_to_rad   = pi / T(180);
    static constexpr T rad_to_deg   = T(180) / pi;
};

/**
 * @brief Numerical tolerances and machine epsilon
 * @tparam T The numeric type
 */
template<typename T>
struct Tolerances {
    static_assert(std::is_floating_point_v<T>,
                  "Tolerances only defined for floating-point types");

    // Machine epsilon
    static constexpr T epsilon      = std::numeric_limits<T>::epsilon();

    // Default tolerance (1000 * machine epsilon)
    static constexpr T tolerance    = T(1000) * epsilon;

    // Strict tolerance (10 * machine epsilon)
    static constexpr T strict       = T(10) * epsilon;

    // Loose tolerance (10000 * machine epsilon)
    static constexpr T loose        = T(10000) * epsilon;

    // Square root of epsilon (useful for finite differences)
    static inline const T sqrt_epsilon = std::sqrt(epsilon);

    // Cube root of epsilon (useful for numerical derivatives)
    static inline const T cbrt_epsilon = std::cbrt(epsilon);

    // Smallest positive normalized value
    static constexpr T min_positive = std::numeric_limits<T>::min();

    // Largest representable value
    static constexpr T max_value    = std::numeric_limits<T>::max();

    // Infinity
    static constexpr T infinity     = std::numeric_limits<T>::infinity();

    // Not-a-Number
    static constexpr T nan          = std::numeric_limits<T>::quiet_NaN();
};

/**
 * @brief Convenient aliases for common types
 */
template<typename T> inline constexpr T pi           = Constants<T>::pi;
template<typename T> inline constexpr T two_pi       = Constants<T>::two_pi;
template<typename T> inline constexpr T half_pi      = Constants<T>::half_pi;
template<typename T> inline constexpr T quarter_pi   = Constants<T>::quarter_pi;
template<typename T> inline constexpr T inv_pi       = Constants<T>::inv_pi;
template<typename T> inline constexpr T inv_two_pi   = Constants<T>::inv_two_pi;

template<typename T> inline constexpr T e            = Constants<T>::e;
template<typename T> inline constexpr T log2e        = Constants<T>::log2e;
template<typename T> inline constexpr T log10e       = Constants<T>::log10e;
template<typename T> inline constexpr T ln2          = Constants<T>::ln2;
template<typename T> inline constexpr T ln10         = Constants<T>::ln10;

template<typename T> inline constexpr T sqrt2        = Constants<T>::sqrt2;
template<typename T> inline constexpr T sqrt3        = Constants<T>::sqrt3;
template<typename T> inline constexpr T inv_sqrt2    = Constants<T>::inv_sqrt2;
template<typename T> inline constexpr T inv_sqrt3    = Constants<T>::inv_sqrt3;

template<typename T> inline constexpr T phi          = Constants<T>::phi;

template<typename T> inline constexpr T deg_to_rad   = Constants<T>::deg_to_rad;
template<typename T> inline constexpr T rad_to_deg   = Constants<T>::rad_to_deg;

template<typename T> inline constexpr T epsilon      = Tolerances<T>::epsilon;
template<typename T> inline constexpr T tolerance    = Tolerances<T>::tolerance;
template<typename T> inline constexpr T strict_tol   = Tolerances<T>::strict;
template<typename T> inline constexpr T loose_tol    = Tolerances<T>::loose;
template<typename T> inline const T sqrt_epsilon = Tolerances<T>::sqrt_epsilon;
template<typename T> inline const T cbrt_epsilon = Tolerances<T>::cbrt_epsilon;
template<typename T> inline constexpr T min_positive = Tolerances<T>::min_positive;
template<typename T> inline constexpr T max_value    = Tolerances<T>::max_value;
template<typename T> inline constexpr T infinity     = Tolerances<T>::infinity;

/**
 * @brief Comparison functions with tolerance
 */

/**
 * @brief Check if two values are approximately equal
 * @param a First value
 * @param b Second value
 * @param tol Tolerance (default: 1000 * epsilon)
 * @return true if |a - b| <= tol * max(|a|, |b|, 1)
 */
template<typename T>
inline constexpr bool approx_equal(T a, T b, T tol = tolerance<T>) {
    static_assert(std::is_floating_point_v<T>,
                  "approx_equal only defined for floating-point types");
    const T scale = std::max({std::abs(a), std::abs(b), T(1)});
    return std::abs(a - b) <= tol * scale;
}

/**
 * @brief Check if a value is approximately zero
 * @param a Value to check
 * @param tol Tolerance (default: 1000 * epsilon)
 * @return true if |a| <= tol
 */
template<typename T>
inline constexpr bool approx_zero(T a, T tol = tolerance<T>) {
    static_assert(std::is_floating_point_v<T>,
                  "approx_zero only defined for floating-point types");
    return std::abs(a) <= tol;
}

/**
 * @brief Check if a value is positive (greater than tolerance)
 * @param a Value to check
 * @param tol Tolerance (default: 1000 * epsilon)
 * @return true if a > tol
 */
template<typename T>
inline constexpr bool is_positive(T a, T tol = tolerance<T>) {
    static_assert(std::is_floating_point_v<T>,
                  "is_positive only defined for floating-point types");
    return a > tol;
}

/**
 * @brief Check if a value is negative (less than -tolerance)
 * @param a Value to check
 * @param tol Tolerance (default: 1000 * epsilon)
 * @return true if a < -tol
 */
template<typename T>
inline constexpr bool is_negative(T a, T tol = tolerance<T>) {
    static_assert(std::is_floating_point_v<T>,
                  "is_negative only defined for floating-point types");
    return a < -tol;
}

/**
 * @brief Check if a value is finite (not infinite or NaN)
 * @param a Value to check
 * @return true if value is finite
 */
template<typename T>
inline constexpr bool is_finite(T a) {
    static_assert(std::is_floating_point_v<T>,
                  "is_finite only defined for floating-point types");
    return std::isfinite(a);
}

/**
 * @brief Degrees to radians conversion
 * @param degrees Angle in degrees
 * @return Angle in radians
 */
template<typename T>
inline constexpr T to_radians(T degrees) {
    static_assert(std::is_floating_point_v<T>,
                  "to_radians only defined for floating-point types");
    return degrees * deg_to_rad<T>;
}

/**
 * @brief Radians to degrees conversion
 * @param radians Angle in radians
 * @return Angle in degrees
 */
template<typename T>
inline constexpr T to_degrees(T radians) {
    static_assert(std::is_floating_point_v<T>,
                  "to_degrees only defined for floating-point types");
    return radians * rad_to_deg<T>;
}

// =============================================================================
// Constants namespace for compatibility with test expectations
// =============================================================================
namespace constants {

// Mathematical constants (double precision defaults)
inline constexpr double PI         = Constants<double>::pi;
inline constexpr double PI_2       = Constants<double>::half_pi;
inline constexpr double PI_4       = Constants<double>::quarter_pi;
inline constexpr double TWO_PI     = Constants<double>::two_pi;
inline constexpr double INV_PI     = Constants<double>::inv_pi;

inline constexpr double E          = Constants<double>::e;
inline constexpr double LN_2       = Constants<double>::ln2;
inline constexpr double LN_10      = Constants<double>::ln10;
inline constexpr double LOG10_E    = Constants<double>::log10e;
inline constexpr double LOG2_E     = Constants<double>::log2e;

inline constexpr double SQRT_2     = Constants<double>::sqrt2;
inline constexpr double SQRT_3     = Constants<double>::sqrt3;
inline constexpr double SQRT_5     = 2.2360679774997896964091736687312L;
inline constexpr double INV_SQRT_2  = Constants<double>::inv_sqrt2;
inline constexpr double INV_SQRT_3  = Constants<double>::inv_sqrt3;

inline constexpr double PHI        = Constants<double>::phi;

// Angle conversion functions
template<typename T>
inline constexpr T deg_to_rad(T degrees) {
    return degrees * Constants<T>::deg_to_rad;
}

template<typename T>
inline constexpr T rad_to_deg(T radians) {
    return radians * Constants<T>::rad_to_deg;
}

// Templated tolerances
template<typename T>
inline constexpr T tolerance() {
    return Tolerances<T>::tolerance;
}

template<typename T>
inline constexpr T machine_epsilon() {
    return Tolerances<T>::epsilon;
}

// Additional constants and utility functions for tests
inline constexpr double DEFAULT_TOLERANCE = Tolerances<double>::tolerance;
inline constexpr double DEFAULT_REL_TOLERANCE = 1e-12;
inline constexpr double GEOMETRY_TOLERANCE = 1e-10;
inline constexpr double SOLVER_TOLERANCE = Tolerances<double>::strict;
inline constexpr double EPSILON = Tolerances<double>::epsilon;
inline constexpr double INF_VALUE = Tolerances<double>::infinity;  // Renamed from INFINITY
inline constexpr double NOT_A_NUMBER = Tolerances<double>::nan;  // Renamed from NAN
inline constexpr double MAX_DOUBLE = Tolerances<double>::max_value;
inline constexpr double MIN_DOUBLE = Tolerances<double>::min_positive;
inline constexpr double LOWEST_DOUBLE = -Tolerances<double>::max_value;

// Physical constants
inline constexpr double SPEED_OF_LIGHT = 299792458.0;         // m/s
inline constexpr double GRAVITATIONAL_CONSTANT = 6.67430e-11;  // m³/(kg·s²)
inline constexpr double PLANCK_CONSTANT = 6.62607015e-34;      // J·s
inline constexpr double AVOGADRO_NUMBER = 6.02214076e23;       // mol⁻¹
inline constexpr double BOLTZMANN_CONSTANT = 1.380649e-23;     // J/K
inline constexpr double STANDARD_GRAVITY = 9.80665;            // m/s²

// Float and long double versions
inline constexpr float PI_F = static_cast<float>(PI);
inline constexpr float E_F = static_cast<float>(E);
inline constexpr float SQRT_2_F = static_cast<float>(SQRT_2);
inline constexpr float EPSILON_F = Tolerances<float>::epsilon;

inline constexpr long double PI_L = static_cast<long double>(PI);
inline constexpr long double E_L = static_cast<long double>(E);
inline constexpr long double SQRT_2_L = static_cast<long double>(SQRT_2);
inline constexpr long double EPSILON_L = Tolerances<long double>::epsilon;

// Additional mathematical constants
inline constexpr double SQRT_PI = 1.7724538509055160272981674833411L;

// Utility functions
template<typename T>
inline constexpr int sign(T value) {
    return (T(0) < value) - (value < T(0));
}

template<typename T>
inline constexpr bool is_zero(T value, T tol = DEFAULT_TOLERANCE) {
    return std::abs(value) <= tol;
}

template<typename T>
inline bool near(T a, T b, T tol = DEFAULT_TOLERANCE) {
    return std::abs(a - b) <= tol;
}

template<typename T>
inline bool near_relative(T a, T b, T rel_tol = DEFAULT_REL_TOLERANCE) {
    T scale = std::max(std::abs(a), std::abs(b));
    return std::abs(a - b) <= rel_tol * scale;
}

template<typename T>
inline constexpr T clamp(T value, T min_val, T max_val) {
    return value < min_val ? min_val : (value > max_val ? max_val : value);
}

template<typename T>
inline constexpr T lerp(T a, T b, T t) {
    return a + t * (b - a);
}

template<typename T>
inline T safe_divide(T numerator, T denominator, T default_val = T(0)) {
    return is_zero(denominator) ? default_val : numerator / denominator;
}

template<typename T>
inline bool isinf(T value) {
    return std::isinf(value);
}

template<typename T>
inline bool isnan(T value) {
    return std::isnan(value);
}

} // namespace constants

// Physical constants for FE analysis
namespace physical_constants {

// Material properties (SI units)
inline constexpr double water_density = 1000.0;         // kg/m³
inline constexpr double steel_density = 7850.0;         // kg/m³
inline constexpr double aluminum_density = 2700.0;      // kg/m³

inline constexpr double water_viscosity = 0.001;        // Pa·s at 20°C
inline constexpr double air_viscosity = 1.81e-5;        // Pa·s at 20°C

inline constexpr double steel_youngs_modulus = 200e9;   // Pa
inline constexpr double aluminum_youngs_modulus = 70e9; // Pa

inline constexpr double steel_poisson_ratio = 0.3;      // dimensionless
inline constexpr double aluminum_poisson_ratio = 0.33;  // dimensionless

// Physical constants
inline constexpr double gravity = 9.80665;              // m/s²
inline constexpr double gas_constant = 8.314462618;     // J/(mol·K)
inline constexpr double boltzmann = 1.380649e-23;       // J/K
inline constexpr double avogadro = 6.02214076e23;       // mol⁻¹

} // namespace physical_constants

} // namespace math
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_MATH_CONSTANTS_H
