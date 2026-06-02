/**
 * @file test_MathConstants.cpp
 * @brief Unit tests for MathConstants.h - mathematical constants and tolerances
 */

#include <gtest/gtest.h>
#include "FE/Math/MathConstants.h"
#include <cmath>
#include <limits>
#include <type_traits>

using namespace svmp::FE::math;

// Test fixture for MathConstants tests
class MathConstantsTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

// =============================================================================
// Mathematical Constants Tests
// =============================================================================

TEST_F(MathConstantsTest, PiConstants) {
    // Test PI value
    EXPECT_NEAR(constants::PI, 3.14159265358979323846, 1e-15);

    // Test PI/2
    EXPECT_NEAR(constants::PI_2, constants::PI / 2.0, 1e-15);

    // Test PI/4
    EXPECT_NEAR(constants::PI_4, constants::PI / 4.0, 1e-15);

    // Test 2*PI
    EXPECT_NEAR(constants::TWO_PI, 2.0 * constants::PI, 1e-15);

    // Test 1/PI
    EXPECT_NEAR(constants::INV_PI, 1.0 / constants::PI, 1e-15);

    // Test sqrt(PI)
    EXPECT_NEAR(constants::SQRT_PI, std::sqrt(constants::PI), 1e-15);
}

TEST_F(MathConstantsTest, EulerConstant) {
    // Test e (Euler's number)
    EXPECT_NEAR(constants::E, std::exp(1.0), 1e-15);

    // Test ln(2)
    EXPECT_NEAR(constants::LN_2, std::log(2.0), 1e-15);

    // Test ln(10)
    EXPECT_NEAR(constants::LN_10, std::log(10.0), 1e-15);

    // Test log10(e)
    EXPECT_NEAR(constants::LOG10_E, std::log10(constants::E), 1e-15);

    // Test log2(e)
    EXPECT_NEAR(constants::LOG2_E, std::log2(constants::E), 1e-15);
}

TEST_F(MathConstantsTest, SquareRootConstants) {
    // Test sqrt(2)
    EXPECT_NEAR(constants::SQRT_2, std::sqrt(2.0), 1e-15);

    // Test sqrt(3)
    EXPECT_NEAR(constants::SQRT_3, std::sqrt(3.0), 1e-15);

    // Test sqrt(5)
    EXPECT_NEAR(constants::SQRT_5, std::sqrt(5.0), 1e-15);

    // Test 1/sqrt(2)
    EXPECT_NEAR(constants::INV_SQRT_2, 1.0 / std::sqrt(2.0), 1e-15);

    // Test 1/sqrt(3)
    EXPECT_NEAR(constants::INV_SQRT_3, 1.0 / std::sqrt(3.0), 1e-15);
}

TEST_F(MathConstantsTest, GoldenRatio) {
    // Test golden ratio φ = (1 + sqrt(5))/2
    EXPECT_NEAR(constants::PHI, (1.0 + std::sqrt(5.0)) / 2.0, 1e-15);

    // Property: φ² = φ + 1
    EXPECT_NEAR(constants::PHI * constants::PHI, constants::PHI + 1.0, 1e-14);

    // Property: 1/φ = φ - 1
    EXPECT_NEAR(1.0 / constants::PHI, constants::PHI - 1.0, 1e-14);
}

// =============================================================================
// Angle Conversion Tests
// =============================================================================

TEST_F(MathConstantsTest, DegreesToRadians) {
    // Test common conversions
    EXPECT_NEAR(constants::deg_to_rad(0.0), 0.0, 1e-15);
    EXPECT_NEAR(constants::deg_to_rad(90.0), constants::PI_2, 1e-15);
    EXPECT_NEAR(constants::deg_to_rad(180.0), constants::PI, 1e-15);
    EXPECT_NEAR(constants::deg_to_rad(270.0), 3.0 * constants::PI_2, 1e-15);
    EXPECT_NEAR(constants::deg_to_rad(360.0), constants::TWO_PI, 1e-15);

    // Test negative angles
    EXPECT_NEAR(constants::deg_to_rad(-90.0), -constants::PI_2, 1e-15);
    EXPECT_NEAR(constants::deg_to_rad(-180.0), -constants::PI, 1e-15);

    // Test arbitrary angle
    EXPECT_NEAR(constants::deg_to_rad(45.0), constants::PI_4, 1e-15);
    EXPECT_NEAR(constants::deg_to_rad(30.0), constants::PI / 6.0, 1e-15);
    EXPECT_NEAR(constants::deg_to_rad(60.0), constants::PI / 3.0, 1e-15);
}

TEST_F(MathConstantsTest, RadiansToDegrees) {
    // Test common conversions
    EXPECT_NEAR(constants::rad_to_deg(0.0), 0.0, 1e-13);
    EXPECT_NEAR(constants::rad_to_deg(constants::PI_2), 90.0, 1e-13);
    EXPECT_NEAR(constants::rad_to_deg(constants::PI), 180.0, 1e-13);
    EXPECT_NEAR(constants::rad_to_deg(constants::TWO_PI), 360.0, 1e-13);

    // Test negative angles
    EXPECT_NEAR(constants::rad_to_deg(-constants::PI), -180.0, 1e-13);

    // Test round-trip conversion
    double angle_deg = 123.456;
    double angle_rad = constants::deg_to_rad(angle_deg);
    double back_to_deg = constants::rad_to_deg(angle_rad);
    EXPECT_NEAR(back_to_deg, angle_deg, 1e-13);
}

// =============================================================================
// Machine Precision Tests
// =============================================================================

TEST_F(MathConstantsTest, MachineEpsilon) {
    // Test double precision epsilon
    EXPECT_EQ(constants::EPSILON, std::numeric_limits<double>::epsilon());

    // Test float precision epsilon
    EXPECT_EQ(constants::EPSILON_F, std::numeric_limits<float>::epsilon());

    // Verify epsilon is the smallest value such that 1.0 + epsilon != 1.0
    double one_plus_eps = 1.0 + constants::EPSILON;
    double one_plus_half_eps = 1.0 + constants::EPSILON / 2.0;

    EXPECT_NE(one_plus_eps, 1.0);
    EXPECT_EQ(one_plus_half_eps, 1.0);
}

TEST_F(MathConstantsTest, NumericalLimits) {
    // Test infinity
    EXPECT_TRUE(std::isinf(constants::INF_VALUE));
    EXPECT_GT(constants::INF_VALUE, std::numeric_limits<double>::max());

    // Test NaN
    EXPECT_TRUE(std::isnan(constants::NOT_A_NUMBER));
    EXPECT_NE(constants::NOT_A_NUMBER, constants::NOT_A_NUMBER);  // NaN != NaN

    // Test max/min values
    EXPECT_EQ(constants::MAX_DOUBLE, std::numeric_limits<double>::max());
    EXPECT_EQ(constants::MIN_DOUBLE, std::numeric_limits<double>::min());
    EXPECT_EQ(constants::LOWEST_DOUBLE, std::numeric_limits<double>::lowest());
}

// =============================================================================
// Tolerance Tests
// =============================================================================

TEST_F(MathConstantsTest, DefaultTolerances) {
    // Test default absolute tolerance
    EXPECT_GT(constants::DEFAULT_TOLERANCE, 0.0);
    EXPECT_LT(constants::DEFAULT_TOLERANCE, 1e-10);

    // Test default relative tolerance
    EXPECT_GT(constants::DEFAULT_REL_TOLERANCE, 0.0);
    EXPECT_LT(constants::DEFAULT_REL_TOLERANCE, 1e-10);

    // Test solver tolerance
    EXPECT_GT(constants::SOLVER_TOLERANCE, 0.0);
    EXPECT_LE(constants::SOLVER_TOLERANCE, constants::DEFAULT_TOLERANCE);

    // Test geometry tolerance (typically larger)
    EXPECT_GT(constants::GEOMETRY_TOLERANCE, 0.0);
    EXPECT_GE(constants::GEOMETRY_TOLERANCE, constants::DEFAULT_TOLERANCE);
}

TEST_F(MathConstantsTest, ToleranceComparison) {
    double a = 1.0;
    double b = 1.0 + constants::DEFAULT_TOLERANCE / 2.0;
    double c = 1.0 + constants::DEFAULT_TOLERANCE * 2.0;

    // Values within tolerance should be considered equal
    EXPECT_TRUE(constants::near(a, b, constants::DEFAULT_TOLERANCE));

    // Values outside tolerance should not be equal
    EXPECT_FALSE(constants::near(a, c, constants::DEFAULT_TOLERANCE));

    // Test relative tolerance
    double large_a = 1e10;
    double large_b = large_a * (1.0 + constants::DEFAULT_REL_TOLERANCE / 2.0);
    double large_c = large_a * (1.0 + constants::DEFAULT_REL_TOLERANCE * 2.0);

    EXPECT_TRUE(constants::near_relative(large_a, large_b, constants::DEFAULT_REL_TOLERANCE));
    EXPECT_FALSE(constants::near_relative(large_a, large_c, constants::DEFAULT_REL_TOLERANCE));
}

TEST_F(MathConstantsTest, ZeroComparison) {
    // Test near zero detection
    EXPECT_TRUE(constants::is_zero(0.0));
    EXPECT_TRUE(constants::is_zero(constants::DEFAULT_TOLERANCE / 2.0));
    EXPECT_FALSE(constants::is_zero(constants::DEFAULT_TOLERANCE * 2.0));

    // Test with negative values
    EXPECT_TRUE(constants::is_zero(-constants::DEFAULT_TOLERANCE / 2.0));
    EXPECT_FALSE(constants::is_zero(-constants::DEFAULT_TOLERANCE * 2.0));
}

// =============================================================================
// Physical Constants Tests
// =============================================================================

TEST_F(MathConstantsTest, PhysicalConstants) {
    // Test speed of light (m/s)
    EXPECT_NEAR(constants::SPEED_OF_LIGHT, 299792458.0, 1.0);

    // Test gravitational constant (m³/kg/s²)
    EXPECT_NEAR(constants::GRAVITATIONAL_CONSTANT, 6.67430e-11, 1e-16);

    // Test standard gravity (m/s²)
    EXPECT_NEAR(constants::STANDARD_GRAVITY, 9.80665, 1e-10);

    // Test Planck constant (J⋅s)
    EXPECT_NEAR(constants::PLANCK_CONSTANT, 6.62607015e-34, 1e-42);

    // Test Boltzmann constant (J/K)
    EXPECT_NEAR(constants::BOLTZMANN_CONSTANT, 1.380649e-23, 1e-29);

    // Test Avogadro's number (1/mol)
    EXPECT_NEAR(constants::AVOGADRO_NUMBER, 6.02214076e23, 1e15);
}

// =============================================================================
// Compile-Time Constants Tests
// =============================================================================

TEST_F(MathConstantsTest, CompileTimeConstants) {
    // Test that constants are constexpr (compile-time)
    constexpr double pi = constants::PI;
    constexpr double e = constants::E;
    constexpr double sqrt2 = constants::SQRT_2;

    EXPECT_EQ(pi, constants::PI);
    EXPECT_EQ(e, constants::E);
    EXPECT_EQ(sqrt2, constants::SQRT_2);

    // Test compile-time functions
    constexpr double angle_rad = constants::deg_to_rad(90.0);
    EXPECT_NEAR(angle_rad, constants::PI_2, 1e-15);

    constexpr double angle_deg = constants::rad_to_deg(constants::PI);
    EXPECT_NEAR(angle_deg, 180.0, 1e-13);
}

// =============================================================================
// Type Traits Tests
// =============================================================================

TEST_F(MathConstantsTest, TypedConstants) {
    // Test float versions
    EXPECT_NEAR(constants::PI_F, static_cast<float>(constants::PI), 1e-7f);
    EXPECT_NEAR(constants::E_F, static_cast<float>(constants::E), 1e-7f);
    EXPECT_NEAR(constants::SQRT_2_F, static_cast<float>(constants::SQRT_2), 1e-7f);

    // Test long double versions
    EXPECT_NEAR(constants::PI_L, static_cast<long double>(constants::PI), 1e-18L);
    EXPECT_NEAR(constants::E_L, static_cast<long double>(constants::E), 1e-18L);
}

// =============================================================================
// Special Functions Tests
// =============================================================================

TEST_F(MathConstantsTest, SignFunction) {
    // Test sign function
    EXPECT_EQ(constants::sign(5.0), 1);
    EXPECT_EQ(constants::sign(-5.0), -1);
    EXPECT_EQ(constants::sign(0.0), 0);

    // Test with very small values
    EXPECT_EQ(constants::sign(constants::EPSILON), 1);
    EXPECT_EQ(constants::sign(-constants::EPSILON), -1);

    // Test with infinity
    EXPECT_EQ(constants::sign(constants::INF_VALUE), 1);
    EXPECT_EQ(constants::sign(-constants::INF_VALUE), -1);
}

TEST_F(MathConstantsTest, SafeDivision) {
    // Test safe division
    EXPECT_NEAR(constants::safe_divide(10.0, 2.0), 5.0, 1e-15);
    EXPECT_NEAR(constants::safe_divide(1.0, 3.0), 1.0/3.0, 1e-15);

    // Test division by zero returns default
    EXPECT_EQ(constants::safe_divide(1.0, 0.0, 999.0), 999.0);
    EXPECT_EQ(constants::safe_divide(1.0, constants::EPSILON/2.0, -1.0), -1.0);

    // Test division by near-zero
    double tiny = constants::DEFAULT_TOLERANCE / 10.0;
    EXPECT_EQ(constants::safe_divide(1.0, tiny, 0.0), 0.0);
}

// =============================================================================
// Utility Functions Tests
// =============================================================================

TEST_F(MathConstantsTest, ClampFunction) {
    // Test clamping
    EXPECT_EQ(constants::clamp(5.0, 0.0, 10.0), 5.0);
    EXPECT_EQ(constants::clamp(-5.0, 0.0, 10.0), 0.0);
    EXPECT_EQ(constants::clamp(15.0, 0.0, 10.0), 10.0);

    // Test with same min/max
    EXPECT_EQ(constants::clamp(5.0, 3.0, 3.0), 3.0);

    // Test with infinity
    EXPECT_EQ(constants::clamp(constants::INF_VALUE, 0.0, 10.0), 10.0);
    EXPECT_EQ(constants::clamp(-constants::INF_VALUE, 0.0, 10.0), 0.0);
}

TEST_F(MathConstantsTest, LerpFunction) {
    // Test linear interpolation
    EXPECT_NEAR(constants::lerp(0.0, 10.0, 0.0), 0.0, 1e-15);
    EXPECT_NEAR(constants::lerp(0.0, 10.0, 1.0), 10.0, 1e-15);
    EXPECT_NEAR(constants::lerp(0.0, 10.0, 0.5), 5.0, 1e-15);
    EXPECT_NEAR(constants::lerp(0.0, 10.0, 0.25), 2.5, 1e-15);

    // Test extrapolation
    EXPECT_NEAR(constants::lerp(0.0, 10.0, -0.5), -5.0, 1e-15);
    EXPECT_NEAR(constants::lerp(0.0, 10.0, 1.5), 15.0, 1e-15);

    // Test with negative range
    EXPECT_NEAR(constants::lerp(-10.0, -5.0, 0.5), -7.5, 1e-15);
}
