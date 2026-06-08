/**
 * @file test_ExpressionOps.cpp
 * @brief Unit tests for ExpressionOps.h - expression template operators
 */

#include <gtest/gtest.h>
#include "FE/Math/ExpressionOps.h"
#include "FE/Math/Vector.h"
#include "FE/Math/Matrix.h"
#include <limits>
#include <cmath>
#include <complex>
#include <type_traits>

using namespace svmp::FE::math;
using namespace svmp::FE::math::detail::ops;

// Test fixture for ExpressionOps tests
class ExpressionOpsTest : public ::testing::Test {
protected:
    static constexpr double tolerance = 1e-14;

    void SetUp() override {}
    void TearDown() override {}

    template<typename T>
    bool approx_equal(T a, T b, T tol = tolerance) {
        return std::abs(a - b) <= tol;
    }
};

// =============================================================================
// Binary Operation Tests
// =============================================================================

TEST_F(ExpressionOpsTest, AddOperator) {
    Add op;

    // Integer addition
    EXPECT_EQ(op(5, 3), 8);
    EXPECT_EQ(op(-5, 3), -2);
    EXPECT_EQ(op(-5, -3), -8);

    // Floating point addition
    EXPECT_DOUBLE_EQ(op(3.14, 2.86), 6.0);
    EXPECT_DOUBLE_EQ(op(-1.5, 2.5), 1.0);

    // Mixed types
    auto result = op(3, 2.5);
    EXPECT_TRUE((std::is_same_v<decltype(result), double>));
    EXPECT_DOUBLE_EQ(result, 5.5);
}

TEST_F(ExpressionOpsTest, SubOperator) {
    Sub op;

    // Integer subtraction
    EXPECT_EQ(op(5, 3), 2);
    EXPECT_EQ(op(3, 5), -2);
    EXPECT_EQ(op(-5, -3), -2);

    // Floating point subtraction
    EXPECT_DOUBLE_EQ(op(5.5, 2.5), 3.0);
    EXPECT_DOUBLE_EQ(op(2.5, 5.5), -3.0);

    // Mixed types
    auto result = op(5.5, 2);
    EXPECT_TRUE((std::is_same_v<decltype(result), double>));
    EXPECT_DOUBLE_EQ(result, 3.5);
}

TEST_F(ExpressionOpsTest, MulOperator) {
    Mul op;

    // Integer multiplication
    EXPECT_EQ(op(5, 3), 15);
    EXPECT_EQ(op(-5, 3), -15);
    EXPECT_EQ(op(-5, -3), 15);

    // Floating point multiplication
    EXPECT_DOUBLE_EQ(op(2.5, 4.0), 10.0);
    EXPECT_DOUBLE_EQ(op(-2.5, 4.0), -10.0);

    // Zero multiplication
    EXPECT_EQ(op(0, 100), 0);
    EXPECT_DOUBLE_EQ(op(0.0, 3.14), 0.0);

    // Mixed types
    auto result = op(3, 2.5);
    EXPECT_TRUE((std::is_same_v<decltype(result), double>));
    EXPECT_DOUBLE_EQ(result, 7.5);
}

TEST_F(ExpressionOpsTest, DivOperator) {
    Div op;

    // Integer division
    EXPECT_EQ(op(10, 2), 5);
    EXPECT_EQ(op(10, 3), 3);  // Integer division truncates
    EXPECT_EQ(op(-10, 2), -5);

    // Floating point division
    EXPECT_DOUBLE_EQ(op(10.0, 2.0), 5.0);
    EXPECT_DOUBLE_EQ(op(10.0, 3.0), 10.0/3.0);
    EXPECT_DOUBLE_EQ(op(-10.0, 2.0), -5.0);

    // Mixed types
    auto result = op(10.0, 3);
    EXPECT_TRUE((std::is_same_v<decltype(result), double>));
    EXPECT_DOUBLE_EQ(result, 10.0/3.0);
}

// =============================================================================
// Unary Operation Tests
// =============================================================================

TEST_F(ExpressionOpsTest, NegateOperator) {
    Negate op;

    // Integer negation
    EXPECT_EQ(op(5), -5);
    EXPECT_EQ(op(-5), 5);
    EXPECT_EQ(op(0), 0);

    // Floating point negation
    EXPECT_DOUBLE_EQ(op(3.14), -3.14);
    EXPECT_DOUBLE_EQ(op(-2.71), 2.71);
    EXPECT_DOUBLE_EQ(op(0.0), 0.0);

    // Type preservation
    auto int_result = op(5);
    EXPECT_TRUE((std::is_same_v<decltype(int_result), int>));

    auto double_result = op(5.0);
    EXPECT_TRUE((std::is_same_v<decltype(double_result), double>));
}

TEST_F(ExpressionOpsTest, AbsOperator) {
    Abs op;

    // Integer absolute value
    EXPECT_EQ(op(5), 5);
    EXPECT_EQ(op(-5), 5);
    EXPECT_EQ(op(0), 0);

    // Floating point absolute value
    EXPECT_DOUBLE_EQ(op(3.14), 3.14);
    EXPECT_DOUBLE_EQ(op(-3.14), 3.14);
    EXPECT_DOUBLE_EQ(op(0.0), 0.0);

    // Special cases
    EXPECT_DOUBLE_EQ(op(-0.0), 0.0);

    // Type preservation
    auto int_result = op(-5);
    EXPECT_TRUE((std::is_same_v<decltype(int_result), int>));

    auto double_result = op(-5.0);
    EXPECT_TRUE((std::is_same_v<decltype(double_result), double>));
}

TEST_F(ExpressionOpsTest, SqrtOperator) {
    Sqrt op;

    // Perfect squares
    EXPECT_DOUBLE_EQ(op(4.0), 2.0);
    EXPECT_DOUBLE_EQ(op(9.0), 3.0);
    EXPECT_DOUBLE_EQ(op(16.0), 4.0);
    EXPECT_DOUBLE_EQ(op(25.0), 5.0);

    // Non-perfect squares
    EXPECT_DOUBLE_EQ(op(2.0), std::sqrt(2.0));
    EXPECT_DOUBLE_EQ(op(3.0), std::sqrt(3.0));

    // Special cases
    EXPECT_DOUBLE_EQ(op(0.0), 0.0);
    EXPECT_DOUBLE_EQ(op(1.0), 1.0);

    // Type conversion
    auto result = op(4);  // Integer input
    EXPECT_DOUBLE_EQ(result, 2.0);
}

// =============================================================================
// Constexpr Tests
// =============================================================================

TEST_F(ExpressionOpsTest, ConstexprOperators) {
    // Test that operators can be used in constexpr contexts
    constexpr Add add_op;
    constexpr Sub sub_op;
    constexpr Mul mul_op;
    constexpr Div div_op;
    constexpr Negate neg_op;

    // Compile-time evaluation
    constexpr auto sum = add_op(3, 4);
    constexpr auto diff = sub_op(7, 3);
    constexpr auto prod = mul_op(3, 4);
    constexpr auto quot = div_op(12, 3);
    constexpr auto neg = neg_op(5);

    EXPECT_EQ(sum, 7);
    EXPECT_EQ(diff, 4);
    EXPECT_EQ(prod, 12);
    EXPECT_EQ(quot, 4);
    EXPECT_EQ(neg, -5);

    // Static assertions to verify compile-time evaluation
    static_assert(add_op(2, 3) == 5);
    static_assert(sub_op(5, 2) == 3);
    static_assert(mul_op(3, 4) == 12);
    static_assert(div_op(10, 2) == 5);
    static_assert(neg_op(3) == -3);
}

// =============================================================================
// Type Deduction Tests
// =============================================================================

TEST_F(ExpressionOpsTest, TypeDeduction) {
    Add add_op;
    Sub sub_op;
    Mul mul_op;
    Div div_op;

    // int + int -> int
    auto int_result = add_op(3, 4);
    EXPECT_TRUE((std::is_same_v<decltype(int_result), int>));

    // double + double -> double
    auto double_result = add_op(3.0, 4.0);
    EXPECT_TRUE((std::is_same_v<decltype(double_result), double>));

    // int + double -> double
    auto mixed_result1 = add_op(3, 4.0);
    EXPECT_TRUE((std::is_same_v<decltype(mixed_result1), double>));

    // double + int -> double
    auto mixed_result2 = add_op(3.0, 4);
    EXPECT_TRUE((std::is_same_v<decltype(mixed_result2), double>));

    // float + double -> double
    auto float_double_result = add_op(3.0f, 4.0);
    EXPECT_TRUE((std::is_same_v<decltype(float_double_result), double>));
}

// =============================================================================
// Complex Expression Tests
// =============================================================================

TEST_F(ExpressionOpsTest, ChainedOperations) {
    Add add_op;
    Sub sub_op;
    Mul mul_op;
    Div div_op;
    Negate neg_op;

    // Simulate complex expression: -(a + b) * c / d
    double a = 2.0, b = 3.0, c = 4.0, d = 2.0;

    auto sum = add_op(a, b);       // 5.0
    auto negated = neg_op(sum);    // -5.0
    auto product = mul_op(negated, c);  // -20.0
    auto result = div_op(product, d);   // -10.0

    EXPECT_DOUBLE_EQ(result, -10.0);
}

TEST_F(ExpressionOpsTest, MixedPrecisionChain) {
    Add add_op;
    Mul mul_op;

    // Mixed precision chain
    int a = 2;
    float b = 3.5f;
    double c = 1.5;

    auto step1 = add_op(a, b);    // int + float -> float (5.5f)
    auto step2 = mul_op(step1, c); // float + double -> double (8.25)

    EXPECT_TRUE((std::is_same_v<decltype(step2), double>));
    EXPECT_DOUBLE_EQ(step2, 8.25);
}

// =============================================================================
// Operator Integration with Vector/Matrix Tests
// =============================================================================

TEST_F(ExpressionOpsTest, VectorIntegration) {
    Vector<double, 3> v1{1.0, 2.0, 3.0};
    Vector<double, 3> v2{4.0, 5.0, 6.0};

    // Test that operators work correctly in vector expressions
    Vector<double, 3> sum = v1 + v2;
    Vector<double, 3> diff = v1 - v2;
    Vector<double, 3> neg = -v1;
    Vector<double, 3> scaled = v1 * 2.0;

    EXPECT_DOUBLE_EQ(sum[0], 5.0);
    EXPECT_DOUBLE_EQ(sum[1], 7.0);
    EXPECT_DOUBLE_EQ(sum[2], 9.0);

    EXPECT_DOUBLE_EQ(diff[0], -3.0);
    EXPECT_DOUBLE_EQ(diff[1], -3.0);
    EXPECT_DOUBLE_EQ(diff[2], -3.0);

    EXPECT_DOUBLE_EQ(neg[0], -1.0);
    EXPECT_DOUBLE_EQ(neg[1], -2.0);
    EXPECT_DOUBLE_EQ(neg[2], -3.0);

    EXPECT_DOUBLE_EQ(scaled[0], 2.0);
    EXPECT_DOUBLE_EQ(scaled[1], 4.0);
    EXPECT_DOUBLE_EQ(scaled[2], 6.0);
}

TEST_F(ExpressionOpsTest, MatrixIntegration) {
    Matrix<double, 2, 2> m1{{1.0, 2.0}, {3.0, 4.0}};
    Matrix<double, 2, 2> m2{{5.0, 6.0}, {7.0, 8.0}};

    // Test that operators work correctly in matrix expressions
    Matrix<double, 2, 2> sum = m1 + m2;
    Matrix<double, 2, 2> diff = m1 - m2;
    Matrix<double, 2, 2> neg = -m1;
    Matrix<double, 2, 2> scaled = m1 * 2.0;

    EXPECT_DOUBLE_EQ(sum(0, 0), 6.0);
    EXPECT_DOUBLE_EQ(sum(0, 1), 8.0);
    EXPECT_DOUBLE_EQ(sum(1, 0), 10.0);
    EXPECT_DOUBLE_EQ(sum(1, 1), 12.0);

    EXPECT_DOUBLE_EQ(diff(0, 0), -4.0);
    EXPECT_DOUBLE_EQ(diff(0, 1), -4.0);
    EXPECT_DOUBLE_EQ(diff(1, 0), -4.0);
    EXPECT_DOUBLE_EQ(diff(1, 1), -4.0);

    EXPECT_DOUBLE_EQ(neg(0, 0), -1.0);
    EXPECT_DOUBLE_EQ(neg(0, 1), -2.0);
    EXPECT_DOUBLE_EQ(neg(1, 0), -3.0);
    EXPECT_DOUBLE_EQ(neg(1, 1), -4.0);

    EXPECT_DOUBLE_EQ(scaled(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(scaled(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(scaled(1, 0), 6.0);
    EXPECT_DOUBLE_EQ(scaled(1, 1), 8.0);
}

// =============================================================================
// Edge Cases and Special Values Tests
// =============================================================================

TEST_F(ExpressionOpsTest, SpecialFloatingPointValues) {
    Add add_op;
    Sub sub_op;
    Mul mul_op;
    Div div_op;
    Abs abs_op;
    Negate neg_op;

    // Infinity handling
    double inf = std::numeric_limits<double>::infinity();
    EXPECT_DOUBLE_EQ(add_op(inf, 1.0), inf);
    EXPECT_DOUBLE_EQ(sub_op(inf, 1.0), inf);
    EXPECT_DOUBLE_EQ(mul_op(inf, 2.0), inf);
    EXPECT_DOUBLE_EQ(div_op(inf, 2.0), inf);
    EXPECT_DOUBLE_EQ(abs_op(inf), inf);
    EXPECT_DOUBLE_EQ(neg_op(inf), -inf);

    // NaN handling
    double nan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_TRUE(std::isnan(add_op(nan, 1.0)));
    EXPECT_TRUE(std::isnan(sub_op(nan, 1.0)));
    EXPECT_TRUE(std::isnan(mul_op(nan, 2.0)));
    EXPECT_TRUE(std::isnan(div_op(nan, 2.0)));
    EXPECT_TRUE(std::isnan(abs_op(nan)));
    EXPECT_TRUE(std::isnan(neg_op(nan)));

    // Division by zero
    EXPECT_DOUBLE_EQ(div_op(1.0, 0.0), inf);
    EXPECT_DOUBLE_EQ(div_op(-1.0, 0.0), -inf);
    EXPECT_TRUE(std::isnan(div_op(0.0, 0.0)));
}

TEST_F(ExpressionOpsTest, LargeAndSmallValues) {
    Add add_op;
    Mul mul_op;

    // Large values
    double large = 1e308;
    double result = add_op(large, large);
    EXPECT_TRUE(std::isinf(result));  // Overflow to infinity

    // Small values
    double tiny = std::numeric_limits<double>::min();
    double tiny_result = mul_op(tiny, 0.5);
    EXPECT_GT(tiny_result, 0.0);  // Should still be positive
    EXPECT_LT(tiny_result, tiny);  // But smaller

    // Denormalized numbers
    double denorm = std::numeric_limits<double>::denorm_min();
    double denorm_result = add_op(denorm, denorm);
    EXPECT_EQ(denorm_result, 2.0 * denorm);
}

// =============================================================================
// SFINAE and Compile-time Constraint Tests
// =============================================================================

TEST_F(ExpressionOpsTest, SFINAECompatibility) {
    // Test that operators work with any arithmetic types
    Add add_op;

    // Various integer types
    EXPECT_EQ(add_op(int8_t(3), int8_t(4)), 7);
    EXPECT_EQ(add_op(int16_t(100), int16_t(200)), 300);
    EXPECT_EQ(add_op(int32_t(1000), int32_t(2000)), 3000);
    EXPECT_EQ(add_op(int64_t(10000), int64_t(20000)), 30000);

    // Unsigned types
    EXPECT_EQ(add_op(uint8_t(3), uint8_t(4)), 7u);
    EXPECT_EQ(add_op(uint16_t(100), uint16_t(200)), 300u);
    EXPECT_EQ(add_op(uint32_t(1000), uint32_t(2000)), 3000u);

    // Floating point types
    EXPECT_FLOAT_EQ(add_op(3.0f, 4.0f), 7.0f);
    EXPECT_DOUBLE_EQ(add_op(3.0, 4.0), 7.0);

    // Long double
    long double ld1 = 3.0L;
    long double ld2 = 4.0L;
    EXPECT_DOUBLE_EQ(add_op(ld1, ld2), 7.0L);
}

// =============================================================================
// Template Instantiation Tests
// =============================================================================

TEST_F(ExpressionOpsTest, TemplateInstantiations) {
    // Test that operators can be instantiated with various types
    Add add_op;
    Sub sub_op;
    Mul mul_op;
    Div div_op;
    Abs abs_op;
    Sqrt sqrt_op;
    Negate neg_op;

    // Custom types that support arithmetic operations
    struct CustomNumber {
        double value;
        CustomNumber(double v) : value(v) {}
        CustomNumber operator+(const CustomNumber& other) const { return CustomNumber(value + other.value); }
        CustomNumber operator-(const CustomNumber& other) const { return CustomNumber(value - other.value); }
        CustomNumber operator*(const CustomNumber& other) const { return CustomNumber(value * other.value); }
        CustomNumber operator/(const CustomNumber& other) const { return CustomNumber(value / other.value); }
        CustomNumber operator-() const { return CustomNumber(-value); }
        bool operator==(const CustomNumber& other) const { return value == other.value; }
    };

    CustomNumber cn1(3.0);
    CustomNumber cn2(4.0);

    auto cn_sum = add_op(cn1, cn2);
    EXPECT_EQ(cn_sum.value, 7.0);

    auto cn_diff = sub_op(cn1, cn2);
    EXPECT_EQ(cn_diff.value, -1.0);

    auto cn_prod = mul_op(cn1, cn2);
    EXPECT_EQ(cn_prod.value, 12.0);

    auto cn_quot = div_op(cn1, cn2);
    EXPECT_EQ(cn_quot.value, 0.75);

    auto cn_neg = neg_op(cn1);
    EXPECT_EQ(cn_neg.value, -3.0);
}

// =============================================================================
// Complex Number Support Tests
// =============================================================================

TEST_F(ExpressionOpsTest, ComplexNumberSupport) {
    Add add_op;
    Sub sub_op;
    Mul mul_op;
    Div div_op;
    Negate neg_op;

    std::complex<double> c1(3.0, 4.0);
    std::complex<double> c2(1.0, 2.0);

    auto c_sum = add_op(c1, c2);
    EXPECT_DOUBLE_EQ(c_sum.real(), 4.0);
    EXPECT_DOUBLE_EQ(c_sum.imag(), 6.0);

    auto c_diff = sub_op(c1, c2);
    EXPECT_DOUBLE_EQ(c_diff.real(), 2.0);
    EXPECT_DOUBLE_EQ(c_diff.imag(), 2.0);

    auto c_prod = mul_op(c1, c2);
    EXPECT_DOUBLE_EQ(c_prod.real(), -5.0);  // (3+4i)(1+2i) = 3+6i+4i+8i² = 3+10i-8 = -5+10i
    EXPECT_DOUBLE_EQ(c_prod.imag(), 10.0);

    auto c_neg = neg_op(c1);
    EXPECT_DOUBLE_EQ(c_neg.real(), -3.0);
    EXPECT_DOUBLE_EQ(c_neg.imag(), -4.0);
}
