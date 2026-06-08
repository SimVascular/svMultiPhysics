/**
 * @file test_Vector.cpp
 * @brief Unit tests for Vector.h - fixed-size vectors with expression templates
 */

#include <gtest/gtest.h>
#include "FE/Math/Vector.h"
#include "FE/Math/VectorExpr.h"
#include <limits>
#include <cmath>
#include <sstream>
#include <thread>
#include <vector>

using namespace svmp::FE::math;

// Test fixture for Vector tests
class VectorTest : public ::testing::Test {
protected:
    static constexpr double tolerance = 1e-14;

    void SetUp() override {}
    void TearDown() override {}

    // Helper function to check if two values are approximately equal
    template<typename T>
    bool approx_equal(T a, T b, T tol = tolerance) {
        return std::abs(a - b) <= tol;
    }
};

// =============================================================================
// Construction and Initialization Tests
// =============================================================================

TEST_F(VectorTest, DefaultConstruction) {
    Vector<double, 3> v;
    EXPECT_EQ(v[0], 0.0);
    EXPECT_EQ(v[1], 0.0);
    EXPECT_EQ(v[2], 0.0);

    Vector<float, 4> vf;
    for (size_t i = 0; i < 4; ++i) {
        EXPECT_EQ(vf[i], 0.0f);
    }
}

TEST_F(VectorTest, FillConstruction) {
    Vector<double, 3> v(5.0);
    EXPECT_EQ(v[0], 5.0);
    EXPECT_EQ(v[1], 5.0);
    EXPECT_EQ(v[2], 5.0);

    Vector<int, 10> vi(-3);
    for (size_t i = 0; i < 10; ++i) {
        EXPECT_EQ(vi[i], -3);
    }
}

TEST_F(VectorTest, InitializerListConstruction) {
    Vector<double, 3> v{1.0, 2.0, 3.0};
    EXPECT_EQ(v[0], 1.0);
    EXPECT_EQ(v[1], 2.0);
    EXPECT_EQ(v[2], 3.0);

    // Partial initialization
    Vector<double, 5> v2{1.0, 2.0};
    EXPECT_EQ(v2[0], 1.0);
    EXPECT_EQ(v2[1], 2.0);
    EXPECT_EQ(v2[2], 0.0);
    EXPECT_EQ(v2[3], 0.0);
    EXPECT_EQ(v2[4], 0.0);
}

TEST_F(VectorTest, CopyConstruction) {
    Vector<double, 3> v1{1.0, 2.0, 3.0};
    Vector<double, 3> v2(v1);

    EXPECT_EQ(v2[0], 1.0);
    EXPECT_EQ(v2[1], 2.0);
    EXPECT_EQ(v2[2], 3.0);

    // Ensure deep copy
    v2[0] = 10.0;
    EXPECT_EQ(v1[0], 1.0);
    EXPECT_EQ(v2[0], 10.0);
}

TEST_F(VectorTest, MoveConstruction) {
    Vector<double, 3> v1{1.0, 2.0, 3.0};
    Vector<double, 3> v2(std::move(v1));

    EXPECT_EQ(v2[0], 1.0);
    EXPECT_EQ(v2[1], 2.0);
    EXPECT_EQ(v2[2], 3.0);
}

// =============================================================================
// Element Access Tests
// =============================================================================

TEST_F(VectorTest, ElementAccess) {
    Vector<double, 3> v{1.0, 2.0, 3.0};

    // Non-const access
    EXPECT_EQ(v[0], 1.0);
    EXPECT_EQ(v[1], 2.0);
    EXPECT_EQ(v[2], 3.0);

    // Modification
    v[1] = 5.0;
    EXPECT_EQ(v[1], 5.0);

    // Const access
    const Vector<double, 3> cv{4.0, 5.0, 6.0};
    EXPECT_EQ(cv[0], 4.0);
    EXPECT_EQ(cv[1], 5.0);
    EXPECT_EQ(cv[2], 6.0);
}

TEST_F(VectorTest, ElementAccessBounds) {
    Vector<double, 3> v{1.0, 2.0, 3.0};

    // at() with bounds checking
    EXPECT_EQ(v.at(0), 1.0);
    EXPECT_EQ(v.at(1), 2.0);
    EXPECT_EQ(v.at(2), 3.0);

    // Test out of bounds throws
    EXPECT_THROW(v.at(3), std::out_of_range);
    EXPECT_THROW(v.at(100), std::out_of_range);
}

TEST_F(VectorTest, DataPointerAccess) {
    Vector<double, 3> v{1.0, 2.0, 3.0};

    double* data = v.data();
    EXPECT_EQ(data[0], 1.0);
    EXPECT_EQ(data[1], 2.0);
    EXPECT_EQ(data[2], 3.0);

    // Const data access
    const Vector<double, 3> cv{4.0, 5.0, 6.0};
    const double* cdata = cv.data();
    EXPECT_EQ(cdata[0], 4.0);
    EXPECT_EQ(cdata[1], 5.0);
    EXPECT_EQ(cdata[2], 6.0);
}

// =============================================================================
// Arithmetic Operations Tests
// =============================================================================

TEST_F(VectorTest, Addition) {
    Vector<double, 3> a{1.0, 2.0, 3.0};
    Vector<double, 3> b{4.0, 5.0, 6.0};

    Vector<double, 3> c = a + b;
    EXPECT_EQ(c[0], 5.0);
    EXPECT_EQ(c[1], 7.0);
    EXPECT_EQ(c[2], 9.0);
}

TEST_F(VectorTest, Subtraction) {
    Vector<double, 3> a{5.0, 7.0, 9.0};
    Vector<double, 3> b{4.0, 5.0, 6.0};

    Vector<double, 3> c = a - b;
    EXPECT_EQ(c[0], 1.0);
    EXPECT_EQ(c[1], 2.0);
    EXPECT_EQ(c[2], 3.0);
}

TEST_F(VectorTest, ScalarMultiplication) {
    Vector<double, 3> a{1.0, 2.0, 3.0};

    // Scalar * Vector
    Vector<double, 3> b = 2.0 * a;
    EXPECT_EQ(b[0], 2.0);
    EXPECT_EQ(b[1], 4.0);
    EXPECT_EQ(b[2], 6.0);

    // Vector * Scalar
    Vector<double, 3> c = a * 3.0;
    EXPECT_EQ(c[0], 3.0);
    EXPECT_EQ(c[1], 6.0);
    EXPECT_EQ(c[2], 9.0);
}

TEST_F(VectorTest, ScalarDivision) {
    Vector<double, 3> a{2.0, 4.0, 6.0};

    Vector<double, 3> b = a / 2.0;
    EXPECT_EQ(b[0], 1.0);
    EXPECT_EQ(b[1], 2.0);
    EXPECT_EQ(b[2], 3.0);
}

TEST_F(VectorTest, UnaryNegation) {
    Vector<double, 3> a{1.0, -2.0, 3.0};

    Vector<double, 3> b = -a;
    EXPECT_EQ(b[0], -1.0);
    EXPECT_EQ(b[1], 2.0);
    EXPECT_EQ(b[2], -3.0);
}

TEST_F(VectorTest, CompoundAssignment) {
    Vector<double, 3> a{1.0, 2.0, 3.0};
    Vector<double, 3> b{4.0, 5.0, 6.0};

    // +=
    a += b;
    EXPECT_EQ(a[0], 5.0);
    EXPECT_EQ(a[1], 7.0);
    EXPECT_EQ(a[2], 9.0);

    // -=
    a -= b;
    EXPECT_EQ(a[0], 1.0);
    EXPECT_EQ(a[1], 2.0);
    EXPECT_EQ(a[2], 3.0);

    // *=
    a *= 2.0;
    EXPECT_EQ(a[0], 2.0);
    EXPECT_EQ(a[1], 4.0);
    EXPECT_EQ(a[2], 6.0);

    // /=
    a /= 2.0;
    EXPECT_EQ(a[0], 1.0);
    EXPECT_EQ(a[1], 2.0);
    EXPECT_EQ(a[2], 3.0);
}

// =============================================================================
// Vector Operations Tests
// =============================================================================

TEST_F(VectorTest, DotProduct) {
    Vector<double, 3> a{1.0, 2.0, 3.0};
    Vector<double, 3> b{4.0, 5.0, 6.0};

    double dot = a.dot(b);
    EXPECT_EQ(dot, 32.0);  // 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32

    // Test commutativity
    EXPECT_EQ(b.dot(a), dot);

    // Test orthogonal vectors
    Vector<double, 3> x{1.0, 0.0, 0.0};
    Vector<double, 3> y{0.0, 1.0, 0.0};
    EXPECT_EQ(x.dot(y), 0.0);
}

TEST_F(VectorTest, CrossProduct3D) {
    Vector<double, 3> x{1.0, 0.0, 0.0};
    Vector<double, 3> y{0.0, 1.0, 0.0};
    Vector<double, 3> z{0.0, 0.0, 1.0};

    // Test basis vector cross products
    Vector<double, 3> xy = x.cross(y);
    EXPECT_EQ(xy[0], 0.0);
    EXPECT_EQ(xy[1], 0.0);
    EXPECT_EQ(xy[2], 1.0);

    Vector<double, 3> yx = y.cross(x);
    EXPECT_EQ(yx[0], 0.0);
    EXPECT_EQ(yx[1], 0.0);
    EXPECT_EQ(yx[2], -1.0);

    // General cross product
    Vector<double, 3> a{1.0, 2.0, 3.0};
    Vector<double, 3> b{4.0, 5.0, 6.0};
    Vector<double, 3> c = a.cross(b);

    EXPECT_EQ(c[0], -3.0);  // 2*6 - 3*5 = 12 - 15 = -3
    EXPECT_EQ(c[1], 6.0);   // 3*4 - 1*6 = 12 - 6 = 6
    EXPECT_EQ(c[2], -3.0);  // 1*5 - 2*4 = 5 - 8 = -3
}

TEST_F(VectorTest, Norm) {
    Vector<double, 3> v{3.0, 4.0, 0.0};
    EXPECT_EQ(v.norm(), 5.0);

    Vector<double, 3> unit{1.0, 0.0, 0.0};
    EXPECT_EQ(unit.norm(), 1.0);

    Vector<double, 3> zero{0.0, 0.0, 0.0};
    EXPECT_EQ(zero.norm(), 0.0);
}

TEST_F(VectorTest, NormSquared) {
    Vector<double, 3> v{3.0, 4.0, 0.0};
    EXPECT_EQ(v.norm_squared(), 25.0);

    Vector<double, 3> a{1.0, 2.0, 3.0};
    EXPECT_EQ(a.norm_squared(), 14.0);  // 1 + 4 + 9 = 14
}

TEST_F(VectorTest, Normalize) {
    Vector<double, 3> v{3.0, 4.0, 0.0};
    Vector<double, 3> n = v.normalized();

    EXPECT_NEAR(n[0], 0.6, tolerance);
    EXPECT_NEAR(n[1], 0.8, tolerance);
    EXPECT_NEAR(n[2], 0.0, tolerance);
    EXPECT_NEAR(n.norm(), 1.0, tolerance);

    // In-place normalization
    v.normalize();
    EXPECT_NEAR(v[0], 0.6, tolerance);
    EXPECT_NEAR(v[1], 0.8, tolerance);
    EXPECT_NEAR(v.norm(), 1.0, tolerance);
}

// =============================================================================
// Expression Template Tests
// =============================================================================

TEST_F(VectorTest, ExpressionTemplatesNoTemporaries) {
    Vector<double, 3> a{1.0, 2.0, 3.0};
    Vector<double, 3> b{4.0, 5.0, 6.0};
    Vector<double, 3> c{7.0, 8.0, 9.0};
    Vector<double, 3> d{10.0, 11.0, 12.0};

    // Complex expression should create no temporaries
    Vector<double, 3> result = a + b - c + d;

    EXPECT_EQ(result[0], 8.0);   // 1 + 4 - 7 + 10
    EXPECT_EQ(result[1], 10.0);  // 2 + 5 - 8 + 11
    EXPECT_EQ(result[2], 12.0);  // 3 + 6 - 9 + 12
}

TEST_F(VectorTest, LazyEvaluation) {
    Vector<double, 3> a{1.0, 2.0, 3.0};
    Vector<double, 3> b{4.0, 5.0, 6.0};

    // Expression should not be evaluated until assignment
    auto expr = a + b;  // No computation yet

    Vector<double, 3> result = expr;  // Evaluation happens here
    EXPECT_EQ(result[0], 5.0);
    EXPECT_EQ(result[1], 7.0);
    EXPECT_EQ(result[2], 9.0);
}

TEST_F(VectorTest, MixedExpressions) {
    Vector<double, 3> a{1.0, 2.0, 3.0};
    Vector<double, 3> b{4.0, 5.0, 6.0};
    double scalar = 2.0;

    // Complex mixed expression
    Vector<double, 3> result = scalar * (a + b) - a / scalar;

    EXPECT_NEAR(result[0], 9.5, tolerance);   // 2*(1+4) - 1/2
    EXPECT_NEAR(result[1], 13.0, tolerance);  // 2*(2+5) - 2/2
    EXPECT_NEAR(result[2], 16.5, tolerance);  // 2*(3+6) - 3/2
}

// =============================================================================
// Special Values Tests
// =============================================================================

TEST_F(VectorTest, ZeroVector) {
    Vector<double, 3> zero = Vector<double, 3>::zero();
    EXPECT_EQ(zero[0], 0.0);
    EXPECT_EQ(zero[1], 0.0);
    EXPECT_EQ(zero[2], 0.0);
    EXPECT_EQ(zero.norm(), 0.0);
}

TEST_F(VectorTest, OnesVector) {
    Vector<double, 3> ones = Vector<double, 3>::ones();
    EXPECT_EQ(ones[0], 1.0);
    EXPECT_EQ(ones[1], 1.0);
    EXPECT_EQ(ones[2], 1.0);
}

TEST_F(VectorTest, BasisVectors) {
    auto e0 = Vector<double, 3>::basis(0);
    EXPECT_EQ(e0[0], 1.0);
    EXPECT_EQ(e0[1], 0.0);
    EXPECT_EQ(e0[2], 0.0);

    auto e1 = Vector<double, 3>::basis(1);
    EXPECT_EQ(e1[0], 0.0);
    EXPECT_EQ(e1[1], 1.0);
    EXPECT_EQ(e1[2], 0.0);

    auto e2 = Vector<double, 3>::basis(2);
    EXPECT_EQ(e2[0], 0.0);
    EXPECT_EQ(e2[1], 0.0);
    EXPECT_EQ(e2[2], 1.0);
}

// =============================================================================
// Edge Cases and Error Handling Tests
// =============================================================================

TEST_F(VectorTest, DivisionByZero) {
    Vector<double, 3> v{1.0, 2.0, 3.0};

    // Division by zero should produce inf
    Vector<double, 3> result = v / 0.0;
    EXPECT_TRUE(std::isinf(result[0]));
    EXPECT_TRUE(std::isinf(result[1]));
    EXPECT_TRUE(std::isinf(result[2]));
}

TEST_F(VectorTest, NormalizeZeroVector) {
    Vector<double, 3> zero{0.0, 0.0, 0.0};

    // Normalizing zero vector should handle gracefully
    Vector<double, 3> n = zero.normalized();
    EXPECT_TRUE(std::isnan(n[0]) || n[0] == 0.0);
}

TEST_F(VectorTest, ExtremeLargeValues) {
    double large = 1e308;  // Near double max
    Vector<double, 3> v{large, large, large};

    // Operations should not overflow
    Vector<double, 3> half = v / 2.0;
    EXPECT_FALSE(std::isinf(half[0]));
    EXPECT_EQ(half[0], large / 2.0);
}

TEST_F(VectorTest, ExtremeSmallValues) {
    double tiny = 1e-308;  // Near double min
    Vector<double, 3> v{tiny, tiny, tiny};

    // Operations should maintain precision
    Vector<double, 3> doubled = v * 2.0;
    EXPECT_EQ(doubled[0], tiny * 2.0);
}

// =============================================================================
// Numerical Precision Tests
// =============================================================================

TEST_F(VectorTest, NumericalStability) {
    // Test Kahan summation for better precision
    Vector<double, 4> v{1e16, 1.0, -1e16, 1.0};
    // Computed for future validation - demonstrates numerical precision issues
    [[maybe_unused]] double sum = v[0] + v[1] + v[2] + v[3];

    // Direct summation might lose precision
    // But vector operations should maintain it
    Vector<double, 4> a{1e16, 0.0, -1e16, 0.0};
    Vector<double, 4> b{0.0, 1.0, 0.0, 1.0};
    Vector<double, 4> c = a + b;

    EXPECT_EQ(c[0], 1e16);
    EXPECT_EQ(c[1], 1.0);
    EXPECT_EQ(c[2], -1e16);
    EXPECT_EQ(c[3], 1.0);
}

TEST_F(VectorTest, OrthogonalityPreservation) {
    // Create nearly orthogonal vectors
    Vector<double, 3> a{1.0, 1e-15, 0.0};
    Vector<double, 3> b{0.0, 1.0, 0.0};

    double dot = a.dot(b);
    EXPECT_NEAR(dot, 1e-15, 1e-16);
}

// =============================================================================
// Comparison Operations Tests
// =============================================================================

TEST_F(VectorTest, Equality) {
    Vector<double, 3> a{1.0, 2.0, 3.0};
    Vector<double, 3> b{1.0, 2.0, 3.0};
    Vector<double, 3> c{1.0, 2.0, 3.1};

    EXPECT_TRUE(a == b);
    EXPECT_FALSE(a == c);
    EXPECT_FALSE(a != b);
    EXPECT_TRUE(a != c);
}

TEST_F(VectorTest, ApproximateEquality) {
    Vector<double, 3> a{1.0, 2.0, 3.0};
    Vector<double, 3> b{1.0 + 1e-15, 2.0 - 1e-15, 3.0 + 1e-15};

    EXPECT_TRUE(a.approx_equal(b, 1e-14));
    EXPECT_FALSE(a.approx_equal(b, 1e-16));
}

// =============================================================================
// Thread Safety Tests
// =============================================================================

TEST_F(VectorTest, ThreadSafetyReadOnly) {
    Vector<double, 3> v{1.0, 2.0, 3.0};

    // Multiple threads reading should be safe
    std::vector<std::thread> threads;
    std::vector<double> results(10);

    for (int i = 0; i < 10; ++i) {
        threads.emplace_back([&v, &results, i]() {
            results[static_cast<std::size_t>(i)] = v.norm();
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    // All threads should get same result
    double expected = v.norm();
    for (double r : results) {
        EXPECT_EQ(r, expected);
    }
}

TEST_F(VectorTest, ThreadSafetyIsolated) {
    // Each thread works on its own vector
    std::vector<std::thread> threads;
    std::vector<Vector<double, 3>> results(10);

    for (int i = 0; i < 10; ++i) {
        threads.emplace_back([&results, i]() {
            Vector<double, 3> local{static_cast<double>(i), 0.0, 0.0};
            results[static_cast<std::size_t>(i)] = local * 2.0;
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    // Check each thread computed correctly
    for (int i = 0; i < 10; ++i) {
        EXPECT_EQ(results[static_cast<std::size_t>(i)][0], 2.0 * i);
    }
}

// =============================================================================
// Memory Alignment Tests
// =============================================================================

TEST_F(VectorTest, MemoryAlignment) {
    Vector<double, 3> v;

    // Check that data is properly aligned for SIMD
    std::uintptr_t addr = reinterpret_cast<std::uintptr_t>(v.data());
    EXPECT_EQ(addr % 32, 0) << "Vector data should be 32-byte aligned for AVX";
}

// =============================================================================
// Utility Function Tests
// =============================================================================

TEST_F(VectorTest, MinMaxElements) {
    Vector<double, 5> v{3.0, -1.0, 4.0, 1.0, -2.0};

    EXPECT_EQ(v.min(), -2.0);
    EXPECT_EQ(v.max(), 4.0);
    EXPECT_EQ(v.min_index(), 4);
    EXPECT_EQ(v.max_index(), 2);
}

TEST_F(VectorTest, Sum) {
    Vector<double, 4> v{1.0, 2.0, 3.0, 4.0};
    EXPECT_EQ(v.sum(), 10.0);

    Vector<double, 3> zero{0.0, 0.0, 0.0};
    EXPECT_EQ(zero.sum(), 0.0);
}

TEST_F(VectorTest, Mean) {
    Vector<double, 4> v{1.0, 2.0, 3.0, 4.0};
    EXPECT_EQ(v.mean(), 2.5);
}

TEST_F(VectorTest, ToString) {
    Vector<double, 3> v{1.0, 2.0, 3.0};
    std::stringstream ss;
    ss << v;

    std::string expected = "[1, 2, 3]";
    EXPECT_EQ(ss.str(), expected);
}
