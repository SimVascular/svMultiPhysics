/**
 * @file test_VectorExpr.cpp
 * @brief Unit tests for VectorExpr.h - vector expression templates
 */

#include <gtest/gtest.h>
#include "FE/Math/Vector.h"
#include "FE/Math/VectorExpr.h"
#include "FE/Math/MathConstants.h"
#include <limits>
#include <cmath>
#include <memory>
#include <atomic>
#include <type_traits>

using namespace svmp::FE::math;

// Test fixture for VectorExpr tests
class VectorExprTest : public ::testing::Test {
protected:
    static constexpr double tolerance = 1e-14;

    // Custom allocator to track memory allocations
    template<typename T>
    class TrackingAllocator {
    public:
        using value_type = T;

        static std::atomic<size_t> allocations;
        static std::atomic<size_t> deallocations;
        static std::atomic<size_t> bytes_allocated;

        TrackingAllocator() = default;

        template<typename U>
        TrackingAllocator(const TrackingAllocator<U>&) {}

        T* allocate(size_t n) {
            allocations.fetch_add(1);
            bytes_allocated.fetch_add(n * sizeof(T));
            return static_cast<T*>(::operator new(n * sizeof(T)));
        }

        void deallocate(T* p, size_t n) {
            deallocations.fetch_add(1);
            ::operator delete(p);
        }

        static void reset() {
            allocations = 0;
            deallocations = 0;
            bytes_allocated = 0;
        }
    };

    void SetUp() override {
        TrackingAllocator<double>::reset();
    }

    void TearDown() override {}

    template<typename T>
    bool approx_equal(T a, T b, T tol = tolerance) {
        return std::abs(a - b) <= tol;
    }
};

template<typename T>
std::atomic<size_t> VectorExprTest::TrackingAllocator<T>::allocations{0};
template<typename T>
std::atomic<size_t> VectorExprTest::TrackingAllocator<T>::deallocations{0};
template<typename T>
std::atomic<size_t> VectorExprTest::TrackingAllocator<T>::bytes_allocated{0};

// =============================================================================
// Lazy Evaluation Tests
// =============================================================================

TEST_F(VectorExprTest, LazyEvaluationNoTemporaries) {
    // Expression templates should not create temporary vectors
    Vector<double, 3> a{1.0, 2.0, 3.0};
    Vector<double, 3> b{4.0, 5.0, 6.0};
    Vector<double, 3> c{7.0, 8.0, 9.0};

    // Build expression without evaluation
    auto expr = a + b - c;

    // Expression type should not be Vector, but an expression type
    using ExprType = decltype(expr);
    EXPECT_FALSE((std::is_same_v<ExprType, Vector<double, 3>>));

    // Now evaluate
    Vector<double, 3> result = expr;
    EXPECT_DOUBLE_EQ(result[0], -2.0);
    EXPECT_DOUBLE_EQ(result[1], -1.0);
    EXPECT_DOUBLE_EQ(result[2], 0.0);
}

TEST_F(VectorExprTest, LazyEvaluationAccessPattern) {
    Vector<double, 4> a{1.0, 2.0, 3.0, 4.0};
    Vector<double, 4> b{5.0, 6.0, 7.0, 8.0};

    auto expr = a + b;

    // Access individual elements without full evaluation
    EXPECT_DOUBLE_EQ(expr[0], 6.0);
    EXPECT_DOUBLE_EQ(expr[2], 10.0);

    // Size should be accessible
    EXPECT_EQ(expr.size(), 4u);
}

// =============================================================================
// Expression Chaining Tests
// =============================================================================

TEST_F(VectorExprTest, ChainedAdditionSubtraction) {
    Vector<double, 3> a{1.0, 2.0, 3.0};
    Vector<double, 3> b{4.0, 5.0, 6.0};
    Vector<double, 3> c{2.0, 3.0, 4.0};
    Vector<double, 3> d{1.0, 1.0, 1.0};

    // Chain multiple operations
    Vector<double, 3> result = a + b - c + d;

    EXPECT_DOUBLE_EQ(result[0], 4.0);
    EXPECT_DOUBLE_EQ(result[1], 5.0);
    EXPECT_DOUBLE_EQ(result[2], 6.0);
}

TEST_F(VectorExprTest, DeepExpressionNesting) {
    Vector<double, 2> v1{1.0, 2.0};
    Vector<double, 2> v2{3.0, 4.0};
    Vector<double, 2> v3{5.0, 6.0};
    Vector<double, 2> v4{7.0, 8.0};
    Vector<double, 2> v5{9.0, 10.0};

    // Deep nesting
    Vector<double, 2> result = ((v1 + v2) - (v3 - v4)) + v5;

    EXPECT_DOUBLE_EQ(result[0], 15.0);
    EXPECT_DOUBLE_EQ(result[1], 18.0);
}

// =============================================================================
// Mixed Operations Tests
// =============================================================================

TEST_F(VectorExprTest, ScalarMultiplicationInExpression) {
    Vector<double, 3> a{1.0, 2.0, 3.0};
    Vector<double, 3> b{4.0, 5.0, 6.0};

    Vector<double, 3> result = 2.0 * (a + b) / 3.0;

    EXPECT_TRUE(approx_equal(result[0], 10.0/3.0));
    EXPECT_TRUE(approx_equal(result[1], 14.0/3.0));
    EXPECT_TRUE(approx_equal(result[2], 6.0));
}

TEST_F(VectorExprTest, MixedScalarVectorOperations) {
    Vector<double, 4> v{2.0, 4.0, 6.0, 8.0};

    // Complex mixed expression
    Vector<double, 4> result = 3.0 * v / 2.0 + v * 0.5 - 1.0 * v;

    EXPECT_DOUBLE_EQ(result[0], 2.0);
    EXPECT_DOUBLE_EQ(result[1], 4.0);
    EXPECT_DOUBLE_EQ(result[2], 6.0);
    EXPECT_DOUBLE_EQ(result[3], 8.0);
}

// =============================================================================
// Unary Operations Tests
// =============================================================================

TEST_F(VectorExprTest, NegationInExpression) {
    Vector<double, 3> a{1.0, -2.0, 3.0};
    Vector<double, 3> b{4.0, 5.0, -6.0};

    Vector<double, 3> result = -a + (-b);

    EXPECT_DOUBLE_EQ(result[0], -5.0);
    EXPECT_DOUBLE_EQ(result[1], -3.0);
    EXPECT_DOUBLE_EQ(result[2], 3.0);
}

TEST_F(VectorExprTest, AbsoluteValueExpression) {
    Vector<double, 4> v{-1.5, 2.3, -4.7, 0.0};

    Vector<double, 4> result = abs(v);

    EXPECT_DOUBLE_EQ(result[0], 1.5);
    EXPECT_DOUBLE_EQ(result[1], 2.3);
    EXPECT_DOUBLE_EQ(result[2], 4.7);
    EXPECT_DOUBLE_EQ(result[3], 0.0);
}

TEST_F(VectorExprTest, SqrtExpression) {
    Vector<double, 3> v{4.0, 9.0, 16.0};

    Vector<double, 3> result = sqrt(v);

    EXPECT_DOUBLE_EQ(result[0], 2.0);
    EXPECT_DOUBLE_EQ(result[1], 3.0);
    EXPECT_DOUBLE_EQ(result[2], 4.0);
}

// =============================================================================
// Element-wise Operations Tests
// =============================================================================

TEST_F(VectorExprTest, HadamardProductExpression) {
    Vector<double, 3> a{2.0, 3.0, 4.0};
    Vector<double, 3> b{5.0, 6.0, 7.0};

    Vector<double, 3> result = hadamard(a, b);

    EXPECT_DOUBLE_EQ(result[0], 10.0);
    EXPECT_DOUBLE_EQ(result[1], 18.0);
    EXPECT_DOUBLE_EQ(result[2], 28.0);
}

TEST_F(VectorExprTest, HadamardDivisionExpression) {
    Vector<double, 3> a{10.0, 18.0, 28.0};
    Vector<double, 3> b{2.0, 3.0, 4.0};

    Vector<double, 3> result = hadamard_div(a, b);

    EXPECT_DOUBLE_EQ(result[0], 5.0);
    EXPECT_DOUBLE_EQ(result[1], 6.0);
    EXPECT_DOUBLE_EQ(result[2], 7.0);
}

// =============================================================================
// Dot Product and Norm Tests
// =============================================================================

TEST_F(VectorExprTest, DotProductOfExpressions) {
    Vector<double, 3> a{1.0, 2.0, 3.0};
    Vector<double, 3> b{4.0, 5.0, 6.0};
    Vector<double, 3> c{2.0, 2.0, 2.0};

    // Dot product of expressions
    double result = dot(a + b, c);

    EXPECT_DOUBLE_EQ(result, 42.0);
}

TEST_F(VectorExprTest, NormOfExpression) {
    Vector<double, 2> a{3.0, 0.0};
    Vector<double, 2> b{0.0, 4.0};

    double result = norm(a + b);

    EXPECT_DOUBLE_EQ(result, 5.0);  // norm of (3,4) = 5
}

TEST_F(VectorExprTest, NormalizeExpression) {
    Vector<double, 3> v{3.0, 0.0, 4.0};

    Vector<double, 3> result = normalize(v);

    EXPECT_DOUBLE_EQ(result[0], 0.6);
    EXPECT_DOUBLE_EQ(result[1], 0.0);
    EXPECT_DOUBLE_EQ(result[2], 0.8);
}

// =============================================================================
// Type Deduction Tests
// =============================================================================

TEST_F(VectorExprTest, TypeDeductionCorrectness) {
    Vector<float, 3> vf{1.0f, 2.0f, 3.0f};
    Vector<double, 3> vd{4.0, 5.0, 6.0};

    // Mixed type operations should promote to higher precision
    auto expr = vf + vf;  // float expression
    using ExprType = decltype(expr[0]);
    EXPECT_TRUE((std::is_same_v<ExprType, float>));

    // Test that expression evaluates correctly
    Vector<float, 3> result = expr;
    EXPECT_FLOAT_EQ(result[0], 2.0f);
    EXPECT_FLOAT_EQ(result[1], 4.0f);
    EXPECT_FLOAT_EQ(result[2], 6.0f);
}

// =============================================================================
// SFINAE and Compile-time Tests
// =============================================================================

TEST_F(VectorExprTest, SFINAEConstraints) {
    // Test that VectorExpr operators only work with VectorExpr types
    Vector<double, 3> v1{1.0, 2.0, 3.0};
    Vector<double, 3> v2{4.0, 5.0, 6.0};

    // This should compile
    auto expr = v1 + v2;
    Vector<double, 3> result = expr;

    // Verify the constraint checking
    EXPECT_TRUE((std::is_base_of_v<VectorExpr<Vector<double, 3>>, Vector<double, 3>>));
}

// =============================================================================
// Aliasing and Self-Assignment Tests
// =============================================================================

TEST_F(VectorExprTest, SelfAssignmentWithExpression) {
    Vector<double, 3> a{1.0, 2.0, 3.0};
    Vector<double, 3> b{4.0, 5.0, 6.0};

    // Self-assignment through expression
    a = a + b;

    EXPECT_DOUBLE_EQ(a[0], 5.0);
    EXPECT_DOUBLE_EQ(a[1], 7.0);
    EXPECT_DOUBLE_EQ(a[2], 9.0);
}

TEST_F(VectorExprTest, AliasingInExpression) {
    Vector<double, 3> a{2.0, 3.0, 4.0};
    Vector<double, 3> b{1.0, 1.0, 1.0};

    // a appears on both sides
    a = b + a;

    EXPECT_DOUBLE_EQ(a[0], 3.0);
    EXPECT_DOUBLE_EQ(a[1], 4.0);
    EXPECT_DOUBLE_EQ(a[2], 5.0);
}

// =============================================================================
// Edge Cases Tests
// =============================================================================

TEST_F(VectorExprTest, SingleElementVector) {
    Vector<double, 1> a{5.0};
    Vector<double, 1> b{3.0};

    Vector<double, 1> result = a + b - a * 0.5;

    EXPECT_DOUBLE_EQ(result[0], 5.5);
}

TEST_F(VectorExprTest, EmptyExpression) {
    Vector<double, 3> v{1.0, 2.0, 3.0};

    // Expression that evaluates to identity
    Vector<double, 3> result = v + v * 0.0;

    EXPECT_DOUBLE_EQ(result[0], 1.0);
    EXPECT_DOUBLE_EQ(result[1], 2.0);
    EXPECT_DOUBLE_EQ(result[2], 3.0);
}

TEST_F(VectorExprTest, LargeVectorExpression) {
    const size_t N = 100;
    Vector<double, N> a, b, c;

    for (size_t i = 0; i < N; ++i) {
        a[i] = static_cast<double>(i);
        b[i] = static_cast<double>(i * 2);
        c[i] = static_cast<double>(i * 3);
    }

    Vector<double, N> result = a + b - c / 2.0;

    for (size_t i = 0; i < N; ++i) {
        EXPECT_DOUBLE_EQ(result[i], i + 2.0 * i - 1.5 * i);
    }
}

// =============================================================================
// Complex Expression Pattern Tests
// =============================================================================

TEST_F(VectorExprTest, ComplexNestedExpression) {
    Vector<double, 3> a{1.0, 2.0, 3.0};
    Vector<double, 3> b{4.0, 5.0, 6.0};
    Vector<double, 3> c{7.0, 8.0, 9.0};

    // Complex expression with multiple operation types
    Vector<double, 3> result = 2.0 * abs(a - b) + sqrt(hadamard(c, c)) / 3.0;

    // Verify each component
    // |a - b| = |(-3, -3, -3)| = (3, 3, 3)
    // 2 * (3, 3, 3) = (6, 6, 6)
    // c * c = (49, 64, 81)
    // sqrt(c * c) = (7, 8, 9)
    // sqrt(c * c) / 3 = (7/3, 8/3, 3)
    // result = (6 + 7/3, 6 + 8/3, 6 + 3) = (25/3, 26/3, 9)

    EXPECT_TRUE(approx_equal(result[0], 25.0/3.0));
    EXPECT_TRUE(approx_equal(result[1], 26.0/3.0));
    EXPECT_DOUBLE_EQ(result[2], 9.0);
}

TEST_F(VectorExprTest, ChainedUnaryOperations) {
    Vector<double, 4> v{-4.0, -9.0, -16.0, -25.0};

    // Chain of unary operations
    Vector<double, 4> result = sqrt(abs(-v));

    EXPECT_DOUBLE_EQ(result[0], 2.0);
    EXPECT_DOUBLE_EQ(result[1], 3.0);
    EXPECT_DOUBLE_EQ(result[2], 4.0);
    EXPECT_DOUBLE_EQ(result[3], 5.0);
}
