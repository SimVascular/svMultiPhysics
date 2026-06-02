/**
 * @file test_MatrixExpr.cpp
 * @brief Unit tests for MatrixExpr.h - matrix expression templates
 */

#include <gtest/gtest.h>
#include "FE/Math/Matrix.h"
#include "FE/Math/MatrixExpr.h"
#include "FE/Math/Vector.h"
#include "FE/Math/MathConstants.h"
#include <limits>
#include <cmath>
#include <memory>
#include <atomic>
#include <type_traits>

using namespace svmp::FE::math;

// Test fixture for MatrixExpr tests
class MatrixExprTest : public ::testing::Test {
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
std::atomic<size_t> MatrixExprTest::TrackingAllocator<T>::allocations{0};
template<typename T>
std::atomic<size_t> MatrixExprTest::TrackingAllocator<T>::deallocations{0};
template<typename T>
std::atomic<size_t> MatrixExprTest::TrackingAllocator<T>::bytes_allocated{0};

// =============================================================================
// Lazy Evaluation Tests
// =============================================================================

TEST_F(MatrixExprTest, LazyEvaluationNoTemporaries) {
    // Expression templates should not create temporary matrices
    Matrix<double, 2, 2> A{{1.0, 2.0}, {3.0, 4.0}};
    Matrix<double, 2, 2> B{{5.0, 6.0}, {7.0, 8.0}};
    Matrix<double, 2, 2> C{{9.0, 10.0}, {11.0, 12.0}};

    // Build expression without evaluation
    auto expr = A + B - C;

    // Expression type should not be Matrix, but an expression type
    using ExprType = decltype(expr);
    EXPECT_FALSE((std::is_same_v<ExprType, Matrix<double, 2, 2>>));

    // Now evaluate
    Matrix<double, 2, 2> result = expr;
    EXPECT_DOUBLE_EQ(result(0, 0), -3.0);
    EXPECT_DOUBLE_EQ(result(0, 1), -2.0);
    EXPECT_DOUBLE_EQ(result(1, 0), -1.0);
    EXPECT_DOUBLE_EQ(result(1, 1), 0.0);
}

TEST_F(MatrixExprTest, LazyEvaluationAccessPattern) {
    Matrix<double, 3, 3> A;
    Matrix<double, 3, 3> B;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            A(i, j) = i * 3 + j + 1;
            B(i, j) = (i * 3 + j + 1) * 2;
        }
    }

    auto expr = A + B;

    // Access individual elements without full evaluation
    EXPECT_DOUBLE_EQ(expr(0, 0), 3.0);
    EXPECT_DOUBLE_EQ(expr(1, 1), 15.0);
    EXPECT_DOUBLE_EQ(expr(2, 2), 27.0);

    // Size should be accessible
    EXPECT_EQ(expr.rows(), 3u);
    EXPECT_EQ(expr.cols(), 3u);
}

// =============================================================================
// Matrix Multiplication Tests
// =============================================================================

TEST_F(MatrixExprTest, MatrixMultiplicationExpression) {
    Matrix<double, 2, 3> A{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
    Matrix<double, 3, 2> B{{7.0, 8.0}, {9.0, 10.0}, {11.0, 12.0}};

    // Matrix multiplication should produce 2x2 result
    Matrix<double, 2, 2> C = A * B;

    // Verify results
    EXPECT_DOUBLE_EQ(C(0, 0), 58.0);   // 1*7 + 2*9 + 3*11
    EXPECT_DOUBLE_EQ(C(0, 1), 64.0);   // 1*8 + 2*10 + 3*12
    EXPECT_DOUBLE_EQ(C(1, 0), 139.0);  // 4*7 + 5*9 + 6*11
    EXPECT_DOUBLE_EQ(C(1, 1), 154.0);  // 4*8 + 5*10 + 6*12
}

TEST_F(MatrixExprTest, ChainedMatrixMultiplication) {
    Matrix<double, 2, 2> A{{1.0, 2.0}, {3.0, 4.0}};
    Matrix<double, 2, 2> B{{5.0, 6.0}, {7.0, 8.0}};
    Matrix<double, 2, 2> C{{9.0, 10.0}, {11.0, 12.0}};

    // Chain matrix multiplications: (A * B) * C
    Matrix<double, 2, 2> result = A * B * C;

    // First compute A * B
    Matrix<double, 2, 2> AB = A * B;
    EXPECT_DOUBLE_EQ(AB(0, 0), 19.0);  // 1*5 + 2*7
    EXPECT_DOUBLE_EQ(AB(0, 1), 22.0);  // 1*6 + 2*8
    EXPECT_DOUBLE_EQ(AB(1, 0), 43.0);  // 3*5 + 4*7
    EXPECT_DOUBLE_EQ(AB(1, 1), 50.0);  // 3*6 + 4*8

    // Then (A * B) * C
    EXPECT_DOUBLE_EQ(result(0, 0), 413.0);  // 19*9 + 22*11
    EXPECT_DOUBLE_EQ(result(0, 1), 454.0);  // 19*10 + 22*12
    EXPECT_DOUBLE_EQ(result(1, 0), 937.0);  // 43*9 + 50*11
    EXPECT_DOUBLE_EQ(result(1, 1), 1030.0); // 43*10 + 50*12
}

// =============================================================================
// Mixed Operations Tests
// =============================================================================

TEST_F(MatrixExprTest, MixedMatrixOperations) {
    Matrix<double, 3, 3> A, B, C, D;

    // Initialize matrices
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            A(i, j) = i + j + 1;
            B(i, j) = (i + 1) * (j + 1);
            C(i, j) = i * j + 1;
            D(i, j) = 1.0;
        }
    }

    // Complex expression: A * B + C * D
    Matrix<double, 3, 3> result = A * B + C * D;

    // Verify a few key elements
    Matrix<double, 3, 3> AB = A * B;
    Matrix<double, 3, 3> CD = C * D;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_DOUBLE_EQ(result(i, j), AB(i, j) + CD(i, j));
        }
    }
}

TEST_F(MatrixExprTest, ScalarMultiplicationInExpression) {
    Matrix<double, 2, 2> A{{1.0, 2.0}, {3.0, 4.0}};
    Matrix<double, 2, 2> B{{5.0, 6.0}, {7.0, 8.0}};

    Matrix<double, 2, 2> result = 2.0 * (A + B) / 3.0;

    EXPECT_TRUE(approx_equal(result(0, 0), 4.0));
    EXPECT_TRUE(approx_equal(result(0, 1), 16.0/3.0));
    EXPECT_TRUE(approx_equal(result(1, 0), 20.0/3.0));
    EXPECT_TRUE(approx_equal(result(1, 1), 8.0));
}

// =============================================================================
// Transpose Tests
// =============================================================================

TEST_F(MatrixExprTest, TransposeExpression) {
    Matrix<double, 2, 3> A{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};

    auto AT = transpose(A);

    // Check dimensions
    EXPECT_EQ(AT.rows(), 3u);
    EXPECT_EQ(AT.cols(), 2u);

    // Check values
    EXPECT_DOUBLE_EQ(AT(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(AT(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(AT(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(AT(1, 1), 5.0);
    EXPECT_DOUBLE_EQ(AT(2, 0), 3.0);
    EXPECT_DOUBLE_EQ(AT(2, 1), 6.0);
}

TEST_F(MatrixExprTest, TransposeInExpression) {
    Matrix<double, 3, 2> A{{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
    Matrix<double, 3, 2> B{{7.0, 8.0}, {9.0, 10.0}, {11.0, 12.0}};

    // Compute A^T * B (should be 2x2)
    Matrix<double, 2, 2> result = transpose(A) * B;

    EXPECT_DOUBLE_EQ(result(0, 0), 89.0);   // 1*7 + 3*9 + 5*11
    EXPECT_DOUBLE_EQ(result(0, 1), 98.0);   // 1*8 + 3*10 + 5*12
    EXPECT_DOUBLE_EQ(result(1, 0), 116.0);  // 2*7 + 4*9 + 6*11
    EXPECT_DOUBLE_EQ(result(1, 1), 128.0);  // 2*8 + 4*10 + 6*12
}

// =============================================================================
// Unary Operations Tests
// =============================================================================

TEST_F(MatrixExprTest, NegationInExpression) {
    Matrix<double, 2, 2> A{{1.0, -2.0}, {3.0, -4.0}};
    Matrix<double, 2, 2> B{{5.0, 6.0}, {-7.0, 8.0}};

    Matrix<double, 2, 2> result = -A + (-B);

    EXPECT_DOUBLE_EQ(result(0, 0), -6.0);
    EXPECT_DOUBLE_EQ(result(0, 1), -4.0);
    EXPECT_DOUBLE_EQ(result(1, 0), 4.0);
    EXPECT_DOUBLE_EQ(result(1, 1), -4.0);
}

TEST_F(MatrixExprTest, AbsoluteValueExpression) {
    Matrix<double, 2, 3> M{{-1.5, 2.3, -4.7}, {0.0, -3.2, 5.1}};

    Matrix<double, 2, 3> result = abs(M);

    EXPECT_DOUBLE_EQ(result(0, 0), 1.5);
    EXPECT_DOUBLE_EQ(result(0, 1), 2.3);
    EXPECT_DOUBLE_EQ(result(0, 2), 4.7);
    EXPECT_DOUBLE_EQ(result(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(result(1, 1), 3.2);
    EXPECT_DOUBLE_EQ(result(1, 2), 5.1);
}

TEST_F(MatrixExprTest, SqrtExpression) {
    Matrix<double, 2, 2> M{{4.0, 9.0}, {16.0, 25.0}};

    Matrix<double, 2, 2> result = sqrt(M);

    EXPECT_DOUBLE_EQ(result(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(result(0, 1), 3.0);
    EXPECT_DOUBLE_EQ(result(1, 0), 4.0);
    EXPECT_DOUBLE_EQ(result(1, 1), 5.0);
}

// =============================================================================
// Element-wise Operations Tests
// =============================================================================

TEST_F(MatrixExprTest, HadamardProductExpression) {
    Matrix<double, 2, 3> A{{2.0, 3.0, 4.0}, {5.0, 6.0, 7.0}};
    Matrix<double, 2, 3> B{{8.0, 9.0, 10.0}, {11.0, 12.0, 13.0}};

    Matrix<double, 2, 3> result = hadamard(A, B);

    EXPECT_DOUBLE_EQ(result(0, 0), 16.0);
    EXPECT_DOUBLE_EQ(result(0, 1), 27.0);
    EXPECT_DOUBLE_EQ(result(0, 2), 40.0);
    EXPECT_DOUBLE_EQ(result(1, 0), 55.0);
    EXPECT_DOUBLE_EQ(result(1, 1), 72.0);
    EXPECT_DOUBLE_EQ(result(1, 2), 91.0);
}

TEST_F(MatrixExprTest, HadamardDivisionExpression) {
    Matrix<double, 2, 2> A{{10.0, 18.0}, {28.0, 36.0}};
    Matrix<double, 2, 2> B{{2.0, 3.0}, {4.0, 6.0}};

    Matrix<double, 2, 2> result = hadamard_div(A, B);

    EXPECT_DOUBLE_EQ(result(0, 0), 5.0);
    EXPECT_DOUBLE_EQ(result(0, 1), 6.0);
    EXPECT_DOUBLE_EQ(result(1, 0), 7.0);
    EXPECT_DOUBLE_EQ(result(1, 1), 6.0);
}

// =============================================================================
// Norm and Trace Tests
// =============================================================================

TEST_F(MatrixExprTest, FrobeniusNormOfExpression) {
    Matrix<double, 2, 2> A{{1.0, 2.0}, {3.0, 4.0}};
    Matrix<double, 2, 2> B{{2.0, 2.0}, {2.0, 2.0}};

    double norm_sq = frobenius_norm_squared(A - B);
    double norm = frobenius_norm(A - B);

    // (A - B) = [[-1, 0], [1, 2]]
    // norm_squared = 1 + 0 + 1 + 4 = 6
    EXPECT_DOUBLE_EQ(norm_sq, 6.0);
    EXPECT_DOUBLE_EQ(norm, std::sqrt(6.0));
}

TEST_F(MatrixExprTest, TraceOfExpression) {
    Matrix<double, 3, 3> A;
    Matrix<double, 3, 3> B;

    // Initialize as diagonal matrices
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            A(i, j) = (i == j) ? (i + 1) : 0.0;  // diag(1, 2, 3)
            B(i, j) = (i == j) ? (i + 4) : 0.0;  // diag(4, 5, 6)
        }
    }

    double tr = trace(A + B);

    // trace(A + B) = trace(diag(5, 7, 9)) = 21
    EXPECT_DOUBLE_EQ(tr, 21.0);
}

// =============================================================================
// Type Deduction Tests
// =============================================================================

TEST_F(MatrixExprTest, TypeDeductionCorrectness) {
    Matrix<float, 2, 2> Mf{{1.0f, 2.0f}, {3.0f, 4.0f}};
    Matrix<double, 2, 2> Md{{5.0, 6.0}, {7.0, 8.0}};

    // Float expression
    auto expr = Mf + Mf;
    using ExprType = decltype(expr(0, 0));
    EXPECT_TRUE((std::is_same_v<ExprType, float>));

    // Test that expression evaluates correctly
    Matrix<float, 2, 2> result = expr;
    EXPECT_FLOAT_EQ(result(0, 0), 2.0f);
    EXPECT_FLOAT_EQ(result(1, 1), 8.0f);
}

// =============================================================================
// SFINAE and Compile-time Tests
// =============================================================================

TEST_F(MatrixExprTest, SFINAEConstraints) {
    // Test that MatrixExpr operators only work with MatrixExpr types
    Matrix<double, 2, 2> M1{{1.0, 2.0}, {3.0, 4.0}};
    Matrix<double, 2, 2> M2{{5.0, 6.0}, {7.0, 8.0}};

    // This should compile
    auto expr = M1 + M2;
    Matrix<double, 2, 2> result = expr;

    // Verify the constraint checking
    EXPECT_TRUE((std::is_base_of_v<MatrixExpr<Matrix<double, 2, 2>>, Matrix<double, 2, 2>>));
}

// =============================================================================
// Aliasing and Self-Assignment Tests
// =============================================================================

TEST_F(MatrixExprTest, SelfAssignmentWithExpression) {
    Matrix<double, 2, 2> A{{1.0, 2.0}, {3.0, 4.0}};
    Matrix<double, 2, 2> B{{5.0, 6.0}, {7.0, 8.0}};

    // Self-assignment through expression
    A = A + B;

    EXPECT_DOUBLE_EQ(A(0, 0), 6.0);
    EXPECT_DOUBLE_EQ(A(0, 1), 8.0);
    EXPECT_DOUBLE_EQ(A(1, 0), 10.0);
    EXPECT_DOUBLE_EQ(A(1, 1), 12.0);
}

TEST_F(MatrixExprTest, AliasingInExpression) {
    Matrix<double, 2, 2> A{{2.0, 3.0}, {4.0, 5.0}};
    Matrix<double, 2, 2> B{{1.0, 1.0}, {1.0, 1.0}};

    // A appears on both sides
    A = B + A;

    EXPECT_DOUBLE_EQ(A(0, 0), 3.0);
    EXPECT_DOUBLE_EQ(A(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(A(1, 0), 5.0);
    EXPECT_DOUBLE_EQ(A(1, 1), 6.0);
}

// =============================================================================
// Edge Cases Tests
// =============================================================================

TEST_F(MatrixExprTest, SingleElementMatrix) {
    Matrix<double, 1, 1> A{{5.0}};
    Matrix<double, 1, 1> B{{3.0}};

    Matrix<double, 1, 1> result = A + B - A * 0.5;

    EXPECT_DOUBLE_EQ(result(0, 0), 5.5);
}

TEST_F(MatrixExprTest, NonSquareMatrixOperations) {
    Matrix<double, 2, 4> A;
    Matrix<double, 2, 4> B;

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 4; ++j) {
            A(i, j) = i * 4 + j + 1;
            B(i, j) = (i * 4 + j + 1) * 2;
        }
    }

    Matrix<double, 2, 4> result = A + B - A * 0.5;

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 4; ++j) {
            double expected = A(i, j) + B(i, j) - A(i, j) * 0.5;
            EXPECT_DOUBLE_EQ(result(i, j), expected);
        }
    }
}

// =============================================================================
// Diagonal Matrix Tests
// =============================================================================

TEST_F(MatrixExprTest, DiagonalMatrixExpression) {
    Vector<double, 3> v{1.0, 2.0, 3.0};

    auto diag = DiagonalExpr<Vector<double, 3>>(v);

    // Check dimensions
    EXPECT_EQ(diag.rows(), 3u);
    EXPECT_EQ(diag.cols(), 3u);

    // Check values
    EXPECT_DOUBLE_EQ(diag(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(diag(1, 1), 2.0);
    EXPECT_DOUBLE_EQ(diag(2, 2), 3.0);

    // Off-diagonal should be zero
    EXPECT_DOUBLE_EQ(diag(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(diag(1, 0), 0.0);
}

TEST_F(MatrixExprTest, DiagonalMatrixInExpression) {
    Vector<double, 2> v{2.0, 3.0};
    Matrix<double, 2, 2> A{{1.0, 1.0}, {1.0, 1.0}};

    auto diag = DiagonalExpr<Vector<double, 2>>(v);
    Matrix<double, 2, 2> result = A + diag;

    EXPECT_DOUBLE_EQ(result(0, 0), 3.0);
    EXPECT_DOUBLE_EQ(result(0, 1), 1.0);
    EXPECT_DOUBLE_EQ(result(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(result(1, 1), 4.0);
}

// =============================================================================
// Complex Expression Pattern Tests
// =============================================================================

TEST_F(MatrixExprTest, ComplexNestedExpression) {
    Matrix<double, 2, 2> A{{1.0, 2.0}, {3.0, 4.0}};
    Matrix<double, 2, 2> B{{5.0, 6.0}, {7.0, 8.0}};
    Matrix<double, 2, 2> C{{9.0, 10.0}, {11.0, 12.0}};

    // Complex expression with multiple operation types
    Matrix<double, 2, 2> result = 2.0 * abs(A - B) + sqrt(hadamard(C, C)) / 3.0;

    // |A - B| = |[-4, -4], [-4, -4]| = [4, 4], [4, 4]
    // 2 * [4, 4], [4, 4] = [8, 8], [8, 8]
    // C * C (element-wise) = [81, 100], [121, 144]
    // sqrt(C * C) = [9, 10], [11, 12]
    // sqrt(C * C) / 3 = [3, 10/3], [11/3, 4]
    // result = [11, 34/3], [35/3, 12]

    EXPECT_DOUBLE_EQ(result(0, 0), 11.0);
    EXPECT_TRUE(approx_equal(result(0, 1), 34.0/3.0));
    EXPECT_TRUE(approx_equal(result(1, 0), 35.0/3.0));
    EXPECT_DOUBLE_EQ(result(1, 1), 12.0);
}

TEST_F(MatrixExprTest, MatrixVectorMixedExpression) {
    Matrix<double, 3, 3> A;
    Vector<double, 3> v{1.0, 2.0, 3.0};

    // Create identity matrix
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            A(i, j) = (i == j) ? 1.0 : 0.0;
        }
    }

    // Create diagonal from vector and add to identity
    auto diag = DiagonalExpr<Vector<double, 3>>(v);
    Matrix<double, 3, 3> result = A + diag;

    // Result should be diag(2, 3, 4)
    EXPECT_DOUBLE_EQ(result(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(result(1, 1), 3.0);
    EXPECT_DOUBLE_EQ(result(2, 2), 4.0);
    EXPECT_DOUBLE_EQ(result(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(result(1, 0), 0.0);
}
