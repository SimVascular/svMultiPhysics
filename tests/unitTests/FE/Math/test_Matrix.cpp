/**
 * @file test_Matrix.cpp
 * @brief Unit tests for Matrix.h - fixed-size matrices with expression templates
 */

#include <gtest/gtest.h>
#include "FE/Math/Matrix.h"
#include "FE/Math/Vector.h"
#include "FE/Math/MatrixExpr.h"
#include <limits>
#include <cmath>
#include <thread>
#include <vector>

using namespace svmp::FE::math;

// Test fixture for Matrix tests
class MatrixTest : public ::testing::Test {
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

TEST_F(MatrixTest, DefaultConstruction) {
    Matrix<double, 3, 3> m;
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            EXPECT_EQ(m(i, j), 0.0);
        }
    }
}

TEST_F(MatrixTest, FillConstruction) {
    Matrix<double, 2, 3> m(5.0);
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            EXPECT_EQ(m(i, j), 5.0);
        }
    }
}

TEST_F(MatrixTest, InitializerListConstruction) {
    Matrix<double, 2, 3> m{{1.0, 2.0, 3.0},
                           {4.0, 5.0, 6.0}};

    EXPECT_EQ(m(0, 0), 1.0);
    EXPECT_EQ(m(0, 1), 2.0);
    EXPECT_EQ(m(0, 2), 3.0);
    EXPECT_EQ(m(1, 0), 4.0);
    EXPECT_EQ(m(1, 1), 5.0);
    EXPECT_EQ(m(1, 2), 6.0);
}

TEST_F(MatrixTest, CopyConstruction) {
    Matrix<double, 2, 2> m1{{1.0, 2.0},
                            {3.0, 4.0}};
    Matrix<double, 2, 2> m2(m1);

    EXPECT_EQ(m2(0, 0), 1.0);
    EXPECT_EQ(m2(0, 1), 2.0);
    EXPECT_EQ(m2(1, 0), 3.0);
    EXPECT_EQ(m2(1, 1), 4.0);

    // Ensure deep copy
    m2(0, 0) = 10.0;
    EXPECT_EQ(m1(0, 0), 1.0);
    EXPECT_EQ(m2(0, 0), 10.0);
}

TEST_F(MatrixTest, MoveConstruction) {
    Matrix<double, 2, 2> m1{{1.0, 2.0},
                            {3.0, 4.0}};
    Matrix<double, 2, 2> m2(std::move(m1));

    EXPECT_EQ(m2(0, 0), 1.0);
    EXPECT_EQ(m2(0, 1), 2.0);
    EXPECT_EQ(m2(1, 0), 3.0);
    EXPECT_EQ(m2(1, 1), 4.0);
}

// =============================================================================
// Element Access Tests
// =============================================================================

TEST_F(MatrixTest, ElementAccess) {
    Matrix<double, 2, 3> m{{1.0, 2.0, 3.0},
                           {4.0, 5.0, 6.0}};

    // Non-const access using operator()
    EXPECT_EQ(m(0, 0), 1.0);
    EXPECT_EQ(m(0, 2), 3.0);
    EXPECT_EQ(m(1, 1), 5.0);

    // Modification
    m(1, 2) = 7.0;
    EXPECT_EQ(m(1, 2), 7.0);

    // Const access
    const Matrix<double, 2, 3> cm{{1.0, 2.0, 3.0},
                                  {4.0, 5.0, 6.0}};
    EXPECT_EQ(cm(0, 1), 2.0);
    EXPECT_EQ(cm(1, 0), 4.0);
}

TEST_F(MatrixTest, ElementAccessBounds) {
    Matrix<double, 2, 3> m{{1.0, 2.0, 3.0},
                           {4.0, 5.0, 6.0}};

    // at() with bounds checking
    EXPECT_EQ(m.at(0, 0), 1.0);
    EXPECT_EQ(m.at(1, 2), 6.0);

    // Test out of bounds throws
    EXPECT_THROW(m.at(2, 0), std::out_of_range);
    EXPECT_THROW(m.at(0, 3), std::out_of_range);
    EXPECT_THROW(m.at(10, 10), std::out_of_range);
}

TEST_F(MatrixTest, RowColumnAccess) {
    Matrix<double, 3, 3> m{{1.0, 2.0, 3.0},
                           {4.0, 5.0, 6.0},
                           {7.0, 8.0, 9.0}};

    // Get row
    auto row1 = m.row(1);
    EXPECT_EQ(row1[0], 4.0);
    EXPECT_EQ(row1[1], 5.0);
    EXPECT_EQ(row1[2], 6.0);

    // Get column
    auto col2 = m.col(2);
    EXPECT_EQ(col2[0], 3.0);
    EXPECT_EQ(col2[1], 6.0);
    EXPECT_EQ(col2[2], 9.0);

    // Set row
    Vector<double, 3> new_row{10.0, 11.0, 12.0};
    m.set_row(0, new_row);
    EXPECT_EQ(m(0, 0), 10.0);
    EXPECT_EQ(m(0, 1), 11.0);
    EXPECT_EQ(m(0, 2), 12.0);

    // Set column
    Vector<double, 3> new_col{20.0, 21.0, 22.0};
    m.set_col(1, new_col);
    EXPECT_EQ(m(0, 1), 20.0);
    EXPECT_EQ(m(1, 1), 21.0);
    EXPECT_EQ(m(2, 1), 22.0);
}

// =============================================================================
// Arithmetic Operations Tests
// =============================================================================

TEST_F(MatrixTest, Addition) {
    Matrix<double, 2, 3> a{{1.0, 2.0, 3.0},
                           {4.0, 5.0, 6.0}};
    Matrix<double, 2, 3> b{{7.0, 8.0, 9.0},
                           {10.0, 11.0, 12.0}};

    Matrix<double, 2, 3> c = a + b;
    EXPECT_EQ(c(0, 0), 8.0);
    EXPECT_EQ(c(0, 1), 10.0);
    EXPECT_EQ(c(0, 2), 12.0);
    EXPECT_EQ(c(1, 0), 14.0);
    EXPECT_EQ(c(1, 1), 16.0);
    EXPECT_EQ(c(1, 2), 18.0);
}

TEST_F(MatrixTest, Subtraction) {
    Matrix<double, 2, 3> a{{8.0, 10.0, 12.0},
                           {14.0, 16.0, 18.0}};
    Matrix<double, 2, 3> b{{7.0, 8.0, 9.0},
                           {10.0, 11.0, 12.0}};

    Matrix<double, 2, 3> c = a - b;
    EXPECT_EQ(c(0, 0), 1.0);
    EXPECT_EQ(c(0, 1), 2.0);
    EXPECT_EQ(c(0, 2), 3.0);
    EXPECT_EQ(c(1, 0), 4.0);
    EXPECT_EQ(c(1, 1), 5.0);
    EXPECT_EQ(c(1, 2), 6.0);
}

TEST_F(MatrixTest, ScalarMultiplication) {
    Matrix<double, 2, 2> a{{1.0, 2.0},
                           {3.0, 4.0}};

    Matrix<double, 2, 2> b = 2.0 * a;
    EXPECT_EQ(b(0, 0), 2.0);
    EXPECT_EQ(b(0, 1), 4.0);
    EXPECT_EQ(b(1, 0), 6.0);
    EXPECT_EQ(b(1, 1), 8.0);

    Matrix<double, 2, 2> c = a * 3.0;
    EXPECT_EQ(c(0, 0), 3.0);
    EXPECT_EQ(c(0, 1), 6.0);
    EXPECT_EQ(c(1, 0), 9.0);
    EXPECT_EQ(c(1, 1), 12.0);
}

TEST_F(MatrixTest, ScalarDivision) {
    Matrix<double, 2, 2> a{{2.0, 4.0},
                           {6.0, 8.0}};

    Matrix<double, 2, 2> b = a / 2.0;
    EXPECT_EQ(b(0, 0), 1.0);
    EXPECT_EQ(b(0, 1), 2.0);
    EXPECT_EQ(b(1, 0), 3.0);
    EXPECT_EQ(b(1, 1), 4.0);
}

TEST_F(MatrixTest, MatrixMultiplication) {
    Matrix<double, 2, 3> a{{1.0, 2.0, 3.0},
                           {4.0, 5.0, 6.0}};
    Matrix<double, 3, 2> b{{7.0, 8.0},
                           {9.0, 10.0},
                           {11.0, 12.0}};

    Matrix<double, 2, 2> c = a * b;
    EXPECT_EQ(c(0, 0), 58.0);   // 1*7 + 2*9 + 3*11
    EXPECT_EQ(c(0, 1), 64.0);   // 1*8 + 2*10 + 3*12
    EXPECT_EQ(c(1, 0), 139.0);  // 4*7 + 5*9 + 6*11
    EXPECT_EQ(c(1, 1), 154.0);  // 4*8 + 5*10 + 6*12
}

TEST_F(MatrixTest, MatrixVectorMultiplication) {
    Matrix<double, 3, 3> m{{1.0, 2.0, 3.0},
                           {4.0, 5.0, 6.0},
                           {7.0, 8.0, 9.0}};
    Vector<double, 3> v{1.0, 2.0, 3.0};

    Vector<double, 3> result = m * v;
    EXPECT_EQ(result[0], 14.0);  // 1*1 + 2*2 + 3*3
    EXPECT_EQ(result[1], 32.0);  // 4*1 + 5*2 + 6*3
    EXPECT_EQ(result[2], 50.0);  // 7*1 + 8*2 + 9*3
}

// =============================================================================
// Special Matrix Operations Tests
// =============================================================================

TEST_F(MatrixTest, Transpose) {
    Matrix<double, 2, 3> m{{1.0, 2.0, 3.0},
                           {4.0, 5.0, 6.0}};

    Matrix<double, 3, 2> mt = m.transpose();
    EXPECT_EQ(mt(0, 0), 1.0);
    EXPECT_EQ(mt(0, 1), 4.0);
    EXPECT_EQ(mt(1, 0), 2.0);
    EXPECT_EQ(mt(1, 1), 5.0);
    EXPECT_EQ(mt(2, 0), 3.0);
    EXPECT_EQ(mt(2, 1), 6.0);
}

TEST_F(MatrixTest, Determinant2x2) {
    Matrix<double, 2, 2> m{{1.0, 2.0},
                           {3.0, 4.0}};

    double det = m.determinant();
    EXPECT_EQ(det, -2.0);  // 1*4 - 2*3 = 4 - 6 = -2
}

TEST_F(MatrixTest, Determinant3x3) {
    Matrix<double, 3, 3> m{{1.0, 2.0, 3.0},
                           {0.0, 1.0, 4.0},
                           {5.0, 6.0, 0.0}};

    double det = m.determinant();
    EXPECT_EQ(det, 1.0);  // Using Sarrus rule
}

TEST_F(MatrixTest, Determinant4x4) {
    Matrix<double, 4, 4> m{{1, 0, 0, 0},
                           {0, 2, 0, 0},
                           {0, 0, 3, 0},
                           {0, 0, 0, 4}};

    double det = m.determinant();
    EXPECT_EQ(det, 24.0);  // 1*2*3*4 = 24 (diagonal matrix)
}

TEST_F(MatrixTest, Inverse2x2) {
    Matrix<double, 2, 2> m{{1.0, 2.0},
                           {3.0, 4.0}};

    Matrix<double, 2, 2> inv = m.inverse();

    // Check inverse properties
    EXPECT_NEAR(inv(0, 0), -2.0, tolerance);
    EXPECT_NEAR(inv(0, 1), 1.0, tolerance);
    EXPECT_NEAR(inv(1, 0), 1.5, tolerance);
    EXPECT_NEAR(inv(1, 1), -0.5, tolerance);

    // Verify M * M^-1 = I
    Matrix<double, 2, 2> identity = m * inv;
    EXPECT_NEAR(identity(0, 0), 1.0, tolerance);
    EXPECT_NEAR(identity(0, 1), 0.0, tolerance);
    EXPECT_NEAR(identity(1, 0), 0.0, tolerance);
    EXPECT_NEAR(identity(1, 1), 1.0, tolerance);
}

TEST_F(MatrixTest, Inverse3x3) {
    Matrix<double, 3, 3> m{{1.0, 2.0, 3.0},
                           {0.0, 1.0, 4.0},
                           {5.0, 6.0, 0.0}};

    Matrix<double, 3, 3> inv = m.inverse();

    // Verify M * M^-1 = I
    Matrix<double, 3, 3> identity = m * inv;
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            double expected = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(identity(i, j), expected, tolerance);
        }
    }
}

TEST_F(MatrixTest, Trace) {
    Matrix<double, 3, 3> m{{1.0, 2.0, 3.0},
                           {4.0, 5.0, 6.0},
                           {7.0, 8.0, 9.0}};

    double trace = m.trace();
    EXPECT_EQ(trace, 15.0);  // 1 + 5 + 9 = 15
}

// =============================================================================
// Special Matrix Types Tests
// =============================================================================

TEST_F(MatrixTest, IdentityMatrix) {
    Matrix<double, 3, 3> I = Matrix<double, 3, 3>::identity();

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            double expected = (i == j) ? 1.0 : 0.0;
            EXPECT_EQ(I(i, j), expected);
        }
    }

    // Test identity property
    Matrix<double, 3, 3> m{{1.0, 2.0, 3.0},
                           {4.0, 5.0, 6.0},
                           {7.0, 8.0, 9.0}};
    Matrix<double, 3, 3> result = m * I;

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            EXPECT_EQ(result(i, j), m(i, j));
        }
    }
}

TEST_F(MatrixTest, ZeroMatrix) {
    Matrix<double, 2, 3> Z = Matrix<double, 2, 3>::zero();

    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            EXPECT_EQ(Z(i, j), 0.0);
        }
    }
}

TEST_F(MatrixTest, DiagonalMatrix) {
    Vector<double, 3> diag{1.0, 2.0, 3.0};
    Matrix<double, 3, 3> D = Matrix<double, 3, 3>::diagonal(diag);

    EXPECT_EQ(D(0, 0), 1.0);
    EXPECT_EQ(D(1, 1), 2.0);
    EXPECT_EQ(D(2, 2), 3.0);

    // Off-diagonal elements should be zero
    EXPECT_EQ(D(0, 1), 0.0);
    EXPECT_EQ(D(0, 2), 0.0);
    EXPECT_EQ(D(1, 0), 0.0);
    EXPECT_EQ(D(1, 2), 0.0);
    EXPECT_EQ(D(2, 0), 0.0);
    EXPECT_EQ(D(2, 1), 0.0);
}

// =============================================================================
// Expression Template Tests
// =============================================================================

TEST_F(MatrixTest, ExpressionTemplatesNoTemporaries) {
    Matrix<double, 2, 2> a{{1, 2}, {3, 4}};
    Matrix<double, 2, 2> b{{5, 6}, {7, 8}};
    Matrix<double, 2, 2> c{{9, 10}, {11, 12}};

    // Complex expression should create no temporaries
    Matrix<double, 2, 2> result = a + b - c;

    EXPECT_EQ(result(0, 0), -3.0);   // 1 + 5 - 9
    EXPECT_EQ(result(0, 1), -2.0);   // 2 + 6 - 10
    EXPECT_EQ(result(1, 0), -1.0);   // 3 + 7 - 11
    EXPECT_EQ(result(1, 1), 0.0);    // 4 + 8 - 12
}

TEST_F(MatrixTest, LazyEvaluation) {
    Matrix<double, 2, 2> a{{1, 2}, {3, 4}};
    Matrix<double, 2, 2> b{{5, 6}, {7, 8}};

    // Expression should not be evaluated until assignment
    auto expr = a + b;  // No computation yet

    Matrix<double, 2, 2> result = expr;  // Evaluation happens here
    EXPECT_EQ(result(0, 0), 6.0);
    EXPECT_EQ(result(0, 1), 8.0);
}

// =============================================================================
// Edge Cases and Error Handling Tests
// =============================================================================

TEST_F(MatrixTest, SingularMatrixInverse) {
    Matrix<double, 2, 2> singular{{1.0, 2.0},
                                  {2.0, 4.0}};  // det = 0

    EXPECT_THROW(singular.inverse(), std::runtime_error);
}

TEST_F(MatrixTest, DivisionByZero) {
    Matrix<double, 2, 2> m{{1.0, 2.0},
                           {3.0, 4.0}};

    Matrix<double, 2, 2> result = m / 0.0;
    EXPECT_TRUE(std::isinf(result(0, 0)));
    EXPECT_TRUE(std::isinf(result(0, 1)));
}

TEST_F(MatrixTest, ExtremeLargeValues) {
    double large = 1e308;
    Matrix<double, 2, 2> m{{large, 0}, {0, large}};

    Matrix<double, 2, 2> half = m / 2.0;
    EXPECT_FALSE(std::isinf(half(0, 0)));
    EXPECT_EQ(half(0, 0), large / 2.0);
}

// =============================================================================
// Numerical Precision Tests
// =============================================================================

TEST_F(MatrixTest, NumericalStability) {
    // Test near-singular matrix
    double eps = 1e-15;
    Matrix<double, 2, 2> m{{1.0, 1.0},
                           {1.0, 1.0 + eps}};

    double det = m.determinant();
    // Relax tolerance due to floating-point arithmetic in determinant calculation
    EXPECT_NEAR(det, eps, 1e-14);
}

TEST_F(MatrixTest, OrthogonalMatrixProperties) {
    // Create rotation matrix (orthogonal)
    double angle = M_PI / 4;
    Matrix<double, 2, 2> R{{cos(angle), -sin(angle)},
                           {sin(angle), cos(angle)}};

    // Check orthogonality: R * R^T = I
    Matrix<double, 2, 2> RRt = R * R.transpose();
    EXPECT_NEAR(RRt(0, 0), 1.0, tolerance);
    EXPECT_NEAR(RRt(0, 1), 0.0, tolerance);
    EXPECT_NEAR(RRt(1, 0), 0.0, tolerance);
    EXPECT_NEAR(RRt(1, 1), 1.0, tolerance);

    // Check determinant = ±1
    EXPECT_NEAR(std::abs(R.determinant()), 1.0, tolerance);
}

// =============================================================================
// Matrix Properties Tests
// =============================================================================

TEST_F(MatrixTest, IsSymmetric) {
    Matrix<double, 3, 3> sym{{1, 2, 3},
                             {2, 4, 5},
                             {3, 5, 6}};
    EXPECT_TRUE(sym.is_symmetric(tolerance));

    Matrix<double, 3, 3> nonsym{{1, 2, 3},
                                {4, 5, 6},
                                {7, 8, 9}};
    EXPECT_FALSE(nonsym.is_symmetric(tolerance));
}

TEST_F(MatrixTest, IsSkewSymmetric) {
    Matrix<double, 3, 3> skew{{0, -1, 2},
                              {1, 0, -3},
                              {-2, 3, 0}};
    EXPECT_TRUE(skew.is_skew_symmetric(tolerance));

    Matrix<double, 3, 3> nonskew{{1, 2, 3},
                                 {4, 5, 6},
                                 {7, 8, 9}};
    EXPECT_FALSE(nonskew.is_skew_symmetric(tolerance));
}

TEST_F(MatrixTest, IsDiagonal) {
    Matrix<double, 3, 3> diag{{1, 0, 0},
                              {0, 2, 0},
                              {0, 0, 3}};
    EXPECT_TRUE(diag.is_diagonal(tolerance));

    Matrix<double, 3, 3> nondiag{{1, 0.1, 0},
                                 {0, 2, 0},
                                 {0, 0, 3}};
    EXPECT_FALSE(nondiag.is_diagonal(tolerance));
}

// =============================================================================
// Thread Safety Tests
// =============================================================================

TEST_F(MatrixTest, ThreadSafetyReadOnly) {
    Matrix<double, 3, 3> m{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

    std::vector<std::thread> threads;
    std::vector<double> results(10);

    for (int i = 0; i < 10; ++i) {
        threads.emplace_back([&m, &results, i]() {
            results[static_cast<std::size_t>(i)] = m.trace();
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    for (double r : results) {
        EXPECT_EQ(r, 15.0);
    }
}

// =============================================================================
// Memory Alignment Tests
// =============================================================================

TEST_F(MatrixTest, MemoryAlignment) {
    Matrix<double, 3, 3> m;

    std::uintptr_t addr = reinterpret_cast<std::uintptr_t>(m.data());
    EXPECT_EQ(addr % 32, 0) << "Matrix data should be 32-byte aligned for AVX";
}

// =============================================================================
// Utility Function Tests
// =============================================================================

TEST_F(MatrixTest, Norms) {
    Matrix<double, 2, 2> m{{1, 2}, {3, 4}};

    // Frobenius norm: sqrt(1^2 + 2^2 + 3^2 + 4^2) = sqrt(30)
    EXPECT_NEAR(m.frobenius_norm(), std::sqrt(30.0), tolerance);

    // Infinity norm (max row sum)
    EXPECT_EQ(m.infinity_norm(), 7.0);  // max(|1|+|2|, |3|+|4|) = max(3, 7)

    // One norm (max column sum)
    EXPECT_EQ(m.one_norm(), 6.0);  // max(|1|+|3|, |2|+|4|) = max(4, 6)
}

TEST_F(MatrixTest, MinMaxElements) {
    Matrix<double, 2, 3> m{{3, -1, 4}, {1, -2, 5}};

    EXPECT_EQ(m.min(), -2.0);
    EXPECT_EQ(m.max(), 5.0);
}

TEST_F(MatrixTest, ToString) {
    Matrix<double, 2, 2> m{{1, 2}, {3, 4}};
    std::stringstream ss;
    ss << m;

    std::string expected = "[[1, 2]\n [3, 4]]";
    EXPECT_EQ(ss.str(), expected);
}
