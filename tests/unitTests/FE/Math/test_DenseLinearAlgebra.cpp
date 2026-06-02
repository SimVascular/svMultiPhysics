/**
 * @file test_DenseLinearAlgebra.cpp
 * @brief Tests for shared dense linear algebra utilities.
 */

#include <gtest/gtest.h>

#include "FE/Common/FEException.h"
#include "FE/Math/DenseLinearAlgebra.h"

#include <cmath>
#include <span>
#include <vector>

using namespace svmp::FE;
using namespace svmp::FE::math;

namespace {

Real multiply_entry(const std::vector<Real>& A,
                    const std::vector<Real>& B,
                    std::size_t n,
                    std::size_t row,
                    std::size_t col) {
    Real sum = Real(0);
    for (std::size_t k = 0; k < n; ++k) {
        sum += A[row * n + k] * B[k * n + col];
    }
    return sum;
}

} // namespace

TEST(DenseLinearAlgebra, InvertsScaledMatrix) {
    const std::vector<Real> A{
        Real(1.0e9), Real(2.0e6),
        Real(3.0e3), Real(4.0)
    };

    const auto inv = invert_dense_matrix(A, 2u, "scaled 2x2");
    for (std::size_t row = 0; row < 2u; ++row) {
        for (std::size_t col = 0; col < 2u; ++col) {
            const Real expected = (row == col) ? Real(1) : Real(0);
            EXPECT_NEAR(multiply_entry(A, inv, 2u, row, col), expected, Real(1.0e-10));
        }
    }
}

TEST(DenseLinearAlgebra, FactorizationSolvesMultipleRightHandSides) {
    const std::vector<Real> A{
        Real(4), Real(2), Real(0),
        Real(2), Real(5), Real(1),
        Real(0), Real(1), Real(3)
    };

    const auto solver = factor_dense_matrix(A, 3u, "symmetric 3x3");
    EXPECT_EQ(solver.diagnostics.rank, 3u);

    const std::vector<Real> rhs{Real(2), Real(4), Real(6)};
    const auto x = solver.solve(std::span<const Real>(rhs.data(), rhs.size()));
    ASSERT_EQ(x.size(), 3u);

    for (std::size_t row = 0; row < 3u; ++row) {
        Real ax = Real(0);
        for (std::size_t col = 0; col < 3u; ++col) {
            ax += A[row * 3u + col] * x[col];
        }
        EXPECT_NEAR(ax, rhs[row], Real(1.0e-12));
    }

    std::vector<Real> second_rhs{Real(1), Real(-2), Real(0.5)};
    const auto original_second_rhs = second_rhs;
    solver.solve_in_place(std::span<Real>(second_rhs.data(), second_rhs.size()));
    for (std::size_t row = 0; row < 3u; ++row) {
        Real ax = Real(0);
        for (std::size_t col = 0; col < 3u; ++col) {
            ax += A[row * 3u + col] * second_rhs[col];
        }
        EXPECT_NEAR(ax, original_second_rhs[row], Real(1.0e-12));
    }
}

TEST(DenseLinearAlgebra, FactorizationSolvesDenseRightHandSideBlock) {
    const std::vector<Real> A{
        Real(4), Real(2), Real(0),
        Real(2), Real(5), Real(1),
        Real(0), Real(1), Real(3)
    };

    const auto solver = factor_dense_matrix(A, 3u, "symmetric 3x3 block");

    std::vector<Real> rhs{
        Real(2), Real(1),
        Real(4), Real(-2),
        Real(6), Real(0.5)
    };
    const auto original_rhs = rhs;
    solver.solve_in_place(std::span<Real>(rhs.data(), rhs.size()), 2u);

    for (std::size_t rhs_col = 0; rhs_col < 2u; ++rhs_col) {
        for (std::size_t row = 0; row < 3u; ++row) {
            Real ax = Real(0);
            for (std::size_t col = 0; col < 3u; ++col) {
                ax += A[row * 3u + col] * rhs[col * 2u + rhs_col];
            }
            EXPECT_NEAR(ax, original_rhs[row * 2u + rhs_col], Real(1.0e-12));
        }
    }
}

TEST(DenseLinearAlgebra, HighConditionInverseUsesSvdFallback) {
    const std::vector<Real> high_condition{
        Real(1), Real(0),
        Real(0), Real(1.0e-13)
    };

    const auto result =
        invert_dense_matrix_with_diagnostics(high_condition, 2u, "high-condition diagonal");
    EXPECT_EQ(result.diagnostics.rank, 2u);
#if defined(FE_HAS_EIGEN) && FE_HAS_EIGEN
    EXPECT_GT(result.diagnostics.condition_estimate,
              dense_matrix_condition_fallback_threshold());
    EXPECT_TRUE(result.used_svd_fallback);
#else
    EXPECT_FALSE(result.used_svd_fallback);
#endif

    for (std::size_t row = 0; row < 2u; ++row) {
        for (std::size_t col = 0; col < 2u; ++col) {
            const Real expected = (row == col) ? Real(1) : Real(0);
            EXPECT_NEAR(multiply_entry(high_condition, result.inverse, 2u, row, col),
                        expected,
                        Real(1.0e-12));
        }
    }
}

TEST(DenseLinearAlgebra, DiagnosticValidationRejectsUnsupportedCondition) {
#if !(defined(FE_HAS_EIGEN) && FE_HAS_EIGEN)
    GTEST_SKIP() << "condition rejection requires FE_ENABLE_EIGEN diagnostics";
#endif
    DenseInverseResult result;
    result.diagnostics.rank = 2u;
    result.diagnostics.condition_estimate =
        dense_matrix_condition_error_threshold() * Real(10);

    EXPECT_GT(result.diagnostics.condition_estimate,
              dense_matrix_condition_error_threshold());
    EXPECT_THROW(validate_dense_inverse_diagnostics(
                     result, 2u, "excessive-condition diagonal"),
                 FEException);
}

TEST(DenseLinearAlgebra, ThrowsForScaleAwareSingularPivot) {
    const std::vector<Real> singular{
        Real(1.0e12), Real(2.0e12),
        Real(0.5e12), Real(1.0e12)
    };

    EXPECT_THROW((void)invert_dense_matrix(singular, 2u, "singular 2x2"),
                 FEException);
}

TEST(DenseLinearAlgebra, FactorizationThrowsForRankDeficientMatrix) {
    const std::vector<Real> singular{
        Real(1), Real(2),
        Real(2), Real(4)
    };

    EXPECT_THROW((void)factor_dense_matrix(singular, 2u, "rank-one 2x2"),
                 FEException);
}

TEST(DenseLinearAlgebra, RankUsesScaleAwareTolerance) {
    const std::vector<Real> rank_one{
        Real(1.0e8), Real(2.0e8),
        Real(3.0e8), Real(6.0e8)
    };
    EXPECT_EQ(dense_matrix_rank(rank_one, 2u, 2u), 1u);

    const std::vector<Real> full_rank{
        Real(1.0e8), Real(2.0e8),
        Real(3.0e8), Real(6.1e8)
    };
    EXPECT_EQ(dense_matrix_rank(full_rank, 2u, 2u), 2u);
}

TEST(DenseLinearAlgebra, DiagnosticsReportRankAndConditionEstimate) {
    const std::vector<Real> diagonal{
        Real(4), Real(0),
        Real(0), Real(0.5)
    };
    const auto full =
        dense_matrix_diagnostics(diagonal, 2u, 2u, "diagonal 2x2");
    EXPECT_EQ(full.rank, 2u);
#if defined(FE_HAS_EIGEN) && FE_HAS_EIGEN
    EXPECT_NEAR(full.largest_singular_value, Real(4), Real(1.0e-14));
    EXPECT_NEAR(full.smallest_retained_singular_value, Real(0.5), Real(1.0e-14));
    EXPECT_NEAR(full.condition_estimate, Real(8), Real(1.0e-14));
#else
    EXPECT_TRUE(std::isinf(full.condition_estimate));
#endif

    const std::vector<Real> rank_one{
        Real(1), Real(2),
        Real(2), Real(4)
    };
    const auto deficient =
        dense_matrix_diagnostics(rank_one, 2u, 2u, "rank-one 2x2");
    EXPECT_EQ(deficient.rank, 1u);
    EXPECT_TRUE(std::isinf(deficient.condition_estimate));
}

TEST(DenseLinearAlgebra, PseudoInverseHandlesSingularMatrixWithoutNormalEquations) {
#if !(defined(FE_HAS_EIGEN) && FE_HAS_EIGEN)
    GTEST_SKIP() << "rank-revealing pseudo-inverse requires FE_ENABLE_EIGEN";
#endif
    const std::vector<Real> rank_one{
        Real(1), Real(2),
        Real(2), Real(4)
    };

    const auto pinv =
        rank_revealing_pseudo_inverse(rank_one, 2u, 2u, "rank-one 2x2");
    EXPECT_EQ(pinv.rank, 1u);
    EXPECT_NEAR(pinv.inverse[0], Real(0.04), Real(1.0e-13));
    EXPECT_NEAR(pinv.inverse[1], Real(0.08), Real(1.0e-13));
    EXPECT_NEAR(pinv.inverse[2], Real(0.08), Real(1.0e-13));
    EXPECT_NEAR(pinv.inverse[3], Real(0.16), Real(1.0e-13));

    std::vector<Real> projection(4u, Real(0));
    for (std::size_t row = 0; row < 2u; ++row) {
        for (std::size_t col = 0; col < 2u; ++col) {
            for (std::size_t a = 0; a < 2u; ++a) {
                for (std::size_t b = 0; b < 2u; ++b) {
                    projection[row * 2u + col] +=
                        rank_one[row * 2u + a] * pinv.inverse[a * 2u + b] *
                        rank_one[b * 2u + col];
                }
            }
            EXPECT_NEAR(projection[row * 2u + col],
                        rank_one[row * 2u + col],
                        Real(1.0e-12));
        }
    }
}

TEST(DenseLinearAlgebra, PseudoInverseDropsNearZeroSingularValues) {
#if !(defined(FE_HAS_EIGEN) && FE_HAS_EIGEN)
    GTEST_SKIP() << "rank-revealing pseudo-inverse requires FE_ENABLE_EIGEN";
#endif
    const std::vector<Real> near_singular{
        Real(1), Real(0),
        Real(0), Real(1.0e-18)
    };

    const auto pinv =
        rank_revealing_pseudo_inverse(near_singular, 2u, 2u, "near-singular 2x2");
    EXPECT_EQ(pinv.rank, 1u);
    EXPECT_GT(pinv.tolerance, Real(1.0e-18));
    EXPECT_NEAR(pinv.inverse[0], Real(1), Real(1.0e-14));
    EXPECT_NEAR(pinv.inverse[1], Real(0), Real(1.0e-14));
    EXPECT_NEAR(pinv.inverse[2], Real(0), Real(1.0e-14));
    EXPECT_NEAR(pinv.inverse[3], Real(0), Real(1.0e-14));
}
