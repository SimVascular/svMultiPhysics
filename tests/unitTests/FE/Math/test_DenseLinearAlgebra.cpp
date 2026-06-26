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

double multiply_entry(const std::vector<double>& A,
                    const std::vector<double>& B,
                    std::size_t n,
                    std::size_t row,
                    std::size_t col) {
    double sum = double(0);
    for (std::size_t k = 0; k < n; ++k) {
        sum += A[row * n + k] * B[k * n + col];
    }
    return sum;
}

} // namespace

TEST(DenseLinearAlgebra, InvertsScaledMatrix) {
    const std::vector<double> A{
        double(1.0e9), double(2.0e6),
        double(3.0e3), double(4.0)
    };

    const auto inv = invert_dense_matrix(A, 2u, "scaled 2x2");
    for (std::size_t row = 0; row < 2u; ++row) {
        for (std::size_t col = 0; col < 2u; ++col) {
            const double expected = (row == col) ? double(1) : double(0);
            EXPECT_NEAR(multiply_entry(A, inv, 2u, row, col), expected, double(1.0e-10));
        }
    }
}

TEST(DenseLinearAlgebra, FactorizationSolvesMultipleRightHandSides) {
    const std::vector<double> A{
        double(4), double(2), double(0),
        double(2), double(5), double(1),
        double(0), double(1), double(3)
    };

    const auto solver = factor_dense_matrix(A, 3u, "symmetric 3x3");
    EXPECT_EQ(solver.diagnostics.rank, 3u);

    const std::vector<double> rhs{double(2), double(4), double(6)};
    const auto x = solver.solve(std::span<const double>(rhs.data(), rhs.size()));
    ASSERT_EQ(x.size(), 3u);

    for (std::size_t row = 0; row < 3u; ++row) {
        double ax = double(0);
        for (std::size_t col = 0; col < 3u; ++col) {
            ax += A[row * 3u + col] * x[col];
        }
        EXPECT_NEAR(ax, rhs[row], double(1.0e-12));
    }

    std::vector<double> second_rhs{double(1), double(-2), double(0.5)};
    const auto original_second_rhs = second_rhs;
    solver.solve_in_place(std::span<double>(second_rhs.data(), second_rhs.size()));
    for (std::size_t row = 0; row < 3u; ++row) {
        double ax = double(0);
        for (std::size_t col = 0; col < 3u; ++col) {
            ax += A[row * 3u + col] * second_rhs[col];
        }
        EXPECT_NEAR(ax, original_second_rhs[row], double(1.0e-12));
    }
}

TEST(DenseLinearAlgebra, FactorizationSolvesDenseRightHandSideBlock) {
    const std::vector<double> A{
        double(4), double(2), double(0),
        double(2), double(5), double(1),
        double(0), double(1), double(3)
    };

    const auto solver = factor_dense_matrix(A, 3u, "symmetric 3x3 block");

    std::vector<double> rhs{
        double(2), double(1),
        double(4), double(-2),
        double(6), double(0.5)
    };
    const auto original_rhs = rhs;
    solver.solve_in_place(std::span<double>(rhs.data(), rhs.size()), 2u);

    for (std::size_t rhs_col = 0; rhs_col < 2u; ++rhs_col) {
        for (std::size_t row = 0; row < 3u; ++row) {
            double ax = double(0);
            for (std::size_t col = 0; col < 3u; ++col) {
                ax += A[row * 3u + col] * rhs[col * 2u + rhs_col];
            }
            EXPECT_NEAR(ax, original_rhs[row * 2u + rhs_col], double(1.0e-12));
        }
    }
}

// Every other matrix in this file already has its largest pivot on the
// diagonal, so these cases cover the row-exchange branch in factor_dense_matrix,
// the inverse path used by SerendipityBasis, and the permutation replay in
// solve_in_place.
TEST(DenseLinearAlgebra, FactorizationPivotsThroughZeroLeadingDiagonal) {
    const std::vector<double> swap_2x2{
        double(0), double(1),
        double(1), double(0)
    };

    const auto solver = factor_dense_matrix(swap_2x2, 2u, "swap 2x2");
    const std::vector<double> rhs{double(3), double(7)};
    const auto x = solver.solve(std::span<const double>(rhs.data(), rhs.size()));
    ASSERT_EQ(x.size(), 2u);
    EXPECT_NEAR(x[0], double(7), double(1.0e-14));
    EXPECT_NEAR(x[1], double(3), double(1.0e-14));

    const auto inv = invert_dense_matrix(swap_2x2, 2u, "swap 2x2");
    for (std::size_t row = 0; row < 2u; ++row) {
        for (std::size_t col = 0; col < 2u; ++col) {
            EXPECT_NEAR(inv[row * 2u + col], swap_2x2[row * 2u + col], double(1.0e-14));
        }
    }

    // Every column requires a row exchange during elimination.
    const std::vector<double> permuted_scaled{
        double(0), double(0), double(1), double(0),
        double(1), double(0), double(0), double(0),
        double(0), double(0), double(0), double(2),
        double(0), double(3), double(0), double(0)
    };

    const auto inv4 = invert_dense_matrix(permuted_scaled, 4u, "permuted scaled 4x4");
    for (std::size_t row = 0; row < 4u; ++row) {
        for (std::size_t col = 0; col < 4u; ++col) {
            const double expected = (row == col) ? double(1) : double(0);
            EXPECT_NEAR(multiply_entry(permuted_scaled, inv4, 4u, row, col),
                        expected,
                        double(1.0e-14));
        }
    }
}

TEST(DenseLinearAlgebra, WideMultiRhsSolveWithPivoting) {
    // Requires a row swap in column 0 and uses a wide right-hand-side block to
    // exercise the row-interleaved multi-RHS layout end to end.
    const std::vector<double> A{
        double(0), double(2), double(1),
        double(4), double(1), double(0),
        double(1), double(0), double(3)
    };
    constexpr std::size_t kRhsCount = 33u;

    const auto solver = factor_dense_matrix(A, 3u, "pivoting 3x3");

    std::vector<double> rhs(3u * kRhsCount, double(0));
    for (std::size_t row = 0; row < 3u; ++row) {
        for (std::size_t r = 0; r < kRhsCount; ++r) {
            rhs[row * kRhsCount + r] =
                double(1) + static_cast<double>(row) - double(0.25) * static_cast<double>(r % 7u);
        }
    }
    const auto original_rhs = rhs;

    solver.solve_in_place(std::span<double>(rhs.data(), rhs.size()), kRhsCount);

    for (std::size_t r = 0; r < kRhsCount; ++r) {
        for (std::size_t row = 0; row < 3u; ++row) {
            double ax = double(0);
            for (std::size_t col = 0; col < 3u; ++col) {
                ax += A[row * 3u + col] * rhs[col * kRhsCount + r];
            }
            EXPECT_NEAR(ax, original_rhs[row * kRhsCount + r], double(1.0e-12))
                << "rhs column " << r << ", row " << row;
        }
    }
}

TEST(DenseLinearAlgebra, SolveInPlaceValidatesInputs) {
    const std::vector<double> identity{
        double(1), double(0),
        double(0), double(1)
    };
    const auto solver = factor_dense_matrix(identity, 2u, "identity 2x2");

    std::vector<double> rhs{double(1), double(2)};
    EXPECT_THROW(solver.solve_in_place(std::span<double>(rhs.data(), rhs.size()), 0u),
                 FEException);

    std::vector<double> wrong_size{double(1), double(2), double(3)};
    EXPECT_THROW(
        solver.solve_in_place(std::span<double>(wrong_size.data(), wrong_size.size()), 1u),
        FEException);

    DenseLUSolver unfactored;
    unfactored.n = 2u;
    unfactored.error_message_label = "unfactored";
    EXPECT_FALSE(unfactored.empty());
    EXPECT_THROW(unfactored.solve_in_place(std::span<double>(rhs.data(), rhs.size()), 1u),
                 FEException);
}

TEST(DenseLinearAlgebra, DiagnosticValidationRejectsRankMismatch) {
    DenseInverseResult result;
    result.diagnostics.rank = 1u;

    EXPECT_THROW(validate_dense_inverse_diagnostics(result, 2u, "rank mismatch"),
                 FEException);
}

TEST(DenseLinearAlgebra, RankHandlesNonSquareMatrices) {
    const std::vector<double> wide_full{
        double(1), double(0), double(2),
        double(0), double(1), double(-1)
    };
    EXPECT_EQ(dense_matrix_rank(wide_full, 2u, 3u), 2u);

    const std::vector<double> tall_rank_one{
        double(1), double(2),
        double(2), double(4),
        double(3), double(6)
    };
    EXPECT_EQ(dense_matrix_rank(tall_rank_one, 3u, 2u), 1u);
}

TEST(DenseLinearAlgebra, HighConditionInverseUsesSvdFallback) {
    const std::vector<double> high_condition{
        double(1), double(0),
        double(0), double(1.0e-13)
    };

    const auto result =
        invert_dense_matrix_with_diagnostics(high_condition, 2u, "high-condition diagonal");
    EXPECT_EQ(result.diagnostics.rank, 2u);
    EXPECT_GT(result.diagnostics.condition_estimate,
              dense_matrix_condition_fallback_threshold());
    EXPECT_TRUE(result.used_svd_fallback);

    for (std::size_t row = 0; row < 2u; ++row) {
        for (std::size_t col = 0; col < 2u; ++col) {
            const double expected = (row == col) ? double(1) : double(0);
            EXPECT_NEAR(multiply_entry(high_condition, result.inverse, 2u, row, col),
                        expected,
                        double(1.0e-12));
        }
    }
}

TEST(DenseLinearAlgebra, DiagnosticValidationRejectsUnsupportedCondition) {
    DenseInverseResult result;
    result.diagnostics.rank = 2u;
    result.diagnostics.condition_estimate =
        dense_matrix_condition_error_threshold() * double(10);

    EXPECT_GT(result.diagnostics.condition_estimate,
              dense_matrix_condition_error_threshold());
    EXPECT_THROW(validate_dense_inverse_diagnostics(
                     result, 2u, "excessive-condition diagonal"),
                 FEException);
}

TEST(DenseLinearAlgebra, ThrowsForScaleAwareSingularPivot) {
    const std::vector<double> singular{
        double(1.0e12), double(2.0e12),
        double(0.5e12), double(1.0e12)
    };

    EXPECT_THROW((void)invert_dense_matrix(singular, 2u, "singular 2x2"),
                 FEException);
}

TEST(DenseLinearAlgebra, FactorizationThrowsForRankDeficientMatrix) {
    const std::vector<double> singular{
        double(1), double(2),
        double(2), double(4)
    };

    EXPECT_THROW((void)factor_dense_matrix(singular, 2u, "rank-one 2x2"),
                 FEException);
}

TEST(DenseLinearAlgebra, RankUsesScaleAwareTolerance) {
    const std::vector<double> rank_one{
        double(1.0e8), double(2.0e8),
        double(3.0e8), double(6.0e8)
    };
    EXPECT_EQ(dense_matrix_rank(rank_one, 2u, 2u), 1u);

    const std::vector<double> full_rank{
        double(1.0e8), double(2.0e8),
        double(3.0e8), double(6.1e8)
    };
    EXPECT_EQ(dense_matrix_rank(full_rank, 2u, 2u), 2u);
}

TEST(DenseLinearAlgebra, DiagnosticsReportRankAndConditionEstimate) {
    const std::vector<double> diagonal{
        double(4), double(0),
        double(0), double(0.5)
    };
    const auto full =
        dense_matrix_diagnostics(diagonal, 2u, 2u, "diagonal 2x2");
    EXPECT_EQ(full.rank, 2u);
    EXPECT_NEAR(full.largest_singular_value, double(4), double(1.0e-14));
    EXPECT_NEAR(full.smallest_retained_singular_value, double(0.5), double(1.0e-14));
    EXPECT_NEAR(full.condition_estimate, double(8), double(1.0e-14));

    const std::vector<double> rank_one{
        double(1), double(2),
        double(2), double(4)
    };
    const auto deficient =
        dense_matrix_diagnostics(rank_one, 2u, 2u, "rank-one 2x2");
    EXPECT_EQ(deficient.rank, 1u);
    EXPECT_TRUE(std::isinf(deficient.condition_estimate));
}

TEST(DenseLinearAlgebra, PseudoInverseHandlesSingularMatrixWithoutNormalEquations) {
    const std::vector<double> rank_one{
        double(1), double(2),
        double(2), double(4)
    };

    const auto pinv =
        rank_revealing_pseudo_inverse(rank_one, 2u, 2u, "rank-one 2x2");
    EXPECT_EQ(pinv.rank, 1u);
    EXPECT_NEAR(pinv.inverse[0], double(0.04), double(1.0e-13));
    EXPECT_NEAR(pinv.inverse[1], double(0.08), double(1.0e-13));
    EXPECT_NEAR(pinv.inverse[2], double(0.08), double(1.0e-13));
    EXPECT_NEAR(pinv.inverse[3], double(0.16), double(1.0e-13));

    std::vector<double> projection(4u, double(0));
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
                        double(1.0e-12));
        }
    }
}

TEST(DenseLinearAlgebra, PseudoInverseDropsNearZeroSingularValues) {
    const std::vector<double> near_singular{
        double(1), double(0),
        double(0), double(1.0e-18)
    };

    const auto pinv =
        rank_revealing_pseudo_inverse(near_singular, 2u, 2u, "near-singular 2x2");
    EXPECT_EQ(pinv.rank, 1u);
    EXPECT_GT(pinv.tolerance, double(1.0e-18));
    EXPECT_NEAR(pinv.inverse[0], double(1), double(1.0e-14));
    EXPECT_NEAR(pinv.inverse[1], double(0), double(1.0e-14));
    EXPECT_NEAR(pinv.inverse[2], double(0), double(1.0e-14));
    EXPECT_NEAR(pinv.inverse[3], double(0), double(1.0e-14));
}
