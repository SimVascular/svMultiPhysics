// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_TEST_QUADRATURE_TEST_HELPERS_H
#define SVMP_TEST_QUADRATURE_TEST_HELPERS_H

#include "FE/Common/FEException.h"
#include "FE/Quadrature/QuadratureRule.h"

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <exception>
#include <string_view>
#include <utility>

namespace svmp::FE::quadrature::test {

struct CompensatedSum {
    long double sum = 0.0L;
    long double correction = 0.0L;

    void add(long double term)
    {
        const long double next = sum + term;
        correction += std::abs(sum) >= std::abs(term)
            ? (sum - next) + term : (term - next) + sum;
        sum = next;
    }

    long double value() const { return sum + correction; }
};

inline void expect_same_ordered_rule(const QuadratureRule& left, const QuadratureRule& right)
{
    ASSERT_EQ(left.num_points(), right.num_points());
    EXPECT_EQ(left.cell_family(), right.cell_family());
    EXPECT_EQ(left.dimension(), right.dimension());
    EXPECT_EQ(left.reference_cell_measure(), right.reference_cell_measure());
    EXPECT_EQ(left.polynomial_exactness(), right.polynomial_exactness());
    for (std::size_t q = 0; q < left.num_points(); ++q) {
        ASSERT_EQ(left.point(q), right.point(q)) << "point at sample " << q;
        ASSERT_EQ(left.weight(q), right.weight(q)) << "weight at sample " << q;
    }
}

template<class Exception, class Function>
void expect_fe_exception(Function&& function, StatusCode status, std::string_view diagnostic)
{
    ASSERT_FALSE(diagnostic.empty());
    try {
        std::forward<Function>(function)();
    } catch (const Exception& error) {
        EXPECT_EQ(error.status_code(), status);
        EXPECT_GT(error.context().line(), 0);
        EXPECT_NE(error.message().find(diagnostic), std::string_view::npos) << error.what();
        return;
    } catch (const std::exception& error) {
        FAIL() << "Unexpected exception: " << error.what();
    } catch (...) {
        FAIL() << "Unexpected non-standard exception";
    }
    FAIL() << "Expected an exception containing " << diagnostic;
}

template<class Function>
void expect_invalid_argument_with_message(Function&& function, std::string_view diagnostic)
{
    expect_fe_exception<InvalidArgumentException>(
        std::forward<Function>(function), StatusCode::InvalidArgument, diagnostic);
}

} // namespace svmp::FE::quadrature::test

#endif // SVMP_TEST_QUADRATURE_TEST_HELPERS_H
