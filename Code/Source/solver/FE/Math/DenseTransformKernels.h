/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See License file.
 */

#ifndef SVMP_FE_MATH_DENSETRANSFORMKERNELS_H
#define SVMP_FE_MATH_DENSETRANSFORMKERNELS_H

#include "Types.h"

#include <algorithm>
#include <array>
#include <cstddef>

namespace svmp {
namespace FE {
namespace math {

constexpr std::size_t dense_transform_blocked_min_rows() noexcept { return 32u; }
constexpr std::size_t dense_transform_blocked_min_rhs() noexcept { return 4u; }

inline void dense_transform_batched_row_major(
    const Real* SVMP_RESTRICT matrix,
    std::size_t rows,
    std::size_t cols,
    const Real* SVMP_RESTRICT input,
    std::size_t input_row_stride,
    Real* SVMP_RESTRICT output,
    std::size_t output_row_stride,
    std::size_t rhs_count) {
    if (rows == 0u || cols == 0u || rhs_count == 0u) {
        return;
    }

    if (rows < dense_transform_blocked_min_rows() ||
        rhs_count < dense_transform_blocked_min_rhs()) {
        for (std::size_t row = 0; row < rows; ++row) {
            const Real* matrix_row = matrix + row * cols;
            Real* output_row = output + row * output_row_stride;
            for (std::size_t rhs = 0; rhs < rhs_count; ++rhs) {
                Real value = Real(0);
                for (std::size_t col = 0; col < cols; ++col) {
                    value += matrix_row[col] * input[col * input_row_stride + rhs];
                }
                output_row[rhs] = value;
            }
        }
        return;
    }

    constexpr std::size_t kRhsBlock = 32u;
    for (std::size_t row = 0; row < rows; ++row) {
        const Real* matrix_row = matrix + row * cols;
        Real* output_row = output + row * output_row_stride;
        for (std::size_t rhs_base = 0; rhs_base < rhs_count; rhs_base += kRhsBlock) {
            const std::size_t block_size = std::min(kRhsBlock, rhs_count - rhs_base);
            std::array<Real, kRhsBlock> accum{};
            for (std::size_t col = 0; col < cols; ++col) {
                const Real coeff = matrix_row[col];
                const Real* input_row = input + col * input_row_stride + rhs_base;
                for (std::size_t rhs = 0; rhs < block_size; ++rhs) {
                    accum[rhs] += coeff * input_row[rhs];
                }
            }
            for (std::size_t rhs = 0; rhs < block_size; ++rhs) {
                output_row[rhs_base + rhs] = accum[rhs];
            }
        }
    }
}

} // namespace math
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_MATH_DENSETRANSFORMKERNELS_H
