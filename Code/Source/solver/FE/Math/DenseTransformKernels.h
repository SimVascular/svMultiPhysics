// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_FE_MATH_DENSETRANSFORMKERNELS_H
#define SVMP_FE_MATH_DENSETRANSFORMKERNELS_H

#include "FEException.h"
#include "Types.h"

#include <Eigen/Core>

#include <cstddef>
#include <span>

namespace svmp {
namespace FE {
namespace math {

/// \brief Apply a row-major dense matrix to a batch of right-hand sides.
///
/// Computes output = matrix * input where matrix is rows-by-cols (row-major),
/// input holds cols rows of rhs_count values each (row stride
/// input_row_stride), and output holds rows rows of rhs_count values each
/// (row stride output_row_stride). Strides may exceed rhs_count for padded
/// layouts; padding entries are left untouched.
inline void dense_transform_batched_row_major(
    std::span<const Real> matrix,
    std::size_t rows,
    std::size_t cols,
    std::span<const Real> input,
    std::size_t input_row_stride,
    std::span<Real> output,
    std::size_t output_row_stride,
    std::size_t rhs_count) {
    if (rows == 0u || cols == 0u || rhs_count == 0u) {
        return;
    }

    svmp::throw_if<FEException>(matrix.size() < rows * cols, SVMP_HERE,
                              "dense_transform_batched_row_major: matrix span is too small");
    svmp::throw_if<FEException>(input_row_stride < rhs_count, SVMP_HERE,
                              "dense_transform_batched_row_major: input stride is smaller than RHS count");
    svmp::throw_if<FEException>(output_row_stride < rhs_count, SVMP_HERE,
                              "dense_transform_batched_row_major: output stride is smaller than RHS count");
    svmp::throw_if<FEException>(
        input.size() < (cols - 1u) * input_row_stride + rhs_count, SVMP_HERE,
        "dense_transform_batched_row_major: input span is too small");
    svmp::throw_if<FEException>(
        output.size() < (rows - 1u) * output_row_stride + rhs_count, SVMP_HERE,
        "dense_transform_batched_row_major: output span is too small");

    using RowMajorMatrix =
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using ConstMap = Eigen::Map<const RowMajorMatrix>;
    using ConstStridedMap =
        Eigen::Map<const RowMajorMatrix, Eigen::Unaligned, Eigen::OuterStride<>>;
    using StridedMap =
        Eigen::Map<RowMajorMatrix, Eigen::Unaligned, Eigen::OuterStride<>>;

    const ConstMap matrix_map(matrix.data(),
                              static_cast<Eigen::Index>(rows),
                              static_cast<Eigen::Index>(cols));
    const ConstStridedMap input_map(
        input.data(),
        static_cast<Eigen::Index>(cols),
        static_cast<Eigen::Index>(rhs_count),
        Eigen::OuterStride<>(static_cast<Eigen::Index>(input_row_stride)));
    StridedMap output_map(
        output.data(),
        static_cast<Eigen::Index>(rows),
        static_cast<Eigen::Index>(rhs_count),
        Eigen::OuterStride<>(static_cast<Eigen::Index>(output_row_stride)));

    output_map.noalias() = matrix_map * input_map;
}

} // namespace math
} // namespace FE
} // namespace svmp

#endif // SVMP_FE_MATH_DENSETRANSFORMKERNELS_H
