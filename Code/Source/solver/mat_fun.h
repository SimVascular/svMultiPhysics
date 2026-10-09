// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef MAT_FUN_H 
#define MAT_FUN_H 
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"
#include "eigen3/unsupported/Eigen/CXX11/Tensor"
#include <stdexcept>

#include "Array.h"
#include "consts.h"
#include "Tensor4.h"
#include "Vector.h"
#include "FE/Common/FEException.h"

/// @brief The classes defined here duplicate the data structures in the 
/// Fortran MATFUN module defined in MATFUN.f. 
///
/// This module defines data structures for generally performed matrix and tensor operations.
///
/// \todo [TODO:DaveP] this should just be a namespace?
//
namespace mat_fun {
    /// @brief A 2nd order tensor, nsd x nsd, fixed size and stack allocated.
    /// Used for the deformation gradient, stresses and similar quantities that
    /// have a known size at compile time.
    template<int nsd>
    using Matrix = Eigen::Matrix<double, nsd, nsd>;

    /// @brief A 4th order tensor, nsd x nsd x nsd x nsd, fixed size and stack
    /// allocated. Used for the material elasticity tensor and other 4th order tensors
    /// that have a known size at compile time.
    template<int nsd>
    using Tensor = Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>>;

    /// @brief One nsd-vector per element node, so nsd x eNoN. Row count is fixed
    /// at compile time while column count is the element's node count, known only
    /// at run time, so it is bounded by consts::maxNoN to stay stack allocated.
    /// Used for shape function gradients and other per-node vector quantities.
    template <int nsd>
    using NodalMatrix = Eigen::Matrix<double, nsd, Eigen::Dynamic, 0, nsd, consts::maxNoN>;

    /// @brief One scalar per element node, so eNoN entries. Dynamic length bounded
    /// by consts::maxNoN to stay stack allocated, as for NodalMatrix. Used for shape
    /// function values and other per-node scalar quantities.
    using NodalVector = Eigen::Matrix<double, Eigen::Dynamic, 1, 0, consts::maxNoN, 1>;

    // The eigen_view overloads below wrap an Array or Vector in an Eigen::Map that
    // shares its storage, so the container must outlive the view.

    /// @brief Read-only Eigen view of an Array, sharing its storage.
    ///
    /// @tparam rows Row count, fixed at compile time; the columns are taken from the Array.
    template <int rows>
    Eigen::Map<const Eigen::Matrix<double, rows, Eigen::Dynamic>>
    eigen_view(const Array<double>& A) {
        if (A.nrows() != rows) {
          svmp::raise<svmp::FE::InvalidArgumentException>(
              "A view of " + std::to_string(rows) + " rows was requested for an array with " +
              std::to_string(A.nrows()) + " rows.");
        }
        return {A.data(), rows, A.ncols()};
    }

    /// @brief Read-only Eigen view of a whole Array, sharing its storage.
    inline Eigen::Map<const Eigen::MatrixXd>
    eigen_view(const Array<double>& A) {
        return {A.data(), A.nrows(), A.ncols()};
    }

    /// @brief Read-only Eigen view of a band of rows of an Array, sharing its storage.
    ///
    /// Used where an Array contains several fields and only a contiguous band is viewed.
    /// For example,
    ///
    ///     A = [ A00, A01, A02 ;         eigen_view_rows<2>(A, 1) = [ A10, A11, A12 ;
    ///           A10, A11, A12 ;                                      A20, A21, A22 ]
    ///           A20, A21, A22 ;
    ///           A30, A31, A32 ]
    ///
    /// @param A The Array to view.
    /// @param first Index of the band's first row within the Array.
    /// @tparam rows Number of rows, fixed at compile time, typically nsd.
    template <int rows>
    Eigen::Map<const Eigen::Matrix<double, rows, Eigen::Dynamic>, 0, Eigen::OuterStride<>>
    eigen_view_rows(const Array<double>& A, const int first) {
        if (first < 0 || first + rows > A.nrows()) {
          svmp::raise<svmp::FE::InvalidArgumentException>(
              "A view of " + std::to_string(rows) + " rows starting at row " + std::to_string(first) +
              " was requested for an array with " + std::to_string(A.nrows()) + " rows.");
        }
        return {A.data() + first, rows, A.ncols(), Eigen::OuterStride<>(A.nrows())};
    }

    /// @brief Writable Eigen view of a whole Array, sharing its storage.
    inline Eigen::Map<Eigen::MatrixXd>
    eigen_view_mutable(Array<double>& A) {
        return {A.data(), A.nrows(), A.ncols()};
    }

    /// @brief Read-only Eigen view of a whole Vector, sharing its storage.
    inline Eigen::Map<const Eigen::VectorXd>
    eigen_view(const Vector<double>& v) {
        return {v.data(), v.size()};
    }

    /// @brief Read-only Eigen view of a Vector, sharing its storage.
    ///
    /// @tparam rows Entry count, fixed at compile time.
    template <int rows>
    Eigen::Map<const Eigen::Matrix<double, rows, 1>>
    eigen_view(const Vector<double>& v) {
        if (v.size() != rows) {
          svmp::raise<svmp::FE::InvalidArgumentException>(
              "A view of " + std::to_string(rows) + " entries was requested for a vector with " +
              std::to_string(v.size()) + " entries.");
        }
        return Eigen::Map<const Eigen::Matrix<double, rows, 1>>(v.data());
    }

    // Function to convert Array<double> to Eigen::Matrix
    template <typename MatrixType>
    MatrixType convert_to_eigen_matrix(const Array<double>& src) {
        MatrixType mat;
        for (int i = 0; i < mat.rows(); ++i)
            for (int j = 0; j < mat.cols(); ++j)
                mat(i, j) = src(i, j);
        return mat;
    }

    // Function to convert Eigen::Matrix to Array<double>
    template <typename MatrixType>
    void convert_to_array(const MatrixType& mat, Array<double>& dest) {
        for (int i = 0; i < mat.rows(); ++i)
            for (int j = 0; j < mat.cols(); ++j)
                dest(i, j) = mat(i, j);
    }

    // Function to convert a higher-dimensional array like Dm
    template <typename MatrixType>
    void copy_Dm(const MatrixType& mat, Array<double>& dest) {
        if ((mat.rows() != dest.nrows()) || (mat.cols() != dest.ncols())) {
          const std::string mat_dims = "(" + std::to_string(mat.rows()) + "x" + std::to_string(mat.cols()) + ")";
          const std::string dest_dims = "(" + std::to_string(dest.nrows()) + "x" + std::to_string(dest.ncols()) + ")";
          svmp::raise<svmp::FE::InvalidArgumentException>(
              "The 'mat" + mat_dims + "' and 'dest" + dest_dims +
              "' arrays have incompatible sizes.");
        }

        for (int i = 0; i < mat.rows(); ++i) {
            for (int j = 0; j < mat.cols(); ++j) {
                dest(i, j) = mat(i, j);
            }
        }
    }

    template <int nsd>
    Eigen::Matrix<double, nsd, 1> cross_product(const Eigen::Matrix<double, nsd, 1>& u, const Eigen::Matrix<double, nsd, 1>& v) {
        if constexpr (nsd == 2) {
            return Eigen::Matrix<double, 2, 1>(v(1), - v(0));
        }
        else if constexpr (nsd == 3) {
            return u.cross(v);
        }
        else {
            throw std::runtime_error("[cross_product] Invalid number of spatial dimensions '" + std::to_string(nsd) + "'. Valid dimensions are 2 or 3.");
        }
    }

    double mat_ddot(const Array<double>& A, const Array<double>& B, const int nd);
    
    template <int nsd>
    double double_dot_product(const Matrix<nsd>& A, const Matrix<nsd>& B) {
        return A.cwiseProduct(B).sum();
    }
    
    double mat_det(const Array<double>& A, const int nd);
    Array<double> mat_dev(const Array<double>& A, const int nd);

    Array<double> mat_dyad_prod(const Vector<double>& u, const Vector<double>& v, const int nd);

    Array<double> mat_id(const int nsd);
    Array<double> mat_inv(const Array<double>& A, const int nd, bool debug = false);
    Array<double> mat_inv_ge(const Array<double>& A, const int nd, bool debug = false);
    Array<double> mat_inv_lp(const Array<double>& A, const int nd);

    /// @brief Multiply a matrix by a vector, returning A*v.
    ///
    /// @param[in] A matrix with as many columns as v has entries.
    /// @param[in] v vector.
    /// @return the product, of size rows(A).
    ///
    /// Throws InvalidArgumentException if the sizes are incompatible.
    Vector<double> mat_mul(const Array<double>& A, const Vector<double>& v);

    /// @brief Multiply two matrices, returning A*B.
    ///
    /// @param[in] A left operand.
    /// @param[in] B right operand, with as many rows as A has columns.
    /// @return the product, of size rows(A) by cols(B).
    ///
    /// Throws InvalidArgumentException if the sizes are incompatible. The
    /// result is freshly allocated, so the operands may alias it, as in
    /// A = mat_mul(A, B).
    Array<double> mat_mul(const Array<double>& A, const Array<double>& B);

    /// @brief Multiply two matrices, writing A*B into an existing result.
    ///
    /// @param[in]  A left operand.
    /// @param[in]  B right operand, with as many rows as A has columns.
    /// @param[out] result the product. The caller sizes it rows(A) by cols(B).
    ///
    /// Throws InvalidArgumentException if the sizes are incompatible. Preferred
    /// in loops, where it reuses the caller's storage instead of allocating a
    /// result on every call. The result must not alias A or B.
    void mat_mul(const Array<double>& A, const Array<double>& B, Array<double>& result);

    /// @brief Matrix product with the operand shape supplied at compile time.
    ///
    /// Overloads of mat_mul rather than differently named helpers, so a call
    /// site states the shape and otherwise reads exactly as before:
    ///
    /// @code
    ///   mat_mul(Dm, Bm.rslice(b), DBm);        // runtime shape check
    ///   mat_mul<6, 6, 3>(Dm, Bm.rslice(b), DBm);  // no check, same arguments
    /// @endcode
    ///
    /// The generic mat_mul overload above will dispatch to this overload
    /// when the shapes are known at compile time.
    ///
    /// @tparam M rows of A and of the result
    /// @tparam K columns of A and rows of B, the contracted dimension
    /// @tparam N columns of B and of the result
    template <int M, int K, int N>
    void mat_mul(const Array<double>& A, const Array<double>& B,
                 Array<double>& C)
    {
      Eigen::Map<const Eigen::Matrix<double, M, K>> a(A.data());
      Eigen::Map<const Eigen::Matrix<double, K, N>> b(B.data());
      Eigen::Map<Eigen::Matrix<double, M, N>>       c(C.data());

      c.noalias() = a * b;
    }

    /// @brief As above, but with the column count known only at run time.
    ///
    /// For operands with one column per element node, where the width depends on
    /// the element type. The row counts are still compile-time, which is where
    /// most of the benefit comes from.
    template <int M, int K>
    void mat_mul(const Array<double>& A, const Array<double>& B,
                 Array<double>& C)
    {
      Eigen::Map<const Eigen::Matrix<double, M, K>> a(A.data());
      Eigen::Map<const Eigen::Matrix<double, K, Eigen::Dynamic>> b(B.data(), K, B.ncols());
      Eigen::Map<Eigen::Matrix<double, M, Eigen::Dynamic>>       c(C.data(), M, C.ncols());

      c.noalias() = a * b;
    }

    Array<double> mat_symm(const Array<double>& A, const int nd);
    Array<double> mat_symm_prod(const Vector<double>& u, const Vector<double>& v, const int nd);

    double mat_trace(const Array<double>& A, const int nd);

    Tensor4<double> ten_asym_prod12(const Array<double>& A, const Array<double>& B, const int nd);
    Tensor4<double> ten_ddot(const Tensor4<double>& A, const Tensor4<double>& B, const int nd);
    Tensor4<double> ten_ddot_2412(const Tensor4<double>& A, const Tensor4<double>& B, const int nd);
    Tensor4<double> ten_ddot_3424(const Tensor4<double>& A, const Tensor4<double>& B, const int nd);

    /**
     * @brief Contracts two 4th order tensors A and B over two dimensions.
     *
     * For example, if dimsA = {0, 1} and dimsB = {2, 3} this is
     *  C_klmn = A_ijkl B_mnij   (sum over i, j)
     *
     * @tparam nsd Number of spatial dimensions; each tensor is nsd^4.
     * @param[in] A,B Fourth order tensors to contract.
     * @param[in] dimsA,dimsB Indices of the contracted dimensions of A and B.
     * @return The contracted tensor.
     */
    template <int nsd>
    Tensor<nsd>
    double_dot_product(const Tensor<nsd>& A, const std::array<int, 2>& dimsA, 
                        const Tensor<nsd>& B, const std::array<int, 2>& dimsB) {
        
        // Fast path for dimsA = dimsB = {2,3}: C_ijmn = A_ijkl * B_mnkl.
        if (dimsA[0] == 2 && dimsA[1] == 3 && dimsB[0] == 2 && dimsB[1] == 3) {
            constexpr int N = nsd * nsd;
            Tensor<nsd> C;
            Eigen::Map<const Eigen::Matrix<double, N, N>> a(A.data());
            Eigen::Map<const Eigen::Matrix<double, N, N>> b(B.data());
            Eigen::Map<Eigen::Matrix<double, N, N>> c(C.data());
            c.noalias() = a * b.transpose();
            return C;
        }

        // Define the contraction dimensions
        Eigen::array<Eigen::IndexPair<int>, 2> contractionDims = {
            Eigen::IndexPair<int>(dimsA[0], dimsB[0]), // Contract A's dimsA[0] with B's dimsB[0]
            Eigen::IndexPair<int>(dimsA[1], dimsB[1])  // Contract A's dimsA[1] with B's dimsB[1]
        };

        // Return the double dot product
        return A.contract(B, contractionDims);
    }

    Tensor4<double> ten_dyad_prod(const Array<double>& A, const Array<double>& B, const int nd);
    
    /**
     * @brief Compute the dyadic product of two 2nd order tensors A and B, C_ijkl = A_ij * B_kl
     * 
     * @tparam nsd, the number of spatial dimensions
     * @param A, the first 2nd order tensor
     * @param B, the second 2nd order tensor
     * @return Tensor<nsd>
     */
    template <int nsd>
    Tensor<nsd> 
    dyadic_product(const Matrix<nsd>& A, const Matrix<nsd>& B) {
        // Initialize the result tensor
        Tensor<nsd> C;
        constexpr int N = nsd * nsd;

        // Column-major storage flattens index pairs: c(ij,kl) = a(ij) * b(kl).
        Eigen::Map<const Eigen::Matrix<double, N, 1>> a(A.data());
        Eigen::Map<const Eigen::Matrix<double, N, 1>> b(B.data());
        Eigen::Map<Eigen::Matrix<double, N, N>> c(C.data());
        c.noalias() = a * b.transpose();
        return C;
    }

    Tensor4<double> ten_ids(const int nd);

    /**
     * @brief Create a 4th order identity tensor:
     * I_ijkl = 0.5 * (δ_ik * δ_jl + δ_il * δ_jk)
     * 
     * @tparam nsd, the number of spatial dimensions
     * @return Tensor<nsd> 
     */
    template <int nsd>
    Tensor<nsd>
    fourth_order_identity() {
        // Initialize as zero
        Tensor<nsd> I;
        I.setZero();

        // Set only non-zero entries
        for (int i = 0; i < nsd; ++i) {
            for (int j = 0; j < nsd; ++j) {
                I(i,j,i,j) += 0.5;
                I(i,j,j,i) += 0.5;
            }
        }

        return I;
    }

    Array<double> ten_mddot(const Tensor4<double>& A, const Array<double>& B, const int nd);

    Tensor4<double> ten_symm_prod(const Array<double>& A, const Array<double>& B, const int nd);
    
    /// @brief Create a 4th order tensor from symmetric outer product of two matrices: C_ijkl = 0.5 * (A_ik * B_jl + A_il * B_jk)
    ///
    /// @tparam nsd Number of spatial dimensions.
    /// @param[in] A,B Second order tensors.
    /// @return The resulting 4th order tensor.
    template <int nsd>
    Tensor<nsd>
    symmetric_dyadic_product(const Matrix<nsd>& A, const Matrix<nsd>& B) {
        
        // Initialize the result tensor
        Tensor<nsd> C;

        // Compute the symmetric product: C_ijkl = 0.5 * (A_ik * B_jl + A_il * B_jk)
        for (int l = 0; l < nsd; ++l) {
            for (int k = 0; k < nsd; ++k) {
                // blk views the (k,l) block of C, so blk(i,j) is C(i,j,k,l).
                Eigen::Map<Eigen::Matrix<double, nsd, nsd>> blk(C.data() + nsd * nsd * (k + nsd * l));
                blk.noalias() = 0.5 * (A.col(k) * B.col(l).transpose()
                                     + A.col(l) * B.col(k).transpose());
            }
        }

        // Return the symmetric product
        return C;
    }

    Tensor4<double> ten_transpose(const Tensor4<double>& A, const int nd);

    /**
     * @brief Performs a tensor transpose operation on a 4th order tensor A, B_ijkl = A_klij
     * 
     * @tparam nsd, the number of spatial dimensions
     * @param A, the input 4th order tensor
     * @return Tensor<nsd>
     */
    template <int nsd>
    Tensor<nsd>
    transpose(const Tensor<nsd>& A) {

        // Initialize the result tensor
        Tensor<nsd> B;

        // Permute the tensor indices to perform the transpose operation
        for (int i = 0; i < nsd; ++i) {
            for (int j = 0; j < nsd; ++j) {
                for (int k = 0; k < nsd; ++k) {
                    for (int l = 0; l < nsd; ++l) {
                        B(i,j,k,l) = A(k,l,i,j);
                    }
                }
            }
        }

        return B;
    }

    Array<double> transpose(const Array<double>& A);

    void ten_init(const int nd);

};

#endif
