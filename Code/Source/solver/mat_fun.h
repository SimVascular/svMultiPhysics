/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef MAT_FUN_H 
#define MAT_FUN_H 
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"
#include "eigen3/unsupported/Eigen/CXX11/Tensor"
#include <stdexcept>

#include "Array.h"
#include "Tensor4.h"
#include "Vector.h"

/// @brief The classes defined here duplicate the data structures in the 
/// Fortran MATFUN module defined in MATFUN.f. 
///
/// This module defines data structures for generally performed matrix and tensor operations.
///
/// \todo [TODO:DaveP] this should just be a namespace?
//
namespace mat_fun {
    // Helper function to convert Array<double> to Eigen::Matrix
    template <typename MatrixType>
    MatrixType toEigenMatrix(const Array<double>& src) {
        MatrixType mat;
        for (int i = 0; i < mat.rows(); ++i)
            for (int j = 0; j < mat.cols(); ++j)
                mat(i, j) = src(i, j);
        return mat;
    }

    // Helper function to convert Eigen::Matrix to Array<double>
    template <typename MatrixType>
    void toArray(const MatrixType& mat, Array<double>& dest) {
        for (int i = 0; i < mat.rows(); ++i)
            for (int j = 0; j < mat.cols(); ++j)
                dest(i, j) = mat(i, j);
    }

    // Helper function to convert a higher-dimensional array like Dm
    template <typename MatrixType>
    void copyDm(const MatrixType& mat, Array<double>& dest, int rows, int cols) {
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                dest(i, j) = mat(i, j);
    }

    template <int nsd>
    Eigen::Matrix<double, nsd, 1> cross_eigen(const Eigen::Matrix<double, nsd, 1>& u, const Eigen::Matrix<double, nsd, 1>& v) {
        if constexpr (nsd == 2) {
            return Eigen::Matrix<double, 2, 1>(v(1), - v(0));
        }
        else if constexpr (nsd == 3) {
            return u.cross(v);
        }
        else {
            throw std::runtime_error("Invalid number of spatial dimensions");
        }
    }

    double mat_ddot(const Array<double>& A, const Array<double>& B, const int nd);
    
    template <int nsd>
    double mat_ddot_eigen(const Eigen::Matrix<double, nsd, nsd>& A, const Eigen::Matrix<double, nsd, nsd>& B) {
        return A.cwiseProduct(B).sum();
    }
    
    double mat_det(const Array<double>& A, const int nd);
    Array<double> mat_dev(const Array<double>& A, const int nd);

    Array<double> mat_dyad_prod(const Vector<double>& u, const Vector<double>& v, const int nd);

    template <int nsd>
    Eigen::Matrix<double, nsd, nsd> mat_dyad_prod_eigen(const Eigen::Matrix<double, nsd, 1>& u, const Eigen::Matrix<double, nsd, 1>& v) {
        return u * v.transpose();
    }

    Array<double> mat_id(const int nsd);
    Array<double> mat_inv(const Array<double>& A, const int nd, bool debug = false);
    Array<double> mat_inv_ge(const Array<double>& A, const int nd, bool debug = false);
    Array<double> mat_inv_lp(const Array<double>& A, const int nd);

    Vector<double> mat_mul(const Array<double>& A, const Vector<double>& v);
    Array<double> mat_mul(const Array<double>& A, const Array<double>& B);
    void mat_mul(const Array<double>& A, const Array<double>& B, Array<double>& result);
    void mat_mul6x3(const Array<double>& A, const Array<double>& B, Array<double>& C);

    Array<double> mat_symm(const Array<double>& A, const int nd);
    Array<double> mat_symm_prod(const Vector<double>& u, const Vector<double>& v, const int nd);

    template <int nsd>
    Eigen::Matrix<double, nsd, nsd> mat_symm_prod_eigen(const Eigen::Matrix<double, nsd, 1>& u, const Eigen::Matrix<double, nsd, 1>& v) {
        return 0.5 * (u * v.transpose() + v * u.transpose());
    }

    double mat_trace(const Array<double>& A, const int nd);

    Tensor4<double> ten_asym_prod12(const Array<double>& A, const Array<double>& B, const int nd);
    Tensor4<double> ten_ddot(const Tensor4<double>& A, const Tensor4<double>& B, const int nd);
    Tensor4<double> ten_ddot_2412(const Tensor4<double>& A, const Tensor4<double>& B, const int nd);
    Tensor4<double> ten_ddot_3424(const Tensor4<double>& A, const Tensor4<double>& B, const int nd);

    /**
     * @brief Contracts two 4th order tensors A and B, C_ijkl = A_ijmn * B_klmn over the last two dimensions
     * 
     * @tparam nsd, the number of spatial dimensions
     * @param A, the first 4th order tensor
     * @param B, the second 4th order tensor
     * @return Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>> 
     */
    template <int nsd>
    Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>>
    ten_ddot_eigen(const Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>>& A, 
                    const Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>>& B) 
    {
        // Define the contraction dimensions
        Eigen::array<Eigen::IndexPair<int>, 2> contractionDims = {
            Eigen::IndexPair<int>(2, 2), // Contract A's 2nd dimension with B's 2nd
            Eigen::IndexPair<int>(3, 3)  // Contract A's 3rd dimension with B's 3rd
        };

        // Return the double dot product
        return A.contract(B, contractionDims);

        // For some reason, the Eigen::Tensor contract function above is faster than the manual implementation below

        // // Initialize the result tensor
        // Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>> C;

        // // Compute the double dot product: C_ijkl = A_ijmn * B_klmn
        // for (int i = 0; i < nsd; ++i) {
        //     for (int j = 0; j < nsd; ++j) {
        //         for (int k = 0; k < nsd; ++k) {
        //             for (int l = 0; l < nsd; ++l) {
        //                 C(i,j,k,l) = 0.0;
        //                 for (int m = 0; m < nsd; ++m) {
        //                     for (int n = 0; n < nsd; ++n) {
        //                         C(i,j,k,l) += A(i,j,m,n) * B(k,l,m,n);
        //                     }
        //                 }
        //             }
        //         }
        //     }
        // }

        // return C;
    }

    /**
     * @brief Contracts two 4th order tensors A and B, C_ijkl = A_imjn * B_mnkl over the 1st and 3rd dimensions of A and 0th and 2nd dimensions of B
     * 
     * @tparam nsd, the number of spatial dimensions
     * @param A, the first 4th order tensor
     * @param B, the second 4th order tensor
     * @return Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>> 
     */
    template <int nsd>
    Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>>
    ten_ddot_2412_eigen(const Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>>& A, 
                         const Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>>& B) 
    {
        // Define the contraction dimensions
        Eigen::array<Eigen::IndexPair<int>, 2> contractionDims = {
            Eigen::IndexPair<int>(1, 0), // Contract A's 1st dimension with B's 0th
            Eigen::IndexPair<int>(3, 1)  // Contract A's 3rd dimension with B's 1st
        };

        // Return the double dot product
        return A.contract(B, contractionDims);

        // // Initialize the result tensor
        // Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>> C;

        // // Compute the double dot product: C_ijkl = A_imjn * B_mnkl
        // for (int i = 0; i < nsd; ++i) {
        //     for (int j = 0; j < nsd; ++j) {
        //         for (int k = 0; k < nsd; ++k) {
        //             for (int l = 0; l < nsd; ++l) {
        //                 C(i,j,k,l) = 0.0;
        //                 for (int m = 0; m < nsd; ++m) {
        //                     for (int n = 0; n < nsd; ++n) {
        //                         C(i,j,k,l) += A(i,m,j,n) * B(m,n,k,l);
        //                     }
        //                 }
        //             }
        //         }
        //     }
        // }

        // return C;
    }

    /**
     * @brief Contracts two 4th order tensors A and B, C_ijkl = A_ijmn * B_kmln over the 2nd and 3rd dimensions of A and 1st and 3rd dimensions of B
     * 
     * @tparam nsd, the number of spatial dimensions
     * @param A, the first 4th order tensor
     * @param B, the second 4th order tensor
     * @return Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>> 
     */
    template <int nsd>
    Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>>
    ten_ddot_3424_eigen(const Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>>& A, 
                         const Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>>& B) 
    {
        // Define the contraction dimensions
        Eigen::array<Eigen::IndexPair<int>, 2> contractionDims = {
            Eigen::IndexPair<int>(2, 1), // Contract A's 2nd dimension with B's 1st
            Eigen::IndexPair<int>(3, 3)  // Contract A's 3rd dimension with B's 3rd
        };

        // Return the double dot product
        return A.contract(B, contractionDims);

        // // Initialize the result tensor
        // Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>> C;

        // // Compute the double dot product: C_ijkl = A_ijmn * B_kmln
        // for (int i = 0; i < nsd; ++i) {
        //     for (int j = 0; j < nsd; ++j) {
        //         for (int k = 0; k < nsd; ++k) {
        //             for (int l = 0; l < nsd; ++l) {
        //                 C(i,j,k,l) = 0.0;
        //                 for (int m = 0; m < nsd; ++m) {
        //                     for (int n = 0; n < nsd; ++n) {
        //                         C(i,j,k,l) += A(i,j,m,n) * B(k,m,l,n);
        //                     }
        //                 }
        //             }
        //         }
        //     }
        // }

        // return C;
    }

    Tensor4<double> ten_dyad_prod(const Array<double>& A, const Array<double>& B, const int nd);
    
    /**
     * @brief Compute the dyadic product of two 2nd order tensors A and B, C_ijkl = A_ij * B_kl
     * 
     * @tparam nsd, the number of spatial dimensions
     * @param A, the first 2nd order tensor
     * @param B, the second 2nd order tensor
     * @return Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>> 
     */
    template <int nsd>
    Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>> 
    ten_dyad_prod_eigen(const Eigen::Matrix<double, nsd, nsd>& A, const Eigen::Matrix<double, nsd, nsd>& B) {
       
        // // Create empty contraction dimensions
        // Eigen::array<Eigen::IndexPair<int>, 0> contractionDims = {};

        // // Convert the input matrices to Eigen::TensorFixedSize so we can use the contract method
        // Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd>> A_tensor = Eigen::TensorMap<const Eigen::Tensor<const double, 2>>(A.data(), nsd, nsd);
        // Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd>> B_tensor = Eigen::TensorMap<const Eigen::Tensor<const double, 2>>(B.data(), nsd, nsd);

        // // The contraction of A and B with no dimensions specified is the outer product
        // // of the two matrices, C_ijkl = A_ij * B_kl
        // return A_tensor.contract(B_tensor, contractionDims);

        // For some reason, the Eigen::Tensor contract function above is slower than the manual implementation below
        
        // Initialize the result tensor
        Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>> C;

        // Compute the dyadic product: C_ijkl = A_ij * B_kl
        for (int i = 0; i < nsd; ++i) {
            for (int j = 0; j < nsd; ++j) {
                for (int k = 0; k < nsd; ++k) {
                    for (int l = 0; l < nsd; ++l) {
                        C(i,j,k,l) = A(i,j) * B(k,l);
                    }
                }
            }
        }

        return C;
    }

    Tensor4<double> ten_ids(const int nd);

    /**
     * @brief Create a 4th order identity tensor:
     * I_ijkl = 0.5 * (δ_ik * δ_jl + δ_il * δ_jk)
     * 
     * @tparam nsd, the number of spatial dimensions
     * @return Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>> 
     */
    template <int nsd>
    Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>> 
    ten_ids_eigen() {
        // Initialize as zero
        Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>> I;
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
    /// Reproduces 'FUNCTION TEN_SYMMPROD(A, B, nd) RESULT(C)'.
    //
    template <int nsd>
    Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>>
    ten_symm_prod_eigen(const Eigen::Matrix<double, nsd, nsd>& A, const Eigen::Matrix<double, nsd, nsd>& B) {
        
        // // Create empty contraction dimensions
        // Eigen::array<Eigen::IndexPair<int>, 0> contractionDims = {};

        // // Convert the input matrices to Eigen::TensorFixedSize so we can use the contract method
        // Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd>> A_tensor = Eigen::TensorMap<const Eigen::Tensor<const double, 2>>(A.data(), nsd, nsd);
        // Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd>> B_tensor = Eigen::TensorMap<const Eigen::Tensor<const double, 2>>(B.data(), nsd, nsd);


        // // Compute the outer product of A and B, A_ij * B_kl
        // // The contraction of A and B with no dimensions specified is the outer product
        // Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>> AB = A_tensor.contract(B_tensor, contractionDims); // A_{ij} * B_{kl}

        // // Compute the product term1 = A_ik * B_jl and term2 = A_il * B_jk by shuffling the indices
        // Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>> term1, term2;
        // term1.shuffle(Eigen::array<int, 4>{0, 2, 1, 3}) = AB; // A_{ik} * B_{jl}
        // term2.shuffle(Eigen::array<int, 4>{0, 3, 1, 2}) = AB; // A_{il} * B_{jk}

        // // Compute the symmetric product: C_{ijkl} = 0.5 * (A_{ik} * B_{jl} + A_{il} * B_{jk})
        // Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>> C = 0.5 * (term1 + term2);

        // For some reason, the Eigen::Tensor contract function above is slower than the manual implementation below
        // Initialize the result tensor
        Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>> C;

        // Compute the symmetric product: C_ijkl = 0.5 * (A_ik * B_jl + A_il * B_jk)
        for (int i = 0; i < nsd; ++i) {
            for (int j = 0; j < nsd; ++j) {
                for (int k = 0; k < nsd; ++k) {
                    for (int l = 0; l < nsd; ++l) {
                        C(i,j,k,l) = 0.5 * (A(i,k) * B(j,l) + A(i,l) * B(j,k));
                    }
                }
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
     * @return Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>> 
     */
    template <int nsd>
    Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>>
    ten_transpose_eigen(const Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>>& A) {

        // Initialize the result tensor
        Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>> B;

        // // Permute the tensor indices to perform the transpose operation
        // B.shuffle(Eigen::array<int, 4>{2, 3, 0, 1}) = A;

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

