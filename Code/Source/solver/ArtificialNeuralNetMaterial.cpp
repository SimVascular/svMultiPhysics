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

/* This material model implementation is based on the following paper: 
Peirlinck, M., Hurtado, J.A., Rausch, M.K. et al. A universal material model subroutine 
for soft matter systems. Engineering with Computers 41, 905–927 (2025). 
https://doi.org/10.1007/s00366-024-02031-w */

// The functions required for CANN material model implementations are defined here

#include "ArtificialNeuralNetMaterial.h"
#include "ComMod.h"
#include "mat_fun.h"
using namespace mat_fun;

/// @brief 0th layer output of CANN for activation func kf, input x
void ArtificialNeuralNetMaterial::uCANN_h0(const double x, const int kf, double &f, double &df, double &ddf) const {
    if (kf == 1) {
        f = x;
        df = 1;
        ddf = 0;
    } else if (kf == 2) {
        f = (std::abs(x) + x) / 2;
        df = 0.5 * (std::abs(x) / x + 1);
        ddf = 0;
    } else if (kf == 3) {
        f = std::abs(x);
        df = std::abs(x) / x;
        ddf = 0;
    }
}

/// @brief 1st layer output of CANN for activation func kf, input x, weight W
void ArtificialNeuralNetMaterial::uCANN_h1(const double x, const int kf, const double W, double &f, double &df, double &ddf) const {
    if (kf == 1) {
        f = W * x;
        df = W;
        ddf = 0;
    } else if (kf == 2) {
        f = W * W * x * x;
        df = 2 * W * W * x;
        ddf = 2 * W * W;
    }
}

/// @brief 2nd layer output of CANN for activation func kf, input x, weight W
void ArtificialNeuralNetMaterial::uCANN_h2(const double x, const int kf, const double W, double &f, double &df, double &ddf) const {
    if (kf == 1) {
        f = W * x;
        df = W;
        ddf = 0;
    } else if (kf == 2) {
        f = std::exp(W * x) - 1;
        df = W * std::exp(W * x);
        ddf = W * W * std::exp(W * x);
    } else if (kf == 3) {
        f = -std::log(1 - W * x);
        df = W / (1 - W * x);
        ddf = -W * W / ((1 - W * x) * (1 - W * x));
    }
}

/// @brief Updates psi and its derivatives
void ArtificialNeuralNetMaterial::uCANN(
    const double xInv, const int kInv,
    const int kf0, const int kf1, const int kf2,
    const double W0, const double W1, const double W2,
    double &psi, double (&dpsi)[9], double (&ddpsi)[9]
) const {
    double f0, df0, ddf0;
    uCANN_h0(xInv, kf0, f0, df0, ddf0);
    double f1, df1, ddf1;
    uCANN_h1(f0, kf1, W0, f1, df1, ddf1);
    double f2, df2, ddf2;
    uCANN_h2(f1, kf2, W1, f2, df2, ddf2);

    psi += W2 * f2;
    dpsi[kInv - 1] += W2 * df2 * df1 * df0;
    ddpsi[kInv - 1] += W2 * ((ddf2 * df1 * df1 + df2 * ddf1) * df0 * df0 + df2 * df1 * ddf0);
}

/// @brief function to build psi and dpsidI1 to 9
void ArtificialNeuralNetMaterial::evaluate(const double aInv[9], double &psi, double (&dpsi)[9], double (&ddpsi)[9]) const {
    // Initializing
    psi = 0;
    for (int i = 0; i < 9; ++i) {
        dpsi[i] = 0;
        ddpsi[i] = 0;
    }

    double ref[9] = {3, 3, 1, 1, 1, 0, 0, 1, 1};

    for (int i = 0; i < num_rows; ++i) {
        int kInv = this->invariant_indices(i);
        int kf0 = this->activation_functions(i, 0);
        int kf1 = this->activation_functions(i, 1);
        int kf2 = this->activation_functions(i, 2);
        double W0 = this->weights(i, 0);
        double W1 = this->weights(i, 1);
        double W2 = this->weights(i, 2);

        double xInv = aInv[kInv - 1] - ref[kInv - 1];
        uCANN(xInv, kInv, kf0, kf1, kf2, W0, W1, W2, psi, dpsi, ddpsi);
    }
}

template<size_t nsd>
void ArtificialNeuralNetMaterial::computeInvariantsAndDerivatives(
const Matrix<nsd>& C, const Matrix<nsd>& fl, int nfd, double J2d, double J4d, const Matrix<nsd>& Ci,
const Matrix<nsd>& Idm, const double Tfa, Matrix<nsd>& N1, double& psi, double (&Inv)[9],std::array<Matrix<nsd>,9>& dInv,
std::array<Tensor<nsd>,9>& ddInv) const {

    // double Inv[9] = {0};
    Matrix<nsd> C2 = C * C;

    Inv[0] = J2d * C.trace();
    Inv[1] = 0.50 * (Inv[0]*Inv[0] - J4d * (C*C).trace());
    Inv[2] = C.determinant();
    Inv[3] = J2d * (fl.col(0).dot(C * fl.col(0)));
    Inv[4] = J4d * (fl.col(0).dot(C2 * fl.col(0)));

    Matrix<nsd> dInv1 = -Inv[0]/3 * Ci + J2d * Idm;
    Matrix<nsd> dInv2 = (C2.trace()/3)*Ci + Inv[0]*dInv1 + J4d*C;
    Matrix<nsd> dInv3 = Inv[2]*Ci;
    N1 = fl.col(0)*fl.col(0).transpose();
    Matrix<nsd> dInv4 = -Inv[3]/3*Ci + J2d*N1;
    Matrix<nsd> dInv5 = J4d*(N1*C + C*N1) - Inv[4]/3*Ci;

    Matrix<nsd> dInv6, dInv7, dInv8, dInv9;
    Tensor<nsd> ddInv6, ddInv7, ddInv8, ddInv9;
    // initialize to 0
    dInv6.setZero();
    dInv7.setZero();
    dInv8.setZero();
    dInv9.setZero();
    ddInv6.setZero();
    ddInv7.setZero();
    ddInv8.setZero();
    ddInv9.setZero();

    Tensor<nsd> dCidC = -symmetric_dyadic_product<nsd>(Ci, Ci);
    Matrix<nsd> dJ4ddC = -2.0/3.0 * J4d * Ci;

    Tensor<nsd> ddInv1 = (-1.0/3.0)*(dyadic_product<nsd>(dInv1,Ci) + Inv[0]*dCidC + J2d*dyadic_product(Ci,Idm));
    Tensor<nsd> ddInv2 = dyadic_product<nsd>(dInv1,dInv1) + Inv[0]*ddInv1 + (1.0/3.0)*C2.trace()*dCidC 
                        + (1.0/3.0)*dyadic_product<nsd>((C2.trace()*dJ4ddC + 2*J4d*C),Ci) 
                        + dyadic_product<nsd>(dJ4ddC,C) - J4d*fourth_order_identity<nsd>();
    Tensor<nsd> ddInv3 = dyadic_product<nsd>(dInv3,Ci) + Inv[2]*dCidC;
    Tensor<nsd> ddInv4 = (-1.0/3.0)*(dyadic_product<nsd>(dInv4,Ci) + J2d*dyadic_product<nsd>(Ci,N1) + Inv[3]*dCidC);
    Matrix<nsd> sum1 = (N1*C + C*N1);
    Tensor<nsd> ddInv5 = (-1.0/3.0)*(dyadic_product<nsd>(dInv5,Ci) + Inv[4]*dCidC + 2*J4d*dyadic_product(Ci,sum1))
                        + J4d*(2*symmetric_dyadic_product(N1,Idm) - dyadic_product(N1,Idm)
                        + 2*symmetric_dyadic_product(Idm,N1) - dyadic_product(Idm,N1));

    if (nfd == 2) {
        Inv[5] = J2d * (fl.col(0).dot(C * fl.col(1)));
        Inv[6] = J4d * (fl.col(0).dot(C2 * fl.col(1)));
        Inv[7] = J2d * (fl.col(1).dot(C * fl.col(1)));
        Inv[8] = J4d * (fl.col(1).dot(C2 * fl.col(1)));

        Matrix<nsd> N2 = fl.col(1)*fl.col(1).transpose();
        Matrix<nsd> N12 = 0.5*(fl.col(0)*fl.col(1).transpose() + fl.col(1)*fl.col(0).transpose());
        

        dInv6 = -Inv[5]/3*Ci + J2d*N12;
        dInv7 = J4d*(N12*C + C*N12) - Inv[6]/3*Ci;
        dInv8 = -Inv[7]/3*Ci + J2d*N2;
        dInv9 = J4d*(N2*C + C*N2) - Inv[8]/3*Ci;

        ddInv6 = -1.0/3.0*(dyadic_product(dInv6,Ci) + J2d*dyadic_product(Ci,N12) + Inv[5]*dCidC);
        Matrix<nsd> sum12 = (N12*C + C*N12);
        ddInv7 = -1.0/3.0*(dyadic_product(dInv7,Ci) + Inv[6]*dCidC + 2*J4d*dyadic_product(Ci,sum12))
                + J4d*(2*symmetric_dyadic_product(N12,Idm) - dyadic_product(N12,Idm)
                + 2*symmetric_dyadic_product(Idm,N12) - dyadic_product(Idm,N12));
        ddInv8 = -1.0/3.0*(dyadic_product(dInv8,Ci) + J2d*dyadic_product(Ci,N2) + Inv[7]*dCidC);
        Matrix<nsd> sum2 = (N2*C + C*N2);
        ddInv9 = -1.0/3.0*(dyadic_product(dInv9,Ci) + Inv[8]*dCidC + 2*J4d*dyadic_product(Ci,sum2))
                + J4d*(2*symmetric_dyadic_product(N2,Idm) - dyadic_product(N2,Idm)
                + 2*symmetric_dyadic_product(Idm,N2) - dyadic_product(Idm,N2));
    }

    dInv = {dInv1, dInv2, dInv3, dInv4, dInv5, dInv6, dInv7, dInv8, dInv9};
    ddInv = {ddInv1, ddInv2, ddInv3, ddInv4, ddInv5, ddInv6, ddInv7, ddInv8, ddInv9};
}


// Template instantiation
template void ArtificialNeuralNetMaterial::computeInvariantsAndDerivatives<2>(
const Matrix<2>& C, const Matrix<2>& fl, int nfd, double J2d, double J4d, const Matrix<2>& Ci,
const Matrix<2>& Idm, const double Tfa, Matrix<2>& N1, double& psi, double (&Inv)[9], std::array<Matrix<2>,9>& dInv,
std::array<Tensor<2>,9>& ddInv) const;

template void ArtificialNeuralNetMaterial::computeInvariantsAndDerivatives<3>(
const Matrix<3>& C, const Matrix<3>& fl, int nfd, double J2d, double J4d, const Matrix<3>& Ci,
const Matrix<3>& Idm, const double Tfa, Matrix<3>& N1, double& psi, double (&Inv)[9], std::array<Matrix<3>,9>& dInv,
std::array<Tensor<3>,9>& ddInv) const;
