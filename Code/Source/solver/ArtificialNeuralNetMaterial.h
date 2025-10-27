// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

/* This material model implementation is based on the following paper: 
Peirlinck, M., Hurtado, J.A., Rausch, M.K. et al. A universal material model subroutine 
for soft matter systems. Engineering with Computers 41, 905–927 (2025). 
https://doi.org/10.1007/s00366-024-02031-w */

#ifndef ArtificialNeuralNet_model_H
#define ArtificialNeuralNet_model_H

#include "mat_fun.h"
#include "utils.h"
#include "Parameters.h"
#include <vector>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"
#include "eigen3/unsupported/Eigen/CXX11/Tensor"

using namespace mat_fun;

// Class for parameter table for material models discovered by constitutive artificial neural network (CANN)

/* This material model implementation is based on the following paper: 
Peirlinck, M., Hurtado, J.A., Rausch, M.K. et al. A universal material model subroutine 
for soft matter systems. Engineering with Computers 41, 905–927 (2025). 
https://doi.org/10.1007/s00366-024-02031-w */

class ArtificialNeuralNetMaterial
{
  public:

    // Invariant indices
    Vector<int> invariant_indices;

    // Activation functions
    Array<int> activation_functions;

    // Weights
    Array<double> weights;

    // Number of rows in parameter table
    int num_rows;

    // Outputs from each layer
    void uCANN_h0(const double x, const int kf, double &f, double &df, double &ddf) const;
    void uCANN_h1(const double x, const int kf, const double W, double &f, double &df, double &ddf) const;
    void uCANN_h2(const double x, const int kf, const double W, double &f, double &df, double &ddf) const;

    // Strain energy and derivatives
    void uCANN(const double xInv, const int kInv,
           const int kf0, const int kf1, const int kf2,
           const double W0, const double W1, const double W2,
           double &psi, double (&dpsi)[9], double (&ddpsi)[9]) const;


    void evaluate(const double aInv[9], double &psi, double (&dpsi)[9], double (&ddpsi)[9]) const;

    // Helper for compute_pk2cc
    template<size_t nsd>
    void computeInvariantsAndDerivatives(
    const Matrix<nsd>& C, const Matrix<nsd>& fl, int nfd, double J2d, double J4d, const Matrix<nsd>& Ci,
    const Matrix<nsd>& Idm, const double Tfa, Matrix<nsd>& N1, double& psi, double (&Inv)[9], std::array<Matrix<nsd>,9>& dInv,
    std::array<Tensor<nsd>,9>& ddInv) const; 
    
};

#endif // ArtificialNeuralNet_model_H