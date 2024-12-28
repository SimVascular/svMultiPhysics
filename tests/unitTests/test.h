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

// --------------------------------------------------------------
// To run the tests in test.cpp
// 1.  Navigate to <svMultiPhysics_root_directory>/build/svMultiPhysics-build/Source/solver
// 2.  Run `make` to build the tests
// 3.  Run `ctest --verbose` to run the tests


// --------------------------------------------------------------
// To add a new material model to test:
// In test.h:
// 1. Create a new material parameters class derived from MatParams (e.g. NeoHookeanParams)
// 2. Create a new material model test class derived from TestMaterialModel (e.g. TestNeoHookean)
// 3. Implement required functions in the material model test class:
//    - Constructor: Sets the material model type and parameters for svMultiPhysics (e.g. TestNeoHookean())
//    - printMaterialParameters(): Prints the material parameters
//    - computeStrainEnergy(): Computes the strain energy density function
// In test.cpp:
// 4. Create a new text fixture class derived from ::testing::Test (e.g. NeoHookeanTest)
//    - In this you set the values of the material parameters for testing
// 5. Add tests for your new material model (e.g. TEST_F(NeoHookeanTest, TestPK2StressIdentityF))
// --------------------------------------------------------------


#include <stdlib.h>
#include <iostream>
#include <random>
#include <chrono>
#include "gtest/gtest.h"   // include GoogleTest
#include "mat_fun.h"
#include "mat_models.h"

// --------------------------------------------------------------
// ---------------------- Helper functions ----------------------
// --------------------------------------------------------------

/**
 * @brief Creates an identity deformation gradient F.
 *
 * @param[out] F The deformation gradient tensor to be set to the identity matrix.
 * @return The deformation gradient tensor F set to the identity matrix.
 */
Array<double> create_identity_F(const int N) {
    Array<double> F(N, N);
    for (int i = 0; i < N; i++) {
        for (int J = 0; J < N; J++) {
            F(i, J) = (i == J);
        }
    }
    return F;
}

/**
 * @brief Create a ones matrix.
 * 
 */
Array<double> create_ones_matrix(const int N) {
    Array<double> A(N, N);
    for (int i = 0; i < N; i++) {
        for (int J = 0; J < N; J++) {
            A(i, J) = 1.0;
        }
    }
    return A;
}

/**
 * @brief Generates a random double value.
 *
 * This function generates a random double value within a specified range.
 *
 * @param[in] min The minimum value of the range.
 * @param[in] max The maximum value of the range.
 * @return A random double value between min and max.
 */
inline double getRandomDouble(double min, double max) {
    // Uncomment to use a random seed
    //unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    unsigned int seed = 42;
    static std::default_random_engine engine(seed);
    std::uniform_real_distribution<double> distribution(min, max);
    return distribution(engine);
}

/**
 * @brief Creates a random deformation gradient F with values between min and max, and det(F) > 0.
 *
 * This function generates a random deformation gradient tensor F such that the determinant of F is greater than 0.
 * @param[in] N The size of the deformation gradient tensor (NxN).
 * @param[in] min The minimum value for the elements of the deformation gradient tensor (default is 0.1).
 * @param[in] max The maximum value for the elements of the deformation gradient tensor (default is 10.0).
 * @return A random deformation gradient tensor F.
 */
Array<double> create_random_F(const int N, double min=0.1, double max=10.0) {
    // Create a random deformation gradient with values between min and max, 
    // and det(F) > 0
    Array<double> F(N, N);
    double J = -1.0;
    while (J < 0) {
        for (int i = 0; i < N; i++) {
            for (int J = 0; J < N; J++) {
                F(i,J) = getRandomDouble(min, max);
            }
        }
        J = mat_fun::mat_det(F, N);
    }
    return F;
}

/**
 * @brief Creates a deformation matrix F of random deviations from the identity matrix.
 *
 * This function generates a random deformation gradient tensor F with values perturbed from the identity matrix,
 * such that the determinant of F is greater than 0.
 *
 * @param[in] N The size of the deformation gradient tensor (NxN).
 * @param[in] max_deviation The maximum deviation from the identity matrix elements.
 * @return A random deformation gradient tensor F.
 */
Array<double> create_random_perturbed_identity_F(const int N, double max_deviation) {
    // Create a random deformation gradient with values perturbed from the identity matrix, 
    // and det(F) > 0
    Array<double> F(N, N);
    double J = -1.0;
    while (J < 0) {
        for (int i = 0; i < N; i++) {
            for (int J = 0; J < N; J++) {
                F(i,J) = (i == J) + max_deviation * getRandomDouble(-1.0, 1.0);
            }
        }
        J = mat_fun::mat_det(F, N);
    }
    return F;
}

/**
 * @brief Perturbs the deformation gradient F by delta times a random number between -1 and 1.
 *
 * This function perturbs the given deformation gradient tensor F by adding delta times a random number 
 * between -1 and 1 to each element, and stores the perturbed deformation gradient in F_tilde.
 *
 * @param[in] F The original deformation gradient tensor.
 * @param[in] delta The perturbation factor.
 * @param[out] F_tilde The perturbed deformation gradient tensor.
 * @param[in] N The size of the deformation gradient tensor (NxN).
 * @return None.
 */
Array<double> perturb_random_F(const Array<double> &F, const double delta) {

    int N = F.nrows(); // Size of the deformation gradient tensor
    assert (N == F.ncols()); // Check that F is square
    
    // Perturb the deformation gradient and store in F_tilde
    Array<double> F_tilde(N, N);
    double dF_iJ;
    for (int i = 0; i < N; i++) {
        for (int J = 0; J < N; J++) {
            dF_iJ = delta * getRandomDouble(-1.0, 1.0);
            F_tilde(i,J) = F(i,J) + dF_iJ; // perturbed deformation gradient
        }
    }
    return F_tilde;
}

/**
 * @brief Computes the Jacobian J, right Cauchy-Green deformation tensor C, and Green-Lagrange strain tensor E from the deformation gradient F.
 *
 * This function computes the Jacobian of the deformation gradient tensor F, the right Cauchy-Green deformation tensor C, 
 * and the Green-Lagrange strain tensor E.
 *
 * @param[in] F The deformation gradient tensor.
 * @param[out] J The computed Jacobian of F.
 * @param[out] C The computed right Cauchy-Green deformation tensor.
 * @param[out] E The computed Green-Lagrange strain tensor.
 * @return None.
 */
void calc_JCE(const Array<double> &F, double &J, Array<double> &C, Array<double> &E) {

    int N = F.nrows(); // Size of the deformation gradient tensor
    assert (N == F.ncols()); // Check that F is square

    // Compute Jacobian of F
    J = mat_fun::mat_det(F, N);

    // Compute transpose of F
    auto F_T = mat_fun::transpose(F);

    // Compute right Cauchy-Green deformation tensor
    C = mat_fun::mat_mul(F_T, F);

    // Compute Green-Lagrange strain tensor
    E = 0.5 * (C - mat_fun::mat_id(N));
}

/**
 * @brief Structure to store solid mechanics terms used to compute strain energy density functions.
 *
 * @tparam N The size of the deformation gradient tensor (NxN).
 */
struct solidMechanicsTerms {
    double J;           /**< Jacobian of the deformation gradient tensor. */
    Array<double> C;    /**< Right Cauchy-Green deformation tensor. */
    Array<double> E;    /**< Green-Lagrange strain tensor. */
    Array<double> E2;   /**< Second-order Green-Lagrange strain tensor. */
    Array<double> C_bar;/**< Modified right Cauchy-Green deformation tensor. */
    double I1;          /**< First invariant of the right Cauchy-Green deformation tensor. */
    double I2;          /**< Second invariant of the right Cauchy-Green deformation tensor. */
    double Ib1;         /**< First invariant of the modified right Cauchy-Green deformation tensor. */
    double Ib2;         /**< Second invariant of the modified right Cauchy-Green deformation tensor. */
};

/**
 * @brief Computes the solid mechanics terms used to compute strain energy density functions.
 *
 * This function computes various solid mechanics terms such as the Jacobian, right Cauchy-Green deformation tensor,
 * Green-Lagrange strain tensor, and their invariants from the given deformation gradient tensor F.
 *
 * @param[in] F The deformation gradient tensor.
 * @return A structure containing the computed solid mechanics terms.
 */
solidMechanicsTerms calcSolidMechanicsTerms(const Array<double> &F) {

    int N = F.nrows(); // Size of the deformation gradient tensor
    assert (N == F.ncols()); // Check that F is square

    solidMechanicsTerms out;

    const double N_d = static_cast<double>(N); // Convert N to double for calculations

    // Jacobian of F
    out.J = mat_fun::mat_det(F, N);

    // Transpose of F
    auto F_T = mat_fun::transpose(F);

    // Right Cauchy-Green deformation tensor
    out.C = mat_fun::mat_mul(F_T, F);

    // Right Cauchy-Green deformation tensor squared
    auto C2 = mat_fun::mat_mul(out.C, out.C);

    // Green-Lagrange strain tensor
    out.E = 0.5 * (out.C - mat_fun::mat_id(N));

    // Green-Lagrange strain tensor squared
    out.E2 = mat_fun::mat_mul(out.E, out.E);

    // Modified right Cauchy-Green deformation tensor
    out.C_bar = pow(out.J, (-2.0/N_d)) * out.C;

    // Modified right Cauchy-Green deformation tensor squared
    auto C_bar2 = mat_fun::mat_mul(out.C_bar, out.C_bar);

    // Invariants of C
    out.I1 = mat_fun::mat_trace(out.C, N);
    out.I2 = 0.5 * (pow(out.I1, 2) - mat_fun::mat_trace(C2, N));

    // Invariants of C_bar
    out.Ib1 = mat_fun::mat_trace(out.C_bar, N);
    out.Ib2 = 0.5 * (pow(out.Ib1, 2) - mat_fun::mat_trace(C_bar2, N));

    // Check that invariants satisfy expected relationship
    EXPECT_NEAR( pow(out.J, (-2.0/3.0)) * out.I1, out.Ib1, 1e-9 * out.Ib1);
    EXPECT_NEAR( pow(out.J, (-4.0/3.0)) * out.I2, out.Ib2, 1e-9 * out.Ib2);
    
    return out;
}

/**
 * @brief Computes a linear regression line y = mx + b for given x and y data.
 * 
 * @param x x data points.
 * @param y y data points.
 * @return std::pair<double, double> A pair containing the slope (m) and the y-intercept (b).
 */
std::pair<double, double> computeLinearRegression(const std::vector<double>& x, const std::vector<double>& y) {
    int n = x.size();
    double sum_x = std::accumulate(x.begin(), x.end(), 0.0);
    double sum_y = std::accumulate(y.begin(), y.end(), 0.0);
    double sum_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
    double sum_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);

    double m = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x * sum_x);
    double b = (sum_y - m * sum_x) / n;

    return std::make_pair(m, b);
}

// --------------------------------------------------------------
// -------------------- Mock svMultiPhysics object -------------------
// --------------------------------------------------------------


class MockCepMod : public CepMod {
public:
    MockCepMod() {
        // initialize if needed 
    }
    // Mock methods if needed
};
class MockdmnType : public dmnType {
public:
    MockdmnType() {
        // initialize if needed 
    }
    // MockstModelType mockStM;
    // Mock methods if needed
};
class MockmshType : public mshType {
public:
    MockmshType() {
        // initialize if needed 
    }
    // Mock methods if needed
};
class MockeqType : public eqType {
public:
    MockeqType() {
        // initialize if needed 
    }
    MockdmnType mockDmn;
    // Mock methods if needed
};
class MockComMod : public ComMod {
public:
    MockComMod() {
        // initialize if needed 
        nsd = 3;
    }
    MockeqType mockEq;
    MockmshType mockMsh;
    // Mock methods if needed
};


// --------------------------------------------------------------
// ------------------ Test Material Model Classes ---------------
// --------------------------------------------------------------

// Class for testing material models in svMultiPhysics
class TestMaterialModel {
public:
    MockComMod com_mod;
    MockCepMod cep_mod;
    int nFn;
    Array<double> fN;
    double ya_g;
    bool ustruct;

    TestMaterialModel(const consts::ConstitutiveModelType matType, const consts::ConstitutiveModelType penType) {
        int nsd = com_mod.nsd;
        mat_fun::ten_init(nsd);                        // initialize tensor index pointer for mat_fun

        // Set material and penalty models
        auto &dmn = com_mod.mockEq.mockDmn;
        dmn.stM.isoType = matType;                            // Mat_model
        dmn.stM.volType = penType;                            // Dilational_penalty_model
        
        // Initialize fibers and other material parameters
        nFn = 2;                          // Number of fiber directions
        fN = Array<double>(nsd, nFn);     // Fiber directions array (initialized to zeros)
        ya_g = 0.0;                       // ?

        // Flag to use struct or ustruct material models
        // If struct, calls compute_pk2cc() and uses strain energy composed of isochoric and volumetric parts
        // If ustruct, calls compute_pk2cc() and uses strain energy composed of isochoric part only
        ustruct = false;

    // Material parameters are set in each derived class
    }

    // Pure virtual method to print material parameters
    virtual void printMaterialParameters() = 0;

    // Pure virtual method for computing Strain Energy
    virtual double computeStrainEnergy(const Array<double> &F) = 0;

    /**
     * @brief Computes the PK2 stress tensor S and material elasticity tensor Dm for a given deformation gradient F.
     *
     * This function computes the PK2 stress tensor S and the material elasticity tensor Dm from the deformation gradient F.
     * If `ustruct` is true, the deviatoric part of the PK2 stress tensor is returned using the `compute_pk2cc` function.
     *
     * @param[in] F The deformation gradient tensor.
     * @param[out] S The computed PK2 stress tensor.
     * @param[out] Dm The computed material elasticity tensor.
     * @return None, but fills S and Dm with the computed values.
     */
    void compute_pk2cc(const Array<double> &F, Array<double> &S,  Array<double> &Dm) {
        auto &dmn = com_mod.mockEq.mockDmn;

        double J = 0; // Jacobian (not used in this testing)
        
        if (ustruct) {
            dmn.phys = consts::EquationType::phys_ustruct;
        } else {
            dmn.phys = consts::EquationType::phys_struct;
        }

        // Call compute_pk2cc to compute S and Dm
        mat_models::compute_pk2cc(com_mod, cep_mod, dmn, F, nFn, fN, ya_g, S, Dm, J);

    }

       /**
     * @brief Computes the solid density, isothermal compressibility coefficient, and their derivatives for a given pressure.
     *
     * This function computes the solid density (rho), isothermal compressibility coefficient (beta), 
     * and their derivatives with respect to pressure (drho and dbeta) for a given pressure (p) using the g_vol_pen() function 
     * from mat_models.h.
     *
     * @param[in] p Pressure.
     * @param[in] rho0 Initial solid density.
     * @param[out] rho Computed solid density.
     * @param[out] beta Computed isothermal compressibility coefficient.
     * @param[out] drho Computed Derivative of solid density with respect to pressure.
     * @param[out] dbeta Computed Derivative of beta with respect to pressure.
     * @param[in] Ja Jacobian (not used in this function).
     * @return None.
     */
    void g_vol_pen(const double p, const double rho0, double &rho, double &beta, double &drho, double &dbeta, const double Ja) {
        auto &dmn = com_mod.mockEq.mockDmn;
        dmn.prop[consts::PhysicalProperyType::solid_density] = rho0; // Set initial solid density

        mat_models::g_vol_pen(com_mod, dmn, p, rho, beta, drho, dbeta, Ja);
    }

    /**
     * @brief Computes the PK2 stress tensor S(F) from the strain energy density Psi(F) using finite differences.
     *
     * Analytically, we should have S = dPsi/dE. Since we have Psi(F), we cannot directly compute S. 
     * Instead, we compute S = F^-1 * P, where P = dPsi/dF is computed using finite differences in each component of F.
     *
     * Pseudocode (for first order finite difference):
     * - Compute strain energy density Psi(F)
     * - For each component of F, F[i][J]
     *      - Perturb F[i][J] by delta to get F_tilde
     *      - Compute Psi(F_tilde)
     *      - Compute dPsi = Psi(F_tilde) - Psi(F)
     *      - Compute P[i][J] = dPsi / delta
     * - Compute S = F^-1 * P
     * 
     * @param[in] F The deformation gradient tensor.
     * @param[in] delta The perturbation scaling factor.
     * @param[in] order The order of the finite difference scheme (1 for first order, 2 for second order, etc.).
     * @param[out] S The computed PK2 stress tensor.
     * @param[in] N The size of the deformation gradient tensor (NxN).
     * @return None, but fills S with the computed values.
     */
    void calcPK2StressFiniteDifference(const Array<double> &F, const double delta, const int order, Array<double> & S) {
        
        int N = F.nrows(); // Size of the deformation gradient tensor
        assert(F.ncols() == N); // Check that F is square

        // Compute strain energy density given F
        double Psi = computeStrainEnergy(F);

        // Compute 1st PK stress P_iJ = dPsi / dF[i][J] using finite difference, component by component
        Array<double> P(N, N);
        if (order == 1){
            Array<double> F_tilde(N, N); // perturbed deformation gradient
            for (int i = 0; i < N; i++) {
                for (int J = 0; J < N; J++) {
                    // Perturb the iJ-th component of F by delta
                    for (int k = 0; k < N; k++) {
                        for (int l = 0; l < N; l++) {
                            F_tilde(k, l) = F(k, l);
                        }
                    }
                    F_tilde(i,J) += delta;

                    // Compute Psi_MR for perturbed deformation gradient
                    double Psi_tilde = computeStrainEnergy(F_tilde);

                    // Compute differences in Psi
                    double dPsi = Psi_tilde - Psi;

                    // Compute P(i,J) = dPsi / dF(i,J)
                    P(i, J) = dPsi / delta;
                }
            }
        }
        else if (order == 2){
            Array<double> F_plus(N, N); // positive perturbed deformation gradient
            Array<double> F_minus(N, N); // negative perturbed deformation gradient
            for (int i = 0; i < N; i++) {
                for (int J = 0; J < N; J++) {
                    // Perturb the iJ-th component of F by +-delta
                    for (int k = 0; k < N; k++) {
                        for (int l = 0; l < N; l++) {
                            F_plus(k,l) = F(k,l);
                            F_minus(k,l) = F(k,l);
                        }
                    }
                    F_plus(i,J) += delta;
                    F_minus(i,J) -= delta;

                    // Compute Psi_MR for perturbed deformation gradient
                    double Psi_plus = computeStrainEnergy(F_plus);
                    double Psi_minus = computeStrainEnergy(F_minus);

                    // Compute differences in Psi
                    double dPsi = Psi_plus - Psi_minus;

                    // Compute P(i,J) = dPsi / dF(i,J)
                    P(i,J) = dPsi / (2.0 * delta);
                }
            }
        }
        

        // Compute S_ref = F^-1 * P_ref
        auto F_inv = mat_fun::mat_inv(F, N);
        S = mat_fun::mat_mul(F_inv, P);
    }

    /**
     * @brief Computes the PK2 stress tensor S(F) from the strain energy density Psi(F) using finite differences and checks the order of convergence.
     * 
     * @param[in] F Deformation gradient.
     * @param[in] delta_min Minimum perturbation scaling factor.
     * @param[in] delta_max Maximum perturbation scaling factor.
     * @param[in] order Order of the finite difference scheme (1 for first order, 2 for second order, etc.).
     * @param[in] convergence_order_tol Tolerance for comparing convergence order with expected value
     * @param[in] verbose Show values error and order of convergence if true.
     */
    void testPK2StressConvergenceOrder(const Array<double> &F, const double delta_max, const double delta_min, const int order, const double convergence_order_tol, const bool verbose = false) {
        // Check delta_max > delta_min
        if (delta_max <= delta_min) {
            std::cerr << "Error: delta_max must be greater than delta_min." << std::endl;
            return;
        }

        // Check order is 1 or 2
        if (order != 1 && order != 2) {
            std::cerr << "Error: order must be 1 or 2." << std::endl;
            return;
        }

        int N = F.ncols();
        assert(F.nrows() == N);

        // Create list of deltas for convergence test (delta = delta_max, delta_max/2, delta_max/4, ...)
        std::vector<double> deltas;
        double delta = delta_max;
        while (delta >= delta_min) {
            deltas.push_back(delta);
            delta /= 2.0;
        }

        // Compute S(F) from compute_pk2cc()
        Array<double> S(3,3), Dm(6,6);
        compute_pk2cc(F, S, Dm);

        // Compute finite difference S for each delta and store error in list
        std::vector<double> errors;
        Array<double> S_fd(3,3);
        for (int i = 0; i < deltas.size(); i++) {
            calcPK2StressFiniteDifference(F, deltas[i], order, S_fd);

            // Compute Frobenius norm of error between S and S_fd
            double error = 0.0;
            for (int I = 0; I < 3; I++) {
                for (int J = 0; J < 3; J++) {
                    error += pow(S(I,J) - S_fd(I,J), 2);
                }
            }
            error = sqrt(error);

            // Store error in list
            errors.push_back(error);
        }

        // Compute order of convergence by fitting a line to log(delta) vs log(error)
        std::vector<double> log_deltas, log_errors;
        for (int i = 0; i < deltas.size(); i++) {
            log_deltas.push_back(log(deltas[i]));
            log_errors.push_back(log(errors[i]));
        }

        // Fit a line to log(delta) vs log(error)
        // m is the slope (order of convergence), b is the intercept
        auto [m, b] = computeLinearRegression(log_deltas, log_errors);

        // Check that order of convergence is > order - convergence_order_tol
        EXPECT_GT(m, order - convergence_order_tol);

        // Print results if verbose
        if (verbose) {
            std::cout << "Slope (order of convergence): " << m << std::endl;
            std::cout << "Intercept: " << b << std::endl;
            std::cout << "Errors: ";
            for (int i = 0; i < errors.size(); i++) {
                std::cout << errors[i] << " ";
            }
            std::cout << std::endl;
            std::cout << std::endl;
            
            std::cout << "F = " << std::endl;
            for (int i = 0; i < 3; i++) {
                for (int J = 0; J < 3; J++) {
                    std::cout << F(i,J) << " ";
                }
                std::cout << std::endl;
            }

            std::cout << "S = " << std::endl;
            for (int i = 0; i < 3; i++) {
                for (int J = 0; J < 3; J++) {
                    std::cout << S(i,J) << " ";
                }
                std::cout << std::endl;
            }

            std::cout << std::endl;
        }
    }

    /**
     * @brief Compute perturbation in strain energy density (dPsi) given perturbation in the deformation gradient (dF).
     * 
     * @param F Deformation gradient
     * @param dF Deformation gradient perturbation shape
     * @param delta Deformation gradient perturbation scaling factor
     * @param order Order of the finite difference scheme (1 for first order, 2 for second order, etc.)
     * @param dPsi Strain energy density perturbation
     */
    void calcdPsiFiniteDifference(const Array<double> &F, const Array<double> &dF, const double delta, const int order, double &dPsi) {

        int N = F.ncols();
        assert(F.nrows() == N);

        // Compute strain energy density given F
        double Psi = computeStrainEnergy(F);

        // Compute dPsi using finite difference, given dF
        if (order == 1){
            Array<double> F_tilde(N, N); // perturbed deformation gradient
            for (int i = 0; i < N; i++) {
                for (int J = 0; J < N; J++) {
                    // Perturb the iJ-th component of F by delta * dF(i,J)
                    for (int k = 0; k < N; k++) {
                        for (int l = 0; l < N; l++) {
                            F_tilde(k,l) = F(k,l) + delta * dF(k,l);
                        }
                    }
                }
            }

            // Compute Psi_tilde for perturbed deformation gradient
            double Psi_tilde = computeStrainEnergy(F_tilde);

            // Compute differences in Psi
            dPsi = Psi_tilde - Psi;
        }
        else if (order == 2){
            Array<double> F_plus(N,N); // positive perturbed deformation gradient
            Array<double> F_minus(N,N); // negative perturbed deformation gradient
            for (int i = 0; i < N; i++) {
                for (int J = 0; J < N; J++) {
                    // Perturb the iJ-th component of F by +-delta * dF(i,J)
                    for (int k = 0; k < N; k++) {
                        for (int l = 0; l < N; l++) {
                            F_plus(k,l) = F(k,l) + delta * dF(k,l);
                            F_minus(k,l) = F(k,l) - delta * dF(k,l);
                        }
                    }
                }
            }

            // Compute Psi_plus and Psi_minus for perturbed deformation gradient
            double Psi_plus = computeStrainEnergy(F_plus);
            double Psi_minus = computeStrainEnergy(F_minus);

            // Compute differences in Psi
            dPsi = Psi_plus - Psi_minus;
        }
    }

    /**
     * @brief Compute perturbed Green-Lagrange strain tensor (dE) given perturbed deformation gradient (dF) using finite differences
     * 
     * @param F Deformation gradient
     * @param dF Deformation gradient perturbation shape
     * @param delta  Deformation gradient perturbation scaling factor
     * @param order Order of the finite difference scheme (1 for first order, 2 for second order, etc.)
     * @param dE  Green-Lagrange strain tensor perturbation
     */
    void calcdEFiniteDifference(const Array<double> &F, const Array<double> &dF, const double delta, const int order, Array<double> &dE) {

        int N = F.ncols();
        assert(F.nrows() == N);

        // Compute E from F
        double J;
        Array<double> C(N, N), E(N, N);
        calc_JCE(F, J, C, E);

        // Compute dE using finite difference, given dF
        if (order == 1){
            Array<double> F_tilde(N,N); // perturbed deformation gradient
            for (int i = 0; i < N; i++) {
                for (int J = 0; J < N; J++) {
                    // Perturb the iJ-th component of F by delta * dF(i,J)
                    for (int k = 0; k < N; k++) {
                        for (int l = 0; l < N; l++) {
                            F_tilde(k,l) = F(k,l) + delta * dF(k,l);
                        }
                    }
                }
            }

            // Compute perturbed E_tilde from F_tilde
            double J_tilde;
            Array<double> C_tilde(N, N), E_tilde(N, N);
            calc_JCE(F_tilde, J_tilde, C_tilde, E_tilde);

            // Compute differences in E
            dE = E_tilde - E;
        }
        else if (order == 2){
            Array<double> F_plus(N,N); // positive perturbed deformation gradient
            Array<double> F_minus(N,N); // negative perturbed deformation gradient
            for (int i = 0; i < N; i++) {
                for (int J = 0; J < N; J++) {
                    // Perturb the iJ-th component of F by +-delta * dF(i,J)
                    for (int k = 0; k < N; k++) {
                        for (int l = 0; l < N; l++) {
                            F_plus(k,l) = F(k,l) + delta * dF(k,l);
                            F_minus(k,l) = F(k,l) - delta * dF(k,l);
                        }
                    }
                }
            }

            // Compute perturbed E_plus and E_minus from F_plus and F_minus
            double J_plus, J_minus;
            Array<double> C_plus(N, N), E_plus(N, N), C_minus(N, N), E_minus(N, N);
            calc_JCE(F_plus, J_plus, C_plus, E_plus);
            calc_JCE(F_minus, J_minus, C_minus, E_minus);

            // Compute differences in E
            dE = (E_plus - E_minus);
        }
    }

    /**
     * @brief Compute contraction of PK2 stress with perturbation in Green-Lagrange strain tensor (S:dE) given perturbation in deformation gradient (dF) using finite differences.
     * 
     * @param F Deformation gradient
     * @param dF Deformation gradient perturbation shape
     * @param delta Deformation gradient perturbation scaling factor
     * @param order Order of the finite difference scheme (1 for first order, 2 for second order, etc.)
     * @param SdE PK2 stress tensor times the perturbation in the Green-Lagrange strain tensor
     */
    void calcSdEFiniteDifference(const Array<double> &F, const Array<double> &dF, const double delta, const int order, double &SdE) {
        
        int N = F.ncols();
        assert(F.nrows() == N);

        // Compute S(F) from compute_pk2cc()
        Array<double> S(N,N), Dm(2*N,2*N);
        compute_pk2cc(F, S, Dm);

        // Compute dE using finite difference, given dF
        Array<double> dE(N,N);
        calcdEFiniteDifference(F, dF, delta, order, dE);

        // Compute S:dE
        SdE = mat_fun::mat_ddot(S, dE, N);
    }

    /**
     * @brief Tests the consistency of the PK2 stress tensor S(F) from compute_pk2cc() with the strain energy density Psi(F) provided by the user.
     *
     * Analytically, we should have S = dPsi/dE. This function checks whether S:dE = dPsi, where dE and dPsi are computed using finite differences in F.
     *
     * Pseudocode:
     * - Compute Psi(F)
     * - Compute S(F) from compute_pk2cc()
     * - For many random dF
     *      - Compute dPsi = Psi(F + dF) - Psi(F)
     *      - Compute dE = E(F + dF) - E(F)
     *      - Check that S:dE = dPsi
     * 
     * @param[in] F Deformation gradient.
     * @param[in] n_iter Number of random perturbations to test.
     * @param[in] rel_tol Relative tolerance for comparing dPsi and S:dE.
     * @param[in] abs_tol Absolute tolerance for comparing dPsi and S:dE.
     * @param[in] delta Perturbation scaling factor.
     * @param[in] verbose Show values of S, dE, SdE, and dPsi if true.
     * @return None.
     *
     */
    void testPK2StressConsistentWithStrainEnergy(const Array<double> &F, int n_iter, double rel_tol, double abs_tol, double delta, bool verbose = false) {
        int order = 2;

        int N = F.ncols();
        assert(F.nrows() == N);

        // Generate many random dF and check that S:dE = dPsi
        // S was obtained from compute_pk2cc(), and dPsi = Psi(F + dF) - Psi(F)
        double dPsi, SdE;
        for (int i = 0; i < n_iter; i++) {
            // Generate random dF
            auto dF = create_random_F(N, 0.0, 1.0);

            // Compute dPsi
            calcdPsiFiniteDifference(F, dF, delta, order, dPsi);

            // Compute SdE
            calcSdEFiniteDifference(F, dF, delta, order, SdE);

            // Check that S:dE = dPsi
            EXPECT_NEAR(SdE, dPsi, fmax(abs_tol, rel_tol * fabs(dPsi)));
            
            // Print results if verbose
            if (verbose) {
                std::cout << "Iteration " << i << ":" << std::endl;

                printMaterialParameters();

                std::cout << "F =" << std::endl;
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        std::cout << F(i,j) << " ";
                    }
                    std::cout << std::endl;
                }

                std::cout << "SdE = " << SdE << ", dPsi = " << dPsi << std::endl;
                std::cout << std::endl;
            }
        }
    }

    /**
     * @brief Tests the order of convergence of the consistency between dPsi and S:dE using finite differences.
     * 
     * Analytically, we should have S = dPsi/dE. This function determines the order of convergence of S:dE = dPsi, where dE and dPsi are computed using finite differences in F.
     *
     * Pseudocode:
     * - Compute Psi(F)
     * - Compute S(F) from compute_pk2cc()
     * - For each component-wise perturbation dF
     *      - Compute dPsi
     *      - Compute dE
     *      - Compute error S:dE - dPsi
     * - Compute order of convergence by fitting a line to log(delta) vs log(error)
     * 
     * Note that the order of convergence should be order + 1, because we are comparing differences (dPsi and S:dE)
     * instead of derivatives (e.g. dPsi/dF and S:dE/dF).
     * @param[in] F Deformation gradient.
     * @param[in] delta_max Maximum perturbation scaling factor.
     * @param[in] delta_min Minimum perturbation scaling factor.
     * @param[in] order Order of the finite difference scheme (1 for first order, 2 for second order, etc.).
     * @param[in] convergence_order_tol Tolerance for comparing convergence order with expected value
     * @param verbose Show values of errors and order of convergence if true.
     */
    void testPK2StressConsistencyConvergenceOrder(const Array<double> &F, double delta_max, double delta_min, int order, const double convergence_order_tol, bool verbose = false) {
        
        int N = F.ncols();
        assert(F.nrows() == N);

        // Check that delta_max > delta_min
        if (delta_max <= delta_min) {
            std::cerr << "Error: delta_max must be greater than delta_min." << std::endl;
            return;
        }

        // Check that order is 1 or 2
        if (order != 1 && order != 2) {
            std::cerr << "Error: order must be 1 or 2." << std::endl;
            return;
        }
        // Loop over perturbations to each component of F, dF
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                // Generate dF with 1.0 in (i,j) component
                Array<double> dF(N, N);
                dF(i,j) = 1.0;

                // Create list of deltas for convergence test (delta = delta_max, delta_max/2, delta_max/4, ...)
                std::vector<double> deltas;
                double delta = delta_max;
                while (delta >= delta_min) {
                    deltas.push_back(delta);
                    delta /= 2.0;
                }

                // Compute dPsi and S:dE for each delta and store error in list
                std::vector<double> errors;
                double dPsi, SdE;

                for (int i = 0; i < deltas.size(); i++) {
                    calcdPsiFiniteDifference(F, dF, deltas[i], order, dPsi);
                    calcSdEFiniteDifference(F, dF, deltas[i], order, SdE);

                    // Compute error between dPsi and S:dE
                    double error = fabs(dPsi - SdE);

                    // Store error in list
                    errors.push_back(error);
                }

                // Compute order of convergence by fitting a line to log(delta) vs log(error)
                std::vector<double> log_deltas, log_errors;
                for (int i = 0; i < deltas.size(); i++) {
                    log_deltas.push_back(log(deltas[i]));
                    log_errors.push_back(log(errors[i]));
                }

                // Fit a line to log(delta) vs log(error)
                // m is the slope (order of convergence), b is the intercept
                auto [m, b] = computeLinearRegression(log_deltas, log_errors);

                // Check that order of convergence is > (order + 1) - convergence_order_tol
                EXPECT_GT(m, order + 1 - convergence_order_tol);

                // Print results if verbose
                if (verbose) {
                    std::cout << "Slope (order of convergence): " << m << std::endl;
                    std::cout << "Intercept: " << b << std::endl;
                    std::cout << "Errors: ";
                    for (int i = 0; i < errors.size(); i++) {
                        std::cout << errors[i] << " ";
                    }
                    std::cout << std::endl;
                    std::cout << std::endl;
                    
                    std::cout << "F = " << std::endl;
                    for (int i = 0; i < 3; i++) {
                        for (int J = 0; J < 3; J++) {
                            std::cout << F(i,J) << " ";
                        }
                        std::cout << std::endl;
                    }
                }
            }
        }
    }

    /**
     * @brief Compute perturbation in PK2 stress (dS) given perturbation in deformation gradient (dF) using finite differences
     * 
     * @param F Deformation gradient
     * @param dF Deformation gradient perturbation shape
     * @param delta Deformation gradient perturbation scaling factor
     * @param order Order of the finite difference scheme (1 for first order, 2 for second order, etc.)
     * @param dS PK2 stress tensor perturbation
     */
    void calcdSFiniteDifference(const Array<double> &F, const Array<double> &dF, const double delta, const int order, Array<double> &dS) {
        
        int N = F.ncols();
        assert(F.nrows() == N);

        // Compute S(F) from compute_pk2cc()
        Array<double> S(N,N), Dm(2*N,2*N);
        compute_pk2cc(F, S, Dm);

        // Compute dS using finite difference, given dF
        if (order == 1){
            Array<double> F_tilde(N,N); // perturbed deformation gradient
            for (int i = 0; i < N; i++) {
                for (int J = 0; J < N; J++) {
                    // Perturb the iJ-th component of F by delta * dF[i][J]
                    for (int k = 0; k < N; k++) {
                        for (int l = 0; l < N; l++) {
                            F_tilde(k,l) = F(k,l) + delta * dF(k,l);
                        }
                    }
                }
            }

            // Compute perturbed S_tilde from F_tilde
            Array<double> S_tilde(N,N), Dm_tilde(2*N,2*N);
            compute_pk2cc(F_tilde, S_tilde, Dm_tilde);

            // Compute differences in S
            dS = S_tilde - S;
        }
        else if (order == 2){
            Array<double> F_plus(N,N); // positive perturbed deformation gradient
            Array<double> F_minus(N,N); // negative perturbed deformation gradient
            for (int i = 0; i < N; i++) {
                for (int J = 0; J < N; J++) {
                    // Perturb the iJ-th component of F by +-delta * dF(i,J)
                    for (int k = 0; k < N; k++) {
                        for (int l = 0; l < N; l++) {
                            F_plus(k,l) = F(k,l) + delta * dF(k,l);
                            F_minus(k,l) = F(k,l) - delta * dF(k,l);
                        }
                    }
                }
            }

            // Compute perturbed S_plus and S_minus from F_plus and F_minus
            Array<double> S_plus(N,N), Dm_plus(2*N,2*N);
            Array<double> S_minus(N,N), Dm_minus(2*N,2*N);

            compute_pk2cc(F_plus, S_plus, Dm_plus);
            compute_pk2cc(F_minus, S_minus, Dm_minus);

            // Compute differences in S
            dS = S_plus - S_minus;
        }
    }

    /**
     * @brief Compute material elasticity tensor contracted with perturbation in Green-Lagrange strain tensor (CC:dE) given perturbation in deformation gradient (dF) using finite differences
     * 
     * @param F Deformation gradient
     * @param dF Deformation gradient perturbation shape
     * @param delta Deformation gradient perturbation scaling factor
     * @param order Order of the finite difference scheme (1 for first order, 2 for second order, etc.)
     * @param CCdE Material elasticity tensor times the perturbation in the Green-Lagrange strain tensor
     */
    void calcCCdEFiniteDifference(const Array<double> &F, const Array<double> &dF, const double delta, const int order, Array<double> &CCdE) {
        
        int N = F.ncols();
        assert(F.nrows() == N);

        // Compute CC(F) from compute_pk2cc()
        Array<double> S(N, N), Dm(2*N, 2*N);
        compute_pk2cc(F, S, Dm);
        
        // Create Tensor for CC
        Tensor4<double> CC(N, N, N, N);

        mat_models::voigt_to_cc(N, Dm, CC);

        // Compute dE using finite difference, given dF
        Array<double> dE(N, N);
        calcdEFiniteDifference(F, dF, delta, order, dE);

        // Compute CC:dE
        CCdE = mat_fun::ten_mddot(CC, dE, N);
    }


    /**
     * @brief Tests the consistency of the material elasticity tensor CC(F) from compute_pk2cc() with the PK2 stress tensor S(F) from compute_pk2cc().
     *
     * Analytically, we should have CC:dE = dS. This function checks whether CC:dE = dS, where dE and dS are computed using finite differences in F.
     *
     * Pseudocode:
     * - Compute S(F) and CC(F) from compute_pk2cc()
     * - For each component-wise perturbation dF
     *      - Compute S(F + dF) from compute_pk2cc()
     *      - Compute dS = S(F + dF) - S(F)
     *      - Compute dE from dF
     *      - Check that CC:dE = dS
     * 
     * @param[in] F Deformation gradient.
     * @param[in] rel_tol Relative tolerance for comparing dS and CC:dE.
     * @param[in] abs_tol Absolute tolerance for comparing dS and CC:dE.
     * @param[in] delta Perturbation scaling factor.
     * @param[in] verbose Show values of CC, dE, CCdE, and dS if true.
     * @return None.
     */
    void testMaterialElasticityConsistentWithPK2Stress(const Array<double> &F, double rel_tol, double abs_tol, double delta, bool verbose = false) {
        int order = 2;

        int N = F.ncols();
        assert(F.nrows() == N);

        // Compute E from F
        double J;
        Array<double> C(N, N), E(N, N);
        calc_JCE(F, J, C, E);

        // Compute S_ij(F)
        // Compute CC_ijkl(F). 
        // CC is provided in Voigt notation as Dm, and we will convert it to CC
        Array<double> S(N, N), Dm(2*N, 2*N);
        compute_pk2cc(F, S, Dm); // S from solver 

        // Calculate CC from Dm
        Tensor4<double> CC(N, N, N, N);
        mat_models::voigt_to_cc(N, Dm, CC);
    
        // ------- Ancillary test ---------
        // Calculate Dm_check from CC
        Array<double> Dm_check(2*N, 2*N);
        mat_models::cc_to_voigt(N, CC, Dm_check);
    
        // Check that Dm_check = Dm, for sanity
        for (int i = 0; i < 2*N; i++) {
            for (int j = 0; j < 2*N; j++) {
                EXPECT_NEAR(Dm_check(i,j), Dm(i,j), abs_tol);
            }
        }
        // -------------------------------
    
        // Generate many random dF and check that CC:dE = dS
        // CC was obtained from compute_pk2cc(), and dS = S(F + dF) - S(F), 
        // where S is also obtained from compute_pk2cc()
        Array<double> dS(N, N), CCdE(N, N);
        
        // Loop over perturbations to each component of F, dF
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                // Generate dF with 1.0 in (i,j) component
                Array<double> dF(N, N);
                dF(i,j) = 1.0;
    
                // Compute dS
                calcdSFiniteDifference(F, dF, delta, order, dS);

                // Compute CC:dE
                calcCCdEFiniteDifference(F, dF, delta, order, CCdE);
        
                // Check that CC_ijkl dE_kl = dS_ij
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        EXPECT_NEAR(CCdE(i,j), dS(i,j), fmax(abs_tol, rel_tol * fabs(dS(i,j))));
                    }
                }
        
                // Print results if verbose
                if (verbose) {
                    std::cout << "Iteration " << i << ":" << std::endl;
        
                    printMaterialParameters();
        
                    std::cout << "F =" << std::endl;
                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < 3; j++) {
                            std::cout << F(i,j) << " ";
                        }
                        std::cout << std::endl;
                    }
        
                    std::cout << "CC =" << std::endl;
                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < 3; j++) {
                            for (int k = 0; k < 3; k++) {
                                for (int l = 0; l < 3; l++) {
                                    std::cout << CC(i,j,k,l) << " ";
                                }
                                std::cout << std::endl;
                            }
                        }
                        std::cout << std::endl;
                    }
        
                    std::cout << "dS =" << std::endl;
                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < 3; j++) {
                            std::cout << dS(i,j) << " ";
                        }
                        std::cout << std::endl;
                    }
        
                    std::cout << "CCdE =" << std::endl;
                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < 3; j++) {
                            std::cout << CCdE(i,j) << " ";
                        }
                        std::cout << std::endl;
                    }
                    std::cout << std::endl;
                }
            }
        }
    }

    /**
     * @brief Tests the order of convergence of the consistency between CC:dE and dS using finite differences.
     * 
     * Analytically, we should have CC:dE = dS. This function determines the order of convergence of CC:dE = dS, where dE and dS are computed using finite differences in F.
     *
     * Pseudocode:
     * - For each component-wise perturbation dF
     *     - For decreasing delta
     *       - Compute dS
     *       - Compute CC:dE 
     *       - Compute error CC:dE - dS
     * - Compute order of convergence by fitting a line to log(delta) vs log(error)
     * 
     * Note that the order of convergence should be order + 1, because we are comparing differences (dS and CC:dE)
     * instead of derivatives (e.g. dS/dF and CC:dE/dF).
     * @param[in] F Deformation gradient.
     * @param[in] delta_max Maximum perturbation scaling factor.
     * @param[in] delta_min Minimum perturbation scaling factor.
     * @param[in] order Order of the finite difference scheme (1 for first order, 2 for second order, etc.).
     * @param[in] convergence_order_tol Tolerance for comparing convergence order with expected value
     * @param[in] verbose Show values of errors and order of convergence if true.
     */
    void testMaterialElasticityConsistencyConvergenceOrder(const Array<double> &F, double delta_max, double delta_min, int order, const double convergence_order_tol, bool verbose = false) {
        
        int N = F.ncols();
        assert(F.nrows() == N);
        
        // Check that delta_max > delta_min
        if (delta_max <= delta_min) {
            std::cerr << "Error: delta_max must be greater than delta_min." << std::endl;
            return;
        }

        // Check that order is 1 or 2
        if (order != 1 && order != 2) {
            std::cerr << "Error: order must be 1 or 2." << std::endl;
            return;
        }

        // Create list of deltas for convergence test (delta = delta_max, delta_max/2, delta_max/4, ...)
        std::vector<double> deltas;
        double delta = delta_max;
        while (delta >= delta_min) {
            deltas.push_back(delta);
            delta /= 2.0;
        }

        // Loop over perturbations to each component of F, dF
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                // Generate dF with 1.0 in (i,j) component
                Array<double> dF(N, N);
                dF(i,j) = 1.0;
            
                // Compute dS and CC:dE for each delta and store error in list
                std::vector<double> errors;
                Array<double> dS(N, N), CCdE(N, N);
                for (int i = 0; i < deltas.size(); i++) {
                    calcdSFiniteDifference(F, dF, deltas[i], order, dS);
                    calcCCdEFiniteDifference(F, dF, deltas[i], order, CCdE);

                    // Compute Frobenius norm of error between dS and CC:dE
                    double error = 0.0;
                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < 3; j++) {
                            error += pow(dS(i,j) - CCdE(i,j), 2);
                        }
                    }
                    error = sqrt(error);

                    // Store error in list
                    errors.push_back(error);
                }

                // Compute order of convergence by fitting a line to log(delta) vs log(error)
                std::vector<double> log_deltas, log_errors;
                for (int i = 0; i < deltas.size(); i++) {
                    log_deltas.push_back(log(deltas[i]));
                    log_errors.push_back(log(errors[i]));
                }

                // Fit a line to log(delta) vs log(error)
                // m is the slope (order of convergence), b is the intercept
                auto [m, b] = computeLinearRegression(log_deltas, log_errors);

                // Check that order of convergence is > (order + 1) - convergence_order_tol
                EXPECT_GT(m, order + 1 - convergence_order_tol);

                // Print results if verbose
                if (verbose) {
                    std::cout << "Iteration " << i << ":" << std::endl;
                    std::cout << "Slope (order of convergence): " << m << std::endl;
                    std::cout << "Intercept: " << b << std::endl;
                    std::cout << "Errors: ";
                    for (int i = 0; i < errors.size(); i++) {
                        std::cout << errors[i] << " ";
                    }
                    std::cout << std::endl;
                    std::cout << std::endl;
                    
                    std::cout << "F = " << std::endl;
                    for (int i = 0; i < 3; i++) {
                        for (int J = 0; J < 3; J++) {
                            std::cout << F(i,J) << " ";
                        }
                        std::cout << std::endl;
                    }

                    std::cout << "dF = " << std::endl;
                    for (int i = 0; i < 3; i++) {
                        for (int J = 0; J < 3; J++) {
                            std::cout << dF(i,J) << " ";
                        }
                        std::cout << std::endl;
                    }
                    std::cout << std::endl;
                }
            }
        }
    }

    /**
     * @brief Compares the PK2 stress tensor S(F) with a reference solution.
     *
     * This function computes the PK2 stress tensor S(F) from the deformation gradient F using compute_pk2cc() 
     * and compares it with a reference solution S_ref. The comparison is done using relative and absolute tolerances.
     *
     * @param[in] F Deformation gradient.
     * @param[in] S_ref Reference solution for PK2 stress.
     * @param[in] rel_tol Relative tolerance for comparing S with S_ref.
     * @param[in] abs_tol Absolute tolerance for comparing S with S_ref.
     * @param[in] verbose Show values of F, S, and S_ref if true.
     * @return None.
     */
    void testPK2StressAgainstReference(const Array<double> &F, const Array<double> &S_ref, double rel_tol, double abs_tol, bool verbose = false) {
        
        int N = F.ncols();
        assert(F.nrows() == N);

        // Compute S(F) from compute_pk2cc()
        Array<double> S(N, N), Dm(2*N, 2*N);
        compute_pk2cc(F, S, Dm);
    
        // Compare S with reference solution
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                EXPECT_NEAR(S(i,j), S_ref(i,j), fmax(abs_tol, rel_tol * fabs(S_ref(i,j))));
            }
        }
    
        // Print results if verbose
        if (verbose) {
            printMaterialParameters();
    
            std::cout << "F =" << std::endl;
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    std::cout << F(i,j) << " ";
                }
                std::cout << std::endl;
            }
    
            std::cout << "S =" << std::endl;
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    std::cout << S(i,j) << " ";
                }
                std::cout << std::endl;
            }
    
            std::cout << "S_ref =" << std::endl;
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    std::cout << S_ref(i,j) << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }

    /**
     * @brief Compares the material elasticity tensor CC(F) with a reference solution.
     *
     * This function computes the material elasticity tensor CC(F) from the deformation gradient F using compute_pk2cc() 
     * and compares it with a reference solution CC_ref. The comparison is done using relative and absolute tolerances.
     *
     * @param[in] F Deformation gradient.
     * @param[in] CC_ref Reference solution for material elasticity tensor.
     * @param[in] rel_tol Relative tolerance for comparing CC with CC_ref.
     * @param[in] abs_tol Absolute tolerance for comparing CC with CC_ref.
     * @param[in] verbose Show values of F, CC, and CC_ref if true.
     * @return None.
     */
    void testMaterialElasticityAgainstReference(const Array<double> &F, const Tensor4<double> &CC_ref, double rel_tol, double abs_tol, bool verbose = false) {
        
        int N = F.ncols();
        assert(F.nrows() == N);

        // Compute CC(F) from compute_pk2cc()
        Array<double> S(N, N), Dm(2*N, 2*N);
        compute_pk2cc(F, S, Dm);

        // Calculate CC from Dm
        Tensor4<double> CC(N, N, N, N);
        mat_models::voigt_to_cc(N, Dm, CC);
    
        // Compare CC with reference solution
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    for (int l = 0; l < 3; l++) {
                        EXPECT_NEAR(CC(i,j,k,l), CC_ref(i,j,k,l),
                        fmax(abs_tol, rel_tol * fabs(CC_ref(i,j,k,l))));  
                    }
                }
            }
        }
    
        // Print results if verbose
        if (verbose) {
            printMaterialParameters();
    
            std::cout << "F =" << std::endl;
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    std::cout << F(i,j) << " ";
                }
                std::cout << std::endl;
            }
    
            std::cout << "CC =" << std::endl;
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    for (int k = 0; k < N; k++) {
                        for (int l = 0; l < N; l++) {
                            std::cout << CC(i,j,k,l) << " ";
                        }
                        std::cout << std::endl;
                    }
                }
                std::cout << std::endl;
            }
    
            std::cout << "CC_ref =" << std::endl;
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    for (int k = 0; k < N; k++) {
                        for (int l = 0; l < N; l++) {
                            std::cout << CC_ref(i,j,k,l) << " ";
                        }
                        std::cout << std::endl;
                    }
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }

       /**
     * @brief Tests rho, beta, drho/dp, and dbeta/dp from g_vol_pen() against reference solutions.
     *
     * This function computes rho, beta, drho/dp, and dbeta/dp from the pressure p using g_vol_pen() 
     * and compares them with reference solutions.
     * These values are required for treating a volumetric penalty term in the ustruct formulation.
     *
     * @param[in] p Pressure.
     * @param[in] rho0 Initial solid density.
     * @param[in] rho_ref Reference solution for rho.
     * @param[in] beta_ref Reference solution for beta.
     * @param[in] drhodp_ref Reference solution for drho/dp.
     * @param[in] dbetadp_ref Reference solution for dbeta/dp.
     * @param[in] rel_tol Relative tolerance for comparing rho and beta with reference solutions.
     * @param[in] abs_tol Absolute tolerance for comparing rho and beta with reference solutions.
     * @param[in] verbose Show values of p, rho, beta, rho_ref, beta_ref if true.
     * @return None.
     */
    void testRhoBetaAgainstReference(double p, double rho0, double rho_ref, double beta_ref, double drhodp_ref, double dbetadp_ref, double rel_tol, double abs_tol, bool verbose = false) {
        double rho, beta, drhodp, dbetadp;
        double Ja = 1.0; // Active strain Jacobian (not used in this function)
    
        // Compute rho, beta, drhodp, dbetadp from g_vol_pen()
        g_vol_pen(p, rho0, rho, beta, drhodp, dbetadp, Ja);
    
        // Compare rho, beta, drho, dbeta with reference solutions
        EXPECT_NEAR(rho, rho_ref, fmax(abs_tol, rel_tol * fabs(rho_ref)));
        EXPECT_NEAR(beta, beta_ref, fmax(abs_tol, rel_tol * fabs(beta_ref)));
        EXPECT_NEAR(drhodp, drhodp_ref, fmax(abs_tol, rel_tol * fabs(drhodp_ref)));
        EXPECT_NEAR(dbetadp, dbetadp_ref, fmax(abs_tol, rel_tol * fabs(dbetadp_ref)));
    
        // Print results if verbose
        if (verbose) {
            printMaterialParameters();
    
            std::cout << "p = " << p << std::endl;
            std::cout << "rho0 = " << rho0 << std::endl;
            std::cout << "rho = " << rho << ", rho_ref = " << rho_ref << std::endl;
            std::cout << "beta = " << beta << ", beta_ref = " << beta_ref << std::endl;
            std::cout << "drhodp = " << drhodp << ", drhodp_ref = " << drhodp_ref << std::endl;
            std::cout << "dbetadp = " << dbetadp << ", dbetadp_ref = " << dbetadp_ref << std::endl;
            std::cout << std::endl;
        }
    }
};



// --------------------------------------------------------------
// --------------------- Material Parameters Classes ------------
// --------------------------------------------------------------

// Class to contain material parameters
class MatParams {
public:
    virtual ~MatParams() {} // Virtual destructor for proper cleanup
};

// Class to contain Neo-Hookean material parameters
class NeoHookeanParams : public MatParams {
public:
    double C10;

    // Default constructor
    NeoHookeanParams() : C10(0.0) {}

    // Constructor with parameters
    NeoHookeanParams(double c10) : C10(c10) {}

};

// Class to contain Mooney-Rivlin material parameters
class MooneyRivlinParams : public MatParams {
public:
    double C01;
    double C10;

    // Default constructor
    MooneyRivlinParams() : C01(0.0), C10(0.0) {}

    // Constructor with parameters
    MooneyRivlinParams(double c01, double c10) : C01(c01), C10(c10) {}

};

// Class to contain Holzapfel-Ogden material parameters
class HolzapfelOgdenParams : public MatParams {
public:
    double a;    
    double b;
    double a_f;
    double b_f;
    double a_s;
    double b_s;
    double a_fs;
    double b_fs;
    double f[3];    // Fiber direction
    double s[3];    // Sheet direction

    double k; // Smoothed Heaviside function parameter

    // Default constructor
    HolzapfelOgdenParams() : a(0.0), b(0.0), a_f(0.0), b_f(0.0), a_s(0.0), b_s(0.0), a_fs(0.0), b_fs(0.0), k(0.0) {
        for (int i = 0; i < 3; i++) {
            f[i] = 0.0;
            s[i] = 0.0;
        }
    }

    // Constructor with parameters
    HolzapfelOgdenParams(double a, double b, double a_f, double b_f, double a_s, double b_s, double a_fs, double b_fs, double k, double f[3], double s[3]) : a(a), b(b), a_f(a_f), b_f(b_f), a_s(a_s), b_s(b_s), a_fs(a_fs), b_fs(b_fs), k(k) {
        for (int i = 0; i < 3; i++) {
            this->f[i] = f[i];
            this->s[i] = s[i];
        }
    }
};

// Class to contain Holzapfel-Ogden (Modified Anisortopy) material parameters
class HolzapfelOgdenMAParams : public MatParams {
public:
    double a;    
    double b;
    double a_f;
    double b_f;
    double a_s;
    double b_s;
    double a_fs;
    double b_fs;
    double f[3];    // Fiber direction
    double s[3];    // Sheet direction

    double k; // Smoothed Heaviside function parameter

    // Default constructor
    HolzapfelOgdenMAParams() : a(0.0), b(0.0), a_f(0.0), b_f(0.0), a_s(0.0), b_s(0.0), a_fs(0.0), b_fs(0.0), k(0.0) {
        for (int i = 0; i < 3; i++) {
            f[i] = 0.0;
            s[i] = 0.0;
        }
    }

    // Constructor with parameters
    HolzapfelOgdenMAParams(double a, double b, double a_f, double b_f, double a_s, double b_s, double a_fs, double b_fs, double k, double f[3], double s[3]) : a(a), b(b), a_f(a_f), b_f(b_f), a_s(a_s), b_s(b_s), a_fs(a_fs), b_fs(b_fs), k(k) {
        for (int i = 0; i < 3; i++) {
            this->f[i] = f[i];
            this->s[i] = s[i];
        }
    }
};

// Class to contain volumetric penalty parameters (just the penalty parameter)
class VolumetricPenaltyParams : public MatParams {
public:
    double kappa;

    // Default constructor
    VolumetricPenaltyParams() : kappa(0.0) {}

    // Constructor with parameters
    VolumetricPenaltyParams(double kappa) : kappa(kappa) {}
};


// --------------------------------------------------------------
// --------------------- Material Model Classes -----------------
// --------------------------------------------------------------

/**
 * @brief Class for testing the Neo-Hookean material model.
 *
 * This class provides methods to set up and test the Neo-Hookean material model, including 
 * computing the strain energy and printing material parameters.
 */
class TestNeoHookean : public TestMaterialModel {
public:

    /**
     * @brief Parameters for the Neo-Hookean material model.
     */
    NeoHookeanParams params;

    /**
     * @brief Constructor for the TestNeoHookean class.
     *
     * Initializes the Neo-Hookean material parameters for svMultiPhysics.
     *
     * @param[in] params_ Parameters for the Neo-Hookean material model.
     */
    TestNeoHookean(const NeoHookeanParams &params_) : TestMaterialModel( consts::ConstitutiveModelType::stIso_nHook, consts::ConstitutiveModelType::stVol_ST91),
        params(params_) 
        {
        // Set Neo-Hookean material parameters for svMultiPhysics
        auto &dmn = com_mod.mockEq.mockDmn;
        dmn.stM.C10 = params.C10;
        dmn.stM.Kpen = 0.0;         // Zero volumetric penalty parameter
    }

    /**
     * @brief Prints the Neo-Hookean material parameters.
     */
    void printMaterialParameters() {
        std::cout << "C10 = " << params.C10 << std::endl;
    }

    /**
     * @brief Computes the strain energy for the Neo-Hookean material model.
     *
     * @param[in] F Deformation gradient.
     * @return Strain energy density for the Neo-Hookean material model.
     */
    double computeStrainEnergy(const Array<double> &F) {
        // Compute solid mechanics terms
        solidMechanicsTerms smTerms = calcSolidMechanicsTerms(F);

        // Strain energy density for Neo-Hookean material model
        // Psi_iso = C10 * (Ib1 - 3)
        double Psi_iso = params.C10 * (smTerms.Ib1 - 3.);

        return Psi_iso;
    }
};

/**
 * @brief Class for testing the Mooney-Rivlin material model.
 *
 * This class provides methods to set up and test the Mooney-Rivlin material model, including 
 * computing the strain energy and printing material parameters.
 */
class TestMooneyRivlin : public TestMaterialModel {
public:

    /**
     * @brief Parameters for the Mooney-Rivlin material model.
     */
    MooneyRivlinParams params;

    /**
     * @brief Constructor for the TestMooneyRivlin class.
     *
     * Initializes the Mooney-Rivlin material parameters for svMultiPhysics.
     *
     * @param[in] params_ Parameters for the Mooney-Rivlin material model.
     */
    TestMooneyRivlin(const MooneyRivlinParams &params_) : TestMaterialModel( consts::ConstitutiveModelType::stIso_MR, consts::ConstitutiveModelType::stVol_ST91),
        params(params_) 
        {
        // Set Mooney-Rivlin material parameters for svMultiPhysics
        auto &dmn = com_mod.mockEq.mockDmn;
        dmn.stM.C01 = params.C01;
        dmn.stM.C10 = params.C10;
        dmn.stM.Kpen = 0.0;         // Zero volumetric penalty parameter
    }

    /**
     * @brief Prints the Mooney-Rivlin material parameters.
     */
    void printMaterialParameters() {
        std::cout << "C01 = " << params.C01 << ", C10 = " << params.C10 << std::endl;
    }

    /**
     * @brief Computes the strain energy for the Mooney-Rivlin material model.
     *
     * @param[in] F Deformation gradient.
     * @return Strain energy density for the Mooney-Rivlin material model.
     */
    double computeStrainEnergy(const Array<double> &F) {
        // Compute solid mechanics terms
        solidMechanicsTerms smTerms = calcSolidMechanicsTerms(F);

        // Strain energy density for Mooney-Rivlin material model
        // Psi_iso = C10 * (Ib1 - 3) + C01 * (Ib2 - 3)
        double Psi_iso = params.C10 * (smTerms.Ib1 - 3.) + params.C01 * (smTerms.Ib2 - 3.);

        return Psi_iso;
    }
};


/**
 * @brief Class for testing the Holzapfel-Ogden material model. 
 * 
 * This class provides methods to set up and test the Holzapfel-Ogden material 
 * model, including computing the strain energy and printing material parameters.
 */
class TestHolzapfelOgden : public TestMaterialModel {
public:

    /**
     * @brief Parameters for the Holzapfel-Ogden material model.
     */
    HolzapfelOgdenParams params;

    /**
     * @brief Constructor for the TestHolzapfelOgden class.
     *
     * Initializes the Holzapfel-Ogden material parameters for svMultiPhysics.
     *
     * @param[in] params_ Parameters for the Holzapfel-Ogden material model.
     */
    TestHolzapfelOgden(const HolzapfelOgdenParams &params_) : TestMaterialModel( consts::ConstitutiveModelType::stIso_HO, consts::ConstitutiveModelType::stVol_ST91),
        params(params_) 
        {
        // Set Holzapfel-Ogden material parameters for svMultiPhysics
        auto &dmn = com_mod.mockEq.mockDmn;
        dmn.stM.a = params.a;
        dmn.stM.b = params.b;
        dmn.stM.aff = params.a_f;
        dmn.stM.bff = params.b_f;
        dmn.stM.ass = params.a_s;
        dmn.stM.bss = params.b_s;
        dmn.stM.afs = params.a_fs;
        dmn.stM.bfs = params.b_fs;
        dmn.stM.khs = params.k;     // Smoothed Heaviside function parameter
        dmn.stM.Kpen = 0.0;         // Zero volumetric penalty parameter

        // Set number of fiber directions and fiber directions
        nFn = 2;
        Vector<double> f = {params.f[0], params.f[1], params.f[2]};
        Vector<double> s = {params.s[0], params.s[1], params.s[2]};
        fN.set_col(0, f);
        fN.set_col(1, s);
    }

    /**
     * @brief Prints the Holzapfel-Ogden material parameters.
     */
    void printMaterialParameters() {
        std::cout << "a = " << params.a << std::endl;
        std::cout << "b = " << params.b << std::endl;
        std::cout << "a_f = " << params.a_f << std::endl;
        std::cout << "b_f = " << params.b_f << std::endl;
        std::cout << "a_s = " << params.a_s << std::endl;
        std::cout << "b_s = " << params.b_s << std::endl;
        std::cout << "a_fs = " << params.a_fs << std::endl;
        std::cout << "b_fs = " << params.b_fs << std::endl;
        std::cout << "k = " << params.k << std::endl;
        std::cout << "f = " << "[" << params.f[0] << " " << params.f[1] << " " << params.f[2] << "]" << std::endl;
        std::cout << "s = " << "[" << params.s[0] << " " << params.s[1] << " " << params.s[2] << "]" << std::endl;
    }

    /**
     * @brief Smoothed Heaviside function centered at x = 1.
     * 
     * @param[in] x Input value.
     * @param[in] k Smoothing parameter.
     * @return Smoothed Heaviside function.
     */
    double chi(const double x, const double k=100) const {
        return 1. / (1. + exp(-k * (x - 1.)));
    }

    /**
     * @brief Computes the strain energy for the Holzapfel-Ogden material model.
     *
     * @param[in] F Deformation gradient.
     * @return Strain energy density for the Holzapfel-Ogden material model.
     */
    double computeStrainEnergy(const Array<double> &F) {
        // Compute solid mechanics terms
        solidMechanicsTerms smTerms = calcSolidMechanicsTerms(F);

        // Material parameters
        double a = params.a;
        double b = params.b;
        double a_f = params.a_f;
        double b_f = params.b_f;
        double a_s = params.a_s;
        double b_s = params.b_s;
        double a_fs = params.a_fs;
        double b_fs = params.b_fs;

        // Smoothed Heaviside parameter
        double k = params.k;

        // Fiber and sheet directions
        Vector<double> f = {params.f[0], params.f[1], params.f[2]};
        Vector<double> s = {params.s[0], params.s[1], params.s[2]};

        // Strain energy density for Holzapfel-Ogden material model

        // Formulation with fully decoupled isochoric-volumetric split
        // Uses I1_bar, I4_bar_f, I4_bar_s, I8_bar_fs (bar = isochoric)
        // Psi = a/2b * exp{b(I1_bar - 3)} 
        //       + a_f/2b_f * chi(I4_bar_f) * (exp{b_f(I4_bar_f - 1)^2} - 1
        //       + a_s/2b_s * chi(I4_bar_s) * (exp{b_s(I4_bar_s - 1)^2} - 1
        //       + a_fs/2b_fs * (exp{b_fs*I8_bar_fs^2} - 1)
        // This corresponds to the HO implementation in svMultiPhysics

        // Invariants
        double I1_bar = smTerms.Ib1;
        // I4_bar_f = f . C_bar . f
        auto C_bar_f = mat_fun::mat_mul(smTerms.C_bar, f);
        double I4_bar_f = f.dot(C_bar_f);
        // I4_bar_s = s . C_bar . s
        auto C_bar_s = mat_fun::mat_mul(smTerms.C_bar, s);
        double I4_bar_s = s.dot(C_bar_s);
        // I8_bar_fs = f . C_bar . s
        double I8_bar_fs = f.dot(C_bar_s);

        // Strain energy density for Holzapfel-Ogden material model with modified anisotropic invariants (bar quantities)
        double Psi = 0.0;
        Psi += a / (2.0 * b) * exp(b * (I1_bar - 3.0));                             // Isotropic term
        Psi += a_f / (2.0 * b_f) * chi(I4_bar_f, k) * (exp(b_f * pow(I4_bar_f - 1.0, 2)) - 1.0);   // Fiber term
        Psi += a_s / (2.0 * b_s) * chi(I4_bar_s, k) * (exp(b_s * pow(I4_bar_s - 1.0, 2)) - 1.0);   // Sheet term
        Psi += a_fs / (2.0 * b_fs) * (exp(b_fs * pow(I8_bar_fs, 2)) - 1.0);                   // Cross-fiber term

        return Psi;
    }
};


/**
 * @brief Class for testing the Holzapfel-Ogden (Modified Anisotropy) material model. 
 * 
 * This class provides methods to set up and test the Holzapfel-Ogden-ma material 
 * model, including computing the strain energy and printing material parameters.
 *
 */
class TestHolzapfelOgdenMA : public TestMaterialModel {
public:

    /**
     * @brief Parameters for the Holzapfel-Ogden ma material model.
     */
    HolzapfelOgdenMAParams params;

    /**
     * @brief Constructor for the TestHolzapfelOgdenMA class.
     *
     * Initializes the Holzapfel-Ogden material parameters for svMultiPhysics.
     *
     * @param[in] params_ Parameters for the Holzapfel-Ogden ma material model.
     */
    TestHolzapfelOgdenMA(const HolzapfelOgdenMAParams &params_) : TestMaterialModel( consts::ConstitutiveModelType::stIso_HO_ma, consts::ConstitutiveModelType::stVol_ST91),
        params(params_) 
        {
        // Set Holzapfel-Ogden material parameters for svMultiPhysics
        auto &dmn = com_mod.mockEq.mockDmn;
        dmn.stM.a = params.a;
        dmn.stM.b = params.b;
        dmn.stM.aff = params.a_f;
        dmn.stM.bff = params.b_f;
        dmn.stM.ass = params.a_s;
        dmn.stM.bss = params.b_s;
        dmn.stM.afs = params.a_fs;
        dmn.stM.bfs = params.b_fs;
        dmn.stM.khs = params.k;     // Smoothed Heaviside function parameter
        dmn.stM.Kpen = 0.0;         // Zero volumetric penalty parameter

        // Set number of fiber directions and fiber directions
        nFn = 2;
        Vector<double> f = {params.f[0], params.f[1], params.f[2]};
        Vector<double> s = {params.s[0], params.s[1], params.s[2]};
        fN.set_col(0, f);
        fN.set_col(1, s);
    }

    /**
     * @brief Prints the Holzapfel-Ogden material parameters.
     */
    void printMaterialParameters() {
        std::cout << "a = " << params.a << std::endl;
        std::cout << "b = " << params.b << std::endl;
        std::cout << "a_f = " << params.a_f << std::endl;
        std::cout << "b_f = " << params.b_f << std::endl;
        std::cout << "a_s = " << params.a_s << std::endl;
        std::cout << "b_s = " << params.b_s << std::endl;
        std::cout << "a_fs = " << params.a_fs << std::endl;
        std::cout << "b_fs = " << params.b_fs << std::endl;
        std::cout << "k = " << params.k << std::endl;
        std::cout << "f = " << "[" << params.f[0] << " " << params.f[1] << " " << params.f[2] << "]" << std::endl;
        std::cout << "s = " << "[" << params.s[0] << " " << params.s[1] << " " << params.s[2] << "]" << std::endl;
    }

    /**
     * @brief Smoothed Heaviside function centered at x = 1.
     * 
     * @param[in] x Input value.
     * @param[in] k Smoothing parameter.
     * @return Smoothed Heaviside function.
     */
    double chi(const double x, const double k=100) const {
        return 1. / (1. + exp(-k * (x - 1.)));
    }

    /**
     * @brief Computes the strain energy for the Holzapfel-Ogden material model.
     *
     * @param[in] F Deformation gradient.
     * @return Strain energy density for the Holzapfel-Ogden material model.
     */
    double computeStrainEnergy(const Array<double> &F) {
        // Compute solid mechanics terms
        solidMechanicsTerms smTerms = calcSolidMechanicsTerms(F);

        // Material parameters
        double a = params.a;
        double b = params.b;
        double a_f = params.a_f;
        double b_f = params.b_f;
        double a_s = params.a_s;
        double b_s = params.b_s;
        double a_fs = params.a_fs;
        double b_fs = params.b_fs;

        // Smoothed Heaviside parameter
        double k = params.k;

        // Fiber and sheet directions
        Vector<double> f = {params.f[0], params.f[1], params.f[2]};
        Vector<double> s = {params.s[0], params.s[1], params.s[2]};

        // Strain energy density for Holzapfel-Ogden material model

        // Formulation used by cardiac mechanics benchmark paper (Arostica et al., 2024)
        // Uses I1_bar (bar = isochoric), but I4_f, I4_s, I8_fs (not bar)
        // Psi = a/2b * exp{b(I1_bar - 3)} 
        //       + a_f/2b_f * chi(I4_f) * (exp{b_f(I4_f - 1)^2} - 1
        //       + a_s/2b_s * chi(I4_s) * (exp{b_s(I4_s - 1)^2} - 1
        //       + a_fs/2b_fs * (exp{b_fs*I8_fs^2} - 1)
        // This corresponds to the HO-ma (modified anisotropy) implementation in svMultiPhysics

        // Invariants
        double I1_bar = smTerms.Ib1;
        // I4_f = f . C . f
        auto C_f = mat_fun::mat_mul(smTerms.C, f);
        double I4_f = f.dot(C_f);
        // I4_s = s . C . s
        auto C_s = mat_fun::mat_mul(smTerms.C, s);
        double I4_s = s.dot(C_s);
        // I8_fs = f . C . s
        double I8_fs = f.dot(C_s);

        // Strain energy density for Holzapfel-Ogden material model with full anisotropic invariants
        double Psi = 0.0;
        Psi += a / (2.0 * b) * exp(b * (I1_bar - 3.0));                             // Isotropic term
        Psi += a_f / (2.0 * b_f) * chi(I4_f, k) * (exp(b_f * pow(I4_f - 1.0, 2)) - 1.0);   // Fiber term
        Psi += a_s / (2.0 * b_s) * chi(I4_s, k) * (exp(b_s * pow(I4_s - 1.0, 2)) - 1.0);   // Sheet term
        Psi += a_fs / (2.0 * b_fs) * (exp(b_fs * pow(I8_fs, 2)) - 1.0);                   // Cross-fiber term

        return Psi;

    }
};


/**
 * @brief Class for testing the quadratic volumetric penalty model.
 *
 * This class provides methods to set up and test the quadratic volumetric penalty model, including 
 * computing the strain energy and printing material parameters.
 */
class TestQuadraticVolumetricPenalty : public TestMaterialModel {
public:

    /**
     * @brief Parameters for the volumetric penalty model.
     */
    VolumetricPenaltyParams params;

    /**
     * @brief Constructor for the TestQuadraticVolumetricPenalty class.
     *
     * Initializes the volumetric penalty parameters for svMultiPhysics.
     *
     * @param[in] params_ Parameters for the volumetric penalty model.
     */
    TestQuadraticVolumetricPenalty(const VolumetricPenaltyParams &params_) : TestMaterialModel( consts::ConstitutiveModelType::stIso_nHook, consts::ConstitutiveModelType::stVol_Quad),
        params(params_) 
        {

        // Set volumetric penalty parameter for svMultiPhysics
        auto &dmn = com_mod.mockEq.mockDmn;
        dmn.stM.Kpen = params.kappa;         // Volumetric penalty parameter

        // Note: Use Neo-Hookean material model for isochoric part, but set parameters to zero
        dmn.stM.C10 = 0.0;         // Zero Neo-Hookean parameter
    }

    /**
     * @brief Prints the volumetric penalty parameters.
     */
    void printMaterialParameters() {
        std::cout << "kappa = " << params.kappa << std::endl;
    }

    /**
     * @brief Computes the strain energy for the quadratic volumetric penalty model.
     *
     * @param[in] F Deformation gradient.
     * @return Strain energy density for the quadratic volumetric penalty model.
     */
    double computeStrainEnergy(const Array<double> &F) {
            
            // Compute solid mechanics terms
            solidMechanicsTerms smTerms = calcSolidMechanicsTerms(F);
    
            // Strain energy density for quadratic volumetric penalty model
            // Psi = kappa/2 * (J - 1)^2  
            double Psi = params.kappa/2.0 * pow(smTerms.J - 1.0, 2);
    
            return Psi;
    }
};

/**
 * @brief Class for testing the Simo-Taylor91 volumetric penalty model.
 *
 * This class provides methods to set up and test the Simo-Taylor91 volumetric penalty model, including 
 * computing the strain energy and printing material parameters.
 */
class TestSimoTaylor91VolumetricPenalty : public TestMaterialModel {
public:

    /**
     * @brief Parameters for the volumetric penalty model.
     */
    VolumetricPenaltyParams params;

    /**
     * @brief Constructor for the TestSimoTaylor91VolumetricPenalty class.
     *
     * Initializes the volumetric penalty parameters for svMultiPhysics.
     *
     * @param[in] params_ Parameters for the volumetric penalty model.
     */
    TestSimoTaylor91VolumetricPenalty(const VolumetricPenaltyParams &params_) : TestMaterialModel( consts::ConstitutiveModelType::stIso_nHook, consts::ConstitutiveModelType::stVol_ST91),
        params(params_) 
        {

        // Set volumetric penalty parameter for svMultiPhysics
        auto &dmn = com_mod.mockEq.mockDmn;
        dmn.stM.Kpen = params.kappa;         // Volumetric penalty parameter

        // Note: Use Neo-Hookean material model for isochoric part, but set parameters to zero
        dmn.stM.C10 = 0.0;         // Zero Neo-Hookean parameter
    }

    /**
     * @brief Prints the volumetric penalty parameters.
     */
    void printMaterialParameters() {
        std::cout << "kappa = " << params.kappa << std::endl;
    }

    /**
     * @brief Computes the strain energy for the Simo-Taylor91 volumetric penalty model.
     *
     * @param[in] F Deformation gradient.
     * @return Strain energy density for the Simo-Taylor91 volumetric penalty model.
     */
    double computeStrainEnergy(const Array<double> &F) {
            
            // Compute solid mechanics terms
            solidMechanicsTerms smTerms = calcSolidMechanicsTerms(F);
    
            // Strain energy density for Simo-Taylor91 volumetric penalty model
            // Psi = kappa/4 * (J^2 - 1 - 2*ln(J))
            double Psi = params.kappa/4.0 * (pow(smTerms.J, 2) - 1.0 - 2.0 * log(smTerms.J));
    
            return Psi;
    }
};

/**
 * @brief Class for testing the Miehe94 volumetric penalty model.
 *
 * This class provides methods to set up and test the Miehe94 volumetric penalty model, including 
 * computing the strain energy and printing material parameters.
 */
class TestMiehe94VolumetricPenalty : public TestMaterialModel {
public:

    /**
     * @brief Parameters for the volumetric penalty model.
     */
    VolumetricPenaltyParams params;

    /**
     * @brief Constructor for the TestMiehe94VolumetricPenalty class.
     *
     * Initializes the volumetric penalty parameters for svMultiPhysics.
     *
     * @param[in] params_ Parameters for the volumetric penalty model.
     */
    TestMiehe94VolumetricPenalty(const VolumetricPenaltyParams &params_) : TestMaterialModel( consts::ConstitutiveModelType::stIso_nHook, consts::ConstitutiveModelType::stVol_M94),
        params(params_) 
        {

        // Set volumetric penalty parameter for svMultiPhysics
        auto &dmn = com_mod.mockEq.mockDmn;
        dmn.stM.Kpen = params.kappa;         // Volumetric penalty parameter

        // Note: Use Neo-Hookean material model for isochoric part, but set parameters to zero
        dmn.stM.C10 = 0.0;         // Zero Neo-Hookean parameter
    }

    /**
     * @brief Prints the volumetric penalty parameters.
     */
    void printMaterialParameters() {
        std::cout << "kappa = " << params.kappa << std::endl;
    }

    /**
     * @brief Computes the strain energy for the Miehe94 volumetric penalty model.
     *
     * @param[in] F Deformation gradient.
     * @return Strain energy density for the Miehe94 volumetric penalty model.
     */
    double computeStrainEnergy(const Array<double> &F) {
            
            // Compute solid mechanics terms
            solidMechanicsTerms smTerms = calcSolidMechanicsTerms(F);
    
            // Strain energy density for Miehe94 volumetric penalty model
            // Psi = kappa * (J - ln(J) - 1)
            double Psi = params.kappa * (smTerms.J - log(smTerms.J) - 1.0);
    
            return Psi;
    }
};
