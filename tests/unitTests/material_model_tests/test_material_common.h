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

#ifndef TEST_MATERIAL_COMMON_H
#define TEST_MATERIAL_COMMON_H

#include "../test_common.h"
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
inline Array<double> create_identity_F(const int N) {
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
inline Array<double> create_ones_matrix(const int N) {
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
inline double getRandomDouble(const double min, const double max) {
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
inline Array<double> create_random_F(const int N, const double min=0.1, const double max=10.0) {
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
inline Array<double> create_random_perturbed_identity_F(const int N, const double max_deviation) {
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
inline Array<double> perturb_random_F(const Array<double> &F, const double delta) {

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
inline void calc_JCE(const Array<double> &F, double &J, Array<double> &C, Array<double> &E) {

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
inline solidMechanicsTerms calcSolidMechanicsTerms(const Array<double> &F) {

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
inline std::pair<double, double> computeLinearRegression(const std::vector<double>& x, const std::vector<double>& y) {
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
// ------------------ Base Material Model Testing Class ---------
// --------------------------------------------------------------

// Class for testing material models in svMultiPhysics
class TestMaterialModel : public TestBase {
public:
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
     * @brief Computes the PK2 stress tensor S(F) from the strain energy density Psi(F) using finite differences and checks the order of convergence using a reference S_ref for exact solution.
     * Using this to compare CANN PK2 stress against other models.
     * @param[in] delta_min Minimum perturbation scaling factor.
     * @param[in] delta_max Maximum perturbation scaling factor.
     * @param[in] order Order of the finite difference scheme (1 for first order, 2 for second order, etc.).
     * @param[in] convergence_order_tol Tolerance for comparing convergence order with expected value
     * @param[in] verbose Show values error and order of convergence if true.
     **/
    
    void testPK2StressConvergenceOrderAgainstReference(const Array<double>& F, const Array<double>& S_ref, const double delta_max, const double delta_min, const int order, const double convergence_order_tol, const bool verbose = false) {
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
        // Create list of deltas for convergence test (delta = delta_max, delta_max/2, delta_max/4, ...)
        std::vector<double> deltas;
        double delta = delta_max;
        while (delta >= delta_min) {
            deltas.push_back(delta);
            delta /= 2.0;
        }
        // Compute finite difference S for each delta and store error in list
        std::vector<double> errors;
        Array<double> S_fd;
        for (int i = 0; i < deltas.size(); i++) {
            calcPK2StressFiniteDifference(F, deltas[i], order, S_fd);

            // Compute Frobenius norm of error between S and S_fd
            double error = 0.0;
            for (int I = 0; I < 3; I++) {
                for (int J = 0; J < 3; J++) {
                    error += pow(S_ref(I,J) - S_fd(I,J), 2);
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

            // std::cout << "F = " << std::endl;
            // for (int i = 0; i < 3; i++) {
            //     for (int J = 0; J < 3; J++) {
            //         std::cout << F[i][J] << " ";
            //     }
            //     std::cout << std::endl;
            // }

            // std::cout << "S = " << std::endl;
            // for (int i = 0; i < 3; i++) {
            //     for (int J = 0; J < 3; J++) {
            //         std::cout << S[i][J] << " ";
            //     }
            //     std::cout << std::endl;
            // }

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
     * @brief Tests the order of convergence of the consistency between CC:dE and dS using finite differences, with CCdE and dS as INPUTS. Used to compare CANN with other models
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
     * @param[in] dS Change in pk2 due to perturbation dF in F
     * @param[in] CCdE CC:dE
     * @param[in] deltas scaling factors for perturbations
     * @param[in] order Order of the finite difference scheme (1 for first order, 2 for second order, etc.).
     * @param[in] convergence_order_tol Tolerance for comparing convergence order with expected value
     * @param[in] verbose Show values of errors and order of convergence if true.
     */
    void testMaterialElasticityConsistencyConvergenceOrderBetweenMaterialModels(const Array<double>& F, Array<double>& dS, Array<double>& CCdE, std::vector<double> deltas, int order, const double convergence_order_tol, bool verbose = false) {

        // Loop over perturbations to each component of F, dF
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                std::vector<double> errors;
                // Compute Frobenius norm of error between dS and CC:dE
                for (int i = 0; i < deltas.size(); i++) {
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
                if (std::isnan(m)) {
                    std::ostringstream oss;
                    oss << "Error: m value nan. "
                        << ", F = [";

                    // Append each element of F to the string stream
                    for (int i = 0; i < 3; ++i) {
                        for (int j = 0; j < 3; ++j) {
                            oss << F(i,j);
                            if (j < 3 - 1) oss << ", ";
                        }
                        if (i < 3 - 1) oss << "; ";
                    }
                    oss << "]";
                    throw std::runtime_error(oss.str());
                }
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
     * @brief Generate perturbation dF in F
     * 
     * @param[in] F Deformation gradient.
     * @param[in] dF Perturbation in deformation gradient.
     * @param[in] delta_max Maximum perturbation scaling factor.
     * @param[in] delta_min Minimum perturbation scaling factor.
     * @param[in] deltas scaling factors for perturbations
     * @param[in] order Order of the finite difference scheme (1 for first order, 2 for second order, etc.).
     * @param[in] verbose Show values of errors and order of convergence if true.
     */
    void generatePerturbationdF(const Array<double>& F, Array<double>& dF, double delta_max, double delta_min, std::vector<double> deltas, int order, bool verbose=false) {
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
        double delta = delta_max;
        while (delta >= delta_min) {
            deltas.push_back(delta);
            delta /= 2.0;
        }
        // Loop over perturbations to each component of F, dF
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                // Generate dF with 1.0 in (i,j) component
                dF(i,j) = 1.0;

                // Print results if verbose
                if (verbose) {
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
     * @brief Calculate the reference material elasticity tensor CC(F) for comparison with Neo-Hookean.
     *
     * This function computes the material elasticity tensor CC(F) from the deformation gradient F using compute_pk2cc() 
     *
     * @param[in] F Deformation gradient.
     * @param[in] CC_ref Reference solution for material elasticity tensor.
     * @param[in] verbose Show values of F, CC, and CC_ref if true.
     * @return None.
     */
    void calcMaterialElasticityReference(const Array<double> &F, Tensor4<double> &CC_ref, bool verbose = false) {
        int N = F.ncols();
        assert(F.nrows() == N);

        // Compute CC(F) from compute_pk2cc()
        Array<double> S(N, N), Dm(2*N, 2*N);
        compute_pk2cc(F, S, Dm);

        // Calculate CC_ref from Dm
        Tensor4<double> CC(N, N, N, N);
        mat_models::voigt_to_cc(N, Dm, CC_ref);

        // mat_fun_carray::print("CC_ref from calc function",CC_ref);

        // Print results if verbose
        if (verbose) {
            printMaterialParameters();

            std::cout << "F =" << std::endl;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    std::cout << F(i,j) << " ";
                }
                std::cout << std::endl;
            }

            std::cout << "CC_ref =" << std::endl;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) {
                        for (int l = 0; l < 3; l++) {
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
// ------------------ Base Material Parameters Class ------------
// --------------------------------------------------------------

// Class to contain material parameters
class MatParams {
public:
    virtual ~MatParams() {} // Virtual destructor for proper cleanup
};

// --------------------------------------------------------------
// ------------------ Material Test Fixture Class ---------------
// --------------------------------------------------------------

/**
 * @brief Test fixture class containing common setup for all material model test
 * 
 */
class MaterialTestFixture : public ::testing::Test {
protected:
    // Variables common across tests
    static constexpr double DEFORMATION_PERTURBATION_SMALL = 0.003; // Small perturbation factor
    static constexpr double DEFORMATION_PERTURBATION_MEDIUM = 0.03; // Medium perturbation factor
    static constexpr double DEFORMATION_PERTURBATION_LARGE = 0.3; // Large perturbation factor
    static constexpr int n_F = 50; // Number of deformation gradients F to test for each small, medium, and large perturbation
    static constexpr double REL_TOL = 1e-3; // relative tolerance for comparing values
    static constexpr double ABS_TOL = 1e-11; // absolute tolerance for comparing values
    //double delta = 1e-7; // perturbation scaling factor
    static constexpr double DELTA_MAX = 1e-4; // maximum perturbation scaling factor
    static constexpr double DELTA_MIN = 1e-6; // minimum perturbation scaling factor
    static constexpr int ORDER = 1; // Order of finite difference method
    static constexpr double CONVERGENCE_ORDER_TOL = 0.02; // Tolerance for comparing convergence order with expected value
    
    bool verbose = false; // Show values of S, dE, SdE and dPsi

    // Vectors to store the Array<double> deformation gradients
    std::vector<Array<double>> F_small_list;
    std::vector<Array<double>> F_medium_list;
    std::vector<Array<double>> F_large_list;

    void SetUp() override {

        // Create random deformation gradients for small perturbations
        for (int i = 0; i < n_F; i++) {
            F_small_list.push_back(create_random_perturbed_identity_F(3, DEFORMATION_PERTURBATION_SMALL));
        }

        // Create random deformation gradients for medium perturbations
        for (int i = 0; i < n_F; i++) {
            F_medium_list.push_back(create_random_perturbed_identity_F(3, DEFORMATION_PERTURBATION_MEDIUM));
        }

        // Create random deformation gradients for large perturbations
        for (int i = 0; i < n_F; i++) {
            F_large_list.push_back(create_random_perturbed_identity_F(3, DEFORMATION_PERTURBATION_LARGE));
        }
    }

    void TearDown() override {}
};


#endif