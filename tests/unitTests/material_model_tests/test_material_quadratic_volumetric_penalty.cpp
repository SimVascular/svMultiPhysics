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

#include "test_material_quadratic_volumetric_penalty.h"

// ----------------------------------------------------------------------------
// ---------------- Quadratic Volumetric Penalty Model ------------------------
// ----------------------------------------------------------------------------
/**
 * @brief Test fixture class for the Quadratic Volumetric penalty model.
 * 
 * This class sets up the necessary parameters and objects for testing the Quadratic Volumetric penalty model.
 */
class QuadraticVolumetricPenaltyTest : public MaterialTestFixture {
    protected:
        // Material parameters object
        QuadraticVolumetricPenaltyParams params;
    
        // Add the test object
        TestQuadraticVolumetricPenalty* TestQVP;
    
        // Setup method to initialize variables before each test
        void SetUp() override {
    
            MaterialTestFixture::SetUp();
    
            // Set random values for the Quadratic penalty parameters between 1000 and 10000
            params.kappa = getRandomDouble(1000.0, 10000.0);
    
            // Initialize the test object
            TestQVP = new TestQuadraticVolumetricPenalty(params);
        }
    
        // TearDown method to clean up after each test, if needed
        void TearDown() override {
            // Clean up the test object
            delete TestQVP;
        }
    };
    
    /**
     * @brief Test fixture class for STRUCT Quadratic penalty model.
     */
    class STRUCT_QuadraticVolumetricPenaltyTest : public QuadraticVolumetricPenaltyTest {
    protected:
        void SetUp() override {
            QuadraticVolumetricPenaltyTest::SetUp();
    
            // Use struct
            //TestQVP->ustruct = false;
        }
    };
    
    /**
     * @brief Test fixture class for USTRUCT Quadratic penalty model.
     */
    class USTRUCT_QuadraticVolumetricPenaltyTest : public QuadraticVolumetricPenaltyTest {
    protected:
        void SetUp() override {
            QuadraticVolumetricPenaltyTest::SetUp();
    
            // Use ustruct
            //TestQVP->ustruct = true;
        }
    };

// ------------------------------ STRUCT Tests --------------------------------

// Test PK2 stress zero for F = I
TEST_F(STRUCT_QuadraticVolumetricPenaltyTest, TestPK2StressIdentityF) {
    //verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.0, 0.0, 0.0},
                        {0.0, 1.0, 0.0},
                        {0.0, 0.0, 1.0}};
    Array<double> S_ref(3, 3); // PK2 stress initialized to zero
    TestQVP->testPK2StressAgainstReference(F, S_ref, REL_TOL, ABS_TOL, verbose);
}

// Test PK2 stress zero for prescribed isochoric deformation
TEST_F(STRUCT_QuadraticVolumetricPenaltyTest, TestPK2StressPrescribedIsochoricDeformation) {
    //verbose = true; // Show values of S and S_ref

    // Check isochoric deformation produces zero PK2 stress
    Array<double> F = {{1.1, 0.0, 0.0},
                        {0.0, 1.2, 0.0},
                        {0.0, 0.0, 1.0/(1.1*1.2)}};
    Array<double> S_ref(3, 3); // PK2 stress initialized to zero
    TestQVP->testPK2StressAgainstReference(F, S_ref, REL_TOL, ABS_TOL, verbose);
}


// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (small)
TEST_F(STRUCT_QuadraticVolumetricPenaltyTest, TestPK2StressConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_small_list) {

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestQVP->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (medium)
TEST_F(STRUCT_QuadraticVolumetricPenaltyTest, TestPK2StressConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_medium_list) {

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestQVP->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (large)
TEST_F(STRUCT_QuadraticVolumetricPenaltyTest, TestPK2StressConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_large_list) {

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestQVP->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (small)
TEST_F(STRUCT_QuadraticVolumetricPenaltyTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_small_list) {

        // Check order of convergence of consistency of material elasticity
        TestQVP->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (medium)
TEST_F(STRUCT_QuadraticVolumetricPenaltyTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_medium_list) {

        // Check order of convergence of consistency of material elasticity
        TestQVP->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (large)
TEST_F(STRUCT_QuadraticVolumetricPenaltyTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_large_list) {

        // Check order of convergence of consistency of material elasticity
        TestQVP->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}


// ------------------------------ USTRUCT Tests --------------------------------

// Test rho, beta, drho/dp and dbeta/dp for random pressure
TEST_F(USTRUCT_QuadraticVolumetricPenaltyTest, TestRhoBeta) {
    //verbose = true; // Show values of rho, beta, drho/dp and dbeta/dp and reference values

    // Generate random pressure p between 100 and 1000
    double p = getRandomDouble(100.0, 1000.0);

    // Other parameters
    double rho0 = 1000.0; // Reference density
    double kappa = TestQVP->params.kappa; // Volumetric penalty parameter

    // Compute reference values for rho, beta, drho/dp and dbeta/dp
    // See ustruct paper (https://doi.org/10.1016/j.cma.2018.03.045) Section 2.4
    double rho_ref = rho0 / (1.0 - p/kappa);
    double beta_ref = 1.0 / (kappa - p);
    double drhodp_ref = -rho0 / pow(1.0 - p/kappa, 2) * (-1.0/kappa); // Derivative of rho with respect to p
    double dbetadp_ref = -1.0 / pow(kappa - p, 2) * (-1.0); // Derivative of beta with respect to p

    // Check rho, beta, drho/dp and dbeta/dp against reference values
    TestQVP->testRhoBetaAgainstReference(p, rho0, rho_ref, beta_ref, drhodp_ref, dbetadp_ref, REL_TOL, ABS_TOL, verbose);
}