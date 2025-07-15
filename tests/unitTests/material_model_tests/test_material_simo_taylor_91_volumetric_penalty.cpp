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

#include "test_material_simo_taylor_91_volumetric_penalty.h"

// ----------------------------------------------------------------------------
// --------------------------- Simo-Taylor 91 Volumetric Penalty --------------
// ----------------------------------------------------------------------------


/**
 * @brief Test fixture class for the Simo-Taylor91 Volumetric penalty model.
 * 
 * This class sets up the necessary parameters and objects for testing the Simo-Taylor91 Volumetric penalty model.
 */
class SimoTaylor91VolumetricPenaltyTest : public MaterialTestFixture {
    protected:
        // Material parameters object
        SimoTaylor91VolumetricPenaltyParams params;
    
        // Add the test object
        TestSimoTaylor91VolumetricPenalty* TestST91;
    
        // Setup method to initialize variables before each test
        void SetUp() override {
    
            MaterialTestFixture::SetUp();
    
            // Set random values for the Simo-Taylor91 penalty parameters between 1000 and 10000
            params.kappa = getRandomDouble(1000.0, 10000.0);
    
            // Initialize the test object
            TestST91 = new TestSimoTaylor91VolumetricPenalty(params);
        }
    
        // TearDown method to clean up after each test, if needed
        void TearDown() override {
            // Clean up the test object
            delete TestST91;
        }
    };
    
    /**
     * @brief Test fixture class for STRUCT Simo-Taylor91 penalty model.
     */
    class STRUCT_SimoTaylor91VolumetricPenaltyTest : public SimoTaylor91VolumetricPenaltyTest {
    protected:
        void SetUp() override {
            SimoTaylor91VolumetricPenaltyTest::SetUp();
    
            // Use struct
            //TestST91->ustruct = false;
        }
    };
    
    /**
     * @brief Test fixture class for USTRUCT Simo-Taylor91 penalty model.
     */
    class USTRUCT_SimoTaylor91VolumetricPenaltyTest : public SimoTaylor91VolumetricPenaltyTest {
    protected:
        void SetUp() override {
            SimoTaylor91VolumetricPenaltyTest::SetUp();
    
            // Use ustruct
            //TestST91->ustruct = true;
        }
    };

// ------------------------------ STRUCT Tests --------------------------------


// Test PK2 stress zero for F = I
TEST_F(STRUCT_SimoTaylor91VolumetricPenaltyTest, TestPK2StressIdentityF) {
    //verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.0, 0.0, 0.0},
                        {0.0, 1.0, 0.0},
                        {0.0, 0.0, 1.0}};
    Array<double> S_ref(3, 3); // PK2 stress initialized to zero
    TestST91->testPK2StressAgainstReference(F, S_ref, REL_TOL, ABS_TOL, verbose);
}

// Test PK2 stress zero for prescribed isochoric deformation
TEST_F(STRUCT_SimoTaylor91VolumetricPenaltyTest, TestPK2StressPrescribedIsochoricDeformation) {
    //verbose = true; // Show values of S and S_ref

    // Check isochoric deformation produces zero PK2 stress
    Array<double> F = {{1.1, 0.0, 0.0},
                        {0.0, 1.2, 0.0},
                        {0.0, 0.0, 1.0/(1.1*1.2)}};
    Array<double> S_ref(3, 3); // PK2 stress initialized to zero
    TestST91->testPK2StressAgainstReference(F, S_ref, REL_TOL, ABS_TOL, verbose);
}


// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (small)
TEST_F(STRUCT_SimoTaylor91VolumetricPenaltyTest, TestPK2StressConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_small_list) {

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestST91->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (medium)
TEST_F(STRUCT_SimoTaylor91VolumetricPenaltyTest, TestPK2StressConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_medium_list) {

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestST91->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (large)
TEST_F(STRUCT_SimoTaylor91VolumetricPenaltyTest, TestPK2StressConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_large_list) {

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestST91->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (small)
TEST_F(STRUCT_SimoTaylor91VolumetricPenaltyTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_small_list) {

        // Check order of convergence of consistency of material elasticity
        TestST91->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (medium)
TEST_F(STRUCT_SimoTaylor91VolumetricPenaltyTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_medium_list) {

        // Check order of convergence of consistency of material elasticity
        TestST91->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (large)
TEST_F(STRUCT_SimoTaylor91VolumetricPenaltyTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_large_list) {

        // Check order of convergence of consistency of material elasticity
        TestST91->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// ------------------------------ USTRUCT Tests --------------------------------

// Test rho and beta values for random p
TEST_F(USTRUCT_SimoTaylor91VolumetricPenaltyTest, TestRhoBeta) {
    //verbose = true; // Show values of rho, beta and rho_ref

    // Generate random pressure p between 100 and 1000
    double p = getRandomDouble(100.0, 1000.0);

    // Other parameters
    double rho0 = 1000.0; // Reference density
    double kappa = TestST91->params.kappa; // Volumetric penalty parameter

    // Compute reference values for rho, beta, drho/dp and dbeta/dp
    // See ustruct paper (https://doi.org/10.1016/j.cma.2018.03.045) Section 2.4
    double rho_ref = rho0 / kappa * (sqrt(p*p + kappa*kappa) + p);
    double beta_ref = 1.0 / sqrt(p*p + kappa*kappa);
    double drhodp_ref = rho0 / kappa * (1.0/2.0 * (1.0/sqrt(p*p + kappa*kappa)) * 2.0*p + 1.0); // Derivative of rho with respect to p
    double dbetadp_ref = -1.0/2.0 * 1.0/pow(p*p + kappa*kappa, 3.0/2.0) * 2.0*p; // Derivative of beta with respect to p

    // Check rho, beta, drho/dp and dbeta/dp against reference values
    TestST91->testRhoBetaAgainstReference(p, rho0, rho_ref, beta_ref, drhodp_ref, dbetadp_ref, REL_TOL, ABS_TOL, verbose);
}



