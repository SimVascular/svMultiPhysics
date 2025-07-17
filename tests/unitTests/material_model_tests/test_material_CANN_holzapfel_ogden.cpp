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

#include "test_material_CANN_holzapfel_ogden.h"
#include "test_material_holzapfel_ogden.h"



// ----------------------------------------------------------------------------
// -------- Compare CANN w/ HO params against HO Model ------------------------
// ----------------------------------------------------------------------------
/**
 * @brief Test fixture class for comparing the CANN model with Holzapfel Ogden mdoel.
 *
 * This class sets up the necessary parameters and objects for comparing CANN model with HO parameters against the Holzapfel-Ogden model. Assuming k = 0, so chi = 0.5.
 */
class HolzapfelOgdenCompareTest : public MaterialTestFixture {
protected:
    // Material parameters objects
    HolzapfelOgdenParams params_HO; // HO parameters
    CANN_HO_Params params_CANN_HO; // CANN with HO parameters

    // Add the test objects
    TestHolzapfelOgden* TestHO;
    TestCANN_HO* TestCANNHO;

    // Setup method to initialize variables before each test
    void SetUp() override {
        MaterialTestFixture::SetUp();

        // Set Holzapfel-Ogden parameters from cardiac benchmark paper
        params_HO.a = 59.0; // Pa
        params_HO.a_f = 18472.0; // Pa
        params_HO.a_s = 2481.0; // Pa
        params_HO.a_fs = 216.0; // Pa
        params_HO.b = 8.023; // no units
        params_HO.b_f = 16.026; // no units
        params_HO.b_s = 11.12; // no units
        params_HO.b_fs = 11.436; // no units
        // setting heaviside parameter to 0
        params_HO.k = 0.0; // no units

        // Set random values for f between 0 and 1 and normalize
        params_HO.f[0] = getRandomDouble(0.0, 1.0);
        params_HO.f[1] = getRandomDouble(0.0, 1.0);
        params_HO.f[2] = getRandomDouble(0.0, 1.0);
        double norm_f = sqrt(params_HO.f.dot(params_HO.f));
        params_HO.f = params_HO.f / norm_f;

        // Create s orthogonal to f
        if (fabs(params_HO.f[0]) < 0.9) { // Check if f[0] is not the dominant component
            params_HO.s[0] = 0;
            params_HO.s[1] = params_HO.f[2];
            params_HO.s[2] = -params_HO.f[1];
        } else { // If f[0] is the dominant component, use another approach
            params_HO.s[0] = -params_HO.f[2];
            params_HO.s[1] = 0;
            params_HO.s[2] = params_HO.f[0];
        }

        // Normalize s
        double norm_s = sqrt(params_HO.s.dot(params_HO.s));
        params_HO.s = params_HO.s / norm_s;

        // Check f.s = 0
        double dot_fs = params_HO.f.dot(params_HO.s);
        if (fabs(dot_fs) > 1e-6) {
            std::cout << "f.s = " << dot_fs << std::endl;
            std::cout << "f = [" << params_HO.f[0] << ", " << params_HO.f[1] << ", " << params_HO.f[2] << "]" << std::endl;
            std::cout << "s = [" << params_HO.s[0] << ", " << params_HO.s[1] << ", " << params_HO.s[2] << "]" << std::endl;
            throw std::runtime_error("f and s are not orthogonal");
        }
        
        // Initializing paramter table
        params_CANN_HO.Table.resize(4);  // Ensure it has 4 entries

        params_CANN_HO.Table[0].invariant_index.value_ = 1;
        params_CANN_HO.Table[0].activation_functions.value_ = {1,1,2};
        params_CANN_HO.Table[0].weights.value_ = {1.0, params_HO.b, params_HO.a / (2 * params_HO.b)};

        params_CANN_HO.Table[1].invariant_index.value_ = 4;
        params_CANN_HO.Table[1].activation_functions.value_ = {1,2,2};
        params_CANN_HO.Table[1].weights.value_ = {1.0, params_HO.b_f, 0.5 * params_HO.a_f / (2 * params_HO.b_f)}; // 0.5 for heaviside func

        params_CANN_HO.Table[2].invariant_index.value_ = 8;
        params_CANN_HO.Table[2].activation_functions.value_ = {1,2,2};
        params_CANN_HO.Table[2].weights.value_ = {1.0, params_HO.b_s, 0.5 * params_HO.a_s / (2 * params_HO.b_s)};

        params_CANN_HO.Table[3].invariant_index.value_ = 6;
        params_CANN_HO.Table[3].activation_functions.value_ = {1,2,2};
        params_CANN_HO.Table[3].weights.value_ = {1.0, params_HO.b_fs, params_HO.a_fs / (2 * params_HO.b_fs)};

        // initializing fiber directions for CANN
        params_CANN_HO.f[0] = params_HO.f[0]; params_CANN_HO.f[1] = params_HO.f[1]; params_CANN_HO.f[2] = params_HO.f[2];
        params_CANN_HO.s[0] = params_HO.s[0]; params_CANN_HO.s[1] = params_HO.s[1]; params_CANN_HO.s[2] = params_HO.s[2];

        // Initialize the test objects
        TestHO = new TestHolzapfelOgden(params_HO);
        TestCANNHO = new TestCANN_HO(params_CANN_HO);
    }

    // TearDown method to clean up after each test, if needed
    void TearDown() override {
        // Clean up the test objects
        delete TestHO;
        delete TestCANNHO;
    }
};

/**
 * @brief Test fixture class for STRUCT Holzapfel-Ogden Compare model.
 */
class STRUCT_CANNHolzapfelOgdenTest : public HolzapfelOgdenCompareTest {
protected:
    void SetUp() override {
        HolzapfelOgdenCompareTest::SetUp();

        // Use struct
        TestHO->ustruct = false;
        TestCANNHO->ustruct = false;
    }
};

/**
 * @brief Test fixture class for USTRUCT Neo-Hookean Compare model.
 */
class USTRUCT_CANNHolzapfelOgdenTest : public HolzapfelOgdenCompareTest {
protected:
    void SetUp() override {
        HolzapfelOgdenCompareTest::SetUp();

        // Use ustruct
        TestHO->ustruct = true;
        TestCANNHO->ustruct = true;
    }
};


// ------------------------------ STRUCT Tests --------------------------------

// Test PK2 stress zero for F = I
TEST_F(STRUCT_CANNHolzapfelOgdenTest, TestPK2StressIdentityF) {
    verbose = false; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.0, 0.0, 0.0},
                       {0.0, 1.0, 0.0},
                       {0.0, 0.0, 1.0}};
    Array<double> S_ref(3,3); // PK2 stress initialized to zero - want to get result from NH and set that to S_ref
    Array<double> Dm(6,6);
    TestHO->compute_pk2cc(F,S_ref,Dm); // Computing S_ref from NH
    TestCANNHO->testPK2StressAgainstReference(F, S_ref, REL_TOL, ABS_TOL, verbose); // Comparing with CANN
}

// Test PK2 stress
TEST_F(STRUCT_CANNHolzapfelOgdenTest, TestPK2StressTriaxialStretch) {
    verbose = false; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.1, 0.0, 0.0},
                       {0.0, 1.2, 0.0},
                       {0.0, 0.0, 0.757}};
    Array<double> S_ref(3,3); // PK2 stress initialized to zero - want to get result from NH and set that to S_ref
    Array<double> Dm(6,6);
    TestHO->compute_pk2cc(F,S_ref,Dm); // Computing S_ref from NH
    TestCANNHO->testPK2StressAgainstReference(F, S_ref, REL_TOL, ABS_TOL, verbose); // Comparing with CANN
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (small)
TEST_F(STRUCT_CANNHolzapfelOgdenTest, TestPK2StressConvergenceOrderAgainstReferenceRandomFSmall) {
    verbose = false; // Show order of convergence, errors, F, S

    for (const auto& F : F_small_list) {

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        Array<double> S_ref(3,3), Dm(6,6);
        TestCANNHO->compute_pk2cc(F,S_ref,Dm); // Computing S_ref from CANN
        TestHO->testPK2StressConvergenceOrderAgainstReference(F, S_ref, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (medium)
TEST_F(STRUCT_CANNHolzapfelOgdenTest, TestPK2StressConvergenceOrderAgainstReferenceRandomFMedium) {
    verbose = false; // Show order of convergence, errors, F, S

    for (const auto& F : F_medium_list) {

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        Array<double> S_ref(3,3), Dm(6,6);
        TestCANNHO->compute_pk2cc(F,S_ref,Dm); // Computing S_ref from CANN
        TestHO->testPK2StressConvergenceOrderAgainstReference(F, S_ref, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (large)
TEST_F(STRUCT_CANNHolzapfelOgdenTest, TestPK2StressConvergenceOrderAgainstReferenceRandomFLarge) {
    verbose = false; // Show order of convergence, errors, F, S

    for (const auto& F : F_large_list) {

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        Array<double> S_ref(3,3), Dm(6,6);
        TestCANNHO->compute_pk2cc(F,S_ref,Dm); // Computing S_ref from CANN
        TestHO->testPK2StressConvergenceOrderAgainstReference(F,S_ref, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (small)
TEST_F(STRUCT_CANNHolzapfelOgdenTest, TestMaterialElasticityConsistencyConvergenceOrderAgainstHORandomFSmall) {
    verbose = false; // Show order of convergence, errors, F, S

    for (const auto& F : F_small_list) {
        // Generating perturbation
        Array<double> dF(3,3);
        std::vector<double> deltas;
        TestHO->generatePerturbationdF(F,dF,DELTA_MAX, DELTA_MIN, deltas,ORDER,verbose);
        // Calculating dS and CCdE
        Array<double> dS(3,3), CCdE(3,3);
        for (int i = 0; i < deltas.size(); i++) {
            TestCANNHO->calcCCdEFiniteDifference(F, dF, deltas[i], ORDER, CCdE);
            TestHO->calcdSFiniteDifference(F, dF, deltas[i], ORDER, dS);

            // Check order of convergence of consistency of material elasticity
            TestHO->testMaterialElasticityConsistencyConvergenceOrderBetweenMaterialModels(F, dS, CCdE, deltas, ORDER, CONVERGENCE_ORDER_TOL, verbose);
        }
    }
}

// Test order of convergence of consistency of material elasticity for random F (medium)
TEST_F(STRUCT_CANNHolzapfelOgdenTest, TestMaterialElasticityConsistencyConvergenceOrderAgainstHORandomFMedium) {
    verbose = false; // Show order of convergence, errors, F, S

    for (const auto& F : F_medium_list) {
        // Generating perturbation
        Array<double> dF(3,3);
        std::vector<double> deltas;
        TestHO->generatePerturbationdF(F,dF,DELTA_MAX, DELTA_MIN, deltas,ORDER,verbose);
        // Calculating dS and CCdE
        Array<double> dS(3,3), CCdE(3,3);
        for (int i = 0; i < deltas.size(); i++) {
            TestCANNHO->calcCCdEFiniteDifference(F, dF, deltas[i], ORDER, CCdE);
            TestHO->calcdSFiniteDifference(F, dF, deltas[i], ORDER, dS);

            // Check order of convergence of consistency of material elasticity
            TestHO->testMaterialElasticityConsistencyConvergenceOrderBetweenMaterialModels(F, dS, CCdE, deltas, ORDER, CONVERGENCE_ORDER_TOL, verbose);
        }
    }
}

// Test order of convergence of consistency of material elasticity for random F (large)
TEST_F(STRUCT_CANNHolzapfelOgdenTest, TestMaterialElasticityConsistencyConvergenceOrderAgainstHORandomFLarge) {
    verbose = false; // Show order of convergence, errors, F, S

    for (const auto& F : F_large_list) {
        // Generating perturbation
        Array<double> dF(3,3);
        std::vector<double> deltas;
        TestHO->generatePerturbationdF(F,dF,DELTA_MAX, DELTA_MIN, deltas,ORDER,verbose);
        // Calculating dS and CCdE
        Array<double> dS(3,3), CCdE(3,3);
        for (int i = 0; i < deltas.size(); i++) {
            TestCANNHO->calcCCdEFiniteDifference(F, dF, deltas[i], ORDER, CCdE);
            TestHO->calcdSFiniteDifference(F, dF, deltas[i], ORDER, dS);

            // Check order of convergence of consistency of material elasticity
            TestHO->testMaterialElasticityConsistencyConvergenceOrderBetweenMaterialModels(F, dS, CCdE, deltas, ORDER, CONVERGENCE_ORDER_TOL, verbose);
        }
    }
}

// Test PK2 stress
TEST_F(STRUCT_CANNHolzapfelOgdenTest, TestMaterialElasticityAgainstReferenceIdentityF) {
    verbose = false; // Show values of CC and CC_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.0, 0.0, 0.0},
                       {0.0, 1.0, 0.0},
                       {0.0, 0.0, 1.0}};
    Tensor4<double> CC_ref(3,3,3,3); // CC_ref initialized to zero
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            for (int k = 0; k < 3; k++){
                for (int l = 0; l < 3; l++){
                    CC_ref(i,j,k,l) = 0;
                } 
            }
        }
    }
    TestHO->calcMaterialElasticityReference(F,CC_ref,verbose);
    TestCANNHO->testMaterialElasticityAgainstReference(F, CC_ref, REL_TOL, ABS_TOL, verbose); // Comparing with CANN
}

// Test PK2 stress
TEST_F(STRUCT_CANNHolzapfelOgdenTest, TestMaterialElasticityAgainstReference) {
    verbose = false; // Show values of CC and CC_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.1, 0.0, 0.0},
                       {0.0, 1.2, 0.0},
                       {0.0, 0.0, 0.757}};
    Tensor4<double> CC_ref(3,3,3,3); // CC_ref initialized to zero
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            for (int k = 0; k < 3; k++){
                for (int l = 0; l < 3; l++){
                    CC_ref(i,j,k,l) = 0;
                } 
            }
        }
    }
    TestHO->calcMaterialElasticityReference(F,CC_ref,verbose);
    TestCANNHO->testMaterialElasticityAgainstReference(F, CC_ref, REL_TOL, ABS_TOL, verbose); // Comparing with CANN
}







