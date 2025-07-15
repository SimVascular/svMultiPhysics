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

#include "test_material_mooney_rivlin.h"

// ----------------------------------------------------------------------------
// --------------------------- Mooney-Rivlin Material -------------------------
// ----------------------------------------------------------------------------

/**
 * @brief  Test fixture class for the Mooney-Rivlin material model.
 * 
 * This class sets up the necessary parameters and objects for testing the Mooney-Rivlin material model.
 */
class MooneyRivlinTest : public MaterialTestFixture {
protected:
    // Material parameters object
    MooneyRivlinParams params;

    // Add the test object
    TestMooneyRivlin* TestMR;

    // Setup method to initialize variables before each test
    void SetUp() override {

        MaterialTestFixture::SetUp();

        // Set random values for the Mooney-Rivlin parameters between 1000 and 10000
        params.C01 = getRandomDouble(1000.0, 10000.0);
        params.C10 = getRandomDouble(1000.0, 10000.0);

        // Initialize the test object
        TestMR = new TestMooneyRivlin(params);
    }

    // TearDown method to clean up after each test, if needed
    void TearDown() override {
        // Clean up the test object
        delete TestMR;
    }
};

/**
 * @brief Test fixture class for STRUCT Mooney-Rivlin material model.
 */
class STRUCT_MooneyRivlinTest : public MooneyRivlinTest {
protected:
    void SetUp() override {
        MooneyRivlinTest::SetUp();

        // Use struct
        TestMR->ustruct = false;
    }
};

/**
 * @brief Test fixture class for USTRUCT Mooney-Rivlin material model.
 */
class USTRUCT_MooneyRivlinTest : public MooneyRivlinTest {
protected:
    void SetUp() override {
        MooneyRivlinTest::SetUp();

        // Use ustruct
        TestMR->ustruct = true;
    }
};

// ------------------------------ STRUCT Tests --------------------------------

// Test PK2 stress zero for F = I
TEST_F(STRUCT_MooneyRivlinTest, TestPK2StressIdentityF) {
    //verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.0, 0.0, 0.0},
                        {0.0, 1.0, 0.0},
                        {0.0, 0.0, 1.0}};
    Array<double> S_ref(3, 3); // PK2 stress initialized to zero
    TestMR->testPK2StressAgainstReference(F, S_ref, REL_TOL, ABS_TOL, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (small)
TEST_F(STRUCT_MooneyRivlinTest, TestPK2StressConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_small_list) {

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestMR->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (medium)
TEST_F(STRUCT_MooneyRivlinTest, TestPK2StressConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_medium_list) {

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestMR->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (large)
TEST_F(STRUCT_MooneyRivlinTest, TestPK2StressConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_large_list) {

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestMR->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (small)
TEST_F(STRUCT_MooneyRivlinTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_small_list) {

        // Check order of convergence of consistency of material elasticity
        TestMR->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (medium)
TEST_F(STRUCT_MooneyRivlinTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_medium_list) {

        // Check order of convergence of consistency of material elasticity
        TestMR->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (large)
TEST_F(STRUCT_MooneyRivlinTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_large_list) {

        // Check order of convergence of consistency of material elasticity
        TestMR->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}


// ------------------------------ USTRUCT Tests --------------------------------

// Test PK2 stress zero for F = I
TEST_F(USTRUCT_MooneyRivlinTest, TestPK2StressIdentityF) {
    //verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.0, 0.0, 0.0},
                        {0.0, 1.0, 0.0},
                        {0.0, 0.0, 1.0}};
    Array<double> S_ref(3, 3); // PK2 stress initialized to zero
    TestMR->testPK2StressAgainstReference(F, S_ref, REL_TOL, ABS_TOL, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (small)
TEST_F(USTRUCT_MooneyRivlinTest, TestPK2StressConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_small_list) {

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestMR->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (medium)
TEST_F(USTRUCT_MooneyRivlinTest, TestPK2StressConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_medium_list) {

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestMR->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (large)
TEST_F(USTRUCT_MooneyRivlinTest, TestPK2StressConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_large_list) {

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestMR->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (small)
TEST_F(USTRUCT_MooneyRivlinTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_small_list) {

        // Check order of convergence of consistency of material elasticity
        TestMR->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (medium)
TEST_F(USTRUCT_MooneyRivlinTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_medium_list) {

        // Check order of convergence of consistency of material elasticity
        TestMR->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (large)
TEST_F(USTRUCT_MooneyRivlinTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_large_list) {

        // Check order of convergence of consistency of material elasticity
        TestMR->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}
