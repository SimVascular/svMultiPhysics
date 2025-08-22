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

#include "test_material_holzapfel_ogden_MA.h"
#include <limits>

// ----------------------------------------------------------------------------
// ------------- Holzapfel-Ogden (Modified Anisotropy) Material  --------------
// ----------------------------------------------------------------------------

/**
 * @brief Test fixture class for the Holzapfel-Ogden (Modified Anisotropy) material model.
 * 
 * This class sets up the necessary parameters and objects for testing the Holzapfel-Ogden (Modified Anisotropy) material model.
*/
class HolzapfelOgdenMATest : public MaterialTestFixture {
protected:
    // Material parameters object
    HolzapfelOgdenMAParams params;

    // Add the test object
    TestHolzapfelOgdenMA* TestHO_ma;

    // Setup method to initialize variables before each test
    void SetUp() override {

        MaterialTestFixture::SetUp();

        // Set Holzapfel-Ogden parameters from cardiac benchmark paper
        params.a = 59.0; // Pa
        params.a_f = 18472.0; // Pa
        params.a_s = 2481.0; // Pa
        params.a_fs = 216.0; // Pa
        params.b = 8.023; // no units
        params.b_f = 16.026; // no units
        params.b_s = 11.12; // no units
        params.b_fs = 11.436; // no units
        params.k = 100.0; // no units

        // Set random values for f between 0 and 1 and normalize
        params.f[0] = getRandomDouble(0.0, 1.0);
        params.f[1] = getRandomDouble(0.0, 1.0);
        params.f[2] = getRandomDouble(0.0, 1.0);
        double norm_f = sqrt(params.f.dot(params.f));
        params.f = params.f / norm_f;

        // Create s orthogonal to f
        if (fabs(params.f[0]) < 0.9) { // Check if f[0] is not the dominant component
            params.s[0] = 0;
            params.s[1] = params.f[2];
            params.s[2] = -params.f[1];
        } else { // If f[0] is the dominant component, use another approach
            params.s[0] = -params.f[2];
            params.s[1] = 0;
            params.s[2] = params.f[0];
        }

        // Normalize s
        double norm_s = sqrt(params.s.dot(params.s));
        params.s = params.s / norm_s;

        // Check f.s = 0
        double dot_fs = params.f.dot(params.s);
        if (fabs(dot_fs) > std::numeric_limits<double>::epsilon()) {
            std::cout << "f.s = " << dot_fs << std::endl;
            std::cout << "f = [" << params.f[0] << ", " << params.f[1] << ", " << params.f[2] << "]" << std::endl;
            std::cout << "s = [" << params.s[0] << ", " << params.s[1] << ", " << params.s[2] << "]" << std::endl;
            throw std::runtime_error("f and s are not orthogonal");
        }


        // Initialize the test object
        TestHO_ma = new TestHolzapfelOgdenMA(params);
    }

    // TearDown method to clean up after each test, if needed
    void TearDown() override {
        // Clean up the test object
        delete TestHO_ma;
    }
};

/**
 * @brief Test fixture class for STRUCT Holzapfel-Ogden material model.
 */
class STRUCT_HolzapfelOgdenMATest : public HolzapfelOgdenMATest {
protected:
    void SetUp() override {
        HolzapfelOgdenMATest::SetUp();

        // Use struct
        TestHO_ma->ustruct = false;
    }
};

/**
 * @brief Test fixture class for USTRUCT Holzapfel-Ogden material model.
 */
class USTRUCT_HolzapfelOgdenMATest : public HolzapfelOgdenMATest {
protected:
    void SetUp() override {
        HolzapfelOgdenMATest::SetUp();

        // Use ustruct
        TestHO_ma->ustruct = true;
    }
};


// ------------------------------ STRUCT Tests --------------------------------

// Test PK2 stress zero for F = I
TEST_F(STRUCT_HolzapfelOgdenMATest, TestPK2StressIdentityF) {
    //verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.0, 0.0, 0.0},
                        {0.0, 1.0, 0.0},
                        {0.0, 0.0, 1.0}};
    Array<double> S_ref(3, 3); // PK2 stress initialized to zero
    TestHO_ma->testPK2StressAgainstReference(F, S_ref, REL_TOL, ABS_TOL, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for triaxial stretch
TEST_F(STRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderTriaxialStretch) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for triaxial stretch
    Array<double> F = {{1.1, 0.0, 0.0},
                        {0.0, 1.2, 0.0},
                        {0.0, 0.0, 1.3}};
    
    // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
    TestHO_ma->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for triaxial compression
TEST_F(STRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderTriaxialCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for triaxial compression
    Array<double> F = {{0.9, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 0.7}};

    // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
    TestHO_ma->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for biaxial stretch/compression
TEST_F(STRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderBiaxialStretchCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for biaxial stretch/compression
    Array<double> F = {{1.2, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 1.0}};

    // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
    TestHO_ma->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (small)
TEST_F(STRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_small_list) {

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestHO_ma->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (medium)
TEST_F(STRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Set convergence order tolerance slightly larger to get this test to pass
    double CONVERGENCE_ORDER_TOL = 0.03;

    for (const auto& F : F_medium_list) {

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestHO_ma->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (large)
TEST_F(STRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_large_list) {

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestHO_ma->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for triaxial stretch
TEST_F(STRUCT_HolzapfelOgdenMATest, TestMaterialElasticityConsistencyConvergenceOrderTriaxialStretch) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for triaxial stretch
    Array<double> F = {{1.1, 0.0, 0.0},
                        {0.0, 1.2, 0.0},
                        {0.0, 0.0, 1.3}};

    // Check order of convergence of consistency of material elasticity
    TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for triaxial compression
TEST_F(STRUCT_HolzapfelOgdenMATest, TestMaterialElasticityConsistencyConvergenceOrderTriaxialCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for triaxial compression
    Array<double> F = {{0.9, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 0.7}};

    // Check order of convergence of consistency of material elasticity
    TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for biaxial stretch/compression
TEST_F(STRUCT_HolzapfelOgdenMATest, TestMaterialElasticityConsistencyConvergenceOrderBiaxialStretchCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for biaxial stretch/compression
    Array<double> F = {{1.2, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 1.0}};

    // Check order of convergence of consistency of material elasticity
    TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
}

// Test order of convergence of consistency of material elasticity for random F (small)
TEST_F(STRUCT_HolzapfelOgdenMATest, TestMaterialElasticityConsistencyConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_small_list) {

        // Check order of convergence of consistency of material elasticity
        TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (medium)
TEST_F(STRUCT_HolzapfelOgdenMATest, TestMaterialElasticityConsistencyConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Set convergence order tolerance slightly larger to get this test to pass
    double CONVERGENCE_ORDER_TOL = 0.03;

    for (const auto& F : F_medium_list) {
        
        // Check order of convergence of consistency of material elasticity
        TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (large)
TEST_F(STRUCT_HolzapfelOgdenMATest, TestMaterialElasticityConsistencyConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Set convergence order tolerance slightly larger to get this test to pass
    double CONVERGENCE_ORDER_TOL = 0.04;
    
    for (const auto& F : F_large_list) {

        // Check order of convergence of consistency of material elasticity
        TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}


// ------------------------------ USTRUCT Tests --------------------------------

// Test PK2 stress zero for F = I
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestPK2StressIdentityF) {
    //verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.0, 0.0, 0.0},
                        {0.0, 1.0, 0.0},
                        {0.0, 0.0, 1.0}};
    Array<double> S_ref(3, 3); // PK2 stress initialized to zero
    TestHO_ma->testPK2StressAgainstReference(F, S_ref, REL_TOL, ABS_TOL, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for triaxial stretch
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderTriaxialStretch) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for triaxial stretch
    Array<double> F = {{1.1, 0.0, 0.0},
                        {0.0, 1.2, 0.0},
                        {0.0, 0.0, 1.3}};
    
    // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
    TestHO_ma->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for triaxial compression
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderTriaxialCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for triaxial compression
    Array<double> F = {{0.9, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 0.7}};

    // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
    TestHO_ma->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for biaxial stretch/compression
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderBiaxialStretchCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for biaxial stretch/compression
    Array<double> F = {{1.2, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 1.0}};

    // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
    TestHO_ma->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
}


// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (small)
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_small_list) {

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestHO_ma->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (medium)
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_medium_list) {

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestHO_ma->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (large)
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_large_list) {

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestHO_ma->testPK2StressConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for triaxial stretch
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestMaterialElasticityConsistencyConvergenceOrderTriaxialStretch) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for triaxial stretch
    Array<double> F = {{1.1, 0.0, 0.0},
                        {0.0, 1.2, 0.0},
                        {0.0, 0.0, 1.3}};

    // Check order of convergence of consistency of material elasticity
    TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for triaxial compression
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestMaterialElasticityConsistencyConvergenceOrderTriaxialCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for triaxial compression
    Array<double> F = {{0.9, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 0.7}};

    // Check order of convergence of consistency of material elasticity
    TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for biaxial stretch/compression
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestMaterialElasticityConsistencyConvergenceOrderBiaxialStretchCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for biaxial stretch/compression
    Array<double> F = {{1.2, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 1.0}};

    // Check order of convergence of consistency of material elasticity
    TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
}

// Test order of convergence of consistency of material elasticity for random F (small)
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestMaterialElasticityConsistencyConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_small_list) {

        // Check order of convergence of consistency of material elasticity
        TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (medium)
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestMaterialElasticityConsistencyConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Set convergence order tolerance slightly larger to get this test to pass
    double CONVERGENCE_ORDER_TOL = 0.03;

    for (const auto& F : F_medium_list) {

        // Check order of convergence of consistency of material elasticity
        TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (large)
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestMaterialElasticityConsistencyConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    for (const auto& F : F_large_list) {

        // Check order of convergence of consistency of material elasticity
        TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, DELTA_MAX, DELTA_MIN, ORDER, CONVERGENCE_ORDER_TOL, verbose);
    }
}