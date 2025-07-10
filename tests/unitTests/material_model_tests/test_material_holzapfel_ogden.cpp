#include "test_material_holzapfel_ogden.h"





// ----------------------------------------------------------------------------
// ----------------------- Holzapfel-Ogden Material ---------------------------
// ----------------------------------------------------------------------------

/**
 * @brief Test fixture class for the Holzapfel-Ogden material model.
 * 
 * This class sets up the necessary parameters and objects for testing the Holzapfel-Ogden material model.
*/
class HolzapfelOgdenTest : public :: MaterialTestFixture {
protected:
    // Material parameters object
    HolzapfelOgdenParams params;

    // Add the test object
    TestHolzapfelOgden* TestHO;

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
        double norm_f = sqrt(params.f[0]*params.f[0] + params.f[1]*params.f[1] + params.f[2]*params.f[2]);
        params.f[0] /= norm_f; params.f[1] /= norm_f; params.f[2] /= norm_f;

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
        double norm_s = sqrt(params.s[0]*params.s[0] + params.s[1]*params.s[1] + params.s[2]*params.s[2]);
        params.s[0] /= norm_s; params.s[1] /= norm_s; params.s[2] /= norm_s;

        // Check f.s = 0
        double dot_fs = params.f[0]*params.s[0] + params.f[1]*params.s[1] + params.f[2]*params.s[2];
        if (fabs(dot_fs) > 1e-6) {
            std::cout << "f.s = " << dot_fs << std::endl;
            std::cout << "f = [" << params.f[0] << ", " << params.f[1] << ", " << params.f[2] << "]" << std::endl;
            std::cout << "s = [" << params.s[0] << ", " << params.s[1] << ", " << params.s[2] << "]" << std::endl;
            throw std::runtime_error("f and s are not orthogonal");
        }


        // Initialize the test object
        TestHO = new TestHolzapfelOgden(params);
    }

    // TearDown method to clean up after each test, if needed
    void TearDown() override {
        // Clean up the test object
        delete TestHO;
        TestHO = nullptr;
    }
};

/**
 * @brief Test fixture class for STRUCT Holzapfel-Ogden material model.
 */
class STRUCT_HolzapfelOgdenTest : public HolzapfelOgdenTest {
protected:
    void SetUp() override {
        HolzapfelOgdenTest::SetUp();

        // Use struct
        TestHO->ustruct = false;
    }
};

/**
 * @brief Test fixture class for USTRUCT Holzapfel-Ogden material model.
 */
class USTRUCT_HolzapfelOgdenTest : public HolzapfelOgdenTest {
protected:
    void SetUp() override {
        HolzapfelOgdenTest::SetUp();

        // Use ustruct
        TestHO->ustruct = true;
    }
};

// ------------------------------ STRUCT Tests --------------------------------

// Test PK2 stress zero for F = I
TEST_F(STRUCT_HolzapfelOgdenTest, TestPK2StressIdentityF) {
    //verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.0, 0.0, 0.0},
                        {0.0, 1.0, 0.0},
                        {0.0, 0.0, 1.0}};
    Array<double> S_ref(3, 3); // PK2 stress initialized to zero
    TestHO->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for triaxial stretch
TEST_F(STRUCT_HolzapfelOgdenTest, TestPK2StressConvergenceOrderTriaxialStretch) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for triaxial stretch
    Array<double> F = {{1.1, 0.0, 0.0},
                        {0.0, 1.2, 0.0},
                        {0.0, 0.0, 1.3}};

    // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
    TestHO->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for triaxial compression
TEST_F(STRUCT_HolzapfelOgdenTest, TestPK2StressConvergenceOrderTriaxialCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for triaxial compression
    Array<double> F = {{0.9, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 0.7}};

    // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
    TestHO->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for biaxial stretch/compression
TEST_F(STRUCT_HolzapfelOgdenTest, TestPK2StressConvergenceOrderBiaxialStretchCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for biaxial stretch/compression
    Array<double> F = {{1.2, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 1.0}};

    // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
    TestHO->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (small)
TEST_F(STRUCT_HolzapfelOgdenTest, TestPK2StressConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestHO->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (medium)
TEST_F(STRUCT_HolzapfelOgdenTest, TestPK2StressConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestHO->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (large)
TEST_F(STRUCT_HolzapfelOgdenTest, TestPK2StressConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestHO->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}
//triaxial extension, compression and biaxial extension - 3 new tests; struct and ustruct HO and HO-ma.
// Test order of convergence of consistency of material elasticity for triaxial stretch
TEST_F(STRUCT_HolzapfelOgdenTest, TestMaterialElasticityConsistencyConvergenceOrderTriaxialStretch) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for triaxial stretch
    Array<double> F = {{1.1, 0.0, 0.0},
                        {0.0, 1.2, 0.0},
                        {0.0, 0.0, 1.3}};

    // Check order of convergence of consistency of material elasticity
    TestHO->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for triaxial compression
TEST_F(STRUCT_HolzapfelOgdenTest, TestMaterialElasticityConsistencyConvergenceOrderTriaxialCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for triaxial compression
    Array<double> F = {{0.9, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 0.7}};

    // Check order of convergence of consistency of material elasticity
    TestHO->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for biaxial stretch/compression
TEST_F(STRUCT_HolzapfelOgdenTest, TestMaterialElasticityConsistencyConvergenceOrderBiaxialStretchCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for biaxial stretch/compression
    Array<double> F = {{1.2, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 1.0}};

    // Check order of convergence of consistency of material elasticity
    TestHO->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence of consistency of material elasticity for random F (small)
TEST_F(STRUCT_HolzapfelOgdenTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestHO->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (medium)
TEST_F(STRUCT_HolzapfelOgdenTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestHO->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (large)
TEST_F(STRUCT_HolzapfelOgdenTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestHO->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// ------------------------------ USTRUCT Tests --------------------------------

// Test PK2 stress zero for F = I
TEST_F(USTRUCT_HolzapfelOgdenTest, TestPK2StressIdentityF) {
    //verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.0, 0.0, 0.0},
                        {0.0, 1.0, 0.0},
                        {0.0, 0.0, 1.0}};
    Array<double> S_ref(3, 3); // PK2 stress initialized to zero
    TestHO->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for triaxial stretch
TEST_F(USTRUCT_HolzapfelOgdenTest, TestPK2StressConvergenceOrderTriaxialStretch) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for triaxial stretch
    Array<double> F = {{1.1, 0.0, 0.0},
                        {0.0, 1.2, 0.0},
                        {0.0, 0.0, 1.3}};

    // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
    TestHO->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for triaxial compression
TEST_F(USTRUCT_HolzapfelOgdenTest, TestPK2StressConvergenceOrderTriaxialCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for triaxial compression
    Array<double> F = {{0.9, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 0.7}};
    
    // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
    TestHO->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for biaxial stretch/compression
TEST_F(USTRUCT_HolzapfelOgdenTest, TestPK2StressConvergenceOrderBiaxialStretchCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for biaxial stretch/compression
    Array<double> F = {{1.2, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 1.0}};

    // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
    TestHO->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (small)
TEST_F(USTRUCT_HolzapfelOgdenTest, TestPK2StressConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestHO->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (medium)
TEST_F(USTRUCT_HolzapfelOgdenTest, TestPK2StressConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestHO->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (large)
TEST_F(USTRUCT_HolzapfelOgdenTest, TestPK2StressConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestHO->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for triaxial stretch
TEST_F(USTRUCT_HolzapfelOgdenTest, TestMaterialElasticityConsistencyConvergenceOrderTriaxialStretch) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for triaxial stretch
    Array<double> F = {{1.1, 0.0, 0.0},
                        {0.0, 1.2, 0.0},
                        {0.0, 0.0, 1.3}};

    // Check order of convergence of consistency of material elasticity
    TestHO->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for triaxial compression
TEST_F(USTRUCT_HolzapfelOgdenTest, TestMaterialElasticityConsistencyConvergenceOrderTriaxialCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for triaxial compression
    Array<double> F = {{0.9, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 0.7}};

    // Check order of convergence of consistency of material elasticity
    TestHO->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for biaxial stretch/compression
TEST_F(USTRUCT_HolzapfelOgdenTest, TestMaterialElasticityConsistencyConvergenceOrderBiaxialStretchCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for biaxial stretch/compression
    Array<double> F = {{1.2, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 1.0}};

    // Check order of convergence of consistency of material elasticity
    TestHO->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence of consistency of material elasticity for random F (small)
TEST_F(USTRUCT_HolzapfelOgdenTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestHO->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (medium)
TEST_F(USTRUCT_HolzapfelOgdenTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestHO->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (large)
TEST_F(USTRUCT_HolzapfelOgdenTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestHO->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}