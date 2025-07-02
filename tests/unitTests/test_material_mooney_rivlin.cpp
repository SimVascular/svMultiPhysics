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
        TestMR = nullptr;
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
    TestMR->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (small)
TEST_F(STRUCT_MooneyRivlinTest, TestPK2StressConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestMR->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (medium)
TEST_F(STRUCT_MooneyRivlinTest, TestPK2StressConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestMR->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (large)
TEST_F(STRUCT_MooneyRivlinTest, TestPK2StressConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestMR->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (small)
TEST_F(STRUCT_MooneyRivlinTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestMR->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (medium)
TEST_F(STRUCT_MooneyRivlinTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestMR->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (large)
TEST_F(STRUCT_MooneyRivlinTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestMR->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
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
    TestMR->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (small)
TEST_F(USTRUCT_MooneyRivlinTest, TestPK2StressConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestMR->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (medium)
TEST_F(USTRUCT_MooneyRivlinTest, TestPK2StressConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestMR->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (large)
TEST_F(USTRUCT_MooneyRivlinTest, TestPK2StressConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestMR->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (small)
TEST_F(USTRUCT_MooneyRivlinTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestMR->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (medium)
TEST_F(USTRUCT_MooneyRivlinTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestMR->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (large)
TEST_F(USTRUCT_MooneyRivlinTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestMR->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}
