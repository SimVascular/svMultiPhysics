#include "test_material_neohookean.h"

// ----------------------------------------------------------------------------
// --------------------------- Neo-Hookean Material ---------------------------
// ----------------------------------------------------------------------------

/**
 * @brief Test fixture class for the Neo-Hookean material model.
 *
 * This class sets up the necessary parameters and objects for testing the Neo-Hookean material model.
 */
class NeoHookeanTest : public MaterialTestFixture {
protected:
    // Material parameters object
    NeoHookeanParams params;

    // Add the test object
    TestNeoHookean* TestNH;

    // Setup method to initialize variables before each test
    void SetUp() override {

        MaterialTestFixture::SetUp();

        // Set random values for the Neo-Hookean parameters between 1000 and 10000
        params.C10 = getRandomDouble(1000.0, 10000.0);

        // Initialize the test object
        TestNH = new TestNeoHookean(params);
    }

    // TearDown method to clean up after each test, if needed
    void TearDown() override {
        // Clean up the test object
        delete TestNH;
        TestNH = nullptr;
    }
};

/**
 * @brief Test fixture class for STRUCT Neo-Hookean material model.
 */
class STRUCT_NeoHookeanTest : public NeoHookeanTest {
protected:
    void SetUp() override {
        NeoHookeanTest::SetUp();

        // Use struct
        TestNH->ustruct = false;
    }
};

/**
 * @brief Test fixture class for USTRUCT Neo-Hookean material model.
 */
class USTRUCT_NeoHookeanTest : public NeoHookeanTest {
protected:
    void SetUp() override {
        NeoHookeanTest::SetUp();

        // Use ustruct
        TestNH->ustruct = true;
    }
};

// ------------------------------ STRUCT Tests --------------------------------

// Test PK2 stress zero for F = I
TEST_F(STRUCT_NeoHookeanTest, TestPK2StressIdentityF) {
    //verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.0, 0.0, 0.0},
                        {0.0, 1.0, 0.0},
                        {0.0, 0.0, 1.0}};
    Array<double> S_ref(3, 3); // PK2 stress initialized to zero
    TestNH->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (small)
TEST_F(STRUCT_NeoHookeanTest, TestPK2StressConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S
    
    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestNH->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (medium)
TEST_F(STRUCT_NeoHookeanTest, TestPK2StressConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) { 
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestNH->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (large)
TEST_F(STRUCT_NeoHookeanTest, TestPK2StressConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestNH->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (small)
TEST_F(STRUCT_NeoHookeanTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestNH->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (medium)
TEST_F(STRUCT_NeoHookeanTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestNH->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (large)
TEST_F(STRUCT_NeoHookeanTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestNH->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// ------------------------------ USTRUCT Tests --------------------------------

// Test PK2 stress zero for F = I
TEST_F(USTRUCT_NeoHookeanTest, TestPK2StressIdentityF) {
    //verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.0, 0.0, 0.0},
                        {0.0, 1.0, 0.0},
                        {0.0, 0.0, 1.0}};
    Array<double> S_ref(3, 3); // PK2 stress initialized to zero
    TestNH->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (small)
TEST_F(USTRUCT_NeoHookeanTest, TestPK2StressConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestNH->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (medium)
TEST_F(USTRUCT_NeoHookeanTest, TestPK2StressConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestNH->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (large)
TEST_F(USTRUCT_NeoHookeanTest, TestPK2StressConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestNH->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (small)
TEST_F(USTRUCT_NeoHookeanTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestNH->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (medium)
TEST_F(USTRUCT_NeoHookeanTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestNH->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (large)
TEST_F(USTRUCT_NeoHookeanTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestNH->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}
