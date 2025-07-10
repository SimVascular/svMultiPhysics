#include "test_material_miehe_94_volumetric_penalty.h"

// ----------------------------------------------------------------------------
// --------------------------- Miehe 94 Volumetric Penalty ---------------------
// ----------------------------------------------------------------------------

/**
 * @brief Test fixture class for the Miehe94 Volumetric penalty model.
 * 
 * This class sets up the necessary parameters and objects for testing the Miehe94 Volumetric penalty model.
 */
class Miehe94VolumetricPenaltyTest : public MaterialTestFixture {
    protected:
        // Material parameters object
        Miehe94VolumetricPenaltyParams params;
    
        // Add the test object
        TestMiehe94VolumetricPenalty* TestM94;
    
        // Setup method to initialize variables before each test
        void SetUp() override {
    
            MaterialTestFixture::SetUp();
    
            // Set random values for the Miehe94 penalty parameters between 1000 and 10000
            params.kappa = getRandomDouble(1000.0, 10000.0);
    
            // Initialize the test object
            TestM94 = new TestMiehe94VolumetricPenalty(params);
        }
    
        // TearDown method to clean up after each test, if needed
        void TearDown() override {
            // Clean up the test object
            delete TestM94;
            TestM94 = nullptr;
        }
    };
    
    /**
     * @brief Test fixture class for STRUCT Miehe94 penalty model.
     */
    class STRUCT_Miehe94VolumetricPenaltyTest : public Miehe94VolumetricPenaltyTest {
    protected:
        void SetUp() override {
            Miehe94VolumetricPenaltyTest::SetUp();
    
            // Use struct
            //TestM94->ustruct = false;
        }
    };
    
    /**
     * @brief Test fixture class for USTRUCT Miehe94 penalty model.
     */
    class USTRUCT_Miehe94VolumetricPenaltyTest : public Miehe94VolumetricPenaltyTest {
    protected:
        void SetUp() override {
            Miehe94VolumetricPenaltyTest::SetUp();
    
            // Use ustruct
            //TestM94->ustruct = true;
        }
    };

// Test PK2 stress zero for F = I
TEST_F(STRUCT_Miehe94VolumetricPenaltyTest, TestPK2StressIdentityF) {
    //verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.0, 0.0, 0.0},
                        {0.0, 1.0, 0.0},
                        {0.0, 0.0, 1.0}};
    Array<double> S_ref(3, 3); // PK2 stress initialized to zero
    TestM94->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress zero for prescribed isochoric deformation
TEST_F(STRUCT_Miehe94VolumetricPenaltyTest, TestPK2StressPrescribedIsochoricDeformation) {
    //verbose = true; // Show values of S and S_ref

    // Check isochoric deformation produces zero PK2 stress
    Array<double> F = {{1.1, 0.0, 0.0},
                        {0.0, 1.2, 0.0},
                        {0.0, 0.0, 1.0/(1.1*1.2)}};
    Array<double> S_ref(3, 3); // PK2 stress initialized to zero
    TestM94->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (small)
TEST_F(STRUCT_Miehe94VolumetricPenaltyTest, TestPK2StressConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestM94->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (medium)
TEST_F(STRUCT_Miehe94VolumetricPenaltyTest, TestPK2StressConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestM94->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (large)
TEST_F(STRUCT_Miehe94VolumetricPenaltyTest, TestPK2StressConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestM94->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (small)
TEST_F(STRUCT_Miehe94VolumetricPenaltyTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestM94->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (medium)
TEST_F(STRUCT_Miehe94VolumetricPenaltyTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestM94->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (large)
TEST_F(STRUCT_Miehe94VolumetricPenaltyTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestM94->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// ------------------------------ USTRUCT Tests --------------------------------

// Test rho and beta values for random p
TEST_F(USTRUCT_Miehe94VolumetricPenaltyTest, TestRhoBeta) {
    //verbose = true; // Show values of rho, beta and rho_ref

    // Generate random pressure p between 100 and 1000
    double p = getRandomDouble(100.0, 1000.0);

    // Other parameters
    double rho0 = 1000.0; // Reference density
    double kappa = TestM94->params.kappa; // Volumetric penalty parameter

    // Compute reference values for rho, beta, drho/dp and dbeta/dp
    // See ustruct paper (https://doi.org/10.1016/j.cma.2018.03.045) Section 2.4
    double rho_ref = rho0 * (1.0 + p/kappa);
    double beta_ref = 1.0 / (p + kappa);
    double drhodp_ref = rho0 / kappa; // Derivative of rho with respect to p
    double dbetadp_ref = -1.0 / pow(p + kappa, 2); // Derivative of beta with respect to p

    // Check rho, beta, drho/dp and dbeta/dp against reference values
    TestM94->testRhoBetaAgainstReference(p, rho0, rho_ref, beta_ref, drhodp_ref, dbetadp_ref, rel_tol, abs_tol, verbose);
}


