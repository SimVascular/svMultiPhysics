#include "test_material_CANN_neohookean.h"
#include "test_material_neohookean.h"


// ----------------------------------------------------------------------------
// -------- Compare CANN w/ NH params against NH Model ------------------------
// ----------------------------------------------------------------------------
/**
 * @brief Test fixture class for comparing the CANN model with Neo-Hookean mdoel.
 *
 * This class sets up the necessary parameters and objects for comparing CANN model with NH parameters against the Neo-Hookean model.
 */
class NeoHookeanCompareTest : public MaterialTestFixture {
protected:
    // Material parameters objects
    NeoHookeanParams params_NH; // Neo-Hookean parameters
    CANN_NH_Params params_CANN_NH; // CANN with neo-hookean parameters

    // Add the test objects
    TestNeoHookean* TestNH;
    TestCANN_NH* TestCANNNH;

    // Setup method to initialize variables before each test
    void SetUp() override {
        MaterialTestFixture::SetUp();

        // Set values for the Neo-Hookean parameters
        double C10 = 4.0094326666666664e+07;

        params_NH.C10 = C10;

        params_CANN_NH.Table.resize(1);  // Ensure it has at least one entry
        params_CANN_NH.Table[0].invariant_index.value_ = 1;
        params_CANN_NH.Table[0].activation_functions.value_ = {1,1,1};
        params_CANN_NH.Table[0].weights.value_ = {1.0,1.0,C10};

        // Initialize the test objects
        TestNH = new TestNeoHookean(params_NH);
        TestCANNNH = new TestCANN_NH(params_CANN_NH);

        if (TestCANNNH == nullptr) {
            throw std::runtime_error("TestCANNNH is not properly initialized");
        }
        if (TestNH == nullptr) {
            throw std::runtime_error("TestNH is not properly initialized");
        }
    }

    // TearDown method to clean up after each test, if needed
    void TearDown() override {
        // Clean up the test objects
        delete TestNH;
        TestNH = nullptr;
        delete TestCANNNH;
        TestCANNNH = nullptr;
    }
};

/**
 * @brief Test fixture class for STRUCT Neo-Hookean Compare model.
 */
class STRUCT_CANNNeoHookeanTest : public NeoHookeanCompareTest {
protected:
    void SetUp() override {
        NeoHookeanCompareTest::SetUp();

        // Use struct
        TestNH->ustruct = false;
        TestCANNNH->ustruct = false;
    }
};

/**
 * @brief Test fixture class for USTRUCT Neo-Hookean Compare model.
 */
class USTRUCT_CANNNeoHookeanTest : public NeoHookeanCompareTest {
protected:
    void SetUp() override {
        NeoHookeanCompareTest::SetUp();

        // Use ustruct
        TestNH->ustruct = true;
        TestCANNNH->ustruct = true;
    }
};


// ----------------------------------------------------------------------------
// ------------- CANN w/ NH param against Neo-Hookean Material ----------------
// ----------------------------------------------------------------------------

// ------------------------------ STRUCT Tests --------------------------------

// Test PK2 stress zero for F = I
TEST_F(STRUCT_CANNNeoHookeanTest, TestPK2StressIdentityF) {
    verbose = false; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.0, 0.0, 0.0},
                       {0.0, 1.0, 0.0},
                       {0.0, 0.0, 1.0}};
    Array<double> S_ref(3,3); // PK2 stress initialized to zero - want to get result from NH and set that to S_ref
    Array<double> Dm(6,6);
    TestNH->compute_pk2cc(F,S_ref,Dm); // Computing S_ref from NH
    TestCANNNH->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose); // Comparing with CANN
}

// Test PK2 stress
TEST_F(STRUCT_CANNNeoHookeanTest, TestPK2StressTriaxialStretch) {
    verbose = false; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.1, 0.0, 0.0},
                       {0.0, 1.2, 0.0},
                       {0.0, 0.0, 0.757}};
    Array<double> S_ref(3,3); // PK2 stress initialized to zero - want to get result from NH and set that to S_ref
    Array<double> Dm(6,6);
    TestNH->compute_pk2cc(F,S_ref,Dm); // Computing S_ref from NH
    TestCANNNH->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose); // Comparing with CANN
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (small)
TEST_F(STRUCT_CANNNeoHookeanTest, TestPK2StressConvergenceOrderAgainstReferenceRandomFSmall) {
    verbose = false; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to C array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        Array<double> S_ref(3,3), Dm(6,6);
        TestCANNNH->compute_pk2cc(F,S_ref,Dm); // Computing S_ref from CANN
        TestNH->testPK2StressConvergenceOrderAgainstReference(F, S_ref, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (medium)
TEST_F(STRUCT_CANNNeoHookeanTest, TestPK2StressConvergenceOrderAgainstReferenceRandomFMedium) {
    verbose = false; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to C array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        Array<double> S_ref(3,3), Dm(6,6);
        TestCANNNH->compute_pk2cc(F,S_ref,Dm); // Computing S_ref from CANN
        TestNH->testPK2StressConvergenceOrderAgainstReference(F, S_ref, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (large)
TEST_F(STRUCT_CANNNeoHookeanTest, TestPK2StressConvergenceOrderAgainstReferenceRandomFLarge) {
    verbose = false; // Show order of convergence, errors, F, S

    // Use struct
    //TestQVP->ustruct = false;
    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to C array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        Array<double> S_ref(3,3), Dm(6,6);
        TestCANNNH->compute_pk2cc(F,S_ref,Dm); // Computing S_ref from CANN
        TestNH->testPK2StressConvergenceOrderAgainstReference(F,S_ref, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (small)
TEST_F(STRUCT_CANNNeoHookeanTest, TestMaterialElasticityConsistencyConvergenceOrderAgainstNHRandomFSmall) {
    verbose = false; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to C array
        convertToArray(F_std, F);
        // Generating perturbation
        Array<double> dF(3,3);
        std::vector<double> deltas;
        TestNH->generatePerturbationdF(F,dF,delta_max, delta_min, deltas,order,verbose);
        // Calculating dS and CCdE
        Array<double> dS(3,3), CCdE(3,3);
        for (int i = 0; i < deltas.size(); i++) {
            TestCANNNH->calcCCdEFiniteDifference(F, dF, deltas[i], order, CCdE);
            TestNH->calcdSFiniteDifference(F, dF, deltas[i], order, dS);

            // Check order of convergence of consistency of material elasticity
            TestNH->testMaterialElasticityConsistencyConvergenceOrderBetweenMaterialModels(F, dS, CCdE, deltas, order, convergence_order_tol, verbose);
        }
    }
}

// Test order of convergence of consistency of material elasticity for random F (medium)
TEST_F(STRUCT_CANNNeoHookeanTest, TestMaterialElasticityConsistencyConvergenceOrderAgainstNHRandomFMedium) {
    verbose = false; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to C array
        convertToArray(F_std, F);
        // Generating perturbation
        Array<double> dF(3,3);
        std::vector<double> deltas;
        TestNH->generatePerturbationdF(F,dF,delta_max, delta_min, deltas,order,verbose);
        // Calculating dS and CCdE
        Array<double> dS(3,3), CCdE(3,3);
        for (int i = 0; i < deltas.size(); i++) {
            TestCANNNH->calcCCdEFiniteDifference(F, dF, deltas[i], order, CCdE);
            TestNH->calcdSFiniteDifference(F, dF, deltas[i], order, dS);

            // Check order of convergence of consistency of material elasticity
            TestNH->testMaterialElasticityConsistencyConvergenceOrderBetweenMaterialModels(F, dS, CCdE, deltas, order, convergence_order_tol, verbose);
        }
    }
}

// Test order of convergence of consistency of material elasticity for random F (large)
TEST_F(STRUCT_CANNNeoHookeanTest, TestMaterialElasticityConsistencyConvergenceOrderAgainstNHRandomFLarge) {
    verbose = false; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to array
        convertToArray(F_std, F);
        // Generating perturbation
        Array<double> dF(3,3);
        std::vector<double> deltas;
        TestNH->generatePerturbationdF(F,dF,delta_max, delta_min, deltas,order,verbose);
        // Calculating dS and CCdE
        Array<double> dS(3,3), CCdE(3,3);
        for (int i = 0; i < deltas.size(); i++) {
            TestCANNNH->calcCCdEFiniteDifference(F, dF, deltas[i], order, CCdE);
            TestNH->calcdSFiniteDifference(F, dF, deltas[i], order, dS);

            // Check order of convergence of consistency of material elasticity
            TestNH->testMaterialElasticityConsistencyConvergenceOrderBetweenMaterialModels(F, dS, CCdE, deltas, order, convergence_order_tol, verbose);
        }
    }
}

// Test PK2 stress
TEST_F(STRUCT_CANNNeoHookeanTest, TestMaterialElasticityAgainstReferenceIdentityF) {
    verbose = false; // Show values of S and S_ref

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
    TestNH->calcMaterialElasticityReference(F,CC_ref,verbose);
    TestCANNNH->testMaterialElasticityAgainstReference(F, CC_ref, rel_tol, abs_tol, verbose); // Comparing with CANN
}

// Test PK2 stress
TEST_F(STRUCT_CANNNeoHookeanTest, TestMaterialElasticityAgainstReference) {
    verbose = false; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.1, 0.0, 0.0},
                       {0.0, 1.2, 0.0},
                       {0.0, 0.0, 1.3}};
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
    TestNH->calcMaterialElasticityReference(F,CC_ref,verbose);
    TestCANNNH->testMaterialElasticityAgainstReference(F, CC_ref, rel_tol, abs_tol, verbose); // Comparing with CANN
}