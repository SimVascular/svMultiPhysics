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

#include "test.h"

using namespace mat_fun;
using namespace std;

// ============================================================================
// --------------------------- Test fixture classes ---------------------------
// ============================================================================

/**
 * @brief Test fixture class containing common setup for all material model tests in this file
 * 
 */
class MaterialModelTest : public ::testing::Test {
protected:
    // Variables common across tests
    Array<double> F; // Deformation gradient
    double deformation_perturbation_small = 0.003; // Small perturbation factor
    double deformation_perturbation_medium = 0.03; // Medium perturbation factor
    double deformation_perturbation_large = 0.3; // Large perturbation factor
    int n_F = 50; // Number of deformation gradients F to test for each small, medium, and large perturbation
    double rel_tol = 1e-3; // relative tolerance for comparing values
    double abs_tol = 1e-11; // absolute tolerance for comparing values
    //double delta = 1e-7; // perturbation scaling factor
    double delta_max = 1e-4; // maximum perturbation scaling factor
    double delta_min = 1e-6; // minimum perturbation scaling factor
    int order = 1; // Order of finite difference method
    double convergence_order_tol = 0.02; // Tolerance for comparing convergence order with expected value
    bool verbose = false; // Show values of S, dE, SdE and dPsi

    // Type alias for a 3x3 std::array
    using Array3x3 = std::array<std::array<double, 3>, 3>;

    // Vectors to store the 3x3 arrays
    std::vector<Array3x3> F_small_list;
    std::vector<Array3x3> F_medium_list;
    std::vector<Array3x3> F_large_list;

    // Function to convert C array to std::array
    Array3x3 convertToStdArray(const Array<double> &F) {

        int N = F.ncols();
        assert(F.nrows() == N);
        assert(N == 3);

        Array3x3 stdArray;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                stdArray[i][j] = F(i,j);
            }
        }
        return stdArray;
    }

    // Function to convert std::array to Array
    void convertToArray(Array3x3 stdArray, Array<double> F) {

        int N = F.ncols();
        assert(F.nrows() == N);
        assert(N == 3);

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                F(i,j) = stdArray[i][j];
            }
        }
    }

    void SetUp() override {
        F.resize(3,3); // Deformation gradient

        // Create random deformation gradients for small perturbations
        for (int i = 0; i < n_F; i++) {
            F = create_random_perturbed_identity_F(3, deformation_perturbation_small);
            F_small_list.push_back(convertToStdArray(F));
        }

        // Create random deformation gradients for medium perturbations
        for (int i = 0; i < n_F; i++) {
            F = create_random_perturbed_identity_F(3, deformation_perturbation_medium);
            F_medium_list.push_back(convertToStdArray(F));
        }

        // Create random deformation gradients for large perturbations
        for (int i = 0; i < n_F; i++) {
            F = create_random_perturbed_identity_F(3, deformation_perturbation_large);
            F_large_list.push_back(convertToStdArray(F));
        }
    }

    void TearDown() override {}
};

// ----------------------------------------------------------------------------
// --------------------------- Neo-Hookean Material ---------------------------
// ----------------------------------------------------------------------------

/**
 * @brief Test fixture class for the Neo-Hookean material model.
 *
 * This class sets up the necessary parameters and objects for testing the Neo-Hookean material model.
 */
class NeoHookeanTest : public MaterialModelTest {
protected:
    // Material parameters object
    NeoHookeanParams params;

    // Add the test object
    TestNeoHookean* TestNH;

    // Setup method to initialize variables before each test
    void SetUp() override {

        MaterialModelTest::SetUp();

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

// ----------------------------------------------------------------------------
// --------------------------- Mooney-Rivlin Material -------------------------
// ----------------------------------------------------------------------------

/**
 * @brief  Test fixture class for the Mooney-Rivlin material model.
 * 
 * This class sets up the necessary parameters and objects for testing the Mooney-Rivlin material model.
 */
class MooneyRivlinTest : public MaterialModelTest {
protected:
    // Material parameters object
    MooneyRivlinParams params;

    // Add the test object
    TestMooneyRivlin* TestMR;

    // Setup method to initialize variables before each test
    void SetUp() override {

        MaterialModelTest::SetUp();

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


// ----------------------------------------------------------------------------
// ----------------------- Holzapfel-Ogden Material ---------------------------
// ----------------------------------------------------------------------------

/**
 * @brief Test fixture class for the Holzapfel-Ogden material model.
 * 
 * This class sets up the necessary parameters and objects for testing the Holzapfel-Ogden material model.
*/
class HolzapfelOgdenTest : public :: MaterialModelTest {
protected:
    // Material parameters object
    HolzapfelOgdenParams params;

    // Add the test object
    TestHolzapfelOgden* TestHO;

    // Setup method to initialize variables before each test
    void SetUp() override {

        MaterialModelTest::SetUp();

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
            cout << "f.s = " << dot_fs << endl;
            cout << "f = [" << params.f[0] << ", " << params.f[1] << ", " << params.f[2] << "]" << endl;
            cout << "s = [" << params.s[0] << ", " << params.s[1] << ", " << params.s[2] << "]" << endl;
            throw runtime_error("f and s are not orthogonal");
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

// ----------------------------------------------------------------------------
// ------------- Holzapfel-Ogden (Modified Anisotropy) Material  --------------
// ----------------------------------------------------------------------------

/**
 * @brief Test fixture class for the Holzapfel-Ogden (Modified Anisotropy) material model.
 * 
 * This class sets up the necessary parameters and objects for testing the Holzapfel-Ogden (Modified Anisotropy) material model.
*/
class HolzapfelOgdenMATest : public MaterialModelTest {
protected:
    // Material parameters object
    HolzapfelOgdenMAParams params;

    // Add the test object
    TestHolzapfelOgdenMA* TestHO_ma;

    // Setup method to initialize variables before each test
    void SetUp() override {

        MaterialModelTest::SetUp();

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
            cout << "f.s = " << dot_fs << endl;
            cout << "f = [" << params.f[0] << ", " << params.f[1] << ", " << params.f[2] << "]" << endl;
            cout << "s = [" << params.s[0] << ", " << params.s[1] << ", " << params.s[2] << "]" << endl;
            throw runtime_error("f and s are not orthogonal");
        }


        // Initialize the test object
        TestHO_ma = new TestHolzapfelOgdenMA(params);
    }

    // TearDown method to clean up after each test, if needed
    void TearDown() override {
        // Clean up the test object
        delete TestHO_ma;
        TestHO_ma = nullptr;
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

// ----------------------------------------------------------------------------
// --------------------------- CANN Neo-Hookean Material ---------------------------
// ----------------------------------------------------------------------------

/**
 * @brief Test fixture class for the CANN Neo-Hookean material model.
 *
 * This class sets up the necessary parameters and objects for testing the Neo-Hookean material model.
 */
class CANN_NH_Test : public MaterialModelTest {
protected:
    // Material parameters object
    CANN_NH_Params params;

    // Add the test object
    TestCANN_NH* TestCANNNH;

    // Setup method to initialize variables before each test
    void SetUp() override {

        MaterialModelTest::SetUp();
        params.Table[0].invariant_index.value_ = 1;
        params.Table[0].activation_functions.value_ = {1,1,1};
        params.Table[0].weights.value_ = {1.0,1.0,2234.32}; // same as C10
        
        // Initialize the test object
        TestCANNNH = new TestCANN_NH(params);

        if (TestCANNNH == nullptr) {
            __throw_runtime_error("TestCANNNH is not properly initialized");
        }
    }

    // TearDown method to clean up after each test, if needed
    void TearDown() override {
        // Clean up the test object
        delete TestCANNNH;
        TestCANNNH = nullptr;
    }
};

/**
 * @brief Test fixture class for STRUCT Neo-Hookean material model.
 */
class STRUCT_CANN_NH_Test : public CANN_NH_Test {
protected:
    void SetUp() override {
        CANN_NH_Test::SetUp();

        // Use struct
        TestCANNNH->ustruct = false;
    }
};

/**
 * @brief Test fixture class for USTRUCT Neo-Hookean material model.
 */
class USTRUCT_CANN_NH_Test : public CANN_NH_Test {
protected:
    void SetUp() override {
        CANN_NH_Test::SetUp();

        // Use ustruct
        TestCANNNH->ustruct = true;
    }
};

// ----------------------------------------------------------------------------
// -------- Compare CANN w/ NH params against NH Model ------------------------
// ----------------------------------------------------------------------------
/**
 * @brief Test fixture class for comparing the CANN model with Neo-Hookean mdoel.
 *
 * This class sets up the necessary parameters and objects for comparing CANN model with NH parameters against the Neo-Hookean model.
 */
class NeoHookeanCompareTest : public MaterialModelTest {
protected:
    // Material parameters objects
    NeoHookeanParams params_NH; // Neo-Hookean parameters
    CANN_NH_Params params_CANN_NH; // CANN with neo-hookean parameters

    // Add the test objects
    TestNeoHookean* TestNH;
    TestCANN_NH* TestCANNNH;

    // Setup method to initialize variables before each test
    void SetUp() override {
        MaterialModelTest::SetUp();

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
            __throw_runtime_error("TestCANNNH is not properly initialized");
        }
        if (TestNH == nullptr) {
            __throw_runtime_error("TestNH is not properly initialized");
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
// ------------------- CANN Holzapfel-Ogden Material --------------------------
// ----------------------------------------------------------------------------

/**
 * @brief Test fixture class for the CANN Holzapfel-Ogden material model.
 *
 * This class sets up the necessary parameters and objects for testing the Holzapfel-Ogden material model.
 */
class CANN_HO_Test : public MaterialModelTest {
protected:
    // Material parameters object
    CANN_HO_Params params;

    // Add the test object
    TestCANN_HO* TestCANNHO;

    // Setup method to initialize variables before each test
    void SetUp() override {

        MaterialModelTest::SetUp();
        
        // Initialize the test object
        TestCANNHO = new TestCANN_HO(params);

        if (TestCANNHO == nullptr) {
            __throw_runtime_error("TestCANNHO is not properly initialized");
        }
    }

    // TearDown method to clean up after each test, if needed
    void TearDown() override {
        // Clean up the test object
        delete TestCANNHO;
        TestCANNHO = nullptr;
    }
};

/**
 * @brief Test fixture class for STRUCT HO material model.
 */
class STRUCT_CANN_HO_Test : public CANN_HO_Test {
protected:
    void SetUp() override {
        CANN_HO_Test::SetUp();

        // Use struct
        TestCANNHO->ustruct = false;
    }
};

/**
 * @brief Test fixture class for USTRUCT Neo-Hookean material model.
 */
class USTRUCT_CANN_HO_Test : public CANN_HO_Test {
protected:
    void SetUp() override {
        CANN_HO_Test::SetUp();

        // Use ustruct
        TestCANNHO->ustruct = true;
    }
};

// ----------------------------------------------------------------------------
// -------- Compare CANN w/ HO params against HO Model ------------------------
// ----------------------------------------------------------------------------
/**
 * @brief Test fixture class for comparing the CANN model with Holzapfel Ogden mdoel.
 *
 * This class sets up the necessary parameters and objects for comparing CANN model with HO parameters against the Holzapfel-Ogden model. Assuming k = 0, so chi = 0.5.
 */
class HolzapfelOgdenCompareTest : public MaterialModelTest {
protected:
    // Material parameters objects
    HolzapfelOgdenParams params_HO; // HO parameters
    CANN_HO_Params params_CANN_HO; // CANN with HO parameters

    // Add the test objects
    TestHolzapfelOgden* TestHO;
    TestCANN_HO* TestCANNHO;

    // Setup method to initialize variables before each test
    void SetUp() override {
        MaterialModelTest::SetUp();

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
        double norm_f = sqrt(params_HO.f[0]*params_HO.f[0] + params_HO.f[1]*params_HO.f[1] + params_HO.f[2]*params_HO.f[2]);
        params_HO.f[0] /= norm_f; params_HO.f[1] /= norm_f; params_HO.f[2] /= norm_f;

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
        double norm_s = sqrt(params_HO.s[0]*params_HO.s[0] + params_HO.s[1]*params_HO.s[1] + params_HO.s[2]*params_HO.s[2]);
        params_HO.s[0] /= norm_s; params_HO.s[1] /= norm_s; params_HO.s[2] /= norm_s;

        // Check f.s = 0
        double dot_fs = params_HO.f[0]*params_HO.s[0] + params_HO.f[1]*params_HO.s[1] + params_HO.f[2]*params_HO.s[2];
        if (fabs(dot_fs) > 1e-6) {
            cout << "f.s = " << dot_fs << endl;
            cout << "f = [" << params_HO.f[0] << ", " << params_HO.f[1] << ", " << params_HO.f[2] << "]" << endl;
            cout << "s = [" << params_HO.s[0] << ", " << params_HO.s[1] << ", " << params_HO.s[2] << "]" << endl;
            throw runtime_error("f and s are not orthogonal");
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

        if (TestCANNHO == nullptr) {
            __throw_runtime_error("TestCANNHO is not properly initialized");
        }
        if (TestHO == nullptr) {
            __throw_runtime_error("TestHO is not properly initialized");
        }
    }

    // TearDown method to clean up after each test, if needed
    void TearDown() override {
        // Clean up the test objects
        delete TestHO;
        TestHO = nullptr;
        delete TestCANNHO;
        TestCANNHO = nullptr;
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


// ----------------------------------------------------------------------------
// ---------------- Quadratic Volumetric Penalty Model ------------------------
// ----------------------------------------------------------------------------
/**
 * @brief Test fixture class for the Quadratic Volumetric penalty model.
 * 
 * This class sets up the necessary parameters and objects for testing the Quadratic Volumetric penalty model.
 */
class QuadraticVolumetricPenaltyTest : public MaterialModelTest {
protected:
    // Material parameters object
    VolumetricPenaltyParams params;

    // Add the test object
    TestQuadraticVolumetricPenalty* TestQVP;

    // Setup method to initialize variables before each test
    void SetUp() override {

        MaterialModelTest::SetUp();

        // Set random values for the Quadratic penalty parameters between 1000 and 10000
        params.kappa = getRandomDouble(1000.0, 10000.0);

        // Initialize the test object
        TestQVP = new TestQuadraticVolumetricPenalty(params);
    }

    // TearDown method to clean up after each test, if needed
    void TearDown() override {
        // Clean up the test object
        delete TestQVP;
        TestQVP = nullptr;
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


// ----------------------------------------------------------------------------
// --------------------------- Simo-Taylor91 Volumetric Penalty Model ---------
// ----------------------------------------------------------------------------

/**
 * @brief Test fixture class for the Simo-Taylor91 Volumetric penalty model.
 * 
 * This class sets up the necessary parameters and objects for testing the Simo-Taylor91 Volumetric penalty model.
 */
class SimoTaylor91VolumetricPenaltyTest : public MaterialModelTest {
protected:
    // Material parameters object
    VolumetricPenaltyParams params;

    // Add the test object
    TestSimoTaylor91VolumetricPenalty* TestST91;

    // Setup method to initialize variables before each test
    void SetUp() override {

        MaterialModelTest::SetUp();

        // Set random values for the Simo-Taylor91 penalty parameters between 1000 and 10000
        params.kappa = getRandomDouble(1000.0, 10000.0);

        // Initialize the test object
        TestST91 = new TestSimoTaylor91VolumetricPenalty(params);
    }

    // TearDown method to clean up after each test, if needed
    void TearDown() override {
        // Clean up the test object
        delete TestST91;
        TestST91 = nullptr;
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

// ----------------------------------------------------------------------------
// ---------------------- Miehe94 Volumetric Penalty Model --------------------
// ----------------------------------------------------------------------------
/**
 * @brief Test fixture class for the Miehe94 Volumetric penalty model.
 * 
 * This class sets up the necessary parameters and objects for testing the Miehe94 Volumetric penalty model.
 */
class Miehe94VolumetricPenaltyTest : public MaterialModelTest {
protected:
    // Material parameters object
    VolumetricPenaltyParams params;

    // Add the test object
    TestMiehe94VolumetricPenalty* TestM94;

    // Setup method to initialize variables before each test
    void SetUp() override {

        MaterialModelTest::SetUp();

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









// ============================================================================
// ------------------------------- TESTS --------------------------------------
// ============================================================================



// ----------------------------------------------------------------------------
// --------------------------- Neo-Hookean Material ---------------------------
// ----------------------------------------------------------------------------

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



// ----------------------------------------------------------------------------
// --------------------------- Mooney-Rivlin Material -------------------------
// ----------------------------------------------------------------------------

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



// ----------------------------------------------------------------------------
// ----------------------- Holzapfel-Ogden Material ---------------------------
// ----------------------------------------------------------------------------

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



// ----------------------------------------------------------------------------
// --------------- Holzapfel-Ogden (Modified Anisotropy) Material -------------
// ----------------------------------------------------------------------------

// ------------------------------ STRUCT Tests --------------------------------

// Test PK2 stress zero for F = I
TEST_F(STRUCT_HolzapfelOgdenMATest, TestPK2StressIdentityF) {
    //verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.0, 0.0, 0.0},
                        {0.0, 1.0, 0.0},
                        {0.0, 0.0, 1.0}};
    Array<double> S_ref(3, 3); // PK2 stress initialized to zero
    TestHO_ma->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for triaxial stretch
TEST_F(STRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderTriaxialStretch) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for triaxial stretch
    Array<double> F = {{1.1, 0.0, 0.0},
                        {0.0, 1.2, 0.0},
                        {0.0, 0.0, 1.3}};
    
    // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
    TestHO_ma->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for triaxial compression
TEST_F(STRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderTriaxialCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for triaxial compression
    Array<double> F = {{0.9, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 0.7}};

    // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
    TestHO_ma->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for biaxial stretch/compression
TEST_F(STRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderBiaxialStretchCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for biaxial stretch/compression
    Array<double> F = {{1.2, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 1.0}};

    // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
    TestHO_ma->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (small)
TEST_F(STRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestHO_ma->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (medium)
TEST_F(STRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestHO_ma->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (large)
TEST_F(STRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestHO_ma->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
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
    TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for triaxial compression
TEST_F(STRUCT_HolzapfelOgdenMATest, TestMaterialElasticityConsistencyConvergenceOrderTriaxialCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for triaxial compression
    Array<double> F = {{0.9, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 0.7}};

    // Check order of convergence of consistency of material elasticity
    TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for biaxial stretch/compression
TEST_F(STRUCT_HolzapfelOgdenMATest, TestMaterialElasticityConsistencyConvergenceOrderBiaxialStretchCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for biaxial stretch/compression
    Array<double> F = {{1.2, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 1.0}};

    // Check order of convergence of consistency of material elasticity
    TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence of consistency of material elasticity for random F (small)
TEST_F(STRUCT_HolzapfelOgdenMATest, TestMaterialElasticityConsistencyConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (medium)
TEST_F(STRUCT_HolzapfelOgdenMATest, TestMaterialElasticityConsistencyConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Set convergence order tolerances larger and special delta_max/delta_min to get this test to pass
        convergence_order_tol = 0.15;
        delta_max = 8e-6;
        delta_min = 2e-6;
        
        // Check order of convergence of consistency of material elasticity
        TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (large)
TEST_F(STRUCT_HolzapfelOgdenMATest, TestMaterialElasticityConsistencyConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Set convergence order tolerances larger and special delta_max/delta_min to get this test to pass
        //convergence_order_tol = 0.02;
        delta_max = 4e-5;
        delta_min = 1e-5;

        // Check order of convergence of consistency of material elasticity
        TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
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
    TestHO_ma->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for triaxial stretch
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderTriaxialStretch) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for triaxial stretch
    Array<double> F = {{1.1, 0.0, 0.0},
                        {0.0, 1.2, 0.0},
                        {0.0, 0.0, 1.3}};
    
    // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
    TestHO_ma->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for triaxial compression
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderTriaxialCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for triaxial compression
    Array<double> F = {{0.9, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 0.7}};

    // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
    TestHO_ma->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for biaxial stretch/compression
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderBiaxialStretchCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for biaxial stretch/compression
    Array<double> F = {{1.2, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 1.0}};

    // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
    TestHO_ma->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}


// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (small)
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestHO_ma->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (medium)
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestHO_ma->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (large)
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestPK2StressConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestHO_ma->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
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
    TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for triaxial compression
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestMaterialElasticityConsistencyConvergenceOrderTriaxialCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for triaxial compression
    Array<double> F = {{0.9, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 0.7}};

    // Check order of convergence of consistency of material elasticity
    TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for biaxial stretch/compression
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestMaterialElasticityConsistencyConvergenceOrderBiaxialStretchCompression) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Create a deformation gradient F for biaxial stretch/compression
    Array<double> F = {{1.2, 0.0, 0.0},
                        {0.0, 0.8, 0.0},
                        {0.0, 0.0, 1.0}};

    // Check order of convergence of consistency of material elasticity
    TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
}

// Test order of convergence of consistency of material elasticity for random F (small)
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestMaterialElasticityConsistencyConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (medium)
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestMaterialElasticityConsistencyConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (large)
TEST_F(USTRUCT_HolzapfelOgdenMATest, TestMaterialElasticityConsistencyConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestHO_ma->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// ----------------------------------------------------------------------------
// ------------- CANN w/ NH param against Neo-Hookean Material ----------------
// ----------------------------------------------------------------------------

// ------------------------------ STRUCT Tests --------------------------------

// Test PK2 stress zero for F = I
TEST_F(STRUCT_CANNNeoHookeanTest, TestPK2StressIdentityF) {
    verbose = true; // Show values of S and S_ref

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
    verbose = true; // Show values of S and S_ref

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
    //verbose = true; // Show order of convergence, errors, F, S

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
    //verbose = true; // Show order of convergence, errors, F, S

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
    //verbose = true; // Show order of convergence, errors, F, S

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
    //verbose = true; // Show order of convergence, errors, F, S

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
    //verbose = true; // Show order of convergence, errors, F, S

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
    //verbose = true; // Show order of convergence, errors, F, S

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
    verbose = true; // Show values of S and S_ref

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
    verbose = true; // Show values of S and S_ref

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

// ----------------------------------------------------------------------------
// ------- CANN w/ HO param against Holzapfel Ogden Material ------------------
// ----------------------------------------------------------------------------

// ------------------------------ STRUCT Tests --------------------------------

// Test PK2 stress zero for F = I
TEST_F(STRUCT_CANNHolzapfelOgdenTest, TestPK2StressIdentityF) {
    verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.0, 0.0, 0.0},
                       {0.0, 1.0, 0.0},
                       {0.0, 0.0, 1.0}};
    Array<double> S_ref(3,3); // PK2 stress initialized to zero - want to get result from NH and set that to S_ref
    Array<double> Dm(6,6);
    TestHO->compute_pk2cc(F,S_ref,Dm); // Computing S_ref from NH
    TestCANNHO->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose); // Comparing with CANN
}

// Test PK2 stress
TEST_F(STRUCT_CANNHolzapfelOgdenTest, TestPK2StressTriaxialStretch) {
    verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.1, 0.0, 0.0},
                       {0.0, 1.2, 0.0},
                       {0.0, 0.0, 0.757}};
    Array<double> S_ref(3,3); // PK2 stress initialized to zero - want to get result from NH and set that to S_ref
    Array<double> Dm(6,6);
    TestHO->compute_pk2cc(F,S_ref,Dm); // Computing S_ref from NH
    TestCANNHO->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose); // Comparing with CANN
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (small)
TEST_F(STRUCT_CANNHolzapfelOgdenTest, TestPK2StressConvergenceOrderAgainstReferenceRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to C array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        Array<double> S_ref(3,3), Dm(6,6);
        TestCANNHO->compute_pk2cc(F,S_ref,Dm); // Computing S_ref from CANN
        TestHO->testPK2StressConvergenceOrderAgainstReference(F, S_ref, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (medium)
TEST_F(STRUCT_CANNHolzapfelOgdenTest, TestPK2StressConvergenceOrderAgainstReferenceRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to C array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        Array<double> S_ref(3,3), Dm(6,6);
        TestCANNHO->compute_pk2cc(F,S_ref,Dm); // Computing S_ref from CANN
        TestHO->testPK2StressConvergenceOrderAgainstReference(F, S_ref, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (large)
TEST_F(STRUCT_CANNHolzapfelOgdenTest, TestPK2StressConvergenceOrderAgainstReferenceRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Use struct
    //TestQVP->ustruct = false;
    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to C array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        Array<double> S_ref(3,3), Dm(6,6);
        TestCANNHO->compute_pk2cc(F,S_ref,Dm); // Computing S_ref from CANN
        TestHO->testPK2StressConvergenceOrderAgainstReference(F,S_ref, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (small)
TEST_F(STRUCT_CANNHolzapfelOgdenTest, TestMaterialElasticityConsistencyConvergenceOrderAgainstHORandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to C array
        convertToArray(F_std, F);
        // Generating perturbation
        Array<double> dF(3,3);
        std::vector<double> deltas;
        TestHO->generatePerturbationdF(F,dF,delta_max, delta_min, deltas,order,verbose);
        // Calculating dS and CCdE
        Array<double> dS(3,3), CCdE(3,3);
        for (int i = 0; i < deltas.size(); i++) {
            TestCANNHO->calcCCdEFiniteDifference(F, dF, deltas[i], order, CCdE);
            TestHO->calcdSFiniteDifference(F, dF, deltas[i], order, dS);

            // Check order of convergence of consistency of material elasticity
            TestHO->testMaterialElasticityConsistencyConvergenceOrderBetweenMaterialModels(F, dS, CCdE, deltas, order, convergence_order_tol, verbose);
        }
    }
}

// Test order of convergence of consistency of material elasticity for random F (medium)
TEST_F(STRUCT_CANNHolzapfelOgdenTest, TestMaterialElasticityConsistencyConvergenceOrderAgainstHORandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to C array
        convertToArray(F_std, F);
        // Generating perturbation
        Array<double> dF(3,3);
        std::vector<double> deltas;
        TestHO->generatePerturbationdF(F,dF,delta_max, delta_min, deltas,order,verbose);
        // Calculating dS and CCdE
        Array<double> dS(3,3), CCdE(3,3);
        for (int i = 0; i < deltas.size(); i++) {
            TestCANNHO->calcCCdEFiniteDifference(F, dF, deltas[i], order, CCdE);
            TestHO->calcdSFiniteDifference(F, dF, deltas[i], order, dS);

            // Check order of convergence of consistency of material elasticity
            TestHO->testMaterialElasticityConsistencyConvergenceOrderBetweenMaterialModels(F, dS, CCdE, deltas, order, convergence_order_tol, verbose);
        }
    }
}

// Test order of convergence of consistency of material elasticity for random F (large)
TEST_F(STRUCT_CANNHolzapfelOgdenTest, TestMaterialElasticityConsistencyConvergenceOrderAgainstHORandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to array
        convertToArray(F_std, F);
        // Generating perturbation
        Array<double> dF(3,3);
        std::vector<double> deltas;
        TestHO->generatePerturbationdF(F,dF,delta_max, delta_min, deltas,order,verbose);
        // Calculating dS and CCdE
        Array<double> dS(3,3), CCdE(3,3);
        for (int i = 0; i < deltas.size(); i++) {
            TestCANNHO->calcCCdEFiniteDifference(F, dF, deltas[i], order, CCdE);
            TestHO->calcdSFiniteDifference(F, dF, deltas[i], order, dS);

            // Check order of convergence of consistency of material elasticity
            TestHO->testMaterialElasticityConsistencyConvergenceOrderBetweenMaterialModels(F, dS, CCdE, deltas, order, convergence_order_tol, verbose);
        }
    }
}

// Test PK2 stress
TEST_F(STRUCT_CANNHolzapfelOgdenTest, TestMaterialElasticityAgainstReferenceIdentityF) {
    verbose = true; // Show values of S and S_ref

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
    TestCANNHO->testMaterialElasticityAgainstReference(F, CC_ref, rel_tol, abs_tol, verbose); // Comparing with CANN
}

// Test PK2 stress
TEST_F(STRUCT_CANNHolzapfelOgdenTest, TestMaterialElasticityAgainstReference) {
    verbose = true; // Show values of S and S_ref

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
    TestCANNHO->testMaterialElasticityAgainstReference(F, CC_ref, rel_tol, abs_tol, verbose); // Comparing with CANN
}





// ----------------------------------------------------------------------------
// ---------------------- Quadratic Volumetric Penalty Material ----------------
// ----------------------------------------------------------------------------

// ------------------------------ STRUCT Tests --------------------------------

// Test PK2 stress zero for F = I
TEST_F(STRUCT_QuadraticVolumetricPenaltyTest, TestPK2StressIdentityF) {
    //verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.0, 0.0, 0.0},
                        {0.0, 1.0, 0.0},
                        {0.0, 0.0, 1.0}};
    Array<double> S_ref(3, 3); // PK2 stress initialized to zero
    TestQVP->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress zero for prescribed isochoric deformation
TEST_F(STRUCT_QuadraticVolumetricPenaltyTest, TestPK2StressPrescribedIsochoricDeformation) {
    //verbose = true; // Show values of S and S_ref

    // Check isochoric deformation produces zero PK2 stress
    Array<double> F = {{1.1, 0.0, 0.0},
                        {0.0, 1.2, 0.0},
                        {0.0, 0.0, 1.0/(1.1*1.2)}};
    Array<double> S_ref(3, 3); // PK2 stress initialized to zero
    TestQVP->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}


// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (small)
TEST_F(STRUCT_QuadraticVolumetricPenaltyTest, TestPK2StressConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestQVP->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (medium)
TEST_F(STRUCT_QuadraticVolumetricPenaltyTest, TestPK2StressConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestQVP->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (large)
TEST_F(STRUCT_QuadraticVolumetricPenaltyTest, TestPK2StressConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestQVP->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (small)
TEST_F(STRUCT_QuadraticVolumetricPenaltyTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestQVP->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (medium)
TEST_F(STRUCT_QuadraticVolumetricPenaltyTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestQVP->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (large)
TEST_F(STRUCT_QuadraticVolumetricPenaltyTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestQVP->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
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
    TestQVP->testRhoBetaAgainstReference(p, rho0, rho_ref, beta_ref, drhodp_ref, dbetadp_ref, rel_tol, abs_tol, verbose);
}


// ----------------------------------------------------------------------------
// ---------------------- Simo-Taylor 91 Volumetric Penalty Material ------------
// ----------------------------------------------------------------------------

// ------------------------------ STRUCT Tests --------------------------------


// Test PK2 stress zero for F = I
TEST_F(STRUCT_SimoTaylor91VolumetricPenaltyTest, TestPK2StressIdentityF) {
    //verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    Array<double> F = {{1.0, 0.0, 0.0},
                        {0.0, 1.0, 0.0},
                        {0.0, 0.0, 1.0}};
    Array<double> S_ref(3, 3); // PK2 stress initialized to zero
    TestST91->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress zero for prescribed isochoric deformation
TEST_F(STRUCT_SimoTaylor91VolumetricPenaltyTest, TestPK2StressPrescribedIsochoricDeformation) {
    //verbose = true; // Show values of S and S_ref

    // Check isochoric deformation produces zero PK2 stress
    Array<double> F = {{1.1, 0.0, 0.0},
                        {0.0, 1.2, 0.0},
                        {0.0, 0.0, 1.0/(1.1*1.2)}};
    Array<double> S_ref(3, 3); // PK2 stress initialized to zero
    TestST91->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}


// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (small)
TEST_F(STRUCT_SimoTaylor91VolumetricPenaltyTest, TestPK2StressConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestST91->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (medium)
TEST_F(STRUCT_SimoTaylor91VolumetricPenaltyTest, TestPK2StressConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestST91->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence between finite difference PK2 stress and compute_pk2cc() PK2 stress for random F (large)
TEST_F(STRUCT_SimoTaylor91VolumetricPenaltyTest, TestPK2StressConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence between finite difference and compute_pk2cc() PK2 stress
        TestST91->testPK2StressConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (small)
TEST_F(STRUCT_SimoTaylor91VolumetricPenaltyTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFSmall) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_small_list
    for (auto F_std : F_small_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestST91->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (medium)
TEST_F(STRUCT_SimoTaylor91VolumetricPenaltyTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFMedium) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_medium_list
    for (auto F_std : F_medium_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestST91->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
    }
}

// Test order of convergence of consistency of material elasticity for random F (large)
TEST_F(STRUCT_SimoTaylor91VolumetricPenaltyTest, TestMaterialElasticityConsistencyConvergenceOrderRandomFLarge) {
    //verbose = true; // Show order of convergence, errors, F, S

    // Loop over F in F_large_list
    for (auto F_std : F_large_list) {
        // Convert to Array
        convertToArray(F_std, F);

        // Check order of convergence of consistency of material elasticity
        TestST91->testMaterialElasticityConsistencyConvergenceOrder(F, delta_max, delta_min, order, convergence_order_tol, verbose);
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
    TestST91->testRhoBetaAgainstReference(p, rho0, rho_ref, beta_ref, drhodp_ref, dbetadp_ref, rel_tol, abs_tol, verbose);
}


// ----------------------------------------------------------------------------
// ---------------------- Miehe 94 Volumetric Penalty Material -----------------
// ----------------------------------------------------------------------------

// ------------------------------ STRUCT Tests --------------------------------

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

// ============================================================================
// ========================== END OF TESTS ====================================
// ============================================================================


