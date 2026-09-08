# Unit Tests for svMultiPhysics

This directory contains unit tests for svMultiPhysics

## Building and Running Tests

### Prerequisites
- svMultiPhysics must be built with unit tests enabled:
  ```bash
  cmake -DENABLE_UNIT_TEST=ON ..
  ```

### Building Tests
1. Navigate to the build directory:
   ```bash
   cd <svMultiPhysics_root_directory>/build/svMultiPhysics-build/Source/solver
   ```

2. Build the tests:
   ```bash
   make
   ```

### Running Tests
Run all unit tests:
```bash
ctest --verbose
```

Or run the test executable directly:
```bash
./run_all_unit_tests
```

## Cellular ODE Trajectory Tests

The `active_stress_tests/` and `ionic_model_tests/` directories contain
standalone trajectory tests for cellular ODE models. Each test advances one
model with prescribed inputs and compares selected trajectory checkpoints with
trusted data in `reference_data/`. Source provenance and regeneration
instructions are provided in `reference_generators/`.

Current coverage includes the Nash-Panfilov and Regazzoni active-stress models,
and the Aliev-Panfilov, FitzHugh-Nagumo, Bueno-Orovio, and ten
Tusscher-Panfilov ionic models. Run all of these trajectory tests from the
repository root with:

```bash
build/svMultiPhysics-build/bin/run_all_unit_tests \
  --gtest_filter='ActiveStressTrajectory.*:IonicModelTrajectory.*'
```

## Common Files
### `test_common.h`
- Contains mock objects for svMultiPhysics components
- Provides `TestBase` class with common test infrastructure
- Includes Google Test framework



## Material Model Testing

Each material model follows this pattern:

1. **Parameter Class** (e.g., `NeoHookeanParams`)
   - Inherits from `MatParams`
   - Contains material-specific parameters

2. **Test Material Class** (e.g., `TestNeoHookean`)
   - Inherits from `TestMaterialModel`
   - Implements required virtual methods:
     - `printMaterialParameters()`
     - `computeStrainEnergy()`
   - Sets up material parameters for svMultiPhysics

3. **Test Fixture Class** (e.g., `NeoHookeanTest`)
   - Inherits from `MaterialTestFixture`
   - Sets up test parameters and objects

4. **Test Cases**
   - Use Google Test macros (`TEST_F`)
   - Test the material model second Piola-Kirchoff stress and material elasticity tensor in a variety of ways

### Material Models Currently Tested

#### Hyperelastic Models
- **Neo-Hookean**: Basic hyperelastic model
- **Mooney-Rivlin**: Two-parameter hyperelastic model
- **Holzapfel-Ogden**: Anisotropic hyperelastic model for soft tissues
- **Holzapfel-Ogden Modified Anisotropy**: Enhanced version with modified anisotropy

#### Artificial Neural Network Models
- **CANN Neo-Hookean**: Neural network approximation of Neo-Hookean behavior
- **CANN Holzapfel-Ogden**: Neural network approximation of Holzapfel-Ogden behavior

#### Volumetric Penalty Models
- **Quadratic**: Simple quadratic volumetric penalty
- **Simo-Taylor 91**: Simo-Taylor volumetric penalty model
- **Miehe 94**: Miehe volumetric penalty model


### Adding a New Material Model

To add a new material model for testing:

#### 1. Create Header File (`test_material_newmodel.h`)
```cpp
#ifndef TEST_MATERIAL_NEWMODEL_H
#define TEST_MATERIAL_NEWMODEL_H

#include "test_material_common.h"

// Parameter class
class NewModelParams : public MatParams {
public:
    double param1, param2;
    // Constructors...
};

// Test material class
class TestNewModel : public TestMaterialModel {
public:
    NewModelParams params;
    
    TestNewModel(const NewModelParams &params_) : 
        TestMaterialModel(consts::ConstitutiveModelType::stIso_NewModel, 
                         consts::ConstitutiveModelType::stVol_ST91),
        params(params_) {
        // Set parameters for svMultiPhysics
    }
    
    void printMaterialParameters() override {
        // Print parameters
    }
    
    double computeStrainEnergy(const Array<double> &F) override {
        // Compute strain energy density
    }
};

#endif
```

#### 2. Create Source File (`test_material_newmodel.cpp`)
```cpp
#include "test_material_newmodel.h"

class NewModelTest : public MaterialTestFixture {
protected:
    NewModelParams params;
    TestNewModel* TestNM;
    
    void SetUp() override {
        MaterialTestFixture::SetUp();
        // Set up parameters
        TestNM = new TestNewModel(params);
    }
    
    void TearDown() override {
        delete TestNM;
    }
};

// Test cases
TEST_F(NewModelTest, TestPK2StressIdentityF) {
    // Test implementation
}
```


#### 3. Automatic Discovery of Test Files
The CMakeLists.txt automatically finds all `.cpp` files under `unitTests/`, so no build configuration changes are needed.

