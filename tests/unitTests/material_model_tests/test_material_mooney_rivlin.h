#ifndef TEST_MATERIAL_MOONEY_RIVLIN_H
#define TEST_MATERIAL_MOONEY_RIVLIN_H

#include "test_material_common.h"



// Class to contain Mooney-Rivlin material parameters
class MooneyRivlinParams : public MatParams {
public:
    double C01;
    double C10;

    // Default constructor
    MooneyRivlinParams() : C01(0.0), C10(0.0) {}

    // Constructor with parameters
    MooneyRivlinParams(double c01, double c10) : C01(c01), C10(c10) {}

};


/**
 * @brief Class for testing the Mooney-Rivlin material model.
 *
 * This class provides methods to set up and test the Mooney-Rivlin material model, including 
 * computing the strain energy and printing material parameters.
 */
class TestMooneyRivlin : public TestMaterialModel {
public:

    /**
     * @brief Parameters for the Mooney-Rivlin material model.
     */
    MooneyRivlinParams params;

    /**
     * @brief Constructor for the TestMooneyRivlin class.
     *
     * Initializes the Mooney-Rivlin material parameters for svMultiPhysics.
     *
     * @param[in] params_ Parameters for the Mooney-Rivlin material model.
     */
    TestMooneyRivlin(const MooneyRivlinParams &params_) : TestMaterialModel( consts::ConstitutiveModelType::stIso_MR, consts::ConstitutiveModelType::stVol_ST91),
        params(params_) 
        {
        // Set Mooney-Rivlin material parameters for svMultiPhysics
        auto &dmn = com_mod.mockEq.mockDmn;
        dmn.stM.C01 = params.C01;
        dmn.stM.C10 = params.C10;
        dmn.stM.Kpen = 0.0;         // Zero volumetric penalty parameter
    }

    /**
     * @brief Prints the Mooney-Rivlin material parameters.
     */
    void printMaterialParameters() {
        std::cout << "C01 = " << params.C01 << ", C10 = " << params.C10 << std::endl;
    }

    /**
     * @brief Computes the strain energy for the Mooney-Rivlin material model.
     *
     * @param[in] F Deformation gradient.
     * @return Strain energy density for the Mooney-Rivlin material model.
     */
    double computeStrainEnergy(const Array<double> &F) {
        // Compute solid mechanics terms
        solidMechanicsTerms smTerms = calcSolidMechanicsTerms(F);

        // Strain energy density for Mooney-Rivlin material model
        // Psi_iso = C10 * (Ib1 - 3) + C01 * (Ib2 - 3)
        double Psi_iso = params.C10 * (smTerms.Ib1 - 3.) + params.C01 * (smTerms.Ib2 - 3.);

        return Psi_iso;
    }
};
#endif // TEST_MATERIAL_MOONEY_RIVLIN_H