#ifndef TEST_MATERIAL_NEOHOOKEAN_H
#define TEST_MATERIAL_NEOHOOKEAN_H

#include "test_material_common.h"




// Class to contain Neo-Hookean material parameters
class NeoHookeanParams : public MatParams {
public:
    double C10;

    // Default constructor
    NeoHookeanParams() : C10(0.0) {}

    // Constructor with parameters
    NeoHookeanParams(double c10) : C10(c10) {}

};

/**
 * @brief Class for testing the Neo-Hookean material model.
 *
 * This class provides methods to set up and test the Neo-Hookean material model, including 
 * computing the strain energy and printing material parameters.
 */
class TestNeoHookean : public TestMaterialModel {
public:

    /**
     * @brief Parameters for the Neo-Hookean material model.
     */
    NeoHookeanParams params;

    /**
     * @brief Constructor for the TestNeoHookean class.
     *
     * Initializes the Neo-Hookean material parameters for svMultiPhysics.
     *
     * @param[in] params_ Parameters for the Neo-Hookean material model.
     */
    TestNeoHookean(const NeoHookeanParams &params_) : TestMaterialModel( consts::ConstitutiveModelType::stIso_nHook, consts::ConstitutiveModelType::stVol_ST91),
        params(params_) 
        {
        // Set Neo-Hookean material parameters for svMultiPhysics
        auto &dmn = com_mod.mockEq.mockDmn;
        dmn.stM.C10 = params.C10;
        dmn.stM.Kpen = 0.0;         // Zero volumetric penalty parameter
    }

    /**
     * @brief Prints the Neo-Hookean material parameters.
     */
    void printMaterialParameters() {
        std::cout << "C10 = " << params.C10 << std::endl;
    }

    /**
     * @brief Computes the strain energy for the Neo-Hookean material model.
     *
     * @param[in] F Deformation gradient.
     * @return Strain energy density for the Neo-Hookean material model.
     */
    double computeStrainEnergy(const Array<double> &F) {
        // Compute solid mechanics terms
        solidMechanicsTerms smTerms = calcSolidMechanicsTerms(F);

        // Strain energy density for Neo-Hookean material model
        // Psi_iso = C10 * (Ib1 - 3)
        double Psi_iso = params.C10 * (smTerms.Ib1 - 3.);

        return Psi_iso;
    }
};
#endif