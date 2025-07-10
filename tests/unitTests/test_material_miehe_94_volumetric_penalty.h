#ifndef TEST_MATERIAL_MIEHE_94_VOLUMETRIC_PENALTY_H
#define TEST_MATERIAL_MIEHE_94_VOLUMETRIC_PENALTY_H

#include "test_material_common.h"


// Class to contain Miehe 94 volumetric penalty parameters (just the penalty parameter)
class Miehe94VolumetricPenaltyParams : public MatParams {
    public:
        double kappa;
    
        // Default constructor
        Miehe94VolumetricPenaltyParams() : kappa(0.0) {}
    
        // Constructor with parameters
        Miehe94VolumetricPenaltyParams(double kappa) : kappa(kappa) {}
    };

/**
 * @brief Class for testing the Miehe94 volumetric penalty model.
 *
 * This class provides methods to set up and test the Miehe94 volumetric penalty model, including 
 * computing the strain energy and printing material parameters.
 */
class TestMiehe94VolumetricPenalty : public TestMaterialModel {
    public:
    
        /**
         * @brief Parameters for the volumetric penalty model.
         */
        Miehe94VolumetricPenaltyParams params;
    
        /**
         * @brief Constructor for the TestMiehe94VolumetricPenalty class.
         *
         * Initializes the volumetric penalty parameters for svMultiPhysics.
         *
         * @param[in] params_ Parameters for the volumetric penalty model.
         */
        TestMiehe94VolumetricPenalty(const Miehe94VolumetricPenaltyParams &params_) : TestMaterialModel( consts::ConstitutiveModelType::stIso_nHook, consts::ConstitutiveModelType::stVol_M94),
            params(params_) 
            {
    
            // Set volumetric penalty parameter for svMultiPhysics
            auto &dmn = com_mod.mockEq.mockDmn;
            dmn.stM.Kpen = params.kappa;         // Volumetric penalty parameter
    
            // Note: Use Neo-Hookean material model for isochoric part, but set parameters to zero
            dmn.stM.C10 = 0.0;         // Zero Neo-Hookean parameter
        }
    
        /**
         * @brief Prints the volumetric penalty parameters.
         */
        void printMaterialParameters() {
            std::cout << "kappa = " << params.kappa << std::endl;
        }
    
        /**
         * @brief Computes the strain energy for the Miehe94 volumetric penalty model.
         *
         * @param[in] F Deformation gradient.
         * @return Strain energy density for the Miehe94 volumetric penalty model.
         */
        double computeStrainEnergy(const Array<double> &F) {
                
                // Compute solid mechanics terms
                solidMechanicsTerms smTerms = calcSolidMechanicsTerms(F);
        
                // Strain energy density for Miehe94 volumetric penalty model
                // Psi = kappa * (J - ln(J) - 1)
                double Psi = params.kappa * (smTerms.J - log(smTerms.J) - 1.0);
        
                return Psi;
        }
    };
#endif