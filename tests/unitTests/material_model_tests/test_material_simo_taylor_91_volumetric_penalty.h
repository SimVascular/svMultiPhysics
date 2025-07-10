#ifndef TEST_MATERIAL_SIMO_TAYLOR_91_VOLUMETRIC_PENALTY_H
#define TEST_MATERIAL_SIMO_TAYLOR_91_VOLUMETRIC_PENALTY_H

#include "test_material_common.h"


// Class to contain Simo-Taylor 91 volumetric penalty parameters (just the penalty parameter)
class SimoTaylor91VolumetricPenaltyParams : public MatParams {
    public:
        double kappa;
    
        // Default constructor
        SimoTaylor91VolumetricPenaltyParams() : kappa(0.0) {}
    
        // Constructor with parameters
        SimoTaylor91VolumetricPenaltyParams(double kappa) : kappa(kappa) {}
    };

/**
 * @brief Class for testing the Simo-Taylor91 volumetric penalty model.
 *
 * This class provides methods to set up and test the Simo-Taylor91 volumetric penalty model, including 
 * computing the strain energy and printing material parameters.
 */
class TestSimoTaylor91VolumetricPenalty : public TestMaterialModel {
    public:
    
        /**
         * @brief Parameters for the volumetric penalty model.
         */
        SimoTaylor91VolumetricPenaltyParams params;
    
        /**
         * @brief Constructor for the TestSimoTaylor91VolumetricPenalty class.
         *
         * Initializes the volumetric penalty parameters for svMultiPhysics.
         *
         * @param[in] params_ Parameters for the volumetric penalty model.
         */
        TestSimoTaylor91VolumetricPenalty(const SimoTaylor91VolumetricPenaltyParams &params_) : TestMaterialModel( consts::ConstitutiveModelType::stIso_nHook, consts::ConstitutiveModelType::stVol_ST91),
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
         * @brief Computes the strain energy for the Simo-Taylor91 volumetric penalty model.
         *
         * @param[in] F Deformation gradient.
         * @return Strain energy density for the Simo-Taylor91 volumetric penalty model.
         */
        double computeStrainEnergy(const Array<double> &F) {
                
                // Compute solid mechanics terms
                solidMechanicsTerms smTerms = calcSolidMechanicsTerms(F);
        
                // Strain energy density for Simo-Taylor91 volumetric penalty model
                // Psi = kappa/4 * (J^2 - 1 - 2*ln(J))
                double Psi = params.kappa/4.0 * (pow(smTerms.J, 2) - 1.0 - 2.0 * log(smTerms.J));
        
                return Psi;
        }
    };
#endif