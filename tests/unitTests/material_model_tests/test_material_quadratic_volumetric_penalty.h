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

#ifndef TEST_MATERIAL_QUADRATIC_VOLUMETRIC_PENALTY_H
#define TEST_MATERIAL_QUADRATIC_VOLUMETRIC_PENALTY_H

#include "test_material_common.h"

// Class to contain quadratic volumetric penalty parameters (just the penalty parameter)
class QuadraticVolumetricPenaltyParams : public MatParams {
    public:
        double kappa;
    
        // Default constructor
        QuadraticVolumetricPenaltyParams() : kappa(0.0) {}
    
        // Constructor with parameters
        QuadraticVolumetricPenaltyParams(double kappa) : kappa(kappa) {}
    };
    

    /**
 * @brief Class for testing the quadratic volumetric penalty model.
 *
 * This class provides methods to set up and test the quadratic volumetric penalty model, including 
 * computing the strain energy and printing material parameters.
 */
class TestQuadraticVolumetricPenalty : public TestMaterialModel {
public:

    /**
     * @brief Parameters for the volumetric penalty model.
     */
    QuadraticVolumetricPenaltyParams params;

    /**
     * @brief Constructor for the TestQuadraticVolumetricPenalty class.
     *
     * Initializes the volumetric penalty parameters for svMultiPhysics.
     *
     * @param[in] params_ Parameters for the volumetric penalty model.
     */
    TestQuadraticVolumetricPenalty(const QuadraticVolumetricPenaltyParams &params_) : TestMaterialModel( consts::ConstitutiveModelType::stIso_nHook, consts::ConstitutiveModelType::stVol_Quad),
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
     * @brief Computes the strain energy for the quadratic volumetric penalty model.
     *
     * @param[in] F Deformation gradient.
     * @return Strain energy density for the quadratic volumetric penalty model.
     */
    double computeStrainEnergy(const Array<double> &F) {
            
            // Compute solid mechanics terms
            solidMechanicsTerms smTerms = calcSolidMechanicsTerms(F);
    
            // Strain energy density for quadratic volumetric penalty model
            // Psi = kappa/2 * (J - 1)^2  
            double Psi = params.kappa/2.0 * pow(smTerms.J - 1.0, 2);
    
            return Psi;
    }
};
#endif