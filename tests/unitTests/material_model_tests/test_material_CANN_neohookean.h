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

#ifndef TEST_MATERIAL_CANN_NEOHOOKEAN_H
#define TEST_MATERIAL_CANN_NEOHOOKEAN_H

#include "test_material_common.h"


// Class to contain CANN model with Neo-Hookean material parameters
class CANN_NH_Params : public MatParams {
public:
    std::vector<CANNRow> Table;

    const int nrows = 1;

    // Default constructor
    CANN_NH_Params() {

        // Resize Table to ensure there's at least 1 element
        Table.resize(nrows);  // Ensure there's space for at least one row

        Table[0].invariant_index.value_ = 1;
        Table[0].activation_functions.value_ = {1,1,1};
        Table[0].weights.value_ = {1.0,1.0,40.0943265e6};
      };

    // Constructor with parameters
    CANN_NH_Params(std::vector<CANNRow> TableValues) {
        Table = TableValues;
    };

};

/**
 * @brief Class for testing the CANN model with Neo-Hookean material model parameters.
 *
 * This class provides methods to set up and test the Neo-Hookean material model, including 
 * computing the strain energy and printing material parameters.
 */
class TestCANN_NH : public TestMaterialModel {
public:

    /**
     * @brief Parameters for the CANN material model.
     */
    CANN_NH_Params params;

    /**
     * @brief Constructor for the TestCANN_NH class.
     *
     * Initializes the CANN - NeoHooke material parameters.
     *
     * @param[in] params_ Parameters for the CANN Neo-Hookean material model.
     */
    TestCANN_NH(const CANN_NH_Params &params_) : TestMaterialModel( consts::ConstitutiveModelType::stArtificialNeuralNet, consts::ConstitutiveModelType::stVol_ST91),
        params(params_) 
        {
        // Set Neo-Hookean material parameters
        auto &dmn = com_mod.mockEq.mockDmn;

        dmn.stM.paramTable.num_rows = params.nrows;
        
        // Resize Arrays and Vectors to ensure there is enough space
        dmn.stM.paramTable.invariant_indices.resize(dmn.stM.paramTable.num_rows);
        dmn.stM.paramTable.activation_functions.resize(dmn.stM.paramTable.num_rows,3);
        dmn.stM.paramTable.weights.resize(dmn.stM.paramTable.num_rows,3);

        // Populate components of the table in stM
        for (size_t i = 0; i < dmn.stM.paramTable.num_rows; i++)
        {
            // Store invariant index
            dmn.stM.paramTable.invariant_indices[i] = params.Table[i].invariant_index.value_;

            // Store activation function values
            dmn.stM.paramTable.activation_functions(i,0) = params.Table[i].activation_functions.value_[0];
            dmn.stM.paramTable.activation_functions(i,1) = params.Table[i].activation_functions.value_[1];
            dmn.stM.paramTable.activation_functions(i,2) = params.Table[i].activation_functions.value_[2];

            // Store weight values
            dmn.stM.paramTable.weights(i,0) = params.Table[i].weights.value_[0];
            dmn.stM.paramTable.weights(i,1) = params.Table[i].weights.value_[1];
            dmn.stM.paramTable.weights(i,2) = params.Table[i].weights.value_[2];

        }

        dmn.stM.Kpen = 0.0;         // Zero volumetric penalty parameter
    }

/**
     * @brief Prints the CANN Neo-Hookean material parameters.
     */
    void printMaterialParameters() {
        for (int i = 0; i < params.nrows; i++){
            std::cout << "ROW: " << i+1 << std::endl;
            std::cout << "Invariant number: " << params.Table[i].invariant_index << std::endl;
            std::cout << "Activation function 0: " << params.Table[i].activation_functions.value()[0] << std::endl;
            std::cout << "Activation function 1: " << params.Table[i].activation_functions.value()[1] << std::endl;
            std::cout << "Activation function 2: " << params.Table[i].activation_functions.value()[2] << std::endl;
            std::cout << "Weight 0: " << params.Table[i].weights[0] << std::endl;
            std::cout << "Weight 1: " << params.Table[i].weights[1] << std::endl;
            std::cout << "Weight 2: " << params.Table[i].weights[2] << std::endl;
        }
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
        double Psi_iso = params.Table[0].weights[2] * (smTerms.Ib1 - 3.); //w[0][6] = C10

        return Psi_iso;
    }
};


#endif // TEST_MATERIAL_CANN_NEOHOOKEAN_H