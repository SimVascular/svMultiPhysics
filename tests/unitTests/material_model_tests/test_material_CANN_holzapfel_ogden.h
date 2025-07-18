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

#ifndef TEST_MATERIAL_CANN_HOLZAPFEL_OGDEN_H
#define TEST_MATERIAL_CANN_HOLZAPFEL_OGDEN_H

#include "test_material_common.h"


// Class to contain CANN model with Holzapfel-Ogden material parameters
class CANN_HO_Params : public MatParams {
public:
    std::vector<CANNRow> Table;
    // Define fiber directions
    Vector<double> f;    // Fiber direction
    Vector<double> s;    // Sheet direction

    const int nrows = 4;

    // Default constructor
    CANN_HO_Params() {

        // Resize Table to ensure there's at least 4 elements
        Table.resize(nrows);  // Ensure there's space for 4 rows

        f.resize(3);
        s.resize(3);
      };

    // Constructor with parameters
    CANN_HO_Params(std::vector<CANNRow> TableValues, Vector<double> f, Vector<double> s) {
        Table = TableValues;

        this -> f = f;
        this -> s = s;
    };

};


/**
 * @brief Class for testing the CANN model with Holzapfel-Ogden material model parameters.
 *
 * This class provides methods to set up and test the Neo-Hookean material model, including 
 * computing the strain energy and printing material parameters.
 */
class TestCANN_HO : public TestMaterialModel {
public:

    /**
     * @brief Parameters for the CANN material model.
     */
    CANN_HO_Params params;

    /**
     * @brief Constructor for the TestCANN_HO class.
     *
     * Initializes the CANN - HO material parameters
     *
     * @param[in] params_ Parameters for the CANN HO material model.
     */
    TestCANN_HO(const CANN_HO_Params &params_) : TestMaterialModel( consts::ConstitutiveModelType::stArtificialNeuralNet, consts::ConstitutiveModelType::stVol_ST91),
        params(params_) 
        {
        // Set HO material parameters for svFSIplus
        auto &dmn = com_mod.mockEq.mockDmn;

        // Set number of rows in the table
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

        // Set number of fiber directions and fiber directions
        nFn = 2;
        fN.set_col(0, params.f);
        fN.set_col(1, params.s);
    }

/**
     * @brief Prints the CANN HO material parameters.
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
     * @brief Computes the strain energy for the Holzapfel-Ogden material model.
     *
     * @param[in] F Deformation gradient.
     * @return Strain energy density for the Neo-Hookean material model.
     */
    double computeStrainEnergy(const Array<double> &F) {
        // Compute solid mechanics terms
        solidMechanicsTerms smTerms = calcSolidMechanicsTerms(F);

        // Fiber and sheet directions
        Vector<double> f = params.f;
        Vector<double> s = params.s;

        // Strain energy density for Holzapfel-Ogden material model

        // Formulation with fully decoupled isochoric-volumetric split
        // Uses I1_bar, I4_bar_f, I4_bar_s, I8_bar_fs (bar = isochoric)
        // Psi = a/2b * exp{b(I1_bar - 3)} 
        //       + a_f/2b_f * chi(I4_bar_f) * (exp{b_f(I4_bar_f - 1)^2} - 1
        //       + a_s/2b_s * chi(I4_bar_s) * (exp{b_s(I4_bar_s - 1)^2} - 1
        //       + a_fs/2b_fs * (exp{b_fs*I8_bar_fs^2} - 1)
        // We set k = 0 which leads to chi = 1/2 for all terms.
        
        // Invariants
        double I1_bar = smTerms.Ib1;
        // I4_bar_f = f . C_bar . f
        auto C_bar_f = mat_fun::mat_mul(smTerms.C_bar, f);
        double I4_bar_f = f.dot(C_bar_f);
        // I4_bar_s = s . C_bar . s
        auto C_bar_s = mat_fun::mat_mul(smTerms.C_bar, s);
        double I4_bar_s = s.dot(C_bar_s);
        // I8_bar_fs = f . C_bar . s
        double I8_bar_fs = f.dot(C_bar_s);

        // Strain energy density for Holzapfel-Ogden material model with modified anisotropic invariants (bar quantities)
        double Psi = 0.0;
        int nterms = params.nrows;
        Psi += params.Table[0].weights[2] * exp(params.Table[0].weights[1] * (I1_bar - 3.0)); // isotropic term
        Psi += params.Table[1].weights[2] * (exp(params.Table[1].weights[1] * pow(I4_bar_f - 1.0, 2)) - 1.0);   // Fiber term; 0.5 included in params
        Psi += params.Table[2].weights[2] * (exp(params.Table[2].weights[1] * pow(I4_bar_s - 1.0, 2)) - 1.0);   // Sheet term
        Psi += params.Table[3].weights[2] * (exp(params.Table[3].weights[1] * pow(I8_bar_fs, 2)) - 1.0);                   // Cross-fiber term


        return Psi;
    }
};

#endif // TEST_MATERIAL_CANN_HOLZAPFEL_OGDEN_H