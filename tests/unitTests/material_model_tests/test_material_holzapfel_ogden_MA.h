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

#ifndef TEST_MATERIAL_HOLOZAPFEL_OGDEN_MA_H
#define TEST_MATERIAL_HOLOZAPFEL_OGDEN_MA_H

#include "test_material_common.h"


// Class to contain Holzapfel-Ogden (Modified Anisotropy) material parameters
class HolzapfelOgdenMAParams : public MatParams {
public:
    double a;    
    double b;
    double a_f;
    double b_f;
    double a_s;
    double b_s;
    double a_fs;
    double b_fs;
    Vector<double> f;    // Fiber direction
    Vector<double> s;    // Sheet direction

    double k; // Smoothed Heaviside function parameter

    // Default constructor
    HolzapfelOgdenMAParams() : a(0.0), b(0.0), a_f(0.0), b_f(0.0), a_s(0.0), b_s(0.0), a_fs(0.0), b_fs(0.0), k(0.0) {
        f.resize(3);
        s.resize(3);
    }

    // Constructor with parameters
    HolzapfelOgdenMAParams(double a, double b, double a_f, double b_f, double a_s, double b_s, double a_fs, double b_fs, double k, Vector<double> f, Vector<double> s) : a(a), b(b), a_f(a_f), b_f(b_f), a_s(a_s), b_s(b_s), a_fs(a_fs), b_fs(b_fs), k(k), f(f), s(s) {
    }
};

/**
 * @brief Class for testing the Holzapfel-Ogden (Modified Anisotropy) material model. 
 * 
 * This class provides methods to set up and test the Holzapfel-Ogden-ma material 
 * model, including computing the strain energy and printing material parameters.
 *
 */
class TestHolzapfelOgdenMA : public TestMaterialModel {
public:

    /**
     * @brief Parameters for the Holzapfel-Ogden ma material model.
     */
    HolzapfelOgdenMAParams params;

    /**
     * @brief Constructor for the TestHolzapfelOgdenMA class.
     *
     * Initializes the Holzapfel-Ogden material parameters for svMultiPhysics.
     *
     * @param[in] params_ Parameters for the Holzapfel-Ogden ma material model.
     */
    TestHolzapfelOgdenMA(const HolzapfelOgdenMAParams &params_) : TestMaterialModel( consts::ConstitutiveModelType::stIso_HO_ma, consts::ConstitutiveModelType::stVol_ST91),
        params(params_) 
        {
        // Set Holzapfel-Ogden material parameters for svMultiPhysics
        auto &dmn = com_mod.mockEq.mockDmn;
        dmn.stM.a = params.a;
        dmn.stM.b = params.b;
        dmn.stM.aff = params.a_f;
        dmn.stM.bff = params.b_f;
        dmn.stM.ass = params.a_s;
        dmn.stM.bss = params.b_s;
        dmn.stM.afs = params.a_fs;
        dmn.stM.bfs = params.b_fs;
        dmn.stM.khs = params.k;     // Smoothed Heaviside function parameter
        dmn.stM.Kpen = 0.0;         // Zero volumetric penalty parameter

        // Set number of fiber directions and fiber directions
        nFn = 2;
        fN.set_col(0, params.f);
        fN.set_col(1, params.s);
    }

    /**
     * @brief Prints the Holzapfel-Ogden material parameters.
     */
    void printMaterialParameters() {
        std::cout << "a = " << params.a << std::endl;
        std::cout << "b = " << params.b << std::endl;
        std::cout << "a_f = " << params.a_f << std::endl;
        std::cout << "b_f = " << params.b_f << std::endl;
        std::cout << "a_s = " << params.a_s << std::endl;
        std::cout << "b_s = " << params.b_s << std::endl;
        std::cout << "a_fs = " << params.a_fs << std::endl;
        std::cout << "b_fs = " << params.b_fs << std::endl;
        std::cout << "k = " << params.k << std::endl;
        std::cout << "f = " << "[" << params.f[0] << " " << params.f[1] << " " << params.f[2] << "]" << std::endl;
        std::cout << "s = " << "[" << params.s[0] << " " << params.s[1] << " " << params.s[2] << "]" << std::endl;
    }

    /**
     * @brief Smoothed Heaviside function centered at x = 1.
     * 
     * @param[in] x Input value.
     * @param[in] k Smoothing parameter.
     * @return Smoothed Heaviside function.
     */
    double chi(const double x, const double k=100) const {
        return 1. / (1. + exp(-k * (x - 1.)));
    }

    /**
     * @brief Computes the strain energy for the Holzapfel-Ogden material model.
     *
     * @param[in] F Deformation gradient.
     * @return Strain energy density for the Holzapfel-Ogden material model.
     */
    double computeStrainEnergy(const Array<double> &F) {
        // Compute solid mechanics terms
        solidMechanicsTerms smTerms = calcSolidMechanicsTerms(F);

        // Material parameters
        double a = params.a;
        double b = params.b;
        double a_f = params.a_f;
        double b_f = params.b_f;
        double a_s = params.a_s;
        double b_s = params.b_s;
        double a_fs = params.a_fs;
        double b_fs = params.b_fs;

        // Smoothed Heaviside parameter
        double k = params.k;

        // Fiber and sheet directions
        Vector<double> f = params.f;
        Vector<double> s = params.s;

        // Strain energy density for Holzapfel-Ogden material model

        // Formulation used by cardiac mechanics benchmark paper (Arostica et al., 2024)
        // Uses I1_bar (bar = isochoric), but I4_f, I4_s, I8_fs (not bar)
        // Psi = a/2b * exp{b(I1_bar - 3)} 
        //       + a_f/2b_f * chi(I4_f) * (exp{b_f(I4_f - 1)^2} - 1
        //       + a_s/2b_s * chi(I4_s) * (exp{b_s(I4_s - 1)^2} - 1
        //       + a_fs/2b_fs * (exp{b_fs*I8_fs^2} - 1)
        // This corresponds to the HO-ma (modified anisotropy) implementation in svMultiPhysics

        // Invariants
        double I1_bar = smTerms.Ib1;
        // I4_f = f . C . f
        auto C_f = mat_fun::mat_mul(smTerms.C, f);
        double I4_f = f.dot(C_f);
        // I4_s = s . C . s
        auto C_s = mat_fun::mat_mul(smTerms.C, s);
        double I4_s = s.dot(C_s);
        // I8_fs = f . C . s
        double I8_fs = f.dot(C_s);

        // Strain energy density for Holzapfel-Ogden material model with full anisotropic invariants
        double Psi = 0.0;
        Psi += a / (2.0 * b) * exp(b * (I1_bar - 3.0));                             // Isotropic term
        Psi += a_f / (2.0 * b_f) * chi(I4_f, k) * (exp(b_f * pow(I4_f - 1.0, 2)) - 1.0);   // Fiber term
        Psi += a_s / (2.0 * b_s) * chi(I4_s, k) * (exp(b_s * pow(I4_s - 1.0, 2)) - 1.0);   // Sheet term
        Psi += a_fs / (2.0 * b_fs) * (exp(b_fs * pow(I8_fs, 2)) - 1.0);                   // Cross-fiber term

        return Psi;

    }
};


#endif // TEST_MATERIAL_HOLOZAPFEL_OGDEN_MA_H