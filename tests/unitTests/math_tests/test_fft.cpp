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

#include "fft.h"
#include "ComMod.h"
#include "../test_common.h"
#include <cmath>

class FFTTest : public ::testing::Test {
protected:
    void SetUp() override {}
    
    void TearDown() override {}

    // Generate temporal values for function f(t) = sin(t) + cos(t) + 0.1*t
    void CreateTemporalValues(int N, double x_start, double x_end, 
                             std::vector<std::vector<double>>& temporal_values) {
        temporal_values.clear();
        temporal_values.reserve(N);
        
        double step = (x_end - x_start) / (N - 1);
        for (int i = 0; i < N; ++i) {
            double t = x_start + i * step;
            double y = std::sin(t) + std::cos(t) + 0.1 * t;
            temporal_values.push_back({t, y});
        }
    }
    
    void InitializeFourierCoefficients(fcType& gt, int d = 1, int n = 16) {
        gt.d = d;
        gt.n = n;
        gt.qi.resize(gt.d);
        gt.qs.resize(gt.d);
        gt.r.resize(gt.d, gt.n);
        gt.i.resize(gt.d, gt.n);
    }
};

TEST_F(FFTTest, SinCosLinearCombination) {
    // Creates a temporal values function: f(t) = sin(t) + cos(t) + 0.1*t 
    // Test finds the interpolated fourier coefficients of this function using fft.cpp
    int N = 100;            // 100 timesteps
    double x_start = 0.0;   // start time
    double x_end = 10.0;    // end time 
    std::vector<std::vector<double>> temporal_values;
    
    CreateTemporalValues(N, x_start, x_end, temporal_values);
    
    fcType gt;
    // Apply interpolation to 1 dimensional data with n = 16 Fourier modes
    InitializeFourierCoefficients(gt, 1, 16);
    
    // Compute the Fourier coefficients
    fft(N, temporal_values, gt);

    // Check the slope (first Fourier coefficient)
    ASSERT_NEAR(gt.qs[0], -0.13830, 1e-2) << "Expected slope ~-0.13830";

    // Check the real and imaginary components of the first three Fourier coefficients
    ASSERT_NEAR(gt.r(0, 0), 0.32094, 1e-2) << "Expected first real coefficient to be close to 0.32094";
    ASSERT_NEAR(gt.i(0, 0), 0.0, 1e-2) << "Expected first imaginary coefficient to be close to 0.0";
    ASSERT_NEAR(gt.r(0, 1), 0.42759, 1e-2) << "Expected second real coefficient to be close to 0.42759";
    ASSERT_NEAR(gt.i(0, 1), 1.25295, 1e-2) << "Expected second imaginary coefficient to be close to 1.25295";
    ASSERT_NEAR(gt.r(0, 2), -0.44685, 1e-2) << "Expected third real coefficient to be close to -0.44685";
    ASSERT_NEAR(gt.i(0, 2), -0.65403, 1e-2) << "Expected third imaginary coefficient to be close to -0.65403";
}