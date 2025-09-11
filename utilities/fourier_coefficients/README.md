# Fourier coefficients

The **fft_temporal_values.py**  script takes the **lumen_inlet.flow** file, computes the Fourier coefficients using a fft python function (which is identical to svMultiphysics' `fft` routine), and returns the **.fcs** file ready for simulation. It also provides a visualization of the fourier interpolation, which is shown below:

<p align="center">
   <img src="./fft_reconstruction.png" width="600">
</p>

An integration test using this **.fcs** file is available in [pipe_RCR_3d_fourier_coeff](../../tests/cases/fluid/pipe_RCR_3d_fourier_coeff/). 


