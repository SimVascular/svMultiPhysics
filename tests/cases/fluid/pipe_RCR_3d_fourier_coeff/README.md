
# **Problem Description**

Simulate unsteady fluid flow in a pipe. This case is identical to <a href="https://github.com/SimVascular/svFSIplus/tree/main/tests/cases/fluid/pipe_RCR_3d"> Fluid RCR 3D Pipe </a>, except that fourier coefficients are read by the solver (instead of the usual temporal values file). For additional validation, the **results_002.vtu** output file used in this case is the same as the `pipe_RCR_3d` case. 

# Inlet flow boundary condition

Interpolated Fourier coefficients are provided for the **lumen_inlet** boundary condition. These coefficients are specified in the **lumen_inlet.fcs** file, which is an alternative to the flow data that is provided in the **lumen_inlet.flow**. Providing the fourier coefficient file skips the fourier interpolation function (**fft.cpp**) in svMultiphysics. Documentation on the file format can be accessed [here](https://simvascular.github.io/documentation/multi_physics.html#data_file_formats_boundary_condition_fourier) and the `lumen_inlet` boundary condition block is provided in the [solver.xml](./solver.xml). 

The **lumen_inlet.fcs** file can be generated using the **fft_temporal_values.py** which is available in the `utilities` directory ([here](../../../../utilities/fourier_coefficients/)). This script takes the lumen_inlet.flow file, computes the Fourier coefficients using an identical fft python function, and returns the **.fcs** file ready for simulation. It also provides a visualization of the fourier interpolation, which is shown below:

<p align="center">
   <img src="./fft_reconstruction.png" width="600">
</p>


# Outlet RCR boundary condition

An RCR boundary condition is defined for the **lumen_outlet** outlet face.
