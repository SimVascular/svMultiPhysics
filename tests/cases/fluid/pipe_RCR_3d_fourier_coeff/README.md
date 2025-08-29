
# **Problem Description**

Simulate unsteady fluid flow in a pipe. This case is identical to <a href="https://github.com/SimVascular/svFSIplus/tree/main/tests/cases/fluid/pipe_RCR_3d"> Fluid RCR 3D Pipe </a>, except that fourier coefficients are read by the solver (instead of the usual temporal values file). 

# Inlet flow boundary condition

Interpolated Fourier coefficients are provided for the **lumen_inlet** boundary condition. These coefficients are specified in the **lumen_inlet.fcs** file, which is an alternative to the flow data that is provided in the **lumen_inlet.flow**. Providing the fourier coefficient file skips the fourier interpolation function (**fft.cpp**) in svMultiphysics. 

```
<Add_BC name="lumen_inlet" > 
    <Type> Dir </Type> 
    <Time_dependence> Unsteady </Time_dependence> 
    <Fourier_coefficients_file_path> lumen_inlet.fcs </Fourier_coefficients_file_path> 
    <Profile> Parabolic </Profile> 
    <Impose_flux> true </Impose_flux> 
</Add_BC> 
```

## File format for Fourier coefficients

```
<Line 1>  <first_time_point>    <last_time_point>
<Line 2>    <initial_value>     <linear_slope>
<Line 3>  <number_of_Fourier_modes>
<Line 4>  <real_Fourier_mode_1>    <imag_Fourier_mode_1>
<Line 5>  <real_Fourier_mode_2>   <imag_Fourier_mode_2>  
...
<Line N+3>  <real_Fourier_mode_N>   <imag_Fourier_mode_N>  
```

This file can be generated using the **fft_temporal_values.py**. This script takes the lumen_inlet.flow file, computes the Fourier coefficients using an identical fft python function, and returns the **.fcs** file ready for simulation. It also provides a visualization of the fourier interpolation, which is shown below:

<p align="center">
   <img src="./fft_reconstruction.png" width="600">
</p>


# Outlet RCR boundary condition

An RCR boundary condition is defined for the **lumen_outlet** outlet face.

```
<Add_BC name="lumen_outlet" > 
  <Type> Neu </Type> 
  <Time_dependence> RCR </Time_dependence> 
  <RCR_values> 
    <Capacitance> 1.5e-5 </Capacitance> 
    <Distal_resistance> 1212 </Distal_resistance> 
    <Proximal_resistance> 121 </Proximal_resistance> 
     <Distal_pressure> 0 </Distal_pressure> 
    <Initial_pressure> 0 </Initial_pressure> 
  </RCR_values> 
</Add_BC>
```
