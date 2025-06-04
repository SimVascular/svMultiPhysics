
# **Problem Description**

Solve the same problem as in [fluid/pipe_RCR_3d](https://github.com/SimVascular/svMultiPhysics/tree/main/tests/cases/fluid/pipe_RCR_3d) replacing the RCR boundary condition with a resistance computed using the SimVascular [svZeroDSolver](https://simvascular.github.io/documentation/rom_simulation.html#0d-solver). 

# Introduction

The svZeroDSolver simulates bulk cardiovascular flow rates and pressures using an arbitrary zero-dimensional (0D) lumped parameter model (LPM) of a discrete network of components analogous to electrical circuits. It provides an Application Programming Interface (API) that allows it to communicate and interact with external software applications directly using function calls to programmatically define custom inflow and outflow boundary conditions for a CFD simulation. The svMultiPhysics solver can directly access the svZeroDSolver API by loading the svZeroDSolver as a shared (dynamic) library available after installing it from [SimTK](https://simtk.org/frs/?group_id=188) or building it from source.


# Enable using the svZeroDSolver

The svMultiPhysics solver XML `svZeroDSolver_interface` subsection under `Add_equation` enables using the svZeroDSolver. The following `svZeroDSolver_interface` parameters are used for this test 
```
<svZeroDSolver_interface> 
  <Coupling_type> semi-implicit </Coupling_type>
  <Configuration_file> svzerod_3Dcoupling.json </Configuration_file>
  <Shared_library> ../../../../svZeroDSolver/build/src/interface/libsvzero_interface.dylib </Shared_library>  
  <Initial_flows> 0.0 </Initial_flows>
  <Initial_pressures> 0.0 </Initial_pressures>
</svZeroDSolver_interface>
```
For more detailed documentation see the [Solver Parameter Input File](https://simvascular.github.io/documentation/multi_physics.html#solver-input-file
) `Equation / svZeroDSolver Interface Subsection`.

## Shared library
The svZeroDSolver application is built from source as part of the continuous integration (CI) testing. The shared library **libsvzero_interface.dylib** is created as part of the build. 

## Configuration file
The **svzerod_3Dcoupling.json** file is used to configure a svZeroDSolver simulation. It contains information describing the blocks (elements) of the LPM used to create the values for the resistance boundary condition communicated to the svMultiPhysics solver.

The `external_solver_coupling_blocks` section of the `svzerod_3Dcoupling.json` file identifies the names of the blocks (`name:` field) used in the svMultiPhysics solver XML file `Add_BC` keyword for boundary conditions coupled to the svZeroDSolver. A single blocka named **RCR_coupling** is used for this test. 

# Interfacing to the svZeroDSolver

An RCR boundary condition is defined for the **lumen_outlet** outlet face using the `Add_BC` keyword

```
<Add_BC name="lumen_outlet" > 
  <Type> Neu </Type> 
  <Time_dependence> Coupled </Time_dependence> 
  <svZeroDSolver_block> RCR_coupling </svZeroDSolver_block>  
</Add_BC> 
```
where 
- `Time_dependence` set to `Coupled` identifies this boundary condition as being coupled to the svZeroDSolver
- `svZeroDSolver_block` identifies the svZeroDSolver block used to generate values for this boiundary condition


