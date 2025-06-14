# Problem Description
This test case simulates an idealized model of a left ventricle (LV) using a NeoHookean material. The LV is inflated at an approximately constant rate of change of volume using a boundary condition coupled to the SimVascular [svZeroDSolver](https://simvascular.github.io/documentation/rom_simulation.html#0d-solver). 

# svZeroDSolver
The svZeroDSolver simulates bulk cardiovascular flow rates and pressures using an arbitrary zero-dimensional (0D) lumped parameter model (LPM) of a discrete network of components analogous to electrical circuits.  The svMultiPhysics solver can directly access the svZeroDSolver by loading the svZeroDSolver as a shared (dynamic) library available after installing it from [SimTK](https://simtk.org/frs/?group_id=188) or building it from source.

The svZeroDSolver application is built from source as part of the continuous integration (CI) testing. The shared library **libsvzero_interface.dylib** is created as part of the build. 

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


# Configuration file
The **svzerod_3Dcoupling.json** file is used to configure a svZeroDSolver simulation. It contains information describing the blocks (elements) of the LPM used to create the values for the resistance boundary condition communicated to the svMultiPhysics solver.

The `external_solver_coupling_blocks` section of the `svzerod_3Dcoupling.json` file identifies the names of the blocks (`name:` field) used in the svMultiPhysics solver XML file `Add_BC` keyword for boundary conditions coupled to the svZeroDSolver. A single block named **LV_IN** is used for this test. 

The `boundary_conditions` section of the **svzerod_3Dcoupling.json** file
```
    "boundary_conditions": [
        {
            "bc_name": "P_SOURCE",
            "bc_type": "PRESSURE",
            "bc_values": {
                "P": [0.0, 1.0E7, 1.0E7],
                "t": [0.0, 0.1, 1.0]
            }
        }
    ],
```
defines a **P_SOURCE** pressure boundary condition that ramps from 0.0 to 1.0E7 and plateaus at 1.0E7 for the time period **t**.

The `vessels` section  
```
    "vessels": [
        {
            "boundary_conditions": {
                "inlet": "P_SOURCE"
            },
            "vessel_id": 0,
            "vessel_length": 10.0,
            "vessel_name": "branch0_seg0",
            "zero_d_element_type": "BloodVessel",
            "zero_d_element_values": {
                "R_poiseuille": 1.0E5
            }
        }
    ]
```
defines a blood vessel block used as a large resistor. It receives pressure from **P_SOURCE** and has a **R_poiseuille** resistance of 1.0E5. 

# Coupling a boundary condition to the svZeroDSolver
The LV is inflated using a boundary condition on the endocardium, the innermost layer of the LV model, defined as an LPM consisting of a large pressure source and large resistor computed using the svZeroDSolver.

The boundary condition is defined for the **endo** endocardium surface using the `Add_BC` keyword
```
<Add_BC name="endo" >
  <Type> Neu </Type>
  <Time_dependence> Coupled </Time_dependence>
  <Follower_pressure_load> true </Follower_pressure_load>
  <svZeroDSolver_block> LV_IN </svZeroDSolver_block>  
</Add_BC>
```
where 
- `Time_dependence` set to `Coupled` identifies this boundary condition as being coupled to the svZeroDSolver
- `svZeroDSolver_block` identifies the svZeroDSolver block used to generate values for this boiundary condition
  



       
