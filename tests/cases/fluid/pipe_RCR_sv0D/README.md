
# **Problem Description**

Solve the same problem as in [fluid/pipe_RCR_3d](../pipe_RCR_3D) replacing the RCR boundary condition with a resistance computed 
using the SimVascular [svZeroDSolver](https://simvascular.github.io/documentation/rom_simulation.html#0d-solver). 

## Introduction

The svZeroDSolver simulates bulk cardiovascular flow rates and pressures using an arbitrary zero-dimensional (0D) lumped parameter model (LPM) of a discrete network of components analogous to electrical circuits. It provides an Application Programming Interface (API) that allows it to communicate and interact with external software applications directly using function calls to programmatically define custom inflow and outflow boundary conditions for a CFD simulation. The svMultiPhysics solver can directly access the svZeroDSolver API by loading the svZeroDSolver as a shared (dynamic) library available after installing the svZeroDSolver.


### Build svZeroDSolver
Importantly, to automatically run test cases with `pytest` (see below), you need to build `svZeroDSolver` in the folder
```
./svZeroDSolver/build
``` 
in the repository root.

To do so, you can run the following in the svMultiPhysics repository root:
```
git clone https://github.com/SimVascular/svZeroDSolver.git
cd svZeroDSolver
mkdir build
cd build
cmake ..
make -j2
``` 

## Configuration of sv0DSolver

The following files require user's attention: [svFSI.xml](./svFSI.xml), [svzerod_3Dcoupling.json](./svzerod_3Dcoupling.json) and [svZeroD_interface.dat](./svZeroD_interface.dat).

### svFSI.xml

The input file [svFSI_genBC.xml](./svFSI.xml) follows the master input file as a template. Some specific input options are discussed below:

```
   <Couple_to_svZeroD type="SI">
   </Couple_to_svZeroD>
```

This tells the solver that the 0d models will be calculated through sv0DSolver. Options to couple 0D codes with svFSI are `N`: none; `I`: implicit; `SI`: semi-implicit; `E`: explicit.

```
   <Add_BC name="lumen_inlet" > 
      <Type> Dir </Type> 
      <Time_dependence> Unsteady </Time_dependence> 
      <Temporal_values_file_path> lumen_inlet.flw</Temporal_values_file_path> 
      <Zero_out_perimeter> true </Zero_out_perimeter> 
      <Impose_flux> true </Impose_flux> 
   </Add_BC> 

   <Add_BC name="lumen_outlet" > 
      <Type> Neu </Type> 
      <Time_dependence> Coupled </Time_dependence> 
   </Add_BC> 
```

In this example, we use the LPN for just the outlet RCR boundary condition and use a file to specify the inlet flow conditions.

### svzerod_3Dcoupling.json

This is the configuration file for sv0DSolver and contains the elements of the 0D model being coupled to the 3D simulation. 

For more information on the available parameters and elements, documentation is available here: [svZeroDSolver](https://github.com/SimVascular/svZeroDSolver)

**The following are necessary in "simulation_parameters" for a coupled simulation:**
"coupled_simulation": true,
"steady_initial": false

The external coupling block is what connects the 3D element to the 0D model. sv0D allows you to create a name for this element and specify its type (in this case, we are interested in **flow** out of the pipe). It is connected at the **inlet** of the block with the name **RCR**. Values of **time** (t) are set to the beginning and end of a cardiac cycle (0.0 to 1.0 s) and the corresponding **flow values** (Q) are set to 1.0, as this flow will be received from the 3D simulation.

The RCR boundary condition block sets up the RCR element with the desired resistance and pressure values.

```
{
    "simulation_parameters": {
        "coupled_simulation": true,
        "number_of_time_pts": 100,
        "output_all_cycles": true,
        "steady_initial": false
    },
    "boundary_conditions": [
        {
            "bc_name": "RCR",
            "bc_type": "RCR",
            "bc_values": {
                "Rp": 121.0,
                "Rd": 1212.0,
                "C": 1.5e-4,
                "Pd": 0.0
            }
        }
    ],
    "external_solver_coupling_blocks": [
        {
            "name": "RCR_coupling",
            "type": "FLOW",
            "location": "inlet",
            "connected_block": "RCR",
            "periodic": false,
            "values": {
                "t": [0.0, 1.0],
                "Q": [1.0, 1.0]
            }
        }
    ],
    "junctions": [],
    "vessels": []
}
```

### svZeroD_interface.dat

This file sets up the interface between svMultiPhysics and sv0DSolver. It requires the path of the dynamic library for svZeroDSolver and the input file (svzerod_3Dcoupling.json) discussed above.

This file also matches the external coupling blocks in the 0D model to the coupled surfaces in svMultiPhysics:
The first element in each line should be the name of the block from the json file and the second element should be the index of the coupled surface in svMultiPhysics. In this case, there is only one coupled surface with index 0.

```
svZeroD external coupling block names to surface IDs (where surface IDs are from *.svpre file):
RCR_coupling 0
```

The next lines initialize the pressure and flow of these coupled surfaces in the 0D model:
0 indicates that the values will not be initialized, and 1 indicates that they will be initialized to the value provided afterwards.

```
Initialize external coupling block flows:
0

External coupling block initial flows (one number is provided, it is applied to all coupling blocks):
0.0

Initialize external coupling block pressures:
1

External coupling block initial pressures (one number is provided, it is applied to all coupling blocks):
0.0
```
