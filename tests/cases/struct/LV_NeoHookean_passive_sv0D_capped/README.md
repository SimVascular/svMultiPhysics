This test case simulates an idealized left ventricle (LV) with a NeoHookean material model
coupled to a lumped-parameter network (LPN), implemented in sv0DSolver. The LPN consists of a large pressure source and large resistor, which together produce an approximately constant flowrate into
the LV. This inflates the LV at an approximately constant rate of change of volume.
This is the same test case as LV_NeoHookean_passive_sv0D except for a virtual cap surface on the endocardial surface.
This is used to correctly calculate the flow rate if the base/top moves.
Dirichlet BC on the base has been removed, and a Robin BC has been added to the epicardial surface.

### Build svZeroDSolver
Importantly, to automatically run test cases with `pytest` (see below), you need to build `svZeroDSolver` in the folder
```
./svZeroDSolver/build
``` 
in the repository root.

To do so, you can run the following in the svFSIplus repository root:
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

The input file [svFSI_genBC.xml](./svFSI.xml) follows the master input file as a template, but with coupling to sv0DSolver specified in the options:

```
   <Couple_to_svZeroD type="SI">
   </Couple_to_svZeroD>
```

This tells the solver that the 0d models will be calculated through sv0DSolver. Options to couple 0D codes with svFSI are `N`: none; `I`: implicit; `SI`: semi-implicit; `E`: explicit.

```
   <Add_BC name="top" > 
      <Type> Dirichlet </Type> 
      <Value> 0.0 </Value>
   </Add_BC> 

   <Add_BC name="endo" > 
      <Type> Neu </Type> 
      <Time_dependence> Coupled </Time_dependence> 
      <Follower_pressure_load> true </Follower_pressure_load> 
   </Add_BC> 
```

In this example, we use the LPN for the flow into the endocardial surface.

### svzerod_3Dcoupling.json

This is the configuration file for sv0DSolver and contains the elements of the 0D model being coupled to the 3D simulation. 

For more information on the available parameters and elements, documentation is available here: [svZeroDSolver](https://github.com/SimVascular/svZeroDSolver)

**The following are necessary in "simulation_parameters" for a coupled simulation:**
"coupled_simulation": true,
"steady_initial": false

The boundary condition "P_SOURCE" defines the pressure source, starting with a ramp from 0.0 to 1.0E7 and then a plateau.

The large resistor is defined as a blood vessel block that receives pressure from "P_SOURCE" and has resistance 1.0E5.

The external coupling block "LV_IN" specifies the connection between the 0D model and the coupled surface from svFSIplus. Here, we indicate that it will receive values of flow (type: FLOW) from the outlet (location: outlet) of the resistor block (connected_block: branch0_seg0).

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
            "bc_name": "P_SOURCE",
            "bc_type": "PRESSURE",
            "bc_values": {
                "P": [0.0, 1.0E7, 1.0E7],
                "t": [0.0, 0.1, 1.0]
            }
        }
    ],
    "external_solver_coupling_blocks": [
        {
            "name": "LV_IN",
            "type": "FLOW",
            "location": "outlet",
            "connected_block": "branch0_seg0",
            "periodic": false,
            "values": {
                "Q": [1.0, 1.0],
                "t": [0.0, 1.0]
            }
        }
    ],
    "junctions": [],
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
}
```

### svZeroD_interface.dat

This file sets up the interface between svFSIplus and sv0DSolver. It requires the path of the dynamic library for svZeroDSolver and the input file (svzerod_3Dcoupling.json) discussed above.

This file also matches the external coupling blocks in the 0D model to the coupled surfaces in svFSIplus:
The first element in each line should be the name of the block from the json file and the second element should be the index of the coupled surface in svFSIplus. In this case, there is only one coupled surface with index 0.

```
svZeroD external coupling block names to surface IDs (where surface IDs are from *.svpre file):
LV_IN 0
```
