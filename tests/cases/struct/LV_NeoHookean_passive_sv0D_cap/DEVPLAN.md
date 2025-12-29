## Development plan
1. Replicate non-capped behavior using object-oriented BC. Use test [LV_NeoHookean_passive_sv0D_cap](../LV_NeoHookean_passive_sv0D) to check this. This will allow to identify all the places where the code needs to be modified.
2. Implement cap using ZeroDBoundaryCondition class.
3. Compare Devin's results (using GenBC) with the object-oriented implementation.
4. (Not a requirement) Refactor/remove the Neumann coupled implementation, such that if using sv0D, the object-oriented class must always be used. 

### Where I am at
* I called the new BC Neu0D (I am up for changing it). When used in the input file, the only thing needed is 
```
   <Add_BC name="endo" > 
      <Type> Neu0D </Type> 
      <Time_dependence> Coupled </Time_dependence> 
      <Follower_pressure_load> true </Follower_pressure_load> 
      <svZeroDSolver_block> LV_IN </svZeroDSolver_block>
      <svZeroDSolver_cap> mesh/mesh-surfaces/endo_cap.vtp </svZeroDSolver_cap>  <!--  (optional) -->
   </Add_BC> 
```
* Created class template and added the lines to read the new BC and call the constructor (the constructor is not doing anything yet).
* Identified all the places in the code that need to be modified. I replicated the behavior of using Neu. 

### Files that need to be modified
Small update, e.g., adding new types/class/etc or if statements to call the right function
* [ComMod.h](../../../../Code/Source/solver/ComMod.h): add new class
* [const.h](../../../../Code/Source/solver/consts.h): add new BC type
* [Parameters.cpp](../../../../Code/Source/solver/Parameters.cpp): add input file option `svZeroDSolver_cap`
* [fils_struct.hpp](../../../../Code/Source/liner_solver/fils_struct.hpp): add new BC type
* [solve.cpp](../../../../Code/Source/liner_solver/solve.cpp): add Neu0D to if statement
* [read_files.cpp](../../../../Code/Source/solver/read_files.cpp): reading options, initializing class, follower pressure

Important updates
* [baf_ini.cpp](../../../../Code/Source/solver/baf_ini.cpp): 
    * `baf_ini` assignment of internal variables
    * `fsi_ls_ini`: compute integral of normal vector over the face
* [set_bc.cpp](../../../../Code/Source/solver/set_bc.cpp): 
   * `calc_der_cpl_bc`: Calculating flowrates ar 3D Neumann BC
   * `set_bc_cpl`: Calculating flowrates ar 3D Neumann BC (initialization?)
   * `set_bc_neu`: Initialize (?) Neumann BC
* [svZeroD_subroutines.cpp](../../../../Code/Source/solver/svZeroD_subroutines.cpp):
   * `get_coupled_QP`: as the name says
   * `init_svZeroD`: initialization
   * `calc_svZeroD`: calc 0D

Note: I left a commented `<<dev_cap>>` marking the places where the code needs to be changed.

### Added files
* [ZeroDBoundaryCondition.cpp](../../../../Code/Source/solver/ZeroDBoundaryCondition.cpp)
* [ZeroDBoundaryCondition.h](../../../../Code/Source/solver/ZeroDBoundaryCondition.h)
