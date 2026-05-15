<h1> svFSI_cinterface </h1>

 This directory contains the C++ solver implementation. The C++ code attempts to replicate the data 
 structures and flow of control of the Fortran implementation. 

 The `Simulation` class stores solver parameters and data, and has methods to read solver parameters 
 from a XML file and create mesh data (e.g. nodes and element connectivity) from VTK .vtu and
 .vtp files. A `Simulation` object can be created and stored in the the Fortran `CMMOD` module. It can then
 be called to replicated the Fortran flow of control and compare data created by C++ and Fortran.
 
 
 <h1> Reading  the solver parameter input XML file </h1>
 
 The `Parameters` class is used to read and store solver parameters from an XML file. The XML file organization and parameter names replicate the old input text file except that parameter names have spaces replaced by underscores.


<h1> Calling C++ from Fortran </h1>

 The `simulation_wrap.h,.cpp` files provide an interface used to call C++ functions from Fortran.
 Data created in the C++ code can be returned to the Fortran code for comparison.
 
 The following interface functions are defiined
 
 - create_simulation() - Creates and returns a Simulation object.

 - get_mesh_data() - Returns mesh data for a specific mesh name.

 - read_files() - Reads in the solver parameter input XML file and creates all mesh data.


<h1> C++ organizatiion to replicate Fortran </h1>

 Most of the Fortran code is replicated in C++ using the same file and subroutine names converted to 
 lower case with underscores. For example
 
```
   ================================================================================================
                Fortran                       |                      C++ 
   ================================================================================================
           SUBROUTINE READFILES               |                  read_files()
   ------------------------------------------------------------------------------------------------
           SUBROUTINE READMSH                 |                  read_msh()
   ------------------------------------------------------------------------------------------------
               LOADMSH.f                      |                  load_msh.cpp
   ------------------------------------------------------------------------------------------------
               VTKXML.f                       |                  vtk_xml.cpp
   ------------------------------------------------------------------------------------------------
```

All Fortan subroutines located in a particular file will typically have a C++ implementation in a similarly named file. 

C++ functions are defined within a `namespace` defined for each Fortran file. For example, the functioins in  `load_msh.cpp` are defined within the `load_msh` `namespace`. Some `namespaces` are named with a `_ns` suffix to prevent conflicts with function names (e.g. `read_files_ns`).

# Fortran Modules 

 C++ classes are used to implement Fortran modules. Fortran variable names are retained to prevent (or maintain) confusion. There are no global variables.
 Modules are accessed from a `Simulation` object.

 ## CmMod 

 The `CmMod` class implements the `CMMOD` module defined in `COMU.f`.  This is used for MPI communication.


 ## ComMod  

 The `ComMod` class implments the `COMOD` module defiined in `MOD.f`.

 All solver data is stored in the `ComMod` class using the same names used in the `COMOD` module.
 
 <h1 id="flow_of_control"> Flow of control </h1>
 
 The following outlines the code flow of control for simulations.

 - [PROGRAM MAIN](#program_main)
   - [sim_interface = simulation()](#function_create_simulation)
   - [call sim_interface % read_files(in_file_name)](#subroutine_read_files)
     - [call read_files_c(this % Simulation_object, c_file_name)](#read_files_c)
       - [<i>read_files_ns::</i><b>read_files(simulation, std::string(file_name))</b>](#read_files_ns_read_files)
         - [simulation-><b>read_parameters()</b>](#simulation_read_parameters) - Read solver parameter XML file
         - [simulation-><b>set_module_parameters()</b>](#simulation_set_module_parameters) - Set ComMod module varliables
         - [<i>read_msh_ns::</i><b>read_msh(simulation)</b>](#read_msh_ns_read_msh) - Read all mesh and BCs data
           - [<i>load_msh::</i><b>read_sv(param, mesh)</b>](#load_msh_read_sv) - Read mesh nodal coordinates and element connectivity
             - [<i>vtk_xml::</i><b>read_vtu(mesh_path, mesh)</b>](#read_vtu) - Read a mesh from a SimVascular .vtu file
               - [<i>vtk_xml_parser::</i><b>load_vtu(file_name, mesh)</b>](#vtk_xml_parser_load_vtk) - Read a mesh from a .vtu file and store its data into mesh
             - [<i>nn::</i><b>select_ele(simulation, mesh)</b>](#nn_select_ele) - Set mesh variables for the input element type
               - [<b>set_3d_element_props\[mesh.eNoN\](mesh)</b>](#set_3d_element_props) - Set element properties based on the number of element nodes
               - [<b>get_gip(simulation, mesh)</b>](#get_gip_mesh) - Set mesh w and xi arrays used for Gauss integration
               - [<b>get_gnn(simulation, g, mesh)</b>](#get_gnn_mesh) - Create mesh N and Nx shape function arrays for each Gaus point g
               - [<b>get_nn_bnds(simulation, mesh)</b>](#get_nn_bnds_mesh) - Create bounds on Gauss integration points and shape functions
             - [<i>read_msh_ns::</i><b>check_ien(simulation, mesh)</b>](#read_msh_ns_check_ien) - Check/change the mesh connectivity and node ordering 
               - [<b>check_element_conn\[eType\](mesh)</b>](#check_element_conn)
             - [<i>vtk_xml::</i><b>read_vtp(face_path, face)</b>](#vtk_xml_read_vtp) - Read a face nodal coordinates and element connectivity from a .vtp file
               - [<i>vtk_xml_parser::</i><b>load_vtp(file_name, face)</b>](#vtk_xml_parser_load_vtp) - Store data read from a .vtp file into a faceType object
             - [<i>nn::</i><b>select_eleb(simulation, mesh, face)</b>](#nn_select_eleb) - Set face properties for the input element type
               - [<b>set_face_element_props\[face.eNoN\](insd, face)</b>](#set_face_element_props)
               - [<b>get_gip(simulation, face)</b>](#get_gip_face)
               - [<b>get_gnn(simulation, g, face)</b>](#get_gnn_face)
           - [<b>set_projector(simulation, avNds)</b>](#set_projector) - Process projection faces 
             - [<b>all_fun::find_face(com_mod.msh, ctmpi, iM, iFa)</b>](#find_face) - Process projection faces 
             - [<b>match_faces(com_mod, face1, face2, tol, lPrj)</b>](#match_faces) - Process projection faces 
               - [<b>find_blk(nsd, nBkd, nFlt, xMin, dx, coord)</b>](#find_blk) - Compute the block ID for the given coordinate 
           - [<b>read_fib_nff(simulation, com_mod.msh\[iM\], cTmp, "FIB_DIR", i)</b>](#read_fib_nff) - Read fiber orientation
             - [<b>vtk_xml_parser::load_fiber_direction_vtu(fName, kwrd, idx, simulation->com_mod.nsd, mesh)</b>](#load_fiber_direction_vtu) - Read fiber orientation from vtu
           - [<b>vtk_xml::read_vtu_pdata(cTmp, "Stress", com_mod.nsd, com_mod.nsymd, 0, com_mod.msh\[iM\])</b>](#read_vtu_pdata) - Read prestress data 
           - [<b>load_var_ini(simulation, com_mod)</b>](#load_var_ini) - Set initial mesh pressure, velocity or displacement from a file 
           - <b>Read contact model parameters (not implmented, no tests)</b>
         - [<b>read_eq(simulation, simulation, eq_params, eq)</b>](#read_eq) - Iterate to read equations
           - [<b>set_equation_properties(simulation, eq_params, lEq, propL, outPuts, nDOP)</b>](#set_equation_properties) - Set equation properties 
             - [<b>set_equation_props[eq_type](simulation, eq_params, lEq, propL, outPuts, nDOP)</b>](#set_equation_props) - Set equation properties (function map)
           - [<b>read_outputs(simulation, eq_params, lEq, nDOP, outPuts)</b>](#read_outputs) - Set output parameters 
           - [<b>read_bc(simulation, Set output parameters)</b>](#read_bc) - Iterate to read boundary conditions 
             - [<b>read_trac_bcff(com_mod, lBc.gm, com_mod.msh\[iM\].fa\[iFa\], file_name)</b>](#read_trac_bcff) - Read traction data 
             - [<b>read_temporal_values_file(file_name, lBc)</b>](#read_temporal_values_file) - Set boundary condition temporal values read in from a file
             - [<b>read_fourier_coeff_values_file(file_name, lBc)</b>](#read_fourier_coeff_values_file) - Set boundary condition Fourier coefficients read in from a file 
             - [<b>read_bct(com_mod, lBc.gm, com_mod.msh\[iM\].fa\[iFa\], file_name)</b>](#read_bct) - Reads general velocity data from bct.vtp 
             - [<b>read_temporal_and_spatial_values_file(com_mod, com_mod.msh\[iM\], com_mod.msh\[iM\].fa\[iFa\], file_name, lBc)</b>](#read_temporal_and_spatial_values_file) - Read in a file containing temporal and spatial values 
             - [<b>read_spatial_values_file(com_mod, com_mod.msh\[iM\], com_mod.msh\[iM\].fa\[iFa\], file_name, lBc)</b>](#read_spatial_values_file) - Read in a file containing spatial values 
             - [<b>vtk_xml::read_vtp_pdata(cTmp, "Stress", com_mod.nsd, com_mod.nsymd, data_series, face)</b>](#read_vtp_pdata) - Read prestress tensor for CMM
             - [<b>vtk_xml::read_vtp_pdata(cTmp, "Displacement", com_mod.nsd, com_mod.nsd, data_series, face)</b>](#read_vtp_pdata) - Read displacement for CMM
   - [<b>call sim_interface % distribute_data()</b>](#subroutine_distribute_data)
     - [<b>call distribute_data_c(this % Simulation_object)</b>](#distribute_data)
       - [<b>distribute(Simulation* simulation))</b>](#distribute) - Partition and distribute data across processors
         - [<b>cm.bcast(cm_mod, &chnl_mod.pClr)</b>](#bcast) - 
         - [<b>cm.bcast(cm_mod, chnl_mod.appPath)</b>](#bcast) - 
         - [<b>cm.bcast(cm_mod, &com_mod.nMsh)</b>](#bcast) - 
         - [<b>cm.bcast(cm_mod, &com_mod.nsd)</b>](#bcast) - 
         - [<b>cm.bcast(cm_mod, wrk)</b>](#bcast) - 
         - [<i>all_fun::</i><b>split_jobs(task_id, nMsh, num_proc, wgt, wrk)</b>](#split_jobs) - Spliting "m" jobs between "n" workers
         - [<b>part_msh(simulation, com_mod.msh\[iM\], gmtl, num_proc, iWgt)</b>](#part_msh) - Partition each mesh iM
           - [<i>nn::</i><b>select_ele(simulation, lM)</b>](#select_ele) - Spliting "m" jobs between "n" workers
           - [<b>split_(&nEl, &eNoN, &eNoNb, lM.IEN.data_, &num_proc, lM.eDist.data_, wgt.data_, part.data_)</b>](#split_) - Partitioning using ParMetis 
         - [<b>part_face(simulation, msh, face, tMs\[iM\].fa\[iFa\], gmtl)</b>](#part_face) - Partitioning iFa faces for each mesh iM  
           - [<i>nn::</i><b>select_eleb(simulation, lM, gFa)</b>](#select_eleb) - Set face properties for the input element type
         - [<i>all_fun::</i><b>local_rv(com_mod, cm_mod, cm, tmpX)</b>](#local_rv) - Send a real vector to all the processors
         - [<b>dist_eq(com_mod, cm_mod, cm, tMs, gmtl, cep_mod, com_mod.eq\[iEq\])</b>](#dist_eq) - Distribute equation iEq to processors 
           - [<b>dist_mat_consts(com_mod, cm_mod, cm, dmn.stM)</b>](#dist_mat_consts) - Distribute material properties to processors 
           - [<b>dist_bc(com_mod, cm_mod, cm, lEq.bc\[iBc\], tMs, gmtl)</b>](#dist_bc) - Distribute boundary condition data to processors 
           - [<b>dist_visc_model(com_mod, cm_mod, cm, dmn.visc)</b>](#dist_visc_model) - Distribute viscosity model to processors 
           - [<b>dist_bf(com_mod, cm_mod, cm, lEq.bf\[iBf\])</b>](#dist_bf) - Distribute boundary fource to processors 
   - [<b>call sim_interface % initialize_solution(timeP)</b>](#subroutine_initialize_solution)
     - [<b>call initialize_solution_c(this % Simulation_object, timeP)</b>](#initialize_solution_c)
       - [<b>initialize_solution(Simulation* simulation)</b>](#initialize_solution) - Initialize or finalize svFSI variables/structures
         - [<b>mat_fun::ten_init(nsd)</b>](#ten_init) - Initialize tensor operations
         - [<b>lhsa_ns::lhsa(simulation, nnz)</b>](#lhsa) - Constructing stiffness matrix
         - [<b>fsi_linear_solver::fsils_commu_create(communicator, cm.com())</b>](#fsils_commu_create) -
         - [<b>fsi_linear_solver::fsils_lhs_create(com_mod.lhs, communicator, com_mod.gtnNo, com_mod.tnNo, nnz,com_mod.ltg, com_mod.rowPtr, com_mod.colPtr, nFacesLS)</b>](#fsils_lhs_create) - Initialize FSILS structures
         - [<b>cep_ion::cep_init(simulation)</b>](#cep_init) - 
           - [<b>IonicModel::init(Vector<double>& X, Vector<double>& Xg) const</b>](#ionic_model_init) -
         - [<b>fs::init_fs_msh(com_mod, mesh)</b>](#init_fs_msh) -
         - [<b>fs::init_fs_face(com_mod, mesh, mesh.fa[iFa])</b>](#init_fs_face) -
         - [<b>all_fun::integ(com_mod, cm_mod, i, s, 0, 0)</b>](#integ) - Calculating the volume of each domain
         - [<b>baf_ini_ns::baf_ini(simulation)</b>](#baf_ini) - Preparing faces
           - [<b>set_bc::rcr_init(com_mod, cm_mod)</b>](#)
           - [<b>set_bc::genBC_Integ_X(com_mod, cm_mod, "I")</b>](#)
           - [<b>set_bc::calc_der_cpl_bc(com_mod, cm_mod)</b>](#)
           - [<b>fsi_ls_ini(com_mod, cm_mod, bc, com_mod.msh[iM].fa[iFa], lsPtr)</b>](#fsi_ls_ini) 
           - [<b>fsils_bc_create(com_mod.lhs, lsPtr, i, nsd, BcType::BC_TYPE_Dir, gNodes)</b>](#fsils_bc_create)
         - [<b>set_bc::set_bc_dir(com_mod, com_mod.Ao, com_mod.Yo, com_mod.Do)</b>](#set_bc_dir) - Preparing BCs
           - [<b>set_bc::set_bc_dir_l(com_mod, bc, com_mod.msh[iM].fa[iFa], tmpA, tmpY, lDof)</b>](#set_bc_dir_l)
         - [<b>txt_ns::txt(simulation, true)</b>]() -
   - [<b>call sim_interface % iterate_solution()</b>](#subroutine_iterate_solution)
     - [<b>call interate_solution_c(this % Simulation_object)</b>](#iterate_solution)
       - [<b> pic::picp(simulation) </b>](#picp) - Predictor
         - [<b> cep_ion::cep_integ(simulation, iEq, e, Do) </b>](#cep_integ)
       - [<b> set_bc::set_bc_dir(com_mod, An, Yn, Dn) </b>](#set_bc_dir) - Apply Dirichlet BCs strongly
       - <b> Inner Loop </b>
         - [<b> set_bc::set_bc_cpl(com_mod, cm_mod) </b>](#set_bc_cpl) - If com_mod.cplBC.coupled
         - [<b> set_bc::set_bc_dir(com_mod, An, Yn, Dn) </b>](#set_bc_dir) - If com_mod.cplBC.coupled 
         - [<b> pic::pici(simulation, Ag, Yg, Dg) </b>](#pici) - Initiator step
         - [<b> ls_ns::ls_alloc(com_mod, eq) </b>](#ls_alloc)
         - [<b> bf::set_bf(com_mod, Dg) </b>](#set_bf)
           - [<b> bf::set_bf_l(com_mod, eq.bf[iBf], com_mod.msh[iM], Dg) </b>](#set_bf_l)
             - [<b> ifft(com_mod, lBf.bt, f, rtmp) </b>](#) - If bfType_ustd
             - [<b> igbc(com_mod, lBf.bm, bfl, xl) </b>](#) - If bfType_gen
             - [<b> bf::bf_construct(com_mod, lM, e, eNoN, idof, xl, dl, bfl, ptr) </b>](#bf_construct) - For shell follower pressre loads or init  pressure for CMM
         - [<b> eq_assem::global_eq_assem(com_mod, com_mod.msh[iM], Ag, Yg, Dg) </b>](#global_eq_assem) - For each mesh iM
           - [<b> fluid::construct_fluid(com_mod, lM, Ag, Yg) </b>](#construct_fluid) -  If EquationType::phys_fluid
             - [<b> fs::get_thood_fs(com_mod, fs, lM, vmsStab, 1) </b>](#)
             - [<b> nn::gnn(fs[1].eNoN, nsd, nsd, Nx, xql, Nqx, Jac, ksix) </b>](#)
             - [<b> nn::gn_nxx(l, fs[0].eNoN, nsd, nsd, Nx, Nxx, xwl, Nwx, Nwxx) </b>](#)
             - [<b> fluid_3d_m(com_mod, vmsStab, fs[0].eNoN, fs[1].eNoN, w, ksix, N0, N1, Nwx, Nqx, Nwxx, al, yl, bfl, lR, lK, K_inverse_darcy_permeability)</b>](#) - If nsd=3
             - [<b> fluid_2d_m(com_mod, vmsStab, fs[0].eNoN, fs[1].eNoN, w, ksix, N0, N1, Nwx, Nqx, Nwxx, al, yl, bfl, lR, lK, K_inverse_darcy_permeability)</b>](#) - If nsd=2
             - [<b> trilinos_doassem_(const_cast<int&>(eNoN), ptr.data(), lK.data(), lR.data())</b>](#) - If using Trilinos
             - [<b> lhsa_ns::do_assem(com_mod, eNoN, ptr, lK, lR)</b>](#do_assem) - If not using Trilinos
         - [<b> set_bc::set_bc_neu(com_mod, cm_mod, Yg, Dg) </b>](#set_bc_neu)
         - [<b> set_bc::set_bc_cmm(com_mod, cm_mod, Ag, Dg) </b>](#set_bc_cmm)
         - [<b> set_bc::set_bc_dir_w(com_mod, Yg, Dg) </b>](#)
         - [<b> contact::contact_forces(com_mod, cm_mod, Dg) </b>](#)
         - [<b> all_fun::commu(com_mod, com_mod.R) </b>](#) - Synchronize residual R across processes
         - [<b> ustruct::ustruct_r(com_mod, Yg) </b>](#)
         - [<b> fs::thood_val_rc(com_mod) </b>](#)
         - [<b> set_bc::set_bc_undef_neu(com_mod) </b>](#)
         - [<b> ls_ns::ls_solve(com_mod, eq, incL, res) </b>](#)
         - [<b> pic::picc(simulation) </b>](#)
         
      

 

<!-- ============================================================================================ -->
<!--                                      Fortran Implementation                                  -->
<!-- ============================================================================================ -->
 
 <h1 id="fortran_implementation_details"> Fortran Implementation </h1>
 
 The following sections provide some implementation details of the Fortran/C++ interface.
 
 <h2 id="program_main"> PROGRAM MAIN </h2>
 
  Fortran main program defined in `MAIN.f`.
  
  - `sim_interface = simulation()` - Create a C++ `Simulation` object, calls `create_simulation()` function.
  - `call sim_interface % read_files(in_file_name)` - Read solver parameters from a solver input XML file and create all of the mesh data. 

<!-- =========================== -->
<!-- create_simulation (Fortran) -->
<!-- =========================== -->

<h2 id="function_create_simulation"> function create_simulation() </h2>

[simulation_interface_mod.f90](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI/simulation_interface_mod.f90)

Fortran / C++ interface function to create a `Simulation` object.

- `create_simulation % Simulation_object = create_simulation_c()` - Calls C++ `create_simulation` functon.


<!-- ======================= -->
<!-- create_simulation (c++) -->
<!-- ======================= -->

<h2 id="create_simulation"> Simulation* create_simulation() </h2>

[simulation_wrap.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/simulation_wrap.cpp)

C++ / Fortran interface function used to create a C++ `Simulation` object and return it to Fortran.

- `Simulation::Simulation()`- Create a C++ Simulation object
  - Create `ComMod` member, set defaults in ctor
  - Create `Parameters` member

<!-- ========== -->
<!-- read_files -->
<!-- ========== -->

<h2 id="subroutine_read_files"> subroutine read_files(this, file_name) </h2>

[simulation_interface_mod.f90](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI/simulation_interface_mod.f90)

Fortran / C++ interface function used to call the C++ `read_files()` function.

- `call read_files_c(this % Simulation_object, c_file_name)` - Call C++ `read_files` functon.


<!-- ========================== -->
<!--  distribute_data (Fortran) -->
<!-- ========================== -->

<h2 id="subroutine_distribute_data"> subroutine distribute_data(this) </h2>

[simulation_interface_mod.f90](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI/simulation_interface_mod.f90)

Fortran / C++ interface function used to call the C++ `()` function.

- `call distribute_data_c(this % Simulation_object)`

<!-- ============================== -->
<!--  initialize_solution (Fortran) -->
<!-- ============================== -->

<h2 id="subroutine_initialize_solution"> subroutine initialize_solution(this, timep) </h2>

[simulation_interface_mod.f90](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI/simulation_interface_mod.f90)

Fortran / C++ interface function used to call the C++ `()` function.

- `call initialize_solution_c(this % Simulation_object, timep_cptr)`
 
<!-- ============================== -->
<!--  iterate_solution (Fortran)    -->
<!-- ============================== -->

<h2 id="subroutine_iterate_solution"> subroutine iterate_solution(this) </h2>

[simulation_interface_mod.f90](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI/simulation_interface_mod.f90)

- `call iterate_solution_c(this % Simulation_object)`



<!-- ============================================================================================ -->
<!--                                      C++ Implementation                                      -->
<!-- ============================================================================================ -->

<h1 id="cpp_implementation_details"> C++ Implementation </h1>

The following sections provide some implementation details of the C++ code replicating Fortran functionality. 


<!-- ========== -->
<!-- read_files -->
<!-- ========== -->

<h2 id="read_files_c">  void read_files(Simulation* simulation, const char* file_name) </h2>

[simulation_wrap.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/simulation_wrap.cpp)

C++ / Fortran interface function used to read in an XML file and all mesh and BC data.

- `read_files_ns::read_files(simulation, std::string(file_name))` 


<h2 id="read_files_ns_read_files"> void read_files_ns::read_files(Simulation* simulation, const std::string& file_name) </h2>

[read_files.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/read_files.cpp)

 Read in an XML file and, all mesh and BC data. Replicates `SUBROUTINE READFILES` in `READFILES.f`
 
 - `Simulation::read_parameters()` - Read solver parameter XML file
 - `Simulation::set_parameters()` - Set module (e.g.`ComMod`) and `Simulation` member data from `Parameters` data
 - `read_msh_ns::read_msh(simulation)` - Read mesh and BCs data

<!-- =============== -->
<!-- read_parameters -->
<!-- =============== -->

<h2 id="simulation_read_parameters"> void Simulation::read_parameters(const std::string& file_name) </h2>

[Simulation.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/Simulation.cpp)

Read solver parameter XML file. 

- `parameters_.read_parameters(file_name)`


<!-- ===================== -->
<!-- set_module_parameters -->
<!-- ===================== -->

<h2 id="simulation_set_module_parameters"> void Simulation::set_module_parameters() </h2>

[Simulation.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/Simulation.cpp)

Set `ComMod` module varliables.

<!-- ======== -->
<!-- read_msh -->
<!-- ======== -->

<h2 id="read_msh_ns_read_msh"> void read_msh_ns::read_msh() </h2>

[read_msh.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/read_msh.cpp)

Read all mesh and BCs data. Replicates `SUBROUTINE READMSH` in `READMSH.f`

- Set the number of meshes`com_mod.nMsh`
- Allocate `com_mod.msh` storing a list of `mshType` objects
- Iterate over each `mesh` in `com_mod.msh` defined by each `Add_mesh` parameter in the XML file with parameter values in `param`
  - `load_msh::read_sv(param, mesh)` - Read in a volume meshes and face meshes from a VTK files
  - `read_msh_ns::check_ien(simulation, mesh)` - Check the mesh element node ordering and change element node ordering if needed
- Re-arranging x and finding the size of the entire domain
- Renumber face node IDs 
- `read_fib_nff(simulation, com_mod.msh[iM], cTmp, "FIB_DIR", i)` - Read fiber orientation
  - `vtk_xml_parser::load_fiber_direction_vtu(fName, kwrd, idx, simulation->com_mod.nsd, mesh)`
- Iterate over each `mesh` iM in `com_mod.msh` to set read prestress data
  - vtk_xml::read_vtu_pdata(cTmp, "Stress", com_mod.nsd, com_mod.nsymd, 0, com_mod.msh[iM]) - Read prestress data
- If have prestress data then set com_mod.pS0() 
- Set initial mesh pressure, velocity or displacement from a file
  - `load_var_ini(simulation, com_mod)`


<!-- ======= -->
<!-- read_sv -->
<!-- ======= -->

<h2 id="load_msh_read_sv"> void load_msh::read_sv(Simulation* simulation, mshType& mesh, const MeshParameters& param) </h2> 

[load_mesh.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/load_msh.cpp)

Read mesh nodal coordinates and element connectivity. Replicates `SUBROUTINE READSV(list, lM)` in `LOADMSH.f`.

- `vtk_xml::read_vtu(mesh_path, mesh)` -  Read in volume mesh
- `nn::select_ele(simulation, mesh)` - Set mesh element properites for the input element type


<!-- ======== -->
<!-- read_vtu -->
<!-- ======== -->

<h2 id="read_vtu"> void read_vtu(const std::string& file_name, mshType& mesh) </h2>

[vtk_xml.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/vtk_xml.cpp)

Read a mesh from a SimVascular `.vtu` file. Replicates Fortran `SUBROUTINE READVTU(lM, fName)` defined in `VTKXML.f`.

- `vtk_xml_parser::load_vtk(VtkFileFormat::VTU, file_name, mesh)` - Read a VTK `.vtu` file and store its data into `mesh`.

<!-- ======== -->
<!-- load_vtu -->
<!-- ======== -->

<h2 id="vtk_xml_parser_load_vtk"> void vtk_xml_parser::load_vtu(const std::string& file_name, mshType& mesh) </h2>

[vtk_xml_parser.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/vtk_xml_parser.cpp)

Read a mesh from a `.vtu` file and store its data into `mesh`. This calls VTK functions.

Replicates 
```
subroutine loadVTK(vtk,fName,istat)
subroutine getVTK_numPoints(vtk,nn,istat)
subroutine getVTK_numElems(vtk,ne,istat)
subroutine getVTK_nodesPerElem(vtk,eNoN,istat)
subroutine getVTK_pointCoords(vtk,x,istat)
subroutine getVTK_elemIEN(vtk,ien,istat)
```
defined in `vtkXMLParser.f90`.

The following `mesh` variables are set
```
 mesh.gnNo - number of nodes
 mesh.x - node coordinates
 mesh.gN - node IDs 
 mesh.gnEl - number of elements
 mesh.eNoN - number of noders per element
 mesh.gIEN - element connectivity (num_nodes_per_elem, num_elems)
```

<!-- =========== -->
<!-- select_ele  -->
<!-- =========== -->

<h2 id="nn_select_ele"> nn::select_ele(Simulation* simulation, mshType& mesh) </h2>

[nn.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/nn.cpp)

Set mesh variables for the input element type.

- `set_3d_element_props[mesh.eNoN](mesh)` - Set element properties based on the number of element nodes
- `get_gip(simulation, mesh)` - Set mesh `w` and `xi` arrays used for Gauss integration
- `get_gnn(simulation, g, mesh)` - Create mesh `N` and `Nx` shape function arrays for each Gaus point `g`
- `get_nn_bnds(simulation, mesh)` - Create bounds on Gauss integration points and shape functions


<!-- ==================== -->
<!-- set_3d_element_props -->
<!-- ==================== -->

<h2 id="set_3d_element_props"> std::map&lt;int, std::function&lt;void(int, mshType&)&gt;&gt; set_3d_element_props\[\] </h2>

[nn_elem_props.h](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/nn_elem_props.h)

A map used to set 3D element properties based on the number of element nodes. This replicates the case statement in the Fortran 'SUBROUTINE SELECTELE(lM)' defined in NN.f. 


<!-- ======= -->
<!-- get_gip -->
<!-- ======= -->

<h2 id="get_gip_mesh"> nn::get_gip(Simulation* simulation, mshType& mesh) </h2>

[nn.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/nn.cpp)

Set mesh `w` and `xi` arrays used for Gauss integration.

- `set_element_gauss_int_data[mesh.eType](mesh)`


<!-- ======= -->
<!-- get_gnn -->
<!-- ======= -->

<h2 id="get_gnn_mesh"> nn:get_gnn(Simulation* simulation, int gaus_pt, mshType& mesh) </h2>

[nn.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/nn.cpp)

Create mesh `N` and `Nx` shape function arrays for a Gauss point `gaus_pt`.

- set_element_shape_data[mesh.eType](g, mesh);

<!-- =========== -->
<!-- get_nn_bnds -->
<!-- =========== -->

<h2 id="get_nn_bnds_mesh"> nn::get_nn_bnds(Simulation* simulation, mshType& mesh) </h2>

[nn.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/nn.cpp)

Create bounds on Gauss integration points and shape functions. Replicates Fortran `SUBROUTINE GETNNBNDSlM%eType, lM%eNoN, lM%xib, lM%Nb)`.


<!-- ========= -->
<!-- check_ien -->
<!-- ========= -->

<h2 id="read_msh_ns_check_ien"> read_msh_ns::check_ien(Simulation* simulation, mshType& mesh) </h2>

[read_msh.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/read_msh.cpp)

Check the mesh connectivity and node ordering. It may reorder element connectivity. Replicates the Fortran `CHECKIEN` subroutine defined in `READMSH.f`.

- `check_element_conn[eType](mesh)`


<!-- ================== -->
<!-- check_element_conn -->
<!-- ================== -->

<h2 id="check_element_conn"> std::map&lt;consts::ElementType, std::function&lt;void(mshType&)&gt;&gt; read_msh_ns::check_element_conn </h2>

[read_msh.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/read_msh.cpp)

A map of function pointers used to check element connecivity.

- check_tet_conn((mshType& mesh)
- check_wedge_conn(mshType& mesh)


<!-- ======== -->
<!-- read_vtp -->
<!-- ======== -->

<h2 id="vtk_xml_read_vtp"> vtk_xml::read_vtp(const std::string& file_name, faceType& face) </h2> 

[vtk_xml.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/vtk_xml.cpp)

Read a face nodal coordinates and element connectivity from a `SimVascular` `.vtp` file. Sets data in `faceType` face: 
- face.eNoN - number of noders per element
- face.gebc - EBC array (gE + gIEN) 
- face.gnEl - globel number of elements
- face.nEl - number of elements
- face.nNo - number of nodes
- face.x - node coordinates

Replicates Fortran READVTP subroutine defined in VTKXML.f.

- `vtk_xml_parser::load_vtp(file_name, face)`


<!-- ======== -->
<!-- load_vtp -->
<!-- ======== -->

<h2 id="vtk_xml_parser_load_vtp"> vtk_xml_parser::load_vtp(const std::string& file_name, faceType& face) </h2>

[vtk_xml_parser.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/vtk_xml_parser.cpp)

Store a surface mesh read in from a VTK `.vtp` file into a faceType object.

- auto vtk_polydata = reader->GetOutput()
- auto points = vtk_polydata->GetPoints()
- store_nodal_coords(points, face)
- `store_nodal_ids(vtk_polydata, face)`
- `store_element_conn(vtk_polydata, face)`
- `store_element_ids(vtk_polydata, face)`


<!-- =========== -->
<!-- select_eleb -->
<!-- =========== -->

<h2 id="nn_select_eleb"> nn::select_eleb(Simulation* simulation, mshType& mesh, faceType& face)) </h2>

[nn.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/nn.cpp)

Set face properties for the input element type. 

- `get_gip(simulation, face)`
- `get_gnn(simulation, g, face)`


<!-- ====================== -->
<!-- set_face_element_props -->
<!-- ====================== -->

<h2 id="set_face_element_props"> std::map&lt;int, std::function&lt;void(int, mshType&)&gt;&gt; set_face_element_props </h2>

[nn_elem_props.h](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/nn_elem_props.h)

A map type used to set element properties.


<!-- ======= -->
<!-- get_gip -->
<!-- ======= -->

<h2 id="get_gip_face"> get_gip(Simulation* simulation, faceType& face) </h2>

[nn.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/nn.cpp)

- `set_face_gauss_int_data[face.eType](face)`


<!-- ============ -->
<!-- get_gnn_face -->
<!-- ============ -->

<h2 id="get_gnn_face"> get_gnn(Simulation* simulation, int gaus_pt, faceType& face) </h2>

[nn.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/nn.cpp)

Computes shape functions and derivatives at given natural coords.

- `set_face_shape_data[face.eType](gaus_pt, face)`


<!-- ============= -->
<!-- set_projector -->
<!-- ============= -->

<h2 id="set_projector"> set_projector(simulation, avNds) </h2>

[read_msh.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/read_msh.cpp)

Associates two faces with each other and sets gN, matches the nodal coordinates for each projection face.

- `all_fun::find_face(com_mod.msh, ctmpi, iM, iFa)`

- `match_faces(com_mod, face1, face2, tol, lPrj)`


<!-- ========= -->
<!-- find_face -->
<!-- ========= -->

<h2 id="find_face"> find_face(com_mod.msh, ctmpi, iM, iFa) </h2>

[all_fun.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/all_fun.cpp)


<!-- =========== -->
<!-- match_faces -->
<!-- =========== -->

<h2 id="match_faces"> match_faces(const ComMod& com_mod, const faceType& lFa, const faceType& pFa, const double ptol, utils::stackType& lPrj) </h2>

[read_msh.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/read_msh.cpp)

Match isoparameteric faces to each other.

- `find_blk(nsd, nBkd, nFlt, xMin, dx, coord)`


<!-- ======== -->
<!-- find_blk -->
<!-- ======== -->

<h2 id="find_blk"> int find_blk(const int nsd, const int nBkd, const std::vector<bool>& nFlt, const Vector<double>&xMin, const Vector<double>&dx, const Vector<double>& x)
</h2>

[read_msh.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/read_msh.cpp)

Compute the block ID for the given coordinate.


<!-- ============ -->
<!-- read_fib_nff -->
<!-- ============ -->

<h2 id="read_fib_nff"> read_fib_nff(Simulation* simulation, mshType& mesh, const std::string& fName, const std::string& kwrd, const int idx) </h2>

[read_msh.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/read_msh.cpp)

Read fiber direction from a vtu file. 

- `vtk_xml_parser::load_fiber_direction_vtu(fName, kwrd, idx, simulation->com_mod.nsd, mesh)`


<!-- ======================== -->
<!-- load_fiber_direction_vtu -->
<!-- ======================== -->

<h2 id="load_fiber_direction_vtu"> load_fiber_direction_vtu(const std::string& file_name, const std::string& data_name, const int idx,
    const int nsd, mshType& mesh" </h2>

[vtk_xml_parser.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/vtk_xml_parser.cpp)

Read fiber direction data from a VTK VTU file and copy it into a mesh. 


<!-- ============== -->
<!-- read_vtu_pdata -->
<!-- ============== -->

<h2 id="read_vtu_pdata">  read_vtu_pdata(const std::string& fName, const std::string& kwrd, const int nsd, const int m, const int idx, mshType& mesh) </h2>

[vtk_xml.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/vtk_xml.cpp)

Read prestress data from a vtu file. 

```
  auto vtk_data = VtkData::create_reader(fName);
  int num_elems = vtk_data->num_elems();
  int num_points = vtk_data->num_points();
```

<!-- ============ -->
<!-- load_var_ini -->
<!-- ============ -->

<h2 id="load_var_ini">  load_var_ini(Simulation* simulation, ComMod& com_mod) </h2>

[read_msh.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/read_msh.cpp)

Read initial field values (pressure, velocity or displacement). 


<!-- ======= -->
<!-- read_eq -->
<!-- ======= -->

<h2 id="read_eq">  read_eq(Simulation* simulation, EquationParameters& eq_params, eqType& lEq) </h2>

[read_files.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/read_files.cpp)

Set equation parameters.

- `set_equation_properties(simulation, eq_params, lEq, propL, outPuts, nDOP)`


<!-- ======================= -->
<!-- set_equation_properties -->
<!-- ======================= -->

<h2 id="set_equation_properties"> set_equation_properties(Simulation* simulation, EquationParameters& eq_params, eqType& lEq, EquationProps& propL, 
EquationOutputs& outPuts, EquationNdop& nDOP) </h2>

[read_files.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/read_files.cpp)

Set equation properties.

- eq_type = equation_name_to_type.at(eq_type_str)

- `set_equation_props[eq_type](simulation, eq_params, lEq, propL, outPuts, nDOP)` - Execute function based on 'eq_type'


<!-- ================== -->
<!-- set_equation_props -->
<!-- ================== -->

<h2 id="set_equation_props"> std::map<consts::EquationType, std::function&lt;void(Simulation*, EquationParameters&, eqType&, EquationProps&, 
EquationOutputs&, EquationNdop&)&gt;&gt;) </h2>

[set_equation_props.h](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/set_equation_props.h)

A map of lambda functions used to set equation properties.


<!-- ============ -->
<!-- read_outputs -->
<!-- ============ -->

<h2 id="read_outputs"> read_outputs(Simulation* simulation, EquationParameters& eq_params, eqType& lEq, EquationNdop& nDOP,  EquationOutputs& outPuts) </h2>

[read_files.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/read_files.cpp)

Set output parameters.


<!-- ======= -->
<!-- read_bc -->
<!-- ======= -->

<h2 id="read_bc">  read_bc(Simulation* simulation, EquationParameters& eq_params, eqType& lEq, BoundaryConditionParameters& bc_params, bcType& lBc) </h2>

[read_files.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/read_files.cpp)

Read boundary condition data.

- `read_trac_bcff(com_mod, lBc.gm, com_mod.msh[iM].fa[iFa], file_name)`
- `read_temporal_values_file(file_name, lBc)` 
- `read_fourier_coeff_values_file(file_name, lBc)`
- `read_bct(com_mod, lBc.gm, com_mod.msh[iM].fa[iFa], file_name)`
- `read_temporal_and_spatial_values_file(com_mod, com_mod.msh[iM], com_mod.msh[iM].fa[iFa], file_name, lBc)`
- `read_spatial_values_file(com_mod, com_mod.msh[iM], com_mod.msh[iM].fa[iFa], file_name, lBc)`
- `vtk_xml::read_vtp_pdata(cTmp, "Stress", com_mod.nsd, com_mod.nsymd, data_series, face)`
- `vtk_xml::read_vtp_pdata(cTmp, "Displacement", com_mod.nsd, com_mod.nsd, data_series, face)`


<!-- ============== -->
<!-- read_trac_bcff -->
<!-- ============== -->

<h2 id="read_trac_bcff"> read_trac_bcff(ComMod& com_mod, MBType& lMB, faceType& lFa, const std::string& fName) </h2>

[read_files.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/read_files.cpp)

Reads pressure/traction data from a vtp file and stores in moving BC data structure.

```
VtkVtpData vtp_data(fName)
int num_points = vtp_data.num_points();
```

<!-- ========================= -->
<!-- read_temporal_values_file -->
<!-- ========================= -->

<h2 id="read_temporal_values_file"> read_temporal_values_file(const std::string& file_name, bcType& lBc)  </h2>

[read_files.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/read_files.cpp)

Set boundary condition temporal values read in from a file. 


<!-- ============================== -->
<!-- read_fourier_coeff_values_file -->
<!-- ============================== -->

<h2 id="read_fourier_coeff_values_file"> read_fourier_coeff_values_file(const std::string& file_name, bcType& lBc) </h2>

[read_files.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/read_files.cpp)

Set boundary condition Fourier coefficients read in from a file.

<!-- ======= -->
<!-- read_bct-->
<!-- ======= -->

<h2 id="read_bct"> read_bct(ComMod& com_mod, MBType& lMB, faceType& lFa, const std::string& fName) </h2>

[read_files.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/read_files.cpp)

Reads general velocity data from bct.vtp.

```
  VtkVtpData vtp_data(fName);
  int num_points = vtp_data.num_points();
  if (num_points == 0) {
```

<!-- ======================== -->
<!-- read_spatial_values_file -->
<!-- ======================== -->

<h2 id="read_spatial_values_file"> read_spatial_values_file(const ComMod& com_mod, const mshType& msh, const faceType& lFa,
    const std::string& file_name, bcType& lBc) </h2>

[read_files.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/read_files.cpp)

Read in a file containing spatial values used for a boundary condition.

<!-- ====================== -->
<!--  distribute_data (c++) -->
<!-- ====================== -->

<h2 id="distribute_data"> distribute_data(Simulation* simulation) </h2>

[simulation_wrap.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/simulation_wrap.cpp)

C++ / Fortran interface function used to distribute data to MPI processes. 

- `distribute(simulation)`

<!-- ========== -->
<!-- distribute -->
<!-- ========== -->

<h2 id="distribute"> distribute(Simulation* simulation)  </h2>

[distribute.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/distribute.cpp)

Partition and distribute data across processors.

- `all_fun::split_jobs(task_id, nMsh, num_proc, wgt, wrk)`

- `part_msh(simulation, com_mod.msh[iM], gmtl, num_proc, iWgt)`

- `part_face(simulation, msh, face, tMs[iM].fa[iFa], gmtl)`


<!-- ===== -->
<!-- bcast -->
<!-- ===== -->

<h2 id="bcast"> cmType::bcast(const CmMod& cm_mod, ...) </h2>

[CmMod.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/CmMod.cpp)

Interfacee to `MPI_Bcast()`.


<!-- ========== -->
<!-- split_jobs -->
<!-- ========== -->

<h2 id="split_jobs"> split_jobs(int tid, int m, int n, Array<double>& A, Vector<double>& b) </h2>

[all_fun.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/all_fun.cpp)

Spliting "m" jobs between "n" workers.


<!-- ======== -->
<!-- part_msh -->
<!-- ======== -->

<h2 id="part_msh"> part_msh(Simulation* simulation, mshType& lM, Vector<int>& gmtl, int nP, Vector<float>& wgt) </h2>

[distribute.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/distribute.cpp)

Partition a mesh amongst N processors.

- Broadcast mesh data.
```
  cm.bcast(cm_mod, &lM.lShpF);
  cm.bcast(cm_mod, &lM.lShl);
  cm.bcast(cm_mod, &lM.lFib);

  cm.bcast(cm_mod, &eType);

  cm.bcast(cm_mod, &lM.eNoN);
  cm.bcast(cm_mod, &lM.nFa);
  cm.bcast(cm_mod, &lM.nFs);
  cm.bcast(cm_mod, &lM.nG);
  cm.bcast(cm_mod, &lM.gnEl);
  cm.bcast(cm_mod, &lM.gnNo);
  cm.bcast(cm_mod, lM.name);
  cm.bcast(cm_mod, &lM.nFn);
  cm.bcast(cm_mod, &lM.scF);
```

- Set face properties for the input element type
```
nn::select_ele(simulation, lM)`
```

- Scattering the lM.gIEN array to all processors.
```
  MPI_Scatterv(lM.gIEN.data_, sCount.data_, disp.data_, cm_mod::mpint, lM.IEN.data_, nEl*eNoN, cm_mod::mpint, cm_mod.master, cm.com())
```

- Doing partitioning, using ParMetis
```
  auto edgecut = split_(&nEl, &eNoN, &eNoNb, lM.IEN.data_, &num_proc, lM.eDist.data_,  wgt.data_, part.data_)`
```

- Gathering the parts inside master
```
  MPI_Gatherv(part.data_, nEl, cm_mod::mpint, gPart.data_, sCount.data_, disp.data_, cm_mod::mpint, cm_mod.master, cm.com())
```

- Communicating eId, if neccessary
```
   MPI_Scatterv(tempIEN.data_, sCount.data_, disp.data_, cm_mod::mpint, lM.IEN.data_, nEl*eNoN, cm_mod::mpint, cm_mod.master, cm.com())
```

- Communicating fN, if neccessary
```
  MPI_Scatterv(tmpFn.data_, sCount.data_, disp.data_, cm_mod::mpreal, lM.fN.data_, nEl*nFn*nsd, cm_mod::mpreal, cm_mod.master, cm.com());
```

- scattering the sorted lM%IEN to all processors.
```
  MPI_Scatterv(tempIEN.data_, sCount.data_, disp.data_, cm_mod::mpint, lM.IEN.data_, nEl*eNoN, cm_mod::mpint, cm_mod.master, cm.com());
```

- Constructing the initial global to local pointer


<!-- ========= -->
<!-- part_face -->
<!-- ========= -->

<h2 id="part_face"> part_face(Simulation* simulation, mshType& lM, faceType& lFa, faceType& gFa, Vector<int>& gmtl) </h2>

[distribute.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/distribute.cpp)

Partition faces across processors.

- Broadcast the number of nodes and elements of to slaves 
```
  cm.bcast(cm_mod, &gFa.d);
  cm.bcast(cm_mod, &gFa.eNoN);
  cm.bcast(cm_mod, &gFa.iM);
  cm.bcast(cm_mod, &gFa.nEl);
  cm.bcast(cm_mod, &gFa.gnEl);
  cm.bcast(cm_mod, &gFa.nNo);
```

- `nn::select_eleb(simulation, lM, gFa)` - Set face properties for the input element type.


<!-- ======= -->
<!-- dist_eq -->
<!-- ======= -->

<h2 id="dist_eq"> dist_eq(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, const std::vector<mshType>& tMs, 
const Vector<int>& gmtl, CepMod& cep_mod, eqType& lEq) </h2>

[distribute.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/distribute.cpp)

Distribute equation data to processors.

- Distribute equation parameters
```
  cm.bcast(cm_mod, &lEq.nOutput);
  cm.bcast(cm_mod, &lEq.coupled);
  cm.bcast(cm_mod, &lEq.maxItr);
  cm.bcast(cm_mod, &lEq.minItr);
  cm.bcast(cm_mod, &lEq.roInf);
  cm.bcast_enum(cm_mod, &lEq.phys);
  ...
```

- Distribute linear solver settings
```
  cm.bcast(cm_mod, &lEq.FSILS.foC);
  cm.bcast_enum(cm_mod, &lEq.FSILS.LS_type);
  cm.bcast(cm_mod, &lEq.FSILS.RI.relTol);
  cm.bcast(cm_mod, &lEq.FSILS.GM.relTol);
  cm.bcast(cm_mod, &lEq.FSILS.CG.relTol);
  ...
```

- Distribute domain properties

- `dist_mat_consts(com_mod, cm_mod, cm, dmn.stM)`


<!-- =============== -->
<!-- dist_mat_consts -->
<!-- =============== -->

<h2 id="dist_mat_consts">  dist_mat_consts(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, stModelType& lStM) </h2>

[distribute.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/distribute.cpp)

Distribute material properties to all processors.


<!-- ======= -->
<!-- dist_bc -->
<!-- ======= -->

<h2 id="dist_bc">  dist_bc(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, bcType& lBc, const std::vector<mshType>& tMs, const Vector<int>& gmtl) </h2>

[distribute.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/distribute.cpp)

Distribute boundary condition data to all processors.


<!-- =============== -->
<!-- dist_visc_model -->
<!-- =============== -->

<h2 id="dist_visc_model">  dist_visc_model(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, fluidViscModelType& lVis) </h2>

[distribute.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/distribute.cpp)

Distribute viscosity model to processors. 


<!-- ======= -->
<!-- dist_bf -->
<!-- ======= -->

<h2 id="dist_bf">  dist_bf(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, bfType& lBf) </h2>

[distribute.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/distribute.cpp)

Distribute boundary fource to processors. 
 
 
<!-- ========================= -->
<!--  initialize_solution(c++) -->
<!-- ========================= -->

<h2 id="initialize_solution"> initialize_solution(Simulation* simulation, double* timep) </h2>

[simulation_wrap.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/simulation_wrap.cpp)

C++ / Fortran interface function used to initialize the solution. 

- `initialize(simulation, timeP)`
 
<!-- =========== -->
<!--  initialize -->
<!-- =========== -->

<h2 id="initialize"> initialize(Simulation* simulation, Vector<double>& timeP) </h2>

[initialize.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/initialize.cpp)

 Initialize or finalize svFSI variables/structures.

 Sets the following for each com_mod.eq[].
 ```
  eq.am
  eq.dof
  eq.pNorm
  eq.af
  eq.beta
  eq.gam
  eq.s
  eq.e
 
  com_mod.Ao.resize(tDof,tnNo);
  com_mod.An.resize(tDof,tnNo);
  com_mod.Yo.resize(tDof,tnNo);
  com_mod.Yn.resize(tDof,tnNo);
  com_mod.Do.resize(tDof,tnNo);
  com_mod.Dn.resize(tDof,tnNo);
  com_mod.Bf.resize(nsd,tnNo);
 
  com_mod.An = com_mod.Ao;
  com_mod.Yn = com_mod.Yo;
  com_mod.Dn = com_mod.Do;

```
 
- `mat_fun::ten_init(nsd)`
 
- `lhsa_ns::lhsa(simulation, nnz)` - Constructing stiffness matrix
 
- `fsi_linear_solver::fsils_lhs_create(com_mod.lhs, communicator, com_mod.gtnNo, com_mod.tnNo, nnz,
      com_mod.ltg, com_mod.rowPtr, com_mod.colPtr, nFacesLS)`
 
- `cep_ion::cep_init(simulation)`
 
- `fs::init_fs_msh(com_mod, mesh)`
 
- `baf_ini_ns::baf_ini(simulation)`
 
- `set_bc::set_bc_dir(com_mod, com_mod.Ao, com_mod.Yo, com_mod.Do)`
 
- `txt_ns::txt(simulation, true)`

<!-- =========== -->
<!--     lsh     -->
<!-- =========== -->

<h2 id="lsh"> void lhsa(Simulation* simulation, int& nnz) </h2>

[lhsa.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/lhsa.cpp)

Create data structure and assembling LHS sparse matrix.

Modifies:
 ```
 com_mod.idMap
 com_mod.colPtr
 com_mod.rowPtr
 ```
 
<!-- ==================== -->
<!--  fsils_commu_create  -->
<!-- ==================== -->

<h2 id="fsils_commu_create">  void fsils_commu_create(FSILS_commuType& commu, cm_mod::MpiCommWorldType commi) </h2>

[commu.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSILS_cinterface/commu.cpp)
 
Modifies:
 ```
  commu.task
  commu.nTasks
  
  commu.foC = true;
  commu.comm = commi;
  commu.nTasks = 1;
  commu.task   = 0;
  commu.master = 0;
 ```

 
<!-- ==================== -->
<!--   fsils_lhs_create   -->
<!-- ==================== -->

<h2 id="fsils_lhs_create"> void fsils_lhs_create(FSILS_lhsType& lhs, FSILS_commuType& commu, int gnNo, int nNo, int nnz, Vector<int>& gNodes, Vector<int> &rowPtr, Vector<int>& colPtr, int nFaces) </h2>

[lhs.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSILS_cinterface/lhs.cpp)

Create data structure and assembling LHS sparse matrix.

Modifies:
```
  lhs.foC 
  lhs.gnNo 
  lhs.nNo 
  lhs.nnz 
  lhs.commu 
  lhs.nFaces 
  lhs.mynNo 

  lhs.colPtr
  lhs.rowPtr
  lhs.diagPtr
  lhs.map
  lhs.face
```
 
<!-- ==================== -->
<!--        cep_init      -->
<!-- ==================== -->

<h2 id="cep_init"> void cep_init(Simulation* simulation) </h2>

[cep_ion.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/cep_ion.cpp)

 Modifies
 ```
 cep_mod.Xion
 ```
 
 - `IonicModel::init(Xl, Xgl)`
 
 - `all_fun::commu(com_mod, sA)`
 - `all_fun::commu(com_mod, sF)`
 
<!-- ==================== -->
<!--     init_fs_msh      -->
<!-- ==================== -->

<h2 id="init_fs_msh"> void init_fs_msh(const ComMod& com_mod, mshType& lM) </h2>

[fs.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/fs.cpp)

 Modifies
 ```
  lM.fs.resize(lM.nFs)
  lM.fs[0].lShpF = lM.lShpF;
  lM.fs[0].eType = lM.eType;
  lM.fs[0].eNoN  = lM.eNoN;
  lM.fs[0].nG    = lM.nG;
 ```
 
 - `nn::get_gn_nxx(insd, ind2, lM.fs[0].eType, lM.fs[0].eNoN, g, lM.fs[0].xi, lM.fs[0].Nxx)`
 
 - `set_thood_fs(lM.fs[1], lM.fs[0].eType)`
 
 - `alloc_fs(lM.fs[1], nsd, insd)`


<!-- ==================== -->
<!--        baf_ini       -->
<!-- ==================== -->

<h2 id="baf_ini"> void baf_ini(Simulation* simulation) </h2>

[baf_ini.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/baf_ini.cpp)

 Modifies
 ```
   com_mod.cplBC.fa
   com_mod.cplBC.xn

   com_mod.cplBC.fa[i].RCR.Rp = bc.RCR.Rp;
   com_mod.cplBC.fa[i].RCR.C  = bc.RCR.C;
   com_mod.cplBC.fa[i].RCR.Rd = bc.RCR.Rd;
   com_mod.cplBC.fa[i].RCR.Pd = bc.RCR.Pd;
   com_mod.cplBC.fa[i].RCR.Xo = bc.RCR.Xo;

 ```
 
 - `face_ini(simulation, msh, face)` - Compute face normals and area
 
 - `shl_ini(com_mod, cm_mod, com_mod.msh[iM])`
 
 - `bc_ini(com_mod, cm_mod, bc, com_mod.msh[iM].fa[iFa])` - Initialize face BC profile
 
 - `shl_bc_ini(com_mod, cm_mod, bc, com_mod.msh[iM].fa[iFa], com_mod.msh[iM])`
 
 - `set_bc::rcr_init(com_mod, cm_mod)`
 
 - `set_bc::genBC_Integ_X(com_mod, cm_mod, "I")`
 
 - `set_bc::calc_der_cpl_bc(com_mod, cm_mod)`
 
 - `fsi_ls_ini(com_mod, cm_mod, bc, com_mod.msh[iM].fa[iFa], lsPtr)`
 
 - `fsils_bc_create(com_mod.lhs, lsPtr, i, nsd, BcType::BC_TYPE_Dir, gNodes)`
 

 
<!-- =================== -->
<!--      fsi_ls_ini     -->
<!-- =================== -->

<h2 id="fsi_ls_ini"> void fsi_ls_ini(ComMod& com_mod, const CmMod& cm_mod, bcType& lBc, const faceType& lFa, int& lsPtr) </h2>

[baf_ini.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/baf_ini.cpp)

 Modifies
 ```
 
 ```
 
 - `fsils_bc_create(com_mod.lhs, lsPtr, lFa.nNo, nsd, BcType::BC_TYPE_Dir, gNodes, sVl)`
 
 
<!-- =================== -->
<!--    fsils_bc_create  -->
<!-- =================== -->

<h2 id="fsils_bc_create"> void fsils_bc_create(FSILS_lhsType& lhs, int faIn, int nNo, int dof, BcType BC_type, const Vector<int>& gNodes,
    Array<double> Val) </h2>

[bc.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSILS_cinterface/bc.cpp)

 Modifies
 ```
 lhs.face[faIn]
 ```
 
 - `MPI_Allreduce(&a, &Ac, 1, cm_mod::mpint, MPI_SUM, lhs.commu.comm)`
 
 - `fsils_commuv(lhs, dof, v)`
 
  
<!-- =================== -->
<!--       set_bc_dir    -->
<!-- =================== -->

<h2 id="set_bc_dir"> void set_bc_dir(ComMod& com_mod, Array<double>& lA, Array<double>& lY, Array<double>& lD) </h2>

[set_bc.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/set_bc.cpp)

 Modifies
 ```
com_mod.Ad
lA(tDof, tnNo)
lY(tDof, tnNo)
lD(tDof, tnNo)
 ```
 
 - `set_bc::set_bc_dir_l(com_mod, bc, com_mod.msh[iM].fa[iFa], tmpA, tmpY, lDof)`
 
<!-- =================== -->
<!--     set_bc_dir_l    -->
<!-- =================== -->

<h2 id="set_bc_dir_l"> void set_bc_dir_l(ComMod& com_mod, const bcType& lBc, const faceType& lFa, Array<double>& lA, Array<double>& lY, int lDof) </h2>

[set_bc.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/set_bc.cpp)

 Modifies
 ```
 ```
 
 - `igbc(com_mod, lBc.gm, lY, lA)` - calculating values by the inverse of general BC
 
 - `ifft(com_mod, lBc.gt, dirY_v, dirA_v)`
 
 
<!-- ==================== -->
<!--   iterate_solution   -->
<!-- ==================== -->

<h2 id="iterate_solution"> void iterate_solution(Simulation* simulation) </h2>

[simulation_wrap.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/simulation_wrap.cpp)
 
 - `read_msh_ns::calc_mesh_props(com_mod, cm_mod, com_mod.nMsh, com_mod.msh)`
 
 - `pic::picp(simulation)`
 
 - `set_bc::set_bc_dir(com_mod, An, Yn, Dn)`
 
 - Inner loop for iteration
   - `set_bc::set_bc_cpl(com_mod, cm_mod)`
   - `set_bc::set_bc_dir(com_mod, An, Yn, Dn)`
   - `pic::pici(simulation, Ag, Yg, Dg)`
   - `ls_ns::ls_alloc(com_mod, eq)`
   - `bf::set_bf(com_mod, Dg)`
   - `eq_assem::global_eq_assem(com_mod, com_mod.msh[iM], Ag, Yg, Dg)` - For each mesh
   - `set_bc::set_bc_neu(com_mod, cm_mod, Yg, Dg)`
   - `set_bc::set_bc_cmm(com_mod, cm_mod, Ag, Dg)`
   - `set_bc::set_bc_dir_w(com_mod, Yg, Dg)` - Apply weakly applied Dirichlet BCs
   - `contact::contact_forces(com_mod, cm_mod, Dg)`
   - `all_fun::commu(com_mod, com_mod.R)` - Synchronize residual across processes
   - `ustruct::ustruct_r(com_mod, Yg)`
   - `fs::thood_val_rc(com_mod)`
   - `set_bc::set_bc_undef_neu(com_mod)`
   - `ls_ns::ls_solve(com_mod, eq, incL, res)`
   - `pic::picc(simulation)`
   
<!-- ==================== -->
<!--          picp        -->
<!-- ==================== -->

<h2 id="picp"> void picp(Simulation* simulation) </h2>

[pic.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/pic.cpp)
 
Modifies
 ```
 pS0 
 Ad 
 Ao 
 Yo 
 Do 
 An 
 Yn 
 Dn 
```
 
 - `cep_ion::cep_integ(simulation, iEq, e, Do)`

 
<!-- ==================== -->
<!--      ls_alloc        -->
<!-- ==================== -->

<h2 id="ls_alloc"> void ls_alloc(ComMod& com_mod, eqType& lEq) </h2>

[ls.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/ls.cpp)
 
Allocate `com_mod.R` and `com_mod.Val` arrays.
 
Modifies
```
  com_mod.R - Residual vector
  com_mod.Val - LHS matrix
```
 
 - if using trilinos
   - `trilinos_lhs_create_(gtnNo, lhs.mynNo, tnNo, lhs.nnz, tls.ltg.data(), com_mod.ltg.data(), com_mod.rowPtr.data(),
        com_mod.colPtr.data(), dof)`

 
<!-- ==================== -->
<!--        set_bf        -->
<!-- ==================== -->

<h2 id="set_bf"> void set_bf(ComMod& com_mod, Array<double>& Dg) </h2>

[bf.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/bf.cpp)
 
 - `set_bf_l(com_mod, eq.bf[iBf], com_mod.msh[iM], Dg)`
 
<!-- ==================== -->
<!--        set_b_l       -->
<!-- ==================== -->

<h2 id="set_bf_l"> void set_bf_l(ComMod& com_mod, bfType& lBf, mshType& lM, Array<double>& Dg) </h2>

[bf.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/bf.cpp)
 
 Modifies
 ```
   com_mod.Bf
 ```
 
- `ifft(com_mod, lBf.bt, f, rtmp)` - If BodyForceType::bfType_ustd
 
- `igbc(com_mod, lBf.bm, bfl, xl)` - If BodyForceType::bfType_gen
 
- `bf_construct(com_mod, lM, e, eNoN, idof, xl, dl, bfl, ptr)`
 
<!-- ==================== -->
<!--     bf_construct     -->
<!-- ==================== -->

<h2 id="bf_construct"> void bf_construct(ComMod& com_mod, const mshType& lM, const int e, const int eNoN, const int idof, Array<double>& xl,
    const Array<double>& dl, const Array<double>& bfl, const Vector<int>& ptr) </h2>

[bf.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/bf.cpp)

 - `cmm::bcmmi(com_mod, eNoN, idof, w, N, Nx, xl, bfl, lR)`
 
 - `trilinos_doassem_(const_cast<int&>(eNoN), ptr.data(), lK.data(), lR.data())` - If using Trilinos
 
 - `lhsa_ns::do_assem(com_mod, eNoN, ptr, lK, lR)` - If not using Trilinos
 
 
<!-- ==================== -->
<!--    global_eq_assem   -->
<!-- ==================== -->

<h2 id="global_eq_assem"> void global_eq_assem(ComMod& com_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Yg, const Array<double>& Dg) </h2>

[eq_assem.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/eq_assem.cpp)

 - `fluid::construct_fluid(com_mod, lM, Ag, Yg) -  If EquationType::phys_fluid
 
<!-- ==================== -->
<!--    construct_fluid   -->
<!-- ==================== -->

<h2 id="construct_fluid"> void construct_fluid(ComMod& com_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Yg) </h2>

[fluid.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/fluid.cpp)
 
This is for solving fluid transport equation solving Navier-Stokes equations. Dirichlet boundary conditions are either treated
strongly or weakly.
 
 - `fs::get_thood_fs(com_mod, fs, lM, vmsStab, 1)`
 
 - `nn::gnn(fs[1].eNoN, nsd, nsd, Nx, xql, Nqx, Jac, ksix)` 
 
 - `nn::gn_nxx(l, fs[0].eNoN, nsd, nsd, Nx, Nxx, xwl, Nwx, Nwxx)`
 
 - `fluid_3d_m(com_mod, vmsStab, fs[0].eNoN, fs[1].eNoN, w, ksix, N0, N1, Nwx, Nqx, Nwxx, al, yl, bfl, lR, lK, K_inverse_darcy_permeability)` - If nsd=3
 
 - `fluid_2d_m(com_mod, vmsStab, fs[0].eNoN, fs[1].eNoN, w, ksix, N0, N1, Nwx, Nqx, Nwxx, al, yl, bfl, lR, lK, K_inverse_darcy_permeability)` - If nsd=2
 
 - `trilinos_doassem_(const_cast<int&>(eNoN), ptr.data(), lK.data(), lR.data())` - If using Trilinos
 
 - `lhsa_ns::do_assem(com_mod, eNoN, ptr, lK, lR)` - If not using Trilinos

<!-- ==================== -->
<!--       do_assem       -->
<!-- ==================== -->

<h2 id="do_assem"> void do_assem(ComMod& com_mod, const int d, const Vector<int>& eqN, const Array3<double>& lK, const Array<double>& lR) </h2>

[lhsa.cpp](https://github.com/ktbolt/svFSI/blob/Implement-svFSI-using-cpp_19/Code/Source/svFSI_cinterface/lhsa.cpp)

 Modifies
 ```
 com_mod.R
 com_mod.Val
 ```
 
 
 

 
