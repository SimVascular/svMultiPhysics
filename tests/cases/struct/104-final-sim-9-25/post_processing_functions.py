# Before running on Sherlock, need to have installed: 
#
# pyvista: pip3 install pyvista
# ffmpy: pip3 install ffmpy
#
# Also, need to load the following required modules on Sherlock: 
#
# ml system mesa ffmpeg
import glob
import os
import pyvista as pv # for VTK mesh manipulations
import numpy as np 
#import ffmpy # to convert .avi movie to .mp4 movie

def get_timestep(result_file):
    '''
    Extracts the timestep of a results_###.vtu file
    '''
    # If file path is provided, get the file name
    file_name = os.path.basename(result_file)

    # Get the ###.vtu part
    s = file_name.split('_')[-1]

    # Get the ### part
    s = s.split('.')[0]

    # Return as integer
    return int(s)

def get_start_end_step(results_folder):
    """
    Automatically determine the start time, end time, and step size based on all
    svFSI results file in results_folder.

    Args:
        results_folder: A string of absolute file path of folder containing results of 
        svFSI simulation. This usually ends with 16-procs/ or other 
        number.
    
    Returns:
        (start_time, end_time, step): A tuple of 3 integers, giving the first 
        time step of results to process, the last time step, and the step size. 

    """

    # Get list of all .vtu files in results_folder sorted by time step
    list_of_files = sorted( filter( os.path.isfile,
                            glob.glob(os.path.join(results_folder, '*.vtu')) ), key = get_timestep)

    # Get start time from the first result file (list_of_files[0])
    start_file_name = os.path.basename(list_of_files[0])
    start_time = int("".join([i for i in start_file_name if i.isdigit()]))
    print('Start time:', start_time)

    # Get end time from the last result file (list_of_files[-1])
    end_file_name = os.path.basename(list_of_files[-1])
    end_time = int("".join([i for i in end_file_name if i.isdigit()]))
    print('End time:', end_time)

    # Get step size by looking at second time step
    start_plus_one_file_name = os.path.basename(list_of_files[1])
    start_time_plus_one = int("".join([i for i in start_plus_one_file_name if i.isdigit()]))
    step = start_time_plus_one - start_time
    print('Step:', step)

    return (start_time, end_time, step)

def calc_volume_struct(start_time, end_time, step, results_folder, reference_surface):
        """
        Calculate the ventricular lumen volume at each time step from the results of 
        an svFSI struct simulation, in which a model of the myocardium is inflated.

        Calculate the volume in the following steps
        1) Sample the result.vtu file onto the reference surface
        2) Warp the samples surface by the Displacement
        3) Flat fill any holes in the warped surface
        4) Calculate the volume of the warped and filled surface

        The units of volume are whatever units used in .vtu files, cubed. For example,
        if units of length in the .vtu files are microns, then the volume calculated
        here is cubic microns. 

        Args:
            start_time: The first svFSI result file to process
            
            end_time: The last svFSI result file to process
            
            step: The step in svFSI result files to process
            
            results_folder: The absolute file path of the svFSI results folder 
            (usually something/something/16-procs/)
            
            reference_surface: The absolute file path of the .vtp file containing 
            the undeformed surface corresponding to the deformed surface of which 
            we want to compute the volume.

            output_folder: The absolute file path of the folder to which we output
            intermediary results of the volume calculation process. 

        Returns: (t, vol), a tuple of lists of length number of time steps. t 
        contains the time step, and vol contains the volume at that time step.
        """
        
        # Create folder to contain intermediary meshes (mostly for checking for errors)
        output_folder = results_folder + '/../' + 'calc_volume_struct_output'
        # checking if the directory exists
        if not os.path.exists(output_folder):
            # if the directory is not present then create it.
            os.makedirs(output_folder)

        print('\n## Calculating volumes ##')

        # Load reference surface onto which we sample
        ref_surface = pv.read(f"{reference_surface}")
        
        # Compute initial lumen volume
        t = []
        vol = []
        k = 0
        # The initial volume is obtain by computing the lumen volume of the reference
        # configuration, which is simply volume of the filled in reference surface.
        ref_lumen = ref_surface.fill_holes(100) # 100 is the largest size of hole to fill

        # --------------------------------------------------------------------
        
        # Recompute normals, incase the normals of the cap are opposite
        ref_lumen.compute_normals(inplace=True)

        # Save warped and filled lumen (to check geometry and normals)
        # (Hopefully the normals on the filled cap will be consistent with the normals
        # on the rest of the surface, but you should check to make sure.)
        ref_lumen.save(f'{output_folder}/resampled_warped_and_filled_{k:03d}.vtp')

        # Add time and volume to arrays
        t.append(k)
        vol.append(ref_lumen.volume)

        print(f"Iteration: {k}, Volume: {ref_lumen.volume}")


        # Loop through results files at each time > 0
        for k in range(start_time, end_time+1, step):
            # Load results VTU mesh
            result = pv.read(os.path.join(results_folder, f"result_{k:03d}.vtu"))

            # Sample result onto ref_surface
            resampled = ref_surface.sample(result)

            # Warp resampled surface by displacement (needed for current configuration 
            # normals, as well volume calculation)
            warped = resampled.warp_by_vector('Displacement')

            ## ------ Compute volume by warping resampled, then filling the hole ----- ##

            #warped.save(f'{output_folder}/resampled_warped_{k:03d}.vtp')

            # Fill the holes
            #                           EITHER
            # -----------------------------------------------------------------------
            # Fix the mesh (i.e. fill the holes where the inlet and outlet caps are)
            #meshfix = mf.MeshFix(warped)
            #meshfix.repair()

            # Convert back to pyvista polydata mesh
            #lumen = meshfix.mesh

            ## PyMeshFix fills the hole, but leads to normals oriented inconsistently ##
            ## The PyVista function .fill_holes() produces a bad cap mesh, but with   ##
            ## consistent normal vectors. The bad mesh is not important if we're only ##
            ## computing volume                                                       ##
            
            # ---------------------------- OR ------------------------------------
            lumen = warped.fill_holes(100) # 100 is the largest size of hole to fill

            # --------------------------------------------------------------------
            
            # Recompute normals, incase the normals of the cap are opposite
            lumen.compute_normals(inplace=True)

            # Save warped and filled lumen (to check geometry and normals)
            # (Hopefully the normals on the filled cap will be consistent with the normals
            # on the rest of the surface, but you should check to make sure.)
            lumen.save(f'{output_folder}/resampled_warped_and_filled_{k:03d}.vtp')

            # Add time and volume to arrays
            t.append(k)
            vol.append(lumen.volume)

            print(f"Iteration: {k}, Volume: {lumen.volume}")
        
        return (t, vol)
def calc_dVdt_struct(start_time, end_time, step, results_folder, reference_surface):
    '''
    Computes the rate of change of volume of a closed (or partially closed) 
    surface. Intended to be used to calculate the rate of change of ventricular
    volume from a struct simulation of an inflating ventricle.

    Compute the rate of change of volume by computing the velocity flux over the
    reference surface, which should be the endo surface of the ventricle.

    dV/dt = int_{Gamma_t} (u . n) dA

    where Gamma_t is the reference surface in the current configuration, u is 
    the velocity on the reference surface, and n is the surface normal vector on
    the reference surface

    Args:
        start_time: The first svFSI result file to process
        
        end_time: The last svFSI result file to process
        
        step: The step in svFSI result files to process
        
        results_folder: The absolute file path of the svFSI results folder 
        (usually something/something/16-procs/)
        
        reference_surface: The absolute file path of the .vtp file containing 
        the undeformed surface corresponding to the deformed surface of which 
        we want to compute the volume.

        output_folder: The absolute file path of the folder to which we output
        intermediary results of the volume calculation process. 

    Returns: (t, dVdt), a tuple of lists of length number of time steps. t 
    contains the time step, and Q contains the rate of change of volume at that 
    time step.
    '''

    # Create folder to contain intermediary meshes (mostly for checking for errors)
    output_folder = results_folder + '/../' + 'calc_volume_struct_output'
    # checking if the directory exists
    if not os.path.exists(output_folder):
        # if the directory is not present then create it.
        os.makedirs(output_folder)
    
    # Initialize volume array and time step lists
    dVdt = [0]
    t = [0]

    # Loop through results files at each time
    for k in range(start_time, end_time+1, step):
        # Load results VTU mesh
        result = pv.read(os.path.join(results_folder, f"result_{k:03d}.vtu"))

        # Load reference surface onto which we sample
        ref_surface = pv.read(f"{reference_surface}")

        # Sample result onto ref_surface
        ref= ref_surface.sample(result)

        # Warp ref surface by displacement (needed for current configuration 
        # normals, as well volume calculation)
        warped = ref.warp_by_vector('Displacement')

        # Compute average velocity for each element
        warped = warped.point_data_to_cell_data()
        ref = ref.point_data_to_cell_data()

        # Recompute normals. Need this after point_data_to_cell_data(), because 
        # the average of point normals is not the same as a cell normal
        warped.compute_normals(inplace=True)
        ref.compute_normals(inplace=True)

        # Compute element areas
        warped = warped.compute_cell_sizes()
        ref = ref.compute_cell_sizes()

        # Save warped and filled lumen (to check geometry and normals)
        # (Hopefully the normals on the filled cap will be consistent with the normals
        # on the rest of the surface, but you should check to make sure.)
        warped.save(f'{output_folder}/warped_{k:03d}.vtp')
        ref.save(f'{output_folder}/ref_{k:03d}.vtp')

        # Compute dVdt from warped by computing velocity flux. Also, compute 
        # velocity flux using quantities from reference surface to compare.
        # 
        # !! Looks like the velocity flux output in B_ST_Velocity_flux.txt is 
        # computed using reference surface normals and reference surface areas !! 
        dVdt_t = 0
        
        cell_normals = warped.cell_data['Normals']
        cell_vels = warped.cell_data['Velocity']
        cell_areas = warped.cell_data['Area']

        cell_normals_ref = ref.cell_data['Normals']
        cell_vels_ref = ref.cell_data['Velocity']
        cell_areas_ref = ref.cell_data['Area']

        u_dot_n = np.sum(cell_vels * cell_normals, axis = 1)
        dVdt_t = np.sum(u_dot_n * cell_areas)

        u_dot_n_ref = np.sum(cell_vels_ref * cell_normals_ref, axis = 1)
        dVdt_t_ref = np.sum(u_dot_n_ref * cell_areas_ref)

        # Add time and volume to arrays
        t.append(k)
        dVdt.append(dVdt_t)

        print(f"Iteration: {k}, dVdt: {dVdt_t}, dVdt_ref: {dVdt_t_ref}")
    
    return (t, dVdt)

def calc_pressure_struct(input_file, pressure_dat_file, t):
    """
    Calculate the ventricular lumen pressure at each time step from the results of 
    an svFSI struct simulation, in which a model of the myocardium is inflated.

    Calculates pressure on the endo surface by reading the input file and 
    pressure load file. The units of pressure are whatever units are used in the
    pressure load file (usually Pa)

    Args:
        input_file: The svFSI input file of the simulation whose results we are
        now processing.

        pressure_dat_file: The file containing the pressure load information for
        the simulation whose results we are now processing.

        t: List of time steps (iterations) at which we want to output the 
        pressure.
    Returns:
        pressure: Pressure evaluated at list of time steps in t
    """
    print('\n## Calculating pressure ##')

    # Compute pressure at each time step by reading input file and pressure file
    with open(input_file) as f:
        lines = f.readlines()
        for line in lines:
            if 'Number of time steps' in line:
                num_time_steps = int(line.split(':')[-1])
        
    with open(pressure_dat_file) as f:
        f.readline()
        line = f.readline()
        vals = line.split()
        t0 = float(vals[0])
        p0 = float(vals[1])
        line = f.readline()
        vals = line.split()
        t1 = float(vals[0])
        p1 = float(vals[1])

    pressure = np.linspace(p0, p1, num_time_steps+1)

    # Sample pressure onto t
    pressure = pressure[t]

    return pressure

def calc_pressure_struct_genBC(genBC_out_file, col, t):
    """
    Calculate the ventricular lumen pressure at each time step from the results of 
    an svFSI struct simulation with genBC, in which a model of the myocardium is 
    inflated.

    Calculates pressure on the endo surface by reading the genBC output file,
    usually called AllData. The units of pressure are whatever units are used in 
    the genBC input file.

    Args:
        genBC_out_file: The file containing the pressure load information for
        the simulation whose results we are now processing. Usually AllData

        col: The column in AllData that contains the pressure at each time step
        col = 1 refers to the first column

        t: List of time steps (iterations) at which we want to output the 
        pressure.
    Returns:
        pressure: Pressure evaluated at list of time steps in t
    """
    print('\n## Calculating pressure ##')

    # Read pressure from AllData
    pressure = [0]
    with open(genBC_out_file) as f:
        lines = f.readlines()
        for line in lines:
            vals = line.split()
            pressure.append(float(vals[col-1]))

    pressure = np.array(pressure)

    # Sample pressure onto t
    pressure = pressure[t]

    return pressure
def calc_volume_struct_genBC(genBC_out_file, col, t):
    """
    Calculate the accumulated lumen volume at each time step from the results of 
    an svFSI struct simulation with genBC, in which a model of the myocardium is 
    inflated.

    Calculates accumulated volumeby reading the genBC output file,
    usually called AllData. In genBC, you must define an additional unknown 
    x(i) such that f(i) = -Q(1). This will compute x(i) as the integral of the
    flow rate, i.e. the accumulated volume

    Args:
        genBC_out_file: The file containing the pressure load information for
        the simulation whose results we are now processing. Usually AllData

        col: The column in AllData that contains the pressure at each time step
        col = 1 refers to the first column

        t: List of time steps (iterations) at which we want to output the 
        pressure.
    Returns:
        volume: Volume evaluated at list of time steps in t
    """
    print('\n## Calculating volume from AllData ##')

    # Read pressure from AllData
    volume = [0]
    with open(genBC_out_file) as f:
        lines = f.readlines()
        for line in lines:
            vals = line.split()
            volume.append(float(vals[col-1]))

    volume = np.array(volume)

    # Sample pressure onto t
    volume = volume[t]

    return volume

def calc_volume_FSI(volume_file):
    """
    Calculate the ventricular lumen volume at each time step from the results of 
    an svFSI FSI simulation, in which a model of the myocardium is inflated.

    Calculate the volume by reading the output file V_FS_Volume_average.txt. This
    is produced if the following is included in svFSI.inp

     Output: Volume_integral {
      Volume: t
    }

    This gives the volume of the fluid domain at each time step.

    The units of volume are whatever units used in svFSI.inp, cubed. For example,
    if units of length in svFSI.inp are meters, then the volume calculated
    here is cubic meters. Note that the mesh scale factor affects these units. If
    the .vtu file has units of microns, but we apply a mesh scale factor of 1e-6
    to convert to meters in svFSI.inp, then the volume output here will be in meters^3 

    Args:
        volume_file: Path to file V_FS_Volume_average.txt from simulation you 
        want to process

    Returns: (t, vol), a tuple of lists of length number of time steps. t 
    contains the time step, and vol contains the volume at that time step.
    """
    print('\n## Calculating volume ##')

    # Initialize time and volume arrays
    t = []
    vol = []

    # Read volume data from V_FS_Volume_average.txt 
    with open(volume_file) as f:
        header = f.readline()               # Read header
        initial_volumes = f.readline()      # Read domain initial volumes
        f.readline()                        # Skip blank line

        # Read volume data
        # (For some reason, the first volume is the initial volume, even though
        # this corresponds to after the first time step)
        lines = f.readlines()
        k = 0
        for line in lines:

            # Add time iteration to t array
            t.append(k)
            k += 1

            # Add fluid domain (Domain 0) volume to volume array
            V_t = float(line.split()[0])
            vol.append(V_t)

    return (t, vol)


def calc_volume_FSI_2(start_time, end_time, step, results_folder):
    """
    Calculate the ventricular lumen volume at each time step from the results of 
    an svFSI FSI simulation, in which a model of the myocardium is inflated.

    Calculate the volume in the following steps
    1) Threshold to extract fluid domain
    2) Warp the fluid domain by the Displacement
    3) Calculate the volume of the warped fluid domain

    The units of volume are whatever units used in .vtu files, cubed. For example,
    if units of length in the .vtu files are microns, then the volume calculated
    here is cubic microns. 

    Args:
        start_time: The first svFSI result file to process
        
        end_time: The last svFSI result file to process
        
        step: The step in svFSI result files to process
        
        results_folder: The absolute file path of the svFSI results folder 
        (usually something/something/16-procs/)

    Returns: (t, vol), a tuple of lists of length number of time steps. t 
    contains the time step, and vol contains the volume at that time step.
    """
    
    # Create folder to contain intermediary meshes (mostly for checking for errors)
    output_folder = results_folder + '/../' + 'calc_volume_FSI_output'
    # checking if the directory exists
    if not os.path.exists(output_folder):
        # if the directory is not present then create it.
        os.makedirs(output_folder)

    print('\n## Calculating volumes ##')

    # Initialize volume array and time step lists
    vol = [0.0]
    t = [0]

    # Loop through results files at each time
    for k in range(start_time, end_time+1, step):
        # Load results VTU mesh
        result = pv.read(os.path.join(results_folder, f"result_{k:03d}.vtu"))

        # Threshold to fluid domain
        fluid = result.threshold(value = (1,1), scalars = 'Domain_ID')

        # Warp fluid domain by displacement (yields current fluid domain)
        warped = fluid.warp_by_vector('Displacement')

        # Save warped fluid domain (to check geometry and normals)
        # (Hopefully the normals on the filled cap will be consistent with the normals
        # on the rest of the surface, but you should check to make sure.)
        #warped.save(f'{output_folder}/fluid_warped_{k:03d}.vtp')
        warped.save(f'{output_folder}/fluid_warped_{k:03d}.vtu')

        

        # Add time and volume to arrays
        t.append(k)
        vol.append(warped.volume)

        print(f"Iteration: {k}, Volume: {warped.volume}")
    
    return (t, vol)


def calc_pressure_FSI(pressure_file):
    """
    Calculate the ventricular lumen pressure at each time step from the results of 
    an svFSI FSI simulation, in which a model of the myocardium is inflated.

    Calculate the pressure by reading the output file V_FS_Pressure_average.txt. This
    is produced if the following is included in svFSI.inp

     Output: Volume_integral {
      Pressure: t
    }
    This gives the volume averaged fluid pressure at each time step.


    The units of pressure are deduced by whatever units system you are using in
    the input file. For SI (m, kg, s), the units of pressure are Pa.

    Args:
        pressure_file: Path to file V_FS_Pressure_average.txt from simulation you 
        want to process

    Returns: (t, pressure), a tuple of lists of length number of time steps. t 
    contains the time step, and pressure contains the pressure at that time step.
    """
    print('\n## Calculating pressure ##')

    # Initialize time and pressure arrays
    t = [0]
    pressure = [0]

    # Read pressure data from V_FS_Pressure_average.txt 
    with open(pressure_file) as f:
        header = f.readline()               # Read header
        initial_volumes = f.readline()      # Read domain initial volumes
        f.readline()                        # Skip blank line

        # Read pressure data
        lines = f.readlines()
        for line in lines:

            # Add next time iteration to t array
            t.append(t[-1] + 1)

            # Add fluid domain (Domain 0) pressure to pressure array
            p_t = float(line.split()[0])
            pressure.append(p_t)

    return (t, pressure)

def calc_pressure_FSI_2(start_time, end_time, step, results_folder):
    """
    Calculate the ventricular lumen pressure at each time step from the results of 
    an svFSI FSI simulation, in which a model of the myocardium is inflated.

    Calculate the pressure in the following steps
    1) Threshold to extract fluid domain
    2) Warp the fluid domain by the Displacement
    3) Compute the volume of each cell
    4) Compute the pressure in each cell
    5) Compute the volume-averaged pressure

    The units of pressure are deduced by whatever units system you are using in
    the input file. For SI (m, kg, s), the units of pressure are Pa.

    Args:
        start_time: The first svFSI result file to process
        
        end_time: The last svFSI result file to process
        
        step: The step in svFSI result files to process
        
        results_folder: The absolute file path of the svFSI results folder 
        (usually something/something/16-procs/)

    Returns: (t, pressure), a tuple of lists of length number of time steps. t 
    contains the time step, and pressure contains the pressure at that time step.
    """
    
    # Create folder to contain intermediary meshes (mostly for checking for errors)
    output_folder = results_folder + '/../' + 'calc_pressure_FSI_output'
    # checking if the directory exists
    if not os.path.exists(output_folder):
        # if the directory is not present then create it.
        os.makedirs(output_folder)

    print('\n## Calculating pressures ##')

    # Initialize volume array and time step lists
    pressure = [0.0]
    t = [0]

    # Loop through results files at each time
    for k in range(start_time, end_time+1, step):
        # Load results VTU mesh
        result = pv.read(os.path.join(results_folder, f"result_{k:03d}.vtu"))

        # Threshold to fluid domain
        fluid = result.threshold(value = (1,1), scalars = 'Domain_ID')

        # Warp fluid domain by displacement (yields current fluid domain)
        warped = fluid.warp_by_vector('Displacement')

        # Compute volumes and areas
        sized = warped.compute_cell_sizes()
    

        # Compute pressure in each cell from pressure at each point
        sized = sized.point_data_to_cell_data()



        # Save processed fluid domain (to check geometry, normals, and values)
        # (Hopefully the normals on the filled cap will be consistent with the normals
        # on the rest of the surface, but you should check to make sure.)
        #warped.save(f'{output_folder}/fluid_warped_{k:03d}.vtp')
        sized.save(f'{output_folder}/fluid_processed_{k:03d}.vtu')

        # Grab volumes for all elements in mesh
        elem_volumes = sized.cell_data['Volume']
        # Grab pressure for all elements in mesh
        elem_pressures = sized.cell_data['Pressure']

        # Compute volume averaged pressure by looping over fluid domain elements
        p_avg = 0
        for i in range(len(elem_volumes)):
            p_avg += elem_pressures[i] * elem_volumes[i]
        p_avg = p_avg / sum(elem_volumes)
        

        # Add time and pressure to arrays
        t.append(k)
        pressure.append(p_avg)

        print(f"Iteration: {k}, Pressure: {p_avg}")
    
    return (t, pressure)

#def avi2mp4(input_movie_file, output_movie_file):
#    ff = ffmpy.FFmpeg(
#        inputs={input_movie_file: None},
#        outputs={output_movie_file: '-pix_fmt yuv420p -y'}
#        # flag -pix_fmt yuv420p needed so that movie is compatible with QuickTime
#        # flag -y needed to overwrite mp4 file if it already exists
#        )
#    ff.run()