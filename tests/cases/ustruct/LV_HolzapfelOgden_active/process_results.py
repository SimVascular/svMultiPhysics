import pyvista as pv
import numpy as np
import os
import matplotlib.pyplot as plt


def get_array_at_point(mesh, point, array_name):
    """
    Get the value of an array at a point in a mesh
    """
    point = pv.PolyData(point)
    sampled_point = point.sample(mesh)
    return sampled_point[array_name]

def compute_volume(ref_mesh_with_displacement, ref_lumen):
    """
    Compute the volume enclosed by a surface mesh warped by a displacement field
    """

    # Sample displacement field at the lumen points
    resampled_lumen = ref_lumen.sample(ref_mesh_with_displacement)

    # Warp the lumen mesh
    warped_lumen = ref_lumen.warp_by_vector('Displacement')

    # Compute the volume enclosed by the warped lumen
    volume = warped_lumen.volume

    return volume


def sort_files(results_path):
    """
    Sort .vtu files by the numeric digits before .vtu and after result_
    """

    files = os.listdir(results_path)
    files = [file for file in files if file.endswith('.vtu')]
    files = sorted(files, key=lambda x: int(x.split('_')[-1].split('.')[0]))

    return files

def get_displacements_at_points(start_timestep, end_timestep, step, timestep_size, results_folder, sample_points):
    """
    Get the displacements at a list of points in a series of .vtu files
    
    Args:
        start_timestep: The first svFSI result file to process
        
        end_timestep: The last svFSI result file to process
        
        step: The step in svFSI result files to process

        timestep_size: The size of the timestep in seconds
        
        results_folder: The absolute file path of the svFSI results folder 
        (usually something/something/16-procs/)
        
        sample_points: List of points at which to sample displacements

    Returns: (t, displacements), a tuple where t contains the time points and 
    displacements is a 3D array of shape (n_time_steps, n_points, 3) containing 
    the displacement vectors at each point and time.
    """
    
    print('\n## Calculating displacements at points ##')

    # Initialize arrays to store time and displacements
    t = []
    displacements_list = []
    
    # If start_timestep is 0, this is a special case, as we need to handle the reference state
    # The smallest first result file is result_001.vtu, so we start at step for actual results
    if start_timestep == 0:
        # For time 0, we assume zero displacement (reference state)
        t.append(0)
        zero_displacements = np.zeros((len(sample_points), 3))
        displacements_list.append(zero_displacements)
        print(f"Time: {0}, Reference state (zero displacement)")
        start_timestep = step
    
    # Loop through results files at each time > 0
    for k in range(start_timestep, end_timestep+1, step):
        # Load results VTU mesh
        result_file = os.path.join(results_folder, f"result_{k:03d}.vtu")
        
        if not os.path.exists(result_file):
            print(f"Warning: File {result_file} not found, skipping...")
            continue
            
        result = pv.read(result_file)
        
        # Initialize array for displacements at this time step
        time_displacements = np.zeros((len(sample_points), 3))
        
        # Get displacement at each sample point
        for j, point in enumerate(sample_points):
            displacement = get_array_at_point(result, point, 'Displacement')
            time_displacements[j, :] = displacement.reshape(3)
        
        # Add time and displacements to arrays
        t.append(k * timestep_size)
        displacements_list.append(time_displacements)
        
        print(f"Time: {k * timestep_size:.3f}s, Processed file: result_{k:03d}.vtu")
    
    # Convert list to numpy array
    displacements = np.array(displacements_list)
    
    return (t, displacements)

def calc_volume_3D(start_timestep, end_timestep, step, timestep_size, results_folder, reference_surface, save_intermediate_data=False, intermediate_output_folder=None):
    """
    Calculate the ventricular lumen volume at each time step from the results of 
    an svFSI struct simulation, in which a model of the myocardium is simulated

    Calculate the volume in the following steps
    1) Sample the result.vtu file onto the reference surface
    2) Warp the samples surface by the Displacement
    3) Flat fill any holes in the warped surface
    4) Calculate the volume of the warped and filled surface

    The units of volume are whatever units used in .vtu files, cubed. For example,
    if units of length in the .vtu files are microns, then the volume calculated
    here is cubic microns. 

    Args:
        start_timestep: The first svFSI result file to process
        
        end_timestep: The last svFSI result file to process
        
        step: The step in svFSI result files to process

        timestep_size: The size of the timestep in seconds
        
        results_folder: The absolute file path of the svFSI results folder 
        (usually something/something/16-procs/)
        
        reference_surface: The absolute file path of the .vtp file containing 
        the undeformed surface corresponding to the deformed surface of which 
        we want to compute the volume.

        save_intermediate_data: Whether to save intermediate data (resampled, warped, and filled surfaces)

        intermediate_output_folder: The folder to save the intermediate data in. If None, the intermediate data is saved in the same folder as the reference_surface.

    Returns: (t, vol), a tuple of lists of length number of time steps. t 
    contains the time, and vol contains the volume at that time.
    """
    
    # Create folder to contain intermediary meshes (mostly for checking for errors)
    if save_intermediate_data:
        assert intermediate_output_folder is not None, "If save_intermediate_data is True, intermediate_output_folder must be provided"
            
        # checking if the directory exists
        if not os.path.exists(intermediate_output_folder):
            # if the directory is not present then create it.
            os.makedirs(intermediate_output_folder)

    print('\n## Calculating volumes ##')

    # Load reference surface onto which we sample
    ref_surface = pv.read(f"{reference_surface}")
    
    # Initialize arrays to store time and volume
    t = []
    vol = []

    # Make reference surface watertight
    ref_lumen = ref_surface.fill_holes(100) # 100 is the largest size of hole to fill
    
    # Recompute normals, incase the normals of the cap are opposite
    ref_lumen.compute_normals(inplace=True)

    # Save filled lumen (to check geometry and normals)
    # (Hopefully the normals on the filled cap will be consistent with the normals
    # on the rest of the surface, but you should check to make sure.)
    if save_intermediate_data:
        ref_lumen.save(os.path.join(intermediate_output_folder,  f'resampled_warped_and_filled_{0:03d}.vtp'))
    
    # Compute volume of ref_lumen
    print(f"Iteration: {0}, Volume: {ref_lumen.volume}")
    t.append(0)
    vol.append(ref_lumen.volume)

    # If start_timestep is 0, this is a special case, as we have already computed the volume
    # of the reference surface. Also, the smallest first result file is result_001.vtu,
    # so we start at step.
    if start_timestep == 0:
        start_timestep = step
    
    # Loop through results files at each time > 0
    for k in range(start_timestep, end_timestep+1, step):
        # Load results VTU mesh
        result = pv.read(os.path.join(results_folder, f"result_{k:03d}.vtu"))

        # Sample result onto ref_lumen
        resampled_lumen = ref_lumen.sample(result)

        # Warp resampled surface by displacement (needed for current configuration 
        # normals, as well volume calculation)
        warped_lumen = resampled_lumen.warp_by_vector('Displacement')

        # Save warped and filled lumen (to check geometry and normals)
        if save_intermediate_data:
            warped_lumen.save(os.path.join(intermediate_output_folder, f'resampled_warped_and_filled_{k:03d}.vtp'))
        
        # Add time and volume to arrays
        t.append(k*timestep_size)
        vol.append(warped_lumen.volume)

        print(f"Iteration: {k}, Volume: {warped_lumen.volume}")
    
    return (t, vol)

# Set points needed to calculate the displacements
p0 = np.array([0.025,0.03,0.0])
p1 = np.array([0.0,0.03,0.0])

sample_points = [p0, p1]

results_path = 'results'

# Get the displacements at the sample points using the new function signature
t_displacements, displacements = get_displacements_at_points(0, 1000, 10, 1e-3, results_path, sample_points)

data = np.load('Step_1_US_P1_h5.npz')
#print(data['time'])

coords = ['x', 'y', 'z']

# Plot the displacements for both points and for each of three coordinates on  2 x 3 grid
plt.figure()
for i in range(3):
    plt.subplot(3,2,2*i+1)
    plt.plot(t_displacements, displacements[:,0,i], label='svMultiPhysics')
    plt.plot(data['time'], data['u_0'][:,i], label = 'benchmark')
    plt.xlabel('Time [s]')
    plt.ylabel(coords[i] + ' Displacement [m]')
    plt.legend()
    plt.tight_layout()

    plt.subplot(3,2,2*i+2)
    plt.plot(t_displacements, displacements[:,1,i], label='svMultiPhysics')
    plt.plot(data['time'], data['u_1'][:,i], label = 'benchmark')
    plt.xlabel('Time [s]')
    plt.ylabel(coords[i] + ' Displacement [m]')
    plt.legend()
    plt.tight_layout()

# label column 1 'p_0' and column 2 'p_1'
plt.subplot(3,2,1)
plt.title(r'$p_0$')
plt.subplot(3,2,2)
plt.title(r'$p_1$')

plt.tight_layout()
plt.savefig('displacements.png')
plt.show()


# Plot P-V loop
t, vol = calc_volume_3D(0, 1000, 10, 1e-3, results_path, 'mesh/mesh-surfaces/endo.vtp')
vol = np.array(vol) * 1e6 # Convert to mL

# Load pressure data from the second column of pressure.dat
pressure_data = np.loadtxt('pressure.dat')
pressure_time = pressure_data[1:, 0]  # First column: time
pressure_values = pressure_data[1:, 1] / 133.322387415 # Second column: pressure

# Interpolate pressure to match volume time points
from scipy.interpolate import interp1d
pressure_interp = interp1d(pressure_time, pressure_values, kind='linear', bounds_error=False, fill_value='extrapolate')
pressure_at_volume_times = pressure_interp(t)

plt.figure(figsize=(10, 8))
plt.plot(vol, pressure_at_volume_times, 'b-', linewidth=2, label='Pressure-Volume Loop')
plt.xlabel('Volume [mL]', fontsize=12)
plt.ylabel('Pressure [mmHg]', fontsize=12)
plt.title('Pressure-Volume Loop', fontsize=14)
plt.grid(True, alpha=0.3)
#plt.legend(fontsize=12)

plt.tight_layout()
plt.savefig('p-v_loop.png', dpi=300, bbox_inches='tight')
plt.show()

# Print some statistics
print(f"Maximum volume: {np.max(vol):.2f} mL")
print(f"Minimum volume: {np.min(vol):.2f} mL")
print(f"Stroke volume: {np.max(vol) - np.min(vol):.2f} mL")
print(f"Maximum pressure: {np.max(pressure_at_volume_times):.2f} mmHg")
print(f"Minimum pressure: {np.min(pressure_at_volume_times):.2f} mmHg")
