# Before running on Sherlock, need to have installed: 
#
# pyvista: pip3 install pyvista
# ffmpy: pip3 install ffmpy
#
# Also, need to load the following required modules on Sherlock: 
#
# ml system mesa ffmpeg
#
# mesa is required just to allow pyvista to run properly. ffmpeg is required for ffpy
# to convert the movie.avi file to movie.mpy. system is required before mesa or ffmpeg
# can be loaded

# USAGE:
# python3 process_results_struct_only_Sherlock.py
# or
# python3 process_results_struct_only_Sherlock.py <simulation_folder>


# Python script to calculate lumen volume and fluid pressure on endo surface from 
# results of svFSI struct, and convert a .avi animation of result to .mp4
# 
# Requires two meshes
# - svFSI result.vtu (with Displacement array)
# - reference surface file (polygonal mesh), representing the undeformed surface
#   corresponding to the deformed surface of which we want to compute the volume
#
# We calculate the volume in the following steps
# 1) Sample the result.vtu file onto the reference surface
# 2) Warp the samples surface by the Displacement
# 3) Flat fill any holes in the warped surface
# 4) Calculate the volume of the warped and filled surface
#
# We calculate the fluid pressure on the endo surface by reading the input file
# and pressure load file.
# 
# This script USED TO use pymeshfix, which is a library for "cleaning up" surface meshes
# that should represent solid objects. This is used to flat fill the holes.
# It also uses pyvista, a Python interface to VTK

import os # for checking and creating directories, and loading modules
#import pymeshfix as mf # for "fixing" meshes (in our case, filling holes of endo surface)
import matplotlib.pyplot as plt
import numpy as np
import sys

# Add post-processing folder to Python path so we can access functions in it
sys.path.insert(0, '/Users/aaronbrown/GitHub/Cardiac_Modeling_Tools/useful-functions/Aaron_post_processing')
sys.path.insert(0, '/scratch/users/abrown97/sandbox/3D0D_coupling_sims')
# Import custom functions useful for post-processing
from post_processing_functions import *

## -------- PARAMETERS TO CHANGE -------------------- ##

# Simulation folder file path
sim_folder = '/Users/dseyler/Documents/Marsden_Lab/svMultiPhysics/tests/cases/struct/104-final-sim-9-25'

# Optionally read sim_folder as command line argument
if len(sys.argv) == 2:
    print(sys.argv)
    sim_folder = os.path.abspath(sys.argv[1])
    print(sim_folder)


# Undeformed surface, corresponding to deformed surface on result.vtu of which we want to compute the volume
reference_surface = os.path.join(sim_folder, 'P1_0.005_Exact/mesh-surfaces/endo_capped.vtp')

# Undeformed vtu file (the original volume mesh)
reference_volume = os.path.join(sim_folder, 'P1_0.005_Exact/mesh-complete.mesh.vtu')

# svFSI results folder, containing results.vtu
results_folder = os.path.join(sim_folder, '1-procs')

# Input file path
input_file = os.path.join(sim_folder, 'solver.xml')

# Pressure file (unused)
pressure_dat_file = os.path.join(sim_folder, 'pressure.dat')


## ----- For Paraview movie------ ##
# Camera Position and View Up. In Paraview, query the current camera position by clicking on
# the Adjust Camera button (one of the camera icons in the row icons containing
# Hover Points On, Hover Cell On, etc.)
camera_position = [-53.38446917875837, -243.23976974151938, 172.3987150979099]
view_up = [0.1822109963898, 0.5416373559188973, 0.5416373559188973]

# Path to pvbatch, required to run Paraview via Python
#pvbatch_path = "/home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch"

# Path to Paraview Python script used to create a paraview animation
#pvpython_script_path = "/scratch/users/abrown97/ZebrafishProjects/01-ventricle_inflation/101-post-processing/create_paraview_movie_struct.py"

## -------------------- END PARAMETERS TO CHANGE ------------------------ ## 


# Automatically determine the start time, end time, and step size based on all
# results file in results_folder
(start_time, end_time, step) = (10, 490, 10)#get_start_end_step(results_folder)
# Option to manually set start time, end time, and time step of results files to process
#start_time = 5
#end_time = 85
#step = 5

# Compute lumen volume from simulation results
(t, vol) = calc_volume_struct(start_time, end_time, step, results_folder, reference_surface)
vol_cm3 = np.array(vol) # cm^3

# Compute lumen volume from AllData file
vol_0D = calc_volume_struct_genBC(os.path.join(sim_folder, "AllData"), 6, t)
# Add on initial volume
vol_0D += vol_cm3[0] # cm^3

# Compute flow rate = dVdt from AllData file
dVdt_0D = calc_volume_struct_genBC(os.path.join(sim_folder, "AllData"), 15, t)

# Compute lumen dVdt from simulation results
(t_dVdt, dVdt) = calc_dVdt_struct(start_time, end_time, step, results_folder, reference_surface)
# Convert volume to cubic meters/s
dVdt = np.array(dVdt) # m/s * cm^2
dVdt_m3 = dVdt #* 10**(-4) # m^3/s

# Compute lumen pressure at iterations in t
pressure = calc_pressure_struct_genBC(os.path.join(sim_folder, "AllData"), 1, t) # dynes/cm^2

# Convert pressure from dynes/cm^2 to mmHg
pressure_mmHg = pressure * 0.000750062



# Combine time step, pressure, and volume into one array
PV = np.column_stack((t, pressure_mmHg, vol_cm3))

print('\n## Outputing pressure-volume data and plot ##')

# Write pressure and volume to file
np.savetxt(results_folder + '/../' + "pv.txt", PV, header = 'Timestep Pressure[mmHg] Volume[m^3]')


# Plot pressure vs. volume
fig, ax = plt.subplots()
ax.plot(vol_cm3, pressure_mmHg, linewidth=2.0, marker = 'o')
ax.set_xlabel('Volume [m^3]')
ax.set_ylabel('Pressure [mmHg]')
#plt.xlim([0,0.4])
#plt.ylim([-2, 14])
plt.savefig(os.path.join(sim_folder, 'pv_plot'))

# Plot dVdt vs. time
fig, ax = plt.subplots()
ax.plot(t_dVdt, dVdt_m3, linewidth=2.0, marker = 'o')
ax.set_xlabel('Time')
ax.set_ylabel('dVdt')
#plt.xlim([0,0.4])
#plt.ylim([-2, 14])
plt.savefig(os.path.join(sim_folder, 'dVdt_plot'))

# Plot dVdt vs. Pressure
fig, ax = plt.subplots()
ax.plot(dVdt_m3, pressure_mmHg, linewidth=2.0, marker = 'o')
ax.set_xlabel('dVdt')
ax.set_ylabel('Pressure')
#plt.xlim([0,0.4])
#plt.ylim([-2, 14])
plt.savefig(os.path.join(sim_folder, 'dVdtvsP_plot'))

# Plot 3D and 0D dVdt
fig, ax = plt.subplots()
ax.plot(dVdt_m3, label = 'dVdt_3D')
ax.plot(dVdt_0D, label = 'dVdt_0D', linestyle = '--')
ax.set_xlabel('Timestep')
ax.set_ylabel('dVdt (m^3/s)')
ax.legend()
#plt.xlim([0,0.4])
#plt.ylim([-2, 14])
plt.savefig(os.path.join(sim_folder, 'dVdt3D_vs_dVdt0D'))

# Plot 3D and 0D volume
fig, ax = plt.subplots()
ax.plot(vol_cm3, label = 'V_3D')
ax.plot(vol_0D, label = 'V_0D', linestyle = '--')
ax.set_xlabel('Timestep')
ax.set_ylabel('Volume (m^3)')
ax.legend()
#plt.xlim([0,0.4])
#plt.ylim([-2, 14])
plt.savefig(os.path.join(sim_folder, 'V3D_vs_V0D'))



# ------------------------- Create movie -----------------------------#
# Create Paraview movie .avi file using create_paraview_movie.py
# Do so via a command line call, since we need to run create_paraview_movie.py 
# using pvbatch, not python3
print('\n## Creating Paraview movie ##')
os.system(f"{pvbatch_path} {pvpython_script_path} {results_folder} " +
        f"{start_time} {end_time} {step} " +
        f"{camera_position[0]} {camera_position[1]} {camera_position[2]} " +
        f"{view_up[0]} {view_up[1]} {view_up[2]}")


print('\n## Converting movie.avi to movie.mp4 ##')
# Convert .avi file produced by Paraview Python script to .mp4 file
input_movie_file = os.path.join(sim_folder, 'movie.avi')
output_movie_file = os.path.join(sim_folder, 'movie.mp4')
avi2mp4(input_movie_file, output_movie_file)
