# Python script to calculate lumen volume and fluid pressure on endo surface from 
# results of svFSI struct
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

import os # for checking and creating directories, and loading modules
import matplotlib.pyplot as plt
import numpy as np
import sys

# Import custom functions useful for post-processing
from post_processing_functions import *

## -------- PARAMETERS TO CHANGE -------------------- ##

# Simulation folder file path
sim_folder = os.path.dirname(os.path.realpath(__file__))

# Optionally read sim_folder as command line argument
if len(sys.argv) == 2:
    print(sys.argv)
    sim_folder = os.path.abspath(sys.argv[1])
    print(sim_folder)


# Undeformed surface, corresponding to deformed surface on result.vtu of which we want to compute the volume
reference_surface = os.path.join(sim_folder, 'mesh/mesh-surfaces/endo.vtp')

# Undeformed vtu file (the original volume mesh)
reference_volume = os.path.join(sim_folder, 'mesh/mesh-complete.mesh.vtu')

# svFSI results folder, containing results.vtu
#results_folder = os.path.join(sim_folder, 'results_svfsi/')
results_folder = os.path.join(sim_folder, '6-procs/')

# File containing genBC output
#alldata_file = 'AllData_svfsi_vvedula22'
alldata_file = 'AllData'

# Input file path
input_file = os.path.join(sim_folder, 'svFSI.inp')

# Pressure file (unused)
pressure_dat_file = os.path.join(sim_folder, 'pressure.dat')


## -------------------- END PARAMETERS TO CHANGE ------------------------ ## 


# Automatically determine the start time, end time, and step size based on all
# results file in results_folder
print(results_folder)
(start_time, end_time, step) = get_start_end_step(results_folder)
# Option to manually set start time, end time, and time step of results files to process
#start_time = 5
#end_time = 85
#step = 5

# Compute lumen volume from simulation results
(t_3D, vol_3D) = calc_volume_struct(start_time, end_time, step, results_folder, reference_surface)
vol_3D_cm3 = np.array(vol_3D) * (100)**3 # cm^3

# Compute lumen dVdt from simulation results
(t_dVdt_3D, dVdt_3D) = calc_dVdt_struct(start_time, end_time, step, results_folder, reference_surface)
# Convert volume to cm^3/s
dVdt_3D = np.array(dVdt_3D) # cm/s * m^2
dVdt_3D_cm3 = dVdt_3D * (100)**2 # cm^3/s

# Compute lumen volume from AllData file
vol_0D_cm3 = calc_volume_struct_genBC(os.path.join(sim_folder, alldata_file), 2, t_3D)
# Add on initial volume
vol_0D_cm3 += vol_3D_cm3[0] # cm^3

# Compute flow rate = dVdt from AllData file
dVdt_0D_cm3 = calc_volume_struct_genBC(os.path.join(sim_folder, alldata_file), 4, t_dVdt_3D)


# Compute lumen pressure at iterations in t
pressure = calc_pressure_struct_genBC(os.path.join(sim_folder, "AllData"), 1, t_3D) # dynes/cm^2

# Convert pressure from dynes/cm^2 to mmHg
#pressure_mmHg = pressure * 0.000750062



# Combine time step, pressure, and volume into one array
#PV = np.column_stack((t_3D, pressure, vol_3D_cm3))

#print('\n## Outputing pressure-volume data and plot ##')

# Write pressure and volume to file
#np.savetxt(results_folder + '/../' + "pv.txt", PV, header = 'Timestep Pressure[mmHg] Volume[m^3]')


# Plot pressure vs. volume
fig, ax = plt.subplots()
ax.plot(vol_3D_cm3, pressure, linewidth=2.0, marker = 'o')
ax.set_xlabel('Volume [cm^3]')
ax.set_ylabel('Pressure [dyne/cm^2]')
#plt.xlim([0,0.4])
#plt.ylim([-2, 14])
plt.savefig(os.path.join(sim_folder, 'pv_plot'))

# Plot dVdt vs. time
#fig, ax = plt.subplots()
#ax.plot(t_dVdt, dVdt_m3, linewidth=2.0, marker = 'o')
#ax.set_xlabel('Time')
#ax.set_ylabel('dVdt')
#plt.xlim([0,0.4])
#plt.ylim([-2, 14])
#plt.savefig(os.path.join(sim_folder, 'dVdt_plot'))

# Plot dVdt vs. Pressure
#fig, ax = plt.subplots()
#ax.plot(dVdt_m3, pressure_mmHg, linewidth=2.0, marker = 'o')
#ax.set_xlabel('dVdt')
#ax.set_ylabel('Pressure')
#plt.xlim([0,0.4])
#plt.ylim([-2, 14])
#plt.savefig(os.path.join(sim_folder, 'dVdtvsP_plot'))

# Plot 3D and 0D dVdt
fig, ax = plt.subplots()
ax.plot(dVdt_3D_cm3, label = 'dVdt_3D')
ax.plot(dVdt_0D_cm3, label = 'dVdt_0D', linestyle = '--')
ax.set_xlabel('Timestep')
ax.set_ylabel('dVdt (cm^3/s)')
ax.legend()
#plt.xlim([0,0.4])
#plt.ylim([-2, 14])
plt.savefig(os.path.join(sim_folder, 'dVdt3D_vs_dVdt0D'))

# Plot 3D and 0D volume
fig, ax = plt.subplots()
ax.plot(vol_3D_cm3, label = 'V_3D')
ax.plot(vol_0D_cm3, label = 'V_0D', linestyle = '--')
ax.set_xlabel('Timestep')
ax.set_ylabel('Volume (cm^3)')
ax.legend()
#plt.xlim([0,0.4])
#plt.ylim([-2, 14])
plt.savefig(os.path.join(sim_folder, 'V3D_vs_V0D'))
