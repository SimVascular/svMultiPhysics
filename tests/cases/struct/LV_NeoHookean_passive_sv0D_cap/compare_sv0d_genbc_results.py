#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
Created on 2026/02/09 14:36:55

@author: Javiera Jilberto Vallejos 
'''

import numpy as np
import meshio as io

genbc_path = 'genbc_results'
sv0d_path = 'results'

ts = 1
import matplotlib.pyplot as plt

results_genbc = io.read(f'{genbc_path}/result_{ts:03d}.vtu')
results_sv0d = io.read(f'{sv0d_path}/result_{ts:03d}.vtu')

variables = ['Displacement', 'Velocity']
labels = ['Displacement (m)', 'Velocity (m/s)']

# Initialize dictionaries to store differences and mean values over time
point_data_diffs = {}
point_data_sv0d = {}
point_data_genbc = {}
cell_data_diffs = {}
cell_data_sv0d = {}
cell_data_genbc = {}

# Loop through timesteps
for ts in range(1, 201):
    results_genbc = io.read(f'{genbc_path}/result_{ts:03d}.vtu')
    results_sv0d = io.read(f'{sv0d_path}/result_{ts:03d}.vtu')
    
    # Process point data
    for key in results_sv0d.point_data.keys():
        if key in results_genbc.point_data and key in variables:
            diff = np.abs(results_sv0d.point_data[key] - results_genbc.point_data[key])
            if key not in point_data_diffs:
                point_data_diffs[key] = []
                point_data_sv0d[key] = []
                point_data_genbc[key] = []
            point_data_diffs[key].append(np.mean(diff))
            point_data_sv0d[key].append(np.mean(results_sv0d.point_data[key]))
            point_data_genbc[key].append(np.mean(results_genbc.point_data[key]))
    
    # Process cell data
    for key in results_sv0d.cell_data.keys():
        if key in results_genbc.cell_data and key in variables:
            diff = np.abs(results_sv0d.cell_data[key][0] - results_genbc.cell_data[key][0])
            if key not in cell_data_diffs:
                cell_data_diffs[key] = []
                cell_data_sv0d[key] = []
                cell_data_genbc[key] = []
            cell_data_diffs[key].append(np.mean(diff))
            cell_data_sv0d[key].append(np.mean(results_sv0d.cell_data[key][0]))
            cell_data_genbc[key].append(np.mean(results_genbc.cell_data[key][0]))

# Convert to numpy arrays
for key in point_data_diffs:
    point_data_diffs[key] = np.array(point_data_diffs[key])
    point_data_sv0d[key] = np.array(point_data_sv0d[key])
    point_data_genbc[key] = np.array(point_data_genbc[key])
    
for key in cell_data_diffs:
    cell_data_diffs[key] = np.array(cell_data_diffs[key])
    cell_data_sv0d[key] = np.array(cell_data_sv0d[key])
    cell_data_genbc[key] = np.array(cell_data_genbc[key])

# Plot traces
fig, axes = plt.subplots(len(point_data_diffs) + len(cell_data_diffs), 1, figsize=(6, 6), sharex=True)
if len(axes.shape) == 0:
    axes = [axes]

idx = 0
for key in point_data_sv0d:
    axes[idx].plot(range(1, 201), point_data_sv0d[key], label='sv0d')
    axes[idx].plot(range(1, 201), point_data_genbc[key], label='genbc')
    axes[idx].set_ylabel(f'{labels[variables.index(key)]}')
    idx += 1

for key in cell_data_sv0d:
    axes[idx].plot(range(1, 201), cell_data_sv0d[key], label='sv0d')
    axes[idx].plot(range(1, 201), cell_data_genbc[key], label='genbc')
    axes[idx].set_ylabel(f'{labels[variables.index(key)]}')
    idx += 1

axes[-1].set_xlabel('Timestep')
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right')
plt.tight_layout()
plt.savefig('comparison.png')
plt.show()
