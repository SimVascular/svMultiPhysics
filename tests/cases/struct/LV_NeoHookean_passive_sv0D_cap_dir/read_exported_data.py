#!/usr/bin/env python3
"""
Example script to read binary data exported from svMultiPhysics using Array::write().

The C++ code exports binary files with format:
1. Size (4 bytes, int32)
2. Data (size * sizeof(double) bytes)

Files are named with pattern: {label}_cm.bin

Usage:
    python read_exported_data.py

Example arrays exported:
    - R_neu_1_1_cm.bin: Residual vector
    - Val_neu_1_1_cm.bin: LHS matrix values (CSR format, blocked)
    - Yg_neu_1_1_cm.bin: Velocity/state variables
    - Dg_neu_1_1_cm.bin: Displacement variables
"""

import sys
import os
import glob
import numpy as np
from scipy import sparse
from scipy.io import mmread
import struct

import re

def read_binary_array(filename):
    """
    Read binary array exported by Array::write() or Array3::write().
    
    The file format is:
    - 4 bytes: size (int32)
    - size * 8 bytes: data (double precision)
    
    Args:
        filename: Path to the binary file (e.g., "R_neu_1_1_cm.bin")
    
    Returns:
        numpy array of shape (size,) - needs to be reshaped based on context
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")
    
    with open(filename, 'rb') as f:
        # Read size (4 bytes, int32)
        size_bytes = f.read(4)
        if len(size_bytes) != 4:
            raise ValueError(f"Could not read size from {filename}")
        size = struct.unpack('i', size_bytes)[0]
        
        # Read data (size * 8 bytes for doubles)
        data = np.fromfile(f, dtype=np.float64, count=size)
        
        if len(data) != size:
            raise ValueError(f"Expected {size} elements, but read {len(data)} from {filename}")
        
        return data

def reshape_array_to_2d(data, nrows, ncols):
    """
    Reshape 1D array to 2D using Fortran (column-major) ordering.
    
    Args:
        data: 1D numpy array
        nrows: number of rows
        ncols: number of columns
    
    Returns:
        2D numpy array of shape (nrows, ncols)
    """
    if len(data) != nrows * ncols:
        raise ValueError(f"Data size {len(data)} does not match shape ({nrows}, {ncols})")
    return data.reshape(nrows, ncols, order='F')

def read_residual(filename):
    """
    Read residual vector from text file (Matrix Market export format).
    
    The file format is:
    # Shape: nrows ncols
    value1
    value2
    ...
    
    Data is stored in column-major order (Fortran-style).
    """
    with open(filename, 'r') as f:
        # Read shape from comment
        first_line = f.readline()
        match = re.search(r'Shape:\s+(\d+)\s+(\d+)', first_line)
        if match:
            nrows, ncols = int(match.group(1)), int(match.group(2))
        else:
            raise ValueError(f"Could not parse shape from file {filename}")
        
        # Read data
        data = []
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                data.append(float(line))
        
        # Reshape to column-major order
        residual = np.array(data).reshape(nrows, ncols, order='F')
        return residual

def read_matrix(filename):
    """
    Read sparse matrix from Matrix Market format file.
    
    Returns a scipy.sparse matrix (CSR format).
    """
    matrix = mmread(filename)
    return matrix

def read_residual_matrix(res_file, mat_file):
    """
    Read residual and matrix from given files (text format).
    
    Returns:
        residual: numpy array
        matrix: scipy.sparse matrix
    """
    residual = read_residual(res_file)
    matrix = read_matrix(mat_file)
    return residual, matrix


# Example usage: Read binary files from Array::write()
# Adjust nrows/ncols based on your problem (e.g., tDof=3, tnNo=618 -> nrows=3, ncols=618)

def read_sv0d_binary_data(timestep=1, iteration=1, base_dir='results_svfsiplus', nrows=3, ncols=618):
    """
    Read binary data exported by Array::write() for sv0D simulation.
    
    Args:
        timestep: Time step number
        iteration: Inner iteration number
        base_dir: Base directory for results
        nrows: Number of rows (usually tDof = degrees of freedom per node)
        ncols: Number of columns (usually tnNo = total number of nodes)
    
    Returns:
        dict with keys for different stages: 'R_neu', 'R_as', 'R_solve', etc.
    """
    suffix = f"_{timestep}_{iteration}_cm.bin"
    results = {}
    
    # Define all files to read at different stages
    files_to_read = {
        # After Neumann BC application
        'R_neu': ('R_neu', True),      # Residual (needs reshape)
        'Yg_neu': ('Yg_neu', True),    # Velocity/state (needs reshape)
        'Dg_neu': ('Dg_neu', True),    # Displacement (needs reshape)
        'Val_neu': ('Val_neu', False), # Matrix values (sparse, no reshape)
        
        'Yg_neu0': ('Yg_neu0', True),    # Velocity/state (needs reshape)
        'Dg_neu0': ('Dg_neu0', True),    # Displacement (needs reshape)
        
        'Yg_neu1': ('Yg_neu1', True),    # Velocity/state (needs reshape)
        'Dg_neu1': ('Dg_neu1', True),    # Displacement (needs reshape)
        
        # After assembly
        'R_as': ('R_as', True),        # Residual after assembly
        'Val_as': ('Val_as', False),   # Matrix after assembly
        
        # After solve (if available)
        'R_solve0': ('R_solve0', True),  # Residual after solve
        'Val_solve0': ('Val_solve0', False), # Matrix after solve
        
        'R_solve': ('R_solve', True),  # Residual after solve
        'Val_solve': ('Val_solve', False), # Matrix after solve
        
        # Other intermediate stages
        'Yg_pic': ('Yg_pic', True),    # After predictor step
        'Dg_pic': ('Dg_pic', True),
        'Ag_pic': ('Ag_pic', True),
        
        'Yg_vor_neu': ('Yg_vor_neu', True),  # Before Neumann BC
        'Dg_vor_neu': ('Dg_vor_neu', True),
        
        'Val_bf': ('Val_bf', False),   # After body forces
        'Val_alloc': ('Val_alloc', False), # After allocation
    }
    
    for key, (prefix, needs_reshape) in files_to_read.items():
        filepath = os.path.join(base_dir, f"{prefix}{suffix}")
        print(filepath)
        if os.path.exists(filepath):
            # try:
                data = read_binary_array(filepath)
                
                if needs_reshape and len(data) == nrows * ncols:
                    results[key] = reshape_array_to_2d(data, nrows, ncols)
                    print(f"Read {key}: shape={results[key].shape}")
                else:
                    results[key] = data
                    print(f"Read {key}: size={len(data)}")
            # except Exception as e:
            #     print(f"Warning: Could not read {key}: {e}")
    
    return results

# Example: Read sv0D data
print("\n" + "="*60)
print("=== Reading sv0D binary data ===")
print("="*60)
sv0d_data = read_sv0d_binary_data(timestep=1, iteration=1, base_dir='.', nrows=3, ncols=618)

# Example: Read genBC data for comparison
print("\n" + "="*60)
print("=== Reading genBC binary data ===")
print("="*60)
genbc_base = '/Users/jjv/software/svfork/cap/svMultiPhysics/tests/cases/struct/LV_NeoHookean_passive_genBC_capped/'
genbc_data = read_sv0d_binary_data(timestep=1, iteration=1, base_dir=genbc_base, nrows=3, ncols=618)

# Print summary of what was read
print("\n" + "="*60)
print("=== Summary ===")
print("="*60)
print(f"sv0D files read: {list(sv0d_data.keys())}")
print(f"genBC files read: {list(genbc_data.keys())}")
print(f"Common files: {set(sv0d_data.keys()) & set(genbc_data.keys())}")

# Compare residuals at different stages
print("\n" + "="*60)
print("=== Residual Comparison (sv0D - genBC) ===")
print("="*60)
stages_to_compare = ['R_neu', 'R_as', 'R_solve']

for stage in stages_to_compare:
    if stage in sv0d_data and stage in genbc_data:
        print(f"\n{stage}:")
        diff = sv0d_data[stage] - genbc_data[stage]
        sv0d_norm = np.linalg.norm(sv0d_data[stage])
        genbc_norm = np.linalg.norm(genbc_data[stage])
        rel_diff = np.linalg.norm(diff) / genbc_norm if genbc_norm > 0 else float('inf')
        
        print(f"  Max abs diff:  {np.max(np.abs(diff)):12.6e}")
        print(f"  RMS diff:      {np.sqrt(np.mean(diff**2)):12.6e}")
        print(f"  sv0D norm:     {sv0d_norm:12.6e}")
        print(f"  genBC norm:    {genbc_norm:12.6e}")
        print(f"  Rel diff:      {rel_diff:12.6e}")
        
        # Show location of maximum difference
        if isinstance(diff, np.ndarray) and diff.ndim == 2:
            max_idx = np.unravel_index(np.argmax(np.abs(diff)), diff.shape)
            print(f"  Max diff at:   ({max_idx[0]}, {max_idx[1]})")
            print(f"    sv0D value:  {sv0d_data[stage][max_idx]:12.6e}")
            print(f"    genBC value: {genbc_data[stage][max_idx]:12.6e}")

# Compare matrix values (Val) at different stages
print("\n" + "="*60)
print("=== Matrix Comparison (sv0D - genBC) ===")
print("="*60)
val_stages = ['Val_neu', 'Val_as', 'Val_solve']

for stage in val_stages:
    if stage in sv0d_data and stage in genbc_data:
        print(f"\n{stage}:")
        diff = sv0d_data[stage] - genbc_data[stage]
        sv0d_norm = np.linalg.norm(sv0d_data[stage])
        genbc_norm = np.linalg.norm(genbc_data[stage])
        rel_diff = np.linalg.norm(diff) / genbc_norm if genbc_norm > 0 else float('inf')
        
        print(f"  Max abs diff:  {np.max(np.abs(diff)):12.6e}")
        print(f"  RMS diff:      {np.sqrt(np.mean(diff**2)):12.6e}")
        print(f"  sv0D norm:     {sv0d_norm:12.6e}")
        print(f"  genBC norm:    {genbc_norm:12.6e}")
        print(f"  Rel diff:      {rel_diff:12.6e}")
