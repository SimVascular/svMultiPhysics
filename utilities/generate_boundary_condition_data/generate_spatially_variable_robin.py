#!/usr/bin/env python3
"""
Script to generate spatially varying Robin boundary condition VTP files.

This script reads a VTP file and creates a new VTP file with spatially varying
stiffness and damping values based on user-defined functions of the coordinates (x, y, z).

To use this script:
1. Edit the configuration section below to set your input/output files and functions
2. Run: python generate_spatially_variable_robin.py

The script supports mathematical expressions using x, y, z coordinates and common
mathematical functions (sin, cos, exp, sqrt, etc.).

Requirements:
    - numpy
    - pyvista
"""

import numpy as np
import pyvista as pv
import os
import math


def safe_eval(expression: str, x: float, y: float, z: float) -> float:
    """
    Safely evaluate a mathematical expression with x, y, z variables.
    
    Args:
        expression: Mathematical expression as string
        x, y, z: Coordinate values
        
    Returns:
        Evaluated result
        
    Raises:
        ValueError: If expression is invalid or contains unsafe operations
    """
    # Define safe namespace for evaluation
    safe_dict = {
        'x': x, 'y': y, 'z': z,
        'sin': math.sin, 'cos': math.cos, 'tan': math.tan,
        'exp': math.exp, 'log': math.log, 'log10': math.log10,
        'sqrt': math.sqrt, 'abs': abs, 'pow': pow,
        'pi': math.pi, 'e': math.e,
        'min': min, 'max': max,
        '__builtins__': {}
    }
    
    try:
        return float(eval(expression, safe_dict))
    except Exception as e:
        raise ValueError(f"Error evaluating expression '{expression}': {e}")


def read_vtp_file(filepath: str) -> pv.PolyData:
    """
    Read a VTP file and return the PyVista polydata object.
    
    Args:
        filepath: Path to the VTP file
        
    Returns:
        PyVista polydata object
        
    Raises:
        FileNotFoundError: If file doesn't exist
        RuntimeError: If file cannot be read
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"VTP file not found: {filepath}")
    
    try:
        return pv.read(filepath)
    except Exception as e:
        raise RuntimeError(f"Failed to read VTP file '{filepath}': {e}")


def write_vtp_file(polydata: pv.PolyData, filepath: str) -> None:
    """
    Write a polydata object to a VTP file.
    
    Args:
        polydata: PyVista polydata object
        filepath: Output file path
        
    Raises:
        RuntimeError: If file cannot be written
    """
    try:
        polydata.save(filepath, binary=True)
    except Exception as e:
        raise RuntimeError(f"Failed to write VTP file '{filepath}': {e}")


def get_coordinates(polydata: pv.PolyData) -> np.ndarray:
    """
    Extract coordinates from polydata.
    
    Args:
        polydata: PyVista polydata object
        
    Returns:
        Array of shape (n_points, 3) with x, y, z coordinates
    """
    if polydata.n_points == 0:
        raise ValueError("No points found in VTP file")
    
    return polydata.points


def add_point_data(polydata: pv.PolyData, name: str, data: np.ndarray) -> None:
    """
    Add point data array to polydata.
    
    Args:
        polydata: PyVista polydata object
        name: Name of the data array
        data: 1D numpy array of data values
    """
    polydata[name] = data


def generate_spatially_varying_robin_bc(
    input_vtp: str,
    output_vtp: str,
    stiffness_func: str,
    damping_func: str,
    stiffness_scale: float = 1.0,
    damping_scale: float = 1.0,
    min_stiffness: float = 0.0,
    min_damping: float = 0.0,
    verbose: bool = False
) -> None:
    """
    Generate spatially varying Robin boundary condition VTP file.
    
    Args:
        input_vtp: Path to input VTP file
        output_vtp: Path to output VTP file
        stiffness_func: Mathematical expression for stiffness as function of x, y, z
        damping_func: Mathematical expression for damping as function of x, y, z
        stiffness_scale: Scaling factor for stiffness values
        damping_scale: Scaling factor for damping values
        min_stiffness: Minimum allowed stiffness value
        min_damping: Minimum allowed damping value
        verbose: Print detailed information
    """
    if verbose:
        print(f"Reading VTP file: {input_vtp}")
    
    # Read input VTP file
    polydata = read_vtp_file(input_vtp)
    coords = get_coordinates(polydata)
    n_points = polydata.n_points
    
    if verbose:
        print(f"Found {n_points} points in the mesh")
        print(f"Coordinate ranges:")
        print(f"  X: [{coords[:, 0].min():.6f}, {coords[:, 0].max():.6f}]")
        print(f"  Y: [{coords[:, 1].min():.6f}, {coords[:, 1].max():.6f}]")
        print(f"  Z: [{coords[:, 2].min():.6f}, {coords[:, 2].max():.6f}]")
    
    # Initialize arrays for stiffness and damping
    stiffness_values = np.zeros(n_points)
    damping_values = np.zeros(n_points)
    
    # Evaluate functions at each point
    if verbose:
        print("Evaluating stiffness and damping functions...")
    
    for i in range(n_points):
        x, y, z = coords[i, :]
        
        try:
            # Evaluate stiffness function
            stiffness_raw = safe_eval(stiffness_func, x, y, z)
            stiffness_values[i] = max(min_stiffness, stiffness_scale * stiffness_raw)
            
            # Evaluate damping function
            damping_raw = safe_eval(damping_func, x, y, z)
            damping_values[i] = max(min_damping, damping_scale * damping_raw)
            
        except ValueError as e:
            raise ValueError(f"Error at point {i} (x={x:.6f}, y={y:.6f}, z={z:.6f}): {e}")
    
    if verbose:
        print(f"Stiffness range: [{stiffness_values.min():.6e}, {stiffness_values.max():.6e}]")
        print(f"Damping range: [{damping_values.min():.6e}, {damping_values.max():.6e}]")
    
    # Create output polydata (copy of input)
    output_polydata = polydata.copy()
    
    # Add stiffness and damping arrays
    add_point_data(output_polydata, "Stiffness", stiffness_values)
    add_point_data(output_polydata, "Damping", damping_values)
    
    # Show mesh with Stiffness
    output_polydata.plot(scalars="Stiffness", show_edges=True, cmap="viridis")
    
    # Write output file
    if verbose:
        print(f"Writing output VTP file: {output_vtp}")
    
    write_vtp_file(output_polydata, output_vtp)
    
    if verbose:
        print("Successfully generated spatially varying Robin BC VTP file!")


# =============================================================================
# CONFIGURATION SECTION - EDIT THESE VALUES
# =============================================================================

# Input and output file paths
os.chdir(os.path.dirname(os.path.realpath(__file__)))
INPUT_VTP_FILE = "mesh/mesh-surfaces/Y0.vtp"  # Path to input VTP file
OUTPUT_VTP_FILE = "Y0_spatially_varying_robin.vtp"  # Path to output VTP file

# Mathematical expressions for stiffness and damping as functions of x, y, z
# Available functions: sin, cos, tan, exp, log, log10, sqrt, abs, pow, min, max
# Available constants: pi, e
STIFFNESS_FUNCTION = "100 * z"  # Example: linear variation in z-direction
DAMPING_FUNCTION = "0"      # Example: no damping

# Scaling factors
STIFFNESS_SCALE = 1.0  # Scaling factor for stiffness values
DAMPING_SCALE = 1.0    # Scaling factor for damping values

# Minimum values (stiffness and damping cannot be negative)
MIN_STIFFNESS = 0.0
MIN_DAMPING = 0.0

# Verbose output
VERBOSE = True

# =============================================================================
# EXAMPLE FUNCTIONS - UNCOMMENT AND MODIFY AS NEEDED
# =============================================================================

# Linear variation in x-direction
# STIFFNESS_FUNCTION = "1000 * x"
# DAMPING_FUNCTION = "10 * x"

# Quadratic variation from center
# STIFFNESS_FUNCTION = "1000 * (x**2 + y**2 + z**2)"
# DAMPING_FUNCTION = "10"

# Sinusoidal variation
# STIFFNESS_FUNCTION = "1000 * (1 + sin(2*pi*x))"
# DAMPING_FUNCTION = "10 * cos(pi*y)"

# Exponential decay
# STIFFNESS_FUNCTION = "1000 * exp(-x)"
# DAMPING_FUNCTION = "10 * exp(-0.5*x)"

# Step function
# STIFFNESS_FUNCTION = "1000 * (x > 2.5)"
# DAMPING_FUNCTION = "10 * (y > 2.5)"

# Gaussian-like distribution
# STIFFNESS_FUNCTION = "1000 * exp(-(x**2 + y**2 + z**2))"
# DAMPING_FUNCTION = "10"

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    try:
        print("Generating spatially varying Robin boundary condition VTP file...")
        print(f"Input file: {INPUT_VTP_FILE}")
        print(f"Output file: {OUTPUT_VTP_FILE}")
        print(f"Stiffness function: {STIFFNESS_FUNCTION}")
        print(f"Damping function: {DAMPING_FUNCTION}")
        print()
        
        generate_spatially_varying_robin_bc(
            input_vtp=INPUT_VTP_FILE,
            output_vtp=OUTPUT_VTP_FILE,
            stiffness_func=STIFFNESS_FUNCTION,
            damping_func=DAMPING_FUNCTION,
            stiffness_scale=STIFFNESS_SCALE,
            damping_scale=DAMPING_SCALE,
            min_stiffness=MIN_STIFFNESS,
            min_damping=MIN_DAMPING,
            verbose=VERBOSE
        )
        
        print("✓ Successfully generated spatially varying Robin BC VTP file!")
        
    except Exception as e:
        print(f"✗ Error: {e}")
        exit(1)
