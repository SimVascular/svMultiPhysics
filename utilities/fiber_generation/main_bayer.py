#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
# SPDX-License-Identifier: BSD-3-Clause
# -*-coding:utf-8 -*-
"""Main script for generating biventricular fibers using the Bayer method.

This module implements fiber generation for biventricular heart models using
the Laplace-Dirichlet rule-based method described in:
Bayer et al. 2012, "A Novel Rule-Based Algorithm for Assigning Myocardial 
Fiber Orientation to Computational Heart Models"
https://doi.org/10.1007/s10439-012-0593-5

The script supports command-line arguments for customization of mesh paths,
output directories, and solver executables.
"""

import argparse
import os
import pyvista as pv
from fiber_generation.laplace_solver import LaplaceSolver
from fiber_generation.fiber_generator import FibGenBayer
from fiber_generation.surface_names import SurfaceName
from fiber_generation.surface_utils import generate_epi_apex
from time import time


if __name__ == "__main__":

    ###########################################################
    ############  USER INPUTS  ################################
    ###########################################################

    run_flag = True
    svmultiphysics_exec = "svmultiphysics "

    mesh_path = "example/biv_truncated/VOLUME.vtu"
    outdir = "example/biv_truncated/output_bayer"
    surfaces_dir = 'example/biv_truncated/mesh-surfaces'
                    
    # Parameters for the Bayer et al. method https://doi.org/10.1007/s10439-012-0593-5. 
    params = {
        "ALFA_END": 60.0,
        "ALFA_EPI": -60.0,
        "BETA_END": -20.0,
        "BETA_EPI": 20.0,
    }


    ###########################################################
    ############  FIBER GENERATION  ###########################
    ###########################################################

    # Optional CLI overrides
    parser = argparse.ArgumentParser(description="Generate fibers using the Bayer method.")
    parser.add_argument("--svmultiphysics-exec", default=svmultiphysics_exec, help="svMultiPhysics executable/command (default: %(default)s)")
    parser.add_argument("--mesh-path", default=mesh_path, help="Path to the volumetric mesh .vtu (default: %(default)s)")
    parser.add_argument(
        "--surfaces-dir",
        default=surfaces_dir,
        help="Directory containing mesh surfaces; default: <parent of mesh_path>/mesh-surfaces",
    )
    parser.add_argument("--outdir", default=outdir, help="Output directory (default: %(default)s)")
    args = parser.parse_args()

    svmultiphysics_exec = args.svmultiphysics_exec
    if not svmultiphysics_exec.endswith(" "):
        svmultiphysics_exec = svmultiphysics_exec + " "

    mesh_path = args.mesh_path
    outdir = args.outdir

    if args.surfaces_dir is None:
        surfaces_dir = os.path.join(os.path.dirname(mesh_path), "mesh-surfaces")
    else:
        surfaces_dir = os.path.abspath(args.surfaces_dir)

    # Make sure the paths are full paths
    mesh_path = os.path.abspath(mesh_path)
    outdir = os.path.abspath(outdir)
    surfaces_dir = os.path.abspath(surfaces_dir)

    # Define surface paths
    surface_paths = {SurfaceName.EPICARDIUM: f'{surfaces_dir}/EPI.vtp',
                    SurfaceName.EPICARDIUM_APEX: f'{surfaces_dir}/EPI_APEX.vtp',
                    SurfaceName.BASE: f'{surfaces_dir}/BASE.vtp',
                    SurfaceName.ENDOCARDIUM_LV: f'{surfaces_dir}/LV.vtp',
                    SurfaceName.ENDOCARDIUM_RV: f'{surfaces_dir}/RV.vtp'}
    
    # Create output directory if needed
    os.makedirs(outdir, exist_ok=True)
    
    # Check if the EPICARDIUM_APEX surface exists; if not create it
    start = time()
    if not os.path.exists(surface_paths[SurfaceName.EPICARDIUM_APEX]):
        print("Generating EPICARDIUM_APEX surface...")
        generate_epi_apex(surface_paths)
        
    # Initialize Laplace solver
    solver = LaplaceSolver(mesh_path, surface_paths, svmultiphysics_exec)

    # Run the Laplace solver
    if run_flag:
        print("Running Laplace solver...")
        laplace_results_file = solver.run("bayer", outdir)
    else:
        laplace_results_file = os.path.join(outdir, 'result_001.vtu')

    # Initialize fiber generator
    print("\nGenerating fibers using Bayer method...")
    fib_gen = FibGenBayer()

    # Load Laplace results
    fib_gen.load_laplace_results(laplace_results_file)

    # Generate fiber directions
    F, S, T = fib_gen.generate_fibers(params)
    print(f"generate fibers (Bayer method) elapsed time: {time() - start:.3f} s")
    
    # Write fibers to output directory
    fib_gen.write_fibers(outdir)

    # Save the result mesh
    result_mesh_path = os.path.join(outdir, "results_bayer.vtu")
    fib_gen.mesh.save(result_mesh_path)
    print(f"\nResults saved to: {result_mesh_path}")
