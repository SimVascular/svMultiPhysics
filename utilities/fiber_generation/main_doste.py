#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
# SPDX-License-Identifier: BSD-3-Clause
# -*-coding:utf-8 -*-
"""Main script for generating biventricular fibers using the Doste method.

This module implements fiber generation for biventricular heart models using
the Laplace-Dirichlet rule-based method described in:
Doste et al. 2019, "A rule-based method to model myocardial fiber orientation
in cardiac biventricular geometries with outflow tracts"
https://doi.org/10.1002/cnm.3185

The script supports command-line arguments for customization of mesh paths,
output directories, and solver executables.
"""

import os
import argparse
import pyvista as pv
from fiber_generation.laplace_solver import LaplaceSolver
from fiber_generation.fiber_generator import FibGenDoste
from fiber_generation.surface_names import SurfaceName
from fiber_generation.surface_utils import generate_epi_apex
from time import time


if __name__ == "__main__":
    ###########################################################
    ############  USER INPUTS  ################################
    ###########################################################

    run_flag = True
    svmultiphysics_exec = "svmultiphysics "

    mesh_path = "example/biv_with_outflow_tracts/mesh-complete.mesh.vtu"
    outdir = "example/biv_with_outflow_tracts/output_doste"
    surfaces_dir = 'example/biv_with_outflow_tracts/mesh-surfaces'


    # Parameters from the Doste paper https://doi.org/10.1002/cnm.3185
    params = {
        'AENDORV': 90,
        'AEPIRV': -25,
        'AENDOLV': 60,
        'AEPILV': -60,

        'AOTENDOLV': 90, 
        'AOTENDORV': 90,
        'AOTEPILV': 0,
        'AOTEPIRV': 0,

        'BENDORV': 0,
        'BEPIRV': 20,
        'BENDOLV': -20,
        'BEPILV': 20,
    }

    ###########################################################
    ############  FIBER GENERATION  ###########################
    ###########################################################

    # Optional CLI overrides
    parser = argparse.ArgumentParser(description="Generate fibers using the Doste method.")
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
    surface_paths = {
        SurfaceName.EPICARDIUM: f'{surfaces_dir}/epi.vtp',
        SurfaceName.EPICARDIUM_APEX: f'{surfaces_dir}/epi_apex.vtp',
        SurfaceName.AORTIC_VALVE: f'{surfaces_dir}/av.vtp',
        SurfaceName.MITRAL_VALVE: f'{surfaces_dir}/mv.vtp',
        SurfaceName.TRICUSPID_VALVE: f'{surfaces_dir}/tv.vtp',
        SurfaceName.PULMONARY_VALVE: f'{surfaces_dir}/pv.vtp',
        SurfaceName.ENDOCARDIUM_LV: f'{surfaces_dir}/endo_lv.vtp',
        SurfaceName.ENDOCARDIUM_RV: f'{surfaces_dir}/endo_rv.vtp',
        SurfaceName.BASE: f'{surfaces_dir}/top.vtp'
    }

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
        laplace_results_file = solver.run("doste", outdir)
    else:
        laplace_results_file = os.path.join(outdir, 'result_001.vtu')

    # Initialize fiber generator
    print("\nGenerating fibers using Doste method...")
    fib_gen = FibGenDoste()

    # Load Laplace results
    fib_gen.load_laplace_results(laplace_results_file)

    # Generate fiber directions
    F, S, T = fib_gen.generate_fibers(params)
    print(f"generate fibers (Doste method) elapsed time: {time() - start:.3f} s")

    # Write fibers to output directory
    fib_gen.write_fibers(outdir)

    # Save the result mesh
    result_mesh_path = os.path.join(outdir, "results_doste.vtu")
    fib_gen.mesh.save(result_mesh_path)
    print(f"\nResults saved to: {result_mesh_path}")
