#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
# SPDX-License-Identifier: BSD-3-Clause
# -*-coding:utf-8 -*-

import os
import numpy as np
import pyvista as pv

from fiber_generation.fiber_generator import FibGenDoste
from scipy import stats
import matplotlib.pyplot as plt


def project_to_plane(fiber, vector1, vector2):
    """
    Project fiber vectors onto the plane formed by vector1 and vector2.
    
    Parameters
    ----------
    fiber : ndarray (n, 3)
        Array of fiber vectors to project
    vector1 : ndarray (n, 3)
        First vector defining the plane for each row
    vector2 : ndarray (n, 3)
        Second vector defining the plane for each row
        
    Returns
    -------
    projected : ndarray (n, 3)
        Fiber vectors projected onto the plane formed by vector1 and vector2
    """
    # Calculate the normal to the plane (cross product of vector1 and vector2)
    normal = np.cross(vector1, vector2)
    
    # Normalize the normal vectors
    normal_norm = np.linalg.norm(normal, axis=1, keepdims=True)
    normal_normalized = normal / normal_norm
    
    # Project fiber onto the plane by removing the component along the normal
    # projected = fiber - (fiber · normal) * normal
    dot_product = np.sum(fiber * normal_normalized, axis=1, keepdims=True)
    projected = fiber - dot_product * normal_normalized
    
    # Normalize the projected vectors
    projected_norm = np.linalg.norm(projected, axis=1, keepdims=True)
    projected_normalized = projected / projected_norm
    
    return projected_normalized

def calculate_alpha_beta_angles(f, eC, eL, eT):
    # Project fiber to the plane formed by eL and eC
    f_projected = project_to_plane(f, eL, eC)

    beta_dot = np.abs(np.sum(f_projected * f, axis=1))
    beta_dot = np.clip(beta_dot, 0, 1)  # Ensure values are within valid range
    abs_beta_angle = np.rad2deg(np.arccos(beta_dot))
    sign_beta = -np.sign(np.sum(f * eT, axis=1))
    beta_angle = abs_beta_angle * sign_beta

    alpha_dot = np.abs(np.sum(eC * f_projected, axis=1))
    alpha_dot = np.clip(alpha_dot, 0, 1)  # Ensure values are within valid range
    abs_alpha_angle = np.rad2deg(np.arccos(alpha_dot))
    sign_alpha = np.sign(np.sum(f_projected * eL, axis=1))
    alpha_angle = abs_alpha_angle * sign_alpha
    
    return alpha_angle, beta_angle, f_projected



if __name__ == "__main__":

    outdir = "example/biv_with_outflow_tracts/output_doste"
    save_vtu = True

    params_zero = {
        'AENDORV': 0,
        'AEPIRV': 0,
        'AENDOLV': 0,
        'AEPILV': 0,

        'AOTENDOLV': 0, 
        'AOTENDORV': 0,
        'AOTEPILV': 0,
        'AOTEPIRV': 0,

        'BENDORV': 0,
        'BEPIRV': 0,
        'BENDOLV': 0,
        'BEPILV': 0,
    }

    params_alpha = {
        'AENDORV': 90,
        'AEPIRV': -25,
        'AENDOLV': 60,
        'AEPILV': -60,

        'AOTENDOLV': 90, 
        'AOTENDORV': 90,
        'AOTEPILV': 0,
        'AOTEPIRV': 0,

        'BENDORV': 0,
        'BEPIRV': 0,
        'BENDOLV': 0,
        'BEPILV': 0,
    }

    params_beta = {
        'AENDORV': 0,
        'AEPIRV': 0,
        'AENDOLV': 0,
        'AEPILV': 0,

        'AOTENDOLV': 0, 
        'AOTENDORV': 0,
        'AOTEPILV': 0,
        'AOTEPIRV': 0,

        'BENDORV': 0,
        'BEPIRV': 20,
        'BENDOLV': -20,
        'BEPILV': 20,
    }

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

    # Read laplace solutions
    laplace_results_file = os.path.join(outdir, 'result_001.vtu')

    # Initialize fiber generator
    fib_gen = FibGenDoste()
    fib_gen.load_laplace_results(laplace_results_file)

    # Calculate orthogonal basis vectors
    eC, eL, eT = fib_gen.generate_fibers(params_zero)

    # Sanity check 1: Only alpha rotation
    f_alpha, n_alpha, s_alpha = fib_gen.generate_fibers(params_alpha)
    ref_alpha_only_a, ref_beta_only_a = fib_gen.get_angle_fields(params_alpha)
    alpha_only_a, beta_only_a, f_projected = calculate_alpha_beta_angles(f_alpha, eC, eL, eT)

    if save_vtu:
        fib_gen.mesh.cell_data.clear()
        fib_gen.mesh.cell_data['f'] = f_alpha
        fib_gen.mesh.cell_data['s'] = s_alpha
        fib_gen.mesh.cell_data['n'] = n_alpha
        fib_gen.mesh.cell_data['alpha_only_a'] = alpha_only_a
        fib_gen.mesh.cell_data['beta_only_a'] = beta_only_a
        fib_gen.mesh.cell_data['alpha_ref_a'] = ref_alpha_only_a
        fib_gen.mesh.cell_data['beta_ref_a'] = ref_beta_only_a
        fib_gen.mesh.cell_data['diff_alpha_a'] = alpha_only_a - ref_alpha_only_a
        fib_gen.mesh.cell_data['diff_beta_a'] = beta_only_a - ref_beta_only_a
        fib_gen.mesh.save('example/biv_with_outflow_tracts/validation_doste_onlyalpha.vtu')

    # Sanity check 2: Only beta rotation
    f_beta, n_beta, s_beta = fib_gen.generate_fibers(params_beta)
    ref_alpha_only_b, ref_beta_only_b = fib_gen.get_angle_fields(params_beta)
    alpha_only_b, beta_only_b, f_projected = calculate_alpha_beta_angles(f_beta, eC, eL, eT)

    if save_vtu:
        fib_gen.mesh.cell_data.clear()
        fib_gen.mesh.cell_data['f'] = f_beta
        fib_gen.mesh.cell_data['s'] = s_beta
        fib_gen.mesh.cell_data['n'] = n_beta
        fib_gen.mesh.cell_data['alpha_only_b'] = alpha_only_b
        fib_gen.mesh.cell_data['beta_only_b'] = beta_only_b
        fib_gen.mesh.cell_data['alpha_ref_b'] = ref_alpha_only_b
        fib_gen.mesh.cell_data['beta_ref_b'] = ref_beta_only_b
        fib_gen.mesh.cell_data['diff_alpha_b'] = alpha_only_b - ref_alpha_only_b
        fib_gen.mesh.cell_data['diff_beta_b'] = beta_only_b - ref_beta_only_b
        fib_gen.mesh.save('example/biv_with_outflow_tracts/validation_doste_onlybeta.vtu')
    
    # Alpha and beta rotation combined
    eC, eL, eT = fib_gen.generate_fibers(params_zero)
    f_combined, n_combined, s_combined = fib_gen.generate_fibers(params)
    ref_alpha_combined, ref_beta_combined = fib_gen.get_angle_fields(params)
    alpha_combined, beta_combined, f_projected = calculate_alpha_beta_angles(f_combined, eC, eL, eT)

    if save_vtu:
        fib_gen.mesh.cell_data.clear()
        fib_gen.mesh.cell_data['f'] = f_combined
        fib_gen.mesh.cell_data['s'] = s_combined
        fib_gen.mesh.cell_data['n'] = n_combined
        fib_gen.mesh.cell_data['alpha_combined'] = alpha_combined
        fib_gen.mesh.cell_data['beta_combined'] = beta_combined
        fib_gen.mesh.cell_data['alpha_ref'] = ref_alpha_combined
        fib_gen.mesh.cell_data['beta_ref'] = ref_beta_combined
        fib_gen.mesh.cell_data['diff_alpha'] = alpha_combined - ref_alpha_combined
        fib_gen.mesh.cell_data['diff_beta'] = beta_combined - ref_beta_combined
        fib_gen.mesh.save('example/biv_with_outflow_tracts/validation_doste_combined.vtu')
        
    # Create figure with correlation plots
    fig, axes = plt.subplots(1, 3, figsize=(8, 3.5), constrained_layout=True)
    fig.suptitle(r'$\alpha$ and $\beta$ angle correlations', fontsize=16)

    # Alpha Only
    axes[0].scatter(ref_alpha_only_a, alpha_only_a, alpha=0.02, s=10, color='blue')
    axes[0].scatter(ref_beta_only_a, beta_only_a, alpha=0.02, s=10, color='red')
    
    slope_a, intercept_a, r_value_a, _, _ = stats.linregress(ref_alpha_only_a, alpha_only_a)
    x_line = np.array([-90, 90])
    axes[0].plot(x_line, slope_a * x_line + intercept_a, 'b-', lw=1, label=f'α: R²={r_value_a**2:.3f}, m={slope_a:.3f}')
    
    axes[0].plot([-90, 90], [-90, 90], 'k--', lw=1, alpha=0.5)
    axes[0].set_title(r'$\alpha=\pm 60, \beta=0$')
    axes[0].set_xlabel('Scalar Interpolation (degrees)')
    axes[0].set_ylabel('Fiber angles (degrees)')
    axes[0].legend(fontsize=8, loc='upper left')

    # Beta Only
    axes[1].scatter(ref_alpha_only_b, alpha_only_b, alpha=0.02, s=10, color='blue')
    axes[1].scatter(ref_beta_only_b, beta_only_b, alpha=0.02, s=10, color='red')
    
    slope_b, intercept_b, r_value_b, _, _ = stats.linregress(ref_beta_only_b, beta_only_b)
    axes[1].plot(x_line, slope_b * x_line + intercept_b, 'r-', lw=1, label=f'β: R²={r_value_b**2:.3f}, m={slope_b:.3f}')
    
    axes[1].plot([-90, 90], [-90, 90], 'k--', lw=1, alpha=0.5)
    axes[1].set_title(r'$\alpha=0, \beta=\pm 20$')
    axes[1].set_xlabel('Scalar Interpolation (degrees)')
    axes[1].set_ylabel('Fiber angles (degrees)')
    axes[1].legend(fontsize=8, loc='upper left')

    # Combined
    axes[2].scatter(ref_alpha_combined, alpha_combined, alpha=0.02, s=10, color='blue')
    axes[2].scatter(ref_beta_combined, beta_combined, alpha=0.02, s=10, color='red')
    
    slope_a, intercept_a, r_value_a, _, _ = stats.linregress(ref_alpha_combined, alpha_combined)
    axes[2].plot(x_line, slope_a * x_line + intercept_a, 'b-', lw=1, label=f'α: R²={r_value_a**2:.3f}, m={slope_a:.3f}')
    
    slope_b, intercept_b, r_value_b, _, _ = stats.linregress(ref_beta_combined, beta_combined)
    axes[2].plot(x_line, slope_b * x_line + intercept_b, 'r-', lw=1, label=f'β: R²={r_value_b**2:.3f}, m={slope_b:.3f}')
    axes[2].plot([-90, 90], [-90, 90], 'k--', lw=1, alpha=0.5)
    axes[2].set_title(r'$\alpha=\pm 60, \beta=\pm 20$')
    axes[2].set_xlabel('Scalar Interpolation (degrees)')
    axes[2].set_ylabel('Fiber angles (degrees)')
    axes[2].legend(fontsize=8, loc='upper left')

    plt.savefig('example/biv_with_outflow_tracts/doste_angle_correlations.png', dpi=150)
