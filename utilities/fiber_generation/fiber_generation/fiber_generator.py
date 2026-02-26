#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
# SPDX-License-Identifier: BSD-3-Clause
# -*-coding:utf-8 -*-
"""Object-oriented fiber generation module for biventricular heart models.

This module provides classes to generate myocardial fiber orientations for
biventricular heart models using Laplace-Dirichlet rule-based methods.

Supports both:
    - Bayer et al. (2012): Truncated BiV geometry
    - Doste et al. (2019): BiV geometry with outflow tracts

References:
    Bayer et al. 2012: https://doi.org/10.1007/s10439-012-0593-5
    Doste et al. 2019: https://doi.org/10.1002/cnm.3185
"""

import os
import numpy as np
import pyvista as pv

import fiber_generation.quat_utils as qu


class FibGen:
    """Base class for fiber generation.
    
    Provides common utilities for computing fiber orientations from
    Laplace field solutions.
    
    Attributes:
        mesh: PyVista mesh with cell-centered data.
        lap: Dictionary of Laplace solution values at cells.
        grad: Dictionary of gradient arrays at cells (N, 3).
    """
    
    def __init__(self):
        """Initialize the FibGen base class."""
        self.mesh = None
        self.lap = None
        self.grad = None
    
    def normalize(self, x):
        """Normalize each row of an (N, 3) array.
        
        Zero-length rows remain zero after normalization.
        
        Args:
            x: Array-like of shape (N, 3).
        
        Returns:
            np.ndarray of shape (N, 3) with row-wise normalized vectors.
        """
        a = np.asarray(x, dtype=float)
        if a.ndim != 2 or a.shape[1] != 3:
            raise ValueError("normalize expects an array of shape (N, 3)")
        norms = np.linalg.norm(a, axis=1, keepdims=True)
        safe_norms = np.where(norms == 0.0, 1.0, norms)
        out = a / safe_norms
        zero_rows = (norms.squeeze() == 0.0)
        if np.any(zero_rows):
            out[zero_rows] = 0.0
        return out
    
    def scale_to_range(self, arr, range=(0.0, 1.0)):
        """Scale array to specified range.
        
        Args:
            arr: Input array to scale.
            range: Tuple (min, max) specifying the target range.
            
        Returns:
            np.ndarray: Scaled array with values in the specified range. If all values are equal,
                returns array filled with the midpoint of the range.
        """
        arr = np.asarray(arr, dtype=float)
        amin = np.min(arr)
        amax = np.max(arr)
        if amax > amin:
            return range[0] + (arr - amin) * (range[1] - range[0]) / (amax - amin)
        else:
            return np.ones_like(arr) * ((range[0] + range[1]) / 2.0)
    
    def compute_gradients(self, mesh, field_names):
        """Compute gradients for specified fields at points.
        
        Args:
            mesh: PyVista mesh with point data.
            field_names: List of field names to compute gradients for.
            
        Returns:
            PyVista mesh with gradient arrays added to point_data.
        """
        for name in field_names:
            if name not in mesh.point_data:
                raise KeyError(f"Field '{name}' not found in mesh point_data")
            
            gmesh = mesh.compute_derivative(scalars=name, gradient=True, preference='point')
            mesh.point_data[name + "_grad"] = np.asarray(gmesh.point_data["gradient"])
        
        return mesh
    
    def axis(self, gL, gT):
        """Construct orthogonal coordinate systems from two gradient fields.
        
        Creates an orthonormal basis [eC, eL, eT] for each element where:
        - eL is aligned with gL (normalized longitudinal)
        - eT is orthogonal to eL and in the plane of gT (transmural)
        - eC is the cross product of eL and eT (circumferential)
        
        Args:
            gL: Array of shape (N, 3) representing the longitudinal gradient.
            gT: Array of shape (N, 3) representing the transmural gradient.
        
        Returns:
            np.ndarray: Array of shape (N, 3, 3) where columns are
                [eC (circumferential), eL (longitudinal), eT (transmural)].
        """
        gL = np.asarray(gL, dtype=float)
        gT = np.asarray(gT, dtype=float)
        ne = gL.shape[0]
        
        # eL = normalized longitudinal
        eL = self.normalize(gL)
        
        # eT = gT - proj_{eL}(gT), orthogonal to eL
        proj = np.sum(eL * gT, axis=1)[:, None] * eL
        eT = gT - proj
        eT = self.normalize(eT)
        
        # eC = cross(eL, eT), circumferential
        eC = np.cross(eL, eT, axisa=1, axisb=1)
        eC = self.normalize(eC)
        
        # Build basis matrix Q = [eC, eL, eT]
        Q = np.zeros((ne, 3, 3), dtype=float)
        Q[:, :, 0] = eC
        Q[:, :, 1] = eL
        Q[:, :, 2] = eT
        
        return Q
    
    def calculate_angle(self, trans, endo_value, epi_value):
        """Compute angle varying linearly from endo to epi.
        
        Args:
            trans: Transmural coordinate array (N,), values in [0, 1].
            endo_value: Angle value at endocardium (scalar).
            epi_value: Angle value at epicardium (scalar).
        
        Returns:
            np.ndarray: Angle values at each point (N,).
        """
        return endo_value * (1 - trans) + epi_value * trans
    
    def orient_matrix(self, Q, alpha, beta):
        """Apply alpha and beta rotations to orthogonal matrices.
        
        Rotates Q by alpha about the z-axis (transmural) and then
        by beta about the y-axis (longitudinal direction).
        
        Args:
            Q: Array of shape (N, 3, 3) containing orthogonal matrices.
            alpha: Array of shape (N,) with rotation angles (radians) about z-axis.
            beta: Array of shape (N,) with rotation angles (radians) about y-axis.
        
        Returns:
            np.ndarray: Array of shape (N, 3, 3) containing rotated matrices.
        """
        Q = np.asarray(Q, dtype=float)
        ne = Q.shape[0]
        
        ca = np.cos(alpha)
        sa = np.sin(alpha)
        cb = np.cos(beta)
        sb = np.sin(beta)
        
        # Rotation about z-axis (Ra)
        Ra = np.zeros((ne, 3, 3), dtype=float)
        Ra[:, 0, 0] = ca
        Ra[:, 0, 1] = -sa
        Ra[:, 1, 0] = sa
        Ra[:, 1, 1] = ca
        Ra[:, 2, 2] = 1.0
        
        # Rotation about y-axis (Rb)
        Rb = np.zeros((ne, 3, 3), dtype=float)
        Rb[:, 0, 0] = cb
        Rb[:, 0, 2] = sb
        Rb[:, 1, 1] = 1.0
        Rb[:, 2, 0] = -sb
        Rb[:, 2, 2] = cb
        
        # Compose rotations and apply to Q
        RaRb = np.einsum('nij,njk->nik', Ra, Rb)
        Qt = np.einsum('nij,njk->nik', Q, RaRb)
        
        return Qt
        
    def orient_rodrigues(self, Q, alpha, beta):
        """Rotate basis using Rodrigues rotation formula (Doste method).
        
        Applies two successive rotations using Rodrigues formula:
        1. Rotate by alpha about the transmural axis (eT)
        2. Rotate by beta about the rotated longitudinal axis
        
        Args:
            Q: Array of shape (N, 3, 3) containing basis matrices.
                Columns are [eC (circumferential), eL (longitudinal), eT (transmural)].
            alpha: Array of shape (N,) with rotation angles (radians) about transmural axis.
            beta: Array of shape (N,) with rotation angles (radians) about rotated longitudinal axis.
        
        Returns:
            np.ndarray: Array of shape (N, 3, 3) containing rotated basis matrices.
        """
        Q = np.asarray(Q, dtype=float)
        alpha = np.asarray(alpha, dtype=float)
        beta = np.asarray(beta, dtype=float)
        
        n = Q.shape[0]
        
        # Extract basis vectors
        eC = Q[:, :, 0]  # Circumferential
        eL = Q[:, :, 1]  # Longitudinal
        eT = Q[:, :, 2]  # Transmural
        
        # Normalize basis vectors
        eC = self.normalize(eC)
        eL = self.normalize(eL)
        eT = self.normalize(eT)

        # Vectorized Rodrigues rotation:
        # v_rot = v*cos(theta) + (k x v)*sin(theta) + k*(k·v)*(1-cos(theta))

        def rot(v, k, theta):
            v = np.asarray(v, dtype=float)
            k = np.asarray(k, dtype=float)
            theta = np.asarray(theta, dtype=float)

            if v.ndim != 2 or v.shape[1] != 3:
                raise ValueError("orient_rodrigues expects vectors of shape (N, 3)")
            if k.shape != v.shape:
                raise ValueError("orient_rodrigues expects axes of shape (N, 3) matching vectors")
                
            # cos(theta) and sin(theta) are broadcasted to (N, 1)
            ct = np.cos(theta)[:, None]
            st = np.sin(theta)[:, None]

            kv = np.cross(k, v, axis=1)
            kdotv = np.einsum('ij,ij->i', k, v)[:, None]
            return v * ct + kv * st + k * kdotv * (1.0 - ct)

        # First rotation: alpha about transmural axis eT
        eC1 = rot(eC, eT, alpha)
        eL1 = rot(eL, eT, alpha)
        eT1 = eT  # unchanged

        # Second rotation: beta about rotated longitudinal axis eL1 (normalized)
        eL1_axis = self.normalize(eL1)
        eC2 = rot(eC1, eL1_axis, beta)
        eL2 = eL1  # unchanged (rotation about itself)
        eT2 = rot(eT1, eL1_axis, beta)

        result = np.zeros((n, 3, 3), dtype=float)
        result[:, :, 0] = eC2
        result[:, :, 1] = eL2
        result[:, :, 2] = eT2

        return result

    
    def interpolate_basis(self, Q1, Q2, t, correct_slerp=False):
        """Spherical linear interpolation between batches of rotation matrices.
        
        Performs SLERP on rotation matrices represented as quaternions internally.
        
        Args:
            Q1: Array of shape (N, 3, 3) containing starting rotation matrices.
            Q2: Array of shape (N, 3, 3) containing ending rotation matrices.
            t: Array of shape (N,) with interpolation values in [0, 1].
            correct_slerp: If True, use quaternion correction to ensure shortest path.
                Defaults to False.
        
        Returns:
            np.ndarray: Array of shape (N, 3, 3) containing interpolated rotation matrices.
        """
        
        # Prepare inputs
        t = np.clip(np.asarray(t, dtype=float), 0.0, 1.0)
        
        # Ensure shortest path on the unit 4-sphere
        if correct_slerp:
            q1 = np.zeros((len(t), 4), dtype=float)
            q2 = np.zeros((len(t), 4), dtype=float)
            q1, q2 = qu.find_best_quaternions(Q1, Q2)
            dot = np.einsum('ni,ni->n', q1, q2)
        else:
            q1 = qu.rotm_to_quat_batch(Q1)
            q2 = qu.rotm_to_quat_batch(Q2)
            dot = np.einsum('ni,ni->n', q1, q2)
        if np.any(dot < 0.0):
            neg_mask = dot < 0.0
            q2[neg_mask] = -q2[neg_mask]
            dot[neg_mask] = -dot[neg_mask]
            
        # SLERP weights
        dot_clipped = np.clip(dot, -1.0, 1.0)
        theta0 = np.arccos(dot_clipped)
        sin_theta0 = np.sin(theta0)
        
        # Threshold for linear interpolation
        lin_mask = sin_theta0 < 1e-6
        q = np.empty_like(q1)
        
        if np.any(~lin_mask):
            theta = theta0[~lin_mask] * t[~lin_mask]
            s0 = np.sin(theta0[~lin_mask] - theta) / sin_theta0[~lin_mask]
            s1 = np.sin(theta) / sin_theta0[~lin_mask]
            q[~lin_mask] = (s0[:, None] * q1[~lin_mask]) + (s1[:, None] * q2[~lin_mask])
        
        if np.any(lin_mask):
            tl = t[lin_mask][:, None]
            q[lin_mask] = (1.0 - tl) * q1[lin_mask] + tl * q2[lin_mask]
        
        # Normalize and convert back to rotation matrices
        q /= np.linalg.norm(q, axis=1, keepdims=True)
        return qu.quat_to_rotm_batch(q)
    
    def generate_fibers(self, params):
        """Generate fiber directions. Override in subclasses."""
        raise NotImplementedError("Subclasses must implement generate_fibers()")
    
    def write_fibers(self, outdir):
        """Write fiber, sheet, and normal directions to VTU files.
        
        Saves three separate files with fiber directions stored in 'FIB_DIR' field:
        - fiber.vtu: Fiber directions
        - sheet.vtu: Sheet normal directions
        - normal.vtu: Sheet-normal directions
        
        Args:
            outdir: Output directory path where files will be saved.
        """
        # Create a copy of the mesh without any data
        mesh_out = self.mesh.copy(deep=True)
        mesh_out.clear_data()

        # Fiber direction
        mesh_out.cell_data['FIB_DIR'] = self.mesh.cell_data['fiber']
        mesh_out.save(os.path.join(outdir, "fiber.vtu"))

        # Sheet direction
        mesh_out.cell_data['FIB_DIR'] = self.mesh.cell_data['sheet']
        mesh_out.save(os.path.join(outdir, "sheet.vtu"))

        # Normal direction
        mesh_out.cell_data['FIB_DIR'] = self.mesh.cell_data['sheet-normal']
        mesh_out.save(os.path.join(outdir, "normal.vtu"))


class FibGenBayer(FibGen):
    """Fiber generator using the Bayer et al. (2012) method.
    
    Suitable for truncated biventricular geometries. Implements the rule-based
    algorithm described in Bayer et al. 2012:
    https://doi.org/10.1007/s10439-012-0593-5
    """
    
    # Field names in Laplace solution
    FIELD_NAMES = ['Trans_EPI', 'Trans_LV', 'Trans_RV', 'Long_AB']
    
    def __init__(self):
        """Initialize the Bayer fiber generator."""
        super().__init__()


    def rescale_fields(self, mesh):
        """Rescale Laplace fields to [0, 1] range.
        
        Args:
            mesh: PyVista mesh with point data containing the fields to rescale.
        
        Returns:
            PyVista mesh with rescaled fields in point_data.
        """
        for name in self.FIELD_NAMES:
            if name not in mesh.point_data:
                raise KeyError(f"Field '{name}' not found in mesh point_data")
            mesh.point_data[name] = self.scale_to_range(mesh.point_data[name], range=(0.0, 1.0))
        return mesh
    
    
    def load_laplace_results(self, file_path):
        """Load Laplace-Dirichlet solution for Bayer method.
        
        Args:
            file_path: Path to the .vtu file with Laplace solution.
        
        Returns:
            tuple: (lap, grad) dictionaries with Laplace values and gradients.
        """
        print(f"   Loading Laplace solution <--- {file_path}")
        result_mesh = pv.read(file_path)

        # Normalize fields to [0, 1]
        result_mesh = self.rescale_fields(result_mesh)

        # Compute gradients for the required fields
        print("   Computing gradients at points")
        result_mesh = self.compute_gradients(result_mesh, self.FIELD_NAMES)
        
        # Convert point-data to cell-data
        mesh_cells = result_mesh.point_data_to_cell_data()
        self.mesh = mesh_cells
        
        # Extract Laplace values and gradients 
        self.lap = {}
        self.grad = {}
        
        for key in self.FIELD_NAMES:
            self.lap[key] = np.asarray(mesh_cells.cell_data[key])
            self.grad[key] = np.asarray(mesh_cells.cell_data[key + "_grad"])
        
        
        return self.lap, self.grad
    

    def generate_fibers(self, params, flip_rv=True, correct_slerp=False):
        """Generate fiber directions using the Bayer method.
        
        Args:
            params: Dictionary with keys:
                - ALFA_END: Endocardial helix angle (degrees)
                - ALFA_EPI: Epicardial helix angle (degrees)
                - BETA_END: Endocardial transverse angle (degrees)
                - BETA_EPI: Epicardial transverse angle (degrees)
            flip_rv: If True, flip circumferential and transmural directions in RV.
                Defaults to True.
            correct_slerp: If True, use quaternion correction for SLERP interpolation.
                Defaults to False.
        
        Returns:
            tuple: (F, S, T) fiber, sheet, and normal directions (N, 3) each.
        """
        if self.lap is None or self.grad is None:
            raise ValueError("Must call load_laplace_results() first")
        
        # Convert parameters to radians (consistent with Doste method)
        params = {k: np.deg2rad(v) for k, v in params.items()}
                
        print("   Computing fiber directions at cells")
        
        # Interpolation factor between LV and RV
        d = self.lap['Trans_RV'] / (self.lap['Trans_LV'] + self.lap['Trans_RV'])
        alfaS = self.calculate_angle(d, params['ALFA_END'], -params['ALFA_END'])
        betaS = self.calculate_angle(d, params['BETA_END'], -params['BETA_END'])
                    
        # Wall angles (interpolated from endo to epi)
        alfaW = self.calculate_angle(self.lap['Trans_EPI'], params['ALFA_END'], params['ALFA_EPI'])
        betaW = self.calculate_angle(self.lap['Trans_EPI'], params['BETA_END'], params['BETA_EPI'])
        
        # Build LV and RV basis
        Q_LV0 = self.axis(self.grad['Long_AB'], -self.grad['Trans_LV'])
        Q_LV = self.orient_matrix(Q_LV0, alfaS, np.sign(params['BETA_END'])*np.abs(betaS))
        
        Q_RV0 = self.axis(self.grad['Long_AB'], self.grad['Trans_RV']) 
        Q_RV = self.orient_matrix(Q_RV0, alfaS, np.sign(params['BETA_END'])*np.abs(betaS))    
        
        # Interpolate between LV and RV (endocardial layer)
        Q_END = self.interpolate_basis(Q_LV, Q_RV, d, correct_slerp=correct_slerp)

        # Flip circumferential and transmural directions in RV
        if flip_rv:
            Q_END[d > 0.5,:,0] = -Q_END[d > 0.5,:,0]
            Q_END[d > 0.5,:,2] = -Q_END[d > 0.5,:,2]
        
        # Build epicardial basis
        Q_EPI0 = self.axis(self.grad['Long_AB'], self.grad['Trans_EPI'])
        Q_EPI = self.orient_matrix(Q_EPI0, alfaW, betaW)
        
        # Interpolate from endo to epi
        FST = self.interpolate_basis(Q_END, Q_EPI, self.lap['Trans_EPI'], correct_slerp=correct_slerp)
        
        F = FST[:, :, 0]  # Fiber direction
        S = FST[:, :, 1]  # Sheet normal
        T = FST[:, :, 2]  # Sheet direction
        
        self.mesh.cell_data['fiber'] = F
        self.mesh.cell_data['sheet-normal'] = S
        self.mesh.cell_data['sheet'] = T
    
        return F, S, T
        
    def get_angle_fields(self, params):
        """Compute global alpha and beta angle fields.
        
        Helper function to compute spatially-varying helix and transverse angle fields
        by interpolating between septum and wall values.
        
        Args:
            params: Dictionary with angle parameters (in degrees or radians).
            
        Returns:
            tuple: (alfa, beta) arrays of helix and transverse angles at each cell.
        """

        # Interpolation factor between LV and RV
        d = self.lap['Trans_RV'] / (self.lap['Trans_LV'] + self.lap['Trans_RV'])
        
        # Septum angles (interpolated between LV and RV)
        alfaS = self.calculate_angle(d, params['ALFA_END'], -params['ALFA_END'])
        betaS = self.calculate_angle(d, params['BETA_END'], -params['BETA_END'])
        alfaS = np.abs(alfaS)   # Note this is doing the same as flipping the sign
        betaS = np.sign(params['BETA_END'])*np.abs(betaS)   # Note this is doing the same as flipping the sign
        
        # Wall angles (interpolated from endo to epi)
        alfaW = self.calculate_angle(self.lap['Trans_EPI'], params['ALFA_END'], params['ALFA_EPI'])
        betaW = self.calculate_angle(self.lap['Trans_EPI'], params['BETA_END'], params['BETA_EPI'])

        alfa = alfaS * (1 - self.lap['Trans_EPI']) + alfaW * self.lap['Trans_EPI']
        beta = betaS * (1 - self.lap['Trans_EPI']) + betaW * self.lap['Trans_EPI']

        return alfa, beta


class FibGenDoste(FibGen):
    """Fiber generator using the Doste et al. (2019) method.
    
    Suitable for biventricular geometries with outflow tracts. Implements
    the algorithm described in Doste et al. 2019:
    https://doi.org/10.1002/cnm.3185
    """
    
    # Field names in Laplace solution
    FIELD_NAMES = ['Trans_BiV', 'Long_AV', 'Long_MV', 'Long_PV', 'Long_TV',
                   'Weight_LV', 'Weight_RV', 'Trans_EPI', 'Trans_LV', 'Trans_RV', 'Trans']
    
    def __init__(self):
        """Initialize the Doste fiber generator."""
        super().__init__()
    

    def rescale_fields(self, mesh):
        """Rescale Laplace fields to [0, 1] range.
        
        Args:
            mesh: PyVista mesh with point data containing the fields to rescale.
        
        Returns:
            PyVista mesh with rescaled fields in point_data.
        """
        for name in self.FIELD_NAMES:
            if name not in mesh.point_data:
                raise KeyError(f"Field '{name}' not found in mesh point_data")
            
            # Set respective range
            if name in ['Trans_BiV', 'Trans']:
                range = (-2.0, 1.0)
            else:
                range = (0.0, 1.0)

            mesh.point_data[name] = self.scale_to_range(mesh.point_data[name], range=range)

        return mesh
    

    def load_laplace_results(self, file_path):
        """Load Laplace-Dirichlet solution for Doste method.
        
        Args:
            file_path: Path to the .vtu file with Laplace solution.
        
        Returns:
            tuple: (lap, grad) dictionaries with Laplace values and gradients.
        """
        print(f"   Loading Laplace solution <--- {file_path}")
        result_mesh = pv.read(file_path)

        # Normalize fields to the original range
        result_mesh = self.rescale_fields(result_mesh)
        
        # Compute gradients for the required fields
        print("   Computing gradients at points")
        result_mesh = self.compute_gradients(result_mesh, self.FIELD_NAMES)
        
        # Convert point-data to cell-data
        mesh_cells = result_mesh.point_data_to_cell_data()
        self.mesh = mesh_cells
        
        # Extract Laplace values and gradients using mapped names
        self.lap = {}
        self.grad = {}
        
        for key in self.FIELD_NAMES:
            self.lap[key] = np.asarray(mesh_cells.cell_data[key])
            self.grad[key] = np.asarray(mesh_cells.cell_data[key + "_grad"])
        
        return self.lap, self.grad
    
    
    def _redistribute_weight(self, weight, up, low):
        """Redistribute weight values to center their distribution.
        
        Args:
            weight: Array of weight values.
            up: Upper quantile threshold.
            low: Lower quantile threshold.
        
        Returns:
            np.ndarray: Redistributed weight values in [0, 1].
        """
        new_weight = weight.copy()
        
        upper_lim = np.quantile(weight, up)
        while upper_lim == 0:
            up += 0.1
            upper_lim = np.quantile(weight, up)
        
        lower_lim = np.quantile(weight, low)
        
        new_weight[new_weight > upper_lim] = upper_lim
        new_weight[new_weight < lower_lim] = lower_lim
        
        return (new_weight - np.min(new_weight)) / (np.max(new_weight) - np.min(new_weight))
    
    def _compute_basis_vectors(self):
        """Compute local orthogonal basis vectors for LV and RV.
        
        Returns:
            dict: Dictionary with basis vectors for LV, RV, and global.
        """
        lap, grad = self.lap, self.grad
        
        # Calculate combined LV longitudinal
        lv_glong = (grad['Long_MV'] * lap['Weight_LV'][:, None] + 
                   grad['Long_AV'] * (1 - lap['Weight_LV'][:, None]))

        # Calculate LV basis
        Q_lv = self.axis(-lv_glong, -grad['Trans_LV'])   # Minus signs to match Bayer convention
        eC_lv = Q_lv[:, :, 0]  # Circumferential
        eL_lv = Q_lv[:, :, 1]  # Longitudinal
        eT_lv = Q_lv[:, :, 2]  # Transmural
        
        # Calculate combined RV longitudinal
        rv_glong = (grad['Long_TV'] * lap['Weight_RV'][:, None] + 
                   grad['Long_PV'] * (1 - lap['Weight_RV'][:, None]))
        Q_rv = self.axis(-rv_glong, -grad['Trans_RV'])   # Minus signs to match Bayer convention
        eC_rv = Q_rv[:, :, 0]  # Circumferential
        eL_rv = Q_rv[:, :, 1]  # Longitudinal
        eT_rv = Q_rv[:, :, 2]  # Transmural
        
        return {
            'eC_lv': eC_lv, 'eT_lv': eT_lv, 'eL_lv': eL_lv,
            'eC_rv': eC_rv, 'eT_rv': eT_rv, 'eL_rv': eL_rv,
        }
    
    def _compute_angles(self, params):
        """Compute spatially-varying alpha and beta angles.
        
        Args:
            params: Dictionary with angle parameters (must be in radians).
        
        Returns:
            dict: Dictionary of angle arrays including:
                - alpha_wall_lv, alpha_wall_rv: Wall helix angles for LV/RV
                - beta_wall_lv, beta_wall_rv: Wall transverse angles for LV/RV
                - alfaS: Septum helix angle
                - beta_septum: Septum transverse angle
        """
        lap = self.lap
        
        # Redistribute weights
        lv_weight = self._redistribute_weight(lap['Weight_LV'], 0.7, 0.01)
        rv_weight = self._redistribute_weight(lap['Weight_RV'], 0.1, 0.001)
        
        # LV angles
        alpha_lv_endo = params['AENDOLV'] * lv_weight + params['AOTENDOLV'] * (1 - lv_weight)
        alpha_lv_epi = params['AEPILV'] * lv_weight + params['AOTEPILV'] * (1 - lv_weight)
        alpha_wall_lv = self.calculate_angle(lap['Trans_EPI'], alpha_lv_endo, alpha_lv_epi)
        beta_wall_lv = self.calculate_angle(lap['Trans_EPI'], params['BENDOLV'], params['BEPILV']) * lv_weight
        
        # RV angles
        alpha_rv_endo = params['AENDORV'] * rv_weight + params['AOTENDORV'] * (1 - rv_weight)
        alpha_rv_epi = params['AEPIRV'] * rv_weight + params['AOTEPIRV'] * (1 - rv_weight)
        alpha_wall_rv = self.calculate_angle(lap['Trans_EPI'], alpha_rv_endo, alpha_rv_epi)
        beta_wall_rv = self.calculate_angle(lap['Trans_EPI'], params['BENDORV'], params['BEPIRV']) * rv_weight
        
        # Septum angles
        sep = lap['Trans'].copy()
        sep[sep < 0.0] = sep[sep < 0.0] / 2
        sep = np.abs(sep)   # This gives a field that is 1 at both endo but that assigns 2/3  of the septum to the lv
        alpha_septum = alpha_lv_endo * sep * (lap['Trans_BiV'] < 0) + alpha_rv_endo * sep * (lap['Trans_BiV'] > 0)
        beta_septum = params['BENDOLV'] * sep * (lap['Trans_BiV'] < 0)  * lv_weight + params['BENDORV'] * sep * (lap['Trans_BiV'] > 0) * rv_weight
        
        return {
            'alpha_wall_lv': alpha_wall_lv, 'beta_wall_lv': beta_wall_lv,
            'alpha_wall_rv': alpha_wall_rv, 'beta_wall_rv': beta_wall_rv,
            'alpha_septum': alpha_septum, 'beta_septum': beta_septum
        }
    
    
    def generate_fibers(self, params):
        """Generate fiber directions using the Doste method.
        
        Args:
            params: Dictionary with angle parameters (in degrees):
                - AENDOLV, AEPILV: LV endo/epi helix angles
                - AENDORV, AEPIRV: RV endo/epi helix angles
                - AOTENDOLV, AOTEPILV: LV outflow tract angles
                - AOTENDORV, AOTEPIRV: RV outflow tract angles
                - BENDOLV, BEPILV: LV endo/epi transverse angles
                - BENDORV, BEPIRV: RV endo/epi transverse angles
        
        Returns:
            tuple: (F, S, T) fiber, sheet, and normal directions (N, 3) each.
        """
        if self.lap is None or self.grad is None:
            raise ValueError("Must call load_laplace_results() first")
        
        # Convert parameters to radians (consistent with Bayer method)
        params_rad = {k: np.deg2rad(v) for k, v in params.items()}
        
        print("   Computing basis vectors")
        basis = self._compute_basis_vectors()
        
        print("   Computing angles")
        angles = self._compute_angles(params_rad)
        
        print("   Computing local basis")
        # Build basis matrices from vectors
        Q_lv = np.stack([basis['eC_lv'], basis['eL_lv'], basis['eT_lv']], axis=-1)
        Q_rv = np.stack([basis['eC_rv'], basis['eL_rv'], basis['eT_rv']], axis=-1)
        
        # Septum basis
        Qlv_sep = self.orient_rodrigues(
            Q_lv, angles['alpha_septum'], angles['beta_septum']
        )
        Qrv_sep = self.orient_rodrigues(
            Q_rv, angles['alpha_septum'], angles['beta_septum']
        )
        
        # Wall basis
        Qlv_wall = self.orient_rodrigues(
            Q_lv, angles['alpha_wall_lv'], angles['beta_wall_lv']
        )
        Qrv_wall = self.orient_rodrigues(
            Q_rv, angles['alpha_wall_rv'], angles['beta_wall_rv']
        )
        
        print("   Interpolating basis")
        # Get interpolation factor between LV and RV
        interp_biv = self.lap['Trans_BiV'].copy()  # Need to normalize this to [0, 1] for interpolation
        interp_biv[interp_biv < 0.0] = interp_biv[interp_biv < 0.0] / 2
        interp_biv = (interp_biv + 1) / 2

        # Get discontinous septal fibers
        Qsep = Qrv_sep.copy()
        Qsep[interp_biv < 0.5] = Qlv_sep[interp_biv < 0.5]
        
        # Interpolate across ventricles
        Qepi = self.interpolate_basis(Qlv_wall, Qrv_wall, interp_biv)
        
        # Interpolate from endo to epi
        Q = self.interpolate_basis(Qsep, Qepi, self.lap['Trans_EPI'])
        
        print("   Done!")
        F = Q[:, :, 0]  # Fiber direction
        S = Q[:, :, 1]  # Sheet normal
        T = Q[:, :, 2]  # Sheet direction
        
        self.mesh.cell_data['fiber'] = F
        self.mesh.cell_data['sheet-normal'] = S
        self.mesh.cell_data['sheet'] = T
        
        return F, S, T


    def get_angle_fields(self, params):
        """Compute global alpha and beta angle fields.
        
        Helper function to compute spatially-varying helix and transverse angle fields
        for the Doste method by interpolating between septum and wall values across
        ventricles and transmural depth.
        
        Args:
            params: Dictionary with angle parameters (in degrees or radians).
            
        Returns:
            tuple: (alfa, beta) arrays of helix and transverse angles at each cell.
        """

        # Compute angles
        angles = self._compute_angles(params)
        
        # Normalize Trans_BiV for interpolation (matching generate_fibers implementation)
        interp_biv = self.lap['Trans_BiV'].copy()
        interp_biv[interp_biv < 0.0] = interp_biv[interp_biv < 0.0] / 2
        interp_biv = (interp_biv + 1) / 2

        # Interpolate wall angles between LV and RV
        alpha_epi = angles['alpha_wall_lv'] * (1 - interp_biv) + angles['alpha_wall_rv'] * interp_biv
        beta_epi = angles['beta_wall_lv'] * (1 - interp_biv) + angles['beta_wall_rv'] * interp_biv

        # Interpolate from septum (endo) to wall (epi) transmurally
        alfa = angles['alpha_septum'] * (1 - self.lap['Trans_EPI']) + alpha_epi * self.lap['Trans_EPI']
        beta = angles['beta_septum'] * (1 - self.lap['Trans_EPI']) + beta_epi * self.lap['Trans_EPI']
        
        return alfa, beta