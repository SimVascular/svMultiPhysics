#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
# SPDX-License-Identifier: BSD-3-Clause
# -*-coding:utf-8 -*-
"""Surface mesh processing utilities for biventricular heart models.

This module provides classes and utilities for processing surface meshes,
including generating epicardial apex surfaces.
"""

import os
import numpy as np
import pyvista as pv
from .surface_names import SurfaceName


def get_normal_plane_svd(points):
    """Find the plane that best fits a set of points using SVD.
    
    Args:
        points: Array of shape (N, 3) representing 3D points.
    
    Returns:
        tuple: A tuple containing:
            - normal (np.ndarray): Unit normal vector to the fitted plane.
            - centroid (np.ndarray): Centroid of the input points.
    """
    centroid = np.mean(points, axis=0)
    svd = np.linalg.svd(points - centroid)
    normal = svd[2][-1]
    normal = normal / np.linalg.norm(normal)
    return normal, centroid


def generate_epi_apex(surface_paths):
    """Generate the epicardial apex surface from the epicardial surface of the BiV.
    
    This function identifies the apex point of the epicardium and creates a surface
    mesh containing elements that include the apex point. The surface is saved with
    global node and element IDs.
    
    Args:
        surface_paths: Dictionary mapping SurfaceName enum values to full file paths.
    
    The function requires:
        - SurfaceName.EPICARDIUM: Epicardial surface
        - SurfaceName.BASE: Base surface (for finding apex)
        - SurfaceName.EPICARDIUM_APEX: Output surface name (will be created)
    """
    # Load the epi surface
    epi_path = surface_paths[SurfaceName.EPICARDIUM]
    epi_mesh = pv.read(epi_path)
    epi_points = epi_mesh.points
    epi_cells = epi_mesh.faces
    epi_eNoN = epi_cells[0]
    epi_cells = epi_cells.reshape((-1, epi_eNoN + 1))
    epi_cells = epi_cells[:, 1:]
    epi_global_node_id = epi_mesh.point_data['GlobalNodeID']
    epi_global_cell_id = epi_mesh.cell_data['GlobalElementID']

    # Load the base surface
    base_path = surface_paths[SurfaceName.BASE]
    base_mesh = pv.read(base_path)
    base_global_node_id = base_mesh.point_data['GlobalNodeID']

    # Extract the boundary of the epi surface (at the top) to find the apex point
    epi_base_global_node_id = np.intersect1d(epi_global_node_id, base_global_node_id)
    epi_base_nodes = np.where(np.isin(epi_global_node_id, epi_base_global_node_id))[0]

    # Get normal
    base_normal, base_centroid = get_normal_plane_svd(epi_points[epi_base_nodes, :])

    # Find the index of the apex point of the epi surface
    distance = np.abs(base_normal @ (epi_points - base_centroid).T)
    epi_apex_point_index = np.argmax(distance)

    # Find elements containing the apex point
    epi_apex_cell_index = np.where(epi_cells == epi_apex_point_index)[0]

    # Create epi_apex mesh
    submesh_cells = epi_cells[epi_apex_cell_index]
    submesh_xyz = np.zeros([len(np.unique(submesh_cells)), epi_points.shape[1]])
    map_mesh_submesh = np.ones(epi_points.shape[0], dtype=int) * -1
    map_submesh_mesh = np.zeros(submesh_xyz.shape[0], dtype=int)
    child_elems_new = np.zeros(submesh_cells.shape, dtype=int)

    cont = 0
    for e in range(submesh_cells.shape[0]):
        for i in range(submesh_cells.shape[1]):
            if map_mesh_submesh[submesh_cells[e, i]] == -1:
                child_elems_new[e, i] = cont
                submesh_xyz[cont] = epi_points[submesh_cells[e, i]]
                map_mesh_submesh[submesh_cells[e, i]] = cont
                map_submesh_mesh[cont] = submesh_cells[e, i]
                cont += 1
            else:
                child_elems_new[e, i] = map_mesh_submesh[submesh_cells[e, i]]

    epi_apex_cells_type = np.full((child_elems_new.shape[0], 1), epi_eNoN)
    epi_apex_cells = np.hstack((epi_apex_cells_type, child_elems_new))
    epi_apex_cells = np.hstack(epi_apex_cells)

    # Get global IDs
    epi_apex_global_node_id = epi_global_node_id[map_submesh_mesh]
    epi_apex_global_cell_id = epi_global_cell_id[epi_apex_cell_index]

    # Create and save mesh
    epi_apex_mesh = pv.PolyData(submesh_xyz, epi_apex_cells)
    epi_apex_mesh.point_data.set_array(epi_apex_global_node_id, 'GlobalNodeID')
    epi_apex_mesh.cell_data.set_array(epi_apex_global_cell_id, 'GlobalElementID')

    epi_apex_path = surface_paths[SurfaceName.EPICARDIUM_APEX]
    epi_apex_mesh.save(epi_apex_path)
