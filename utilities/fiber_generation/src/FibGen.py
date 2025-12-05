#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
Created on 2025/11/21 20:43:23

@author: Javiera Jilberto Vallejos 
'''

import os
import re
import numpy as np
import pyvista as pv
import time
import copy


def normalize(x):
    """
    Normalize each row of an (N, 3) array. Zero rows remain zero.
    
    Args:
        x: array-like of shape (N, 3)
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



def get_normal_plane_svd(points):   # Find the plane that minimizes the distance given N points
    centroid = np.mean(points, axis=0)
    svd = np.linalg.svd(points - centroid)
    normal = svd[2][-1]
    normal = normal/np.linalg.norm(normal)
    return normal, centroid


def generate_epi_apex(mesh_path, surfaces_dir, surface_names):
    '''
    Generate the epi apex and epi mid surfaces from the epi surface of the BiV.
    
    Parameters:
    -----------
    surfaces_dir : str
        Directory containing surface meshes
    surface_names : list of str
        List of surface mesh filenames
    '''
    
    # Load the epi surface
    epi_name = os.path.join(surfaces_dir, surface_names['epi'])
    epi_mesh = pv.read(epi_name)
    epi_points = epi_mesh.points
    epi_cells = epi_mesh.faces
    epi_eNoN = epi_cells[0]
    epi_cells = epi_cells.reshape((-1, epi_eNoN + 1))
    epi_cells = epi_cells[:, 1:]
    epi_global_node_id = epi_mesh.point_data['GlobalNodeID']
    epi_global_cell_id = epi_mesh.cell_data['GlobalElementID']

    # Load the base surface
    base_name = os.path.join(surfaces_dir, surface_names['base'])
    base_mesh = pv.read(base_name)
    base_global_node_id = base_mesh.point_data['GlobalNodeID']

    # Extract the boundary of the epi surface (at the top) to find the apex point
    epi_base_global_node_id = np.intersect1d(epi_global_node_id, base_global_node_id)
    epi_base_nodes = np.where(np.isin(epi_global_node_id, epi_base_global_node_id))[0]

    # # Get normal
    base_normal, base_centroid = get_normal_plane_svd(epi_points[epi_base_nodes, :])

    # Find the index of the apex point of the epi surface
    distance = np.abs(base_normal@(epi_points - base_centroid).T)
    epi_apex_point_index = np.argmax(distance)

    # Find elements containing the apex point
    epi_apex_cell_index = np.where(epi_cells == epi_apex_point_index)[0]

    # Create epi_apex mesh
    submesh_cells = epi_cells[epi_apex_cell_index]
    submesh_xyz = np.zeros([len(np.unique(submesh_cells)), epi_points.shape[1]])
    map_mesh_submesh = np.ones(epi_points.shape[0], dtype=int)*-1
    map_submesh_mesh = np.zeros(submesh_xyz.shape[0], dtype=int)
    child_elems_new = np.zeros(submesh_cells.shape, dtype=int)

    cont = 0
    for e in range(submesh_cells.shape[0]):
        for i in range(submesh_cells.shape[1]):
            if map_mesh_submesh[submesh_cells[e,i]] == -1:
                child_elems_new[e,i] = cont
                submesh_xyz[cont] = epi_points[submesh_cells[e,i]]
                map_mesh_submesh[submesh_cells[e,i]] = cont
                map_submesh_mesh[cont] = submesh_cells[e,i]
                cont += 1
            else:
                child_elems_new[e,i] = map_mesh_submesh[submesh_cells[e,i]]

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

    epi_apex_name = os.path.join(surfaces_dir, surface_names['epi_apex'])
    epi_apex_mesh.save(epi_apex_name)



def runLaplaceSolver(mesh_dir, surfaces_dir, mesh_file, exec_svmultiphysics, template_file, outdir, surface_names):
    xml_template_path = template_file
    out_name = os.path.join(surfaces_dir, "../svFibers_BiV.xml")
    
    with open(xml_template_path, 'r') as svFile:
        xml_content = svFile.read()
    
    # Update mesh file path using regex
    mesh_pattern = r'(<Mesh_file_path>)\s+[^\s<]+[^<]*(</Mesh_file_path>)'
    xml_content = re.sub(mesh_pattern, r'\1 ' + mesh_file + r' \2', xml_content)
    
    # Update face file paths - need to identify which face by checking context
    # Read lines to determine context
    lines = xml_content.split('\n')
    updated_lines = []
    i = 0
    while i < len(lines):
        line = lines[i]
        
        # Check if this line has a face name
        face_match = re.search(r'name="([^"]+)"', line)
        face_name = face_match.group(1) if face_match else None
        
        # Look ahead for Face_file_path
        if face_name and i + 1 < len(lines) and "<Face_file_path>" in lines[i + 1]:
            # Determine which file to use based on face name
            if face_name == "epi":
                new_path = os.path.join(surfaces_dir, surface_names['epi'])
            elif face_name == "epi_top":
                new_path = os.path.join(surfaces_dir, surface_names['base'])
            elif face_name == "epi_apex":
                new_path = os.path.join(surfaces_dir, surface_names['epi_apex'])
            elif face_name == "endo_lv":
                new_path = os.path.join(surfaces_dir, surface_names['endo_lv'])
            elif face_name == "endo_rv":
                new_path = os.path.join(surfaces_dir, surface_names['endo_rv'])
            elif face_name == "mv":
                new_path = os.path.join(surfaces_dir, surface_names['mv'])
            elif face_name == "tv":
                new_path = os.path.join(surfaces_dir, surface_names['tv'])
            elif face_name == "av":
                new_path = os.path.join(surfaces_dir, surface_names['av'])
            elif face_name == "pv":
                new_path = os.path.join(surfaces_dir, surface_names['pv'])
            else:
                new_path = None
            
            if new_path:
                # Add current line
                updated_lines.append(line)
                # Replace the path in the next line
                i += 1
                face_pattern = r'(<Face_file_path>)\s+[^\s<]+[^<]*(</Face_file_path>)'
                updated_line = re.sub(face_pattern, r'\1 ' + new_path + r' \2', lines[i])
                updated_lines.append(updated_line)
                i += 1
                continue
        
        # Add line as-is
        updated_lines.append(line)
        i += 1
    
    xml_content = '\n'.join(updated_lines)
    
    # Update save results folder using regex
    save_pattern = r'(<Save_results_in_folder>)\s+[^\s<]+[^<]*(</Save_results_in_folder>)'
    xml_content = re.sub(save_pattern, r'\1 ' + outdir + r' \2', xml_content)

    with open(out_name, 'w') as svFileNew:
        svFileNew.write(xml_content)

    print("   Running svMultiPhysics solver")
    print(f"   {exec_svmultiphysics + out_name}")
    os.system(exec_svmultiphysics + out_name)

    return outdir + '/results_001.vtu'


def loadLaplaceSolnBayer(fileName):
    '''
    Load a solution to a Laplace-Dirichlet problem from a .vtu file and extract
    the solution and its gradients at the cells.

    ARGS:
    fileName : str
        Path to the .vtu file with the Laplace solution. The solution should be
        defined at the nodes. The Laplace fields should be named as follows:
        - Phi_BiV_EPI: Laplace field for the endocardium
        - Phi_BiV_LV: Laplace field for the left ventricle
        - Phi_BiV_RV: Laplace field for the right ventricle
        - Phi_BiV_AB: Laplace field for the apex to base direction
    '''

    DATASTR1 = 'Phi_BiV_EPI'
    DATASTR2 = 'Phi_BiV_LV'
    DATASTR3 = 'Phi_BiV_RV'
    DATASTR4 = 'Phi_BiV_AB'

    print("   Loading Laplace solution <---   %s" % (fileName))

    # Read mesh with pyvista
    result_mesh = pv.read(fileName)

    # Convert point-data to cell-data (keep point data passed to cells)
    mesh_cells = result_mesh.point_data_to_cell_data()

    print("   Extracting solution and estimating gradients at cells")

    # Get cell centers (Nx3) and scalar cell arrays (N,)
    cPhiEP = np.asarray(mesh_cells.cell_data[DATASTR1])
    cPhiLV = np.asarray(mesh_cells.cell_data[DATASTR2])
    cPhiRV = np.asarray(mesh_cells.cell_data[DATASTR3])
    cPhiAB = np.asarray(mesh_cells.cell_data[DATASTR4])

    # Use pyvista's compute_derivative to get cell gradients
    # compute_derivative will add arrays named '<DATA>_grad' to the cell_data
    grad_mesh = mesh_cells.compute_derivative(scalars=DATASTR1, gradient=True, preference='cell')
    cGPhiEP = np.asarray(grad_mesh.cell_data['gradient'])

    grad_mesh = mesh_cells.compute_derivative(scalars=DATASTR2, gradient=True, preference='cell')
    cGPhiLV = np.asarray(grad_mesh.cell_data['gradient'])

    grad_mesh = mesh_cells.compute_derivative(scalars=DATASTR3, gradient=True, preference='cell')
    cGPhiRV = np.asarray(grad_mesh.cell_data['gradient'])

    grad_mesh = mesh_cells.compute_derivative(scalars=DATASTR4, gradient=True, preference='cell')
    cGPhiAB = np.asarray(grad_mesh.cell_data['gradient'])

    # Use the mesh with cell-data (but without the large scalar arrays) as result_mesh
    mesh_cells.cell_data[DATASTR1 + '_grad'] = cGPhiEP
    mesh_cells.cell_data[DATASTR2 + '_grad'] = cGPhiLV
    mesh_cells.cell_data[DATASTR3 + '_grad'] = cGPhiRV
    mesh_cells.cell_data[DATASTR4 + '_grad'] = cGPhiAB

    return mesh_cells, cPhiEP, cPhiLV, cPhiRV, cPhiAB, \
        cGPhiEP, cGPhiLV, cGPhiRV, cGPhiAB


def bislerp(Q1, Q2, interp_func):
    """
    Vectorized spherical interpolation between batches of rotation matrices.
    Q1, Q2: (N, 3, 3)
    interp_func: (N,) values in [0,1]
    Returns Q: (N, 3, 3)
    Notes:
      - Uses wxyz quaternion convention internally, matching quat2rot below.
      - Avoids per-element Python/Scipy objects for performance.
    """
    def rotm_to_quat_batch(R):
        # R: (N,3,3) -> q: (N,4) [w,x,y,z]
        t = np.einsum('nii->n', R)  # trace
        q = np.zeros((R.shape[0], 4), dtype=float)
        # Branch where trace is positive
        mask_t = t > 0.0
        if np.any(mask_t):
            S = np.sqrt(t[mask_t] + 1.0) * 2.0
            q[mask_t, 0] = 0.25 * S
            q[mask_t, 1] = (R[mask_t, 2, 1] - R[mask_t, 1, 2]) / S
            q[mask_t, 2] = (R[mask_t, 0, 2] - R[mask_t, 2, 0]) / S
            q[mask_t, 3] = (R[mask_t, 1, 0] - R[mask_t, 0, 1]) / S
        # For remaining, choose major diagonal
        mask_f = ~mask_t
        if np.any(mask_f):
            Rf = R[mask_f]
            m00 = Rf[:, 0, 0]
            m11 = Rf[:, 1, 1]
            m22 = Rf[:, 2, 2]
            idx = np.argmax(np.stack([m00, m11, m22], axis=1), axis=1)
            mf_idx = np.nonzero(mask_f)[0]
            # Case idx==0
            m0 = idx == 0
            if np.any(m0):
                S = np.sqrt(1.0 + Rf[m0, 0, 0] - Rf[m0, 1, 1] - Rf[m0, 2, 2]) * 2.0
                rows = mf_idx[m0]
                q[rows, 0] = (Rf[m0, 2, 1] - Rf[m0, 1, 2]) / S
                q[rows, 1] = 0.25 * S
                q[rows, 2] = (Rf[m0, 0, 1] + Rf[m0, 1, 0]) / S
                q[rows, 3] = (Rf[m0, 0, 2] + Rf[m0, 2, 0]) / S
            # Case idx==1
            m1 = idx == 1
            if np.any(m1):
                S = np.sqrt(1.0 + Rf[m1, 1, 1] - Rf[m1, 0, 0] - Rf[m1, 2, 2]) * 2.0
                rows = mf_idx[m1]
                q[rows, 0] = (Rf[m1, 0, 2] - Rf[m1, 2, 0]) / S
                q[rows, 1] = (Rf[m1, 0, 1] + Rf[m1, 1, 0]) / S
                q[rows, 2] = 0.25 * S
                q[rows, 3] = (Rf[m1, 1, 2] + Rf[m1, 2, 1]) / S
            # Case idx==2
            m2 = idx == 2
            if np.any(m2):
                S = np.sqrt(1.0 + Rf[m2, 2, 2] - Rf[m2, 0, 0] - Rf[m2, 1, 1]) * 2.0
                rows = mf_idx[m2]
                q[rows, 0] = (Rf[m2, 1, 0] - Rf[m2, 0, 1]) / S
                q[rows, 1] = (Rf[m2, 0, 2] + Rf[m2, 2, 0]) / S
                q[rows, 2] = (Rf[m2, 1, 2] + Rf[m2, 2, 1]) / S
                q[rows, 3] = 0.25 * S

        # Normalize for numerical safety
        q /= np.linalg.norm(q, axis=1, keepdims=True)
        return q

    def quat_to_rotm_batch(q):
        # q: (N,4) [w,x,y,z] -> R: (N,3,3)
        w = q[:, 0]
        x = q[:, 1]
        y = q[:, 2]
        z = q[:, 3]
        x2 = x * x
        y2 = y * y
        z2 = z * z
        wx = w * x
        wy = w * y
        wz = w * z
        xy = x * y
        xz = x * z
        yz = y * z
        R = np.zeros((q.shape[0], 3, 3), dtype=float)
        R[:, 0, 0] = 1.0 - 2.0 * y2 - 2.0 * z2
        R[:, 1, 0] = 2.0 * xy + 2.0 * wz
        R[:, 2, 0] = 2.0 * xz - 2.0 * wy
        R[:, 0, 1] = 2.0 * xy - 2.0 * wz
        R[:, 1, 1] = 1.0 - 2.0 * x2 - 2.0 * z2
        R[:, 2, 1] = 2.0 * yz + 2.0 * wx
        R[:, 0, 2] = 2.0 * xz + 2.0 * wy
        R[:, 1, 2] = 2.0 * yz - 2.0 * wx
        R[:, 2, 2] = 1.0 - 2.0 * x2 - 2.0 * y2
        return R

    # Prepare inputs
    t = np.clip(np.asarray(interp_func, dtype=float), 0.0, 1.0)
    q1 = rotm_to_quat_batch(np.asarray(Q1, dtype=float))
    q2 = rotm_to_quat_batch(np.asarray(Q2, dtype=float))

    # Ensure shortest path on the unit 4-sphere
    dot = np.sum(q1 * q2, axis=1)
    neg_mask = dot < 0.0
    if np.any(neg_mask):
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
        # Nearly identical orientations: perform linear interpolation and normalize later
        tl = t[lin_mask][:, None]
        q[lin_mask] = (1.0 - tl) * q1[lin_mask] + tl * q2[lin_mask]

    # Normalize and convert back to rotation matrices
    q /= np.linalg.norm(q, axis=1, keepdims=True)
    return quat_to_rotm_batch(q)

def axis(u, v):
    """
    u, v: (nelems, 3)
    return Q: (nelems, 3, 3) where columns are [e0 (circ), e1 (long), e2 (trans)] per element
    """
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    ne = u.shape[0]

    # e1 = normalize rows of u
    e1 = normalize(u)

    # e2 = v - proj_{e1}(v)
    proj = np.sum(e1 * v, axis=1)[:, None] * e1
    e2 = v - proj
    e2 = normalize(e2)

    # e0 = cross(e1, e2) normalized
    e0 = np.cross(e1, e2, axisa=1, axisb=1)
    e0 = normalize(e0)

    Q = np.zeros((ne, 3, 3), dtype=float)
    Q[:, :, 0] = e0
    Q[:, :, 1] = e1
    Q[:, :, 2] = e2

    return Q

def orient(Q, alpha, beta):
    """
    Given an orthogonal matrix Q (ne,3,3), rotate each Q[i] by alpha[i] about
    the z-axis and then by beta[i] about the x-axis. alpha and beta are arrays 
    of length ne.
    """
    Q = np.asarray(Q, dtype=float)
    ne = Q.shape[0]

    ca = np.cos(alpha)
    sa = np.sin(alpha)
    cb = np.cos(beta)
    sb = np.sin(beta)

    # Rotation about z (Ra) and x (Rb) for each element
    Ra = np.zeros((ne, 3, 3), dtype=float)
    Ra[:, 0, 0] = ca
    Ra[:, 0, 1] = -sa
    Ra[:, 1, 0] = sa
    Ra[:, 1, 1] = ca
    Ra[:, 2, 2] = 1.0

    # Rb = np.zeros((ne, 3, 3), dtype=float)
    # Rb[:, 0, 0] = 1.0
    # Rb[:, 1, 1] = cb
    # Rb[:, 1, 2] = sb
    # Rb[:, 2, 1] = -sb
    # Rb[:, 2, 2] = cb

    # Rb = np.zeros((ne, 3, 3), dtype=float)
    # Rb[:, 0, 0] = 1.0
    # Rb[:, 1, 1] = cb
    # Rb[:, 1, 2] = -sb
    # Rb[:, 2, 1] = sb
    # Rb[:, 2, 2] = cb

    Rb = np.zeros((ne, 3, 3), dtype=float)
    Rb[:, 0, 0] = cb
    Rb[:, 0, 2] = sb
    Rb[:, 1, 1] = 1.0
    Rb[:, 2, 0] = -sb
    Rb[:, 2, 2] = cb

    # Compose rotations and apply to Q per element
    RaRb = np.einsum('nij,njk->nik', Ra, Rb)
    Qt = np.einsum('nij,njk->nik', Q, RaRb)

    return Qt

def getFiberDirectionsBayer(Phi_EPI, Phi_LV, Phi_RV,
                                 gPhi_EPI, gPhi_LV, gPhi_RV, gPhi_AB, 
                                 params, intermediate=False):
    '''
    Compute the fiber directions at the center of each cell
    '''

    # Unpack parameters
    ALFA_END = np.deg2rad(params["ALFA_END"])
    ALFA_EPI = np.deg2rad(params["ALFA_EPI"])
    BETA_END = np.deg2rad(params["BETA_END"])
    BETA_EPI = np.deg2rad(params["BETA_EPI"])

    print("   Computing fiber directions at cells")

    d = Phi_RV / (Phi_LV + Phi_RV)
    alfaS = ALFA_END * (1 - d) - ALFA_END * d
    betaS = BETA_END * (1 - d) - BETA_END * d
    alfaW = ALFA_END * (1 - Phi_EPI) + ALFA_EPI * Phi_EPI
    betaW = BETA_END * (1 - Phi_EPI) + BETA_EPI * Phi_EPI

    Q_LV0 = axis(gPhi_AB, -gPhi_LV)
    Q_LV = orient(Q_LV0, alfaS, betaS)
    Q_RV0 = axis(gPhi_AB, gPhi_RV)  # Note that gPhi_RV points the other way
    Q_RV = orient(Q_RV0, alfaS, -betaS)  # Therefore, we need a minus in betaS

    Q_END = bislerp(Q_LV, Q_RV, d)
    Q_END[d > 0.5,:,0] = -Q_END[d > 0.5,:,0]
    Q_END[d > 0.5,:,2] = -Q_END[d > 0.5,:,2]

    Q_EPI0 = axis(gPhi_AB, gPhi_EPI)
    Q_EPI = orient(Q_EPI0, alfaW, betaW)

    FST = bislerp(Q_END, Q_EPI, Phi_EPI)

    F = FST[:, :, 0]
    S = FST[:, :, 1]
    T = FST[:, :, 2]

    if intermediate:
        return F, S, T, Q_LV[:,:,0], Q_RV[:,:,0], Q_END[:,:,0], Q_EPI[:,:,0]

    return F, S, T


def get_alpha_beta_angles_Bayer(F, Phi_EPI, Phi_LV, Phi_RV,
                                 gPhi_EPI, gPhi_LV, gPhi_RV, gPhi_AB, 
                                 params):
    '''
    Sanity check routine
    Compute alpha and beta angles at cells given fiber directions F and Laplace gradients

    '''

    ALFA_END = np.deg2rad(params["ALFA_END"])
    ALFA_EPI = np.deg2rad(params["ALFA_EPI"])

    d = Phi_RV / (Phi_LV + Phi_RV)
    alfaS = ALFA_END * (1 - d) - ALFA_END * d
    alfaW = ALFA_END * (1 - Phi_EPI) + ALFA_EPI * Phi_EPI

    # Alpha angle 
    Q_LV = axis(gPhi_AB, -gPhi_LV)          # (N,3,3)
    Q_RV = axis(gPhi_AB, gPhi_RV)            # (N,3,3)
    Q = np.copy(Q_LV)
    Q[d > 0.5,:,0] = -Q_RV[d > 0.5,:,0]
    Q[d > 0.5,:,2] = -Q_RV[d > 0.5,:,2]
    C = Q[:, :, 0]                       # (N,3)
    L = Q[:, :, 1]                       # (N,3)

    # Angle in radians between F and circumferential vector
    # alpha is positive in the direction of the longitudinal vector
    cosang = np.clip(np.sum(F * C, axis=1), -1.0, 1.0)
    sinang = np.clip(np.sum(F * L, axis=1), -1.0, 1.0)
    alpha_angle = np.sign(sinang) * np.arccos(np.abs(cosang))

    # Beta angle
    Q_LV = orient(axis(gPhi_AB, -gPhi_LV), 
                    alfaS,
                    0)
    Q_RV = orient(axis(gPhi_AB, gPhi_RV), 
                    alfaS,
                    0)
    
    Q_END = bislerp(Q_LV, Q_RV, d)
    Q_END[d > 0.5,:,0] = -Q_END[d > 0.5,:,0]
    Q_END[d > 0.5,:,2] = -Q_END[d > 0.5,:,2]

    Q_EPI0 = axis(gPhi_AB, gPhi_EPI)
    Q_EPI = orient(Q_EPI0, alfaW, 0.0)

    Q = bislerp(Q_END, Q_EPI, Phi_EPI)
    Cr = Q[:, :, 0]
    Tr = Q[:, :, 2]

    # Angle in radians between F and rotated circumferential vector
    # beta is negative in the direction of the transmural vector
    cosang = np.clip(np.sum(F * Cr, axis=1), -1.0, 1.0)
    sinang = np.clip(np.sum(F * Tr, axis=1), -1.0, 1.0)
    beta_angle = - np.sign(sinang) * np.arccos(np.abs(cosang)) # Note the minus sign to match definition

    return np.rad2deg(alpha_angle), np.rad2deg(beta_angle), C, Cr

def generate_fibers_BiV_Bayer_cells(outdir, laplace_results_file, params, return_angles=False, return_intermediate=False):
    '''
    Generate fiber directions on a truncated BiV ventricular geometry using the
    Laplace-Dirichlet rule-based method of Bayer et al. 2012

    ARGS:
    laplace_results_file : str
        Path to the .vtu mesh with Laplace fields defined at nodes
    params : dict
        Dictionary of parameters for fiber generation
    '''
    
    t1 = time.time()
    print("========================================================")

    # Load Laplace solution    
    result_mesh, Phi_EPI, Phi_LV, Phi_RV, Phi_AB, \
    gPhi_EPI, gPhi_LV, gPhi_RV, gPhi_AB = loadLaplaceSolnBayer(laplace_results_file)


    # Write the fiber directions to a vtu files
    output_mesh = copy.deepcopy(result_mesh)
    # Ensure only FIB_DIR is present
    for k in list(output_mesh.cell_data.keys()):
        output_mesh.cell_data.remove(k)
    for k in list(output_mesh.point_data.keys()):
        output_mesh.point_data.remove(k)
        
    # Generate fiber directions
    out = getFiberDirectionsBayer(Phi_EPI, Phi_LV, Phi_RV,
                                 gPhi_EPI, gPhi_LV, gPhi_RV, gPhi_AB, 
                                 params, intermediate=return_intermediate)
    
    if return_intermediate:
        F, S, T, eC_LV, eC_RV, eC_END, eC_EPI = out
        result_mesh.cell_data['F'] = F
        result_mesh.cell_data['S'] = S
        result_mesh.cell_data['T'] = T
        result_mesh.cell_data['eC_LV'] = eC_LV
        result_mesh.cell_data['eC_RV'] = eC_RV
        result_mesh.cell_data['eC_END'] = eC_END
        result_mesh.cell_data['eC_EPI'] = eC_EPI
    else:
        F, S, T = out
        result_mesh.cell_data['F'] = F
        result_mesh.cell_data['S'] = S
        result_mesh.cell_data['T'] = T

    print("   Writing domains and fibers to VTK data structure")


    fname1 = os.path.join(outdir, "fibersLong.vtu")
    print("   Writing to vtu file   --->   %s" % (fname1))
    output_mesh.cell_data.set_array(F, 'FIB_DIR')
    output_mesh.save(fname1)

    fname1 = os.path.join(outdir, "fibersSheet.vtu")
    print("   Writing to vtu file   --->   %s" % (fname1))
    output_mesh.cell_data.set_array(T, 'FIB_DIR')
    output_mesh.save(fname1)

    fname1 = os.path.join(outdir, "fibersNormal.vtu")
    print("   Writing to vtu file   --->   %s" % (fname1))
    output_mesh.cell_data.set_array(S, 'FIB_DIR')
    output_mesh.save(fname1)

    t2 = time.time()
    print('\n   Total time: %.3fs' % (t2-t1))
    print("========================================================")

    if return_angles:
        alpha_angle, beta_angle, eC, eCr = get_alpha_beta_angles_Bayer(F, Phi_EPI, Phi_LV, Phi_RV,
                                 gPhi_EPI, gPhi_LV, gPhi_RV, gPhi_AB, 
                                 params)
        result_mesh.cell_data['Alpha_Angle'] = alpha_angle
        result_mesh.cell_data['Beta_Angle'] = beta_angle
        result_mesh.cell_data['eC'] = eC
        result_mesh.cell_data['eCr'] = eCr
    

    return result_mesh



def loadLaplaceSolnDoste(fileName):
    '''
    Load a solution to a Laplace-Dirichlet problem from a .vtu file and extract
    the solution and its gradients at the cells.

    ARGS:
    fileName : str
        Path to the .vtu file with the Laplace solution. The solution should be
        defined at the nodes. 

    Returns:
    lap : dict
        Dictionary of Laplace solution at cells
    grad : dict
        Dictionary of gradients at cells
    '''

    varnames = ['Trans_BiV', 'Long_AV', 'Long_MV', 'Long_PV', 'Long_TV', 'Weight_LV', 
                'Weight_RV', 'Trans_EPI', 'Trans_LV', 'Trans_RV']

    print("   Loading Laplace solution <---   %s" % (fileName))

    # Read mesh with pyvista
    result_mesh = pv.read(fileName)

    # Convert point-data to cell-data (keep point data passed to cells)
    mesh_cells = result_mesh.point_data_to_cell_data()

    print("   Extracting solution and estimating gradients at cells")

    # Map VTU array names to internal keys expected by downstream code
    name_map = {
        'Trans_BiV': 'ven_trans',
        'Long_AV': 'lv_av_long',
        'Long_MV': 'lv_mv_long',
        'Long_PV': 'rv_pv_long',
        'Long_TV': 'rv_tv_long',
        'Weight_LV': 'lv_weight',
        'Weight_RV': 'rv_weight',
        'Weight_RV_OP': 'rv_op_weight',
        'Trans_EPI': 'epi_trans',
        'Trans_LV': 'lv_trans',
        'Trans_RV': 'rv_trans',
    }

    lap = {}
    grad = {}

    for vname in varnames:
        if vname not in mesh_cells.cell_data:
            print(f"   Warning: '{vname}' not found in cell_data; skipping")
            continue
        key = name_map[vname]

        # Cell-centered Laplace values
        lap[key] = np.asarray(mesh_cells.cell_data[vname])

        # Cell-centered gradients via PyVista
        gmesh = mesh_cells.compute_derivative(scalars=vname, gradient=True, preference='cell')
        grad[key] = np.asarray(gmesh.cell_data['gradient'])

    return result_mesh, lap, grad


def compute_basis_vectors(lap, grad):
    # LV
    # longitudinal
    lv_glong = grad['lv_mv_long']*lap['lv_weight'][:,None] + grad['lv_av_long']*(1 - lap['lv_weight'][:,None])
    eL_lv = normalize(lv_glong)

    # transmural
    lv_gtrans = grad['lv_trans'] - (eL_lv*grad['lv_trans'])*eL_lv
    eT_lv = normalize(lv_gtrans)

    # circumferential
    eC_lv = np.cross(eL_lv, eT_lv, axisa=1, axisb=1)
    eC_lv = normalize(eC_lv)

    # Ensuring orthogonality
    eT_lv = np.cross(eC_lv, eL_lv, axisa=1, axisb=1)
    eT_lv = normalize(eT_lv)

    # RV
    # longitudinal
    rv_glong = grad['rv_tv_long']*lap['rv_weight'][:,None]  + grad['rv_pv_long']*(1 - lap['rv_weight'][:,None] )
    eL_rv = normalize(rv_glong)

    # transmural
    rv_gtrans = grad['rv_trans'] - (eL_rv*grad['rv_trans'])*eL_rv
    eT_rv = normalize(rv_gtrans)

    # circumferential
    eC_rv = np.cross(eL_rv, eT_rv, axisa=1, axisb=1)
    eC_rv = normalize(eC_rv)

    # Ensuring orthogonality
    eT_rv = np.cross(eC_rv, eL_rv, axisa=1, axisb=1)
    eT_rv = normalize(eT_rv)

    # Write out global circumferential vector
    eC = eC_rv*(1-lap['ven_trans'][:,None]) + eC_lv*lap['ven_trans'][:,None]
    eC = normalize(eC)

    basis = {'eC_lv': eC_lv,
                    'eT_lv': eT_lv,
                    'eL_lv': eL_lv,
                    'eC_rv': eC_rv,
                    'eT_rv': eT_rv,
                    'eL_rv': eL_rv,
                    'eC': eC}

    return basis


def redistribute_weight(weight, up, low, strategy='centre'):
    new_weight = weight.copy()

    if strategy == 'flip':
        # Shift all weights
        new_mean = 1 - np.mean(weight)
        shift = new_mean - np.mean(weight)
        new_weight = new_weight + shift

        # Cut off values outside of range 0 - 1
        new_weight[new_weight > 1] = 1
        new_weight[new_weight < 0] = 0

        # Redistribute new tail
        new_weight = (new_weight - np.min(new_weight)) / (np.max(new_weight) - np.min(new_weight))
        tmp = new_weight.copy()

        if shift > 0:
            tmp[tmp >= new_mean] = np.nan
            tmp = (tmp - np.nanmin(tmp)) / (new_mean - np.nanmin(tmp))
        elif shift < 0:
            tmp[tmp <= new_mean] = np.nan
            tmp = (tmp - new_mean) / (np.nanmax(tmp) - new_mean)

        tmp[np.isnan(tmp)] = new_weight[np.isnan(tmp)]
        new_weight = tmp

    else:  # cut off tails so that weights are centered
        # Find upper and lower limits
        upper_lim = np.quantile(weight, up)
        while upper_lim == 0:
            print('Upper limit is 0, increasing upper limit')
            up += 0.1
            upper_lim = np.quantile(weight, up)
        
        lower_lim = np.quantile(weight, low)

        # Set upper and lower values to limits
        new_weight[new_weight > upper_lim] = upper_lim
        new_weight[new_weight < lower_lim] = lower_lim

        # Redistribute/normalize values
        new_weight = (new_weight - np.min(new_weight)) / (np.max(new_weight) - np.min(new_weight))

    return new_weight


def compute_alpha_beta_angles(lap, params):
    # Modify weights so the effect of outflow tracts is localized
    lv_weight = redistribute_weight(lap['lv_weight'], 0.7, 0.01)
    rv_weight = redistribute_weight(lap['rv_weight'], 0.1, 0.001)

    # LV
    alpha_lv_endo_long = params['AENDOLV'] * lv_weight + params['AOTENDOLV'] * (1 - lv_weight)  # Endo
    alpha_lv_epi_long = params['AEPILV'] * lv_weight + params['AOTEPILV'] * (1 - lv_weight)
    
    alpha_wall_lv = alpha_lv_endo_long * (1 - lap['epi_trans']) + alpha_lv_epi_long * lap['epi_trans']
    beta_wall_lv = (params['BENDOLV'] * (1 - lap['epi_trans']) + params['BEPILV'] * lap['epi_trans']) * lv_weight

    # RV
    alpha_rv_endo_long = params['AENDORV'] * rv_weight + params['AOTENDORV'] * (1 - rv_weight)
    alpha_rv_epi_long = params['AEPIRV'] * rv_weight + params['AOTEPIRV'] * (1 - rv_weight)

    alpha_wall_rv = alpha_rv_endo_long * (1 - lap['epi_trans']) + alpha_rv_epi_long * lap['epi_trans']
    beta_wall_rv = (params['BENDORV'] * (1 - lap['epi_trans']) + params['BEPIRV'] * lap['epi_trans']) * rv_weight

    # Septum
    sep = np.abs(lap['ven_trans'] - 0.5)
    sep = (sep - np.min(sep)) / (np.max(sep) - np.min(sep))
    alpha_septum = (alpha_lv_endo_long * sep * lap['lv_trans']) + (alpha_rv_endo_long * sep * lap['rv_trans'])
    beta_septum = (params['BENDOLV'] * lap['lv_trans'] * lv_weight) + (params['BENDORV'] * lap['rv_trans'] * rv_weight)

    angles = {'alpha_lv_endo_long': alpha_lv_endo_long,
            'alpha_lv_epi_long': alpha_lv_epi_long,
            'alpha_wall_lv': alpha_wall_lv,
            'beta_wall_lv': beta_wall_lv,
            'alpha_rv_endo_long': alpha_rv_endo_long,
            'alpha_rv_epi_long': alpha_rv_epi_long,
            'alpha_wall_rv': alpha_wall_rv,
            'beta_wall_rv': beta_wall_rv,
            'alpha_septum': alpha_septum,
            'beta_septum': beta_septum
            }

    return angles


def rotate_basis(eC, eL, eT, alpha, beta):
    eC = normalize(eC)
    eT = normalize(eT)
    eL = normalize(eL)

    # Matrix of directional vectors
    Q = np.stack([eC, eL, eT], axis=-1)
    Q = np.transpose(Q, (2, 1, 0))

    # Create rotation matrix - from Doste code
    axis = eT
    R = np.array([[np.cos(alpha) + (axis[:, 0]**2)*(1 - np.cos(alpha)), axis[:,0] * axis[:,1]*(1 - np.cos(alpha)) - axis[:,2]*np.sin(alpha), axis[:,0]*axis[:,2]*(1 - np.cos(alpha)) + axis[:,1]*np.sin(alpha)],
                            [axis[:,1]*axis[:,0]*(1 - np.cos(alpha)) + axis[:,2]*np.sin(alpha), np.cos(alpha) + (axis[:,1]**2)*(1 - np.cos(alpha)), axis[:,1]*axis[:, 2]*(1 - np.cos(alpha)) - axis[:, 0]*np.sin(alpha)],
                            [axis[:,2]*axis[:,0]*(1 - np.cos(alpha)) - axis[:,1]*np.sin(alpha), axis[:,2]*axis[:,1]*(1 - np.cos(alpha)) + axis[:, 0]*np.sin(alpha), np.cos(alpha)+(axis[:, 2]**2)*(1 - np.cos(alpha))]])

    # Rotate the circumferential direction around the transmural direction
    QX = np.zeros_like(R)
    for i in range(len(eC)):
        QX[:, :, i] = np.matmul(Q[:, :, i], R[:, :, i])

    # Second rotation (beta) about QX
    axis2 = QX[1, :, :].T
    R2 = np.array([
        [np.cos(beta) + (axis2[:,0]**2)*(1 - np.cos(beta)), axis2[:,0]*axis2[:, 1]*(1 - np.cos(beta)) - axis2[:,2] * np.sin(beta), axis2[:,0] * axis2[:,2] * (1 - np.cos(beta)) + axis2[:,1] * np.sin(beta)],
        [axis2[:,1] * axis2[:,0]*(1 - np.cos(beta)) + axis2[:,2]*np.sin(beta), np.cos(beta) + (axis2[:,1]**2)*(1 - np.cos(beta)), axis2[:,1] * axis2[:,2] * (1 - np.cos(beta)) - axis2[:,0] * np.sin(beta)],
        [axis2[:,2] * axis2[:,0]*(1 - np.cos(beta)) - axis2[:,1]*np.sin(beta), axis2[:, 2] * axis2[:,1] * (1 - np.cos(beta)) + axis2[:,0] * np.sin(beta), np.cos(beta) + (axis2[:,2]**2) * (1 - np.cos(beta))]
    ])

    QX2 = np.zeros((R.shape[2], 3, 3), dtype=float)
    for i in range(len(eC)):
        QX2[i] = np.matmul(QX[:, :, i], R2[:, :, i]).T

    return QX2


def compute_local_basis(basis, angles):
    Qlv_septum = rotate_basis(basis['eC_lv'], basis['eL_lv'], basis['eT_lv'], angles['alpha_septum'], angles['beta_septum'])
    Qrv_septum = rotate_basis(basis['eC_rv'], basis['eL_rv'], basis['eT_rv'], angles['alpha_septum'], angles['beta_septum'])
    Qlv_epi = rotate_basis(basis['eC_lv'], basis['eL_lv'], basis['eT_lv'], angles['alpha_wall_lv'], angles['beta_wall_lv'])
    Qrv_epi = rotate_basis(basis['eC_rv'], basis['eL_rv'], basis['eT_rv'], angles['alpha_wall_rv'], angles['beta_wall_rv'])

    local_basis = {'Qlv_septum': Qlv_septum,
                    'Qrv_septum': Qrv_septum,
                    'Qlv_epi': Qlv_epi,
                    'Qrv_epi': Qrv_epi,
                    }

    return local_basis


def interpolate_local_basis(lap, local_basis):

    epi_trans = lap['epi_trans']

    Qrv_septum = local_basis['Qrv_septum']
    Qlv_septum = local_basis['Qlv_septum']
    Qrv_epi = local_basis['Qrv_epi']
    Qlv_epi = local_basis['Qlv_epi']

    Qepi = bislerp(Qrv_epi, Qlv_epi, lap['ven_trans'])
    Qendo = bislerp(Qrv_septum, Qlv_septum, lap['ven_trans'])
    Q = bislerp(Qendo, Qepi, epi_trans)

    return Q, Qepi


def getFiberDirectionsDoste(lap, grad, params, intermediate=False):
    # Convert parameters from degrees to radians
    for key in params:
        params[key] = np.deg2rad(params[key])

    print('Computing basis vectors')
    basis = compute_basis_vectors(lap, grad)

    print('Computing angles')
    angles = compute_alpha_beta_angles(lap, params)

    print('Computing local basis')
    local_basis = compute_local_basis(basis, angles)

    print('Interpolating basis')
    Q, Qepi = interpolate_local_basis(lap, local_basis)

    print('Done!')
    f = Q[:, :, 0]
    s = Q[:, :, 1]
    n = Q[:, :, 2]

    if intermediate:
        return f, s, n, basis, angles, local_basis, Qepi[:,:,0]

    return f, s, n


def generate_fibers_BiV_Doste_cells(outdir, laplace_results_file, params, return_angles=False, return_intermediate=False):
    '''
    Generate fiber directions on a truncated BiV ventricular geometry using the
    Laplace-Dirichlet rule-based method of Bayer et al. 2012

    ARGS:
    laplace_results_file : str
        Path to the .vtu mesh with Laplace fields defined at nodes
    params : dict
        Dictionary of parameters for fiber generation
    '''
    
    t1 = time.time()
    print("========================================================")

    # Load Laplace solution    
    result_mesh,lap, grad = loadLaplaceSolnDoste(laplace_results_file)

    # Write the fiber directions to a vtu files
    output_mesh = copy.deepcopy(result_mesh)
    # Ensure only FIB_DIR is present
    for k in list(output_mesh.cell_data.keys()):
        output_mesh.cell_data.remove(k)
    for k in list(output_mesh.point_data.keys()):
        output_mesh.point_data.remove(k)

    # Generate fiber directions
    out = getFiberDirectionsDoste(lap, grad, params, intermediate=return_intermediate)

    if return_intermediate:
        F, S, T, basis, angles, local_basis, Qepi = out
        result_mesh.cell_data['F'] = F
        result_mesh.cell_data['S'] = S
        result_mesh.cell_data['T'] = T
        result_mesh.cell_data['Qepi'] = Qepi
        for k, v in basis.items():
            result_mesh.cell_data[k] = v
        for k, v in angles.items():
            result_mesh.cell_data[k] = v
        for k, v in local_basis.items():
            # Flatten local basis matrices to store
            flattened = v[:, :, 0]
            result_mesh.cell_data[k] = flattened
    else:
        F, S, T = out
        result_mesh.cell_data['F'] = F
        result_mesh.cell_data['S'] = S
        result_mesh.cell_data['T'] = T

    print("   Writing domains and fibers to VTK data structure")


    fname1 = os.path.join(outdir, "fibersLong.vtu")
    print("   Writing to vtu file   --->   %s" % (fname1))
    output_mesh.cell_data.set_array(F, 'FIB_DIR')
    output_mesh.save(fname1)

    fname1 = os.path.join(outdir, "fibersSheet.vtu")
    print("   Writing to vtu file   --->   %s" % (fname1))
    output_mesh.cell_data.set_array(T, 'FIB_DIR')
    output_mesh.save(fname1)

    fname1 = os.path.join(outdir, "fibersNormal.vtu")
    print("   Writing to vtu file   --->   %s" % (fname1))
    output_mesh.cell_data.set_array(S, 'FIB_DIR')
    output_mesh.save(fname1)

    t2 = time.time()
    print('\n   Total time: %.3fs' % (t2-t1))
    print("========================================================")

    if return_angles:
        alpha_angle, beta_angle, eC, eCr = get_alpha_beta_angles_Doste(F, lap, grad, params)
        result_mesh.cell_data['Alpha_Angle'] = alpha_angle
        result_mesh.cell_data['Beta_Angle'] = beta_angle
        result_mesh.cell_data['eC'] = eC
        result_mesh.cell_data['eCr'] = eCr


    return result_mesh



def get_alpha_beta_angles_Doste(F, lap, grad, params):
    '''
    Sanity check routine for Doste-based fibers.
    Compute alpha and beta angles at cells given fiber directions F and the
    Laplace/basis fields used by the Doste method.

    Returns:
      alpha_angle_deg, beta_angle_deg, eC_ref, Cr_ref
        - alpha, beta in degrees
        - eC_ref: reference circumferential vector (before rotations)
        - Cr_ref: circumferential vector after applying only alpha rotation
    '''

    # Reconstruct base vectors used by Doste
    basis = compute_basis_vectors(lap, grad)
    eC_global = basis['eC']  # blended circumferential

    # Build a blended longitudinal direction for alpha sign (LV/RV mix)
    ven = lap['ven_trans'][:, None]
    eL_blend = normalize(basis['eL_rv'] * (1.0 - ven) + basis['eL_lv'] * ven)

    # Alpha: signed angle between F and circumferential in the tangent plane,
    # sign taken along longitudinal direction (consistent with Bayer routine, 
    # but negative because the longitudinal direction is opposite to the Doste paper)
    cos_a = np.clip(np.sum(F * eC_global, axis=1), -1.0, 1.0)
    sin_a = np.clip(np.sum(F * eL_blend, axis=1), -1.0, 1.0)
    alpha_angle = -np.sign(sin_a) * np.arccos(np.abs(cos_a))

    # Build reference frame after ONLY alpha rotation (beta = 0)
    angles = compute_alpha_beta_angles(lap, params)  # radians
    Qlv_septum_a = rotate_basis(basis['eC_lv'], basis['eL_lv'], basis['eT_lv'],
                                angles['alpha_septum'], 0.0)
    Qrv_septum_a = rotate_basis(basis['eC_rv'], basis['eL_rv'], basis['eT_rv'],
                                angles['alpha_septum'], 0.0)
    Qlv_epi_a = rotate_basis(basis['eC_lv'], basis['eL_lv'], basis['eT_lv'],
                             angles['alpha_wall_lv'], 0.0)
    Qrv_epi_a = rotate_basis(basis['eC_rv'], basis['eL_rv'], basis['eT_rv'],
                             angles['alpha_wall_rv'], 0.0)

    Qepi_a = bislerp(Qrv_epi_a, Qlv_epi_a, lap['ven_trans'])
    Qendo_a = bislerp(Qrv_septum_a, Qlv_septum_a, lap['ven_trans'])
    Qa = bislerp(Qendo_a, Qepi_a, lap['epi_trans'])

    Cr = Qa[:, :, 0]  # circumferential after alpha-only rotation
    Tr = Qa[:, :, 2]  # transmural after alpha-only rotation

    # Beta: signed angle between F and Cr, sign w.r.t. Tr (negative by convention)
    cos_b = np.clip(np.sum(F * Cr, axis=1), -1.0, 1.0)
    sin_b = np.clip(np.sum(F * Tr, axis=1), -1.0, 1.0)
    beta_angle = np.sign(sin_b) * np.arccos(np.abs(cos_b))

    return np.rad2deg(alpha_angle), np.rad2deg(beta_angle), eC_global, Cr

