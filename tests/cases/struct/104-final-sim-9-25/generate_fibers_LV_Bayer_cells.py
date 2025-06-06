'''
Generates myofiber directions on a truncated left ventricular geometry using 
the Laplace-Dirichlet rule-based method of Bayer et al. 2012
(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3518842/), modified for an LV only
geometry

The fiber directions are computed at the center of each cell of the mesh. The 
Laplace solution is computed using the svFSI solver. The fiber directions are 
written to a vtu file.

Written by Lei Shi. Modifications by Aaron Brown (2024)

Input:
    1. Mesh file: mesh-complete.mesh.vtu
    2. Mesh surfaces: top.vtp, epi.vtp, endo.vtp
    3. svFSI input file: svFSI.inp

Output:
    1. Fiber directions: fibersLong.vtu, fibersSheet.vtu, fibersNormal.vtu
       in same folder as mesh files
    2. Laplace solution: result_020.vtu in lps_thickness folder
    3. Epi apex and epi mid surfaces: epi_apex.vtp, epi_mid.vtp in mesh-surfaces folder
       Epi apex contains all cells that contain the apex point of the epi surface,
       which is the point farthest from the top surface. Epi mid contains the 
       remaining cells of the epi surface.
'''

import os
import sys
import vtk
import time
import numpy as np
from vtkmodules.util import numpy_support as vtknp
import copy
import pyvista as pv

#----------------------------------------------------------------------
def distance_to_surface(point, surface_normal, surface_point):
    '''
    Calculate the distance from a point to a point on a surface
    '''
    # Ensure that the input vectors are NumPy arrays
    point = np.array(point)
    surface_normal = np.array(surface_normal)
    surface_point = np.array(surface_point)

    # Calculate the vector from the surface point to the target point
    v = point - surface_point

    # Calculate the dot product between the vector v and the surface normal
    dot_product = np.dot(v, surface_normal)

    # Calculate the absolute value of the dot product to get the distance
    distance = abs(dot_product)

    return distance

#----------------------------------------------------------------------

#---------------------------------------------------------------------
def normalize(u):
    '''
    Calculate the normalized vector of a given vector
    '''
    u_norm = np.linalg.norm(u)
    if u_norm > 0.0:
        return u / u_norm
    return u

#----------------------------------------------------------------------
def rot2quat(R):
    """
    ROT2QUAT - Transform Rotation matrix into normalized quaternion.
    Usage: q = rot2quat(R)
    Input:
    R - 3-by-3 Rotation matrix
    Output:
    q - 4-by-1 quaternion, with form [w x y z], where w is the scalar term.
    """
    tr = R[0, 0] + R[1, 1] + R[2, 2]

    if tr > 0:
        S = np.sqrt(tr + 1.0) * 2  # S=4*qw
        qw = 0.25 * S
        qx = (R[2, 1] - R[1, 2]) / S
        qy = (R[0, 2] - R[2, 0]) / S
        qz = (R[1, 0] - R[0, 1]) / S
    elif (R[0, 0] > R[1, 1]) and (R[0, 0] > R[2, 2]):
        S = np.sqrt(1.0 + R[0, 0] - R[1, 1] - R[2, 2]) * 2  # S=4*qx
        qw = (R[2, 1] - R[1, 2]) / S
        qx = 0.25 * S
        qy = (R[0, 1] + R[1, 0]) / S
        qz = (R[0, 2] + R[2, 0]) / S
    elif R[1, 1] > R[2, 2]:
        S = np.sqrt(1.0 + R[1, 1] - R[0, 0] - R[2, 2]) * 2  # S=4*qy
        qw = (R[0, 2] - R[2, 0]) / S
        qx = (R[0, 1] + R[1, 0]) / S
        qy = 0.25 * S
        qz = (R[1, 2] + R[2, 1]) / S
    else:
        S = np.sqrt(1.0 + R[2, 2] - R[0, 0] - R[1, 1]) * 2  # S=4*qz
        qw = (R[1, 0] - R[0, 1]) / S
        qx = (R[0, 2] + R[2, 0]) / S
        qy = (R[1, 2] + R[2, 1]) / S
        qz = 0.25 * S

    return normalize(np.array([qw, qx, qy, qz]))
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def quat2rot(q: np.ndarray) -> np.ndarray:
    """Convert quaternion to rotation matrix

    Parameters
    ----------
    q : np.ndarray
        Quaternion

    Returns
    -------
    np.ndarray
        Rotation matrix
    """
    R = np.zeros((3, 3))
    w = q[0]
    x = q[1]
    y = q[2]
    z = q[3]

    x2 = x * x
    y2 = y * y
    z2 = z * z

    wx = w * x
    wy = w * y
    wz = w * z

    xy = x * y
    xz = x * z

    yz = y * z

    R[0][0] = 1.0 - 2.0 * y2 - 2.0 * z2
    R[1][0] = 2.0 * xy + 2.0 * wz
    R[2][0] = 2.0 * xz - 2.0 * wy
    R[0][1] = 2.0 * xy - 2.0 * wz
    R[1][1] = 1.0 - 2.0 * x2 - 2.0 * z2
    R[2][1] = 2.0 * yz + 2.0 * wx
    R[0][2] = 2.0 * xz + 2.0 * wy
    R[1][2] = 2.0 * yz - 2.0 * wx
    R[2][2] = 1.0 - 2.0 * x2 - 2.0 * y2

    return R
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def slerp(q1: np.ndarray, q2: np.ndarray, t: float) -> np.ndarray:
    """Spherical linear interpolation from `q1` to `q2` at `t`

    Parameters
    ----------
    q1 : np.ndarray
        Source quaternion
    q2 : np.ndarray
        Target quaternion
    t : float
        Interpolation factor, between 0 and 1

    Returns
    -------
    np.ndarray
        The spherical linear interpolation between `q1` and `q2` at `t`
    """
    dot = q1.dot(q2)
    q3 = q2
    if dot < 0.0:
        dot = -dot
        q3 = -q2

    if dot < 0.9999:
        angle = np.arccos(dot)
        a = np.sin(angle * (1 - t)) / np.sin(angle)
        b = np.sin(angle * t) / np.sin(angle)
        return a * q1 + b * q3

    # Angle is close to zero - do linear interpolation
    return q1 * (1 - t) + q3 * t
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def loadLaplaceSoln(fileName):
    '''
    Load a solution to a Laplace-Dirichlet problem from a .vtu file and extract
    the solution and its gradients at the cells.

    ARGS:
    fileName : str
        Path to the .vtu file with the Laplace solution. The solution should be
        defined at the nodes. The Laplace fields should be named as follows:
        - Phi_LV_EPI: Laplace field for the endocardium
        - Phi_LV_LV: Laplace field for the left ventricle
        - Phi_LV_RV: Laplace field for the right ventricle
        - Phi_LV_AB: Laplace field for the apex to base direction
    '''

    DATASTR1 = 'Phi_LV_EPI'
    DATASTR2 = 'Phi_LV_LV'
    DATASTR3 = 'Phi_LV_RV'
    DATASTR4 = 'Phi_LV_AB'

    print("   Loading Laplace solution   <---   %s" % (fileName))
    vtuReader = vtk.vtkXMLUnstructuredGridReader()
    vtuReader.SetFileName(fileName)
    vtuReader.Update()

    result_mesh = pv.read(fileName)

    print("   Extracting solution and its gradients at cells")

    pt2Cell = vtk.vtkPointDataToCellData()
    pt2Cell.SetInputConnection(vtuReader.GetOutputPort())
    pt2Cell.PassPointDataOn()

    gradFilter = vtk.vtkGradientFilter()
    gradFilter.SetInputConnection(pt2Cell.GetOutputPort())

    print(f"      Reading {DATASTR1} into cPhiEP") 
    gradFilter.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,\
        DATASTR1)
    gradFilter.SetResultArrayName(DATASTR1 + '_grad')
    gradFilter.Update()
    vtuMesh = gradFilter.GetOutput()
    cPhiEP  = vtuMesh.GetCellData().GetArray(DATASTR1)
    cGPhiEP = vtuMesh.GetCellData().GetArray(DATASTR1+'_grad')

    print(f"      Reading {DATASTR2} into cPhiLV")
    gradFilter.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,\
        DATASTR2)
    gradFilter.SetResultArrayName(DATASTR2 + '_grad')
    gradFilter.Update()
    vtuMesh = gradFilter.GetOutput()
    cPhiLV  = vtuMesh.GetCellData().GetArray(DATASTR2)
    cGPhiLV = vtuMesh.GetCellData().GetArray(DATASTR2+'_grad')

    print(f"      Reading {DATASTR3} into cPhiRV")
    gradFilter.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,\
        DATASTR3)
    gradFilter.SetResultArrayName(DATASTR3 + '_grad')
    gradFilter.Update()
    vtuMesh = gradFilter.GetOutput()
    cPhiRV  = vtuMesh.GetCellData().GetArray(DATASTR3)
    cGPhiRV = vtuMesh.GetCellData().GetArray(DATASTR3+'_grad')

    print(f"      Reading {DATASTR4} into cPhiAB")
    gradFilter.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,\
        DATASTR4)
    gradFilter.SetResultArrayName(DATASTR4 + '_grad')
    gradFilter.Update()
    vtuMesh = gradFilter.GetOutput()
    cPhiAB  = vtuMesh.GetCellData().GetArray(DATASTR4)
    cGPhiAB = vtuMesh.GetCellData().GetArray(DATASTR4+'_grad')

    # Clean unnecessary arrays
    vtuMesh.GetPointData().RemoveArray(DATASTR1)
    vtuMesh.GetCellData().RemoveArray(DATASTR1)
    vtuMesh.GetCellData().RemoveArray(DATASTR1+'_grad')

    vtuMesh.GetPointData().RemoveArray(DATASTR2)
    vtuMesh.GetCellData().RemoveArray(DATASTR2)
    vtuMesh.GetCellData().RemoveArray(DATASTR2+'_grad')

    vtuMesh.GetPointData().RemoveArray(DATASTR3)
    vtuMesh.GetCellData().RemoveArray(DATASTR3)
    vtuMesh.GetCellData().RemoveArray(DATASTR3+'_grad')

    vtuMesh.GetPointData().RemoveArray(DATASTR4)
    vtuMesh.GetCellData().RemoveArray(DATASTR4)
    vtuMesh.GetCellData().RemoveArray(DATASTR4+'_grad')

    cPhiEP = vtknp.vtk_to_numpy(cPhiEP)
    cPhiLV = vtknp.vtk_to_numpy(cPhiLV)
    cPhiRV = vtknp.vtk_to_numpy(cPhiRV)
    cPhiAB = vtknp.vtk_to_numpy(cPhiAB)
    cGPhiEP = vtknp.vtk_to_numpy(cGPhiEP)
    cGPhiLV = vtknp.vtk_to_numpy(cGPhiLV)
    cGPhiRV = vtknp.vtk_to_numpy(cGPhiRV)
    cGPhiAB = vtknp.vtk_to_numpy(cGPhiAB)

    return result_mesh, cPhiEP, cPhiLV, cPhiRV, cPhiAB, \
        cGPhiEP, cGPhiLV, cGPhiRV, cGPhiAB
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def axis (u, v):
    '''
    Given two vectors u and v, compute an orthogonal matrix Q whose first
    column is u, second column is othogonal to u in the direction of v, and
    third column is orthogonal to both u and v.
    '''

    e1 = normalize(u)

    e2 = v - (e1.dot(v)) * e1
    e2 = normalize(e2)

    e0 = np.cross(e1, e2)
    e0 = normalize(e0)

    Q  = np.zeros((3,3))
    Q[:,0] = e0
    Q[:,1] = e1
    Q[:,2] = e2

    return Q
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def orient(Q, alpha, beta):
    '''
    Given an orthogonal matrix Q, rotate it by alpha about the z-axis and
    then by beta about the x-axis.
    '''
    ca = np.cos(alpha)
    sa = np.sin(alpha)
    cb = np.cos(beta)
    sb = np.sin(beta)

    Ra = np.array([ [ ca,  -sa,  0.0],
                    [ sa,   ca,  0.0],
                    [0.0,  0.0,  1.0]])

    # Rb = np.array([ [1.0,  0.0,  0.0],
    #                 [0.0,   cb,   sb],
    #                 [0.0,  -sb,   cb]])
    
    Rb = np.array([ [1.0,  0.0,  0.0],
                    [0.0,   cb,   -sb],
                    [0.0,  sb,   cb]])

    Qt = np.matmul(Q, np.matmul(Ra, Rb))

    return Qt
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def bislerp(Q_A, Q_B, t):
    '''
    :param Q_A: ndarray
    :param Q_B: ndarray
    :param t: float
    :return: ndarray
    Linear interpolation of two orthogonal matrices.
    '''
    qa = rot2quat(Q_A)
    qb = rot2quat(Q_B)

    quat_array = np.array([
        [-qa[1], qa[0], qa[3], -qa[2]],
        [-qa[2], -qa[3], qa[0], qa[1]],
        [-qa[3], qa[2], -qa[1], qa[0]],
    ])

    qm = qa
    max_dot = abs(qm.dot(qb))

    for v in quat_array[0:]:
        dot = abs(v.dot(qb))
        if dot > max_dot:
            max_dot = dot
            qm = v

    qm_slerp = slerp(qm, qb, t)

    return quat2rot(qm_slerp)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def getFiberDirections(vtuMesh, Phi_EP, Phi_LV, Phi_RV, \
                        gPhi_EP, gPhi_LV, gPhi_RV, gPhi_AB, \
                        ALFA_END, ALFA_EPI, BETA_END, BETA_EPI):
    '''
    Compute the fiber directions at the center of each cell
    '''
    EPS = sys.float_info.epsilon

    numCells = vtuMesh.GetNumberOfCells()

    print("   Computing fiber directions at cells")
    F = np.zeros((numCells, 3))
    S = np.zeros((numCells, 3))
    T = np.zeros((numCells, 3))

    j = 1
    k = 1
    print ("      Progress "),
    sys.stdout.flush()
    for iCell in range(0, numCells):
        phiEP = Phi_EP[iCell]
        phiLV = Phi_LV[iCell]
        phiRV = Phi_RV[iCell]

        gPhiEP = gPhi_EP[iCell, :]
        gPhiLV = gPhi_LV[iCell, :]
        gPhiRV = gPhi_RV[iCell, :]
        gPhiAB = gPhi_AB[iCell, :]

        d = phiRV / max(EPS, phiLV + phiRV)
        alfaS = ALFA_END * (1 - d) - ALFA_END * d
        betaS = BETA_END * (1 - d) - BETA_END * d
        alfaW = ALFA_END * (1 - phiEP) + ALFA_EPI * phiEP
        betaW = BETA_END * (1 - phiEP) + BETA_EPI * phiEP

        Q_LV = axis(gPhiAB, - gPhiLV)
        Q_LV = orient(Q_LV, alfaS, betaS)

        Q_RV = axis(gPhiAB, gPhiRV)
        Q_RV = orient(Q_RV, alfaS, betaS)

        Q_END = bislerp(Q_LV, Q_RV, d)

        Q_EPI = axis(gPhiAB, gPhiEP)
        Q_EPI = orient(Q_EPI, alfaW, betaW)

        FST = bislerp(Q_END, Q_EPI, phiEP)

        F[iCell, :] = np.array([FST[0, 0], FST[1, 0], FST[2, 0]])
        S[iCell, :] = np.array([FST[0, 1], FST[1, 1], FST[2, 1]])
        T[iCell, :] = np.array([FST[0, 2], FST[1, 2], FST[2, 2]])
        if iCell==j:
            print ("%d%%  " % ((k-1)*10)),
            sys.stdout.flush()
            k = k + 1
            j = int(float((k-1)*numCells)/10.0)
    print("[Done!]")

    return F, S, T

#----------------------------------------------------------------------

#----------------------------------------------------------------------
def runLaplaceSolver(fd_name, fname, exec_svfsi):
    in_name = "svFSI_LV.xml"
    out_name = os.path.join(fd_name, in_name)
    with open(in_name, 'r') as svFile:
        svRead = svFile.readlines()
    for item in range(len(svRead)):
        # if "Constitutive model" in svRead[item]:
        if "Mesh file path: ./mesh-complete.vtu" in svRead[item]:
            svRead[item] = "   Mesh file path: " + fname + "\n"
        if "Face file path: ./top.vtp" in svRead[item]:
            svRead[item] = "      Face file path: " + fd_name + "mesh-surfaces/top.vtp" + "\n"
        if "Face file path: ./epi.vtp" in svRead[item]:
            svRead[item] = "      Face file path: " + fd_name + "mesh-surfaces/epi.vtp" + "\n"
        if "Face file path: ./epi_apex.vtp" in svRead[item]:
            svRead[item] = "      Face file path: " + fd_name + "mesh-surfaces/epi_apex.vtp" + "\n"
        if "Face file path: ./endo.vtp" in svRead[item]:
            svRead[item] = "      Face file path: " + fd_name + "mesh-surfaces/endo.vtp" + "\n"
        if "Save results in folder:" in svRead[item]:
            svRead[item] = "Save results in folder: " + fd_name + "lps_thickness" + "\n"

    with open(out_name, 'w') as svFileNew:
        svFileNew.writelines(svRead)

    print("   Running svFSI solver")
    print(f"   {exec_svfsi + out_name}")
    os.system(exec_svfsi + out_name)

    return
#----------------------------------------------------------------------

#----------------------------------------------------------------------

def generate_epi_apex_old(fd_name):

    # Generate epi_apex and epi_mid from epi surface
    top_name = fd_name + "mesh-surfaces/top.vtp"
    top_mesh = pv.read(top_name)
    top_mesh = top_mesh.compute_normals()
    #top_normal_rep = top_mesh['Normals'][0, :]
    top_normal_rep = np.mean(top_mesh['Normals'], axis=0)
    top_normal_rep = top_normal_rep / np.linalg.norm(top_normal_rep)
    #top_point_rep = top_mesh.points[0, :]
    top_point_rep = np.mean(top_mesh.points, axis=0)
    

    epi_name = fd_name + "mesh-surfaces/epi.vtp"
    epi_mesh = pv.read(epi_name)
    epi_points = epi_mesh.points
    epi_cells = epi_mesh.faces
    epi_global_node_id = epi_mesh.point_data['GlobalNodeID']
    epi_global_cell_id = epi_mesh.cell_data['GlobalElementID']
    epi_eNoN = epi_cells[0]
    epi_cells = epi_cells.reshape((-1, epi_eNoN + 1))
    epi_cells = epi_cells[:, 1:]

    # Find the index of the apex point of the epi surface
    epi_apex_point_index = 0
    for i in range(epi_points.shape[0]):
        epi_point_dis_to_top = distance_to_surface(epi_points[i, :], top_normal_rep, top_point_rep)
        if i == 0:
            epi_apex_point_dis = epi_point_dis_to_top
        elif epi_point_dis_to_top > epi_apex_point_dis:
            epi_apex_point_index = i
            epi_apex_point_dis = epi_point_dis_to_top

    # Generate epi_apex mesh
    epi_apex_cell_index = np.where(epi_cells == epi_apex_point_index)[0]
    epi_apex_name = fd_name + "mesh-surfaces/epi_apex.vtp"
    epi_apex_cells_surf = epi_cells[epi_apex_cell_index, :]
    epi_apex_cells = np.zeros_like(epi_apex_cells_surf)
    epi_apex_points_index = np.unique(np.hstack(epi_apex_cells_surf))
    epi_apex_points = epi_points[epi_apex_points_index, :]
    for j in range(epi_apex_cells_surf.shape[0]):
        for k in range(epi_apex_cells_surf.shape[1]):
            temp_index = np.where(epi_apex_points_index == epi_apex_cells_surf[j, k])
            temp_index = temp_index[0]
            epi_apex_cells[j, k] = temp_index

    epi_apex_cells_type = np.full((epi_apex_cells.shape[0], 1), epi_eNoN)
    epi_apex_cells = np.hstack((epi_apex_cells_type, epi_apex_cells))
    epi_apex_cells = np.hstack(epi_apex_cells)

    epi_apex_global_node_id = epi_global_node_id[epi_apex_points_index]
    epi_apex_global_cell_id = epi_global_cell_id[epi_apex_cell_index]

    epi_apex_mesh = pv.PolyData(epi_apex_points, epi_apex_cells)

    epi_apex_mesh.point_data.set_array(epi_apex_global_node_id, 'GlobalNodeID')
    epi_apex_mesh.cell_data.set_array(epi_apex_global_cell_id, 'GlobalElementID')

    epi_apex_mesh.save(epi_apex_name)

    # Generate epi_mid mesh (epi_mid = epi - epi_apex)
    epi_mid_cell_index = np.where(~np.any(epi_cells == epi_apex_point_index, axis=1))[0]
    epi_mid_name = fd_name + "mesh-surfaces/epi_mid.vtp"
    epi_mid_cells_surf = epi_cells[epi_mid_cell_index, :]
    epi_mid_cells = np.zeros_like(epi_mid_cells_surf)
    epi_mid_points_index = np.unique(np.hstack(epi_mid_cells_surf))
    epi_mid_points = epi_points[epi_mid_points_index, :]
    for j in range(epi_mid_cells_surf.shape[0]):
        for k in range(epi_mid_cells_surf.shape[1]):
            temp_index = np.where(epi_mid_points_index == epi_mid_cells_surf[j, k])
            temp_index = temp_index[0]
            epi_mid_cells[j, k] = temp_index

    epi_mid_cells_type = np.full((epi_mid_cells.shape[0], 1), epi_eNoN)
    epi_mid_cells = np.hstack((epi_mid_cells_type, epi_mid_cells))
    epi_mid_cells = np.hstack(epi_mid_cells)

    epi_mid_global_node_id = epi_global_node_id[epi_mid_points_index]
    epi_mid_global_cell_id = epi_global_cell_id[epi_mid_cell_index]

    epi_mid_mesh = pv.PolyData(epi_mid_points, epi_mid_cells)

    epi_mid_mesh.point_data.set_array(epi_mid_global_node_id, 'GlobalNodeID')
    epi_mid_mesh.cell_data.set_array(epi_mid_global_cell_id, 'GlobalElementID')

    epi_mid_mesh.save(epi_mid_name)

    return

#----------------------------------------------------------------------

#----------------------------------------------------------------------

def generate_epi_apex(fd_name):
    '''
    Generate the epi apex and epi mid surfaces from the epi surface of the LV.
    '''

    # Generate epi_apex and epi_mid from epi surface
    # top_name = fd_name + "mesh-surfaces/top.vtp"
    # top_mesh = pv.read(top_name)
    # top_mesh = top_mesh.compute_normals()
    # #top_normal_rep = top_mesh['Normals'][0, :]
    # top_normal_rep = np.mean(top_mesh['Normals'], axis=0)
    # top_normal_rep = top_normal_rep / np.linalg.norm(top_normal_rep)
    # #top_point_rep = top_mesh.points[0, :]
    # top_point_rep = np.mean(top_mesh.points, axis=0)
    
    # Load the epi surface
    epi_name = fd_name + "mesh-surfaces/epi.vtp"
    epi_mesh = pv.read(epi_name)
    epi_points = epi_mesh.points
    epi_cells = epi_mesh.faces
    epi_global_node_id = epi_mesh.point_data['GlobalNodeID']
    epi_global_cell_id = epi_mesh.cell_data['GlobalElementID']
    epi_eNoN = epi_cells[0]
    epi_cells = epi_cells.reshape((-1, epi_eNoN + 1))
    epi_cells = epi_cells[:, 1:]

    # Extract the boundary of the epi surface (at the top) to find the apex point
    epi_boundary = epi_mesh.extract_feature_edges(non_manifold_edges=False, feature_edges=False, manifold_edges=False)

    # Triangulate the boundary
    epi_boundary_triangulated = epi_boundary.delaunay_2d()

    # Compute the center and average normal vector of the triangulated boundary
    top_point_rep = np.mean(epi_boundary_triangulated.points, axis=0)
    epi_boundary_triangulated = epi_boundary_triangulated.compute_normals()
    top_normal_rep = np.mean(epi_boundary_triangulated['Normals'], axis=0)
    top_normal_rep = top_normal_rep / np.linalg.norm(top_normal_rep)

    # Show the normal vector and centroid
    #p = pv.Plotter()
    #p.add_mesh(epi_boundary_triangulated, color='white', show_edges=True)
    #p.add_arrows(top_point_rep, top_normal_rep, mag=1, color='r')
    #p.show()

    # Find the index of the apex point of the epi surface
    epi_apex_point_index = 0
    for i in range(epi_points.shape[0]):
        epi_point_dis_to_top = distance_to_surface(epi_points[i, :], top_normal_rep, top_point_rep)
        if i == 0:
            epi_apex_point_dis = epi_point_dis_to_top
        elif epi_point_dis_to_top > epi_apex_point_dis:
            epi_apex_point_index = i
            epi_apex_point_dis = epi_point_dis_to_top

    # Generate epi_apex mesh
    epi_apex_cell_index = np.where(epi_cells == epi_apex_point_index)[0]
    epi_apex_name = fd_name + "mesh-surfaces/epi_apex.vtp"
    epi_apex_cells_surf = epi_cells[epi_apex_cell_index, :]
    epi_apex_cells = np.zeros_like(epi_apex_cells_surf)
    epi_apex_points_index = np.unique(np.hstack(epi_apex_cells_surf))
    epi_apex_points = epi_points[epi_apex_points_index, :]
    for j in range(epi_apex_cells_surf.shape[0]):
        for k in range(epi_apex_cells_surf.shape[1]):
            temp_index = np.where(epi_apex_points_index == epi_apex_cells_surf[j, k])
            temp_index = temp_index[0]
            epi_apex_cells[j, k] = temp_index

    epi_apex_cells_type = np.full((epi_apex_cells.shape[0], 1), epi_eNoN)
    epi_apex_cells = np.hstack((epi_apex_cells_type, epi_apex_cells))
    epi_apex_cells = np.hstack(epi_apex_cells)

    epi_apex_global_node_id = epi_global_node_id[epi_apex_points_index]
    epi_apex_global_cell_id = epi_global_cell_id[epi_apex_cell_index]

    epi_apex_mesh = pv.PolyData(epi_apex_points, epi_apex_cells)

    epi_apex_mesh.point_data.set_array(epi_apex_global_node_id, 'GlobalNodeID')
    epi_apex_mesh.cell_data.set_array(epi_apex_global_cell_id, 'GlobalElementID')

    epi_apex_mesh.save(epi_apex_name)

    # Generate epi_mid mesh (epi_mid = epi - epi_apex)
    epi_mid_cell_index = np.where(~np.any(epi_cells == epi_apex_point_index, axis=1))[0]
    epi_mid_name = fd_name + "mesh-surfaces/epi_mid.vtp"
    epi_mid_cells_surf = epi_cells[epi_mid_cell_index, :]
    epi_mid_cells = np.zeros_like(epi_mid_cells_surf)
    epi_mid_points_index = np.unique(np.hstack(epi_mid_cells_surf))
    epi_mid_points = epi_points[epi_mid_points_index, :]
    for j in range(epi_mid_cells_surf.shape[0]):
        for k in range(epi_mid_cells_surf.shape[1]):
            temp_index = np.where(epi_mid_points_index == epi_mid_cells_surf[j, k])
            temp_index = temp_index[0]
            epi_mid_cells[j, k] = temp_index

    epi_mid_cells_type = np.full((epi_mid_cells.shape[0], 1), epi_eNoN)
    epi_mid_cells = np.hstack((epi_mid_cells_type, epi_mid_cells))
    epi_mid_cells = np.hstack(epi_mid_cells)

    epi_mid_global_node_id = epi_global_node_id[epi_mid_points_index]
    epi_mid_global_cell_id = epi_global_cell_id[epi_mid_cell_index]

    epi_mid_mesh = pv.PolyData(epi_mid_points, epi_mid_cells)

    epi_mid_mesh.point_data.set_array(epi_mid_global_node_id, 'GlobalNodeID')
    epi_mid_mesh.cell_data.set_array(epi_mid_global_cell_id, 'GlobalElementID')

    epi_mid_mesh.save(epi_mid_name)




def generate_fibers_LV_Bayer_cells(mesh_with_Laplace_fields_path, params):
    '''
    Generate fiber directions on a truncated left ventricular geometry using the
    Laplace-Dirichlet rule-based method of Bayer et al. 2012

    ARGS:
    mesh_with_Laplace_fields_path : str
        Path to the .vtu mesh with Laplace fields defined at nodes
    params : dict
        Dictionary of parameters for fiber generation
    '''
    
    t1 = time.time()
    print("========================================================")

    # Get directory of mesh
    fd_name = os.path.dirname(os.path.dirname(mesh_with_Laplace_fields_path))
    
    result_mesh, Phi_EPI, Phi_LV, Phi_RV, Phi_AB, \
    gPhi_EPI, gPhi_LV, gPhi_RV, gPhi_AB = loadLaplaceSoln(mesh_with_Laplace_fields_path)

    # Unpack the input data
    ALFA_END = params['ALFA_END']
    ALFA_EPI = params['ALFA_EPI']
    BETA_END = params['BETA_END']
    BETA_EPI = params['BETA_EPI']

    # Generate fiber directions
    F, S, T = getFiberDirections(result_mesh, Phi_EPI, Phi_LV, Phi_RV,
                                 gPhi_EPI, gPhi_LV, gPhi_RV, gPhi_AB, 
                                 ALFA_END, ALFA_EPI, BETA_END, BETA_EPI)

    print("   Writing domains and fibers to VTK data structure")

    # Write the fiber directions to a vtu files
    output_mesh = copy.deepcopy(result_mesh)
    #output_mesh.cell_data.remove('Proc_ID')
    #output_mesh.point_data.remove(DATASTR1)
    #output_mesh.point_data.remove(DATASTR2)
    #output_mesh.point_data.remove(DATASTR3)
    #output_mesh.point_data.remove(DATASTR4)

    fname1 = os.path.join(fd_name, "fibersLong.vtu")
    print("   Writing to vtu file   --->   %s" % (fname1))
    output_mesh.cell_data.set_array(F, 'FIB_DIR')
    output_mesh.save(fname1)

    fname1 = os.path.join(fd_name, "fibersSheet.vtu")
    print("   Writing to vtu file   --->   %s" % (fname1))
    output_mesh.cell_data.remove('FIB_DIR')
    output_mesh.cell_data.set_array(T, 'FIB_DIR')
    output_mesh.save(fname1)

    fname1 = os.path.join(fd_name, "fibersNormal.vtu")
    print("   Writing to vtu file   --->   %s" % (fname1))
    output_mesh.cell_data.remove('FIB_DIR')
    output_mesh.cell_data.set_array(S, 'FIB_DIR')
    output_mesh.save(fname1)

    t2 = time.time()
    print('\n   Total time: %.3fs' % (t2-t1))
    print("========================================================")
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Main function
if __name__ == '__main__':

    # Whether to solve the Laplace equation or not
    run_flag = True

    # Define fiber and sheet angles
    # Values for human LV from Piersanti_2021 Eq. 8
    # Note: assuming beta is associated with CCW rotation about the fiber axis (new convention)
    params = {}
    # Fiber angle at endocardium
    params['ALFA_END'] = np.radians(60.0)
    # Fiber angle at epicardium
    params['ALFA_EPI'] = np.radians(-60.0)
    # Sheet angle at endocardium
    params['BETA_END'] = np.radians(20.0)
    # Sheet angle at epicardium
    params['BETA_EPI'] = np.radians(-20.0)

    #----------------------------------------------------------------------
    # Navigate to script directory
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

     # Define the path to the mesh and fibers
    fd_name = "./P1_0.005_Exact/"

    # Define the path to the svFSI executable
    exec_svfsi = "../../../../build/svMultiphysics-build/bin/svmultiphysics "

    # Generate the apex and mid epi surfaces
    generate_epi_apex(fd_name)

    # Run the Laplace solver
    if run_flag:
        runLaplaceSolver(fd_name, fd_name + "mesh-complete.mesh.vtu", exec_svfsi)

    # Get the mesh with Laplace fields
    mesh_with_Laplace_fields_path = fd_name + "lps_thickness/result_020.vtu"

    # Generate the fiber directions
    generate_fibers_LV_Bayer_cells(mesh_with_Laplace_fields_path, params)
    

#----------------------------------------------------------------------