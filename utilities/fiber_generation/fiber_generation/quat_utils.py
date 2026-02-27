#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

def rotm_to_quat_batch(R):
    # R: (N,3,3) -> q: (N,4) [w,x,y,z]
    trace = np.einsum('nii->n', R)
    q = np.zeros((R.shape[0], 4), dtype=float)
    
    # Branch where trace is positive
    mask_t = trace > 0.0
    if np.any(mask_t):
        S = np.sqrt(trace[mask_t] + 1.0) * 2.0
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
        
        for case_idx, (i, j, k) in enumerate([(0, 1, 2), (1, 0, 2), (2, 0, 1)]):
            mask_case = idx == case_idx
            if np.any(mask_case):
                S = np.sqrt(1.0 + Rf[mask_case, i, i] - Rf[mask_case, j, j] - Rf[mask_case, k, k]) * 2.0
                rows = mf_idx[mask_case]
                if case_idx == 0:
                    q[rows, 0] = (Rf[mask_case, 2, 1] - Rf[mask_case, 1, 2]) / S
                    q[rows, 1] = 0.25 * S
                    q[rows, 2] = (Rf[mask_case, 0, 1] + Rf[mask_case, 1, 0]) / S
                    q[rows, 3] = (Rf[mask_case, 0, 2] + Rf[mask_case, 2, 0]) / S
                elif case_idx == 1:
                    q[rows, 0] = (Rf[mask_case, 0, 2] - Rf[mask_case, 2, 0]) / S
                    q[rows, 1] = (Rf[mask_case, 0, 1] + Rf[mask_case, 1, 0]) / S
                    q[rows, 2] = 0.25 * S
                    q[rows, 3] = (Rf[mask_case, 1, 2] + Rf[mask_case, 2, 1]) / S
                else:
                    q[rows, 0] = (Rf[mask_case, 1, 0] - Rf[mask_case, 0, 1]) / S
                    q[rows, 1] = (Rf[mask_case, 0, 2] + Rf[mask_case, 2, 0]) / S
                    q[rows, 2] = (Rf[mask_case, 1, 2] + Rf[mask_case, 2, 1]) / S
                    q[rows, 3] = 0.25 * S
    
    # Normalize for numerical safety
    q /= np.linalg.norm(q, axis=1, keepdims=True)
    return q

def quat_to_rotm_batch(q):
    # q: (N,4) [w,x,y,z] -> R: (N,3,3)
    w, x, y, z = q[:, 0], q[:, 1], q[:, 2], q[:, 3]
    x2, y2, z2 = x * x, y * y, z * z
    wx, wy, wz = w * x, w * y, w * z
    xy, xz, yz = x * y, x * z, y * z
    
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

def quat_multiply_batch(q1, q2):
    """Multiply two quaternions: q1 * q2 (batch operation).
    
    Args:
        q1: First quaternion (N, 4) [w, x, y, z].
        q2: Second quaternion (N, 4) [w, x, y, z].
    
    Returns:
        np.ndarray: Product quaternion (N, 4) [w, x, y, z].
    """
    w1, x1, y1, z1 = q1[:, 0], q1[:, 1], q1[:, 2], q1[:, 3]
    w2, x2, y2, z2 = q2[:, 0], q2[:, 1], q2[:, 2], q2[:, 3]
    
    result = np.zeros_like(q1)
    result[:, 0] = w1*w2 - x1*x2 - y1*y2 - z1*z2  # w
    result[:, 1] = w1*x2 + x1*w2 + y1*z2 - z1*y2  # x
    result[:, 2] = w1*y2 - x1*z2 + y1*w2 + z1*x2  # y
    result[:, 3] = w1*z2 + x1*y2 - y1*x2 + z1*w2  # z
    
    return result

def find_best_quaternions(Q1, Q2):
    """Find the best quaternion representations from rotation matrices (vectorized).
    
    Args:
        Q1: Reference rotation matrix to rotate (N, 3, 3).
        Q2: Target rotation matrix (N, 3, 3).
    
    Returns:
        tuple: (q1_best, q2) where both are (N, 4) quaternions [w, x, y, z].
    """
    n = Q1.shape[0]
    
    # Convert rotation matrices to quaternions
    q1 = rotm_to_quat_batch(Q1)
    q2 = rotm_to_quat_batch(Q2)

    # Create 4 candidate quaternions for each element (original + 3 variations)
    candidates = np.zeros((n, 4, 4), dtype=float)
    
    # Option 0: Original
    candidates[:, 0] = q1
    
    # Option 1: [-q1[1], q1[0], q1[3], -q1[2]]
    candidates[:, 1, 0] = -q1[:, 1]
    candidates[:, 1, 1] = q1[:, 0]
    candidates[:, 1, 2] = q1[:, 3]
    candidates[:, 1, 3] = -q1[:, 2]
    
    # Option 2: [-q1[2], -q1[3], q1[0], q1[1]]
    candidates[:, 2, 0] = -q1[:, 2]
    candidates[:, 2, 1] = -q1[:, 3]
    candidates[:, 2, 2] = q1[:, 0]
    candidates[:, 2, 3] = q1[:, 1]
    
    # Option 3: [-q1[3], q1[2], -q1[1], q1[0]]
    candidates[:, 3, 0] = -q1[:, 3]
    candidates[:, 3, 1] = q1[:, 2]
    candidates[:, 3, 2] = -q1[:, 1]
    candidates[:, 3, 3] = q1[:, 0]
    
    # Compute absolute dot products with q2 for all candidates
    dots = np.abs(np.einsum('ni,nci->nc', q2, candidates))
    
    # Find the candidate with the largest dot product
    best_idx = np.argmax(dots, axis=1)
    
    # Select best quaternion for each element
    q1_best = candidates[np.arange(n), best_idx]
    
    return q1_best, q2
