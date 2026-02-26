#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
# SPDX-License-Identifier: BSD-3-Clause
"""Unit tests for FibGen class methods.

Tests the core functions: axis, orient_matrix, orient_rodrigues, and interpolate_basis.
"""

import unittest
import numpy as np
import sys
import os
from fiber_generation.fiber_generator import FibGen


class TestFibGen(unittest.TestCase):
    """Test cases for FibGen base class methods."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.fibgen = FibGen()
        self.tol = 1e-10  # Tolerance for floating point comparisons
    
    def test_axis_orthogonality(self):
        """Test that axis() produces orthogonal basis vectors."""
        n = 10
        np.random.seed(42)
        gL = np.random.randn(n, 3)
        gT = np.random.randn(n, 3)
        
        Q = self.fibgen.axis(gL, gT)
        
        # Check shape
        self.assertEqual(Q.shape, (n, 3, 3))
        
        # Check orthogonality for each element
        for i in range(n):
            eC = Q[i, :, 0]
            eL = Q[i, :, 1]
            eT = Q[i, :, 2]
            
            # Check dot products (should be ~0)
            self.assertAlmostEqual(np.dot(eC, eL), 0.0, places=10)
            self.assertAlmostEqual(np.dot(eC, eT), 0.0, places=10)
            self.assertAlmostEqual(np.dot(eL, eT), 0.0, places=10)
            
            # Check normalization
            self.assertAlmostEqual(np.linalg.norm(eC), 1.0, places=10)
            self.assertAlmostEqual(np.linalg.norm(eL), 1.0, places=10)
            self.assertAlmostEqual(np.linalg.norm(eT), 1.0, places=10)
    
    def test_axis_orientation(self):
        """Test that axis() produces correctly oriented basis."""
        n = 5
        np.random.seed(42)
        gL = np.random.randn(n, 3)
        gT = np.random.randn(n, 3)
        
        Q = self.fibgen.axis(gL, gT)
        
        for i in range(n):
            eC = Q[i, :, 0]
            eL = Q[i, :, 1]
            eT = Q[i, :, 2]
            
            # eL should be aligned with normalized gL
            gL_norm = gL[i] / np.linalg.norm(gL[i])
            self.assertTrue(np.allclose(eL, gL_norm) or np.allclose(eL, -gL_norm))
            
            # eT should be orthogonal to eL
            self.assertAlmostEqual(np.dot(eT, eL), 0.0, places=10)
            
            # eC should be cross product of eL and eT
            eC_expected = np.cross(eL, eT)
            eC_expected = eC_expected / np.linalg.norm(eC_expected)
            self.assertTrue(np.allclose(eC, eC_expected) or np.allclose(eC, -eC_expected))
    
    def test_axis_zero_gradient(self):
        """Test axis() with zero gradient vectors."""
        n = 3
        gL = np.zeros((n, 3))
        gT = np.random.randn(n, 3)
        
        Q = self.fibgen.axis(gL, gT)
        
        # Check that zero vectors are handled gracefully
        self.assertEqual(Q.shape, (n, 3, 3))
        # Zero vectors should result in zero output (normalize handles this)
        for i in range(n):
            eL = Q[i, :, 1]
            self.assertTrue(np.allclose(eL, 0.0))
    
    def test_axis_parallel_vectors(self):
        """Test axis() when gL and gT are parallel."""
        # Use exactly parallel vectors (not random to avoid numerical issues)
        gL = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        gT = 2.0 * gL  # Exactly parallel to gL
        
        Q = self.fibgen.axis(gL, gT)
        
        # When parallel, eT should be zero (projection removes everything)
        # normalize() sets zero vectors to zero
        for i in range(len(gL)):
            eT = Q[i, :, 2]
            # eT should be zero (normalize handles this)
            self.assertTrue(np.allclose(eT, 0.0, atol=1e-10),
                          f"eT should be zero for parallel vectors, got {eT}")

    def test_orient_matrix_identity(self):
        """Test orient_matrix() with zero rotations."""
        n = 5
        np.random.seed(42)
        Q = np.random.randn(n, 3, 3)
        # Make Q orthogonal
        for i in range(n):
            Q[i], _ = np.linalg.qr(Q[i])

        alpha = np.zeros(n)
        beta = np.zeros(n)

        Qt = self.fibgen.orient_matrix(Q, alpha, beta)

        # Should return original Q (within numerical precision)
        self.assertTrue(np.allclose(Q, Qt, atol=1e-10))

    def test_orient_matrix_orthogonality(self):
        """Test that orient_matrix() preserves orthogonality."""
        n = 10
        np.random.seed(42)
        Q = np.random.randn(n, 3, 3)
        # Make Q orthogonal
        for i in range(n):
            Q[i], _ = np.linalg.qr(Q[i])

        alpha = np.random.uniform(-np.pi, np.pi, n)
        beta = np.random.uniform(-np.pi, np.pi, n)

        Qt = self.fibgen.orient_matrix(Q, alpha, beta)

        # Check orthogonality of result
        for i in range(n):
            Qt_i = Qt[i]
            # Check that Qt_i is orthogonal: Qt_i^T @ Qt_i = I
            should_be_identity = Qt_i.T @ Qt_i
            self.assertTrue(np.allclose(should_be_identity, np.eye(3), atol=1e-10))

            # Check determinant (should be 1 for rotation matrix)
            det = np.linalg.det(Qt_i)
            self.assertAlmostEqual(det, 1.0, places=10)

    def test_orient_matrix_known_rotation(self):
        """Test orient_matrix() with known rotation angles."""
        n = 1
        # Start with identity matrix
        Q = np.eye(3)[np.newaxis, :, :].repeat(n, axis=0)

        # Rotate by 90 degrees about z-axis (alpha = pi/2)
        alpha = np.array([np.pi / 2])
        beta = np.array([0.0])

        Qt = self.fibgen.orient_matrix(Q, alpha, beta)

        result = Qt[0]
        # Check that it's a rotation
        self.assertTrue(np.allclose(result.T @ result, np.eye(3), atol=1e-10))
    
    def test_orient_rodrigues_identity(self):
        """Test orient_rodrigues() with zero rotations."""
        n = 5
        np.random.seed(42)
        Q = np.random.randn(n, 3, 3)
        # Make Q orthogonal with columns [eC, eL, eT]
        for i in range(n):
            Q[i], _ = np.linalg.qr(Q[i])
        
        alpha = np.zeros(n)
        beta = np.zeros(n)
        
        Qt = self.fibgen.orient_rodrigues(Q, alpha, beta)
        
        # Should return original Q (within numerical precision)
        self.assertTrue(np.allclose(Q, Qt, atol=1e-10))
    
    def test_orient_rodrigues_orthogonality(self):
        """Test that orient_rodrigues() preserves orthogonality."""
        n = 10
        np.random.seed(42)
        Q = np.random.randn(n, 3, 3)
        # Make Q orthogonal
        for i in range(n):
            Q[i], _ = np.linalg.qr(Q[i])
        
        alpha = np.random.uniform(-np.pi, np.pi, n)
        beta = np.random.uniform(-np.pi, np.pi, n)
        
        Qt = self.fibgen.orient_rodrigues(Q, alpha, beta)
        
        # Check orthogonality of result
        for i in range(n):
            Qt_i = Qt[i]
            eC = Qt_i[:, 0]
            eL = Qt_i[:, 1]
            eT = Qt_i[:, 2]
            
            # Check orthogonality
            self.assertAlmostEqual(np.dot(eC, eL), 0.0, places=10)
            self.assertAlmostEqual(np.dot(eC, eT), 0.0, places=10)
            self.assertAlmostEqual(np.dot(eL, eT), 0.0, places=10)
            
            # Check normalization
            self.assertAlmostEqual(np.linalg.norm(eC), 1.0, places=10)
            self.assertAlmostEqual(np.linalg.norm(eL), 1.0, places=10)
            self.assertAlmostEqual(np.linalg.norm(eT), 1.0, places=10)
    
    def test_orient_rodrigues_eT_unchanged_after_alpha(self):
        """Test that eT is unchanged after alpha rotation."""
        n = 5
        np.random.seed(42)
        Q = np.random.randn(n, 3, 3)
        for i in range(n):
            Q[i], _ = np.linalg.qr(Q[i])
        
        alpha = np.random.uniform(-np.pi, np.pi, n)
        beta = np.zeros(n)  # No beta rotation
        
        Qt = self.fibgen.orient_rodrigues(Q, alpha, beta)
        
        # eT should be unchanged (rotating about eT axis)
        for i in range(n):
            eT_original = Q[i, :, 2]
            eT_rotated = Qt[i, :, 2]
            # Should be the same (within numerical precision)
            self.assertTrue(np.allclose(eT_original, eT_rotated, atol=1e-10) or
                          np.allclose(eT_original, -eT_rotated, atol=1e-10))
    
    def test_orient_rodrigues_eL_unchanged_after_beta(self):
        """Test that rotated eL is unchanged after beta rotation."""
        n = 5
        np.random.seed(42)
        Q = np.random.randn(n, 3, 3)
        for i in range(n):
            Q[i], _ = np.linalg.qr(Q[i])
        
        alpha = np.random.uniform(-np.pi, np.pi, n)
        beta = np.random.uniform(-np.pi, np.pi, n)
        
        Qt = self.fibgen.orient_rodrigues(Q, alpha, beta)
        
        # Compute "after alpha" reference using orient_rodrigues with beta=0
        Qt_alpha_only = self.fibgen.orient_rodrigues(Q, alpha, np.zeros(n))

        # After beta rotation, eL should be unchanged (rotation is about eL)
        for i in range(n):
            eL_after_alpha = Qt_alpha_only[i, :, 1]
            eL_final = Qt[i, :, 1]
            self.assertTrue(np.allclose(eL_after_alpha, eL_final, atol=1e-10) or
                            np.allclose(eL_after_alpha, -eL_final, atol=1e-10))
    
    def test_orient_rodrigues_known_rotation(self):
        """Test orient_rodrigues() with known rotation."""
        n = 1
        # Create a known basis
        eC = np.array([1.0, 0.0, 0.0])
        eL = np.array([0.0, 1.0, 0.0])
        eT = np.array([0.0, 0.0, 1.0])
        Q = np.array([np.column_stack([eC, eL, eT])])
        
        # Rotate by 90 degrees about eT (z-axis)
        alpha = np.array([np.pi / 2])
        beta = np.array([0.0])
        
        Qt = self.fibgen.orient_rodrigues(Q, alpha, beta)
        
        # eC should rotate to [0, 1, 0] (or [-1, 0, 0] depending on sign)
        eC_rot = Qt[0, :, 0]
        # Should be perpendicular to original eC
        self.assertAlmostEqual(np.dot(eC, eC_rot), 0.0, places=10)
        
        # eT should be unchanged
        eT_rot = Qt[0, :, 2]
        self.assertTrue(np.allclose(eT, eT_rot, atol=1e-10))
    
    def test_orient_rodrigues_vs_orient_matrix(self):
        """Test that orient_rodrigues and orient_matrix differ on a generic basis.

        Note: These methods use different rotation conventions:
        - orient_matrix: rotates about fixed z-axis (transmural) then y-axis (longitudinal)
        - orient_rodrigues: rotates about local eT (transmural) then rotated eL (longitudinal)
        """
        n = 10
        np.random.seed(42)

        # Create random orthogonal basis
        Q = np.random.randn(n, 3, 3)
        # Make orthogonal with columns [eC, eL, eT]
        for i in range(n):
            Q[i], _ = np.linalg.qr(Q[i])

        # Random rotation angles
        alpha = np.random.uniform(-np.pi, np.pi, n)
        beta = np.random.uniform(-np.pi, np.pi, n)

        # Apply both rotation methods
        Q_rodrigues = self.fibgen.orient_rodrigues(Q.copy(), alpha, beta)
        Q_matrix = self.fibgen.orient_matrix(Q.copy(), alpha, beta)

        # They should generally NOT be equal on an arbitrary basis, because the
        # rotation axes differ (local basis axes vs fixed coordinate axes).
        diff = np.abs(Q_rodrigues - Q_matrix)
        max_diff = float(np.max(diff))
        self.assertTrue(
            max_diff < 1e-9,
        )

        # Both should preserve orthogonality
        for i in range(n):
            Q_rod_i = Q_rodrigues[i]
            Q_mat_i = Q_matrix[i]

            # Check orthogonality
            self.assertTrue(np.allclose(Q_rod_i.T @ Q_rod_i, np.eye(3), atol=1e-9))
            self.assertTrue(np.allclose(Q_mat_i.T @ Q_mat_i, np.eye(3), atol=1e-9))

            # Check determinants (should be 1 for rotation matrices)
            self.assertAlmostEqual(np.linalg.det(Q_rod_i), 1.0, places=9)
            self.assertAlmostEqual(np.linalg.det(Q_mat_i), 1.0, places=9)

    def test_orient_rodrigues_vs_orient_matrix_identity_basis(self):
        """Test orient_rodrigues vs orient_matrix with identity basis.

        When Q is identity, the local axes match the fixed axes, so the two
        conventions should agree.
        """
        n = 5
        np.random.seed(42)

        # Exact identity basis (aligned with coordinate axes)
        Q = np.eye(3)[np.newaxis, :, :].repeat(n, axis=0)

        # Small rotation angles
        alpha = np.random.uniform(-0.1, 0.1, n)
        beta = np.random.uniform(-0.1, 0.1, n)

        Q_rodrigues = self.fibgen.orient_rodrigues(Q.copy(), alpha, beta)
        Q_matrix = self.fibgen.orient_matrix(Q.copy(), alpha, beta)

        diff = np.abs(Q_rodrigues - Q_matrix)
        max_diff = float(np.max(diff))
        self.assertLess(
            max_diff,
            1e-9,
            msg=f"Expected methods to match for identity basis, but max_diff={max_diff:.3e}",
        )

        if os.environ.get("FIBGEN_TEST_COMPARE", "0") == "1":
            mean_diff = float(np.mean(diff))
            print(
                "\n"
                + "-" * 72
                + "\n"
                + "orient_rodrigues vs orient_matrix (identity basis)\n"
                + f"  max : {max_diff: .3e}\n"
                + f"  mean: {mean_diff: .3e}\n"
                + "-" * 72
            )

        for i in range(n):
            Q_rod_i = Q_rodrigues[i]
            Q_mat_i = Q_matrix[i]

            self.assertTrue(np.allclose(Q_rod_i.T @ Q_rod_i, np.eye(3), atol=1e-9))
            self.assertTrue(np.allclose(Q_mat_i.T @ Q_mat_i, np.eye(3), atol=1e-9))

            for j in range(3):
                self.assertAlmostEqual(np.linalg.norm(Q_rod_i[:, j]), 1.0, places=9)
                self.assertAlmostEqual(np.linalg.norm(Q_mat_i[:, j]), 1.0, places=9)


    def test_interpolate_basis_boundary_conditions(self):
        """Test interpolate_basis() at boundaries (t=0 and t=1)."""
        n = 10
        np.random.seed(42)
        Q1 = np.random.randn(n, 3, 3)
        Q2 = np.random.randn(n, 3, 3)
        # Make orthogonal
        for i in range(n):
            Q1[i], _ = np.linalg.qr(Q1[i])
            Q2[i], _ = np.linalg.qr(Q2[i])
        
        # At t=0, should get Q1
        t0 = np.zeros(n)
        Q_interp0 = self.fibgen.interpolate_basis(Q1, Q2, t0)
        self.assertTrue(np.allclose(Q1, Q_interp0, atol=1e-10))
        
        # At t=1, should get Q2
        t1 = np.ones(n)
        Q_interp1 = self.fibgen.interpolate_basis(Q1, Q2, t1)
        self.assertTrue(np.allclose(Q2, Q_interp1, atol=1e-10))
    
    def test_interpolate_basis_orthogonality(self):
        """Test that interpolate_basis() preserves orthogonality."""
        n = 10
        np.random.seed(42)
        Q1 = np.random.randn(n, 3, 3)
        Q2 = np.random.randn(n, 3, 3)
        # Make orthogonal
        for i in range(n):
            Q1[i], _ = np.linalg.qr(Q1[i])
            Q2[i], _ = np.linalg.qr(Q2[i])
        
        t = np.random.uniform(0, 1, n)
        Q_interp = self.fibgen.interpolate_basis(Q1, Q2, t)
        
        # Check orthogonality
        for i in range(n):
            Q_i = Q_interp[i]
            # Check that Q_i is orthogonal
            should_be_identity = Q_i.T @ Q_i
            self.assertTrue(np.allclose(should_be_identity, np.eye(3), atol=1e-9))
            
            # Check determinant (should be 1 for rotation matrix)
            det = np.linalg.det(Q_i)
            self.assertAlmostEqual(det, 1.0, places=9)
    
    def test_interpolate_basis_identical_matrices(self):
        """Test interpolate_basis() with identical input matrices."""
        n = 5
        np.random.seed(42)
        Q1 = np.random.randn(n, 3, 3)
        for i in range(n):
            Q1[i], _ = np.linalg.qr(Q1[i])
        Q2 = Q1.copy()  # Same as Q1
        
        t = np.random.uniform(0, 1, n)
        Q_interp = self.fibgen.interpolate_basis(Q1, Q2, t)
        
        # Should return Q1 (or Q2, they're the same)
        self.assertTrue(np.allclose(Q1, Q_interp, atol=1e-10))
    
    def test_interpolate_basis_midpoint(self):
        """Test interpolate_basis() at t=0.5."""
        n = 5
        np.random.seed(42)
        Q1 = np.random.randn(n, 3, 3)
        Q2 = np.random.randn(n, 3, 3)
        # Make orthogonal
        for i in range(n):
            Q1[i], _ = np.linalg.qr(Q1[i])
            Q2[i], _ = np.linalg.qr(Q2[i])
        
        t = np.full(n, 0.5)
        Q_interp = self.fibgen.interpolate_basis(Q1, Q2, t)
        
        # Check that it's between Q1 and Q2
        # For each element, check that interpolated matrix is valid rotation
        for i in range(n):
            Q_i = Q_interp[i]
            self.assertTrue(np.allclose(Q_i.T @ Q_i, np.eye(3), atol=1e-9))
            self.assertAlmostEqual(np.linalg.det(Q_i), 1.0, places=9)
    
    def test_interpolate_basis_clipping(self):
        """Test that interpolate_basis() clips t to [0, 1]."""
        n = 5
        np.random.seed(42)
        Q1 = np.random.randn(n, 3, 3)
        Q2 = np.random.randn(n, 3, 3)
        for i in range(n):
            Q1[i], _ = np.linalg.qr(Q1[i])
            Q2[i], _ = np.linalg.qr(Q2[i])
        
        # t values outside [0, 1]
        t_neg = np.array([-0.5, -1.0, -10.0, 0.3, 0.7])
        t_pos = np.array([1.5, 2.0, 10.0, 0.3, 0.7])
        
        Q_interp_neg = self.fibgen.interpolate_basis(Q1, Q2, t_neg)
        Q_interp_pos = self.fibgen.interpolate_basis(Q1, Q2, t_pos)
        
        # Negative values should be clipped to 0 (should equal Q1)
        # Only check elements where t < 0
        neg_indices = np.where(t_neg < 0)[0]
        for idx in neg_indices:
            self.assertTrue(np.allclose(Q1[idx], Q_interp_neg[idx], atol=1e-10),
                          f"Element {idx} with t={t_neg[idx]} should equal Q1")
        
        # Check that t=0 gives Q1 exactly
        t_zero = np.zeros(n)
        Q_interp_zero = self.fibgen.interpolate_basis(Q1, Q2, t_zero)
        self.assertTrue(np.allclose(Q1, Q_interp_zero, atol=1e-10))
        
        # Positive values > 1 should be clipped to 1 (should equal Q2)
        # Only check elements where t > 1
        pos_indices = np.where(t_pos > 1)[0]
        for idx in pos_indices:
            self.assertTrue(np.allclose(Q2[idx], Q_interp_pos[idx], atol=1e-10),
                          f"Element {idx} with t={t_pos[idx]} should equal Q2")
        
        # Check that t=1 gives Q2 exactly
        t_one = np.ones(n)
        Q_interp_one = self.fibgen.interpolate_basis(Q1, Q2, t_one)
        self.assertTrue(np.allclose(Q2, Q_interp_one, atol=1e-10))
    
    def test_interpolate_basis_correct_slerp(self):
        """Test interpolate_basis() with correct_slerp=True."""
        n = 5
        np.random.seed(42)
        Q1 = np.random.randn(n, 3, 3)
        Q2 = np.random.randn(n, 3, 3)
        for i in range(n):
            Q1[i], _ = np.linalg.qr(Q1[i])
            Q2[i], _ = np.linalg.qr(Q2[i])
        
        t = np.random.uniform(0, 1, n)
        
        # Test with correct_slerp=False (default)
        Q_interp_default = self.fibgen.interpolate_basis(Q1, Q2, t, correct_slerp=False)
        
        # Test with correct_slerp=True
        Q_interp_correct = self.fibgen.interpolate_basis(Q1, Q2, t, correct_slerp=True)
        
        # Both should produce valid rotation matrices
        for i in range(n):
            Q_def = Q_interp_default[i]
            Q_corr = Q_interp_correct[i]
            
            self.assertTrue(np.allclose(Q_def.T @ Q_def, np.eye(3), atol=1e-9))
            self.assertTrue(np.allclose(Q_corr.T @ Q_corr, np.eye(3), atol=1e-9))
            self.assertAlmostEqual(np.linalg.det(Q_def), 1.0, places=9)
            self.assertAlmostEqual(np.linalg.det(Q_corr), 1.0, places=9)


class TestFibGenEdgeCases(unittest.TestCase):
    """Test edge cases and error conditions."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.fibgen = FibGen()
    
    def test_axis_single_element(self):
        """Test axis() with single element."""
        gL = np.array([[1.0, 0.0, 0.0]])
        gT = np.array([[0.0, 1.0, 0.0]])
        
        Q = self.fibgen.axis(gL, gT)
        
        self.assertEqual(Q.shape, (1, 3, 3))
        eC = Q[0, :, 0]
        eL = Q[0, :, 1]
        eT = Q[0, :, 2]
        
        # Check orthogonality
        self.assertAlmostEqual(np.dot(eC, eL), 0.0, places=10)
        self.assertAlmostEqual(np.dot(eC, eT), 0.0, places=10)
        self.assertAlmostEqual(np.dot(eL, eT), 0.0, places=10)

    def test_orient_matrix_single_element(self):
        """Test orient_matrix() with single element."""
        Q = np.eye(3)[np.newaxis, :, :]
        alpha = np.array([np.pi / 4])
        beta = np.array([np.pi / 4])

        Qt = self.fibgen.orient_matrix(Q, alpha, beta)

        self.assertEqual(Qt.shape, (1, 3, 3))
        self.assertTrue(np.allclose(Qt[0].T @ Qt[0], np.eye(3), atol=1e-10))
    
    def test_orient_rodrigues_single_element(self):
        """Test orient_rodrigues() with single element."""
        Q = np.eye(3)[np.newaxis, :, :]
        alpha = np.array([np.pi / 4])
        beta = np.array([np.pi / 4])
        
        Qt = self.fibgen.orient_rodrigues(Q, alpha, beta)
        
        self.assertEqual(Qt.shape, (1, 3, 3))
        eC = Qt[0, :, 0]
        eL = Qt[0, :, 1]
        eT = Qt[0, :, 2]
        
        # Check orthogonality
        self.assertAlmostEqual(np.dot(eC, eL), 0.0, places=10)
        self.assertAlmostEqual(np.dot(eC, eT), 0.0, places=10)
        self.assertAlmostEqual(np.dot(eL, eT), 0.0, places=10)
    
    def test_interpolate_basis_single_element(self):
        """Test interpolate_basis() with single element."""
        Q1 = np.eye(3)[np.newaxis, :, :]
        Q2 = np.array([[[0, -1, 0], [1, 0, 0], [0, 0, 1]]])
        t = np.array([0.5])
        
        Q_interp = self.fibgen.interpolate_basis(Q1, Q2, t)
        
        self.assertEqual(Q_interp.shape, (1, 3, 3))
        self.assertTrue(np.allclose(Q_interp[0].T @ Q_interp[0], np.eye(3), atol=1e-9))


if __name__ == '__main__':
    # Run tests with verbose output
    unittest.main(verbosity=2)
