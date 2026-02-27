#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
# SPDX-License-Identifier: BSD-3-Clause
"""Fiber generation package for biventricular heart models.

This package provides classes and utilities to generate myocardial fiber
orientations for biventricular heart models using Laplace-Dirichlet
rule-based methods.

Modules:
    fiber_generator: Core fiber generation classes (FibGen, FibGenBayer, FibGenDoste)
    laplace_solver: Laplace-Dirichlet equation solver
    surface_names: Surface name definitions for heart geometry
    surface_utils: Utility functions for surface operations
    quat_utils: Quaternion utilities for rotation operations
"""

from fiber_generation.fiber_generator import FibGen, FibGenBayer, FibGenDoste
from fiber_generation.laplace_solver import LaplaceSolver
from fiber_generation.surface_names import SurfaceName

__all__ = [
    'FibGen',
    'FibGenBayer',
    'FibGenDoste',
    'LaplaceSolver',
    'SurfaceName',
]

__version__ = '0.1.0'
