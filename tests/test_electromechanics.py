from .conftest import run_with_reference
import os
import subprocess

# Common folder for all tests in this file
base_folder = "electromechanics"

# Fields to test
fields = [
    "Membrane_potential",
    "Calcium",
    "Cauchy_stress",
    "Def_grad",
    "Displacement",
    "Jacobian",
    "Stress",
    "Strain",
    "Velocity",
    "VonMises_stress",
    "Active_tension_fibers",
    "Active_tension_sheets",
    "Active_tension_normal",
]


def test_slab(n_proc):
    test_folder = "slab"
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max=1)
