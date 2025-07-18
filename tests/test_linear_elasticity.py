from .conftest import run_with_reference
import os
import subprocess

# Common folder for all tests in this file
base_folder = "linear-elasticity"

# Fields to test
fields = [
    "Displacement",
    "Jacobian",
    "Strain",
    "Stress", 
    "VonMises_stress",
]

def test_beam(n_proc):
    test_folder = "beam"
    t_max = 1
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)


