from .conftest import run_with_reference
import os
import subprocess

# Common folder for all tests in this file
base_folder = "ris"

# Fields to test
fields = ["Displacement", "Pressure", "Velocity", "Traction", "WSS"]

def test_pipe_ris_3d(n_proc):
    test_folder = "pipe_ris_3d"
    t_max = 5
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)