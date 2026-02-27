from .conftest import run_with_reference
import os
import subprocess

# Common folder for all tests in this file
base_folder = "uris"

# Fields to test
fields_cfd = ["Velocity", "Pressure", "Traction", "WSS", "Vorticity", "Divergence"]

def test_pipe_uris_cfd(n_proc):
    test_folder = "pipe_uris_cfd"
    t_max = 5
    run_with_reference(base_folder, test_folder, fields_cfd, n_proc, t_max)

fields_fsi = ["Displacement", "Pressure", "Velocity"]

def test_pipe_uris_fsi(n_proc):
    test_folder = "pipe_uris_fsi"
    t_max = 5
    run_with_reference(base_folder, test_folder, fields_fsi, n_proc, t_max)