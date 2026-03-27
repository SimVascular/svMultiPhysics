from .conftest import run_with_reference

# Common folder for all tests in this file
base_folder = "fsi"

# Fields to test
fields = ["Displacement", "Pressure", "Velocity"]


def test_pipe_3d(n_proc):
    test_folder = "pipe_3d"
    t_max = 5
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)

def test_pipe_3d_petsc(n_proc):
    test_folder = "pipe_3d_petsc"
    t_max = 5
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)

def test_pipe_3d_trilinos_bj(n_proc):
    test_folder = "pipe_3d_trilinos_bj"
    t_max = 5
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)

def test_pipe_3d_trilinos_ml(n_proc):
    test_folder = "pipe_3d_trilinos_ml"
    t_max = 5
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)

def test_pipe_RCR_3d(n_proc):
    test_folder = "pipe_RCR_3d"
    t_max = 5
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)