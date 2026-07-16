from .conftest import run_with_reference, skip_if_no_petsc, skip_if_no_trilinos

# Common folder for all tests in this file
base_folder = "fsi"

# Fields to test
fields = ["Displacement", "Pressure", "Velocity"]


def test_pipe_3d(n_proc):
    test_folder = "pipe_3d"
    t_max = 5
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)

@skip_if_no_petsc
def test_pipe_3d_petsc(n_proc):
    test_folder = "pipe_3d_petsc"
    t_max = 5
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)

@skip_if_no_trilinos
def test_pipe_3d_trilinos_bj(n_proc):
    test_folder = "pipe_3d_trilinos_bj"
    t_max = 5
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)

@skip_if_no_trilinos
def test_pipe_3d_trilinos_ml(n_proc):
    test_folder = "pipe_3d_trilinos_ml"
    t_max = 5
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)

def test_pipe_RCR_3d(n_proc):
    test_folder = "pipe_RCR_3d"
    t_max = 5
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)

def test_pipe_3d_partitioned(n_proc):
    test_folder = "pipe_3d_partitioned"
    t_max = 1
    # A single 1-proc result is the reference for all processor counts, compared
    # at the global RTOL. Each sub-mesh (fluid/solid/mesh) is partitioned
    # independently, so multi-proc runs only match the 1-proc run once the
    # Dirichlet-Neumann coupling has converged tightly. Two settings in
    # solver.xml make that possible at the original tolerances:
    #   - Time_step_size = 1e-2: a larger step weakens the added-mass coupling
    #     stiffness (which grows ~1/dt for partitioned FSI), so the Aitken
    #     iteration converges robustly and identically at 1, 3 and 4 procs
    #     instead of stagnating at 4 procs.
    #   - Coupling_tolerance = 1e-6: at this dt all proc counts reach the same
    #     coupling fixed point in ~10 iterations, leaving the partition-count
    #     round-off 5-10 orders below the global RTOL for every field.
    run_with_reference(base_folder, test_folder, fields=[], n_proc=n_proc, t_max=t_max,
                       comparisons=[
                           {"fields": ["Velocity", "Pressure"],
                            "name_ref": "result_fluid_001.vtu",
                            "name_result": "result_fluid_001.vtu"},
                           {"fields": ["Displacement", "VonMises_stress"],
                            "name_ref": "result_solid_001.vtu",
                            "name_result": "result_solid_001.vtu"},
                       ])