from .conftest import run_with_reference

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
    run_with_reference(base_folder, "slab", fields, n_proc, t_max=1,
                       name_inp="solver_NashPanfilov.xml",
                       name_ref="result_NashPanfilov_001.vtu")


def test_slab_regazzoni(n_proc):
    run_with_reference(base_folder, "slab", fields, n_proc, t_max=1,
                       name_inp="solver_Regazzoni.xml",
                       name_ref="result_Regazzoni_001.vtu")
