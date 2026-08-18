import pytest

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


@pytest.mark.parametrize("model", ["NashPanfilov", "Regazzoni"])
def test_slab(model, n_proc):
    run_with_reference(base_folder, "slab", fields, n_proc, t_max=1,
                       name_inp=f"solver_{model}.xml",
                       name_ref=f"result_{model}_001.vtu")
