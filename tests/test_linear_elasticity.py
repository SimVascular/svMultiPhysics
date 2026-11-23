import shutil
import xml.etree.ElementTree as ET
from pathlib import Path

import meshio
import numpy as np
import pytest

from .conftest import RTOL, run_by_name, run_with_reference

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


@pytest.fixture
def beam_case(tmp_path):
    source = Path(__file__).parent / "cases" / base_folder / "beam"
    return shutil.copytree(source, tmp_path / "beam",
                           ignore=shutil.ignore_patterns("*-procs", "result_*.vtu"))


def test_beam_prescribed_prestress(beam_case, n_proc):
    tree = ET.parse(beam_case / "solver.xml")
    tree.find("GeneralSimulationParameters/Start_saving_after_time_step").text = "0"
    mesh_input = tree.find("Add_mesh")
    mesh = meshio.read(beam_case / mesh_input.findtext("Mesh_file_path").strip())
    prestress = np.zeros((len(mesh.points), 6))
    prestress[:, 0] = 100.0
    mesh.point_data["Stress"] = prestress
    mesh.write(beam_case / "prestress.vtu")
    ET.SubElement(mesh_input, "Prestress_file_path").text = "prestress.vtu"
    tree.write(beam_case / "solver.xml")

    # At zero displacement the prescribed uniaxial stress is the entire stress.
    result = run_by_name(str(beam_case), "solver.xml", 0, n_proc)
    np.testing.assert_allclose(result.point_data["Stress"], prestress)
    np.testing.assert_allclose(result.point_data["VonMises_stress"], 100.0)


def test_beam_prestress_initialization(beam_case, n_proc):
    tree = ET.parse(beam_case / "solver.xml")
    general = tree.find("GeneralSimulationParameters")
    general.find("Number_of_time_steps").text = "2"
    # With alpha_f = 1, the full first-step stress carries over to the second step.
    general.find("Spectral_radius_of_infinite_time_step").text = "0"
    equation = tree.find("Add_equation")
    ET.SubElement(equation, "Prestress").text = "true"
    output = equation.find("Output")
    for field in ("Jacobian", "VonMises_stress"):
        output.remove(output.find(field))
    tree.write(beam_case / "solver.xml")

    second = run_by_name(str(beam_case), "solver.xml", 2, n_proc)
    first = meshio.read(beam_case / f"{n_proc}-procs/result_001.vtu")
    young = float(equation.findtext("Elasticity_modulus"))
    poisson = float(equation.findtext("Poisson_ratio"))
    shear_modulus = young / (2 * (1 + poisson))
    prestress = shear_modulus * first.point_data["Strain"][:, 3:]
    peak = np.max(np.abs(prestress))
    assert peak > 0
    expected = prestress + shear_modulus * second.point_data["Strain"][:, 3:]
    np.testing.assert_allclose(second.point_data["Stress"][:, 3:], expected,
                               rtol=0, atol=RTOL["Stress"] * peak)

