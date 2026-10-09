from .conftest import make_interior_node_case, run_by_name, run_with_reference
import copy
import numpy as np
import meshio
import pytest
import xml.etree.ElementTree as ET
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


def make_interior_structural_case(tmp_path):
    root, points, cells, ids = make_interior_node_case(tmp_path, 3)
    root.find("GeneralSimulationParameters/Spectral_radius_of_infinite_time_step").text = "0.5"
    equation = root.find("Add_equation")
    equation.set("type", "lElas")
    for child in list(equation):
        if child.tag.startswith("Darcy_") or child.tag in ("Fluid_density", "Source_term"):
            equation.remove(child)
    for tag, value in [("Density", "1"), ("Elasticity_modulus", "10"),
                       ("Poisson_ratio", "0.3"), ("Force_y", "1")]:
        ET.SubElement(equation, tag).text = value
    equation.find("Max_iterations").text = "20"
    equation.find("LS").set("type", "GMRES")
    equation.find("LS/Linear_algebra/Preconditioner").text = "fsils"
    ET.SubElement(equation.find("LS"), "Absolute_tolerance").text = "1e-13"
    equation.find("Output").clear()
    equation.find("Output").set("type", "Spatial")
    for field in ["Displacement", "Velocity"]:
        ET.SubElement(equation.find("Output"), field).text = "true"
    return root, points, cells, ids


@pytest.mark.parametrize("integral", [False, True])
@pytest.mark.parametrize("prescription", ["components", "history", "masked_history"])
def test_interior_node_components(tmp_path, n_proc, integral, prescription):
    root, points, _, ids = make_interior_structural_case(tmp_path)
    equation = root.find("Add_equation")
    bc = equation.find("Add_BC")
    ET.SubElement(bc, "Impose_on_state_variable_integral").text = str(integral).lower()
    if prescription == "components":
        # A second prescription on the same nodes constrains a disjoint component.
        slope = 0
        bc.remove(bc.find("Temporal_and_spatial_values_file_path"))
        bc.find("Time_dependence").text = "Steady"
        ET.SubElement(bc, "Value").text = "0.02"
        ET.SubElement(bc, "Effective_direction").text = "1 0 0"
        other = copy.deepcopy(bc)
        other.set("name", "z_component")
        other.find("Effective_direction").text = "0 0 1"
        other.find("Value").text = "0.03"
        equation.append(other)
        expected = np.array([[0.02, 0.03]] * 2)
        components = [0, 2]
    else:
        # Distinct linear histories, listed in reverse node-set order.
        slope = np.array([[0.01, -0.02, 0.03], [-0.04, 0.05, -0.06]])
        bc.find("Temporal_and_spatial_values_file_path").text = "vector.dat"
        expected = np.array([[0.01, 0.02, 0.03], [0.04, 0.05, 0.06]])
        components = [0, 2] if prescription == "masked_history" else [0, 1, 2]
        if prescription == "masked_history":
            ET.SubElement(bc, "Effective_direction").text = "1 0 1"
        expected = expected[:, components]
        slope = slope[:, components]
        (tmp_path / "vector.dat").write_text(
            f"{len(components)} 2 2\n0 1\n" + "".join(
                f"{ids[i]} " + " ".join(map(str, np.r_[expected[i], expected[i] + slope[i]])) + "\n"
                for i in [1, 0]))
    ET.ElementTree(root).write(tmp_path / "solver.xml")
    result = run_by_name(tmp_path, "solver.xml", 2, n_proc)
    order = [np.argmin(np.linalg.norm(result.points - points[i-1], axis=1)) for i in ids]
    for step in [0, 1, 2]:
        sample = meshio.read(tmp_path / f"{n_proc}-procs/result_{step:03d}.vtu")
        selected = [np.argmin(np.linalg.norm(sample.points - points[i-1], axis=1)) for i in ids]
        constrained = sample.point_data["Displacement" if integral else "Velocity"][selected]
        np.testing.assert_allclose(constrained[:, components], expected + slope * step * 0.1, atol=1e-10, rtol=0)
        if integral:
            np.testing.assert_allclose(sample.point_data["Velocity"][selected][:, components], slope, atol=1e-10)
    # A mask that accidentally fixes every component suppresses this response.
    if len(components) < 3:
        assert np.all(result.point_data["Displacement"][order, 1] > 1e-5)
    assert np.isfinite(result.point_data["Displacement"]).all()
    if n_proc > 1:
        serial = run_by_name(tmp_path, "solver.xml", 2)
        serial_order = np.lexsort(serial.points.T)
        parallel_order = np.lexsort(result.points.T)
        np.testing.assert_allclose(result.points[parallel_order], serial.points[serial_order], atol=1e-12)
        for field in ("Displacement", "Velocity"):
            np.testing.assert_allclose(result.point_data[field][parallel_order],
                                       serial.point_data[field][serial_order], atol=1e-9, rtol=1e-7)


@pytest.mark.parametrize("last_component", ["1", "1.5", "2"])
@pytest.mark.parametrize("node_set_last", [False, True])
@pytest.mark.parametrize("first_mask,last_mask_format", [("1 0", "{}"), ("(1, 0)", "({})")])
def test_interior_node_repeated_component_mask(
        tmp_path, last_component, node_set_last, first_mask, last_mask_format):
    """Validate every mask tag, regardless of the node-set reference's position."""
    root, points, _, ids = make_interior_structural_case(tmp_path)
    bc = root.find("Add_equation/Add_BC")
    bc.remove(bc.find("Temporal_and_spatial_values_file_path"))
    bc.find("Time_dependence").text = "Steady"
    ET.SubElement(bc, "Value").text = "0.02"
    ET.SubElement(bc, "Impose_on_state_variable_integral").text = "false"
    ET.SubElement(bc, "Effective_direction").text = first_mask
    ET.SubElement(bc, "Effective_direction").text = last_mask_format.format(last_component)
    if node_set_last:
        node_set = bc.find("Node_set")
        bc.remove(node_set)
        bc.append(node_set)
    ET.ElementTree(root).write(tmp_path / "solver.xml")
    error = None if last_component == "1" else "node set.*component mask entries must be 0 or 1"
    result = run_by_name(tmp_path, "solver.xml", 2, expected_error=error)
    if error is None:
        selected = [np.argmin(np.linalg.norm(result.points - points[i-1], axis=1)) for i in ids]
        np.testing.assert_allclose(result.point_data["Velocity"][selected][:, [0, 2]], 0.02, atol=1e-10)
