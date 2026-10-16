import copy
import math

import meshio
import numpy as np
import pytest
import xml.etree.ElementTree as ET

from .conftest import add_test_boundary_face, make_interior_node_case, run_by_name, run_with_reference


def test_interior_node_pressure(n_proc):
    run_with_reference("darcy", "interior_node_pressure", ["Darcy_pressure"], n_proc, 2)


@pytest.mark.parametrize("nsd", [2, 3])
def test_interior_node_scalar(tmp_path, n_proc, nsd):
    """Catch ignored constraints, incorrect ID/value mapping, and free-row errors."""
    root, points, cells, ids = make_interior_node_case(tmp_path, nsd)
    ET.ElementTree(root).write(tmp_path / "solver.xml")
    result = run_by_name(tmp_path, "solver.xml", 2, n_proc)

    # Assemble an independent linear simplex diffusion problem with natural
    # exterior boundaries and eliminate the prescribed degrees of freedom.
    stiffness = np.zeros((len(points), len(points)))
    load = np.zeros(len(points))
    for cell in cells:
        coordinates = points[cell, :nsd]
        shape = np.column_stack((np.ones(nsd + 1), coordinates))
        gradients = np.linalg.inv(shape)[1:, :]
        volume = abs(np.linalg.det(shape)) / math.factorial(nsd)
        stiffness[np.ix_(cell, cell)] += volume * gradients.T @ gradients
        load[cell] += 0.25 * volume / (nsd + 1)
    fixed = ids - 1
    free = np.setdiff1d(np.arange(len(points)), fixed)
    expected = np.zeros(len(points))
    expected[fixed] = [2.0, 7.0]
    expected[free] = np.linalg.solve(
        stiffness[np.ix_(free, free)],
        load[free] - stiffness[np.ix_(free, fixed)] @ expected[fixed],
    )
    # Match physical coordinates, independently of the solver's output order.
    order = [np.argmin(np.linalg.norm(result.points - point, axis=1)) for point in points]
    np.testing.assert_allclose(result.points[order], points, atol=1e-12)
    actual = result.point_data["Darcy_pressure"].reshape(-1)[order]
    np.testing.assert_allclose(actual[fixed], [2.0, 7.0], atol=1e-10, rtol=0)
    np.testing.assert_allclose(actual, expected, atol=1e-8, rtol=1e-8)
    np.testing.assert_allclose((stiffness @ actual - load)[free], 0, atol=1e-8)


@pytest.mark.parametrize("ids", ["", "0", "-1", "17", "2.5", "3x", "6 6", "999999999999999999999"])
def test_interior_node_invalid_ids(tmp_path, ids):
    root, _, _, _ = make_interior_node_case(tmp_path)
    root.find("Add_mesh/Add_node_set/Node_IDs").text = ids
    ET.ElementTree(root).write(tmp_path / "solver.xml")
    run_by_name(tmp_path, "solver.xml", 2, expected_error="node.*(ID|empty|value)|Node_IDs")


@pytest.mark.parametrize("tag,value,error", [
    ("Mesh_name", "missing", "Unknown mesh"),
    ("Node_set", "missing", "Unknown node set"),
    ("Type", "Neumann", "node.*Dirichlet"),
    ("Weakly_applied", "true", "node.*Weakly_applied"),
    ("Profile", "Parabolic", "node.*Profile"),
    ("Bct_file_path", "unused.vtp", "node.*Bct_file_path"),
    ("Time_dependence", "Coupled", "node.*Time_dependence"),
    ("Effective_direction", "0 1", "node.*component"),
    ("Effective_direction", "1 0.5", "node.*component"),
    ("Initial_displacements_file_path", "unused.vtu", "node.*Initial_displacements_file_path"),
    ("Prestress_file_path", "unused.vtu", "node.*Prestress_file_path"),
    ("Penalty_parameter_normal", "1", "node.*Penalty_parameter_normal"),
    ("Penalty_parameter_tangential", "1", "node.*Penalty_parameter_tangential"),
    ("Impose_on_state_variable_integral", "true", "node.*Impose_on_state_variable_integral"),
])
def test_interior_node_invalid_options(tmp_path, tag, value, error):
    root, _, _, _ = make_interior_node_case(tmp_path)
    bc = root.find("Add_equation/Add_BC")
    element = bc.find(tag)
    if element is None:
        element = ET.SubElement(bc, tag)
    element.text = value
    ET.ElementTree(root).write(tmp_path / "solver.xml")
    run_by_name(tmp_path, "solver.xml", 2, expected_error=error)


@pytest.mark.parametrize("value", ["2.5", "2 3", "2junk", "nan", "inf"])
@pytest.mark.parametrize("repeated", [False, True])
def test_interior_node_value_literals(tmp_path, value, repeated):
    root, points, _, ids = make_interior_node_case(tmp_path)
    bc = root.find("Add_equation/Add_BC")
    bc.remove(bc.find("Temporal_and_spatial_values_file_path"))
    bc.find("Time_dependence").text = "Steady"
    if repeated:
        ET.SubElement(bc, "Value").text = "1"
    ET.SubElement(bc, "Value").text = value
    ET.ElementTree(root).write(tmp_path / "solver.xml")
    error = None if value == "2.5" else "node set.*Value must be one finite number"
    result = run_by_name(tmp_path, "solver.xml", 2, expected_error=error)
    if error is None:
        order = [np.argmin(np.linalg.norm(result.points - points[i-1], axis=1)) for i in ids]
        np.testing.assert_allclose(result.point_data["Darcy_pressure"].reshape(-1)[order], 2.5, atol=1e-10)


@pytest.mark.parametrize("data", [
    "1 2 1\n0 1\n6 7 7\n",  # Missing selected node.
    "1 2 2\n0 1\n6 7 7\n6 2 2\n",  # Duplicate record.
    "1 2 2\n0 1\n6 7 7\n1 2 2\n",  # Outside the node set.
    "1 2 2\n0 1\n6 7 7\n11 2\n",  # Truncated data.
    "1 2 2\n0 1\n6 7 7\n11 2 2\nextra\n",
    "1 2 2\n0 1\n6 nan 7\n11 2 2\n",
    "1 2 2\n1 2\n6 7 7\n11 2 2\n",
    "1 2 2\n0 0\n6 7 7\n11 2 2\n",
    "1 1 2\n0\n6 7\n11 2\n",
    "2 2 2\n0 1\n6 7 8 7 8\n11 2 3 2 3\n",
    "1 2 2\n0 1\n6.5 7\n11 2 2\n",  # A fractional ID must not become a value.
    "1 2 2.0 1\n6 7 7\n11 2 2\n",  # A fractional count must not become time zero.
])
def test_interior_node_invalid_values(tmp_path, data):
    root, _, _, _ = make_interior_node_case(tmp_path)
    (tmp_path / "values.dat").write_text(data)
    ET.ElementTree(root).write(tmp_path / "solver.xml")
    run_by_name(tmp_path, "solver.xml", 2, expected_error="temporal and spatial values|node.*component")


@pytest.mark.parametrize("mode", ["steady", "temporal", "fourier", "general", "file"])
def test_interior_node_value_modes(tmp_path, mode):
    root, points, _, ids = make_interior_node_case(tmp_path)
    bc = root.find("Add_equation/Add_BC")
    if mode in ("steady", "temporal", "fourier"):
        bc.remove(bc.find("Temporal_and_spatial_values_file_path"))
        bc.find("Time_dependence").text = "Steady" if mode == "steady" else "Unsteady"
        if mode == "steady":
            ET.SubElement(bc, "Value").text = "3"
        elif mode == "temporal":
            ET.SubElement(bc, "Temporal_values_file_path").text = "history.dat"
            ET.SubElement(bc, "Ramp_function").text = "true"
            (tmp_path / "history.dat").write_text("2 1\n0 3\n1 5\n")
        else:
            ET.SubElement(bc, "Fourier_coefficients_file_path").text = "history.fcs"
            (tmp_path / "history.fcs").write_text("0 1\n3 0\n2\n0 0\n1 0\n")
    elif mode == "general":
        (tmp_path / "values.dat").write_text(f"1 2 2\n0 1\n{ids[1]} 7 5\n{ids[0]} 2 3\n")
        root.find("GeneralSimulationParameters/Number_of_time_steps").text = "11"
    else:
        node_set = root.find("Add_mesh/Add_node_set")
        node_set.remove(node_set.find("Node_IDs"))
        ET.SubElement(node_set, "Node_IDs_file_path").text = "nodes.dat"
        (tmp_path / "nodes.dat").write_text(f"{ids[0]}\n{ids[1]}\n")
    ET.ElementTree(root).write(tmp_path / "solver.xml")
    last = 11 if mode == "general" else 2
    run_by_name(tmp_path, "solver.xml", last)
    for step in [0, 1, last]:
        result = meshio.read(tmp_path / f"1-procs/result_{step:03d}.vtu")
        order = [np.argmin(np.linalg.norm(result.points - points[i-1], axis=1)) for i in ids]
        time = step * 0.1
        expected = {"steady": [3, 3], "temporal": [3 + 2*time] * 2,
                    "fourier": [3 + np.cos(2*np.pi*time)] * 2,
                    "general": [2 + time % 1, 7 - 2*(time % 1)], "file": [2, 7]}[mode]
        np.testing.assert_allclose(result.point_data["Darcy_pressure"].reshape(-1)[order], expected, atol=1e-9, rtol=0)


@pytest.mark.parametrize("other", ["node_set", "face"])
def test_interior_node_overlap(tmp_path, n_proc, other):
    """Reject any other Dirichlet prescription of a node-set component, even an identical one."""
    root, points, cells, ids = make_interior_node_case(tmp_path, 3)
    equation = root.find("Add_equation")
    node_bc = equation.find("Add_BC")
    if other == "node_set":
        second = copy.deepcopy(node_bc)
        second.set("name", "second_pressure")
    else:
        ids[0] = add_test_boundary_face(root, tmp_path, points, cells)[0]
        root.find("Add_mesh/Add_node_set/Node_IDs").text = " ".join(map(str, ids))
        (tmp_path / "values.dat").write_text(f"1 2 2\n0 1\n{ids[0]} 2 2\n{ids[1]} 7 7\n")
        second = ET.Element("Add_BC", name="surface")
        for tag, value in [("Type", "Dirichlet"), ("Value", "2"), ("Zero_out_perimeter", "false")]:
            ET.SubElement(second, tag).text = value
    # Read the other condition first, so the check cannot rely on input order.
    equation.insert(list(equation).index(node_bc), second)
    ET.ElementTree(root).write(tmp_path / "solver.xml")
    run_by_name(tmp_path, "solver.xml", 2, n_proc, expected_error="node set 'interior'.*also prescribe")


def test_interior_node_restart(tmp_path):
    root, points, cells, ids = make_interior_node_case(tmp_path)
    meshio.write_points_cells(tmp_path / "initial.vtu", points, [("triangle", cells)],
                              point_data={"Pressure": np.full(len(points), 17.0)})
    ET.SubElement(root.find("Add_mesh"), "Initial_pressures_file_path").text = "initial.vtu"
    ET.ElementTree(root).write(tmp_path / "solver.xml")
    first = run_by_name(tmp_path, "solver.xml", 2)
    initial = meshio.read(tmp_path / "1-procs/result_000.vtu")
    order = [np.argmin(np.linalg.norm(initial.points - points[i-1], axis=1)) for i in ids]
    np.testing.assert_allclose(initial.point_data["Darcy_pressure"].reshape(-1)[order], [2, 7], atol=1e-10)
    free = np.setdiff1d(np.arange(len(points)), order)
    np.testing.assert_allclose(initial.point_data["Darcy_pressure"].reshape(-1)[free], 17, atol=1e-10)
    root.find("GeneralSimulationParameters/Continue_previous_simulation").text = "true"
    root.find("GeneralSimulationParameters/Number_of_time_steps").text = "4"
    ET.ElementTree(root).write(tmp_path / "solver.xml")
    restarted = run_by_name(tmp_path, "solver.xml", 4, clean=False)
    np.testing.assert_allclose(restarted.point_data["Darcy_pressure"], first.point_data["Darcy_pressure"], atol=1e-8)
    np.testing.assert_allclose(restarted.field_data["TimeValue"], 0.4, atol=1e-12)


def test_interior_node_multiple_meshes(tmp_path):
    root, points, cells, _ = make_interior_node_case(tmp_path)
    other_points = points + [2, 0, 0]
    meshio.write_points_cells(tmp_path / "other.vtu", other_points, [("triangle", cells)])
    other_mesh = copy.deepcopy(root.find("Add_mesh"))
    other_mesh.set("name", "other")
    other_mesh.find("Mesh_file_path").text = "other.vtu"
    root.insert(2, other_mesh)
    equation = root.find("Add_equation")
    equation.find("Source_term").text = "0"
    bc = equation.find("Add_BC")
    bc.remove(bc.find("Temporal_and_spatial_values_file_path"))
    bc.find("Time_dependence").text = "Steady"
    ET.SubElement(bc, "Value").text = "2"
    other_bc = copy.deepcopy(bc)
    other_bc.set("name", "other_pressure")
    other_bc.find("Mesh_name").text = "other"
    other_bc.find("Value").text = "11"
    equation.append(other_bc)
    ET.ElementTree(root).write(tmp_path / "solver.xml")
    result = run_by_name(tmp_path, "solver.xml", 2)
    pressure = result.point_data["Darcy_pressure"].reshape(-1)
    np.testing.assert_allclose(pressure[result.points[:, 0] < 1.5], 2, atol=1e-9)
    np.testing.assert_allclose(pressure[result.points[:, 0] > 1.5], 11, atol=1e-9)


def test_interior_node_equation_offset(tmp_path):
    root, points, _, ids = make_interior_node_case(tmp_path)
    heat = copy.deepcopy(root.find("Add_equation"))
    heat.set("type", "heatS")
    for child in list(heat):
        if child.tag.startswith("Darcy_") or child.tag == "Fluid_density":
            heat.remove(child)
    ET.SubElement(heat, "Density").text = "0"
    ET.SubElement(heat, "Conductivity").text = "1"
    heat.find("Output/Darcy_pressure").tag = "Temperature"
    bc = heat.find("Add_BC")
    bc.find("Time_dependence").text = "Steady"
    bc.remove(bc.find("Temporal_and_spatial_values_file_path"))
    ET.SubElement(bc, "Value").text = "13"
    root.insert(2, heat)
    ET.ElementTree(root).write(tmp_path / "solver.xml")
    result = run_by_name(tmp_path, "solver.xml", 2)
    order = [np.argmin(np.linalg.norm(result.points - points[i-1], axis=1)) for i in ids]
    np.testing.assert_allclose(result.point_data["Darcy_pressure"].reshape(-1)[order], [2, 7], atol=1e-10)
    np.testing.assert_allclose(result.point_data["Temperature"], 13, atol=1e-8)


@pytest.mark.parametrize("kind", ["both_sources", "neither_source", "missing_file", "duplicate_set",
                                 "missing_mesh_reference", "missing_set_reference", "inactive_domain",
                                 "empty_name", "ambiguous_mesh", "coupling_interface"])
def test_interior_node_invalid_definition(tmp_path, kind):
    root, _, _, _ = make_interior_node_case(tmp_path)
    mesh = root.find("Add_mesh")
    node_set = mesh.find("Add_node_set")
    bc = root.find("Add_equation/Add_BC")
    if kind in ("both_sources", "missing_file"):
        ET.SubElement(node_set, "Node_IDs_file_path").text = "missing.dat"
    if kind in ("neither_source", "missing_file"):
        node_set.remove(node_set.find("Node_IDs"))
    if kind == "duplicate_set":
        mesh.append(copy.deepcopy(node_set))
    if kind == "empty_name":
        node_set.set("name", " ")
    if kind == "ambiguous_mesh":
        root.insert(2, copy.deepcopy(mesh))
        root.find("GeneralSimulationParameters/Number_of_time_steps").text = "0"
    if kind == "coupling_interface":
        ET.SubElement(ET.SubElement(bc, "Coupling_interface"), "svZeroDSolver_block").text = "unused"
    if kind == "missing_mesh_reference":
        bc.remove(bc.find("Mesh_name"))
    if kind == "missing_set_reference":
        bc.remove(bc.find("Node_set"))
    if kind == "inactive_domain":
        ET.SubElement(mesh, "Domain").text = "1"
        equation = root.find("Add_equation")
        domain = ET.SubElement(equation, "Domain", id="2")
        for child in list(equation):
            if child.tag.startswith("Darcy_") or child.tag in ("Fluid_density", "Source_term"):
                equation.remove(child)
                domain.append(child)
    ET.ElementTree(root).write(tmp_path / "solver.xml")
    run_by_name(tmp_path, "solver.xml", 2, expected_error=(
        "node set.*empty name" if kind == "empty_name" else "Ambiguous mesh" if kind == "ambiguous_mesh" else
        "exactly one|cannot open node-ID|Duplicate node set|requires Mesh_name|without active|node.*Coupling_interface"))


@pytest.mark.parametrize("kind", ["remeshing", "merged_ids"])
def test_interior_node_topology_validation(tmp_path, kind):
    root, points, cells, _ = make_interior_node_case(tmp_path, 3)
    if kind == "remeshing":
        ET.SubElement(root.find("GeneralSimulationParameters"), "Simulation_requires_remeshing").text = "true"
    else:
        first = add_test_boundary_face(root, tmp_path, points, cells, "first", 0)
        second = add_test_boundary_face(root, tmp_path, points, cells, "second", 1)
        root.find("Add_mesh/Add_node_set/Node_IDs").text = " ".join(map(str, np.r_[first, second]))
        projection = ET.SubElement(root, "Add_projection", name="first")
        ET.SubElement(projection, "Project_from_face").text = "second"
        ET.SubElement(projection, "Projection_tolerance").text = "2"
    ET.ElementTree(root).write(tmp_path / "solver.xml")
    run_by_name(tmp_path, "solver.xml", 2, expected_error="node.*(remeshing|merged|duplicate)")
