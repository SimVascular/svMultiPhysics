from .conftest import (add_test_boundary_face, make_interior_node_case, run_by_name,
                       run_with_reference, skip_if_no_petsc, skip_if_no_trilinos)
import numpy as np
import os
import subprocess
import xml.etree.ElementTree as ET

# Common folder for all tests in this file
base_folder = "fluid"

# Fields to test
fields = ["Velocity", "Pressure", "Traction", "WSS", "Vorticity", "Divergence"]


def test_pipe_RCR_3d(n_proc):
    test_folder = "pipe_RCR_3d"
    t_max = 2
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)

def test_pipe_RCR_3d_fourier_coeff(n_proc):
    test_folder = "pipe_RCR_3d_fourier_coeff"
    t_max = 2
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)

@skip_if_no_petsc
def test_pipe_RCR_3d_petsc(n_proc):
    test_folder = "pipe_RCR_3d_petsc"
    t_max = 2
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)

@skip_if_no_trilinos
def test_pipe_RCR_3d_trilinos_ilut(n_proc):
    test_folder = "pipe_RCR_3d_trilinos_ilut"
    t_max = 2
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)

@skip_if_no_trilinos
def test_pipe_RCR_3d_trilinos_bj(n_proc):
    test_folder = "pipe_RCR_3d_trilinos_bj"
    t_max = 2
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)

def test_pipe_RCR_weak_dir_3d(n_proc):
    test_folder = "pipe_RCR_weak_dir_3d"
    t_max = 2
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)

def test_pipe_RCR_genBC(n_proc):
    test_folder = "pipe_RCR_genBC"
    t_max = 2

    # Remove old genBC output
    os.chdir(os.path.join("cases", base_folder, test_folder))
    for name in ["AllData", "InitialData", "GenBC.int"]:
        if os.path.isfile(name):
            os.remove(name)

    # Compile genBC
    os.chdir("genBC")
    subprocess.run(["make", "clean"], check=True)
    subprocess.run(["make"], check=True)

    # Change back to original directory
    os.chdir("../../../..")

    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)

def test_pipe_RCR_sv0D(n_proc):
    test_folder = "pipe_RCR_sv0D"
    t_max = 2
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)


def test_driven_cavity_2d(n_proc):
    test_folder = "driven_cavity_2d"
    t_max = 2
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)


def test_driven_cavity_2d_porous(n_proc):
    test_folder = "driven_cavity_2d_porous"
    t_max = 2
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)


def test_dye_AD(n_proc):
    test_folder = "dye_AD"
    run_with_reference(base_folder, test_folder, ['Pressure', 'Velocity', 'Concentration'], n_proc)


def test_precomputed_dye_AD(n_proc):
    test_folder = "precomputed_dye_AD"
    run_with_reference(base_folder, test_folder, ['Velocity', 'Concentration'], n_proc)


def test_newtonian(n_proc):
    test_folder = "newtonian"
    run_with_reference(base_folder, test_folder, fields, n_proc)


def test_casson(n_proc):
    test_folder = "casson"
    run_with_reference(base_folder, test_folder, fields, n_proc)


def test_carreau_yasuda(n_proc):
    test_folder = "carreau_yasuda"
    run_with_reference(base_folder, test_folder, fields, n_proc)


def test_iliac_artery(n_proc):
    test_folder = "iliac_artery"
    run_with_reference(base_folder, test_folder, fields, n_proc)

@skip_if_no_trilinos
def test_iliac_artery_trilinos_gmres_ilut(n_proc):
    test_folder = "iliac_artery_trilinos_gmres_ilut"
    run_with_reference(base_folder, test_folder, fields, n_proc)

def test_quadratic_tet10(n_proc):
    test_folder = "quadratic_tet10"
    t_max = 3
    fields = ["Velocity", "Pressure", "Vorticity", "Divergence"]
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max) 

def test_interior_node_with_coupled_face(tmp_path):
    root, points, cells, ids = make_interior_node_case(tmp_path, 3)
    add_test_boundary_face(root, tmp_path, points, cells)
    equation = root.find("Add_equation")
    equation.set("type", "fluid")
    for child in list(equation):
        if child.tag.startswith("Darcy_") or child.tag in ("Fluid_density", "Source_term"):
            equation.remove(child)
    ET.SubElement(equation, "Density").text = "1"
    ET.SubElement(ET.SubElement(equation, "Viscosity", model="Constant"), "Value").text = "1"
    equation.find("LS").set("type", "GMRES")
    equation.find("LS/Linear_algebra/Preconditioner").text = "fsils"
    equation.find("Output/Darcy_pressure").tag = "Velocity"
    bc = equation.find("Add_BC")
    bc.remove(bc.find("Temporal_and_spatial_values_file_path"))
    bc.find("Time_dependence").text = "Steady"
    ET.SubElement(bc, "Value").text = "0.02"
    outlet = ET.SubElement(equation, "Add_BC", name="surface")
    ET.SubElement(outlet, "Type").text = "Neumann"
    ET.SubElement(outlet, "Time_dependence").text = "RCR"
    rcr = ET.SubElement(outlet, "RCR_values")
    for tag, value in [("Capacitance", "1"), ("Proximal_resistance", "1"),
                       ("Distal_resistance", "1"), ("Distal_pressure", "0"), ("Initial_pressure", "0")]:
        ET.SubElement(rcr, tag).text = value
    ET.ElementTree(root).write(tmp_path / "solver.xml")
    result = run_by_name(tmp_path, "solver.xml", 2)
    order = [np.argmin(np.linalg.norm(result.points - points[i-1], axis=1)) for i in ids]
    np.testing.assert_allclose(result.point_data["Velocity"][order], 0.02, atol=1e-10)
    assert np.isfinite(result.point_data["Velocity"]).all()
