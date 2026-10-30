from .conftest import run_by_name, run_with_reference, skip_if_no_petsc, skip_if_no_trilinos
import os
import subprocess
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path
import pytest

# Common folder for all tests in this file
base_folder = "fluid"

# Fields to test
fields = ["Velocity", "Pressure", "Traction", "WSS", "Vorticity", "Divergence"]


def test_pipe_RCR_3d(n_proc):
    test_folder = "pipe_RCR_3d"
    t_max = 2
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)

def test_pipe_restart_with_spaces(n_proc, tmp_path):
    case = tmp_path / "pipe with spaces"
    shutil.copytree(Path(__file__).parent / "cases/fluid/pipe_RCR_3d", case,
                    ignore=shutil.ignore_patterns("*-procs", "histor.dat"))
    tree = ET.parse(case / "solver.xml")
    general = tree.getroot().find("GeneralSimulationParameters")
    general.find("Number_of_time_steps").text = "1"
    general.find("Increment_in_saving_VTK_files").text = "1"
    ET.SubElement(general, "Restart_file_name").text = "restart snapshot"
    tree.write(case / "first step.xml")
    run_by_name(str(case), "first step.xml", 1, n_proc)
    assert (case / f"{n_proc}-procs/restart snapshot_last.bin").is_file()

    general.find("Continue_previous_simulation").text = "true"
    general.find("Number_of_time_steps").text = "2"
    tree.write(case / "second step.xml")
    run_with_reference("", str(case), fields, n_proc, 2, name_inp="second step.xml", clean=False)

def test_pipe_output_error(n_proc, tmp_path):
    case = tmp_path / "unwritable output"
    shutil.copytree(Path(__file__).parent / "cases/fluid/pipe_RCR_3d", case,
                    ignore=shutil.ignore_patterns("*-procs", "histor.dat"))
    (case / f"{n_proc}-procs").write_text("preserve this file")
    with pytest.raises(RuntimeError, match=f"{n_proc}-procs"):
        run_by_name(str(case), "solver.xml", 2, n_proc, clean=False, timeout=30)
    assert (case / f"{n_proc}-procs").read_text() == "preserve this file"

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

def test_pipe_RCR_genBC(n_proc, tmp_path):
    case = tmp_path / "genBC with spaces"
    shutil.copytree(Path(__file__).parent / "cases/fluid/pipe_RCR_genBC", case,
                    ignore=shutil.ignore_patterns("*-procs", "obj", "genBC.exe", "AllData", "InitialData", "GenBC.int"))
    subprocess.run(["make"], cwd=case / "genBC", check=True)
    executable = case / "genBC/coupled solver.exe"
    (case / "genBC/genBC.exe").rename(executable)
    tree = ET.parse(case / "solver.xml")
    tree.find(".//ZeroD_code_file_path").text = str(executable)
    tree.write(case / "solver.xml")
    run_with_reference("", str(case), fields, n_proc, 2)

def test_pipe_RCR_genBC_missing_executable(n_proc, tmp_path):
    case = tmp_path / "missing genBC"
    shutil.copytree(Path(__file__).parent / "cases/fluid/pipe_RCR_genBC", case,
                    ignore=shutil.ignore_patterns("*-procs", "genBC"))
    with pytest.raises(RuntimeError, match="Could not start genBC"):
        run_by_name(str(case), "solver.xml", 2, n_proc)

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
