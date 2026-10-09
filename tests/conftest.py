import numpy as np

import math
import itertools
import pytest
import os
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
import meshio
import vtk
from vtk.util.numpy_support import numpy_to_vtk

this_file_dir = os.path.abspath(os.path.dirname(__file__))
cpp_exec = os.path.join(this_file_dir, "..", "build", "svMultiPhysics-build", "bin", "svmultiphysics")
cpp_exec_p = os.path.join(this_file_dir, "..", "build-petsc", "svMultiPhysics-build", "bin", "svmultiphysics")


def read_cmake_cache_variable(cache_path, key):
    """Return the value of `key` in a CMakeCache.txt (lines are KEY:TYPE=VALUE),
    or None if the cache file or the key cannot be found."""
    try:
        with open(cache_path) as cache:
            for line in cache:
                line = line.strip()
                if line.startswith(key + ":"):
                    return line.split("=", 1)[1].strip() if "=" in line else ""
    except OSError:
        return None
    return None


def cmake_cache_path_for(exe_path):
    """CMakeCache.txt of the build that produced `exe_path`, which lives at
    <build>/svMultiPhysics-build/bin/svmultiphysics."""
    return os.path.join(os.path.dirname(os.path.dirname(exe_path)), "CMakeCache.txt")


# Whether svMultiPhysics was built with PETSc / Trilinos, read from the
# CMakeCache.txt of the corresponding build (PETSc tests use the separate build
# at cpp_exec_p; Trilinos is linked into the main build at cpp_exec). PETSc is
# enabled when SV_PETSC_DIR is a non-empty path, Trilinos when SV_USE_TRILINOS
# is ON. A missing cache (e.g. the build does not exist) means "not available".
HAS_PETSC = bool(read_cmake_cache_variable(cmake_cache_path_for(cpp_exec_p), "SV_PETSC_DIR"))
HAS_TRILINOS = (
    read_cmake_cache_variable(cmake_cache_path_for(cpp_exec), "SV_USE_TRILINOS") or ""
).upper() in ("ON", "1", "TRUE", "YES")

# Reusable markers to decorate PETSc/Trilinos tests at their definition site.
skip_if_no_petsc = pytest.mark.skipif(
    not HAS_PETSC,
    reason="svMultiPhysics not built with PETSc (SV_PETSC_DIR empty in CMakeCache.txt)",
)
skip_if_no_trilinos = pytest.mark.skipif(
    not HAS_TRILINOS,
    reason="svMultiPhysics not built with Trilinos (SV_USE_TRILINOS=OFF in CMakeCache.txt)",
)


def _detect_oversubscribe_flag():
    """Return the mpirun flag needed to allow more ranks than physical cores.

    Open MPI requires ``--oversubscribe``; Intel MPI / MPICH (Hydra) allow
    oversubscription by default and reject the unknown flag, so no flag is used.
    """
    try:
        proc = subprocess.run(
            ["mpirun", "--version"], capture_output=True, text=True, check=False
        )
        version = (proc.stdout or "") + (proc.stderr or "")
    except FileNotFoundError:
        version = ""

    if "Open MPI" in version or "OpenRTE" in version:
        return "--oversubscribe"
    return ""


# Detected once at import; empty string for Intel MPI / MPICH.
OVERSUBSCRIBE_FLAG = _detect_oversubscribe_flag()

# Relative tolerances for each tested field
RTOL = {
    "Membrane_potential": 1.0e-10,
    "Calcium": 1.0e-10,
    "Cauchy_stress": 1.0e-4,
    "Concentration": 1.0e-10,
    "Darcy_pressure": 1.0e-10,
    "Def_grad": 1.0e-10,
    "Divergence": 1.0e-9,
    "Displacement": 1.0e-10,
    "Jacobian": 1.0e-10,
    "Pressure": 1.0e-6,
    "Stress": 1.0e-4,
    "Strain": 1.0e-10,
    "Temperature": 1.0e-10,
    "Traction": 1.0e-6,
    "Velocity": 1.0e-7,
    "VonMises_stress": 1.0e-3,
    "Vorticity": 1.0e-7,
    "WSS": 1.0e-8,
    "Fiber_stretch": 1.0e-10,
    "Fiber_stretch_rate": 1.0e-10,
    "Active_tension_fibers": 1.0e-10,
    "Active_tension_sheets": 1.0e-10,
    "Active_tension_normal": 1.0e-10,
}

# Relative tolerance for the TimeValue field data. The solver accumulates the
# time as time += dt while the test computes t_max * dt, so the two differ only
# by floating point round-off.
RTOL_TIME_VALUE = 1.0e-12

# Number of processors to test
PROCS = [1, 3, 4]


def read_time_step_size(name_inp):
    """
    Read the time step size from a svMultiPhysics input file
    Args:
        name_inp: path to the svMultiPhysics input file (.xml)

    Returns:
    Time step size
    """
    general = ET.parse(name_inp).getroot().find("GeneralSimulationParameters")
    if general is None:
        raise ValueError("No GeneralSimulationParameters in " + name_inp)

    time_step_size = general.find("Time_step_size")
    if time_step_size is None:
        raise ValueError("No Time_step_size in " + name_inp)

    return float(time_step_size.text)


# Fixture to parametrize the number of processors for all tests
@pytest.fixture(params=PROCS)
def n_proc(request):
    return request.param


def add_test_boundary_face(root, folder, points, cells, name="surface", x=0):
    """Write one exterior triangle of a generated tetrahedral mesh as a VTP."""
    for element, cell in enumerate(cells):
        nodes = cell[np.isclose(points[cell, 0], x)]
        if len(nodes) == 3:
            break
    else:
        raise ValueError("No boundary triangle found")
    vtk_points = vtk.vtkPoints()
    vtk_points.SetData(numpy_to_vtk(points[nodes], deep=True))
    polygons = vtk.vtkCellArray()
    polygons.InsertNextCell(3, (0, 1, 2))
    surface = vtk.vtkPolyData()
    surface.SetPoints(vtk_points)
    surface.SetPolys(polygons)
    for attributes, label, values in [
        (surface.GetPointData(), "GlobalNodeID", nodes + 1),
        (surface.GetCellData(), "GlobalElementID", [element + 1]),
    ]:
        array = numpy_to_vtk(np.asarray(values, dtype=np.int32), deep=True)
        array.SetName(label)
        attributes.AddArray(array)
    writer = vtk.vtkXMLPolyDataWriter()
    path = folder / f"{name}.vtp"
    writer.SetFileName(str(path))
    writer.SetInputData(surface)
    if writer.Write() != 1:
        raise OSError(f"Could not write {path}")
    face = ET.SubElement(root.find("Add_mesh"), "Add_face", name=name)
    ET.SubElement(face, "Face_file_path").text = f"{name}.vtp"
    return nodes + 1


def run_by_name(folder, name, t_max, n_proc=1, *, expected_error=None, clean=True, timeout=120):
    """
    Run a test case and return results
    Args:
        folder: location from which test will be executed
        name: name of svMultiPhysics input file (.xml)
        t_max: time step to compare
        n_proc: number of processors

    Returns:
    Simulation results
    """

    # remove old results folders if they exist
    dir_path = os.path.join(folder, str(n_proc) + "-procs")
    if clean and os.path.exists(dir_path):
        shutil.rmtree(dir_path)

    # run simulation (PETSc tests use a dedicated build; see cpp_exec_p)
    exe = cpp_exec_p if "petsc" in str(folder) else cpp_exec
    cmd = ["mpirun"]
    if n_proc > 1 and OVERSUBSCRIBE_FLAG:
        cmd.append(OVERSUBSCRIBE_FLAG)
    cmd.extend(["-np", str(n_proc), exe, name])

    # Run the command while capturing the return code and stderr output. This
    # way, if something goes wrong, we can raise an appropriate error message.
    completed = subprocess.run(
        cmd, cwd=folder, capture_output=True, text=True, timeout=timeout
    )
    if completed.stdout:
        print(completed.stdout, end="")

    # Print the captured stderr to console, so it is visible. Notice that this
    # will print stderr after stdout, so they might be out of order (printing
    # them in order while capturing is apparently not easy through subprocess).
    if completed.stderr:
        print(completed.stderr, end="", file=sys.stderr)

    if expected_error is not None:
        assert completed.returncode != 0, "Invalid input unexpectedly succeeded"
        assert re.search(expected_error, completed.stdout + completed.stderr, re.I)
        return

    # If something went wrong, raise an error with the captured stderr output in
    # the message.
    if completed.returncode != 0:
        raise RuntimeError(
            "Exit code {}: {}\n".format(completed.returncode, completed.stderr)
        )

    # read results
    fname = os.path.join(
        folder, str(n_proc) + "-procs", "result_" + str(t_max).zfill(3) + ".vtu"
    )
    if not os.path.exists(fname):
        raise RuntimeError("No svMultiPhysics output: " + fname)
    return meshio.read(fname)


def run_with_reference(
    base_folder,
    test_folder,
    fields,
    n_proc=1,
    t_max=1,
    name_ref=None,
    name_inp="solver.xml",
    check_time_value=True,
):
    """
    Run a test case and compare it to a stored reference solution
    Args:
        folder: location from which test will be executed
        fields: array fields to compare (e.g. ["Pressure", "Velocity"])
        n_proc: number of processors
        t_max: time step to compare
        name_inp: name of svMultiPhysics input file (.xml)
        name_ref: name of reference file (.vtu)
        check_time_value: whether to compare the TimeValue field data against
            the time reached at time step t_max
    """
    # default reference name
    if not name_ref:
        name_ref = "result_" + str(t_max).zfill(3) + ".vtu"

    # run simulation
    folder = os.path.join(this_file_dir, "cases", base_folder, test_folder)
    res = run_by_name(folder, name_inp, t_max, n_proc, timeout=None)

    # read reference
    fname = os.path.join(folder, name_ref)
    ref = meshio.read(fname)

    # check results
    msg = ""

    # check the time attached to the result as field data. This assumes a
    # constant time step size, which does not hold if the case sets
    # Number_of_initialization_time_steps.
    if check_time_value:
        if "TimeValue" not in res.field_data.keys():
            raise ValueError("Field data TimeValue not in simulation result")

        time_value = res.field_data["TimeValue"][0]
        time_expected = t_max * read_time_step_size(os.path.join(folder, name_inp))

        if not math.isclose(time_value, time_expected, rel_tol=RTOL_TIME_VALUE):
            msg += "Test failed in field data TimeValue."
            msg += " Result is " + str(time_value)
            msg += " instead of " + str(time_expected) + ".\n"

    for f in fields:
        # extract field
        if f not in res.point_data.keys():
            raise ValueError("Field " + f + " not in simulation result")
        a = res.point_data[f]

        if f not in ref.point_data.keys():
            raise ValueError("Field " + f + " not in reference result")
        b = ref.point_data[f]

        # truncate last dimension if solution is 2D but reference is 3D
        if len(a.shape) == 2:
            if a.shape[1] == 2 and b.shape[1] == 3:
                assert not np.any(b[:, 2])
                b = b[:, :2]

        # pick tolerance for current field
        if f not in RTOL:
            raise ValueError("No tolerance defined for field " + f)
        rtol = RTOL[f]

        # relative difference (as computed in np.isclose)
        # note that we consider rtol as absolute zero (and as relative tolerance)
        a_fl = a.flatten()
        b_fl = b.flatten()
        rel_diff = np.abs(a_fl - b_fl) - rtol - rtol * np.abs(b_fl)

        # throw error if not all results are within relative tolerance
        close = rel_diff <= 0.0
        if not np.all(close):
            # portion of individual results that are above the tolerance
            wrong = 1 - np.sum(close) / close.size

            # location of maximum relative difference
            i_max = rel_diff.argmax()

            # maximum relative difference
            max_rel = rel_diff[i_max]

            # maximum absolute difference at same location
            max_abs = np.abs(a_fl[i_max] - b_fl[i_max])

            # throw error message for pytest
            msg += "Test failed in field " + f + "."
            msg += " Results differ by more than rtol=" + str(rtol)
            msg += " in {:.1%}".format(wrong)
            msg += " of results."
            msg += " Max. rel. difference is"
            msg += " {:.1e}".format(max_rel)
            msg += " (abs. {:.1e}".format(max_abs) + ")\n"
    # check all fields first and then throw error if any failed
    if msg:
        raise AssertionError(msg)


def make_interior_node_case(folder, nsd=2):
    """Generate a small simplex mesh and editable input without surface files."""
    grid = list(itertools.product(range(4), repeat=nsd))
    lookup = {point: i for i, point in enumerate(grid)}
    points = np.zeros((len(grid), 3))
    points[:, :nsd] = np.asarray(grid) / 3.0
    cells = []
    for corner in itertools.product(range(3), repeat=nsd):
        for axes in itertools.permutations(range(nsd)):
            vertex = list(corner)
            simplex = [lookup[tuple(vertex)]]
            for axis in axes:
                vertex[axis] += 1
                simplex.append(lookup[tuple(vertex)])
            if np.linalg.det((points[simplex[1:], :nsd] - points[simplex[0], :nsd]).T) < 0:
                simplex[0], simplex[1] = simplex[1], simplex[0]
            cells.append(simplex)
    cells = np.asarray(cells)
    ids = np.array([lookup[(2,) * nsd] + 1, lookup[(1,) * nsd] + 1])
    vtk_points = vtk.vtkPoints()
    vtk_points.SetData(numpy_to_vtk(points, deep=True))
    vtk_cells = vtk.vtkCellArray()
    for cell in cells:
        vtk_cells.InsertNextCell(len(cell), cell)
    volume = vtk.vtkUnstructuredGrid()
    volume.SetPoints(vtk_points)
    volume.SetCells(vtk.VTK_TRIANGLE if nsd == 2 else vtk.VTK_TETRA, vtk_cells)
    node_ids = numpy_to_vtk(np.arange(len(points)) + 1001, deep=True)
    node_ids.SetName("GlobalNodeID")
    volume.GetPointData().AddArray(node_ids)
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(str(folder / "mesh.vtu"))
    writer.SetInputData(volume)
    if writer.Write() != 1:
        raise OSError(f"Could not write {folder / 'mesh.vtu'}")
    root = ET.Element("svMultiPhysicsFile", version="0.1")
    general = ET.SubElement(root, "GeneralSimulationParameters")
    for tag, value in {
        "Continue_previous_simulation": "false",
        "Number_of_spatial_dimensions": nsd,
        "Number_of_time_steps": 2,
        "Time_step_size": 0.1,
        "Spectral_radius_of_infinite_time_step": 0,
        "Save_results_to_VTK_format": "true",
        "Name_prefix_of_saved_VTK_files": "result",
        "Increment_in_saving_VTK_files": 1,
        "Start_saving_after_time_step": 0,
        "Increment_in_saving_restart_files": 1,
    }.items():
        ET.SubElement(general, tag).text = str(value)
    mesh = ET.SubElement(root, "Add_mesh", name="volume")
    ET.SubElement(mesh, "Mesh_file_path").text = "mesh.vtu"
    node_set = ET.SubElement(mesh, "Add_node_set", name="interior")
    ET.SubElement(node_set, "Node_IDs").text = " ".join(map(str, ids))
    equation = ET.SubElement(root, "Add_equation", type="darcy")
    for tag, value in {
        "Coupled": "true", "Min_iterations": 1, "Max_iterations": 6, "Tolerance": "1e-11",
        "Darcy_permeability": 1, "Darcy_fluid_viscosity": 1, "Darcy_compressibility": 0,
        "Fluid_density": 1, "Source_term": 0.25,
    }.items():
        ET.SubElement(equation, tag).text = str(value)
    output = ET.SubElement(equation, "Output", type="Spatial")
    ET.SubElement(output, "Darcy_pressure").text = "true"
    ls = ET.SubElement(equation, "LS", type="CG")
    algebra = ET.SubElement(ls, "Linear_algebra", type="fsils")
    ET.SubElement(algebra, "Preconditioner").text = "rcs"
    ET.SubElement(ls, "Tolerance").text = "1e-12"
    ET.SubElement(ls, "Max_iterations").text = "500"
    bc = ET.SubElement(equation, "Add_BC", name="interior_pressure")
    for tag, value in {
        "Mesh_name": "volume", "Node_set": "interior", "Type": "Dirichlet",
        "Time_dependence": "General",
        "Temporal_and_spatial_values_file_path": "values.dat",
    }.items():
        ET.SubElement(bc, tag).text = value
    # Records intentionally reverse the node-set order.
    (folder / "values.dat").write_text(f"1 2 2\n0 1\n{ids[1]} 7 7\n{ids[0]} 2 2\n")
    return root, points, cells, ids
