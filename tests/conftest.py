import numpy as np

import pytest
import os
import shutil
import platform
import subprocess
import meshio

is_not_Darwin = True
if platform.system() == "Darwin": is_not_Darwin = False

this_file_dir = os.path.abspath(os.path.dirname(__file__))
cpp_exec = os.path.join(this_file_dir, "..", "build", "svMultiPhysics-build", "bin", "svmultiphysics")
cpp_exec_p = os.path.join(this_file_dir, "..", "build-petsc", "svMultiPhysics-build", "bin", "svmultiphysics")

# Relative tolerances for each tested field
RTOL = {
    "Membrane_potential": 1.0e-10,
    "Calcium": 1.0e-10,
    "Cauchy_stress": 1.0e-4,
    "Concentration": 1.0e-10,
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
}

# Number of processors to test
PROCS = [1, 3, 4]


# Fixture to parametrize the number of processors for all tests
@pytest.fixture(params=PROCS)
def n_proc(request):
    return request.param


def _run_simulation(folder, name_inp, n_proc):
    """Run the solver once, removing any previous output directory first."""
    dir_path = os.path.join(folder, str(n_proc) + "-procs")
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)

    exec_path = cpp_exec_p if "petsc" in folder else cpp_exec
    cmd = " ".join([
        "mpirun",
        "--oversubscribe" if n_proc > 1 else "",
        "-np", str(n_proc),
        exec_path,
        name_inp,
    ])
    subprocess.call(cmd, cwd=folder, shell=True)


def _read_result(folder, n_proc, name_result, t_max):
    """Read a single VTU result file from the n_proc output directory."""
    result_file = name_result if name_result else "result_" + str(t_max).zfill(3) + ".vtu"
    fname = os.path.join(folder, str(n_proc) + "-procs", result_file)
    if not os.path.exists(fname):
        raise RuntimeError("No svMultiPhysics output: " + fname)
    return meshio.read(fname)


def _compare_fields(res, ref, fields):
    """Compare fields between a result and reference mesh. Returns an error string."""
    msg = ""
    for f in fields:
        if f not in res.point_data:
            raise ValueError("Field " + f + " not in simulation result")
        a = res.point_data[f]

        if f not in ref.point_data:
            raise ValueError("Field " + f + " not in reference result")
        b = ref.point_data[f]

        if len(a.shape) == 2:
            if a.shape[1] == 2 and b.shape[1] == 3:
                assert not np.any(b[:, 2])
                b = b[:, :2]

        if f not in RTOL:
            raise ValueError("No tolerance defined for field " + f)
        rtol = RTOL[f]

        a_fl = a.flatten()
        b_fl = b.flatten()
        rel_diff = np.abs(a_fl - b_fl) - rtol - rtol * np.abs(b_fl)

        close = rel_diff <= 0.0
        if not np.all(close):
            wrong = 1 - np.sum(close) / close.size
            i_max = rel_diff.argmax()
            max_rel = rel_diff[i_max]
            max_abs = np.abs(a_fl[i_max] - b_fl[i_max])

            msg += "Test failed in field " + f + "."
            msg += " Results differ by more than rtol=" + str(rtol)
            msg += " in {:.1%}".format(wrong)
            msg += " of results."
            msg += " Max. rel. difference is"
            msg += " {:.1e}".format(max_rel)
            msg += " (abs. {:.1e}".format(max_abs) + ")\n"
    return msg


def run_by_name(folder, name, t_max, n_proc=1, name_result=None):
    """Run a test case and return results (legacy single-file interface)."""
    _run_simulation(folder, name, n_proc)
    return _read_result(folder, n_proc, name_result, t_max)


def run_with_reference(
    base_folder,
    test_folder,
    fields,
    n_proc=1,
    t_max=1,
    name_ref=None,
    name_inp="solver.xml",
    name_result=None,
    comparisons=None,
):
    """
    Run a test case once and compare one or more output files to stored references.

    Args:
        fields:      fields to compare for the primary output file
        n_proc:      number of processors
        t_max:       timestep index used for default file naming
        name_inp:    solver input XML filename
        name_ref:    reference VTU filename (default: result_{t_max:03d}.vtu)
        name_result: result VTU filename inside {n_proc}-procs/ (default: same as name_ref)
        comparisons: list of dicts for multi-file comparison, each with keys:
                       "fields"      — list of field names to compare
                       "name_ref"    — reference VTU filename
                       "name_result" — result VTU filename (optional, defaults to name_ref)
                     When provided, fields/name_ref/name_result are ignored.
    """
    folder = os.path.join("cases", base_folder, test_folder)

    if not is_not_Darwin and ("petsc" in folder or "trilinos" in folder):
        return

    # Build the comparison list from either the new or legacy parameters
    if comparisons is None:
        if not name_ref:
            name_ref = "result_" + str(t_max).zfill(3) + ".vtu"
        comparisons = [{"fields": fields, "name_ref": name_ref, "name_result": name_result}]

    # Run the simulation once
    _run_simulation(folder, name_inp, n_proc)

    # Compare each requested output file against its reference
    msg = ""
    for comp in comparisons:
        res = _read_result(folder, n_proc, comp.get("name_result"), t_max)
        ref = meshio.read(os.path.join(folder, comp["name_ref"]))
        msg += _compare_fields(res, ref, comp["fields"])

    if msg:
        raise AssertionError(msg)
