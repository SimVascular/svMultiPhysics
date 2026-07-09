import numpy as np

import pytest
import os
import shutil
import subprocess
import sys
import meshio

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


def run_by_name(folder, name, t_max, n_proc=1):
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
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)

    # run simulation (PETSc tests use a dedicated build; see cpp_exec_p)
    exe = cpp_exec_p if "petsc" in folder else cpp_exec
    cmd = " ".join(
        [
            "mpirun",
            OVERSUBSCRIBE_FLAG if n_proc > 1 else "",
            "-np",
            str(n_proc),
            exe,
            name,
        ]
    )

    # Run the command while capturing the return code and stderr output. This
    # way, if something goes wrong, we can raise an appropriate error message.
    completed = subprocess.run(
        cmd, cwd=folder, shell=True, stderr=subprocess.PIPE, text=True
    )

    # Print the captured stderr to console, so it is visible. Notice that this
    # will print stderr after stdout, so they might be out of order (printing
    # them in order while capturing is apparently not easy through subprocess).
    if completed.stderr:
        print(completed.stderr, end="", file=sys.stderr)

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
):
    """
    Run a test case and compare it to a stored reference solution
    Args:
        folder: location from which test will be executed
        fields: array fields to compare (e.g. ["Pressure", "Velocity"])
        n_proc: number of processors
        t_max: time step to compare
        name_inp: name of svMultiPhysics input file (.xml)
        name_ref: name of refence file (.vtu)
    """
    # default reference name
    if not name_ref:
        name_ref = "result_" + str(t_max).zfill(3) + ".vtu"

    # run simulation
    folder = os.path.join("cases", base_folder, test_folder)
    res = run_by_name(folder, name_inp, t_max, n_proc)

    # read reference
    fname = os.path.join(folder, name_ref)
    ref = meshio.read(fname)

    # check results
    msg = ""
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
