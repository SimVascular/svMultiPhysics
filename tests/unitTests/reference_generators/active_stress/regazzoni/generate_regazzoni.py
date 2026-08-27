#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause

"""Generate the canonical Regazzoni reference from pinned authors' C++ code."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
from typing import Iterable


PINNED_COMMIT = "26f05df28891df7b3c69f16bb136cdced6b63c4d"
DT_MS = 1.0
DT_S = DT_MS / 1000.0
UPDATE_COUNT = 600
CHECKPOINT_STEPS = (0, 30, 99, 157, 299, 599)
STATE_COUNT = 20

BASELINE_CALCIUM_MM = 1.0e-4
PEAK_CALCIUM_MM = 9.0e-4
CALCIUM_RISE_MS = 20.0
CALCIUM_DECAY_MS = 50.0
CALCIUM_ONSET_MS = 10.0

RESTING_SL_UM = 2.2
MINIMUM_SL_UM = 2.134
SHORTENING_ONSET_MS = 30.0
MINIMUM_TIME_MS = 150.0
RECOVERY_TIME_MS = 350.0

REFERENCE_SOURCE_PATHS = (
    "models_cpp/model_RDQ20_MF.cpp",
    "models_cpp/model_RDQ20_MF.hpp",
    "models_cpp/sarcomere.cpp",
    "models_cpp/sarcomere.hpp",
    "params/params_RDQ20-MF_human_body-temperature.json",
)

COMPATIBLE_PREAMBLE = """# formulation source: Regazzoni, Dede', Quarteroni (2020), PLOS Comput. Biol.,
#   doi:10.1371/journal.pcbi.1008294. Supporting Information S3 uses Forward
#   Euler for RU dynamics and an exponential integrator for XB dynamics.
# oracle discretization: authors' C++ reference implementation,
#   https://github.com/FrancescoRegazzoni/cardiac-activation,
#   commit 26f05df28891df7b3c69f16bb136cdced6b63c4d. That implementation uses
#   Forward-Euler RU substeps and implicit Euler for XB dynamics, which is the
#   discretization followed by svMultiPhysics. The following code patch was
#   applied before running the oracle:
#   - XB_A initialization: the reference code declares
#       Eigen::Matrix<double,4,4> XB_A;
#     and sets only 10 of 16 entries, leaving 6 off-diagonal entries as
#     uninitialized memory (undefined behavior). The oracle was generated with
#     XB_A.setZero() inserted before entry population, matching the
#     svMultiPhysics implementation which uses Eigen::Matrix<double,4,4>::Zero().
# SL protocol: the raised-cosine ramp below is a custom oracle driver protocol;
#   the reference main.cpp uses an exponential SL sigmoid. At outer-step start
#   t, the driver supplies SL(t) and the forward-difference velocity
#     dSL/dt = [SL(t+dt) - SL(t)] / dt.
# State ordering: s0..s15 are RU probabilities P(TL,TC,TR,CC), with
#   index = 8*TL + 4*TC + 2*TR + CC (TL outermost, CC innermost); s16..s19 are
#   [mu_P^0, mu_P^1, mu_N^0, mu_N^1], matching the reference serialization.
# Parameters: Regazzoni 2020 human body-temperature calibration defaults.
# Ca transient: double-exponential, c0=1e-4 mM, cmax=9e-4 mM,
#               tau_rise=20ms, tau_decay=50ms, onset=10ms.
# SL: raised-cosine ramp, SL0=2.2um, SL_min=2.134um,
#     shortening onset=30ms, minimum at 150ms, recovery complete at 350ms.
# dt=1ms, 600 steps.
# Active tension units: reference uses kPa; values in this file are in MPa
#   (divided by 1000; a_XB=22.894 MPa in svMP corresponds to 22894 kPa in ref).
# Values extracted from the existing test without modification.
step,Ta,s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19
"""

DRIVER_SOURCE = r'''// Narrow driver for the pinned authors' RDQ20-MF model.
#include "model_RDQ20_MF.hpp"

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  if (argc != 4) {
    std::fprintf(stderr,
                 "usage: regazzoni_reference_driver <params.json> <inputs.csv> <raw.csv>\n");
    return 2;
  }

  model_RDQ20_MF model(argv[1]);
  std::ifstream input(argv[2]);
  std::ofstream output(argv[3]);
  if (!input || !output) {
    std::fprintf(stderr, "failed to open an input or output file\n");
    return 3;
  }

  output << std::setprecision(17);
  output << "step,t_start_s,t_end_s,dt_s,Ca_uM,SL_um,dSLdt_um_per_s";
  for (int i = 0; i < 20; ++i)
    output << ",s" << i;
  output << ",Ta_kPa\n";

  std::string line;
  std::getline(input, line);  // header
  std::vector<double> state(20, 0.0);
  state[0] = 1.0;

  while (std::getline(input, line)) {
    if (line.empty())
      continue;
    std::stringstream stream(line);
    std::string cell;
    std::vector<double> values;
    while (std::getline(stream, cell, ','))
      values.push_back(std::stod(cell));
    if (values.size() != 7) {
      std::fprintf(stderr, "expected seven input columns\n");
      return 4;
    }

    const long step = static_cast<long>(values[0]);
    const double t_start_s = values[1];
    const double t_end_s = values[2];
    const double dt_s = values[3];
    const double calcium_uM = values[4];
    const double sarcomere_length_um = values[5];
    const double dSLdt_um_per_s = values[6];

    model.solve_time_step(state, calcium_uM, sarcomere_length_um,
                          dSLdt_um_per_s, dt_s);
    const double tension_kPa =
        model.get_active_tension(state, sarcomere_length_um);

    output << step << ',' << t_start_s << ',' << t_end_s << ',' << dt_s
           << ',' << calcium_uM << ',' << sarcomere_length_um << ','
           << dSLdt_um_per_s;
    for (double value : state)
      output << ',' << value;
    output << ',' << tension_kPa << '\n';
  }
  return 0;
}
'''


def run(command: list[str], *, cwd: Path | None = None) -> subprocess.CompletedProcess[bytes]:
    printable = " ".join(shlex.quote(part) for part in command)
    print(f"+ {printable}")
    return subprocess.run(command, cwd=cwd, check=True)


def git_object(repo: Path, path: str) -> bytes:
    result = subprocess.run(
        ["git", "-C", str(repo), "show", f"{PINNED_COMMIT}:{path}"],
        check=True,
        stdout=subprocess.PIPE,
    )
    return result.stdout


def stage_pinned_sources(repo: Path, source_dir: Path, patch_path: Path) -> Path:
    if not repo.is_dir():
        raise FileNotFoundError(
            f"cardiac-activation checkout not found: {repo}\n"
            "Clone https://github.com/FrancescoRegazzoni/cardiac-activation "
            "or pass --reference-repo."
        )
    subprocess.run(
        ["git", "-C", str(repo), "cat-file", "-e", f"{PINNED_COMMIT}^{{commit}}"],
        check=True,
    )

    source_dir.mkdir(parents=True)
    for repository_path in REFERENCE_SOURCE_PATHS:
        destination = source_dir / Path(repository_path).name
        destination.write_bytes(git_object(repo, repository_path))

    run(["patch", "--silent", "-p1", "-i", str(patch_path)], cwd=source_dir)
    patched = (source_dir / "model_RDQ20_MF.cpp").read_text(encoding="utf-8")
    marker = "XB_A.setZero();"
    if patched.count(marker) != 1:
        raise RuntimeError("the deterministic XB_A patch was not applied exactly once")

    driver_path = source_dir / "regazzoni_reference_driver.cpp"
    driver_path.write_text(DRIVER_SOURCE, encoding="utf-8")
    return driver_path


def include_with_header(
    explicit: Path | None,
    environment_name: str,
    header: Path,
    candidates: Iterable[Path],
) -> Path:
    possible: list[Path] = []
    if explicit is not None:
        possible.append(explicit.expanduser())
    configured = os.environ.get(environment_name)
    if configured:
        possible.append(Path(configured).expanduser())
    possible.extend(candidates)
    for candidate in possible:
        if (candidate / header).is_file():
            return candidate.resolve()
    raise FileNotFoundError(
        f"could not locate {header}; pass the corresponding command-line option "
        f"or set {environment_name}"
    )


def eigen_candidates() -> list[Path]:
    candidates = [
        Path("/opt/homebrew/include/eigen3"),
        Path("/usr/local/include/eigen3"),
        Path("/usr/include/eigen3"),
    ]
    spack_root = Path.home() / "spack" / "opt" / "spack"
    if spack_root.is_dir():
        candidates.extend(
            path.parent.parent
            for path in spack_root.glob("**/include/eigen3/Eigen/Dense")
        )
    return candidates


def boost_candidates() -> list[Path]:
    return [
        Path("/opt/homebrew/include"),
        Path("/usr/local/include"),
        Path("/usr/include"),
    ]


def compile_driver(
    source_dir: Path,
    executable: Path,
    cxx: str,
    eigen_include: Path,
    boost_include: Path,
) -> None:
    compiler = shutil.which(cxx)
    if compiler is None:
        raise FileNotFoundError(f"C++ compiler not found: {cxx}")
    run(
        [
            compiler,
            "-std=c++17",
            "-O2",
            "-I",
            str(eigen_include),
            "-I",
            str(boost_include),
            "-I",
            str(source_dir),
            str(source_dir / "regazzoni_reference_driver.cpp"),
            str(source_dir / "model_RDQ20_MF.cpp"),
            str(source_dir / "sarcomere.cpp"),
            "-o",
            str(executable),
        ]
    )


def calcium_at(t_ms: float) -> float:
    """Calcium in mM at the left endpoint of an outer step."""
    if t_ms < CALCIUM_ONSET_MS:
        return BASELINE_CALCIUM_MM
    peak_time = (
        math.log(CALCIUM_RISE_MS / CALCIUM_DECAY_MS)
        * CALCIUM_RISE_MS
        * CALCIUM_DECAY_MS
        / (CALCIUM_RISE_MS - CALCIUM_DECAY_MS)
    )
    peak_raw = math.exp(-peak_time / CALCIUM_DECAY_MS) - math.exp(
        -peak_time / CALCIUM_RISE_MS
    )
    elapsed = t_ms - CALCIUM_ONSET_MS
    raw = math.exp(-elapsed / CALCIUM_DECAY_MS) - math.exp(
        -elapsed / CALCIUM_RISE_MS
    )
    return BASELINE_CALCIUM_MM + (
        PEAK_CALCIUM_MM - BASELINE_CALCIUM_MM
    ) * raw / peak_raw


def sarcomere_length_at(t_ms: float) -> float:
    """Raised-cosine sarcomere length in micrometers."""
    if SHORTENING_ONSET_MS <= t_ms <= MINIMUM_TIME_MS:
        fraction = (t_ms - SHORTENING_ONSET_MS) / (
            MINIMUM_TIME_MS - SHORTENING_ONSET_MS
        )
        return RESTING_SL_UM - (RESTING_SL_UM - MINIMUM_SL_UM) * (
            1.0 - math.cos(math.pi * fraction)
        ) / 2.0
    if MINIMUM_TIME_MS < t_ms <= RECOVERY_TIME_MS:
        fraction = (t_ms - MINIMUM_TIME_MS) / (
            RECOVERY_TIME_MS - MINIMUM_TIME_MS
        )
        return MINIMUM_SL_UM + (RESTING_SL_UM - MINIMUM_SL_UM) * (
            1.0 - math.cos(math.pi * fraction)
        ) / 2.0
    return RESTING_SL_UM


def protocol_rows() -> list[tuple[int, float, float, float, float, float, float]]:
    rows = []
    for step in range(UPDATE_COUNT):
        t_start_ms = step * DT_MS
        t_end_ms = t_start_ms + DT_MS
        sl_start = sarcomere_length_at(t_start_ms)
        sl_end = sarcomere_length_at(t_end_ms)
        dSLdt_um_per_s = (sl_end - sl_start) / DT_S
        rows.append(
            (
                step,
                t_start_ms / 1000.0,
                t_end_ms / 1000.0,
                DT_S,
                calcium_at(t_start_ms) * 1000.0,
                sl_start,
                dSLdt_um_per_s,
            )
        )
    return rows


def write_inputs(path: Path, rows: list[tuple[float, ...]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(
            [
                "step",
                "t_start_s",
                "t_end_s",
                "dt_s",
                "Ca_uM",
                "SL_um",
                "dSLdt_um_per_s",
            ]
        )
        for row in rows:
            writer.writerow([row[0], *(format(value, ".17g") for value in row[1:])])


def load_raw_trajectory(path: Path) -> list[dict[str, float]]:
    expected_header = [
        "step",
        "t_start_s",
        "t_end_s",
        "dt_s",
        "Ca_uM",
        "SL_um",
        "dSLdt_um_per_s",
        *(f"s{i}" for i in range(STATE_COUNT)),
        "Ta_kPa",
    ]
    rows: list[dict[str, float]] = []
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != expected_header:
            raise ValueError(f"unexpected raw trajectory columns: {reader.fieldnames}")
        for raw in reader:
            rows.append({key: float(value) for key, value in raw.items()})

    if len(rows) != UPDATE_COUNT:
        raise ValueError(f"expected {UPDATE_COUNT} trajectory rows, found {len(rows)}")
    for expected_step, row in enumerate(rows):
        if row["step"] != expected_step:
            raise ValueError("raw trajectory steps are not consecutive")
        values = list(row.values())
        if not all(math.isfinite(value) for value in values):
            raise ValueError(f"non-finite value at step {expected_step}")
    return rows


def compatible_csv(rows: list[dict[str, float]]) -> bytes:
    by_step = {int(row["step"]): row for row in rows}
    if tuple(step for step in CHECKPOINT_STEPS if step in by_step) != CHECKPOINT_STEPS:
        raise ValueError("one or more required checkpoints are missing")

    lines = [COMPATIBLE_PREAMBLE]
    for step in CHECKPOINT_STEPS:
        row = by_step[step]
        values = [row["Ta_kPa"] / 1000.0]
        values.extend(row[f"s{i}"] for i in range(STATE_COUNT))
        lines.append(f"{step}," + ",".join(format(value, ".16e") for value in values) + "\n")
    return "".join(lines).encode("utf-8")


def atomic_write(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_bytes(content)
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Regenerate the patched authors-C++ Regazzoni oracle."
    )
    parser.add_argument(
        "--reference-repo",
        type=Path,
        required=True,
        help="cardiac-activation git checkout containing the pinned commit",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="svMultiPhysics-compatible checkpoint CSV",
    )
    parser.add_argument("--cxx", default=os.environ.get("CXX", "c++"))
    parser.add_argument("--eigen-include", type=Path)
    parser.add_argument("--boost-include", type=Path)
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    patch_path = script_dir / "zero_initialize_xb_matrix.patch"
    if not patch_path.is_file():
        raise FileNotFoundError(f"required patch not found: {patch_path}")
    eigen_include = include_with_header(
        args.eigen_include,
        "EIGEN3_INCLUDE_DIR",
        Path("Eigen/Dense"),
        eigen_candidates(),
    )
    boost_include = include_with_header(
        args.boost_include,
        "BOOST_INCLUDE_DIR",
        Path("boost/property_tree/json_parser.hpp"),
        boost_candidates(),
    )

    with tempfile.TemporaryDirectory(prefix="regazzoni-oracle-") as temporary:
        work = Path(temporary)
        source_dir = work / "source"
        driver_source = stage_pinned_sources(
            args.reference_repo.resolve(), source_dir, patch_path
        )
        if not driver_source.is_file():
            raise RuntimeError("reference driver source was not staged")
        executable = work / "regazzoni_reference_driver"
        compile_driver(source_dir, executable, args.cxx, eigen_include, boost_include)

        inputs = work / "inputs.csv"
        raw_output = work / "raw_trajectory.csv"
        inputs_rows = protocol_rows()
        write_inputs(inputs, inputs_rows)
        run(
            [
                str(executable),
                str(source_dir / "params_RDQ20-MF_human_body-temperature.json"),
                str(inputs),
                str(raw_output),
            ]
        )
        trajectory = load_raw_trajectory(raw_output)
        computed_content = compatible_csv(trajectory)
        atomic_write(args.output.resolve(), computed_content)

    print(f"wrote compatible checkpoint CSV: {args.output.resolve()}")
    print("PASS: generated a fresh authors-C++ compatible trajectory")


if __name__ == "__main__":
    main()
