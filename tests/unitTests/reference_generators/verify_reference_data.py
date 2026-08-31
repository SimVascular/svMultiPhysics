#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause

"""Regenerate and verify all unit-test reference trajectories."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import math
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parent
REGAZZONI_ABS_TOLERANCE = 5.0e-15
REGAZZONI_REL_TOLERANCE = 5.0e-13


@dataclass(frozen=True)
class ReferenceSpecification:
    label: str
    generator: Path
    filename: str
    arguments: tuple[str, ...] = ()


PURE_PYTHON_REFERENCES = (
    ReferenceSpecification(
        "Nash-Panfilov",
        ROOT / "active_stress/nash_panfilov/generate_nash_panfilov.py",
        "active_stress_nash_panfilov_twitch.csv",
    ),
    ReferenceSpecification(
        "Aliev-Panfilov",
        ROOT / "ionic_model/aliev_panfilov/generate_aliev_panfilov.py",
        "ionic_aliev_panfilov_stimulated_trajectory.csv",
    ),
    ReferenceSpecification(
        "FitzHugh-Nagumo",
        ROOT / "ionic_model/fitzhugh_nagumo/generate_fitzhugh_nagumo.py",
        "ionic_fitzhugh_nagumo_fe_trajectory.csv",
    ),
    ReferenceSpecification(
        "Bueno-Orovio EPI",
        ROOT / "ionic_model/bueno_orovio/generate_bueno_orovio.py",
        "ionic_bueno_orovio_epi_trajectory.csv",
        ("--profile", "epi"),
    ),
    ReferenceSpecification(
        "Bueno-Orovio ENDO",
        ROOT / "ionic_model/bueno_orovio/generate_bueno_orovio.py",
        "ionic_bueno_orovio_endo_trajectory.csv",
        ("--profile", "endo"),
    ),
    ReferenceSpecification(
        "Bueno-Orovio M",
        ROOT / "ionic_model/bueno_orovio/generate_bueno_orovio.py",
        "ionic_bueno_orovio_m_trajectory.csv",
        ("--profile", "m"),
    ),
    ReferenceSpecification(
        "TP06 EPI",
        ROOT / "ionic_model/ten_tusscher_panfilov/generate_ttp.py",
        "ionic_ttp_epi_trajectory.csv",
        ("--profile", "epi"),
    ),
    ReferenceSpecification(
        "TP06 ENDO",
        ROOT / "ionic_model/ten_tusscher_panfilov/generate_ttp.py",
        "ionic_ttp_endo_trajectory.csv",
        ("--profile", "endo"),
    ),
    ReferenceSpecification(
        "TP06 M",
        ROOT / "ionic_model/ten_tusscher_panfilov/generate_ttp.py",
        "ionic_ttp_m_trajectory.csv",
        ("--profile", "m"),
    ),
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_generator(
    specification: ReferenceSpecification,
    output: Path,
    extra_arguments: tuple[str, ...] = (),
) -> None:
    command = [
        sys.executable,
        str(specification.generator),
        *specification.arguments,
        *extra_arguments,
        "--output",
        str(output),
    ]
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != 0:
        detail = result.stderr.strip() or result.stdout.strip()
        raise RuntimeError(f"generator failed: {' '.join(command)}\n{detail}")
    if not output.is_file():
        raise RuntimeError(f"generator did not create {output}")


def generate_twice(
    specification: ReferenceSpecification,
    directory: Path,
    extra_arguments: tuple[str, ...] = (),
) -> Path:
    first = directory / f"first-{specification.filename}"
    second = directory / f"second-{specification.filename}"
    run_generator(specification, first, extra_arguments)
    run_generator(specification, second, extra_arguments)
    if first.read_bytes() != second.read_bytes():
        raise RuntimeError(
            f"nondeterministic output: {sha256(first)} != {sha256(second)}"
        )
    return first


def compare_bytes(generated: Path, reference: Path) -> None:
    if generated.read_bytes() != reference.read_bytes():
        raise RuntimeError(
            "byte mismatch: "
            f"generated={sha256(generated)} reference={sha256(reference)}"
        )


def csv_structure(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as stream:
        lines = [line for line in stream if line.strip() and not line.startswith("#")]
    reader = csv.DictReader(lines)
    if reader.fieldnames is None:
        raise RuntimeError(f"CSV has no header: {path}")
    return tuple(reader.fieldnames), list(reader)


def compare_regazzoni(generated: Path, reference: Path) -> tuple[float, float]:
    generated_header, generated_rows = csv_structure(generated)
    reference_header, reference_rows = csv_structure(reference)
    if generated_header != reference_header:
        raise RuntimeError(
            f"Regazzoni schema mismatch: {generated_header} != {reference_header}"
        )
    if len(generated_rows) != len(reference_rows):
        raise RuntimeError("Regazzoni checkpoint counts differ")

    maximum_absolute = 0.0
    maximum_scaled = 0.0
    for generated_row, reference_row in zip(
        generated_rows, reference_rows, strict=True
    ):
        if generated_row["step"] != reference_row["step"]:
            raise RuntimeError("Regazzoni checkpoint ordering differs")
        for column in generated_header[1:]:
            actual = float(generated_row[column])
            expected = float(reference_row[column])
            absolute = abs(actual - expected)
            scaled = absolute / max(1.0, abs(expected))
            maximum_absolute = max(maximum_absolute, absolute)
            maximum_scaled = max(maximum_scaled, scaled)
            if not math.isclose(
                actual,
                expected,
                rel_tol=REGAZZONI_REL_TOLERANCE,
                abs_tol=REGAZZONI_ABS_TOLERANCE,
            ):
                raise RuntimeError(
                    f"Regazzoni mismatch at step {generated_row['step']}, {column}: "
                    f"generated={actual:.17g}, reference={expected:.17g}"
                )
    return maximum_absolute, maximum_scaled


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--regazzoni-reference-repo", type=Path)
    parser.add_argument("--skip-regazzoni", action="store_true")
    parser.add_argument("--regazzoni-eigen-include", type=Path)
    parser.add_argument("--regazzoni-boost-include", type=Path)
    parser.add_argument("--cxx")
    args = parser.parse_args()
    if not args.skip_regazzoni and args.regazzoni_reference_repo is None:
        parser.error(
            "--regazzoni-reference-repo is required unless --skip-regazzoni is used"
        )
    return args


def main() -> int:
    args = parse_args()
    reference_directory = args.repo.resolve() / "tests/unitTests/reference_data"
    if not reference_directory.is_dir():
        raise SystemExit(f"reference_data directory not found: {reference_directory}")

    failures = 0
    with tempfile.TemporaryDirectory(prefix="svmp-reference-verify-") as temporary:
        work = Path(temporary)
        for specification in PURE_PYTHON_REFERENCES:
            try:
                generated = generate_twice(specification, work)
                reference = reference_directory / specification.filename
                compare_bytes(generated, reference)
                print(f"{specification.label} PASS sha256={sha256(generated)}")
            except Exception as error:  # noqa: BLE001 - report every reference
                failures += 1
                print(f"{specification.label} FAIL {error}")

        if args.skip_regazzoni:
            print("Regazzoni SKIP")
        else:
            specification = ReferenceSpecification(
                "Regazzoni",
                ROOT / "active_stress/regazzoni/generate_regazzoni.py",
                "active_stress_regazzoni_twitch.csv",
            )
            extra = (
                "--reference-repo",
                str(args.regazzoni_reference_repo.resolve()),
            )
            if args.regazzoni_eigen_include is not None:
                extra += ("--eigen-include", str(args.regazzoni_eigen_include))
            if args.regazzoni_boost_include is not None:
                extra += ("--boost-include", str(args.regazzoni_boost_include))
            if args.cxx is not None:
                extra += ("--cxx", args.cxx)
            try:
                generated = generate_twice(specification, work, extra)
                maximum_absolute, maximum_scaled = compare_regazzoni(
                    generated, reference_directory / specification.filename
                )
                print(
                    "Regazzoni PASS "
                    f"sha256={sha256(generated)} "
                    f"max_abs={maximum_absolute:.3e} "
                    f"max_scaled={maximum_scaled:.3e}"
                )
            except Exception as error:  # noqa: BLE001 - concise aggregate report
                failures += 1
                print(f"Regazzoni FAIL {error}")

    print(f"summary: {len(PURE_PYTHON_REFERENCES) + (not args.skip_regazzoni) - failures} PASS, {failures} FAIL")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
