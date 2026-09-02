#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause

"""Regenerate and verify the in-repository unit-test reference generators."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parent
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
) -> None:
    command = [
        sys.executable,
        str(specification.generator),
        *specification.arguments,
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
) -> Path:
    first = directory / f"first-{specification.filename}"
    second = directory / f"second-{specification.filename}"
    run_generator(specification, first)
    run_generator(specification, second)
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    return parser.parse_args()


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

        print(
            "Regazzoni NOT CHECKED -- no in-repository generator; see "
            "active_stress/regazzoni/README.md"
        )

    print(f"summary: {len(PURE_PYTHON_REFERENCES) - failures} PASS, {failures} FAIL")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
