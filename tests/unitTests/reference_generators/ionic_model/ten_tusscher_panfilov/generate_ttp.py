#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause

"""Generate independent TP06 EPI/ENDO/M FE/Rush--Larsen trajectory oracles."""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
from typing import Iterator, TextIO

from ttp_model import (
    Gates,
    MainState,
    Parameters,
    default_endo_gates,
    default_endo_state,
    default_epi_gates,
    default_epi_state,
    default_m_gates,
    default_m_state,
    endo_parameters,
    hybrid_forward_euler_step,
    named_values,
    m_parameters,
    stimulus_current,
)


DT_MS = 0.005
UPDATE_COUNT = 120_000
EPI_CHECKPOINTS = (
    0,
    2_000,
    2_001,
    2_200,
    2_263,
    3_193,
    3_833,
    5_672,
    6_965,
    7_217,
    12_456,
    15_090,
    20_000,
    50_852,
    54_710,
    60_472,
    80_000,
    104_054,
    120_000,
)
ENDO_CHECKPOINTS = (
    0,
    2_000,
    2_001,
    2_200,
    2_374,
    3_520,
    6_003,
    7_354,
    7_622,
    20_000,
    20_464,
    37_064,
    49_329,
    53_028,
    53_338,
    58_781,
    60_000,
    80_000,
    106_073,
    120_000,
)
M_CHECKPOINTS = (
    0,
    2_000,
    2_001,
    2_200,
    2_269,
    3_144,
    5_755,
    7_395,
    7_622,
    16_162,
    20_000,
    46_611,
    59_665,
    66_376,
    66_838,
    72_997,
    80_000,
    100_000,
    114_081,
    120_000,
)
STATE_COLUMNS = (
    "V_mV",
    "Ki",
    "Nai",
    "Cai",
    "Ca_ss",
    "Ca_SR",
    "R_prime",
    "Xr1",
    "Xr2",
    "Xs",
    "m",
    "h",
    "j",
    "d",
    "f",
    "f2",
    "fCass",
    "s",
    "r",
)


@dataclass(frozen=True)
class CanonicalTrajectory:
    rows: tuple[tuple[int, tuple[float, ...]], ...]


@dataclass(frozen=True)
class ProfileConfiguration:
    name: str
    zone_id: int
    parameters: Parameters
    initial_state: MainState
    initial_gates: Gates
    checkpoints: tuple[int, ...]


PROFILE_CONFIGURATIONS = {
    "epi": ProfileConfiguration(
        "EPI", 1, Parameters(), default_epi_state(), default_epi_gates(), EPI_CHECKPOINTS
    ),
    "endo": ProfileConfiguration(
        "ENDO",
        2,
        endo_parameters(),
        default_endo_state(),
        default_endo_gates(),
        ENDO_CHECKPOINTS,
    ),
    "m": ProfileConfiguration(
        "M",
        3,
        m_parameters(),
        default_m_state(),
        default_m_gates(),
        M_CHECKPOINTS,
    ),
}


def iter_full_trajectory(
    profile: str = "epi",
) -> Iterator[tuple[int, tuple[float, ...]]]:
    """Yield every completed-update state from the canonical TP06 protocol."""

    configuration = PROFILE_CONFIGURATIONS[profile]
    state = configuration.initial_state
    gates = configuration.initial_gates

    for step in range(UPDATE_COUNT + 1):
        values = tuple(value for _, value in named_values(state, gates))
        if len(values) != len(STATE_COLUMNS):
            raise AssertionError(f"expected 19 states, found {len(values)}")
        if not all(isfinite(value) for value in values):
            raise AssertionError(f"non-finite state after update {step}")
        yield step, values

        if step == UPDATE_COUNT:
            continue
        old_time_ms = step * DT_MS
        result = hybrid_forward_euler_step(
            state,
            gates,
            dt=DT_MS,
            parameters=configuration.parameters,
            phenotype=profile,
            stimulus_current=stimulus_current(old_time_ms),
            sac_coefficient=0.0,
        )
        state, gates = result.state, result.gates


def generate_canonical_trajectory(
    profile: str = "epi",
) -> CanonicalTrajectory:
    """Run one frozen 600 ms protocol and collect canonical checkpoints."""

    configuration = PROFILE_CONFIGURATIONS[profile]
    checkpoints = configuration.checkpoints

    if tuple(sorted(set(checkpoints))) != checkpoints:
        raise AssertionError("checkpoints must be unique and strictly increasing")
    if checkpoints[0] != 0 or checkpoints[-1] != UPDATE_COUNT:
        raise AssertionError("checkpoints must include the initial and final states")

    checkpoint_set = set(checkpoints)
    rows: list[tuple[int, tuple[float, ...]]] = []
    stimulated_update_count = 0

    for step, values in iter_full_trajectory(profile):
        if step > 0:
            old_time_ms = (step - 1) * DT_MS
            stimulated_update_count += stimulus_current(old_time_ms) != 0.0
        if step in checkpoint_set:
            rows.append((step, values))

    if len(rows) != len(checkpoints):
        raise AssertionError("not every canonical checkpoint was recorded")
    if tuple(step for step, _ in rows) != checkpoints:
        raise AssertionError("recorded checkpoint ordering changed")
    if stimulated_update_count != 200:
        raise AssertionError(
            f"expected 200 stimulated updates, found {stimulated_update_count}"
        )
    return CanonicalTrajectory(tuple(rows))


def write_csv(
    trajectory: CanonicalTrajectory, profile: str, stream: TextIO
) -> None:
    configuration = PROFILE_CONFIGURATIONS[profile]
    stream.write(
        f"# Provenance: independent TP06 {configuration.name} "
        "FE/Rush-Larsen implementation; not generated by svMP.\n"
    )
    stream.write(
        f"# Protocol: zone_id={configuration.zone_id}, dt={DT_MS:g} ms, "
        f"{UPDATE_COUNT} updates, Ksac=0, no pre-pacing.\n"
    )
    stream.write("# Stimulus: Istim=-52 pA/pF for 10 <= t < 11 ms; zero otherwise.\n")
    stream.write(
        f"# Initial state: published/CellML TP06 {configuration.name} defaults.\n"
    )
    stream.write("# Checkpoint N is the state after exactly N completed updates.\n")
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(("step", *STATE_COLUMNS))
    for step, values in trajectory.rows:
        writer.writerow((step, *values))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, help="write canonical CSV to this path instead of stdout"
    )
    parser.add_argument(
        "--profile", choices=tuple(PROFILE_CONFIGURATIONS), default="epi"
    )
    args = parser.parse_args()

    trajectory = generate_canonical_trajectory(args.profile)
    if args.output is None:
        write_csv(trajectory, args.profile, sys.stdout)
    else:
        with args.output.open("w", encoding="utf-8", newline="") as stream:
            write_csv(trajectory, args.profile, stream)


if __name__ == "__main__":
    main()
