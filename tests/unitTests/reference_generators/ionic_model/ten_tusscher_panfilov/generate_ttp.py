#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause

"""Generate canonical TP06 EPI/ENDO/M FE/Rush--Larsen trajectory CSVs."""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from math import exp, isfinite
from pathlib import Path
from types import ModuleType
from typing import Iterator, TextIO

import ten_tusscher_2006_endo as endo_model
import ten_tusscher_2006_epi as epi_model
import ten_tusscher_2006_m as m_model


DT_MS = 0.005
UPDATE_COUNT = 120_000
STIMULUS_AMPLITUDE = -52.0
STIMULUS_START_MS = 10.0
STIMULUS_END_MS = 11.0
EPI_CHECKPOINTS = (
    0, 2_000, 2_001, 2_200, 2_263, 3_193, 3_833, 5_672, 6_965, 7_217,
    12_456, 15_090, 20_000, 50_852, 54_710, 60_472, 80_000, 104_054,
    120_000,
)
ENDO_CHECKPOINTS = (
    0, 2_000, 2_001, 2_200, 2_374, 3_520, 6_003, 7_354, 7_622, 20_000,
    20_464, 37_064, 49_329, 53_028, 53_338, 58_781, 60_000, 80_000,
    106_073, 120_000,
)
M_CHECKPOINTS = (
    0, 2_000, 2_001, 2_200, 2_269, 3_144, 5_755, 7_395, 7_622, 16_162,
    20_000, 46_611, 59_665, 66_376, 66_838, 72_997, 80_000, 100_000,
    114_081, 120_000,
)
STATE_COLUMNS = (
    "V_mV", "Ki", "Nai", "Cai", "Ca_ss", "Ca_SR", "R_prime", "Xr1",
    "Xr2", "Xs", "m", "h", "j", "d", "f", "f2", "fCass", "s", "r",
)

# Generated state order is V, Ki, Nai, Cai, Xr1, Xr2, Xs, m, h, j, Ca_ss,
# d, f, f2, fCass, s, r, Ca_SR, R_prime. The CSV stores svMP's seven X
# states followed by its twelve Xg states.
OUTPUT_STATE_INDICES = (0, 1, 2, 3, 10, 17, 18, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16)
MAIN_STATE_INDICES = (0, 1, 2, 3, 10, 17, 18)
# (generated state index, steady-state algebraic index, time-constant index)
GATE_KINETICS = (
    (4, 0, 34), (5, 1, 35), (6, 2, 36), (7, 3, 37), (8, 4, 38),
    (9, 5, 39), (11, 6, 42), (12, 7, 20), (13, 8, 21), (14, 9, 22),
    (15, 10, 23), (16, 11, 24),
)


@dataclass(frozen=True)
class CanonicalTrajectory:
    rows: tuple[tuple[int, tuple[float, ...]], ...]


@dataclass(frozen=True)
class ProfileConfiguration:
    model: ModuleType
    checkpoints: tuple[int, ...]


PROFILE_CONFIGURATIONS = {
    "epi": ProfileConfiguration(epi_model, EPI_CHECKPOINTS),
    "endo": ProfileConfiguration(endo_model, ENDO_CHECKPOINTS),
    "m": ProfileConfiguration(m_model, M_CHECKPOINTS),
}


def stimulus_current(time_ms: float) -> float:
    if not isfinite(time_ms):
        raise ValueError("stimulus time must be finite")
    if STIMULUS_START_MS <= time_ms < STIMULUS_END_MS:
        return STIMULUS_AMPLITUDE
    return 0.0


def iter_full_trajectory(profile: str = "epi") -> Iterator[tuple[int, tuple[float, ...]]]:
    """Yield every completed-update state from the canonical TP06 protocol."""

    configuration = PROFILE_CONFIGURATIONS[profile]
    states, constants = configuration.model.initConsts()
    if len(states) != 19:
        raise AssertionError(f"expected 19 states, found {len(states)}")

    for step in range(UPDATE_COUNT + 1):
        values = tuple(states[index] for index in OUTPUT_STATE_INDICES)
        if not all(isfinite(value) for value in values):
            raise AssertionError(f"non-finite state after update {step}")
        yield step, values

        if step == UPDATE_COUNT:
            continue

        old_time_ms = step * DT_MS
        rates, algebraic = configuration.model.computeRates(
            old_time_ms, states, constants, stimulus_current(old_time_ms)
        )
        updated_states = list(states)
        for state_index in MAIN_STATE_INDICES:
            updated_states[state_index] = states[state_index] + DT_MS * rates[state_index]
        for state_index, steady_state_index, time_constant_index in GATE_KINETICS:
            steady_state = algebraic[steady_state_index]
            time_constant = algebraic[time_constant_index]
            updated_states[state_index] = steady_state - (
                steady_state - states[state_index]
            ) * exp(-DT_MS / time_constant)
        states = updated_states


def generate_canonical_trajectory(profile: str = "epi") -> CanonicalTrajectory:
    """Run one frozen 600 ms protocol and collect canonical checkpoints."""

    checkpoints = PROFILE_CONFIGURATIONS[profile].checkpoints
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

    if tuple(step for step, _ in rows) != checkpoints:
        raise AssertionError("not every canonical checkpoint was recorded")
    if stimulated_update_count != 200:
        raise AssertionError(
            f"expected 200 stimulated updates, found {stimulated_update_count}"
        )
    return CanonicalTrajectory(tuple(rows))


def write_csv(trajectory: CanonicalTrajectory, stream: TextIO) -> None:
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
        write_csv(trajectory, sys.stdout)
    else:
        with args.output.open("w", encoding="utf-8", newline="") as stream:
            write_csv(trajectory, stream)


if __name__ == "__main__":
    main()
