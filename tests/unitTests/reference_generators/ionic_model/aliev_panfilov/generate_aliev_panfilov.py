#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause

"""Generate the canonical stimulated Aliev--Panfilov reference."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path
from typing import Iterable

from aliev_panfilov_1996 import (
    compute_computed_constants,
    compute_rates,
    create_states_array,
    create_variables_array,
    initialise_variables,
)


BASE_DT_MS = 0.1
BASE_STEPS = 6000
INITIAL_V_MV = -80.0
INITIAL_W = 0.001
STIMULUS_START_MS = 10.0
STIMULUS_END_MS = 12.0
STIMULUS_AMPLITUDE = -35.714
CHECKPOINT_STEPS = (
    0, 100, 101, 118, 120, 257, 1000, 3000, 3782, 3874, 4082, 4500, 6000
)


def public_stimulus_at_time(time_ms: float) -> float:
    """Return public stimulus; negative current is depolarizing for AP."""
    if STIMULUS_START_MS <= time_ms < STIMULUS_END_MS:
        return STIMULUS_AMPLITUDE
    return 0.0


def cellml_stimulus(
    time_ms: float,
    states: list[float],
    rates: list[float],
    variables: list[float],
    variable_index: int,
) -> float:
    """Map svMP public stimulus to the CellML current convention."""
    return -variables[0] * public_stimulus_at_time(time_ms)


def fe_solution(
    dt_ms: float,
    sample_times_ms: Iterable[float],
    final_time_ms: float | None = None,
) -> dict[float, tuple[float, float]]:
    sample_times = tuple(sample_times_ms)
    sample_indices = {round(time / dt_ms): time for time in sample_times}
    for index, time in sample_indices.items():
        if not math.isclose(index * dt_ms, time, rel_tol=0.0, abs_tol=1.0e-12):
            raise ValueError(f"sample time {time} is not aligned with dt={dt_ms}")

    states = create_states_array()
    rates = create_states_array()
    variables = create_variables_array()
    initialise_variables(0.0, states, rates, variables, cellml_stimulus)
    states[0] = INITIAL_V_MV
    states[1] = INITIAL_W
    compute_computed_constants(variables)

    values: dict[float, tuple[float, float]] = {}
    if 0 in sample_indices:
        values[sample_indices[0]] = (states[0], states[1])
    final_index = max(sample_indices)
    if final_time_ms is not None:
        requested_final_index = round(final_time_ms / dt_ms)
        if not math.isclose(
            requested_final_index * dt_ms,
            final_time_ms,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ):
            raise ValueError(f"final time {final_time_ms} is not aligned with dt={dt_ms}")
        final_index = max(final_index, requested_final_index)

    for index in range(final_index):
        compute_rates(index * dt_ms, states, rates, variables, cellml_stimulus)
        states[0] += dt_ms * rates[0]
        states[1] += dt_ms * rates[1]
        if not all(math.isfinite(value) for value in states):
            raise FloatingPointError(f"non-finite AP state at step {index + 1}")
        completed_steps = index + 1
        if completed_steps in sample_indices:
            values[sample_indices[completed_steps]] = (states[0], states[1])
    return values


def write_oracle(output: Path, samples: dict[float, tuple[float, float]]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("step", "V_mV", "w"))
        for step in CHECKPOINT_STEPS:
            time_ms = step * BASE_DT_MS
            voltage, w = samples[time_ms]
            writer.writerow((step, f"{voltage:.16e}", f"{w:.16e}"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    sample_times = tuple(step * BASE_DT_MS for step in CHECKPOINT_STEPS)
    samples = fe_solution(BASE_DT_MS, sample_times, BASE_STEPS * BASE_DT_MS)
    write_oracle(args.output, samples)
    digest = hashlib.sha256(args.output.read_bytes()).hexdigest()
    print(f"wrote {args.output}")
    print(f"sha256={digest}")
    print("checkpoints=" + ",".join(map(str, CHECKPOINT_STEPS)))


if __name__ == "__main__":
    main()
