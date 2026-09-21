#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause

"""Generate the canonical Nash--Panfilov active-tension reference."""

from __future__ import annotations

import argparse
import csv
import io
import math
from pathlib import Path
import sys
from typing import Iterable

import nash_panfilov_2004 as model

# Canonical calcium protocol.
C0 = 1.0e-4
CMAX = 9.0e-4
TAU_RISE_MS = 20.0
TAU_DECAY_MS = 50.0
ONSET_MS = 10.0
BASE_DT_MS = 1.0
BASE_STEPS = 200
CHECKPOINT_STEPS = (0, 10, 30, 60, 99, 149, 199)


def peak_raw_calcium_factor() -> float:
    peak_time = (
        math.log(TAU_RISE_MS / TAU_DECAY_MS)
        * TAU_RISE_MS
        * TAU_DECAY_MS
        / (TAU_RISE_MS - TAU_DECAY_MS)
    )
    return math.exp(-peak_time / TAU_DECAY_MS) - math.exp(
        -peak_time / TAU_RISE_MS
    )


PEAK_RAW_CA_FACTOR = peak_raw_calcium_factor()


def calcium_at(t_ms: float) -> float:
    """Raised double-exponential calcium transient used by the GTest."""
    if t_ms < ONSET_MS:
        return C0
    shifted_time = t_ms - ONSET_MS
    raw = math.exp(-shifted_time / TAU_DECAY_MS) - math.exp(
        -shifted_time / TAU_RISE_MS
    )
    return C0 + (CMAX - C0) * raw / PEAK_RAW_CA_FACTOR


def fe_solution(dt_ms: float, sample_times_ms: Iterable[float]) -> dict[float, float]:
    """Forward Euler with left-endpoint calcium, matching the GTest protocol."""
    sample_times = tuple(sample_times_ms)
    sample_indices = {round(time / dt_ms): time for time in sample_times}
    for index, time in sample_indices.items():
        if not math.isclose(index * dt_ms, time, rel_tol=0.0, abs_tol=1.0e-12):
            raise ValueError(f"sample time {time} is not aligned with dt={dt_ms}")

    states = model.create_states_array()
    rates = [math.nan] * model.STATE_COUNT
    variables = model.create_variables_array()
    states[0] = 0.0

    def prescribed_calcium(voi, states, rates, variables, variable_index):
        return calcium_at(voi)

    model.initialise_variables(
        0.0, states, rates, variables, prescribed_calcium
    )

    values: dict[float, float] = {}
    final_index = max(sample_indices)
    for index in range(final_index):
        t_ms = index * dt_ms
        model.compute_rates(
            t_ms, states, rates, variables, prescribed_calcium
        )
        states[0] += dt_ms * rates[0]
        if not math.isfinite(states[0]):
            raise FloatingPointError(
                f"non-finite Nash-Panfilov state at step {index + 1}"
            )
        completed_steps = index + 1
        if completed_steps in sample_indices:
            values[sample_indices[completed_steps]] = states[0]
    return values


def canonical_csv() -> str:
    """Return the deterministic CSV consumed by the ActiveStress unit test.

    ActiveStress checkpoint label N denotes the state after outer update N,
    whose left-endpoint time is N*dt. Thus label 0 is the state after the first
    update, unlike the IonicModel convention where checkpoint 0 is initial.
    """

    sample_times = tuple((step + 1) * BASE_DT_MS for step in CHECKPOINT_STEPS)
    values = fe_solution(BASE_DT_MS, sample_times)
    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(("step", "Ta"))
    for step, time in zip(CHECKPOINT_STEPS, sample_times, strict=True):
        writer.writerow((step, f"{values[time]:.16e}"))
    return output.getvalue()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, help="write canonical CSV instead of stdout"
    )
    args = parser.parse_args()

    text = canonical_csv()
    if args.output is None:
        sys.stdout.write(text)
    else:
        args.output.write_text(text, encoding="utf-8")



if __name__ == "__main__":
    main()
