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


ALPHA = 0.01
GAMMA = 0.002  # svMultiPhysics parameter name: a
B = 0.15
C = 8.0
MU1 = 0.2
MU2 = 0.3

VOLTAGE_SCALE_MV = 100.0
VOLTAGE_OFFSET_MV = -80.0
TIME_SCALE_MS = 12.90
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


def ap_rhs(
    u: float, w: float, stimulus: float = 0.0, stretch_current: float = 0.0
) -> tuple[float, float]:
    """Göktepe--Kuhl split AP equations with svMP current-sign convention."""
    du = C * u * (u - ALPHA) * (1.0 - u) - u * w
    du += -stimulus + stretch_current
    dw = (GAMMA + MU1 * w / (MU2 + u)) * (
        -w - C * u * (u - B - 1.0)
    )
    return du, dw


def voltage_from_u(u: float) -> float:
    return VOLTAGE_SCALE_MV * u + VOLTAGE_OFFSET_MV


def public_stimulus_at_time(time_ms: float) -> float:
    """Return public stimulus; negative current is depolarizing for AP."""
    if STIMULUS_START_MS <= time_ms < STIMULUS_END_MS:
        return STIMULUS_AMPLITUDE
    return 0.0


def internal_stimulus_at_time(time_ms: float) -> float:
    return public_stimulus_at_time(time_ms) * TIME_SCALE_MS / VOLTAGE_SCALE_MV


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

    u = (INITIAL_V_MV - VOLTAGE_OFFSET_MV) / VOLTAGE_SCALE_MV
    w = INITIAL_W
    values: dict[float, tuple[float, float]] = {}
    if 0 in sample_indices:
        values[sample_indices[0]] = (u, w)
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
    dt_tau = dt_ms / TIME_SCALE_MS
    for index in range(final_index):
        stimulus = internal_stimulus_at_time(index * dt_ms)
        du, dw = ap_rhs(u, w, stimulus=stimulus)
        u += dt_tau * du
        w += dt_tau * dw
        if not (math.isfinite(u) and math.isfinite(w)):
            raise FloatingPointError(f"non-finite AP state at step {index + 1}")
        completed_steps = index + 1
        if completed_steps in sample_indices:
            values[sample_indices[completed_steps]] = (u, w)
    return values


def write_oracle(output: Path, samples: dict[float, tuple[float, float]]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("step", "V_mV", "w"))
        for step in CHECKPOINT_STEPS:
            time_ms = step * BASE_DT_MS
            u, w = samples[time_ms]
            writer.writerow((step, f"{voltage_from_u(u):.16e}", f"{w:.16e}"))


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
