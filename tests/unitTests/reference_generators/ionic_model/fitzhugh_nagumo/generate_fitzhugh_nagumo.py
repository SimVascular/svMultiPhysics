#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause

"""Generate the canonical stimulated FitzHugh--Nagumo reference."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys


ALPHA = -0.50
A = 0.0
B = -0.60
C = 50.0

State = tuple[float, float]
INITIAL_STATE: State = (0.0, 0.0)
STIMULUS_AMPLITUDE = 0.5
STIMULUS_START = 0.10
STIMULUS_END = 0.12
DT = 5.0e-4
UPDATE_COUNT = 3000
END_TIME = DT * UPDATE_COUNT

CHECKPOINT_STEPS = (
    0,
    200,
    201,
    240,
    409,
    551,
    680,
    1000,
    1255,
    1379,
    1502,
    2000,
    2800,
    3000,
)


def stimulus_at_time(time: float) -> float:
    return STIMULUS_AMPLITUDE if STIMULUS_START <= time < STIMULUS_END else 0.0


def intrinsic_rhs(state: State) -> State:
    u, w = state
    return (
        C * (u * (u - ALPHA) * (1.0 - u) - w),
        u - B * w + A,
    )


def stimulated_rhs(time: float, state: State) -> State:
    du, dw = intrinsic_rhs(state)
    return du + stimulus_at_time(time), dw


def forward_euler(dt: float, end_time: float) -> list[State]:
    step_count_float = end_time / dt
    step_count = int(round(step_count_float))
    if not math.isclose(step_count_float, step_count, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("end_time must be an integer multiple of dt")

    states = [INITIAL_STATE]
    for step in range(step_count):
        old_state = states[step]
        du, dw = stimulated_rhs(step * dt, old_state)
        new_state = old_state[0] + dt * du, old_state[1] + dt * dw
        if not all(math.isfinite(value) for value in new_state):
            raise RuntimeError(
                f"Forward-Euler trajectory contains non-finite values at step {step + 1}"
            )
        states.append(new_state)
    return states


def csv_text(states: list[State]) -> str:
    lines = ["step,u,w"]
    for step in CHECKPOINT_STEPS:
        u, w = states[step]
        lines.append(f"{step},{u:.16e},{w:.16e}")
    return "\n".join(lines) + "\n"


def write_csv(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write the canonical CSV")
    args = parser.parse_args()

    states = forward_euler(DT, END_TIME)
    text = csv_text(states)
    if args.output:
        write_csv(args.output, text)
    else:
        sys.stdout.write(text)


if __name__ == "__main__":
    main()
