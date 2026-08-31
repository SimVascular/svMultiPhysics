#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause

"""Generate independent Bueno-Orovio EPI/ENDO/M forward-Euler oracles."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence, TextIO


@dataclass(frozen=True)
class Parameters:
    """One phenotype parameter set for the final 2008 equations."""

    u_o: float
    u_u: float
    theta_v: float
    theta_w: float
    theta_v_minus: float
    theta_o: float
    tau_v1_minus: float
    tau_v2_minus: float
    tau_v_plus: float
    tau_w1_minus: float
    tau_w2_minus: float
    k_w_minus: float
    u_w_minus: float
    tau_w_plus: float
    tau_fi: float
    tau_o1: float
    tau_o2: float
    tau_so1: float
    tau_so2: float
    k_so: float
    u_so: float
    tau_s1: float
    tau_s2: float
    k_s: float
    u_s: float
    tau_si: float
    tau_w_inf: float
    w_inf_star: float


@dataclass(frozen=True)
class Profile:
    zone_id: int
    parameters: Parameters
    update_count: int
    checkpoints: tuple[int, ...]


EPI_PARAMETERS = Parameters(
    u_o=0.0, u_u=1.55, theta_v=0.30, theta_w=0.13,
    theta_v_minus=0.006, theta_o=0.006,
    tau_v1_minus=60.0, tau_v2_minus=1150.0, tau_v_plus=1.4506,
    tau_w1_minus=60.0, tau_w2_minus=15.0, k_w_minus=65.0,
    u_w_minus=0.03, tau_w_plus=200.0, tau_fi=0.11,
    tau_o1=400.0, tau_o2=6.0, tau_so1=30.0181, tau_so2=0.9957,
    k_so=2.0458, u_so=0.65, tau_s1=2.7342, tau_s2=16.0,
    k_s=2.0994, u_s=0.9087, tau_si=1.8875, tau_w_inf=0.07,
    w_inf_star=0.94,
)

ENDO_PARAMETERS = Parameters(
    u_o=0.0, u_u=1.56, theta_v=0.30, theta_w=0.13,
    theta_v_minus=0.20, theta_o=0.006,
    tau_v1_minus=75.0, tau_v2_minus=10.0, tau_v_plus=1.4506,
    tau_w1_minus=6.0, tau_w2_minus=140.0, k_w_minus=200.0,
    u_w_minus=0.016, tau_w_plus=280.0, tau_fi=0.10,
    tau_o1=470.0, tau_o2=6.0, tau_so1=40.0, tau_so2=1.2,
    k_so=2.0, u_so=0.65, tau_s1=2.7342, tau_s2=2.0,
    k_s=2.0994, u_s=0.9087, tau_si=2.9013, tau_w_inf=0.0273,
    w_inf_star=0.78,
)

# Intentional regression value: the published M-cell tau_s2 is 4 ms.
M_PARAMETERS = Parameters(
    u_o=0.0, u_u=1.61, theta_v=0.30, theta_w=0.13,
    theta_v_minus=0.10, theta_o=0.005,
    tau_v1_minus=80.0, tau_v2_minus=1.4506, tau_v_plus=1.4506,
    tau_w1_minus=70.0, tau_w2_minus=8.0, k_w_minus=200.0,
    u_w_minus=0.016, tau_w_plus=280.0, tau_fi=0.078,
    tau_o1=410.0, tau_o2=7.0, tau_so1=91.0, tau_so2=0.8,
    k_so=2.1, u_so=0.6, tau_s1=2.7342, tau_s2=2.0,
    k_s=2.0994, u_s=0.9087, tau_si=3.3849, tau_w_inf=0.01,
    w_inf_star=0.5,
)

V_OFFSET_MV = -84.0
V_SCALE_MV = 85.7
INITIAL_PUBLIC_STATE = (-84.0, 1.0, 1.0, 0.0)
STIMULUS_AMPLITUDE_PUBLIC = -35.714
STIMULUS_START_MS = 10.0
STIMULUS_END_MS = 12.0
ZERO_STIMULUS = 0.0
ZERO_SAC = 0.0
DT_MS = 0.01
PROFILES = {
    "epi": Profile(
        zone_id=1,
        parameters=EPI_PARAMETERS,
        update_count=60_000,
        checkpoints=(
            0, 1000, 1001, 1109, 1198, 1200, 2083, 6018, 8589,
            10000, 20000, 23059, 27518, 28294, 28403, 30247,
            40000, 60000,
        ),
    ),
    "endo": Profile(
        zone_id=2,
        parameters=ENDO_PARAMETERS,
        update_count=120_000,
        checkpoints=(
            0, 1000, 1001, 1107, 1200, 1972, 2500, 4944,
            10000, 20000, 24067, 28075, 28557, 28755, 28876,
            30720, 40000, 60000, 120000,
        ),
    ),
    "m": Profile(
        zone_id=3,
        parameters=M_PARAMETERS,
        update_count=120_000,
        checkpoints=(
            0, 1000, 1001, 1106, 1200, 2021, 2500, 4651,
            10000, 20000, 30000, 36369, 44724, 46147, 46456,
            46640, 48735, 60000, 120000,
        ),
    ),
}

State = tuple[float, float, float, float]
Row = tuple[int, float, float, float, float]


def heaviside(x: float) -> float:
    """Paper convention: H(x)=0 for x<0 and H(x)=1 otherwise."""
    return 0.0 if x < 0.0 else 1.0


def voltage_to_u(voltage_mv: float) -> float:
    return (voltage_mv - V_OFFSET_MV) / V_SCALE_MV


def u_to_voltage(u: float) -> float:
    return V_OFFSET_MV + V_SCALE_MV * u


def stimulus_public_at_time(time_ms: float) -> float:
    """Return the half-open public stimulus pulse at an old-state time."""
    if STIMULUS_START_MS <= time_ms < STIMULUS_END_MS:
        return STIMULUS_AMPLITUDE_PUBLIC
    return ZERO_STIMULUS


def stimulus_internal_at_time(time_ms: float) -> float:
    """Scale the public pA/pF stimulus as IonicModel::integ does."""
    return stimulus_public_at_time(time_ms) / V_SCALE_MV


def rhs(
    state: State,
    parameters: Parameters,
    i_stim: float = ZERO_STIMULUS,
    i_sac: float = ZERO_SAC,
) -> State:
    """Return (du/dt,dv/dt,dw/dt,ds/dt), with time measured in ms."""
    u, v, w, s = state
    p = parameters
    h_v = heaviside(u - p.theta_v)
    h_w = heaviside(u - p.theta_w)
    h_v_minus = heaviside(u - p.theta_v_minus)
    h_o = heaviside(u - p.theta_o)

    tau_v_minus = (
        (1.0 - h_v_minus) * p.tau_v1_minus
        + h_v_minus * p.tau_v2_minus
    )
    tau_w_minus = p.tau_w1_minus + 0.5 * (
        p.tau_w2_minus - p.tau_w1_minus
    ) * (
        1.0 + math.tanh(p.k_w_minus * (u - p.u_w_minus))
    )
    tau_so = p.tau_so1 + 0.5 * (p.tau_so2 - p.tau_so1) * (
        1.0 + math.tanh(p.k_so * (u - p.u_so))
    )
    tau_s = (1.0 - h_w) * p.tau_s1 + h_w * p.tau_s2
    tau_o = (1.0 - h_o) * p.tau_o1 + h_o * p.tau_o2
    v_inf = 1.0 - h_v_minus
    w_inf = (1.0 - h_o) * (1.0 - u / p.tau_w_inf) + h_o * p.w_inf_star

    i_fi = -v * h_v * (u - p.theta_v) * (p.u_u - u) / p.tau_fi
    i_so = (u - p.u_o) * (1.0 - h_w) / tau_o + h_w / tau_so
    i_si = -h_w * w * s / p.tau_si

    du_dt = -(i_fi + i_so + i_si + i_stim) + i_sac
    dv_dt = (
        (1.0 - h_v) * (v_inf - v) / tau_v_minus
        - h_v * v / p.tau_v_plus
    )
    dw_dt = (
        (1.0 - h_w) * (w_inf - w) / tau_w_minus
        - h_w * w / p.tau_w_plus
    )
    s_inf = 0.5 * (1.0 + math.tanh(p.k_s * (u - p.u_s)))
    ds_dt = (s_inf - s) / tau_s
    return du_dt, dv_dt, dw_dt, ds_dt


def forward_euler_public(
    state: State,
    profile: Profile,
    time_ms: float,
    dt_ms: float = DT_MS,
) -> State:
    """One literal simultaneous FE update using the public voltage state."""
    voltage_mv, v, w, s = state
    internal_state = (voltage_to_u(voltage_mv), v, w, s)
    du_dt, dv_dt, dw_dt, ds_dt = rhs(
        internal_state,
        profile.parameters,
        i_stim=stimulus_internal_at_time(time_ms),
    )
    next_u = internal_state[0] + dt_ms * du_dt
    return (
        u_to_voltage(next_u),
        v + dt_ms * dv_dt,
        w + dt_ms * dw_dt,
        s + dt_ms * ds_dt,
    )


def generate_full_trajectory(profile: Profile) -> list[State]:
    trajectory = [INITIAL_PUBLIC_STATE]
    state = INITIAL_PUBLIC_STATE
    for completed_updates in range(profile.update_count):
        old_state_time_ms = completed_updates * DT_MS
        state = forward_euler_public(state, profile, old_state_time_ms)
        trajectory.append(state)
    return trajectory


def selected_rows(trajectory: Sequence[State], profile: Profile) -> list[Row]:
    return [(step, *trajectory[step]) for step in profile.checkpoints]


def validate_canonical(
    trajectory: Sequence[State], rows: Sequence[Row], profile: Profile
) -> None:
    if len(trajectory) != profile.update_count + 1:
        raise AssertionError("trajectory has the wrong number of completed updates")
    if len(rows) != len(profile.checkpoints):
        raise AssertionError("wrong CSV row count")
    if tuple(row[0] for row in rows) != profile.checkpoints:
        raise AssertionError("CSV checkpoints are missing or out of order")

    for step, state in enumerate(trajectory):
        if not all(math.isfinite(value) for value in state):
            raise AssertionError(f"non-finite state at step {step}")


def write_rows(rows: Iterable[Row], destination: TextIO) -> None:
    writer = csv.writer(destination, lineterminator="\n")
    writer.writerow(("step", "V_mV", "v", "w", "s"))
    # csv.writer uses Python's shortest round-trip representation for floats.
    writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--profile",
        choices=tuple(PROFILES),
        default="epi",
        help="phenotype to generate (default: epi)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="write canonical CSV to this path instead of stdout",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    profile = PROFILES[args.profile]
    trajectory = generate_full_trajectory(profile)
    rows = selected_rows(trajectory, profile)
    validate_canonical(trajectory, rows, profile)
    if args.output is None:
        write_rows(rows, sys.stdout)
    else:
        with args.output.open("w", encoding="utf-8", newline="") as destination:
            write_rows(rows, destination)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
