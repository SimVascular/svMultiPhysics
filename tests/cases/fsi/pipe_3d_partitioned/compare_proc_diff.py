#!/usr/bin/env python3
"""Plot absolute difference of centerline pressure and velocity between
1-proc and 3-proc partitioned FSI results along z, for each coupling tolerance."""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

from compare_fsi import read_mesh, extract_centerline

CASES = [
    {
        "label":  "cp=1e-6 (baseline)",
        "dir_1":  "results_cp1e-6_1proc",
        "dir_3":  "results_cp1e-6_3proc",
    },
    {
        "label":  "cp=1e-7",
        "dir_1":  "results_cp1e-7_1proc",
        "dir_3":  "results_cp1e-7_3proc",
    },
    {
        "label":  "cp=1e-6 tight nltol",
        "dir_1":  "results_cp1e-6_tight_nltol_1proc",
        "dir_3":  "results_cp1e-6_tight_nltol_3proc",
    },
    {
        "label":  "cp=1e-6 tight lstol",
        "dir_1":  "results_cp1e-6_tight_lstol_1proc",
        "dir_3":  "results_cp1e-6_tight_lstol_3proc",
    },
]

STYLES = ["b-o", "r--s", "g-.^", "m:D"]


def load_centerline(result_dir, prefix, step, field):
    f = os.path.join(result_dir, f"{prefix}_{step:03d}.vtu")
    if not os.path.exists(f):
        raise FileNotFoundError(f)
    m = read_mesh(f)
    z, vals = extract_centerline(m, field)
    return z, vals


def abs_diff_on_common_z(z1, v1, z3, v3):
    """Interpolate both onto the union of z points and return |v1 - v3|."""
    z_common = np.union1d(np.round(z1, 10), np.round(z3, 10))
    v1_interp = np.interp(z_common, z1, v1)
    v3_interp = np.interp(z_common, z3, v3)
    return z_common, np.abs(v1_interp - v3_interp)


def main():
    parser = argparse.ArgumentParser(
        description="Absolute 1proc vs 3proc difference along centerline"
    )
    parser.add_argument("--step", type=int, default=3, help="Time step to compare")
    args = parser.parse_args()

    plt.rcParams.update({"font.size": 10, "figure.figsize": (8, 5)})
    step = args.step

    # ---- pressure difference ----
    fig, ax = plt.subplots()
    for ci, case in enumerate(CASES):
        z1, p1 = load_centerline(case["dir_1"], "result_fluid", step, "Pressure")
        z3, p3 = load_centerline(case["dir_3"], "result_fluid", step, "Pressure")
        if len(z1) == 0 or len(z3) == 0:
            print(f"WARNING: empty centerline pressure for {case['label']}, skipping")
            continue
        z, diff = abs_diff_on_common_z(z1, p1, z3, p3)
        ax.plot(z, diff, STYLES[ci], ms=3, label=case["label"])
    ax.set_xlabel("z")
    ax.set_ylabel("|pressure 1proc − 3proc|")
    ax.set_title(f"Centerline pressure difference (step {step})")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig("diff_pressure_z.pdf")
    print("Saved diff_pressure_z.pdf")

    # ---- axial velocity difference ----
    fig, ax = plt.subplots()
    for ci, case in enumerate(CASES):
        z1, v1 = load_centerline(case["dir_1"], "result_fluid", step, "Velocity")
        z3, v3 = load_centerline(case["dir_3"], "result_fluid", step, "Velocity")
        if len(z1) == 0 or len(z3) == 0:
            print(f"WARNING: empty centerline velocity for {case['label']}, skipping")
            continue
        v1_ax = v1[:, 2] if v1.ndim > 1 else v1
        v3_ax = v3[:, 2] if v3.ndim > 1 else v3
        z, diff = abs_diff_on_common_z(z1, v1_ax, z3, v3_ax)
        ax.plot(z, diff, STYLES[ci], ms=3, label=case["label"])
    ax.set_xlabel("z")
    ax.set_ylabel("|velocity 1proc − 3proc|")
    ax.set_title(f"Centerline axial velocity difference (step {step})")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig("diff_velocity_z.pdf")
    print("Saved diff_velocity_z.pdf")

    plt.close("all")
    print("Done.")


if __name__ == "__main__":
    main()
