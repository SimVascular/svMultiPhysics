#!/usr/bin/env python3
"""
Calculate a common phase-aware RCR initial pressure from:

  1) One *.flow file in the current directory
     - the filename can be anything, e.g.:
           inflow.flow
           random.flow
           shifted.flow
     - first line is ignored
     - subsequent columns: time, stored flow
     - stored flow = - physical inflow

  2) solver.xml
     - all boundary conditions with Time_dependence == "RCR" are detected
     - Rp, Rd, C, Pd are read automatically
     - the number and names of RCR outlets may differ between models

Method
------
1. Estimate each outlet flow fraction:

       f_i = [1 / (Rp_i + Rd_i + Rg_i)]
             --------------------------------
             sum_j [1 / (Rp_j + Rd_j + Rg_j)]

   By default, Rg_i = 0.

2. Estimate outlet flow:

       Q_i(t) = f_i Q_in(t)

3. Calculate the periodic RCR capacitor state at the waveform start:

                    integral_0^T exp[-(T-s)/(Rd_i C_i)] Q_i(s) ds
       Pc_i*(0) = Pd_i + ------------------------------------------------
                              C_i [1 - exp(-T/(Rd_i C_i))]

4. Use the arithmetic mean of all outlet Pc_i*(0) values as one common
   RCR Initial_pressure and update every RCR outlet in solver.xml.

Default behavior
----------------
Running:

    python calculate_rcr_initial_pressure.py

expects exactly:

    - one *.flow file
    - solver.xml

in the current directory.

The *.flow filename can be anything.

Examples:

    inflow.flow
    random.flow
    shifted_001.flow

If exactly one *.flow file exists, it is automatically detected.

The script:

  - creates a backup: solver.xml.before_rcr_pc_update.bak
  - updates solver.xml in place
  - writes rcr_initial_pressure_summary.csv

Optional arguments
------------------
    --flow FILE
        Explicitly specify a flow file.
        Useful if multiple *.flow files exist in the directory.

    --solver FILE
        Explicitly specify the solver XML file.

    --output FILE
        Write an updated XML to another file instead of replacing
        the original solver.xml.

    --summary FILE
        Specify the CSV summary filename.

    --dry-run
        Calculate and print only; do not write XML.
"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np


MMHG_IN_CGS = 1333.22

# Optional effective 3D geometric resistance by outlet face.
# Leave empty to use exactly Rp + Rd, as in the proposed method.
#
# Example:
# GEOMETRIC_RESISTANCE_BY_FACE = {
#     "outflow": 120.0,
# }
GEOMETRIC_RESISTANCE_BY_FACE: Dict[str, float] = {}


@dataclass
class RCROutlet:
    name: str
    rp: float
    rd: float
    capacitance: float
    distal_pressure: float
    old_initial_pressure: float
    initial_pressure_element: ET.Element
    geometric_resistance: float = 0.0

    @property
    def total_resistance_for_split(self) -> float:
        return self.rp + self.rd + self.geometric_resistance

    @property
    def tau(self) -> float:
        return self.rd * self.capacitance


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate one common phase-aware RCR Initial_pressure from "
            "one *.flow file and solver.xml, then update the RCR entries."
        )
    )

    parser.add_argument(
        "--flow",
        type=Path,
        default=None,
        help=(
            "Optional input flow file. When omitted, the script "
            "automatically detects the only *.flow file in the "
            "current directory."
        ),
    )

    parser.add_argument(
        "--solver",
        type=Path,
        default=Path("solver.xml"),
        help="Input solver XML. Default: solver.xml",
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "Optional output XML path. When omitted, the input solver.xml "
            "is backed up and updated in place."
        ),
    )

    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("rcr_initial_pressure_summary.csv"),
        help=(
            "Output CSV summary. "
            "Default: rcr_initial_pressure_summary.csv"
        ),
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Calculate and print results without writing an XML file.",
    )

    return parser.parse_args()


def find_flow_file(flow_argument: Path | None) -> Path:
    """
    Determine the input .flow file.

    If --flow is explicitly provided:
        use that file.

    Otherwise:
        search the current directory for *.flow.

    Exactly one *.flow file must exist for automatic detection.
    """

    # ---------------------------------------------------------
    # Case 1: user explicitly supplied --flow FILE
    # ---------------------------------------------------------
    if flow_argument is not None:

        if not flow_argument.is_file():
            raise FileNotFoundError(
                f"Specified flow file not found: {flow_argument}"
            )

        return flow_argument

    # ---------------------------------------------------------
    # Case 2: automatically search current directory
    # ---------------------------------------------------------
    flow_files = sorted(Path(".").glob("*.flow"))

    if len(flow_files) == 0:
        raise FileNotFoundError(
            "No .flow file was found in the current directory.\n"
            "Expected exactly one *.flow file together with solver.xml."
        )

    if len(flow_files) > 1:
        filenames = "\n".join(
            f"  - {path.name}" for path in flow_files
        )

        raise ValueError(
            "Multiple .flow files were found in the current directory:\n"
            f"{filenames}\n\n"
            "Exactly one .flow file is expected for automatic detection.\n"
            "Remove the extra .flow files or explicitly select one using:\n\n"
            "    python calculate_rcr_initial_pressure.py "
            "--flow FILE.flow"
        )

    return flow_files[0]


def read_inflow_file(
    path: Path,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Read the detected *.flow file.

    The first line is ignored.

    Each remaining non-empty line must contain:

        time  stored_flow

    The physical inflow used in the calculation is:

        Q_in = -stored_flow
    """

    if not path.is_file():
        raise FileNotFoundError(
            f"Flow file not found: {path}"
        )

    times: List[float] = []
    stored_flows: List[float] = []

    with path.open("r") as stream:

        # The first line is intentionally ignored.
        first_line = stream.readline()

        if first_line == "":
            raise ValueError(
                f"Flow file is empty: {path}"
            )

        for line_number, line in enumerate(
            stream,
            start=2,
        ):
            stripped = line.strip()

            if not stripped:
                continue

            tokens = stripped.split()

            if len(tokens) < 2:
                raise ValueError(
                    f"{path}, line {line_number}: "
                    "expected at least two columns."
                )

            try:
                time_value = float(tokens[0])
                stored_flow_value = float(tokens[1])

            except ValueError as exc:
                raise ValueError(
                    f"{path}, line {line_number}: "
                    "non-numeric time or flow."
                ) from exc

            times.append(time_value)
            stored_flows.append(stored_flow_value)

    if len(times) < 2:
        raise ValueError(
            "The flow file must contain at least two data rows."
        )

    time = np.asarray(
        times,
        dtype=float,
    )

    stored_flow = np.asarray(
        stored_flows,
        dtype=float,
    )

    if (
        not np.all(np.isfinite(time))
        or not np.all(np.isfinite(stored_flow))
    ):
        raise ValueError(
            "The flow file contains NaN or infinite values."
        )

    if not np.all(np.diff(time) > 0.0):
        raise ValueError(
            "Flow-file time values must be strictly increasing."
        )

    # Normalize the cycle time to begin at zero.
    time = time - time[0]

    # Stored flow follows the SimVascular sign convention.
    physical_inflow = -stored_flow

    if np.mean(physical_inflow) <= 0.0:
        print(
            "WARNING: the mean physical inflow (-stored flow) "
            "is not positive. "
            "Confirm the sign convention in the .flow file.",
            file=sys.stderr,
        )

    # A periodic input normally repeats its first point
    # at the final time.
    scale = max(
        1.0,
        float(np.max(np.abs(physical_inflow))),
    )

    endpoint_mismatch = (
        abs(
            physical_inflow[-1]
            - physical_inflow[0]
        )
        / scale
    )

    if endpoint_mismatch > 1.0e-3:
        print(
            "WARNING: first and last physical-flow values differ "
            "by more than 0.1% of the waveform scale. "
            "The code will still use the supplied interval "
            "as one cardiac cycle.",
            file=sys.stderr,
        )

    return time, physical_inflow


def require_float(
    parent: ET.Element,
    tag: str,
    outlet_name: str,
) -> float:

    element = parent.find(tag)

    if (
        element is None
        or element.text is None
    ):
        raise ValueError(
            f"RCR outlet '{outlet_name}' "
            f"is missing <{tag}> in solver.xml."
        )

    try:
        value = float(
            element.text.strip()
        )

    except ValueError as exc:
        raise ValueError(
            f"RCR outlet '{outlet_name}' "
            f"has a non-numeric <{tag}> value."
        ) from exc

    if not math.isfinite(value):
        raise ValueError(
            f"RCR outlet '{outlet_name}' "
            f"has a non-finite <{tag}> value."
        )

    return value


def read_solver_rcr_outlets(
    solver_path: Path,
) -> Tuple[
    ET.ElementTree,
    List[RCROutlet],
]:
    """
    Find every <Add_BC> whose
    <Time_dependence> text is RCR.
    """

    if not solver_path.is_file():
        raise FileNotFoundError(
            f"Solver XML not found: {solver_path}"
        )

    try:
        tree = ET.parse(
            solver_path
        )

    except ET.ParseError as exc:
        raise ValueError(
            f"Could not parse XML file: {solver_path}"
        ) from exc

    root = tree.getroot()

    outlets: List[RCROutlet] = []

    for boundary_condition in root.findall(
        ".//Add_BC"
    ):

        time_dependence = (
            boundary_condition.findtext(
                "Time_dependence"
            )
        )

        if (
            time_dependence is None
            or time_dependence.strip().upper()
            != "RCR"
        ):
            continue

        outlet_name = boundary_condition.get(
            "name",
            "<unnamed RCR outlet>",
        )

        values = boundary_condition.find(
            "RCR_values"
        )

        if values is None:
            raise ValueError(
                f"RCR outlet '{outlet_name}' "
                "is missing <RCR_values>."
            )

        initial_pressure_element = values.find(
            "Initial_pressure"
        )

        if (
            initial_pressure_element is None
            or initial_pressure_element.text is None
        ):
            raise ValueError(
                f"RCR outlet '{outlet_name}' "
                "is missing <Initial_pressure>."
            )

        outlet = RCROutlet(

            name=outlet_name,

            rp=require_float(
                values,
                "Proximal_resistance",
                outlet_name,
            ),

            rd=require_float(
                values,
                "Distal_resistance",
                outlet_name,
            ),

            capacitance=require_float(
                values,
                "Capacitance",
                outlet_name,
            ),

            distal_pressure=require_float(
                values,
                "Distal_pressure",
                outlet_name,
            ),

            old_initial_pressure=float(
                initial_pressure_element.text.strip()
            ),

            initial_pressure_element=(
                initial_pressure_element
            ),

            geometric_resistance=float(
                GEOMETRIC_RESISTANCE_BY_FACE.get(
                    outlet_name,
                    0.0,
                )
            ),
        )

        if outlet.rp < 0.0:
            raise ValueError(
                f"RCR outlet '{outlet_name}' "
                "has negative Rp."
            )

        if outlet.rd <= 0.0:
            raise ValueError(
                f"RCR outlet '{outlet_name}' "
                "must have Rd > 0."
            )

        if outlet.capacitance <= 0.0:
            raise ValueError(
                f"RCR outlet '{outlet_name}' "
                "must have C > 0."
            )

        if (
            outlet.total_resistance_for_split
            <= 0.0
        ):
            raise ValueError(
                f"RCR outlet '{outlet_name}' "
                "has non-positive Rp + Rd + Rg."
            )

        outlets.append(
            outlet
        )

    if not outlets:
        raise ValueError(
            "No RCR boundary conditions "
            "were found in solver.xml."
        )

    return tree, outlets


def calculate_flow_fractions(
    outlets: Sequence[RCROutlet],
) -> np.ndarray:
    """
    Resistance-based flow split:

        f_i = [1 / (Rp_i + Rd_i + Rg_i)]
              / sum_j [1 / (Rp_j + Rd_j + Rg_j)]
    """

    total_resistances = np.asarray(
        [
            outlet.total_resistance_for_split
            for outlet in outlets
        ],
        dtype=float,
    )

    conductances = (
        1.0
        / total_resistances
    )

    fractions = (
        conductances
        / np.sum(conductances)
    )

    if not np.isclose(
        np.sum(fractions),
        1.0,
        atol=1.0e-12,
    ):
        raise RuntimeError(
            "Calculated flow fractions "
            "do not sum to one."
        )

    return fractions


def calculate_periodic_pc0(
    time: np.ndarray,
    outlet_flow: np.ndarray,
    outlet: RCROutlet,
) -> float:
    """
    Calculate the periodic capacitor state Pc_i*(0):

                         integral_0^T exp[-(T-s)/tau_i] Q_i(s) ds
        Pc_i*(0) = Pd_i + -----------------------------------------
                              C_i [1 - exp(-T/tau_i)]

        tau_i = Rd_i C_i
    """

    if time.shape != outlet_flow.shape:
        raise ValueError(
            "Time and outlet-flow arrays "
            "have different lengths."
        )

    cycle_period = float(
        time[-1]
    )

    if cycle_period <= 0.0:
        raise ValueError(
            "The cardiac-cycle period "
            "must be positive."
        )

    tau = outlet.tau

    exponential_weight = np.exp(
        -(cycle_period - time)
        / tau
    )

    integrand = (
        exponential_weight
        * outlet_flow
    )

    if hasattr(
        np,
        "trapezoid",
    ):
        weighted_flow_integral = float(
            np.trapezoid(
                integrand,
                time,
            )
        )

    else:
        weighted_flow_integral = float(
            np.trapz(
                integrand,
                time,
            )
        )

    denominator = (
        outlet.capacitance
        * (
            1.0
            - math.exp(
                -cycle_period
                / tau
            )
        )
    )

    if (
        denominator <= 0.0
        or not math.isfinite(
            denominator
        )
    ):
        raise ValueError(
            "Invalid periodic-state denominator "
            f"for outlet '{outlet.name}'."
        )

    pc0 = (
        outlet.distal_pressure
        + weighted_flow_integral
        / denominator
    )

    if not math.isfinite(
        pc0
    ):
        raise ValueError(
            "Calculated non-finite Pc*(0) "
            f"for outlet '{outlet.name}'."
        )

    return pc0


def write_summary_csv(
    path: Path,
    outlets: Sequence[RCROutlet],
    fractions: np.ndarray,
    outlet_pc0: np.ndarray,
    common_pc0: float,
    cycle_period: float,
) -> None:

    fieldnames = [
        "face",
        "Rp",
        "Rd",
        "C",
        "Pd",
        "Rg_optional",
        "Rp_plus_Rd_plus_Rg",
        "flow_fraction",
        "tau_RdC_sec",
        "tau_over_T",
        "old_Initial_pressure_cgs",
        "old_Initial_pressure_mmHg",
        "predicted_outlet_Pc0_cgs",
        "predicted_outlet_Pc0_mmHg",
        "common_Pc0_cgs",
        "common_Pc0_mmHg",
    ]

    with path.open(
        "w",
        newline="",
    ) as stream:

        writer = csv.DictWriter(
            stream,
            fieldnames=fieldnames,
        )

        writer.writeheader()

        for (
            outlet,
            fraction,
            pc0,
        ) in zip(
            outlets,
            fractions,
            outlet_pc0,
        ):

            writer.writerow(
                {
                    "face":
                        outlet.name,

                    "Rp":
                        outlet.rp,

                    "Rd":
                        outlet.rd,

                    "C":
                        outlet.capacitance,

                    "Pd":
                        outlet.distal_pressure,

                    "Rg_optional":
                        outlet.geometric_resistance,

                    "Rp_plus_Rd_plus_Rg":
                        outlet.total_resistance_for_split,

                    "flow_fraction":
                        fraction,

                    "tau_RdC_sec":
                        outlet.tau,

                    "tau_over_T":
                        outlet.tau
                        / cycle_period,

                    "old_Initial_pressure_cgs":
                        outlet.old_initial_pressure,

                    "old_Initial_pressure_mmHg":
                        outlet.old_initial_pressure
                        / MMHG_IN_CGS,

                    "predicted_outlet_Pc0_cgs":
                        pc0,

                    "predicted_outlet_Pc0_mmHg":
                        pc0
                        / MMHG_IN_CGS,

                    "common_Pc0_cgs":
                        common_pc0,

                    "common_Pc0_mmHg":
                        common_pc0
                        / MMHG_IN_CGS,
                }
            )


def update_solver_xml(
    tree: ET.ElementTree,
    outlets: Sequence[RCROutlet],
    common_pc0: float,
    solver_path: Path,
    output_path: Path | None,
) -> Path:
    """
    When output_path is None:

      - back up solver.xml
      - replace solver.xml in place

    Otherwise:

      - leave the input solver.xml unchanged
      - write to output_path
    """

    for outlet in outlets:

        outlet.initial_pressure_element.text = (
            f"{common_pc0:.10f}"
        )

    # Improve readability when supported
    # by the Python version.
    if hasattr(
        ET,
        "indent",
    ):
        ET.indent(
            tree,
            space="    ",
        )

    if output_path is None:

        backup_path = (
            solver_path.with_name(
                solver_path.name
                + ".before_rcr_pc_update.bak"
            )
        )

        if backup_path.exists():

            counter = 1

            while True:

                candidate = (
                    solver_path.with_name(
                        solver_path.name
                        + (
                            ".before_rcr_pc_update_"
                            f"{counter}.bak"
                        )
                    )
                )

                if not candidate.exists():
                    backup_path = candidate
                    break

                counter += 1

        shutil.copy2(
            solver_path,
            backup_path,
        )

        destination = (
            solver_path
        )

        print(
            f"\nBackup created: {backup_path}"
        )

    else:
        destination = (
            output_path
        )

    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    tree.write(
        destination,
        encoding="utf-8",
        xml_declaration=True,
    )

    return destination


def main() -> int:

    args = parse_arguments()

    # ---------------------------------------------------------
    # Automatically determine the *.flow input file.
    # ---------------------------------------------------------
    flow_path = find_flow_file(
        args.flow
    )

    # ---------------------------------------------------------
    # Read input files.
    # ---------------------------------------------------------
    time, physical_inflow = (
        read_inflow_file(
            flow_path
        )
    )

    tree, outlets = (
        read_solver_rcr_outlets(
            args.solver
        )
    )

    cycle_period = float(
        time[-1]
    )

    fractions = (
        calculate_flow_fractions(
            outlets
        )
    )

    outlet_pc0_values: List[float] = []

    for (
        outlet,
        fraction,
    ) in zip(
        outlets,
        fractions,
    ):

        estimated_outlet_flow = (
            fraction
            * physical_inflow
        )

        pc0 = calculate_periodic_pc0(
            time=time,
            outlet_flow=estimated_outlet_flow,
            outlet=outlet,
        )

        outlet_pc0_values.append(
            pc0
        )

    outlet_pc0 = np.asarray(
        outlet_pc0_values,
        dtype=float,
    )

    # ---------------------------------------------------------
    # One identical Initial_pressure
    # for every RCR outlet.
    # ---------------------------------------------------------
    common_pc0 = float(
        np.mean(
            outlet_pc0
        )
    )

    # ---------------------------------------------------------
    # Warn when distal pressures differ.
    #
    # The simple conductance-based flow split assumes
    # approximately the same upstream-to-distal
    # pressure difference.
    # ---------------------------------------------------------
    distal_pressures = np.asarray(
        [
            outlet.distal_pressure
            for outlet in outlets
        ],
        dtype=float,
    )

    if (
        np.ptp(
            distal_pressures
        )
        > 1.0e-10
    ):
        print(
            "\nWARNING: distal pressures differ among outlets. "
            "The requested resistance-only flow-fraction "
            "approximation may be less accurate in this case.",
            file=sys.stderr,
        )

    # ---------------------------------------------------------
    # Calculate mean inflow.
    # ---------------------------------------------------------
    if hasattr(
        np,
        "trapezoid",
    ):

        mean_flow = (
            float(
                np.trapezoid(
                    physical_inflow,
                    time,
                )
            )
            / cycle_period
        )

    else:

        mean_flow = (
            float(
                np.trapz(
                    physical_inflow,
                    time,
                )
            )
            / cycle_period
        )

    # ---------------------------------------------------------
    # Print summary.
    # ---------------------------------------------------------
    print()
    print(
        "============================================================"
    )
    print(
        "Phase-aware common RCR initial-pressure calculation"
    )
    print(
        "============================================================"
    )

    print(
        f"Flow file                 : {flow_path}"
    )

    print(
        f"Solver XML                : {args.solver}"
    )

    print(
        f"Cardiac-cycle period      : "
        f"{cycle_period:.8f} s"
    )

    print(
        f"Detected RCR outlet count : "
        f"{len(outlets)}"
    )

    print(
        f"Mean physical inflow      : "
        f"{mean_flow:.8f}"
    )

    print()

    header = (
        f"{'Outlet':<24}"
        f"{'f_i':>12}"
        f"{'Rp':>16}"
        f"{'Rd':>16}"
        f"{'C':>16}"
        f"{'Pc_i*(0) [mmHg]':>20}"
    )

    print(
        header
    )

    print(
        "-"
        * len(header)
    )

    for (
        outlet,
        fraction,
        pc0,
    ) in zip(
        outlets,
        fractions,
        outlet_pc0,
    ):

        print(
            f"{outlet.name:<24}"
            f"{fraction:>12.6f}"
            f"{outlet.rp:>16.6g}"
            f"{outlet.rd:>16.6g}"
            f"{outlet.capacitance:>16.6g}"
            f"{pc0 / MMHG_IN_CGS:>20.6f}"
        )

    print()

    print(
        "------------------------------------------------------------"
    )

    print(
        "Common Initial_pressure   : "
        f"{common_pc0:.10f} "
        "[CGS pressure]"
    )

    print(
        "Common Initial_pressure   : "
        f"{common_pc0 / MMHG_IN_CGS:.6f} "
        "mmHg"
    )

    print(
        "Rounded integer value     : "
        f"{int(round(common_pc0))}"
    )

    print(
        "------------------------------------------------------------"
    )

    # ---------------------------------------------------------
    # Write CSV summary.
    # ---------------------------------------------------------
    write_summary_csv(
        path=args.summary,
        outlets=outlets,
        fractions=fractions,
        outlet_pc0=outlet_pc0,
        common_pc0=common_pc0,
        cycle_period=cycle_period,
    )

    print(
        f"\nSummary CSV written: "
        f"{args.summary}"
    )

    # ---------------------------------------------------------
    # Dry-run mode.
    # ---------------------------------------------------------
    if args.dry_run:

        print(
            "Dry run: solver.xml was not changed."
        )

        return 0

    # ---------------------------------------------------------
    # Update solver.xml.
    # ---------------------------------------------------------
    updated_path = update_solver_xml(
        tree=tree,
        outlets=outlets,
        common_pc0=common_pc0,
        solver_path=args.solver,
        output_path=args.output,
    )

    print(
        f"Updated XML written: "
        f"{updated_path}"
    )

    print(
        f"All {len(outlets)} RCR "
        "<Initial_pressure> values were set "
        f"to {common_pc0:.10f}."
    )

    return 0


if __name__ == "__main__":

    try:
        raise SystemExit(
            main()
        )

    except Exception as error:

        print(
            f"\nERROR: {error}",
            file=sys.stderr,
        )

        raise SystemExit(1)