# RCR Initial Pressure Utility

This utility calculates a phase-aware initial capacitor pressure for RCR boundary conditions and updates the `Initial_pressure` values of all RCR outlets in an svMultiPhysics `solver.xml` file.

## Requirements

- Python 3
- NumPy
- One svMultiPhysics `solver.xml` file
- One inflow waveform file with the `.flow` extension

## Usage

Place `solver.xml` and one `.flow` file in the simulation directory, then run:

```bash
python /path/to/calculate_rcr_initial_pressure.py
```

The `.flow` file is detected automatically when exactly one is present. A specific file can also be selected:

```bash
python /path/to/calculate_rcr_initial_pressure.py --flow inflow.flow
```

To calculate the pressure without modifying `solver.xml`, use:

```bash
python /path/to/calculate_rcr_initial_pressure.py --dry-run
```

## Output

By default, the utility:

- creates `solver.xml.before_rcr_pc_update.bak`;
- updates `solver.xml` in place; and
- writes `rcr_initial_pressure_summary.csv`.

Run the following command for all available options:

```bash
python /path/to/calculate_rcr_initial_pressure.py --help
```
