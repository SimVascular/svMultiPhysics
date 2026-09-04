# Initialize svZeroD boundary conditions

A Python script used to initialize svZeroD boundary conditions from the VTU file used by svMultiPhysics for velocity and pressure initial conditions. 

## Purpose

This helper reads an svMultiPhysics XML configuration, locates the VTU restart file
containing point pressures and velocities, integrates per-face average pressure and
surface flux for coupled boundary faces, and produces `initial_condition` entries
that match the svZeroDSolver coupling configuration in the JSON configuration file.

## Key behavior

- Parses the XML for `Add_mesh` and `svZeroDSolver_interface/Configuration_file`.
- Loads the VTU (unstructured grid) and the face surface files (VTP).
- Integrates pressure and oriented flux on each coupled face and computes:
  - `pressure_c:<rcr_block>` (pressure corrected by Rp*flow)
  - `flow:<coupling_block>:<rcr_block>`
  - `pressure:<coupling_block>:<rcr_block>`
- Prints a JSON block with the `initial_condition` mapping; optionally writes it
  into the svZeroDSolver JSON configuration file.

## Requirements

- Python 3
- NumPy
- One svMultiPhysics `solver.xml` file
- One coupled svZeroDSolver `svzerod_3Dcoupling.json` file
- One VTU restart file containing pressure and velocity point data
- VTP surface files for each coupled boundary face
- vtk with Python bindings (VTK Python package must provide `vtkXMLUnstructuredGridReader`,
  `vtkXMLPolyDataReader`, and `vtk.util.numpy_support.vtk_to_numpy`).

Install (example):

```bash
python3 -m pip install numpy vtk
```

## Usage

Basic (print-only):

```bash
python3 update_coupled_json_initial_conditions.py path/to/solver_vtu_and_bcs.xml
```

Write changes into the JSON configuration file referenced by the XML:

```bash
python3 update_coupled_json_initial_conditions.py path/to/solver_vtu_and_bcs.xml --write
```

Options
- `--no-backup` : when used with `--write` do not create a `.bak` copy of the JSON file.
- `--flow-mode {abs,oriented}` : use absolute flow values (`abs`, default) or oriented
  (signed) flux (`oriented`).
- `--pressure-array` : VTU point-data array name for pressure (default: `Pressure`).
- `--velocity-array` : VTU point-data array name for velocity (default: `Velocity`).

## Compact example template

A small, portable configuration-only example is available in:

```./example_template/
```

It contains:

- `solver.xml`: svMultiPhysics setup with one inlet, one wall, and four coupled RCR outlets.
- `svzerod_3Dcoupling.json`: matching svZeroD coupling blocks and example
  `initial_condition` values generated from a longer simulation.
- `inlet.flow`: a tiny dummy inlet waveform.

The compact example intentionally does not include mesh, surface, or restart
files. The paths in `./example_template/solver.xml` show the expected layout:

```
mesh-complete/mesh-complete.mesh.vtu
mesh-complete/mesh-surfaces/*.vtp
num-procs/result_6000.vtu
```

To make the compact example runnable, add those files or update the paths in
`solver.xml`, then run:

```bash
python3 ../update_coupled_json_initial_conditions.py solver.xml --write
```

from inside `./example_template`.

## Notes & assumptions

- The script expects the VTU referenced by `Initial_pressures_file_path` and
  `Initial_velocities_file_path` to be the same file.
- Each coupled `Add_BC` in the XML must reference an `Add_face` (surface VTP) with a
  matching name.
- Surface files must be readable VTP files describing the face geometry.
- The JSON configuration must contain the expected `external_solver_coupling_blocks`
  and RCR boundary conditions with `Rp` values.
- Each XML `svZeroDSolver_block` name must match a JSON
  `external_solver_coupling_blocks[].name`.
- Each JSON coupling block `connected_block` must match an RCR
  `boundary_conditions[].bc_name`.

## Output

The script prints a summary of computed face pressure/flow/area and creates a JSON object of the form:

```json
{
  "initial_condition": { "...": ... }
}
```

When `--write` is used the script will update the JSON file in-place (and create
a `.bak` copy unless `--no-backup` is specified).

## Example

```bash
python3 update_coupled_json_initial_conditions.py solver_vtu_and_bcs.xml --write
```

## Troubleshooting

- If the script errors about missing point arrays, confirm the array names with a
  VTU inspector and pass `--pressure-array`/`--velocity-array` accordingly.
- If surfaces report zero area, inspect the VTP files for degenerate cells.
