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
python3 initialize_0D_boundary_conditions.py path/to/solver_vtu_and_bcs.xml
```

Write changes into the JSON configuration file referenced by the XML:

```bash
python3 initialize_0D_boundary_conditions.py path/to/solver_vtu_and_bcs.xml --write
```

Options
- `--no-backup` : when used with `--write` do not create a `.bak` copy of the JSON file.
- `--flow-mode {abs,oriented}` : use absolute flow values (`abs`, default) or oriented
  (signed) flux (`oriented`).
- `--pressure-array` : VTU point-data array name for pressure (default: `Pressure`).
- `--velocity-array` : VTU point-data array name for velocity (default: `Velocity`).

## Example using `tests/cases/fluid/pipe_RCR_sv0D`

This directory includes a lightweight example XML that reuses the mesh, surface,
restart, and svZeroD JSON files from:

```bash
./tests/cases/fluid/pipe_RCR_sv0D
```

If those test assets are present as Git LFS pointer files in your checkout, fetch
the LFS contents before running the example.

Run it from the repository root in print-only mode:

```bash
python3 utilities/initialize_0D_boundary_conditions/initialize_0D_boundary_conditions.py \
  utilities/initialize_0D_boundary_conditions/pipe_RCR_sv0D_example/solver_initial_conditions.xml
```

The example XML adds the `Initial_pressures_file_path` and
`Initial_velocities_file_path` entries needed by this utility and points both to
`tests/cases/fluid/pipe_RCR_sv0D/result_002.vtu`. It also reuses the test case's
`lumen_outlet` coupled boundary and `RCR_coupling` svZeroD block.

To update the test case JSON in place, rerun with `--write`:

```bash
python3 utilities/initialize_0D_boundary_conditions/initialize_0D_boundary_conditions.py \
  utilities/initialize_0D_boundary_conditions/pipe_RCR_sv0D_example/solver_initial_conditions.xml \
  --write
```

By default, `--write` creates
`tests/cases/fluid/pipe_RCR_sv0D/svzerod_3Dcoupling.json.bak`.

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
python3 initialize_0D_boundary_conditions.py solver_vtu_and_bcs.xml --write
```

## Troubleshooting

- If the script errors about missing point arrays, confirm the array names with a
  VTU inspector and pass `--pressure-array`/`--velocity-array` accordingly.
- If the script says a VTU or VTP is a Git LFS pointer file, install/enable Git
  LFS and fetch the real test assets before running the example.
- If surfaces report zero area, inspect the VTP files for degenerate cells.
