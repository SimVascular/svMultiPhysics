# Example 0D Coupled Boundary Conditions Example

This directory is a small, portable example of the files needed to connect
svMultiPhysics coupled outlet boundary conditions to svZeroDSolver RCR blocks.

It intentionally does not include the 3D mesh, face surfaces, or VTU restart
data. The paths in `solver.xml` show the expected layout, but the referenced
files are placeholders.

## Files

- `solver.xml`: minimal svMultiPhysics configuration showing coupled outlet BCs.
- `svzerod_3Dcoupling.json`: matching svZeroD coupling configuration.
- `inlet.flow`: tiny example inlet waveform referenced by `solver.xml`.

## Expected External Data

To run the initial-condition update script against this example, provide:

- `mesh-complete/mesh-complete.mesh.vtu`
- `mesh-complete/mesh-surfaces/*.vtp` for each face listed in `solver.xml`
- `num-procs/result_6000.vtu` with point-data arrays named `Pressure` and `Velocity`
- a valid path to `libsvzero_interface.so`

Then run from this directory or adjust paths accordingly:

```bash
python3 ../../update_coupled_json_initial_conditions.py solver.xml --write
```

The important naming contract is:

1. Each coupled XML `Add_BC name` must match an XML `Add_face name`.
2. Each XML `svZeroDSolver_block` must match a JSON
   `external_solver_coupling_blocks[].name`.
3. Each JSON coupling block `connected_block` must match an RCR
   `boundary_conditions[].bc_name`.

