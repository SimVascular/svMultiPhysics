# pipe_RCR_sv0D_units

This test is identical to `fluid/pipe_RCR_sv0D` but uses mixed units:

- solver.xml: cgs units (dynes/cm^2, cm^3/s)
- svZeroD JSON: clinical units (mmHg, mL/s)

The `svZeroDSolver_interface` specifies conversion factors so that values are converted automatically:

- `Pressure_conversion_factor =  0.00075006157584566` (dynes/cm^2 -> mmHg)
- `Flowrate_conversion_factor = 1.0` (cm^3/s -> mL/s)

Mesh and inflow files are referenced from the original test directory `pipe_RCR_sv0D/` to avoid duplication.
Also, the reference result `result_002.vtu` in `pipe_RCR_sv0D/` is used for comparison
