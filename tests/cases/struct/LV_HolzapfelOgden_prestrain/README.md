# Prestraining the LV

Two-stage example of `<Prestrain>` on the `LV_HolzapfelOgden_passive`
geometry at an endocardial pressure of 10 mmHg (13332.2 dyne/cm²). Both stages
read the mesh and fibers in `mesh/`.

The imaged geometry is taken to be subject to the load. The prestrain is
a deformation gradient `F0` at every Gauss point, from an unknown stress-free
configuration to the mesh, and a displacement `u` from the mesh gives
`F = (I + Grad u) F0`. The coordinates never change; instead every Newton
iteration adds the displacement it reached to the prestrain,
`F0 <- (I + Grad u) F0`, and the next iteration starts again from rest. This is the incremental scheme of Gee,
Förster & Wall (Int. J. Numer. Meth. Biomed. Engng. 26, 2010, section 3) with
the increments composed on the imaged mesh, so the state converges to one that
carries the load on the imaged geometry without deforming, and rigid motions
of the deformed body leave the stress unchanged.

| Folder | Input | Reference |
|---|---|---|
| `1-prestrain/` | `solver.xml`: prestrain run, 3 steps | `result_003.vtu` |
| `2-forward/` | `solver.xml`: forward run, 3 steps, reads `prestrain.vtu` | `result_003.vtu` |

`2-forward/prestrain.vtu` is the converged prestrain: the result of the
`1-prestrain` input run for 20 steps. The regression test
(`test_LV_HolzapfelOgden_prestrain` in `tests/test_struct.py`) runs each stage
from its folder and compares its last result with the reference.

```
cd 1-prestrain && mpirun -np 1 <build>/bin/svmultiphysics solver.xml   # writes 1-procs/result_NNN.vtu
cd 2-forward   && mpirun -np 1 <build>/bin/svmultiphysics solver.xml   # reads prestrain.vtu
```

To regenerate `2-forward/prestrain.vtu`, set `Number_of_time_steps` to 20 in
`1-prestrain/solver.xml`, run it, and copy `1-procs/result_020.vtu`.

`1-prestrain` applies the full load from the first step as a dead load
(`Follower_pressure_load false`, as in the paper). Every Newton iteration
starts from rest and composes the displacement it reaches into the prestrain.
With `Adaptive_time_step` in the general simulation parameters the time step grows
as the state approaches equilibrium, so once it is large the residual falls
quadratically and reaches round-off by step 14 of a 20-step run. Every result
carries the prestrain as the cell arrays `Prestrain_F_g<g>`, one per Gauss
point.

`2-forward` loads the converged prestrain through `<Prestrain_file_path>`
and holds the same 10 mmHg as a follower load, so the geometry should not move:
its displacement stays at round-off. Raise the pressure in
