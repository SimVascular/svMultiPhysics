# Prestraining the LV

Two-phase example of `<Prestrain>` on the `LV_HolzapfelOgden_passive`
geometry at an endocardial pressure of 10 mmHg (13332.2 dyne/cm²).

The imaged geometry is taken to be the loaded configuration. The prestrain is
a deformation gradient `F0` at every Gauss point, from an unknown stress-free
configuration to the mesh, and a displacement `u` from the mesh gives
`F = (I + Grad u) F0`. The coordinates never change; instead every step adds
the displacement it reached to the prestrain, `F0 <- (I + Grad u) F0`, and the
next step starts again from rest. This is the incremental scheme of Gee,
Förster & Wall (Int. J. Numer. Meth. Biomed. Engng. 26, 2010, section 3) with
the increments composed on the imaged mesh, so the state converges to one that
carries the load on the imaged geometry without deforming, and rigid motions
of the deformed body leave the stress unchanged.

```
mpirun -np 1 <build>/bin/svmultiphysics prestrain.xml    # writes prestrain/result_NNN.vtu
mpirun -np 1 <build>/bin/svmultiphysics forward.xml      # reads prestrain/result_020.vtu
python check_hold.py                                     # equilibrium hold test
```

`prestrain.xml` applies the full load from the first step as a dead load
(`Follower_pressure_load false`, as in the paper). With `Pseudo_transient`
the run is a pseudo-transient continuation: every step starts from rest and
prints `Pseudo-transient: max update`, and `Max_iterations 1` makes each step
a single Newton iteration. Without them every step is an ordinary time step,
iterated to convergence. With
`Pseudo_transient_adaptive_dt` the pseudo time step grows as the state
approaches equilibrium. With one Newton iteration per step the tangent is that
of the composed update, so once the time step is large the residual falls
quadratically and reaches round-off by step 10 of the 20 here. The last
result carries the prestrain as
the cell arrays `Prestrain_F_g<g>`, one per Gauss point. A prestrain run
solves the struct (or ustruct) equation alone and does not accept
`<Line_search_experimental>`; a forward run that reads the prestrain may use
both.

`forward.xml` loads that prestrain through `<Prestrain_file_path>`
and holds the same 10 mmHg as a follower load. `check_hold.py` asserts that the
geometry does not move and that the reported stress reproduces the prestrain
stress. Raise the pressure in `forward.xml` to simulate the prestrained LV; run
it without the file path to see the naive inflation for contrast.

The forward run may be started on any number of ranks regardless of how many
the prestrain run used: the prestrain is written in the mesh's original
element order and partitioned on input like the fibre directions.
