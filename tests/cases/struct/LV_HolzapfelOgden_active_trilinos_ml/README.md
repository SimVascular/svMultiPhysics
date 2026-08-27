
# **Problem Description**

Left-ventricle contraction with a Holzapfel-Ogden material and active stress, solved with the
Trilinos linear algebra interface and the MueLu algebraic multigrid preconditioner
(`trilinos-ml`).

This is the `trilinos-ml` variant of [`LV_HolzapfelOgden_active`](../LV_HolzapfelOgden_active).
The physics, mesh, loading and reference solution are identical; only the `<LS>` block differs.

## Why this case exists

MueLu is told how many degrees of freedom each node carries through its `number of equations`
parameter, and uses it to amalgamate the matrix node-by-node before aggregating. Every other test
using MueLu is `fluid` or `FSI`, where that value is 4. This is the only MueLu case with `struct`
physics, where it is 3.

With a wrong value, MueLu's `CoalesceDropFactory` walks off the end of the row arrays during
aggregation and the run aborts. The failure needs a mesh large enough that MueLu actually builds a
multigrid hierarchy — smaller `struct` cases such as `block_compression` (193 nodes) converge in a
single Krylov iteration and never reach the aggregation path, so they do not exercise this at all.
This case (4577 nodes) aborts at 2 or more MPI ranks when `number of equations` is wrong, and
passes at 1 rank either way.

The reference `result_001.vtu` is the same file as in `LV_HolzapfelOgden_active`, so this case also
checks that the Trilinos/MueLu path reproduces the FSILS path to within the standard field
tolerances.

## References

Holzapfel, Gerhard A., and Ray W. Ogden. Constitutive Modelling of Passive Myocardium: A
Structurally Based Framework for Material Characterization. *Philosophical Transactions of the
Royal Society A* 367, no. 1902 (2009): 3445–75. https://doi.org/10.1098/rsta.2009.0091.
