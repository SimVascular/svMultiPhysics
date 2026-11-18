
Solve for a solid beam subject to a transverse load at the tip.

The `Stress` reference uses the stored displacements and the material constants in
`solver.xml`. For each tetrahedron, compute the symmetric displacement gradient
and apply Hooke's law. At each node, average the resulting stress tensors from
the adjacent tetrahedra, weighted by their volumes.

