# **Problem Description**

Inflation of a thick-walled cylinder discretized with quadratic tetrahedra
(TET10), solved with the displacement-based solid equation.

The purpose of this case is to exercise the **stress post-processing on a mesh
with curved element faces**. The mid-edge nodes of the elements on the
`endo` and `epi` faces sit on the exact cylindrical surfaces, so they lie
*outside* the straight-sided tetrahedron spanned by their four corner nodes
(336 of the 864 mid-edge nodes of this mesh do, by up to 13% of a shape
function). `post::tpost` computes the stress-type outputs at the corner nodes
and interpolates them to the mid-edge nodes with the lower-order (TET4) basis;
it used to recover the position of a mid-edge node by inverse-mapping its
physical coordinates into that straight-sided tetrahedron, which fails on a
curved mesh and aborted the run with

```
Error in computing shape functions
```

The case therefore enables every stress-type output. A linear mesh would not
catch the problem, and neither would a quadratic mesh with straight edges.

## Geometry and model

A 90 degree sector of the wall of a cylinder, `Ri` = 10 mm, `Re` = 12 mm,
`L` = 4 mm, with 2 elements through the wall (144 elements, 325 nodes).

- `base`, `top`: `u_z` = 0 (fixed length, free to slide in the plane)
- `y0`, `x0`: symmetry planes, `u_y` = 0 and `u_x` = 0
- `endo`: follower pressure load of 2 kPa

The material is neo-Hookean with `mu` = 100 kPa and a Simo-Taylor 91
dilational penalty with `kappa` = 2 `mu`. The penalty is deliberately soft:
the mesh is coarse and a nearly incompressible penalty makes the
displacement-based element lock, which has nothing to do with what this case
tests.

## Mesh

The mesh is structured: each cell of a (radius, angle, axis) grid is split into
6 tetrahedra, and the mid-edge nodes are placed on the exact cylinder, so the
`endo` and `epi` faces are curved to second order.
