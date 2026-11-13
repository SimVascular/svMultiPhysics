# **Problem Description**

Inflation of a thick-walled cylinder discretized with quadratic tetrahedra
(TET10), solved with the displacement-based solid equation.

The case covers the stress post-processing on a mesh with curved element faces.
The mid-edge nodes of the elements on the `endo` and `epi` faces sit on the
exact cylindrical surfaces, so 336 of the 864 mid-edge nodes of this mesh lie
outside the straight-sided tetrahedron spanned by their four corner nodes, by
up to 13% of a shape function. `post::tensor_post` computes the stress-type outputs
at the corner nodes and interpolates them to the mid-edge nodes, so the case
enables every stress-type output. A linear mesh, or a quadratic mesh with
straight edges, leaves that interpolation untested.

## Geometry and model

A 90 degree sector of the wall of a cylinder, `Ri` = 10 mm, `Re` = 12 mm,
`L` = 4 mm, with 2 elements through the wall (144 elements, 325 nodes).

- `base`, `top`: `u_z` = 0 (fixed length, free to slide in the plane)
- `y0`, `x0`: symmetry planes, `u_y` = 0 and `u_x` = 0
- `endo`: follower pressure load of 2 kPa

The material is neo-Hookean with `mu` = 100 kPa and a Simo-Taylor 91
dilational penalty with `kappa` = 2 `mu`. The penalty is soft on purpose: the
mesh is coarse, and a nearly incompressible penalty makes the
displacement-based element lock, which is not what this case covers.

## Mesh

The mesh is structured: each cell of a (radius, angle, axis) grid is split into
6 tetrahedra, and the mid-edge nodes lie on the exact cylinder, so the `endo`
and `epi` faces are curved to second order.
