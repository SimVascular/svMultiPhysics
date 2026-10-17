# **Problem Description**

Simulate the block compression problem with a prestress field read from a VTK file.

This case uses the mesh, material model, boundary conditions, and load from
`block_compression`. The additional `Prestress_file_path` entry points to a VTK
file containing a nodal `Stress` field. The test checks that the prestress read
from that file contributes to the reported stress outputs.

## Prestress Field

The prestress input is stored in `prestress.vtu` as a constant six-component
nodal `Stress` field. The only nonzero component is `Sxx = -48000`.

![Prestress Sxx](prestress_sxx.png)

Constant prestress `Sxx` component.

## Stress Output Comparison

The images below show the reported `Stress` `Sxx` component for the original
block compression case and for this case with the prestress field included.

![Stress Sxx without prestress](result_without_prestress_sxx.png)

Stress `Sxx` without prestress.

![Stress Sxx with prestress](result_with_prestress_sxx.png)

Stress `Sxx` with prestress.

## References

Liu, Ju, and Alison L. Marsden. A Unified Continuum and Variational Multiscale Formulation for Fluids, Solids, and Fluid Structure Interaction. *Computer Methods in Applied Mechanics and Engineering* 337 (August 2018): 549 97. https://doi.org/10.1016/j.cma.2018.03.045.
