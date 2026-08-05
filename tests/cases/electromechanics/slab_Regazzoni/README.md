
# **Problem Description**

Simulate cardiac electromechanics on a slab of myocardial tissue. This test
couples cardiac electrophysiology (`CEP`) to solid mechanics (`struct`),
reproducing the geometry and stimulation setting of the Niederer electrophysiology
benchmark [1] with the addition of active contraction and finite-strain
mechanics.

## Electrophysiology

The propagation of the transmembrane potential is modeled with the
ten-Tusscher-Panfilov (`TTP`) cell activation model [2, 3], using epicardial
parameters (included from `../../cep/ttp_parameters/ttp_epicardium_parameters.xml`)
and anisotropic conductivity aligned with the fiber direction. The domain is split
into two `Domain`s: an unstimulated region (`domain 1`) and a stimulated region
(`domain 2`) where an external `Istim` stimulus initiates depolarization.

```
<Stimulus type="Istim" >
  <Amplitude> -35.714 </Amplitude>
  <Start_time> 0.0 </Start_time>
  <Duration> 2.0 </Duration>
  <Cycle_length> 10000.0 </Cycle_length>
</Stimulus>
```

## Mechanics

The tissue is modeled as a nearly incompressible Holzapfel-Ogden material with
modified anisotropy (`HolzapfelOgden-ModifiedAnisotropy`) [4]. Active contraction
is driven by the calcium concentration computed by the electrophysiology model,
through the RDQ20-MF mean-field active-stress model [5], configured with the
published human body-temperature calibration. The scalar active tension is
distributed along the fiber, sheet, and sheet-normal directions using the same
directional weights as the NashPanfilov slab case.

```
<Active_stress>
  <Model>Regazzoni</Model>
  <Directional_distribution>
    <Fiber_direction> 0.7 </Fiber_direction>
    <Sheet_direction> 0.2 </Sheet_direction>
    <Sheet_normal_direction> 0.1 </Sheet_normal_direction>
  </Directional_distribution>
  ...
</Active_stress>
```

The slab is fixed with a zero-displacement Dirichlet boundary condition on the
`X1` face, and contracts as the depolarization wave propagates through the tissue.

The (0.7, 0.2, 0.1) fiber/sheet/sheet-normal directional weights are an
svMultiPhysics extension of the paper's fiber-only active stress formulation; they
are not prescribed by the RDQ20-MF model itself.

## Regression Reference

svMultiPhysics stores `T_act = a_XB * (μ_P^1 + μ_N^1) * φ(SL)` — the scalar
RDQ20-MF active tension — and distributes it into per-direction fields
(`Active_tension_fibers`, `Active_tension_sheets`, `Active_tension_normal`) using
the directional weights `η`. Each per-direction field stores `η · T_act`; their
sum recovers `T_act` because the directional weights sum to one. Assembly of
this scalar into the continuum active stress tensor follows the existing
svMultiPhysics mechanics convention. The formulation of that assembly will be
addressed separately.

The active tension fields in `result_001.vtu` were validated node-by-node against
the C++ reference implementation at commit
[`26f05df`](https://github.com/FrancescoRegazzoni/cardiac-activation/commit/26f05df28891df7b3c69f16bb136cdced6b63c4d).
Both implementations use the same implicit-Euler XB scheme, so agreement is
to machine precision (~1e-16 relative error). The comparison evaluates `T_act`
from the svMultiPhysics output directly against the reference C++ active tension,
using the calcium and sarcomere-length inputs from this one-step test. The remaining
fields in the VTU serve as integrated svMultiPhysics regression references and were
not independently validated by the RDQ20-MF reference code.

## References

[1] S. A. Niederer, E. Kerfoot, A. P. Benson, et al. Verification of cardiac tissue
electrophysiology simulators using an N-version benchmark. Philosophical Transactions
of the Royal Society A, 369(1954):4331–4351, 2011.

[2] K. H. W. J. ten Tusscher, D. Noble, P. J. Noble, and A. V. Panfilov. A model for
human ventricular tissue. American Journal of Physiology-Heart and Circulatory
Physiology, 286(4):H1573–H1589, apr 2004.

[3] K. H. W. J. ten Tusscher and A. V. Panfilov. Alternans and spiral breakup in a
human ventricular tissue model. American Journal of Physiology-Heart and Circulatory
Physiology, 291(3):H1088–H1100, sep 2006.

[4] G. A. Holzapfel and R. W. Ogden. Constitutive modelling of passive myocardium: a
structurally based framework for material characterization. Philosophical Transactions
of the Royal Society A, 367(1902):3445–3475, 2009.

[5] F. Regazzoni, L. Dede', and A. Quarteroni. Biophysically detailed mathematical
models of multiscale cardiac active mechanics. PLOS Computational Biology,
16(10):e1008294, 2020.
