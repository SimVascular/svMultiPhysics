
# **Problem Description**

Simulate cardiac electromechanics on a slab of myocardial tissue. This
directory contains two solver configurations that share the same geometry and
electrophysiology setup but differ in the active-stress model:

| Configuration file           | Active-stress model |
|------------------------------|---------------------|
| `solver_NashPanfilov.xml`    | Nash-Panfilov       |
| `solver_Regazzoni.xml`       | RDQ20-MF (Regazzoni)|

Both configurations couple cardiac electrophysiology (`CEP`) to solid mechanics
(`struct`), reproducing the geometry and stimulation setting of the Niederer
electrophysiology benchmark [1] with the addition of active contraction and
finite-strain mechanics.

## Shared Geometry and Electrophysiology

The mesh is a rectangular slab (`mesh/`) with two boundary faces `X0` and `X1`.

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

The tissue is modeled as a nearly incompressible Holzapfel-Ogden material with
modified anisotropy (`HolzapfelOgden-ModifiedAnisotropy`) [4]. The slab is fixed
with a zero-displacement Dirichlet boundary condition on the `X1` face, and
contracts as the depolarization wave propagates through the tissue.

## Nash-Panfilov variant (`solver_NashPanfilov.xml`)

Active contraction is driven by the calcium concentration computed by the
electrophysiology model, through the Nash-Panfilov active-stress model [5] with a
directional distribution along the fiber, sheet, and sheet-normal directions.

```
<Active_stress>
  <Model>NashPanfilov</Model>
  <Directional_distribution>
    <Fiber_direction> 0.7 </Fiber_direction>
    <Sheet_direction> 0.2 </Sheet_direction>
    <Sheet_normal_direction> 0.1 </Sheet_normal_direction>
  </Directional_distribution>
  ...
</Active_stress>
```

**Regression reference:** `result_NashPanfilov_001.vtu`

## Regazzoni variant (`solver_Regazzoni.xml`)

Active contraction is driven by the calcium concentration computed by the
electrophysiology model, through the RDQ20-MF mean-field active-stress model [6],
configured with the published human body-temperature calibration expressed in the
solver's unit system (time in ms, calcium in mM, length in µm). The scalar active
tension is distributed along the fiber, sheet, and sheet-normal directions using
the same directional weights as the Nash-Panfilov variant.

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

The (0.7, 0.2, 0.1) fiber/sheet/sheet-normal directional weights are an
svMultiPhysics extension of the paper's fiber-only active stress formulation; they
are not prescribed by the RDQ20-MF model itself.

**Regression reference:** `result_Regazzoni_001.vtu`

### Validation

svMultiPhysics stores `T_act = a_XB * (μ_P^1 + μ_N^1) * φ(SL)` — the scalar
RDQ20-MF active tension — and distributes it into per-direction fields
(`Active_tension_fibers`, `Active_tension_sheets`, `Active_tension_normal`) using
the directional weights `η`. Each per-direction field stores `η · T_act`; their
sum recovers `T_act` because the directional weights sum to one. Assembly of
this scalar into the continuum active stress tensor follows the existing
svMultiPhysics mechanics convention. The formulation of that assembly will be
addressed separately.

The active tension fields in `result_Regazzoni_001.vtu` were validated
node-by-node against the C++ reference implementation at commit
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

[5] M. P. Nash and A. V. Panfilov. Electromechanical model of excitable tissue to
study reentrant cardiac arrhythmias. Progress in Biophysics and Molecular Biology,
85(2-3):501–522, 2004.

[6] F. Regazzoni, L. Dede', and A. Quarteroni. Biophysically detailed mathematical
models of multiscale cardiac active mechanics. PLOS Computational Biology,
16(10):e1008294, 2020.
