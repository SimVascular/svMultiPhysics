# Aliev--Panfilov ionic reference

`aliev_panfilov_1996.py` is adapted from Python generated with libCellML 0.6.3
from Physiome's [Aliev--Panfilov 1996 CellML model][cellml]. The model is in
exposure [`afd`][exposure], from workspace `a1` at pinned
[revision `1c3a018574af68610e2b95973f32fa831ea3096f`][revision]. The CellML
source has SHA-256
`995f3315fc4dbb5173ee03d7e7a10e8b6c3b033e2ec6e21a6f9e8a49a1c47ba0`.
Its metadata names David Nickerson as the CellML document creator and cites
R. R. Aliev and A. V. Panfilov as the model authors.

The AP artifact has no model-specific license page. The
[Physiome repository-wide policy][citation] licenses publicly accessible
content under CC BY 3.0 Unported ([terms][license]); attribution and adaptation
notices are retained in `aliev_panfilov_1996.py`.

## Generation, adaptation, and protocol

To reproduce the intermediate code, parse the CellML 1.1 source with
non-strict `Parser(False)`. Add the global unit `pms_per_uApmmsq`, defined as
`pms * uApmmsq^-1`. Change `units="dimensionless"` to `units="uApmmsq"` on
`interface/v`, `ionic_current/v`, and `recovery_variable/v`. Change
`units="pms"` to `units="pms_per_uApmmsq"` on `interface/d` and
`recovery_variable/d`. Register `interface/Istim` as the
`AnalyserExternalVariable` and analyze the processing model with libCellML
0.6.3. These are CellML interface and unit-annotation repairs; no MathML
numerical equation is changed.

Generate `Generator.implementationCode()` with
`GeneratorProfile.Profile.PYTHON` and the default libCellML 0.6.3 Python-profile
settings. Write the returned source unchanged as UTF-8 with LF line endings
and its existing final newline. The generated Python has SHA-256
`2023bad501b460aaa0ada70795d2659976c48b3ba7f9a0d98b765d3bd835a498`.

The derived core applies the [Goktepe--Kuhl][goktepe-kuhl] formulation used by
svMultiPhysics. Physiome's common threshold is split into `alpha=0.01` in the
voltage equation and `b=0.15` in the recovery equation. The voltage convention
is shifted to `V=100*u-80 mV`; `Cm=0.129` gives
`Cm*(Vp-Vr)=12.9`, encoding the 12.9 ms time scale. Recovery rates are
`0.002/12.9` and `0.2/12.9`, and the current mapping is
`I_CellML=-Cm*I_svMP`. The generated recovery state maps numerically to svMP
`w`.

The thin driver starts from `(V,w)=(-80 mV,0.001)`, applies
`Istim=-35.714 pA/pF` for `10 <= t < 12 ms`, and uses simultaneous Forward
Euler with `dt=0.1 ms` for 6000 updates, with zero SAC. Stimulus is evaluated
at old-state time. Output contains only the `step,V_mV,w` header and selected
checkpoint rows.

```bash
python3 generate_aliev_panfilov.py \
  --output ionic_aliev_panfilov_stimulated_trajectory.csv
```

[exposure]: https://models.physiomeproject.org/e/afd
[cellml]: https://models.physiomeproject.org/workspace/a1/rawfile/1c3a018574af68610e2b95973f32fa831ea3096f/models/1996_aliev/model.xml
[revision]: https://models.physiomeproject.org/workspace/a1/@@file/1c3a018574af68610e2b95973f32fa831ea3096f/models/1996_aliev/
[citation]: https://models.cellml.org/e/49f/experiments/cell-model.xml/license_citation
[license]: https://creativecommons.org/licenses/by/3.0/legalcode
[goktepe-kuhl]: https://doi.org/10.1002/nme.2571
