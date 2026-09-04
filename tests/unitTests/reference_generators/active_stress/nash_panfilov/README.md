# Nash--Panfilov active-stress reference

`nash_panfilov_2004.py` is adapted from Python generated with libCellML 0.6.3
from the `active_tension` component of Physiome's
[Nash--Panfilov 2004 CellML model][cellml]. The exposure is
`d96a64b94d824692955e06ec878a2d09`; the pinned workspace
[revision is `40f16c45ea89e60ca78582b18c7e93d40752073a`][revision], and the
downloaded CellML has SHA-256
`a2ae1974d0b434d1ed01d87ecf5048d0aec508408dce6117a0b3bc58e19b2967`.
The source is [CC BY 3.0 Unported][citation] ([terms][license]). Attribution
from the [model metadata][metadata] and raw CellML is retained in
`nash_panfilov_2004.py`: Martyn Nash is the model author, and David Nickerson
is also named as a document creator.

## Extraction, adaptation, and protocol

The Physiome web generator reports the complete model as underconstrained
because `interface/Istim` is externally supplied. The processing model instead
isolates `active_tension` and supplies its input `active_tension/u` externally.

To reproduce the intermediate code, parse a CellML 1.0 processing copy with
non-strict `Parser(False)` and retain only the global `ms` and `kPa` units and
the `active_tension` component. Set `active_tension/e0` to
`initial_value="1.0"` and `active_tension/kTa` to `initial_value="47.9"`.
Remove `public_interface="in"` from `t`, `u`, `e0`, and `kTa`; remove both
`public_interface="out"` and `private_interface="out"` from `Ta`. Remove only
the root MathML attributes `cmeta:id="Ta_deriv_eqn"` and
`cmeta:id="e_calc_eqn"`; the inner MathML IDs `Ta_deriv` and `e_calc` remain
unchanged. Register `active_tension/u` as the `AnalyserExternalVariable` and
analyze with libCellML 0.6.3. The existing dimensional warning for `e` is
retained rather than repaired.

Generate `Generator.implementationCode()` with
`GeneratorProfile.Profile.PYTHON` and the default libCellML 0.6.3 Python-profile
settings. Write the returned source unchanged as UTF-8 with LF line endings
and its existing final newline. The generated Python has SHA-256
`2cef2e3ebed9a4d8b1fbcac1c4a48de21d86bdfed24ea6b36a6fc5a57b849d32`.

The adapted core replaces normalized-voltage activation with prescribed
intracellular calcium, the original piecewise rate with the
[Göktepe--Kuhl smooth rate][goktepe-kuhl], and the target tension with the
svMultiPhysics calcium-dependent target. Its six parameter values match those
configured by the test. The thin driver owns the normalized double-exponential
calcium transient, initial tension, Forward-Euler updates, checkpoints, and CSV
output.

The canonical protocol uses `dt=1 ms` for 200 updates. ActiveStress checkpoint
label `N` is the state after outer update `N`, so label 0 follows the first
update. Output contains only the `step,Ta` header and checkpoint rows.

```bash
python3 generate_nash_panfilov.py \
  --output active_stress_nash_panfilov_twitch.csv
```

[cellml]: https://models.physiomeproject.org/exposure/d96a64b94d824692955e06ec878a2d09/nash_panfilov_2004.cellml/source_text
[revision]: https://models.physiomeproject.org/workspace/nash_panfilov_2004/rawfile/40f16c45ea89e60ca78582b18c7e93d40752073a/nash_panfilov_2004.cellml
[metadata]: https://models.physiomeproject.org/exposure/d96a64b94d824692955e06ec878a2d09/nash_panfilov_2004.cellml/cmeta
[citation]: https://models.physiomeproject.org/exposure/d96a64b94d824692955e06ec878a2d09/nash_panfilov_2004.cellml/license_citation
[license]: https://creativecommons.org/licenses/by/3.0/legalcode
[goktepe-kuhl]: https://doi.org/10.1007/s00466-009-0434-z
