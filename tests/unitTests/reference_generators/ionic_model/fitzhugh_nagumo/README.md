# FitzHugh--Nagumo ionic reference

`fitzhugh_1961.py` is adapted from the [CellML API-generated Python][python]
for Physiome's [FitzHugh (1961) model][cellml], *Impulses and Physiological
States in Theoretical Models of Nerve Membrane*. The exposure is
`cf32346a9e5c4b2cdb559b11da5f1ae1`, from workspace `fitzhugh_1961` at
[revision `afcf226d5ada1c3b419b515f19a18095fba28e01`][revision].
The generated source has SHA-256
`2cf433fe361e7c4a865c756c925b18b2222bd327d762eb05d5273a9348649c10`.

The source is [CC BY 3.0 Unported][citation] ([terms][license]);
[model attribution][metadata] is retained in `fitzhugh_1961.py`.

## Adaptation and protocol

The two generated rate assignments and `initConsts()` are unchanged.
`computeRates()` takes a caller-supplied stimulus instead of the built-in
pulse. Unused legends, algebraic-output and piecewise helpers, demonstration
integration/plotting, and their imports are removed; only Python's standard
library is needed.

The driver sets CellML `alpha=-0.5`, `gamma=-0.6`, and `epsilon=0.02` for the
tested svMP parameters `alpha=-0.5`, `a=0`, `b=-0.6`, and `c=50`.
It uses `T=50*t` (CellML milliseconds versus svMP model-time units),
`dT=50*dt`, and `I_CellML=Istim/50`; this mapping is only claimed for `a=0`.
CellML states `(v,w)` map directly to CSV columns `(u,w)`.

The trajectory starts at the exact unstable equilibrium `(u,w)=(0,0)`, applies
`Istim=0.5` for `0.10 <= t < 0.12`, and uses Forward Euler with `dt=0.0005`
for 3000 updates, with zero SAC. Stimulus is evaluated at old-state time
`n*dt`; the CellML-coordinate pulse is `+0.01` on `[5,6)` ms, and the FE step
is `0.025` ms. Both states advance simultaneously. The trajectory stops during
recovery from the first triggered cycle and before the next autonomous upstroke.

Output contains only the `step,u,w` header and the selected checkpoint rows.

```bash
python3 generate_fitzhugh_nagumo.py \
  --output ionic_fitzhugh_nagumo_fe_trajectory.csv
```

[python]: https://models.physiomeproject.org/exposure/cf32346a9e5c4b2cdb559b11da5f1ae1/fitzhugh_1961.cellml/@@cellml_codegen/Python/raw
[cellml]: https://models.physiomeproject.org/exposure/cf32346a9e5c4b2cdb559b11da5f1ae1/fitzhugh_1961.cellml/source_text
[revision]: https://models.physiomeproject.org/workspace/fitzhugh_1961/rawfile/afcf226d5ada1c3b419b515f19a18095fba28e01/fitzhugh_1961.cellml
[metadata]: https://models.physiomeproject.org/exposure/cf32346a9e5c4b2cdb559b11da5f1ae1/fitzhugh_1961.cellml/cmeta
[citation]: https://models.physiomeproject.org/exposure/cf32346a9e5c4b2cdb559b11da5f1ae1/fitzhugh_1961.cellml/license_citation
[license]: https://creativecommons.org/licenses/by/3.0/legalcode
