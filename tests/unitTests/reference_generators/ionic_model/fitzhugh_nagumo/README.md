# FitzHugh--Nagumo ionic reference

The generator implements the FitzHugh--Nagumo parameterization
`alpha=-0.5`, `a=0`, `b=-0.6`, and `c=50`. Provenance was cross-checked against
Physiome workspace `fitzhugh_1961`, revision
`8f99d77593097ffe84e6936425a0d8d2aaa8a880`, using
`t_CellML = 50 * t_svMP`; generated CellML code is not required.

The trajectory starts at the exact unstable equilibrium `(u,w)=(0,0)`, applies
`Istim=0.5` for `0.10 <= t < 0.12`, and uses Forward Euler with `dt=0.0005`
for 3000 updates. It stops during recovery from the first triggered cycle and
before the next autonomous upstroke.

```bash
python3 generate_fitzhugh_nagumo.py \
  --output ionic_fitzhugh_nagumo_fe_trajectory.csv
```
