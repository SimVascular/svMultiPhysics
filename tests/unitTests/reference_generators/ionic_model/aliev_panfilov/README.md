# Aliev--Panfilov ionic reference

The generator follows the split-parameter Aliev--Panfilov formulation used by
svMultiPhysics, following [Göktepe and Kuhl (2009)](https://doi.org/10.1002/nme.2571)
and based on the original [Aliev and Panfilov (1996)](https://doi.org/10.1016/0960-0779\(95\)00089-5)
model.

The canonical trajectory starts from `(-80 mV, 0.001)`, applies
`Istim=-35.714 pA/pF` for `10 <= t < 12 ms`, and uses Forward Euler with
`dt=0.1 ms` for 6000 updates. `Ksac=0`. Checkpoint `N` is the state after
exactly `N` completed updates.

```bash
python3 generate_aliev_panfilov.py \
  --output ionic_aliev_panfilov_stimulated_trajectory.csv
```
