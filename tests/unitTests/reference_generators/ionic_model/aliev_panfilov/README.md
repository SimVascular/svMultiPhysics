# Aliev--Panfilov ionic reference

The generator implements the split Aliev--Panfilov formulation and the
voltage/time mapping documented by Goktepe and Kuhl.

The canonical trajectory starts from `(-80 mV, 0.001)`, applies
`Istim=-35.714 pA/pF` for `10 <= t < 12 ms`, and uses Forward Euler with
`dt=0.1 ms` for 6000 updates. `Ksac=0`. Checkpoint `N` is the state after
exactly `N` completed updates.

```bash
python3 generate_aliev_panfilov.py \
  --output ionic_aliev_panfilov_stimulated_trajectory.csv
```
