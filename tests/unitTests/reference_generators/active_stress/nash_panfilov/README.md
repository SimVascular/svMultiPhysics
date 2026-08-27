# Nash--Panfilov active-stress reference

The generator implements the Nash--Panfilov active-tension equation in the
smooth Goktepe--Kuhl form. The unit-test model replaces transmembrane potential
with intracellular calcium; this adaptation is specific to svMultiPhysics and
is not part of either cited formulation.

The canonical protocol uses the slab-calibration parameters, a normalized
double-exponential calcium transient, `dt=1 ms`, and 200 Forward-Euler updates.
ActiveStress checkpoint label `N` is the state after outer update `N`, so label
0 follows the first update.

```bash
python3 generate_nash_panfilov.py \
  --output active_stress_nash_panfilov_twitch.csv
```
