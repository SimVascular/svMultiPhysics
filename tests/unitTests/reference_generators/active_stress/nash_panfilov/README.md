# Nash--Panfilov active-stress reference

The generator follows the [Nash and Panfilov (2004)](https://doi.org/10.1016/j.pbiomolbio.2004.01.016)
active-tension model as reformulated by
[Göktepe and Kuhl (2009)](https://doi.org/10.1007/s00466-009-0434-z).
svMultiPhysics uses intracellular calcium in place of the transmembrane-potential
activation variable.

The canonical protocol uses the slab-calibration parameters, a normalized
double-exponential calcium transient, `dt=1 ms`, and 200 Forward-Euler updates.
ActiveStress checkpoint label `N` is the state after outer update `N`, so label
0 follows the first update.

```bash
python3 generate_nash_panfilov.py \
  --output active_stress_nash_panfilov_twitch.csv
```
