# Bueno--Orovio ionic references

The generator implements the final Bueno-Orovio, Cherry, and Fenton 2008
minimal ventricular model for epicardial (EPI), endocardial (ENDO), and M-cell
phenotypes. The older CellML transcription represents an earlier formulation
and is not used.

The public state is `(V_mV,v,w,s)`, with `V_mV=-84+85.7*u`. All profiles start
from `(-84 mV,1,1,0)`, apply `Istim=-35.714 pA/pF` for
`10 <= t < 12 ms`, and use simultaneous Forward Euler with `dt=0.01 ms` and
`Ksac=0`. EPI runs for 600 ms; ENDO and M run for 1200 ms. Checkpoint `N` is
the state after exactly `N` completed updates.

The published 2008 table gives M-cell `tau_s2=4 ms`; svMultiPhysics currently
uses `2 ms`. The canonical M reference intentionally uses 2 ms to test the
shipped parameter set.

Generate any phenotype with:

```bash
python3 generate_bueno_orovio.py --profile epi \
  --output ionic_bueno_orovio_epi_trajectory.csv
```

Valid profiles are `epi`, `endo`, and `m`.
